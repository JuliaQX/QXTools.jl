include("tensor_cache.jl")

using ITensors
using YAML
using JLD2
using Random
using QXTn

const DSL_VERSION = VersionNumber("0.2")

export generate_dsl_files, generate_parameter_file

"""
    write_version_header(io::IO)

Function to write version head to DSL file with current version constant
"""
function write_version_header(io::IO)
    write(io, "# version: $(DSL_VERSION)\n")
end

"""
    write_dsl_load_header(tnc::TensorNetworkCircuit, dsl_io::IO, data_io::JLD2.JldFile)

Write the header part of the dsl file which loads tensors from their corresponding data symbols
and also prepares the data file
"""
function write_dsl_load_header(tnc::TensorNetworkCircuit, dsl_io::IO, data_io::JLD2.JLDFile)
    # create a tensor cache for saving intitial tensor data
    tc = TensorCache()

    write_version_header(dsl_io)

    # iterate over tensors and add to TensorCache
    write(dsl_io, "outputs $(qubits(tnc))\n")
    for (tensor_symbol, tensor) in pairs(tnc)
        if !(tensor_symbol in output_tensors(tnc))
            data_symbol = push!(tc, tensor_data(tensor, consider_hyperindices=true))
            write(dsl_io, "load $tensor_symbol $data_symbol\n")
        end
    end
    # save tensor data labels
    save_cache(tc, data_io)
    write(data_io, "output_0", [1., 0.])
    write(data_io, "output_1", [0., 1.])
end

"""
    generate_dsl_files(tnc::TensorNetworkCircuit,
                       prefix::String,
                       plan::Array{NTuple{3, Symbol}, 1},
                       sliced_bond_groups::Array{<:Array{<:Index, 1}, 1},
                       force::Bool=false)

Function to create a dsl and data files to contracting the given tensor network circuit
with the plan provided
"""
function generate_dsl_files(tnc::TensorNetworkCircuit,
                            prefix::String,
                            plan::Array{NTuple{3, Symbol}, 1},
                            sliced_bond_groups::Array{<:Array{<:Index, 1}, 1};
                            force::Bool=true)

    dsl_filename = "$(prefix).qx"
    data_filename = "$(prefix).jld2"
    param_filename = "$(prefix).yml"

    @assert force || !isfile(dsl_filename) "Error $(dsl_filename) already exists"
    @assert force || !isfile(data_filename) "Error $(data_filename) already exists"
    open(dsl_filename, "w") do dsl_io; jldopen(data_filename, "w") do data_io
        # write header of DSL file with data loading
        write_dsl_load_header(tnc, dsl_io, data_io)

        # we create a network with mocked tensors
        tnc_copy = copy(tnc)
        tn_copy = tnc_copy.tn

        # for the output tensors we replace them with output specific symbols
        # and create a symbol to symbol map to keep track of the mapping
        changed_ids = Dict{Symbol, Symbol}()
        for (i, tensor_sym) in enumerate(output_tensors(tnc_copy))
            new_sym = Symbol("\$o$i")
            replace_tensor_symbol!(tn_copy, tensor_sym, new_sym)
            changed_ids[tensor_sym] = new_sym
        end

        # create views on tensors
        for (i, slice_bond_group) in enumerate(sliced_bond_groups)
            related_tensors = union([tn_copy[slice_bond] for slice_bond in slice_bond_group]...)
            # slice_tensors = tn_copy[slice_bond]
            for tensor_sym in related_tensors
                new_sym = Symbol("$(tensor_sym)_\$v$i")
                tensor_bonds_in_group = intersect(inds(tn_copy[tensor_sym]), slice_bond_group)
                write_view_command(dsl_io, tn_copy, tensor_sym, new_sym, tensor_bonds_in_group[1], "\$v$(i)")
                replace_tensor_symbol!(tn_copy, tensor_sym, new_sym)
                changed_ids[tensor_sym] = new_sym
                write(dsl_io, "del $tensor_sym\n")
            end
        end

        # perform contraction
        for (A_sym, B_sym, C_sym) in plan
            while haskey(changed_ids, A_sym) A_sym = changed_ids[A_sym] end
            while haskey(changed_ids, B_sym) B_sym = changed_ids[B_sym] end
            write_ncon_command(dsl_io, tn_copy, A_sym, B_sym, C_sym)
            write(dsl_io, "del $A_sym\n")
            write(dsl_io, "del $B_sym\n")
            contract_pair!(tn_copy, A_sym, B_sym, C_sym; mock=true)
        end

        # if there is more than one tensor remaining we assume these are scalers and
        # contract in order which should amount to a product over the scalars
        if length(tn_copy) > 1
            tensors = collect(keys(tn_copy))
            A_sym = tensors[1]
            for j in 2:length(tensors)
                B_sym = tensors[j]
                C_sym = next_tensor_id!(tn_copy)
                write_ncon_command(dsl_io, tn_copy, A_sym, B_sym, C_sym)
                write(dsl_io, "del $A_sym\n")
                write(dsl_io, "del $B_sym\n")
                contract_pair!(tn_copy, A_sym, B_sym, C_sym; mock=true)
                A_sym = C_sym
            end
        end
        output_tensor = first(keys(tn_copy))
        write(dsl_io, "save $output_tensor output\n")
    end; end
    nothing
end

"""
    write_ncon_command(io::IO, tn::TensorNetwork, A_sym::Symbol, B_sym::Symbol, C_sym::Symbol)

Function that constructs and writes command describing pairwise contraction of tensors pointed to
by symbols A_sym and B_sym to a tensor pointed to by symbol C_sym to the given IO stream in the format

```
ncon C_sym C_labels A_sym A_labels B_sym B_labels
````

where the labels are comma separated lists integers which follow the Einstein summation convention. Where
one of the tensors is rank 0 (a scalar) a '0' is used as a placeholder for the labels to aid parsing.

"""
function write_ncon_command(io::IO, tn::TensorNetwork, A_sym::Symbol, B_sym::Symbol, C_sym::Symbol)
    A_indices = copy(QXTn.inds(tn[A_sym]))
    B_indices = copy(QXTn.inds(tn[B_sym]))
    common_indices = intersect(A_indices, B_indices)
    all_indices = union(A_indices, B_indices)
    C_indices = setdiff(all_indices, common_indices)

    # assign unique number to each index across all indices of both tensors
    all_index_map = Dict{Index, Int64}((y => x for (x, y) in enumerate(all_indices)))

    # update index map to map indices which belong to the same hyper index group
    # to the same number
    function update_all_index_map(sym)
        for group in hyperindices(tn[sym])
            @assert length(group) > 1 "Expect hyperedge groups to have at least 2 indices for $sym"
            if length(intersect(group, common_indices)) > 0
                ref_index = sort(intersect(group, common_indices), by=x->all_index_map[x])[1]
            else
                ref_index = group[1]
            end
            ref_pos = all_index_map[ref_index]
            for r in setdiff(group, ref_index)
                all_index_map[r] = ref_pos
            end
        end
    end

    update_all_index_map(A_sym)
    update_all_index_map(B_sym)

    C_labels = length(C_indices) == 0 ? "0" : join(unique(getindex.([all_index_map], C_indices)), ",")

    # must update A_indices (and B_indices) so only have a single index for each hyper edge group of A (B)
    function update_indices(indices, hyper_index_groups)
        indices = copy(indices)
        for group in hyper_index_groups
            ref_index = sort(group, by=x -> all_index_map[x])[1]
            for other in setdiff(group, ref_index)
                replace!(indices, other => ref_index)
            end
        end
        unique(indices)
    end

    A_labels = length(A_indices) == 0 ? "0" : join(getindex.([all_index_map], update_indices(A_indices, hyperindices(tn[A_sym]))), ",")
    B_labels = length(B_indices) == 0 ? "0" : join(getindex.([all_index_map], update_indices(B_indices, hyperindices(tn[B_sym]))), ",")
    write(io, "ncon $(C_sym) $(C_labels) $(A_sym) $(A_labels) $(B_sym) $(B_labels)\n")
    nothing
end

"""
    write_view_command(dsl_io::IO, tn::TensorNetwork, tensor_sym::Symbol, new_sym::Symbol, slice_bond::Index, bond_label::String)

Write command to create a view on an existing tensor to the dsl_io stream that is passed.
"""
function write_view_command(dsl_io::IO, tn::TensorNetwork, tensor_sym::Symbol, new_sym::Symbol, slice_bond::Index, bond_label::String)
    indices = QXTn.inds(tn[tensor_sym])
    index_map = OrderedDict{Index, Int64}((y => x for (x, y) in enumerate(indices)))

    for group in hyperindices(tn[tensor_sym])
        @assert length(group) >= 2 "Expect all hyperindex groups to have 2 or more elements"
        ref_index = index_map[group[1]]
        for r in group[2:end]
            index_map[r] = ref_index
        end
    end
    index_pos = Dict{Int64, Int64}()
    pos = 0
    for v in values(index_map)
        if !haskey(index_pos, v)
            pos += 1
            index_pos[v] = pos
        end
    end
    position_of_index = index_pos[index_map[slice_bond]]
    write(dsl_io, "view $new_sym $tensor_sym $position_of_index $bond_label\n")
end

"""
    generate_parameter_file(filename::String
                            sliced_bond_groups::Array{<:Array{<:Index, 1}, 1},
                            amplitudes::Union{Base.Generator, <: AbstractArray})

Generate a yml file with parameters corresponding to partitions and dimensions of each
along with the amplitudes to contract for.
"""
function generate_parameter_file(filename_prefix::String,
                                 sliced_bond_groups::Array{<:Array{<:Index, 1}, 1},
                                 amplitudes::Union{Base.Generator, <: AbstractArray})
    partition_dims = Dict{String, Int64}()
    for (i, sliced_bond_group) in enumerate(sliced_bond_groups)
        partition_dims["v$i"] = dim(sliced_bond_group[1])
    end
    partition_parameters = Dict("parameters" => partition_dims)
    config = Dict("partitions" => partition_parameters, "amplitudes" => collect(amplitudes))
    YAML.write_file("$(filename_prefix).yml", config)
end
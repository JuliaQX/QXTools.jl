include("tensor_cache.jl")
include("cmds.jl")
include("compute_tree.jl")
include("tree_opt.jl")
include("tree_stats.jl")

using YAML
using JLD2
using Random
using QXTns

const DSL_VERSION = VersionNumber("0.4")

export generate_dsl_files, generate_parameter_file

"""
    write_version_header(io::IO)

Function to write version head to DSL file with current version constant
"""
function write_version_header(io::IO)
    write(io, "# version: $(DSL_VERSION)\n")
end

"""
    write_metadata_header(io::IO, metadata::OrderedDict, indent::String="")

Function to write contraction plan metadata to DSL file.
"""
function write_metadata_header(io::IO, metadata::OrderedDict, indent::String="")
    for (k, v) in metadata
        if typeof(v) <: OrderedDict
            write(io, "# $(indent)$k :\n")
            write_metadata_header(io, v, indent * "  ")
        else
            write(io, "# $(indent)$k : $(v)\n")
        end
    end
end

"""
    write_dsl_load_header(tnc::TensorNetworkCircuit, dsl_io::IO, data_io::JLD2.JldFile)

Write the header part of the dsl file which loads tensors from their corresponding data symbols
and also prepares the data file
"""
function write_dsl_load_header(tnc::TensorNetworkCircuit, dsl_io::IO, data_io::JLD2.JLDFile;
                               metadata::OrderedDict=OrderedDict())
    # create a tensor cache for saving intitial tensor data
    tc = TensorCache()

    write_version_header(dsl_io)
    write_metadata_header(dsl_io, metadata)

    # iterate over tensors and add to TensorCache
    write(dsl_io, "outputs $(qubits(tnc))\n")
    for tensor_symbol in keys(tnc)
        if !(tensor_symbol in output_tensors(tnc))
            data_symbol = push!(tc, tensor_data(tnc, tensor_symbol))
            write(dsl_io, "load $tensor_symbol $data_symbol\n")
        end
    end
    # save tensor data labels
    save_cache(tc, data_io)
end

"""
    generate_dsl_files(tnc::TensorNetworkCircuit,
                       prefix::String,
                       plan::Array{NTuple{3, Symbol}, 1},
                       sliced_bond_groups::Array{<:Array{<:Index, 1}, 1};
                       force::Bool=true,
                       metadata::Dict{String, Any}=Dict{String, Any}())

Function to create a dsl and data files to contracting the given tensor network circuit
with the plan provided
"""
function generate_dsl_files(tnc::TensorNetworkCircuit,
                            prefix::String,
                            plan::Array{NTuple{3, Symbol}, 1},
                            sliced_bond_groups::Array{<:Array{<:Index, 1}, 1};
                            force::Bool=true,
                            metadata::OrderedDict=OrderedDict())

    dsl_filename = "$(prefix).qx"
    data_filename = "$(prefix).jld2"

    @assert force || !isfile(dsl_filename) "Error $(dsl_filename) already exists"
    @assert force || !isfile(data_filename) "Error $(data_filename) already exists"
    open(dsl_filename, "w") do dsl_io; jldopen(data_filename, "w") do data_io
        # write header of DSL file with data loading
        write_dsl_load_header(tnc, dsl_io, data_io; metadata=metadata)

        # we create a network with mocked tensors
        tnc_copy = copy(tnc)
        tn_copy = tnc_copy.tn

        # for the output tensors we replace them with output specific symbols
        # and create a symbol to symbol map to keep track of the mapping
        changed_ids = Dict{Symbol, Symbol}()
        for (i, tensor_sym) in enumerate(output_tensors(tnc_copy))
            new_sym = Symbol("o$i")
            replace_tensor_symbol!(tn_copy, tensor_sym, new_sym)
            changed_ids[tensor_sym] = new_sym
        end

        # create views on tensors
        for (i, slice_bond_group) in enumerate(sliced_bond_groups)
            related_tensors = union([tn_copy[slice_bond] for slice_bond in slice_bond_group]...)
            for tensor_sym in related_tensors
                new_sym = Symbol("$(tensor_sym)_s")
                write_view_command(dsl_io, tn_copy, tensor_sym, new_sym, slice_bond_group, "v$(i)")
                replace_tensor_symbol!(tn_copy, tensor_sym, new_sym)
                changed_ids[tensor_sym] = new_sym
            end
        end

        # perform contraction
        contract_cmds = ContractCommand[]
        for (A_sym, B_sym, C_sym) in plan
            while haskey(changed_ids, A_sym) A_sym = changed_ids[A_sym] end
            while haskey(changed_ids, B_sym) B_sym = changed_ids[B_sym] end
            push!(contract_cmds, gen_ncon_command(tn_copy, A_sym, B_sym, C_sym))
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
                push!(contract_cmds, gen_ncon_command(tn_copy, A_sym, B_sym, C_sym))
                contract_pair!(tn_copy, A_sym, B_sym, C_sym; mock=true)
                A_sym = C_sym
            end
        end
        contraction_tree = build_tree(contract_cmds)
        # remove_repeated!(contraction_tree)
        # permute_and_merge!(contraction_tree)
        # join_remaining!(contraction_tree)
        write(dsl_io, contraction_tree)

        output_tensor = first(keys(tn_copy))
        write(dsl_io, "save $output_tensor output\n")
    end; end
    nothing
end

"""
    gen_ncon_command(tn::TensorNetwork, A_sym::Symbol, B_sym::Symbol, C_sym::Symbol)

Function that constructs and writes command describing pairwise contraction of tensors pointed to
by symbols A_sym and B_sym to a tensor pointed to by symbol C_sym to the given IO stream in the format

```
ncon C_sym C_labels A_sym A_labels B_sym B_labels
````

where the labels are comma separated lists integers which follow the Einstein summation convention. Where
one of the tensors is rank 0 (a scalar) a '0' is used as a placeholder for the labels to aid parsing.

"""
function gen_ncon_command(tn::TensorNetwork, A_sym::Symbol, B_sym::Symbol, C_sym::Symbol)
    r = contraction_indices(tn, A_sym, B_sym)
    ContractCommand(C_sym, r.c_labels, A_sym, r.a_labels, B_sym, r.b_labels)
end

"""
    write_view_command(dsl_io::IO, tn::TensorNetwork, tensor_sym::Symbol, new_sym::Symbol, slice_bonds::Vector{<:Index}, bond_label::String)

Write command to create a view on an existing tensor to the dsl_io stream that is passed.
"""
function write_view_command(dsl_io::IO, tn::TensorNetwork, tensor_sym::Symbol, new_sym::Symbol, slice_bonds::Vector{<:Index}, bond_label::String)
    indices = hyperindices(tn, tensor_sym, all_indices=true)
    # find which group overlaps with slice_bonds
    index_position = findfirst(x -> length(intersect(slice_bonds, x)) > 0, indices)
    write(dsl_io, "view $new_sym $tensor_sym $index_position $bond_label\n")
end

"""
    generate_parameter_file(filename::String,
                            sliced_bond_groups::Array{<:Array{<:Index, 1}, 1},
                            amplitudes::Union{Base.Generator, <: AbstractArray})

Generate a yml file with parameters corresponding to partitions and dimensions of each
along with the amplitudes to contract for.

output:
    fix_M: false
    M: 0.001
    num_samples: 10
    output_method: rejection
    seed: ~
partitions:
    parameters:
        v1: 2
        v2: 2
"""
function generate_parameter_file(filename_prefix::String,
                                 sliced_bond_groups::Array{<:Array{<:Index, 1}, 1},
                                 output_parameters)
    partition_dims = OrderedDict{String, Int64}()
    for (i, sliced_bond_group) in enumerate(sliced_bond_groups)
        partition_dims["v$i"] = QXTns.dim(sliced_bond_group[1])
    end
    partition_parameters = Dict("parameters" => partition_dims)

    config = Dict()
    config["output"] = output_parameters
    config["partitions"] = partition_parameters

    YAML.write_file("$(filename_prefix).yml", config)
end
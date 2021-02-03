module DSL

include("tensor_cache.jl")

using ITensors
using YAML
using JLD
using Random
using ..TN

export generate_dsl_files, generate_parameter_file


"""
    write_dsl_load_header(tnc::TensorNetworkCircuit, dsl_filename::String, data_filename::String)

Write the header part of the dsl file which loads tensors from their corresponding data symbols
and also prepares the data file
"""
function write_dsl_load_header(tnc::TensorNetworkCircuit, dsl_filename::String, data_filename::String)
    # create a tensor cache for saving intitial tensor data
    tc = TensorCache()
    # iterate over tensors and add to TensorCache
    open(dsl_filename, "w") do io
        write(io, "outputs $(qubits(tnc))\n")
        for (tensor_symbol, tensor) in pairs(tnc)
            if !(tensor_symbol in output_tensors(tnc))
                data_symbol = push!(tc, tensor_data(tensor))
                write(io, "load $tensor_symbol $data_symbol\n")
            end
        end
    end
    # save tensor data labels
    save_cache(tc, data_filename)
    f = jldopen(data_filename, "r+")
    write(f, "output_0", [1., 0.])
    write(f, "output_1", [0., 1.])
    close(f)
end

"""
    generate_dsl_files(tnc::TensorNetworkCircuit,
                       prefix::String
                       plan::Array{<:Index, 1},
                       bonds::Array{<:Index, 1})

Function to create a dsl and data files to contracting the given tensor network circuit
with the plan provided
"""
function generate_dsl_files(tnc::TensorNetworkCircuit,
                            prefix::String,
                            plan::Array{<:Index, 1},
                            sliced_bonds::Array{<:Index, 1})

    dsl_filename = "$(prefix).tl"
    data_filename = "$(prefix).jld"
    param_filename = "$(prefix).yml"

    write_dsl_load_header(tnc, dsl_filename, data_filename)


    # wr create a network with mocked tensors
    tn_mock = use_mock_tensors(tnc.tn)

    for (i, tensor_sym) in enumerate(output_tensors(tnc))
        replace_tensor_symbol!(tn_mock, tensor_sym, Symbol("\$o$i"))
    end

    # create views on tensors
    # @assert all([!(x in output_indices(tnc)) for x in sliced_bonds]) "Cannot slice output indices"
    open(dsl_filename, "a") do io
        for (i, slice_bond) in enumerate(sliced_bonds)
            slice_tensors = tn_mock[slice_bond]
            for tensor_sym in slice_tensors
                position_of_index =  findfirst(x -> x == slice_bond, inds(tn_mock[tensor_sym]))
                new_sym = Symbol("$(tensor_sym)_\$v$i")
                replace_tensor_symbol!(tn_mock, tensor_sym, new_sym)
                write(io, "view $new_sym $tensor_sym $position_of_index \$v$(i)\n")
                write(io, "del $tensor_sym\n")
            end
        end
    end

    # perform contraction
    open(dsl_filename, "a") do io
        for index in plan
            if haskey(tn_mock, index)
                tensor_pair = tn_mock[index]
                if length(tensor_pair) == 2
                    A_sym = tensor_pair[1]
                    B_sym = tensor_pair[2]
                    ncon_labels = contract_ncon_indices(tn_mock, A_sym, B_sym)
                    C_sym = contract_pair!(tn_mock, A_sym, B_sym)
                    A_labels = join(ncon_labels[1], ",")
                    B_labels = join(ncon_labels[2], ",")
                    write(io, "ncon $C_sym $A_sym $A_labels $B_sym $B_labels\n")
                    write(io, "del $A_sym\n")
                    write(io, "del $B_sym\n")
                end
            end
        end
        if length(tn_mock) > 1
            tensors = collect(keys(tn_mock))
            A_sym = tensors[1]
            for j in 2:length(tensors)
                B_sym = tensors[j]
                ncon_labels = contract_ncon_indices(tn_mock, A_sym, B_sym)
                C_sym = contract_pair!(tn_mock, A_sym, B_sym)
                parse_labels = x -> length(x) == 0 ? "0" : join(x, ",")
                A_labels = parse_labels(ncon_labels[1])
                B_labels = parse_labels(ncon_labels[2])
                write(io, "ncon $C_sym $A_sym $A_labels $B_sym $B_labels\n")
                write(io, "del $A_sym\n")
                write(io, "del $B_sym\n")
                A_sym = C_sym
            end
        end
        output_tensor = first(keys(tn_mock))
        write(io, "save $output_tensor output\n")
    end
    nothing
end


"""
    generate_parameter_file(tnc::TensorNetworkCircuit,
                            filename::String
                            sliced_bonds::Array{<:Index, 1},
                            number_amplitudes::Int64,
                            seed::Union{Nothing, Int64}=nothing)

Generate a yml file with parameters corresponding to partitions and dimensions of each
along with the amplitudes to contract for.
"""
function generate_parameter_file(tnc::TensorNetworkCircuit,
                                 filename_prefix::String,
                                 sliced_bonds::Array{<:Index, 1},
                                 number_amplitudes::Int64,
                                 seed::Union{Nothing, Int64}=nothing)
    rng = MersenneTwister(seed)
    bitstring = n -> join(rand(rng, ['0', '1'], n))
    amplitude_strs = [bitstring(qubits(tnc)) for _ in 1:number_amplitudes]

    partition_dims = Dict{String, Int64}()
    for (i, sliced_bond) in enumerate(sliced_bonds)
        partition_dims["v$i"] = dim(sliced_bond)
    end
    partition_parameters = Dict("parameters" => partition_dims)
    config = Dict("partitions" => partition_parameters, "amplitudes" => unique(amplitude_strs))
    YAML.write_file("$(filename_prefix).yml", config)
end

end # module
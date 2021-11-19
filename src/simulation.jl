using Random
using YAML

export single_amplitude
export amplitudes_uniform, amplitudes_all
export generate_simulation_files, run_simulation
export generate_parameter_file, generate_dsl_files

"""
    amplitudes_all(qubits::Int)

Return generator for all the amplitudes of 'qubits' qubits
"""
function amplitudes_all(qubits::Int)
    (bitstring(x-1)[end-qubits+1:end] for x in 1:2^qubits)
end

"""
    amplitudes_uniform(qubits::Int, seed::Union{Int, Nothing}, number_amplitudes::Int)

Return generator to generate uniformly distributed bit strings of 'qubits' bits
(with replacement).
"""
function amplitudes_uniform(qubits::Int, seed::Union{Int, Nothing}, number_amplitudes::Int)
    rng = MersenneTwister(seed)
    bitstring = n -> join(rand(rng, ['0', '1'], n))
    (bitstring(qubits) for _ in 1:number_amplitudes)
end

"""
    generate_simulation_files(circ::QXZoo.Circuit.Circ,
                              output_prefix::String="simulation_input",
                              number_bonds_to_slice::Int=2;
                              decompose::Bool=true,
                              seed::Union{Int64, Nothing}=nothing,
                              output_args::Union{OrderedDict, Nothing}=nothing,
                              kwargs...)

Function to generate files required by qxrun to simulate the given circuit. This includes
a .tl file with a list of the operations, a .jld file with the initial tensors and a .yml
file with the parameters to use during the simulation.

# Keywords
- `number_bonds_to_slice::Int=2`: the number of edges to slice.
- `output_prefix::String="simulation_input"`: the prefix to be used for the simulation files.
- `num_amplitudes::Union{Int64, Nothing}=nothing`: the number of amplitudes to compute.
- `seed::Union{Int64, Nothing}=nothing`: the seed to be used by flow cutter to find a tree
                                         decomposition and for randomly selecting amplitudes to compute..
- `decompose::Bool=true`: set if two qubit gates should be decompoed when the circuit is converted to a tensor network.
- `kwargs`: all other kwargs are passed to `contraction_scheme` when it is called.
"""
function generate_simulation_files(circ::QXZoo.Circuit.Circ,
                                   output_prefix::String="simulation_input",
                                   number_bonds_to_slice::Int=2;
                                   decompose::Bool=true,
                                   seed::Union{Int64, Nothing}=nothing,
                                   output_args::Union{OrderedDict, Nothing}=nothing,
                                   kwargs...)

    @info("Convert circuit to tensor network")
    tnc = convert_to_tnc(circ; decompose=decompose)
    @info("Tensor network created with $(length(tnc)) tensors and $(length(bonds(tnc))) bonds")

    @info("Get contraction plan and edges to slice using QXGraphDecompositions")
    fc_seed = (seed === nothing) ? -1 : seed # flow cutter seed expects -1 in place of nothing.
    bond_groups, plan, metadata = contraction_scheme(tnc.tn, number_bonds_to_slice;
                                                     seed=fc_seed,
                                                     kwargs...)

    compute_tree = build_compute_graph(tnc, plan, bond_groups)

    @info("Prepare DSL and data files")
    generate_dsl_files(compute_tree, output_prefix; force=true, metadata=metadata)

    if output_args === nothing output_args = output_params_dict(qubits(tnc)) end
    generate_parameter_file(output_prefix, output_args)
end

"""
    single_amplitude(tnc::TensorNetworkCircuit, plan::Array{<:Index, 1}, amplitude::Union{String, Nothing}=nothing)

Contract the given tensor network using the given plan to calculate the amplitude of the given bitstring.
Creates a copy of the tensor network, replaces the outputs with those corresponding to the bitstring and then
contracts the network and returns the scalar amplitude.
"""
function single_amplitude(tnc::TensorNetworkCircuit, plan::Array{NTuple{3, Symbol}, 1}, amplitude::Union{String, Nothing}=nothing)
    sim_tnc = copy(tnc)
    add_output!(sim_tnc, amplitude)
    output = contract_tn!(sim_tnc, plan)
    output[1]
end

"""
    run_simulation(circ::QXZoo.Circuit.Circ;
                   num_amplitudes::Union{Int64, Nothing}=nothing,
                   seed::Union{Int64, Nothing}=nothing)

Function to run a simulation for the given circuit. If no number of amplitudes are given
then values for all possible amplitudes are found (only feasible for small systems). Returns
a dictionary of bitstring, amplitude pairs.
"""
function run_simulation(circ::QXZoo.Circuit.Circ;
                        num_amplitudes::Union{Int64, Nothing}=nothing,
                        seed::Union{Int64, Nothing}=nothing)
    @info("Convert circuit to tensor network")
    tnc = convert_to_tnc(circ)
    @info("Tensor network created with $(length(tnc)) tensors and $(length(bonds(tnc))) bonds")

    @info("Convert tensor network to graph")
    g = convert_to_graph(tnc)
    @info("Graph created: $(g)")

    @info("Get contraction plan and edges to slice using QXGraphDecompositions")
    plan = flow_cutter_contraction_plan(tnc; hypergraph=true)

    # to prevent memory issues don't attempt to get all amplitudes for large numbers of qubits
    # unless specifically demanded
    if num_amplitudes === nothing && qubits(tnc) > 30
        num_amplitudes = 1000
    end

    if num_amplitudes === nothing
        amplitudes = amplitudes_all(qubits(tnc))
    else
        amplitudes = amplitudes_uniform(qubits(tnc), seed, num_amplitudes)
    end

    results = OrderedDict{String, ComplexF64}()
    for amplitude in unique(amplitudes)
        results[amplitude] = single_amplitude(tnc, plan, amplitude)
    end
    results
end

"""
    generate_parameter_file(filename_prefix::String,
                            output_parameters)

Generate a yml file with details of how outputs are sampled

output:
    output_method: rejection
    params:
        fix_M: false
        M: 0.001
        num_samples: 10
        seed: ~
"""
function generate_parameter_file(filename_prefix::String,
                                 output_parameters)
    config = Dict()
    config["output"] = output_parameters

    YAML.write_file("$(filename_prefix).yml", config)
end

export amplitudes_uniform, amplitudes_all
export generate_simulation_files, run_simulation

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
    generate_simulation_files(circ::QXZoo.Circuit.Circ;
                              number_bonds_to_slice::Int=2,
                              output_prefix::String="simulation_input",
                              num_amplitudes::Union{Int64, Nothing}=nothing,
                              seed::Union{Int64, Nothing}=nothing)

Function to generate files required by qxrun to simulate the given circuit. This includes
a .tl file with a list of the operations, a .jld file with the initial tensors and a .yml
file with the parameters to use during the simulation.
"""
function generate_simulation_files(circ::QXZoo.Circuit.Circ;
                                   number_bonds_to_slice::Int=2,
                                   output_prefix::String="simulation_input",
                                   num_amplitudes::Union{Int64, Nothing}=nothing,
                                   seed::Union{Int64, Nothing}=nothing)
    @info("Convert circuit to tensor network")
    tnc = convert_to_tnc(circ)
    @info("Tensor network created with $(length(tnc)) tensors and $(length(bonds(tnc))) bonds")

    @info("Convert tensor network to graph")
    g = convert_to_graph(tnc)
    @info("Graph created: $(g)")

    @info("Get contraction plan and edges to slice using qxgraph")
    bonds_to_slice, plan = contraction_scheme(tnc.tn, number_bonds_to_slice)

    if num_amplitudes === nothing
        amplitudes = amplitudes_all(qubits(tnc))
    else
        amplitudes = amplitudes_uniform(qubits(tnc), seed, num_amplitudes)
    end

    @info("Write parameter file for retrieving $num_amplitudes amplitudes")
    generate_parameter_file(tnc, output_prefix, bonds_to_slice, amplitudes, seed)

    @info("Prepare DSL and data files")
    generate_dsl_files(tnc, output_prefix, plan, bonds_to_slice)
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

    @info("Get contraction plan and edges to slice using qxgraph")
    plan = quickbb_contraction_plan(tnc)

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



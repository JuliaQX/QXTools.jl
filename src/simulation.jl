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
                              seed::Union{Int64, Nothing}=nothing,
                              decompose::Bool=true,
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
function generate_simulation_files(circ::QXZoo.Circuit.Circ;
                                   number_bonds_to_slice::Int=2,
                                   decompose::Bool=true,
                                   output_method::Symbol=:list,
                                   output_prefix::String="simulation_input",
                                   seed::Union{Int64, Nothing}=nothing,
                                   num_outputs::Union{Int64, Nothing}=nothing,
                                   M::Float64=0.0001,
                                   fix_M::Bool=false,
                                   bitstrings::Union{Vector{String}, Nothing}=nothing,
                                   kwargs...)

    @info("Convert circuit to tensor network")
    tnc = convert_to_tnc(circ; decompose=decompose)
    @info("Tensor network created with $(length(tnc)) tensors and $(length(bonds(tnc))) bonds")

    @info("Get contraction plan and edges to slice using QXGraphDecompositions")
    fc_seed = (seed === nothing) ? -1 : seed # flow cutter seed expects -1 in place of nothing.
    bonds_to_slice, plan, metadata = contraction_scheme(tnc.tn, number_bonds_to_slice;
                                                        seed=fc_seed,
                                                        kwargs...)

    bond_groups_to_slice = expand_slice_bonds_to_hyperindices(tnc.tn, bonds_to_slice)

    @info("Write parameter file for retrieving $num_outputs amplitudes")
    output_args = OrderedDict()
    output_args[:method] = output_method
    output_params = OrderedDict()
    if output_method == :rejection
        output_params[:num_qubits] = tnc.qubits
        output_params[:M] = M
        output_params[:fix_M] = fix_M
        output_params[:seed] = seed
        output_params[:num_samples] = num_outputs === nothing ? 10 : num_outputs

    elseif output_method == :List
        if bitstrings === nothing
            bitstrings = amplitudes_all(qubits(tnc))
        end
        if num_outputs === nothing
            num_outputs = length(bitstrings)
        end
        output_params[:num_samples] = num_outputs
        output_params[:bitstrings] = collect(bitstrings)[1:num_outputs]

    elseif output_method == :uniform
        output_params[:num_qubits] = tnc.qubits
        output_params[:num_samples] = num_outputs === nothing ? 10 : num_outputs
        output_params[:seed] = seed
    end
    output_args[:params] = output_params
    generate_parameter_file(output_prefix, bond_groups_to_slice, output_args)

    @info("Prepare DSL and data files")
    generate_dsl_files(tnc, output_prefix, plan, bond_groups_to_slice; force=true, metadata=metadata)
end

"""
    expand_slice_bonds_to_hyperindices(tn::TensorNetwork, bonds_to_slice::Array{<: Index, 1})

Given a list of bonds to slice it is necessary to expand these to include bonds in the same hyper edge group. In
this function for each of the bonds to slice we identify their hyper indices if any and return an array of
groups of hyper edges to slice.
"""
function expand_slice_bonds_to_hyperindices(tn::TensorNetwork, bonds_to_slice::Array{<: Index, 1})
    bond_groups = Array{Array{<:Index, 1}, 1}()
    if length(bonds_to_slice) > 0
        push!(bond_groups, find_connected_indices(tn, bonds_to_slice[1]))
        if length(bonds_to_slice) > 1
            for b in bonds_to_slice[2:end]
                if !any([b in g for g in bond_groups])
                    push!(bond_groups, find_connected_indices(tn, b))
                end
            end
        end
    end
    bond_groups
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
    plan = quickbb_contraction_plan(tnc)

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

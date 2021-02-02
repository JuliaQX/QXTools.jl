module QXSim
using Logging
import QXZoo

# utilities circuit manipulation
export create_qft_circuit
include("circuits/circuits.jl")
using .Circuits

# data structions and functions tensor networks
export convert_to_tnc, next_tensor_id, convert_to_graph
export TensorNetwork, bonds, simple_contraction, tensor_data, neighbours
export TensorNetworkCircuit, qubits, add_input!, add_output!
export TensorCache, save_cache
export use_mock_tensors, contract_tensors
export quickbb_contraction_plan, contraction_scheme
include("tn/tn.jl")
using .TN

# data structures and functions for dealing with DSL and parameter files
export generate_dsl_files, generate_parameter_file
include("dsl/dsl.jl")
using .DSL

export generate_simulation_files

function generate_simulation_files(circ::QXZoo.Circuit.Circ;
                                   number_bonds_to_slice::Int=2,
                                   output_prefix::String="simulation_input",
                                   num_amplitudes::Int=100,
                                   seed::Union{Int64, Nothing}=nothing)
    @info("Convert circuit to tensor network")
    tnc = convert_to_tnc(circ)
    @info("Tensor network created with $(length(tnc)) tensors and $(length(bonds(tnc))) bonds")

    @info("Convert tensor network to graph")
    g = convert_to_graph(tnc)
    @info("Graph created: $(g)")

    @info("Get contraction plan and edges to slice using qxgraph")
    bonds_to_slice, plan = contraction_scheme(tnc.tn, number_bonds_to_slice)

    @info("Write parameter file for retrieving $num_amplitudes amplitudes")
    generate_parameter_file(tnc, output_prefix, bonds_to_slice, num_amplitudes, seed)

    @info("Prepare DSL and data files")
    generate_dsl_files(tnc, output_prefix, plan, bonds_to_slice)
end

end

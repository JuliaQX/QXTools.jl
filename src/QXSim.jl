module QXSim

# utilities circuit manipulation
export create_qft_circuit
include("circuits/circuits.jl")
using .Circuits

# data structions and functions tensor networks
export convert_to_tnc, next_tensor_id, convert_to_graph
export TensorNetwork, bonds, simple_contraction, tensor_data, neighbours
export TensorNetworkCircuit, qubits, add_input!, add_output!
export TensorCache, save_cache
include("tn/tn.jl")
using .TN

export generate_dsl_file
include("dsl/dsl.jl")
using .DSL
end

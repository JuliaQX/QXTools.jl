module QXSim

# utilities circuit manipulation
export create_qft_circuit
include("circuits/circuits.jl")
using .Circuits

# data structions and functions tensor networks
export convert_to_tnc
export convert_to_graph
include("tn/tn.jl")
using .TN

export generate_dsl_file
include("dsl/dsl.jl")
using .DSL

end

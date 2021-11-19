module QXTools
using Logging
using Reexport
using QXZoo
using QXTns
import DataStructures: OrderedDict

export TensorNetworkCircuit, add_input!, add_output!

# circuit manipulation
include("circuits/circuits.jl")
# circuit to tn conversion
include("tn_conversion.jl")
# contraction planning
include("contraction_planning.jl")
# dsl and parameter files
include("compute_graph/compute_graph.jl")
# functions for output arguments and parsing
include("outputs.jl")
# simulation utilities
include("simulation.jl")

@reexport using QXTools.Circuits

end

module QXSim

# utilities for handling circuits and conversion to tensor networks
# uses QXZOo
export convert_to_tnc
export convert_to_graph

include("tn/tn.jl")
include("tn/graph.jl")

end

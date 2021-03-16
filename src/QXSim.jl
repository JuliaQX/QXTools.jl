module QXSim
using Logging
import QXZoo
import DataStructures: OrderedDict

# utilities circuit manipulation
include("circuits/circuits.jl")

include("tn_conversion.jl")

include("contraction_planning.jl")

# data structures and functions for dealing with DSL and parameter files
include("dsl/dsl.jl")

include("simulation.jl")

end

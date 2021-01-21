using QXSim
using Test

@testset "QXSim.jl" begin
    include("test_tn.jl")
    include("test_dsl.jl")
end

using JLD

@testset "Test tensor cache" begin

tc = QXSim.TensorCache()
a = rand(Float64, 3, 4, 5)
sym = QXSim.add_tensor!(tc, a)
@test sym == QXSim.add_tensor!(tc, a)
@test sym != QXSim.add_tensor!(tc, rand(Float64, 3, 4, 5))
# test approximate matching
@test sym == QXSim.add_tensor!(tc, a .+ 0.1*eps(Float64))
@test sym != QXSim.add_tensor!(tc, a .+ 1.1*eps(Float64))

@test tc[sym] == a
@test length(tc) == 3

# test saving tensor cache
try
    QXSim.save("tmp.jld", tc)

    loaded_data = JLD.load("tmp.jld")
    @test all([(loaded_data[x] == tc[Symbol(x)]) for x in keys(loaded_data)])
finally
    rm("tmp.jld", force=true)
end

end
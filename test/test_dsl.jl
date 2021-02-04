using JLD

@testset "Test tensor cache" begin
    tc = QXSim.TensorCache()
    a = rand(Float64, 3, 4, 5)
    sym = push!(tc, a)
    @test sym == push!(tc, a)
    @test sym != push!(tc, rand(Float64, 3, 4, 5))
    # test approximate matching
    @test sym == push!(tc, a .+ 0.1*eps(Float64))
    @test sym != push!(tc, a .+ 1.1*eps(Float64))

    @test tc[sym] == a
    @test length(tc) == 3

    # test saving tensor cache
    mktempdir() do path
        save_cache(tc, joinpath(path, "tmp.jld"))
        loaded_data = JLD.load(joinpath(path, "tmp.jld"))
        @test all([(loaded_data[x] == tc[Symbol(x)]) for x in keys(loaded_data)])
    end
end
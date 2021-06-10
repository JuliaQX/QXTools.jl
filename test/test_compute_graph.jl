using JLD2
using FileIO
using QXTns

@testset "Test tensor cache" begin
    tc = QXTools.TensorCache()
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
        save_cache(tc, joinpath(path, "tmp.jld2"))
        loaded_data = load(joinpath(path, "tmp.jld2"))
        @test all([(loaded_data[x] == tc[Symbol(x)]) for x in keys(loaded_data)])
    end
end

@testset "Test build contraction tree" begin
    circ = create_test_circuit()
    tnc = convert_to_tnc(circ)
    plan = min_fill_contraction_plan(tnc)
    cg = build_compute_graph(tnc, plan)
    # test that number of nodes in compute graph is
    # tensors in network + length of plan + 1 for save node
    @test length(cg.root) == length(tnc) + length(plan) + 1
end
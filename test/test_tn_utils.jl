using QXSim.TN
using LinearAlgebra

@testset "Test the tensor network utlities" begin
    # test single qubit gates
    h = [[1., 1.] [1., -1.]]
    @test QXSim.TN.find_hyper_edges(h) == []
    z = [[1., 0.] [0., -1.]]
    @test QXSim.TN.find_hyper_edges(z) == [(1, 2)]

    # test 2 qubit identity
    id = collect(Diagonal(ones(4)))
    @test QXSim.TN.find_hyper_edges(id) == [(1, 2)]
    @test QXSim.TN.find_hyper_edges(reshape(id, (2, 2, 2, 2))) == [(1, 3), (2, 4)]
    B, C = QXSim.TN.decompose_gate(reshape(id, (2, 2, 2, 2)))
    @test size(B)[3] == 1
    @test size(C)[1] == 1
    @test QXSim.TN.find_hyper_edges(B) == [(1, 2)]
    @test QXSim.TN.find_hyper_edges(C) == [(2, 3)]

    # test 2 qubit cz
    id = collect(Diagonal(ones(4)))
    @test QXSim.TN.find_hyper_edges(id) == [(1, 2)]
    @test QXSim.TN.find_hyper_edges(reshape(id, (2, 2, 2, 2))) == [(1, 3), (2, 4)]
    B, C = QXSim.TN.decompose_gate(reshape(id, (2, 2, 2, 2)))
    @test size(B)[3] == 1
    @test size(C)[1] == 1
    @test QXSim.TN.find_hyper_edges(B) == [(1, 2)]
    @test QXSim.TN.find_hyper_edges(C) == [(2, 3)]
end
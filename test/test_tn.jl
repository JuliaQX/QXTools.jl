
using LightGraphs

@testset "Test circuit to Tensor Network Circuit conversion" begin
    circ = QXSim.create_test_circuit()
    tnc = convert_to_tnc(circ)

    @test tnc.qubits == 3
    @test all(tnc.input_indices .!= tnc.output_indices)
end

@testset "Test tensor network to graph conversion" begin
    # prepare circuit and network
    circ = QXSim.create_test_circuit()
    tnc = convert_to_tnc(circ)

    # convert to graph without inputs or outputs
    g = convert_to_graph(tnc)
    @test nv(g) == 3 # 3 vertices for gates
    @test ne(g) == 2 # 2 edges between gates on 2nd qubit

    # convert to graph with both inputs and outputs
    QXSim.add_input!(tnc)
    QXSim.add_output!(tnc)
    g = convert_to_graph(tnc)
    @test nv(g) == 9 # 3 gate, 3 input and 3 output
    @test ne(g) == 8 # 2 between gates and 6 to inputs and outputs

end
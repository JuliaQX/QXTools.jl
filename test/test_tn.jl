
using LightGraphs

@testset "Test circuit to Tensor Network Circuit conversion" begin
    circ = QXSim.create_test_circuit()
    tnc = convert_to_tnc(circ, no_input=true, no_output=true)

    @test tnc.qubits == 3
    @test length(tnc.tn.data) == 3
    @test all(tnc.input_indices .!= tnc.output_indices)

    tnc = convert_to_tnc(circ, no_input=false, no_output=true)
    @test length(tnc.tn.data) == 6
end

@testset "Test simple circuit contraction to verify convertion to circuit working" begin
    circ = QXSim.create_test_circuit()
    # we add input but no output to get full output vector
    tnc = convert_to_tnc(circ, no_input=false, no_output=true)

    output = QXSim.simple_contraction(tnc.tn)
    ref = zeros(8)
    ref[[1,8]] .= 1/sqrt(2)
    @test all(store(output) .≈ ref)
end

@testset "Test tensor network to graph conversion" begin
    # prepare circuit and network
    circ = QXSim.create_test_circuit()
    tnc = convert_to_tnc(circ, no_input=true, no_output=true)

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
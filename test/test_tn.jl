
using LightGraphs
using ITensors

@testset "Test circuit to Tensor Network Circuit conversion" begin
    circ = QXSim.create_test_circuit()

    # check tensornetwork size when no input or output
    tnc = convert_to_tnc(circ, no_input=true, no_output=true)
    @test tnc.qubits == 3
    @test length(tnc) == 3
    # make sure there are gates on all qubits as expected
    @test all(tnc.input_indices .!= tnc.output_indices)

    # convert again adding input
    tnc = convert_to_tnc(circ, no_input=false, no_output=true)
    @test length(tnc) == 6

    # check data of first tensor matches that of first gate
    @test QXSim.tensor_data(tnc, 1) == QXSim.gate_matrix(first(circ.circ_ops))
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

@testset "Test mock tensors" begin
    # create a mock tensor and check dimension reported correctly
    dim = [2 << 20, 2<< 20]
    a_store = QXSim.MockTensor{Complex{Float64}}(dim)
    @test length(a_store) == prod(dim)

    # create two tensors using mock tensor as storage and contract
    b_store = QXSim.MockTensor{Complex{Float64}}([dim[1], 2, 4]);
    a_inds = Index.(dim)
    a = ITensor(a_store, a_inds)
    b = ITensor(b_store, [a_inds[1], Index(2), Index(4)])
    c = a * b
    @test length(store(c)) == dim[2] * 2 * 4

    # test replacing all stores with mock tensor in a tensor network
    circ = QXSim.create_test_circuit()
    tnc = convert_to_tnc(circ)
    new_tn = QXSim.use_mock_tensors(tnc.tn)
    @test all([typeof(store(x)) <: QXSim.MockTensor for x in new_tn.data])
end
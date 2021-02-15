
using QXSim
using QXSim.TN
using QXSim.Circuits

import QXGraph
import LightGraphs
using ITensors

@testset "Test TensorNetwork struct and interface functions" begin
    # create empty tensor network
    tn = TensorNetwork()

    # create two tensors with one common index
    i1 = Index(2)
    i2 = Index(3)
    a_inds = [i1, Index(4)]
    b_inds = [i1, i2]
    tensor_a = push!(tn, a_inds, rand(2, 4))
    tensor_b = push!(tn, b_inds, rand(2, 3))
    @test length(tn) == 2
    @test length(collect(tn)) == 2

    # check the number of bonds 2
    # and that i1 connects to two tensors
    @test length(bonds(tn)) == 3
    @test length(tn[i1]) == 2
    @test all(tn[i1] .== [tensor_a, tensor_b])

    @test length(keys(tn)) == 2
    @test length(tn) == 2

    @test size(tensor_data(tn, tensor_a)) == (2, 4)

    tn2 = TensorNetwork()
    tensor_c = push!(tn2, [i2], rand(3))
    tn3 = merge(tn, tn2)
    @test length(tn3) == 3
    @test length(tn3[i2]) == 2

    @test length(neighbours(tn, tensor_b)) == 1
    @test length(neighbours(tn3, tensor_b)) == 2

    # Test delete!
    delete!(tn, tensor_a)
    @test length(tn) == 1
    @test haskey(tn.tensor_map, tensor_a) == false
    @test length(neighbours(tn, tensor_b)) == 0

    # Test TN copy
    a = TensorNetwork()
    tensor_a = push!(a, a_inds, rand(2, 4))
    tensor_b = push!(a, b_inds, rand(2, 3))
    b = copy(a)
    @test length(a) == length(b)
    @test !haskey(b, tensor_a) # test that symbols have changed
    delete!(a, tensor_a)
    @test length(a) != length(b)
end

@testset "Test TensorNetworkCircuit struct and interface functions" begin
    # create empty tensor network
    tnc = TensorNetworkCircuit(2)

    @test qubits(tnc) == 2
    @test length(tnc) == 0

    # add single qubit gate
    push!(tnc, [1], rand(2, 2))
    @test length(tnc) == 1

    # add two qubit gate
    push!(tnc, [1, 2], rand(4, 4), decompose=false)
    @test length(tnc) == 2

    # adding input
    add_input!(tnc, "01")
    @test length(tnc) == 4
    @test tensor_data(tnc, input_tensors(tnc)[1]) == [1., 0.]
    @test tensor_data(tnc, input_tensors(tnc)[2]) == [0., 1.]

    # adding output
    add_output!(tnc, "+-")
    @test length(tnc) == 6
    @test tensor_data(tnc, output_tensors(tnc)[1]) == [1., 1.]/sqrt(2)
    @test tensor_data(tnc, output_tensors(tnc)[2]) == [1., -1.]/sqrt(2)

    # test Base.copy on tnc
    tnc2 = copy(tnc)
    @test tensor_data(tnc2, output_tensors(tnc2)[1]) == [1., 1.]/sqrt(2)
    @test tensor_data(tnc2, output_tensors(tnc2)[2]) == [1., -1.]/sqrt(2)
    @test length(tnc) == length(tnc2)
end

@testset "Test circuit to Tensor Network Circuit conversion" begin
    circ = create_test_circuit()

    # check tensornetwork size when no input or output
    tnc = convert_to_tnc(circ, no_input=true, no_output=true, decompose=false)
    @test qubits(tnc) == 3
    @test length(tnc) == 3
    # make sure there are gates on all qubits as expected
    @test all(tnc.input_indices .!= tnc.output_indices)

    tnc = convert_to_tnc(circ, no_input=true, no_output=true, decompose=true)
    @test qubits(tnc) == 3
    @test length(tnc) == 5
    # make sure there are gates on all qubits as expected
    @test all(tnc.input_indices .!= tnc.output_indices)

    # convert again adding input
    tnc = convert_to_tnc(circ, no_input=false, no_output=true, decompose=false)
    @test length(tnc) == 6

    # check data of first tensor matches that of first gate
    @test tensor_data(tnc, first(keys(tnc))) == gate_matrix(first(circ.circ_ops))
end

@testset "Test simple circuit contraction to verify convertion to circuit working" begin
    circ = create_test_circuit()
    # we add input but no output to get full output vector
    tnc = convert_to_tnc(circ, no_input=false, no_output=true)

    output = simple_contraction(tnc)
    ref = zeros(8)
    ref[[1,8]] .= 1/sqrt(2)
    @test all(output .≈ ref)

    tnc = convert_to_tnc(circ, no_input=false, no_output=true, decompose=false)

    output = simple_contraction(tnc)
    ref = zeros(8)
    ref[[1,8]] .= 1/sqrt(2)
    @test all(output .≈ ref)
end


@testset "Test mock tensors" begin
    # create a mock tensor and check dimension reported correctly
    dim = [2 << 20, 2<< 20]
    a_store = MockTensor{Complex{Float64}}(dim)
    @test length(a_store) == prod(dim)

    # create two tensors using mock tensor as storage and contract
    b_store = MockTensor{Complex{Float64}}([dim[1], 2, 4]);
    a_inds = Index.(dim)
    a = ITensor(a_store, a_inds)
    b = ITensor(b_store, [a_inds[1], Index(2), Index(4)])
    c = contract_tensors(a, b)
    @test length(store(c)) == dim[2] * 2 * 4

    # test replacing all stores with mock tensor in a tensor network
    circ = create_test_circuit()
    tnc = convert_to_tnc(circ)
    new_tn = use_mock_tensors(tnc.tn)
    @test all([typeof(store(x)) <: MockTensor for x in new_tn])
end

@testset "Test tensor decomnposition" begin
    # prepare the circuit.
    circ = create_test_circuit()
    tnc = convert_to_tnc(circ, no_input=false, no_output=false, decompose=false)
    @assert length(tnc.tn.tensor_map) == 9

    # Find a tensor with four indices to decompose.
    t_id = :_
    for id in keys(tnc.tn.tensor_map)
        if length(inds(tnc.tn.tensor_map[id])) == 4
            t_id = id
            break
        end
    end
    @assert length(inds(tnc.tn.tensor_map[t_id])) == 4
    left_inds = collect(inds(tnc.tn.tensor_map[t_id]))
    left_inds = left_inds[1:2]

    # test svd
    Uid, Sid, Vid = replace_with_svd!(tnc.tn, t_id, left_inds;
                                            maxdim=2,
                                            cutoff=1e-13)
    @test length(tnc.tn.tensor_map) == 9 + 2

    SVid = contract_pair!(tnc.tn, Sid, Vid)
    USVid = contract_pair!(tnc.tn, Uid, SVid)
    @assert length(tnc.tn.tensor_map) == 9

    Uid, Vid = decompose_tensor!(tnc, USVid, left_inds;
                                        maxdim=2,
                                        cutoff=1e-13)
    @test length(tnc.tn.tensor_map) == 9 + 1
end
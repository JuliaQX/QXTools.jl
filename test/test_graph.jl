using QXSim.Circuits
using QXSim.TN
using QXGraph

using LightGraphs
using ITensors

@testset "Test tensor network to graph conversion" begin
    # prepare circuit and network
    circ = create_test_circuit()
    tnc = convert_to_tnc(circ, no_input=true, no_output=true)

    # convert to graph without inputs or outputs
    g = convert_to_graph(tnc)
    @test LightGraphs.nv(g) == 3 # 3 vertices for gates
    @test LightGraphs.ne(g) == 2 # 2 edges between gates on 2nd qubit

    # convert to graph with both inputs and outputs
    add_input!(tnc)
    add_output!(tnc)
    g = convert_to_graph(tnc)
    @test LightGraphs.nv(g) == 9 # 3 gate, 3 input and 3 output
    @test LightGraphs.ne(g) == 8 # 2 between gates and 6 to inputs and outputs


    # prepare the circuit.
    circ = create_test_circuit()
    tnc = convert_to_tnc(circ, no_input=true, no_output=true)
    add_input!(tnc)
    add_output!(tnc)

    # create the line graph for the circuit.
    g, symbol_map = convert_to_line_graph(tnc.tn)
    @test QXGraph.nv(g) == 8 # 2 between gates and 6 to inputs and outputs
    @test QXGraph.ne(g) == 13 # 1-qubit gate -> 1 edge, 2-qubit gate -> 6 edges. 1+6+6.

    # Check conversion to linegraph of network's hypergraph.
    tnc = TensorNetworkCircuit(2)
    push!(tnc, [1], rand(2, 2))
    push!(tnc, [1, 2], rand(4, 4); diagonal=true)
    push!(tnc, [1], rand(2, 2); diagonal=true)
    g, symbol_map = convert_to_line_graph(tnc.tn; use_tags=true)
    @test QXGraph.nv(g) == 3
    @test QXGraph.ne(g) == 2
    @test length(symbol_map) == 3
    @test sum([typeof(inds) <:Index ? 1 : length(inds) for inds in values(symbol_map)]) == 6
end


@testset "Test QXGraph contraction and slicing" begin
    # prepare the circuit.
    circ = create_test_circuit()
    tnc = convert_to_tnc(circ, no_input=false, no_output=true)

    # test contraction plan
    plan = quickbb_contraction_plan(tnc)
    @test length(plan) == 8 # This includes open indices.
    @test length(unique(plan)) == 8

    # test contracting the network
    output = contract_tn!(tnc, plan)
    ref = zeros(8); ref[[1,8]] .= 1/sqrt(2)
    @test all(store(output) .≈ ref)

    # prepare another, un-contracted circuit.
    circ = create_test_circuit()
    tnc = convert_to_tnc(circ, no_input=false, no_output=false)

    # Test contraction scheme function
    edges_to_slice, plan = contraction_scheme(tnc, 3)
    @test length(edges_to_slice) == 3 # Should have 3 edges to slice
    @test length(plan) == length(tnc.tn.bond_map) - 3 # modified plan should be smaller.
end

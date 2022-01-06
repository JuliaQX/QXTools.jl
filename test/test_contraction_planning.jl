using QXTools.Circuits
using QXGraphDecompositions
using QXZoo
using QXTns

using LightGraphs

@testset "Test tensor network to graph conversion" begin
    # prepare circuit and network
    circ = create_test_circuit()
    tnc = convert_to_tnc(circ, no_input=true, no_output=true, decompose=false)

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
    tnc = convert_to_tnc(circ, no_input=true, no_output=true, decompose=false)
    add_input!(tnc)
    add_output!(tnc)

    # create the line graph for the circuit.
    g, symbol_map = convert_to_line_graph(tnc.tn)
    @test QXGraphDecompositions.nv(g) == 8 # 2 between gates and 6 to inputs and outputs
    @test QXGraphDecompositions.ne(g) == 13 # 1-qubit gate -> 1 edge, 2-qubit gate -> 6 edges. 1+6+6.

    # Check conversion to linegraph of network's hypergraph.
    circ = create_test_circuit()
    tnc = convert_to_tnc(circ, no_input=false, no_output=false, decompose=true)
    g, symbol_map = convert_to_line_graph(tnc; use_hyperedges=true)
    @test QXGraphDecompositions.nv(g) == 6 # circuit has 6 hyperedges
    @test QXGraphDecompositions.ne(g) == 7
    @test length(symbol_map) == 6
    @test QXTns.counter(length.(values(symbol_map))) == Dict(4=>2, 2=>4)
end


@testset "Test QXGraphDecompositions contraction and slicing" begin
    # prepare the circuit.
    circ = create_test_circuit()
    tnc = convert_to_tnc(circ, no_input=false, no_output=true, decompose=false)

    # test contraction plan
    plan = flow_cutter_contraction_plan(tnc; time=15)
    @test length(plan) == 5

    # test contracting the network
    output = contract_tn!(tnc, plan)
    output = reshape(output, prod(size(output)))
    ref = zeros(8); ref[[1,8]] .= 1/sqrt(2)
    @test all(output .≈ ref)

    # prepare another, un-contracted circuit.
    circ = create_test_circuit()
    tnc = convert_to_tnc(circ, no_input=false, no_output=false, decompose=false)

    # Test contraction scheme function
    edges_to_slice, plan = contraction_scheme(tnc, 3; time = 15, hypergraph=false)
    @test length(vcat(edges_to_slice...)) >= 3 # Should have 3 edges to slice
    @test length(plan) == length(tnc) - 1 - 3 # modified plan should be smaller.

    # Test contraction scheme function
    edges_to_slice, plan = contraction_scheme(tnc, 0; time = 15)
    @test length(edges_to_slice) == 0 # Should have 0 edges to slice
    @test length(plan) == length(tnc) - 1

    # Test contraction plan creation with hypergraph.
    circ = create_test_circuit()
    tnc = convert_to_tnc(circ, no_input=false, no_output=false, decompose=false)
    plan = flow_cutter_contraction_plan(tnc; time=30, hypergraph=true)
    @test length(plan) == 8

    # Test contraction scheme function with hypergraph.
    edges_to_slice, plan = contraction_scheme(tnc, 3; time = 15, hypergraph=true)
    # Should have 5 indices to slice (1 hyperedge with 2 indices)
    @test length(vcat(edges_to_slice...)) >= 5
    @test length(plan) == 3 # modified plan should be smaller.

    # Create a network consisting of one hyperedge which is too large for netcon to
    # contract. Test if fallback method for contracting large hyperedges produces a valid
    # contraction plan.
    circ = Circuit.Circ(1)
    for i in 1:40
        Circuit.add_gatecall!(circ, DefaultGates.z(1))
    end
    tnc = convert_to_tnc(circ, no_input=false, no_output=false, decompose=false)
    plan = flow_cutter_contraction_plan(tnc; time = 15, hypergraph=true)
    @test length(plan) == 41
end

@testset "Test flow cutter contraction" begin
    # A test circuit to test flow cutter on.
    circ = create_test_circuit()
    tnc = convert_to_tnc(circ, no_input=false, no_output=true, decompose=false)

    # test contraction plan
    plan = flow_cutter_contraction_plan(tnc; time=15)
    @test length(plan) == 5

    # test contracting the network
    output = contract_tn!(tnc, plan)
    output = reshape(output, prod(size(output)))
    ref = zeros(8); ref[[1,8]] .= 1/sqrt(2)
    @test all(output .≈ ref)

    # A test circuit to test flow cutter on.
    circ = create_test_circuit()
    tnc = convert_to_tnc(circ, no_input=false, no_output=true, decompose=false)

    # test contraction plan with 0 seconds so the min fill heursitic is used as a fall back
    # when flow cutter can't find a tree decomposition.
    plan = flow_cutter_contraction_plan(tnc; time=0)
    @test length(plan) == 5

    # test contracting the network
    output = contract_tn!(tnc, plan)
    output = reshape(output, prod(size(output)))
    ref = zeros(8); ref[[1,8]] .= 1/sqrt(2)
    @test all(output .≈ ref)
end

@testset "Test min fill contraction" begin
    # A test circuit to test the min fill heuristic on.
    circ = create_test_circuit()
    tnc = convert_to_tnc(circ, no_input=false, no_output=true, decompose=false)

    # test contraction plan
    plan = min_fill_contraction_plan(tnc)
    @test length(plan) == 5

    # test contracting the network
    output = contract_tn!(tnc, plan)
    output = reshape(output, prod(size(output)))
    ref = zeros(8); ref[[1,8]] .= 1/sqrt(2)
    @test all(output .≈ ref)
end

@testset "Test netcon contraction" begin
    # prepare the circuit.
    circ = create_test_circuit()
    tnc = convert_to_tnc(circ, no_input=false, no_output=true, decompose=false)

    # test contraction plan
    plan = netcon(tnc)
    @test length(plan) == 5

    # test contracting the network
    output = contract_tn!(tnc, plan)
    output = reshape(output, prod(size(output)))
    ref = zeros(8); ref[[1,8]] .= 1/sqrt(2)
    @test all(output .≈ ref)
end

@testset "Test elimination order to contraction plan conversion" begin
    # Prepare the circuit.
    circ = create_test_circuit()
    tnc = convert_to_tnc(circ, no_input=false, no_output=false, decompose=true)

    # Create the line graphs for the network's graph and hypergraph.
    lg, symbol_map = convert_to_line_graph(tnc)
    hyper_lg, hyper_symbol_map = convert_to_line_graph(tnc; use_hyperedges=true)

    # Get an elimination order for the line graphs.
    tw, order = QXTools.qxg.min_fill(lg)
    order = [symbol_map[edge] for edge in order]
    hyper_tw, hyper_order = QXTools.qxg.min_fill(hyper_lg)
    hyper_order = [hyper_symbol_map[edge] for edge in hyper_order]

    # Try converting the orders to contraction plans
    plan = QXTools.order_to_contraction_plan(order, tnc.tn)
    hyper_plan = QXTools.order_to_contraction_plan(hyper_order, tnc.tn)

    # Test contraction plan has correct length.
    @test length(plan) == 10
    @test length(hyper_plan) == 10
end
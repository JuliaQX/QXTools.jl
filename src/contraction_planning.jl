import QXGraph; qxg = QXGraph
import LightGraphs
using QXTn

export convert_to_graph, convert_to_line_graph
export quickbb_contraction_plan, contraction_scheme
export netcon

# **************************************************************************************** #
#                             Converting tensor networks to graphs
# **************************************************************************************** #

"""
    convert_to_graph(tn::TensorNetwork)

Return a SimpleGraph representing the graph structure of the tensor network `tn`.
"""
function convert_to_graph(tn::TensorNetwork)
    g = LightGraphs.SimpleGraph(length(tn))
    tensor_map = Dict{Symbol, Int64}([k => i for (i, k) in enumerate(keys(tn))])
    # TODO: add structures to tensor network to avoid N^2 complexity here
    for bond in bonds(tn)
        tensors = tn[bond]
        if length(tensors) == 2
            LightGraphs.add_edge!(g, tensor_map[tensors[1]], tensor_map[tensors[2]])
        end
    end
    g
end

convert_to_graph(tnc::TensorNetworkCircuit) = convert_to_graph(tnc.tn)

"""
    convert_to_line_graph(tn::TensorNetwork; use_tags::Bool=false)

Create a labeled graph representing the line graph of 'tn'.

If use_tags is false, each Index struct in 'tn' is converted to a symbol and used as a label 
for the corresponding vertex in the line graph. A dictionary is also returned which maps 
these symbols to the original Index structs.

# Keywords
- `use_tags::Bool=false`: set it index tags should be used as line graph vertex labels.
"""
function convert_to_line_graph(tn::TensorNetwork; use_tags::Bool=false)
    # Create the line graph.
    g = qxg.LabeledGraph()

    # Add a vertex for each edge/hyperedge in tn.
    symbol_map = Dict{Symbol, Array{<:Index, 1}}()
    for index in bonds(tn)
        # Create the symbol to label the vertex representing the index
        index_symbol = index_to_symbol(index; use_tags=use_tags)

        # If the index belongs to a hyperedge, the corresponding vertex may already exist.
        v = qxg.get_vertex(g, index_symbol)
        if v === nothing
            qxg.add_vertex!(g, index_symbol)
        end

        # Add index to dictionary keeping track of how symbols represent sets of indices.
        push!(get!(symbol_map, index_symbol, Array{Index, 1}()), index)
    end

    # Connect all vertices whose indices belong to the same tensor in tn. Each tensor in tn 
    # should then be represented as a clique in the line graph.
    for tensor in values(tn)
        indices = inds(tensor)
        for i = 1:length(indices)-1
            for j = i+1:length(indices)
                vi = index_to_symbol(indices[i]; use_tags=use_tags)
                vj = index_to_symbol(indices[j]; use_tags=use_tags)
                if !(vi==vj)
                    qxg.add_edge!(g, vi, vj)
                end
            end
        end
    end
    g, symbol_map
end

"""
    index_to_symbol(ind::Index; use_tags::Bool=false)

Convert the Index `ind` to a symbol.
"""
function index_to_symbol(ind::QXTn.Index; use_tags::Bool=false)
    if use_tags
        symb = String(ind.tags[1])*"_"*String(ind.tags[2])
        return Symbol(symb)
    else
        return Symbol(ind)
    end
end



# **************************************************************************************** #
#                                QuickBB contraction
# **************************************************************************************** #

"""
    quickbb_contraction_plan(tn::TensorNetwork)

Use QuickBB to create a contraction plan for 'tn'.

# Keywords
- `time::Integer=0`: the number of second to run the quickbb binary for.
- `order::Symbol=:_`: the branching order to be used by quickbb (:random or :min_fill).
- `hypergraph::Bool=false`: set if hyperedges exist in `tn` and shoudl be accounted for.
"""
function quickbb_contraction_plan(tn::TensorNetwork; 
                                  time::Integer=120, 
                                  order::Symbol=:_,
                                  hypergraph::Bool=false)
    # Convert tn to a line graph and pass it to quickbb to find an elimination order.
    lg, symbol_map = convert_to_line_graph(tn, use_tags=hypergraph)
    tw, order = qxg.quickbb(lg; time=time, order=order)

    # Convert the elimination order to an array of Index structs in tn.
    order = [symbol_map[index_symbol] for index_symbol in order]

    # Convert the elimination order into a contraction plan.
    order_to_contraction_plan(order, tn)
end

quickbb_contraction_plan(tnc::TensorNetworkCircuit; kwargs...) = quickbb_contraction_plan(tnc.tn; kwargs...)



# **************************************************************************************** #
#                         Contraction with automatic slicing
# **************************************************************************************** #

"""
    contraction_scheme(tn::TensorNetwork, num::Integer)

Return an array of 'num' indices to slice in the tensor network 'tn' and a contraction plan 
for the remaining tensor network.

# Keywords
- `time::Integer=120`: number of seconds to run quickbb when looking for contraction plan.
- `order::Symbol=:_`: the branching order to be used by quickbb (:random or :min_fill).
- `score_function::Symbol=:direct_treewidth`: function to maximise when selecting vertices 
                                            to remove. (:degree, :direct_treewidth)
- `hypergraph::Bool=false`: set if hyperedges exist in `tn` and shoudl be accounted for.
"""
function contraction_scheme(tn::TensorNetwork, num::Integer; 
                            time::Integer=120,
                            order::Symbol=:min_fill,
                            score_function::Symbol=:direct_treewidth,
                            hypergraph::Bool=false)
    # Create the line graph for the given tn and pass it to quickbb to find a contraction 
    # plan.
    lg, symbol_map = convert_to_line_graph(tn; use_tags=hypergraph)
    tw, order = qxg.quickbb(lg; time=time, order=order)

    # Use the greedy treewidth deletion algorithm to select indices in tn to slice.
    scheme = qxg.greedy_treewidth_deletion(lg, num; 
                                           score_function=:direct_treewidth, 
                                           elim_order=order)
    sliced_lg, edges_to_slice, modified_order, new_tw = scheme

    # Convert the contraction plan to an array of Index structs before returning.
    modified_order = [symbol_map[index_symbol] for index_symbol in modified_order]
    edges_to_slice = [symbol_map[index_symbol] for index_symbol in edges_to_slice]
    contraction_plan = order_to_contraction_plan(modified_order, tn)

    # TODO: concatenating index arrays discards distinguishability of edges and hyperedges
    # which may not be desirable. We may not want to slice hyperedges in the same way we 
    # slice regular edges. 
    edges_to_slice = vcat(edges_to_slice...)
    edges_to_slice, contraction_plan
end

contraction_scheme(tnc::TensorNetworkCircuit, args...; kwargs...) = contraction_scheme(tnc.tn, args...; kwargs...)



# **************************************************************************************** #
#                                   Netcon functions
# **************************************************************************************** #
import TensorOperations.optimaltree

"""
    netcon(tensor_map::OrderedDict{Symbol, QXTensor}

Call the TensorOperations netcon implementation on the tensors in `tensor_map`.
"""
function netcon(tensor_map::OrderedDict{Symbol, QXTensor})
    # Generate the arguments for the netcon method from the given network.
    tensor_indices, index_dims, tensor_labels = create_netcon_input(tensor_map)

    # Run the netcon method to find the optimal contraction plan.
    contraction_tree, cost = optimaltree(tensor_indices, index_dims)

    # Replace integers in contraction tree with node labels before returning
    # _replace_leaf_integers_with_labels(contraction_tree, tensor_labels)
    plan, final_tensor_id = tree_to_contraction_plan(contraction_tree, tensor_labels)
    plan
end

netcon(tn::TensorNetwork) = netcon(OrderedDict(keys(tn) .=> values(tn)))
netcon(tnc::TensorNetworkCircuit) = netcon(tnc.tn)

"""
    create_netcon_input(tensor_map::OrderedDict{Symbol, QXTensor})

Create arguments for netcon to find the optimal contraction plan for the tensors
contatined in `tensor_map`.
"""
function create_netcon_input(tensor_map::OrderedDict{Symbol, QXTensor})
    # Collect the indices and corresponding dimensions of the network into arrays.
    num_tensors = length(tensor_map)

    tensor_indices = Array{Array{Symbol, 1}, 1}(undef, num_tensors)
    tensor_dims = Array{Array{Int, 1}, 1}(undef, num_tensors)
    tensor_labels = Array{Symbol, 1}(undef, num_tensors)

    for (i, (tensor_id, tensor)) in enumerate(tensor_map)
        tensor_indices[i] = [index_to_symbol(ind) for ind in inds(tensor)]
        tensor_dims[i] = [dim(ind) for ind in inds(tensor)]
        tensor_labels[i] = tensor_id
    end

    # Make a dictionary of index dimensions for netcon.
    labels = reduce(vcat, tensor_indices); dims = reduce(vcat, tensor_dims)
    index_dims = Dict{Symbol, Int}(labels .=> dims)

    # The netcon implementation only needs tensor_indices and index_dims as
    # input but it returns a nested array of integers which need to be replaced
    # by node labels for contraction. To this end, tensor_labels is also
    # returned.
    tensor_indices, index_dims, tensor_labels
end



# **************************************************************************************** #
#       Converting edge elimination orders & contraction trees into contraction plans
# **************************************************************************************** #

"""
    order_to_contraction_plan(elimination_order::Array{<:Array{<:Index, 1}, 1}, 
                              tn::TensorNetwork)

Convert the given edge elimination order into a contraction plan for `tn`.
"""
function order_to_contraction_plan(elimination_order::Array{<:Array{<:Index, 1}, 1}, 
                                   tn::TensorNetwork)
    plan = Array{NTuple{3, Symbol}, 1}()
    # TODO: Might be best to have an option in copy to replace tensors with mocktensors here.
    new_tn = copy(tn)

    # Render each edge in the order into a contraction of tensors in tn.
    for edge in elimination_order
        # If the edge is a regular edge with two tensors attached, contract those tensors.
        if length(edge) == 1
            index = edge[1]
            if haskey(new_tn.bond_map, index) && (length(new_tn[index]) == 2)
                A_id, B_id = new_tn[index]
                C_id = contract_pair!(new_tn, A_id, B_id)
                append!(plan, [(A_id, B_id, C_id)])
            end

        # If the edge is a hyperedge, use netcon to find a contraction plan for the
        # corresponding tensors and append it to the contraction plan for tn.
        else
            tensors = OrderedDict{Symbol, QXTensor}()
            for index in edge
                if haskey(new_tn, index)
                    tensor_pair = new_tn[index]
                    for t in tensor_pair
                        tensors[t] = new_tn[t]
                    end
                end
            end
            # TODO: should have a strategy available if netcon can't be used here.
            @assert length(tensors) <= 36 "Netcon can only contract upto 36 tensors"
            local_contraction_plan = netcon(tensors)
            for (a, b, c) in local_contraction_plan
                contract_pair!(new_tn, a, b, c)
            end
            append!(plan, local_contraction_plan)
        end
    end
    plan
end

"""
    tree_to_contraction_plan(tree::Array{Any, 1}, labels::Array{Symbol, 1})::Tuple{Array{NTuple{3, Symbol}, 1}, Symbol}

Return the contraction plan described by the netcon contraction tree stored in 'tree'. The
tensors labeled by integers in `tree` are labeled by the corresponding symbol in `labels` in
the returned contraction plan.
"""
function tree_to_contraction_plan(tree::Array{Any, 1}, labels::Array{Symbol, 1})::Tuple{Array{NTuple{3, Symbol}, 1}, Symbol}
    plan_A, A_id = tree_to_contraction_plan(tree[1], labels)
    plan_B, B_id = tree_to_contraction_plan(tree[2], labels)
    C_id = next_tensor_id()
    return [plan_A; plan_B; NTuple{3, Symbol}[(A_id, B_id, C_id)]], C_id
end

function tree_to_contraction_plan(tree::Int, labels::Array{Symbol, 1})::Tuple{Array{NTuple{3, Symbol}, 1}, Symbol}
    [], labels[tree]
end
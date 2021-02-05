import QXGraph; qxg = QXGraph
import LightGraphs
import ITensors

export convert_to_graph, convert_to_line_graph
export quickbb_contraction_plan, contraction_scheme

convert_to_graph(tnc::TensorNetworkCircuit) = convert_to_graph(tnc.tn)

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
    # Create the graph which will be the line graph.
    g = qxg.LabeledGraph()

    # Add a vertex for each edge/hyperedge in tn.
    symbol_map = Dict{Symbol, Union{Index, Array{<:Index, 1}}}()
    for index in bonds(tn)
        index_symbol = index_to_symbol(index; use_tags=use_tags)
        v = qxg.get_vertex(g, index_symbol)
        if v === nothing
            qxg.add_vertex!(g, index_symbol)
        end

        if haskey(symbol_map, index_symbol)
            symbol_map[index_symbol] = vcat(symbol_map[index_symbol], index)
        else
            symbol_map[index_symbol] = index
        end
    end

    # Connect all vertices whose indices belong to the same tensor in tn.
    # Each tensor in tn should then be represented as a clique in the line graph.
    for tensor in values(tn.tensor_map)
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

Convert the ITensors Index `ind` to a symbol.
"""
function index_to_symbol(ind::Index; use_tags::Bool=false)
    if use_tags
        symb = String(ind.tags[1])*"_"*String(ind.tags[2])
        return Symbol(symb)
    else
        return Symbol(ind)
    end
end


"""
    quickbb_contraction_plan(tn::TensorNetwork)

Use QuickBB to create a contraction plan for 'tn'.

# Keywords
- `time::Integer=0`: the number of second to run the quickbb binary for.
- `order::Symbol=:_`: the branching order to be used by quickbb (:random or :min_fill).
"""
function quickbb_contraction_plan(tn::TensorNetwork; time::Integer=120, order::Symbol=:_)
    # Convert tn to a line graph and pass it to quickbb to find a contraction plan.
    lg, symbol_map = convert_to_line_graph(tn)
    tw, plan = qxg.quickbb(lg; time=time, order=order)

    # Convert the contraction plan to an array of Index structs in tn.
    [symbol_map[index_symbol] for index_symbol in plan]
end


"""
    contraction_scheme(tn::TensorNetwork, num::Integer)

Return an array of 'num' indices to slice in the tensor network 'tn' and a contraction plan 
for the remaining tensor network.

# Keywords
- `time::Integer=120`: number of seconds to run quickbb when looking for contraction plan.
- `order::Symbol=:_`: the branching order to be used by quickbb (:random or :min_fill).
- `score_function::Symbol=:direct_treewidth`: function to maximise when selecting vertices 
                                            to remove. (:degree, :direct_treewidth)
"""
function contraction_scheme(tn::TensorNetwork, num::Integer; 
                            time::Integer=120,
                            order::Symbol=:min_fill,
                            score_function::Symbol=:direct_treewidth)
    # Create the line graph for the given tn and pass it to quickbb to find a contraction 
    # plan.
    lg, symbol_map = convert_to_line_graph(tn)
    tw, plan = qxg.quickbb(lg; time=time, order=order)

    # Use the greedy treewidth deletion algorithm to select indices in tn to slice.
    scheme = qxg.greedy_treewidth_deletion(lg, num; 
                                           score_function=:direct_treewidth, 
                                           elim_order=plan)
    sliced_lg, edges_to_slice, modified_plan, new_tw = scheme

    # Convert the contraction plan to an array of Index structs before returning.
    modified_plan = [symbol_map[index_symbol] for index_symbol in modified_plan]
    edges_to_slice = [symbol_map[index_symbol] for index_symbol in edges_to_slice]
    edges_to_slice, modified_plan
end
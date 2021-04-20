import QXGraphDecompositions; qxg = QXGraphDecompositions
import LightGraphs
import QXTns

export convert_to_graph, convert_to_line_graph
export quickbb_contraction_plan, contraction_scheme
export flow_cutter_contraction_plan, min_fill_contraction_plan
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
    convert_to_line_graph(tn::TensorNetwork; use_hyperedges::Bool=false)

Create a labeled graph representing the line graph of 'tn'.

If `use_hyperedges` is false, each edge in `tn` is assigned a symbol which labels a unique
vertex in the line graph. If it is true, the hyperedges of `tn` are assigned symbols instead.

A dictionary is also returned which maps these symbols to the corresponding edges which are
represented as tensor sets.

# Keywords
- `use_hyperedges::Bool=false`: set if hyperedges should be used to create the line graph.
"""
function convert_to_line_graph(tn::TensorNetwork; use_hyperedges::Bool=false)
    edges = use_hyperedges ? get_hyperedges(tn) : [tn[ind] for ind in bonds(tn)]

    # Create the line graph and dict mapping labels to hyperedges.
    g = qxg.LabeledGraph(length(edges))
    symbol_map = Dict{Symbol, Array{Symbol, 1}}(qxg.labels(g) .=> edges)

    # Connect the vertices of g if corresponding hyperedges overlap.
    for i in 1:length(edges)-1
        for j in i+1:length(edges)
            if !isempty(intersect(edges[i], edges[j]))
                qxg.add_edge!(g, i, j)
            end
        end
    end
    g, symbol_map
end

convert_to_line_graph(tnc::TensorNetworkCircuit; kwargs...) = convert_to_line_graph(tnc.tn; kwargs...)



# **************************************************************************************** #
#                              Contraction Planning
# **************************************************************************************** #

"""
    quickbb_contraction_plan(tn::TensorNetwork)

Use QuickBB to create a contraction plan for 'tn'.

# Keywords
- `time::Integer=120`: the number of second to run the quickbb binary for.
- `order::Symbol=:min_fill`: the branching order to be used by quickbb (:random or :min_fill).
- `hypergraph::Bool=false`: set if hyperedges exist in `tn` and should be accounted for.
"""
function quickbb_contraction_plan(tn::TensorNetwork;
                                  time::Integer=120,
                                  order::Symbol=:min_fill,
                                  hypergraph::Bool=false)
    # Convert tn to a line graph and pass it to quickbb to find an elimination order.
    lg, symbol_map = convert_to_line_graph(tn, use_hyperedges=hypergraph)
    if qxg.nv(lg) > 1 # quickkbb fails on trivial graphs with 1 vertex.
        order, metadata = qxg.quickbb(lg; time=time, order=order)
    elseif qxg.nv(lg) == 1
        order = qxg.labels(lg)
    end

    # Convert the elimination order to an array of Index structs in tn.
    order = [symbol_map[index_symbol] for index_symbol in order]

    # Convert the elimination order into a contraction plan.
    order_to_contraction_plan(order, tn)
end

quickbb_contraction_plan(tnc::TensorNetworkCircuit; kwargs...) = quickbb_contraction_plan(tnc.tn; kwargs...)

"""
    quickbb_contraction_plan(tensors::OrderedDict{Symbol, Array{Index, 1}})

Use QuickBB to create a contraction plan for the set of tensors described by the OrderedDict
`tensors`.

The keys of `tensors` are assumed to be ids/names of tensors and the values are arrays of
indices belonging to the corrseponding tensor.

# Keywords
- `time::Integer=0`: the number of second to run the quickbb binary for.
- `order::Symbol=:_`: the branching order to be used by quickbb (:random or :min_fill).
"""
function quickbb_contraction_plan(tensors::OrderedDict{Symbol, Array{Index, 1}};
                                  time::Integer=0,
                                  order::Symbol=:min_fill)
    # Create a graph for the tensors and convert it to to a line graph.
    tensor_ids = collect(keys(tensors))
    g = qxg.LabeledGraph(length(tensor_ids))
    for i = 1:length(tensor_ids)-1
        for j = i+1:length(tensor_ids)
            if !isempty(intersect(tensors[tensor_ids[i]], tensors[tensor_ids[j]]))
                qxg.add_edge!(g, i, j)
            end
        end
    end
    lg = qxg.line_graph(g)

    # Call quickbb to find a vertex elimination order for the line graph.
    order, metadata = qxg.quickbb(lg; time=time, order=order)

    # Convert the elimination order to an array of tensor symbol pairs.
    order = [parse.(Int, split(String(lg_vertex), '_')) for lg_vertex in order]
    order = [[tensor_ids[edge[1]], tensor_ids[edge[2]]] for edge in order]

    # Convert the elimination order into a contraction plan.
    order_to_contraction_plan(order, tensors)
end


"""
    flow_cutter_contraction_plan(tn::TensorNetwork;
                                 time::Integer=60,
                                 seed::Integer=-1,
                                 hypergraph::Bool=false)

Use flow cutter to create a contraction plan for 'tn'.

If flow cutter does not find a tree decomposotion to turn into a contraction plan within
the specified amount of time then the min fill heuristic is used to find a contraction plan.

# Keywords
- `time::Integer=60`: The number of seconds to run the flow cutter for.
- `seed::Integer=-1`: Sets the seed used by flow cutter. If negative then flow cutter will 
                      choose a seed.
- `hypergraph::Bool=false`: Sets if hyperedges in `tn` should be accounted for.
"""
function flow_cutter_contraction_plan(tn::TensorNetwork; 
                                      time::Integer=60,
                                      seed::Integer=-1,
                                      hypergraph::Bool=false)
    # Convert tn to a line graph and pass it to flow cutter to find an tree decomposition.
    lg, symbol_map = convert_to_line_graph(tn, use_hyperedges=hypergraph)

    # Use flow cutter to try find a tree decomposition of the line graph.
    td = qxg.flow_cutter(lg, time; seed=seed)

    # If a tree decomposition was found, convert it into a vertex elimination order for lg,
    # otherwise use the min fill heuristic to find an elimination order.
    if haskey(td, :treewidth)
        order = Symbol.(qxg.order_from_tree_decomposition(td))
    else
        tw, order = qxg.min_fill(lg)
    end

    # Convert the elimination order to an array of Index structs in tn.
    order = [symbol_map[lg_vertex] for lg_vertex in order]

    # Convert the elimination order into a contraction plan.
    order_to_contraction_plan(order, tn)
end

flow_cutter_contraction_plan(tnc::TensorNetworkCircuit; kwargs...) = flow_cutter_contraction_plan(tnc.tn; kwargs...)


"""
    min_fill_contraction_plan(tn::TensorNetwork;
                              hypergraph::Bool=false)

Use the min fill heuristic to create a contraction plan for 'tn'.

# Keywords
- `hypergraph::Bool=false`: set if hyperedges in `tn` should be accounted for.
"""
function min_fill_contraction_plan(tn::TensorNetwork;
                                   hypergraph::Bool=false)
    # Convert tn to a line graph and pass it to flow cutter to find an tree decomposition.
    lg, symbol_map = convert_to_line_graph(tn, use_hyperedges=hypergraph)

    # Use flow cutter to try find a tree decomposition of the line graph.
    tw, order = qxg.min_fill(lg)

    # Convert the elimination order to an array of Index structs in tn.
    order = [symbol_map[lg_vertex] for lg_vertex in order]

    # Convert the elimination order into a contraction plan.
    order_to_contraction_plan(order, tn)
end

min_fill_contraction_plan(tnc::TensorNetworkCircuit; kwargs...) = min_fill_contraction_plan(tnc.tn; kwargs...)


"""
    min_fill_contraction_plan(tensors::OrderedDict{Symbol, Array{Index, 1}})

Use the min fill heuristic to create a contraction plan for the set of tensors described by 
the OrderedDict `tensors`.

The keys of `tensors` are assumed to be ids/names of tensors and the values are arrays of
indices belonging to the corrseponding tensor.
"""
function min_fill_contraction_plan(tensors::OrderedDict{Symbol, Array{Index, 1}})
    # Create a graph for the tensors and convert it to to a line graph.
    tensor_ids = collect(keys(tensors))
    g = qxg.LabeledGraph(length(tensor_ids))
    for i = 1:length(tensor_ids)-1
        for j = i+1:length(tensor_ids)
            if !isempty(intersect(tensors[tensor_ids[i]], tensors[tensor_ids[j]]))
                qxg.add_edge!(g, i, j)
            end
        end
    end
    lg = qxg.line_graph(g)

    # Use the min fill heuristic to find a vertex elimination order for the line graph.
    tw, order = qxg.min_fill(lg)

    # Convert the elimination order to an array of tensor symbol pairs.
    order = [parse.(Int, split(String(lg_vertex), '_')) for lg_vertex in order]
    order = [[tensor_ids[edge[1]], tensor_ids[edge[2]]] for edge in order]

    # Convert the elimination order into a contraction plan.
    order_to_contraction_plan(order, tensors)
end


# **************************************************************************************** #
#                         Contraction with automatic slicing
# **************************************************************************************** #

"""
    contraction_scheme(tn::TensorNetwork, num::Integer)

Return an array of 'num' indices to slice in the tensor network 'tn' and a contraction plan
for the remaining tensor network.

# Keywords
- `time::Integer=120`: number of seconds to run quickbb when looking for contraction plan.
- `qbb_order::Symbol=:min_fill`: the branching order to be used by quickbb (:random or :min_fill).
- `lb::Bool=false`: set if a lowerbound for the treewidth should be computed.
- `score_function::Symbol=:direct_treewidth`: function to maximise when selecting vertices
                                            to remove. (:degree, :direct_treewidth)
- `hypergraph::Bool=true`: set if hyperedges exist in `tn` and should be accounted for.
"""
function contraction_scheme(tn::TensorNetwork, num::Integer;
                            time::Integer=120,
                            seed::Integer=-1,
                            score_function::Symbol=:direct_treewidth,
                            hypergraph::Bool=true)
    # Create the line graph for the given tn.
    lg, symbol_map = convert_to_line_graph(tn; use_hyperedges=hypergraph)

    # Use flow cutter to try find a tree decomposition of the line graph.
    td = qxg.flow_cutter(lg, time; seed=seed)
    flow_cutter_metadata = OrderedDict(1:length(td[:comments]) .=> td[:comments])

    # If a tree decomposition was found, convert it into a vertex elimination order for lg,
    # otherwise use the min fill heuristic to find an elimination order.
    if haskey(td, :treewidth)
        order = Symbol.(qxg.order_from_tree_decomposition(td))
        method_used = "flow cutter"
    else
        tw, order = qxg.min_fill(lg)
        method_used = "min fill heuristic"
    end

    # Create a dictionary for metadata regarding the contraction plan.
    contraction_metadata = OrderedDict{String, Any}()
    contraction_metadata["Method used"] = method_used
    contraction_metadata["Time allocated"] = time
    contraction_metadata["Seed used"] = seed
    contraction_metadata["Returned metadata"] = flow_cutter_metadata
    contraction_metadata["Hypergraph used"] = hypergraph
    contraction_metadata["Hyperedge contraction method"] = "Netcon where possible, min fill heuristic otherwise."

    # Use the greedy treewidth deletion algorithm to select indices in tn to slice.
    scheme = qxg.greedy_treewidth_deletion(lg, num;
                                           score_function=score_function,
                                           elim_order=order)
    sliced_lg, edges_to_slice, modified_orders, treewidths = scheme

    # Create a dictionary for metadata regarding sliced edges.
    slicing_metadata = OrderedDict{String, Any}()
    slicing_metadata["Method used"] = "greedy treewidth deletion"
    slicing_metadata["Edges sliced"] = num
    slicing_metadata["Score fucntion used"] = score_function
    slicing_metadata["Treewidths after slicing consecutive edges"] = treewidths

    # Convert the contraction plan to an array of Index structs before returning.
    contraction_order = !isempty(modified_orders) ? modified_orders[end] : order
    contraction_order = [symbol_map[edge_symbol] for edge_symbol in contraction_order]
    edges_to_slice = [symbol_map[edge_symbol] for edge_symbol in edges_to_slice]
    contraction_plan = order_to_contraction_plan(contraction_order, tn)

    # Convert edges_to_slice to an array of Index structs to slice.
    # TODO: I'm not confident there are no edge cases were this method for converting
    # edges to index arrays doesn't work. For example, when there are two consecutive
    # 2-qubit gates sharing both a hyperedge and a normal edge. Would be nice if it could
    # be made a bit cleaner also.
    indices_to_slice = Array{Array{<:Index, 1}, 1}(undef, length(edges_to_slice))
    for (i, edge) in enumerate(edges_to_slice)
        indices = Array{Index, 1}()
        for i in 1:length(edge)-1
            for j in i+1:length(edge)
                union!(indices, intersect(tn[edge[i]].indices, tn[edge[j]].indices))
            end
        end
        indices_to_slice[i] = indices
    end
    if length(indices_to_slice) > 0
        indices_to_slice = vcat(indices_to_slice...)
    else # in case there are no edges to slice ensure return array has correct type
        indices_to_slice = Index[]
    end

    metadata = OrderedDict{String, Any}()
    metadata["Determination of contraction plan"] = contraction_metadata
    metadata["Slicing"] = slicing_metadata

    indices_to_slice, contraction_plan, metadata
end

contraction_scheme(tnc::TensorNetworkCircuit, args...; kwargs...) = contraction_scheme(tnc.tn, args...; kwargs...)



# **************************************************************************************** #
#                                   Netcon functions
# **************************************************************************************** #
import TensorOperations.optimaltree

"""
    netcon(tn::TensorNetwork}

Call the TensorOperations netcon implementation on the tensors in `tensor_map`.

`tensor_map` is assumed to be a dictionary mapping tensor ids to arrays of indices connected
to the corresponding tensor.
"""
function netcon(tn::TensorNetwork, tensor_map::Union{OrderedDict{Symbol, <:Array{<:Index, 1}}, Nothing}=nothing)
    if tensor_map === nothing
        tensor_map = OrderedDict(keys(tn) .=> QXTns.inds.(values(tn)))
    end
    # Generate the arguments for the netcon method from the given network.
    tensor_indices, index_dims, tensor_labels = create_netcon_input(tensor_map)

    # Run the netcon method to find the optimal contraction plan.
    contraction_tree, cost = optimaltree(tensor_indices, index_dims)

    # Replace integers in contraction tree with node labels before returning
    plan, final_tensor_id = tree_to_contraction_plan(tn, contraction_tree, tensor_labels)
    plan
end

netcon(tnc::TensorNetworkCircuit) = netcon(tnc.tn)

"""
    create_netcon_input(tensor_map::OrderedDict{Symbol, <:Array{<:Index, 1}})

Create arguments for netcon to find the optimal contraction plan for the tensors
contatined in `tensor_map`.
"""
function create_netcon_input(tensor_map::OrderedDict{Symbol, <:Array{<:Index, 1}})
    # Collect the indices and corresponding dimensions of the network into arrays.
    num_tensors = length(tensor_map)

    # Allocate space for indices, dims and labels of tensors to contract.
    tensor_indices = Array{Array{Symbol, 1}, 1}(undef, num_tensors)
    tensor_dims = Array{Array{Int, 1}, 1}(undef, num_tensors)
    tensor_labels = Array{Symbol, 1}(undef, num_tensors)

    index_ids = Dict{Index, Symbol}()
    for (i, (tensor_id, tensor_inds)) in enumerate(tensor_map)
        tensor_indices[i] = [get!(index_ids, ind, Symbol("ind_$(i)_$(j)"))
                            for (j, ind) in enumerate(tensor_inds)]
        tensor_dims[i] = [QXTns.dim(ind) for ind in tensor_inds]
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
    order_to_contraction_plan(elimination_order::Array{<:Array{Symbol, 1}, 1},
                              tn::Union{TensorNetwork, OrderedDict{Symbol, Array{Index, 1}}}
                              )::Array{NTuple{3, Symbol}, 1}

Convert the given edge elimination order into a contraction plan for `tn`.

`tn` can be a TensorNetwork or an OrderedDict describing a set of tensors, in which case
the keys of `tn` are assumed to be tensor ids/names and the values are arrays containing the
indices of the corresponding tensor.
"""
function order_to_contraction_plan(elimination_order::Array{Array{Symbol, 1}, 1},
                                   tn::Union{TensorNetwork, OrderedDict{Symbol, Array{Index, 1}}}
                                   )::Array{NTuple{3, Symbol}, 1}
    # An array to hold the contraction plan.
    plan = Array{NTuple{3, Symbol}, 1}()

    # A dictionary to keep track of which tensors are replaced by which intermediate tensors
    # at different stages of the contraction process. Initially, before any pairwise
    # contractions, none of the tensors are replaced by intermediates, so all tensor ids are
    # mapped to themselves.
    intermediate_tensor_map = Dict{Symbol, Symbol}(keys(tn) .=> keys(tn))

    # Convert each edge in the elimination order to a set of pairwise contractions and
    # append it to plan.
    for (i, edge) in enumerate(elimination_order)
        if length(edge) == 2
            # For edges consisting of just two tensors, find the intermediates they belong
            # to and add their pairwise contraction to the plan.
            A_id = _get_intermediate_tensor(intermediate_tensor_map, edge[1])
            B_id = _get_intermediate_tensor(intermediate_tensor_map, edge[2])
            if !(A_id == B_id)
                # id for the new intermediate created by contracting A_id and B_id.
                I_id = Symbol("I$i")
                append!(plan, [(A_id, B_id, I_id)])

                # Update intermediate_tensor_map with new intermediate.
                intermediate_tensor_map[A_id] = I_id
                intermediate_tensor_map[B_id] = I_id
                intermediate_tensor_map[I_id] = I_id
            end

        elseif length(edge) > 2
            # For edges with more than 2 tensors, collect all of the tensors and
            # intermediate tensors that belong to the edge and find a contraction plan for
            # them. Append the contraction this plan to plan.
            tensors_to_contract = OrderedDict{Symbol, Array{Index, 1}}()
            for t_id in edge
                I_id = _get_intermediate_tensor(intermediate_tensor_map, t_id)
                inds = typeof(tn) <: TensorNetwork ? QXTns.inds(tn[t_id]) : tn[t_id]
                tensors_to_contract[I_id] = symdiff(get(tensors_to_contract, I_id, []), inds)
            end
            if length(tensors_to_contract) > 1
                # Check if netcon can be used on the given set of tensors. If not, use
                # a fallback method to find a contraction plan.
                tensor_sizes = prod.([QXTns.dim.(inds) for inds in values(tensors_to_contract)])
                if length(tensors_to_contract) < 37 && any(tensor_sizes .< 2^62)
                    local_contraction_plan = netcon(tn, tensors_to_contract)
                else
                    local_contraction_plan = min_fill_contraction_plan(tensors_to_contract)
                end
                append!(plan, local_contraction_plan)

                # Update intermediate_tensor_map with new intermediates.
                for (A_id, B_id, I_id) in local_contraction_plan
                    intermediate_tensor_map[A_id] = I_id
                    intermediate_tensor_map[B_id] = I_id
                    intermediate_tensor_map[I_id] = I_id
                end
            end
        end
    end
    plan
end

"""
    _get_intermediate_tensor(intermediate_tensor_map::Dict{Symbol, Symbol}, t_id::Symbol)::Symbol

Return the id of the intermediate tensor of which the tensor `t_id` is a factor.

During a tensor network contraction, tensors and intermediate tensors are contracted
together and replaced by the result. The dictionary `intermediate_tensor_map` is
assumed to describe which intermediate tensor has replaced a given tensor during a
contraction of a tensor network.
"""
function _get_intermediate_tensor(intermediate_tensor_map::Dict{Symbol, Symbol}, t_id::Symbol)::Symbol
    while true
        I_id = intermediate_tensor_map[t_id]
        # If t_id = I_id then t_id has not been replaced by an intermediate tensor yet.
        (t_id == I_id) && return I_id
        t_id = I_id
    end
end

"""
    tree_to_contraction_plan(tn::TensorNetwork, tree::Array{Any, 1}, labels::Array{Symbol, 1})::Tuple{Array{NTuple{3, Symbol}, 1}, Symbol}

Return the contraction plan described by the netcon contraction tree stored in 'tree'. The
tensors labeled by integers in `tree` are labeled by the corresponding symbol in `labels` in
the returned contraction plan.
"""
function tree_to_contraction_plan(tn::TensorNetwork, tree::Array{Any, 1}, labels::Array{Symbol, 1})::Tuple{Array{NTuple{3, Symbol}, 1}, Symbol}
    plan_A, A_id = tree_to_contraction_plan(tn, tree[1], labels)
    plan_B, B_id = tree_to_contraction_plan(tn, tree[2], labels)
    C_id = next_tensor_id!(tn)
    return [plan_A; plan_B; NTuple{3, Symbol}[(A_id, B_id, C_id)]], C_id
end

function tree_to_contraction_plan(::TensorNetwork, tree::Int, labels::Array{Symbol, 1})::Tuple{Array{NTuple{3, Symbol}, 1}, Symbol}
    [], labels[tree]
end
using AbstractTrees

using QXContexts

export build_compute_tree

"""
build_compute_tree(tnc::TensorNetworkCircuit,
                   plan::Vector{<:Tuple}),
                   bond_groups::Union{Nothing, Vector{<:Vector{<:Index}}}=nothing)

Function to build a compute tree from a tensor network and plan
"""
function build_compute_tree(tnc::TensorNetworkCircuit,
                            plan::Vector{<:Tuple},
                            bond_groups::Union{Nothing, Vector{<:Vector{<:Index}}}=nothing)
    tnc = copy(tnc)
    tn = tnc.tn
    # create load and outputs nodes
    nodes = Dict{Symbol, ComputeNode}()
    tc = TensorCache()

    # add a load node for each tensor
    for t in keys(tnc)
        data_symbol = push!(tc, tensor_data(tnc, t))
        op = LoadCommand(t, data_symbol, collect(size(tnc[t])))
        nodes[t] = ComputeNode{LoadCommand}(op)
    end

    # overwrite tensors corresponding to output tensors with output comamnds
    for (i, o) in enumerate(output_tensors(tnc))
        op = OutputCommand(o, i, size(tnc[o])[1])
        nodes[o] = ComputeNode{OutputCommand}(op)
    end

    changed_ids = Dict{Symbol, Symbol}()
    # add view commands to slice tensors
    if bond_groups !== nothing
        for (i, bg) in enumerate(bond_groups)
            related_tensors = union([tnc[b] for b in bg]...)
            slice_sym = Symbol("v$(i)")
            slice_dim = QXTns.dim(bg[1])
            for t in related_tensors
                new_sym = Symbol("$(t)_s")
                indices = hyperindices(tnc, t, all_indices=true)
                # find which group overlaps with slice_bonds
                index_position = findfirst(x -> length(intersect(bg, x)) > 0, indices)
                op = ViewCommand(new_sym, t, slice_sym, index_position, slice_dim)
                replace_tensor_symbol!(tn, t, new_sym)
                changed_ids[t] = new_sym
                node = ComputeNode{ViewCommand}(op)
                push!(node.children, nodes[t])
                nodes[t].parent = node
                nodes[new_sym] = node
            end
        end
    end

    # now add contraction nodes in
    for c in plan
        A_sym, B_sym, C_sym = c
        while haskey(changed_ids, A_sym) A_sym = changed_ids[A_sym] end
        while haskey(changed_ids, B_sym) B_sym = changed_ids[B_sym] end
        r = contraction_indices(tn, A_sym, B_sym)
        op = ContractCommand(C_sym, r.c_labels, A_sym, r.a_labels, B_sym, r.b_labels)
        node = ComputeNode{ContractCommand}(op)
        for s in [A_sym, B_sym]
            nodes[s].parent = node
            push!(node.children, nodes[s])
        end
        nodes[C_sym] = node
        contract_pair!(tn, A_sym, B_sym, C_sym; mock=true)
    end

    parentless = collect(keys(filter(x -> !isdefined(x[2], :parent), nodes)))
    reduce(parentless[2:end], init=parentless[1]) do x, y
        r = contraction_indices(tn, x, y)
        op = ContractCommand(:dummy, r.c_labels, x, r.a_labels, y, r.b_labels)
        # cmd = gen_ncon_command(tn, x, y, :dummy)
        sym = contract_pair!(tn, x, y, mock=true)
        op.output_name = sym
        node = ComputeNode{ContractCommand}(op)
        nodes[sym] = node
        for s in [x, y]
            push!(node.children, nodes[s])
            nodes[s].parent = node
        end
        sym # return sym to contract with next parentless node
    end
    parentless = collect(keys(filter(x -> !isdefined(x[2], :parent), nodes)))
    @assert length(parentless) == 1 "Only root node should have no parent"
    root = parentless[1]
    node = ComputeNode(SaveCommand(:output, root))
    push!(node.children, nodes[root])
    nodes[root].parent = node
    ComputeTree(node, convert(Dict, tc))
end
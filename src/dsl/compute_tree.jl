using AbstractTrees

using QXContexts

"""
    build_tree(tn::TensorNetwork, plan::Vector{<:Tuple})

Function to build a compute tree from a tensor network and plan
"""
function QXContexts.build_tree(tn::TensorNetwork, plan::Vector{<:Tuple})
    tn = copy(tn)
    cmd_map = Dict(x[3] => (x[1] => x[2]) for x in plan)
    all_outputs = Set(keys(cmd_map))
    all_inputs = vcat(collect.(values(cmd_map))...)
    parentless = collect(setdiff(all_outputs, all_inputs))
    # @assert length(parentless) == 1 "Plan leads to disconnected graph"
    # if length(parentless) > 1
        # for i = 1:length(parentless)-1
            # cmd_map[]
        # end
    # end
    tree_nodes = Dict{Symbol, ComputeNode}()
    for sym in parentless
        tree_nodes[sym] = _create_node(tn, cmd_map, sym)
    end
    reduce(parentless[2:end], init=parentless[1]) do x, y
        cmd = gen_ncon_command(tn, x, y, :dummy)
        sym = contract_pair!(tn, x, y, mock=true)
        cmd.output_name = sym
        node = ComputeNode{ContractCommand}(cmd)
        tree_nodes[sym] = node
        if haskey(tree_nodes, x)
            node.left = tree_nodes[x]
            node.left.parent = node
            delete!(tree_nodes, x)
        end
        if haskey(tree_nodes, y)
            node.right = tree_nodes[y]
            node.right.parent = node
            delete!(tree_nodes, y)
        end
        sym
    end
    first(values(tree_nodes))
end

QXContexts.build_tree(tnc::TensorNetworkCircuit, plan::Vector{<:Tuple}) = build_tree(tnc.tn, plan)

function _create_node(tn::TensorNetwork, cmd_map, output_sym)
    inputs = cmd_map[output_sym]
    node = ComputeNode{ContractCommand}()
    if haskey(cmd_map, inputs[1])
        node.left = _create_node(tn, cmd_map, inputs[1])
        node.left.parent = node
    end
    if haskey(cmd_map, inputs[2])
        node.right = _create_node(tn, cmd_map, inputs[2])
        node.right.parent = node
    end
    cmd = gen_ncon_command(tn, inputs..., output_sym)
    contract_pair!(tn, inputs..., output_sym, mock=true)
    node.data = cmd
    node
end
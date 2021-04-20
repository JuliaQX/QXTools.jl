using AbstractTrees

# In this file we define a tree data structure which can be used for optimisation passes
# over contraction commands

"""Generic binary tree node data structure"""
mutable struct BinaryNode{T}
    data::T
    parent::BinaryNode{T}
    left::BinaryNode{T}
    right::BinaryNode{T}

    # Root constructor
    BinaryNode{T}(data) where T = new{T}(data)
    # Child node constructor
    BinaryNode{T}(data, parent::BinaryNode{T}) where T = new{T}(data, parent)
end
BinaryNode(data) = BinaryNode{typeof(data)}(data)

"""
    leftchild(data, parent::BinaryNode)

Returns left child node or nothing
"""
function leftchild(data, parent::BinaryNode)
    !isdefined(parent, :left) || error("left child is already assigned")
    node = typeof(parent)(data, parent)
    parent.left = node
end

"""
    rightchild(data, parent::BinaryNode)

Returns right child node or nothing
"""
function rightchild(data, parent::BinaryNode)
    !isdefined(parent, :right) || error("right child is already assigned")
    node = typeof(parent)(data, parent)
    parent.right = node
end

"""
    AbstractTrees.children(node::BinaryNode)

Implement children function from AbstractTrees package
"""
function AbstractTrees.children(node::BinaryNode)
    if isdefined(node, :left)
        if isdefined(node, :right)
            return (node.left, node.right)
        end
        return (node.left,)
    end
    isdefined(node, :right) && return (node.right,)
    return ()
end

## functions and data structures for trees of contraction commands

"""Structure to represent a contraction command"""
struct ContractCommand
    output_name::Symbol
    output_idxs::Vector{Int}
    left_name::Symbol
    left_idxs::Vector{Int}
    right_name::Symbol
    right_idxs::Vector{Int}
end

"""
    Base.write(io::IO, cmd::ContractCommand)

Function to serialise command to the given IO stream
"""
function Base.write(io::IO, cmd::ContractCommand)
    j = x -> length(x) == 0 ? "0" : join(x, ",")
    write(io, "ncon $(cmd.output_name) $(j(cmd.output_idxs)) $(cmd.left_name) $(j(cmd.left_idxs)) $(cmd.right_name) $(j(cmd.right_idxs))\n")
end

"""
    Base.write(io::IO, tree::tree::BinaryNode{ContractCommand})

Function to serialise tree to the given IO stream
"""
function Base.write(io::IO, tree::BinaryNode{ContractCommand})
    for cmd in PostOrderDFS(tree)
        write(io, cmd.data)
    end
end

"""
    build_tree(cmds::Vector{ContractCommand})

Function to construct a tree from a list of contraction commands
"""
function build_tree(cmds::Vector{ContractCommand})
    cmd_map = Dict(x.output_name => deepcopy(x) for x in cmds)
    all_outputs = Set(keys(cmd_map))
    get_names = x -> [x.left_name, x.right_name]
    all_inputs = vcat(get_names.(values(cmd_map))...)
    root_sym = collect(setdiff(all_outputs, all_inputs))[1]
    root_node = _create_node(cmd_map, root_sym)
    root_node
end

function _create_node(cmd_map, cmd_sym)
    cmd = cmd_map[cmd_sym]
    node = BinaryNode{ContractCommand}(cmd)
    if haskey(cmd_map, cmd.left_name)
        node.left = _create_node(cmd_map, cmd.left_name)
        node.left.parent = node
    end
    if haskey(cmd_map, cmd.right_name)
        node.right = _create_node(cmd_map, cmd.right_name)
        node.right.parent = node
    end
    node
end

AbstractTrees.printnode(io::IO, node::BinaryNode{ContractCommand}) = print(io, node.data.output_name)

"""
    remove_repeated!(indices::Vector{Int64})

Given a list of integers, only keep first occurence of each
"""
function remove_repeated!(indices::Vector{Int64})
    to_remove = Int64[]
    for i in 1:length(indices)
        if findfirst(x -> x == indices[i], indices) != i
            push!(to_remove, i)
        end
    end
    for i in sort(to_remove, rev=true) popat!(indices, i) end
    indices
end

function _remove_repeated_in_output!(cmd::ContractCommand, indices::Vector{Int64})
    # first we create a map of all indices and start with them mapping to themselves
    all_idx = Set([cmd.output_idxs..., cmd.left_idxs..., cmd.right_idxs...])
    all_idx_map = Dict(x => x for x in all_idx)

    # go through indices and update index map with first occurence
    for i in 1:length(indices)
        first_pos = findfirst(x -> x == indices[i], indices)
        all_idx_map[cmd.output_idxs[i]] = cmd.output_idxs[first_pos]
    end

    # replace indices with those they map to in the map
    replace!(cmd.left_idxs, all_idx_map...)
    replace!(cmd.right_idxs, all_idx_map...)

    # we can now remove duplicate indices from output_idxs
    filter!(x -> all_idx_map[x] == x, cmd.output_idxs)
    nothing
end

"""
    remove_repeated!(node::BinaryNode{ContractCommand}, indices=nothing)

Descend contraction tree simplifying contraction commands by removing repeated indices.
Repeated indices occur where two indices of a tensor are hyper indices of another tensor.
The process of removing these given a node consists of:
1. Check if there are repeated indices in lists of left or right indices
2. If there are then
"""
function remove_repeated!(node::BinaryNode{ContractCommand}, indices=nothing)
    if indices !== nothing
        _remove_repeated_in_output!(node.data, indices)
    end
    is_hyper = x -> length(x) > length(Set(x))
    if is_hyper(node.data.left_idxs) && isdefined(node, :left)
        remove_repeated!(node.left, node.data.left_idxs)
        remove_repeated!(node.data.left_idxs)
    end
    if is_hyper(node.data.right_idxs) && isdefined(node, :right)
        remove_repeated!(node.right, node.data.right_idxs)
        remove_repeated!(node.data.right_idxs)
    end
end
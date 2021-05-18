using AbstractTrees

# In this file we define a tree data structure which can be used for optimisation passes
# over contraction commands

"""Generic binary tree node data structure"""
mutable struct ComputeNode{T}
    data::T
    parent::ComputeNode{T}
    left::ComputeNode{T}
    right::ComputeNode{T}
    # Root constructor
    ComputeNode{T}(data) where T = new{T}(data)
    # Child node constructor
    ComputeNode{T}(data, parent::ComputeNode{T}) where T = new{T}(data, parent)
end
ComputeNode(data) = ComputeNode{typeof(data)}(data)

"""
    AbstractTrees.children(node::ComputeNode)

Implement children function from AbstractTrees package
"""
function AbstractTrees.children(node::ComputeNode)
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
    reshape_groups::Vector{Int}
    function ContractCommand(output_name::Symbol,
                             output_idxs::Vector{Int},
                             left_name::Symbol,
                             left_idxs::Vector{Int},
                             right_name::Symbol,
                             right_idxs::Vector{Int})
        new(output_name, output_idxs, left_name, left_idxs, right_name, right_idxs, ones(Int, length(output_idxs)))
    end
end

function ContractCommand(cmd::AbstractString)
    @assert match(r"^ncon", cmd) !== nothing "Command must begin with \"ncon\""
    p = split(cmd, " ")[2:end]
    s = x -> Symbol(x)
    l = x -> x == "0" ? Int[] : map(y -> parse(Int, y), split(x, ","))
    ContractCommand(s(p[1]), l(p[2]), s(p[3]), l(p[4]), s(p[5]), l(p[6]))
end


"""
    Base.write(io::IO, cmd::ContractCommand)

Function to serialise command to the given IO stream
"""
function Base.write(io::IO, cmd::ContractCommand)
    j = x -> length(x) == 0 ? "0" : join(x, ",")
    write(io, "ncon $(cmd.output_name) $(j(cmd.output_idxs)) $(cmd.left_name) $(j(cmd.left_idxs)) $(cmd.right_name) $(j(cmd.right_idxs))\n")
    if length(cmd.reshape_groups) < length(cmd.output_idxs)
        write(io, "reshape $(cmd.output_name) $(j(cmd.reshape_groups))\n")
    end
end

"""
    Base.write(io::IO, tree::tree::ComputeNode{ContractCommand})

Function to serialise tree to the given IO stream
"""
function Base.write(io::IO, tree::ComputeNode{ContractCommand})
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
    node = ComputeNode{ContractCommand}(cmd)
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

AbstractTrees.printnode(io::IO, node::ComputeNode{ContractCommand}) = print(io, node.data.output_name)

"""
    _remove_positions_in_reshape!(positions::Array{Int}, reshape_groups::Vector{Int})

Remove the given positions from the list of groups

Example:
reshape_groups: [1,3,1]
positions: [1,3]
output: [2,1]
"""
function _remove_positions_in_reshape(positions::Array{Int}, reshape_groups::Vector{Int})
    expanded_positions = vcat([fill(i, x) for (i, x) in enumerate(reshape_groups)]...)
    expanded_positions = expanded_positions[setdiff(1:length(expanded_positions), positions)]
    final_groups = Int[]
    if length(expanded_positions) == 0 return final_groups end

    val = expanded_positions[1]
    count = 1
    for i in 2:length(expanded_positions)
        if expanded_positions[i] != val
            push!(final_groups, count)
            count = 1
            val = expanded_positions[i]
        else
            count += 1
        end
    end
    push!(final_groups, count)
    final_groups
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

    # we can now remove duplicate indices from output_idxs and also update reshape indices to be consistent
    # find positions of values that will be remoeved
    remove_positions = findall(x -> all_idx_map[x] == x, cmd.output_idxs)
    # remove repeated indices in output_indices
    filter!(x -> all_idx_map[x] == x, cmd.output_idxs)
    cmd.reshape_groups = _remove_positions_in_reshape(remove_positions, cmd.reshape_groups)
    nothing
end

"""
    remove_repeated!(node::ComputeNode{ContractCommand}, indices=nothing)

Descend contraction tree simplifying contraction commands by removing repeated indices.
Repeated indices occur where two indices of a tensor are hyper indices of another tensor.
The process of removing these given a node consists of:

1. Check if there are repeated indices in lists of left or right indices
2. If there are repeated indices for either set of indices, remove any after the first
occurence and call this function on the relevant child passing the original index list

For example, when the node has the format
````
ncon d 1 c 2,2,3 e 1,3
```
where `c` is created with
```
ncon c 1,3,4 a 1,2 b 2,3,4
```
we first update initial command to
```
ncon d 1 c 2,3 e 1,3
```
and then update the command which gives c to
```
ncon c 1,4 a 1,2 b 2,1,4
```
In this way we can decreast the rank of the largest tensors
"""
function remove_repeated!(node::ComputeNode{ContractCommand}, indices=nothing)
    if indices !== nothing
        _remove_repeated_in_output!(node.data, indices)
    end
    is_hyper = x -> length(x) > length(Set(x))
    if is_hyper(node.data.left_idxs) && isdefined(node, :left)
        remove_repeated!(node.left, node.data.left_idxs)
        unique!(node.data.left_idxs)
    end
    if is_hyper(node.data.right_idxs) && isdefined(node, :right)
        remove_repeated!(node.right, node.data.right_idxs)
        unique!(node.data.right_idxs)
    end
end

"""
    permute_and_merge!(node::ComputeNode{ContractCommand}, new_index_order=nothing)

Descend contraction tree simplifying contraction commands by permuting and merging together
indices so that contractions operations are more efficient.

For example, when there are a sequence of contractions like
````
ncon c 1,2,3,4 a 1,2,3,4 b 0
ncon e 1,2,3,4 d 1,2,3,4 e 0
ncon f 1,4,5,6 c 1,2,3,4 e 5,6,3,2
```
in the last contraction indices 2,3 are contracted and thus can be merged. We first permute as
```
ncon d 1,4,5,6 c 2,3,1,4 e 2,3,5,6
```
and then joining (2,3) => 2 which gives
```
ncon ncon d 1,4,5,6 c 2,1,4 e 2,5,6
```
This would require previous commands to be changed so full set would be
````
ncon c (2,3),1,4 a 1,2,3,4 b 0
ncon e (4,3),1,2 d 1,2,3,4 e 0
ncon f 1,4,5,6 c 2,1,4 e 2,5,6
```
where the parentheses indicate that these indices should be merged.

Additionally, any indices appearing in all three tuples should be moved to the end.
In this case we would have a recursive process which passes down indices list
"""
function permute_and_merge!(node::ComputeNode{ContractCommand}, new_index_order=nothing)
    _permute_and_merge!(node.data, new_index_order)

    # identify batched, common and remaining indices
    batch_idxs = batched_indices(node.data) # these go at the end
    common_idxs = setdiff(intersect(node.data.left_idxs, node.data.right_idxs), batch_idxs) # these go at the start in order appear in left
    common_idxs = filter(x -> x in common_idxs, node.data.left_idxs) # ensure it's ordered by left idxs
    remaining_idxs = setdiff(union(node.data.left_idxs, node.data.right_idxs), union(common_idxs, batch_idxs)) # these go in the middle
    remaining_idxs = filter(x -> x in remaining_idxs, node.data.output_idxs) # reorder by output

    if isdefined(node, :left)
        l_index_map = Dict(x => i for (i, x) in enumerate(node.data.left_idxs))
        m = x -> [l_index_map[y] for y in x]
        left_remaining_idxs = filter(x -> x in node.data.left_idxs, remaining_idxs)
        # left_remaining_idx_groups = find_overlaps(node.data.output_idxs, left_remaining_idxs)
        new_index_order = length(common_idxs) == 0 ? Vector{Vector{Int}}() : [m(common_idxs)]
        append!(new_index_order, [map(x -> [x], m(left_remaining_idxs))..., map(x -> [x], m(batch_idxs))...])
        # append!(new_index_order, [m.(left_remaining_idx_groups)..., map(x -> [x], m(batch_idxs))...])
        permute_and_merge!(node.left, new_index_order)
        new_left_idxs = length(common_idxs) > 0 ? [common_idxs[1]] : Int[]
        append!(new_left_idxs, [left_remaining_idxs..., batch_idxs...])
        empty!(node.data.left_idxs)
        append!(node.data.left_idxs, new_left_idxs)
    end
    if isdefined(node, :right)
        r_index_map = Dict(x => i for (i, x) in enumerate(node.data.right_idxs))
        m = x -> [r_index_map[y] for y in x]
        right_remaining_idxs = filter(x -> x in node.data.right_idxs, remaining_idxs)
        # right_remaining_idx_groups = find_overlaps(node.data.output_idxs, right_remaining_idxs)
        new_index_order = length(common_idxs) == 0 ? Vector{Vector{Int}}() : [m(common_idxs)]
        append!(new_index_order, [map(x -> [x], m(right_remaining_idxs))..., map(x -> [x], m(batch_idxs))...])
        # append!(new_index_order, [m.(right_remaining_idx_groups)..., map(x -> [x], m(batch_idxs))...])
        permute_and_merge!(node.right, new_index_order)
        new_right_idxs = length(common_idxs) > 0 ? [common_idxs[1]] : Int[]
        append!(new_right_idxs, [right_remaining_idxs..., batch_idxs...])
        empty!(node.data.right_idxs)
        append!(node.data.right_idxs, new_right_idxs)
    end
end

batched_indices(cmd) = intersect(cmd.output_idxs, cmd.left_idxs, cmd.right_idxs)

function _permute_and_merge!(cmd::ContractCommand, new_index_order)
    @show "Updating $(cmd.output_name) with $new_index_order"
    if new_index_order !== nothing
        # get mapping between indices as seen by next command, ordered 1,2..n where n is rank
        # and grouped according to reshape_groups
        current_idxs = Dict{Int, Vector{Int}}()
        pos = 1
        for (i, group_size) in enumerate(cmd.reshape_groups)
            current_idxs[i] = cmd.output_idxs[pos:pos+group_size-1]
            pos += group_size
        end
        flat_order = vcat(new_index_order...)
        cmd.output_idxs[:] = vcat(map(x -> current_idxs[x], flat_order)...)
        # update reshape groups to reflect any new merges
        new_reshape_groups = []
        pos = 1
        for l in length.(new_index_order)
            push!(new_reshape_groups, sum(cmd.reshape_groups[pos:pos+l-1]))
            pos += l
        end
        empty!(cmd.reshape_groups)
        append!(cmd.reshape_groups, new_reshape_groups)
    end
end

function find_overlaps(array::Vector{T}, sub_array::Vector{T}) where T
    sub_pos = 1
    groups = Vector{Vector{T}}()
    while sub_pos <= length(sub_array)
        pos = findfirst(x -> x == sub_array[sub_pos], array)
        i = 1
        while pos + i <= length(array) &&
              sub_pos + i <= length(sub_array) &&
              array[pos+i] == sub_array[sub_pos+i]
            i += 1
        end
        push!(groups, sub_array[sub_pos:sub_pos+i-1])
        pos += i
        sub_pos += i
    end
    groups
end
"""Structure to represent a contraction command"""
mutable struct ContractCommand
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
        new(output_name,
            output_idxs,
            left_name,
            left_idxs,
            right_name,
            right_idxs,
            ones(Int, length(output_idxs)))
    end

    function ContractCommand(output_name::Symbol,
                             output_idxs::Vector{Int},
                             left_name::Symbol,
                             left_idxs::Vector{Int},
                             right_name::Symbol,
                             right_idxs::Vector{Int},
                             reshape_groups::Vector{Int})
        new(output_name,
        output_idxs,
        left_name,
        left_idxs,
        right_name,
        right_idxs,
        reshape_groups)
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
        groups = Vector{Vector{Int}}()
        pos = 1
        for j in cmd.reshape_groups
            push!(groups, collect(pos:pos+j-1))
            pos += j
        end
        groups_str = join(map(x -> join(x, ","), groups), ";")
        write(io, "reshape $(cmd.output_name) $(groups_str)\n")
    end
end

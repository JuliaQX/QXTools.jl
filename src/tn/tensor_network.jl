using DataStructures
using ITensors

# TensorNetwork struct and public functions
export next_tensor_id
export TensorNetwork, bonds, simple_contraction, tensor_data, neighbours

const qxsim_ids = Dict{Symbol, Int64}(:tensor_id => 0)
@noinline next_tensor_id() = begin qxsim_ids[:tensor_id] += 1; Symbol("t$(qxsim_ids[:tensor_id])") end

"""Tensor network data-structure"""
mutable struct TensorNetwork
    tensor_map::OrderedDict{Symbol, ITensor}
    bond_map::OrderedDict{Index, Vector{Symbol}}

    # constructors
    TensorNetwork(tensor_map::OrderedDict{Symbol, ITensor},
                  bond_map::OrderedDict{<:Index, Vector{Symbol}}) = new(tensor_map, bond_map)
    TensorNetwork() = new(OrderedDict{Symbol, ITensor}(), OrderedDict{Index, Vector{Symbol}}())
end

"""
    TensorNetwork(array::Vector{<: ITensor})

Outer consturctor to create a tensor network object from an array of
ITensor objects
"""
function TensorNetwork(array::Vector{<: ITensor})
    tensor_map = OrderedDict{Symbol, ITensor}()
    bond_map = OrderedDict{Index, Vector{Symbol}}()
    for (i, tensor) in enumerate(array)
        tensor_id = next_tensor_id()
        tensor_map[tensor_id] = tensor
        for bond in inds(tensor)
            if haskey(bond_map, bond)
                push!(bond_map[bond], tensor_id)
            else
                bond_map[bond] = [tensor_id]
            end
        end
    end
    TensorNetwork(tensor_map, bond_map)
end

Base.length(tn::TensorNetwork) = length(tn.tensor_map)
tensor_data(tn::TensorNetwork, i::Symbol) = tensor_data(tn.tensor_map[i])
Base.values(tn::TensorNetwork) = values(tn.tensor_map)
Base.iterate(tn::TensorNetwork) = iterate(values(tn))
Base.iterate(tn::TensorNetwork, state) = iterate(values(tn), state)
Base.eltype(tn::TensorNetwork) = ITensor
Base.keys(tn::TensorNetwork) = keys(tn.tensor_map)
Base.getindex(tn::TensorNetwork, i::Symbol) = tn.tensor_map[i]
Base.getindex(tn::TensorNetwork, i::T) where T <: Index = tn.bond_map[i]
bonds(tn::TensorNetwork) = keys(tn.bond_map)

"""
    neighbours(tn::TensorNetwork, tensor::Symbol)

Function get the symbols of the neighbouring tensors
"""
function neighbours(tn::TensorNetwork, tensor::Symbol)
    tensor_indices = inds(tn[tensor])
    connected_tensors = unique(vcat([tn[x] for x in tensor_indices]...))
    setdiff(connected_tensors, [tensor])
end

"""
    merge(a::TensorNetwork, b::TensorNetwork)

Join two networks together
"""
function Base.merge(a::TensorNetwork, b::TensorNetwork)
    tensor_map = merge(a.tensor_map, b.tensor_map)
    bond_map = OrderedDict{Index, Vector{Symbol}}()
    for (tid, tensor) in pairs(tensor_map)
        for bond in inds(tensor)
            if !haskey(bond_map, bond)
                bond_map[bond] = [tid]
            else
                push!(bond_map[bond], tid)
            end
        end
    end
    c = TensorNetwork(tensor_map, bond_map)

end

"""
    tensor_data(tensor::ITensor)

Get the data associated with given tensor
"""
function tensor_data(tensor::ITensor)
    reshape(convert(Array, store(tensor)), Tuple([dim(x) for x in inds(tensor)]))
end

"""
    push!(tn::TensorNetwork,
          indices::Vector{Index},
          data::Array{T, N}) where {T, N}

Function to add a tensor to the tensor network tensor data and indices
"""
function Base.push!(tn::TensorNetwork,
                    indices::Vector{<:Index},
                    data::Array{T, N}) where {T, N}
    @assert size(data) == Tuple(dim.(indices))
    data_store = NDTensors.Dense(reshape(data, prod(size(data))))
    tensor = ITensor(data_store, indices)
    tid = next_tensor_id()
    tn.tensor_map[tid] = tensor
    for bond in indices
        if haskey(tn.bond_map, bond)
            push!(tn.bond_map[bond], tid)
        else
            tn.bond_map[bond] = [tid]
        end
    end
    tid
end

"""
    push!(tn::TensorNetwork,
          tensor::ITensor{N}) where {N}

Function to add a tensor
"""
function Base.push!(tn::TensorNetwork,
               tensor::ITensor{N}) where {N}
    @assert size(data) == Tuple(dim.(indices))
    data_store = NDTensors.Dense(reshape(data, prod(size(data))))
    tensor = ITensor(data_store, indices)
    tid = next_tensor_id()
    tn.tensor_map[tid] = tensor
    for bond in indices
        if haskey(tn.bond_map, bond)
            push!(tn.bond_map[bond], tid)
        else
            tn.bond_map[bond] = [tid]
        end
    end
    tid
end


"""
    simple_contraction(tn::TensorNetwork)

Function to perfrom a simple contraction, contracting all tensors in order.
Only useful for very small networks for testing.
"""
function simple_contraction(tn::TensorNetwork)
    reduce(*, tn, init=ITensor(1.))
end
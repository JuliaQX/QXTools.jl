using DataStructures
using ITensors

# TensorNetwork struct and public functions
export next_tensor_id
export TensorNetwork, bonds, simple_contraction, simple_contraction!, tensor_data, neighbours
export decompose_tensor!, replace_with_svd!
export contract_tn!, contract_pair!, replace_tensor_symbol!, contract_ncon_indices
export use_mock_tensors

const qxsim_ids = Dict{Symbol, Int64}(:tensor_id => 0)
@noinline next_tensor_id() = begin qxsim_ids[:tensor_id] += 1; Symbol("t$(qxsim_ids[:tensor_id])") end

"""Tensor network data-structure"""
mutable struct TensorNetwork
    tensor_map::OrderedDict{Symbol, ITensor}
    bond_map::OrderedDict{Index, Vector{Symbol}}
end

# constructors
TensorNetwork() = TensorNetwork(OrderedDict{Symbol, ITensor}(), OrderedDict{Index, Vector{Symbol}}())

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

Base.copy(tn::TensorNetwork) = TensorNetwork(collect(values(tn)))
Base.length(tn::TensorNetwork) = length(tn.tensor_map)
Base.values(tn::TensorNetwork) = values(tn.tensor_map)
Base.iterate(tn::TensorNetwork) = iterate(values(tn))
Base.iterate(tn::TensorNetwork, state) = iterate(values(tn), state)
Base.eltype(tn::TensorNetwork) = ITensor
Base.keys(tn::TensorNetwork) = keys(tn.tensor_map)
Base.getindex(tn::TensorNetwork, i::Symbol) = tn.tensor_map[i]
Base.haskey(tn::TensorNetwork, i::Symbol) = haskey(tn.tensor_map, i)
Base.getindex(tn::TensorNetwork, i::T) where T <: Index = tn.bond_map[i]
Base.haskey(tn::TensorNetwork, i::T) where T <: Index = haskey(tn.bond_map, i)
Base.show(io::IO, ::MIME"text/plain", tn::TensorNetwork) = print(io, "TensorNetwork(tensors => $(length(tn)), bonds => $(length(bonds(tn))))")

bonds(tn::TensorNetwork) = keys(tn.bond_map)
tensor_data(tn::TensorNetwork, i::Symbol) = tensor_data(tn.tensor_map[i])

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
    tid = next_tensor_id()
    tn.tensor_map[tid] = tensor
    for bond in inds(tensor)
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
    store(reduce(contract_tensors, tn, init=ITensor(1.)))
end

"""
    simple_contraction!(tn::TensorNetwork)

Function to perfrom a simple contraction, contracting all tensors in order.
Only useful for very small networks for testing.
"""
function simple_contraction!(tn::TensorNetwork)
    tensor_syms = collect(keys(tn))
    A = tensor_syms[1]
    for B in tensor_syms[2:end]
        A = contract_pair!(tn, A, B)
    end
    store(tn[A])
end

"""
    contract_pair!(tn::TensorNetwork, A_id::Symbol, B_id::Symbol)

Contract the tensors in 'tn' with ids 'A_id' and 'B_id'.
"""
function contract_pair!(tn::TensorNetwork, A_id::Symbol, B_id::Symbol)
    # Get and contract the tensors A and B to create tensor C.
    A = tn.tensor_map[A_id]
    B = tn.tensor_map[B_id]
    C_id = next_tensor_id()
    C = contract_tensors(A, B)

    # Remove the contracted indices from the bond map in tn. Also, replace all references
    # in tn to tensors A and B with a reference to tensor C.
    for ind in commoninds(A, B)
        delete!(tn.bond_map, ind)
    end
    for ind in noncommoninds(A, B)
        tn.bond_map[ind] = replace(tn.bond_map[ind], A_id=>C_id, B_id=>C_id)
    end

    # Add tensor C to the tn and remove both A and B.
    tn.tensor_map[C_id] = C
    delete!(tn.tensor_map, A_id); delete!(tn.tensor_map, B_id)
    C_id
end


"""
    contract_tn!(tn::TensorNetwork, plan)

Contract the indices of 'tn' according to 'plan'.
Contract each index in the contraction plan while skipping indices that are not shared
by exactly two tensors in tn. (i.e. open indices and hyper edges are skipped.)
"""
function contract_tn!(tn::TensorNetwork, plan::Array{<:Index, 1})
    for index in plan
        if haskey(tn, index)
            tensor_pair = tn[index]
            if length(tensor_pair) == 2
                contract_pair!(tn, tensor_pair...)
            end
        end
    end
    simple_contraction!(tn)
end

"""
    decompose_tensor!(tn::TensorNetwork,
                      tensor_id::Symbol,
                      left_indices::Array{<:Index, 1};
                      contract_S_with::Symbol=:V,
                      kwargs...)

Function to decompose a tensor in a tensor network using svd.

# Keywords
- `contract_S_with::Symbol=:V`: the maxtrix which should absorb the matrix of singular values
- `maxdim::Int`: the maximum number of singular values to keep.
- `mindim::Int`: the minimum number of singular values to keep.
- `cutoff::Float64`: set the desired truncation error of the SVD.
"""
function decompose_tensor!(tn::TensorNetwork,
                           tensor_id::Symbol,
                           left_indices::Array{<:Index, 1};
                           contract_S_with::Symbol=:V,
                           kwargs...)
    # Decompose the tensor into its svd factors.
    U_id, S_id, V_id = replace_with_svd!(tn, tensor_id, left_indices; kwargs)

    # Absorb the singular values tensor into eith U or V.
    if contract_S_with == :V_id
        return U_id, contract_pair!(tn, S_id, V_id)
    else
        return contract_pair!(tn, S_id, U_id), V_id
    end
end


"""
    replace_with_svd!(tn::TensorNetwork, 
                      tensor_id::Symbol,
                      left_indices::Array{<:Index, 1};
                      kwargs...)

Function to replace a tensor in a tensor network with its svd.

The indices contained in 'left_indices' are considered the row indices of the tensor when
the svd is performed.

# Keywords
- `maxdim::Int`: the maximum number of singular values to keep.
- `mindim::Int`: the minimum number of singular values to keep.
- `cutoff::Float64`: set the desired truncation error of the SVD.
"""
function replace_with_svd!(tn::TensorNetwork, 
                            tensor_id::Symbol,
                            left_indices::Array{<:Index, 1};
                            kwargs...)
    # Get the tensor and decompose it.
    tensor = tn.tensor_map[tensor_id]
    U, S, V = svd(tensor, left_indices; use_absolute_cutoff=true, kwargs...)

    # Remove the original tensor and add its svd factors to the network.
    delete!(tn, tensor_id)
    U_id = push!(tn, U)
    S_id = push!(tn, S)
    V_id = push!(tn, V)
    U_id, S_id, V_id
end


"""
    delete!(tn::TensorNetwork, tensor_id::Symbol)

Function to remove a tensor from a tensor network.
"""
function Base.delete!(tn::TensorNetwork, tensor_id::Symbol)
    tensor = tn.tensor_map[tensor_id]
    for index in inds(tensor)
        # remove tensor_id from tn.bond_map[index]
        filter!(id -> id ≠ tensor_id, tn.bond_map[index])        
    end
    delete!(tn.tensor_map, tensor_id)
end

"""
    contract_ncon_indices(tn::TensorNetwork, A_sym::Symbol, B_sym)

Function return indices in ncon format for contraction of tensors with given symbols.
Returns two tuples for indices in each with convention that negative values are remaining
indices and positive values are indices being contracted over.

For example if (1, -1), (-2, 1) is returned, this menas that the first index of tensor A
A is contracted with the second index of  tensor B and the resulting tensor will have
indices corresponding to the second index of the first tensor and first index of the second
tensor.
"""
function contract_ncon_indices(tn::TensorNetwork, A_sym::Symbol, B_sym::Symbol)
    _contract_ncon_indices(inds(tn[A_sym]), inds(tn[B_sym]))
end

"""
    _contract_ncon_indices(A_inds::IndexSet{M}, B_inds::IndexSet{N}) where {M, N}

Function return indices in ncon format for contraction of tensors with given index sets.
Returns two tuples for indices in each with convention that negative values are remaining
indices and positive values are indices being contracted over.

For example if (1, -1), (-2, 1) is returned, this menas that the first index of tensor A
A is contracted with the second index of  tensor B and the resulting tensor will have
indices corresponding to the second index of the first tensor and first index of the second
tensor.
"""
function _contract_ncon_indices(A_inds::IndexSet{M}, B_inds::IndexSet{N}) where {M, N}
    labels = ITensors.compute_contraction_labels(A_inds, B_inds)
    # ITensors uses a different convention with negatie and positive reversed and
    # find lowest positive
    all_positive_labels = [x for x in vcat(collect.(labels)...) if x > 0]
    offset = length(all_positive_labels) > 0 ? minimum(all_positive_labels) - 1 : 0
    [Tuple([x > 0 ? -(x - offset) : -x for x in ls]) for ls in labels]
end

"""
    replace_tensor_symbol!(tn::TensorNetwork, orig_sym::Symbol, new_sym::Symbol)

Replace the given symbol with the given new symbol
"""
function replace_tensor_symbol!(tn::TensorNetwork, orig_sym::Symbol, new_sym::Symbol)
    tensor = tn[orig_sym]
    for ind in inds(tensor)
        replace!(tn[ind], orig_sym => new_sym)
    end
    tn.tensor_map[new_sym] = tensor
    delete!(tn.tensor_map, orig_sym)
end

"""
    use_mock_tensors(tn::TensorNetwork)

Function to create a copy of the given tensor network with the storage replaced by
mock tensors
"""
function use_mock_tensors(tn::TensorNetwork)
    new_tensor_map = OrderedDict{Symbol, ITensor}()
    bond_map = OrderedDict{Index, Vector{Symbol}}()
    for (sym, tensor) in pairs(tn)
        mock_tensor = MockTensor{ComplexF64}(collect(size(tensor_data(tn, sym))))
        new_tensor_map[sym] = ITensor(mock_tensor, inds(tensor))
    end
    for i in bonds(tn)
        bond_map[i] = tn[i]
    end
    TensorNetwork(new_tensor_map, bond_map)
end
"""
In this module we define a MockTensor which inherits from the NDTensors.TensorStorage.
This can be used as the store for ITensor objects, but only stores the dimensions and so
can be used for tracking the dimension of contractions which would not be possible if
storing all the data.
"""

using NDTensors
using ITensors
import Base: length, size, getindex, copy
import NDTensors: tensor, contract

export MockTensor, size, length, copy, tensor, contract
export use_mock_tensors

""" Tensor store struct that just tracks tensor dimensions"""
struct MockTensor{Elt} <: NDTensors.TensorStorage{Elt}
    size::Array{Int64, 1}
end

"""Overload functions from base to make MockTensor usable"""
copy(a::MockTensor{Elt}) where Elt = MockTensor{Elt}(a.size)
length(a::MockTensor{Elt}) where Elt = prod(size(a))
size(a::MockTensor{Elt}) where Elt = a.size
getindex(a::MockTensor{Elt}, i::Int64) where Elt = NaN
getindex(a::MockTensor{Elt}, i...) where Elt = NaN
tensor(a::MockTensor{Elt}, inds) where Elt = NDTensors.Tensor{Elt, length(inds), MockTensor{Elt}, ITensors.IndexSet}(inds, a)
getindex(a::NDTensors.Tensor{Elt, N, StoreT, IndsT}, i) where {Elt, N , StoreT <: MockTensor, IndsT} = NaN
getindex(a::NDTensors.Tensor{Elt, N, StoreT, IndsT}, i...) where {Elt, N , StoreT <: MockTensor, IndsT} = NaN


"""
contract(T1::NDTensors.Tensor{Elt, N, StoreT, IndsT},
         labelsT1,
         T2::NDTensors.Tensor{Elt, M, StoreT, IndsT},
         labelsT2,
         labelsR = NDTensors.contract_labels(labelsT1, labelsT2)) where {Elt, N, M, StoreT <: MockTensor, IndsT}

Overloaded contract function from NDTensors which implements
contraction for tensors using MockTensor objects as storage.
"""
function contract(T1::NDTensors.Tensor{Elt, N, StoreT, IndsT},
                  labelsT1,
                  T2::NDTensors.Tensor{Elt, M, StoreT, IndsT},
                  labelsT2,
                  labelsR = NDTensors.contract_labels(labelsT1, labelsT2)) where {Elt, N, M, StoreT <: MockTensor, IndsT}
    final_dims = zeros(Int64, length(labelsR))
    T1_dims, T2_dims = size(T1), size(T2)
    final_inds = Array{ITensors.Index}(undef, length(labelsR))

    for (dims, labels, inds) in zip([T1_dims, T2_dims], [labelsT1, labelsT2], [T1.inds, T2.inds])
        for (i, li) in enumerate(labels)
            pos = findfirst(x -> x == li, labelsR)
            if pos !== nothing
                final_dims[pos] = dims[i]
                final_inds[pos] = inds[i]
            end
        end
    end
    tensor_store = MockTensor{Complex{Float64}}(final_dims)

    return NDTensors.Tensor(tensor_store, final_inds)
end

"""
    use_mock_tensors(tn::TensorNetwork)

Function to create a copy of the given tensor network with the storage replaced by
mock tensors
"""
function use_mock_tensors(tn::TensorNetwork)
    newtn = TensorNetwork()
    for tensor in tn.data
        mock_tensor = MockTensor{ComplexF64}(collect(size(tensor)))
        push!(newtn.data, ITensor(mock_tensor, inds(tensor)))
    end
    newtn
end
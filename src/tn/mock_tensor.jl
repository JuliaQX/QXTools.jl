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

export MockTensor, size, length, copy, tensor, contract_tensors


""" Tensor store struct that just tracks tensor dimensions"""
struct MockTensor{Elt} <: NDTensors.TensorStorage{Elt}
    size::Array{Int64, 1}
end

"""Overload functions from base to make MockTensor usable"""
copy(a::MockTensor{Elt}) where Elt = MockTensor{Elt}(a.size)
length(a::MockTensor{Elt}) where Elt = prod(size(a))
size(a::MockTensor{Elt}) where Elt = a.size
getindex(::MockTensor{Elt}, ::Int64) where Elt = NaN
getindex(::MockTensor{Elt}, i...) where Elt = NaN
tensor(a::MockTensor{Elt}, inds) where Elt = NDTensors.Tensor{Elt, length(inds), MockTensor{Elt}, ITensors.IndexSet}(inds, a)
getindex(::NDTensors.Tensor{Elt, N, StoreT, IndsT}, ::Any) where {Elt, N , StoreT <: MockTensor, IndsT} = NaN
getindex(::NDTensors.Tensor{Elt, N, StoreT, IndsT}, i...) where {Elt, N , StoreT <: MockTensor, IndsT} = NaN


"""
contract(T1::NDTensors.Tensor{Elt, N, StoreT, IndsT},
         labelsT1,
         T2::NDTensors.Tensor{Elt, M, StoreT, IndsT},
         labelsT2,
         labelsR = NDTensors.contract_labels(labelsT1, labelsT2)) where {Elt, N, M, StoreT <: MockTensor, IndsT}

Overloaded contract function from NDTensors which implements
contraction for tensors using MockTensor objects as storage.
"""
function mock_contract(T1::NDTensors.Tensor{Elt1, N, StoreT, IndsT},
                       labelsT1,
                       T2::NDTensors.Tensor{Elt2, M, StoreT, IndsT},
                       labelsT2,
                       labelsR = NDTensors.contract_labels(labelsT1, labelsT2)) where {Elt1, Elt2, N, M, StoreT <: MockTensor, IndsT}

    ITensors.disable_warn_order()
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


function contract_tensors(A::ITensors.ITensor, B::ITensors.ITensor)
    (labelsA,labelsB) = ITensors.compute_contraction_labels(inds(A),inds(B))
    if typeof(store(A)) <: MockTensor && typeof(store(B)) <: MockTensor
        CT = mock_contract(tensor(A), labelsA, tensor(B), labelsB)
    else
        CT = contract(tensor(A), labelsA, tensor(B), labelsB)
    end
    C = itensor(CT)
    warnTensorOrder = ITensors.get_warn_order()
    if !isnothing(warnTensorOrder) > 0 &&
        order(C) >= warnTensorOrder
        #@warn "Contraction resulted in ITensor with $(order(C)) indices, which is greater than or equal to the ITensor order warning threshold $warnTensorOrder. You can modify the threshold with functions like `set_warn_order!(::Int)`, `reset_warn_order!()`, and `disable_warn_order!()`."
        println("Contraction resulted in ITensor with $(order(C)) indices, which is greater than or equal to the ITensor order warning threshold $warnTensorOrder. You can modify the threshold with functions like `set_warn_order!(::Int)`, `reset_warn_order!()`, and `disable_warn_order!()`.")
        show(stdout, MIME"text/plain"(), stacktrace())
        println()
        end
    return C
end


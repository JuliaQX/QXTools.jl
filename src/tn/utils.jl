using LinearAlgebra

"""
    function decompose_gate!(gate_data::Array{<:Number, 4},
                             threshold::AbstractFloat=1e-15)

Function to decompose a tensor into two smaller tensors
"""
function decompose_gate(gate_data::Array{<:Number, 4},
                        threshold::AbstractFloat=1e-15)
    left_positions = [1, 3]
    right_positions = [2, 4]
    dims = size(gate_data)
    left_dims = [dims[x] for x in left_positions]
    right_dims = [dims[x] for x in right_positions]

    A = permutedims(gate_data, vcat(left_positions, right_positions))
    A = reshape(A, Tuple([prod(left_dims), prod(right_dims)]))

    # Use SVD here but QR could also be used
    F = svd(A)

    # find number of singular values above the threshold
    chi = sum(F.S .> threshold)
    s = sqrt.(F.S[1:chi])

    # assume that singular values and basis of U and V matrices are sorted
    # in descending order of singular value
    B = reshape(F.U[:, 1:chi] * Diagonal(s), Tuple(vcat(left_dims, [chi,])))
    C = reshape(Diagonal(s) * F.Vt[1:chi, :], Tuple(vcat([chi,], right_dims)))

    return B, C
end

"""
    function find_hyper_edges(A::AbstractArray{Elt, N}) where {Elt, N}

Function to identify hyper edges of tensors. Returns an array of tuples of indices
of the original tensor which can be identified
"""
function find_hyper_edges(A::AbstractArray{Elt, N}) where {Elt, N}
    if N == 2
        if isdiagonal(A)
            return [(1, 2)]
        else
            return []
        end
    elseif N > 2
        groups = []
        for i in 1:N-1
            for j in i+1:N
                if isdiagonal(A, Pair(i, j))
                    # check if there is overlap iwth any existing groups
                    if any(x -> length(intersect(x, (i, j))) > 0, groups)                        
                        # groups with intersections
                        igroups = [x for x in groups if length(intersect(x, (i, j))) > 0]
                        # filter these from main groups list
                        filter!(x -> length(intersect(x, (i, j))) == 0, groups)
                        merged_group = Tuple(reduce(union, igroups, init=[i, j]))
                        push!(groups, merged_group)
                    else
                        push!(groups, (i, j))
                    end
                end
            end
        end
        return groups
    end
end


"""
    isdiagonal(A::AbstractArray{Elt, N}, Pair{Int64, Int64}) where {Elt, N}

Function to check if the given matrix is diagonal along given axes
"""
function isdiagonal(A::AbstractArray{Elt, N}, p::Pair{Int64, Int64}) where {Elt, N}
    po = collect(1:N)
    function exchange!(a, e1, e2)  a[e1], a[e2] = a[e2], a[e1] end
    exchange!(po, 1, p[1])
    exchange!(po, 2, p[2])
    A = permutedims(A, po)
    d = size(A)
    A = reshape(A, (d[1], d[2], prod(d[3:end])))
    all(x -> isdiagonal(A[:, :, x]), 1:prod(d[3:end]))    
end


"""
    isdiagonal(A::AbstractArray{Elt, 2}) where Elt

Function to check if the given matrix is diagonal
"""
function isdiagonal(A::AbstractArray{Elt, 2}) where Elt
    if real(eltype(A)) <: Integer
        all([A[i, j] == 0 for i in 1:size(A)[1], j in 1:size(A)[2] if i != j])
    else
        all([abs(A[i, j]) < eps(real(eltype(A))) for i in 1:size(A)[1], j in 1:size(A)[2] if i != j])
    end
end
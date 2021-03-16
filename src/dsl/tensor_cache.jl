import JLD2

export TensorCache, save_cache

"""TensorCache structure for storing unique tensors"""
mutable struct TensorCache{Elt}
    tensors::Dict{Vector{Int64}, Array{Tuple{Symbol, Vector{Elt}}, 1}}
    key_dim_map::Dict{Symbol, Vector{Int64}}
    id_val::Int64
    label::String
    TensorCache{T}(label::String="data_") where T = new{T}(Dict{Vector{Int64},
                                                           Array{Tuple{Symbol, Vector{T}}, 1}}(),
                                                           Dict{Symbol, Vector{Int64}}(), 1, label)
    TensorCache() = TensorCache{ComplexF64}()
end

"""
    next_symbol!(tc::TensorCache{T}) where T

Returns the next symbol with format :label_i where i is incremented each time. With default
"data_" label, first symbols will be :data_1, :data_2 etc...
"""
function next_symbol!(tc::TensorCache{T}) where T
    sym = Symbol("$(tc.label)$(tc.id_val)")
    tc.id_val += 1
    sym
end

"""
    push!(tc::TensorCache{T1}, data::Array{T2, N}) where {T1, T2, N}

Function to add the provided tensor to the tensor cache and return the symbol matching the
tensor. If there is already a tensor in the cache which matches (up to numerical precision)
the symbol of the pre-existing tensor is returned.

# Examples
```
julia> tc = TensorCache();
julia> push!(tc, [1.0, π])
:data_1
julia> push!(tc, [1.1, π])
:data_2
julia> push!(tc, [1., π])
:data_1
```
"""
function Base.push!(tc::TensorCache{T1}, data::Array{T2, N}) where {T1, T2, N}
    dim = collect(size(data))
    flat_data = convert(Array{T1,1}, reshape(data, prod(dim)))
    if !haskey(tc.tensors, dim)
        sym = next_symbol!(tc)
        tc.tensors[dim] = [(sym, flat_data)]
        tc.key_dim_map[sym] = dim
    else
        match = findfirst(x -> x < eps(real(T1)), [maximum(abs.(flat_data .- x[2])) for x in tc.tensors[dim]])
        if match !== nothing
            sym = tc.tensors[dim][match][1]
        else
            sym = next_symbol!(tc)
            push!(tc.tensors[dim], (sym, flat_data))
            tc.key_dim_map[sym] = dim
        end
    end
    return sym
end

"""
    getindex(tc::TensorCache{T}, sym::Symbol) where T

Overloaded getindex for TensorCache which returns tensor data given a symbol
"""
function Base.getindex(tc::TensorCache{T}, sym::Symbol) where T
    if !haskey(tc.key_dim_map, sym) throw(error("No symbol $sym in cache")) end
    dim = tc.key_dim_map[sym]
    pos = findfirst(x -> x[1] == sym, tc.tensors[dim])
    reshape(tc.tensors[dim][pos][2], Tuple(dim))
end

Base.length(tc::TensorCache{T}) where T = length(tc.key_dim_map)

"""
    save_cache(tc::TensorCache{T}, io::JLD2.JldFile) where T

Function to save tensor data to a JLD file using symbol names as keys
"""
function save_cache(tc::TensorCache{T}, io::JLD2.JLDFile) where T
    for key in keys(tc.key_dim_map)
        write(io, "$(key)", tc[key])
    end
end

"""
    save_cache(tc::TensorCache{T}, filename::String) where T

Function to save tensor data to a JLD file using symbol names as keys
"""
function save_cache(tc::TensorCache{T}, filename::String) where T
    if splitext(filename)[2] != ".jld2" throw(error("Filename must have suffix \".jld2\"")) end
    jldopen(filename, "a+") do io
        save_cache(tc, io)
    end
end

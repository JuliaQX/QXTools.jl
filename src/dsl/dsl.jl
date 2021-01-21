module DSL

using JLD
using ITensors
using ..TN

import Base: length, getindex

export TensorCache, add_tensor!, generate_dsl_file
export length, getindex, save

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
    add_tensor!(tc::TensorCache{T1}, data::Array{T2, N}) where {T1, T2, N}

Function to add the provided tensor to the tensor cache and return the symbol matching the
tensor. If there is already a tensor in the cache which matches (up to numerical precision)
the symbol of the pre-existing tensor is returned.

# Examples
```jldoctest
julia> tc = TensorCache();
julia> add_tensor!(tc, [1.0, π])
:data_1
julia> add_tensor!(tc, [1.1, π])
:data_2
julia> add_tensor!(tc, [1., π])
:data_1
```
"""
function add_tensor!(tc::TensorCache{T1}, data::Array{T2, N}) where {T1, T2, N}
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
function getindex(tc::TensorCache{T}, sym::Symbol) where T
    if !haskey(tc.key_dim_map, sym) throw(error("No symbol $sym in cache")) end
    dim = tc.key_dim_map[sym]
    pos = findfirst(x -> x[1] == sym, tc.tensors[dim])
    reshape(tc.tensors[dim][pos][2], Tuple(dim))
end

length(tc::TensorCache{T}) where T = length(tc.key_dim_map)

"""
    save(filename::String, tc::TensorCache{T}) where T

Function to save tensor data to a JLD file using symbol names as keys
"""
function save(filename::String, tc::TensorCache{T}) where T
    if splitext(filename)[2] != ".jld" throw(error("Filename must have suffix \".jld\"")) end
    JLD.save(filename, vcat([[String(key), tc[key]] for key in keys(tc.key_dim_map)]...)...)
end

"""
    generate_dsl_file(tn::TensorNetwork,
                      plan::Array{Index, 1},
                      dsl_filename::String,
                      data_filename::String)

Function to create a dsl file and data from from the given tensor network and contraction
plan
"""
function generate_dsl_file(tn::TensorNetwork,
                           plan::Array{Index, 1},
                           dsl_filename::String,
                           data_filename::String)
    # create a tensor cache for saving tensor data
    tc = TensorCache()

    tensor_symbols = [Symbol("tensor_$i" for i in 1:length(tn.data))]
    data_symbols = Dict{Symbol, Symbol}()
    # iterate over tensors and add to TensorCache
    open(dsl_filename, "w") do io
        for (tensor_symbol, tensor) in zip(tensor_symbols, tn.data)
            data_symbols[tensor_symbol] = add_tensor!(tc, store(tensor))
            writeln(io, "load $tensor_symbol $data_symbol")
        end
    end
    # save tensor data labels
    save(data_filename, tc)

    # tnc_mock = use_mock_tensors(tnc)

    # for index in plan

end

generate_dsl_file(tnc::TensorNetworkCircuit, args...) = generate_dsl_file(tnc.tn, args...)

end # module
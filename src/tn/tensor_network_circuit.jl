using ITensors

# TensorNetworkCircuit struct and public functions
export TensorNetworkCircuit, add_input!, add_output!, qubits
export contract_tn!, decompose_tensor!
export input_tensors, output_tensors, input_indices, output_indices, single_amplitude

const input_output_tensors = Dict{Char,Array{Float64, 1}}('0' => [1., 0.],
                                                          '1' => [0., 1.],
                                                          '+' => [1., 1.]./sqrt(2),
                                                          '-' => [1., -1.]./sqrt(2))

"""Tensor network circuit data-structure"""
mutable struct TensorNetworkCircuit
    qubits::Int64
    input_indices::Array{Index, 1}
    output_indices::Array{Index, 1}
    input_tensors::Array{Union{Symbol, Nothing}, 1}
    output_tensors::Array{Union{Symbol, Nothing}, 1}
    tn::TensorNetwork
end

function TensorNetworkCircuit(qubits::Int64)
    input_indices = [Index(2, tags="Qubit $i") for i in 1:qubits]
    output_indices = copy(input_indices)
    input_tensors = Array{Union{Symbol, Nothing}, 1}(nothing, qubits)
    output_tensors = Array{Union{Symbol, Nothing}, 1}(nothing, qubits)
    TensorNetworkCircuit(qubits, input_indices, output_indices, input_tensors, output_tensors, TensorNetwork())
end

function Base.copy(tnc::TensorNetworkCircuit)
    tn = TensorNetwork()    
    f = x -> x === nothing ? nothing : push!(tn, tnc[x])
    my_input_tensors = convert(Array{Union{Symbol, Nothing}}, f.(input_tensors(tnc)))
    my_output_tensors = convert(Array{Union{Symbol, Nothing}}, f.(output_tensors(tnc)))
    for x in keys(tnc)
        if !(x in input_tensors(tnc)) && !(x in output_tensors(tnc))
            push!(tn, tnc[x])
        end
    end

    TensorNetworkCircuit(qubits(tnc),
                         copy(input_indices(tnc)),
                         copy(output_indices(tnc)),
                         my_input_tensors,
                         my_output_tensors,
                         tn)
end

Base.eltype(::TensorNetworkCircuit) = ITensor
Base.getindex(tnc::TensorNetworkCircuit, i::Symbol) = tnc.tn[i]
Base.getindex(tnc::TensorNetworkCircuit, i::T) where T <: Index = tnc.tn[i]
Base.iterate(tnc::TensorNetworkCircuit) = iterate(values(tnc))
Base.iterate(tnc::TensorNetworkCircuit, i) = iterate(values(tnc), i)
Base.keys(tnc::TensorNetworkCircuit) = keys(tnc.tn)
Base.length(tnc::TensorNetworkCircuit) = length(tnc.tn)
Base.show(io::IO, ::MIME"text/plain", tnc::TensorNetworkCircuit) = print(io, "TensorNetworkCircuit(qubits => $(qubits(tnc)), gates => $(length(tnc)))")
Base.values(tnc::TensorNetworkCircuit) = values(tnc.tn)

input_tensors(tnc::TensorNetworkCircuit) = tnc.input_tensors
output_tensors(tnc::TensorNetworkCircuit) = tnc.output_tensors
input_indices(tnc::TensorNetworkCircuit) = tnc.input_indices
output_indices(tnc::TensorNetworkCircuit) = tnc.output_indices
bonds(tnc::TensorNetworkCircuit) = bonds(tnc.tn)
qubits(tnc::TensorNetworkCircuit) = tnc.qubits
tensor_data(tnc::TensorNetworkCircuit, i) = tensor_data(tnc.tn, i)


"""
    push!(tnc::TensorNetworkCircuit,
          qubits::Vector{Int64},
          data::Array{T, 2}) where T

Function to add a gate to the tensor network circuit given
the qubits it acts on and an array of the matrix elements
"""
function Base.push!(tnc::TensorNetworkCircuit,
                    qubits::Vector{Int64},
                    data::Array{T, 2};
                    diagonal_check::Bool=false,
                    decompose::Bool=true) where T                    
    input_indices = tnc.output_indices[qubits]
    tnc.output_indices[qubits] = output_indices = [prime(x) for x in input_indices]
    if length(qubits) == 1 || !decompose         
        indices = [output_indices..., input_indices...]
        @assert prod(size(data)) == prod(dim.(indices)) "Data matrix dimension does not match indices"
        data = reshape(data, Tuple(dim.(indices)))
        push!(tnc, indices, data)
    elseif length(qubits) == 2 && decompose         
        data = reshape(data, Tuple([dim.(output_indices)..., dim.(input_indices)...]))
        A, B = decompose_gate(data)
        if size(A)[3] == size(B)[1] == 1
            A = reshape(A, Tuple(size(A)[1:2]))
            B = reshape(A, Tuple(size(B)[2:3]))
            a_indices = [output_indices[1], input_indices[1]]
            b_indices = [output_indices[2], input_indices[2]]
        else
            virtual_index = Index(size(A)[3])
            a_indices = [output_indices[1], input_indices[1], virtual_index]
            b_indices = [virtual_index, output_indices[2], input_indices[2]]
        end
        push!(tnc, a_indices, A)
        push!(tnc, b_indices, B)    
    else

        @error("Gates spanning more than two qubits not supported with decomposition yet.")
    end
end

Base.push!(tnc::TensorNetworkCircuit, args...) = Base.push!(tnc.tn, args...)

"""
    delete!(tnc::TensorNetworkCircuit, tensor_id::Symbol)

Function to remove a tensor from a tensor network circuit.
"""
function Base.delete!(tnc::TensorNetworkCircuit, tensor_id::Symbol)
    replace!(tnc.input_tensors, tensor_id => nothing)
    replace!(tnc.output_tensors, tensor_id => nothing)
    delete!(tnc.tn, tensor_id)
end

"""
    push_input!(tnc::TensorNetworkCircuit, tensor::Array{Elt, 1}, pos::Int64) where Elt

Function to add a single input tensor to the tensor network circuit at the given position
"""
function push_input!(tnc::TensorNetworkCircuit, tensor::Array{Elt, 1}, pos::Int64) where Elt
    index = tnc.input_indices[pos]
    tensor_sym = tnc.input_tensors[pos]
    if tensor_sym !== nothing
        delete!(tnc, tensor_sym)
    end
    tnc.input_tensors[pos] = push!(tnc, [index], tensor)
end

"""
    push_output!(tnc::TensorNetworkCircuit, tensor::Array{Elt, 1}, pos::Int64) where Elt

Function to add a single output tensor to the tensor network circuit at the given position
"""
function push_output!(tnc::TensorNetworkCircuit, tensor::Array{Elt, 1}, pos::Int64) where Elt
    index = tnc.output_indices[pos]
    tensor_sym = tnc.output_tensors[pos]
    if tensor_sym !== nothing
        delete!(tnc, tensor_sym)
    end
    tnc.output_tensors[pos] = push!(tnc, [index], tensor)    
end

"""
    add_input!(tnc::TensorNetworkCircuit; input::Union{String, Nothing}=nothing)

Function to add input tensors to the circuit
"""
function add_input!(tnc::TensorNetworkCircuit, input::Union{String, Nothing}=nothing)
    if input === nothing input = "0"^qubits(tnc) end
    [push_input!(tnc, input_output_tensors[input[pos]], pos) for pos in 1:qubits(tnc)]
end

"""
    add_output!(tnc::TensorNetworkCircuit; output::Union{String, Nothing}=nothing)

Function to add output tensors to the circuit
"""
function add_output!(tnc::TensorNetworkCircuit, output::Union{String, Nothing}=nothing)    
    if output === nothing output = "0"^qubits(tnc) end
    [push_output!(tnc, input_output_tensors[output[pos]], pos) for pos in 1:qubits(tnc)]
end

contract_tn!(tnc::TensorNetworkCircuit, plan) = contract_tn!(tnc.tn, plan)

contract_tn(tnc::TensorNetworkCircuit, plan) = contract_tn!(copy(tnc), plan)

function decompose_tensor!(tnc::TensorNetworkCircuit, args...; kwargs...)                            
    decompose_tensor!(tnc.tn, args...; kwargs...)
end

simple_contraction(tnc::TensorNetworkCircuit) = simple_contraction(tnc.tn)


function single_amplitude(tnc::TensorNetworkCircuit, plan::Array{<:Index, 1}, amplitude::Union{String, Nothing}=nothing)
    sim_tnc = copy(tnc)
    add_output!(sim_tnc, amplitude)
    output = contract_tn!(sim_tnc, plan)
    output[1]
end


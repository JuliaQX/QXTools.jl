using ITensors

import QXZoo
using ..Circuits

# TensorNetworkCircuit struct and public functions
export TensorNetworkCircuit, add_input!, add_output!, convert_to_tnc, qubits
export quickbb_contraction_plan, contract_tn!, contraction_scheme, decompose_tensor!
export input_tensors, output_tensors, input_indices, output_indices

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

    function TensorNetworkCircuit(qubits::Int64)
        input_indices = [Index(2, tags="Qubit $i") for i in 1:qubits]
        output_indices = copy(input_indices)
        input_tensors = Array{Union{Symbol, Nothing}, 1}(nothing, qubits)
        output_tensors = Array{Union{Symbol, Nothing}, 1}(nothing, qubits)
        new(qubits, input_indices, output_indices, input_tensors, output_tensors, TensorNetwork())
    end
end

Base.eltype(tnc::TensorNetworkCircuit) = ITensor
Base.getindex(tnc::TensorNetworkCircuit, i::Symbol) = tnc.tn[i]
Base.getindex(tnc::TensorNetworkCircuit, i::T) where T <: Index = tnc.tn[i]
Base.iterate(tnc::TensorNetworkCircuit) = iterate(values(tnc))
Base.iterate(tnc::TensorNetworkCircuit, i) = iterate(values(tnc), i)
Base.keys(tnc::TensorNetworkCircuit) = keys(tnc.tn)
Base.length(tnc::TensorNetworkCircuit) = length(tnc.tn)
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
                    data::Array{T, 2}) where T
    input_indices = tnc.output_indices[qubits]
    output_indices = [prime(x) for x in input_indices]
    tnc.output_indices[qubits] = output_indices
    indices = [output_indices..., input_indices...]
    @assert prod(size(data)) == prod(dim.(indices)) "Data matrix dimension does not match indices"
    data = reshape(data, Tuple(dim.(indices)))
    push!(tnc, indices, data)
end

Base.push!(tnc::TensorNetworkCircuit, args...) = Base.push!(tnc.tn, args...)

"""
_convert_to_tnc(circ::QXZoo.Circuit.Circ)

Function to convert a QXZoo circuit to a QXSim tensor network circuit
"""
function _convert_to_tnc(circ::QXZoo.Circuit.Circ)
    tnc = TensorNetworkCircuit(circ.num_qubits)
    for gate in circ.circ_ops
        mat_elements = gate_matrix(gate)
        qubits = gate_qubits(gate)
        push!(tnc, qubits, mat_elements)
    end
    tnc
end

"""
convert_to_tnc(circ::QXZoo.Circuit.Circ;
               input::Union{String, Nothing}=nothing,
               output::Union{String, Nothing}=nothing,
               no_input::Bool=false,
               no_output::Bool=false)

Function to convert a QXZoo circuit to a QXSim tensor network circuit
"""
function convert_to_tnc(circ::QXZoo.Circuit.Circ;
                        input::Union{String, Nothing}=nothing,
                        output::Union{String, Nothing}=nothing,
                        no_input::Bool=false,
                        no_output::Bool=false)
    tnc = _convert_to_tnc(circ)
    if !no_input add_input!(tnc, input) end
    if !no_output add_output!(tnc, output) end
    tnc
end

"""
    push_input!(tnc::TensorNetworkCircuit, tensor::Array{Elt, 1}, pos::Int64) where Elt

Function to add a single input tensor to the tensor network circuit at the given position
"""
function push_input!(tnc::TensorNetworkCircuit, tensor::Array{Elt, 1}, pos::Int64) where Elt
    @assert tnc.input_tensors[pos] === nothing "There is already an input tensor at position $pos"
    index = tnc.input_indices[pos]
    tnc.input_tensors[pos] = push!(tnc, [index], tensor)
end

"""
    push_output!(tnc::TensorNetworkCircuit, tensor::Array{Elt, 1}, pos::Int64) where Elt

Function to add a single output tensor to the tensor network circuit at the given position
"""
function push_output!(tnc::TensorNetworkCircuit, tensor::Array{Elt, 1}, pos::Int64) where Elt
    @assert tnc.output_tensors[pos] === nothing "There is already an input tensor at position $pos"
    index = tnc.output_indices[pos]
    tnc.output_tensors[pos] = push!(tnc, [index], tensor)
end

"""
    add_input!(tnc::TensorNetworkCircuit; input::Union{String, Nothing}=nothing)

Function to add input tensors to the circuit
"""
function add_input!(tnc::TensorNetworkCircuit, input::Union{String, Nothing}=nothing)
    @assert all(tnc.input_tensors .=== nothing) "Circuit already has input tensors"
    if input === nothing input = "0"^qubits(tnc) end
    [push_input!(tnc, input_output_tensors[input[pos]], pos) for pos in 1:qubits(tnc)]
end

"""
    add_output!(tnc::TensorNetworkCircuit; output::Union{String, Nothing}=nothing)

Function to add output tensors to the circuit
"""
function add_output!(tnc::TensorNetworkCircuit, output::Union{String, Nothing}=nothing)
    @assert all(tnc.output_tensors .=== nothing) "Circuit already has output tensors"
    if output === nothing output = "0"^qubits(tnc) end
    [push_output!(tnc, input_output_tensors[output[pos]], pos) for pos in 1:qubits(tnc)]
end

contract_tn!(tnc::TensorNetworkCircuit, plan) = contract_tn!(tnc.tn, plan)

function quickbb_contraction_plan(tnc::TensorNetworkCircuit; 
                                  time::Integer=120, 
                                  order::Symbol=:_)
    quickbb_contraction_plan(tnc.tn; time=time, order=order)
end

function contraction_scheme(tnc::TensorNetworkCircuit, num::Integer; 
                            time::Int=120,
                            order::Symbol=:_,
                            score_function::Symbol=:direct_treewidth)
    contraction_scheme(tnc.tn, num; time=time, order=order, score_function=score_function)
end

function decompose_tensor!(tnc::TensorNetworkCircuit, args...; kwargs...)                            
    decompose_tensor!(tnc.tn, args...; kwargs...)
end

simple_contraction(tnc::TensorNetworkCircuit) = simple_contraction(tnc.tn)

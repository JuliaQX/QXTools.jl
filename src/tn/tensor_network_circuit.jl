using ITensors

import QXZoo
using ..Circuits

# TensorNetworkCircuit struct and public functions
export TensorNetworkCircuit, add_input!, add_output!, convert_to_tnc, qubits
export quickbb_contraction_plan, contract_tn!, contraction_scheme

"""Tensor network circuit data-structure"""
mutable struct TensorNetworkCircuit
    qubits::Int64
    input_indices::Array{Index, 1}
    output_indices::Array{Index, 1}
    tn::TensorNetwork

    function TensorNetworkCircuit(qubits::Int64)
        input_indices = [Index(2, tags="Qubit $i") for i in 1:qubits]
        output_indices = copy(input_indices)
        new(qubits, input_indices, output_indices, TensorNetwork())
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


bonds(tnc::TensorNetworkCircuit) = bonds(tnc.tn)
qubits(tnc::TensorNetworkCircuit) = tnc.qubits
tensor_data(tnc::TensorNetworkCircuit, i) = tensor_data(tnc.tn, i)

"""
    push!(tnc::TensorNetworkCircuit,
          qubits::Vector{Int64},
          data::Array{T, N}) where {T, N}

Function to add a gate to the tensor network circuit given
the qubits it acts on and an array of the matrix elements
"""
function Base.push!(tnc::TensorNetworkCircuit,
                    qubits::Vector{Int64},
                    data::Array{T, N}) where {T, N}
    input_indices = tnc.output_indices[qubits]
    output_indices = [prime(x) for x in input_indices]
    tnc.output_indices[qubits] = output_indices
    push!(tnc.tn, [output_indices..., input_indices...], data)
end

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
    create_input_output_tensors(indices::Array{Index, 1}, input::String)

Function to create array of tensors corresponding to the given indices and config string
"""
function create_input_output_tensors(indices::Array{Index, 1}, input::String)
    input_data = Dict{Char, NDTensors.TensorStorage}('0' => NDTensors.Dense([1., 0.]),
                                                     '1' => NDTensors.Dense([0., 1.]),
                                                     '+' => NDTensors.Dense([1., 1.]./sqrt(2)),
                                                     '-' => NDTensors.Dense([1., -1.]./sqrt(2)))
    @assert length(input) == length(indices) "Input and index set must have the same length"

    TensorNetwork([ITensor(input_data[input[i]], [index]) for (i, index) in enumerate(indices)])
end

"""
    add_input!(tnc::TensorNetworkCircuit; input::Union{String, Nothing}=nothing)

Function to add input tensors to the circuit
"""
function add_input!(tnc::TensorNetworkCircuit, input::Union{String, Nothing}=nothing)
    qubits = tnc.qubits

    # TODO: check if there are already inputs present
    if input === nothing
        input = "0"^qubits
    end

    input_tensors = create_input_output_tensors(tnc.input_indices, input)

    tnc.tn = merge(tnc.tn, input_tensors)
end

"""
    add_output!(tnc::TensorNetworkCircuit; output::Union{String, Nothing}=nothing)

Function to add output tensors to the circuit
"""
function add_output!(tnc::TensorNetworkCircuit, output::Union{String, Nothing}=nothing)
    qubits = tnc.qubits

    # TODO: check if there are already inputs present
    if output === nothing
        output = "0"^qubits
    end

    output_tensors = create_input_output_tensors(tnc.output_indices, output)

    tnc.tn = merge(tnc.tn, output_tensors)
end

simple_contraction(tnc::TensorNetworkCircuit) = simple_contraction(tnc.tn)
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
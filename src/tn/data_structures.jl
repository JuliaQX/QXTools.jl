using ITensors
using QXZoo
import Base: length
using ..Circuits

export TensorNetwork, TensorNetworkCircuit
export compose, add_input!, add_output!, simple_contraction
export convert_to_tnc, add_gate!, tensor_data
export length

"""Tensor network data-structure"""
struct TensorNetwork
    data::Array{ITensor, 1}
    TensorNetwork() = new(Array{ITensor, 1}())
    TensorNetwork(array::Array{<: ITensor, 1}) = new(array)
end

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

"""
    compose(a::TensorNetwork, b::TensorNetwork)

Join the networks together
"""
function compose(a::TensorNetwork, b::TensorNetwork)
    TensorNetwork(vcat(a.data, b.data))
end

"""
    tensor_data(tensor::ITensor)

Get the data associated with given tensor
"""
function tensor_data(tensor::ITensor)
    reshape(convert(Array, store(tensor)), Tuple([dim(x) for x in inds(tensor)]))
end

length(tn::TensorNetwork) = length(tn.data)
tensor_data(tn::TensorNetwork, i::Int64) = tensor_data(tn.data[i])


length(tnc::TensorNetworkCircuit) = length(tnc.tn)
tensor_data(tnc::TensorNetworkCircuit, i::Int64) = tensor_data(tnc.tn, i)


"""
    add_gate!(tnc::TensorNetworkCircuit,
              qubits::Vector{Int64},
              data::Array{T, N}) where {T, N}

Function to add a gate to the tensor network circuit given
the qubits it acts on and an array of the matrix elements
"""
function add_gate!(tnc::TensorNetworkCircuit,
                   qubits::Vector{Int64},
                   data::Array{T, N}) where {T, N}
    input_indices = tnc.output_indices[qubits]
    output_indices = [prime(x) for x in input_indices]
    tnc.output_indices[qubits] = output_indices
    data_store = NDTensors.Dense(reshape(data, prod(size(data))))
    push!(tnc.tn.data, ITensor(data_store, [output_indices..., input_indices...]))
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
        add_gate!(tnc, qubits, mat_elements)
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

    tnc.tn = compose(tnc.tn, input_tensors)
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

    tnc.tn = compose(tnc.tn, output_tensors)
end

"""
    simple_contraction(tn::TensorNetwork)

Function to perfrom a simple contraction, contracting all tensors in order.
Only useful for very small networks for testing.
"""
function simple_contraction(tn::TensorNetwork)
    a = tn.data[1]
    for b in tn.data[2:end]
        a = a * b
    end
    a
end

simple_contraction(tnc::TensorNetworkCircuit) = simple_contraction(tnc.tn)
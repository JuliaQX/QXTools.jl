using ITensors
using QXZoo

"""Tensor network data-structure"""
struct TensorNetwork
    data::Array{ITensor, 1}
    TensorNetwork() = new(Array{ITensor, 1}())
    TensorNetwork(array::Array{<: ITensor, 1}) = new(array)
end

"""
    compose(a::TensorNetwork, b::TensorNetwork)

Join the networks together
"""
function compose(a::TensorNetwork, b::TensorNetwork)
    TensorNetwork(vcat(a.data, b.data))
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
_convert_to_tnc(circ::QXZoo.Circuit.Circ)

Function to convert a QXZoo circuit to a QXSim tensor network circuit
"""
function _convert_to_tnc(circ::QXZoo.Circuit.Circ)
    tnc = TensorNetworkCircuit(circ.num_qubits)

    for gate in circ.circ_ops
        mat_elements = convert(Array{ComplexF64}, QXZoo.GateMap.gates[gate.gate_symbol]())
        mat_elements = reshape(mat_elements, prod(size(mat_elements)))
        if typeof(gate) <: QXZoo.GateOps.GateCall1
            qubit = gate.target + 1
            input_index = tnc.output_indices[qubit]
            output_index = prime(input_index)
            tnc.output_indices[qubit] = output_index
            push!(tnc.tn.data, ITensor(NDTensors.Dense(mat_elements), [output_index, input_index]))
        elseif typeof(gate) <: QXZoo.GateOps.GateCall2
            qubit1 = gate.target + 1
            qubit2 = gate.ctrl + 1
            input_indices = tnc.output_indices[[qubit1, qubit2]]
            output_indices = [prime(x) for x in input_indices]
            tnc.output_indices[[qubit1, qubit2]] = output_indices
            push!(tnc.tn.data, ITensor(NDTensors.Dense(mat_elements), [output_indices..., input_indices...]))
        end
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
    create_test_circuit()

Function to create a 3 qubit QXZoo ghz state preparation circuit useful for testing
"""
function create_test_circuit()
    circ = QXZoo.Circuit.Circ(3)
    QXZoo.Circuit.add_gatecall!(circ, QXZoo.DefaultGates.h(0))
    QXZoo.Circuit.add_gatecall!(circ, QXZoo.DefaultGates.c_x(0, 1))
    QXZoo.Circuit.add_gatecall!(circ, QXZoo.DefaultGates.c_x(1, 2))
    circ
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

    [ITensor(input_data[input[i]], [index]) for (i, index) in enumerate(indices)]
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

    tnc.tn = compose(tnc.tn, TensorNetwork(input_tensors))
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

    tnc.tn = compose(tnc.tn, TensorNetwork(output_tensors))
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

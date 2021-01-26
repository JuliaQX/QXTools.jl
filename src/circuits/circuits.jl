module Circuits

using QXZoo

export create_test_circuit, create_qft_circuit
export gate_matrix, gate_qubits

"""
Functions to create circuits
"""

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
    create_qft_circuit(n::Integer)

Function to create a qft circuit with n qubits
"""
function create_qft_circuit(n::Integer)
    circ = QXZoo.Circuit.Circ(n)
    QXZoo.QFT.apply_qft!(circ, collect(0:n-1))
    QXZoo.Circuit.add_gatecall!(circ, QXZoo.DefaultGates.c_x(0, 1))
    QXZoo.Circuit.add_gatecall!(circ, QXZoo.DefaultGates.c_x(1, 2))
    circ
end

"""
    create_rqc_circuit(num_qubits::Integer, rows::Int, cols::Int, depth::Int)

Function to create a qft circuit with n qubits
"""
function create_rqc_circuit(rows::Int, cols::Int, depth::Int; final_h::Bool=false)
    return QXZoo.RQC.create_RQC(rows, cols, depth; final_Hadamard_layer=final_h)
end

"""
    gate_matrix(gate::QXZoo.GateOps.AGateCall)

Function to get the matrix for a gate
"""
function gate_matrix(gate::QXZoo.GateOps.AGateCall)
    mat_elements = convert(Array{ComplexF64}, QXZoo.GateMap.gates[gate.gate_symbol]())
    # TODO: find better way of identifying multi qubit gates
    dims = size(mat_elements)
    new_dims = Tuple(vcat([ones(Int64, convert(Int64, log2(x))) * 2 for x in dims]...))
    reshape(mat_elements, new_dims)
end

"""
    gate_qubits(gate::QXZoo.GateOps.AGateCall)

Function to get the qubit indices that the gate acts on
"""
function gate_qubits(gate::QXZoo.GateOps.AGateCall)
    if typeof(gate) <: QXZoo.GateOps.GateCall1
        qubits = [gate.target + 1]
    elseif typeof(gate) <: QXZoo.GateOps.GateCall2
        qubits = [gate.target + 1, gate.ctrl + 1]
    end
    qubits
end

end  # end module
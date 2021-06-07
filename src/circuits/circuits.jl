module Circuits

using QXZoo

export create_test_circuit, create_qft_circuit, create_ghz_circuit
export create_rqc_circuit, gate_matrix, gate_qubits

"""
Functions to create circuits
"""

"""
    create_test_circuit()

Function to create a 3 qubit QXZoo ghz state preparation circuit useful for testing
"""
function create_test_circuit()
    circ = QXZoo.Circuit.Circ(3)
    QXZoo.Circuit.add_gatecall!(circ, QXZoo.DefaultGates.h(1))
    QXZoo.Circuit.add_gatecall!(circ, QXZoo.DefaultGates.c_x(1, 2))
    QXZoo.Circuit.add_gatecall!(circ, QXZoo.DefaultGates.c_x(2, 3))
    circ
end

"""
    create_tescreate_ghz_circuitt_circuit(n::Int64)

Function to create an n qubit QXZoo ghz state preparation circuit
"""
function create_ghz_circuit(n::Int64)
    circ = QXZoo.Circuit.Circ(n)
    QXZoo.Circuit.add_gatecall!(circ, QXZoo.DefaultGates.h(1))
    for i in 1:n-1
        QXZoo.Circuit.add_gatecall!(circ, QXZoo.DefaultGates.c_x(i, i+1))
    end
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
    create_rqc_circuit(rows::Int,
                       cols::Int,
                       depth::Int,
                       seed::Union{Int, Nothing}=nothing;
                       final_h::Bool=false)

Function to create a qft circuit with n qubits
"""
function create_rqc_circuit(rows::Int,
                            cols::Int,
                            depth::Int,
                            seed::Union{Int, Nothing}=nothing;
                            final_h::Bool=false)
    return QXZoo.RQC.create_RQC(rows, cols, depth, seed; final_Hadamard_layer=final_h)
end

"""
    gate_matrix(gate::QXZoo.GateOps.AGateCall)

Function to get the matrix for a gate
"""
function gate_matrix(gate::QXZoo.GateOps.AGateCall)
    convert(Array{ComplexF64}, QXZoo.GateMap.gates[gate.gate_symbol]())
end

"""
    gate_qubits(gate::QXZoo.GateOps.AGateCall)

Function to get the qubit indices that the gate acts on
"""
function gate_qubits(gate::QXZoo.GateOps.AGateCall)
    if typeof(gate) <: QXZoo.GateOps.GateCall1
        qubits = [gate.target]
    elseif typeof(gate) <: QXZoo.GateOps.GateCall2
        qubits = [gate.target, gate.ctrl]
    end
    qubits
end

end  # end module
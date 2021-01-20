using ITensors
using QXZoo


"""Tensor network data-structure"""
struct TensorNetwork
    data::Array{ITensor, 1}
    TensorNetwork() = new(Array{ITensor, 1}())
end

"""Tensor network circuit data-structure"""
struct TensorNetworkCircuit
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
    circuit_to_tn(circ::QXZoo.Circuit.Circ)

Function to convert a QXZoo circuit to a QXSim tensor network
"""
function circuit_to_tn(circ::QXZoo.Circuit.Circ)
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
            qubit1 = gate.ctrl + 1
            qubit2 = gate.target + 1
            input_indices = tnc.output_indices[[qubit1, qubit2]]
            output_indices = [prime(x) for x in input_indices]
            tnc.output_indices[[qubit1, qubit2]] = output_indices
            push!(tnc.tn.data, ITensor(NDTensors.Dense(mat_elements), [output_indices..., input_indices...]))
        end
    end
    tnc
end

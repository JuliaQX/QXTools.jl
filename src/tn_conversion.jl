using .Circuits
using .TN
using QXZoo

export convert_to_tnc

"""
_convert_to_tnc(circ::QXZoo.Circuit.Circ; kwargs...)

Function to convert a QXZoo circuit to a QXSim tensor network circuit
"""
function _convert_to_tnc(circ::QXZoo.Circuit.Circ; kwargs...)
    tnc = TensorNetworkCircuit(circ.num_qubits)
    for gate in circ.circ_ops
        mat_elements = gate_matrix(gate)
        qubits = gate_qubits(gate)
        push!(tnc, qubits, mat_elements; kwargs...)
    end
    tnc
end

"""
convert_to_tnc(circ::QXZoo.Circuit.Circ;
               input::Union{String, Nothing}=nothing,
               output::Union{String, Nothing}=nothing,
               no_input::Bool=false,
               no_output::Bool=false,
               kwargs...)

Function to convert a QXZoo circuit to a QXSim tensor network circuit
"""
function convert_to_tnc(circ::QXZoo.Circuit.Circ;
                        input::Union{String, Nothing}=nothing,
                        output::Union{String, Nothing}=nothing,
                        no_input::Bool=false,
                        no_output::Bool=false,
                        kwargs...)
    tnc = _convert_to_tnc(circ; kwargs...)
    if !no_input add_input!(tnc, input) end
    if !no_output add_output!(tnc, output) end
    tnc
end
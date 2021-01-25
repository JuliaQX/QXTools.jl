using ITensors
import QXZoo
import Base: length, iterate, keys, push!, getindex, merge, eltype
using ..Circuits

# TensorNetwork struct and public functions
export next_tensor_id
export TensorNetwork, bonds, simple_contraction, tensor_data, neighbours
# TensorNetworkCircuit struct and public functions
export TensorNetworkCircuit, add_input!, add_output!, convert_to_tnc

# Overloaded functions for dealing with TensorNetworks
export length, iterate, keys, push!, getindex, merge, eltype

const qxsim_ids = Dict{Symbol, Int64}(:tensor_id => 0)
@noinline next_tensor_id() = begin qxsim_ids[:tensor_id] += 1; Symbol("t$(qxsim_ids[:tensor_id])") end

"""Tensor network data-structure"""
mutable struct TensorNetwork
    tensor_map::Dict{Symbol, ITensor}
    bond_map::Dict{Index, Vector{Symbol}}

    # constructors
    TensorNetwork(tensor_map::Dict{Symbol, ITensor},
                  bond_map::Dict{<:Index, Vector{Symbol}}) = new(tensor_map, bond_map)
    TensorNetwork() = new(Dict{Symbol, ITensor}(), Dict{Index, Vector{Symbol}}())
end

"""
    TensorNetwork(array::Vector{<: ITensor})

Outer consturctor to create a tensor network object from an array of
ITensor objects
"""
function TensorNetwork(array::Vector{<: ITensor})
    tensor_map = Dict{Symbol, ITensor}()
    bond_map = Dict{Index, Vector{Symbol}}()
    for (i, tensor) in enumerate(array)
        tensor_id = next_tensor_id()
        tensor_map[tensor_id] = tensor
        for bond in inds(tensor)
            if haskey(bond_map, bond)
                push!(bond_map[bond], tensor_id)
            else
                bond_map[bond] = [tensor_id]
            end
        end
    end
    TensorNetwork(tensor_map, bond_map)
end

length(tn::TensorNetwork) = length(tn.tensor_map)
tensor_data(tn::TensorNetwork, i::Symbol) = tensor_data(tn.tensor_map[i])
iterate(tn::TensorNetwork) = iterate(values(tn.tensor_map))
iterate(tn::TensorNetwork, state) = iterate(values(tn.tensor_map), state)
eltype(tn::TensorNetwork) = ITensor
keys(tn::TensorNetwork) = keys(tn.tensor_map)
getindex(tn::TensorNetwork, i::Symbol) = tn.tensor_map[i]
getindex(tn::TensorNetwork, i::T) where T <: Index = tn.bond_map[i]
bonds(tn::TensorNetwork) = keys(tn.bond_map)

"""
    neighbours(tn::TensorNetwork, tensor::Symbol)

Function get the symbols of the neighbouring tensors
"""
function neighbours(tn::TensorNetwork, tensor::Symbol)
    tensor_indices = inds(tn[tensor])
    connected_tensors = unique(vcat([tn[x] for x in tensor_indices]...))
    setdiff(connected_tensors, [tensor])
end

"""
    merge(a::TensorNetwork, b::TensorNetwork)

Join two networks together
"""
function merge(a::TensorNetwork, b::TensorNetwork)
    tensor_map = merge(a.tensor_map, b.tensor_map)
    bond_map = Dict{Index, Vector{Symbol}}()
    for (tid, tensor) in pairs(tensor_map)
        for bond in inds(tensor)
            if !haskey(bond_map, bond)
                bond_map[bond] = [tid]
            else
                push!(bond_map[bond], tid)
            end
        end
    end
    c = TensorNetwork(tensor_map, bond_map)

end

"""
    tensor_data(tensor::ITensor)

Get the data associated with given tensor
"""
function tensor_data(tensor::ITensor)
    reshape(convert(Array, store(tensor)), Tuple([dim(x) for x in inds(tensor)]))
end

"""
    push!(tn::TensorNetwork,
          indices::Vector{Index},
          data::Array{T, N}) where {T, N}

Function to add a tensor to the tensor network tensor data and indices
"""
function push!(tn::TensorNetwork,
               indices::Vector{<:Index},
               data::Array{T, N}) where {T, N}
    @assert size(data) == Tuple(dim.(indices))
    data_store = NDTensors.Dense(reshape(data, prod(size(data))))
    tensor = ITensor(data_store, indices)
    tid = next_tensor_id()
    tn.tensor_map[tid] = tensor
    for bond in indices
        if haskey(tn.bond_map, bond)
            push!(tn.bond_map[bond], tid)
        else
            tn.bond_map[bond] = [tid]
        end
    end
    tid
end

"""
    simple_contraction(tn::TensorNetwork)

Function to perfrom a simple contraction, contracting all tensors in order.
Only useful for very small networks for testing.
"""
function simple_contraction(tn::TensorNetwork)
    a, state = iterate(tn)
    for b in iterate(tn, state)
        a = a * b
    end
    a
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

length(tnc::TensorNetworkCircuit) = length(tnc.tn)
tensor_data(tnc::TensorNetworkCircuit, i::Int64) = tensor_data(tnc.tn, i)
iterate(tnc::TensorNetworkCircuit, i) = iterate(tnc.tn::TensorNetwork, i)

"""
    push!(tnc::TensorNetworkCircuit,
          qubits::Vector{Int64},
          data::Array{T, N}) where {T, N}

Function to add a gate to the tensor network circuit given
the qubits it acts on and an array of the matrix elements
"""
function push!(tnc::TensorNetworkCircuit,
               qubits::Vector{Int64},
               data::Array{T, N}) where {T, N}
    input_indices = tnc.output_indices[qubits]
    output_indices = [prime(x) for x in input_indices]
    tnc.output_indices[qubits] = output_indices
    data_store = NDTensors.Dense(reshape(data, prod(size(data))))
    push!(tnc.tn, data, [output_indices..., input_indices...])
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

simple_contraction(tnc::TensorNetworkCircuit) = simple_contraction(tnc.tn)
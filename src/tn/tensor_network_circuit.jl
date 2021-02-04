using ITensors

import QXZoo
using ..Circuits

# TensorNetworkCircuit struct and public functions
export TensorNetworkCircuit, add_input!, add_output!, convert_to_tnc, qubits
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
    input_indices = [Index(2, tags="Qubit $i,Hyp 1") for i in 1:qubits]
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
                    diagonal::Bool=false) where T
    input_indices = tnc.output_indices[qubits]
    if diagonal
        output_indices = [prime(x) for x in input_indices]
    else
        tsold = [ind.tags[1] for ind in input_indices]
        tsnew = [parse(Int, String(tag)[4:end]) + 1 for tag in tsold]
        tsnew = ["Hyp $num" for num in tsnew]
        output_indices = [prime(replacetags(x, tsold[i] => tsnew[i])) 
                            for (i, x) in enumerate(input_indices)]
    end
    tnc.output_indices[qubits] = output_indices
    indices = [output_indices..., input_indices...]
    @assert prod(size(data)) == prod(dim.(indices)) "Data matrix dimension does not match indices"
    data = reshape(data, Tuple(dim.(indices)))
    push!(tnc, indices, data)
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


function single_amplitude(tnc::TensorNetworkCircuit, plan::Array{<:Index, 1}, amplitude::Union{String, Nothing}=nothing)
    sim_tnc = copy(tnc)
    add_output!(sim_tnc, amplitude)
    output = contract_tn!(sim_tnc, plan)
    
    store(output)[1]
end


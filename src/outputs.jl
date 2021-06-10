using DataStructures
using ArgParse

export create_output_parser, output_params_dict, process_output_args

function create_output_parser()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--output_method"
            help = "Method for selecting amplitudes to output. (uniform, rejection, list)"
            default = :List
            arg_type = Symbol
        "--num_outputs", "-n"
            help = "Number of amplitudes or bitstring samples to output"
            default = nothing
        "--M"
            help = "Constant to use for rejection sampling"
            default = 0.001
            arg_type = Float64
        "--fix_M"
            help = "Fix rejection sampling constant for frugal sampling"
            default = false
            arg_type = Bool
        "--bitstrings"
            help = "The bitstrings to compute amplitudes for if list output method is selected"
            default = nothing
        "--output_seed"
            help = "Seed to use for sampling output bitstrings"
            default = nothing
    end

    return s
end

"""
function output_params_dict(qubits::Int,
                            num_outputs::Int64=10;
                            output_method::Symbol=:List,
                            seed::Union{Int64, Nothing}=nothing,
                            M::Float64=0.0001,
                            fix_M::Bool=false,
                            bitstrings::Union{Vector{String}, Nothing}=nothing)

Function to construct dictionary with appropriate paramters to describe sampling approach.
"""
function output_params_dict(num_qubits::Int,
                            num_outputs::Int=10;
                            output_method::Symbol=:List,
                            seed::Union{Int, Nothing}=nothing,
                            M::Float64=0.0001,
                            fix_M::Bool=false,
                            bitstrings::Union{Vector{String}, Nothing}=nothing)
    output_args = OrderedDict()
    output_args[:method] = output_method
    output_params = OrderedDict()
    if output_method == :Rejection
        output_params[:num_qubits] = num_qubits
        output_params[:M] = M
        output_params[:fix_M] = fix_M
        output_params[:seed] = seed
        output_params[:num_samples] = num_outputs
    elseif output_method == :List
        if bitstrings === nothing
            bitstrings = collect(amplitudes_uniform(num_qubits, seed, num_outputs))
        else num_outputs = length(bitstrings) end
        output_params[:num_samples] = num_outputs
        output_params[:bitstrings] = bitstrings
    elseif output_method == :Uniform
        output_params[:num_qubits] = num_qubits
        output_params[:num_samples] = num_outputs
        output_params[:seed] = seed
    else
        @error("Output method \"$(output_method)\" not supported")
    end
    output_args[:params] = output_params
    output_args
end

"""
function parse_output_args(qubits, parsed_args)

Process parsed output arguments and construct output parameters dictionary
"""
function process_output_args(num_qubits, parsed_args)
    # Output parameters.
    output_method = parsed_args["output_method"]
    num_outputs = parsed_args["num_outputs"]
    if num_outputs === nothing
        num_outputs = 10
    else num_outputs = parse(Int64, num_outputs) end
    M = parsed_args["M"]
    fix_M = parsed_args["fix_M"]
    bitstrings = parsed_args["bitstrings"]
    if bitstrings !== nothing
        if match(r"[01|,]*", bitstrings) === nothing
            @error("Bistrings must only contains '0','1',','")
        end
        bitstrings = split(bitstrings, ",")
    end
    seed = parsed_args["output_seed"]
    seed === nothing || (seed = parse(Int64, seed))

    output_params_dict(num_qubits,
                       num_outputs;
                       output_method=output_method,
                       seed=seed,
                       M=M,
                       fix_M=fix_M,
                       bitstrings=bitstrings)
end
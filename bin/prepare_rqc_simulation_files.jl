using QXSim
using QXSim.Circuits
using ArgParse
using Logging

function parse_commandline(ARGS)
    s = ArgParseSettings()
    @add_arg_table! s begin
        "--rows", "-r"
            help = "Number of rows"
            default = 3
            arg_type = Int64
        "--cols", "-c"
            help = "Number of columns"
            default = 3
            arg_type = Int64
        "--depth", "-d"
            help = "Number layers"
            default = 8
            arg_type = Int64
        "--sliced_bonds", "-s"
            help = "Number bonds to slice"
            default = 2
            arg_type = Int64
        "--amplitudes", "-a"
            help = "Number of amplitudes"
            default = nothing
        "--seed"
            help = "Seed to use for both circuit and amplitude selection"
            default = nothing
        "--decompose"
            help = "Set if two qubit gates should be decomposed."
            default = true
            arg_type = Bool
        "--hypergraph"
            help = "Set if the hypergraph structure of the tensor network should be used."
            default = true
            arg_type = Bool
        "--time"
            help = "The number of seconds to run quickbb for."
            default = 120
            arg_type = Int64
        "--qbb_order"
            help = "The branching order to be used by quickbb."
            default = :min_fill
            arg_type = Symbol
        "--lb"
            help = "Set if a lowerbound for the treewidth should be computed."
            default = false
            arg_type = Bool
        "--score_function"
            help = "Function to maximise when selecting vertices to remove."
            default = :direct_treewidth
            arg_type = Symbol
        "--prefix", "-p"
            help = "Prefix to use for output files"
            required = true
            arg_type = String
        "--verbose", "-v"
            help = "Verbose output"
            action = :store_true
    end
    return parse_args(ARGS, s)
end

function main(ARGS)
    parsed_args = parse_commandline(ARGS)

    rows = parsed_args["rows"]
    cols = parsed_args["cols"]
    depth = parsed_args["depth"]
    number_bonds_to_slice = parsed_args["sliced_bonds"]
    num_amplitudes = parsed_args["amplitudes"]
    if num_amplitudes !== nothing
        num_amplitudes = parse(Int64, num_amplitudes)
    end
    output_prefix = parsed_args["prefix"]
    seed = parsed_args["seed"]
    if seed !== nothing
        seed = parse(Int64, seed)
    end
    decompose = parsed_args["decompose"]
    hypergraph = parsed_args["hypergraph"]
    time = parsed_args["time"]
    qbb_order = parsed_args["qbb_order"]
    lb = parsed_args["lb"]
    score_function = parsed_args["score_function"]
    verbose = parsed_args["verbose"]

    @info("Create circuit with $(rows * cols) qubits")
    circ = create_rqc_circuit(rows, cols, depth, seed)
    @info("Circuit created with $(circ.circ_ops.len) gates")

    generate_simulation_files(circ;
                              number_bonds_to_slice=number_bonds_to_slice,
                              output_prefix=output_prefix,
                              num_amplitudes=num_amplitudes,
                              seed=seed,
                              decompose=decompose,
                              hypergraph=hypergraph,
                              time=time,
                              qbb_order=qbb_order,
                              lb=lb,
                              score_function=score_function)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
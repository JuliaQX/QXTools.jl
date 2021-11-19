using QXTools
using QXTools.Circuits
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
        "--seed"
            help = "Seed to use for circuit generation, contraction planning and selecting amplitudes to compute."
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
            help = "The number of seconds to run flow cutter for."
            default = 30
            arg_type = Int64
        "--sliced_bonds", "-s"
            help = "Number bonds to slice"
            default = 2
            arg_type = Int64
        "--score_function"
            help = "Function to maximise when selecting vertices to remove."
            default = :tree_trimming
            arg_type = Symbol
        "--prefix", "-p"
            help = "Prefix to use for output files"
            required = true
            arg_type = String
        "--verbose", "-v"
            help = "Verbose output"
            action = :store_true
    end
    import_settings!(s, create_output_parser())

    return parse_args(ARGS, s)
end

function main(ARGS)
    parsed_args = parse_commandline(ARGS)

    # Citcuit parameters.
    rows = parsed_args["rows"]
    cols = parsed_args["cols"]
    depth = parsed_args["depth"]

    # Tensor network parameters.
    decompose = parsed_args["decompose"]
    hypergraph = parsed_args["hypergraph"]

    # Contraction and slicing parameters.
    time = parsed_args["time"]
    score_function = parsed_args["score_function"]
    number_bonds_to_slice = parsed_args["sliced_bonds"]

    # Other parameters.
    seed = parsed_args["seed"]
    if seed !== nothing
        seed = parse(Int64, seed)
    end
    output_prefix = parsed_args["prefix"]
    verbose = parsed_args["verbose"]

    output_args = process_output_args(rows*cols, parsed_args)

    @info("Create circuit with $(rows * cols) qubits")
    circ = create_rqc_circuit(rows, cols, depth, seed)
    @info("Circuit created with $(circ.circ_ops.len) gates")

    generate_simulation_files(circ,
                              output_prefix,
                              number_bonds_to_slice;
                              decompose=decompose,
                              seed=seed,
                              output_args=output_args,
                              hypergraph=hypergraph,
                              time=time,
                              score_function=score_function)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
using QXSim
using ArgParse
using Logging

function parse_commandline(ARGS)
    s = ArgParseSettings()
    @add_arg_table s begin
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
            default = 4
            arg_type = Int64
        "--seed"
            help = "Seed to use for both circuit and amplitude selection"
            default = nothing
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
    output_prefix = parsed_args["prefix"]
    seed = parsed_args["seed"]
    if seed !== nothing
        seed = parse(Int64, seed)
    end
    verbose = parsed_args["verbose"]

    @info("Create circuit with $(rows * cols) qubits")
    circ = QXSim.create_rqc_circuit(rows, cols, depth, seed)
    @info("Circuit created with $(circ.circ_ops.len) gates")

    generate_simulation_files(circ,
                              number_bonds_to_slice=number_bonds_to_slice,
                              output_prefix=output_prefix,
                              num_amplitudes=num_amplitudes,
                              seed=seed)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
using QXSim
using ArgParse
using Logging

function parse_commandline(ARGS)
    s = ArgParseSettings()
        @add_arg_table! s begin
            "--rows", "-r"
                help = "Number of rows"
                default = 3
            "--cols", "-c"
                help = "Number of columns"
                default = 3
            "--depth", "-d"
                help = "Number layers"
                default = 8
            "--sliced_bonds", "-s"
                help = "Number bonds to slice"
                default = 2
            "--amplitudes", "-a"
                help = "Number of amplitudes"
                default = 4
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
    number_sliced_bonds = parsed_args["sliced_bonds"]
    amplitudes = parsed_args["amplitudes"]
    verbose = parsed_args["verbose"]

    @info("Create circuit with $qubits qubits")
    # TODO: replace with call to rqc circuit creation functions
    # circ = QXSim.create_rqc_circuit(rows, cols, depth)
    circ = QXSim.create_ghz_circuit(rows)

    @info("Convert circuit to tensor network")
    tnc = convert_to_tnc(circ)

    @info("Convert tensor network to graph")
    g = convert_to_graph(tnc)

    @info("Get contraction plan using qxgraph")
    # fill from qxgraph
    # TODO: fill in call
    #plan = get_contraction_plan(g)
    #bonds, dims = get_bonds_to_slice(g, plan, sliced_bonds)

    @info("Prepare DSL and data files")
    # using plan generate the DSL and data files
    # TODO: fill in call
    #generate_dsl_files(tnc, g, plan, bonds)

    @info("Write parameter file for retrieving $amplitudes amplitudes")
    # Write parameter file detailing partitions and amplitudes
    # TODO: fill in call
    #write_parameter_file(tnc, amplitudes, (bonds, dims))
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
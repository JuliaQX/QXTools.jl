# using MPI
using QXContexts
using ArgParse
using TimerOutputs
using LinearAlgebra

"""
    parse_commandline(ARGS)

Parse command line arguments and return argument dictionary
"""
function parse_commandline(ARGS)
    s = ArgParseSettings("QXContexts")

    @add_arg_table! s begin
        "--dsl", "-d"
            help = "DSL file path"
            required = true
            arg_type = String
        "--parameter-file", "-p"
            help = "Parameter file path, default is to use dsl filename with .yml suffix"
            default = nothing
            arg_type = String
        "--input-file", "-i"
            help = "Input data file path, default is to use dsl filename with .jld2 suffix"
            default = nothing
            arg_type = String
        "--output-file", "-o"
            help = "Output data file path"
            required = true
            arg_type = String
        "--number-amplitudes", "-a"
            help = "Number of amplitudes to calculate out of number in parameter file"
            default = nothing
            arg_type = Int
        "--number-slices", "-n"
            help = "The number of slices to use out of number given in parameter file"
            default = nothing
            arg_type = Int
        "--sub-comm-size", "-s"
            help = "The number of ranks to assign to each sub-communicator for partitions"
            default = 1
            arg_type = Int
        "--mpi", "-m"
            help = "Use MPI"
            action = :store_true
        "--gpu", "-g"
            help = "Use GPU if available"
            action = :store_true
        "--timings", "-t"
            help = "Enable output of timings with warmup run"
            action = :store_true
        "--blas-threads", "-b"
            help = "Set number of blas threads"
            default = 8
            arg_type = Int
    end
    return parse_args(ARGS, s)
end

"""
    main(ARGS)

QXContexts entry point
"""
function main(ARGS)
    args = parse_commandline(ARGS)

    dsl_file       = args["dsl"]
    input_file     = args["input-file"]
    param_file = args["parameter-file"]
    output_file    = args["output-file"]
    number_amplitudes = args["number-amplitudes"]
    number_slices  = args["number-slices"]
    sub_comm_size  = args["sub-comm-size"]
    use_mpi        = args["mpi"]
    use_gpu        = args["gpu"]
    timings        = args["timings"]
    blas_threads   = args["blas-threads"]

    BLAS.set_num_threads(blas_threads)

    results = execute(dsl_file, input_file, param_file, output_file;
                      use_mpi=use_mpi, sub_comm_size=sub_comm_size,
                      use_gpu=use_gpu, max_amplitudes=number_amplitudes,
                      max_slices=number_slices,
                      timings=timings)

    if timings
        reset_timer!()
        results = execute(dsl_file, input_file, param_file, output_file;
                          use_mpi=use_mpi, sub_comm_size=sub_comm_size,
                          use_gpu=use_gpu, max_amplitudes=number_amplitudes,
                          max_slices=number_slices,
                          timings=timings)
    end
end


if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
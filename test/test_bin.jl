module TestCLI
using Base: decompose
import Logging
using Test
using YAML

# include source of bin file here to avoid world age issues
include("../bin/prepare_rqc_simulation_files.jl")

@testset "Test prepare rqc input cli script" begin
    # create empty temporary directory
    mktempdir() do path
        prefix = joinpath(path, "rqc_3_3_8")
        args = ["-p", prefix, "--time", "10", "-n", "15"]
        Logging.with_logger(Logging.NullLogger()) do # suppress logging
            main(args)
        end
        @test all([isfile(prefix * suffix) for suffix in [".qx", ".jld2", ".yml"]])

        params = YAML.load_file(prefix * ".yml")
        @test length(params["output"]["params"]["bitstrings"]) == 15
    end

    # create empty temporary directory
    mktempdir() do path
        prefix = joinpath(path, "rqc_3_3_8")
        N = 20
        args = ["-p", prefix, "-n", "$N", "--time", "15", "--output_method", "Rejection"]
        Logging.with_logger(Logging.NullLogger()) do # suppress logging
            main(args)
        end
        @test all([isfile(prefix * suffix) for suffix in [".qx", ".jld2", ".yml"]])

        params = YAML.load_file(prefix * ".yml")
        @test params["output"]["method"] == "Rejection"
        @test params["output"]["params"]["num_samples"] == N
    end

    # create empty temporary directory
    mktempdir() do path
        prefix = joinpath(path, "rqc_3_3_8")
        N = 20
        args = ["-p", prefix, "-n", "$N", "--time", "15"]
        Logging.with_logger(Logging.NullLogger()) do # suppress logging
            main(args)
        end
        @test all([isfile(prefix * suffix) for suffix in [".qx", ".jld2", ".yml"]])

        params = YAML.load_file(prefix * ".yml")
        @test length(params["output"]["params"]["bitstrings"]) <= N
    end
end
end
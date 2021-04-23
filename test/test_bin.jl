module TestCLI
import Logging
using Test
using YAML

# include source of bin file here to avoid world age issues
include("../bin/prepare_rqc_simulation_files.jl")

@testset "Test prepare rqc input cli script" begin
    # create empty temporary directory
    mktempdir() do path
        prefix = joinpath(path, "rqc_3_3_8")
        args = ["-p", prefix, "--time", "30", "-a", "15"]
        Logging.with_logger(Logging.NullLogger()) do # suppress logging
            main(args)
        end
        @test all([isfile(prefix * suffix) for suffix in [".qx", ".jld2", ".yml"]])

        params = YAML.load_file(prefix * ".yml")
        @test length(params["amplitudes"]) == 15
        @test length(params["output"]["bitstrings"]) == 2^9
    end

    # create empty temporary directory
    mktempdir() do path
        prefix = joinpath(path, "rqc_3_3_8")
        N = 20
        args = ["-p", prefix, "-a", "$N", "--time", "30", "--output_method", "rejection"]
        Logging.with_logger(Logging.NullLogger()) do # suppress logging
            main(args)
        end
        @test all([isfile(prefix * suffix) for suffix in [".qx", ".jld2", ".yml"]])

        params = YAML.load_file(prefix * ".yml")
        @test params["output"]["output_method"] == "rejection"
        @test params["output"]["num_samples"] == N
    end

    # create empty temporary directory
    mktempdir() do path
        prefix = joinpath(path, "rqc_3_3_8")
        N = 20
        args = ["-p", prefix, "-a", "$N", "--time", "30"]
        Logging.with_logger(Logging.NullLogger()) do # suppress logging
            main(args)
        end
        @test all([isfile(prefix * suffix) for suffix in [".qx", ".jld2", ".yml"]])

        params = YAML.load_file(prefix * ".yml")
        @test length(params["output"]["bitstrings"]) <= N
    end
end
end
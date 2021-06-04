include("tensor_cache.jl")
include("compute_tree.jl")

using YAML
using FileIO

export generate_dsl_files, generate_parameter_file

"""
    generate_dsl_files(compute_tree::ComputeTree,
                       prefix::String;
                       force::Bool=true,
                       metadata=nothing)

Function to create a dsl and data files to contracting the given tensor network circuit
with the plan provided
"""
function generate_dsl_files(compute_tree::ComputeTree,
                            prefix::String;
                            force::Bool=true,
                            metadata=nothing)

    dsl_filename = "$(prefix).qx"
    data_filename = "$(prefix).jld2"

    @assert force || !isfile(dsl_filename) "Error $(dsl_filename) already exists"
    @assert force || !isfile(data_filename) "Error $(data_filename) already exists"

    open(dsl_filename, "w") do dsl_io
        write(dsl_io, compute_tree; metadata=metadata)
    end

    save(data_filename, Dict(String(x) => y for (x, y) in pairs(compute_tree.tensors)))
    nothing
end

"""
    generate_parameter_file(filename_prefix::String,
                            output_parameters)

Generate a yml file with details of how outputs are sampled

output:
    output_method: rejection
    params:
        fix_M: false
        M: 0.001
        num_samples: 10
        seed: ~
"""
function generate_parameter_file(filename_prefix::String,
                                 output_parameters)
    config = Dict()
    config["output"] = output_parameters

    YAML.write_file("$(filename_prefix).yml", config)
end
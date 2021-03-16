push!(LOAD_PATH, "../")
using QXSim, Documenter

makedocs(;
    modules=[QXSim],
    authors="QuantEx team",
    repo="https://github.com/JuliaQX/QXSim.jl/blob/{commit}{path}#L{line}",
    sitename="QXSim.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaQX.github.io/QXSim.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "User's Guide" => "users_guide.md",
        "Tutorials" => [],
        "LICENSE" => "license.md"
    ],
)

deploydocs(;
    repo="github.com/JuliaQX/QXSim.jl",
)

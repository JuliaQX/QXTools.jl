using QXSim
using Documenter

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
    ],
)

deploydocs(;
    repo="github.com/JuliaQX/QXSim.jl",
)

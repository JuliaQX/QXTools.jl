push!(LOAD_PATH, "../")
using QXTools, Documenter

makedocs(;
    modules=[QXTools],
    authors="QuantEx team",
    repo="https://github.com/JuliaQX/QXTools.jl/blob/{commit}{path}#L{line}",
    sitename="QXTools.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuliaQX.github.io/QXTools.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "User's Guide" => "users_guide.md",
        "Tutorials" => ["basics.md","features.md"],
        "LICENSE" => "license.md"
    ],
)

deploydocs(;
    repo="github.com/JuliaQX/QXTools.jl",
)

using FerriteNeumann
using Documenter

DocMeta.setdocmeta!(FerriteNeumann, :DocTestSetup, :(using FerriteNeumann); recursive=true)

makedocs(;
    modules=[FerriteNeumann],
    authors="Knut Andreas Meyer and contributors",
    repo="https://github.com/KnutAM/FerriteNeumann.jl/blob/{commit}{path}#{line}",
    sitename="FerriteNeumann.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://KnutAM.github.io/FerriteNeumann.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/KnutAM/FerriteNeumann.jl",
    devbranch="main",
)

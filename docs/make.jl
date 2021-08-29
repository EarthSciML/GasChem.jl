using GasChemMTK
using Documenter

DocMeta.setdocmeta!(GasChemMTK, :DocTestSetup, :(using DepositionMTK); recursive=true)

makedocs(;
    modules=[GasChemMTK],
    authors="EarthSciML authors and contributors",
    repo="https://github.com/EarthSciML/GasChemMTK.jl/blob/{commit}{path}#{line}",
    sitename="GasChemMTK.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://EarthSciML.github.io/GasChemMTK.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/EarthSciML/GasChemMTK.jl",
)

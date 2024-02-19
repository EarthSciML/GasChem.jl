using GasChem
using Documenter

DocMeta.setdocmeta!(GasChem, :DocTestSetup, :(using GasChem); recursive = true)

makedocs(;
    modules = [GasChem],
    authors = "EarthSciML authors and contributors",
    repo = "https://github.com/EarthSciML/GasChem.jl/blob/{commit}{path}#{line}",
    sitename = "GasChem.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://EarthSciML.github.io/GasChem.jl",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Composing models" => "composing_models.md",
    ],
)

deploydocs(; repo = "github.com/EarthSciML/GasChem.jl")

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
        "SuperFast" => "superfast.md",
        "GEOS-Chem" => [
            "Overview" => "geoschem/overview.md",
            "State Variables" => "geoschem/states.md",
            "Parameters" => "geoschem/params.md",
        ],
        "Composing models" => "composing_models.md",
        "API" => "api.md",
    ],
)

deploydocs(; repo = "github.com/EarthSciML/GasChem.jl.git")

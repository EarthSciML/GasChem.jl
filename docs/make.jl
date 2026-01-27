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
        canonical = "https://gaschem.earthsci.dev",
        assets = String[],
        repolink = "https://github.com/EarthSciML/GasChem.jl"
    ),
    pages = [
        "Home" => "index.md",
        "SuperFast" => "superfast.md",
        "GEOS-Chem" => [
            "Overview" => "geoschem/overview.md",
            "State Variables" => "geoschem/states.md",
            "Parameters" => "geoschem/params.md"
        ],
        "Seinfeld & Pandis Ch. 6" => [
            "OH Production" => "oh_production.md",
            "NOx Photochemistry" => "nox_photochemistry.md",
            "CO Oxidation" => "co_oxidation.md",
            "Methane Oxidation" => "methane_oxidation.md",
            "Combined System" => "combined_system.md"
        ],
        "Composing models" => "composing_models.md",
        "API" => "api.md"
    ]
)

deploydocs(; repo = "github.com/EarthSciML/GasChem.jl.git")

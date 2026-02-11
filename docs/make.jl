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
        repolink = "https://github.com/EarthSciML/GasChem.jl",
        size_threshold = 2 * 1024 * 1024,  # 2 MiB (params page is large due to many GEOS-Chem parameters)
        size_threshold_warn = 500 * 1024    # 500 KiB
    ),
    pages = [
        "Home" => "index.md",
        "Atmospheric Lifetime" => "AtmosphericLifetime.md",
        "SuperFast" => "superfast.md",
        "GEOS-Chem" => [
            "Overview" => "geoschem/overview.md",
            "State Variables" => "geoschem/states.md",
            "Parameters" => "geoschem/params.md",
        ],
        "Seinfeld & Pandis Ch. 6" => [
            "OH Production" => "oh_production.md",
            "NOx Photochemistry" => "nox_photochemistry.md",
            "CO Oxidation" => "co_oxidation.md",
            "Methane Oxidation" => "methane_oxidation.md",
            "Combined System" => "combined_system.md",
        ],
        "Radiation Fundamentals" => "radiation_fundamentals.md",
        "Stratospheric Chemistry" => "StratosphericChemistry.md",
        "Climate Forcing" => "climate_forcing.md",
        "Composing models" => "composing_models.md",
        "API" => "api.md",
    ]
)

deploydocs(; repo = "github.com/EarthSciML/GasChem.jl.git")

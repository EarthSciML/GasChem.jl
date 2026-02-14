```@meta
CurrentModule = GasChem
```

# GasChem: Gas-Phase Atmospheric Chemical Mechanisms

## Installation

To install GasChem.jl, use the Julia package manager:

```julia
using Pkg
Pkg.add("GasChem")
```

## Features

Currently, we have implemented versions of the SuperFast and GEOS-Chem chemical mechanisms, which can be optionally coupled with the Fast-JX photolysis model, which is also implemented here.

## Seinfeld & Pandis Chapter 6: Chemistry of the Troposphere

GasChem.jl also includes implementations of tropospheric chemistry systems from
Seinfeld & Pandis "Atmospheric Chemistry and Physics" (2nd Edition), Chapter 6:

  - **[OH Production](@ref oh_production)**: Production of OH radicals from O3 photolysis (Section 6.1)
  - **[NOx Photochemistry](@ref nox_photochemistry)**: The NOx photochemical cycle and photostationary state (Section 6.2)
  - **[CO Oxidation](@ref co_oxidation)**: CO oxidation and HOx cycling (Section 6.3)
  - **[Methane Oxidation](@ref methane_oxidation)**: Complete methane oxidation mechanism (Section 6.4)
  - **[Combined System](@ref combined_system)**: Integrated tropospheric chemistry model

### Key Concepts

Ozone in the troposphere is produced through the photochemical oxidation of
volatile organic compounds (VOCs) and CO in the presence of nitrogen oxides (NOx).
The OH and HO2 radicals (collectively HOx) cycle rapidly, allowing catalytic
ozone production before termination (mainly OH + NO2 -> HNO3).

Ozone production depends strongly on NOx levels:

  - **High NOx (VOC-limited)**: O3 increases with VOC, decreases with NOx
  - **Low NOx (NOx-limited)**: O3 increases with NOx

**Reference**: Seinfeld, J.H. and Pandis, S.N. (2006). *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change*, 2nd Edition. John Wiley & Sons. Chapter 6.

## Contributing

...coming soon

## Reproducibility

```@raw html
<details><summary>The documentation of this EarthSciML package was built using these direct dependencies,</summary>
```

```@example
using Pkg # hide
Pkg.status() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>and using this machine and Julia version.</summary>
```

```@example
using InteractiveUtils # hide
versioninfo() # hide
```

```@raw html
</details>
```

```@raw html
<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>
```

```@example
using Pkg # hide
Pkg.status(; mode = PKGMODE_MANIFEST) # hide
```

```@raw html
</details>
```

```@eval
using TOML
using Markdown
version = TOML.parse(read("../../Project.toml", String))["version"]
name = TOML.parse(read("../../Project.toml", String))["name"]
link_manifest = "https://github.com/EarthSciML/" * name * ".jl/tree/gh-pages/v" * version *
                "/assets/Manifest.toml"
link_project = "https://github.com/EarthSciML/" * name * ".jl/tree/gh-pages/v" * version *
               "/assets/Project.toml"
Markdown.parse("""You can also download the
[manifest]($link_manifest)
file and the
[project]($link_project)
file.
""")
```

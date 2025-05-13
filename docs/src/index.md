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

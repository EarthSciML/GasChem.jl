# TroposphericChemistry.jl

ModelingToolkit.jl implementations of tropospheric chemistry systems from
Seinfeld & Pandis "Atmospheric Chemistry and Physics" (2nd Edition), Chapter 6.

## Overview

This package implements the major chemical systems governing tropospheric ozone
formation as described in Chapter 6 of Seinfeld & Pandis (2006):

- **[OH Production](@ref oh_production)**: Production of OH radicals from O₃ photolysis (Section 6.1)
- **[NOx Photochemistry](@ref nox_photochemistry)**: The NOx photochemical cycle and photostationary state (Section 6.2)
- **[CO Oxidation](@ref co_oxidation)**: CO oxidation and HOx cycling (Section 6.3)
- **[Methane Oxidation](@ref methane_oxidation)**: Complete methane oxidation mechanism (Section 6.4)
- **[Combined System](@ref combined_system)**: Integrated tropospheric chemistry model

## Key Concepts

### Tropospheric Ozone Production

Ozone in the troposphere is produced through the photochemical oxidation of
volatile organic compounds (VOCs) and CO in the presence of nitrogen oxides (NOx).
The basic mechanism involves:

1. **OH Production**: OH radicals are produced primarily through O₃ photolysis
   at wavelengths < 320 nm, producing O(¹D) which reacts with H₂O.

2. **VOC/CO Oxidation**: OH initiates oxidation of CO and hydrocarbons,
   producing peroxy radicals (HO₂, RO₂).

3. **NO-NO₂ Conversion**: Peroxy radicals convert NO to NO₂ without consuming O₃.

4. **O₃ Production**: NO₂ photolyzes to produce O atoms, which combine with O₂
   to form O₃.

Net reaction: CO + 2O₂ + hν → CO₂ + O₃

### HOx Cycling

The OH and HO₂ radicals (collectively HOx) cycle rapidly:

```
OH → (+ CO, CH₄, VOC) → HO₂ → (+ NO) → OH
```

This catalytic cycle allows a single OH radical to produce multiple O₃ molecules
before termination (mainly OH + NO₂ → HNO₃).

### NOx Regimes

Ozone production depends strongly on NOx levels:

- **High NOx (VOC-limited)**: O₃ increases with VOC, decreases with NOx
- **Low NOx (NOx-limited)**: O₃ increases with NOx

The Ozone Production Efficiency (OPE = ΔO₃/ΔNOx) varies from ~1-3 at high NOx
to ~10-30 at low NOx.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/EarthSciML/TroposphericChemistry.jl")
```

## Quick Start

```julia
using TroposphericChemistry
using ModelingToolkit

# Create the combined system
@named sys = TroposphericChemistrySystem()

# Simplify the system
sys_simplified = structural_simplify(sys)

# Get typical atmospheric conditions
conditions = get_typical_conditions()
```

## Reference

Seinfeld, J.H. and Pandis, S.N. (2006). Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change, 2nd Edition. John Wiley & Sons.
Chapter 6: Chemistry of the Troposphere.

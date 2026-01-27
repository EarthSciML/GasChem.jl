# [Combined System](@id combined_system)

The combined tropospheric chemistry system integrates all individual mechanisms
into a comprehensive model.

## Overview

The `TroposphericChemistrySystem` couples:
- OH production from O₃ photolysis (Section 6.1)
- NOx photochemical cycle (Section 6.2)
- CO oxidation and HOx cycling (Section 6.3)
- Methane oxidation (Section 6.4)

## Coupling Points

### OH-HOx Coupling
- OH is produced from O₃ photolysis and HO₂ + NO
- OH is consumed by CO, CH₄, and other species
- HO₂ is produced from CO/VOC oxidation
- HO₂ is lost to NO (producing OH) and self-reaction

### NOx-HOx Coupling
- HO₂ + NO → OH + NO₂ converts NO to NO₂
- RO₂ + NO → RO + NO₂ converts NO to NO₂
- OH + NO₂ → HNO₃ removes both OH and NOx

### O₃ Budget
- Production: peroxy + NO reactions followed by NO₂ photolysis
- Loss: NO + O₃, OH + O₃, HO₂ + O₃

## Diagnostic Variables

| Variable | Description |
|----------|-------------|
| `NOx` | Total nitrogen oxides (NO + NO₂) |
| `HOx` | Total odd hydrogen (OH + HO₂) |
| `P_O3_total` | Total O₃ production rate |
| `L_O3_total` | Total O₃ loss rate |
| `P_O3_net` | Net O₃ tendency |
| `OPE` | Ozone production efficiency |
| `chain_length` | HOx chain length |
| `Φ` | Photostationary state parameter |

## Atmospheric Conditions

The module provides helper functions for typical conditions:

### Typical Background
```julia
conditions = get_typical_conditions()
# O₃ ~ 40 ppb, NO ~ 0.1 ppb, NO₂ ~ 1 ppb
```

### Urban
```julia
conditions = get_urban_conditions()
# O₃ ~ 80 ppb, NO ~ 10 ppb, NO₂ ~ 30 ppb
```

### Remote/Clean
```julia
conditions = get_remote_conditions()
# O₃ ~ 30 ppb, NO ~ 10 ppt, NO₂ ~ 20 ppt
```

## NOx-VOC Regimes

The combined system captures the NOx-VOC sensitivity of O₃ production:

### High NOx (VOC-limited)
- O₃ increases with VOC, decreases with NOx
- HOx termination: OH + NO₂ → HNO₃
- Short chain length (~1-5)
- Low OPE (~1-3)

### Low NOx (NOx-limited)
- O₃ increases with NOx
- HOx termination: HO₂ + HO₂ → H₂O₂
- Long chain length (~100+)
- High OPE (~10-30)

### O₃ Isopleth

The relationship between O₃, NOx, and VOC can be visualized as an isopleth
diagram showing constant O₃ contours in NOx-VOC space.

```
     VOC ↑
         │    ╱ High O₃
         │   ╱
         │  ╱ Ridge line
         │ ╱
    NOx- │╱   VOC-limited
    limited  ╱
         │  ╱
         └──────────→ NOx
```

## Usage

```julia
using TroposphericChemistry
using ModelingToolkit
using DifferentialEquations

# Create the combined system
@named sys = TroposphericChemistrySystem()

# Simplify
sys_simplified = structural_simplify(sys)

# Get typical conditions
conditions = get_typical_conditions()
```

## API

```@docs
TroposphericChemistrySystem
get_typical_conditions
get_urban_conditions
get_remote_conditions
```

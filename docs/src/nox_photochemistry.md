# [NOx Photochemistry](@id nox_photochemistry)

Implementation of the NOx photochemical cycle, based on Seinfeld & Pandis
Chapter 6, Section 6.2 (pp. 207-212).

## Background

The NOx (NO + NO₂) photochemical cycle is fundamental to tropospheric ozone
chemistry. In the absence of other species, the cycle involves:

```
NO₂ + hν → NO + O        (λ < 424 nm)     [R1]
O + O₂ + M → O₃ + M                        [R2]
NO + O₃ → NO₂ + O₂                         [R3]
```

These reactions cycle rapidly during the day, establishing a "photostationary
state" relationship between O₃, NO, and NO₂.

## Equations

### Equation 6.5: Ground-State Oxygen Atom Steady-State

```math
[O] = \frac{j_{NO_2}[NO_2]}{k_2[O_2][M]}
```

The concentration of ground-state O atoms is determined by the balance between
production from NO₂ photolysis and loss by reaction with O₂.

### Equation 6.6: Photostationary State (Leighton Relationship)

```math
[O_3] = \frac{j_{NO_2}[NO_2]}{k_{NO+O_3}[NO]}
```

This is the fundamental photostationary state relationship, also known as the
Leighton relationship. It predicts the ozone concentration when the NO-NO₂-O₃
system is in photochemical equilibrium.

### Equation 6.7: Photostationary State Parameter (Φ)

```math
\Phi = \frac{j_{NO_2}[NO_2]}{k_{NO+O_3}[NO][O_3]}
```

In perfect photostationary state, Φ = 1. Deviations indicate:
- **Φ > 1**: Net O₃ production (presence of HO₂, RO₂ converting NO to NO₂)
- **Φ < 1**: Net O₃ loss (additional NO sources, or nighttime)

Typical urban measurements show Φ = 1.5-3 during daytime.

### Equation 6.8: Net O₃ Production Rate

```math
\frac{d[O_3]}{dt} = j_{NO_2}[NO_2] - k_{NO+O_3}[NO][O_3]
```

When the system is not in photostationary state, this gives the net rate of
ozone change.

## Rate Constants at 298 K

| Reaction | Rate Constant | Units |
|----------|---------------|-------|
| NO₂ + hν → NO + O | ~8 × 10⁻³ | s⁻¹ (midday) |
| O + O₂ + M → O₃ + M | 6.0 × 10⁻³⁴ | cm⁶ molecule⁻² s⁻¹ |
| NO + O₃ → NO₂ + O₂ | 1.8 × 10⁻¹⁴ | cm³ molecule⁻¹ s⁻¹ |

## Photostationary State Example

For typical urban conditions at midday:
- j_NO₂ = 8 × 10⁻³ s⁻¹
- [NO₂] = 2 ppb = 5 × 10¹⁰ molecules cm⁻³
- [NO] = 1 ppb = 2.5 × 10¹⁰ molecules cm⁻³
- k_NO+O₃ = 1.8 × 10⁻¹⁴ cm³ s⁻¹

```math
[O_3]_{pss} = \frac{8 \times 10^{-3} \times 5 \times 10^{10}}{1.8 \times 10^{-14} \times 2.5 \times 10^{10}} \approx 9 \times 10^{11} \text{ molecules cm}^{-3} \approx 36 \text{ ppb}
```

## Usage

```julia
using TroposphericChemistry
using ModelingToolkit

# Create the NOx photochemistry system
@named nox_sys = NOxPhotochemistry()

# Or use the simplified photostationary state model
@named pss_sys = PhotostationaryState()
```

## API

```@docs
NOxPhotochemistry
PhotostationaryState
```

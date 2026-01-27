# [OH Production](@id oh_production)

Implementation of OH radical production from O₃ photolysis, based on
Seinfeld & Pandis Chapter 6, Section 6.1 (pp. 204-207).

## Background

The hydroxyl radical (OH) is the primary oxidant in the troposphere, initiating
the oxidation of most trace gases. The dominant source of OH in the troposphere
is photolysis of ozone at wavelengths below 320 nm:

```
O₃ + hν → O(¹D) + O₂     (λ < 320 nm)
```

The electronically excited oxygen atom O(¹D) can either be quenched back to
ground state O(³P) by collision with N₂ or O₂, or react with water vapor to
produce OH:

```
O(¹D) + M → O(³P) + M    (M = N₂, O₂)
O(¹D) + H₂O → 2OH
```

## Equations

### Equation 6.1: O(¹D) Steady-State Concentration

```math
[O(^1D)] = \frac{j_{O_3 \to O(^1D)} [O_3]}{k_3[M] + k_4[H_2O]}
```

where:
- ``j_{O_3 \to O(^1D)}`` = photolysis rate of O₃ producing O(¹D) [s⁻¹]
- ``k_3`` = quenching rate constant for O(¹D) + M [cm³ molecule⁻¹ s⁻¹]
- ``k_4`` = rate constant for O(¹D) + H₂O → 2OH [cm³ molecule⁻¹ s⁻¹]
- ``[M]`` = total air concentration [molecules cm⁻³]

### Equation 6.3: OH Production Rate

```math
P_{OH} = 2 j_{O_3 \to O(^1D)} \frac{k_4[H_2O]}{k_3[M]} [O_3]
```

This is the approximate form valid when ``k_3[M] \gg k_4[H_2O]``, which is
typically the case in the troposphere.

### Equation 6.4: OH Yield

```math
\varepsilon_{OH} = \frac{k_4[H_2O]}{k_3[M] + k_4[H_2O]}
```

The OH yield represents the fraction of O(¹D) atoms that react with H₂O to
produce OH rather than being quenched back to O(³P).

## Rate Constants at 298 K

| Reaction | Rate Constant | Units |
|----------|---------------|-------|
| O(¹D) + N₂ → O + N₂ | 2.6 × 10⁻¹¹ | cm³ molecule⁻¹ s⁻¹ |
| O(¹D) + O₂ → O + O₂ | 4.0 × 10⁻¹¹ | cm³ molecule⁻¹ s⁻¹ |
| O(¹D) + H₂O → 2OH | 2.2 × 10⁻¹⁰ | cm³ molecule⁻¹ s⁻¹ |

## Typical Values

For typical tropospheric conditions:
- O₃ photolysis rate: ``j_{O_3} \approx 10^{-5}`` s⁻¹
- OH yield: ``\varepsilon_{OH} \approx 0.05-0.15`` (5-15%)
- OH production rate: ``P_{OH} \approx 10^6`` molecules cm⁻³ s⁻¹

## Usage

```julia
using TroposphericChemistry
using ModelingToolkit

# Create the OH production system
@named oh_sys = OHProduction()

# Simplify
oh_simplified = structural_simplify(oh_sys)
```

## API

```@docs
OHProduction
```

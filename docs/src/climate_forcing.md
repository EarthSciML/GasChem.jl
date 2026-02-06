# Climate Forcing and Global Warming Potentials

## Overview

This module implements the fundamental equations for climate sensitivity with feedbacks, radiative forcing from greenhouse gases, and Global Warming Potentials (GWPs) as presented in Chapter 23 "Climate and Chemical Composition of the Atmosphere" of Seinfeld & Pandis (2006).

**Reference**: Seinfeld, J.H. and Pandis, S.N. (2006) *Atmospheric Chemistry and Physics: From Air Pollution to Climate Change*, 2nd Edition, Chapter 23, John Wiley & Sons.

```@docs
ClimateFeedback
GHGForcing
GlobalWarmingPotential
GWP_exponential
```

## Implementation

### Climate Feedback

The climate feedback equations relate radiative forcing perturbations to equilibrium temperature changes with feedbacks (Eqs. 23.1-23.4, and the unrealized warming relation from p. 1045). Note: The basic no-feedback climate sensitivity is implemented in [`ClimateSensitivity`](@ref) (Chapter 4); this component extends it by incorporating feedback processes.

**Eq. 23.1** - Climate sensitivity with feedbacks:

```math
\Delta T_s = \lambda \Delta F
```

**Eq. 23.2** - No-feedback temperature response:

```math
\Delta T_0 = \lambda_0 \Delta F
```

**Eq. 23.3** - Efficacy of forcing agents:

```math
E_i = \frac{\lambda_i}{\lambda_{CO_2}}
```

**Eq. 23.4** - Effective forcing (incorporating efficacy):

```math
\Delta F_e = \Delta F_i \cdot E_i
```

**Unrealized (committed) warming** (p. 1045):

```math
\Delta T_{\text{unrealized}} = (\Delta F - \Delta F_r) \cdot \lambda
```

#### State Variables

```@example climate
using GasChem, ModelingToolkit, Symbolics
using DynamicQuantities
using DataFrames

sys = ClimateFeedback()

vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars])
```

#### Parameters

```@example climate
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in params],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in params],
    :Description => [ModelingToolkit.getdescription(v) for v in params],
    :Default => [ModelingToolkit.hasdefault(v) ? ModelingToolkit.getdefault(v) : missing for v in params])
```

#### Equations

```@example climate
eqs = equations(sys)
```

### Greenhouse Gas Radiative Forcing

The `GHGForcing` component computes radiative forcing contributions from the major well-mixed greenhouse gases using reference values from IPCC (2001), as cited in Seinfeld & Pandis (2006), p. 1039.

| Species         | Forcing (W/m²) |
|:--------------- |:-------------- |
| CO₂             | 1.46           |
| CH₄             | 0.48           |
| N₂O             | 0.15           |
| Tropospheric O₃ | 0.40           |
| Halocarbons     | 0.34           |
| **Total**       | **2.83**       |

Note: The well-mixed GHG total (CO₂ + CH₄ + N₂O + halocarbons) is 2.43 W/m² as stated in the text; the total of 2.83 W/m² includes tropospheric O₃ forcing.

#### State Variables

```@example climate
ghg = GHGForcing()
vars = unknowns(ghg)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars])
```

#### Parameters

```@example climate
params = parameters(ghg)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in params],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in params],
    :Description => [ModelingToolkit.getdescription(v) for v in params],
    :Default => [ModelingToolkit.hasdefault(v) ? ModelingToolkit.getdefault(v) : missing for v in params])
```

#### Equations

```@example climate
eqs = equations(ghg)
```

### Global Warming Potential

The Global Warming Potential (GWP) is defined as the time-integrated radiative forcing from a pulse emission of 1 kg of a compound relative to 1 kg of CO₂ over a specified time horizon.

**Eq. 23.5** - GWP definition:

```math
\text{GWP} = \frac{\int_0^{t_f} a_A [A(t)] dt}{\int_0^{t_f} a_R [R(t)] dt}
```

**Eq. 23.6** - Absolute GWP:

```math
\text{AGWP}_A = \int_0^{t_f} a_A [A(t)] dt
```

For exponentially decaying species (Problem 23.3):

```math
\text{GWP} = \frac{a_A \tau_A (1 - e^{-t_f/\tau_A})}{a_{CO_2} \tau_{CO_2} (1 - e^{-t_f/\tau_{CO_2}})}
```

#### State Variables

```@example climate
gwp_sys = GlobalWarmingPotential()
vars = unknowns(gwp_sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars])
```

#### Parameters

```@example climate
params = parameters(gwp_sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in params],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in params],
    :Description => [ModelingToolkit.getdescription(v) for v in params],
    :Default => [ModelingToolkit.hasdefault(v) ? ModelingToolkit.getdefault(v) : missing for v in params])
```

#### Equations

```@example climate
eqs = equations(gwp_sys)
```

## Analysis

### GWP Calculations

Using the `GWP_exponential` function to calculate GWPs for different species:

```@example climate
using GasChem, Plots

# Species parameters from Table 23.1
# τ = atmospheric lifetime (years)
# a = radiative efficiency relative to CO₂

# Calculate GWP for CH₄ (τ = 12 yr)
# Using fitted a value to match tabulated GWP₁₀₀ = 23
gwp_ch4_100 = GWP_exponential(12.0, 140.0, 100.0)
println("CH₄ GWP₁₀₀ = $(round(gwp_ch4_100, digits=1))")

# N₂O (τ = 114 yr)
gwp_n2o_100 = GWP_exponential(114.0, 326.0, 100.0)
println("N₂O GWP₁₀₀ = $(round(gwp_n2o_100, digits=1))")

# CO₂ relative to itself
gwp_co2 = GWP_exponential(150.0, 1.0, 100.0)
println("CO₂ GWP₁₀₀ = $(round(gwp_co2, digits=1))")
```

### GWP Time Horizon Dependence (Figure 23.15)

GWP values depend strongly on the time horizon chosen. Short-lived species have higher GWPs at short time horizons, while long-lived species have relatively higher GWPs at longer horizons. This reproduces the behavior shown in Figure 23.15 of Seinfeld & Pandis (2006), which uses log-log axes.

```@example climate
# Calculate GWP vs time horizon for different species
t_horizons = 1:1:500

# HCFC-225ca-like (τ = 2.5 yr)
gwp_hcfc225 = [GWP_exponential(2.5, 100.0, Float64(t)) for t in t_horizons]

# HFC-134a-like (τ = 1.4 yr, high radiative efficiency)
gwp_hfc134a = [GWP_exponential(1.4, 3400.0, Float64(t)) for t in t_horizons]

# N₂O-like (τ = 120 yr)
gwp_n2o = [GWP_exponential(120.0, 326.0, Float64(t)) for t in t_horizons]

# C₂F₆-like (τ = 10000 yr)
gwp_c2f6 = [GWP_exponential(10000.0, 1000.0, Float64(t)) for t in t_horizons]

p = plot(t_horizons, gwp_c2f6, label = "C₂F₆ (τ~10000 yr)", lw = 2,
    xscale = :log10, yscale = :log10)
plot!(p, t_horizons, gwp_hfc134a, label = "HFC-134a (τ~1.4 yr)", lw = 2)
plot!(p, t_horizons, gwp_n2o, label = "N₂O (τ~120 yr)", lw = 2)
plot!(p, t_horizons, gwp_hcfc225, label = "HCFC-225ca (τ~2.5 yr)", lw = 2)
xlabel!(p, "Time Horizon (years)")
ylabel!(p, "Global Warming Potential")
title!(p, "GWP vs Time Horizon (cf. Figure 23.15)")
p
```

### Climate Feedback Factor

The climate feedback factor is defined as the ratio of the climate sensitivity with feedbacks to the no-feedback sensitivity:

```math
\text{Feedback Factor} = \frac{\lambda}{\lambda_0}
```

From page 1040 of Seinfeld & Pandis (2006), typical values range from 1.2 to 3.75, corresponding to equilibrium climate sensitivities for CO₂ doubling (ΔT₂ₓCO₂) of approximately 1.5 to 4.5 K.

### Unrealized Warming

The Earth system has thermal inertia, primarily due to the oceans, which means that the full warming response to current forcing levels has not yet been realized. The unrealized (or committed) warming represents the additional temperature increase that will occur even if forcing is held constant.

From page 1045 of Seinfeld & Pandis (2006):

  - Current forcing: ΔF ≈ 1.7 W/m²
  - Realized warming: ΔT_r ≈ 0.7 K
  - With λ ≈ 0.7 K/(W/m²): unrealized warming ≈ 0.5 K

This means even with no additional emissions, the Earth is committed to approximately 0.5 K of additional warming.

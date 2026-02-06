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

The climate feedback equations relate radiative forcing perturbations to equilibrium temperature changes with feedbacks (Eqs. 23.1-23.4 and 23.7). Note: The basic no-feedback climate sensitivity is implemented in [`ClimateSensitivity`](@ref) (Chapter 4); this component extends it by incorporating feedback processes.

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

**Eq. 23.7** - Unrealized (committed) warming (p. 1045):

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
    :Default => [ModelingToolkit.hasdefault(v) ? ModelingToolkit.getdefault(v) : missing
                 for v in params])
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
    :Default => [ModelingToolkit.hasdefault(v) ? ModelingToolkit.getdefault(v) : missing
                 for v in params])
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
    :Default => [ModelingToolkit.hasdefault(v) ? ModelingToolkit.getdefault(v) : missing
                 for v in params])
```

#### Equations

```@example climate
eqs = equations(gwp_sys)
```

## Analysis

### Table 23.1: Greenhouse Gas Properties and GWPs

Table 23.1 from Seinfeld & Pandis (2006) lists key greenhouse gases with their abundances, lifetimes, and 100-year GWPs. Using the exponential decay GWP formula (Problem 23.3), we can reproduce the tabulated GWP values by fitting the relative radiative efficiency parameter `a`. The table below compares computed GWPs against the textbook values for selected species.

```@example climate
using GasChem, DataFrames

# Species data from Table 23.1 (Seinfeld & Pandis 2006, p. 1044)
# τ = perturbation lifetime (yr), GWP_ref = tabulated 100-yr GWP
# a = relative radiative efficiency (fitted to reproduce GWP_ref)
species_data = [
    ("CH₄", 12.0, 23, 140.0),
    ("N₂O", 114.0, 296, 326.0),
    ("CF₄", 50000.0, 5700, 390.0),
    ("C₂F₆", 10000.0, 11900, 816.0),
    ("SF₆", 3200.0, 22200, 1540.0),
    ("HFC-134a", 13.8, 1300, 692.0),
    ("CFC-11", 45.0, 4600, 753.0),
    ("CFC-12", 100.0, 10600, 1162.0),
    ("CCl₄", 35.0, 1800, 385.0)
]

names_col = [s[1] for s in species_data]
τ_col = [s[2] for s in species_data]
gwp_ref_col = [s[3] for s in species_data]
gwp_calc_col = [round(GWP_exponential(s[2], s[4], 100.0), digits = 0)
                for s in species_data]

DataFrame(
    :Species => names_col,
    Symbol("Lifetime (yr)") => τ_col,
    Symbol("GWP₁₀₀ (Table 23.1)") => gwp_ref_col,
    Symbol("GWP₁₀₀ (Computed)") => gwp_calc_col)
```

### GWP Calculations

Using the `GWP_exponential` function to calculate GWPs for different species:

```@example climate
using Plots

# Species parameters from Table 23.1
# τ = atmospheric lifetime (years)
# a = radiative efficiency relative to CO₂

# Calculate GWP for CH₄ (τ = 12 yr)
# Using fitted a value to match tabulated GWP₁₀₀ = 23
gwp_ch4_100 = GWP_exponential(12.0, 140.0, 100.0)
println("CH₄ GWP₁₀₀ = $(round(gwp_ch4_100, digits = 1))")

# N₂O (τ = 114 yr)
gwp_n2o_100 = GWP_exponential(114.0, 326.0, 100.0)
println("N₂O GWP₁₀₀ = $(round(gwp_n2o_100, digits = 1))")

# CO₂ relative to itself
gwp_co2 = GWP_exponential(150.0, 1.0, 100.0)
println("CO₂ GWP₁₀₀ = $(round(gwp_co2, digits = 1))")
```

### Radiative Forcing After Pulse Emission (Figure 23.13)

Figure 23.13 shows the instantaneous radiative forcing (W m⁻² kg⁻¹) as a function of time after a pulse release of 1 kg of several greenhouse gases (IPCC 1995). For exponentially decaying species, the forcing decays as ``a_A \exp(-t/\tau_A)``. Note that after 1 year, instantaneous forcings can differ by four orders of magnitude.

```@example climate
# Reproduce Figure 23.13: Radiative forcing vs time after pulse release
# For exponentially decaying species: F(t) = a_A × exp(-t/τ_A)
# Using species lifetimes as labeled in Figure 23.13
t_years = 1:1:1000

# Radiative forcing decay: a_A × exp(-t/τ_A)
# a values are in arbitrary units relative to each other (qualitative reproduction)
forcing_c2f6 = [816.0 * exp(-Float64(t) / 10000.0) for t in t_years]
forcing_n2o = [326.0 * exp(-Float64(t) / 120.0) for t in t_years]
forcing_co2 = [1.0 * exp(-Float64(t) / 150.0) for t in t_years]
forcing_hcfc225 = [100.0 * exp(-Float64(t) / 2.5) for t in t_years]
forcing_hfc134a = [692.0 * exp(-Float64(t) / 1.4) for t in t_years]

p = plot(t_years, forcing_c2f6, label = "C₂F₆ (τ≈10000 yr)", lw = 2,
    xscale = :log10, yscale = :log10,
    ylim = (1e-8, 1e4))
plot!(p, t_years, forcing_n2o, label = "N₂O (τ≈120 yr)", lw = 2)
plot!(p, t_years, forcing_co2, label = "CO₂ (τ≈150 yr)", lw = 2)
plot!(p, t_years, forcing_hcfc225, label = "HCFC-225ca (τ≈2.5 yr)", lw = 2)
plot!(p, t_years, forcing_hfc134a, label = "HFC-134a (τ≈1.4 yr)", lw = 2)
xlabel!(p, "Time (years)")
ylabel!(p, "Instantaneous Radiative Forcing (relative)")
title!(p, "Radiative Forcing After Pulse Emission (cf. Figure 23.13)")
p
```

### Absolute GWP vs Time Horizon (Figure 23.14)

Figure 23.14 shows the time-integrated radiative forcing (absolute GWP, W m⁻² kg⁻¹ yr) for a range of greenhouse gases, computed from the pulse-decay curves of Figure 23.13 (IPCC 1995). The AGWP is the numerator of the GWP expression (Eq. 23.6).

```@example climate
# Reproduce Figure 23.14: Absolute GWP (time-integrated forcing) vs time horizon
# AGWP_A = a_A × τ_A × (1 - exp(-t_f/τ_A))
t_horizons_long = 1:1:1000

agwp_c2f6 = [816.0 * 10000.0 * (1 - exp(-Float64(t) / 10000.0)) for t in t_horizons_long]
agwp_hfc134a = [692.0 * 1.4 * (1 - exp(-Float64(t) / 1.4)) for t in t_horizons_long]
agwp_hcfc225 = [100.0 * 2.5 * (1 - exp(-Float64(t) / 2.5)) for t in t_horizons_long]
agwp_n2o = [326.0 * 120.0 * (1 - exp(-Float64(t) / 120.0)) for t in t_horizons_long]
agwp_co2 = [1.0 * 150.0 * (1 - exp(-Float64(t) / 150.0)) for t in t_horizons_long]

p = plot(t_horizons_long, agwp_c2f6, label = "C₂F₆ (τ≈10000 yr)", lw = 2,
    xscale = :log10, yscale = :log10)
plot!(p, t_horizons_long, agwp_hfc134a, label = "HFC-134a (τ≈1.4 yr)", lw = 2)
plot!(p, t_horizons_long, agwp_hcfc225, label = "HCFC-225ca (τ≈2.5 yr)", lw = 2)
plot!(p, t_horizons_long, agwp_n2o, label = "N₂O (τ≈120 yr)", lw = 2)
plot!(p, t_horizons_long, agwp_co2, label = "CO₂ (τ≈150 yr)", lw = 2)
xlabel!(p, "Time Horizon (years)")
ylabel!(p, "Absolute Global Warming Potential (relative)")
title!(p, "Absolute GWP vs Time Horizon (cf. Figure 23.14)")
p
```

### GWP Time Horizon Dependence (Figure 23.15)

GWP values depend strongly on the time horizon chosen. Short-lived species have higher GWPs at short time horizons, while long-lived species have relatively higher GWPs at longer horizons. This reproduces the behavior shown in Figure 23.15 of Seinfeld & Pandis (2006), which uses log-log axes.

```@example climate
# Calculate GWP vs time horizon for species from Figure 23.15 (IPCC 1995)
# Using lifetimes as labeled in the figure
t_horizons = 1:1:500

# C₂F₆ (τ ~ 10,000 yr, as labeled in Figure 23.15)
gwp_c2f6 = [GWP_exponential(10000.0, 816.0, Float64(t)) for t in t_horizons]

# HFC-134a (τ ~ 1.4 yr, as labeled in Figure 23.15)
gwp_hfc134a = [GWP_exponential(1.4, 692.0, Float64(t)) for t in t_horizons]

# N₂O (τ ~ 120 yr, as labeled in Figure 23.15)
gwp_n2o = [GWP_exponential(120.0, 326.0, Float64(t)) for t in t_horizons]

# HCFC-225ca (τ ~ 2.5 yr, as labeled in Figure 23.15)
gwp_hcfc225 = [GWP_exponential(2.5, 100.0, Float64(t)) for t in t_horizons]

p = plot(t_horizons, gwp_c2f6, label = "C₂F₆ (τ≈10000 yr)", lw = 2,
    xscale = :log10, yscale = :log10)
plot!(p, t_horizons, gwp_hfc134a, label = "HFC-134a (τ≈1.4 yr)", lw = 2)
plot!(p, t_horizons, gwp_n2o, label = "N₂O (τ≈120 yr)", lw = 2)
plot!(p, t_horizons, gwp_hcfc225, label = "HCFC-225ca (τ≈2.5 yr)", lw = 2)
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

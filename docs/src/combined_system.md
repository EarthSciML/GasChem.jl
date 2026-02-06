# [Combined System](@id combined_system)

## Overview

The `TroposphericChemistrySystem` integrates the individual chemistry mechanisms
from Chapter 6 into a comprehensive tropospheric chemistry diagnostic model.
It couples:

  - OH production from O3 photolysis (Section 6.1, via `OHProduction`)
  - NOx photochemical cycle (Section 6.2, via `NOxPhotochemistry`)
  - CO oxidation and HOx cycling (Section 6.3, via `COOxidation`)

The systems are coupled through shared species (OH, HO2, NO, NO2, O3) and
the combined system computes aggregate diagnostics including net O3 production,
OPE, and HOx chain length.

Condition systems provide typical atmospheric conditions for different environments.

**Reference**: Seinfeld, J.H. and Pandis, S.N. (2006). *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change*, 2nd Edition. John Wiley & Sons. Chapter 6.

```@docs
TroposphericChemistrySystem
```

```@docs
TypicalConditions
```

```@docs
UrbanConditions
```

```@docs
RemoteConditions
```

```@docs
get_conditions_dict
```

## Implementation

### State Variables

```@example combined
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities, GasChem

sys = TroposphericChemistrySystem()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example combined
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example combined
eqs = equations(sys)
```

### Typical Atmospheric Conditions

The module provides condition systems for three atmospheric environments:

```@example combined
using DataFrames

conditions = [
    ("Background", get_conditions_dict(TypicalConditions())),
    ("Urban", get_conditions_dict(UrbanConditions())),
    ("Remote", get_conditions_dict(RemoteConditions()))
]

# Show key species for each condition
species = [:O3, :NO, :NO2, :CO, :OH, :HO2]
header = vcat([:Species, :Units], [Symbol(c[1]) for c in conditions])

rows = []
for s in species
    ppb_factor = 2.5e16  # m⁻³ per ppb at STP
    vals = [c[2][s] for c in conditions]
    ppb_vals = vals ./ ppb_factor
    push!(rows,
        (
            Species = string(s),
            Units = "ppb",
            Background = round(ppb_vals[1], sigdigits = 3),
            Urban = round(ppb_vals[2], sigdigits = 3),
            Remote = round(ppb_vals[3], sigdigits = 3)
        ))
end
DataFrame(rows)
```

## Analysis

### O3 Production in Different NOx Regimes

The combined system captures the transition between NOx-limited and
VOC-limited regimes. In the NOx-limited regime (remote conditions),
O3 production increases linearly with NOx. In the VOC-limited regime
(urban conditions), adding more NOx can actually decrease O3.

This analysis uses the `TroposphericChemistrySystem` to compute how
net O3 production and OPE vary across a range of NOx levels.

```@example combined
using Plots

# Rate constants at 298 K (SI: m³/s) — same defaults as TroposphericChemistrySystem parameters
k_HO2_NO = 8.1e-12 * 1e-6
k_HO2_HO2 = 2.9e-12 * 1e-6
k_CO_OH = 2.4e-13 * 1e-6
k_CH3O2_NO = 7.7e-12 * 1e-6
k_OH_NO2 = 1.0e-11 * 1e-6
k_NO_O3 = 1.9e-14 * 1e-6
k_OH_O3 = 7.3e-14 * 1e-6
k_HO2_O3 = 2.0e-15 * 1e-6

# Fixed conditions (background, SI: m⁻³)
CO_val = 2.5e18     # 100 ppb
O3_val = 1e18       # 40 ppb
OH_val = 1e12       # typical daytime
CH3O2_val = 1e14    # typical
CH4_val = 4.5e19    # ~1800 ppb
H2O_val = 4e23
M_val = 2.5e25
O2_val = 5.25e24
P_OH = 1e12         # m⁻³/s (for HO2 estimation)

# Vary NO from 10 ppt to 100 ppb
NO_ppb = 10 .^ range(-2, 2, length = 300)
NO_vals = NO_ppb .* 2.5e16  # m⁻³
NO2_vals = 2 .* NO_vals     # assume NO2/NO ratio ~ 2

# Estimate HO2 from steady state in two regimes
HO2_high_NOx = k_CO_OH .* CO_val .* OH_val ./ (k_HO2_NO .* NO_vals)
HO2_low_NOx = sqrt(P_OH / (2 * k_HO2_HO2))
HO2_vals = min.(HO2_high_NOx, HO2_low_NOx)

# Compute diagnostics using rate constants from the combined system
P_O3_total = k_HO2_NO .* HO2_vals .* NO_vals .+ k_CH3O2_NO .* CH3O2_val .* NO_vals
L_O3_total = k_NO_O3 .* NO_vals .* O3_val .+ k_OH_O3 .* OH_val .* O3_val .+ k_HO2_O3 .* HO2_vals .* O3_val
P_O3_net_vals = P_O3_total .- L_O3_total
L_NOx = k_OH_NO2 .* OH_val .* NO2_vals
OPE_vals = P_O3_total ./ L_NOx

p1 = plot(NO_ppb, P_O3_net_vals ./ 1e12,
    xlabel = "NO (ppb)",
    ylabel = "Net P(O₃) (10¹² m⁻³ s⁻¹)",
    title = "Net O₃ Production vs NOx",
    xscale = :log10, linewidth = 2,
    label = "P(O₃)_net", legend = :topleft)
vline!([0.1], linestyle = :dash, color = :gray, label = "Background NO")
vline!([10.0], linestyle = :dash, color = :red, label = "Urban NO")

p2 = plot(NO_ppb, OPE_vals,
    xlabel = "NO (ppb)",
    ylabel = "OPE (mol O₃ / mol NOx)",
    title = "Ozone Production Efficiency",
    xscale = :log10, yscale = :log10, linewidth = 2,
    label = "OPE", legend = :topright, ylims = (0.1, 1000))
hline!([3.0], linestyle = :dot, color = :red, label = "Urban OPE ~ 1-3")
hline!([15.0], linestyle = :dot, color = :blue, label = "Remote OPE ~ 10-30")

plot(p1, p2, layout = (1, 2), size = (900, 400))
savefig("combined_o3_regimes.svg") # hide
```

![O3 production in different NOx regimes](combined_o3_regimes.svg)

The left panel shows that net O3 production peaks at intermediate NOx levels
and decreases at very high NOx due to O3 titration by NO and reduced HO2
concentrations. The right panel shows that OPE decreases monotonically with
increasing NOx, from OPE > 10 in remote conditions to OPE of 1-3 in urban
environments, consistent with Section 6.3 of Seinfeld & Pandis.

### Comparison of Atmospheric Conditions

This table computes key diagnostics for each of the three standard atmospheric
environments using the rate constants from the `TroposphericChemistrySystem`.

```@example combined
environments = [
    ("Background", get_conditions_dict(TypicalConditions())),
    ("Urban", get_conditions_dict(UrbanConditions())),
    ("Remote", get_conditions_dict(RemoteConditions()))
]

results = []
for (name, cond) in environments
    NO_v = cond[:NO]
    NO2_v = cond[:NO2]
    OH_v = cond[:OH]
    HO2_v = cond[:HO2]
    O3_v = cond[:O3]
    CH3O2_v = cond[:CH3O2]

    p_o3 = k_HO2_NO * HO2_v * NO_v + k_CH3O2_NO * CH3O2_v * NO_v
    l_nox = k_OH_NO2 * OH_v * NO2_v
    ope = p_o3 / l_nox
    l_hox = k_OH_NO2 * OH_v * NO2_v + 2 * k_HO2_HO2 * HO2_v^2
    chain = (k_HO2_NO * HO2_v * NO_v) / l_hox

    push!(results,
        (
            Environment = name,
            NO_ppb = round(NO_v / 2.5e16, sigdigits = 3),
            O3_ppb = round(O3_v / 2.5e16, sigdigits = 3),
            P_O3 = round(p_o3, sigdigits = 3),
            OPE = round(ope, sigdigits = 3),
            Chain_Length = round(chain, sigdigits = 3)
        ))
end

DataFrame(results)
```

This table shows how the key photochemical diagnostics vary across
different atmospheric environments, illustrating the transition from
NOx-limited (remote, high OPE and long chain length) to VOC-limited
(urban, low OPE and short chain length) regimes discussed in Chapter 6.

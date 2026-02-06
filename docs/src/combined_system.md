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

Helper functions provide typical atmospheric conditions for different environments.

**Reference**: Seinfeld, J.H. and Pandis, S.N. (2006). *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change*, 2nd Edition. John Wiley & Sons. Chapter 6.

```@docs
TroposphericChemistrySystem
```

```@docs
get_typical_conditions
```

```@docs
get_urban_conditions
```

```@docs
get_remote_conditions
```

## Implementation

### State Variables

```@example combined
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities, GasChem

sys = TroposphericChemistrySystem()
vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape=false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

### Parameters

```@example combined
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape=false)) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

### Equations

```@example combined
eqs = equations(sys)
```

### Typical Atmospheric Conditions

The module provides helper functions for three atmospheric environments:

```@example combined
using DataFrames

conditions = [
    ("Background", get_typical_conditions()),
    ("Urban", get_urban_conditions()),
    ("Remote", get_remote_conditions()),
]

# Show key species for each condition
species = [:O3, :NO, :NO2, :CO, :OH, :HO2]
header = vcat([:Species, :Units], [Symbol(c[1]) for c in conditions])

rows = []
for s in species
    ppb_factor = 2.5e16  # m⁻³ per ppb at STP
    vals = [c[2][s] for c in conditions]
    ppb_vals = vals ./ ppb_factor
    push!(rows, (
        Species = string(s),
        Units = "ppb",
        Background = round(ppb_vals[1], sigdigits=3),
        Urban = round(ppb_vals[2], sigdigits=3),
        Remote = round(ppb_vals[3], sigdigits=3),
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

This analysis uses the rate constants from the combined system to show
how net O3 production and OPE vary across a range of NOx levels.

```@example combined
using Plots

# Rate constants at 298 K (SI: m³/s)
k_HO2_NO = 8.1e-12 * 1e-6   # HO2 + NO [m³/s]
k_CH3O2_NO = 7.7e-12 * 1e-6 # CH3O2 + NO [m³/s]
k_OH_NO2 = 1.0e-11 * 1e-6   # OH + NO2 [m³/s]
k_NO_O3 = 1.8e-14 * 1e-6    # NO + O3 [m³/s]
k_OH_O3 = 7.3e-14 * 1e-6    # OH + O3 [m³/s]
k_HO2_O3 = 2.0e-15 * 1e-6   # HO2 + O3 [m³/s]
k_HO2_HO2 = 2.9e-12 * 1e-6  # HO2 + HO2 [m³/s]
k_CO_OH = 2.4e-13 * 1e-6     # CO + OH [m³/s]

# Fixed conditions (background, SI: m⁻³)
CO = 2.5e18     # 100 ppb
O3 = 1e18       # 40 ppb
OH = 1e12       # typical daytime
CH3O2 = 1e14    # typical
P_OH = 1e12     # m⁻³/s

# Vary NO from 10 ppt to 100 ppb
NO_ppb = 10 .^ range(-2, 2, length=300)
NO = NO_ppb .* 2.5e16  # m⁻³
NO2 = 2 .* NO  # assume NO2/NO ratio ~ 2

# Estimate HO2 from steady state in two regimes
HO2_high_NOx = k_CO_OH .* CO .* OH ./ (k_HO2_NO .* NO)
HO2_low_NOx = sqrt(P_OH / (2 * k_HO2_HO2))
HO2 = min.(HO2_high_NOx, HO2_low_NOx)

# Total O3 production (from HO2 + NO and CH3O2 + NO)
P_O3_total = k_HO2_NO .* HO2 .* NO .+ k_CH3O2_NO .* CH3O2 .* NO

# Total O3 loss
L_O3_total = k_NO_O3 .* NO .* O3 .+ k_OH_O3 .* OH .* O3 .+ k_HO2_O3 .* HO2 .* O3

# Net O3 production
P_O3_net = P_O3_total .- L_O3_total

# NOx loss and OPE
L_NOx = k_OH_NO2 .* OH .* NO2
OPE = P_O3_total ./ L_NOx

p1 = plot(NO_ppb, P_O3_net ./ 1e12,
    xlabel="NO (ppb)",
    ylabel="Net P(O₃) (10¹² m⁻³ s⁻¹)",
    title="Net O₃ Production vs NOx",
    xscale=:log10, linewidth=2,
    label="P(O₃)_net", legend=:topleft)
vline!([0.1], linestyle=:dash, color=:gray, label="Background NO")
vline!([10.0], linestyle=:dash, color=:red, label="Urban NO")

p2 = plot(NO_ppb, OPE,
    xlabel="NO (ppb)",
    ylabel="OPE (mol O₃ / mol NOx)",
    title="Ozone Production Efficiency",
    xscale=:log10, yscale=:log10, linewidth=2,
    label="OPE", legend=:topright, ylims=(0.1, 1000))
hline!([3.0], linestyle=:dot, color=:red, label="Urban OPE ~ 1-3")
hline!([15.0], linestyle=:dot, color=:blue, label="Remote OPE ~ 10-30")

plot(p1, p2, layout=(1, 2), size=(900, 400))
savefig("combined_o3_regimes.svg") # hide
```

![O3 production in different NOx regimes](combined_o3_regimes.svg)

The left panel shows that net O3 production peaks at intermediate NOx levels
and decreases at very high NOx due to O3 titration by NO and reduced HO2
concentrations. The right panel shows that OPE decreases monotonically with
increasing NOx, from OPE > 10 in remote conditions to OPE of 1-3 in urban
environments, consistent with Section 6.3 of Seinfeld & Pandis.

### Comparison of Atmospheric Conditions

```@example combined
# Compare key diagnostics across the three environments
environments = [
    ("Background", get_typical_conditions()),
    ("Urban", get_urban_conditions()),
    ("Remote", get_remote_conditions()),
]

# Compute diagnostics for each environment
results = []
for (name, cond) in environments
    NO_val = cond[:NO]
    NO2_val = cond[:NO2]
    OH_val = cond[:OH]
    HO2_val = cond[:HO2]
    O3_val = cond[:O3]
    CH3O2_val = cond[:CH3O2]

    p_o3 = k_HO2_NO * HO2_val * NO_val + k_CH3O2_NO * CH3O2_val * NO_val
    l_nox = k_OH_NO2 * OH_val * NO2_val
    ope = p_o3 / l_nox
    l_hox = k_OH_NO2 * OH_val * NO2_val + 2 * k_HO2_HO2 * HO2_val^2
    chain = (k_HO2_NO * HO2_val * NO_val) / l_hox

    push!(results, (
        Environment = name,
        NO_ppb = round(NO_val / 2.5e16, sigdigits=3),
        O3_ppb = round(O3_val / 2.5e16, sigdigits=3),
        P_O3 = round(p_o3, sigdigits=3),
        OPE = round(ope, sigdigits=3),
        Chain_Length = round(chain, sigdigits=3),
    ))
end

DataFrame(results)
```

This table shows how the key photochemical diagnostics vary across
different atmospheric environments, illustrating the transition from
NOx-limited (remote, high OPE and long chain length) to VOC-limited
(urban, low OPE and short chain length) regimes discussed in Chapter 6.

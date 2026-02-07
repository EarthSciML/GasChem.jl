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
net O3 production and OPE vary across a range of NOx levels. The HO2
concentration is estimated from steady-state approximations (Eqs. 6.13
and 6.18), and the resulting concentrations are fed into the compiled
system to compute the diagnostics.

```@example combined
using Plots, NonlinearSolve

# Compile the TroposphericChemistrySystem with all species as inputs
sys_nns = ModelingToolkit.toggle_namespacing(sys, false)
input_vars = [sys_nns.O3, sys_nns.NO, sys_nns.NO2, sys_nns.OH, sys_nns.HO2,
    sys_nns.CO, sys_nns.CH3O2, sys_nns.H2O, sys_nns.M, sys_nns.O2]
compiled = mtkcompile(sys; inputs = input_vars)

# Extract rate constants from the system parameters for HO2 estimation
k_HO2_NO_val = Float64(ModelingToolkit.getdefault(sys.co.k_HO2_NO))
k_HO2_HO2_val = Float64(ModelingToolkit.getdefault(sys.co.k_HO2_HO2))
k_CO_OH_val = Float64(ModelingToolkit.getdefault(sys.co.k_CO_OH))

# Fixed conditions (background, SI: m⁻³)
CO_val = 2.5e18     # 100 ppb
O3_val = 1e18       # 40 ppb
OH_val = 1e12       # typical daytime
CH3O2_val = 1e14    # typical
H2O_val = 4e23
M_val = 2.5e25
O2_val = 5.25e24
P_HOx_est = 1e12    # m⁻³/s (for HO2 estimation from Eq. 6.13)

# Vary NO from 10 ppt to 100 ppb
NO_ppb = 10 .^ range(-2, 2, length = 300)
NO_vals = NO_ppb .* 2.5e16  # m⁻³
NO2_vals = 2 .* NO_vals     # assume NO2/NO ratio ~ 2

# Estimate HO2 from steady state (Eqs. 6.13 and 6.18)
HO2_high_NOx = k_CO_OH_val .* CO_val .* OH_val ./ (k_HO2_NO_val .* NO_vals)
HO2_low_NOx = sqrt(P_HOx_est / (2 * k_HO2_HO2_val))
HO2_vals = min.(HO2_high_NOx, HO2_low_NOx)

# Solve the combined system for each NO level
prob = NonlinearProblem(compiled,
    Dict(compiled.O3 => O3_val, compiled.NO => NO_vals[1], compiled.NO2 => NO2_vals[1],
        compiled.OH => OH_val, compiled.HO2 => HO2_vals[1], compiled.CO => CO_val,
        compiled.CH3O2 => CH3O2_val, compiled.H2O => H2O_val, compiled.M => M_val,
        compiled.O2 => O2_val);
    build_initializeprob = false)

P_O3_net_vals = Float64[]
OPE_result = Float64[]
for i in eachindex(NO_ppb)
    newprob = remake(prob,
        p = [compiled.NO => NO_vals[i], compiled.NO2 => NO2_vals[i],
            compiled.HO2 => HO2_vals[i]])
    sol = solve(newprob)
    push!(P_O3_net_vals, sol[compiled.P_O3_net])
    push!(OPE_result, sol[compiled.OPE])
end

p1 = plot(NO_ppb, P_O3_net_vals ./ 1e12,
    xlabel = "NO (ppb)",
    ylabel = "Net P(O₃) (10¹² m⁻³ s⁻¹)",
    title = "Net O₃ Production vs NOx",
    xscale = :log10, linewidth = 2,
    label = "P(O₃)_net", legend = :topleft)
vline!([0.1], linestyle = :dash, color = :gray, label = "Background NO")
vline!([10.0], linestyle = :dash, color = :red, label = "Urban NO")

p2 = plot(NO_ppb, OPE_result,
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
environments by compiling and solving the `TroposphericChemistrySystem` with
the conditions from each environment.

```@example combined
environments = [
    ("Background", get_conditions_dict(TypicalConditions())),
    ("Urban", get_conditions_dict(UrbanConditions())),
    ("Remote", get_conditions_dict(RemoteConditions()))
]

results = []
for (name, cond) in environments
    env_prob = remake(prob,
        p = [compiled.O3 => cond[:O3], compiled.NO => cond[:NO],
            compiled.NO2 => cond[:NO2], compiled.OH => cond[:OH],
            compiled.HO2 => cond[:HO2], compiled.CO => cond[:CO],
            compiled.CH3O2 => cond[:CH3O2], compiled.H2O => cond[:H2O],
            compiled.M => cond[:M], compiled.O2 => cond[:O2]])
    sol = solve(env_prob)

    push!(results,
        (
            Environment = name,
            NO_ppb = round(cond[:NO] / 2.5e16, sigdigits = 3),
            O3_ppb = round(cond[:O3] / 2.5e16, sigdigits = 3),
            P_O3 = round(sol[compiled.P_O3_total], sigdigits = 3),
            OPE = round(sol[compiled.OPE], sigdigits = 3),
            Chain_Length = round(sol[compiled.chain_length], sigdigits = 3)
        ))
end

DataFrame(results)
```

This table shows how the key photochemical diagnostics vary across
different atmospheric environments, illustrating the transition from
NOx-limited (remote, high OPE and long chain length) to VOC-limited
(urban, low OPE and short chain length) regimes discussed in Chapter 6.

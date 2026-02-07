# SuperFast Gas-Phase Atmospheric Chemical Mechanism

## Overview

The Super Fast Chemical Mechanism is one of the simplest representations of atmospheric chemistry. It can efficiently simulate background tropospheric ozone chemistry and perform well for those species included in the mechanism. The mechanism includes 15 tracers and 30 reactions covering methane oxidation, oxidant chemistry, sulfur chemistry, and simplified isoprene chemistry.

**Reference**: Brown-Steiner, B., Selin, N. E., Prinn, R. G., Tilmes, S., Emmons, L., Lamarque, J.-F., and Cameron-Smith, P.: Evaluating simplified chemical mechanisms within present-day simulations of the Community Earth System Model version 1.2 with CAM4 (CESM1.2 CAM-chem): MOZART-4 vs. Reduced Hydrocarbon vs. Super-Fast chemistry, Geosci. Model Dev., 11, 4155-4174, [https://doi.org/10.5194/gmd-11-4155-2018](https://gmd.copernicus.org/articles/11/4155/2018/), 2018.

```@docs
SuperFast
```

## Implementation

The chemical equations are included in the supporting Table S2 of the paper. Rate constants use Arrhenius, Troe, and custom forms implemented as sub-systems that compute temperature- and pressure-dependent rate coefficients.

```@example superfast
using GasChem, EarthSciMLBase, ModelingToolkit, Symbolics
using DynamicQuantities, OrdinaryDiffEqRosenbrock
using Catalyst
using Plots
using DataFrames
using ModelingToolkit: t

model = SuperFast()
```

### State Variables

```@example superfast
vars = unknowns(model)[1:12]
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars],
    :Default => [ModelingToolkit.getdefault(v) for v in vars])
```

### Parameters

```@example superfast
params = parameters(model)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in params],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in params],
    :Description => [ModelingToolkit.getdescription(v) for v in params],
    :Default => [ModelingToolkit.getdefault(v) for v in params])
```

### Equations

```@example superfast
eqs = equations(model)
```

## Analysis

### Base Case Simulation

We can run simulations with the model, optionally changing the initial conditions and parameters. Here we change the initial concentration of O₃ to 15 ppb and the temperature to 293K:

```@example superfast
sys = mtkcompile(model)

tspan = (0.0, 3600*24)
prob = ODEProblem(sys, [sys.O3 => 15], tspan, [sys.T => 293])
sol = solve(prob, Rosenbrock23(), saveat = 10.0)

plot(sol, ylim = (0, 50), xlabel = "Time (s)",
    ylabel = "Concentration (ppb)", legend = :outerright)
```

### Temperature Sensitivity

The following simulations show how ozone concentrations respond to different temperatures, consistent with the temperature dependence of the Arrhenius rate expressions.

```@example superfast
sol1 = solve(ODEProblem(sys, [], tspan, [sys.T => 273]),
    Rosenbrock23(), saveat = 10.0)
sol2 = solve(ODEProblem(sys, [], tspan, [sys.T => 298]),
    Rosenbrock23(), saveat = 10.0)

plot([sol1[sys.O3], sol2[sys.O3]], label = ["T=273K" "T=298K"],
    title = "Change of O3 concentration at different temperatures",
    xlabel = "Time (s)", ylabel = "Concentration (ppb)")
```

### NO₂ Sensitivity

NO₂ is a key precursor for ozone production through photolysis (NO₂ + hv → NO + O₃). Higher NO₂ concentrations lead to increased ozone formation.

```@example superfast
sol_base = solve(ODEProblem(sys, [], tspan),
    Rosenbrock23(), saveat = 10.0)
sol_highNO2 = solve(ODEProblem(sys, [sys.NO2 => 100.0], tspan),
    Rosenbrock23(), saveat = 10.0)

p = plot(sol_base.t, sol_base[sys.O3], label = "Base case (NO₂=0.4 ppt)",
    xlabel = "Time (s)", ylabel = "O₃ concentration (ppb)")
plot!(p, sol_highNO2.t, sol_highNO2[sys.O3], label = "High NO₂ (100 ppb)")
title!(p, "Ozone response to NO₂ perturbation")
p
```

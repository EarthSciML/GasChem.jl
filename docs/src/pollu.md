# POLLU Atmospheric Chemistry Benchmark

## Overview

The POLLU problem is a 20-species, 25-reaction atmospheric chemistry benchmark problem
originally developed at the Dutch National Institute of Public Health and Environmental
Protection (RIVM). It is widely used for testing stiff ODE solvers due to its characteristic
atmospheric chemistry timescale separation.

**Reference**: J. G. Verwer, "Gauss-Seidel iteration for stiff ODEs from chemical kinetics,"
SIAM J. Sci. Comput., 15(5):1243-1259, 1994.

```@docs
Pollu
```

## Implementation

The original problem uses concentrations in ppm and time in minutes. The implementation
converts all rate constants to SI-compatible units (ppb for concentration, seconds for time):

  - First-order rates: divide by 60 (min⁻¹ to s⁻¹)
  - Second-order rates: divide by 60000 (ppm⁻¹ min⁻¹ to ppb⁻¹ s⁻¹)

### State Variables

```@example pollu
using GasChem, ModelingToolkit, DataFrames, Symbolics, DynamicQuantities

model = Pollu()
vars = unknowns(model)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars],
    :Default => [ModelingToolkit.getdefault(v) for v in vars]
)
```

### Parameters

```@example pollu
params = parameters(model)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in params],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in params],
    :Description => [ModelingToolkit.getdescription(v) for v in params],
    :Default => [ModelingToolkit.getdefault(v) for v in params]
)
```

### Equations

```@example pollu
eqs = equations(model)
```

## Analysis

### Reference Solution Comparison

The Verwer (1994) benchmark includes a high-precision reference solution at t = 60 minutes,
computed using RADAU5 with true double precision (uround = 1.01e-19). We reproduce this here
and compare against the reference values.

```@example pollu
using OrdinaryDiffEqRosenbrock, Plots

sys = mtkcompile(model)

# Initial conditions from the reference (converted ppm to ppb)
ic = [
    sys.NO2 => 0.0,
    sys.NO => 200.0,
    sys.O3P => 0.0,
    sys.O3 => 40.0,
    sys.HO2 => 0.0,
    sys.OH => 0.0,
    sys.CH2O => 100.0,
    sys.CO => 300.0,
    sys.ALD => 10.0,
    sys.MEO2 => 0.0,
    sys.C2O3 => 0.0,
    sys.CO2 => 0.0,
    sys.PAN => 0.0,
    sys.CH3O => 0.0,
    sys.HNO3 => 0.0,
    sys.O1D => 0.0,
    sys.SO2 => 7.0,
    sys.SO4 => 0.0,
    sys.NO3 => 0.0,
    sys.N2O5 => 0.0
]

tspan = (0.0, 3600.0)  # 60 minutes in seconds

sol = solve(
    ODEProblem(sys, ic, tspan),
    Rosenbrock23(),
    saveat = 1.0,
    abstol = 1e-12,
    reltol = 1e-12
)

p1 = plot(
    sol.t, sol[sys.NO2], label = "NO₂", ylabel = "Concentration (ppb)", xlabel = "Time (s)")
plot!(p1, sol.t, sol[sys.NO], label = "NO")
plot!(p1, sol.t, sol[sys.O3], label = "O₃")
plot!(p1, sol.t, sol[sys.HNO3], label = "HNO₃")
title!(p1, "Major Nitrogen and Ozone Species")
p1
```

### Reference Solution Table

The following table compares our computed solution at t = 60 minutes against the reference values from Verwer (1994).

```@example pollu
ref_species = [sys.NO2, sys.NO, sys.O3P, sys.O3, sys.HO2, sys.OH, sys.CH2O,
    sys.CO, sys.ALD, sys.MEO2, sys.C2O3, sys.CO2, sys.PAN, sys.CH3O,
    sys.HNO3, sys.O1D, sys.SO2, sys.SO4, sys.NO3, sys.N2O5]

# Reference values at t=60min from Verwer (1994), converted ppm to ppb
ref_vals = [0.5646255480022769e-01, 0.1342484130422339e+00, 0.4139734331099427e-08,
    0.5523140207484359e-02, 0.2018977262302196e-06, 0.1464541863493966e-06,
    0.7784249118997964e-01, 0.3245075353396018e+00, 0.7494013383880406e-02,
    0.1622293157301561e-07, 0.1135863833257075e-07, 0.2230505975721359e-02,
    0.2087162882798630e-03, 0.1396921016840158e-04, 0.8964884856898295e-02,
    0.4352846369330103e-17, 0.6899219696263405e-02, 0.1007803037365946e-03,
    0.1772146513969984e-05, 0.5682943292316392e-04] .* 1e3  # ppm to ppb

computed = [sol[s][end] for s in ref_species]
rel_err = abs.(computed .- ref_vals) ./ max.(abs.(ref_vals), 1e-20)

DataFrame(
    :Species => [string(Symbolics.tosymbol(v, escape=false)) for v in ref_species],
    Symbol("Reference (ppb)") => ref_vals,
    Symbol("Computed (ppb)") => computed,
    Symbol("Relative Error") => rel_err
)
```

### Formaldehyde and CO Evolution

```@example pollu
p2 = plot(sol.t, sol[sys.CH2O], label = "CH₂O",
    ylabel = "Concentration (ppb)", xlabel = "Time (s)")
plot!(p2, sol.t, sol[sys.CO], label = "CO")
plot!(p2, sol.t, sol[sys.ALD], label = "ALD")
title!(p2, "Organic Species")
p2
```

### Sulfur Conservation

The sulfur budget (SO₂ + SO₄) should be conserved since the only sulfur reaction is SO₂ + OH → SO₄ + HO₂.

```@example pollu
sulfur_total = sol[sys.SO2] .+ sol[sys.SO4]
p3 = plot(
    sol.t, sol[sys.SO2], label = "SO₂", ylabel = "Concentration (ppb)", xlabel = "Time (s)")
plot!(p3, sol.t, sol[sys.SO4], label = "SO₄")
plot!(p3, sol.t, sulfur_total, label = "Total S", linestyle = :dash)
title!(p3, "Sulfur Budget Conservation")
p3
```

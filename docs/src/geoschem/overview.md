# GEOS-Chem "fullchem" Gas-Phase Mechanism

This is an implementation of the GEOS-Chem "fullchem" gas-phase mechanism in Julia. 
The original version is used in the [GEOS-Chem](https://geoschem.github.io/) global 3-D atmospheric chemistry transport model.
The version here is adapted from GEOS-Chem version 14.1.1.
This mechanism is the result of many journal articles which are cited in API documentation for the [`GEOSChemGasPhase`](@ref) type.

!!! warning 
    This implementation is a work in progress.
    In particular, it does not yet include heterogeneous chemistry.

## System overview

First, let's initialize the model and we can also look at the first few ODE equations of the reaction network:

```@example 1
using GasChem, EarthSciMLBase
using DifferentialEquations, ModelingToolkit
using DynamicQuantities, Plots
using ModelingToolkit:t

tspan = (0.0, 360.0)
gc = GEOSChemGasPhase()
equations(gc)[1:5]
```

You can see that each reaction has a rate constant; rate constants are specified at the end of the list of equations:

```@example 1
equations(gc)[end-3:end]
```

## Simulation

Now, let's run a simulation and plot the results:

```@example 1
sys = structural_simplify(gc)
vals = ModelingToolkit.get_defaults(sys)
for k in setdiff(unknowns(sys),keys(vals))
    vals[k] = 0 # Set variables with no default to zero.
end
prob = ODEProblem(sys, vals, tspan, vals)
sol = solve(prob, AutoTsit5(Rosenbrock23()))
plot(sol, legend = :outertopright, xlabel = "Time (s)", 
        ylabel = "Concentration (ppb)")
```
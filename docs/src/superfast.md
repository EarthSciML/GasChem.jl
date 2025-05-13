# SuperFast Gas-Phase Atmospheric Chemical Mechanism

The Super Fast Chemical Mechanism is one of the simplest representations of atmospheric chemistry. It can efficiently simulate background tropospheric ozone chemistry and perform well for those species included in the mechanism. The chemical equations used are included in the supporting table S2 of the paper,
["Evaluating simplified chemical mechanisms within present-day simulations of the Community Earth System Model version 1.2 with CAM4 (CESM1.2 CAM-chem):
MOZART-4 vs. Reduced Hydrocarbon vs. Super-Fast chemistry" (2018), Benjamin Brown-Steiner, Noelle E. Selin, Ronald G. Prinn, Simone Tilmes, Louisa Emmons, Jean-François Lamarque, and Philip Cameron-Smith.](https://gmd.copernicus.org/articles/11/4155/2018/)

## Illustrative Example

Here is a simple example of initializing the SuperFast model and running a simulation.
First, we can look at the reaction equations:

```@example 1
using GasChem, EarthSciMLBase, ModelingToolkit
using DynamicQuantities, DifferentialEquations
using Catalyst
using Plots
using ModelingToolkit: t

model = SuperFast()
```

We can also look at them as a graph:

```@example 1
Graph(SuperFast(; name = :SuperFast, rxn_sys = true))
```

## Variables and parameters

The chemical species included in the superfast model are:

```@example 1
vars = unknowns(model)[1:12]
using DataFrames
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [ModelingToolkit.get_unit(v) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars],
    :Default => [ModelingToolkit.getdefault(v) for v in vars])
```

And here are the parameters:

```@example 1
vars = parameters(model)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [ModelingToolkit.get_unit(v) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars],
    :Default => [ModelingToolkit.getdefault(v) for v in vars])
```

## Running simulations

We can run simulations with the model, optionally changing the initial conditions and parameters. For example, we can change the initial concentration of O₃ to 15 ppb and the temperature to 293K:

```@example 1
sys = structural_simplify(model)

tspan = (0.0, 3600*24)
# Change the initial concentration of O₃ to 15 ppb and the temperature to 293K.
prob = ODEProblem(sys, [sys.O3 => 15], tspan, [sys.T => 293])
```

Now we can solve the system and plot the result:

```@example 1
sol = solve(prob, AutoTsit5(Rosenbrock23()), saveat = 10.0)

plot(sol, ylim = (0, 50), xlabel = "Time",
    ylabel = "Concentration (ppb)", legend = :outerright)
```

Finally let's run some simulations with different values for parameter `T`.

```@example 1
sol1 = solve(ODEProblem(sys, [], tspan, [sys.T => 273]),
    AutoTsit5(Rosenbrock23()), saveat = 10.0)
sol2 = solve(ODEProblem(sys, [], tspan, [sys.T => 298]),
    AutoTsit5(Rosenbrock23()), saveat = 10.0)

plot([sol1[sys.O3], sol2[sys.O3]], label = ["T=273K" "T=298K"],
    title = "Change of O3 concentration at different temperatures",
    xlabel = "Time (second)", ylabel = "concentration (ppb)")
```

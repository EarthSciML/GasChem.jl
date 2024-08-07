# SuperFast Gas-Phase Atmospheric Chemical Mechanism

The Super Fast Chemical Mechanism is one of the simplest representations of atmospheric chemistry. It can efficiently simulate background tropheric ozone chemistry and perform well for those species included in the mechanism. The chemical equations used is included in the supporting table S2 of the paper,
["Evaluating simplified chemical mechanisms within present-day simulations of the Community Earth System Model version 1.2 with CAM4 (CESM1.2 CAM-chem):
MOZART-4 vs. Reduced Hydrocarbon vs. Super-Fast chemistry" (2018), Benjamin Brown-Steiner, Noelle E. Selin, Ronald G. Prinn, Simone Tilmes, Louisa Emmons, Jean-François Lamarque, and Philip Cameron-Smith.](https://gmd.copernicus.org/articles/11/4155/2018/)


## Illustrative Example
Here is a simple example of initializing the SuperFast model and running a simulation.
First, we can look at the reaction equations:

```@example 1
using GasChem, EarthSciMLBase, ModelingToolkit, Unitful, DifferentialEquations
using Catalyst

@parameters t [unit = u"s", description = "Time"]
model = SuperFast(t)

model
```

We can also look at them as a graph:

```@example 1
Graph(model)
```


Before running any simulations with the model, we need to convert it into a system of differential equations. We can solve it using the default values for variables and parameters. However, by using the ```@unpack``` command, we can assign new values to specific variables and parameters, allowing for simulations under varied conditions.

We can visualize the differential equation version of the system as follows:
```@example 1
model = couple(model) # Convert the ReactionSystem in Catalyst.jl to CoupledSystem in EarthSciMLBase.jl
sys = structural_simplify(get_mtk(model))

vars = states(sys)  # Give you the variables in the system
var_dict = Dict(string(var) => var for var in vars)
pars = parameters(sys) # Give you the parameters in the system
par_dict = Dict(string(par) => par for par in pars)

tspan = (0.0, 3600*24)
u0 = [var_dict["superfast₊O3(t)"] => 15] # Change the initial concentration of O₃ to 15 ppb
p0 = [par_dict["superfast₊T"] => 293] # temperature = 293K
prob = ODEProblem(sys, u0, tspan, p0)

equations(sys)
```

We can finally solve the system and plot the result as

```@example 1
sol = solve(prob,AutoTsit5(Rosenbrock23()), saveat=10.0)

using Plots
plot(sol, ylim = (0,50), xlabel = "Time", ylabel = "Concentration (ppb)", legend=:outerright)
```

## Variables and parameters
The species included in the superfast model are: O₃, OH, HO₂, O₂, NO, NO₂, CH₄, CH₃O₂, H₂O, CH₂O, CO, CH₃OOH, CH₃O, DMS, SO₂, ISOP, O₁d, H₂O₂.

The parameters in the model that are not constant are the photolysis reaction rates ```jO31D```, ```j2OH```, ```jH2O2```, ```jNO2```, ```jCH2Oa```, ```jCH3OOH``` and temperature ```T```

Let's run some simulation with different values for parameter ```T```.
```@example 1
p1 = [par_dict["superfast₊T"] => 273]
p2 = [par_dict["superfast₊T"] => 298]
sol1 = solve(ODEProblem(sys, [], tspan, p1),AutoTsit5(Rosenbrock23()), saveat=10.0)
sol2 = solve(ODEProblem(sys, [], tspan, p2),AutoTsit5(Rosenbrock23()), saveat=10.0)

plot([sol1[var_dict["superfast₊O3(t)"]],sol2[var_dict["superfast₊O3(t)"]]], label = ["T=273K" "T=298K"], title = "Change of O3 concentration at different temperatures", xlabel="Time (second)", ylabel="concentration (ppb)")
```
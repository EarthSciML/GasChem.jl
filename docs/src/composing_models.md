```@meta
CurrentModule = GasChem
```

# Composing Models

## Illustrative Example
Here is the complete example of composing, visualizing and solving the SuperFast
model and the Fast-JX model, with explanation to follow:

```julia 
using EarthSciMLBase, GasChem, ModelingToolkit, OrdinaryDiffEq, Dates, Unitful, DifferentialEquations


@parameters t
composed_ode = SuperFast(t) + FastJX(t) # Compose two models simply use the "+" operator

start = Dates.datetime2unix(Dates.DateTime(2024, 2, 29))
tspan = (start, start+3600*24*3)
sys = structural_simplify(get_mtk(composed_ode)) # Define the coupled system  
sol = solve(ODEProblem(sys, [], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0) # Solve the coupled system
```

In the composed system, the variable name for O₃ is not ```O3``` but ```superfast₊O3(t)```. So we need some preparation of the result before visualizing. 

```julia
vars = states(sys)  # Get the variables in the composed system
var_dict = Dict(string(var) => var for var in vars)
pols = ["O3", "OH", "NO", "NO2", "CH4", "CH3O2", "CO","CH3OOH", "CH3O", "DMS", "SO2", "ISOP"]
var_names_p = ["superfast₊$(v)(t)" for v in pols]

x_t = unix2datetime.(sol[t]) # Convert from unixtime to date time for visualizing 
```
Then, we could plot the results as:
```@setup 1
using EarthSciMLBase, GasChem, ModelingToolkit, OrdinaryDiffEq, Dates, Unitful, DifferentialEquations


@parameters t
composed_ode = SuperFast(t) + FastJX(t) # Compose two models simply use the "+" operator

start = Dates.datetime2unix(Dates.DateTime(2024, 2, 29))
tspan = (start, start+3600*24*3)
sys = structural_simplify(get_mtk(composed_ode)) # Define the coupled system  
sol = solve(ODEProblem(sys, [], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0) # Solve the coupled system

vars = states(sys)  # Get the variables in the composed system
var_dict = Dict(string(var) => var for var in vars)
pols = ["O3", "OH", "NO", "NO2", "CH4", "CH3O2", "CO","CH3OOH", "CH3O", "DMS", "SO2", "ISOP"]
var_names_p = ["superfast₊$(v)(t)" for v in pols]

x_t = unix2datetime.(sol[t]) # Convert from unixtime to date time for visualizing 
```

```@example 1
using Plots
pp = []
for (i, v) in enumerate(var_names_p)
    name = pols[i]
    push!(pp, Plots.plot(x_t,sol[var_dict[v]],label = "$name", size = (1000, 600), xrotation=45))
end
Plots.plot(pp..., layout=(3, 4))
```

## Adding Emission Data
GasChem.jl incorporates an emissions model that utilizes data from the [US National Emissions Inventory for the year 2016](https://gaftp.epa.gov/Air/emismod/2016/v1/gridded/monthly_netCDF/). This model is activated as an extension when the ```EarthSciData``` package is used.
Here's a simple example:

```@example 2 
using GasChem, EarthSciData # This will trigger the emission extension
using Dates, ModelingToolkit, OrdinaryDiffEq, DifferentialEquations, EarthSciMLBase
using Plots

ModelingToolkit.check_units(eqs...) = nothing 
@parameters t
model_noemis = SuperFast(t)+FastJX(t) # A model with chemistry and photolysis, but no emissions.
model_withemis = SuperFast(t)+FastJX(t)+ 
    NEI2016MonthlyEmis{Float64}("mrggrid_withbeis_withrwc", t, -100.0, 30.0, 1.0, Δz) # The same model with emissions.

sys = structural_simplify(get_mtk(composed_ode))

start = Dates.datetime2unix(Dates.DateTime(2016, 5, 1))
tspan = (start, start+3*24*3600)
sol_noemis = solve(ODEProblem(structural_simplify(get_mtk(model_noemis)), [], tspan, []),
    AutoTsit5(Rosenbrock23()))
sol_withemis = solve(ODEProblem(structural_simplify(get_mtk(model_withemis)), [], tspan, []),
    AutoTsit5(Rosenbrock23()))

plot(
    plot(sol_noemis, title="Model without emissions"),
    plot!(sol_withemis, title="Model with emissions"),
)
```

Here is a plot that makes it easier to see what's going on for each species:
```@example 2
vars = states(sys)  # Get the variables names
var_dict = Dict(string(var) => var for var in vars)

x_t = unix2datetime.(sol[t])

pp = []
for (i, v) in enumerate(var_names_p)
    p = plot(x_t,sol_noemis[var_dict[v]], label="No Emissions", title = v, 
        size = (1000, 600), xrotation=45)
    plot!(x_t,sol_withemis[var_dict[v]], label="With Emissions", )
    push!(pp, p)
end
Plots.plot(pp..., layout=(3, 4))
```
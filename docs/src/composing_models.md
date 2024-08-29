```@meta
CurrentModule = GasChem
```

# Composing Models

## Illustrative Example
Here is the complete example of composing, visualizing and solving the SuperFast
model and the Fast-JX model, with explanation to follow:

```julia 
using EarthSciMLBase, GasChem, ModelingToolkit, Dates, DynamicQuantities, DifferentialEquations
using ModelingToolkit:t

composed_ode = couple(SuperFast(), FastJX()) # Compose two models use the "couple" function

start = Dates.datetime2unix(Dates.DateTime(2024, 2, 29))
tspan = (start, start+3600*24*3)
sys = structural_simplify(convert(ODESystem, composed_ode)) # Define the coupled system  
sol = solve(ODEProblem(sys, [], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0) # Solve the coupled system
```

In the composed system, the variable name for O₃ is not ```O3``` but ```superfast₊O3(t)```. So we need some preparation of the result before visualizing. 

```julia
vars = unknowns(sys)  # Get the variables in the composed system
var_dict = Dict(string(var) => var for var in vars)
pols = ["O3", "OH", "NO", "NO2", "CH4", "CH3O2", "CO","CH3OOH", "CH3O", "DMS", "SO2", "ISOP"]
var_names_p = ["superfast₊$(v)(t)" for v in pols]

x_t = unix2datetime.(sol[t]) # Convert from unixtime to date time for visualizing 
```
Then, we could plot the results as:
```@setup 1
using EarthSciMLBase, GasChem, ModelingToolkit, Dates, DynamicQuantities, DifferentialEquations
using ModelingToolkit:t

composed_ode = couple(SuperFast(), FastJX()) # Compose two models use the "couple" function

start = Dates.datetime2unix(Dates.DateTime(2024, 2, 29))
tspan = (start, start+3600*24*3)
sys = structural_simplify(convert(ODESystem, composed_ode)) # Define the coupled system  
sol = solve(ODEProblem(sys, [], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0) # Solve the coupled system

vars = unknowns(sys)  # Get the variables in the composed system
var_dict = Dict(string(var) => var for var in vars)
pols = ["O3", "OH", "NO", "NO2", "CH4", "CH3O2", "CO","CH3OOH", "CH3O", "DMS", "SO2", "ISOP"]
var_names_p = ["SuperFast₊$(v)(t)" for v in pols]

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
using Dates, ModelingToolkit, DifferentialEquations, EarthSciMLBase
using Plots, DynamicQuantities
using ModelingToolkit:t

@parameters lat = deg2rad(40.0f0) [unit=u"rad"]
@parameters lon = deg2rad(-97.0f0) [unit=u"rad"]
@parameters lev = 1
emis, emis_updater = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", lon, lat, lev; dtype=Float64)

model_noemis = couple(SuperFast(),FastJX()) # A model with chemistry and photolysis, but no emissions.
model_withemis = couple(SuperFast(), FastJX(), emis) # The same model with emissions.

sys_noemis = structural_simplify(convert(ODESystem, model_noemis))
sys_withemis = structural_simplify(convert(ODESystem, model_withemis))

start = Dates.datetime2unix(Dates.DateTime(2016, 5, 1))
tspan = (start, start+3*24*3600)
sol_noemis = solve(ODEProblem(sys_noemis, [], tspan, []), AutoTsit5(Rosenbrock23()))
sol_withemis = solve(ODEProblem(sys_withemis, [], tspan, []), AutoTsit5(Rosenbrock23()))

vars = unknowns(sys_noemis)  # Get the variables in the composed system
var_dict = Dict(string(var) => var for var in vars)
pols = ["O3", "OH", "NO", "NO2", "CH4", "CH3O2", "CO","CH3OOH", "CH3O", "DMS", "SO2", "ISOP"]
var_names_p = ["SuperFast₊$(v)(t)" for v in pols]

pp = []
for (i, v) in enumerate(var_names_p)
    p = plot(unix2datetime.(sol_noemis[t]),sol_noemis[var_dict[v]], label="No Emissions", title = v, 
        size = (1000, 600), xrotation=45)
    plot!(unix2datetime.(sol_withemis[t]),sol_withemis[var_dict[v]], label="With Emissions", )
    push!(pp, p)
end
Plots.plot(pp..., layout=(3, 4))
```
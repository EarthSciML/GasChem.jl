```@meta
CurrentModule = GasChem
```

# Composing Models

## Illustrative Example
Here is the complete example of composing, visualizing and solving the SuperFast
model and the Fast-JX model, with explanation to follow:

```@example 1
using EarthSciMLBase, GasChem, ModelingToolkit, Dates, DynamicQuantities, DifferentialEquations
using ModelingToolkit:t

composed_ode = couple(SuperFast(), FastJX()) # Compose two models use the "couple" function

start = Dates.datetime2unix(Dates.DateTime(2024, 2, 29))
tspan = (start, start+3600*24*3)
sys = convert(ODESystem, composed_ode) # Define the coupled system  
sol = solve(ODEProblem(sys, [], tspan, []), AutoTsit5(Rosenbrock23()), saveat=10.0) # Solve the coupled system
```

In the composed system, the variable name for O₃ is not ```O3``` but ```superfast₊O3(t)```. So we need some preparation of the result before visualizing. 

```@example 1
vars = unknowns(sys)  # Get the variables in the composed system
var_dict = Dict(string(var) => var for var in vars)
pols = ["O3", "OH", "NO", "NO2", "CH3O2", "CO","CH3OOH", "CH2O", "ISOP"]
var_names_p = ["SuperFast₊$(v)(t)" for v in pols]

x_t = unix2datetime.(sol[t]) # Convert from unixtime to date time for visualizing 
```
Then, we could plot the results as:

```@example 1
using Plots
pp = []
for (i, v) in enumerate(var_names_p)
    name = pols[i]
    push!(pp, Plots.plot(x_t,sol[var_dict[v]],label = "$name", size = (1000, 600), xrotation=45))
end
Plots.plot(pp..., layout=(3, 3))
```

## Adding Emission Data
GasChem.jl incorporates an emissions model that utilizes data from the [US National Emissions Inventory for the year 2016](https://gaftp.epa.gov/Air/emismod/2016/v1/gridded/monthly_netCDF/). This model is activated as an extension when the ```EarthSciData``` package is used.
Here's a simple example:

```@example 2 
using GasChem, EarthSciData # This will trigger the emission extension
using Dates, ModelingToolkit, DifferentialEquations, EarthSciMLBase
using Plots, DynamicQuantities
using ModelingToolkit:t

domain = DomainInfo(
    Dates.DateTime(2016, 5, 1), Dates.DateTime(2016, 5, 4);
    lonrange = deg2rad(-129):deg2rad(1):deg2rad(-61),
    latrange = deg2rad(11):deg2rad(1):deg2rad(59),
    levrange = 1:1:3,
    dtype = Float64)

emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain)

model_noemis = couple(SuperFast(),FastJX()) # A model with chemistry and photolysis, but no emissions.
model_withemis = couple(SuperFast(), FastJX(), emis) # The same model with emissions.

sys_noemis = convert(ODESystem, model_noemis)
sys_withemis = convert(ODESystem, model_withemis)

tspan = EarthSciMLBase.get_tspan(domain)
sol_noemis = solve(ODEProblem(sys_noemis, [], tspan, []), AutoTsit5(Rosenbrock23()))
sol_withemis = solve(ODEProblem(sys_withemis, [], tspan, []), AutoTsit5(Rosenbrock23()))

vars = unknowns(sys_noemis)  # Get the variables in the composed system
var_dict = Dict(string(var) => var for var in vars)
pols = ["O3", "OH", "NO", "NO2", "CH3O2", "CO","CH3OOH", "CH2O", "ISOP"]
var_names_p = ["SuperFast₊$(v)(t)" for v in pols]

pp = []
for (i, v) in enumerate(var_names_p)
    p = plot(unix2datetime.(sol_noemis[t]),sol_noemis[var_dict[v]], label="No Emissions", title = v, 
        size = (1000, 600), xrotation=45)
    plot!(unix2datetime.(sol_withemis[t]),sol_withemis[var_dict[v]], label="With Emissions", )
    push!(pp, p)
end
Plots.plot(pp..., layout=(3, 3))
```
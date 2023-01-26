```@meta
CurrentModule = GasChem
```
# Composing models

## Illustrative Example
Here is the complete example of composing, visualizing and solving the SuperFast
model and the Fast-JX model, with explanation to follow:

## Generating ODESystems
First, we need to build the two different components we want to compose together.

```julia
@parameters t 

sf = superfast(t) 
fj = fast_jx(t) 
```

Compose the SuperFast system and the Fast-JX system together.

```julia
connect = compose_fastjx_superfast(fj,sf)
```

This is now a differential-algebraic equation (DAE) of 18 variables. We can then define the resulting ODEProblem and send it over to DifferentialEquations.jl:
```julia
tspan = (0.0, 3600*24*2) #simulating for 2 days
sol = solve(ODEProblem(connect, [], tspan, [], combinatoric_ratelaws=false),Tsit5(), saveat=10.0)

using Plots
plot(sol,ylims=(0,20),xlabel="Time (second)", ylabel="concentration (ppb)",legend=:outertopright)
```
![Example1 Graph](compose_example.svg)
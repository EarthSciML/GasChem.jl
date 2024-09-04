# Model Parameters

## Temperature and Pressure Example

The two main parameters are temperature (`T` (K)) and pressure (`num_density` (molecules/cm³)).
We can explore what happens when we change them:

```@example
using GasChem, EarthSciMLBase
using DifferentialEquations, ModelingToolkit
using DynamicQuantities, Plots
using ModelingToolkit:t

tspan = (0.0, 60.0*60*24*4) # 4 day simulation

# Run a simulation with constant temperature and pressure.
sys = structural_simplify(GEOSChemGasPhase())
vals = ModelingToolkit.get_defaults(sys)
for k in setdiff(unknowns(sys),keys(vals))
    vals[k] = 0 # Set variables with no default to zero.
end
prob = ODEProblem(sys, vals, tspan, vals)
sol1 = solve(prob, AutoTsit5(Rosenbrock23()))

# Now, convert parameters to variables so we can change them over time.
sys2 = param_to_var(GEOSChemGasPhase(), :T, :num_density)

# Vary temperature and pressure over time.
@unpack T, num_density = sys2
@constants T_0 = 300 [unit=u"K"]
@constants t_0 = 1 [unit=u"s"]
eqs = [
    T ~ T_0 + T_0 / 1.5 * sin(2π*t/t_0/(60*60*24)),
    num_density ~ 2.7e19 - 2.5e19*t/t_0/(60*60*24*4),
]
sys2 = extend(sys2,ODESystem(eqs, t; name=:var_T))

# Run the simulation again.
sys2 = structural_simplify(sys2)
vals = ModelingToolkit.get_defaults(sys2)
for k in setdiff(unknowns(sys2),keys(vals))
    vals[k] = 0 # Set variables with no default to zero.
end
prob = ODEProblem(sys2, vals, tspan, vals)
sol2 = solve(prob, AutoTsit5(Rosenbrock23()))

# Plot the results
p1 = plot(sol1.t, sol1[sys2.O3], xticks=:none, label="Constant T and P",
        ylabel = "Concentration (ppb)")
plot!(p1, sol2.t, sol2[sys2.O3], label="Varying T and P")

p2 = plot(sol2.t, sol2[sys2.T], label = "T (K)", xticks=:none)
p3 = plot(sol2.t, sol2[sys2.num_density], label = "num_density (molec/cm³)", xlabel = "Time (s)")

plot(p1, p2, p3, layout=grid(3, 1, heights=[0.7, 0.15, 0.15]))
```

## Model Parameter Information

Here is a list of all of the model parameters:

```@example 1
using GasChem, DataFrames, EarthSciMLBase, ModelingToolkit, DynamicQuantities
@variables t [unit = u"s", description = "Time"]
gc = structural_simplify(GEOSChemGasPhase())
vars = parameters(gc)
DataFrame(
        :Name => [string(Symbolics.tosymbol(v, escape=false)) for v ∈ vars],
        :Units => [ModelingToolkit.get_unit(v) for v ∈ vars],
        :Description => [ModelingToolkit.getdescription(v) for v ∈ vars],
        :Default => [ModelingToolkit.getdefault(v) for v ∈ vars])
```
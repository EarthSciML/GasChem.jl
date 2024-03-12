# Model Parameters

The two main parameters are temperature (`T` (K)) and pressure (`num_density` (molecules/cm³)).
We can explore what happens when we change them:

!!! warning
    This demonstration does not current work.

```julia
using GasChem, EarthSciMLBase
using DifferentialEquations, ModelingToolkit
using Unitful, Plots

tspan = (0.0, 60.0*60*24*4) # 4 day simulation
@variables t [unit = u"s", description = "Time"]
sys = get_mtk(GEOSChemGasPhase(t))

# Convert parameters to variables so we can change them over time.
sys = param_to_var(sys, :T, :num_density)

# Vary temperature and pressure over time.
@constants T_0 = 300 [unit=u"K"]
@constants t_0 = 1 [unit=u"s"]
eqs = [
    sys.T ~ T_0 + T_0 / 150 * sin(2π*t/t_0/(60*60*24)),
    sys.num_density ~ 2.7e19 - 1e19*t/t_0/(60*60*24*4),
]
sys = extend(sys,ODESystem(eqs, t; name=:var_T))

# Run the simulation.
sys = structural_simplify(sys)
vals = ModelingToolkit.get_defaults(sys)
for k in setdiff(states(sys),keys(vals))
    vals[k] = 0 # Set variables with no default to zero.
end
prob = ODEProblem(sys, vals, tspan, vals)
sol = solve(prob, AutoTsit5(Rosenbrock23()))
p1 = plot(sol, legend = :outertopright, xticks=:none, 
        ylabel = "Concentration (nmol/mol)")

p2 = plot(sol.t, sol[sys.T], label = "T", xticks=:none)
p3 = plot(sol.t, sol[sys.numden], label = "numden", xlabel = "Time (s)")

plot(p1, p2, p3, layout=grid(3, 1, heights=[0.7, 0.15, 0.15]))
```

Here is a list of all of the model parameters:

```@example 1
using GasChem, DataFrames, EarthSciMLBase, ModelingToolkit, Unitful
@variables t [unit = u"s", description = "Time"]
gc = structural_simplify(get_mtk(GEOSChemGasPhase(t)))
vars = parameters(gc)
DataFrame(
        :Name => [string(Symbolics.tosymbol(v, escape=false)) for v ∈ vars],
        :Units => [ModelingToolkit.get_unit(v) for v ∈ vars],
        :Description => [ModelingToolkit.getdescription(v) for v ∈ vars],
        :Default => [ModelingToolkit.getdefault(v) for v ∈ vars])
```
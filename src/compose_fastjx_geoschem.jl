"""
Compose GEOSChemGasPhase and FastJX models together.

# Example
```
using GasChem, EarthSciMLBase, DifferentialEquations, Unitful

@parameters t [unit = u"s"]
gf = Couple(GEOSChemGasPhase(t), FastJX(t))
tspan = (0.0, 3600*24*2)
sys = structural_simplify(get_mtk(gf))
sol = solve(ODEProblem(sys, [], tspan, []),Rosenbrock23(), saveat=10.0)
using Plots
plot(sol,ylims=(0,20),xlabel="Time (second)", ylabel="concentration (ppb)",legend=:outertopright)
```
"""

@parameters t [unit = u"s"]

register_coupling(GEOSChemGasPhase(t), FastJX(t)) do c, p
    c = param_to_var(c, :j_3, :j_9, :j_11, :j_7, :j_10)
    @constants uconv = 1 [unit = u"s"]
    @constants c_fixme1 = 10^(-21) [unit = u"s"] # FIXME: Suspicious constant
    ConnectorSystem([
        c.j_9 ~ uconv * p.j_h2o2
        c.j_7 ~ uconv * p.j_CH2Oa
        c.j_10 ~ uconv * p.j_CH3OOH
        c.j_11 ~ uconv * p.j_NO2
        c.j_3 ~ c_fixme1 * p.j_o31D
    ],c , p)
end
"""
Compose GEOSChemGasPhase and FastJX models together.

# Example
```
using GasChem, EarthSciMLBase, DifferentialEquations

@parameters t 
gf = GEOSChemGasPhase(t) + FastJX(t)
tspan = (0.0, 3600*24*2)
sys = structural_simplify(get_mtk(gf))
sol = solve(ODEProblem(sys, [], tspan, []),Rosenbrock23(), saveat=10.0)
using Plots
plot(sol,ylims=(0,20),xlabel="Time (second)", ylabel="concentration (ppb)",legend=:outertopright)
```
"""
function Base.:(+)(g::GEOSChemGasPhase, f::FastJX)::ComposedEarthSciMLSystem
    sys = param_to_var(g.sys, :j_3, :j_9, :j_11, :j_7, :j_10)
    g = GEOSChemGasPhase(sys, g.rxn_sys)
    @constants uconv = 1 [unit = u"s"]
    @constants c_fixme1 = 10^(-21) [unit = u"s"] # FIXME: Suspicious constant
    ComposedEarthSciMLSystem(
        ConnectorSystem([
                g.sys.j_9 ~ uconv * f.sys.j_h2o2
                g.sys.j_7 ~ uconv * f.sys.j_CH2Oa
                g.sys.j_10 ~ uconv * f.sys.j_CH3OOH
                g.sys.j_11 ~ uconv * f.sys.j_NO2
                g.sys.j_3 ~ c_fixme1 * f.sys.j_o31D
            ], g, f),
        g, f,
    )
end
Base.:(+)(f::FastJX, s::GEOSChemGasPhase)::ComposedEarthSciMLSystem = s + f
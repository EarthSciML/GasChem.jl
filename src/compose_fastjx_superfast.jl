"""
Compose superfast and fast-jx models together.

# Example
```
using GasChem, EarthSciMLBase

@parameters t 
sf = SuperFast(t) + FastJX(t)
tspan = (0.0, 3600*24*2)
sol = solve(ODEProblem(structural_simplify(get_mtk(sf)), [], tspan, []),Tsit5(), saveat=10.0)
using Plots
plot(sol,ylims=(0,20),xlabel="Time (second)", ylabel="concentration (ppb)",legend=:outertopright)
```
"""
function Base.:(+)(s::SuperFast, f::FastJX)::ComposedEarthSciMLSystem
    sys = param_to_var(s.sys, :jO31D, :jH2O2, :jNO2, :jCH2Oa, :jCH3OOH)
    s = SuperFast(sys, s.rxn_sys)
    ComposedEarthSciMLSystem(
        ConnectorSystem([
            s.sys.jH2O2 ~ f.sys.j_h2o2
            s.sys.jCH2Oa ~ f.sys.j_CH2Oa
            s.sys.jCH3OOH ~ f.sys.j_CH3OOH
            s.sys.jNO2 ~ f.sys.j_NO2
            s.sys.jO31D ~ f.sys.j_o31D
        ], s, f),
        s, f,
    )
end
Base.:(+)(f::FastJX, s::SuperFast)::ComposedEarthSciMLSystem = s + f
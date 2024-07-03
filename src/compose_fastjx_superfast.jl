"""
Compose superfast and fast-jx models together.

# Example
```
using GasChem, EarthSciMLBase

@parameters t [unit = u"s"]
sf = couple(FastJX(t), SuperFast(t))
tspan = (0.0, 3600*24*2)
sol = solve(ODEProblem(structural_simplify(get_mtk(sf)), [], tspan, []),Tsit5(), saveat=10.0)
using Plots
plot(sol,ylims=(0,20),xlabel="Time (second)", ylabel="concentration (ppb)",legend=:outertopright)
```
"""

@parameters t [unit = u"s"]

register_coupling(SuperFast(t), FastJX(t)) do c, p
    c = param_to_var(convert(ODESystem, c), :jO31D, :jH2O2, :jNO2, :jCH2Oa, :jCH3OOH)
    ConnectorSystem([
        c.jH2O2 ~ p.j_h2o2
        c.jCH2Oa ~ p.j_CH2Oa
        c.jCH3OOH ~ p.j_CH3OOH
        c.jNO2 ~ p.j_NO2
        c.jO31D ~ p.j_o31D
    ], c, p)
end
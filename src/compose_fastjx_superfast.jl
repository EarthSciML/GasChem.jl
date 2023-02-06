"""
Converting the ReactionSystem to an ODESystem adds extra terms for the derivatives of the photolysis rate constants. We need to remove these terms before adding the fast-jx equations.
"""
function remove_D!(superfast)
    @parameters t
    @unpack jO31D, jH2O2, jNO2, jCH2Oa, jCH3OOH = superfast

    D = Differential(t)
    terms_to_remove = [D(jH2O2), D(jCH2Oa), D(jO31D), D(jCH3OOH), D(jNO2)]
    i_remove = []
    for (i, eq) in enumerate(equations(superfast))
        should_remove = false
        for term in terms_to_remove
            if isequal(term, eq.lhs)
                should_remove = true
            end
        end
        should_remove ? push!(i_remove, i) : continue
    end
    deleteat!(equations(superfast), i_remove)
    superfast
end

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
    remove_D!(s.sys)
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
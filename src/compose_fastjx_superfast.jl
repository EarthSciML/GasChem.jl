export remove_D, compose_fastjx_superfast
"""
Converting the ReactionSystem to an ODESystem adds extra terms for the derivatives of the photolysis rate constants. We need to remove these terms before adding the fast-jx equations.
"""
function remove_D(superfast)

    @unpack jO31D, jH2O2, jNO2, jCH2Oa, jCH3OOH = superfast

    sf = convert(ODESystem, superfast, combinatoric_ratelaws = false)
    D = Differential(t)
    terms_to_remove = [D(jH2O2), D(jCH2Oa), D(jO31D), D(jCH3OOH), D(jNO2)]
    sf_eqs = Equation[]
    for eq in equations(sf)
        should_remove = false
        for term in terms_to_remove
            if isequal(term, eq.lhs)
                should_remove = true
            end
        end
        should_remove ? continue : push!(sf_eqs, eq)
    end
    ODESystem(sf_eqs, t, states(sf), parameters(sf); name = nameof(superfast))
end

"""
Compose superfast and fast-jx models together.

# Example
```
using GasChem

@parameters t 
sf = superfast(t)
fj = fast_jx(t)
connect = compose_fastjx_superfast(fj,sf)
tspan = (0.0, 3600*24*2)
sol = solve(ODEProblem(connect, [], tspan, [], combinatoric_ratelaws=false),Tsit5(), saveat=10.0)
using Plots
plot(sol,ylims=(0,20),xlabel="Time (second)", ylabel="concentration (ppb)",legend=:outertopright)
```
"""
function compose_fastjx_superfast(fastjx, superfast)
    r_sf = remove_D(superfast)
    @named connected = ODESystem(
        [
            r_sf.jH2O2 ~ fastjx.j_h2o2
            r_sf.jCH2Oa ~ fastjx.j_CH2Oa
            r_sf.jCH3OOH ~ fastjx.j_CH3OOH
            r_sf.jNO2 ~ fastjx.j_NO2
            r_sf.jO31D ~ fastjx.j_o31D
        ],
        t,
        systems = [r_sf, fastjx],
    )
    simplified_sys = structural_simplify(connected)
    return simplified_sys
end

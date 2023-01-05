using ModelingToolkit
include("SuperFast.jl")
include("Fast-JX.jl")

#sf = superfast()
@named sf = ReactionSystem(superfast(), t, combinatoric_ratelaws=false)
odesys_sf = convert(ODESystem, sf)
fj = fast_jx()

@parameters t
connected = compose(ODESystem([
                        sf.jH2O2 ~ fj.j_h2o2
                        sf.jCH20a ~ fj.j_CH2Oa
                        sf.jCH20b ~ fj.j_CH2Ob
                        sf.jCH3OOH ~ fj.j_CH3OOH
                        sf.jNO2 ~ fj.j_NO2
                        sf.jO31D ~ fj.j_o31D
                      ], t; name=:connected), sf, fj)

states(connected)
equations(connected)

simplified_sys = structural_simplify(connected)

tspan = (0.0, 36000.0)
sol = solve(ODEProblem(simplified_sys, [], tspan, [], combinatoric_ratelaws=false, check_length=false),Tsit5(), saveat=10.0)

using Plots 
plot(sol)


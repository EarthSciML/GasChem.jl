using ModelingToolkit
using DifferentialEquations
include("SuperFast.jl")
include("Fast-JX.jl")

function remove_D(superfast)

  @unpack jO31D, jH2O2, jNO2, jCH2Oa, jCH3OOH = superfast

  sf = convert(ODESystem, superfast, combinatoric_ratelaws=false)

  # Converting the ReactionSystem to an ODESystem adds extra terms for the
  # derivatives of the photolysis rate constants. We need to remove these terms before
  # adding the fast-jx equations.
  D = Differential(t)
  terms_to_remove = [D(jH2O2), D(jCH2Oa), D(jO31D), D(jCH3OOH), D(jNO2)]
  sf_eqs = Equation[]
  for eq in equations(sf)
      should_remove = false
      for term in terms_to_remove
          @info term, eq.lhs, isequal(term, eq.lhs)
          if isequal(term, eq.lhs)
              should_remove = true
          end
      end
      should_remove ? continue : push!(sf_eqs, eq)
  end
  @named superfast = ODESystem(sf_eqs, t, states(sf), parameters(sf))
end

@parameters t 

sf = superfast(t)
r_sf = remove_D(sf)
fj = fast_jx(t)


@named connected = ODESystem([
  r_sf.jH2O2 ~ fj.j_h2o2
  r_sf.jCH2Oa ~ fj.j_CH2Oa
  r_sf.jCH3OOH ~ fj.j_CH3OOH
  r_sf.jNO2 ~ fj.j_NO2
  r_sf.jO31D ~ fj.j_o31D
                      ], t, systems=[r_sf, fj])

simplified_sys = structural_simplify(connected)

tspan = (0.0, 3600*24*2)
sol = solve(ODEProblem(simplified_sys, [], tspan, [], combinatoric_ratelaws=false),Tsit5(), saveat=10.0)
using Plots
plot(sol,ylims=(0,20),xlabel="Time (second)", ylabel="concentration (ppb)",legend=:outertopright)

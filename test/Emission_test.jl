using GasChem, EarthSciData
using Test, Dates, ModelingToolkit, OrdinaryDiffEq, DifferentialEquations, EarthSciMLBase

@testset "emission" begin
    @parameters t
    ModelingToolkit.check_units(eqs...) = nothing
    composed_ode = SuperFast(t) + Emission(t)
    sys = structural_simplify(get_mtk(composed_ode))
    @test length(states(sys)) ≈ 18
end

@testset "Nei2016extension" begin
    @parameters t
    ModelingToolkit.check_units(eqs...) = nothing
    start = Dates.datetime2unix(Dates.DateTime(2016, 5, 1))
    tspan = (start, start+3600*24*3)

    composed_ode = SuperFast(t) + Emission(t)
    sf = SuperFast(t)
    sys = structural_simplify(get_mtk(composed_ode))
    sys2 = structural_simplify(get_mtk(sf))
    
    sol = solve(ODEProblem(structural_simplify(get_mtk(composed_ode)), [], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)
    sol2 = solve(ODEProblem(structural_simplify(get_mtk(sf)), [], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)
    @test sol != sol2
end

    
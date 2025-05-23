using GasChem
using EarthSciMLBase
using DifferentialEquations, ModelingToolkit, DynamicQuantities
using ModelingToolkit: t
using Test

@testset "2wayCoupling" begin
    sol_middle = 10.054760758144594

    sf = couple(SuperFast(), FastJX())
    sys = convert(ODESystem, sf)
    tspan = (0.0, 3600 * 24)
    sol = solve(ODEProblem(sys, [], tspan, []), Rosenbrock23(), saveat = 10.0)
    @test sol[sys.SuperFast₊O3][4320] ≈ sol_middle rtol=1e-4
end

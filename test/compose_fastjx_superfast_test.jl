using GasChem
using EarthSciMLBase
using OrdinaryDiffEqRosenbrock, ModelingToolkit, DynamicQuantities
using ModelingToolkit: t
using Dates
using Test

@testset "2wayCoupling" begin
    sol_middle = 10.054760758144594

    sf = couple(SuperFast(), FastJX(0.0))
    sys = convert(ODESystem, sf)
    tspan = (0.0, 3600 * 24)
    prob = ODEProblem(sys, [], tspan, [])
    sol = solve(prob, Rosenbrock23(), saveat = 10.0)
    @test sol[sys.SuperFast₊O3][4320] ≈ sol_middle rtol=1e-4
end

@testset "2wayCoupling" begin
    sol_middle = 10.054760758144594

    sf = couple(SuperFast(), FastJX_interpolation_troposphere(0.0))
    sys = convert(ODESystem, sf)
    tspan = (0.0, 3600 * 24)
    prob = ODEProblem(sys, [], tspan, [])
    sol = solve(prob, Rosenbrock23(), saveat = 10.0)
    @test sol[sys.SuperFast₊O3][4320] ≈ sol_middle rtol=1e-4
end

@testitem "2wayCoupling" begin
    using EarthSciMLBase
    using OrdinaryDiffEqRosenbrock
    sol_middle = 10.054760758144594
    sf = couple(SuperFast(), FastJX(0.0))
    sys = convert(System, sf)
    tspan = (0.0, 3600 * 24)
    prob = ODEProblem(sys, [], tspan)
    sol = solve(prob, Rosenbrock23(), saveat = 10.0)
    @test sol[sys.SuperFast₊O3][4320]≈sol_middle rtol=1e-4
end

@testitem "2wayCoupling" begin
    using EarthSciMLBase
    using OrdinaryDiffEqRosenbrock
    sol_middle = 10.054760758144594

    sf = couple(SuperFast(), FastJX_interpolation_troposphere(0.0))
    sys = convert(System, sf)
    tspan = (0.0, 3600 * 24)
    prob = ODEProblem(sys, [], tspan)
    sol = solve(prob, Rosenbrock23(), saveat = 10.0)
    @test sol[sys.SuperFast₊O3][4320]≈sol_middle rtol=1e-4
end

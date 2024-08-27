using GasChem
using EarthSciMLBase
using DifferentialEquations, ModelingToolkit, DynamicQuantities
using ModelingToolkit:t

#   Unit Test
@testset "2wayCoupling" begin
    sol_middle = 9.948004877573444

    sf = couple(SuperFast(), FastJX())
    combined_mtk = convert(ODESystem, sf)
    sys = structural_simplify(combined_mtk)
    tspan = (0.0, 3600 * 24)
    sol = solve(
        ODEProblem(sys, [], tspan, []),
        Tsit5(),
        saveat = 10.0,
        abstol = 1e-8,
        reltol = 1e-8,
    )
    @test sol[sys.SuperFast.O3][4320] ≈ sol_middle
end

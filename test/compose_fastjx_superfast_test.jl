using GasChem
using EarthSciMLBase
using DifferentialEquations, ModelingToolkit, Unitful

#   Unit Test
@testset "2wayCoupling" begin
    sol_middle = 9.948004877573444
    @parameters t [unit = u"s"]

    sf = couple(SuperFast(t), FastJX(t))
    sys = structural_simplify(get_mtk(sf))
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

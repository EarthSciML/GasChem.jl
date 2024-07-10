using GasChem
using EarthSciMLBase
using DifferentialEquations, ModelingToolkit, Unitful

#   Unit Test
@testset "2wayCoupling" begin
    sol_middle = [
        9.948004797464273,
        -1.120813966829142e-6,
        4.517836132103682e-6,
        2.1000001328504878e8,
        0.0,
        10.0,
        1699.9493259138553,
        0.010708899259973205,
        450.1742122596299,
        0.1467289371492774,
        274.724815098401,
        1.7953005691445345,
        0.024279753277211612,
        48.41715760586734,
        3.5828423941319696,
        0.06476153381355561,
        6.519766848676772e-6,
        7.312607510153822,
    ]
    @parameters t [unit = u"s"]

    sf = couple(SuperFast(t), FastJX(t))
    tspan = (0.0, 3600 * 24)
    sol = solve(
        ODEProblem(structural_simplify(get_mtk(sf)), [], tspan, []),
        Tsit5(),
        saveat = 10.0,
    )
    @test sol[:, 4320] ≈ sol_middle
end

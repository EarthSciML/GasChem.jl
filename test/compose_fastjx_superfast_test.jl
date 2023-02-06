using GasChem
using EarthSciMLBase
using DifferentialEquations, ModelingToolkit

#   Unit Test
@test begin
    sol_middle = [
        13.391937312188814,
        1.66248182252514e-5,
        6.310062567496124e-5,
        2.10000233322078e8,
        3.1812457381754067,
        6.818754261824493,
        1699.9047452476004,
        1.570733237961087e-5,
        450.3894168691452,
        0.37418117772901327,
        274.5033352141585,
        1.634703057483833,
        0.0476341350477081,
        47.06609495457132,
        4.933905045429342,
        0.03313517657787631,
        0.5352891668840724,
        6.676368920250184,
    ]
    @parameters t
    sf = SuperFast(t) + FastJX(t)
    tspan = (0.0, 3600 * 24)
    sol = solve(
        ODEProblem(structural_simplify(get_mtk(sf)), [], tspan, []),
        Tsit5(),
        saveat = 10.0,
    )
    sol[4320] ≈ sol_middle
end

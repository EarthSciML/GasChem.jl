using GasChem
using EarthSciMLBase
using DifferentialEquations, ModelingToolkit

#   Unit Test
@test begin
    sol_middle = [
        19.077532687377207,
        0.00012384836570105673,
        0.00033229569773400683,
        2.1000037846421498e8,
        4.503643704452021,
        5.496356295548036,
        1699.6284053412858,
        4.093792287799892e-5,
        451.5375535458708,
        0.7970083235178884,
        273.5336586450974,
        0.8735935366030294,
        0.1514726826919632,
        39.4923169140811,
        12.507683085919219,
        0.0005514935961935434,
        0.543857042446943,
        3.3684450021603682,
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

using GasChem, EarthSciData
using Test, Dates, ModelingToolkit, OrdinaryDiffEq, DifferentialEquations, EarthSciMLBase, Unitful

@testset "Nei2016extension" begin
    @parameters t
    ModelingToolkit.check_units(eqs...) = nothing
    @parameters t [unit = u"s"]
    @parameters Δz = 60 [unit = u"m"]
    emis = NEI2016MonthlyEmis{Float64}("mrggrid_withbeis_withrwc", t, -100.0, 30.0, 1.0, Δz)

    sf = SuperFast(t)
    composed_ode = sf+emis
    sys = structural_simplify(get_mtk(composed_ode))
    @test length(states(sys)) ≈ 18
    
    start = Dates.datetime2unix(Dates.DateTime(2016, 5, 1))
    tspan = (start, start+3600*24*3)
    sol = solve(ODEProblem(structural_simplify(get_mtk(composed_ode)), [], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0)
    #@test isapprox(sol[end][end], 1.3654736679144752, atol = 0.001)
end
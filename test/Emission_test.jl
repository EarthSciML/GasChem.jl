using GasChem, EarthSciData
using Test, Dates, ModelingToolkit, DifferentialEquations, EarthSciMLBase, Unitful

@testset "NEI2016Extension3way" begin
    @parameters t [unit = u"s"]

    @parameters lat = 40 
    @parameters lon = -97 
    @parameters lev = 1
    @parameters Δz = 60 [unit = u"m"]
    emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", t, lon, lat, lev, Δz; dtype=Float64)

    model_3way = couple(FastJX(t), SuperFast(t), emis)

    sys = structural_simplify(get_mtk(model_3way))
    @test length(states(sys)) ≈ 18
    
    start = Dates.datetime2unix(Dates.DateTime(2016, 5, 1))
    tspan = (start, start+24*3600)
    sol_3way = solve(ODEProblem(sys, [], tspan, []), AutoTsit5(Rosenbrock23()))
    @test sol_3way[1,end] ≈ 2479.8538498918574
end
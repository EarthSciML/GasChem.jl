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

    eqs = string(equations(sys))
    wanteqs = ["Differential(t)(superfast₊CH2O(t)) ~ superfast₊NEI2016MonthlyEmis_mrggrid_withbeis_withrwc₊FORM(t)"]
    @test contains(string(eqs), wanteqs[1])
end
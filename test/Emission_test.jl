using GasChem, EarthSciData
using Test, Dates, ModelingToolkit, DifferentialEquations, EarthSciMLBase, Unitful

@testset "NEI2016Extension3way" begin
    @parameters t [unit = u"s"]

    @parameters lat = 40
    @parameters lon = -97
    @parameters lev = 1
    emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", t, lon, lat, lev; dtype=Float64)

    model_3way = couple(FastJX(t), SuperFast(t), emis)

    sys = structural_simplify(get_mtk(model_3way))
    @test length(states(sys)) ≈ 18

    eqs = string(equations(sys))
    wanteq = "Differential(t)(SuperFast₊CH2O(t)) ~ SuperFast₊NEI2016MonthlyEmis_FORM(t)"
    @test contains(string(eqs), wanteq)
end
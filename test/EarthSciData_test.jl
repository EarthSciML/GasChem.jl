using Main.GasChem, EarthSciData
using Test, Dates, ModelingToolkit, DifferentialEquations, EarthSciMLBase, DynamicQuantities

@testset "NEI2016Extension3way" begin

    @parameters lat = deg2rad(40.0f0) [unit=u"rad"]
    @parameters lon = deg2rad(-97.0f0) [unit=u"rad"]
    @parameters lev = 1
    emis, emis_updater = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", lon, lat, lev)

    model_3way = couple(FastJX(), SuperFast(), emis)

    sys = structural_simplify(convert(ODESystem, model_3way))
    @test length(unknowns(sys)) ≈ 18

    eqs = string(equations(sys))
    wanteq = "Differential(t)(SuperFast₊CH2O(t)) ~ SuperFast₊NEI2016MonthlyEmis_FORM(t)"
    @test contains(string(eqs), wanteq)
end


@testset "GEOS-FP" begin

    @parameters lat = deg2rad(40.0f0) [unit=u"rad"]
    @parameters lon = deg2rad(-97.0f0) [unit=u"rad"]
    @parameters lev = 1
    geosfp = GEOSFP("4x5")

    model_3way = couple(FastJX(), SuperFast(), geosfp)

    sys = structural_simplify(convert(ODESystem, model_3way))
    @test length(unknowns(sys)) ≈ 18

    eqs = string(equations(convert(ODESystem, model_3way)))
    wanteq = "SuperFast₊T(t) ~ GEOSFP₊I3₊T(t)"
    @test contains(eqs, wanteq)
    wanteq = "FastJX₊T(t) ~ GEOSFP₊I3₊T(t)"
    @test contains(eqs, wanteq)
    wanteq = "SuperFast₊jH2O2(t) ~ FastJX₊j_h2o2(t)"
    @test contains(eqs, wanteq)
end
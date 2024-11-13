using GasChem, EarthSciData
using Test, Dates, ModelingToolkit, DifferentialEquations, EarthSciMLBase, DynamicQuantities

# @testset "NEI2016Extension3way" begin
#     domain = DomainInfo(
#     DateTime(2016, 5, 1),
#     DateTime(2016, 5, 4);
#     lonrange = deg2rad(-115):deg2rad(2.5):deg2rad(-68.75),
#     latrange = deg2rad(25):deg2rad(2):deg2rad(53.7),
#     levrange = 1:15,
#     dtype = Float64)

#     emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain)

#     model_3way = couple(FastJX(), SuperFast(), emis)

#     sys = structural_simplify(convert(ODESystem, model_3way))
#     @test length(unknowns(sys)) ≈ 15

#     eqs = string(equations(sys))
#     wanteq = "Differential(t)(SuperFast₊CH2O(t)) ~ SuperFast₊NEI2016MonthlyEmis_FORM(t)"
#     @test contains(string(eqs), wanteq)
# end


@testset "GEOS-FP" begin

    domain = DomainInfo(
        DateTime(2016, 5, 1),
        DateTime(2016, 5, 4);
        lonrange = deg2rad(-115):deg2rad(2.5):deg2rad(-68.75),
        latrange = deg2rad(25):deg2rad(2):deg2rad(53.7),
        levrange = 1:15,
        dtype = Float64)

    geosfp = GEOSFP("4x5", domain)

    model_3way = couple(FastJX(), SuperFast(), geosfp)

    sys = structural_simplify(convert(ODESystem, model_3way))
    @test length(unknowns(sys)) ≈ 15

    eqs = string(equations(convert(ODESystem, model_3way)))
    wanteq = "SuperFast₊T(t) ~ GEOSFP₊I3₊T(t)"
    @test contains(eqs, wanteq)
    wanteq = "FastJX₊T(t) ~ GEOSFP₊I3₊T(t)"
    @test contains(eqs, wanteq)
    wanteq = "FastJX₊lat(t) ~ rad2deg(GEOSFP₊lat)"
    @test contains(eqs, wanteq)
    wanteq = "SuperFast₊jH2O2(t) ~ FastJX₊j_h2o2(t)"
    @test contains(eqs, wanteq)
end

using GasChem, EarthSciData
using Test, Dates, ModelingToolkit, EarthSciMLBase

@testset "NEI2016Extension3way" begin
    domain = DomainInfo(
    DateTime(2016, 5, 1),
    DateTime(2016, 5, 4);
    lonrange = deg2rad(-115):deg2rad(2.5):deg2rad(-68.75),
    latrange = deg2rad(25):deg2rad(2):deg2rad(53.7),
    levrange = 1:15,
    dtype = Float64)

    emis = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain)

    model_3way = couple(FastJX(), SuperFast(), emis)

    sys = convert(ODESystem, model_3way)
    @test length(unknowns(sys)) ≈ 12

    eqs = string(equations(sys))

    wanteq = "Differential(t)(SuperFast₊CH2O(t)) ~ SuperFast₊NEI2016MonthlyEmis_FORM(t)"
    @test contains(string(eqs), wanteq)
end


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

    sys = convert(ODESystem, model_3way)

    @test length(unknowns(sys)) ≈ 12

    eqs = string(observed(sys))
    wanteq = "SuperFast₊T(t) ~ GEOSFP₊I3₊T(t)"
    @test contains(eqs, wanteq) || contains(eqs, "SuperFast₊T(t) ~ FastJX₊T(t)")
    wanteq = "FastJX₊T(t) ~ GEOSFP₊I3₊T(t)"
    @test contains(eqs, wanteq)
    wanteq = "SuperFast₊jH2O2(t) ~ FastJX₊j_h2o2(t)"
    @test contains(eqs, wanteq)
    wanteq = "FastJX₊lat(t) ~ rad2deg(GEOSFP₊lat)"
    @test contains(eqs, wanteq)
    wanteq = "SuperFast₊P(t) ~ FastJX₊P(t)"
    @test contains(eqs, wanteq)
    wanteq = "FastJX₊P(t) ~ GEOSFP₊P(t)"
    @test contains(eqs, wanteq)
end

@testset "GEOSChemGasPhase couplings" begin
    domain = DomainInfo(DateTime(2016, 5, 1), DateTime(2016, 5, 2);
        latrange=deg2rad(-85.0f0):deg2rad(2):deg2rad(85.0f0),
        lonrange=deg2rad(-180.0f0):deg2rad(2.5):deg2rad(175.0f0),
        levrange=1:10, dtype=Float64)

    csys = couple(
        GEOSChemGasPhase(),
        NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain),
        GEOSFP("4x5", domain),
        FastJX(),
        domain,
    )
    sys = convert(ODESystem, csys)
    eqs = string(observed(sys))
    @test contains(eqs, "GEOSChemGasPhase₊T(t) ~ GEOSFP₊I3₊T(t)") ||
        contains(eqs, "GEOSChemGasPhase₊T(t) ~ FastJX₊T(t)")
    @test contains(eqs, "FastJX₊T(t) ~ GEOSFP₊I3₊T(t)") ||
        contains(eqs, "FastJX₊T(t) ~ GEOSChemGasPhase₊T(t)")
    @test contains(eqs, "GEOSChemGasPhase₊j_11(t) ~ FastJX₊j_NO2(t)")
end

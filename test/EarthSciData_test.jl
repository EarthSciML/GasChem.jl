@testitem "NEI2016Extension3way" begin
    using EarthSciData
    using Dates, ModelingToolkit, EarthSciMLBase

    domain = DomainInfo(
        DateTime(2016, 5, 1),
        DateTime(2016, 5, 4);
        lonrange = deg2rad(-115):deg2rad(2.5):deg2rad(-68.75),
        latrange = deg2rad(25):deg2rad(2):deg2rad(53.7),
        levrange = 1:15
    )

    model_3way = couple(
        FastJX(get_tref(domain)),
        SuperFast(),
        NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain)
    )

    sys = convert(System, model_3way)
    @test length(unknowns(sys)) ≈ 12

    eqs = string(equations(sys))
    # Use regex pattern that handles both ₊ and . namespace separators
    ns = "[₊.]"  # namespace separator pattern
    @test occursin(
        Regex("Differential\\(t\\)\\(SuperFast$(ns)CH2O\\(t\\)\\) ~ SuperFast$(ns)NEI2016MonthlyEmis_FORM\\(t\\)"),
        eqs)
end

@testitem "GEOS-FP" begin
    using EarthSciData
    using Dates, ModelingToolkit, EarthSciMLBase

    domain = DomainInfo(
        DateTime(2016, 5, 1),
        DateTime(2016, 5, 4);
        lonrange = deg2rad(-115):deg2rad(2.5):deg2rad(-68.75),
        latrange = deg2rad(25):deg2rad(2):deg2rad(53.7),
        levrange = 1:15
    )

    model_3way = couple(
        FastJX(domain),
        SuperFast(),
        GEOSFP("4x5", domain)
    )

    sys = convert(System, model_3way)

    @test length(unknowns(sys)) ≈ 12

    eqs = string(observed(sys))
    # Use regex patterns that handle both ₊ and . namespace separators
    ns = "[₊.]"  # namespace separator pattern
    # Test temperature coupling: SuperFast.T should be coupled to either GEOSFP.I3.T or FastJX.T
    @test occursin(Regex("SuperFast$(ns)T\\(t\\) ~ (GEOSFP$(ns)I3$(ns)T|FastJX$(ns)T)\\(t\\)"), eqs)
    # Test FastJX.T coupling to either GEOSFP.I3.T or SuperFast.T
    @test occursin(Regex("FastJX$(ns)T\\(t\\) ~ (GEOSFP$(ns)I3$(ns)T|SuperFast$(ns)T)\\(t\\)"), eqs)
    # Test jH2O2 coupling
    @test occursin(Regex("SuperFast$(ns)jH2O2\\(t\\) ~ FastJX$(ns)j_H2O2\\(t\\)"), eqs)
    # Test latitude coupling - note GEOSFP.lat may not have (t)
    @test occursin(Regex("FastJX$(ns)lat\\(t\\) ~ rad2deg\\(GEOSFP$(ns)lat"), eqs)
    # Test pressure couplings
    @test occursin(Regex("SuperFast$(ns)P\\(t\\) ~ FastJX$(ns)P\\(t\\)"), eqs)
    @test occursin(Regex("FastJX$(ns)P\\(t\\) ~ GEOSFP$(ns)P\\(t\\)"), eqs)
end

@testitem "GEOSChemGasPhase couplings" begin
    using EarthSciData
    using Dates, ModelingToolkit, EarthSciMLBase

    domain = DomainInfo(
        DateTime(2016, 5, 1),
        DateTime(2016, 5, 2);
        latrange = deg2rad(-85.0f0):deg2rad(2):deg2rad(85.0f0),
        lonrange = deg2rad(-180.0f0):deg2rad(2.5):deg2rad(175.0f0),
        levrange = 1:10
    )

    csys = couple(
        GEOSChemGasPhase(),
        NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain),
        GEOSFP("4x5", domain),
        FastJX(get_tref(domain)),
        domain
    )
    sys = convert(System, csys)
    eqs = string(observed(sys))
    # Use regex patterns that handle both ₊ and . namespace separators
    ns = "[₊.]"  # namespace separator pattern
    @test occursin(Regex("GEOSChemGasPhase$(ns)T\\(t\\) ~ (GEOSFP$(ns)I3$(ns)T|FastJX$(ns)T)\\(t\\)"), eqs)
    @test occursin(Regex("FastJX$(ns)T\\(t\\) ~ (GEOSFP$(ns)I3$(ns)T|GEOSChemGasPhase$(ns)T)\\(t\\)"), eqs)
    @test occursin(Regex("GEOSChemGasPhase$(ns)j_11\\(t\\) ~ FastJX$(ns)j_NO2\\(t\\)"), eqs)
end

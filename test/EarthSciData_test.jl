@testitem "NEI2016Extension3way" begin
    using GasChem, EarthSciData
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
    # Check that formaldehyde (CH2O) equation includes NEI2016 emissions term
    @test occursin(r"Differential.*SuperFast.*CH2O.*~.*SuperFast.*NEI2016MonthlyEmis_FORM"i, eqs) ||
        occursin(r"SuperFast.*NEI2016MonthlyEmis_FORM.*~.*Differential.*SuperFast.*CH2O"i, eqs)
end

@testitem "GEOS-FP" begin
    using GasChem, EarthSciData
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
    # Check that expected couplings exist using lenient patterns
    # The exact format depends on ModelingToolkit version, so check for key terms
    # Test temperature coupling: SuperFast T should be coupled somewhere
    @test occursin(r"SuperFast.*T.*~.*(GEOSFP|FastJX).*T"i, eqs) ||
        occursin(r"(GEOSFP|FastJX).*T.*~.*SuperFast.*T"i, eqs)
    # Test FastJX T coupling
    @test occursin(r"FastJX.*T.*~.*(GEOSFP|SuperFast).*T"i, eqs) ||
        occursin(r"(GEOSFP|SuperFast).*T.*~.*FastJX.*T"i, eqs)
    # Test jH2O2 coupling
    @test occursin(r"SuperFast.*jH2O2.*~.*FastJX.*j_H2O2"i, eqs) ||
        occursin(r"FastJX.*j_H2O2.*~.*SuperFast.*jH2O2"i, eqs)
    # Test latitude coupling
    @test occursin(r"FastJX.*lat.*~.*rad2deg.*GEOSFP.*lat"i, eqs) ||
        occursin(r"rad2deg.*GEOSFP.*lat.*~.*FastJX.*lat"i, eqs)
    # Test pressure couplings
    @test occursin(r"SuperFast.*P.*~.*FastJX.*P"i, eqs) ||
        occursin(r"FastJX.*P.*~.*SuperFast.*P"i, eqs)
    @test occursin(r"FastJX.*P.*~.*GEOSFP.*P"i, eqs) ||
        occursin(r"GEOSFP.*P.*~.*FastJX.*P"i, eqs)
end

@testitem "GEOSChemGasPhase couplings" begin
    using GasChem, EarthSciData
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
    # Check that expected couplings exist using lenient patterns
    @test occursin(r"GEOSChemGasPhase.*T.*~.*(GEOSFP|FastJX).*T"i, eqs) ||
        occursin(r"(GEOSFP|FastJX).*T.*~.*GEOSChemGasPhase.*T"i, eqs)
    @test occursin(r"FastJX.*T.*~.*(GEOSFP|GEOSChemGasPhase).*T"i, eqs) ||
        occursin(r"(GEOSFP|GEOSChemGasPhase).*T.*~.*FastJX.*T"i, eqs)
    @test occursin(r"GEOSChemGasPhase.*j_11.*~.*FastJX.*j_NO2"i, eqs) ||
        occursin(r"FastJX.*j_NO2.*~.*GEOSChemGasPhase.*j_11"i, eqs)
end

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
    @test length(unknowns(sys)) ≈ 19

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

    @test length(unknowns(sys)) ≈ 19

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

    # NEI ammonia emission feeds GEOS-Chem NH3 (anthropogenic + fertilizer sectors), the
    # source NH3 for the ISORROPIA aerosol partition (IsorropiaOp).
    @test occursin("NEI2016MonthlyEmis_NH3", string(equations(sys)))
end


@testitem "GEOS-Chem<->NEI: both NH3 sectors survive the translate table" begin
    using EarthSciMLBase, EarthSciData, Dates, ModelingToolkit
    using EarthSciMLBase: DomainInfo
    t0 = DateTime(2016, 3, 10)
    dom = DomainInfo(
        t0, t0 + Hour(3) + Hour(36);
        lonrange = deg2rad(-88):deg2rad(0.625):deg2rad(-86),
        latrange = deg2rad(32):deg2rad(0.5):deg2rad(33.5),
        levrange = 1:2, u_proto = zeros(Float64, 1, 1, 1, 1)
    )
    gfp = EarthSciData.GEOSFP("0.25x0.3125", dom)
    nei = EarthSciData.NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", dom)
    csys = convert(System, couple(dom, GEOSChemGasPhase(), nei, gfp))
    eqs = ModelingToolkit.equations(csys)
    isdiff(e, sp) = occursin("Differential", string(e.lhs)) &&
        occursin(Regex("GEOSChemGasPhase.$(sp)\\("), string(e.lhs))
    nh3eq = string(only(filter(e -> isdiff(e, "NH3"), eqs)))
    @test occursin("NH3_FERT", nh3eq) && occursin("₊NH3(t)", nh3eq)  # both NH3 sectors
end


@testitem "GEOS-Chem coupled invariants: met-driven H2O" begin
    using EarthSciMLBase, EarthSciData, Dates, ModelingToolkit
    using EarthSciMLBase: DomainInfo
    t0 = DateTime(2016, 3, 10)
    dom = DomainInfo(
        t0, t0 + Hour(3) + Hour(36);
        lonrange = deg2rad(-88):deg2rad(0.625):deg2rad(-86),
        latrange = deg2rad(32):deg2rad(0.5):deg2rad(33.5),
        levrange = 1:2, u_proto = zeros(Float64, 1, 1, 1, 1)
    )
    gfp = EarthSciData.GEOSFP("0.25x0.3125", dom)
    nei = EarthSciData.NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", dom)
    csys = convert(System, couple(dom, GEOSChemGasPhase(), nei, gfp))
    us = string.(unknowns(csys))
    # H2O is a met-driven constant species (parameter -> observed), not a state
    @test !any(u -> u == "H2O(t)" || endswith(u, "₊H2O(t)"), us)
    obs = string.([e.lhs for e in ModelingToolkit.observed(csys)])
    @test any(o -> endswith(o, "₊H2O(t)"), obs)
end

@testitem "GEOS-Chem<->NEI: SULF->SO4 and gas SO2 restored" begin
    using EarthSciMLBase, EarthSciData, Dates, ModelingToolkit
    using EarthSciMLBase: DomainInfo
    t0 = DateTime(2016, 3, 10)
    dom = DomainInfo(
        t0, t0 + Hour(3) + Hour(36);
        lonrange = deg2rad(-88):deg2rad(0.625):deg2rad(-86),
        latrange = deg2rad(32):deg2rad(0.5):deg2rad(33.5),
        levrange = 1:2, u_proto = zeros(Float64, 1, 1, 1, 1)
    )
    gfp = EarthSciData.GEOSFP("0.25x0.3125", dom)
    nei = EarthSciData.NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", dom)
    csys = convert(System, couple(dom, GEOSChemGasPhase(), nei, gfp))
    eqs = ModelingToolkit.equations(csys)
    isdiff(e, sp) = occursin("Differential", string(e.lhs)) &&
        occursin(Regex("GEOSChemGasPhase.$(sp)\\("), string(e.lhs))
    so2eq = string(only(filter(e -> isdiff(e, "SO2"), eqs)))
    so4eq = string(only(filter(e -> isdiff(e, "SO4"), eqs)))
    @test occursin("SO2", so2eq)                     # gas SO2 emission arrives
    @test occursin("SULF", so4eq)                    # direct sulfate -> SO4 ...
    @test !occursin("SULF", so2eq)                   # ... and NOT wrongly mapped into SO2
end

@testitem "GEOSChemGasPhase NEI2016 emission mappings" begin
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
        domain
    )
    sys = convert(System, csys)

    # Differential equation(s) for GEOS-Chem species `sp`.
    diffeqs(sp) = filter(
        e -> occursin("Differential", string(e.lhs)) &&
            occursin(Regex("GEOSChemGasPhase.$(sp)\\("), string(e.lhs)),
        equations(sys)
    )
    # Whether species `sp` receives the emission of NEI variable `nei`.
    emits(sp, nei) = any(
        e -> occursin("NEI2016MonthlyEmis_$(nei)", string(e)), diffeqs(sp)
    )

    # Every mapped species must survive as a state, otherwise `emits` below is vacuous.
    for sp in
        ["SO2", "SO4", "ALK4", "HNO4", "PRPE", "C2H6", "TOLU", "XYLE", "MTPA", "MEK"]
        @test length(diffeqs(sp)) == 1
    end

    # SULF is direct sulfate emission: it belongs on SO4, not gas-phase SO2, and the gas
    # SO2 emission must survive (a Dict translate table dropped it via duplicate keys).
    @test emits("SO2", "SO2")
    @test emits("SO4", "SULF")
    @test !emits("SO2", "SULF")

    # Mappings GEOS-Chem makes from this product that were previously missing.
    @test emits("ALK4", "PAR")
    @test emits("HNO4", "PNA")
    # Both OLE and IOLE map to PRPE, so PRPE appears twice in the translate table.
    @test emits("PRPE", "OLE")
    @test emits("PRPE", "IOLE")

    # A sample of the remaining newly mapped species.
    @test emits("C2H6", "ETHA")
    @test emits("TOLU", "TOL")
    @test emits("XYLE", "XYLMN")
    @test emits("MTPA", "TERP")
    @test emits("MEK", "KET")
end

@testitem "NEI2016 NO and HONO are NOx-as-NO2 mass" begin
    using GasChem, EarthSciData
    using Dates, ModelingToolkit, EarthSciMLBase  # ModelingToolkit re-exports Symbolics

    # Molar masses (kg/mol) as used by the coupling.
    MW_NO2 = 46.0055e-3
    MW_NO = 30.01e-3
    MW_HNO2 = 47.014e-3
    MW_CO = 28.0101e-3
    MW_air = 28.97e-3
    nmolpermol = 1.0e9

    t0, t1 = DateTime(2016, 5, 1), DateTime(2016, 5, 2)

    # 1. NO + NO2 + HONO == NOX in every emitting cell: the fields are NOx-as-NO2 mass.
    fs = EarthSciData.NEI2016MonthlyEmisFileSet("mrggrid_withbeis_withrwc", t0, t1)
    readvar(v) = Float64.(fs.ds[v][:, :, 1, :]) # (COL, ROW, month)
    NO, NO2, HONO, NOX = readvar("NO"), readvar("NO2"), readvar("HONO"), readvar("NOX")
    close(fs)

    emitting = NOX .> 0
    @test count(emitting) > 10_000
    @test all(isapprox.((NO .+ NO2 .+ HONO)[emitting], NOX[emitting]; rtol = 1.0e-5))
    @test all(iszero, (NO .+ NO2 .+ HONO)[.!emitting])
    # In their own molar masses the three would sum to ~1.48x NOX.
    @test sum(NO .* (MW_NO2 / MW_NO) .+ NO2 .+ HONO .* (MW_NO2 / MW_HNO2)) / sum(NOX) > 1.4

    # 2. The coupled NO and HONO emissions are the file mass divided by MW_NO2.
    domain = DomainInfo(
        t0,
        t1;
        lonrange = deg2rad(-115):deg2rad(2.5):deg2rad(-68.75),
        latrange = deg2rad(25):deg2rad(2):deg2rad(53.7),
        levrange = 1:15
    )
    nei = NEI2016MonthlyEmis("mrggrid_withbeis_withrwc", domain)

    # Coefficient of NEI2016MonthlyEmis₊neivar in model₊NEI2016MonthlyEmis_neivar (ppb/s per kg/kg/s).
    function emis_factor(sys, model, neivar)
        lhs = "$(model)₊NEI2016MonthlyEmis_$(neivar)(t)"
        eq = only(filter(e -> string(e.lhs) == lhs, observed(sys)))
        subs = Dict()
        for v in Symbolics.get_variables(eq.rhs)
            if occursin("NEI2016MonthlyEmis₊$(neivar)(t)", string(v))
                subs[v] = 1.0
            elseif ModelingToolkit.isconstant(v)
                subs[v] = ModelingToolkit.getdefault(v)
            end
        end
        return Float64(Symbolics.value(Symbolics.substitute(eq.rhs, subs)))
    end
    # Molar mass (kg/mol) by which the coupling divides the NEI mass field.
    divisor_mw(f) = nmolpermol * MW_air / f

    # SuperFast: NO divided by MW_NO2 like NO2; CO is the control.
    sf = convert(System, couple(SuperFast(), nei, domain))
    f_NO = emis_factor(sf, "SuperFast", "NO")
    @test divisor_mw(f_NO) ≈ 46.0055e-3 rtol = 1.0e-6
    @test f_NO ≈ emis_factor(sf, "SuperFast", "NO2") rtol = 1.0e-12
    @test !isapprox(divisor_mw(f_NO), MW_NO; rtol = 0.1)
    @test divisor_mw(emis_factor(sf, "SuperFast", "CO")) ≈ MW_CO rtol = 1.0e-6

    # GEOS-Chem fullchem: NO and HONO both divided by MW_NO2.
    gc = convert(System, couple(GEOSChemGasPhase(), nei, domain))
    f_NO = emis_factor(gc, "GEOSChemGasPhase", "NO")
    f_HONO = emis_factor(gc, "GEOSChemGasPhase", "HONO")
    @test divisor_mw(f_NO) ≈ 46.0055e-3 rtol = 1.0e-6
    @test divisor_mw(f_HONO) ≈ 46.0055e-3 rtol = 1.0e-6
    @test f_NO ≈ emis_factor(gc, "GEOSChemGasPhase", "NO2") rtol = 1.0e-12
    @test !isapprox(divisor_mw(f_NO), MW_NO; rtol = 0.1)
    @test !isapprox(divisor_mw(f_HONO), MW_HNO2; rtol = 0.01)
    @test divisor_mw(emis_factor(gc, "GEOSChemGasPhase", "CO")) ≈ MW_CO rtol = 1.0e-6
end

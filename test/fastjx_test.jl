@testsnippet FastJXSetup begin
    using GasChem, ModelingToolkit
    using SymbolicIndexingInterface: setp, getsym, parameter_values

    function get_fluxes(t, lat, lon, P)
        cos_sza = GasChem.cos_solar_zenith_angle(t, lat, lon)
        return GasChem.calc_direct_fluxes(cos_sza, P)
    end

    fj = mtkcompile(FastJX(0.0))
    # 18:00 UTC is local noon at 97W (cos(SZA) = 0.45).  Do not move this back to 12:00
    # UTC: the sun is then 20 degrees below the horizon there, every actinic flux is
    # exactly zero, and every assertion below that compares the compiled system against
    # `j_mean_*` passes by comparing 0.0 to 0.0.
    test_time = 3600 * 18.0
    # FastJX scales every band's flux by the Earth-Sun distance factor (SOLF); the bare
    # `j_mean_*` helpers below take a flux array and do not. Comparisons between the two
    # have to carry it explicitly.
    solf = GasChem.solar_flux_factor(test_time)
    # [t_ref, lat, long, T, P, H2O]
    p = [0.0, 40.0, -97.0, 298.0, 101325.0, 450.0]
    prob = ODEProblem(
        fj,
        [
            fj.t_ref => 0.0, fj.lat => 40.0, fj.long => -97.0,
            fj.T => 298.0, fj.P => 101325.0, fj.H2O => 450.0,
        ],
        (test_time, test_time + 1)
    )
end

#   Unit Test 0: O3 -> O2 + O(1D)

@testitem "o31D" setup = [FastJXSetup] begin
    u_0 = [
        0.007356510224173006,
        0.007361736924955544,
        0.007570811355827869,
        0.007570811355827869,
    ]
    fluxes = get_fluxes(3600 * 12.0, 30.0, 0.0, 0.9)
    test_0 = [
        GasChem.j_mean_O31D(100.0, fluxes),
        GasChem.j_mean_O31D(220.0, fluxes),
        GasChem.j_mean_O31D(300.0, fluxes),
        GasChem.j_mean_O31D(400.0, fluxes),
    ]
    @test test_0 ≈ u_0 rtol = 1.0e-6
end

#   Unit Test 1: H2O2 -> OH + OH
@testitem "H2O2" setup = [FastJXSetup] begin
    u_1 = [9.556109440917478e-5, 9.784078586085366e-5, 0.00010012047731253253]

    fluxes = get_fluxes(3600 * 12.0, 30.0, 0.0, 0.9)
    test_1 = [
        GasChem.j_mean_H2O2(150.0, fluxes),
        GasChem.j_mean_H2O2(250.0, fluxes),
        GasChem.j_mean_H2O2(350.0, fluxes),
    ]

    @test test_1 ≈ u_1

    j_H2O2_func = getsym(prob, fj.j_H2O2)
    j_H2O2_value = j_H2O2_func(prob)
    j_want = solf * GasChem.j_mean_H2O2(298.0, get_fluxes(test_time, 40.0, -97.0, 101325))
    @test j_H2O2_value ≈ j_want rtol = 0.004
end

# Unit Test 2: CH2O -> H + HO2 + CO
@testitem "H2COa" setup = [FastJXSetup] begin
    u_2 = [8.642676785125311e-5, 8.64227156795509e-5, 8.641551181874694e-5]

    fluxes = get_fluxes(3600 * 12.0, 30.0, 0.0, 0.9)
    test_2 = [
        GasChem.j_mean_H2COa(200.0, fluxes),
        GasChem.j_mean_H2COa(250.0, fluxes),
        GasChem.j_mean_H2COa(300.0, fluxes),
    ]

    @test test_2 ≈ u_2

    j_H2COa_func = getsym(prob, fj.j_H2COa)
    j_H2COa_value = j_H2COa_func(prob)
    @test j_H2COa_value ≈ solf * GasChem.j_mean_H2COa(
        298.0,
        get_fluxes(test_time, 40.0, -97.0, 101325)
    ) rtol = 1.0e-6
end

@testitem "H2COb" setup = [FastJXSetup] begin
    u_2 = [7.379813688974829e-5, 7.383806831956884e-5, 7.39090575281387e-5]

    fluxes = get_fluxes(3600 * 12.0, 30.0, 0.0, 0.9)
    test_2 = [
        GasChem.j_mean_H2COb(200.0, fluxes),
        GasChem.j_mean_H2COb(250.0, fluxes),
        GasChem.j_mean_H2COb(300.0, fluxes),
    ]

    @test test_2 ≈ u_2

    j_H2COb_func = getsym(prob, fj.j_H2COb)
    j_H2COb_value = j_H2COb_func(prob)
    @test j_H2COb_value ≈
        solf * GasChem.j_mean_H2COb(
        298.0,
        get_fluxes(test_time, 40.0, -97.0, 101325)
    ) rtol = 1.0e-6
end

# Unit Test 3: CH3OOH -> OH + HO2 + CH2O
@testitem "CH3OOH" setup = [FastJXSetup] begin
    u_3 = [5.479266685458071e-5, 5.479266685458071e-5, 5.479266685458071e-5]

    test_3 = [
        GasChem.j_mean_CH3OOH(200.0, get_fluxes(3600 * 6.0, 30.0, 0.0, 0.9)),
        GasChem.j_mean_CH3OOH(200.0, get_fluxes(3600 * 12.0, 30.0, 0.0, 0.9)),
        GasChem.j_mean_CH3OOH(200.0, get_fluxes(3600 * 18.0, 30.0, 0.0, 0.9)),
    ]

    @test test_3 ≈ u_3

    j_CH3OOH_func = getsym(prob, fj.j_CH3OOH)
    j_CH3OOH_value = j_CH3OOH_func(prob)
    @test j_CH3OOH_value ≈
        solf * GasChem.j_mean_CH3OOH(
        298.0,
        get_fluxes(test_time, 40.0, -97.0, 101325)
    ) rtol = 1.0e-6
end

# Unit Test 4: NO2 -> NO + O
@testitem "NO2" setup = [FastJXSetup] begin
    u_4 = [0.008441071700500404, 0.00881068096075025, 0.009135937109770115]

    fluxes = get_fluxes(3600 * 12.0, 30.0, 0.0, 0.9)
    test_4 = [
        GasChem.j_mean_NO2(150.0, fluxes),
        GasChem.j_mean_NO2(250.0, fluxes),
        GasChem.j_mean_NO2(300.0, fluxes),
    ]

    @test test_4 ≈ u_4 rtol = 1.0e-6

    j_NO2_func = getsym(prob, fj.j_NO2)
    j_NO2_value = j_NO2_func(prob)
    # guards the assertion below against silently degenerating to 0.0 == 0.0
    @test j_NO2_value > 1.0e-4

    @test j_NO2_value ≈
        solf * GasChem.j_mean_NO2(298.0, get_fluxes(test_time, 40.0, -97.0, 101325)) rtol = 1.0e-6
end

@testitem "GEOS-Chem: CFCl3, H1301, Glyxlc" setup = [FastJXSetup] begin
    # [t_ref, lat, long, T, P, H2O]
    p = [0.0, 40.0, -97.0, 298.0, 101325.0, 450.0]
    j_CFCl3_func = getsym(prob, fj.j_CFCl3)
    j_CFCl3_value = j_CFCl3_func(prob)

    @test j_CFCl3_value ≈
        solf * GasChem.j_mean_CFCl3(298.0, get_fluxes(test_time, 40.0, -97.0, 101325)) rtol = 1.0e-6

    j_H1301_func = getsym(prob, fj.j_H1301)
    j_H1301_value = j_H1301_func(prob)

    @test j_H1301_value ≈
        solf * GasChem.j_mean_H1301(298.0, get_fluxes(test_time, 40.0, -97.0, 101325)) rtol = 1.0e-6

    j_Glyxlc_func = getsym(prob, fj.j_Glyxlc)
    j_Glyxlc_value = j_Glyxlc_func(prob)

    @test j_Glyxlc_value ≈
        solf * GasChem.j_mean_Glyxlc(298.0, get_fluxes(test_time, 40.0, -97.0, 101325)) rtol = 1.0e-6
end

@testitem "Ensure Cos SZA is non-allocating" begin
    using GasChem, AllocCheck
    @check_allocs checkcos(lat, t, long) = GasChem.cos_solar_zenith_angle(t, lat, long)
    checkcos(0.0, 0.0, 0.0)
    checkcos(0.0f0, 0.0f0, 0.0f0)
end

@testitem "FastJX Initialization" begin
    using GasChem, ModelingToolkit
    @test_nowarn mtkcompile(GasChem.FastJX(0.0))
end

@testitem "FastJX_interpolation_troposphere Initialization" begin
    using GasChem, ModelingToolkit
    @test_nowarn mtkcompile(GasChem.FastJX_interpolation_troposphere(0.0))                    # mech=:all (default)
    @test_nowarn mtkcompile(GasChem.FastJX_interpolation_troposphere(0.0; mech = :superfast))
    @test_throws ArgumentError GasChem.FastJX_interpolation_troposphere(0.0; mech = :nope)
end

@testitem "Direct Flux" begin
    using GasChem
    @test GasChem.calc_direct_fluxes(0.42255961917649837, 1013525) ≈ [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.008445690580780573,
        4.368635862759404e6,
        3.861091113407282e11,
        7.162608365510285e12,
        4.381532157326822e13,
        5.5736656209660256e14,
        5.378272528516121e15,
        1.6978055609093792e17,
    ]
end

@testitem "Direct Flux 2" begin
    using GasChem
    t, lat, lon, P = 3600 * 12.0, 30.0, 0.0, 0.9

    cos_sza = GasChem.cos_solar_zenith_angle(t, lat, lon)

    @test GasChem.calc_direct_flux(cos_sza, P, 1) ≈ 1.391000027136e12
    @test GasChem.calc_direct_flux(cos_sza, P, 2) ≈ 1.6270000128e12
    @test GasChem.calc_direct_flux(cos_sza, P, 4) ≈ 9.27799967744e11
    @test GasChem.calc_direct_flux(cos_sza, P, 6) ≈ 4.680000208896e12
    @test GasChem.calc_direct_flux(cos_sza, P, 8) ≈ 1.219000008704e13
    @test GasChem.calc_direct_flux(cos_sza, P, 10) ≈ 4.0489998876672e14
    @test GasChem.calc_direct_flux(cos_sza, P, 12) ≈ 5.88900011606016e14
    @test GasChem.calc_direct_flux(cos_sza, P, 14) ≈ 5.04500011925504e14
    @test GasChem.calc_direct_flux(cos_sza, P, 16) ≈ 3.853000128856064e15
    @test GasChem.calc_direct_flux(cos_sza, P, 18) ≈ 2.1310000789140275e17

    P = 100
    @test GasChem.calc_direct_flux(cos_sza, P, 1) ≈ 3.0643508503689505e6
    P = 500
    @test GasChem.calc_direct_flux(cos_sza, P, 1) ≈ 8.931704683157367e-17
end

@testitem "Direct Flux Twilight" begin
    using GasChem
    P = 1
    cos_sza = -0.1
    @test GasChem.calc_direct_flux(cos_sza, P, 1) ≈ 0.0
    @test GasChem.calc_direct_flux(cos_sza, P, 3) ≈ 890275.9088383563
    @test GasChem.calc_direct_flux(cos_sza, P, 5) ≈ 4.250677812848234e10
    @test GasChem.calc_direct_flux(cos_sza, P, 7) ≈ 1.6784210577555624e10
    @test GasChem.calc_direct_flux(cos_sza, P, 9) ≈ 8.920108354738617e-7
    @test GasChem.calc_direct_flux(cos_sza, P, 11) ≈ 5.4169480803700356e10
    @test GasChem.calc_direct_flux(cos_sza, P, 13) ≈ 1.697289300309447e14
    @test GasChem.calc_direct_flux(cos_sza, P, 15) ≈ 6.569697111320194e14
    @test GasChem.calc_direct_flux(cos_sza, P, 18) ≈ 2.0972515867967904e17

    P = 100
    @test GasChem.calc_direct_flux(cos_sza, P, 3) ≈ 0.0
    @test GasChem.calc_direct_flux(cos_sza, P, 9) ≈ 0.0
    @test GasChem.calc_direct_flux(cos_sza, P, 18) ≈ 4.908683888514731e16
end

#   Dedicated cross-sections for the GEOS-Chem photolysis-completion channels, transcribed
#   from ExtData/CHEM_INPUTS/CLOUD_J/v2024-09/FJX_spec.dat.
@testitem "GEOSChem dedicated photolysis cross-sections" setup = [FastJXSetup] begin
    # Reference j = Σ flux·σ(T)·ϕ at the top-of-atmosphere flux used by the other unit tests.
    fluxes = get_fluxes(3600 * 12.0, 30.0, 0.0, 0.9)
    expected = Dict(
        :ONIT1 => 1.2879551602032073e-5, :ONIT2 => 2.5759103204064146e-5,
        :ETNO3 => 0.00035574512324203873, :IPRNO3 => 0.000426503946987196,
        :NPRNO3 => 0.0004035466117803549, :MVKN => 0.00016795574750408425,
        :MACRN => 0.0005696433077568555, :MACRNP => 0.00015608234645796674,
        :ICN => 0.0003326633053236463, :ETHLN => 0.0002499367255213189,
        :NITP => 6.766986252962223e-5, :HMHP => 3.835302279140898e-5,
        :HP2 => 0.00010958712986581257, :ENOL => 0.00033785045579737445,
        :PROPNN => 0.00011333939, :HPALD1 => 0.0003434080894,
        :HPALD2 => 0.0003256792373, :PrAldP => 0.00023912282840808,
        :BALD => 0.0006601164,
    )
    for (s, val) in expected
        jf = getfield(GasChem, Symbol("j_mean_", s))
        @test jf(298.0, fluxes) ≈ val rtol = 1.0e-3
    end
    # Two-temperature cross-sections (ETNO3, IPRNO3) must actually interpolate in T.
    @test GasChem.j_mean_ETNO3(240.0, fluxes) != GasChem.j_mean_ETNO3(298.0, fluxes)
    @test GasChem.j_mean_IPRNO3(240.0, fluxes) != GasChem.j_mean_IPRNO3(298.0, fluxes)
    # Diurnal: at surface pressure the slant-path airmass zeroes the direct beam at night.
    night = get_fluxes(3600 * 0.0, 30.0, 0.0, 101325.0)
    @test GasChem.j_mean_ONIT1(298.0, night) == 0.0
    @test GasChem.j_mean_HP2(298.0, night) == 0.0
end

#   Five of GEOS-Chem's cross-sections are exact linear combinations of others, and are
#   written that way in Fast-JX.jl rather than copied out.  These identities must hold.
@testitem "GEOSChem combination cross-sections" setup = [FastJXSetup] begin
    fluxes = get_fluxes(3600 * 12.0, 30.0, 0.0, 0.9)
    j(sp) = getfield(GasChem, Symbol("j_mean_", sp))(298.0, fluxes)
    @test j(:CH3OOH) > 0 && j(:ONIT1) > 0 && j(:MACRN) > 0
    @test j(:HP2) ≈ 2 * j(:CH3OOH) rtol = 1.0e-6
    @test j(:HMHP) ≈ 0.7 * j(:CH3OOH) rtol = 1.0e-6
    @test j(:ONIT2) ≈ 2 * j(:ONIT1) rtol = 1.0e-6
    @test j(:NITP) ≈ j(:CH3OOH) + j(:ONIT1) rtol = 1.0e-6
    # bin 17 of the tabulated MACRNP reads 1.734e-23 where the combination gives 1.743e-23
    @test j(:MACRNP) ≈ 0.25 * (j(:CH3OOH) + j(:MACRN)) rtol = 1.0e-4
end

#   Halogen / iodine / remaining inorganic cross-sections, same source file.
@testitem "GEOSChem halogen and inorganic cross-sections" setup = [FastJXSetup] begin
    fluxes = get_fluxes(3600 * 12.0, 30.0, 0.0, 0.9)
    expected = Dict(
        :ClNO2 => 0.001552311385, :Br2 => 0.030246335083853996,
        :BrNO2 => 0.0098252416, :HAC => 8.002663374e-5,
        :H2SO4 => 5.417002e-8, :ClOO => 0.481128504278,
        :MPN => 2.97824521e-5, :I2 => 0.178869914386,
        :HOI => 0.00854073857, :IO => 0.1903769419,
        :OIO => 0.20600377, :INO => 0.05162421918,
        :IONO => 0.005615942886, :IONO2 => 0.013487755090000001,
        :I2O2 => 0.051839342, :CH2I2 => 0.0185509959228,
        :CH2ICl => 0.0017879211685999999, :CH2IBr => 0.004009730029999999,
        :I2O3 => 0.042764693289999996, :IBr => 0.074549024996,
        :ICl => 0.02410865173396,
    )
    for (s, val) in expected
        jf = getfield(GasChem, Symbol("j_mean_", s))
        @test jf(298.0, fluxes) ≈ val rtol = 1.0e-3
    end
    # Three of them are temperature-interpolated and must actually vary with T.
    for s in (:ClNO2, :CH2I2, :CH2ICl, :CH2IBr)
        jf = getfield(GasChem, Symbol("j_mean_", s))
        @test jf(210.0, fluxes) != jf(298.0, fluxes)
    end
    # N2O absorbs only below 214 nm, so it photolyzes in the stratosphere and not at the
    # surface.  Its j-variable had no defining equation before; this pins that it does now.
    @test GasChem.j_mean_N2O(298.0, fluxes) > 1.0e-7
    @test GasChem.j_mean_N2O(298.0, get_fluxes(3600 * 12.0, 30.0, 0.0, 101325.0)) < 1.0e-15
end

#   The cross-section tables and the flux table have to come from the same 18-bin grid.
#   GEOS-Chem's tables put bins 17 and 18 at 380 nm and 574 nm; the example tables shipped
#   in the geoschem/Cloud-J repository re-bin them to 429 nm and 631 nm with a 3.05x larger
#   bin-17 flux.  Mixing the two silently rescales bin 17 by 0.328 and bin 18 by 1.432.
@testitem "Fast-JX bins 17/18 are on GEOS-Chem's grid" begin
    using GasChem
    @test GasChem.WL[17] == 380 && GasChem.WL[18] == 574
    @test GasChem.top_flux[17] ≈ 1.547e16 && GasChem.top_flux[18] ≈ 2.131e17
    # spot values straight out of ExtData/CHEM_INPUTS/CLOUD_J/v2024-09/FJX_spec.dat
    @test GasChem.σ_NO2_interp[17](300.0) ≈ 4.643e-19 rtol = 1.0e-3
    @test GasChem.σ_NO2_interp[18](300.0) ≈ 4.345e-22 rtol = 1.0e-3
    @test GasChem.σ_CH3OOH[17] ≈ 6.973e-23 rtol = 1.0e-3
    @test GasChem.σ_HOCl[17] ≈ 6.529e-21 rtol = 1.0e-3
    @test GasChem.σ_O3_interp[18](298.0) ≈ 1.666e-21 rtol = 1.0e-3
    @test GasChem.σ_HNO4_interp[18](300.0) ≈ 4.694e-23 rtol = 1.0e-3
end

@testitem "solar_flux_factor matches GEOS-Chem SOLFX" begin
    using Dates
    # SOLF = 1 - 0.034*cos((DOY-172)*2pi/365)  (fast_jx_mod.F90 SOLAR_JX):
    # minimum at the aphelion-side solstice (DOY 172), maximum near perihelion.
    t172 = datetime2unix(DateTime(2016, 6, 20, 12))
    t355 = datetime2unix(DateTime(2016, 12, 20, 12))
    @test GasChem.solar_flux_factor(t172) ≈ 0.966 atol = 1.0e-3
    @test GasChem.solar_flux_factor(t355) ≈ 1.0339 atol = 1.0e-3
    # energy-neutral over a full year
    days = [datetime2unix(DateTime(2016, 1, 1) + Day(d)) for d in 0:364]
    @test sum(GasChem.solar_flux_factor.(days)) / 365 ≈ 1.0 atol = 2.0e-3
end

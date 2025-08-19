@testsnippet FastJXSetup begin
    using ModelingToolkit
    using SymbolicIndexingInterface: setp, getsym, parameter_values

    function get_fluxes(t, lat, lon, P)
        cos_sza = GasChem.cos_solar_zenith_angle(t, lat, lon)
        return GasChem.calc_direct_fluxes(cos_sza, P)
    end

    fj = mtkcompile(FastJX(0.0))
    test_time = 3600 * 12.0
    # [t_ref, lat, long, T, P, H2O]
    p = [0.0, 40.0, -97.0, 298.0, 101325.0, 450.0]
    prob = ODEProblem(fj,
        [fj.t_ref => 0.0, fj.lat => 40.0, fj.long => -97.0,
            fj.T => 298.0, fj.P => 101325.0, fj.H2O => 450.0],
        (test_time, test_time + 1))
end

#   Unit Test 0: O3 -> O2 + O(1D)

@testitem "o31D" setup=[FastJXSetup] begin
    u_0 = [
        0.007369280845991884, 
        0.007374722096550664, 
        0.007583971820745277, 
        0.007583971820745277
    ]
    fluxes = get_fluxes(3600 * 12.0, 30.0, 0.0, 0.9)
    test_0 = [
        GasChem.j_mean_O31D(100.0, fluxes),
        GasChem.j_mean_O31D(220.0, fluxes),
        GasChem.j_mean_O31D(300.0, fluxes),
        GasChem.j_mean_O31D(400.0, fluxes)
    ]
    @test test_0≈u_0 rtol=1e-6
end

#   Unit Test 1: H2O2 -> OH + OH
@testitem "H2O2" setup=[FastJXSetup] begin
    u_1 = [9.537314938390102e-5, 9.762101904345852e-5, 9.986888870301603e-5]

    fluxes = get_fluxes(3600 * 12.0, 30.0, 0.0, 0.9)
    test_1 = [
        GasChem.j_mean_H2O2(150.0, fluxes),
        GasChem.j_mean_H2O2(250.0, fluxes),
        GasChem.j_mean_H2O2(350.0, fluxes)
    ]

    @test test_1 ≈ u_1

    j_H2O2_func = getsym(prob, fj.j_h2o2)
    j_H2O2_value = j_H2O2_func(prob)
    j_want = GasChem.j_mean_H2O2(298.0, get_fluxes(3600 * 12.0, 40.0, -97.0, 101325)) 
    @test j_H2O2_value≈j_want rtol=0.004
end

# Unit Test 2: CH2O -> H + HO2 + CO
@testitem "CH2Oa" setup=[FastJXSetup] begin
    u_2 = [8.642676785125311e-5, 8.64227156795509e-5, 8.641551181874694e-5]

    fluxes = get_fluxes(3600 * 12.0, 30.0, 0.0, 0.9)
    test_2 = [
        GasChem.j_mean_H2COa(200.0, fluxes),
        GasChem.j_mean_H2COa(250.0, fluxes),
        GasChem.j_mean_H2COa(300.0, fluxes)
    ]

    @test test_2 ≈ u_2


    j_CH2Oa_func = getsym(prob, fj.j_CH2Oa)
    j_CH2Oa_value = j_CH2Oa_func(prob)
    @test j_CH2Oa_value≈GasChem.j_mean_CH2Oa(298.0,
        get_fluxes(3600 * 12.0, 40.0, -97.0, 101325)
    ) rtol=1e-6
end

@testitem "CH2Ob" setup=[FastJXSetup] begin
    u_2 = [7.16264584129105e-5, 7.166900731492678e-5, 7.174464980740016e-5]

    fluxes = get_fluxes(3600 * 12.0, 30.0, 0.0, 0.9)
    test_2 = [
        GasChem.j_mean_H2COb(200.0, fluxes),
        GasChem.j_mean_H2COb(250.0, fluxes),
        GasChem.j_mean_H2COb(300.0, fluxes)
    ]

    @test test_2 ≈ u_2

    j_CH2Ob_func = getsym(prob, fj.j_CH2Ob)
    j_CH2Ob_value = j_CH2Ob_func(prob)
    @test j_CH2Ob_value≈
    GasChem.j_mean_CH2Ob(
        298.0,
        get_fluxes(3600 * 12.0, 40.0, -97.0, 101325)
    ) rtol=1e-6
end

# Unit Test 3: CH3OOH -> OH + HO2 + CH2O
@testitem "CH3OOH" setup=[FastJXSetup] begin
    u_3 = [5.406743321900099e-5, 5.406743321900099e-5, 5.406743321900099e-5]

    test_3 = [
        GasChem.j_mean_CH3OOH(200.0, get_fluxes(3600 * 6.0, 30.0, 0.0, 0.9)),
        GasChem.j_mean_CH3OOH(200.0, get_fluxes(3600 * 12.0, 30.0, 0.0, 0.9)),
        GasChem.j_mean_CH3OOH(200.0, get_fluxes(3600 * 18.0, 30.0, 0.0, 0.9))
    ]

    @test test_3 ≈ u_3

    j_CH3OOH_func = getsym(prob, fj.j_CH3OOH)
    j_CH3OOH_value = j_CH3OOH_func(prob)
    @test j_CH3OOH_value≈
    GasChem.j_mean_CH3OOH(
        298.0,
        get_fluxes(3600 * 12.0, 40.0, -97.0, 101325)
    ) rtol=1e-6
end

# Unit Test 4: NO2 -> NO + O
@testitem "NO2" setup=[FastJXSetup] begin
    u_4 = [0.003926795211288372, 0.004094143741247612, 0.004261492271206852]

    fluxes = get_fluxes(3600 * 12.0, 30.0, 0.0, 0.9)
    test_4 = [
        GasChem.j_mean_NO2(150.0, fluxes),
        GasChem.j_mean_NO2(250.0, fluxes),
        GasChem.j_mean_NO2(300.0, fluxes)
    ]

    @test test_4≈u_4 rtol=1e-6

    j_NO2_func = getsym(prob, fj.j_NO2)
    j_NO2_value = j_NO2_func(prob)

    @test j_NO2_value≈
    GasChem.j_mean_NO2(298.0, get_fluxes(3600 * 12.0, 40.0, -97.0, 101325)) rtol=1e-6
end

@testitem "GEOS-Chem: CFCl3, H1301, Glyxlc" setup=[FastJXSetup] begin
    # [t_ref, lat, long, T, P, H2O]
    p = [0.0, 40.0, -97.0, 298.0, 101325.0, 450.0]
    j_CFCl3_func = getsym(prob, fj.j_CFCl3)
    j_CFCl3_value = j_CFCl3_func(prob)

    @test j_CFCl3_value≈
    GasChem.j_mean_CFCl3(298.0, get_fluxes(3600 * 12.0, 40.0, -97.0, 101325)) rtol=1e-6

    j_H1301_func = getsym(prob, fj.j_H1301)
    j_H1301_value = j_H1301_func(prob)

    @test j_H1301_value≈
    GasChem.j_mean_H1301(298.0, get_fluxes(3600 * 12.0, 40.0, -97.0, 101325)) rtol=1e-6

    j_Glyxlc_func = getsym(prob, fj.j_Glyxlc)
    j_Glyxlc_value = j_Glyxlc_func(prob)

    @test j_Glyxlc_value≈
    GasChem.j_mean_Glyxlc(298.0, get_fluxes(3600 * 12.0, 40.0, -97.0, 101325)) rtol=1e-6
end

@testitem "Ensure Cos SZA is non-allocating" begin
    using AllocCheck
    @check_allocs checkcos(lat, t, long) = GasChem.cos_solar_zenith_angle(t, lat, long)
    checkcos(0.0, 0.0, 0.0)
    checkcos(0.0f0, 0.0f0, 0.0f0)
end

@testitem "FastJX Initialization" begin
    using ModelingToolkit
    @test_nowarn mtkcompile(GasChem.FastJX(0.0))
end

@testitem "Direct Flux" begin
    using ModelingToolkit
    @parameters P, csa
    x = GasChem.calc_direct_fluxes(csa, P)
    @test substitute(x, Dict(P => 1013525, csa => 0.42255961917649837)) ≈ [
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
        1.6978055609093792e17
    ]
end

@testitem "Direct Flux 2" begin
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

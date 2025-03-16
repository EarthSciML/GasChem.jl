using GasChem
using ModelingToolkit
using Test
using AllocCheck

function get_fluxes(t, lat, lon, P)
    cos_sza = GasChem.cos_solar_zenith_angle(t, lat, lon)
    return GasChem.calc_direct_fluxes(cos_sza, P)
end

#   Unit Test 0: O3 -> O2 + O(1D)
@testset "o31D" begin
    u_0 = [
        3.008127306705564e15,
        3.0354776971213985e15,
        5.082768676900143e15,
        6.079063776373681e15,
    ]
    fluxes = get_fluxes(3600 * 12.0, 30.0, 0.0, 0.9)
    test_0 = [
        GasChem.j_mean_o31D(100.0, fluxes),
        GasChem.j_mean_o31D(220.0, fluxes),
        GasChem.j_mean_o31D(300.0, fluxes),
        GasChem.j_mean_o31D(400.0, fluxes),
    ]
    @test test_0 ≈ u_0 rtol = 1e-6
end


#   Unit Test 1: H2O2 -> OH + OH
@testset "H2O2" begin
    u_1 = [9.55610917750852e-5, 9.78407829679127e-5, 0.0001001204751748322]

    fluxes = get_fluxes(3600 * 12.0, 30.0, 0.0, 0.9)
    test_1 = [
        GasChem.j_mean_H2O2(150.0f0, fluxes),
        GasChem.j_mean_H2O2(250.0f0, fluxes),
        GasChem.j_mean_H2O2(350.0f0, fluxes),
    ]

    @test test_1 ≈ u_1
end

# Unit Test 2: CH2O -> H + HO2 + CO
@testset "CH2Oa" begin
    u_2 = [8.642676369563976e-5, 8.642271238446067e-5, 8.641551005347563e-5]

    fluxes = get_fluxes(3600 * 12.0, 30.0, 0.0, 0.9)
    test_2 = [
        GasChem.j_mean_CH2Oa(200.0, fluxes),
        GasChem.j_mean_CH2Oa(250.0, fluxes),
        GasChem.j_mean_CH2Oa(300.0, fluxes),
    ]

    @test test_2 ≈ u_2
end

@testset "CH2Ob" begin
    u_2 = [4.440814277628268e-5, 4.443217170426024e-5, 4.447488979844258e-5]

    cos_sza = GasChem.cos_solar_zenith_angle(3600 * 12.0, 30.0, 0.0)
    fluxes = GasChem.calc_fluxes(cos_sza)
    test_2 = [
        GasChem.j_mean_CH2Ob(200.0, fluxes),
        GasChem.j_mean_CH2Ob(250.0, fluxes),
        GasChem.j_mean_CH2Ob(300.0, fluxes),
    ]

    @test test_2 ≈ u_2
end

# Unit Test 3: CH3OOH -> OH + HO2 + CH2O
@testset "CH3OOH" begin
    u_3 = [5.479266623743904e-5, 5.479266623743904e-5, 5.479266623743904e-5]

    test_3 = [
        GasChem.j_mean_CH3OOH(200.0, get_fluxes(3600 * 6.0, 30.0, 0.0, 0.9)),
        GasChem.j_mean_CH3OOH(200.0, get_fluxes(3600 * 12.0, 30.0, 0.0, 0.9)),
        GasChem.j_mean_CH3OOH(200.0, get_fluxes(3600 * 18.0, 30.0, 0.0, 0.9)),
    ]

    @test test_3 ≈ u_3
end


# Unit Test 4: NO2 -> NO + O
@testset "NO2" begin
    u_4 = [0.008441071595788706, 0.00881068300155889, 0.00913594103863665]

    fluxes = get_fluxes(3600 * 12.0, 30.0, 0.0, 0.9)
    test_4 = [
        GasChem.j_mean_NO2(150.0, fluxes),
        GasChem.j_mean_NO2(250.0, fluxes),
        GasChem.j_mean_NO2(300.0, fluxes),
    ]

    @test test_4 ≈ u_4 rtol = 1e-6
end

@testset "Ensure Cos SZA is non-allocating" begin
    @check_allocs checkcos(lat, t, long) = GasChem.cos_solar_zenith_angle(t, lat, long)
    checkcos(0.0, 0.0, 0.0)
    checkcos(0.0f0, 0.0f0, 0.0f0)
end

@testset "FastJX Initialization" begin
    fj = GasChem.FastJX()
    @test_nowarn structural_simplify(fj)
end

@testset "Direct Flux" begin
    @parameters P, csa
    x = GasChem.calc_direct_fluxes(csa, P)
    @test substitute(x, Dict(P => 1013525, csa => 0.42255961917649837)) ≈ [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 6.942926208560752e-20, 9.324385539750345e11, 1.6175240854575052e10, 6.451788954313698e10, 2.65230857303719e13, 4.949172766139693e13, 4.1853350675689734e13, 9.0818932182367e13, 6.106135960927899e14, 5.384781489955934e15, 1.75419654839834e17]
end

@testset "Direct Flux 2" begin
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
    @test GasChem.calc_direct_flux(cos_sza, P, 1) ≈ 3.1945754293116e6
    P = 500
    @test GasChem.calc_direct_flux(cos_sza, P, 1) ≈ 1.681154736397082e-16
end

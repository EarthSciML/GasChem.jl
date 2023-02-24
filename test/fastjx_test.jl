using GasChem
using Test

#   Unit Test 0: O3 -> O2 + O(1D)
@testset "o31D" begin
    u_0 = [
        1.4073075332238303e14,
        1.416409373022159e14,
        3.070717954604601e14,
        3.888770405597495e14,
    ]

    test_0 = [
        GasChem.j_mean_o31D(3600 * 12.0, 30.0, 100.0),
        GasChem.j_mean_o31D(3600 * 12.0, 30.0, 220.0),
        GasChem.j_mean_o31D(3600 * 12.0, 30.0, 300.0),
        GasChem.j_mean_o31D(3600 * 12.0, 30.0, 400.0),
    ]

    @test test_0 ≈ u_0 rtol = 1e-6
end

#   Unit Test 1: H2O2 -> OH + OH
@testset "H2O2" begin
    u_1 = [5.455259236723579e-6, 5.589614033487323e-6, 5.723968830251066e-6]

    test_1 = [
        GasChem.j_mean_H2O2(3600 * 12.f0, 30.0f0, 150.f0),
        GasChem.j_mean_H2O2(3600 * 12.f0, 30.f0, 250.f0),
        GasChem.j_mean_H2O2(3600 * 12.f0, 30.f0, 350.f0),
    ]

    @test test_1 ≈ u_1
end

# Unit Test 2: CH2O -> H + HO2 + CO
@testset "CH2O" begin
    u_2 = [3.7478813392330445e-6, 3.7481177225578875e-6, 3.7485379595798317e-6]

    test_2 = [
        GasChem.j_mean_CH2Oa(3600 * 12., 30., 200.),
        GasChem.j_mean_CH2Oa(3600 * 12., 30., 250.),
        GasChem.j_mean_CH2Oa(3600 * 12., 30., 300.),
    ]

    @test test_2 ≈ u_2
end

# Unit Test 3: CH3OOH -> OH + HO2 + CH2O
@testset "CH3OOH" begin
    u_3 = [1.1100414728144508e-7, 3.3305186094197368e-6, 1.1100414728144508e-7]

    test_3 = [
        GasChem.j_mean_CH3OOH(3600 * 6., 30., 200.),
        GasChem.j_mean_CH3OOH(3600 * 12., 30., 200.),
        GasChem.j_mean_CH3OOH(3600 * 18., 30., 200.),
    ]

    @test test_3 ≈ u_3
end

# Unit Test 4: NO2 -> NO + O
@testset "NO2" begin
    u_4 = [0.0020197327512508217, 0.0021057190609635253, 0.0021813870135107055]

    test_4 = [
        GasChem.j_mean_NO2(3600 * 12., 30., 150.),
        GasChem.j_mean_NO2(3600 * 12., 30., 250.),
        GasChem.j_mean_NO2(3600 * 12., 30., 300.),
    ]

    @test test_4 ≈ u_4 rtol = 1e-6
end

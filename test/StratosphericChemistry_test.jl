"""
    Stratospheric Chemistry Tests

Comprehensive tests for the stratospheric chemistry implementation, verifying
against Chapter 5 of Seinfeld & Pandis (2006) "Atmospheric Chemistry and
Physics: From Air Pollution to Climate Change", 2nd Edition.
"""

@testsnippet StratSetup begin
    using Test
    using ModelingToolkit
    using OrdinaryDiffEqDefault
    using OrdinaryDiffEqRosenbrock
    using GasChem

    # CGS to SI conversion factors
    const CGS_TO_SI_CONC = 1e6  # molec/cm^3 → molec/m^3
    const CGS_TO_SI_K2 = 1e-12  # cm^6/molec^2/s → m^6/s
    const CGS_TO_SI_K = 1e-6    # cm^3/molec/s → m^3/s
end

# =============================================================================
# Rate Coefficient Functions — now returning SI values
# =============================================================================

@testitem "Rate Coefficient k_O_O2_M" setup=[StratSetup] tags=[:stratospheric] begin
    # Reference: Table B.2, k = 6.0 x 10^-34 (T/300)^-2.4 cm^6/molec^2/s
    # Function now returns SI: × 1e-12

    T = 300.0
    @test isapprox(GasChem.k_O_O2_M(T), 6.0e-34 * CGS_TO_SI_K2, rtol = 1e-10)

    T = 227.0
    @test isapprox(GasChem.k_O_O2_M(T), 6.0e-34 * (227.0 / 300.0)^(-2.4) * CGS_TO_SI_K2, rtol = 1e-10)

    # Rate increases at lower temperatures (negative exponent)
    @test GasChem.k_O_O2_M(200.0) > GasChem.k_O_O2_M(300.0)
end

@testitem "Rate Coefficient k_O_O3" setup=[StratSetup] tags=[:stratospheric] begin
    # Reference: Table B.1, k = 8.0 x 10^-12 exp(-2060/T) cm^3/molec/s
    # Function now returns SI: × 1e-6

    T = 227.0
    @test isapprox(GasChem.k_O_O3(T), 8.0e-12 * exp(-2060.0 / T) * CGS_TO_SI_K, rtol = 1e-10)

    T = 300.0
    @test isapprox(GasChem.k_O_O3(T), 8.0e-12 * exp(-2060.0 / T) * CGS_TO_SI_K, rtol = 1e-10)

    # Rate increases with temperature (Arrhenius behavior)
    @test GasChem.k_O_O3(300.0) > GasChem.k_O_O3(200.0)
end

@testitem "Rate Coefficient k_NO_O3" setup=[StratSetup] tags=[:stratospheric] begin
    # Reference: Page 154, k = 3.0 × 10⁻¹² exp(-1500/T) cm³/molec/s
    T = 222.0
    @test isapprox(GasChem.k_NO_O3(T), 3.0e-12 * exp(-1500.0 / T) * CGS_TO_SI_K, rtol = 1e-10)
    @test GasChem.k_NO_O3(300.0) > GasChem.k_NO_O3(200.0)
end

@testitem "Rate Coefficient k_NO2_O" setup=[StratSetup] tags=[:stratospheric] begin
    # Reference: Page 154, k = 5.6 × 10⁻¹² exp(180/T) cm³/molec/s
    T = 237.0
    @test isapprox(GasChem.k_NO2_O(T), 5.6e-12 * exp(180.0 / T) * CGS_TO_SI_K, rtol = 1e-10)
end

@testitem "Rate Coefficient k_OH_O3" setup=[StratSetup] tags=[:stratospheric] begin
    # Reference: Page 161, k = 1.7 × 10⁻¹² exp(-940/T) cm³/molec/s
    T = 227.0
    @test isapprox(GasChem.k_OH_O3(T), 1.7e-12 * exp(-940.0 / T) * CGS_TO_SI_K, rtol = 1e-10)
end

@testitem "Rate Coefficient k_Cl_O3" setup=[StratSetup] tags=[:stratospheric] begin
    # Reference: Page 162, k = 2.3 × 10⁻¹¹ exp(-200/T) cm³/molec/s
    T = 251.0
    @test isapprox(GasChem.k_Cl_O3(T), 2.3e-11 * exp(-200.0 / T) * CGS_TO_SI_K, rtol = 1e-10)
end

@testitem "Rate Coefficient k_ClO_O" setup=[StratSetup] tags=[:stratospheric] begin
    # Reference: Page 162, k = 3.0 × 10⁻¹¹ exp(70/T) cm³/molec/s
    T = 251.0
    @test isapprox(GasChem.k_ClO_O(T), 3.0e-11 * exp(70.0 / T) * CGS_TO_SI_K, rtol = 1e-10)
end

@testitem "Rate Coefficient k_OH_HCl" setup=[StratSetup] tags=[:stratospheric] begin
    # Reference: Page 168, k = 2.6 × 10⁻¹² exp(-350/T) cm³/molec/s
    T = 237.0
    @test isapprox(GasChem.k_OH_HCl(T), 2.6e-12 * exp(-350.0 / T) * CGS_TO_SI_K, rtol = 1e-10)
end

@testitem "Rate Coefficient k_HO2_O3" setup=[StratSetup] tags=[:stratospheric] begin
    # Reference: Table B.1, k = 1.0 × 10⁻¹⁴ exp(-490/T) cm³/molec/s
    T = 227.0
    @test isapprox(GasChem.k_HO2_O3(T), 1.0e-14 * exp(-490.0 / T) * CGS_TO_SI_K, rtol = 1e-10)
    @test GasChem.k_HO2_O3(300.0) > GasChem.k_HO2_O3(200.0)
end

@testitem "Rate Coefficient k_Cl_CH4" setup=[StratSetup] tags=[:stratospheric] begin
    # Reference: Table B.1, k = 9.6 × 10⁻¹² exp(-1360/T) cm³/molec/s
    T = 227.0
    @test isapprox(GasChem.k_Cl_CH4(T), 9.6e-12 * exp(-1360.0 / T) * CGS_TO_SI_K, rtol = 1e-10)
    @test GasChem.k_Cl_CH4(300.0) > GasChem.k_Cl_CH4(200.0)
end

@testitem "Rate Coefficient k_Br_O3" setup=[StratSetup] tags=[:stratospheric] begin
    # Reference: Table B.1, k = 1.7 × 10⁻¹¹ exp(-800/T) cm³/molec/s
    T = 227.0
    @test isapprox(GasChem.k_Br_O3(T), 1.7e-11 * exp(-800.0 / T) * CGS_TO_SI_K, rtol = 1e-10)
    @test GasChem.k_Br_O3(300.0) > GasChem.k_Br_O3(200.0)
end

@testitem "Rate Coefficient k_HO2_O" setup=[StratSetup] tags=[:stratospheric] begin
    # Reference: Page 161, k = 3.0 × 10⁻¹¹ exp(200/T) cm³/molec/s
    T = 227.0
    @test isapprox(GasChem.k_HO2_O(T), 3.0e-11 * exp(200.0 / T) * CGS_TO_SI_K, rtol = 1e-10)
end

@testitem "Rate Coefficient k_HO2_NO" setup=[StratSetup] tags=[:stratospheric] begin
    # Reference: Page 158, k = 3.5 × 10⁻¹² exp(250/T) cm³/molec/s
    T = 227.0
    @test isapprox(GasChem.k_HO2_NO(T), 3.5e-12 * exp(250.0 / T) * CGS_TO_SI_K, rtol = 1e-10)
end

@testitem "Rate Coefficient k_ClO_NO" setup=[StratSetup] tags=[:stratospheric] begin
    # Reference: Page 163, k = 6.4 × 10⁻¹² exp(290/T) cm³/molec/s
    T = 227.0
    @test isapprox(GasChem.k_ClO_NO(T), 6.4e-12 * exp(290.0 / T) * CGS_TO_SI_K, rtol = 1e-10)
end

@testitem "Rate Coefficient k_BrO_ClO" setup=[StratSetup] tags=[:stratospheric] begin
    # Reference: Table B.1, BrO + ClO channels
    T = 227.0
    @test isapprox(GasChem.k_BrO_ClO_BrCl(T), 5.8e-13 * exp(170.0 / T) * CGS_TO_SI_K, rtol = 1e-10)
    @test isapprox(GasChem.k_BrO_ClO_ClOO(T), 2.3e-12 * exp(260.0 / T) * CGS_TO_SI_K, rtol = 1e-10)
    # Temperature dependence (positive activation → rate decreases with T)
    @test GasChem.k_BrO_ClO_BrCl(300.0) < GasChem.k_BrO_ClO_BrCl(200.0)
    @test GasChem.k_BrO_ClO_ClOO(300.0) < GasChem.k_BrO_ClO_ClOO(200.0)
end

@testitem "Rate Coefficient k_O1D_M" setup=[StratSetup] tags=[:stratospheric] begin
    # Reference: Page 143
    T = 227.0
    @test isapprox(GasChem.k_O1D_M(T, :O2), 3.2e-11 * exp(70.0 / T) * CGS_TO_SI_K, rtol = 1e-10)
    @test isapprox(GasChem.k_O1D_M(T, :N2), 1.8e-11 * exp(110.0 / T) * CGS_TO_SI_K, rtol = 1e-10)
    # Air weighted average
    k_air = (0.21 * 3.2e-11 * exp(70.0 / T) + 0.79 * 1.8e-11 * exp(110.0 / T)) * CGS_TO_SI_K
    @test isapprox(GasChem.k_O1D_M(T, :air), k_air, rtol = 1e-10)
end

# =============================================================================
# Structural Verification Tests
# =============================================================================

@testitem "Chapman Mechanism Structure" setup=[StratSetup] tags=[:stratospheric] begin
    sys = ChapmanMechanism()

    # 2 ODEs (O, O3) + 1 algebraic (Ox)
    @test length(equations(sys)) == 3

    states = unknowns(sys)
    state_names = [string(s) for s in states]
    @test any(occursin("O(t)", n) for n in state_names)
    @test any(occursin("O3(t)", n) for n in state_names)
    @test any(occursin("Ox(t)", n) for n in state_names)
end

@testitem "NOx Cycle Structure" setup=[StratSetup] tags=[:stratospheric] begin
    sys = NOxCycle()
    @test length(equations(sys)) == 3
    states = unknowns(sys)
    state_names = [string(s) for s in states]
    @test any(occursin("NO(t)", n) for n in state_names)
    @test any(occursin("NO2(t)", n) for n in state_names)
    @test any(occursin("NOx(t)", n) for n in state_names)
end

@testitem "HOx Cycle Structure" setup=[StratSetup] tags=[:stratospheric] begin
    sys = HOxCycle()
    @test length(equations(sys)) == 3
    states = unknowns(sys)
    state_names = [string(s) for s in states]
    @test any(occursin("OH(t)", n) for n in state_names)
    @test any(occursin("HO2(t)", n) for n in state_names)
    @test any(occursin("HOx(t)", n) for n in state_names)
end

@testitem "ClOx Cycle Structure" setup=[StratSetup] tags=[:stratospheric] begin
    sys = ClOxCycle()
    # 4 ODEs (Cl, ClO, HCl, ClONO2) + 2 algebraic (ClOx, Cly)
    @test length(equations(sys)) == 6
    states = unknowns(sys)
    state_names = [string(s) for s in states]
    @test any(occursin("Cl(t)", n) for n in state_names)
    @test any(occursin("ClO(t)", n) for n in state_names)
    @test any(occursin("ClONO2(t)", n) for n in state_names)
    @test any(occursin("ClOx(t)", n) for n in state_names)
end

@testitem "BrOx Cycle Structure" setup=[StratSetup] tags=[:stratospheric] begin
    sys = BrOxCycle()
    # 3 ODEs (Br, BrO, HOBr) + 2 algebraic (BrOx, Bry)
    @test length(equations(sys)) == 5
    states = unknowns(sys)
    state_names = [string(s) for s in states]
    @test any(occursin("Br(t)", n) for n in state_names)
    @test any(occursin("BrO(t)", n) for n in state_names)
    @test any(occursin("BrOx(t)", n) for n in state_names)
end

@testitem "Comprehensive System Structure" setup=[StratSetup] tags=[:stratospheric] begin
    sys = StratosphericOzoneSystem()

    # 14 ODEs + 5 algebraic = 19 equations, 19 unknowns
    @test length(equations(sys)) == 19

    states = unknowns(sys)
    state_names = [string(s) for s in states]

    # Odd oxygen
    @test any(occursin("O(t)", n) for n in state_names)
    @test any(occursin("O3(t)", n) for n in state_names)
    @test any(occursin("O1D(t)", n) for n in state_names)

    # NOx
    @test any(occursin("NO(t)", n) for n in state_names)
    @test any(occursin("NO2(t)", n) for n in state_names)

    # HOx, ClOx, BrOx
    @test any(occursin("OH(t)", n) for n in state_names)
    @test any(occursin("HO2(t)", n) for n in state_names)
    @test any(occursin("Cl(t)", n) for n in state_names)
    @test any(occursin("ClO(t)", n) for n in state_names)
    @test any(occursin("Br(t)", n) for n in state_names)
    @test any(occursin("BrO(t)", n) for n in state_names)
end

# =============================================================================
# Integration Tests — values in SI (m^-3)
# =============================================================================

@testitem "Chapman Mechanism Integration" setup=[StratSetup] tags=[:stratospheric] begin
    sys = ChapmanMechanism()
    compiled_sys = mtkcompile(sys)

    # Initial conditions in SI units (m^-3)
    u0 = [compiled_sys.O => 1e13, compiled_sys.O3 => 3e18]
    tspan = (0.0, 86400.0)  # 1 day

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol = 1e-8, reltol = 1e-8)

    @test sol.retcode == ReturnCode.Success
    @test all(sol[compiled_sys.O] .> 0)
    @test all(sol[compiled_sys.O3] .> 0)
end

# =============================================================================
# Positivity Preservation Tests
# =============================================================================

@testitem "Positivity Preservation" setup=[StratSetup] tags=[:stratospheric] begin
    sys = ChapmanMechanism()
    compiled_sys = mtkcompile(sys)

    # SI units (m^-3)
    u0 = [compiled_sys.O => 1e13, compiled_sys.O3 => 3e18]
    tspan = (0.0, 3600.0)  # 1 hour

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol = 1e-10, reltol = 1e-10)

    @test sol.retcode == ReturnCode.Success
    @test all(sol[compiled_sys.O] .>= 0)
    @test all(sol[compiled_sys.O3] .>= 0)
end

# =============================================================================
# Conservation Law Tests
# =============================================================================

@testitem "Odd Oxygen Conservation (Chapman)" setup=[StratSetup] tags=[:stratospheric] begin
    sys = ChapmanMechanism()
    compiled_sys = mtkcompile(sys)

    # SI units (m^-3)
    u0 = [compiled_sys.O => 1e13, compiled_sys.O3 => 3e18]
    tspan = (0.0, 10.0)  # Very short time

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol = 1e-12, reltol = 1e-12)

    @test sol.retcode == ReturnCode.Success

    Ox_initial = sol[compiled_sys.O][1] + sol[compiled_sys.O3][1]
    Ox_final = sol[compiled_sys.O][end] + sol[compiled_sys.O3][end]
    @test isapprox(Ox_initial, Ox_final, rtol = 1e-4)
end

@testitem "Cly Conservation (ClOx Cycle)" setup=[StratSetup] tags=[:stratospheric] begin
    sys = ClOxCycle()
    compiled_sys = mtkcompile(sys)

    # SI units (m^-3)
    u0 = [
        compiled_sys.Cl => 1e10,
        compiled_sys.ClO => 1e13,
        compiled_sys.HCl => 1e15,
        compiled_sys.ClONO2 => 1e15
    ]
    tspan = (0.0, 1.0)  # Very short time

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol = 1e-12, reltol = 1e-12)

    @test sol.retcode == ReturnCode.Success

    Cly_init = sol[compiled_sys.Cl][1] + sol[compiled_sys.ClO][1] +
               sol[compiled_sys.HCl][1] + sol[compiled_sys.ClONO2][1]
    Cly_final = sol[compiled_sys.Cl][end] + sol[compiled_sys.ClO][end] +
                sol[compiled_sys.HCl][end] + sol[compiled_sys.ClONO2][end]

    @test isapprox(Cly_init, Cly_final, rtol = 1e-6)
end

# =============================================================================
# Steady-State Tests
# =============================================================================

@testitem "Chapman Steady-State [O]/[O3] Ratio (Eq. 5.7)" setup=[StratSetup] tags=[:stratospheric] begin
    # At steady state, [O]/[O3] = j_O3 / (k2 * [O2] * [M]) = j_O3 / (k2 * 0.21 * M^2)
    # Eq. 5.7 from Seinfeld & Pandis (2006)

    sys = ChapmanMechanism()
    compiled_sys = mtkcompile(sys)

    j_O3_val = 4e-4  # s^-1
    T_val = 227.0  # K
    # Rate coefficients now returned in SI directly
    k2_si = GasChem.k_O_O2_M(T_val)
    k4_si = GasChem.k_O_O3(T_val)
    M_si = 3.1e23  # m^-3 (3.1e17 cm^-3 × 1e6)
    O2_val = 0.21

    tspan = (0.0, 3600.0 * 24 * 30)  # 30 days

    prob = ODEProblem(compiled_sys,
        [compiled_sys.O => 1e11, compiled_sys.O3 => 1e16,
            compiled_sys.j_O2 => 1e-11, compiled_sys.j_O3 => j_O3_val,
            compiled_sys.k2 => k2_si, compiled_sys.k4 => k4_si,
            compiled_sys.M => M_si, compiled_sys.O2_mix => O2_val],
        tspan)
    sol = solve(prob, abstol = 1e-8, reltol = 1e-8)

    @test sol.retcode == ReturnCode.Success

    # Check steady-state ratio (Eq. 5.7): [O]/[O3] = j_O3 / (k2 * O2 * M^2)
    O_ss = sol[compiled_sys.O][end]
    O3_ss = sol[compiled_sys.O3][end]
    expected_ratio = j_O3_val / (k2_si * O2_val * M_si * M_si)
    actual_ratio = O_ss / O3_ss

    @test isapprox(actual_ratio, expected_ratio, rtol = 0.1)
end

# =============================================================================
# Analytical Solution Tests
# =============================================================================

@testitem "Chapman Steady-State Ozone (Eq. 5.13)" setup=[StratSetup] tags=[:stratospheric] begin
    # Equation 5.13: [O3]_ss = 0.21 * sqrt(k2 * j_O2 / (k4 * j_O3)) * [M]^(3/2)
    # Using SI units throughout

    T = 227.0
    M_si = 3.1e23  # m^-3
    k2_si = GasChem.k_O_O2_M(T)  # m^6/s (SI)
    k4_si = GasChem.k_O_O3(T)    # m^3/s (SI)
    j_O2 = 1e-11  # s^-1
    j_O3 = 4e-4   # s^-1

    O3_si = 0.21 * sqrt(k2_si * j_O2 / (k4_si * j_O3)) * M_si^1.5

    # The steady-state O3 should be on the order of 10^18 m^-3 at 30 km
    # (equivalent to ~10^12 molec/cm^3 × 1e6)
    @test O3_si > 1e16
    @test O3_si < 1e20
end

@testitem "Time to Steady State (Eq. 5.17)" setup=[StratSetup] tags=[:stratospheric] begin
    # Equation 5.17: tau_O3_ss = (1/4) * sqrt(k2*[M] / (k4 * j_O2 * j_O3))
    # Using SI units: k2 in m^6/s, M in m^-3, k4 in m^3/s

    # 30 km conditions
    T = 227.0
    M_si = 3.1e23  # m^-3
    k2 = GasChem.k_O_O2_M(T)  # m^6/s
    k4 = GasChem.k_O_O3(T)    # m^3/s
    j_O2 = 6e-11
    j_O3 = 1.2e-3

    tau = 0.25 * sqrt(k2 * M_si / (k4 * j_O2 * j_O3))
    tau_hours = tau / 3600.0

    # At 30 km, the reference value is ~160 hours
    @test tau_hours > 50   # order of magnitude check
    @test tau_hours < 500

    # 40 km conditions
    T40 = 251.0
    M40_si = 7.1e22  # m^-3 (7.1e16 cm^-3 × 1e6)
    k2_40 = GasChem.k_O_O2_M(T40)
    k4_40 = GasChem.k_O_O3(T40)
    j_O2_40 = 5e-10
    j_O3_40 = 1.9e-3

    tau_40 = 0.25 * sqrt(k2_40 * M40_si / (k4_40 * j_O2_40 * j_O3_40))
    tau_hours_40 = tau_40 / 3600.0

    # At 40 km, the reference value is ~12-40 hours
    @test tau_hours_40 > 1
    @test tau_hours_40 < 100

    # Timescale should decrease with altitude (faster at higher altitudes)
    @test tau_hours_40 < tau_hours
end

@testitem "[O]/[O3] Ratio Reference Values (Page 145)" setup=[StratSetup] tags=[:stratospheric] begin
    # Reference: Page 145, Seinfeld & Pandis (2006)
    # At z=30km: [O]/[O3] = 3.0 × 10^-5
    # At z=40km: [O]/[O3] = 9.4 × 10^-4

    # Using Eq. 5.7: [O]/[O3] = j_O3 / (k2 * [O2] * [M])
    # All in SI: k2 in m^6/s, M in m^-3

    # 30 km
    T30 = 227.0
    M30_si = 3.1e23  # m^-3
    k2_30 = GasChem.k_O_O2_M(T30)
    j_O3_30 = 1.2e-3  # Total O3 photolysis rate at 30 km

    ratio_30 = j_O3_30 / (k2_30 * 0.21 * M30_si * M30_si)

    # Should be on the order of 10^-5
    @test ratio_30 > 1e-6
    @test ratio_30 < 1e-3

    # 40 km
    T40 = 251.0
    M40_si = 7.1e22  # m^-3
    k2_40 = GasChem.k_O_O2_M(T40)
    j_O3_40 = 1.9e-3

    ratio_40 = j_O3_40 / (k2_40 * 0.21 * M40_si * M40_si)

    # Should be larger than at 30 km (lower M means higher ratio)
    @test ratio_40 > ratio_30
    @test ratio_40 > 1e-5
    @test ratio_40 < 1e-2
end

# =============================================================================
# Limiting Behavior Tests
# =============================================================================

@testitem "Zero Photolysis Rates" setup=[StratSetup] tags=[:stratospheric] begin
    sys = ChapmanMechanism()
    compiled_sys = mtkcompile(sys)

    # SI units (m^-3)
    tspan = (0.0, 3600.0)

    prob = ODEProblem(compiled_sys,
        [compiled_sys.O => 1e13, compiled_sys.O3 => 3e18,
            compiled_sys.j_O2 => 0.0, compiled_sys.j_O3 => 0.0],
        tspan)
    sol = solve(prob, abstol = 1e-10, reltol = 1e-10)

    @test sol.retcode == ReturnCode.Success

    # With no photolysis, O should decrease (consumed by O+O2+M and O+O3)
    @test sol[compiled_sys.O][end] < sol[compiled_sys.O][1]
end

# =============================================================================
# Rate Coefficients Dictionary
# =============================================================================

@testitem "Stratospheric Rate Coefficients Dictionary" setup=[StratSetup] tags=[:stratospheric] begin
    T = 227.0
    M_si = 3.1e23  # m^-3

    coeffs = stratospheric_rate_coefficients(T, M_si)

    expected_keys = [
        :k_O_O2_M, :k_O_O3, :k_O1D_M, :k_NO_O3, :k_NO2_O,
        :k_OH_O3, :k_HO2_O3, :k_HO2_O, :k_HO2_NO,
        :k_Cl_O3, :k_ClO_O, :k_ClO_NO, :k_Cl_CH4, :k_OH_HCl,
        :k_Br_O3, :k_BrO_ClO_BrCl, :k_BrO_ClO_ClOO, :k_ClO_ClO_M
    ]

    for key in expected_keys
        @test haskey(coeffs, key)
        @test coeffs[key] > 0
    end

    @test isapprox(coeffs[:k_O_O2_M], GasChem.k_O_O2_M(T), rtol = 1e-10)
    @test isapprox(coeffs[:k_O_O3], GasChem.k_O_O3(T), rtol = 1e-10)
end

# =============================================================================
# Comprehensive System Integration
# =============================================================================

@testitem "Comprehensive System Integration" setup=[StratSetup] tags=[:stratospheric] begin
    sys = StratosphericOzoneSystem()
    compiled_sys = mtkcompile(sys)

    prob = ODEProblem(compiled_sys, [], (0.0, 3600.0))  # 1 hour with defaults
    sol = solve(prob, Rodas5P(), abstol = 1e-8, reltol = 1e-8)

    @test sol.retcode == ReturnCode.Success

    # All concentrations should remain positive
    @test all(sol[compiled_sys.O] .>= 0)
    @test all(sol[compiled_sys.O3] .> 0)
    @test all(sol[compiled_sys.NO] .>= 0)
    @test all(sol[compiled_sys.NO2] .>= 0)
end

# =============================================================================
# Bry Conservation (BrOx Cycle)
# =============================================================================

@testitem "Bry Conservation (BrOx Cycle)" setup=[StratSetup] tags=[:stratospheric] begin
    sys = BrOxCycle()
    compiled_sys = mtkcompile(sys)

    # SI units (m^-3)
    u0 = [
        compiled_sys.Br => 1e11,
        compiled_sys.BrO => 1e12,
        compiled_sys.HOBr => 1e12
    ]
    tspan = (0.0, 1.0)  # Very short time

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol = 1e-12, reltol = 1e-12)

    @test sol.retcode == ReturnCode.Success

    Bry_init = sol[compiled_sys.Br][1] + sol[compiled_sys.BrO][1] +
               sol[compiled_sys.HOBr][1]
    Bry_final = sol[compiled_sys.Br][end] + sol[compiled_sys.BrO][end] +
                sol[compiled_sys.HOBr][end]

    @test isapprox(Bry_init, Bry_final, rtol = 1e-6)
end

# =============================================================================
# HOx Conservation (HOx Cycle)
# =============================================================================

@testitem "HOx Conservation (HOx Cycle)" setup=[StratSetup] tags=[:stratospheric] begin
    sys = HOxCycle()
    compiled_sys = mtkcompile(sys)

    # SI units (m^-3)
    u0 = [
        compiled_sys.OH => 1e12,
        compiled_sys.HO2 => 1e13
    ]
    tspan = (0.0, 1.0)  # Very short time

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol = 1e-12, reltol = 1e-12)

    @test sol.retcode == ReturnCode.Success

    HOx_init = sol[compiled_sys.OH][1] + sol[compiled_sys.HO2][1]
    HOx_final = sol[compiled_sys.OH][end] + sol[compiled_sys.HO2][end]

    @test isapprox(HOx_init, HOx_final, rtol = 1e-6)
end

# =============================================================================
# NOx Conservation (NOx Cycle)
# =============================================================================

@testitem "NOx Conservation (NOx Cycle)" setup=[StratSetup] tags=[:stratospheric] begin
    sys = NOxCycle()
    compiled_sys = mtkcompile(sys)

    # SI units (m^-3)
    u0 = [
        compiled_sys.NO => 1e15,
        compiled_sys.NO2 => 1e15
    ]
    tspan = (0.0, 1.0)  # Very short time

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol = 1e-12, reltol = 1e-12)

    @test sol.retcode == ReturnCode.Success

    NOx_init = sol[compiled_sys.NO][1] + sol[compiled_sys.NO2][1]
    NOx_final = sol[compiled_sys.NO][end] + sol[compiled_sys.NO2][end]

    @test isapprox(NOx_init, NOx_final, rtol = 1e-4)
end

# =============================================================================
# Heterogeneous N2O5 Rate Coefficient (Eq. 5.31) — now fully SI
# =============================================================================

@testitem "N2O5 Heterogeneous Rate (Eq. 5.31)" setup=[StratSetup] tags=[:stratospheric] begin
    # Reference: Eq. 5.31, Page 180
    # k = (gamma/4) * mean_speed * Ap
    # At T=220K, gamma=0.06, Ap=1e-8 cm^2/cm^3 = 1e-8 * 1e-4 m^2 / 1e-6 m^3 = 1e-6 m^-1
    # (Ap in CGS: cm^2/cm^3 = 1e-4 m^2 / 1e-6 m^3 = 100 m^-1 per unit of cm^2/cm^3)
    # So Ap = 1e-8 cm^2/cm^3 = 1e-8 * 100 m^-1 = 1e-6 m^-1

    gamma = 0.06
    T = 220.0
    Ap_si = 1e-6  # m^-1 (converted from 1e-8 cm^2/cm^3)

    k = GasChem.k_N2O5_H2O_het(gamma, T, Ap_si)

    @test k > 0
    # Expected: ~3.1e-6 s^-1 (from page 184)
    @test isapprox(k, 3.1e-6, rtol = 0.1)
end

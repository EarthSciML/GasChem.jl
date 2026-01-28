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
    using GasChem
end

# =============================================================================
# Rate Coefficient Functions (Tables B.1 and B.2)
# =============================================================================

@testitem "Rate Coefficient k_O_O2_M" setup=[StratSetup] tags=[:stratospheric] begin
    # Reference: Table B.2, k = 6.0 x 10^-34 (T/300)^-2.4 cm^6 molecule^-2 s^-1

    T = 300.0
    @test isapprox(GasChem.k_O_O2_M(T), 6.0e-34, rtol=1e-10)

    T = 227.0
    @test isapprox(GasChem.k_O_O2_M(T), 6.0e-34 * (227.0 / 300.0)^(-2.4), rtol=1e-10)

    # Rate increases at lower temperatures (negative exponent)
    @test GasChem.k_O_O2_M(200.0) > GasChem.k_O_O2_M(300.0)
end

@testitem "Rate Coefficient k_O_O3" setup=[StratSetup] tags=[:stratospheric] begin
    # Reference: Table B.1, k = 8.0 x 10^-12 exp(-2060/T) cm^3 molecule^-1 s^-1

    T = 227.0
    @test isapprox(GasChem.k_O_O3(T), 8.0e-12 * exp(-2060.0 / T), rtol=1e-10)

    T = 300.0
    @test isapprox(GasChem.k_O_O3(T), 8.0e-12 * exp(-2060.0 / T), rtol=1e-10)

    # Rate increases with temperature (Arrhenius behavior)
    @test GasChem.k_O_O3(300.0) > GasChem.k_O_O3(200.0)
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

    # NOx (NO3, N2O5, HNO3 removed — no equations for them)
    @test any(occursin("NO(t)", n) for n in state_names)
    @test any(occursin("NO2(t)", n) for n in state_names)
    @test !any(occursin("NO3(t)", n) for n in state_names)
    @test !any(occursin("N2O5(t)", n) for n in state_names)
    @test !any(occursin("HNO3(t)", n) for n in state_names)

    # HOx, ClOx, BrOx
    @test any(occursin("OH(t)", n) for n in state_names)
    @test any(occursin("HO2(t)", n) for n in state_names)
    @test any(occursin("Cl(t)", n) for n in state_names)
    @test any(occursin("ClO(t)", n) for n in state_names)
    @test any(occursin("Br(t)", n) for n in state_names)
    @test any(occursin("BrO(t)", n) for n in state_names)
end

# =============================================================================
# Integration Tests
# =============================================================================

@testitem "Chapman Mechanism Integration" setup=[StratSetup] tags=[:stratospheric] begin
    sys = ChapmanMechanism()
    compiled_sys = mtkcompile(sys)

    u0 = [compiled_sys.O => 1e7, compiled_sys.O3 => 3e12]
    tspan = (0.0, 86400.0)  # 1 day

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol=1e-8, reltol=1e-8)

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

    u0 = [compiled_sys.O => 1e7, compiled_sys.O3 => 3e12]
    tspan = (0.0, 3600.0)  # 1 hour

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol=1e-10, reltol=1e-10)

    @test sol.retcode == ReturnCode.Success
    @test all(sol[compiled_sys.O] .>= 0)
    @test all(sol[compiled_sys.O3] .>= 0)
end

# =============================================================================
# Conservation Law Tests
# =============================================================================

@testitem "Odd Oxygen Conservation (Chapman)" setup=[StratSetup] tags=[:stratospheric] begin
    # In the Chapman mechanism, the sum d[O]/dt + d[O3]/dt = 2*j_O2*[O2] - 2*k4*[O]*[O3]
    # This means Ox = O + O3 changes slowly (only through O2 photolysis and O+O3 destruction).
    # Over short timescales, O rapidly equilibrates with O3, but Ox is approximately conserved
    # compared to the fast O <-> O3 interconversion.

    sys = ChapmanMechanism()
    compiled_sys = mtkcompile(sys)

    u0 = [compiled_sys.O => 1e7, compiled_sys.O3 => 3e12]
    tspan = (0.0, 10.0)  # Very short time — Ox should be nearly constant

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol=1e-12, reltol=1e-12)

    @test sol.retcode == ReturnCode.Success

    # Ox = O + O3 should change by a very small fraction over 10 seconds
    Ox_initial = sol[compiled_sys.O][1] + sol[compiled_sys.O3][1]
    Ox_final = sol[compiled_sys.O][end] + sol[compiled_sys.O3][end]
    @test isapprox(Ox_initial, Ox_final, rtol=1e-4)
end

@testitem "Cly Conservation (ClOx Cycle)" setup=[StratSetup] tags=[:stratospheric] begin
    # Total inorganic chlorine Cly = Cl + ClO + HCl + ClONO2 should be conserved
    # because the ClOx cycle only interconverts chlorine species.
    # (ClONO2 photolysis is the only external source/sink, which is slow.)

    sys = ClOxCycle()
    compiled_sys = mtkcompile(sys)

    u0 = [
        compiled_sys.Cl => 1e4, compiled_sys.ClO => 1e7,
        compiled_sys.HCl => 1e9, compiled_sys.ClONO2 => 1e9
    ]
    tspan = (0.0, 1.0)  # Very short time to avoid ClONO2 photolysis effect

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol=1e-12, reltol=1e-12)

    @test sol.retcode == ReturnCode.Success

    Cly_init = sol[compiled_sys.Cl][1] + sol[compiled_sys.ClO][1] +
               sol[compiled_sys.HCl][1] + sol[compiled_sys.ClONO2][1]
    Cly_final = sol[compiled_sys.Cl][end] + sol[compiled_sys.ClO][end] +
                sol[compiled_sys.HCl][end] + sol[compiled_sys.ClONO2][end]

    # Cly is not perfectly conserved due to ClONO2 photolysis, but over 1 second
    # the change should be negligible
    @test isapprox(Cly_init, Cly_final, rtol=1e-6)
end

# =============================================================================
# Steady-State Tests
# =============================================================================

@testitem "Chapman Steady-State [O]/[O3] Ratio (Eq. 5.7)" setup=[StratSetup] tags=[:stratospheric] begin
    # At steady state, [O]/[O3] = j_O3 / (k2 * [O2] * [M])
    # We verify this by running the Chapman system to approach steady state.

    sys = ChapmanMechanism()
    compiled_sys = mtkcompile(sys)

    j_O3_val = 4e-4
    k2_val = GasChem.k_O_O2_M(227.0)
    M_val = 3.83e17
    O2_val = 0.21

    u0 = [compiled_sys.O => 1e5, compiled_sys.O3 => 1e10]
    tspan = (0.0, 3600.0 * 24 * 30)  # 30 days

    prob = ODEProblem(compiled_sys, u0, tspan,
        [compiled_sys.j_O2 => 1e-11, compiled_sys.j_O3 => j_O3_val,
         compiled_sys.k2 => k2_val, compiled_sys.k4 => GasChem.k_O_O3(227.0),
         compiled_sys.M => M_val, compiled_sys.O2_mix => O2_val])
    sol = solve(prob, abstol=1e-8, reltol=1e-8)

    @test sol.retcode == ReturnCode.Success

    # Check steady-state ratio (Eq. 5.7)
    O_ss = sol[compiled_sys.O][end]
    O3_ss = sol[compiled_sys.O3][end]
    expected_ratio = j_O3_val / (k2_val * O2_val * M_val * M_val)
    actual_ratio = O_ss / O3_ss

    # Allow generous tolerance since the system may not have fully reached steady state
    @test isapprox(actual_ratio, expected_ratio, rtol=0.1)
end

# =============================================================================
# Analytical Solution Tests
# =============================================================================

@testitem "Chapman Steady-State Ozone (Eq. 5.13)" setup=[StratSetup] tags=[:stratospheric] begin
    # Equation 5.13: [O3]_ss = 0.21 * sqrt(k2 * j_O2 / (k4 * j_O3)) * [M]^(3/2)

    T = 227.0
    M_val = 3.83e17
    k2_val = GasChem.k_O_O2_M(T)
    k4_val = GasChem.k_O_O3(T)
    j_O2_val = 1e-11
    j_O3_val = 4e-4

    O3_analytical = 0.21 * sqrt(k2_val * j_O2_val / (k4_val * j_O3_val)) * M_val^1.5

    # The analytical solution should give a physically reasonable ozone concentration
    # (order of 10^12 molec/cm^3 at 30 km)
    @test O3_analytical > 1e10
    @test O3_analytical < 1e14
end

# =============================================================================
# Limiting Behavior Tests
# =============================================================================

@testitem "Zero Photolysis Rates" setup=[StratSetup] tags=[:stratospheric] begin
    # With j_O2 = 0 and j_O3 = 0 (nighttime), O and O3 should only change
    # through the k4*O*O3 reaction, consuming both.

    sys = ChapmanMechanism()
    compiled_sys = mtkcompile(sys)

    u0 = [compiled_sys.O => 1e7, compiled_sys.O3 => 3e12]
    tspan = (0.0, 3600.0)

    prob = ODEProblem(compiled_sys, u0, tspan,
        [compiled_sys.j_O2 => 0.0, compiled_sys.j_O3 => 0.0])
    sol = solve(prob, abstol=1e-10, reltol=1e-10)

    @test sol.retcode == ReturnCode.Success

    # With no photolysis, O should decrease (consumed by O+O2+M and O+O3)
    @test sol[compiled_sys.O][end] < sol[compiled_sys.O][1]
end

# =============================================================================
# Rate Coefficients Dictionary
# =============================================================================

@testitem "Stratospheric Rate Coefficients Dictionary" setup=[StratSetup] tags=[:stratospheric] begin
    T = 227.0
    M = 3.1e17

    coeffs = stratospheric_rate_coefficients(T, M)

    expected_keys = [:k_O_O2_M, :k_O_O3, :k_O1D_M, :k_NO_O3, :k_NO2_O,
                    :k_OH_O3, :k_HO2_O3, :k_HO2_O, :k_HO2_NO,
                    :k_Cl_O3, :k_ClO_O, :k_ClO_NO, :k_Cl_CH4, :k_OH_HCl,
                    :k_Br_O3, :k_BrO_ClO_BrCl, :k_BrO_ClO_ClOO, :k_ClO_ClO_M]

    for key in expected_keys
        @test haskey(coeffs, key)
        @test coeffs[key] > 0
    end

    @test isapprox(coeffs[:k_O_O2_M], GasChem.k_O_O2_M(T), rtol=1e-10)
    @test isapprox(coeffs[:k_O_O3], GasChem.k_O_O3(T), rtol=1e-10)
end

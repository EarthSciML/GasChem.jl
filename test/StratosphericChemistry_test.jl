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
    const CGS_TO_SI_CONC = 1.0e6  # molec/cm^3 → molec/m^3
    const CGS_TO_SI_K2 = 1.0e-12  # cm^6/molec^2/s → m^6/s
    const CGS_TO_SI_K = 1.0e-6    # cm^3/molec/s → m^3/s

    # Inline rate coefficient computations matching the @constants in the source
    # These are used to verify analytical predictions against the source material
    function k_O_O2_M_si(T)
        # k = 6e-34 (T/300)^-2.4 cm^6/molec^2/s → × 1e-12 for SI
        return 6.0e-34 * (T / 300.0)^(-2.4) * CGS_TO_SI_K2
    end

    function k_O_O3_si(T)
        # k = 8e-12 exp(-2060/T) cm^3/molec/s → × 1e-6 for SI
        return 8.0e-12 * exp(-2060.0 / T) * CGS_TO_SI_K
    end
end

# =============================================================================
# Structural Verification Tests
# =============================================================================

@testitem "Chapman Mechanism Structure" setup = [StratSetup] tags = [:stratospheric] begin
    sys = ChapmanMechanism()

    # 2 ODEs (O, O3) + 1 algebraic (Ox)
    @test length(equations(sys)) == 3

    states = unknowns(sys)
    state_names = [string(s) for s in states]
    @test any(occursin("O(t)", n) for n in state_names)
    @test any(occursin("O3(t)", n) for n in state_names)
    @test any(occursin("Ox(t)", n) for n in state_names)
end

@testitem "NOx Cycle Structure" setup = [StratSetup] tags = [:stratospheric] begin
    sys = NOxCycle()
    @test length(equations(sys)) == 3
    states = unknowns(sys)
    state_names = [string(s) for s in states]
    @test any(occursin("NO(t)", n) for n in state_names)
    @test any(occursin("NO2(t)", n) for n in state_names)
    @test any(occursin("NOx(t)", n) for n in state_names)
end

@testitem "HOx Cycle Structure" setup = [StratSetup] tags = [:stratospheric] begin
    sys = HOxCycle()
    @test length(equations(sys)) == 3
    states = unknowns(sys)
    state_names = [string(s) for s in states]
    @test any(occursin("OH(t)", n) for n in state_names)
    @test any(occursin("HO2(t)", n) for n in state_names)
    @test any(occursin("HOx(t)", n) for n in state_names)
end

@testitem "ClOx Cycle Structure" setup = [StratSetup] tags = [:stratospheric] begin
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

@testitem "BrOx Cycle Structure" setup = [StratSetup] tags = [:stratospheric] begin
    sys = BrOxCycle()
    # 3 ODEs (Br, BrO, HOBr) + 2 algebraic (BrOx, Bry)
    @test length(equations(sys)) == 5
    states = unknowns(sys)
    state_names = [string(s) for s in states]
    @test any(occursin("Br(t)", n) for n in state_names)
    @test any(occursin("BrO(t)", n) for n in state_names)
    @test any(occursin("BrOx(t)", n) for n in state_names)
end

@testitem "Comprehensive System Structure" setup = [StratSetup] tags = [:stratospheric] begin
    sys = StratosphericOzoneSystem()

    # 15 ODEs (O, O1D, O3, N2O, NO, NO2, OH, HO2, Cl, ClO, HCl, ClONO2, Br, BrO, HOBr)
    # + 5 algebraic (Ox, NOx, HOx, ClOx, BrOx) = 20 equations
    @test length(equations(sys)) == 20

    states = unknowns(sys)
    state_names = [string(s) for s in states]

    # Odd oxygen
    @test any(occursin("O(t)", n) for n in state_names)
    @test any(occursin("O3(t)", n) for n in state_names)
    @test any(occursin("O1D(t)", n) for n in state_names)

    # N2O
    @test any(occursin("N2O(t)", n) for n in state_names)

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

    # Families
    @test any(occursin("Ox(t)", n) for n in state_names)
    @test any(occursin("NOx(t)", n) for n in state_names)
    @test any(occursin("HOx(t)", n) for n in state_names)
    @test any(occursin("ClOx(t)", n) for n in state_names)
    @test any(occursin("BrOx(t)", n) for n in state_names)
end

# =============================================================================
# Integration Tests — values in SI (m^-3)
# =============================================================================

@testitem "Chapman Mechanism Integration" setup = [StratSetup] tags = [:stratospheric] begin
    sys = ChapmanMechanism()
    compiled_sys = mtkcompile(sys)

    # Initial conditions in SI units (m^-3)
    u0 = [compiled_sys.O => 1.0e13, compiled_sys.O3 => 3.0e18]
    tspan = (0.0, 86400.0)  # 1 day

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol = 1.0e-8, reltol = 1.0e-8)

    @test sol.retcode == ReturnCode.Success
    @test all(sol[compiled_sys.O] .> 0)
    @test all(sol[compiled_sys.O3] .> 0)
end

@testitem "NOx Cycle Integration" setup = [StratSetup] tags = [:stratospheric] begin
    sys = NOxCycle()
    compiled_sys = mtkcompile(sys)

    u0 = [compiled_sys.NO => 1.0e15, compiled_sys.NO2 => 1.0e15]
    tspan = (0.0, 3600.0)  # 1 hour

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol = 1.0e-8, reltol = 1.0e-8)

    @test sol.retcode == ReturnCode.Success
    @test all(sol[compiled_sys.NO] .> 0)
    @test all(sol[compiled_sys.NO2] .> 0)
end

@testitem "HOx Cycle Integration" setup = [StratSetup] tags = [:stratospheric] begin
    sys = HOxCycle()
    compiled_sys = mtkcompile(sys)

    u0 = [compiled_sys.OH => 1.0e12, compiled_sys.HO2 => 1.0e13]
    tspan = (0.0, 3600.0)  # 1 hour

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol = 1.0e-8, reltol = 1.0e-8)

    @test sol.retcode == ReturnCode.Success
    @test all(sol[compiled_sys.OH] .> 0)
    @test all(sol[compiled_sys.HO2] .> 0)
end

@testitem "ClOx Cycle Integration" setup = [StratSetup] tags = [:stratospheric] begin
    sys = ClOxCycle()
    compiled_sys = mtkcompile(sys)

    u0 = [
        compiled_sys.Cl => 1.0e10,
        compiled_sys.ClO => 1.0e13,
        compiled_sys.HCl => 1.0e15,
        compiled_sys.ClONO2 => 1.0e15,
    ]
    tspan = (0.0, 3600.0)  # 1 hour

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol = 1.0e-8, reltol = 1.0e-8)

    @test sol.retcode == ReturnCode.Success
    @test all(sol[compiled_sys.Cl] .>= 0)
    @test all(sol[compiled_sys.ClO] .>= 0)
    @test all(sol[compiled_sys.HCl] .>= 0)
    @test all(sol[compiled_sys.ClONO2] .>= 0)
end

@testitem "BrOx Cycle Integration" setup = [StratSetup] tags = [:stratospheric] begin
    sys = BrOxCycle()
    compiled_sys = mtkcompile(sys)

    u0 = [
        compiled_sys.Br => 1.0e11,
        compiled_sys.BrO => 1.0e12,
        compiled_sys.HOBr => 1.0e12,
    ]
    tspan = (0.0, 3600.0)  # 1 hour

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol = 1.0e-8, reltol = 1.0e-8)

    @test sol.retcode == ReturnCode.Success
    @test all(sol[compiled_sys.Br] .>= 0)
    @test all(sol[compiled_sys.BrO] .>= 0)
    @test all(sol[compiled_sys.HOBr] .>= 0)
end

# =============================================================================
# Positivity Preservation Tests
# =============================================================================

@testitem "Positivity Preservation" setup = [StratSetup] tags = [:stratospheric] begin
    sys = ChapmanMechanism()
    compiled_sys = mtkcompile(sys)

    # SI units (m^-3)
    u0 = [compiled_sys.O => 1.0e13, compiled_sys.O3 => 3.0e18]
    tspan = (0.0, 3600.0)  # 1 hour

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol = 1.0e-10, reltol = 1.0e-10)

    @test sol.retcode == ReturnCode.Success
    @test all(sol[compiled_sys.O] .>= 0)
    @test all(sol[compiled_sys.O3] .>= 0)
end

# =============================================================================
# Conservation Law Tests
# =============================================================================

@testitem "Odd Oxygen Conservation (Chapman)" setup = [StratSetup] tags = [:stratospheric] begin
    sys = ChapmanMechanism()
    compiled_sys = mtkcompile(sys)

    # SI units (m^-3)
    u0 = [compiled_sys.O => 1.0e13, compiled_sys.O3 => 3.0e18]
    tspan = (0.0, 10.0)  # Very short time

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol = 1.0e-12, reltol = 1.0e-12)

    @test sol.retcode == ReturnCode.Success

    Ox_initial = sol[compiled_sys.O][1] + sol[compiled_sys.O3][1]
    Ox_final = sol[compiled_sys.O][end] + sol[compiled_sys.O3][end]
    @test isapprox(Ox_initial, Ox_final, rtol = 1.0e-4)
end

@testitem "Cly Conservation (ClOx Cycle)" setup = [StratSetup] tags = [:stratospheric] begin
    sys = ClOxCycle()
    compiled_sys = mtkcompile(sys)

    # SI units (m^-3)
    u0 = [
        compiled_sys.Cl => 1.0e10,
        compiled_sys.ClO => 1.0e13,
        compiled_sys.HCl => 1.0e15,
        compiled_sys.ClONO2 => 1.0e15,
    ]
    tspan = (0.0, 1.0)  # Very short time

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol = 1.0e-12, reltol = 1.0e-12)

    @test sol.retcode == ReturnCode.Success

    Cly_init = sol[compiled_sys.Cl][1] + sol[compiled_sys.ClO][1] +
        sol[compiled_sys.HCl][1] + sol[compiled_sys.ClONO2][1]
    Cly_final = sol[compiled_sys.Cl][end] + sol[compiled_sys.ClO][end] +
        sol[compiled_sys.HCl][end] + sol[compiled_sys.ClONO2][end]

    @test isapprox(Cly_init, Cly_final, rtol = 1.0e-6)
end

@testitem "Bry Conservation (BrOx Cycle)" setup = [StratSetup] tags = [:stratospheric] begin
    sys = BrOxCycle()
    compiled_sys = mtkcompile(sys)

    # SI units (m^-3)
    u0 = [
        compiled_sys.Br => 1.0e11,
        compiled_sys.BrO => 1.0e12,
        compiled_sys.HOBr => 1.0e12,
    ]
    tspan = (0.0, 1.0)  # Very short time

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol = 1.0e-12, reltol = 1.0e-12)

    @test sol.retcode == ReturnCode.Success

    Bry_init = sol[compiled_sys.Br][1] + sol[compiled_sys.BrO][1] +
        sol[compiled_sys.HOBr][1]
    Bry_final = sol[compiled_sys.Br][end] + sol[compiled_sys.BrO][end] +
        sol[compiled_sys.HOBr][end]

    @test isapprox(Bry_init, Bry_final, rtol = 1.0e-6)
end

@testitem "HOx Conservation (HOx Cycle)" setup = [StratSetup] tags = [:stratospheric] begin
    sys = HOxCycle()
    compiled_sys = mtkcompile(sys)

    # SI units (m^-3)
    u0 = [
        compiled_sys.OH => 1.0e12,
        compiled_sys.HO2 => 1.0e13,
    ]
    tspan = (0.0, 1.0)  # Very short time

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol = 1.0e-12, reltol = 1.0e-12)

    @test sol.retcode == ReturnCode.Success

    HOx_init = sol[compiled_sys.OH][1] + sol[compiled_sys.HO2][1]
    HOx_final = sol[compiled_sys.OH][end] + sol[compiled_sys.HO2][end]

    @test isapprox(HOx_init, HOx_final, rtol = 1.0e-6)
end

@testitem "NOx Conservation (NOx Cycle)" setup = [StratSetup] tags = [:stratospheric] begin
    sys = NOxCycle()
    compiled_sys = mtkcompile(sys)

    # SI units (m^-3)
    u0 = [
        compiled_sys.NO => 1.0e15,
        compiled_sys.NO2 => 1.0e15,
    ]
    tspan = (0.0, 1.0)  # Very short time

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, abstol = 1.0e-12, reltol = 1.0e-12)

    @test sol.retcode == ReturnCode.Success

    NOx_init = sol[compiled_sys.NO][1] + sol[compiled_sys.NO2][1]
    NOx_final = sol[compiled_sys.NO][end] + sol[compiled_sys.NO2][end]

    @test isapprox(NOx_init, NOx_final, rtol = 1.0e-4)
end

# =============================================================================
# Steady-State Tests
# =============================================================================

@testitem "Chapman Steady-State [O]/[O3] Ratio (Eq. 5.7)" setup = [StratSetup] tags = [:stratospheric] begin
    # At steady state, [O]/[O3] = j_O3 / (k2 * [O2] * [M]) = j_O3 / (k2 * 0.21 * M^2)
    # Eq. 5.7 from Seinfeld & Pandis (2006)

    sys = ChapmanMechanism()
    compiled_sys = mtkcompile(sys)

    j_O3_val = 4.0e-4  # s^-1
    T_val = 227.0  # K
    k2_si = k_O_O2_M_si(T_val)
    M_si = 3.1e23  # m^-3 (3.1e17 cm^-3 × 1e6)
    O2_val = 0.21

    tspan = (0.0, 3600.0 * 24 * 30)  # 30 days

    prob = ODEProblem(
        compiled_sys,
        [
            compiled_sys.O => 1.0e11, compiled_sys.O3 => 1.0e16,
            compiled_sys.j_O2 => 1.0e-11, compiled_sys.j_O3 => j_O3_val,
            compiled_sys.T => T_val,
            compiled_sys.M => M_si, compiled_sys.O2_mix => O2_val,
        ],
        tspan
    )
    sol = solve(prob, abstol = 1.0e-8, reltol = 1.0e-8)

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

@testitem "Chapman Steady-State Ozone (Eq. 5.13)" setup = [StratSetup] tags = [:stratospheric] begin
    # Equation 5.13: [O3]_ss = 0.21 * sqrt(k2 * j_O2 / (k4 * j_O3)) * [M]^(3/2)
    # Using SI units throughout

    T = 227.0
    M_si = 3.1e23  # m^-3
    k2_si = k_O_O2_M_si(T)  # m^6/s (SI)
    k4_si = k_O_O3_si(T)    # m^3/s (SI)
    j_O2 = 1.0e-11  # s^-1
    j_O3 = 4.0e-4   # s^-1

    O3_si = 0.21 * sqrt(k2_si * j_O2 / (k4_si * j_O3)) * M_si^1.5

    # The steady-state O3 should be on the order of 10^18 m^-3 at 30 km
    # (equivalent to ~10^12 molec/cm^3 × 1e6)
    @test O3_si > 1.0e16
    @test O3_si < 1.0e20
end

@testitem "Time to Steady State (Eq. 5.17)" setup = [StratSetup] tags = [:stratospheric] begin
    # Equation 5.17: tau_O3_ss = (1/4) * sqrt(k2*[M] / (k4 * j_O2 * j_O3))
    # Using SI units: k2 in m^6/s, M in m^-3, k4 in m^3/s

    # 30 km conditions
    T = 227.0
    M_si = 3.1e23  # m^-3
    k2 = k_O_O2_M_si(T)  # m^6/s
    k4 = k_O_O3_si(T)    # m^3/s
    j_O2 = 6.0e-11
    j_O3 = 1.2e-3

    tau = 0.25 * sqrt(k2 * M_si / (k4 * j_O2 * j_O3))
    tau_hours = tau / 3600.0

    # At 30 km, the reference value is ~160 hours
    @test tau_hours > 50   # order of magnitude check
    @test tau_hours < 500

    # 40 km conditions
    T40 = 251.0
    M40_si = 7.1e22  # m^-3 (7.1e16 cm^-3 × 1e6)
    k2_40 = k_O_O2_M_si(T40)
    k4_40 = k_O_O3_si(T40)
    j_O2_40 = 5.0e-10
    j_O3_40 = 1.9e-3

    tau_40 = 0.25 * sqrt(k2_40 * M40_si / (k4_40 * j_O2_40 * j_O3_40))
    tau_hours_40 = tau_40 / 3600.0

    # At 40 km, the reference value is ~12-40 hours
    @test tau_hours_40 > 1
    @test tau_hours_40 < 100

    # Timescale should decrease with altitude (faster at higher altitudes)
    @test tau_hours_40 < tau_hours
end

@testitem "[O]/[O3] Ratio Reference Values (Page 145)" setup = [StratSetup] tags = [:stratospheric] begin
    # Reference: Page 145, Seinfeld & Pandis (2006)
    # At z=30km: [O]/[O3] = 3.0 × 10^-5
    # At z=40km: [O]/[O3] = 9.4 × 10^-4

    # Using Eq. 5.7: [O]/[O3] = j_O3 / (k2 * [O2] * [M])
    # All in SI: k2 in m^6/s, M in m^-3

    # 30 km
    T30 = 227.0
    M30_si = 3.1e23  # m^-3
    k2_30 = k_O_O2_M_si(T30)
    j_O3_30 = 1.2e-3  # Total O3 photolysis rate at 30 km

    ratio_30 = j_O3_30 / (k2_30 * 0.21 * M30_si * M30_si)

    # Should be on the order of 10^-5
    @test ratio_30 > 1.0e-6
    @test ratio_30 < 1.0e-3

    # 40 km
    T40 = 251.0
    M40_si = 7.1e22  # m^-3
    k2_40 = k_O_O2_M_si(T40)
    j_O3_40 = 1.9e-3

    ratio_40 = j_O3_40 / (k2_40 * 0.21 * M40_si * M40_si)

    # Should be larger than at 30 km (lower M means higher ratio)
    @test ratio_40 > ratio_30
    @test ratio_40 > 1.0e-5
    @test ratio_40 < 1.0e-2
end

# =============================================================================
# Limiting Behavior Tests
# =============================================================================

@testitem "Zero Photolysis Rates" setup = [StratSetup] tags = [:stratospheric] begin
    sys = ChapmanMechanism()
    compiled_sys = mtkcompile(sys)

    # SI units (m^-3)
    tspan = (0.0, 3600.0)

    prob = ODEProblem(
        compiled_sys,
        [
            compiled_sys.O => 1.0e13, compiled_sys.O3 => 3.0e18,
            compiled_sys.j_O2 => 0.0, compiled_sys.j_O3 => 0.0,
        ],
        tspan
    )
    sol = solve(prob, abstol = 1.0e-10, reltol = 1.0e-10)

    @test sol.retcode == ReturnCode.Success

    # With no photolysis, O should decrease (consumed by O+O2+M and O+O3)
    @test sol[compiled_sys.O][end] < sol[compiled_sys.O][1]
end

# =============================================================================
# Comprehensive System Integration
# =============================================================================

@testitem "Comprehensive System Integration" setup = [StratSetup] tags = [:stratospheric] begin
    sys = StratosphericOzoneSystem()
    compiled_sys = mtkcompile(sys)

    prob = ODEProblem(compiled_sys, [], (0.0, 3600.0))  # 1 hour with defaults
    sol = solve(prob, Rodas5P(), abstol = 1.0e-8, reltol = 1.0e-8)

    @test sol.retcode == ReturnCode.Success

    # All concentrations should remain positive
    @test all(sol[compiled_sys.O] .>= 0)
    @test all(sol[compiled_sys.O3] .> 0)
    @test all(sol[compiled_sys.NO] .>= 0)
    @test all(sol[compiled_sys.NO2] .>= 0)
    @test all(sol[compiled_sys.N2O] .>= 0)
end

# =============================================================================
# N2O Chemistry Tests
# =============================================================================

@testitem "N2O Destruction in Comprehensive System" setup = [StratSetup] tags = [:stratospheric] begin
    # N2O is destroyed by photolysis and reaction with O(1D)
    # It should decrease over time when no replenishment source exists
    sys = StratosphericOzoneSystem()
    compiled_sys = mtkcompile(sys)

    prob = ODEProblem(compiled_sys, [], (0.0, 86400.0))  # 1 day with defaults
    sol = solve(prob, Rodas5P(), abstol = 1.0e-8, reltol = 1.0e-8)

    @test sol.retcode == ReturnCode.Success

    # N2O should decrease over time (no source, only sinks)
    @test sol[compiled_sys.N2O][end] <= sol[compiled_sys.N2O][1]
end

@testitem "N2O as NOx Source (Section 5.3.1)" setup = [StratSetup] tags = [:stratospheric] begin
    # The primary source of stratospheric NOx is N2O + O(1D) → 2NO
    # Verify NO equation includes this source term by checking that
    # starting with zero NO/NO2 but nonzero N2O and O1D, NO is produced

    sys = StratosphericOzoneSystem()
    compiled_sys = mtkcompile(sys)

    # Start with NO=0, NO2=0 but nonzero N2O and O1D
    prob = ODEProblem(
        compiled_sys,
        [
            compiled_sys.NO => 0.0, compiled_sys.NO2 => 0.0,
            compiled_sys.N2O => 9.3e16, compiled_sys.O1D => 5.0e7,
        ],
        (0.0, 100.0)
    )
    sol = solve(prob, Rodas5P(), abstol = 1.0e-8, reltol = 1.0e-8)

    @test sol.retcode == ReturnCode.Success

    # NO should be produced from N2O + O(1D) → 2NO
    @test sol[compiled_sys.NO][end] > 0
end

# =============================================================================
# Rate Coefficient Reference Value Tests
# =============================================================================

@testitem "Rate Coefficients at Reference Temperatures" setup = [StratSetup] tags = [:stratospheric] begin
    # Verify rate coefficients match textbook values at specific temperatures
    # All CGS values from Seinfeld & Pandis Chapter 5

    # k2 at 30 km (T=227K): 1.15e-33 cm^6/molec^2/s (Page 144)
    k2_cgs_227 = 6.0e-34 * (227.0 / 300.0)^(-2.4)
    @test isapprox(k2_cgs_227, 1.15e-33, rtol = 0.05)

    # k2 at 40 km (T=251K): 9.1e-34 cm^6/molec^2/s (Page 144)
    k2_cgs_251 = 6.0e-34 * (251.0 / 300.0)^(-2.4)
    @test isapprox(k2_cgs_251, 9.1e-34, rtol = 0.05)

    # k4 at 30 km (T=227K): 9.2e-16 cm^3/molec/s (Page 145)
    k4_cgs_227 = 8.0e-12 * exp(-2060.0 / 227.0)
    @test isapprox(k4_cgs_227, 9.2e-16, rtol = 0.1)

    # k4 at 40 km (T=251K): 2.2e-15 cm^3/molec/s (Page 145)
    k4_cgs_251 = 8.0e-12 * exp(-2060.0 / 251.0)
    @test isapprox(k4_cgs_251, 2.2e-15, rtol = 0.1)

    # Verify SI conversions
    k2_si_227 = k_O_O2_M_si(227.0)
    @test isapprox(k2_si_227, k2_cgs_227 * CGS_TO_SI_K2, rtol = 1.0e-10)

    k4_si_227 = k_O_O3_si(227.0)
    @test isapprox(k4_si_227, k4_cgs_227 * CGS_TO_SI_K, rtol = 1.0e-10)
end

@testitem "Characteristic Time tau_2 (Eq. 5.4)" setup = [StratSetup] tags = [:stratospheric] begin
    # tau_2 = 1 / (0.21 * k2 * [M]^2) — Eq. 5.4, Page 143
    # At 30 km (T=227K, M=3.1e17 cm^-3): tau_2 ≈ 0.04 s
    # At 40 km (T=251K, M=7.1e16 cm^-3): tau_2 ≈ 1.04 s

    M_30_cgs = 3.1e17  # molec/cm^3
    k2_30_cgs = 6.0e-34 * (227.0 / 300.0)^(-2.4)
    tau2_30 = 1.0 / (0.21 * k2_30_cgs * M_30_cgs^2)
    @test isapprox(tau2_30, 0.04, rtol = 0.2)

    M_40_cgs = 7.1e16
    k2_40_cgs = 6.0e-34 * (251.0 / 300.0)^(-2.4)
    tau2_40 = 1.0 / (0.21 * k2_40_cgs * M_40_cgs^2)
    @test isapprox(tau2_40, 1.04, rtol = 0.2)
end

@testitem "O(1D) Steady-State Concentration (Eq. 5.19)" setup = [StratSetup] tags = [:stratospheric] begin
    # Eq. 5.19: [O(1D)] = j_{O3→O(1D)} * [O3] / (k4 * [M])
    # At 30 km, theta=45°: j_N2O ≈ 5e-8 s^-1, j_{O3→O(1D)} ≈ 15e-5 s^-1
    # k4 = 3.2e-11 cm^3/molec/s (weighted quenching rate)
    # [M] = 3.1e17 cm^-3, [O3] = 3e12 cm^-3
    # [O(1D)] ≈ 45 molec/cm^3 (Page 153)

    T = 227.0
    M_cgs = 3.1e17  # molec/cm^3
    O3_cgs = 3.0e12  # molec/cm^3
    j_O3_O1D = 15.0e-5  # s^-1

    # Weighted quenching rate: k = 0.21 * k_O2 + 0.79 * k_N2
    k_O1D_M_cgs = 0.21 * 3.2e-11 * exp(70.0 / T) + 0.79 * 1.8e-11 * exp(110.0 / T)

    O1D_ss = j_O3_O1D * O3_cgs / (k_O1D_M_cgs * M_cgs)

    # Reference: ~45 molec/cm^3 at 30 km (Page 153)
    @test O1D_ss > 10
    @test O1D_ss < 200
end

@testitem "Ox Lifetime (Eq. 5.8)" setup = [StratSetup] tags = [:stratospheric] begin
    # tau_Ox ≈ 0.21 * k2 * [M]^2 / (k4 * j_O3 * [O3]) — Eq. 5.8, Page 145
    # At 30 km: tau_Ox ≈ 1.2e7 s (~140 days)
    # At 40 km: tau_Ox ≈ 1e6 s (~12 days)

    # 30 km
    M_30_cgs = 3.1e17
    O3_30_cgs = 3.0e12
    k2_30 = 6.0e-34 * (227.0 / 300.0)^(-2.4)
    k4_30 = 8.0e-12 * exp(-2060.0 / 227.0)
    j_O3_30 = 1.2e-3

    tau_Ox_30 = 0.21 * k2_30 * M_30_cgs^2 / (k4_30 * j_O3_30 * O3_30_cgs)
    tau_Ox_30_days = tau_Ox_30 / 86400.0

    # Reference: ~140 days at 30 km (Page 145)
    @test tau_Ox_30_days > 50
    @test tau_Ox_30_days < 500

    # 40 km
    M_40_cgs = 7.1e16
    O3_40_cgs = 0.5e12
    k2_40 = 6.0e-34 * (251.0 / 300.0)^(-2.4)
    k4_40 = 8.0e-12 * exp(-2060.0 / 251.0)
    j_O3_40 = 1.9e-3

    tau_Ox_40 = 0.21 * k2_40 * M_40_cgs^2 / (k4_40 * j_O3_40 * O3_40_cgs)
    tau_Ox_40_days = tau_Ox_40 / 86400.0

    # Reference: ~12 days at 40 km (Page 145)
    @test tau_Ox_40_days > 3
    @test tau_Ox_40_days < 50

    # Lifetime should decrease with altitude
    @test tau_Ox_40_days < tau_Ox_30_days
end

@testitem "Time to Steady State Table Values (Eq. 5.17)" setup = [StratSetup] tags = [:stratospheric] begin
    # Table on Page 147 gives tau_O3^ss at several altitudes
    # All in CGS units for comparison

    # z=20km: T=217K, k4=6e-16, j_O2=1e-11, j_O3=0.7e-3, tau≈1400h
    # z=25km: T=222K, k4=7.5e-16, j_O2=2e-11, j_O3=0.7e-3, tau≈600h
    # z=30km: T=227K, k4=9.2e-16, j_O2=6e-11, j_O3=1.2e-3, tau≈160h
    # z=40km: T=251K, k4=2.2e-15, j_O2=5e-10, j_O3=1.9e-3, tau≈12h
    # z=45km: T=265K, k4=3.4e-15, j_O2=8e-10, j_O3=6e-3, tau≈3h

    altitudes = [20, 25, 30, 40, 45]
    temps = [217.0, 222.0, 227.0, 251.0, 265.0]
    M_cgs = [1.4e18, 6.4e17, 3.1e17, 7.1e16, 3.6e16]
    j_O2_vals = [1.0e-11, 2.0e-11, 6.0e-11, 5.0e-10, 8.0e-10]
    j_O3_vals = [0.7e-3, 0.7e-3, 1.2e-3, 1.9e-3, 6.0e-3]
    tau_ref_hours = [1400.0, 600.0, 160.0, 12.0, 3.0]

    for (i, z) in enumerate(altitudes)
        T = temps[i]
        M = M_cgs[i]
        k2 = 6.0e-34 * (T / 300.0)^(-2.4)
        k4 = 8.0e-12 * exp(-2060.0 / T)
        j_O2 = j_O2_vals[i]
        j_O3 = j_O3_vals[i]

        tau = 0.25 * sqrt(k2 * M / (k4 * j_O2 * j_O3))
        tau_h = tau / 3600.0

        # Allow factor of 2 tolerance due to approximate photolysis rate values
        @test isapprox(tau_h, tau_ref_hours[i], rtol = 0.5)
    end
end

@testitem "NOx Cycle Rate Ratio (Page 155)" setup = [StratSetup] tags = [:stratospheric] begin
    # The ratio k_{NO2+O}[NO2] / (k_{O+O3}[O3]) at 35 km ≈ 4.5
    # This demonstrates NOx cycle is ~5x more effective than Chapman (Page 155)

    T = 237.0  # 35 km
    k_NO2_O = 5.6e-12 * exp(180.0 / T)   # cm^3/molec/s
    k_O_O3 = 8.0e-12 * exp(-2060.0 / T)  # cm^3/molec/s

    rate_coeff_ratio = k_NO2_O / k_O_O3
    # Reference: ~9000 (Page 155)
    @test rate_coeff_ratio > 5000
    @test rate_coeff_ratio < 15000

    # With [NO2] ≈ 1e9 cm^-3 and [O3] ≈ 2e12 cm^-3 at 35 km
    NO2_conc = 1.0e9   # cm^-3
    O3_conc = 2.0e12    # cm^-3
    full_ratio = rate_coeff_ratio * NO2_conc / O3_conc
    # Reference: ≈ 4.5 (Page 155)
    @test isapprox(full_ratio, 4.5, rtol = 0.5)
end

@testitem "ClOx Lifetime Estimates (Page 163)" setup = [StratSetup] tags = [:stratospheric] begin
    # At 40 km (T=251K), [O3] ≈ 0.5e12 cm^-3, [O]/[O3] ≈ 9.4e-4
    # tau_Cl = 1/(k1*[O3]) ≈ 0.2 s
    # tau_ClO = 1/(k2*[O]) ≈ 53 s

    T = 251.0
    O3_cgs = 0.5e12  # cm^-3
    O_O3_ratio = 9.4e-4
    O_cgs = O_O3_ratio * O3_cgs  # ≈ 4.7e8 cm^-3

    k_Cl_O3 = 2.3e-11 * exp(-200.0 / T)
    k_ClO_O = 3.0e-11 * exp(70.0 / T)

    tau_Cl = 1.0 / (k_Cl_O3 * O3_cgs)
    tau_ClO = 1.0 / (k_ClO_O * O_cgs)

    # Reference: tau_Cl ≈ 0.2 s, tau_ClO ≈ 53 s (Page 163)
    @test isapprox(tau_Cl, 0.2, rtol = 0.3)
    @test isapprox(tau_ClO, 53.0, rtol = 0.3)
end

@testitem "HO2/OH Ratio (Eq. 5.28)" setup = [StratSetup] tags = [:stratospheric] begin
    # Equation 5.28: [HO2]/[OH] = k_{OH+O3}[O3] / (k_{HO2+NO}[NO])
    # At 30 km (T=227K): k_{OH+O3} = 1.7e-12 exp(-940/T), k_{HO2+NO} = 3.5e-12 exp(250/T)
    # With [O3] ≈ 2e12 cm^-3 and NO mixing ratio of 3 ppb at 30 km (Page 161)

    T = 227.0
    k_OH_O3 = 1.7e-12 * exp(-940.0 / T)  # cm^3/molec/s
    k_HO2_NO = 3.5e-12 * exp(250.0 / T)  # cm^3/molec/s

    O3_cgs = 2.0e12  # molec/cm^3
    M_cgs = 3.1e17  # molec/cm^3
    NO_cgs = 3.0e-9 * M_cgs  # 3 ppb → molec/cm^3 ≈ 9.3e8

    ratio = k_OH_O3 * O3_cgs / (k_HO2_NO * NO_cgs)

    # Reference: ~4.4 at 30 km (Page 161)
    @test ratio > 1.0
    @test ratio < 20.0
end

@testitem "Cl/ClO Steady-State Ratio (Eq. 5.30)" setup = [StratSetup] tags = [:stratospheric] begin
    # Equation 5.30: [Cl]/[ClO] = (k_{ClO+O}[O] + k_{ClO+NO}[NO]) / (k_{Cl+O3}[O3])
    # At 40 km (T=251K): [Cl]/[ClO] ≈ 0.008 (Page 164)

    T = 251.0
    O3_cgs = 0.5e12  # cm^-3
    O_O3_ratio = 9.4e-4
    O_cgs = O_O3_ratio * O3_cgs  # ≈ 4.7e8 cm^-3
    NO_cgs = 1.0e9  # cm^-3 (Page 164)

    k_Cl_O3 = 2.3e-11 * exp(-200.0 / T)  # cm^3/molec/s
    k_ClO_O = 3.0e-11 * exp(70.0 / T)    # cm^3/molec/s
    k_ClO_NO = 6.4e-12 * exp(290.0 / T)  # cm^3/molec/s

    Cl_ClO_ratio = (k_ClO_O * O_cgs + k_ClO_NO * NO_cgs) / (k_Cl_O3 * O3_cgs)

    # Reference: ≈ 0.008 at 40 km (Page 164)
    @test Cl_ClO_ratio > 0.001
    @test Cl_ClO_ratio < 0.1
end

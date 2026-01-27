"""
    StratosphericChemistry Tests

Comprehensive unit tests for the StratosphericChemistry module, verifying the
ModelingToolkit.jl implementation against Chapter 5 of Seinfeld & Pandis (2006)
"Atmospheric Chemistry and Physics: From Air Pollution to Climate Change", 2nd Edition.

Reference: Seinfeld, J.H. and Pandis, S.N. (2006), Chapter 5, pp. 138-203.
"""

# =============================================================================
# Test Suite 1: Rate Coefficient Functions (Tables B.1 and B.2)
# =============================================================================

@testitem "Rate Coefficient k_O_O2_M" begin
    using GasChem

    # Reference: Table B.2, k = 6.0 x 10^-34 (T/300)^-2.4 cm^6 molecule^-2 s^-1

    # Test at 300 K (reference temperature)
    T = 300.0
    k_expected = 6.0e-34
    k_actual = GasChem.k_O_O2_M(T)
    @test isapprox(k_actual, k_expected, rtol=1e-10)

    # Test at 227 K (typical 30 km altitude, Page 140)
    T = 227.0
    k_expected = 6.0e-34 * (227.0 / 300.0)^(-2.4)
    k_actual = GasChem.k_O_O2_M(T)
    @test isapprox(k_actual, k_expected, rtol=1e-10)

    # Verify rate increases at lower temperatures (negative exponent)
    @test GasChem.k_O_O2_M(200.0) > GasChem.k_O_O2_M(300.0)
end

@testitem "Rate Coefficient k_O_O3" begin
    using GasChem

    # Reference: Table B.1, k = 8.0 x 10^-12 exp(-2060/T) cm^3 molecule^-1 s^-1

    # Test at 227 K
    T = 227.0
    k_expected = 8.0e-12 * exp(-2060.0 / T)
    k_actual = GasChem.k_O_O3(T)
    @test isapprox(k_actual, k_expected, rtol=1e-10)

    # Test at 300 K
    T = 300.0
    k_expected = 8.0e-12 * exp(-2060.0 / T)
    k_actual = GasChem.k_O_O3(T)
    @test isapprox(k_actual, k_expected, rtol=1e-10)

    # Verify rate increases with temperature (Arrhenius behavior with positive Ea)
    @test GasChem.k_O_O3(300.0) > GasChem.k_O_O3(200.0)
end

@testitem "Chapman Mechanism Model Structure" begin
    using GasChem
    using ModelingToolkit

    sys = ChapmanMechanism()

    # Verify the system has the correct number of equations
    # 2 ODEs (O, O3) + 1 algebraic (Ox)
    @test length(equations(sys)) == 3

    # Verify state variables
    states = unknowns(sys)
    state_names = [string(s) for s in states]
    @test any(occursin("O(t)", n) for n in state_names)
    @test any(occursin("O3(t)", n) for n in state_names)
    @test any(occursin("Ox(t)", n) for n in state_names)
end

@testitem "Chapman Mechanism Integration" begin
    using GasChem
    using ModelingToolkit
    using OrdinaryDiffEq

    # Test that the Chapman mechanism can be integrated without errors

    sys = ChapmanMechanism()
    compiled_sys = structural_simplify(sys)

    # Initial conditions
    u0 = [
        compiled_sys.O => 1e7,
        compiled_sys.O3 => 3e12
    ]

    # Time span (1 day in seconds)
    tspan = (0.0, 86400.0)

    # Create and solve problem
    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-8)

    # Verify solution succeeded
    @test sol.retcode == ReturnCode.Success

    # Verify concentrations remain positive
    @test all(sol[compiled_sys.O] .> 0)
    @test all(sol[compiled_sys.O3] .> 0)
end

@testitem "NOx Cycle Model Structure" begin
    using GasChem
    using ModelingToolkit

    sys = NOxCycle()

    # 2 ODEs (NO, NO2) + 1 algebraic (NOx)
    @test length(equations(sys)) == 3

    # Verify state variables
    states = unknowns(sys)
    state_names = [string(s) for s in states]
    @test any(occursin("NO(t)", n) for n in state_names)
    @test any(occursin("NO2(t)", n) for n in state_names)
    @test any(occursin("NOx(t)", n) for n in state_names)
end

@testitem "HOx Cycle Model Structure" begin
    using GasChem
    using ModelingToolkit

    sys = HOxCycle()

    # 2 ODEs (OH, HO2) + 1 algebraic (HOx)
    @test length(equations(sys)) == 3

    # Verify state variables
    states = unknowns(sys)
    state_names = [string(s) for s in states]
    @test any(occursin("OH(t)", n) for n in state_names)
    @test any(occursin("HO2(t)", n) for n in state_names)
    @test any(occursin("HOx(t)", n) for n in state_names)
end

@testitem "ClOx Cycle Model Structure" begin
    using GasChem
    using ModelingToolkit

    sys = ClOxCycle()

    # 4 ODEs (Cl, ClO, HCl, ClONO2) + 2 algebraic (ClOx, Cly)
    @test length(equations(sys)) == 6

    # Verify state variables
    states = unknowns(sys)
    state_names = [string(s) for s in states]
    @test any(occursin("Cl(t)", n) for n in state_names)
    @test any(occursin("ClO(t)", n) for n in state_names)
    @test any(occursin("ClOx(t)", n) for n in state_names)
end

@testitem "BrOx Cycle Model Structure" begin
    using GasChem
    using ModelingToolkit

    sys = BrOxCycle()

    # 3 ODEs (Br, BrO, HOBr) + 2 algebraic (BrOx, Bry)
    @test length(equations(sys)) == 5

    # Verify state variables
    states = unknowns(sys)
    state_names = [string(s) for s in states]
    @test any(occursin("Br(t)", n) for n in state_names)
    @test any(occursin("BrO(t)", n) for n in state_names)
    @test any(occursin("BrOx(t)", n) for n in state_names)
end

@testitem "Comprehensive Stratospheric Ozone System" begin
    using GasChem
    using ModelingToolkit

    sys = StratosphericOzoneSystem()

    # Verify the system can be created without errors
    @test sys isa ODESystem

    # The full system should have many equations (19 currently)
    @test length(equations(sys)) >= 19

    # Verify key species are present
    states = unknowns(sys)
    state_names = [string(s) for s in states]

    # Odd oxygen
    @test any(occursin("O(t)", n) for n in state_names)
    @test any(occursin("O3(t)", n) for n in state_names)

    # NOx
    @test any(occursin("NO(t)", n) for n in state_names)
    @test any(occursin("NO2(t)", n) for n in state_names)

    # HOx
    @test any(occursin("OH(t)", n) for n in state_names)
    @test any(occursin("HO2(t)", n) for n in state_names)

    # ClOx
    @test any(occursin("Cl(t)", n) for n in state_names)
    @test any(occursin("ClO(t)", n) for n in state_names)

    # BrOx
    @test any(occursin("Br(t)", n) for n in state_names)
    @test any(occursin("BrO(t)", n) for n in state_names)
end

@testitem "Stratospheric Rate Coefficients Dictionary" begin
    using GasChem

    T = 227.0
    M = 3.1e17

    coeffs = stratospheric_rate_coefficients(T, M)

    # Verify all expected keys are present
    expected_keys = [:k_O_O2_M, :k_O_O3, :k_O1D_M, :k_NO_O3, :k_NO2_O,
                    :k_OH_O3, :k_HO2_O3, :k_HO2_O, :k_HO2_NO,
                    :k_Cl_O3, :k_ClO_O, :k_ClO_NO, :k_Cl_CH4, :k_OH_HCl,
                    :k_Br_O3, :k_BrO_ClO_BrCl, :k_BrO_ClO_ClOO, :k_ClO_ClO_M]

    for key in expected_keys
        @test haskey(coeffs, key)
        @test coeffs[key] > 0  # All rate coefficients should be positive
    end

    # Verify some specific values
    @test isapprox(coeffs[:k_O_O2_M], GasChem.k_O_O2_M(T), rtol=1e-10)
    @test isapprox(coeffs[:k_O_O3], GasChem.k_O_O3(T), rtol=1e-10)
end

@testitem "Positivity Preservation" begin
    using GasChem
    using ModelingToolkit
    using OrdinaryDiffEq

    # Chemical concentrations must remain non-negative
    # Use ChapmanMechanism which is self-contained (no external inputs)

    sys = ChapmanMechanism()
    compiled_sys = structural_simplify(sys)

    # Initial conditions
    u0 = [
        compiled_sys.O => 1e7,
        compiled_sys.O3 => 3e12
    ]

    # Short integration
    tspan = (0.0, 3600.0)  # 1 hour

    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, Tsit5(), abstol=1e-10, reltol=1e-10)

    @test sol.retcode == ReturnCode.Success
    @test all(sol[compiled_sys.O] .>= 0)
    @test all(sol[compiled_sys.O3] .>= 0)
end

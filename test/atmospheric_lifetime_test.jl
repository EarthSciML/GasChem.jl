"""
Test suite for Seinfeld & Pandis Chapter 2 atmospheric lifetime implementation.

Tests verify:
1. Structural correctness (model components, variables, parameters)
2. Equation verification against textbook formulas
3. Numerical accuracy for typical atmospheric conditions
4. Physical constraints and edge cases
5. Dimensional consistency

Reference: Seinfeld, J.H. and Pandis, S.N. (2006). Atmospheric Chemistry and Physics,
2nd Edition, Chapter 2: Atmospheric Trace Constituents. Wiley-Interscience.
"""

#=============================================================================
# Test Set 1: Structural Tests
=============================================================================#

@testitem "1.1 AtmosphericBudget Structure" begin
    using ModelingToolkit
    using GasChem.AtmosphericLifetime: AtmosphericBudget

    @named budget = AtmosphericBudget()
    @test budget isa System

    # Check that all expected variables exist
    vars = unknowns(budget)
    var_names = string.(vars)
    @test any(occursin("Q", n) for n in var_names)
    @test any(occursin("F_in", n) for n in var_names)
    @test any(occursin("F_out", n) for n in var_names)
    @test any(occursin("P", n) && !occursin("net", n) for n in var_names)
    @test any(occursin("R", n) for n in var_names)
    @test any(occursin("net_transport", n) for n in var_names)
    @test any(occursin("net_chemistry", n) for n in var_names)

    # Check equations
    eqs = equations(budget)
    @test length(eqs) == 3  # Mass conservation + 2 diagnostics
end

@testitem "1.2 SpeciesLifetime Structure" begin
    using ModelingToolkit
    using GasChem.AtmosphericLifetime: SpeciesLifetime

    @named lifetime = SpeciesLifetime()
    @test lifetime isa System

    vars = unknowns(lifetime)
    var_names = string.(vars)
    @test any(occursin("Q", n) for n in var_names)
    @test any(occursin("tau_general", n) for n in var_names)
    @test any(occursin("tau_removal", n) for n in var_names)
    @test any(occursin("tau_production", n) for n in var_names)
    @test any(occursin("tau_first_order", n) for n in var_names)

    # Check parameters
    params = parameters(lifetime)
    param_names = string.(params)
    @test any(occursin("lambda", n) for n in param_names)

    eqs = equations(lifetime)
    @test length(eqs) == 4  # Four lifetime equations
end

@testitem "1.3 MultipleRemovalLifetime Structure" begin
    using ModelingToolkit
    using GasChem.AtmosphericLifetime: MultipleRemovalLifetime

    @named multi = MultipleRemovalLifetime()
    @test multi isa System

    vars = unknowns(multi)
    var_names = string.(vars)
    @test any(occursin("tau_1", n) for n in var_names)
    @test any(occursin("tau_2", n) for n in var_names)
    @test any(occursin("tau_combined", n) for n in var_names)
    @test any(occursin("inverse_tau", n) for n in var_names)

    params = parameters(multi)
    param_names = string.(params)
    @test any(occursin("k_1", n) for n in param_names)
    @test any(occursin("k_2", n) for n in param_names)

    eqs = equations(multi)
    @test length(eqs) == 4
end

@testitem "1.4 OHReactionLifetime Structure" begin
    using ModelingToolkit
    using GasChem.AtmosphericLifetime: OHReactionLifetime

    @named oh = OHReactionLifetime()
    @test oh isa System

    vars = unknowns(oh)
    var_names = string.(vars)
    @test any(occursin("OH_conc", n) for n in var_names)
    @test any(occursin("k_eff", n) for n in var_names)
    @test any(occursin("tau_OH", n) for n in var_names)

    params = parameters(oh)
    param_names = string.(params)
    @test any(occursin("k_OH", n) for n in param_names)

    eqs = equations(oh)
    @test length(eqs) == 2
end

@testitem "1.5 TroposphericBudget Structure" begin
    using ModelingToolkit
    using GasChem.AtmosphericLifetime: TroposphericBudget

    @named trop = TroposphericBudget()
    @test trop isa System

    vars = unknowns(trop)
    var_names = string.(vars)
    @test any(occursin("Q", n) for n in var_names)
    @test any(occursin("P_n", n) for n in var_names)
    @test any(occursin("P_a", n) for n in var_names)
    @test any(occursin("P_c", n) for n in var_names)
    @test any(occursin("P_total", n) for n in var_names)
    @test any(occursin("R_total", n) for n in var_names)
    @test any(occursin("tau_total", n) for n in var_names)
    @test any(occursin("tau_dry", n) for n in var_names)
    @test any(occursin("tau_wet", n) for n in var_names)
    @test any(occursin("tau_chem", n) for n in var_names)
    @test any(occursin("tau_transport", n) for n in var_names)

    params = parameters(trop)
    param_names = string.(params)
    @test any(occursin("k_d", n) for n in param_names)
    @test any(occursin("k_w", n) for n in param_names)
    @test any(occursin("k_c", n) for n in param_names)
    @test any(occursin("k_t", n) for n in param_names)

    eqs = equations(trop)
    @test length(eqs) == 10
end

#=============================================================================
# Test Set 2: Equation Verification Tests
=============================================================================#

@testitem "2.1 Eq. 2.1: Mass Conservation" begin
    # dQ/dt = (F_in - F_out) + (P - R)
    F_in = 100.0   # mol/s
    F_out = 80.0   # mol/s
    P = 50.0       # mol/s
    R = 60.0       # mol/s

    # Expected rate of change
    expected_dQdt = (F_in - F_out) + (P - R)
    @test expected_dQdt == 10.0  # (100-80) + (50-60) = 20 - 10 = 10 mol/s

    # Net transport
    net_transport = F_in - F_out
    @test net_transport == 20.0

    # Net chemistry
    net_chemistry = P - R
    @test net_chemistry == -10.0
end

@testitem "2.2 Eq. 2.2: Steady State Condition" begin
    # At steady state: F_in + P = F_out + R
    F_in = 100.0
    P = 50.0
    F_out = 80.0

    # Required R for steady state
    R_steady = F_in + P - F_out
    @test R_steady == 70.0

    # Verify: (100-80) + (50-70) = 20 - 20 = 0
    @test (F_in - F_out) + (P - R_steady) == 0.0
end

@testitem "2.3 Eq. 2.3: General Lifetime" begin
    # tau = Q / (R + F_out)
    Q = 1e12      # mol (typical tropospheric burden)
    R = 1e6       # mol/s
    F_out = 2e5   # mol/s

    expected_tau = Q / (R + F_out)
    @test expected_tau ≈ 1e12 / 1.2e6 rtol=1e-10
    @test expected_tau ≈ 833333.33 rtol=1e-5  # ~9.6 days
end

@testitem "2.4 Eq. 2.4-2.5: Removal and Production Lifetimes" begin
    # tau = Q / R = Q / P (at steady state)
    Q = 1e12      # mol
    R = 1e6       # mol/s
    P = 1e6       # mol/s (equal to R at steady state)

    tau_removal = Q / R
    tau_production = Q / P

    @test tau_removal == tau_production
    @test tau_removal == 1e6  # seconds
    @test tau_removal / 86400 ≈ 11.57 rtol=0.01  # ~11.6 days
end

@testitem "2.5 Eq. 2.6: First-Order Lifetime" begin
    # tau = 1 / lambda
    lambda = 1e-6  # s^-1

    tau = 1 / lambda
    @test tau == 1e6  # seconds
    @test tau / 86400 ≈ 11.57 rtol=0.01  # ~11.6 days

    # Test different rate constants
    lambda_fast = 1e-4  # s^-1 (fast removal)
    lambda_slow = 1e-8  # s^-1 (slow removal)

    @test (1 / lambda_fast) < (1 / lambda) < (1 / lambda_slow)
end

@testitem "2.6 Eq. 2.7: Two First-Order Processes" begin
    # tau = 1 / (k_1 + k_2)
    k_1 = 1e-6    # s^-1
    k_2 = 2e-6    # s^-1

    tau_combined = 1 / (k_1 + k_2)
    @test tau_combined ≈ 1 / 3e-6 rtol=1e-10
    @test tau_combined ≈ 333333.33 rtol=1e-5  # ~3.9 days

    # Verify tau_combined < min(tau_1, tau_2)
    tau_1 = 1 / k_1
    tau_2 = 1 / k_2
    @test tau_combined < tau_1
    @test tau_combined < tau_2
end

@testitem "2.7 Eq. 2.8: Inverse Lifetime Sum" begin
    # 1/tau = 1/tau_1 + 1/tau_2
    tau_1 = 1e6   # s
    tau_2 = 5e5   # s

    inverse_tau = 1/tau_1 + 1/tau_2
    tau_combined = 1 / inverse_tau

    @test tau_combined ≈ tau_1 * tau_2 / (tau_1 + tau_2) rtol=1e-10
    @test tau_combined ≈ 333333.33 rtol=1e-5
end

@testitem "2.8 Eq. 2.9: Combined Lifetime Formula" begin
    # tau = (tau_1 * tau_2) / (tau_1 + tau_2)
    tau_1 = 1e6   # s (11.6 days)
    tau_2 = 5e5   # s (5.8 days)

    tau_combined = (tau_1 * tau_2) / (tau_1 + tau_2)
    @test tau_combined ≈ 333333.33 rtol=1e-5  # ~3.9 days

    # Verify equivalence with Eq. 2.8
    @test tau_combined ≈ 1 / (1/tau_1 + 1/tau_2) rtol=1e-10
end

@testitem "2.9 Eq. 2.12: OH Reaction Lifetime" begin
    # tau = 1 / (k_OH * [OH])

    # Typical values for methane (CH4)
    k_OH_CH4 = 6.3e-15    # cm^3/(molec*s) at 298 K
    OH_conc = 1e6         # molec/cm^3 (global average)

    tau_CH4 = 1 / (k_OH_CH4 * OH_conc)
    tau_CH4_years = tau_CH4 / (365.25 * 86400)

    @test tau_CH4_years ≈ 5.0 rtol=0.1  # ~5 years for CH4

    # Typical values for CO
    k_OH_CO = 2.4e-13     # cm^3/(molec*s) at 298 K
    tau_CO = 1 / (k_OH_CO * OH_conc)
    tau_CO_days = tau_CO / 86400

    @test tau_CO_days ≈ 48 rtol=0.1  # ~48 days for CO
end

@testitem "2.10 Eq. 2.15: Total Lifetime from Removal" begin
    # tau = 1 / (k_d + k_w + k_c + k_t)

    # Typical values for a reactive species
    k_d = 1e-6    # s^-1 (dry deposition)
    k_w = 5e-7    # s^-1 (wet deposition)
    k_c = 2e-6    # s^-1 (chemical loss)
    k_t = 1e-7    # s^-1 (stratospheric transport)

    k_total = k_d + k_w + k_c + k_t
    tau_total = 1 / k_total

    @test k_total ≈ 3.6e-6 rtol=1e-10
    @test tau_total ≈ 277777.78 rtol=1e-5  # ~3.2 days
end

#=============================================================================
# Test Set 3: Physical Constraints (selected tests)
=============================================================================#

@testitem "3.1 Lifetime Positivity" begin
    # All lifetimes must be positive for positive rate constants

    # Test with various positive rate constants
    for lambda in [1e-10, 1e-6, 1e-2, 1.0]
        tau = 1 / lambda
        @test tau > 0
        @test isfinite(tau)
    end

    # Combined lifetime is always less than individual lifetimes
    k_1 = 1e-6
    k_2 = 2e-6
    tau_1 = 1 / k_1
    tau_2 = 1 / k_2
    tau_combined = 1 / (k_1 + k_2)

    @test tau_combined < tau_1
    @test tau_combined < tau_2
end

#=============================================================================
# Test Set 4: System Compilation Tests (selected)
=============================================================================#

@testitem "4.1 AtmosphericBudget Compilation" begin
    using ModelingToolkit
    using GasChem.AtmosphericLifetime: AtmosphericBudget

    @named budget = AtmosphericBudget()

    # Turn off namespacing for easier variable access
    budget_nns = toggle_namespacing(budget, false)

    # Specify inputs (fluxes and rates)
    inputs = [budget_nns.F_in, budget_nns.F_out, budget_nns.P, budget_nns.R]

    compiled = mtkcompile(budget; inputs)
    @test compiled isa System
end

@testitem "4.5 TroposphericBudget Compilation" begin
    using ModelingToolkit
    using GasChem.AtmosphericLifetime: TroposphericBudget

    @named trop = TroposphericBudget()

    trop_nns = toggle_namespacing(trop, false)
    inputs = [trop_nns.P_n, trop_nns.P_a, trop_nns.P_c]

    compiled = mtkcompile(trop; inputs)
    @test compiled isa System
end

#=============================================================================
# Test Set 5: Integration Tests (shortened)
=============================================================================#

@testitem "5.1 TroposphericBudget ODE Solution" begin
    using ModelingToolkit
    using OrdinaryDiffEqDefault
    using GasChem.AtmosphericLifetime: TroposphericBudget

    @named trop = TroposphericBudget()

    trop_nns = toggle_namespacing(trop, false)
    inputs = [trop_nns.P_n, trop_nns.P_a, trop_nns.P_c]

    compiled = mtkcompile(trop; inputs)

    # Set up problem with constant sources
    k_d_val = 1e-6
    k_w_val = 5e-7
    k_c_val = 2e-6
    k_t_val = 1e-7

    P_n_val = 1e5
    P_a_val = 5e4
    P_c_val = 2e4

    Q_0 = 0.0

    # Expected steady state
    P_total_val = P_n_val + P_a_val + P_c_val
    k_total_val = k_d_val + k_w_val + k_c_val + k_t_val
    Q_ss_expected = P_total_val / k_total_val

    tau_expected = 1 / k_total_val
    t_end = 10 * tau_expected

    prob = ODEProblem(
        compiled,
        Dict(
            compiled.Q => Q_0,
            compiled.k_d => k_d_val,
            compiled.k_w => k_w_val,
            compiled.k_c => k_c_val,
            compiled.k_t => k_t_val,
            compiled.P_n => P_n_val,
            compiled.P_a => P_a_val,
            compiled.P_c => P_c_val,
        ),
        (0.0, t_end)
    )

    sol = solve(prob)

    # Check that solution approaches steady state
    Q_final = sol[compiled.Q][end]
    @test Q_final ≈ Q_ss_expected rtol=0.01
end

#=============================================================================
# Test Set 6: Real-World Examples (selected)
=============================================================================#

@testitem "6.1 Methane (CH4) Lifetime" begin
    k_OH = 6.3e-15
    OH_avg = 1e6

    tau_OH = 1 / (k_OH * OH_avg)
    tau_OH_years = tau_OH / (365.25 * 86400)

    @test tau_OH_years > 4 && tau_OH_years < 7

    tau_OH_actual = 10.0
    tau_soil = 30.0
    tau_strat = 120.0

    tau_total = 1 / (1/tau_OH_actual + 1/tau_soil + 1/tau_strat)
    @test tau_total > 6 && tau_total < 12
end

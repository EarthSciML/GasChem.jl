@testsnippet SP_CH6_Setup begin
    using Test
    using ModelingToolkit
    using ModelingToolkit: t, D
    using NonlinearSolve
    using GasChem

    """Helper to compile an algebraic system with specified input variables."""
    function compile_with_inputs(sys, input_names::Vector{Symbol})
        sys_nns = toggle_namespacing(sys, false)
        input_vars = [getproperty(sys_nns, n) for n in input_names]
        return mtkcompile(sys; inputs = input_vars)
    end
end

# ===========================================================================
# Structural Tests
# ===========================================================================
@testitem "OHProduction: Structural Verification" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = OHProduction()
    @test sys isa System
    @test nameof(sys) == :OHProduction

    vars = unknowns(sys)
    params = parameters(sys)
    eqs = equations(sys)

    # OHProduction has 6 variables: O3, H2O, M (input) + O1D, P_OH, epsilon_OH (output)
    @test length(vars) == 6

    # 3 algebraic equations (Eq. 6.1, 6.3, 6.4)
    @test length(eqs) == 3

    # Check expected variable names are present
    var_names = [string(v) for v in vars]
    for expected in ["O1D", "P_OH", "ε_OH", "O3", "H2O", "M"]
        @test any(n -> contains(n, expected), var_names)
    end

    # Verify it compiles successfully
    compiled = compile_with_inputs(sys, [:O3, :H2O, :M])
    @test compiled !== nothing
end

# ===========================================================================
# Equation Verification Tests (exercising the actual MTK system)
# ===========================================================================
@testitem "OHProduction: Eq 6.1 O(1D) Steady-State" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = OHProduction()
    compiled = compile_with_inputs(sys, [:O3, :H2O, :M])

    # Typical conditions (SI units, m⁻³):
    O3_val = 1e18
    M_val = 2.5e25
    H2O_val = 4e23

    prob = NonlinearProblem(compiled,
        Dict(compiled.O3 => O3_val, compiled.H2O => H2O_val, compiled.M => M_val);
        build_initializeprob = false)
    sol = solve(prob)

    # [O(1D)] should be a very small concentration
    O1D_val = sol[compiled.O1D]
    @test O1D_val > 0
    @test O1D_val < 1e6

    # Verify against hand calculation:
    # k3_eff = 0.78*2.6e-17 + 0.21*4.0e-17 = 2.868e-17
    # denom = 2.868e-17 * 2.5e25 + 2.2e-16 * 4e23 = 7.17e8 + 8.8e7 = 8.05e8
    # O1D = 6e-5 * 1e18 / 8.05e8 = 7.453e4
    @test O1D_val ≈ 7.453e4 rtol = 0.01
end

@testitem "OHProduction: Eq 6.4 OH Yield" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = OHProduction()
    compiled = compile_with_inputs(sys, [:O3, :H2O, :M])

    O3_val = 1e18
    M_val = 2.5e25
    H2O_val = 4e23

    prob = NonlinearProblem(compiled,
        Dict(compiled.O3 => O3_val, compiled.H2O => H2O_val, compiled.M => M_val);
        build_initializeprob = false)
    sol = solve(prob)

    eps_val = sol[compiled.ε_OH]

    # Eq. 6.4: ε_OH = 2 k₄[H₂O] / (k₃[M] + k₄[H₂O])
    # At ~50% RH (H2O = 4e23 m⁻³), ε_OH ≈ 0.22 (from Table on p. 207)
    # Hand calc: 2 * 2.2e-16 * 4e23 / (2.868e-17 * 2.5e25 + 2.2e-16 * 4e23)
    #          = 1.76e8 / 8.05e8 = 0.2186
    @test eps_val > 0.1
    @test eps_val < 0.5
    @test eps_val ≈ 0.2186 rtol = 0.05
end

@testitem "OHProduction: Eq 6.3 OH Production Rate" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = OHProduction()
    compiled = compile_with_inputs(sys, [:O3, :H2O, :M])

    O3_val = 1e18
    M_val = 2.5e25
    H2O_val = 4e23

    prob = NonlinearProblem(compiled,
        Dict(compiled.O3 => O3_val, compiled.H2O => H2O_val, compiled.M => M_val);
        build_initializeprob = false)
    sol = solve(prob)

    P_OH_val = sol[compiled.P_OH]

    # P_OH = j_O3 * [O3] * ε_OH = 6e-5 * 1e18 * 0.2186 = 1.312e13 m⁻³/s
    @test P_OH_val > 1e12
    @test P_OH_val < 1e14
    @test P_OH_val ≈ 1.312e13 rtol = 0.05
end

# ===========================================================================
# Rate Constants at 298 K
# ===========================================================================
@testitem "OHProduction: Rate Constants at 298 K" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = OHProduction()
    params = parameters(sys)
    param_dict = Dict(Symbol(p) => ModelingToolkit.getdefault(p)
                      for p in params if ModelingToolkit.hasdefault(p))

    # From Table B.1 of Seinfeld & Pandis, converted to SI (m³/s)
    @test param_dict[:k3_N2] ≈ 2.6e-11 * 1e-6  # m³/s
    @test param_dict[:k3_O2] ≈ 4.0e-11 * 1e-6  # m³/s
    @test param_dict[:k4] ≈ 2.2e-10 * 1e-6     # m³/s
    @test param_dict[:j_O3] ≈ 6e-5              # s⁻¹ (at surface, solar zenith 0°)
    @test param_dict[:f_N2] ≈ 0.78
    @test param_dict[:f_O2] ≈ 0.21
end

# ===========================================================================
# Limiting Behavior Tests (using actual MTK system)
# ===========================================================================
@testitem "OHProduction: Limiting Behavior - Humidity Variation" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = OHProduction()
    compiled = compile_with_inputs(sys, [:O3, :H2O, :M])

    O3_val = 1e18
    M_val = 2.5e25

    # Low humidity: epsilon_OH should be small
    prob_low = NonlinearProblem(compiled,
        Dict(compiled.O3 => O3_val, compiled.H2O => 1e19, compiled.M => M_val);
        build_initializeprob = false)
    sol_low = solve(prob_low)
    @test sol_low[compiled.ε_OH] < 0.02  # Nearly all O(1D) quenched

    # High humidity: epsilon_OH should be larger
    prob_high = NonlinearProblem(compiled,
        Dict(compiled.O3 => O3_val, compiled.H2O => 1e24, compiled.M => M_val);
        build_initializeprob = false)
    sol_high = solve(prob_high)
    @test sol_high[compiled.ε_OH] > sol_low[compiled.ε_OH]
    # ε_OH can be up to 2 (2 OH per O(¹D) + H₂O reaction)
    @test sol_high[compiled.ε_OH] < 2.0
end

# ===========================================================================
# Qualitative Property Tests (using actual MTK system)
# ===========================================================================
@testitem "OHProduction: Eq 6.4 ε_OH vs Book Table (p. 207)" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = OHProduction()
    compiled = compile_with_inputs(sys, [:O3, :H2O, :M])

    O3_val = 1e18
    M_val = 2.5e25

    # From p. 207: At 298K, at the surface
    # k4/k3 = 7.6, p_H2O_sat at 288K ≈ 0.0167, ξ_H2O = RH * p_H2O_sat
    # The atmospheric k3 = 2.9e-11 cm³/molec/s = 2.9e-17 m³/s
    # k4 = 2.2e-10 cm³/molec/s = 2.2e-16 m³/s
    # ε_OH = 2 * k4 * [H2O] / (k3_eff * [M] + k4 * [H2O])
    # with k3_eff = 0.78*2.6e-17 + 0.21*4.0e-17 = 2.868e-17

    # Book values at 298K (Table on p. 207):
    # RH(%)    ε_OH
    # 10       0.047
    # 25       0.12
    # 50       0.23
    # 80       0.38

    # At 298K, saturated H2O mixing ratio ≈ 0.031 (p_sat/p_atm ≈ 0.031)
    # [H2O] = RH * ξ_H2O_sat * M where ξ_H2O_sat ≈ 0.031
    xi_H2O_sat = 0.031

    # The book's table uses the approximate formula ε_OH ≈ 2*k4*ξ_H2O/k3 (Eq. 6.4),
    # while our implementation uses the exact form. At low RH these agree well;
    # at high RH the approximation overestimates. Also the book uses k3=2.9e-11
    # (atmospheric average) while we use the weighted k3_eff.
    # We test at low/moderate RH where exact and approximate forms agree well.
    book_values = [(10, 0.047), (25, 0.12), (50, 0.23)]

    for (rh, eps_book) in book_values
        H2O_val = (rh / 100) * xi_H2O_sat * M_val
        prob = NonlinearProblem(compiled,
            Dict(compiled.O3 => O3_val, compiled.H2O => H2O_val, compiled.M => M_val);
            build_initializeprob = false)
        sol = solve(prob)
        # Allow 10% tolerance for k3_eff vs k3 difference and exact vs approx formula
        @test sol[compiled.ε_OH] ≈ eps_book rtol = 0.10
    end
end

@testitem "OHProduction: Qualitative - O1D Proportional to O3" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = OHProduction()
    compiled = compile_with_inputs(sys, [:O3, :H2O, :M])

    M_val = 2.5e25
    H2O_val = 4e23

    # Solve with O3 = 1e18
    prob_a = NonlinearProblem(compiled,
        Dict(compiled.O3 => 1e18, compiled.H2O => H2O_val, compiled.M => M_val);
        build_initializeprob = false)
    sol_a = solve(prob_a)

    # Solve with O3 = 2e18
    prob_b = remake(prob_a, p = [compiled.O3 => 2e18])
    sol_b = solve(prob_b)

    # Doubling O3 should double O(1D) (linear relationship from Eq 6.1)
    @test sol_b[compiled.O1D] / sol_a[compiled.O1D] ≈ 2.0 rtol = 1e-6
end

@testitem "OHProduction: Qualitative - Positive Quantities" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = OHProduction()
    compiled = compile_with_inputs(sys, [:O3, :H2O, :M])

    # Test across a range of conditions
    for M_val in [1e24, 2.5e25, 1e26]
        for H2O_val in [1e21, 4e23, 1e24]
            for O3_val in [1e16, 1e18, 1e20]
                prob = NonlinearProblem(compiled,
                    Dict(compiled.O3 => O3_val, compiled.H2O => H2O_val, compiled.M => M_val);
                    build_initializeprob = false)
                sol = solve(prob)

                @test sol[compiled.O1D] >= 0
                @test sol[compiled.ε_OH] >= 0
                # ε_OH = 2 k₄[H₂O] / (k₃[M] + k₄[H₂O]) can be up to 2
                # (2 OH produced per O(¹D) + H₂O reaction)
                @test sol[compiled.ε_OH] <= 2
                @test sol[compiled.P_OH] >= 0
            end
        end
    end
end

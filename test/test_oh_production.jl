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
    # O1D = 1e-5 * 1e18 / 8.05e8 = 1.242e4
    @test O1D_val ≈ 1.242e4 rtol = 0.01
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

    # Typical tropospheric epsilon_OH is 0.05-0.15
    @test eps_val > 0.05
    @test eps_val < 0.15
    @test eps_val ≈ 0.109 rtol = 0.05
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

    # Typical OH production rate ~2e12 m⁻³/s
    @test P_OH_val > 1e11
    @test P_OH_val < 1e13
    @test P_OH_val ≈ 2.18e12 rtol = 0.05
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
    @test param_dict[:j_O3] ≈ 1e-5              # s⁻¹
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
    @test sol_low[compiled.ε_OH] < 0.01  # Nearly all O(1D) quenched

    # High humidity: epsilon_OH should be larger
    prob_high = NonlinearProblem(compiled,
        Dict(compiled.O3 => O3_val, compiled.H2O => 1e24, compiled.M => M_val);
        build_initializeprob = false)
    sol_high = solve(prob_high)
    @test sol_high[compiled.ε_OH] > sol_low[compiled.ε_OH]
    @test sol_high[compiled.ε_OH] < 1.0
end

# ===========================================================================
# Qualitative Property Tests (using actual MTK system)
# ===========================================================================
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
                @test sol[compiled.ε_OH] <= 1
                @test sol[compiled.P_OH] >= 0
            end
        end
    end
end

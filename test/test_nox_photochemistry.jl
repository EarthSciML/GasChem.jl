# ===========================================================================
# Structural Tests
# ===========================================================================
@testitem "NOxPhotochemistry: Structural Verification" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=NOxPhotochemistry()
    @test sys isa System
    @test nameof(sys) == :NOxPhotochemistry

    vars=unknowns(sys)
    params=parameters(sys)
    eqs=equations(sys)

    # 9 variables: NO, NO2, O3, O2, M (input) + O, O3_pss, Phi, P_O3 (output)
    @test length(vars) == 9

    # 4 algebraic equations (Eq. 6.5, 6.6, 6.7, 6.8)
    @test length(eqs) == 4

    # 3 parameters: j_NO2, k_O_O2_M, k_NO_O3
    @test length(params) == 3

    # Check output variable names
    var_names=[string(v) for v in vars]
    for expected in ["O3_pss", "P_O3"]
        @test any(n -> contains(n, expected), var_names)
    end

    # Verify it compiles
    compiled=compile_with_inputs(sys, [:NO, :NO2, :O3, :O2, :M])
    @test compiled !== nothing
end

@testitem "PhotostationaryState: Structural Verification" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=PhotostationaryState()
    @test sys isa System
    @test nameof(sys) == :PhotostationaryState

    vars=unknowns(sys)
    eqs=equations(sys)

    # 5 variables: NO, NO2, O3 (input) + Phi, Phi_deviation (output)
    @test length(vars) == 5

    # 2 equations
    @test length(eqs) == 2

    compiled=compile_with_inputs(sys, [:NO, :NO2, :O3])
    @test compiled !== nothing
end

# ===========================================================================
# Equation Verification Tests (exercising actual MTK system)
# ===========================================================================
@testitem "NOxPhotochemistry: Eq 6.5 O Atom Steady-State" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=NOxPhotochemistry()
    compiled=compile_with_inputs(sys, [:NO, :NO2, :O3, :O2, :M])

    # Typical conditions (SI)
    NO_val=2.5e16    # m⁻³ (1 ppb)
    NO2_val=5e16     # m⁻³ (2 ppb)
    O3_val=1e18      # m⁻³ (~40 ppb)
    O2_val=5.25e24   # m⁻³
    M_val=2.5e25     # m⁻³

    prob=NonlinearProblem(compiled,
        Dict(compiled.NO=>NO_val, compiled.NO2=>NO2_val, compiled.O3=>O3_val,
            compiled.O2=>O2_val, compiled.M=>M_val);
        build_initializeprob = false)
    sol=solve(prob)

    # O atoms are very short-lived, concentration should be very low
    O_val=sol[compiled.O]
    @test O_val > 0
    @test O_val < 1e12
end

@testitem "NOxPhotochemistry: Eq 6.6 Leighton Relationship" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=NOxPhotochemistry()
    compiled=compile_with_inputs(sys, [:NO, :NO2, :O3, :O2, :M])

    # Test at typical urban conditions
    NO_val=2.5e16    # 1 ppb
    NO2_val=5e16     # 2 ppb
    O3_val=1e18
    O2_val=5.25e24
    M_val=2.5e25

    prob=NonlinearProblem(compiled,
        Dict(compiled.NO=>NO_val, compiled.NO2=>NO2_val, compiled.O3=>O3_val,
            compiled.O2=>O2_val, compiled.M=>M_val);
        build_initializeprob = false)
    sol=solve(prob)

    O3_pss=sol[compiled.O3_pss]

    # Photostationary state O3 should be a reasonable tropospheric value
    @test O3_pss > 1e17
    @test O3_pss < 1e19

    # O3_pss = j_NO2 * [NO2] / (k_NO_O3 * [NO])
    # = 8e-3 * 5e16 / (1.9e-20 * 2.5e16) = 4e14 / 4.75e-4 = 8.42e17
    @test O3_pss ≈ 8.42e17 rtol = 0.05
end

@testitem "NOxPhotochemistry: Eq 6.7 Photostationary State Parameter" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=NOxPhotochemistry()
    compiled=compile_with_inputs(sys, [:NO, :NO2, :O3, :O2, :M])

    NO_val=2.5e16
    NO2_val=5e16
    O2_val=5.25e24
    M_val=2.5e25

    # Compute O3_pss for these conditions
    j_NO2=8e-3
    k_NO_O3=1.9e-14*1e-6
    O3_pss_val=j_NO2*NO2_val/(k_NO_O3*NO_val)

    # At PSS value of O3: Phi should be 1
    prob=NonlinearProblem(compiled,
        Dict(compiled.NO=>NO_val, compiled.NO2=>NO2_val, compiled.O3=>O3_pss_val,
            compiled.O2=>O2_val, compiled.M=>M_val);
        build_initializeprob = false)
    sol=solve(prob)
    @test sol[compiled.Φ] ≈ 1.0 rtol = 1e-6

    # Below PSS O3: Phi > 1
    prob_low=remake(prob, p = [compiled.O3=>O3_pss_val*0.5])
    sol_low=solve(prob_low)
    @test sol_low[compiled.Φ] > 1.0

    # Above PSS O3: Phi < 1
    prob_high=remake(prob, p = [compiled.O3=>O3_pss_val*2.0])
    sol_high=solve(prob_high)
    @test sol_high[compiled.Φ] < 1.0
end

@testitem "NOxPhotochemistry: Eq 6.8 Net O3 Production" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=NOxPhotochemistry()
    compiled=compile_with_inputs(sys, [:NO, :NO2, :O3, :O2, :M])

    NO_val=2.5e16
    NO2_val=5e16
    O2_val=5.25e24
    M_val=2.5e25

    j_NO2=8e-3
    k_NO_O3=1.9e-14*1e-6
    O3_pss_val=j_NO2*NO2_val/(k_NO_O3*NO_val)

    # At photostationary state: P_O3 ≈ 0
    prob=NonlinearProblem(compiled,
        Dict(compiled.NO=>NO_val, compiled.NO2=>NO2_val, compiled.O3=>O3_pss_val,
            compiled.O2=>O2_val, compiled.M=>M_val);
        build_initializeprob = false)
    sol=solve(prob)
    @test abs(sol[compiled.P_O3]) / (j_NO2 * NO2_val) < 1e-6

    # Below PSS: net production
    prob_low=remake(prob, p = [compiled.O3=>O3_pss_val*0.5])
    sol_low=solve(prob_low)
    @test sol_low[compiled.P_O3] > 0

    # Above PSS: net loss
    prob_high=remake(prob, p = [compiled.O3=>O3_pss_val*2.0])
    sol_high=solve(prob_high)
    @test sol_high[compiled.P_O3] < 0
end

# ===========================================================================
# Rate Constants at 298 K
# ===========================================================================
@testitem "NOxPhotochemistry: Rate Constants at 298 K" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=NOxPhotochemistry()
    params=parameters(sys)
    param_dict=Dict(Symbol(p)=>ModelingToolkit.getdefault(p)
    for p in params if ModelingToolkit.hasdefault(p))

    @test param_dict[:j_NO2] ≈ 8e-3
    @test param_dict[:k_O_O2_M] ≈ 6.0e-34 * 1e-12
    @test param_dict[:k_NO_O3] ≈ 1.9e-14 * 1e-6
end

# ===========================================================================
# Steady-State and Limiting Behavior Tests
# ===========================================================================
@testitem "NOxPhotochemistry: O3_pss Scales with NO2/NO Ratio" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=NOxPhotochemistry()
    compiled=compile_with_inputs(sys, [:NO, :NO2, :O3, :O2, :M])

    O2_val=5.25e24
    M_val=2.5e25
    NO_val=1e16

    prob=NonlinearProblem(compiled,
        Dict(compiled.NO=>NO_val, compiled.NO2=>1e16, compiled.O3=>1e18,
            compiled.O2=>O2_val, compiled.M=>M_val);
        build_initializeprob = false)
    sol_a=solve(prob)

    prob_b=remake(prob, p = [compiled.NO2=>2e16])
    sol_b=solve(prob_b)

    # O3_pss should scale linearly with NO2
    @test sol_b[compiled.O3_pss] / sol_a[compiled.O3_pss] ≈ 2.0 rtol = 1e-6
end

@testitem "PhotostationaryState: Deviation from PSS" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=PhotostationaryState()
    compiled=compile_with_inputs(sys, [:NO, :NO2, :O3])

    j_NO2=8e-3
    k_NO_O3=1.9e-14*1e-6
    NO_val=2.5e16
    NO2_val=5e16
    O3_pss_val=j_NO2*NO2_val/(k_NO_O3*NO_val)

    # At PSS: deviation = 0
    prob=NonlinearProblem(compiled,
        Dict(compiled.NO=>NO_val, compiled.NO2=>NO2_val, compiled.O3=>O3_pss_val);
        build_initializeprob = false)
    sol=solve(prob)
    @test sol[compiled.Φ_deviation] ≈ 0.0 atol = 1e-6

    # Below PSS O3: positive deviation
    prob_low=remake(prob, p = [compiled.O3=>O3_pss_val*0.5])
    sol_low=solve(prob_low)
    @test sol_low[compiled.Φ_deviation] > 0

    # Above PSS O3: negative deviation
    prob_high=remake(prob, p = [compiled.O3=>O3_pss_val*2.0])
    sol_high=solve(prob_high)
    @test sol_high[compiled.Φ_deviation] < 0
end

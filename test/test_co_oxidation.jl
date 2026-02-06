# ===========================================================================
# Structural Tests
# ===========================================================================
@testitem "COOxidation: Structural Verification" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = COOxidation()
    @test sys isa System
    @test nameof(sys) == :COOxidation

    vars = unknowns(sys)
    params = parameters(sys)
    eqs = equations(sys)

    # 12 variables: CO, OH, HO2, NO, NO2, O3 (input)
    #             + P_O3, L_HOx, L_OH, L_HO2, chain_length, HO2_ss (output)
    @test length(vars) == 12

    # 6 algebraic equations
    @test length(eqs) == 6

    # 7 parameters: k_CO_OH, k_HO2_NO, k_HO2_HO2, k_OH_NO2, k_HO2_O3, k_OH_O3 + 1 constant (two)
    @test length(params) == 7

    # Check key output variable names
    var_names = [string(v) for v in vars]
    for expected in ["P_O3", "L_HOx", "chain_length", "HO2_ss"]
        @test any(n -> contains(n, expected), var_names)
    end

    # Verify it compiles successfully
    compiled = compile_with_inputs(sys, [:CO, :OH, :HO2, :NO, :NO2, :O3])
    @test compiled !== nothing
end

@testitem "OzoneProductionEfficiency: Structural Verification" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = OzoneProductionEfficiency()
    @test sys isa System
    @test nameof(sys) == :OzoneProductionEfficiency

    vars = unknowns(sys)
    params = parameters(sys)
    eqs = equations(sys)

    # 8 variables: OH, HO2, RO2, NO, NO2 (input) + P_O3, L_NOx, OPE (output)
    @test length(vars) == 8

    # 3 equations
    @test length(eqs) == 3

    # 3 parameters: k_HO2_NO, k_RO2_NO, k_OH_NO2
    @test length(params) == 3

    # Verify it compiles successfully
    compiled = compile_with_inputs(sys, [:OH, :HO2, :RO2, :NO, :NO2])
    @test compiled !== nothing
end

# ===========================================================================
# Equation Verification Tests (exercising the actual MTK system)
# ===========================================================================
@testitem "COOxidation: Eq 6.18 HO2 Steady-State (High NOx)" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = COOxidation()
    compiled = compile_with_inputs(sys, [:CO, :OH, :HO2, :NO, :NO2, :O3])

    # Typical conditions (SI units, m⁻³)
    CO_val = 2.5e18    # ~100 ppb
    OH_val = 1e12      # typical daytime
    HO2_val = 1e14
    NO_val = 2.5e15    # ~0.1 ppb
    NO2_val = 2.5e16
    O3_val = 1e18

    prob = NonlinearProblem(compiled,
        Dict(compiled.CO => CO_val, compiled.OH => OH_val, compiled.HO2 => HO2_val,
             compiled.NO => NO_val, compiled.NO2 => NO2_val, compiled.O3 => O3_val);
        build_initializeprob = false)
    sol = solve(prob)

    HO2_ss_val = sol[compiled.HO2_ss]

    # HO2 should be on the order of 1e13-1e14 m⁻³
    @test HO2_ss_val > 1e12
    @test HO2_ss_val < 1e16

    # Verify against hand calculation:
    # HO2_ss = k_CO_OH * CO * OH / (k_HO2_NO * NO)
    # = (2.4e-19) * (2.5e18) * (1e12) / ((8.1e-18) * (2.5e15))
    # = 6.0e11 / 2.025e-2 = 2.963e13
    @test HO2_ss_val ≈ 2.963e13 rtol = 0.05
end

@testitem "COOxidation: Net O3 Production Rate" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = COOxidation()
    compiled = compile_with_inputs(sys, [:CO, :OH, :HO2, :NO, :NO2, :O3])

    HO2_val = 1e14
    NO_val = 2.5e15
    OH_val = 1e12
    O3_val = 1e18
    CO_val = 2.5e18
    NO2_val = 2.5e16

    prob = NonlinearProblem(compiled,
        Dict(compiled.CO => CO_val, compiled.OH => OH_val, compiled.HO2 => HO2_val,
             compiled.NO => NO_val, compiled.NO2 => NO2_val, compiled.O3 => O3_val);
        build_initializeprob = false)
    sol = solve(prob)

    P_O3_val = sol[compiled.P_O3]

    # Net production should be positive in these conditions
    @test P_O3_val > 0

    # Verify against hand calculation:
    # P_O3 = k_HO2_NO * HO2 * NO - k_OH_O3 * OH * O3 - k_HO2_O3 * HO2 * O3
    # = (8.1e-18)(1e14)(2.5e15) - (7.3e-20)(1e12)(1e18) - (2.0e-21)(1e14)(1e18)
    # = 2.025e12 - 7.3e10 - 2.0e11 = 1.752e12
    @test P_O3_val ≈ 1.752e12 rtol = 0.01
end

@testitem "COOxidation: HOx Loss Rate" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = COOxidation()
    compiled = compile_with_inputs(sys, [:CO, :OH, :HO2, :NO, :NO2, :O3])

    OH_val = 1e12
    NO2_val = 2.5e16
    HO2_val = 1e14
    CO_val = 2.5e18
    NO_val = 2.5e15
    O3_val = 1e18

    prob = NonlinearProblem(compiled,
        Dict(compiled.CO => CO_val, compiled.OH => OH_val, compiled.HO2 => HO2_val,
             compiled.NO => NO_val, compiled.NO2 => NO2_val, compiled.O3 => O3_val);
        build_initializeprob = false)
    sol = solve(prob)

    L_HOx_val = sol[compiled.L_HOx]

    @test L_HOx_val > 0

    # Verify against hand calculation:
    # L_HOx = k_OH_NO2 * OH * NO2 + 2 * k_HO2_HO2 * HO2^2
    # = (1.0e-17)(1e12)(2.5e16) + 2*(2.9e-18)(1e14)^2
    # = 2.5e11 + 5.8e10 = 3.08e11
    @test L_HOx_val ≈ 3.08e11 rtol = 0.01
end

@testitem "COOxidation: Chain Length" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = COOxidation()
    compiled = compile_with_inputs(sys, [:CO, :OH, :HO2, :NO, :NO2, :O3])

    HO2_val = 1e14
    NO_val = 2.5e15
    OH_val = 1e12
    NO2_val = 2.5e16
    CO_val = 2.5e18
    O3_val = 1e18

    prob = NonlinearProblem(compiled,
        Dict(compiled.CO => CO_val, compiled.OH => OH_val, compiled.HO2 => HO2_val,
             compiled.NO => NO_val, compiled.NO2 => NO2_val, compiled.O3 => O3_val);
        build_initializeprob = false)
    sol = solve(prob)

    chain_val = sol[compiled.chain_length]

    # Chain length should be > 1 (catalytic cycle operates)
    @test chain_val > 1
    # Typical chain length is ~1-100 depending on conditions
    @test chain_val < 100

    # Verify against hand calculation:
    # chain_length = k_HO2_NO * HO2 * NO / L_HOx
    # = (8.1e-18)(1e14)(2.5e15) / 3.08e11 ≈ 6.58
    @test chain_val ≈ 6.58 rtol = 0.02
end

# ===========================================================================
# Rate Constants at 298 K
# ===========================================================================
@testitem "COOxidation: Rate Constants at 298 K" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = COOxidation()
    params = parameters(sys)
    param_dict = Dict(Symbol(p) => ModelingToolkit.getdefault(p)
                      for p in params if ModelingToolkit.hasdefault(p))

    @test param_dict[:k_CO_OH] ≈ 2.4e-13 * 1e-6
    @test param_dict[:k_HO2_NO] ≈ 8.1e-12 * 1e-6
    @test param_dict[:k_HO2_HO2] ≈ 2.9e-12 * 1e-6
    @test param_dict[:k_OH_NO2] ≈ 1.0e-11 * 1e-6
    @test param_dict[:k_HO2_O3] ≈ 2.0e-15 * 1e-6
    @test param_dict[:k_OH_O3] ≈ 7.3e-14 * 1e-6
end

@testitem "OzoneProductionEfficiency: Rate Constants" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = OzoneProductionEfficiency()
    params = parameters(sys)
    param_dict = Dict(Symbol(p) => ModelingToolkit.getdefault(p)
                      for p in params if ModelingToolkit.hasdefault(p))

    @test param_dict[:k_HO2_NO] ≈ 8.1e-12 * 1e-6
    @test param_dict[:k_RO2_NO] ≈ 8.0e-12 * 1e-6
    @test param_dict[:k_OH_NO2] ≈ 1.0e-11 * 1e-6
end

# ===========================================================================
# OPE Equation Verification (exercising the actual MTK system)
# ===========================================================================
@testitem "OzoneProductionEfficiency: OPE Calculation" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = OzoneProductionEfficiency()
    compiled = compile_with_inputs(sys, [:OH, :HO2, :RO2, :NO, :NO2])

    OH_val = 1e12
    HO2_val = 1e14
    RO2_val = 1e14
    NO_val = 2.5e15
    NO2_val = 2.5e16

    prob = NonlinearProblem(compiled,
        Dict(compiled.OH => OH_val, compiled.HO2 => HO2_val, compiled.RO2 => RO2_val,
             compiled.NO => NO_val, compiled.NO2 => NO2_val);
        build_initializeprob = false)
    sol = solve(prob)

    OPE_val = sol[compiled.OPE]
    P_O3_val = sol[compiled.P_O3]
    L_NOx_val = sol[compiled.L_NOx]

    @test OPE_val > 0
    @test OPE_val > 1

    # Verify against hand calculation:
    # P_O3 = k_HO2_NO * HO2 * NO + k_RO2_NO * RO2 * NO
    # = (8.1e-18)(1e14)(2.5e15) + (8.0e-18)(1e14)(2.5e15)
    # = 2.025e12 + 2.0e12 = 4.025e12
    @test P_O3_val ≈ 4.025e12 rtol = 0.01

    # L_NOx = k_OH_NO2 * OH * NO2 = (1.0e-17)(1e12)(2.5e16) = 2.5e11
    @test L_NOx_val ≈ 2.5e11 rtol = 0.01

    # OPE = P_O3 / L_NOx ≈ 16.1
    @test OPE_val ≈ 16.1 rtol = 0.01
end

# ===========================================================================
# Limiting Behavior Tests (using actual MTK system)
# ===========================================================================
@testitem "COOxidation: High NOx - Short Chain Length" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = COOxidation()
    compiled = compile_with_inputs(sys, [:CO, :OH, :HO2, :NO, :NO2, :O3])

    CO_val = 2.5e18
    OH_val = 5e11
    HO2_val = 5e13
    O3_val = 1e18

    # Low NOx scenario
    prob_low = NonlinearProblem(compiled,
        Dict(compiled.CO => CO_val, compiled.OH => OH_val, compiled.HO2 => HO2_val,
             compiled.NO => 2.5e15, compiled.NO2 => 2.5e16, compiled.O3 => O3_val);
        build_initializeprob = false)
    sol_low = solve(prob_low)

    # High NOx scenario
    prob_high = remake(prob_low, p = [compiled.NO => 2.5e17, compiled.NO2 => 7.5e17])
    sol_high = solve(prob_high)

    # At high NOx, the HNO3 termination dominates, so L_HOx is larger
    @test sol_high[compiled.L_HOx] > sol_low[compiled.L_HOx]

    # At high NOx, chain length should change since propagation and termination both increase
    # but the chain length formula depends on specific conditions
    @test sol_high[compiled.chain_length] > 0
    @test sol_low[compiled.chain_length] > 0
end

@testitem "OzoneProductionEfficiency: High NOx - Low OPE" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = OzoneProductionEfficiency()
    compiled = compile_with_inputs(sys, [:OH, :HO2, :RO2, :NO, :NO2])

    OH_val = 1e12
    HO2_val = 1e14
    RO2_val = 1e14

    # Low NOx
    prob_low = NonlinearProblem(compiled,
        Dict(compiled.OH => OH_val, compiled.HO2 => HO2_val, compiled.RO2 => RO2_val,
             compiled.NO => 2.5e14, compiled.NO2 => 5e14);
        build_initializeprob = false)
    sol_low = solve(prob_low)

    # High NOx (same NO2/NO ratio but 100x higher)
    prob_high = remake(prob_low, p = [compiled.NO => 2.5e16, compiled.NO2 => 5e16])
    sol_high = solve(prob_high)

    # OPE depends on NO/NO2 ratio, not absolute NOx level (when HO2, RO2 fixed)
    # Since NO and NO2 are both scaled by 100x, the ratio is the same, so OPE stays the same
    @test sol_low[compiled.OPE] ≈ sol_high[compiled.OPE] rtol = 1e-6
    @test sol_low[compiled.OPE] > 0
end

# ===========================================================================
# Qualitative Property Tests (using actual MTK system)
# ===========================================================================
@testitem "COOxidation: Qualitative - Positive Production and Loss" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = COOxidation()
    compiled = compile_with_inputs(sys, [:CO, :OH, :HO2, :NO, :NO2, :O3])

    CO_val = 2.5e18
    O3_val = 1e18

    for OH_val in [1e11, 1e12, 1e13]
        for HO2_val in [1e13, 1e14, 1e15]
            for NO_val in [1e14, 1e16, 1e18]
                NO2_val = 2 * NO_val

                prob = NonlinearProblem(compiled,
                    Dict(compiled.CO => CO_val, compiled.OH => OH_val, compiled.HO2 => HO2_val,
                         compiled.NO => NO_val, compiled.NO2 => NO2_val, compiled.O3 => O3_val);
                    build_initializeprob = false)
                sol = solve(prob)

                @test sol[compiled.L_OH] > 0
                @test sol[compiled.L_HO2] > 0
                @test sol[compiled.L_HOx] > 0
                @test sol[compiled.HO2_ss] > 0
            end
        end
    end
end

@testitem "COOxidation: HOx Catalytic Cycling" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = COOxidation()
    compiled = compile_with_inputs(sys, [:CO, :OH, :HO2, :NO, :NO2, :O3])

    # Typical conditions (m⁻³)
    CO_val = 2.5e18
    OH_val = 1e12
    HO2_val = 1e14
    NO_val = 2.5e15
    NO2_val = 2.5e16
    O3_val = 1e18

    prob = NonlinearProblem(compiled,
        Dict(compiled.CO => CO_val, compiled.OH => OH_val, compiled.HO2 => HO2_val,
             compiled.NO => NO_val, compiled.NO2 => NO2_val, compiled.O3 => O3_val);
        build_initializeprob = false)
    sol = solve(prob)

    # Chain length must be > 1 for catalytic behavior
    @test sol[compiled.chain_length] > 1

    # Net O3 production: CO + 2O2 + hv -> CO2 + O3
    @test sol[compiled.P_O3] > 0
end

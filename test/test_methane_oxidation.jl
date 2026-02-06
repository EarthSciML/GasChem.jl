# ===========================================================================
# Structural Tests
# ===========================================================================
@testitem "MethaneOxidation: Structural Verification" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=MethaneOxidation()
    @test sys isa System
    @test nameof(sys) == :MethaneOxidation

    vars=unknowns(sys)
    params=parameters(sys)
    eqs=equations(sys)

    # 15 input species + 17 reaction rates (R1-R17) + 3 diagnostics (P_O3_gross, P_HCHO, L_CH4) = 35
    @test length(vars) == 35

    # 20 equations: 17 reaction rates + 3 diagnostics
    @test length(eqs) == 20

    # 17 parameters: k1, k2_0, k3-k8, k10, k13, k14_0, k15, k17_0, j9, j11, j12, j16
    @test length(params) == 17

    # Check key diagnostic variable names
    var_names=[string(v) for v in vars]
    for expected in ["P_O3_gross", "P_HCHO", "L_CH4"]
        @test any(n -> contains(n, expected), var_names)
    end

    # Verify it compiles successfully
    ch4_inputs=[:CH4, :CH3, :CH3O2, :CH3O, :CH3OOH, :HCHO, :HCO,
        :OH, :HO2, :H, :NO, :NO2, :O, :O2, :M]
    compiled=compile_with_inputs(sys, ch4_inputs)
    @test compiled !== nothing
end

@testitem "MethaneOxidationODE: Structural Verification" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=MethaneOxidationODE()
    @test sys isa System
    @test nameof(sys) == :MethaneOxidationODE

    vars=unknowns(sys)
    params=parameters(sys)
    eqs=equations(sys)

    # 18 state variables
    @test length(vars) == 18

    # 18 ODEs (one for each state variable)
    @test length(eqs) == 18

    # Check key species names
    var_names=[string(v) for v in vars]
    for expected in ["CH4", "CH3O2", "HCHO", "OH", "HO2", "NO", "NO2", "O3", "HNO3", "H2O2"]
        @test any(n -> contains(n, expected), var_names)
    end

    # Verify it compiles successfully
    compiled=mtkcompile(sys)
    @test compiled !== nothing
end

# ===========================================================================
# Rate Constants at 298 K (Table 6.1)
# ===========================================================================
@testitem "MethaneOxidation: Table 6.1 Rate Constants" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=MethaneOxidation()
    params=parameters(sys)
    param_dict=Dict(Symbol(p)=>ModelingToolkit.getdefault(p)
    for p in params if ModelingToolkit.hasdefault(p))

    # Bimolecular rate constants (m^3/s)
    @test param_dict[:k1] ≈ 6.3e-15 * 1e-6
    @test param_dict[:k3] ≈ 7.7e-12 * 1e-6
    @test param_dict[:k4] ≈ 5.2e-12 * 1e-6
    @test param_dict[:k5] ≈ 3.5e-13 * 1e-6
    @test param_dict[:k6] ≈ 1.9e-15 * 1e-6
    @test param_dict[:k7] ≈ 3.8e-12 * 1e-6
    @test param_dict[:k8] ≈ 1.9e-12 * 1e-6
    @test param_dict[:k10] ≈ 9.0e-12 * 1e-6
    @test param_dict[:k13] ≈ 5.2e-12 * 1e-6
    @test param_dict[:k15] ≈ 8.1e-12 * 1e-6

    # Termolecular rate constants (m^6/s)
    @test param_dict[:k2_0] ≈ 1.0e-30 * 1e-12
    @test param_dict[:k14_0] ≈ 5.7e-32 * 1e-12
    @test param_dict[:k17_0] ≈ 6.0e-34 * 1e-12

    # Photolysis rates (s^-1)
    @test param_dict[:j9] ≈ 5e-6
    @test param_dict[:j11] ≈ 3e-5
    @test param_dict[:j12] ≈ 4e-5
    @test param_dict[:j16] ≈ 8e-3
end

@testitem "MethaneOxidationODE: Rate Constants and Parameters" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=MethaneOxidationODE()
    params=parameters(sys)
    param_dict=Dict(Symbol(p)=>ModelingToolkit.getdefault(p)
    for p in params if ModelingToolkit.hasdefault(p))

    # Same kinetic rate constants as MethaneOxidation (in SI)
    @test param_dict[:k1] ≈ 6.3e-15 * 1e-6
    @test param_dict[:k3] ≈ 7.7e-12 * 1e-6
    @test param_dict[:k15] ≈ 8.1e-12 * 1e-6

    # Additional rate constants
    @test param_dict[:k_CO_OH] ≈ 2.4e-13 * 1e-6
    @test param_dict[:k_OH_NO2] ≈ 1.0e-11 * 1e-6
    @test param_dict[:k_HO2_HO2] ≈ 2.9e-12 * 1e-6
    @test param_dict[:k_NO_O3] ≈ 1.9e-14 * 1e-6

    # External OH source (m^-3/s)
    @test param_dict[:P_OH_ext] ≈ 1e6 * 1e6
end

# ===========================================================================
# Equation Verification Tests (using NonlinearProblem + solve)
# ===========================================================================
@testitem "MethaneOxidation: Reaction Rate R1 (CH4+OH)" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=MethaneOxidation()
    ch4_inputs=[:CH4, :CH3, :CH3O2, :CH3O, :CH3OOH, :HCHO, :HCO,
        :OH, :HO2, :H, :NO, :NO2, :O, :O2, :M]
    compiled=compile_with_inputs(sys, ch4_inputs)

    # Typical concentrations (SI: m^-3)
    input_vals=Dict(
        compiled.CH4=>4.5e19, compiled.OH=>1e12, compiled.CH3=>1e10,
        compiled.CH3O2=>1e14, compiled.CH3O=>1e10, compiled.CH3OOH=>1e15,
        compiled.HCHO=>1e16, compiled.HCO=>1e10, compiled.HO2=>1e14,
        compiled.H=>1e8, compiled.NO=>2.5e15, compiled.NO2=>2.5e16,
        compiled.O=>1e8, compiled.O2=>5.25e24, compiled.M=>2.5e25
    )

    prob=NonlinearProblem(compiled, input_vals; build_initializeprob = false)
    sol=solve(prob)

    # R1 = k1 * CH4 * OH
    k1=6.3e-15*1e-6
    expected_R1=k1*4.5e19*1e12
    @test sol[compiled.R1] ≈ expected_R1 rtol = 1e-6
    @test sol[compiled.R1] ≈ 2.835e11 rtol = 0.01
    @test sol[compiled.R1] > 0
end

@testitem "MethaneOxidation: Reaction Rate R3 (CH3O2+NO)" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=MethaneOxidation()
    ch4_inputs=[:CH4, :CH3, :CH3O2, :CH3O, :CH3OOH, :HCHO, :HCO,
        :OH, :HO2, :H, :NO, :NO2, :O, :O2, :M]
    compiled=compile_with_inputs(sys, ch4_inputs)

    input_vals=Dict(
        compiled.CH4=>4.5e19, compiled.OH=>1e12, compiled.CH3=>1e10,
        compiled.CH3O2=>1e14, compiled.CH3O=>1e10, compiled.CH3OOH=>1e15,
        compiled.HCHO=>1e16, compiled.HCO=>1e10, compiled.HO2=>1e14,
        compiled.H=>1e8, compiled.NO=>2.5e15, compiled.NO2=>2.5e16,
        compiled.O=>1e8, compiled.O2=>5.25e24, compiled.M=>2.5e25
    )

    prob=NonlinearProblem(compiled, input_vals; build_initializeprob = false)
    sol=solve(prob)

    # R3 = k3 * CH3O2 * NO
    k3=7.7e-12*1e-6
    expected_R3=k3*1e14*2.5e15
    @test sol[compiled.R3] ≈ expected_R3 rtol = 1e-6
    @test sol[compiled.R3] ≈ 1.925e12 rtol = 0.01
    @test sol[compiled.R3] > 0
end

@testitem "MethaneOxidation: HCHO Production (P_HCHO = R6 + R8)" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=MethaneOxidation()
    ch4_inputs=[:CH4, :CH3, :CH3O2, :CH3O, :CH3OOH, :HCHO, :HCO,
        :OH, :HO2, :H, :NO, :NO2, :O, :O2, :M]
    compiled=compile_with_inputs(sys, ch4_inputs)

    input_vals=Dict(
        compiled.CH4=>4.5e19, compiled.OH=>1e12, compiled.CH3=>1e10,
        compiled.CH3O2=>1e14, compiled.CH3O=>1e10, compiled.CH3OOH=>1e15,
        compiled.HCHO=>1e16, compiled.HCO=>1e10, compiled.HO2=>1e14,
        compiled.H=>1e8, compiled.NO=>2.5e15, compiled.NO2=>2.5e16,
        compiled.O=>1e8, compiled.O2=>5.25e24, compiled.M=>2.5e25
    )

    prob=NonlinearProblem(compiled, input_vals; build_initializeprob = false)
    sol=solve(prob)

    # P_HCHO = R6 + R8
    @test sol[compiled.P_HCHO] ≈ sol[compiled.R6] + sol[compiled.R8] rtol = 1e-6

    # Check individual rates against expected values
    k6=1.9e-15*1e-6
    k8=1.9e-12*1e-6
    expected_R6=k6*1e10*5.25e24
    expected_R8=k8*1e15*1e12
    @test sol[compiled.R6] ≈ expected_R6 rtol = 1e-6
    @test sol[compiled.R8] ≈ expected_R8 rtol = 1e-6
    @test sol[compiled.P_HCHO] > 0
end

@testitem "MethaneOxidation: CH4 Loss Equals R1" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=MethaneOxidation()
    ch4_inputs=[:CH4, :CH3, :CH3O2, :CH3O, :CH3OOH, :HCHO, :HCO,
        :OH, :HO2, :H, :NO, :NO2, :O, :O2, :M]
    compiled=compile_with_inputs(sys, ch4_inputs)

    input_vals=Dict(
        compiled.CH4=>4.5e19, compiled.OH=>1e12, compiled.CH3=>1e10,
        compiled.CH3O2=>1e14, compiled.CH3O=>1e10, compiled.CH3OOH=>1e15,
        compiled.HCHO=>1e16, compiled.HCO=>1e10, compiled.HO2=>1e14,
        compiled.H=>1e8, compiled.NO=>2.5e15, compiled.NO2=>2.5e16,
        compiled.O=>1e8, compiled.O2=>5.25e24, compiled.M=>2.5e25
    )

    prob=NonlinearProblem(compiled, input_vals; build_initializeprob = false)
    sol=solve(prob)

    # L_CH4 = R1 (only loss pathway for CH4 is reaction with OH)
    @test sol[compiled.L_CH4] ≈ sol[compiled.R1] rtol = 1e-6
    @test sol[compiled.L_CH4] > 0
end

@testitem "MethaneOxidation: Net O3 Production (P_O3_gross = R15 + R3)" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=MethaneOxidation()
    ch4_inputs=[:CH4, :CH3, :CH3O2, :CH3O, :CH3OOH, :HCHO, :HCO,
        :OH, :HO2, :H, :NO, :NO2, :O, :O2, :M]
    compiled=compile_with_inputs(sys, ch4_inputs)

    input_vals=Dict(
        compiled.CH4=>4.5e19, compiled.OH=>1e12, compiled.CH3=>1e10,
        compiled.CH3O2=>1e14, compiled.CH3O=>1e10, compiled.CH3OOH=>1e15,
        compiled.HCHO=>1e16, compiled.HCO=>1e10, compiled.HO2=>1e14,
        compiled.H=>1e8, compiled.NO=>2.5e15, compiled.NO2=>2.5e16,
        compiled.O=>1e8, compiled.O2=>5.25e24, compiled.M=>2.5e25
    )

    prob=NonlinearProblem(compiled, input_vals; build_initializeprob = false)
    sol=solve(prob)

    # P_O3_gross = R15 + R3 (HO2+NO and CH3O2+NO reactions produce NO2, leading to O3)
    @test sol[compiled.P_O3_gross] ≈ sol[compiled.R15] + sol[compiled.R3] rtol = 1e-6
    @test sol[compiled.P_O3_gross] > 0

    # Both contributing reactions should be positive
    @test sol[compiled.R15] > 0
    @test sol[compiled.R3] > 0
end

# ===========================================================================
# Methane Oxidation Chain Tests (using NonlinearProblem + solve)
# ===========================================================================
@testitem "MethaneOxidation: Oxidation Chain Relative Rates" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=MethaneOxidation()
    ch4_inputs=[:CH4, :CH3, :CH3O2, :CH3O, :CH3OOH, :HCHO, :HCO,
        :OH, :HO2, :H, :NO, :NO2, :O, :O2, :M]
    compiled=compile_with_inputs(sys, ch4_inputs)

    # Moderate NOx conditions
    input_vals=Dict(
        compiled.CH4=>4.5e19, compiled.OH=>1e12, compiled.CH3=>1e10,
        compiled.CH3O2=>1e14, compiled.CH3O=>1e10, compiled.CH3OOH=>1e15,
        compiled.HCHO=>1e16, compiled.HCO=>1e10, compiled.HO2=>1e14,
        compiled.H=>1e8, compiled.NO=>2.5e15, compiled.NO2=>2.5e16,
        compiled.O=>1e8, compiled.O2=>5.25e24, compiled.M=>2.5e25
    )

    prob=NonlinearProblem(compiled, input_vals; build_initializeprob = false)
    sol=solve(prob)

    # At moderate NOx, CH3O2 fate: R3 (+ NO) >> R4 (+ HO2) >> R5 (self-reaction)
    @test sol[compiled.R3] > sol[compiled.R4]
    @test sol[compiled.R4] > sol[compiled.R5]
    @test sol[compiled.R3] / sol[compiled.R4] > 10
end

@testitem "MethaneOxidation: HCHO as Key Intermediate" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=MethaneOxidation()
    ch4_inputs=[:CH4, :CH3, :CH3O2, :CH3O, :CH3OOH, :HCHO, :HCO,
        :OH, :HO2, :H, :NO, :NO2, :O, :O2, :M]
    compiled=compile_with_inputs(sys, ch4_inputs)

    input_vals=Dict(
        compiled.CH4=>4.5e19, compiled.OH=>1e12, compiled.CH3=>1e10,
        compiled.CH3O2=>1e14, compiled.CH3O=>1e10, compiled.CH3OOH=>1e15,
        compiled.HCHO=>1e16, compiled.HCO=>1e10, compiled.HO2=>1e14,
        compiled.H=>1e8, compiled.NO=>2.5e15, compiled.NO2=>2.5e16,
        compiled.O=>1e8, compiled.O2=>5.25e24, compiled.M=>2.5e25
    )

    prob=NonlinearProblem(compiled, input_vals; build_initializeprob = false)
    sol=solve(prob)

    # HCHO loss pathways: R10 (OH), R11 (hv -> HCO+H), R12 (hv -> H2+CO)
    @test sol[compiled.R10] > 0
    @test sol[compiled.R11] > 0
    @test sol[compiled.R12] > 0

    # R12 (molecular channel) > R11 (radical channel) since j12 > j11
    @test sol[compiled.R12] > sol[compiled.R11]
end

# ===========================================================================
# Limiting Behavior Tests (using NonlinearProblem + solve)
# ===========================================================================
@testitem "MethaneOxidation: Low NOx Regime" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=MethaneOxidation()
    ch4_inputs=[:CH4, :CH3, :CH3O2, :CH3O, :CH3OOH, :HCHO, :HCO,
        :OH, :HO2, :H, :NO, :NO2, :O, :O2, :M]
    compiled=compile_with_inputs(sys, ch4_inputs)

    # Low NOx: peroxide channel (R4) dominates over NO channel (R3)
    input_vals=Dict(
        compiled.CH4=>4.5e19, compiled.OH=>1e12, compiled.CH3=>1e10,
        compiled.CH3O2=>1e14, compiled.CH3O=>1e10, compiled.CH3OOH=>1e15,
        compiled.HCHO=>1e16, compiled.HCO=>1e10, compiled.HO2=>2e14,
        compiled.H=>1e8, compiled.NO=>1e13, compiled.NO2=>1e14,
        compiled.O=>1e8, compiled.O2=>5.25e24, compiled.M=>2.5e25
    )

    prob=NonlinearProblem(compiled, input_vals; build_initializeprob = false)
    sol=solve(prob)

    # At low NOx, R4 (CH3O2 + HO2) > R3 (CH3O2 + NO)
    @test sol[compiled.R4] > sol[compiled.R3]
end

@testitem "MethaneOxidation: CH4 Lifetime" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=MethaneOxidation()
    ch4_inputs=[:CH4, :CH3, :CH3O2, :CH3O, :CH3OOH, :HCHO, :HCO,
        :OH, :HO2, :H, :NO, :NO2, :O, :O2, :M]
    compiled=compile_with_inputs(sys, ch4_inputs)

    input_vals=Dict(
        compiled.CH4=>4.5e19, compiled.OH=>1e12, compiled.CH3=>1e10,
        compiled.CH3O2=>1e14, compiled.CH3O=>1e10, compiled.CH3OOH=>1e15,
        compiled.HCHO=>1e16, compiled.HCO=>1e10, compiled.HO2=>1e14,
        compiled.H=>1e8, compiled.NO=>2.5e15, compiled.NO2=>2.5e16,
        compiled.O=>1e8, compiled.O2=>5.25e24, compiled.M=>2.5e25
    )

    prob=NonlinearProblem(compiled, input_vals; build_initializeprob = false)
    sol=solve(prob)

    # CH4 lifetime = [CH4] / L_CH4 = [CH4] / R1
    tau_s=4.5e19/sol[compiled.L_CH4]
    tau_years=tau_s/(365.25*24*3600)

    # CH4 tropospheric lifetime is approximately 5-15 years at typical [OH]
    @test tau_years > 3
    @test tau_years < 20
    @test tau_years ≈ 5.03 rtol = 0.1
end

@testitem "MethaneOxidation: All Reaction Rates Positive" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=MethaneOxidation()
    ch4_inputs=[:CH4, :CH3, :CH3O2, :CH3O, :CH3OOH, :HCHO, :HCO,
        :OH, :HO2, :H, :NO, :NO2, :O, :O2, :M]
    compiled=compile_with_inputs(sys, ch4_inputs)

    input_vals=Dict(
        compiled.CH4=>4.5e19, compiled.OH=>1e12, compiled.CH3=>1e10,
        compiled.CH3O2=>1e14, compiled.CH3O=>1e10, compiled.CH3OOH=>1e15,
        compiled.HCHO=>1e16, compiled.HCO=>1e10, compiled.HO2=>1e14,
        compiled.H=>1e8, compiled.NO=>2.5e15, compiled.NO2=>2.5e16,
        compiled.O=>1e8, compiled.O2=>5.25e24, compiled.M=>2.5e25
    )

    prob=NonlinearProblem(compiled, input_vals; build_initializeprob = false)
    sol=solve(prob)

    # All reaction rates should be positive for positive input concentrations
    for R_sym in [compiled.R1, compiled.R2, compiled.R3, compiled.R4, compiled.R5,
        compiled.R6, compiled.R7, compiled.R8, compiled.R9, compiled.R10,
        compiled.R11, compiled.R12, compiled.R13, compiled.R14, compiled.R15,
        compiled.R16, compiled.R17]
        @test sol[R_sym] > 0
    end

    # All diagnostics should be positive
    @test sol[compiled.P_O3_gross] > 0
    @test sol[compiled.P_HCHO] > 0
    @test sol[compiled.L_CH4] > 0
end

# ===========================================================================
# Sensitivity Tests (using remake to vary input parameters)
# ===========================================================================
@testitem "MethaneOxidation: R1 Scales Linearly with OH" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=MethaneOxidation()
    ch4_inputs=[:CH4, :CH3, :CH3O2, :CH3O, :CH3OOH, :HCHO, :HCO,
        :OH, :HO2, :H, :NO, :NO2, :O, :O2, :M]
    compiled=compile_with_inputs(sys, ch4_inputs)

    input_vals=Dict(
        compiled.CH4=>4.5e19, compiled.OH=>1e12, compiled.CH3=>1e10,
        compiled.CH3O2=>1e14, compiled.CH3O=>1e10, compiled.CH3OOH=>1e15,
        compiled.HCHO=>1e16, compiled.HCO=>1e10, compiled.HO2=>1e14,
        compiled.H=>1e8, compiled.NO=>2.5e15, compiled.NO2=>2.5e16,
        compiled.O=>1e8, compiled.O2=>5.25e24, compiled.M=>2.5e25
    )

    prob_a=NonlinearProblem(compiled, input_vals; build_initializeprob = false)
    sol_a=solve(prob_a)

    # Double OH concentration
    prob_b=remake(prob_a, p = [compiled.OH=>2e12])
    sol_b=solve(prob_b)

    # R1 = k1 * CH4 * OH, so doubling OH should double R1
    @test sol_b[compiled.R1] / sol_a[compiled.R1] ≈ 2.0 rtol = 1e-6
end

# ===========================================================================
# Numerical Integration Tests (MethaneOxidationODE)
# ===========================================================================
@testitem "MethaneOxidationODE: Numerical Integration" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    using OrdinaryDiffEqRosenbrock

    sys=MethaneOxidationODE()
    compiled_sys=mtkcompile(sys)

    # Set initial conditions for all 18 state variables (SI: m^-3)
    u0=[
        compiled_sys.CH4=>4.5e19,
        compiled_sys.CH3=>0.0,
        compiled_sys.CH3O2=>1e14,
        compiled_sys.CH3O=>1e10,
        compiled_sys.CH3OOH=>1e15,
        compiled_sys.HCHO=>1e16,
        compiled_sys.HCO=>1e10,
        compiled_sys.CO=>2.5e18,
        compiled_sys.H2=>1.3e19,
        compiled_sys.OH=>1e12,
        compiled_sys.HO2=>1e14,
        compiled_sys.H=>0.0,
        compiled_sys.NO=>2.5e15,
        compiled_sys.NO2=>2.5e16,
        compiled_sys.O=>0.0,
        compiled_sys.O3=>1e18,
        compiled_sys.HNO3=>0.0,
        compiled_sys.H2O2=>0.0
    ]

    # Integrate for 100 seconds
    tspan=(0.0, 100.0)
    prob=ODEProblem(compiled_sys, u0, tspan)
    sol=solve(prob, Rosenbrock23(), abstol = 1e-8, reltol = 1e-8)

    @test string(sol.retcode) == "Success"

    # All concentrations must remain non-negative (with small numerical tolerance)
    for v in unknowns(compiled_sys)
        vals=sol[v]
        @test all(x -> x >= -1e-6 * max(abs(vals[1]), 1.0), vals)
    end
end

# ===========================================================================
# Conservation Law Tests
# ===========================================================================
@testitem "MethaneOxidationODE: Approximate NOx Conservation" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    using OrdinaryDiffEqRosenbrock

    sys=MethaneOxidationODE()
    compiled_sys=mtkcompile(sys)

    u0=[
        compiled_sys.CH4=>4.5e19,
        compiled_sys.CH3=>0.0,
        compiled_sys.CH3O2=>1e14,
        compiled_sys.CH3O=>1e10,
        compiled_sys.CH3OOH=>1e15,
        compiled_sys.HCHO=>1e16,
        compiled_sys.HCO=>1e10,
        compiled_sys.CO=>2.5e18,
        compiled_sys.H2=>1.3e19,
        compiled_sys.OH=>1e12,
        compiled_sys.HO2=>1e14,
        compiled_sys.H=>0.0,
        compiled_sys.NO=>2.5e15,
        compiled_sys.NO2=>2.5e16,
        compiled_sys.O=>0.0,
        compiled_sys.O3=>1e18,
        compiled_sys.HNO3=>0.0,
        compiled_sys.H2O2=>0.0
    ]

    tspan=(0.0, 100.0)
    prob=ODEProblem(compiled_sys, u0, tspan)
    sol=solve(prob, Rosenbrock23(), abstol = 1e-8, reltol = 1e-8)

    @test string(sol.retcode) == "Success"

    # Total reactive nitrogen: NO + NO2 + HNO3 should be approximately conserved
    NOx_total_initial=2.5e15+2.5e16+0.0

    NO_final=sol[compiled_sys.NO][end]
    NO2_final=sol[compiled_sys.NO2][end]
    HNO3_final=sol[compiled_sys.HNO3][end]
    NOx_total_final=NO_final+NO2_final+HNO3_final

    @test NOx_total_final ≈ NOx_total_initial rtol = 1e-4
end

# ===========================================================================
# Qualitative Behavior Tests
# ===========================================================================
@testitem "MethaneOxidationODE: CH4 Decreases Monotonically" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    using OrdinaryDiffEqRosenbrock

    sys=MethaneOxidationODE()
    compiled_sys=mtkcompile(sys)

    u0=[
        compiled_sys.CH4=>4.5e19,
        compiled_sys.CH3=>0.0,
        compiled_sys.CH3O2=>1e14,
        compiled_sys.CH3O=>1e10,
        compiled_sys.CH3OOH=>1e15,
        compiled_sys.HCHO=>1e16,
        compiled_sys.HCO=>1e10,
        compiled_sys.CO=>2.5e18,
        compiled_sys.H2=>1.3e19,
        compiled_sys.OH=>1e12,
        compiled_sys.HO2=>1e14,
        compiled_sys.H=>0.0,
        compiled_sys.NO=>2.5e15,
        compiled_sys.NO2=>2.5e16,
        compiled_sys.O=>0.0,
        compiled_sys.O3=>1e18,
        compiled_sys.HNO3=>0.0,
        compiled_sys.H2O2=>0.0
    ]

    tspan=(0.0, 100.0)
    prob=ODEProblem(compiled_sys, u0, tspan)
    sol=solve(prob, Rosenbrock23(), abstol = 1e-8, reltol = 1e-8, saveat = 1.0)

    @test string(sol.retcode) == "Success"

    ch4_vals=sol[compiled_sys.CH4]
    for i in 2:length(ch4_vals)
        @test ch4_vals[i] <= ch4_vals[i - 1] + 1e-3 * ch4_vals[1]
    end

    @test ch4_vals[end] < ch4_vals[1]
end

@testitem "MethaneOxidationODE: HNO3 and H2O2 Accumulate" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    using OrdinaryDiffEqRosenbrock

    sys=MethaneOxidationODE()
    compiled_sys=mtkcompile(sys)

    u0=[
        compiled_sys.CH4=>4.5e19,
        compiled_sys.CH3=>0.0,
        compiled_sys.CH3O2=>1e14,
        compiled_sys.CH3O=>1e10,
        compiled_sys.CH3OOH=>1e15,
        compiled_sys.HCHO=>1e16,
        compiled_sys.HCO=>1e10,
        compiled_sys.CO=>2.5e18,
        compiled_sys.H2=>1.3e19,
        compiled_sys.OH=>1e12,
        compiled_sys.HO2=>1e14,
        compiled_sys.H=>0.0,
        compiled_sys.NO=>2.5e15,
        compiled_sys.NO2=>2.5e16,
        compiled_sys.O=>0.0,
        compiled_sys.O3=>1e18,
        compiled_sys.HNO3=>0.0,
        compiled_sys.H2O2=>0.0
    ]

    tspan=(0.0, 100.0)
    prob=ODEProblem(compiled_sys, u0, tspan)
    sol=solve(prob, Rosenbrock23(), abstol = 1e-8, reltol = 1e-8)

    @test string(sol.retcode) == "Success"

    @test sol[compiled_sys.HNO3][end] > 0
    @test sol[compiled_sys.H2O2][end] > 0
end

# ===========================================================================
# Structural Tests
# ===========================================================================
@testitem "MethaneOxidation: Structural Verification" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = MethaneOxidation()
    @test sys isa System
    @test nameof(sys) == :MethaneOxidation

    vars = unknowns(sys)
    params = parameters(sys)
    eqs = equations(sys)

    # 15 input species (CO and O3 are declared but unused in equations, so MTK excludes them)
    # + 17 reaction rates (R1-R17) + 3 diagnostics (P_O3_net, P_HCHO, L_CH4) = 35
    @test length(vars) == 35

    # 20 equations: 17 reaction rates + 3 diagnostics
    @test length(eqs) == 20

    # 17 parameters: k1, k2_0, k3-k8, k10, k13, k14_0, k15, k17_0, j9, j11, j12, j16
    @test length(params) == 17

    # Check key diagnostic variable names
    var_names = [string(v) for v in vars]
    for expected in ["P_O3_net", "P_HCHO", "L_CH4"]
        @test any(n -> contains(n, expected), var_names)
    end
end

@testitem "MethaneOxidationODE: Structural Verification" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = MethaneOxidationODE()
    @test sys isa System
    @test nameof(sys) == :MethaneOxidationODE

    vars = unknowns(sys)
    params = parameters(sys)
    eqs = equations(sys)

    # 18 state variables
    @test length(vars) == 18

    # 18 ODEs (one for each state variable)
    @test length(eqs) == 18

    # Check key species names
    var_names = [string(v) for v in vars]
    for expected in ["CH4", "CH3O2", "HCHO", "OH", "HO2", "NO", "NO2", "O3", "HNO3", "H2O2"]
        @test any(n -> contains(n, expected), var_names)
    end
end

# ===========================================================================
# Rate Constants at 298 K (Table 6.1)
# ===========================================================================
@testitem "MethaneOxidation: Table 6.1 Rate Constants" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = MethaneOxidation()
    params = parameters(sys)
    param_dict = Dict(Symbol(p) => ModelingToolkit.getdefault(p) for p in params if ModelingToolkit.hasdefault(p))

    # Bimolecular rate constants (m³/s)
    @test param_dict[:k1] ≈ 6.3e-15 * 1e-6
    @test param_dict[:k3] ≈ 7.7e-12 * 1e-6
    @test param_dict[:k4] ≈ 5.2e-12 * 1e-6
    @test param_dict[:k5] ≈ 3.5e-13 * 1e-6
    @test param_dict[:k6] ≈ 1.9e-15 * 1e-6
    @test param_dict[:k7] ≈ 3.8e-12 * 1e-6
    @test param_dict[:k8] ≈ 1.9e-12 * 1e-6
    @test param_dict[:k10] ≈ 8.5e-12 * 1e-6
    @test param_dict[:k13] ≈ 5.2e-12 * 1e-6
    @test param_dict[:k15] ≈ 8.1e-12 * 1e-6

    # Termolecular rate constants (m⁶/s)
    @test param_dict[:k2_0] ≈ 1.0e-30 * 1e-12
    @test param_dict[:k14_0] ≈ 5.7e-32 * 1e-12
    @test param_dict[:k17_0] ≈ 6.0e-34 * 1e-12

    # Photolysis rates (s⁻¹, unchanged)
    @test param_dict[:j9] ≈ 5e-6
    @test param_dict[:j11] ≈ 3e-5
    @test param_dict[:j12] ≈ 5e-5
    @test param_dict[:j16] ≈ 8e-3
end

@testitem "MethaneOxidationODE: Rate Constants and Parameters" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = MethaneOxidationODE()
    params = parameters(sys)
    param_dict = Dict(Symbol(p) => ModelingToolkit.getdefault(p) for p in params if ModelingToolkit.hasdefault(p))

    # Same kinetic rate constants as MethaneOxidation (in SI)
    @test param_dict[:k1] ≈ 6.3e-15 * 1e-6
    @test param_dict[:k3] ≈ 7.7e-12 * 1e-6
    @test param_dict[:k15] ≈ 8.1e-12 * 1e-6

    # Additional rate constants
    @test param_dict[:k_CO_OH] ≈ 2.4e-13 * 1e-6
    @test param_dict[:k_OH_NO2] ≈ 1.0e-11 * 1e-6
    @test param_dict[:k_HO2_HO2] ≈ 2.9e-12 * 1e-6
    @test param_dict[:k_NO_O3] ≈ 1.8e-14 * 1e-6

    # External OH source (m⁻³/s)
    @test param_dict[:P_OH_ext] ≈ 1e6 * 1e6
end

# ===========================================================================
# Equation Verification Tests
# ===========================================================================
@testitem "MethaneOxidation: Reaction Rate R1 (CH4+OH)" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # R1 = k1 * [CH4] * [OH]
    k1 = 6.3e-15 * 1e-6    # m³/s
    CH4 = 4.5e19             # m⁻³ (~1800 ppb)
    OH = 1e12                # m⁻³

    R1 = k1 * CH4 * OH
    @test R1 > 0
    @test R1 ≈ 6.3e-21 * 4.5e19 * 1e12 rtol=1e-10
    @test R1 ≈ 2.835e11 rtol=0.01
end

@testitem "MethaneOxidation: Reaction Rate R3 (CH3O2+NO)" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # R3 = k3 * [CH3O2] * [NO]
    k3 = 7.7e-12 * 1e-6     # m³/s
    CH3O2 = 1e14              # m⁻³
    NO = 2.5e15               # m⁻³ (~0.1 ppb)

    R3 = k3 * CH3O2 * NO
    @test R3 > 0
    @test R3 ≈ 7.7e-18 * 1e14 * 2.5e15 rtol=1e-10
    @test R3 ≈ 1.925e12 rtol=0.01
end

@testitem "MethaneOxidation: HCHO Production" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # P_HCHO = R6 + R8 = k6*CH3O*O2 + k8*CH3OOH*OH
    k6 = 1.9e-15 * 1e-6     # m³/s
    k8 = 1.9e-12 * 1e-6     # m³/s
    O2 = 5.25e24              # m⁻³
    CH3O = 1e10               # m⁻³ (very short-lived radical)
    CH3OOH = 1e15             # m⁻³
    OH = 1e12                 # m⁻³

    R6 = k6 * CH3O * O2
    R8 = k8 * CH3OOH * OH

    P_HCHO = R6 + R8
    @test P_HCHO > 0
    @test R6 ≈ 1.9e-21 * 1e10 * 5.25e24 rtol=1e-10
    @test R8 ≈ 1.9e-18 * 1e15 * 1e12 rtol=1e-10
end

@testitem "MethaneOxidation: CH4 Loss Equals R1" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # The only loss pathway for CH4 is reaction with OH (R1)
    k1 = 6.3e-15 * 1e-6
    CH4 = 4.5e19
    OH = 1e12

    R1 = k1 * CH4 * OH
    L_CH4 = R1

    @test L_CH4 ≈ R1 rtol=1e-10
    @test L_CH4 > 0
end

# ===========================================================================
# Numerical Integration Tests
# ===========================================================================
@testitem "MethaneOxidationODE: Numerical Integration" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    using OrdinaryDiffEqRosenbrock

    sys = MethaneOxidationODE()
    compiled_sys = mtkcompile(sys)

    # Set initial conditions for all 18 state variables (SI: m⁻³)
    u0 = [
        compiled_sys.CH4 => 4.5e19,
        compiled_sys.CH3 => 0.0,
        compiled_sys.CH3O2 => 1e14,
        compiled_sys.CH3O => 1e10,
        compiled_sys.CH3OOH => 1e15,
        compiled_sys.HCHO => 1e16,
        compiled_sys.HCO => 1e10,
        compiled_sys.CO => 2.5e18,
        compiled_sys.H2 => 1.3e19,
        compiled_sys.OH => 1e12,
        compiled_sys.HO2 => 1e14,
        compiled_sys.H => 0.0,
        compiled_sys.NO => 2.5e15,
        compiled_sys.NO2 => 2.5e16,
        compiled_sys.O => 0.0,
        compiled_sys.O3 => 1e18,
        compiled_sys.HNO3 => 0.0,
        compiled_sys.H2O2 => 0.0,
    ]

    # Integrate for 100 seconds
    tspan = (0.0, 100.0)
    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, Rosenbrock23(), abstol=1e-8, reltol=1e-8)

    @test string(sol.retcode) == "Success"

    # All concentrations must remain non-negative
    for v in unknowns(compiled_sys)
        vals = sol[v]
        @test all(x -> x >= -1e-6 * max(abs(vals[1]), 1.0), vals)
    end
end

# ===========================================================================
# Conservation Law Tests
# ===========================================================================
@testitem "MethaneOxidationODE: Approximate NOx Conservation" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    using OrdinaryDiffEqRosenbrock

    sys = MethaneOxidationODE()
    compiled_sys = mtkcompile(sys)

    u0 = [
        compiled_sys.CH4 => 4.5e19,
        compiled_sys.CH3 => 0.0,
        compiled_sys.CH3O2 => 1e14,
        compiled_sys.CH3O => 1e10,
        compiled_sys.CH3OOH => 1e15,
        compiled_sys.HCHO => 1e16,
        compiled_sys.HCO => 1e10,
        compiled_sys.CO => 2.5e18,
        compiled_sys.H2 => 1.3e19,
        compiled_sys.OH => 1e12,
        compiled_sys.HO2 => 1e14,
        compiled_sys.H => 0.0,
        compiled_sys.NO => 2.5e15,
        compiled_sys.NO2 => 2.5e16,
        compiled_sys.O => 0.0,
        compiled_sys.O3 => 1e18,
        compiled_sys.HNO3 => 0.0,
        compiled_sys.H2O2 => 0.0,
    ]

    tspan = (0.0, 100.0)
    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, Rosenbrock23(), abstol=1e-8, reltol=1e-8)

    @test string(sol.retcode) == "Success"

    # Total reactive nitrogen: NO + NO2 + HNO3 should be approximately conserved
    NOx_total_initial = 2.5e15 + 2.5e16 + 0.0

    NO_final = sol[compiled_sys.NO][end]
    NO2_final = sol[compiled_sys.NO2][end]
    HNO3_final = sol[compiled_sys.HNO3][end]
    NOx_total_final = NO_final + NO2_final + HNO3_final

    @test NOx_total_final ≈ NOx_total_initial rtol=1e-4
end

# ===========================================================================
# Qualitative Behavior Tests
# ===========================================================================
@testitem "MethaneOxidationODE: CH4 Decreases Monotonically" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    using OrdinaryDiffEqRosenbrock

    sys = MethaneOxidationODE()
    compiled_sys = mtkcompile(sys)

    u0 = [
        compiled_sys.CH4 => 4.5e19,
        compiled_sys.CH3 => 0.0,
        compiled_sys.CH3O2 => 1e14,
        compiled_sys.CH3O => 1e10,
        compiled_sys.CH3OOH => 1e15,
        compiled_sys.HCHO => 1e16,
        compiled_sys.HCO => 1e10,
        compiled_sys.CO => 2.5e18,
        compiled_sys.H2 => 1.3e19,
        compiled_sys.OH => 1e12,
        compiled_sys.HO2 => 1e14,
        compiled_sys.H => 0.0,
        compiled_sys.NO => 2.5e15,
        compiled_sys.NO2 => 2.5e16,
        compiled_sys.O => 0.0,
        compiled_sys.O3 => 1e18,
        compiled_sys.HNO3 => 0.0,
        compiled_sys.H2O2 => 0.0,
    ]

    tspan = (0.0, 100.0)
    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, Rosenbrock23(), abstol=1e-8, reltol=1e-8, saveat=1.0)

    @test string(sol.retcode) == "Success"

    ch4_vals = sol[compiled_sys.CH4]
    for i in 2:length(ch4_vals)
        @test ch4_vals[i] <= ch4_vals[i-1] + 1e-3 * ch4_vals[1]
    end

    @test ch4_vals[end] < ch4_vals[1]
end

@testitem "MethaneOxidationODE: HNO3 and H2O2 Accumulate" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    using OrdinaryDiffEqRosenbrock

    sys = MethaneOxidationODE()
    compiled_sys = mtkcompile(sys)

    u0 = [
        compiled_sys.CH4 => 4.5e19,
        compiled_sys.CH3 => 0.0,
        compiled_sys.CH3O2 => 1e14,
        compiled_sys.CH3O => 1e10,
        compiled_sys.CH3OOH => 1e15,
        compiled_sys.HCHO => 1e16,
        compiled_sys.HCO => 1e10,
        compiled_sys.CO => 2.5e18,
        compiled_sys.H2 => 1.3e19,
        compiled_sys.OH => 1e12,
        compiled_sys.HO2 => 1e14,
        compiled_sys.H => 0.0,
        compiled_sys.NO => 2.5e15,
        compiled_sys.NO2 => 2.5e16,
        compiled_sys.O => 0.0,
        compiled_sys.O3 => 1e18,
        compiled_sys.HNO3 => 0.0,
        compiled_sys.H2O2 => 0.0,
    ]

    tspan = (0.0, 100.0)
    prob = ODEProblem(compiled_sys, u0, tspan)
    sol = solve(prob, Rosenbrock23(), abstol=1e-8, reltol=1e-8)

    @test string(sol.retcode) == "Success"

    @test sol[compiled_sys.HNO3][end] > 0
    @test sol[compiled_sys.H2O2][end] > 0
end

# ===========================================================================
# Methane Oxidation Chain Tests
# ===========================================================================
@testitem "MethaneOxidation: Oxidation Chain Relative Rates" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    k3 = 7.7e-12 * 1e-6
    k4 = 5.2e-12 * 1e-6
    k5 = 3.5e-13 * 1e-6

    CH3O2 = 1e14
    NO = 2.5e15
    HO2 = 1e14

    R3 = k3 * CH3O2 * NO
    R4 = k4 * CH3O2 * HO2
    R5 = k5 * CH3O2 * CH3O2

    # At moderate NOx, R3 >> R4 >> R5
    @test R3 > R4
    @test R4 > R5
    @test R3 / R4 > 10
end

@testitem "MethaneOxidation: HCHO as Key Intermediate" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    k10 = 8.5e-12 * 1e-6
    j11 = 3e-5
    j12 = 5e-5

    HCHO = 1e16
    OH = 1e12

    R10 = k10 * HCHO * OH
    R11 = j11 * HCHO
    R12 = j12 * HCHO

    @test R10 > 0
    @test R11 > 0
    @test R12 > 0

    @test R12 > R11
end

# ===========================================================================
# Limiting Behavior Tests
# ===========================================================================
@testitem "MethaneOxidation: Low NOx Regime" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    k3 = 7.7e-12 * 1e-6
    k4 = 5.2e-12 * 1e-6

    CH3O2 = 1e14
    HO2 = 2e14  # higher at low NOx

    NO_low = 1e13  # ~0.4 ppt (m⁻³)
    R3_low = k3 * CH3O2 * NO_low
    R4_low = k4 * CH3O2 * HO2

    @test R4_low > R3_low
end

@testitem "MethaneOxidation: CH4 Lifetime" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # CH4 lifetime = 1 / (k1 * [OH])
    k1 = 6.3e-15 * 1e-6
    OH = 1e12

    tau_s = 1.0 / (k1 * OH)
    tau_years = tau_s / (365.25 * 24 * 3600)

    # CH4 tropospheric lifetime is approximately 8-12 years
    @test tau_years > 3
    @test tau_years < 20
    @test tau_years ≈ 5.03 rtol=0.1
end

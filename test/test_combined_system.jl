# ===========================================================================
# Structural Tests
# ===========================================================================
@testitem "TroposphericChemistrySystem: Structural Verification" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=TroposphericChemistrySystem()
    @test sys isa System
    @test nameof(sys) == :TroposphericChemistrySystem

    vars=unknowns(sys)
    params=parameters(sys)
    eqs=equations(sys)

    @test length(eqs) > 0
    @test length(vars) > 0
    @test length(params) > 0

    # Check key diagnostic variable names at the top level
    var_names=[string(v) for v in vars]
    for expected in ["NOx", "HOx", "P_O3_net", "OPE", "chain_length", "L_NOx"]
        @test any(n -> contains(n, expected), var_names)
    end
end

@testitem "TroposphericChemistrySystem: Subsystem Composition" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=TroposphericChemistrySystem()

    subsystems=ModelingToolkit.get_systems(sys)
    subsys_names=[nameof(s) for s in subsystems]

    @test :oh in subsys_names
    @test :nox in subsys_names
    @test :co in subsys_names
    @test length(subsystems) == 3
end

@testitem "TroposphericChemistrySystem: Subsystem Variable Access" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=TroposphericChemistrySystem()

    all_vars=unknowns(sys)
    var_strs=[string(v) for v in all_vars]

    @test any(s -> contains(s, "oh") && contains(s, "O1D"), var_strs)
    @test any(s -> contains(s, "oh") && contains(s, "P_OH"), var_strs)
    @test any(s -> contains(s, "nox") && contains(s, "O3_pss"), var_strs)
    @test any(s -> contains(s, "co") && contains(s, "HO2_ss"), var_strs)
end

# ===========================================================================
# Condition Functions Tests
# ===========================================================================
@testitem "TroposphericChemistrySystem: TypicalConditions" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    cond=get_conditions_dict(TypicalConditions())

    required_keys=[:M, :O2, :H2O, :O3, :NO, :NO2, :CO, :CH4, :OH, :HO2, :CH3O2]
    for key in required_keys
        @test haskey(cond, key)
        @test cond[key] > 0
    end

    # Verify specific values (SI: m^-3)
    @test cond[:M] ≈ 2.5e25
    @test cond[:O2] ≈ 5.25e24
    @test cond[:O3] ≈ 1e18
    @test cond[:OH] ≈ 1e12
    @test cond[:HO2] ≈ 1e14
    @test cond[:CO] ≈ 2.5e18
    @test cond[:CH4] ≈ 4.5e19

    # Physical consistency: O2 is ~21% of M
    @test cond[:O2] / cond[:M] ≈ 0.21 rtol = 0.01
end

@testitem "TroposphericChemistrySystem: UrbanConditions vs TypicalConditions" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    typical=get_conditions_dict(TypicalConditions())
    urban=get_conditions_dict(UrbanConditions())

    @test urban[:NO] > typical[:NO]
    @test urban[:NO2] > typical[:NO2]
    @test urban[:CO] > typical[:CO]
    @test urban[:OH] <= typical[:OH]
    @test urban[:M] ≈ typical[:M]
    @test urban[:CH4] ≈ typical[:CH4]
end

@testitem "TroposphericChemistrySystem: RemoteConditions vs TypicalConditions" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    typical=get_conditions_dict(TypicalConditions())
    remote=get_conditions_dict(RemoteConditions())

    @test remote[:NO] < typical[:NO]
    @test remote[:NO2] < typical[:NO2]
    @test remote[:HO2] > typical[:HO2]
    @test remote[:M] ≈ typical[:M]
end

# ===========================================================================
# Equation Verification: Combined Diagnostics (using NonlinearProblem)
# ===========================================================================
@testitem "TroposphericChemistrySystem: NOx = NO + NO2" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=TroposphericChemistrySystem()
    sys_nns=toggle_namespacing(sys, false)
    inputs=[sys_nns.O3, sys_nns.NO, sys_nns.NO2, sys_nns.OH, sys_nns.HO2,
        sys_nns.CO, sys_nns.CH3O2, sys_nns.H2O, sys_nns.M, sys_nns.O2]
    compiled=mtkcompile(sys; inputs = inputs)

    cond=get_conditions_dict(TypicalConditions())
    u0p=Dict(
        compiled.O3=>cond[:O3], compiled.NO=>cond[:NO], compiled.NO2=>cond[:NO2],
        compiled.OH=>cond[:OH], compiled.HO2=>cond[:HO2], compiled.CO=>cond[:CO],
        compiled.CH3O2=>cond[:CH3O2], compiled.H2O=>cond[:H2O],
        compiled.M=>cond[:M], compiled.O2=>cond[:O2]
    )
    prob=NonlinearProblem(compiled, u0p; build_initializeprob = false)
    sol=solve(prob)

    # NOx should equal NO + NO2
    @test sol[compiled.NOx] ≈ cond[:NO] + cond[:NO2] rtol = 1e-6
end

@testitem "TroposphericChemistrySystem: HOx = OH + HO2" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=TroposphericChemistrySystem()
    sys_nns=toggle_namespacing(sys, false)
    inputs=[sys_nns.O3, sys_nns.NO, sys_nns.NO2, sys_nns.OH, sys_nns.HO2,
        sys_nns.CO, sys_nns.CH3O2, sys_nns.H2O, sys_nns.M, sys_nns.O2]
    compiled=mtkcompile(sys; inputs = inputs)

    cond=get_conditions_dict(TypicalConditions())
    u0p=Dict(
        compiled.O3=>cond[:O3], compiled.NO=>cond[:NO], compiled.NO2=>cond[:NO2],
        compiled.OH=>cond[:OH], compiled.HO2=>cond[:HO2], compiled.CO=>cond[:CO],
        compiled.CH3O2=>cond[:CH3O2], compiled.H2O=>cond[:H2O],
        compiled.M=>cond[:M], compiled.O2=>cond[:O2]
    )
    prob=NonlinearProblem(compiled, u0p; build_initializeprob = false)
    sol=solve(prob)

    # HOx should equal OH + HO2
    @test sol[compiled.HOx] ≈ cond[:OH] + cond[:HO2] rtol = 1e-6

    # HOx is dominated by HO2
    @test sol[compiled.HOx] ≈ cond[:HO2] rtol = 0.01
    @test cond[:HO2] / cond[:OH] > 10
end

@testitem "TroposphericChemistrySystem: O3 Production Diagnostics" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=TroposphericChemistrySystem()
    sys_nns=toggle_namespacing(sys, false)
    inputs=[sys_nns.O3, sys_nns.NO, sys_nns.NO2, sys_nns.OH, sys_nns.HO2,
        sys_nns.CO, sys_nns.CH3O2, sys_nns.H2O, sys_nns.M, sys_nns.O2]
    compiled=mtkcompile(sys; inputs = inputs)

    cond=get_conditions_dict(TypicalConditions())
    u0p=Dict(
        compiled.O3=>cond[:O3], compiled.NO=>cond[:NO], compiled.NO2=>cond[:NO2],
        compiled.OH=>cond[:OH], compiled.HO2=>cond[:HO2], compiled.CO=>cond[:CO],
        compiled.CH3O2=>cond[:CH3O2], compiled.H2O=>cond[:H2O],
        compiled.M=>cond[:M], compiled.O2=>cond[:O2]
    )
    prob=NonlinearProblem(compiled, u0p; build_initializeprob = false)
    sol=solve(prob)

    P_O3=sol[compiled.P_O3_total]
    @test P_O3 > 0
    @test P_O3 > 1e11
    @test P_O3 < 1e14

    # Verify the budget is consistent: P_O3_net = P_O3_total - L_O3_total
    @test sol[compiled.P_O3_net] ≈ sol[compiled.P_O3_total] - sol[compiled.L_O3_total] rtol = 1e-6
end

@testitem "TroposphericChemistrySystem: NOx Loss Rate" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=TroposphericChemistrySystem()
    sys_nns=toggle_namespacing(sys, false)
    inputs=[sys_nns.O3, sys_nns.NO, sys_nns.NO2, sys_nns.OH, sys_nns.HO2,
        sys_nns.CO, sys_nns.CH3O2, sys_nns.H2O, sys_nns.M, sys_nns.O2]
    compiled=mtkcompile(sys; inputs = inputs)

    cond=get_conditions_dict(TypicalConditions())
    u0p=Dict(
        compiled.O3=>cond[:O3], compiled.NO=>cond[:NO], compiled.NO2=>cond[:NO2],
        compiled.OH=>cond[:OH], compiled.HO2=>cond[:HO2], compiled.CO=>cond[:CO],
        compiled.CH3O2=>cond[:CH3O2], compiled.H2O=>cond[:H2O],
        compiled.M=>cond[:M], compiled.O2=>cond[:O2]
    )
    prob=NonlinearProblem(compiled, u0p; build_initializeprob = false)
    sol=solve(prob)

    L_NOx_val=sol[compiled.L_NOx]
    @test L_NOx_val > 0

    # L_NOx = k_OH_NO2 * OH * NO2 (from the co subsystem rate constant)
    # k_OH_NO2 = 1.0e-11 * 1e-6 = 1.0e-17 m^3/s
    # L_NOx = 1.0e-17 * 1e12 * 2.5e16 = 2.5e11 m^-3/s
    @test L_NOx_val ≈ 1.0e-17 * 1e12 * 2.5e16 rtol = 1e-6
end

@testitem "TroposphericChemistrySystem: OPE Calculation" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=TroposphericChemistrySystem()
    sys_nns=toggle_namespacing(sys, false)
    inputs=[sys_nns.O3, sys_nns.NO, sys_nns.NO2, sys_nns.OH, sys_nns.HO2,
        sys_nns.CO, sys_nns.CH3O2, sys_nns.H2O, sys_nns.M, sys_nns.O2]
    compiled=mtkcompile(sys; inputs = inputs)

    cond=get_conditions_dict(TypicalConditions())
    u0p=Dict(
        compiled.O3=>cond[:O3], compiled.NO=>cond[:NO], compiled.NO2=>cond[:NO2],
        compiled.OH=>cond[:OH], compiled.HO2=>cond[:HO2], compiled.CO=>cond[:CO],
        compiled.CH3O2=>cond[:CH3O2], compiled.H2O=>cond[:H2O],
        compiled.M=>cond[:M], compiled.O2=>cond[:O2]
    )
    prob=NonlinearProblem(compiled, u0p; build_initializeprob = false)
    sol=solve(prob)

    OPE_val=sol[compiled.OPE]
    @test OPE_val > 0
    @test OPE_val > 1
    @test OPE_val < 100

    # OPE should be consistent: P_O3_total / L_NOx
    @test OPE_val ≈ sol[compiled.P_O3_total] / sol[compiled.L_NOx] rtol = 1e-6
end

# ===========================================================================
# Subsystem Coupling Tests (using NonlinearProblem)
# ===========================================================================
@testitem "TroposphericChemistrySystem: OH Production Coupling" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=TroposphericChemistrySystem()
    sys_nns=toggle_namespacing(sys, false)
    inputs=[sys_nns.O3, sys_nns.NO, sys_nns.NO2, sys_nns.OH, sys_nns.HO2,
        sys_nns.CO, sys_nns.CH3O2, sys_nns.H2O, sys_nns.M, sys_nns.O2]
    compiled=mtkcompile(sys; inputs = inputs)

    cond=get_conditions_dict(TypicalConditions())
    u0p=Dict(
        compiled.O3=>cond[:O3], compiled.NO=>cond[:NO], compiled.NO2=>cond[:NO2],
        compiled.OH=>cond[:OH], compiled.HO2=>cond[:HO2], compiled.CO=>cond[:CO],
        compiled.CH3O2=>cond[:CH3O2], compiled.H2O=>cond[:H2O],
        compiled.M=>cond[:M], compiled.O2=>cond[:O2]
    )
    prob=NonlinearProblem(compiled, u0p; build_initializeprob = false)
    sol=solve(prob)

    # OH production subsystem should produce valid outputs
    eps_val=sol[compiled.oh.ε_OH]
    @test eps_val > 0 && eps_val < 1

    P_OH_val=sol[compiled.oh.P_OH]
    @test P_OH_val > 0
    @test P_OH_val < 1e14

    O1D_val=sol[compiled.oh.O1D]
    @test O1D_val > 0
end

@testitem "TroposphericChemistrySystem: NOx Cycling Coupling" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=TroposphericChemistrySystem()
    sys_nns=toggle_namespacing(sys, false)
    inputs=[sys_nns.O3, sys_nns.NO, sys_nns.NO2, sys_nns.OH, sys_nns.HO2,
        sys_nns.CO, sys_nns.CH3O2, sys_nns.H2O, sys_nns.M, sys_nns.O2]
    compiled=mtkcompile(sys; inputs = inputs)

    cond=get_conditions_dict(TypicalConditions())
    u0p=Dict(
        compiled.O3=>cond[:O3], compiled.NO=>cond[:NO], compiled.NO2=>cond[:NO2],
        compiled.OH=>cond[:OH], compiled.HO2=>cond[:HO2], compiled.CO=>cond[:CO],
        compiled.CH3O2=>cond[:CH3O2], compiled.H2O=>cond[:H2O],
        compiled.M=>cond[:M], compiled.O2=>cond[:O2]
    )
    prob=NonlinearProblem(compiled, u0p; build_initializeprob = false)
    sol=solve(prob)

    # NOx subsystem outputs
    Phi_val=sol[compiled.nox.Φ]
    @test Phi_val > 0

    O3_pss_val=sol[compiled.nox.O3_pss]
    @test O3_pss_val > 0

    # The NOx subsystem net O3 production rate
    P_O3_nox=sol[compiled.nox.P_O3]
    # This can be positive or negative depending on conditions (Eq. 6.8)
    @test isfinite(P_O3_nox)
end

@testitem "TroposphericChemistrySystem: CO Oxidation Coupling" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=TroposphericChemistrySystem()
    sys_nns=toggle_namespacing(sys, false)
    inputs=[sys_nns.O3, sys_nns.NO, sys_nns.NO2, sys_nns.OH, sys_nns.HO2,
        sys_nns.CO, sys_nns.CH3O2, sys_nns.H2O, sys_nns.M, sys_nns.O2]
    compiled=mtkcompile(sys; inputs = inputs)

    cond=get_conditions_dict(TypicalConditions())
    u0p=Dict(
        compiled.O3=>cond[:O3], compiled.NO=>cond[:NO], compiled.NO2=>cond[:NO2],
        compiled.OH=>cond[:OH], compiled.HO2=>cond[:HO2], compiled.CO=>cond[:CO],
        compiled.CH3O2=>cond[:CH3O2], compiled.H2O=>cond[:H2O],
        compiled.M=>cond[:M], compiled.O2=>cond[:O2]
    )
    prob=NonlinearProblem(compiled, u0p; build_initializeprob = false)
    sol=solve(prob)

    # CO oxidation subsystem: net O3 production from CO cycle
    P_O3_co=sol[compiled.co.P_O3]
    @test P_O3_co > 0  # Should be positive in typical conditions

    # HOx loss rate should be positive
    L_HOx_val=sol[compiled.co.L_HOx]
    @test L_HOx_val > 0

    # Chain length from CO subsystem should match top-level chain_length
    @test sol[compiled.chain_length] ≈ sol[compiled.co.chain_length] rtol = 1e-6

    # Chain length should be > 1 (catalytic cycling)
    @test sol[compiled.chain_length] > 1
end

# ===========================================================================
# Regime Comparison Tests (using NonlinearProblem)
# ===========================================================================
@testitem "TroposphericChemistrySystem: NOx-VOC Regime Differences" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=TroposphericChemistrySystem()
    sys_nns=toggle_namespacing(sys, false)
    inputs=[sys_nns.O3, sys_nns.NO, sys_nns.NO2, sys_nns.OH, sys_nns.HO2,
        sys_nns.CO, sys_nns.CH3O2, sys_nns.H2O, sys_nns.M, sys_nns.O2]
    compiled=mtkcompile(sys; inputs = inputs)

    function solve_system(compiled, cond)
        u0p=Dict(
            compiled.O3=>cond[:O3], compiled.NO=>cond[:NO], compiled.NO2=>cond[:NO2],
            compiled.OH=>cond[:OH], compiled.HO2=>cond[:HO2], compiled.CO=>cond[:CO],
            compiled.CH3O2=>cond[:CH3O2], compiled.H2O=>cond[:H2O],
            compiled.M=>cond[:M], compiled.O2=>cond[:O2]
        )
        prob=NonlinearProblem(compiled, u0p; build_initializeprob = false)
        return solve(prob)
    end

    sol_typical=solve_system(compiled, get_conditions_dict(TypicalConditions()))
    sol_urban=solve_system(compiled, get_conditions_dict(UrbanConditions()))
    sol_remote=solve_system(compiled, get_conditions_dict(RemoteConditions()))

    # All OPE values should be positive
    @test sol_typical[compiled.OPE] > 0
    @test sol_urban[compiled.OPE] > 0
    @test sol_remote[compiled.OPE] > 0

    # All P_O3_total values should be positive
    @test sol_typical[compiled.P_O3_total] > 0
    @test sol_urban[compiled.P_O3_total] > 0
    @test sol_remote[compiled.P_O3_total] > 0
end

@testitem "TroposphericChemistrySystem: O3 Production Budget Consistency" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=TroposphericChemistrySystem()
    sys_nns=toggle_namespacing(sys, false)
    inputs=[sys_nns.O3, sys_nns.NO, sys_nns.NO2, sys_nns.OH, sys_nns.HO2,
        sys_nns.CO, sys_nns.CH3O2, sys_nns.H2O, sys_nns.M, sys_nns.O2]
    compiled=mtkcompile(sys; inputs = inputs)

    cond=get_conditions_dict(TypicalConditions())
    u0p=Dict(
        compiled.O3=>cond[:O3], compiled.NO=>cond[:NO], compiled.NO2=>cond[:NO2],
        compiled.OH=>cond[:OH], compiled.HO2=>cond[:HO2], compiled.CO=>cond[:CO],
        compiled.CH3O2=>cond[:CH3O2], compiled.H2O=>cond[:H2O],
        compiled.M=>cond[:M], compiled.O2=>cond[:O2]
    )
    prob=NonlinearProblem(compiled, u0p; build_initializeprob = false)
    sol=solve(prob)

    P_O3_total=sol[compiled.P_O3_total]
    L_O3_total=sol[compiled.L_O3_total]
    P_O3_net=sol[compiled.P_O3_net]

    @test P_O3_total > 0
    @test L_O3_total > 0
    @test P_O3_net ≈ P_O3_total - L_O3_total rtol = 1e-6
end

# ===========================================================================
# Qualitative Property Tests
# ===========================================================================
@testitem "TroposphericChemistrySystem: Physical Bounds" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    for cond_sys in [TypicalConditions(), UrbanConditions(), RemoteConditions()]
        cond=get_conditions_dict(cond_sys)
        for (key, val) in cond
            @test val > 0
        end
    end
end

@testitem "TroposphericChemistrySystem: Atmospheric Composition Consistency" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=TroposphericChemistrySystem()
    sys_nns=toggle_namespacing(sys, false)
    inputs=[sys_nns.O3, sys_nns.NO, sys_nns.NO2, sys_nns.OH, sys_nns.HO2,
        sys_nns.CO, sys_nns.CH3O2, sys_nns.H2O, sys_nns.M, sys_nns.O2]
    compiled=mtkcompile(sys; inputs = inputs)

    function solve_system(compiled, cond)
        u0p=Dict(
            compiled.O3=>cond[:O3], compiled.NO=>cond[:NO], compiled.NO2=>cond[:NO2],
            compiled.OH=>cond[:OH], compiled.HO2=>cond[:HO2], compiled.CO=>cond[:CO],
            compiled.CH3O2=>cond[:CH3O2], compiled.H2O=>cond[:H2O],
            compiled.M=>cond[:M], compiled.O2=>cond[:O2]
        )
        prob=NonlinearProblem(compiled, u0p; build_initializeprob = false)
        return solve(prob)
    end

    for (i,
        cond_sys) in enumerate([TypicalConditions(), UrbanConditions(), RemoteConditions()])
        cond=get_conditions_dict(cond_sys)

        # O2 should be ~21% of M
        @test cond[:O2] / cond[:M] ≈ 0.21 rtol = 0.01

        # O3 should be much less than O2
        @test cond[:O3] < cond[:O2] * 1e-6

        # NO2 should generally be larger than NO
        @test cond[:NO2] >= cond[:NO] || cond[:NO] < 1e15

        # HO2 >> OH
        @test cond[:HO2] / cond[:OH] > 1

        # CH4 is the most abundant reactive hydrocarbon (except urban CO)
        if i!=2  # not UrbanConditions
            @test cond[:CH4] > cond[:CO]
        end

        # Solve the system and verify all computed outputs are positive and finite
        sol=solve_system(compiled, cond)

        @test sol[compiled.NOx] > 0
        @test sol[compiled.HOx] > 0
        @test sol[compiled.P_O3_total] > 0
        @test sol[compiled.L_O3_total] > 0
        @test sol[compiled.OPE] > 0
        @test sol[compiled.chain_length] > 0
        @test sol[compiled.L_NOx] > 0
        @test isfinite(sol[compiled.P_O3_net])
    end
end

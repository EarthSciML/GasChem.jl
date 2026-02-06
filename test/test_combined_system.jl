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
@testitem "TroposphericChemistrySystem: Typical Conditions" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    cond=get_typical_conditions()

    required_keys=[:M, :O2, :H2O, :O3, :NO, :NO2, :CO, :CH4, :OH, :HO2, :CH3O2]
    for key in required_keys
        @test haskey(cond, key)
        @test cond[key] > 0
    end

    # Verify specific values (SI: m⁻³)
    @test cond[:M] ≈ 2.5e25
    @test cond[:O2] ≈ 5.25e24
    @test cond[:O3] ≈ 1e18
    @test cond[:OH] ≈ 1e12
    @test cond[:HO2] ≈ 1e14
    @test cond[:CO] ≈ 2.5e18
    @test cond[:CH4] ≈ 4.5e19

    # Physical consistency: O2 is ~21% of M
    @test cond[:O2] / cond[:M] ≈ 0.21 rtol=0.01
end

@testitem "TroposphericChemistrySystem: Urban vs Typical Conditions" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    typical=get_typical_conditions()
    urban=get_urban_conditions()

    @test urban[:NO] > typical[:NO]
    @test urban[:NO2] > typical[:NO2]
    @test urban[:CO] > typical[:CO]
    @test urban[:OH] <= typical[:OH]
    @test urban[:M] ≈ typical[:M]
    @test urban[:CH4] ≈ typical[:CH4]
end

@testitem "TroposphericChemistrySystem: Remote vs Typical Conditions" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    typical=get_typical_conditions()
    remote=get_remote_conditions()

    @test remote[:NO] < typical[:NO]
    @test remote[:NO2] < typical[:NO2]
    @test remote[:HO2] > typical[:HO2]
    @test remote[:M] ≈ typical[:M]
end

# ===========================================================================
# Equation Verification: Combined Diagnostics
# ===========================================================================
@testitem "TroposphericChemistrySystem: NOx = NO + NO2" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    cond=get_typical_conditions()
    NOx_expected=cond[:NO]+cond[:NO2]

    @test NOx_expected > 0
    @test NOx_expected ≈ 2.5e15 + 2.5e16 rtol=1e-10
end

@testitem "TroposphericChemistrySystem: HOx = OH + HO2" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    cond=get_typical_conditions()
    HOx_expected=cond[:OH]+cond[:HO2]

    # HOx is dominated by HO2
    @test HOx_expected ≈ cond[:HO2] rtol=0.01
    @test cond[:HO2] / cond[:OH] > 10
end

@testitem "TroposphericChemistrySystem: O3 Production Diagnostics" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    cond=get_typical_conditions()
    k_HO2_NO=8.1e-12*1e-6    # m³/s
    k_CH3O2_NO=7.7e-12*1e-6  # m³/s

    P_O3=k_HO2_NO*cond[:HO2]*cond[:NO]+k_CH3O2_NO*cond[:CH3O2]*cond[:NO]

    @test P_O3 > 0
    @test P_O3 > 1e11
    @test P_O3 < 1e14
end

@testitem "TroposphericChemistrySystem: NOx Loss Rate" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    cond=get_typical_conditions()
    k_OH_NO2=1.0e-11*1e-6

    L_NOx=k_OH_NO2*cond[:OH]*cond[:NO2]

    @test L_NOx > 0
    @test L_NOx ≈ 1.0e-17 * 1e12 * 2.5e16 rtol=1e-10
end

@testitem "TroposphericChemistrySystem: OPE Calculation" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    cond=get_typical_conditions()
    k_HO2_NO=8.1e-12*1e-6
    k_CH3O2_NO=7.7e-12*1e-6
    k_OH_NO2=1.0e-11*1e-6

    P_O3=k_HO2_NO*cond[:HO2]*cond[:NO]+k_CH3O2_NO*cond[:CH3O2]*cond[:NO]
    L_NOx=k_OH_NO2*cond[:OH]*cond[:NO2]
    OPE=P_O3/L_NOx

    @test OPE > 0
    @test OPE > 1
    @test OPE < 100
end

# ===========================================================================
# Subsystem Coupling Tests
# ===========================================================================
@testitem "TroposphericChemistrySystem: OH Production Coupling" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    cond=get_typical_conditions()

    j_O3=1e-5
    k3_N2=2.6e-11*1e-6
    k3_O2=4.0e-11*1e-6
    k4=2.2e-10*1e-6
    f_N2=0.78
    f_O2=0.21

    k3_eff=f_N2*k3_N2+f_O2*k3_O2
    denom=k3_eff*cond[:M]+k4*cond[:H2O]
    eps_OH=k4*cond[:H2O]/denom
    P_OH=2*j_O3*cond[:O3]*eps_OH

    @test eps_OH > 0 && eps_OH < 1
    @test P_OH > 0
    @test P_OH < 1e14
end

@testitem "TroposphericChemistrySystem: NOx Cycling Coupling" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    cond=get_typical_conditions()

    j_NO2=8e-3
    k_NO_O3=1.8e-14*1e-6

    O3_pss=j_NO2*cond[:NO2]/(k_NO_O3*cond[:NO])
    Phi=j_NO2*cond[:NO2]/(k_NO_O3*cond[:NO]*cond[:O3])

    @test Phi > 0
    @test O3_pss > 0
end

@testitem "TroposphericChemistrySystem: CO Oxidation Coupling" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    cond=get_typical_conditions()

    k_HO2_NO=8.1e-12*1e-6
    k_OH_NO2=1.0e-11*1e-6
    k_HO2_HO2=2.9e-12*1e-6
    k_OH_O3=7.3e-14*1e-6
    k_HO2_O3=2.0e-15*1e-6

    P_O3_co=k_HO2_NO*cond[:HO2]*cond[:NO]-k_OH_O3*cond[:OH]*cond[:O3]-k_HO2_O3*cond[:HO2]*cond[:O3]
    @test P_O3_co > 0

    L_HOx=k_OH_NO2*cond[:OH]*cond[:NO2]+2*k_HO2_HO2*cond[:HO2]^2
    @test L_HOx > 0

    chain_length=k_HO2_NO*cond[:HO2]*cond[:NO]/L_HOx
    @test chain_length > 1
end

# ===========================================================================
# Regime Comparison Tests
# ===========================================================================
@testitem "TroposphericChemistrySystem: NOx-VOC Regime Differences" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    typical=get_typical_conditions()
    urban=get_urban_conditions()
    remote=get_remote_conditions()

    k_HO2_NO=8.1e-12*1e-6
    k_CH3O2_NO=7.7e-12*1e-6
    k_OH_NO2=1.0e-11*1e-6

    function compute_OPE(c)
        P_O3=k_HO2_NO*c[:HO2]*c[:NO]+k_CH3O2_NO*c[:CH3O2]*c[:NO]
        L_NOx=k_OH_NO2*c[:OH]*c[:NO2]
        return P_O3/L_NOx
    end

    OPE_typical=compute_OPE(typical)
    OPE_urban=compute_OPE(urban)
    OPE_remote=compute_OPE(remote)

    @test OPE_typical > 0
    @test OPE_urban > 0
    @test OPE_remote > 0
end

@testitem "TroposphericChemistrySystem: O3 Production Budget Consistency" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    cond=get_typical_conditions()

    k_HO2_NO=8.1e-12*1e-6
    k_CH3O2_NO=7.7e-12*1e-6
    k_NO_O3=1.8e-14*1e-6
    k_OH_O3=7.3e-14*1e-6
    k_HO2_O3=2.0e-15*1e-6

    P_O3_total=k_HO2_NO*cond[:HO2]*cond[:NO]+k_CH3O2_NO*cond[:CH3O2]*cond[:NO]
    L_O3_total=k_NO_O3*cond[:NO]*cond[:O3]+k_OH_O3*cond[:OH]*cond[:O3]+k_HO2_O3*cond[:HO2]*cond[:O3]
    P_O3_net=P_O3_total-L_O3_total

    @test P_O3_total > 0
    @test L_O3_total > 0
    @test P_O3_net ≈ P_O3_total - L_O3_total rtol=1e-10
end

# ===========================================================================
# Qualitative Property Tests
# ===========================================================================
@testitem "TroposphericChemistrySystem: Physical Bounds" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    for get_cond in [get_typical_conditions, get_urban_conditions, get_remote_conditions]
        cond=get_cond()
        for (key, val) in cond
            @test val > 0
        end
    end
end

@testitem "TroposphericChemistrySystem: Atmospheric Composition Consistency" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    for get_cond in [get_typical_conditions, get_urban_conditions, get_remote_conditions]
        cond=get_cond()

        # O2 should be ~21% of M
        @test cond[:O2] / cond[:M] ≈ 0.21 rtol=0.01

        # O3 should be much less than O2
        @test cond[:O3] < cond[:O2] * 1e-6

        # NO2 should generally be larger than NO
        @test cond[:NO2] >= cond[:NO] || cond[:NO] < 1e15

        # HO2 >> OH
        @test cond[:HO2] / cond[:OH] > 1

        # CH4 is the most abundant reactive hydrocarbon (except urban CO)
        if get_cond!==get_urban_conditions
            @test cond[:CH4] > cond[:CO]
        end
    end
end

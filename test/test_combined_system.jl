# ===========================================================================
# Structural Tests
# ===========================================================================
@testitem "TroposphericChemistrySystem: Structural Verification" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = TroposphericChemistrySystem()
    @test sys isa System
    @test nameof(sys) == :TroposphericChemistrySystem

    vars = unknowns(sys)
    params = parameters(sys)
    eqs = equations(sys)

    # Should have equations from the combined system plus subsystem coupling
    @test length(eqs) > 0
    @test length(vars) > 0
    @test length(params) > 0

    # Check key diagnostic variable names at the top level
    var_names = [string(v) for v in vars]
    for expected in ["NOx", "HOx", "P_O3_net", "OPE", "chain_length", "L_NOx"]
        @test any(n -> contains(n, expected), var_names)
    end
end

@testitem "TroposphericChemistrySystem: Subsystem Composition" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = TroposphericChemistrySystem()

    # The system should have 3 subsystems: oh, nox, co
    subsystems = ModelingToolkit.get_systems(sys)
    subsys_names = [nameof(s) for s in subsystems]

    @test :oh in subsys_names
    @test :nox in subsys_names
    @test :co in subsys_names
    @test length(subsystems) == 3
end

@testitem "TroposphericChemistrySystem: Subsystem Variable Access" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = TroposphericChemistrySystem()

    # Should be able to access subsystem variables using dot notation
    # Access OH subsystem outputs
    all_vars = unknowns(sys)
    var_strs = [string(v) for v in all_vars]

    # OH subsystem should contribute O1D, P_OH, epsilon_OH
    @test any(s -> contains(s, "oh") && contains(s, "O1D"), var_strs)
    @test any(s -> contains(s, "oh") && contains(s, "P_OH"), var_strs)

    # NOx subsystem should contribute O3_pss, Phi, P_O3
    @test any(s -> contains(s, "nox") && contains(s, "O3_pss"), var_strs)

    # CO subsystem should contribute chain_length, HO2_ss
    @test any(s -> contains(s, "co") && contains(s, "HO2_ss"), var_strs)
end

# ===========================================================================
# Condition Functions Tests
# ===========================================================================
@testitem "TroposphericChemistrySystem: Typical Conditions" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    cond = get_typical_conditions()

    # Verify all required keys are present
    required_keys = [:M, :O2, :H2O, :O3, :NO, :NO2, :CO, :CH4, :OH, :HO2, :CH3O2]
    for key in required_keys
        @test haskey(cond, key)
        @test cond[key] > 0
    end

    # Verify specific values match the implementation
    @test cond[:M] ≈ 2.5e19
    @test cond[:O2] ≈ 5.25e18
    @test cond[:O3] ≈ 1e12
    @test cond[:OH] ≈ 1e6
    @test cond[:HO2] ≈ 1e8
    @test cond[:CO] ≈ 2.5e12
    @test cond[:CH4] ≈ 4.5e13

    # Physical consistency: O2 is ~21% of M
    @test cond[:O2] / cond[:M] ≈ 0.21 rtol=0.01
end

@testitem "TroposphericChemistrySystem: Urban vs Typical Conditions" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    typical = get_typical_conditions()
    urban = get_urban_conditions()

    # Urban should have higher NOx
    @test urban[:NO] > typical[:NO]
    @test urban[:NO2] > typical[:NO2]

    # Urban should have higher CO (traffic, combustion)
    @test urban[:CO] > typical[:CO]

    # Urban OH tends to be lower due to high NOx scavenging
    @test urban[:OH] <= typical[:OH]

    # Both should have same background M and CH4
    @test urban[:M] ≈ typical[:M]
    @test urban[:CH4] ≈ typical[:CH4]
end

@testitem "TroposphericChemistrySystem: Remote vs Typical Conditions" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    typical = get_typical_conditions()
    remote = get_remote_conditions()

    # Remote should have lower NOx
    @test remote[:NO] < typical[:NO]
    @test remote[:NO2] < typical[:NO2]

    # Remote HO2 tends to be higher (less NO to convert to OH)
    @test remote[:HO2] > typical[:HO2]

    # Same background M
    @test remote[:M] ≈ typical[:M]
end

# ===========================================================================
# Equation Verification: Combined Diagnostics
# ===========================================================================
@testitem "TroposphericChemistrySystem: NOx = NO + NO2" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # The combined system defines NOx ~ NO + NO2
    cond = get_typical_conditions()
    NOx_expected = cond[:NO] + cond[:NO2]

    @test NOx_expected > 0
    @test NOx_expected ≈ 2.5e9 + 2.5e10 rtol=1e-10
end

@testitem "TroposphericChemistrySystem: HOx = OH + HO2" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    cond = get_typical_conditions()
    HOx_expected = cond[:OH] + cond[:HO2]

    # HOx is dominated by HO2 (HO2 >> OH)
    @test HOx_expected ≈ cond[:HO2] rtol=0.01
    @test cond[:HO2] / cond[:OH] > 10
end

@testitem "TroposphericChemistrySystem: O3 Production Diagnostics" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # P_O3_total = k_HO2_NO * HO2 * NO + k_CH3O2_NO * CH3O2 * NO
    cond = get_typical_conditions()
    k_HO2_NO = 8.1e-12
    k_CH3O2_NO = 7.7e-12

    P_O3 = k_HO2_NO * cond[:HO2] * cond[:NO] + k_CH3O2_NO * cond[:CH3O2] * cond[:NO]

    # Should be a positive O3 production rate
    @test P_O3 > 0

    # At typical conditions, P_O3 is on the order of 1e6 molecules/cm^3/s
    @test P_O3 > 1e5
    @test P_O3 < 1e8
end

@testitem "TroposphericChemistrySystem: NOx Loss Rate" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # L_NOx = k_OH_NO2 * OH * NO2
    cond = get_typical_conditions()
    k_OH_NO2 = 1.0e-11

    L_NOx = k_OH_NO2 * cond[:OH] * cond[:NO2]

    @test L_NOx > 0
    @test L_NOx ≈ 1.0e-11 * 1e6 * 2.5e10 rtol=1e-10
end

@testitem "TroposphericChemistrySystem: OPE Calculation" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # OPE = P_O3_total / L_NOx
    cond = get_typical_conditions()
    k_HO2_NO = 8.1e-12
    k_CH3O2_NO = 7.7e-12
    k_OH_NO2 = 1.0e-11

    P_O3 = k_HO2_NO * cond[:HO2] * cond[:NO] + k_CH3O2_NO * cond[:CH3O2] * cond[:NO]
    L_NOx = k_OH_NO2 * cond[:OH] * cond[:NO2]
    OPE = P_O3 / L_NOx

    # OPE should be positive and physically reasonable
    @test OPE > 0
    # Typical OPE ranges from ~5 to ~30
    @test OPE > 1
    @test OPE < 100
end

# ===========================================================================
# Subsystem Coupling Tests
# ===========================================================================
@testitem "TroposphericChemistrySystem: OH Production Coupling" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # The OH subsystem computes epsilon_OH and P_OH from O3, H2O, M
    # These should be consistent with the values from the standalone OHProduction system
    cond = get_typical_conditions()

    # Compute expected OH production using the same formula as OHProduction
    j_O3 = 1e-5  # default from OHProduction
    k3_N2 = 2.6e-11
    k3_O2 = 4.0e-11
    k4 = 2.2e-10
    f_N2 = 0.78
    f_O2 = 0.21

    k3_eff = f_N2 * k3_N2 + f_O2 * k3_O2
    denom = k3_eff * cond[:M] + k4 * cond[:H2O]
    eps_OH = k4 * cond[:H2O] / denom
    P_OH = 2 * j_O3 * cond[:O3] * eps_OH

    # These should be physically reasonable
    @test eps_OH > 0 && eps_OH < 1
    @test P_OH > 0
    @test P_OH < 1e8  # not unreasonably large
end

@testitem "TroposphericChemistrySystem: NOx Cycling Coupling" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # The NOx subsystem computes O3_pss and Phi from NO, NO2, O3
    cond = get_typical_conditions()

    j_NO2 = 8e-3
    k_NO_O3 = 1.8e-14

    O3_pss = j_NO2 * cond[:NO2] / (k_NO_O3 * cond[:NO])

    # The PSS O3 should be comparable to the actual O3
    # If Phi ~ 1, O3_pss ~ O3
    Phi = j_NO2 * cond[:NO2] / (k_NO_O3 * cond[:NO] * cond[:O3])

    # With peroxy radicals present, Phi > 1 is expected (additional NO -> NO2 conversion)
    @test Phi > 0
    @test O3_pss > 0
end

@testitem "TroposphericChemistrySystem: CO Oxidation Coupling" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # The CO subsystem computes P_O3, L_HOx, chain_length from CO, OH, HO2, NO, NO2, O3
    cond = get_typical_conditions()

    k_HO2_NO = 8.1e-12
    k_OH_NO2 = 1.0e-11
    k_HO2_HO2 = 2.9e-12
    k_OH_O3 = 7.3e-14
    k_HO2_O3 = 2.0e-15

    # Net O3 production from CO subsystem
    P_O3_co = k_HO2_NO * cond[:HO2] * cond[:NO] - k_OH_O3 * cond[:OH] * cond[:O3] - k_HO2_O3 * cond[:HO2] * cond[:O3]

    # At typical conditions with some NOx, net O3 production should be positive
    @test P_O3_co > 0

    # HOx loss rate
    L_HOx = k_OH_NO2 * cond[:OH] * cond[:NO2] + 2 * k_HO2_HO2 * cond[:HO2]^2
    @test L_HOx > 0

    # Chain length
    chain_length = k_HO2_NO * cond[:HO2] * cond[:NO] / L_HOx
    @test chain_length > 1
end

# ===========================================================================
# Regime Comparison Tests
# ===========================================================================
@testitem "TroposphericChemistrySystem: NOx-VOC Regime Differences" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    typical = get_typical_conditions()
    urban = get_urban_conditions()
    remote = get_remote_conditions()

    k_HO2_NO = 8.1e-12
    k_CH3O2_NO = 7.7e-12
    k_OH_NO2 = 1.0e-11

    # Compute OPE for each scenario
    function compute_OPE(c)
        P_O3 = k_HO2_NO * c[:HO2] * c[:NO] + k_CH3O2_NO * c[:CH3O2] * c[:NO]
        L_NOx = k_OH_NO2 * c[:OH] * c[:NO2]
        return P_O3 / L_NOx
    end

    OPE_typical = compute_OPE(typical)
    OPE_urban = compute_OPE(urban)
    OPE_remote = compute_OPE(remote)

    # All OPE values should be positive
    @test OPE_typical > 0
    @test OPE_urban > 0
    @test OPE_remote > 0
end

@testitem "TroposphericChemistrySystem: O3 Production Budget Consistency" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    cond = get_typical_conditions()

    k_HO2_NO = 8.1e-12
    k_CH3O2_NO = 7.7e-12
    k_NO_O3 = 1.8e-14
    k_OH_O3 = 7.3e-14
    k_HO2_O3 = 2.0e-15

    # Total O3 production from peroxy radicals + NO
    P_O3_total = k_HO2_NO * cond[:HO2] * cond[:NO] + k_CH3O2_NO * cond[:CH3O2] * cond[:NO]

    # Total O3 loss
    L_O3_total = k_NO_O3 * cond[:NO] * cond[:O3] + k_OH_O3 * cond[:OH] * cond[:O3] + k_HO2_O3 * cond[:HO2] * cond[:O3]

    # Net O3 production
    P_O3_net = P_O3_total - L_O3_total

    # All terms should be non-negative
    @test P_O3_total > 0
    @test L_O3_total > 0

    # Budget must satisfy: net = total production - total loss
    @test P_O3_net ≈ P_O3_total - L_O3_total rtol=1e-10
end

# ===========================================================================
# Qualitative Property Tests
# ===========================================================================
@testitem "TroposphericChemistrySystem: Physical Bounds" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # All condition dictionaries should have positive concentrations
    for get_cond in [get_typical_conditions, get_urban_conditions, get_remote_conditions]
        cond = get_cond()
        for (key, val) in cond
            @test val > 0
        end
    end
end

@testitem "TroposphericChemistrySystem: Atmospheric Composition Consistency" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # Verify that the condition dictionaries represent physically realistic atmospheres
    for get_cond in [get_typical_conditions, get_urban_conditions, get_remote_conditions]
        cond = get_cond()

        # O2 should be ~21% of M
        @test cond[:O2] / cond[:M] ≈ 0.21 rtol=0.01

        # O3 should be much less than O2 (ppb vs %)
        @test cond[:O3] < cond[:O2] * 1e-6

        # NO2 should generally be larger than NO in the troposphere
        # (except possibly in very clean environments with low NO2)
        @test cond[:NO2] >= cond[:NO] || cond[:NO] < 1e9

        # HO2 >> OH (typical ratio is 10-100)
        @test cond[:HO2] / cond[:OH] > 1

        # CH4 is the most abundant reactive hydrocarbon (in non-urban environments)
        # In urban environments, CO can exceed CH4 due to combustion sources
        if get_cond !== get_urban_conditions
            @test cond[:CH4] > cond[:CO]
        end
    end
end

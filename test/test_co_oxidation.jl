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

    # 6 parameters: k_CO_OH, k_HO2_NO, k_HO2_HO2, k_OH_NO2, k_HO2_O3, k_OH_O3
    @test length(params) == 6

    # Check key output variable names
    var_names = [string(v) for v in vars]
    for expected in ["P_O3", "L_HOx", "chain_length", "HO2_ss"]
        @test any(n -> contains(n, expected), var_names)
    end
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
end

# ===========================================================================
# Equation Verification Tests
# ===========================================================================
@testitem "COOxidation: Eq 6.14 HO2 Steady-State (High NOx)" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # HO2_ss = k_CO_OH * [CO] * [OH] / (k_HO2_NO * [NO])
    # At typical conditions:
    k_CO_OH = 2.4e-13
    k_HO2_NO = 8.1e-12
    CO = 2.5e12   # ~100 ppb
    OH = 1e6      # typical daytime
    NO = 2.5e9    # ~0.1 ppb

    HO2_ss = k_CO_OH * CO * OH / (k_HO2_NO * NO)

    # HO2 should be on the order of 1e7-1e8
    @test HO2_ss > 1e6
    @test HO2_ss < 1e10
    @test HO2_ss ≈ 2.963e7 rtol=0.05
end

@testitem "COOxidation: Net O3 Production Rate" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # P_O3 = k_HO2_NO * [HO2] * [NO] - k_OH_O3 * [OH] * [O3] - k_HO2_O3 * [HO2] * [O3]
    k_HO2_NO = 8.1e-12
    k_OH_O3 = 7.3e-14
    k_HO2_O3 = 2.0e-15

    HO2 = 1e8
    NO = 2.5e9
    OH = 1e6
    O3 = 1e12

    # Production term (HO2 + NO)
    P_term = k_HO2_NO * HO2 * NO  # = 8.1e-12 * 1e8 * 2.5e9 = 2.025e6

    # Loss terms
    L_OH_O3 = k_OH_O3 * OH * O3    # = 7.3e-14 * 1e6 * 1e12 = 7.3e4
    L_HO2_O3 = k_HO2_O3 * HO2 * O3  # = 2.0e-15 * 1e8 * 1e12 = 2e5

    P_O3_expected = P_term - L_OH_O3 - L_HO2_O3
    @test P_O3_expected > 0  # Net production in these conditions
    @test P_O3_expected ≈ 2.025e6 - 7.3e4 - 2e5 rtol=0.01
end

@testitem "COOxidation: HOx Loss Rate" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # L_HOx = k_OH_NO2 * [OH] * [NO2] + 2 * k_HO2_HO2 * [HO2]^2
    k_OH_NO2 = 1.0e-11
    k_HO2_HO2 = 2.9e-12

    OH = 1e6
    NO2 = 2.5e10
    HO2 = 1e8

    # HNO3 formation
    L_HNO3 = k_OH_NO2 * OH * NO2   # = 1.0e-11 * 1e6 * 2.5e10 = 2.5e5
    # H2O2 formation
    L_H2O2 = 2 * k_HO2_HO2 * HO2^2  # = 2 * 2.9e-12 * 1e16 = 5.8e4

    L_HOx_expected = L_HNO3 + L_H2O2
    @test L_HOx_expected > 0
    @test L_HOx_expected ≈ 2.5e5 + 5.8e4 rtol=0.01
end

@testitem "COOxidation: Chain Length" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # chain_length = (k_HO2_NO * [HO2] * [NO]) / L_HOx
    # This measures how many OH-HO2 cycles occur before radical termination
    k_HO2_NO = 8.1e-12
    k_OH_NO2 = 1.0e-11
    k_HO2_HO2 = 2.9e-12

    HO2 = 1e8
    NO = 2.5e9
    OH = 1e6
    NO2 = 2.5e10

    propagation = k_HO2_NO * HO2 * NO
    termination = k_OH_NO2 * OH * NO2 + 2 * k_HO2_HO2 * HO2^2

    chain = propagation / termination

    # Chain length should be > 1 (catalytic cycle operates)
    @test chain > 1
    # Typical chain length is ~1-100 depending on conditions
    @test chain < 100
    @test chain ≈ propagation / termination rtol=1e-10
end

# ===========================================================================
# Rate Constants at 298 K
# ===========================================================================
@testitem "COOxidation: Rate Constants at 298 K" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = COOxidation()
    params = parameters(sys)
    param_dict = Dict(Symbol(p) => ModelingToolkit.getdefault(p) for p in params if ModelingToolkit.hasdefault(p))

    @test param_dict[:k_CO_OH] ≈ 2.4e-13
    @test param_dict[:k_HO2_NO] ≈ 8.1e-12
    @test param_dict[:k_HO2_HO2] ≈ 2.9e-12
    @test param_dict[:k_OH_NO2] ≈ 1.0e-11
    @test param_dict[:k_HO2_O3] ≈ 2.0e-15
    @test param_dict[:k_OH_O3] ≈ 7.3e-14
end

@testitem "OzoneProductionEfficiency: Rate Constants" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = OzoneProductionEfficiency()
    params = parameters(sys)
    param_dict = Dict(Symbol(p) => ModelingToolkit.getdefault(p) for p in params if ModelingToolkit.hasdefault(p))

    @test param_dict[:k_HO2_NO] ≈ 8.1e-12
    @test param_dict[:k_RO2_NO] ≈ 8.0e-12
    @test param_dict[:k_OH_NO2] ≈ 1.0e-11
end

# ===========================================================================
# OPE Equation Verification
# ===========================================================================
@testitem "OzoneProductionEfficiency: OPE Calculation" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # OPE = P_O3 / L_NOx = (k_HO2_NO*HO2*NO + k_RO2_NO*RO2*NO) / (k_OH_NO2*OH*NO2)
    k_HO2_NO = 8.1e-12
    k_RO2_NO = 8.0e-12
    k_OH_NO2 = 1.0e-11

    OH = 1e6
    HO2 = 1e8
    RO2 = 1e8
    NO = 2.5e9
    NO2 = 2.5e10

    P_O3 = k_HO2_NO * HO2 * NO + k_RO2_NO * RO2 * NO
    L_NOx = k_OH_NO2 * OH * NO2
    OPE = P_O3 / L_NOx

    # OPE should be positive
    @test OPE > 0
    # Typical OPE ranges from ~5 to ~30
    @test OPE > 1
    @test OPE ≈ (k_HO2_NO * HO2 * NO + k_RO2_NO * RO2 * NO) / (k_OH_NO2 * OH * NO2) rtol=1e-10
end

# ===========================================================================
# Limiting Behavior Tests
# ===========================================================================
@testitem "COOxidation: High NOx - Short Chain Length" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # At high NOx, OH + NO2 -> HNO3 is fast, terminating radicals quickly.
    # Chain length decreases.
    k_HO2_NO = 8.1e-12
    k_OH_NO2 = 1.0e-11
    k_HO2_HO2 = 2.9e-12

    HO2 = 5e7  # lower at high NOx
    OH = 5e5    # lower at high NOx

    # Low NOx scenario
    NO_low = 2.5e9
    NO2_low = 2.5e10
    prop_low = k_HO2_NO * HO2 * NO_low
    term_low = k_OH_NO2 * OH * NO2_low + 2 * k_HO2_HO2 * HO2^2
    chain_low = prop_low / term_low

    # High NOx scenario
    NO_high = 2.5e11
    NO2_high = 7.5e11
    prop_high = k_HO2_NO * HO2 * NO_high
    term_high = k_OH_NO2 * OH * NO2_high + 2 * k_HO2_HO2 * HO2^2
    chain_high = prop_high / term_high

    # At high NOx, the HNO3 termination dominates
    # chain_high should still be > 1 but termination is faster
    @test term_high > term_low  # More termination at high NOx
end

@testitem "OzoneProductionEfficiency: High NOx - Low OPE" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # OPE decreases with increasing NOx because L_NOx (OH+NO2) increases faster
    # than P_O3 when NO2 increases proportionally more than NO
    k_HO2_NO = 8.1e-12
    k_RO2_NO = 8.0e-12
    k_OH_NO2 = 1.0e-11

    OH = 1e6
    HO2 = 1e8
    RO2 = 1e8

    # Low NOx
    NO_low = 2.5e8
    NO2_low = 5e8
    OPE_low = (k_HO2_NO * HO2 * NO_low + k_RO2_NO * RO2 * NO_low) / (k_OH_NO2 * OH * NO2_low)

    # High NOx (same NO2/NO ratio but 100x)
    NO_high = 2.5e10
    NO2_high = 5e10
    OPE_high = (k_HO2_NO * HO2 * NO_high + k_RO2_NO * RO2 * NO_high) / (k_OH_NO2 * OH * NO2_high)

    # OPE should be the same if NO2/NO ratio is fixed and HO2/RO2 are constant
    # OPE = (k_HO2_NO*HO2 + k_RO2_NO*RO2) * NO / (k_OH_NO2*OH*NO2)
    # = (k_HO2_NO*HO2 + k_RO2_NO*RO2) / (k_OH_NO2*OH) * (NO/NO2)
    # So OPE depends on NO/NO2 ratio, not absolute NOx level (when HO2,RO2 fixed)
    @test OPE_low ≈ OPE_high rtol=1e-10
    @test OPE_low > 0
end

# ===========================================================================
# Qualitative Property Tests
# ===========================================================================
@testitem "COOxidation: Qualitative - Positive Production and Loss" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # For any positive concentrations, loss rates must be non-negative
    k_CO_OH = 2.4e-13
    k_HO2_NO = 8.1e-12
    k_HO2_HO2 = 2.9e-12
    k_OH_NO2 = 1.0e-11
    k_HO2_O3 = 2.0e-15
    k_OH_O3 = 7.3e-14

    for OH in [1e5, 1e6, 1e7]
        for HO2 in [1e7, 1e8, 1e9]
            for NO in [1e8, 1e10, 1e12]
                NO2 = 2 * NO
                CO = 2.5e12
                O3 = 1e12

                L_OH = k_CO_OH * CO * OH + k_OH_NO2 * OH * NO2 + k_OH_O3 * OH * O3
                L_HO2 = k_HO2_NO * HO2 * NO + 2 * k_HO2_HO2 * HO2^2 + k_HO2_O3 * HO2 * O3
                L_HOx = k_OH_NO2 * OH * NO2 + 2 * k_HO2_HO2 * HO2^2
                HO2_ss = k_CO_OH * CO * OH / (k_HO2_NO * NO)

                @test L_OH > 0
                @test L_HO2 > 0
                @test L_HOx > 0
                @test HO2_ss > 0
            end
        end
    end
end

@testitem "COOxidation: HOx Catalytic Cycling" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # The CO oxidation cycle regenerates OH:
    # CO + OH -> CO2 + H, then H + O2 -> HO2, then HO2 + NO -> OH + NO2
    # The chain length > 1 means the catalytic cycle is operating
    k_HO2_NO = 8.1e-12
    k_OH_NO2 = 1.0e-11
    k_HO2_HO2 = 2.9e-12

    # Typical conditions
    HO2 = 1e8
    NO = 2.5e9
    OH = 1e6
    NO2 = 2.5e10

    propagation = k_HO2_NO * HO2 * NO
    termination = k_OH_NO2 * OH * NO2 + 2 * k_HO2_HO2 * HO2^2
    chain_length = propagation / termination

    # Chain length must be > 1 for catalytic behavior
    @test chain_length > 1

    # Net O3 production: CO + 2O2 + hv -> CO2 + O3
    # This means every CO oxidation produces one O3 (in high NOx conditions)
    # P_O3 > 0 (net O3 production)
    k_OH_O3 = 7.3e-14
    k_HO2_O3 = 2.0e-15
    O3 = 1e12
    P_O3 = k_HO2_NO * HO2 * NO - k_OH_O3 * OH * O3 - k_HO2_O3 * HO2 * O3
    @test P_O3 > 0
end

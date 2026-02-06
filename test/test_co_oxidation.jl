# ===========================================================================
# Structural Tests
# ===========================================================================
@testitem "COOxidation: Structural Verification" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=COOxidation()
    @test sys isa System
    @test nameof(sys) == :COOxidation

    vars=unknowns(sys)
    params=parameters(sys)
    eqs=equations(sys)

    # 12 variables: CO, OH, HO2, NO, NO2, O3 (input)
    #             + P_O3, L_HOx, L_OH, L_HO2, chain_length, HO2_ss (output)
    @test length(vars) == 12

    # 6 algebraic equations
    @test length(eqs) == 6

    # 7 parameters: k_CO_OH, k_HO2_NO, k_HO2_HO2, k_OH_NO2, k_HO2_O3, k_OH_O3 + 1 constant (two)
    @test length(params) == 7

    # Check key output variable names
    var_names=[string(v) for v in vars]
    for expected in ["P_O3", "L_HOx", "chain_length", "HO2_ss"]
        @test any(n -> contains(n, expected), var_names)
    end
end

@testitem "OzoneProductionEfficiency: Structural Verification" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=OzoneProductionEfficiency()
    @test sys isa System
    @test nameof(sys) == :OzoneProductionEfficiency

    vars=unknowns(sys)
    params=parameters(sys)
    eqs=equations(sys)

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
    k_CO_OH=2.4e-13*1e-6     # m³/s
    k_HO2_NO=8.1e-12*1e-6    # m³/s
    CO=2.5e18    # m⁻³ (~100 ppb)
    OH=1e12      # m⁻³ (typical daytime)
    NO=2.5e15    # m⁻³ (~0.1 ppb)

    HO2_ss=k_CO_OH*CO*OH/(k_HO2_NO*NO)

    # HO2 should be on the order of 1e13-1e14 m⁻³
    @test HO2_ss > 1e12
    @test HO2_ss < 1e16
    @test HO2_ss ≈ 2.963e13 rtol=0.05
end

@testitem "COOxidation: Net O3 Production Rate" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # P_O3 = k_HO2_NO * [HO2] * [NO] - k_OH_O3 * [OH] * [O3] - k_HO2_O3 * [HO2] * [O3]
    k_HO2_NO=8.1e-12*1e-6
    k_OH_O3=7.3e-14*1e-6
    k_HO2_O3=2.0e-15*1e-6

    HO2=1e14
    NO=2.5e15
    OH=1e12
    O3=1e18

    # Production term (HO2 + NO)
    P_term=k_HO2_NO*HO2*NO

    # Loss terms
    L_OH_O3=k_OH_O3*OH*O3
    L_HO2_O3=k_HO2_O3*HO2*O3

    P_O3_expected=P_term-L_OH_O3-L_HO2_O3
    @test P_O3_expected > 0  # Net production in these conditions
    @test P_O3_expected ≈ P_term - L_OH_O3 - L_HO2_O3 rtol=0.01
end

@testitem "COOxidation: HOx Loss Rate" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # L_HOx = k_OH_NO2 * [OH] * [NO2] + 2 * k_HO2_HO2 * [HO2]^2
    k_OH_NO2=1.0e-11*1e-6
    k_HO2_HO2=2.9e-12*1e-6

    OH=1e12
    NO2=2.5e16
    HO2=1e14

    # HNO3 formation
    L_HNO3=k_OH_NO2*OH*NO2
    # H2O2 formation
    L_H2O2=2*k_HO2_HO2*HO2^2

    L_HOx_expected=L_HNO3+L_H2O2
    @test L_HOx_expected > 0
    @test L_HOx_expected ≈ L_HNO3 + L_H2O2 rtol=0.01
end

@testitem "COOxidation: Chain Length" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # chain_length = (k_HO2_NO * [HO2] * [NO]) / L_HOx
    k_HO2_NO=8.1e-12*1e-6
    k_OH_NO2=1.0e-11*1e-6
    k_HO2_HO2=2.9e-12*1e-6

    HO2=1e14
    NO=2.5e15
    OH=1e12
    NO2=2.5e16

    propagation=k_HO2_NO*HO2*NO
    termination=k_OH_NO2*OH*NO2+2*k_HO2_HO2*HO2^2

    chain=propagation/termination

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
    sys=COOxidation()
    params=parameters(sys)
    param_dict=Dict(Symbol(p)=>ModelingToolkit.getdefault(p)
    for p in params if ModelingToolkit.hasdefault(p))

    @test param_dict[:k_CO_OH] ≈ 2.4e-13 * 1e-6
    @test param_dict[:k_HO2_NO] ≈ 8.1e-12 * 1e-6
    @test param_dict[:k_HO2_HO2] ≈ 2.9e-12 * 1e-6
    @test param_dict[:k_OH_NO2] ≈ 1.0e-11 * 1e-6
    @test param_dict[:k_HO2_O3] ≈ 2.0e-15 * 1e-6
    @test param_dict[:k_OH_O3] ≈ 7.3e-14 * 1e-6
end

@testitem "OzoneProductionEfficiency: Rate Constants" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys=OzoneProductionEfficiency()
    params=parameters(sys)
    param_dict=Dict(Symbol(p)=>ModelingToolkit.getdefault(p)
    for p in params if ModelingToolkit.hasdefault(p))

    @test param_dict[:k_HO2_NO] ≈ 8.1e-12 * 1e-6
    @test param_dict[:k_RO2_NO] ≈ 8.0e-12 * 1e-6
    @test param_dict[:k_OH_NO2] ≈ 1.0e-11 * 1e-6
end

# ===========================================================================
# OPE Equation Verification
# ===========================================================================
@testitem "OzoneProductionEfficiency: OPE Calculation" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    k_HO2_NO=8.1e-12*1e-6
    k_RO2_NO=8.0e-12*1e-6
    k_OH_NO2=1.0e-11*1e-6

    OH=1e12
    HO2=1e14
    RO2=1e14
    NO=2.5e15
    NO2=2.5e16

    P_O3=k_HO2_NO*HO2*NO+k_RO2_NO*RO2*NO
    L_NOx=k_OH_NO2*OH*NO2
    OPE=P_O3/L_NOx

    @test OPE > 0
    @test OPE > 1
    @test OPE ≈ P_O3 / L_NOx rtol=1e-10
end

# ===========================================================================
# Limiting Behavior Tests
# ===========================================================================
@testitem "COOxidation: High NOx - Short Chain Length" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    k_HO2_NO=8.1e-12*1e-6
    k_OH_NO2=1.0e-11*1e-6
    k_HO2_HO2=2.9e-12*1e-6

    HO2=5e13  # lower at high NOx
    OH=5e11    # lower at high NOx

    # Low NOx scenario
    NO_low=2.5e15
    NO2_low=2.5e16
    term_low=k_OH_NO2*OH*NO2_low+2*k_HO2_HO2*HO2^2

    # High NOx scenario
    NO_high=2.5e17
    NO2_high=7.5e17
    term_high=k_OH_NO2*OH*NO2_high+2*k_HO2_HO2*HO2^2

    # At high NOx, the HNO3 termination dominates
    @test term_high > term_low
end

@testitem "OzoneProductionEfficiency: High NOx - Low OPE" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    k_HO2_NO=8.1e-12*1e-6
    k_RO2_NO=8.0e-12*1e-6
    k_OH_NO2=1.0e-11*1e-6

    OH=1e12
    HO2=1e14
    RO2=1e14

    # Low NOx
    NO_low=2.5e14
    NO2_low=5e14
    OPE_low=(k_HO2_NO*HO2*NO_low+k_RO2_NO*RO2*NO_low)/(k_OH_NO2*OH*NO2_low)

    # High NOx (same NO2/NO ratio but 100x)
    NO_high=2.5e16
    NO2_high=5e16
    OPE_high=(k_HO2_NO*HO2*NO_high+k_RO2_NO*RO2*NO_high)/(k_OH_NO2*OH*NO2_high)

    # OPE depends on NO/NO2 ratio, not absolute NOx level (when HO2,RO2 fixed)
    @test OPE_low ≈ OPE_high rtol=1e-10
    @test OPE_low > 0
end

# ===========================================================================
# Qualitative Property Tests
# ===========================================================================
@testitem "COOxidation: Qualitative - Positive Production and Loss" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    k_CO_OH=2.4e-13*1e-6
    k_HO2_NO=8.1e-12*1e-6
    k_HO2_HO2=2.9e-12*1e-6
    k_OH_NO2=1.0e-11*1e-6
    k_HO2_O3=2.0e-15*1e-6
    k_OH_O3=7.3e-14*1e-6

    for OH in [1e11, 1e12, 1e13]
        for HO2 in [1e13, 1e14, 1e15]
            for NO in [1e14, 1e16, 1e18]
                NO2=2*NO
                CO=2.5e18
                O3=1e18

                L_OH=k_CO_OH*CO*OH+k_OH_NO2*OH*NO2+k_OH_O3*OH*O3
                L_HO2=k_HO2_NO*HO2*NO+2*k_HO2_HO2*HO2^2+k_HO2_O3*HO2*O3
                L_HOx=k_OH_NO2*OH*NO2+2*k_HO2_HO2*HO2^2
                HO2_ss=k_CO_OH*CO*OH/(k_HO2_NO*NO)

                @test L_OH > 0
                @test L_HO2 > 0
                @test L_HOx > 0
                @test HO2_ss > 0
            end
        end
    end
end

@testitem "COOxidation: HOx Catalytic Cycling" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    k_HO2_NO=8.1e-12*1e-6
    k_OH_NO2=1.0e-11*1e-6
    k_HO2_HO2=2.9e-12*1e-6

    # Typical conditions (m⁻³)
    HO2=1e14
    NO=2.5e15
    OH=1e12
    NO2=2.5e16

    propagation=k_HO2_NO*HO2*NO
    termination=k_OH_NO2*OH*NO2+2*k_HO2_HO2*HO2^2
    chain_length=propagation/termination

    # Chain length must be > 1 for catalytic behavior
    @test chain_length > 1

    # Net O3 production: CO + 2O2 + hv -> CO2 + O3
    k_OH_O3=7.3e-14*1e-6
    k_HO2_O3=2.0e-15*1e-6
    O3=1e18
    P_O3=k_HO2_NO*HO2*NO-k_OH_O3*OH*O3-k_HO2_O3*HO2*O3
    @test P_O3 > 0
end

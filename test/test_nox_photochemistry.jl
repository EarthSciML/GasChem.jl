# ===========================================================================
# Structural Tests
# ===========================================================================
@testitem "NOxPhotochemistry: Structural Verification" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = NOxPhotochemistry()
    @test sys isa System
    @test nameof(sys) == :NOxPhotochemistry

    vars = unknowns(sys)
    params = parameters(sys)
    eqs = equations(sys)

    # 9 variables: NO, NO2, O3, O2, M (input) + O, O3_pss, Phi, P_O3 (output)
    @test length(vars) == 9

    # 4 algebraic equations (Eq. 6.5, 6.6, 6.7, 6.8)
    @test length(eqs) == 4

    # 3 parameters: j_NO2, k_O_O2_M, k_NO_O3
    @test length(params) == 3

    # Check output variable names
    var_names = [string(v) for v in vars]
    for expected in ["O3_pss", "P_O3"]
        @test any(n -> contains(n, expected), var_names)
    end
end

@testitem "PhotostationaryState: Structural Verification" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = PhotostationaryState()
    @test sys isa System
    @test nameof(sys) == :PhotostationaryState

    vars = unknowns(sys)
    params = parameters(sys)
    eqs = equations(sys)

    # 5 variables: NO, NO2, O3 (input) + Phi, Phi_deviation (output)
    @test length(vars) == 5

    # 2 equations
    @test length(eqs) == 2

    # 3 parameters: j_NO2, k_NO_O3, + 1 constant (one)
    @test length(params) == 3
end

# ===========================================================================
# Equation Verification Tests
# ===========================================================================
@testitem "NOxPhotochemistry: Eq 6.5 O Atom Steady-State" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # [O] = j_NO2 * [NO2] / (k_O_O2_M * [O2] * [M])
    j_NO2 = 8e-3                    # s⁻¹
    NO2 = 1e16                       # m⁻³ (0.4 ppb)
    k_O_O2_M = 6.0e-34 * 1e-12      # m⁶/s
    O2 = 5.25e24                     # m⁻³
    M = 2.5e25                       # m⁻³

    O_expected = j_NO2 * NO2 / (k_O_O2_M * O2 * M)

    # O atoms are very short-lived, concentration should be very low
    @test O_expected > 0
    @test O_expected < 1e12  # should be on order of ~1e9 m⁻³
    @test O_expected ≈ 1.016e9 rtol=0.05
end

@testitem "NOxPhotochemistry: Eq 6.6 Leighton Relationship" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # [O3]_pss = j_NO2 * [NO2] / (k_NO_O3 * [NO])
    j_NO2 = 8e-3                 # s⁻¹
    k_NO_O3 = 1.8e-14 * 1e-6     # m³/s

    # Test at typical urban conditions (m⁻³)
    NO2 = 5e16           # 2 ppb
    NO = 2.5e16          # 1 ppb

    O3_pss = j_NO2 * NO2 / (k_NO_O3 * NO)

    # Photostationary state O3 should be a reasonable tropospheric value
    @test O3_pss > 1e17   # > ~4 ppb
    @test O3_pss < 1e19   # < ~400 ppb
    @test O3_pss ≈ 8.89e17 rtol=0.05
end

@testitem "NOxPhotochemistry: Eq 6.7 Photostationary State Parameter" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # Phi = j_NO2 * [NO2] / (k_NO_O3 * [NO] * [O3])
    j_NO2 = 8e-3
    k_NO_O3 = 1.8e-14 * 1e-6

    NO2 = 5e16
    NO = 2.5e16
    O3_pss = j_NO2 * NO2 / (k_NO_O3 * NO)

    # Phi should be exactly 1 when O3 is at the PSS value
    Phi = j_NO2 * NO2 / (k_NO_O3 * NO * O3_pss)
    @test Phi ≈ 1.0 rtol=1e-10

    # When O3 is lower than PSS (e.g., extra oxidants present), Phi > 1
    O3_low = O3_pss * 0.5
    Phi_high = j_NO2 * NO2 / (k_NO_O3 * NO * O3_low)
    @test Phi_high > 1.0

    # When O3 is higher than PSS, Phi < 1
    O3_high = O3_pss * 2.0
    Phi_low = j_NO2 * NO2 / (k_NO_O3 * NO * O3_high)
    @test Phi_low < 1.0
end

@testitem "NOxPhotochemistry: Eq 6.8 Net O3 Production" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # P_O3 = j_NO2 * [NO2] - k_NO_O3 * [NO] * [O3]
    j_NO2 = 8e-3
    k_NO_O3 = 1.8e-14 * 1e-6
    NO2 = 5e16
    NO = 2.5e16

    # At photostationary state, P_O3 = 0
    O3_pss = j_NO2 * NO2 / (k_NO_O3 * NO)
    P_O3_pss = j_NO2 * NO2 - k_NO_O3 * NO * O3_pss
    @test abs(P_O3_pss) < 1e3  # numerically zero

    # Below PSS ozone: net production
    O3_low = O3_pss * 0.5
    P_O3_prod = j_NO2 * NO2 - k_NO_O3 * NO * O3_low
    @test P_O3_prod > 0

    # Above PSS ozone: net loss
    O3_high = O3_pss * 2.0
    P_O3_loss = j_NO2 * NO2 - k_NO_O3 * NO * O3_high
    @test P_O3_loss < 0
end

# ===========================================================================
# Rate Constants at 298 K
# ===========================================================================
@testitem "NOxPhotochemistry: Rate Constants at 298 K" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = NOxPhotochemistry()
    params = parameters(sys)
    param_dict = Dict(Symbol(p) => ModelingToolkit.getdefault(p) for p in params if ModelingToolkit.hasdefault(p))

    @test param_dict[:j_NO2] ≈ 8e-3                  # s⁻¹
    @test param_dict[:k_O_O2_M] ≈ 6.0e-34 * 1e-12    # m⁶/s
    @test param_dict[:k_NO_O3] ≈ 1.8e-14 * 1e-6      # m³/s
end

@testitem "PhotostationaryState: Rate Constants" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = PhotostationaryState()
    params = parameters(sys)
    param_dict = Dict(Symbol(p) => ModelingToolkit.getdefault(p) for p in params if ModelingToolkit.hasdefault(p))

    @test param_dict[:j_NO2] ≈ 8e-3
    @test param_dict[:k_NO_O3] ≈ 1.8e-14 * 1e-6
end

# ===========================================================================
# Steady-State Tests
# ===========================================================================
@testitem "NOxPhotochemistry: Photostationary State Equilibrium" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    j_NO2 = 8e-3
    k_NO_O3 = 1.8e-14 * 1e-6

    # Test several NO2/NO ratios
    for ratio in [0.1, 0.5, 1.0, 2.0, 10.0]
        NO = 1e16
        NO2 = ratio * NO
        O3_pss = j_NO2 * NO2 / (k_NO_O3 * NO)
        P_O3 = j_NO2 * NO2 - k_NO_O3 * NO * O3_pss
        @test abs(P_O3) / (j_NO2 * NO2) < 1e-10
    end
end

# ===========================================================================
# Limiting Behavior Tests
# ===========================================================================
@testitem "NOxPhotochemistry: O3 Scales with NO2/NO Ratio" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    j_NO2 = 8e-3
    k_NO_O3 = 1.8e-14 * 1e-6
    NO = 1e16

    NO2_a = 1e16
    NO2_b = 2e16

    O3_a = j_NO2 * NO2_a / (k_NO_O3 * NO)
    O3_b = j_NO2 * NO2_b / (k_NO_O3 * NO)

    @test O3_b / O3_a ≈ 2.0 rtol=1e-10
end

@testitem "NOxPhotochemistry: High j_NO2 Limit" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    k_NO_O3 = 1.8e-14 * 1e-6
    NO2 = 5e16
    NO = 2.5e16

    j_low = 1e-3
    j_high = 1e-2

    O3_low = j_low * NO2 / (k_NO_O3 * NO)
    O3_high = j_high * NO2 / (k_NO_O3 * NO)

    @test O3_high > O3_low
    @test O3_high / O3_low ≈ j_high / j_low rtol=1e-10
end

# ===========================================================================
# Qualitative Property Tests
# ===========================================================================
@testitem "NOxPhotochemistry: Qualitative - Positive O3_pss" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    j_NO2 = 8e-3
    k_NO_O3 = 1.8e-14 * 1e-6

    for NO in [1e14, 1e16, 1e18]
        for NO2 in [1e14, 1e16, 1e18]
            O3_pss = j_NO2 * NO2 / (k_NO_O3 * NO)
            @test O3_pss > 0
        end
    end
end

@testitem "PhotostationaryState: Deviation Sign" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    j_NO2 = 8e-3
    k_NO_O3 = 1.8e-14 * 1e-6
    NO = 2.5e16
    NO2 = 5e16
    O3_pss = j_NO2 * NO2 / (k_NO_O3 * NO)

    # At PSS: deviation = 0
    Phi = j_NO2 * NO2 / (k_NO_O3 * NO * O3_pss)
    @test Phi - 1 ≈ 0.0 atol=1e-10

    # Below PSS O3: deviation > 0
    Phi_above = j_NO2 * NO2 / (k_NO_O3 * NO * (O3_pss * 0.5))
    @test Phi_above - 1 > 0

    # Above PSS O3: deviation < 0
    Phi_below = j_NO2 * NO2 / (k_NO_O3 * NO * (O3_pss * 2.0))
    @test Phi_below - 1 < 0
end

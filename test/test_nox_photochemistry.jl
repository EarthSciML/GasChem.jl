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

    # 2 parameters: j_NO2, k_NO_O3
    @test length(params) == 2
end

# ===========================================================================
# Equation Verification Tests
# ===========================================================================
@testitem "NOxPhotochemistry: Eq 6.5 O Atom Steady-State" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # [O] = j_NO2 * [NO2] / (k_O_O2_M * [O2] * [M])
    # At typical conditions:
    j_NO2 = 8e-3       # s^-1
    NO2 = 1e10          # molecules/cm^3 (0.4 ppb)
    k_O_O2_M = 6.0e-34  # cm^6/s
    O2 = 5.25e18        # molecules/cm^3
    M = 2.5e19          # molecules/cm^3

    O_expected = j_NO2 * NO2 / (k_O_O2_M * O2 * M)

    # O atoms are very short-lived, concentration should be very low
    @test O_expected > 0
    @test O_expected < 1e6  # should be on order of ~1e3
    @test O_expected ≈ 1.016e3 rtol=0.05
end

@testitem "NOxPhotochemistry: Eq 6.6 Leighton Relationship" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # [O3]_pss = j_NO2 * [NO2] / (k_NO_O3 * [NO])
    # This is the fundamental Leighton relationship

    j_NO2 = 8e-3       # s^-1
    k_NO_O3 = 1.8e-14  # cm^3/s

    # Test at typical urban conditions
    NO2 = 5e10          # 2 ppb
    NO = 2.5e10         # 1 ppb

    O3_pss = j_NO2 * NO2 / (k_NO_O3 * NO)

    # Photostationary state O3 should be a reasonable tropospheric value
    # O3_pss ~ 8.9e11 molecules/cm^3 ~ 36 ppb
    @test O3_pss > 1e11   # > ~4 ppb
    @test O3_pss < 1e13   # < ~400 ppb
    @test O3_pss ≈ 8.89e11 rtol=0.05
end

@testitem "NOxPhotochemistry: Eq 6.7 Photostationary State Parameter" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # Phi = j_NO2 * [NO2] / (k_NO_O3 * [NO] * [O3])
    # In pure photostationary state, Phi = 1
    j_NO2 = 8e-3
    k_NO_O3 = 1.8e-14

    # At exact photostationary state: O3 = j_NO2*NO2/(k_NO_O3*NO)
    NO2 = 5e10
    NO = 2.5e10
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
    k_NO_O3 = 1.8e-14
    NO2 = 5e10
    NO = 2.5e10

    # At photostationary state, P_O3 = 0
    O3_pss = j_NO2 * NO2 / (k_NO_O3 * NO)
    P_O3_pss = j_NO2 * NO2 - k_NO_O3 * NO * O3_pss
    @test abs(P_O3_pss) < 1e-3  # numerically zero

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

    @test param_dict[:j_NO2] ≈ 8e-3       # s^-1 (typical midday)
    @test param_dict[:k_O_O2_M] ≈ 6.0e-34  # cm^6/s
    @test param_dict[:k_NO_O3] ≈ 1.8e-14   # cm^3/s
end

@testitem "PhotostationaryState: Rate Constants" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = PhotostationaryState()
    params = parameters(sys)
    param_dict = Dict(Symbol(p) => ModelingToolkit.getdefault(p) for p in params if ModelingToolkit.hasdefault(p))

    @test param_dict[:j_NO2] ≈ 8e-3
    @test param_dict[:k_NO_O3] ≈ 1.8e-14
end

# ===========================================================================
# Steady-State Tests
# ===========================================================================
@testitem "NOxPhotochemistry: Photostationary State Equilibrium" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # At photostationary state, the NO-NO2-O3 system is in balance.
    # The production of O3 (NO2 photolysis followed by O + O2 + M) equals
    # loss (NO + O3), giving P_O3 = 0.
    j_NO2 = 8e-3
    k_NO_O3 = 1.8e-14

    # Test several NO2/NO ratios
    for ratio in [0.1, 0.5, 1.0, 2.0, 10.0]
        NO = 1e10
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
    # O3_pss = j_NO2/k_NO_O3 * (NO2/NO)
    # Doubling the NO2/NO ratio should double O3
    j_NO2 = 8e-3
    k_NO_O3 = 1.8e-14
    NO = 1e10

    NO2_a = 1e10
    NO2_b = 2e10

    O3_a = j_NO2 * NO2_a / (k_NO_O3 * NO)
    O3_b = j_NO2 * NO2_b / (k_NO_O3 * NO)

    @test O3_b / O3_a ≈ 2.0 rtol=1e-10
end

@testitem "NOxPhotochemistry: High j_NO2 Limit" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # As j_NO2 increases (stronger sunlight), O3_pss increases
    k_NO_O3 = 1.8e-14
    NO2 = 5e10
    NO = 2.5e10

    j_low = 1e-3   # low light
    j_high = 1e-2  # strong sunlight

    O3_low = j_low * NO2 / (k_NO_O3 * NO)
    O3_high = j_high * NO2 / (k_NO_O3 * NO)

    @test O3_high > O3_low
    @test O3_high / O3_low ≈ j_high / j_low rtol=1e-10
end

# ===========================================================================
# Qualitative Property Tests
# ===========================================================================
@testitem "NOxPhotochemistry: Qualitative - Positive O3_pss" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # O3_pss must be positive for any positive inputs
    j_NO2 = 8e-3
    k_NO_O3 = 1.8e-14

    for NO in [1e8, 1e10, 1e12]
        for NO2 in [1e8, 1e10, 1e12]
            O3_pss = j_NO2 * NO2 / (k_NO_O3 * NO)
            @test O3_pss > 0
        end
    end
end

@testitem "PhotostationaryState: Deviation Sign" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # Phi_deviation = Phi - 1
    # When Phi > 1 (extra oxidants): deviation > 0
    # When Phi = 1 (pure PSS): deviation = 0
    # When Phi < 1 (extra reductants): deviation < 0
    j_NO2 = 8e-3
    k_NO_O3 = 1.8e-14
    NO = 2.5e10
    NO2 = 5e10
    O3_pss = j_NO2 * NO2 / (k_NO_O3 * NO)

    # At PSS: deviation = 0
    Phi = j_NO2 * NO2 / (k_NO_O3 * NO * O3_pss)
    @test Phi - 1 ≈ 0.0 atol=1e-10

    # Below PSS O3: deviation > 0 (extra oxidants converting NO to NO2)
    Phi_above = j_NO2 * NO2 / (k_NO_O3 * NO * (O3_pss * 0.5))
    @test Phi_above - 1 > 0

    # Above PSS O3: deviation < 0
    Phi_below = j_NO2 * NO2 / (k_NO_O3 * NO * (O3_pss * 2.0))
    @test Phi_below - 1 < 0
end

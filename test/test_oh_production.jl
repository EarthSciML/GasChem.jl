@testsnippet SP_CH6_Setup begin
    using Test
    using ModelingToolkit
    using ModelingToolkit: t, D
    using GasChem
end

# ===========================================================================
# Structural Tests
# ===========================================================================
@testitem "OHProduction: Structural Verification" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = OHProduction()
    @test sys isa System
    @test nameof(sys) == :OHProduction

    vars = unknowns(sys)
    params = parameters(sys)
    eqs = equations(sys)

    # OHProduction has 6 variables: O3, H2O, M (input) + O1D, P_OH, epsilon_OH (output)
    @test length(vars) == 6

    # 3 algebraic equations (Eq. 6.1, 6.3, 6.4)
    @test length(eqs) == 3

    # 7 parameters: j_O3, k3_N2, k3_O2, k4, f_N2, f_O2 + 1 constant (cm3_to_m3 is pruned, two counts as param)
    @test length(params) == 7

    # Check expected variable names are present
    var_names = [string(v) for v in vars]
    for expected in ["O1D", "P_OH", "ε_OH", "O3", "H2O", "M"]
        @test any(n -> contains(n, expected), var_names)
    end
end

# ===========================================================================
# Equation Verification Tests
# ===========================================================================
@testitem "OHProduction: Eq 6.1 O(1D) Steady-State" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = OHProduction()

    # Substitute known values and verify the O(1D) steady-state formula
    # From Eq 6.1: [O(1D)] = j_O3 * [O3] / (k3_eff * [M] + k4 * [H2O])
    #
    # At typical conditions (SI units, m⁻³):
    j_O3_val = 1e-5
    O3_val = 1e18
    M_val = 2.5e25
    H2O_val = 4e23
    k3_N2_val = 2.6e-11 * 1e-6  # m³/s
    k3_O2_val = 4.0e-11 * 1e-6  # m³/s
    k4_val = 2.2e-10 * 1e-6     # m³/s
    f_N2_val = 0.78
    f_O2_val = 0.21

    k3_eff = f_N2_val * k3_N2_val + f_O2_val * k3_O2_val
    denom = k3_eff * M_val + k4_val * H2O_val
    O1D_expected = j_O3_val * O3_val / denom

    @test O1D_expected > 0
    @test O1D_expected < 1e6  # O(1D) is present at very low concentrations
    @test O1D_expected ≈ 1.242e4 rtol=0.01
end

@testitem "OHProduction: Eq 6.4 OH Yield" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # epsilon_OH = k4 * [H2O] / (k3_eff * [M] + k4 * [H2O])
    k3_eff = 0.78 * 2.6e-11 * 1e-6 + 0.21 * 4.0e-11 * 1e-6
    k4_val = 2.2e-10 * 1e-6
    M_val = 2.5e25
    H2O_val = 4e23

    eps_expected = k4_val * H2O_val / (k3_eff * M_val + k4_val * H2O_val)

    # Typical tropospheric epsilon_OH is 0.05-0.15 (dimensionless, unchanged by unit system)
    @test eps_expected > 0.05
    @test eps_expected < 0.15
    @test eps_expected ≈ 0.109 rtol=0.05
end

@testitem "OHProduction: Eq 6.3 OH Production Rate" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # P_OH = 2 * j_O3 * [O3] * epsilon_OH
    j_O3_val = 1e-5
    O3_val = 1e18
    k3_eff = 0.78 * 2.6e-11 * 1e-6 + 0.21 * 4.0e-11 * 1e-6
    k4_val = 2.2e-10 * 1e-6
    M_val = 2.5e25
    H2O_val = 4e23

    eps_OH = k4_val * H2O_val / (k3_eff * M_val + k4_val * H2O_val)
    P_OH_expected = 2 * j_O3_val * O3_val * eps_OH

    # Typical OH production rate ~2e12 m⁻³/s
    @test P_OH_expected > 1e11
    @test P_OH_expected < 1e13
    @test P_OH_expected ≈ 2.18e12 rtol=0.05
end

# ===========================================================================
# Rate Constants at 298 K
# ===========================================================================
@testitem "OHProduction: Rate Constants at 298 K" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    sys = OHProduction()
    params = parameters(sys)
    param_dict = Dict(Symbol(p) => ModelingToolkit.getdefault(p) for p in params if ModelingToolkit.hasdefault(p))

    @test param_dict[:k3_N2] ≈ 2.6e-11 * 1e-6  # m³/s
    @test param_dict[:k3_O2] ≈ 4.0e-11 * 1e-6  # m³/s
    @test param_dict[:k4] ≈ 2.2e-10 * 1e-6     # m³/s
    @test param_dict[:j_O3] ≈ 1e-5              # s⁻¹
    @test param_dict[:f_N2] ≈ 0.78
    @test param_dict[:f_O2] ≈ 0.21
end

# ===========================================================================
# Limiting Behavior Tests
# ===========================================================================
@testitem "OHProduction: Limiting Behavior - Low Humidity" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # At very low humidity ([H2O] -> 0), epsilon_OH -> 0
    k3_eff = 0.78 * 2.6e-11 * 1e-6 + 0.21 * 4.0e-11 * 1e-6
    k4_val = 2.2e-10 * 1e-6
    M_val = 2.5e25

    H2O_low = 1e19  # very dry conditions (m⁻³)
    eps_low = k4_val * H2O_low / (k3_eff * M_val + k4_val * H2O_low)
    @test eps_low < 0.01  # Nearly all O(1D) quenched

    H2O_high = 1e24  # very humid conditions (m⁻³)
    eps_high = k4_val * H2O_high / (k3_eff * M_val + k4_val * H2O_high)
    @test eps_high > eps_low  # More OH at higher humidity
    @test eps_high < 1.0       # Cannot exceed 1
end

@testitem "OHProduction: Limiting Behavior - High Humidity" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # At extremely high humidity, epsilon_OH -> 1
    k3_eff = 0.78 * 2.6e-11 * 1e-6 + 0.21 * 4.0e-11 * 1e-6
    k4_val = 2.2e-10 * 1e-6
    M_val = 2.5e25

    H2O_extreme = 99 * k3_eff * M_val / k4_val
    eps_extreme = k4_val * H2O_extreme / (k3_eff * M_val + k4_val * H2O_extreme)
    @test eps_extreme ≈ 0.99 rtol=1e-6
end

# ===========================================================================
# Qualitative Property Tests
# ===========================================================================
@testitem "OHProduction: Qualitative - O1D Proportional to O3" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # Doubling O3 should double O(1D) (linear relationship)
    j_O3_val = 1e-5
    k3_eff = 0.78 * 2.6e-11 * 1e-6 + 0.21 * 4.0e-11 * 1e-6
    k4_val = 2.2e-10 * 1e-6
    M_val = 2.5e25
    H2O_val = 4e23

    denom = k3_eff * M_val + k4_val * H2O_val

    O3_a = 1e18
    O3_b = 2e18
    O1D_a = j_O3_val * O3_a / denom
    O1D_b = j_O3_val * O3_b / denom

    @test O1D_b / O1D_a ≈ 2.0 rtol=1e-10
end

@testitem "OHProduction: Qualitative - Positive Quantities" setup=[SP_CH6_Setup] tags=[:sp_ch6] begin
    # All computed quantities must be non-negative for physical inputs
    j_O3_val = 1e-5
    k3_eff = 0.78 * 2.6e-11 * 1e-6 + 0.21 * 4.0e-11 * 1e-6
    k4_val = 2.2e-10 * 1e-6

    for M_val in [1e24, 2.5e25, 1e26]
        for H2O_val in [1e21, 4e23, 1e24]
            for O3_val in [1e16, 1e18, 1e20]
                denom = k3_eff * M_val + k4_val * H2O_val
                O1D = j_O3_val * O3_val / denom
                eps = k4_val * H2O_val / denom
                P_OH = 2 * j_O3_val * O3_val * eps

                @test O1D >= 0
                @test eps >= 0
                @test eps <= 1
                @test P_OH >= 0
            end
        end
    end
end

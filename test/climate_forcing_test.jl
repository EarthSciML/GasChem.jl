@testsnippet ClimateSetup begin
    using Test
    using ModelingToolkit
    using NonlinearSolve
    using GasChem
    using DynamicQuantities
end

@testitem "ClimateFeedback structural verification" setup = [ClimateSetup] tags = [:climate] begin
    sys = ClimateFeedback()

    # Check number of equations
    eqs = equations(sys)
    @test length(eqs) == 6

    # Check that key variables exist
    vars = unknowns(sys)
    var_names = [string(Symbolics.tosymbol(v, escape = false)) for v in vars]
    @test "ΔT_s" in var_names
    @test "ΔT_0" in var_names
    @test "λ_0" in var_names
    @test "ΔF_e" in var_names
    @test "ΔT_unrealized" in var_names
    @test "ΔF_r" in var_names
end

@testitem "ClimateFeedback equation verification" setup = [ClimateSetup] tags = [:climate] begin
    # Test Eq. 23.1: ΔTs = λ ΔF
    # For λ = 0.8 K/(W/m²) and ΔF = 3.7 W/m² (2×CO₂), expect ΔTs ≈ 2.96 K
    sys = ClimateFeedback()
    csys = mtkcompile(sys)

    # Set parameters for 2×CO₂ scenario
    prob = NonlinearProblem(
        csys, [csys.λ => 0.8, csys.ΔF => 3.7];
        build_initializeprob = false
    )
    sol = solve(prob)

    ΔTs_computed = sol[csys.ΔT_s]
    @test isapprox(ΔTs_computed, 0.8 * 3.7, rtol = 1.0e-6)
end

@testitem "ClimateFeedback no-feedback response" setup = [ClimateSetup] tags = [:climate] begin
    # Test Eq. 23.2: ΔT₀ = λ₀ ΔF
    # For 2×CO₂ (ΔF ≈ 3.7 W/m²), no-feedback response should be ~1.2-1.3 K
    # With λ₀ = ΔT0_2xCO2/ΔF_2xCO2 = 1.25/3.7 ≈ 0.338 K/(W/m²)
    sys = ClimateFeedback()
    csys = mtkcompile(sys)

    prob = NonlinearProblem(csys, [csys.ΔF => 3.7]; build_initializeprob = false)
    sol = solve(prob)

    λ0_computed = sol[csys.λ_0]
    ΔT0_computed = sol[csys.ΔT_0]

    # λ₀ should be approximately 1.25/3.7 = 0.338
    @test isapprox(λ0_computed, 1.25 / 3.7, rtol = 1.0e-6)

    # ΔT₀ should be λ₀ × ΔF
    @test isapprox(ΔT0_computed, λ0_computed * 3.7, rtol = 1.0e-6)
end

@testitem "ClimateFeedback feedback factor" setup = [ClimateSetup] tags = [:climate] begin
    # Climate feedback factor = λ/λ₀ should be in range 1.2-3.75 (p. 1040)
    # Test that with λ = 0.8 and λ₀ ≈ 0.338, feedback factor ≈ 2.4
    sys = ClimateFeedback()
    csys = mtkcompile(sys)

    prob = NonlinearProblem(
        csys, [csys.λ => 0.8, csys.ΔF => 3.7];
        build_initializeprob = false
    )
    sol = solve(prob)

    λ_value = 0.8
    λ0_computed = sol[csys.λ_0]
    feedback_factor = λ_value / λ0_computed

    # Feedback factor should be in reasonable range
    @test 1.2 < feedback_factor < 3.75
end

@testitem "ClimateFeedback unrealized warming" setup = [ClimateSetup] tags = [:climate] begin
    # Test unrealized warming calculation (p. 1045)
    # Using values from p. 1045: ΔF ≈ 1.7 W/m², ΔTr ≈ 0.7°C, λ ≈ 0.7 K/(W/m²)
    # Expected unrealized warming ≈ 0.5°C

    sys = ClimateFeedback()
    csys = mtkcompile(sys)

    # Parameters from text
    λ_val = 0.7  # K/(W/m²), corresponds to ΔT_2xCO2 ≈ 2.6 K
    ΔF_val = 1.7  # W/m²
    ΔTr_val = 0.7  # K

    prob = NonlinearProblem(
        csys,
        [csys.λ => λ_val, csys.ΔF => ΔF_val, csys.ΔT_r => ΔTr_val];
        build_initializeprob = false
    )
    sol = solve(prob)

    ΔT_unrealized = sol[csys.ΔT_unrealized]

    # ΔT_unrealized = (ΔF - ΔF_r) × λ
    # where ΔF_r = ΔT_r / λ = 0.7 / 0.7 = 1.0 W/m²
    # So ΔT_unrealized = (1.7 - 1.0) × 0.7 = 0.49 K
    expected_unrealized = (ΔF_val - ΔTr_val / λ_val) * λ_val
    @test isapprox(ΔT_unrealized, expected_unrealized, rtol = 1.0e-6)

    # Should be approximately 0.5 K as stated in text
    @test isapprox(ΔT_unrealized, 0.49, rtol = 0.05)
end

@testitem "GHGForcing structural verification" setup = [ClimateSetup] tags = [:climate] begin
    sys = GHGForcing()

    eqs = equations(sys)
    @test length(eqs) == 6

    vars = unknowns(sys)
    var_names = [string(Symbolics.tosymbol(v, escape = false)) for v in vars]
    @test "ΔF_total" in var_names
    @test "ΔF_CO2" in var_names
    @test "ΔF_CH4" in var_names
end

@testitem "GHGForcing reference values" setup = [ClimateSetup] tags = [:climate] begin
    # Test that default forcings match IPCC 2001 values (p. 1039)
    sys = GHGForcing()
    csys = mtkcompile(sys)

    prob = NonlinearProblem(csys, Dict(); build_initializeprob = false)
    sol = solve(prob)

    # Individual forcings from p. 1039
    @test isapprox(sol[csys.ΔF_CO2], 1.46, rtol = 1.0e-6)
    @test isapprox(sol[csys.ΔF_CH4], 0.48, rtol = 1.0e-6)
    @test isapprox(sol[csys.ΔF_N2O], 0.15, rtol = 1.0e-6)
    @test isapprox(sol[csys.ΔF_O3], 0.4, rtol = 1.0e-6)

    # Total forcing (sum of well-mixed GHGs + O3)
    total = 1.46 + 0.48 + 0.15 + 0.4 + 0.34
    @test isapprox(sol[csys.ΔF_total], total, rtol = 1.0e-6)
end

@testitem "GlobalWarmingPotential structural verification" setup = [ClimateSetup] tags = [:climate] begin
    sys = GlobalWarmingPotential()

    eqs = equations(sys)
    @test length(eqs) == 5

    vars = unknowns(sys)
    var_names = [string(Symbolics.tosymbol(v, escape = false)) for v in vars]
    @test "GWP" in var_names
    @test "AGWP_A" in var_names
end

@testitem "GlobalWarmingPotential CO2 self-reference" setup = [ClimateSetup] tags = [:climate] begin
    # GWP of CO₂ relative to itself should be 1.0 for any time horizon
    sys = GlobalWarmingPotential()
    csys = mtkcompile(sys)

    yr_to_s = 365.25 * 24 * 3600  # seconds per year
    # CO₂: τ = 150 yr (effective), a = 1 (reference)
    for t_f_yr in [20.0, 100.0, 500.0]
        prob = NonlinearProblem(
            csys,
            [
                csys.τ_A => 150.0 * yr_to_s, csys.a_A => 1.0,
                csys.t_f => t_f_yr * yr_to_s,
            ];
            build_initializeprob = false
        )
        sol = solve(prob)
        @test isapprox(sol[csys.GWP], 1.0, rtol = 1.0e-6)
    end
end

@testitem "GlobalWarmingPotential CH4 estimate" setup = [ClimateSetup] tags = [:climate] begin
    # Test CH₄ GWP at 100-year horizon
    # From Table 23.1: CH₄ has τ = 12 yr (perturbation lifetime) and GWP₁₀₀ = 23

    # To get GWP = 23 with our formula, we need the radiative efficiency a_CH4
    # GWP = a_CH4 × τ_CH4 × (1 - exp(-100/12)) / (1 × 150 × (1 - exp(-100/150)))
    # 23 = a_CH4 × 12 × 0.9997 / (150 × 0.4866)
    # 23 = a_CH4 × 11.996 / 72.99
    # a_CH4 ≈ 140

    # Using the pure function to verify
    gwp_ch4 = GWP_exponential(12.0, 140.0, 100.0, τ_CO2 = 150.0)
    @test isapprox(gwp_ch4, 23.0, rtol = 0.1)  # Within 10% of tabulated value
end

@testitem "GlobalWarmingPotential time horizon dependence" setup = [ClimateSetup] tags = [:climate] begin
    # Test that GWP varies with time horizon as shown in Figure 23.15
    # For short-lived species, GWP decreases with longer time horizons

    # Short-lived species (τ = 2 yr) should have decreasing GWP with time
    τ_short = 2.0
    a_short = 100.0

    gwp_20 = GWP_exponential(τ_short, a_short, 20.0)
    gwp_100 = GWP_exponential(τ_short, a_short, 100.0)
    gwp_500 = GWP_exponential(τ_short, a_short, 500.0)

    @test gwp_20 > gwp_100
    @test gwp_100 > gwp_500

    # Long-lived species (τ = 10000 yr) should have increasing GWP with time
    τ_long = 10000.0
    a_long = 100.0

    gwp_20_long = GWP_exponential(τ_long, a_long, 20.0)
    gwp_100_long = GWP_exponential(τ_long, a_long, 100.0)
    gwp_500_long = GWP_exponential(τ_long, a_long, 500.0)

    @test gwp_20_long < gwp_100_long
    @test gwp_100_long < gwp_500_long
end

@testitem "GWP_exponential function" setup = [ClimateSetup] tags = [:climate] begin
    # Test the pure function implementation

    # For τ_A = τ_CO2 and a_A = 1, GWP should equal 1
    @test isapprox(GWP_exponential(150.0, 1.0, 100.0), 1.0, rtol = 1.0e-10)

    # For very short-lived species (τ → 0), GWP → 0
    @test GWP_exponential(0.001, 1.0, 100.0) < 0.01

    # For very long-lived species (τ → ∞), GWP → a_A × t_f / (τ_CO2 × (1 - exp(-t_f/τ_CO2)))
    # As τ → ∞: τ(1 - exp(-t_f/τ)) → t_f
    gwp_long = GWP_exponential(1.0e6, 1.0, 100.0)
    expected_long = 100.0 / (150.0 * (1 - exp(-100.0 / 150.0)))
    @test isapprox(gwp_long, expected_long, rtol = 0.01)
end

@testitem "N2O GWP estimate" setup = [ClimateSetup] tags = [:climate] begin
    # From Table 23.1: N₂O has τ = 114 yr (perturbation lifetime) and GWP₁₀₀ = 296
    # Back-calculate radiative efficiency: need a_N2O such that GWP = 296
    # 296 = a_N2O × 114 × (1 - exp(-100/114)) / (150 × (1 - exp(-100/150)))
    # 296 = a_N2O × 114 × 0.582 / (150 × 0.487)
    # 296 = a_N2O × 66.35 / 73.05
    # a_N2O ≈ 326

    gwp_n2o = GWP_exponential(114.0, 326.0, 100.0, τ_CO2 = 150.0)
    @test isapprox(gwp_n2o, 296.0, rtol = 0.1)  # Within 10% of tabulated value
end

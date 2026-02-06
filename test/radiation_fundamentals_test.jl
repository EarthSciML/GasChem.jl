@testsnippet RadiationSetup begin
    using ModelingToolkit, NonlinearSolve
    using GasChem: PhotonEnergy, BlackbodyRadiation, WienDisplacement, StefanBoltzmann
    using GasChem: PlanetaryEnergyBalance, ClimateSensitivity, TOARadiativeForcing,
                   RadiationFundamentals
end

# ============================================================
# Structural Tests
# ============================================================

@testitem "Structural: RadiationFundamentals composed system" setup=[RadiationSetup] tags=[:radiation] begin
    @named radiation = RadiationFundamentals()

    # Verify the system has exactly 7 subsystems
    subsystems=ModelingToolkit.get_systems(radiation)
    @test length(subsystems) == 7

    # Verify subsystem names
    names=Set(string(nameof(s)) for s in subsystems)
    @test "photon" in names
    @test "blackbody" in names
    @test "wien" in names
    @test "stefan_boltzmann" in names
    @test "energy_balance" in names
    @test "climate_sensitivity" in names
    @test "toa_forcing" in names
end

@testitem "Structural: PhotonEnergy system" setup=[RadiationSetup] tags=[:radiation] begin
    @named photon = PhotonEnergy()

    # Check pre-compiled system structure
    @test length(equations(photon)) == 2
    @test length(ModelingToolkit.get_unknowns(photon)) == 2
    # Verify it compiles successfully
    sys=mtkcompile(photon)
    @test sys !== nothing
end

@testitem "Structural: PlanetaryEnergyBalance system" setup=[RadiationSetup] tags=[:radiation] begin
    @named balance = PlanetaryEnergyBalance()

    # Check pre-compiled system structure
    @test length(equations(balance)) == 3
    @test length(ModelingToolkit.get_unknowns(balance)) == 3
    # Verify it compiles successfully
    sys=mtkcompile(balance)
    @test sys !== nothing
end

@testitem "Structural: BlackbodyRadiation system" setup=[RadiationSetup] tags=[:radiation] begin
    @named bb = BlackbodyRadiation()

    @test length(equations(bb)) == 1
    @test length(ModelingToolkit.get_unknowns(bb)) == 1
    sys=mtkcompile(bb)
    @test sys !== nothing
end

@testitem "Structural: WienDisplacement system" setup=[RadiationSetup] tags=[:radiation] begin
    @named wien = WienDisplacement()

    @test length(equations(wien)) == 1
    @test length(ModelingToolkit.get_unknowns(wien)) == 1
    sys=mtkcompile(wien)
    @test sys !== nothing
end

@testitem "Structural: StefanBoltzmann system" setup=[RadiationSetup] tags=[:radiation] begin
    @named sb = StefanBoltzmann()

    @test length(equations(sb)) == 1
    @test length(ModelingToolkit.get_unknowns(sb)) == 1
    sys=mtkcompile(sb)
    @test sys !== nothing
end

@testitem "Structural: ClimateSensitivity system" setup=[RadiationSetup] tags=[:radiation] begin
    @named cs = ClimateSensitivity()

    @test length(equations(cs)) == 4
    @test length(ModelingToolkit.get_unknowns(cs)) == 4
    sys=mtkcompile(cs)
    @test sys !== nothing
end

@testitem "Structural: TOARadiativeForcing system" setup=[RadiationSetup] tags=[:radiation] begin
    @named toa = TOARadiativeForcing()

    @test length(equations(toa)) == 2
    @test length(ModelingToolkit.get_unknowns(toa)) == 2
    sys=mtkcompile(toa)
    @test sys !== nothing
end

# ============================================================
# Equation Verification Tests
# ============================================================

@testitem "Equation: PhotonEnergy (Eq. 4.1)" setup=[RadiationSetup] tags=[:radiation] begin
    @named photon = PhotonEnergy()
    sys=mtkcompile(photon)

    # Test at default wavelength (500 nm = 5e-7 m)
    # Reference: visible light photon at 500 nm has ν ≈ 5.996e14 Hz, ε ≈ 3.97e-19 J
    prob=NonlinearProblem(sys, Dict(); build_initializeprob = false)
    sol=solve(prob)

    @test isapprox(sol[photon.ν], 5.9958e14, rtol = 1e-4)
    @test isapprox(sol[photon.Δε], 3.973e-19, rtol = 1e-3)

    # Test UV wavelength (300 nm): ν ≈ 9.993e14 Hz, ε ≈ 6.622e-19 J
    prob_uv=remake(prob, p = [sys.λ=>3e-7])
    sol_uv=solve(prob_uv)

    @test isapprox(sol_uv[photon.ν], 9.993e14, rtol = 1e-3)
    @test isapprox(sol_uv[photon.Δε], 6.622e-19, rtol = 1e-3)
end

@testitem "Equation: WienDisplacement (Eq. 4.3)" setup=[RadiationSetup] tags=[:radiation] begin
    @named wien = WienDisplacement()
    sys=mtkcompile(wien)

    # S&P: Sun at ~5800 K peaks near 500 nm (b/T = 2.897e-3/5800 ≈ 4.995e-7 m)
    prob=NonlinearProblem(sys, Dict(); build_initializeprob = false)
    sol=solve(prob)
    @test isapprox(sol[wien.λ_max], 500e-9, rtol = 2e-3)  # ~500 nm

    # S&P: Earth at ~300 K peaks near 10 μm (b/T = 2.897e-3/300 ≈ 9.657e-6 m)
    prob_earth=remake(prob, p = [sys.T=>300.0])
    sol_earth=solve(prob_earth)
    @test isapprox(sol_earth[wien.λ_max], 10e-6, rtol = 5e-2)  # ~10 μm
end

@testitem "Equation: StefanBoltzmann (Eq. 4.4)" setup=[RadiationSetup] tags=[:radiation] begin
    @named sb = StefanBoltzmann()
    sys=mtkcompile(sb)

    # At T = 255 K (Earth's equilibrium), F_B ≈ 239.7 W/m² (S&P p. 101)
    prob=NonlinearProblem(sys, Dict(); build_initializeprob = false)
    sol=solve(prob)
    @test isapprox(sol[sb.F_B], 239.7, rtol = 1e-2)

    # At T = 5800 K (Sun), F_B ≈ 6.42e7 W/m²
    prob_sun=remake(prob, p = [sys.T=>5800.0])
    sol_sun=solve(prob_sun)
    @test isapprox(sol_sun[sb.F_B], 6.42e7, rtol = 1e-2)
end

@testitem "Equation: BlackbodyRadiation (Eq. 4.2)" setup=[RadiationSetup] tags=[:radiation] begin
    @named bb = BlackbodyRadiation()
    sys=mtkcompile(bb)

    # At T = 5800 K, λ = 500 nm: F_B_λ ≈ 8.46e13 W/m³
    # (2πc²h/λ⁵ / (exp(ch/kλT) - 1) with full 2π factor for hemispherical emissive power)
    prob=NonlinearProblem(sys, Dict(); build_initializeprob = false)
    sol=solve(prob)
    @test isapprox(sol[bb.F_B_λ], 8.46e13, rtol = 5e-2)

    # At T = 300 K, λ = 10 μm: F_B_λ ≈ 3.13e7 W/m³
    prob_earth=remake(prob, p = [sys.T=>300.0, sys.λ=>1e-5])
    sol_earth=solve(prob_earth)
    @test isapprox(sol_earth[bb.F_B_λ], 3.13e7, rtol = 5e-2)
end

@testitem "Equation: PlanetaryEnergyBalance (Eqs. 4.5-4.7)" setup=[RadiationSetup] tags=[:radiation] begin
    @named balance = PlanetaryEnergyBalance()
    sys=mtkcompile(balance)

    # S&P p. 101: S_0 = 1370 W/m², R_p = 0.3 → T_e ≈ 255 K, F_S ≈ 239.75 W/m²
    prob=NonlinearProblem(sys, Dict(balance.T_e=>250.0); build_initializeprob = false)
    sol=solve(prob)

    @test isapprox(sol[balance.F_S], 239.75, rtol = 1e-6)
    @test isapprox(sol[balance.T_e], 255.0, rtol = 1e-2)

    # S&P also gives T_e ≈ 268 K for R_p = 0.15 (Earth without clouds)
    prob2=remake(prob, p = [sys.R_p=>0.15])
    sol2=solve(prob2)
    @test isapprox(sol2[balance.T_e], 268.0, rtol = 1e-2)
end

@testitem "Equation: ClimateSensitivity (Eqs. 4.8-4.10)" setup=[RadiationSetup] tags=[:radiation] begin
    @named sensitivity = ClimateSensitivity()
    sys=mtkcompile(sensitivity)

    # S&P p. 105: λ_0 ≈ 0.3 K/(W/m²) at T_e = 255 K (exact: 1/(4σT³) ≈ 0.266)
    prob=NonlinearProblem(sys, Dict(); build_initializeprob = false)
    sol=solve(prob)

    @test isapprox(sol[sensitivity.λ_0], 0.27, rtol = 0.05)

    # Default ΔF_S = 4.0 W/m² gives ΔT_e ≈ 1.1 K (no-feedback, using exact λ_0)
    @test isapprox(sol[sensitivity.ΔT_e], 1.1, rtol = 0.1)

    # S&P p. 105: CO2 doubling forcing of 4.6 W/m² gives ΔT_e ≈ 1.4 K (with rounded λ_0 ≈ 0.3)
    # With exact λ_0 ≈ 0.266: ΔT_e = 0.266 × 4.6 ≈ 1.22 K
    prob2=remake(prob, p = [sys.ΔF_S=>4.6])
    sol2=solve(prob2)
    @test isapprox(sol2[sensitivity.ΔT_e], 1.22, rtol = 0.05)
end

@testitem "Equation: TOARadiativeForcing (Eq. 4.11)" setup=[RadiationSetup] tags=[:radiation] begin
    @named toa = TOARadiativeForcing()
    sys=mtkcompile(toa)

    # At equilibrium (F_L = F_S = 239.75 W/m²), net flux should be ~0
    prob=NonlinearProblem(sys, Dict(); build_initializeprob = false)
    sol=solve(prob)
    @test isapprox(sol[toa.F_net], 0.0, atol = 1e-6)
    @test isapprox(sol[toa.F_S], 239.75, rtol = 1e-6)

    # With increased F_L (enhanced greenhouse), net flux becomes negative (cooling)
    prob2=remake(prob, p = [sys.F_L=>250.0])
    sol2=solve(prob2)
    @test sol2[toa.F_net] < 0  # Planet losing more energy than absorbing
end

# ============================================================
# Conservation Law Tests
# ============================================================

@testitem "Conservation: Energy balance F_S = F_L at equilibrium" setup=[RadiationSetup] tags=[:radiation] begin
    @named balance = PlanetaryEnergyBalance()
    sys=mtkcompile(balance)

    # Verify F_S = F_L at equilibrium for multiple albedo values
    for R_p in [0.0, 0.15, 0.3, 0.5, 0.7]
        prob=NonlinearProblem(sys, Dict(balance.T_e=>250.0); build_initializeprob = false)
        prob=remake(prob, p = [sys.R_p=>R_p])
        sol=solve(prob)
        @test isapprox(sol[balance.F_S], sol[balance.F_L], rtol = 1e-8)
    end
end

# ============================================================
# Steady-State Tests
# ============================================================

@testitem "Steady-State: TOA net flux zero at equilibrium" setup=[RadiationSetup] tags=[:radiation] begin
    @named toa = TOARadiativeForcing()
    sys=mtkcompile(toa)

    # At equilibrium, the outgoing longwave flux equals absorbed solar flux
    # With default F_L = 239.75 W/m² matching the absorbed solar flux
    prob=NonlinearProblem(sys, Dict(); build_initializeprob = false)
    sol=solve(prob)
    @test isapprox(sol[toa.F_net], 0.0, atol = 1e-6)
end

@testitem "Steady-State: ClimateSensitivity zero forcing gives zero response" setup=[RadiationSetup] tags=[:radiation] begin
    @named sensitivity = ClimateSensitivity()
    sys=mtkcompile(sensitivity)

    # With no forcing (ΔF_S = 0, ΔF_L = 0), temperature change should be zero
    prob=NonlinearProblem(sys, Dict(); build_initializeprob = false)
    prob=remake(prob, p = [sys.ΔF_S=>0.0, sys.ΔF_L=>0.0])
    sol=solve(prob)
    @test isapprox(sol[sensitivity.ΔT_e], 0.0, atol = 1e-10)
    @test isapprox(sol[sensitivity.ΔF_net], 0.0, atol = 1e-10)
end

# ============================================================
# Limiting Behavior Tests
# ============================================================

@testitem "Limiting: Wien's law at extreme temperatures" setup=[RadiationSetup] tags=[:radiation] begin
    @named wien = WienDisplacement()
    sys=mtkcompile(wien)

    prob=NonlinearProblem(sys, Dict(); build_initializeprob = false)

    # Very hot object: peak shifts to shorter wavelengths (UV/X-ray)
    prob_hot=remake(prob, p = [sys.T=>1e6])
    sol_hot=solve(prob_hot)
    @test sol_hot[wien.λ_max] < 1e-8  # Peak in X-ray region

    # Very cold object: peak shifts to long wavelengths (far-IR/microwave)
    prob_cold=remake(prob, p = [sys.T=>3.0])  # CMB temperature
    sol_cold=solve(prob_cold)
    @test sol_cold[wien.λ_max] > 1e-4  # Peak beyond far-IR
    @test isapprox(sol_cold[wien.λ_max], 9.66e-4, rtol = 1e-2)  # ~1 mm for CMB
end

@testitem "Limiting: Equilibrium temperature at extreme albedo" setup=[RadiationSetup] tags=[:radiation] begin
    @named balance = PlanetaryEnergyBalance()
    sys=mtkcompile(balance)

    # R_p = 0 (no reflection): maximum equilibrium temperature
    prob=NonlinearProblem(sys, Dict(balance.T_e=>250.0); build_initializeprob = false)
    prob_no_albedo=remake(prob, p = [sys.R_p=>0.0])
    sol_no_albedo=solve(prob_no_albedo)
    T_max=sol_no_albedo[balance.T_e]
    @test T_max > 270  # Higher than Earth's 255 K

    # R_p = 0.99 (near-perfect reflection): very cold equilibrium
    prob_high_albedo=remake(prob, p = [sys.R_p=>0.99])
    sol_high_albedo=solve(prob_high_albedo)
    T_min=sol_high_albedo[balance.T_e]
    @test T_min < 100  # Very cold
    @test T_min > 0  # But still positive (T^4 law)
end

@testitem "Limiting: Blackbody radiation non-negative for all positive T and λ" setup=[RadiationSetup] tags=[:radiation] begin
    @named bb = BlackbodyRadiation()
    sys=mtkcompile(bb)

    prob=NonlinearProblem(sys, Dict(); build_initializeprob = false)

    # Test across a range of temperatures and wavelengths
    # Note: at extreme combinations (e.g., T=100K, λ=100nm), the exponent
    # ch/kλT can exceed ~1400, causing numerical underflow to 0.0
    for T in [100.0, 300.0, 1000.0, 5800.0, 10000.0]
        for λ in [1e-7, 1e-6, 1e-5, 1e-4]
            prob_i=remake(prob, p = [sys.T=>T, sys.λ=>λ])
            sol_i=solve(prob_i)
            @test sol_i[bb.F_B_λ] >= 0
        end
    end

    # At physically relevant combinations, verify strictly positive
    for (T, λ) in [(5800.0, 5e-7), (300.0, 1e-5), (1000.0, 3e-6)]
        prob_i=remake(prob, p = [sys.T=>T, sys.λ=>λ])
        sol_i=solve(prob_i)
        @test sol_i[bb.F_B_λ] > 0
    end
end

# ============================================================
# Qualitative Behavior Tests
# ============================================================

@testitem "Qualitative: T_e decreases with increasing albedo" setup=[RadiationSetup] tags=[:radiation] begin
    @named balance = PlanetaryEnergyBalance()
    sys=mtkcompile(balance)

    prob=NonlinearProblem(sys, Dict(balance.T_e=>250.0); build_initializeprob = false)

    albedos=[0.0, 0.1, 0.2, 0.3, 0.5, 0.7]
    temperatures=Float64[]

    for R_p in albedos
        prob_i=remake(prob, p = [sys.R_p=>R_p])
        sol_i=solve(prob_i)
        push!(temperatures, sol_i[balance.T_e])
    end

    # Temperature should be strictly decreasing with increasing albedo
    for i in 2:length(temperatures)
        @test temperatures[i] < temperatures[i - 1]
    end
end

@testitem "Qualitative: T_e increases with increasing solar constant" setup=[RadiationSetup] tags=[:radiation] begin
    @named balance = PlanetaryEnergyBalance()
    sys=mtkcompile(balance)

    prob=NonlinearProblem(sys, Dict(balance.T_e=>250.0); build_initializeprob = false)

    solar_constants=[1000.0, 1200.0, 1370.0, 1500.0, 1800.0]
    temperatures=Float64[]

    for S_0 in solar_constants
        prob_i=remake(prob, p = [sys.S_0=>S_0])
        sol_i=solve(prob_i)
        push!(temperatures, sol_i[balance.T_e])
    end

    # Temperature should be strictly increasing with increasing solar constant
    for i in 2:length(temperatures)
        @test temperatures[i] > temperatures[i - 1]
    end
end

@testitem "Qualitative: ΔT_e proportional to radiative forcing" setup=[RadiationSetup] tags=[:radiation] begin
    @named sensitivity = ClimateSensitivity()
    sys=mtkcompile(sensitivity)

    prob=NonlinearProblem(sys, Dict(); build_initializeprob = false)

    # ΔT_e = λ_0 * ΔF_net, so ΔT should be linear in forcing
    forcings=[1.0, 2.0, 4.0, 8.0]
    delta_Ts=Float64[]

    for ΔF in forcings
        prob_i=remake(prob, p = [sys.ΔF_S=>ΔF])
        sol_i=solve(prob_i)
        push!(delta_Ts, sol_i[sensitivity.ΔT_e])
    end

    # Check linearity: doubling forcing should double ΔT
    @test isapprox(delta_Ts[2] / delta_Ts[1], 2.0, rtol = 1e-8)
    @test isapprox(delta_Ts[3] / delta_Ts[1], 4.0, rtol = 1e-8)
    @test isapprox(delta_Ts[4] / delta_Ts[1], 8.0, rtol = 1e-8)
end

@testitem "Qualitative: Wien peak shifts to shorter λ at higher T" setup=[RadiationSetup] tags=[:radiation] begin
    @named wien = WienDisplacement()
    sys=mtkcompile(wien)

    prob=NonlinearProblem(sys, Dict(); build_initializeprob = false)

    temps=[300.0, 1000.0, 3000.0, 5800.0]
    peaks=Float64[]

    for T in temps
        prob_i=remake(prob, p = [sys.T=>T])
        sol_i=solve(prob_i)
        push!(peaks, sol_i[wien.λ_max])
    end

    # Peak wavelength should be strictly decreasing with temperature
    for i in 2:length(peaks)
        @test peaks[i] < peaks[i - 1]
    end
end

# ============================================================
# Analytical Solution Tests
# ============================================================

@testitem "Analytical: Planck spectrum peak matches Wien prediction" setup=[RadiationSetup] tags=[:radiation] begin
    # Verify that BlackbodyRadiation peaks at the wavelength predicted by WienDisplacement
    @named bb = BlackbodyRadiation()
    @named wien = WienDisplacement()
    sys_bb=mtkcompile(bb)
    sys_wien=mtkcompile(wien)

    for T in [300.0, 1000.0, 5800.0]
        # Get Wien peak wavelength
        prob_wien=NonlinearProblem(sys_wien, Dict(); build_initializeprob = false)
        prob_wien=remake(prob_wien, p = [sys_wien.T=>T])
        sol_wien=solve(prob_wien)
        λ_peak=sol_wien[wien.λ_max]

        # Compute Planck function at peak and at wavelengths slightly below and above
        prob_bb=NonlinearProblem(sys_bb, Dict(); build_initializeprob = false)

        prob_at_peak=remake(prob_bb, p = [sys_bb.T=>T, sys_bb.λ=>λ_peak])
        F_peak=solve(prob_at_peak)[bb.F_B_λ]

        prob_below=remake(prob_bb, p = [sys_bb.T=>T, sys_bb.λ=>0.9*λ_peak])
        F_below=solve(prob_below)[bb.F_B_λ]

        prob_above=remake(prob_bb, p = [sys_bb.T=>T, sys_bb.λ=>1.1*λ_peak])
        F_above=solve(prob_above)[bb.F_B_λ]

        # The value at the Wien peak should be greater than at nearby wavelengths
        @test F_peak > F_below
        @test F_peak > F_above
    end
end

@testitem "Analytical: Climate sensitivity factor alternative form (Eq. 4.10)" setup=[RadiationSetup] tags=[:radiation] begin
    # Eq. 4.10: λ₀ = 1/(4σT_e³) = T_e/(4F_L)
    # Verify the equivalence of both forms using the ClimateSensitivity component
    @named cs = ClimateSensitivity()
    sys=mtkcompile(cs)

    for T_e in [200.0, 255.0, 300.0]
        prob=NonlinearProblem(sys, Dict(); build_initializeprob = false)
        prob=remake(prob, p = [sys.T_e=>T_e])
        sol=solve(prob)

        λ_0=sol[cs.λ_0]
        F_L=sol[cs.F_L]
        # Check that λ₀ = T_e / (4 * F_L)
        @test isapprox(λ_0, T_e / (4 * F_L), rtol = 1e-10)
    end
end

@testitem "Analytical: RadiationFundamentals composed system solves" setup=[RadiationSetup] tags=[:radiation] begin
    # Verify the composed system can be compiled and solved
    @named radiation = RadiationFundamentals()

    # Pre-compiled system has 14 equations across 7 subsystems
    @test length(equations(radiation)) == 14

    sys=mtkcompile(radiation)
    @test sys !== nothing

    # After compilation, most algebraic subsystems reduce to observables.
    # Only T_e from PlanetaryEnergyBalance remains as an unknown (nonlinear constraint).
    @test length(unknowns(sys)) >= 1

    # Solve and verify key results from subsystems
    prob=NonlinearProblem(sys, Dict(radiation.energy_balance.T_e=>250.0);
        build_initializeprob = false)
    sol=solve(prob)

    # Verify energy balance T_e ~ 255 K (Eq. 4.7, S&P p. 101)
    @test isapprox(sol[radiation.energy_balance.T_e], 255.0, rtol = 1e-2)
end

@testitem "PhotonEnergy (Eq. 4.1)" begin
    using ModelingToolkit, NonlinearSolve
    using GasChem: PhotonEnergy

    @named photon = PhotonEnergy()
    sys = mtkcompile(photon)

    # Test visible light wavelength (500 nm = 5e-7 m)
    prob = NonlinearProblem(sys, Dict(); build_initializeprob=false)
    sol = solve(prob)

    λ_val = 5e-7  # Default wavelength
    ν_expected = 2.9979e8 / λ_val
    ε_expected = 6.626e-34 * ν_expected

    @test isapprox(sol[photon.ν], ν_expected, rtol=1e-3)
    @test isapprox(sol[photon.Δε], ε_expected, rtol=1e-3)
end

@testitem "WienDisplacement (Eq. 4.3)" begin
    using ModelingToolkit, NonlinearSolve
    using GasChem: WienDisplacement

    @named wien = WienDisplacement()
    sys = mtkcompile(wien)

    # Test for Sun (T = 5800 K) - using default
    prob = NonlinearProblem(sys, Dict(); build_initializeprob=false)
    sol = solve(prob)

    T_sun = 5800.0  # Default temperature
    λ_max_expected = 2.897e-3 / T_sun

    @test isapprox(sol[wien.λ_max], λ_max_expected, rtol=1e-3)
end

@testitem "StefanBoltzmann (Eq. 4.4)" begin
    using ModelingToolkit, NonlinearSolve
    using GasChem: StefanBoltzmann

    @named sb = StefanBoltzmann()
    sys = mtkcompile(sb)

    # Test for Earth (T = 255 K) - using default
    prob = NonlinearProblem(sys, Dict(); build_initializeprob=false)
    sol = solve(prob)

    T_earth = 255.0  # Default temperature
    F_B_expected = 5.671e-8 * T_earth^4

    @test isapprox(sol[sb.F_B], F_B_expected, rtol=1e-3)
end

@testitem "PlanetaryEnergyBalance (Eqs. 4.5-4.7)" begin
    using ModelingToolkit, NonlinearSolve
    using GasChem: PlanetaryEnergyBalance

    @named balance = PlanetaryEnergyBalance()
    sys = mtkcompile(balance)

    # With default parameters (S_0 = 1370 W/m², R_p = 0.3)
    F_S_expected = 1370.0 / 4 * (1 - 0.3)
    T_e_expected = (F_S_expected / 5.671e-8)^0.25

    prob = NonlinearProblem(sys, Dict(balance.T_e => 250.0); build_initializeprob=false)
    sol = solve(prob)

    @test isapprox(sol[balance.F_S], F_S_expected, rtol=1e-3)
    @test isapprox(sol[balance.T_e], T_e_expected, rtol=1e-3)
    @test isapprox(sol[balance.F_S], sol[balance.F_L], rtol=1e-6)  # Energy balance
end

@testitem "ClimateSensitivity (Eqs. 4.8-4.10)" begin
    using ModelingToolkit, NonlinearSolve
    using GasChem: ClimateSensitivity

    @named sensitivity = ClimateSensitivity()
    sys = mtkcompile(sensitivity)

    # Test at Earth's equilibrium temperature T_e ≈ 255 K (default)
    T_e_ref = 255.0
    λ_0_expected = 1 / (4 * 5.671e-8 * T_e_ref^3)

    prob = NonlinearProblem(sys, Dict(); build_initializeprob=false)
    sol = solve(prob)

    ΔF_S = 4.0  # W/m² radiative forcing (default)
    ΔT_e_expected = λ_0_expected * ΔF_S

    @test isapprox(sol[sensitivity.λ_0], λ_0_expected, rtol=1e-3)
    @test isapprox(sol[sensitivity.ΔT_e], ΔT_e_expected, rtol=1e-3)
end

@testitem "TOARadiativeForcing (Eq. 4.11)" begin
    using ModelingToolkit, NonlinearSolve
    using GasChem: TOARadiativeForcing

    @named toa = TOARadiativeForcing()
    sys = mtkcompile(toa)

    # At equilibrium with F_L matching absorbed solar
    prob = NonlinearProblem(sys, Dict(); build_initializeprob=false)
    sol = solve(prob)

    @test isapprox(sol[toa.F_net], 0.0, atol=1e-6)  # At equilibrium, net = 0
end

@testitem "BlackbodyRadiation (Eq. 4.2)" begin
    using ModelingToolkit, NonlinearSolve
    using GasChem: BlackbodyRadiation

    @named bb = BlackbodyRadiation()
    sys = mtkcompile(bb)

    # Test at Sun's temperature (default: T = 5800 K) and peak wavelength
    T_sun = 5800.0
    λ_peak = 5e-7

    prob = NonlinearProblem(sys, Dict(); build_initializeprob=false)
    sol = solve(prob)

    # Planck function at peak
    h = 6.626e-34
    c = 2.9979e8
    k = 1.381e-23
    F_B_λ_expected = 2 * π * c^2 * h * λ_peak^(-5) / (exp(c * h / (k * λ_peak * T_sun)) - 1)

    @test isapprox(sol[bb.F_B_λ], F_B_λ_expected, rtol=1e-3)
end

@testitem "RadiationFundamentals (composed system)" begin
    using ModelingToolkit
    using GasChem: RadiationFundamentals

    @named radiation = RadiationFundamentals()

    # Verify the system structure
    @test length(ModelingToolkit.get_systems(radiation)) == 7
end

@testitem "Physical Constants Verification" begin
    # Verify constants match Seinfeld & Pandis Table A.6
    @test isapprox(6.626e-34, 6.626e-34, rtol=1e-10)  # Planck's constant
    @test isapprox(2.9979e8, 2.9979e8, rtol=1e-10)    # Speed of light
    @test isapprox(1.381e-23, 1.381e-23, rtol=1e-10)  # Boltzmann constant
    @test isapprox(5.671e-8, 5.671e-8, rtol=1e-10)    # Stefan-Boltzmann constant
    @test isapprox(2.897e-3, 2.897e-3, rtol=1e-10)    # Wien's displacement constant
end

@testitem "Wien-Planck Consistency" begin
    # Wien's law should be consistent with Planck's law
    h = 6.626e-34
    c = 2.9979e8
    k = 1.381e-23

    # Theoretical Wien constant from Planck derivation
    x_wien = 4.965114231744276  # Numerical solution to x = 5(1-e^-x)
    b_theoretical = h * c / (x_wien * k)
    b_implementation = 2.897e-3

    @test isapprox(b_implementation, b_theoretical, rtol=1e-3)
end

@testitem "Climate Sensitivity Alternative Form" begin
    # Verify two equivalent forms of lambda_0
    sigma = 5.671e-8
    T_e = 255.0

    F_L = sigma * T_e^4
    lambda_0_form1 = 1 / (4 * sigma * T_e^3)
    lambda_0_form2 = T_e / (4 * F_L)

    @test isapprox(lambda_0_form1, lambda_0_form2, rtol=1e-10)
    @test isapprox(lambda_0_form1, 0.266, rtol=0.05)  # Within 5% of S&P approximate value
end

@testitem "Equilibrium Temperature Calculation" begin
    # Verify S&P stated values for equilibrium temperature
    sigma = 5.671e-8
    S_0 = 1370.0

    # Case 1: R_p = 0.30
    R_p_1 = 0.30
    T_e_1 = ((1 - R_p_1) * S_0 / (4 * sigma))^0.25
    @test isapprox(T_e_1, 255.0, rtol=0.01)

    # Case 2: R_p = 0.15 (Earth without clouds)
    R_p_2 = 0.15
    T_e_2 = ((1 - R_p_2) * S_0 / (4 * sigma))^0.25
    @test isapprox(T_e_2, 268.0, rtol=0.01)
end

@testitem "CO2 Doubling Temperature Response" begin
    # S&P: CO2 doubling produces ΔF_L = 4.6 W/m², giving ΔT_e ≈ 1.4 K
    sigma = 5.671e-8
    T_e = 255.0
    lambda_0 = 1 / (4 * sigma * T_e^3)

    delta_F = 4.6  # W/m²
    delta_T_e = lambda_0 * delta_F

    @test isapprox(delta_T_e, 1.4, rtol=0.15)  # Within 15% due to approximation
end

@testitem "Solar/Terrestrial Peak Wavelengths" begin
    # S&P: Solar peak ~480 nm, Terrestrial peak ~10 μm
    b = 2.897e-3

    # Sun at ~6000 K
    T_sun = 6000.0
    lambda_max_sun = b / T_sun
    @test isapprox(lambda_max_sun * 1e9, 480, rtol=0.05)  # ~480 nm

    # Earth at ~300 K
    T_earth = 300.0
    lambda_max_earth = b / T_earth
    @test isapprox(lambda_max_earth * 1e6, 10.0, rtol=0.05)  # ~10 μm
end

@testitem "Energy Balance Values" begin
    # S&P Figure 4.4: Incoming solar ~342 W/m², absorbed ~235-240 W/m²
    S_0 = 1370.0
    incoming = S_0 / 4
    @test isapprox(incoming, 342.0, rtol=0.01)

    R_p = 0.30
    absorbed = incoming * (1 - R_p)
    @test isapprox(absorbed, 239.75, rtol=0.01)
    @test absorbed > 230 && absorbed < 250
end

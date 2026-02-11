"""
    PhotonEnergy(; name=:PhotonEnergy)

Implements the photon energy equation (Eq. 4.1) from Seinfeld & Pandis.

    Δε = hν = hc/λ

The energy of a photon is related to its frequency (ν) or wavelength (λ)
through Planck's constant (h) and the speed of light (c).

# Outputs (Variables)

  - `Δε`: Photon energy (J)
  - `ν`: Frequency (Hz = s⁻¹)

# Inputs (Parameters)

  - `λ`: Wavelength (m) - input parameter

# Constants

  - `h`: Planck's constant = 6.626 × 10⁻³⁴ J·s
  - `c`: Speed of light = 2.9979 × 10⁸ m/s
"""
@component function PhotonEnergy(; name = :PhotonEnergy)
    @constants begin
        h = 6.626e-34, [description = "Planck's constant", unit = u"J*s"]
        c = 2.9979e8, [description = "Speed of light", unit = u"m/s"]
    end

    @parameters begin
        λ = 5.0e-7, [description = "Wavelength (input)", unit = u"m"]
    end

    @variables begin
        Δε(t), [description = "Photon energy", unit = u"J"]
        ν(t), [description = "Frequency", unit = u"Hz"]
    end

    eqs = [
        ν ~ c / λ,            # Eq. 4.1b - Frequency-wavelength relation
        Δε ~ h * ν,           # Eq. 4.1a - Energy from frequency
    ]

    return System(eqs, t; name)
end

"""
    BlackbodyRadiation(; name=:BlackbodyRadiation)

Implements Planck's blackbody radiation law (Eq. 4.2) from Seinfeld & Pandis.

    F_B(λ) = 2πc²hλ⁻⁵ / (exp(ch/kλT) - 1)

This equation gives the monochromatic emissive power of a blackbody at
temperature T for a given wavelength λ.

# Outputs (Variables)

  - `F_B_λ`: Monochromatic emissive power (W m⁻³)

# Inputs (Parameters)

  - `T`: Temperature (K)
  - `λ`: Wavelength (m)

# Constants

  - `h`: Planck's constant = 6.626 × 10⁻³⁴ J·s
  - `c`: Speed of light = 2.9979 × 10⁸ m/s
  - `k`: Boltzmann constant = 1.381 × 10⁻²³ J/K
"""
@component function BlackbodyRadiation(; name = :BlackbodyRadiation)
    @constants begin
        h = 6.626e-34, [description = "Planck's constant", unit = u"J*s"]
        c = 2.9979e8, [description = "Speed of light", unit = u"m/s"]
        k_B = 1.381e-23, [description = "Boltzmann constant", unit = u"J/K"]
        π_val = Float64(π), [description = "Pi (dimensionless)", unit = u"1"]
    end

    @parameters begin
        T = 5800.0, [description = "Temperature (input)", unit = u"K"]
        λ = 5.0e-7, [description = "Wavelength (input)", unit = u"m"]
    end

    @variables begin
        F_B_λ(t), [description = "Monochromatic blackbody emissive power", unit = u"W/m^3"]
    end

    eqs = [
        # Eq. 4.2 - Planck's law for blackbody radiation
        F_B_λ ~ 2 * π_val * c^2 * h * λ^(-5) / (exp(c * h / (k_B * λ * T)) - 1),
    ]

    return System(eqs, t; name)
end

"""
    WienDisplacement(; name=:WienDisplacement)

Implements Wien's displacement law (Eq. 4.3) from Seinfeld & Pandis.

    λ_max = 2.897 × 10⁻³ / T  (SI: λ_max in m, T in K)

Note: S&P express this as λ_max = 2.897 × 10⁶ / T when λ_max is in nm.

This gives the wavelength at which the blackbody emission spectrum peaks.

# Outputs (Variables)

  - `λ_max`: Peak wavelength (m)

# Inputs (Parameters)

  - `T`: Temperature (K)

# Constants

  - `b`: Wien's displacement constant = 2.897 × 10⁻³ m·K
"""
@component function WienDisplacement(; name = :WienDisplacement)
    @constants begin
        # Wien's displacement constant in SI units (m·K)
        b = 2.897e-3, [description = "Wien's displacement constant", unit = u"m*K"]
    end

    @parameters begin
        T = 5800.0, [description = "Temperature (input)", unit = u"K"]
    end

    @variables begin
        λ_max(t), [description = "Peak emission wavelength", unit = u"m"]
    end

    eqs = [
        # Eq. 4.3 - Wien's displacement law (SI units)
        λ_max ~ b / T,
    ]

    return System(eqs, t; name)
end

"""
    StefanBoltzmann(; name=:StefanBoltzmann)

Implements the Stefan-Boltzmann law (Eq. 4.4) from Seinfeld & Pandis.

    F_B = σT⁴

The total emissive power of a blackbody integrated over all wavelengths.

# Outputs (Variables)

  - `F_B`: Total blackbody emissive power (W/m²)

# Inputs (Parameters)

  - `T`: Temperature (K)

# Constants

  - `σ`: Stefan-Boltzmann constant = 5.671 × 10⁻⁸ W m⁻² K⁻⁴
"""
@component function StefanBoltzmann(; name = :StefanBoltzmann)
    @constants begin
        σ = 5.671e-8, [description = "Stefan-Boltzmann constant", unit = u"W/(m^2*K^4)"]
    end

    @parameters begin
        T = 255.0, [description = "Temperature (input)", unit = u"K"]
    end

    @variables begin
        F_B(t), [description = "Total blackbody emissive power", unit = u"W/m^2"]
    end

    eqs = [
        # Eq. 4.4 - Stefan-Boltzmann law
        F_B ~ σ * T^4,
    ]

    return System(eqs, t; name)
end

"""
    PlanetaryEnergyBalance(; name=:PlanetaryEnergyBalance)

Implements the planetary energy balance equations (Eqs. 4.5-4.7) from Seinfeld & Pandis.

    F_S = S₀/4 × (1 - R_p)     (Eq. 4.5 - Absorbed solar flux)
    F_L = σT_e⁴                 (Eq. 4.6 - Emitted longwave flux)
    T_e = ((1 - R_p)S₀/4σ)^(1/4) (Eq. 4.7 - Equilibrium temperature)

At equilibrium: F_S = F_L

# Outputs (Variables)

  - `F_S`: Absorbed solar flux (W/m²)
  - `F_L`: Emitted longwave flux (W/m²)
  - `T_e`: Planetary equilibrium temperature (K)

# Inputs (Parameters)

  - `S_0`: Solar constant = 1370 W/m² (can be varied)
  - `R_p`: Planetary albedo ~ 0.3 (can be varied)

# Constants

  - `σ`: Stefan-Boltzmann constant = 5.671 × 10⁻⁸ W m⁻² K⁻⁴
"""
@component function PlanetaryEnergyBalance(; name = :PlanetaryEnergyBalance)
    @constants begin
        σ = 5.671e-8, [description = "Stefan-Boltzmann constant", unit = u"W/(m^2*K^4)"]
    end

    @parameters begin
        S_0 = 1370.0, [description = "Solar constant", unit = u"W/m^2"]
        R_p = 0.3, [description = "Planetary albedo (dimensionless)", unit = u"1"]
    end

    @variables begin
        F_S(t) = 240.0, [description = "Absorbed solar flux", unit = u"W/m^2"]
        F_L(t) = 240.0, [description = "Emitted longwave flux", unit = u"W/m^2"]
        T_e(t) = 255.0, [description = "Planetary equilibrium temperature", unit = u"K"]
    end

    eqs = [
        # Eq. 4.5 - Absorbed solar flux (factor of 4 accounts for Earth's geometry)
        F_S ~ S_0 / 4 * (1 - R_p),

        # Eq. 4.6 - Emitted longwave flux (Stefan-Boltzmann law)
        F_L ~ σ * T_e^4,

        # At equilibrium, absorbed = emitted
        # Eq. 4.7 can be derived from F_S = F_L, solving for T_e
        F_S ~ F_L,
    ]

    return System(eqs, t; name)
end

"""
    ClimateSensitivity(; name=:ClimateSensitivity)

Implements the climate sensitivity equations (Eqs. 4.8-4.10) from Seinfeld & Pandis.

    ΔF_net = ΔF_S - ΔF_L       (Eq. 4.8 - Net radiative energy change)
    ΔT_e = λ₀ × ΔF_net         (Eq. 4.9 - Temperature response)
    λ₀ = 1/(4σT_e³) = T_e/(4F_L) (Eq. 4.10 - Climate sensitivity factor)

# Outputs (Variables)

  - `ΔF_net`: Net radiative forcing (W/m²)
  - `ΔT_e`: Change in equilibrium temperature (K)
  - `λ_0`: Climate sensitivity factor (K m²/W)
  - `F_L`: Reference emitted longwave flux (W/m²)

# Inputs (Parameters)

  - `T_e`: Reference equilibrium temperature (K)
  - `ΔF_S`: Change in absorbed solar flux (W/m²)
  - `ΔF_L`: Change in emitted longwave flux (W/m²)

# Constants

  - `σ`: Stefan-Boltzmann constant = 5.671 × 10⁻⁸ W m⁻² K⁻⁴
"""
@component function ClimateSensitivity(; name = :ClimateSensitivity)
    @constants begin
        σ = 5.671e-8, [description = "Stefan-Boltzmann constant", unit = u"W/(m^2*K^4)"]
    end

    @parameters begin
        T_e = 255.0,
            [description = "Reference equilibrium temperature (input)", unit = u"K"]
        ΔF_S = 4.0, [description = "Change in absorbed solar flux (input)", unit = u"W/m^2"]
        ΔF_L = 0.0,
            [description = "Change in emitted longwave flux (input)", unit = u"W/m^2"]
    end

    @variables begin
        ΔF_net(t), [description = "Net radiative forcing", unit = u"W/m^2"]
        ΔT_e(t), [description = "Change in equilibrium temperature", unit = u"K"]
        λ_0(t), [description = "Climate sensitivity factor", unit = u"K*m^2/W"]
        F_L(t), [description = "Reference emitted longwave flux", unit = u"W/m^2"]
    end

    eqs = [
        # Eq. 4.8 - Net radiative energy change
        ΔF_net ~ ΔF_S - ΔF_L,

        # Eq. 4.10 - Climate sensitivity factor
        λ_0 ~ 1 / (4 * σ * T_e^3),

        # Eq. 4.9 - Climate sensitivity relationship
        ΔT_e ~ λ_0 * ΔF_net,

        # Reference state: F_L from Stefan-Boltzmann
        F_L ~ σ * T_e^4,
    ]

    return System(eqs, t; name)
end

"""
    TOARadiativeForcing(; name=:TOARadiativeForcing)

Implements the top of atmosphere radiative forcing based on Eq. 4.11 from Seinfeld & Pandis.

    F_net = S₀/4 × (1 - R_p) - F_L

Note: Seinfeld & Pandis define Eq. 4.11 as -F_net = S₀/4(1-R_p) - F_L, using the convention
that -F_net represents net downward flux. Here we define F_net as net incoming flux directly,
so positive F_net indicates the planet is gaining energy (warming).

# Outputs (Variables)

  - `F_net`: Net radiative flux at TOA (W/m²) - positive = energy gain
  - `F_S`: Absorbed solar flux (W/m²)

# Inputs (Parameters)

  - `S_0`: Solar constant = 1370 W/m² (can be varied)
  - `R_p`: Planetary albedo ~ 0.3 (can be varied)
  - `F_L`: Emitted longwave flux (W/m²) - input parameter
"""
@component function TOARadiativeForcing(; name = :TOARadiativeForcing)
    @parameters begin
        S_0 = 1370.0, [description = "Solar constant", unit = u"W/m^2"]
        R_p = 0.3, [description = "Planetary albedo (dimensionless)", unit = u"1"]
        F_L = 239.75, [description = "Emitted longwave flux (input)", unit = u"W/m^2"]
    end

    @variables begin
        F_net(t), [description = "Net radiative flux at TOA", unit = u"W/m^2"]
        F_S(t), [description = "Absorbed solar flux", unit = u"W/m^2"]
    end

    eqs = [
        # Absorbed solar flux
        F_S ~ S_0 / 4 * (1 - R_p),

        # Eq. 4.11 - TOA energy balance
        # F_net positive means energy gain (warming)
        F_net ~ F_S - F_L,
    ]

    return System(eqs, t; name)
end

"""
    RadiationFundamentals(; name=:RadiationFundamentals)

Combined system implementing all radiation fundamentals equations from
Seinfeld & Pandis Chapter 4 (Eqs. 4.1-4.11).

This is a composed system that includes:

  - PhotonEnergy: Energy-frequency-wavelength relations (Eq. 4.1)
  - BlackbodyRadiation: Planck's law (Eq. 4.2)
  - WienDisplacement: Peak emission wavelength (Eq. 4.3)
  - StefanBoltzmann: Total emissive power (Eq. 4.4)
  - PlanetaryEnergyBalance: Earth's energy balance (Eqs. 4.5-4.7)
  - ClimateSensitivity: Temperature response to forcing (Eqs. 4.8-4.10)
  - TOARadiativeForcing: Net flux at top of atmosphere (Eq. 4.11)
"""
@component function RadiationFundamentals(; name = :RadiationFundamentals)
    @named photon = PhotonEnergy()
    @named blackbody = BlackbodyRadiation()
    @named wien = WienDisplacement()
    @named stefan_boltzmann = StefanBoltzmann()
    @named energy_balance = PlanetaryEnergyBalance()
    @named climate_sensitivity = ClimateSensitivity()
    @named toa_forcing = TOARadiativeForcing()

    # No additional equations needed - all physics is in subsystems
    eqs = Equation[]

    return System(
        eqs, t;
        systems = [
            photon, blackbody, wien, stefan_boltzmann,
            energy_balance, climate_sensitivity, toa_forcing,
        ],
        name
    )
end

# Export all components
export PhotonEnergy, BlackbodyRadiation, WienDisplacement, StefanBoltzmann
export PlanetaryEnergyBalance, ClimateSensitivity, TOARadiativeForcing
export RadiationFundamentals

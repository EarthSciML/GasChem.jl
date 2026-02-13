export ClimateFeedback, GlobalWarmingPotential, GWP_exponential, GHGForcing

"""
$(TYPEDSIGNATURES)

Climate feedback model relating radiative forcing to equilibrium temperature change with feedbacks.

This model implements the climate sensitivity equations with feedback from Seinfeld & Pandis (2006),
Chapter 23 "Climate and Chemical Composition of the Atmosphere", Section 23.4.

Note: This extends the basic no-feedback climate sensitivity from Chapter 4 (see
[`ClimateSensitivity`](@ref)) by incorporating climate feedback processes.

The model includes:

  - Eq. 23.1: ΔTs = λ ΔF (climate sensitivity with feedbacks)
  - Eq. 23.2: ΔT₀ = λ₀ ΔF (no-feedback temperature response)
  - Eq. 23.3: Eᵢ = λᵢ / λ_CO₂ (efficacy of forcing agents)
  - Eq. 23.4: ΔFₑ = ΔFᵢ Eᵢ (effective forcing)
  - Eq. 23.7: ΔT_unrealized = (ΔF - ΔF_r) λ (unrealized warming, p. 1045)

**Reference**: Seinfeld, J.H. and Pandis, S.N. (2006) *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change*, 2nd Edition, Chapter 23, John Wiley & Sons.

# Example

```julia
using GasChem, ModelingToolkit
sys = ClimateFeedback()
```
"""
@component function ClimateFeedback(; name = :ClimateFeedback)
    @constants begin
        # Radiative forcing for CO₂ doubling (IPCC 2001 best estimate, p. 1037)
        ΔF_2xCO2 = 3.7,
            [description = "Radiative forcing for CO₂ doubling", unit = u"W/m^2"]
        # No-feedback temperature increase for 2×CO₂ (p. 1037-1038)
        ΔT0_2xCO2 = 1.25,
            [description = "No-feedback temperature change for CO₂ doubling", unit = u"K"]
    end

    @parameters begin
        ΔF, [description = "Radiative forcing perturbation", unit = u"W/m^2"]
        λ = 0.8,
            [description = "Climate sensitivity parameter with feedbacks", unit = u"K*m^2/W"]
        λ_CO2 = 0.8,
            [description = "Climate sensitivity for CO₂ doubling", unit = u"K*m^2/W"]
        E_i = 1.0, [description = "Efficacy of forcing agent (dimensionless)", unit = u"1"]
        ΔT_r = 0.7, [description = "Realized warming to date", unit = u"K"]
    end

    @variables begin
        ΔT_s(t),
            [description = "Equilibrium surface temperature change with feedbacks", unit = u"K"]
        ΔT_0(t), [description = "Surface temperature change without feedbacks", unit = u"K"]
        λ_0(t), [description = "No-feedback climate sensitivity", unit = u"K*m^2/W"]
        ΔF_e(t), [description = "Effective radiative forcing", unit = u"W/m^2"]
        ΔT_unrealized(t), [description = "Unrealized (committed) warming", unit = u"K"]
        ΔF_r(t),
            [description = "Forcing corresponding to realized warming", unit = u"W/m^2"]
    end

    eqs = [
        # Eq. 23.1 - Climate sensitivity relation with feedbacks
        ΔT_s ~ λ * ΔF,

        # No-feedback climate sensitivity (derived from ΔT₀ = 1.2-1.3 K for ΔF ≈ 3.7 W/m², p. 1040)
        λ_0 ~ ΔT0_2xCO2 / ΔF_2xCO2,

        # Eq. 23.2 - No-feedback temperature response
        ΔT_0 ~ λ_0 * ΔF,

        # Eq. 23.4 - Effective forcing (incorporating efficacy, Eq. 23.3)
        ΔF_e ~ ΔF * E_i,

        # Forcing for realized warming (inverse of Eq. 23.1)
        ΔF_r ~ ΔT_r / λ,

        # Eq. 23.7 - Unrealized warming (warming commitment, p. 1045)
        ΔT_unrealized ~ (ΔF - ΔF_r) * λ,
    ]

    return System(eqs, t; name)
end

"""
$(TYPEDSIGNATURES)

Greenhouse gas radiative forcing model.

Computes the radiative forcing from well-mixed greenhouse gases using values from
IPCC (2001) as presented in Seinfeld & Pandis (2006), Chapter 23, Section 23.3.

**Reference**: Seinfeld, J.H. and Pandis, S.N. (2006) *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change*, 2nd Edition, Chapter 23, John Wiley & Sons.

# Example

```julia
using GasChem, ModelingToolkit
sys = GHGForcing()
```
"""
@component function GHGForcing(; name = :GHGForcing)
    @constants begin
        # Individual forcings (IPCC 2001, presented in S&P p. 1039)
        ΔF_CO2_ref = 1.46,
            [description = "CO₂ forcing since preindustrial", unit = u"W/m^2"]
        ΔF_CH4_ref = 0.48,
            [description = "CH₄ forcing since preindustrial", unit = u"W/m^2"]
        ΔF_N2O_ref = 0.15,
            [description = "N₂O forcing since preindustrial", unit = u"W/m^2"]
        ΔF_O3_trop_ref = 0.4, [description = "Tropospheric O₃ forcing", unit = u"W/m^2"]
        ΔF_halo_ref = 0.34, [description = "Halocarbon forcing", unit = u"W/m^2"]
    end

    @parameters begin
        f_CO2 = 1.0,
            [description = "Scaling factor for CO₂ forcing (dimensionless)", unit = u"1"]
        f_CH4 = 1.0,
            [description = "Scaling factor for CH₄ forcing (dimensionless)", unit = u"1"]
        f_N2O = 1.0,
            [description = "Scaling factor for N₂O forcing (dimensionless)", unit = u"1"]
        f_O3 = 1.0,
            [
                description = "Scaling factor for tropospheric O₃ forcing (dimensionless)",
                unit = u"1",
            ]
        f_halo = 1.0,
            [description = "Scaling factor for halocarbon forcing (dimensionless)", unit = u"1"]
    end

    @variables begin
        ΔF_CO2(t), [description = "Radiative forcing from CO₂", unit = u"W/m^2"]
        ΔF_CH4(t), [description = "Radiative forcing from CH₄", unit = u"W/m^2"]
        ΔF_N2O(t), [description = "Radiative forcing from N₂O", unit = u"W/m^2"]
        ΔF_O3(t), [description = "Radiative forcing from tropospheric O₃", unit = u"W/m^2"]
        ΔF_halo(t), [description = "Radiative forcing from halocarbons", unit = u"W/m^2"]
        ΔF_total(t), [description = "Total radiative forcing from GHGs", unit = u"W/m^2"]
    end

    eqs = [
        ΔF_CO2 ~ f_CO2 * ΔF_CO2_ref,
        ΔF_CH4 ~ f_CH4 * ΔF_CH4_ref,
        ΔF_N2O ~ f_N2O * ΔF_N2O_ref,
        ΔF_O3 ~ f_O3 * ΔF_O3_trop_ref,
        ΔF_halo ~ f_halo * ΔF_halo_ref,
        ΔF_total ~ ΔF_CO2 + ΔF_CH4 + ΔF_N2O + ΔF_O3 + ΔF_halo,
    ]

    return System(eqs, t; name)
end

"""
$(TYPEDSIGNATURES)

Global Warming Potential (GWP) calculation for exponentially decaying species.

Implements Eq. 23.5 and 23.6 from Seinfeld & Pandis (2006), Chapter 23, for the
special case of exponential decay of atmospheric perturbations. Also implements
the closed-form GWP formula from Problem 23.3 (p. 1050).

The GWP is defined as the ratio of the time-integrated radiative forcing from a
pulse emission of 1 kg of compound A to that of 1 kg of CO₂ over a specified time horizon.

**Equations implemented:**

  - Eq. 23.5: GWP = ∫₀^{tf} aₐ[A(t)]dt / ∫₀^{tf} a_R[R(t)]dt
  - Eq. 23.6: AGWP_A = ∫₀^{tf} aₐ[A(t)]dt
  - Problem 23.3: GWP = [aₐ τₐ (1 - e^{-tf/τₐ})] / [a_CO₂ τ_CO₂ (1 - e^{-tf/τ_CO₂})]

**Reference**: Seinfeld, J.H. and Pandis, S.N. (2006) *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change*, 2nd Edition, Chapter 23, John Wiley & Sons.

# Example

```julia
using GasChem, ModelingToolkit
sys = GlobalWarmingPotential()
```
"""
@component function GlobalWarmingPotential(; name = :GlobalWarmingPotential)
    @constants begin
        # Effective CO₂ lifetime for GWP calculations (assumed, see p. 1043)
        # 150 years converted to seconds
        τ_CO2 = 150.0 * 365.25 * 24 * 3600,
            [description = "Effective CO₂ lifetime for GWP calculations (150 yr)", unit = u"s"]
        # Radiative efficiency of CO₂ (approximate)
        a_CO2 = 1.0,
            [description = "Reference radiative efficiency for CO₂ (normalized)", unit = u"1"]
    end

    @parameters begin
        τ_A, [description = "Atmospheric lifetime of species A", unit = u"s"]
        a_A,
            [
                description = "Radiative efficiency of species A relative to CO₂ (dimensionless)",
                unit = u"1",
            ]
        t_f = 100.0 * 365.25 * 24 * 3600,
            [description = "Time horizon for GWP calculation (default 100 yr)", unit = u"s"]
    end

    @variables begin
        GWP(t),
            [
                description = "Global Warming Potential relative to CO₂ (dimensionless)",
                unit = u"1",
            ]
        AGWP_A(t),
            [description = "Absolute GWP of species A", unit = u"s"]
        AGWP_CO2(t), [description = "Absolute GWP of CO₂", unit = u"s"]
        decay_integral_A(t),
            [description = "Time integral of decay function for species A", unit = u"s"]
        decay_integral_CO2(t),
            [description = "Time integral of decay function for CO₂", unit = u"s"]
    end

    eqs = [
        # Time integral of exponential decay: τ(1 - exp(-tf/τ))
        # This comes from ∫₀^{tf} exp(-t/τ) dt = τ(1 - exp(-tf/τ))
        decay_integral_A ~ τ_A * (1 - exp(-t_f / τ_A)),
        decay_integral_CO2 ~ τ_CO2 * (1 - exp(-t_f / τ_CO2)),

        # Eq. 23.6 - Absolute GWP (proportional to a × decay_integral)
        AGWP_A ~ a_A * decay_integral_A,
        AGWP_CO2 ~ a_CO2 * decay_integral_CO2,

        # Eq. 23.5 / Problem 23.3 - GWP as ratio of AGWPs
        GWP ~ AGWP_A / AGWP_CO2,
    ]

    return System(eqs, t; name)
end

"""
    GWP_exponential(τ_A, a_A, t_f; τ_CO2=150.0)

Calculate the Global Warming Potential using the closed-form solution for exponential decay.

This implements the formula from Problem 23.3 (p. 1050) of Seinfeld & Pandis (2006):

```math
GWP = \\frac{a_A \\tau_A (1 - e^{-t_f/\\tau_A})}{a_{CO_2} \\tau_{CO_2} (1 - e^{-t_f/\\tau_{CO_2}})}
```

# Arguments

  - `τ_A`: Atmospheric lifetime of species A [years]
  - `a_A`: Radiative efficiency of species A relative to CO₂ (dimensionless)
  - `t_f`: Time horizon for GWP calculation [years]
  - `τ_CO2`: Effective CO₂ lifetime (default: 150 years)

# Returns

  - `GWP`: Global Warming Potential relative to CO₂ (dimensionless)

# Example

```julia
# Calculate 100-year GWP for CH₄ (τ = 12 yr, fitted relative efficiency ≈ 140)
gwp_ch4 = GWP_exponential(12.0, 140.0, 100.0)
```

**Reference**: Seinfeld, J.H. and Pandis, S.N. (2006) *Atmospheric Chemistry and Physics:
From Air Pollution to Climate Change*, 2nd Edition, Chapter 23, Problem 23.3.
"""
function GWP_exponential(τ_A, a_A, t_f; τ_CO2 = 150.0)
    # Decay integrals
    decay_A = τ_A * (1 - exp(-t_f / τ_A))
    decay_CO2 = τ_CO2 * (1 - exp(-t_f / τ_CO2))

    # GWP calculation
    return a_A * decay_A / decay_CO2
end

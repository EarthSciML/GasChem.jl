"""
    CO Oxidation Chemistry

Implements equations 6.9-6.24 from Seinfeld & Pandis Chapter 6, Section 6.3.

CO oxidation is the simplest hydrocarbon oxidation mechanism and illustrates
the fundamental HOx cycling that drives tropospheric ozone production:

    CO + OH → CO₂ + H         (R6.9)
    H + O₂ + M → HO₂ + M      (R6.10)
    HO₂ + NO → OH + NO₂       (R6.11)
    NO₂ + hν → NO + O         (R6.12)
    O + O₂ + M → O₃ + M       (R6.13)

Net: CO + 2O₂ + hν → CO₂ + O₃

The key feature is that OH is regenerated through the HO₂ + NO reaction,
allowing catalytic ozone production.

Reference: Seinfeld & Pandis (2006), Section 6.3, pp. 212-219
"""

using ModelingToolkit: t, D

# ============================================================================
# Equations 6.9-6.13: Basic CO Oxidation Cycle
# ============================================================================
# These reactions form the basic cycle for CO oxidation with O₃ production.

# ============================================================================
# Equation 6.18: HO₂ Steady-State (High NOx)
# ============================================================================
# [HO₂] = k_{CO+OH}[CO][OH] / (k_{HO₂+NO}[NO])
#
# At high NOx, HO₂ is removed primarily by reaction with NO.

# ============================================================================
# Equation 6.15: HO₂ Self-Reaction (Low NOx)
# ============================================================================
# HO₂ + HO₂ → H₂O₂ + O₂
#
# At low NOx, HO₂ self-reaction becomes the dominant HOx loss.

# ============================================================================
# Equation 6.17: OH Steady-State
# ============================================================================
# [OH] = P_OH / (k₉[CO] + k_{OH+NO₂}[NO₂] + ...)
#
# OH concentration is determined by balance of production (mainly O₃ photolysis)
# and loss to various reactants.

# ============================================================================
# Equations 6.21-6.24: Ozone Production Efficiency
# ============================================================================
# OPE = Δ[O₃] / Δ[NOx] = moles O₃ produced per mole NOx consumed
#
# This is a key metric for understanding the efficiency of ozone production
# in different NOx regimes.

"""
    COOxidation(; name)

ModelingToolkit System for CO oxidation chemistry diagnostics.

Implements diagnostic calculations for the CO oxidation mechanism (Eqs. 6.9-6.17)
from Seinfeld & Pandis Chapter 6, including HOx cycling and steady-state relationships.

# Input Variables

  - `CO`, `OH`, `HO2`, `NO`, `NO2`, `O3`: Species concentrations [m⁻³]

# Output Variables

  - `P_O3`: Net ozone production rate [m⁻³ s⁻¹]
  - `L_HOx`: Total HOx loss rate [m⁻³ s⁻¹]
  - `chain_length`: HOx chain length (cycles before termination)
  - `HO2_ss`: HO₂ steady-state concentration (high NOx limit) [m⁻³]

# Rate Constants at 298 K

  - k_CO_OH = 2.4 × 10⁻¹³ cm³ molecule⁻¹ s⁻¹ = 2.4 × 10⁻¹⁹ m³ s⁻¹
  - k_HO2_NO = 8.1 × 10⁻¹² cm³ molecule⁻¹ s⁻¹ = 8.1 × 10⁻¹⁸ m³ s⁻¹
  - k_HO2_HO2 = 2.9 × 10⁻¹² cm³ molecule⁻¹ s⁻¹ = 2.9 × 10⁻¹⁸ m³ s⁻¹
  - k_OH_NO2 = 1.0 × 10⁻¹¹ cm³ molecule⁻¹ s⁻¹ = 1.0 × 10⁻¹⁷ m³ s⁻¹
  - k_HO2_O3 = 2.0 × 10⁻¹⁵ cm³ molecule⁻¹ s⁻¹ = 2.0 × 10⁻²¹ m³ s⁻¹
  - k_OH_O3 = 7.3 × 10⁻¹⁴ cm³ molecule⁻¹ s⁻¹ = 7.3 × 10⁻²⁰ m³ s⁻¹
"""
@component function COOxidation(; name = :COOxidation)
    @constants begin
        two = 2,
            [
                description = "Stoichiometric factor for HO₂ self-reaction (dimensionless)", unit = u"1",
            ]
    end

    # Parameters - Rate constants at 298 K (converted to SI)
    @parameters begin
        k_CO_OH = 2.4e-13 * 1.0e-6,
            [description = "CO + OH rate constant (2.4e-13 cm³/molec/s)", unit = u"m^3/s"]
        k_HO2_NO = 8.1e-12 * 1.0e-6,
            [description = "HO₂ + NO rate constant (8.1e-12 cm³/molec/s)", unit = u"m^3/s"]
        k_HO2_HO2 = 2.9e-12 * 1.0e-6,
            [description = "HO₂ + HO₂ rate constant (2.9e-12 cm³/molec/s)", unit = u"m^3/s"]
        k_OH_NO2 = 1.0e-11 * 1.0e-6,
            [description = "OH + NO₂ + M rate constant (1.0e-11 cm³/molec/s)", unit = u"m^3/s"]
        k_HO2_O3 = 2.0e-15 * 1.0e-6,
            [description = "HO₂ + O₃ rate constant (2.0e-15 cm³/molec/s)", unit = u"m^3/s"]
        k_OH_O3 = 7.3e-14 * 1.0e-6,
            [description = "OH + O₃ rate constant (7.3e-14 cm³/molec/s)", unit = u"m^3/s"]
    end

    # Input variables (concentrations in SI: m⁻³)
    @variables begin
        CO(t), [description = "CO concentration", unit = u"m^-3"]
        OH(t), [description = "OH concentration", unit = u"m^-3"]
        HO2(t), [description = "HO₂ concentration", unit = u"m^-3"]
        NO(t), [description = "NO concentration", unit = u"m^-3"]
        NO2(t), [description = "NO₂ concentration", unit = u"m^-3"]
        O3(t), [description = "O₃ concentration", unit = u"m^-3"]
    end

    # Output variables
    @variables begin
        P_O3(t), [description = "Net O₃ production rate", unit = u"m^-3*s^-1"]
        L_HOx(t), [description = "HOx loss rate", unit = u"m^-3*s^-1"]
        L_OH(t), [description = "OH loss rate", unit = u"m^-3*s^-1"]
        L_HO2(t), [description = "HO₂ loss rate", unit = u"m^-3*s^-1"]
        chain_length(t), [description = "HOx chain length (dimensionless)", unit = u"1"]
        HO2_ss(t), [description = "HO₂ steady-state (high NOx)", unit = u"m^-3"]
    end

    # Equations
    eqs = [
        # OH loss rate (to CO, NO₂, O₃)
        L_OH ~ k_CO_OH * CO * OH + k_OH_NO2 * OH * NO2 + k_OH_O3 * OH * O3,

        # HO₂ loss rate (to NO, HO₂, O₃)
        L_HO2 ~ k_HO2_NO * HO2 * NO + two * k_HO2_HO2 * HO2^2 + k_HO2_O3 * HO2 * O3,

        # HOx loss rate (radical termination reactions)
        # Main termination: OH + NO₂ → HNO₃ and HO₂ + HO₂ → H₂O₂
        L_HOx ~ k_OH_NO2 * OH * NO2 + two * k_HO2_HO2 * HO2^2,

        # Equation 6.18: HO₂ steady-state (high NOx limit)
        HO2_ss ~ k_CO_OH * CO * OH / (k_HO2_NO * NO),

        # Net ozone production rate (Eq. 6.9 production minus O₃ loss terms)
        # O₃ production from HO₂ + NO minus loss from OH + O₃ and HO₂ + O₃
        P_O3 ~ k_HO2_NO * HO2 * NO - k_OH_O3 * OH * O3 - k_HO2_O3 * HO2 * O3,

        # HOx chain length: number of OH-HO₂ cycles before termination
        chain_length ~ (k_HO2_NO * HO2 * NO) / L_HOx,
    ]

    return System(eqs, t; name)
end

"""
    OzoneProductionEfficiency(; name)

ModelingToolkit System for calculating Ozone Production Efficiency (OPE).

OPE = Δ[O₃] / Δ[NOx] represents the number of ozone molecules produced
per NOx molecule consumed. This is a key metric for understanding
ozone-NOx-VOC chemistry in different regimes.

From Equations 6.21-6.24:

  - At high NOx (VOC-limited): OPE is low
  - At low NOx (NOx-limited): OPE is high
  - The transition occurs near the "ridge line" in O₃-NOx-VOC space

# Input Variables

  - `OH`, `HO2`, `RO2`, `NO`, `NO2`: Species concentrations [m⁻³]

# Output Variables

  - `P_O3`: Gross ozone production rate [m⁻³ s⁻¹]
  - `L_NOx`: NOx loss rate [m⁻³ s⁻¹]
  - `OPE`: Ozone production efficiency [dimensionless]
"""
@component function OzoneProductionEfficiency(; name = :OzoneProductionEfficiency)
    @parameters begin
        k_HO2_NO = 8.1e-12 * 1.0e-6,
            [description = "HO₂ + NO rate constant (8.1e-12 cm³/molec/s)", unit = u"m^3/s"]
        k_RO2_NO = 8.0e-12 * 1.0e-6,
            [description = "RO₂ + NO rate constant (8.0e-12 cm³/molec/s)", unit = u"m^3/s"]
        k_OH_NO2 = 1.0e-11 * 1.0e-6,
            [description = "OH + NO₂ + M rate constant (1.0e-11 cm³/molec/s)", unit = u"m^3/s"]
    end

    # Input variables
    @variables begin
        OH(t), [description = "OH concentration", unit = u"m^-3"]
        HO2(t), [description = "HO₂ concentration", unit = u"m^-3"]
        RO2(t), [description = "Organic peroxy radical concentration", unit = u"m^-3"]
        NO(t), [description = "NO concentration", unit = u"m^-3"]
        NO2(t), [description = "NO₂ concentration", unit = u"m^-3"]
    end

    # Output variables
    @variables begin
        P_O3(t), [description = "Gross O₃ production rate", unit = u"m^-3*s^-1"]
        L_NOx(t), [description = "NOx loss rate", unit = u"m^-3*s^-1"]
        OPE(t), [description = "Ozone Production Efficiency (dimensionless)", unit = u"1"]
    end

    eqs = [
        # Gross O₃ production from peroxy radical + NO reactions
        P_O3 ~ k_HO2_NO * HO2 * NO + k_RO2_NO * RO2 * NO,

        # NOx loss (mainly OH + NO₂ → HNO₃)
        L_NOx ~ k_OH_NO2 * OH * NO2,

        # OPE = P(O₃) / L(NOx)
        OPE ~ P_O3 / L_NOx,
    ]

    return System(eqs, t; name)
end

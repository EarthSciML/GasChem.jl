"""
    OH Radical Production from O₃ Photolysis

Implements equations 6.1-6.4 from Seinfeld & Pandis Chapter 6, Section 6.1.

The primary source of OH radicals in the troposphere is the photolysis of O₃
at wavelengths below 320 nm, producing electronically excited oxygen atoms O(¹D):

    O₃ + hν → O(¹D) + O₂     (λ < 320 nm)

O(¹D) can either be quenched back to ground state O(³P) by collision with N₂ or O₂,
or react with H₂O to produce OH:

    O(¹D) + M → O(³P) + M    (M = N₂, O₂)
    O(¹D) + H₂O → 2OH

Reference: Seinfeld & Pandis (2006), Section 6.1, pp. 204-207
"""

using ModelingToolkit: t, D

# ============================================================================
# Equation 6.1: O(¹D) Steady-State Concentration
# ============================================================================
# [O(¹D)] = j_{O₃→O(¹D)} [O₃] / (k₃[M] + k₄[H₂O])
#
# where:
#   j_{O₃→O(¹D)} = photolysis rate of O₃ producing O(¹D) [s⁻¹]
#   k₃ = quenching rate constant for O(¹D) + M [cm³ molecule⁻¹ s⁻¹]
#   k₄ = rate constant for O(¹D) + H₂O → 2OH [cm³ molecule⁻¹ s⁻¹]
#   [M] = total air concentration (N₂ + O₂) [molecules cm⁻³]
#   [H₂O] = water vapor concentration [molecules cm⁻³]

# ============================================================================
# Equation 6.3: OH Production Rate
# ============================================================================
# P_OH = 2 * j_{O₃→O(¹D)} * k₄[H₂O] / (k₃[M]) * [O₃]
#
# This is the approximate form valid when k₃[M] >> k₄[H₂O], which is typically
# the case in the troposphere.

# ============================================================================
# Equation 6.4: OH Yield (Fraction of O(¹D) producing OH)
# ============================================================================
# ε_OH = k₄[H₂O] / (k₃[M] + k₄[H₂O])
#
# This represents the fraction of O(¹D) atoms that react with H₂O to produce OH
# rather than being quenched back to O(³P).

"""
    OHProduction(; name)

ModelingToolkit System for OH radical production from O₃ photolysis.

Implements the steady-state O(¹D) concentration (Eq. 6.1), OH production rate (Eq. 6.3),
and OH yield (Eq. 6.4) from Seinfeld & Pandis Chapter 6.

This is an algebraic system that computes diagnostic quantities from input concentrations.

# Input Variables (must be provided)
- `O3`: Ozone concentration [molecules cm⁻³]
- `H2O`: Water vapor concentration [molecules cm⁻³]
- `M`: Total air number density [molecules cm⁻³]

# Output Variables (computed)
- `O1D`: Excited oxygen O(¹D) steady-state concentration [molecules cm⁻³]
- `P_OH`: OH production rate [molecules cm⁻³ s⁻¹]
- `ε_OH`: OH yield (fraction of O(¹D) producing OH) [dimensionless]

# Parameters
- `j_O3`: O₃ photolysis rate producing O(¹D) [s⁻¹]
- `k3_N2`: Rate constant for O(¹D) + N₂ quenching [cm³ molecule⁻¹ s⁻¹]
- `k3_O2`: Rate constant for O(¹D) + O₂ quenching [cm³ molecule⁻¹ s⁻¹]
- `k4`: Rate constant for O(¹D) + H₂O → 2OH [cm³ molecule⁻¹ s⁻¹]
- `f_N2`: Fraction of M that is N₂ [dimensionless]
- `f_O2`: Fraction of M that is O₂ [dimensionless]

# Rate Constants at 298 K (from Table B.1, Seinfeld & Pandis)
- k3_N2 = 2.6 × 10⁻¹¹ cm³ molecule⁻¹ s⁻¹
- k3_O2 = 4.0 × 10⁻¹¹ cm³ molecule⁻¹ s⁻¹
- k4 = 2.2 × 10⁻¹⁰ cm³ molecule⁻¹ s⁻¹
"""
@component function OHProduction(; name=:OHProduction)
    # Parameters
    @parameters begin
        j_O3 = 1e-5, [description = "O₃ photolysis rate producing O(¹D)", unit = u"s^-1"]
        k3_N2 = 2.6e-11, [description = "O(¹D) + N₂ quenching rate", unit = u"cm^3/molec/s"]
        k3_O2 = 4.0e-11, [description = "O(¹D) + O₂ quenching rate", unit = u"cm^3/molec/s"]
        k4 = 2.2e-10, [description = "O(¹D) + H₂O → 2OH rate", unit = u"cm^3/molec/s"]
        f_N2 = 0.78, [description = "Fraction of air that is N₂ (dimensionless)", unit = u"1"]
        f_O2 = 0.21, [description = "Fraction of air that is O₂ (dimensionless)", unit = u"1"]
    end

    # Input variables (concentrations)
    @variables begin
        O3(t), [description = "Ozone concentration", unit = u"molec/cm^3"]
        H2O(t), [description = "Water vapor concentration", unit = u"molec/cm^3"]
        M(t), [description = "Total air number density", unit = u"molec/cm^3"]
    end

    # Computed intermediate: effective quenching rate constant.
    # Note: k3_eff depends on symbolic parameters f_N2, k3_N2, f_O2, k3_O2,
    # so it cannot be placed in a @constants block. It remains a symbolic expression.
    k3_eff = f_N2 * k3_N2 + f_O2 * k3_O2

    # Output variables (diagnostics computed from steady-state relations)
    @variables begin
        O1D(t), [description = "O(¹D) steady-state concentration", unit = u"molec/cm^3"]
        P_OH(t), [description = "OH production rate", unit = u"molec/cm^3/s"]
        ε_OH(t), [description = "OH yield fraction (dimensionless)", unit = u"1"]
    end

    # Equations (algebraic relationships)
    eqs = [
        # Equation 6.1: O(¹D) steady-state concentration
        O1D ~ j_O3 * O3 / (k3_eff * M + k4 * H2O),

        # Equation 6.4: OH yield (fraction of O(¹D) producing OH)
        ε_OH ~ k4 * H2O / (k3_eff * M + k4 * H2O),

        # Equation 6.3: OH production rate (2 OH per O(¹D) + H₂O reaction)
        P_OH ~ 2 * j_O3 * O3 * ε_OH,
    ]

    return System(eqs, t; name, checks=false)
end

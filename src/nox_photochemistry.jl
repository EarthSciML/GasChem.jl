"""
    NOx Photochemical Cycle

Implements equations 6.5-6.8 from Seinfeld & Pandis Chapter 6, Section 6.2.

The NOx (NO + NO₂) photochemical cycle is fundamental to tropospheric ozone chemistry:

    NO₂ + hν → NO + O        (λ < 424 nm)     [R1]
    O + O₂ + M → O₃ + M                        [R2]
    NO + O₃ → NO₂ + O₂                         [R3]

In the absence of other reactions, this cycle leads to the photostationary state
relationship between O₃, NO, and NO₂.

Reference: Seinfeld & Pandis (2006), Section 6.2, pp. 207-212
"""

using ModelingToolkit
using Unitful
using ModelingToolkit: t_nounits as t, D_nounits as D

# ============================================================================
# Equation 6.5: Ground-State Oxygen Atom Steady-State
# ============================================================================
# [O] = j_{NO₂} [NO₂] / (k₂[O₂][M])
#
# The concentration of ground-state O atoms is determined by the balance between
# production from NO₂ photolysis and loss by reaction with O₂.

# ============================================================================
# Equation 6.6: Photostationary State (Leighton Relationship)
# ============================================================================
# [O₃] = j_{NO₂} [NO₂] / (k_{NO+O₃} [NO])
#
# This is the fundamental photostationary state relationship, also known as
# the Leighton relationship. It predicts the ozone concentration when the
# NO-NO₂-O₃ system is in photochemical equilibrium.

# ============================================================================
# Equation 6.7: Photostationary State Parameter (Φ)
# ============================================================================
# Φ = j_{NO₂} [NO₂] / (k_{NO+O₃} [NO] [O₃])
#
# In photostationary state, Φ = 1. Deviations from unity indicate the presence
# of additional oxidants (Φ > 1) or reductants (Φ < 1) perturbing the cycle.

# ============================================================================
# Equation 6.8: Net O₃ Production Rate
# ============================================================================
# d[O₃]/dt = j_{NO₂}[NO₂] - k_{NO+O₃}[NO][O₃]
#
# When the system is not in photostationary state, this gives the net rate
# of ozone change.

"""
    NOxPhotochemistry(; name)

ModelingToolkit System for the NOx photochemical cycle.

Implements the photostationary state relationships (Eqs. 6.5-6.8) from
Seinfeld & Pandis Chapter 6.

# Input Variables
- `NO`: Nitric oxide concentration [molecules cm⁻³]
- `NO2`: Nitrogen dioxide concentration [molecules cm⁻³]
- `O3`: Ozone concentration [molecules cm⁻³]
- `O2`: Molecular oxygen concentration [molecules cm⁻³]
- `M`: Total air number density [molecules cm⁻³]

# Output Variables
- `O`: Ground-state oxygen atom concentration [molecules cm⁻³]
- `O3_pss`: Photostationary state ozone concentration [molecules cm⁻³]
- `Φ`: Photostationary state parameter [dimensionless]
- `P_O3`: Net ozone production rate [molecules cm⁻³ s⁻¹]

# Parameters
- `j_NO2`: NO₂ photolysis rate [s⁻¹]
- `k_O_O2_M`: Rate constant for O + O₂ + M → O₃ [cm⁶ molecule⁻² s⁻¹]
- `k_NO_O3`: Rate constant for NO + O₃ → NO₂ + O₂ [cm³ molecule⁻¹ s⁻¹]

# Rate Constants at 298 K (from Table B.1)
- j_NO2 ≈ 8 × 10⁻³ s⁻¹ (typical midday value)
- k_O_O2_M = 6.0 × 10⁻³⁴ cm⁶ molecule⁻² s⁻¹
- k_NO_O3 = 1.8 × 10⁻¹⁴ cm³ molecule⁻¹ s⁻¹
"""
function NOxPhotochemistry(; name=:NOxPhotochemistry)
    # Parameters
    @parameters begin
        j_NO2 = 8e-3, [description = "NO₂ photolysis rate [s⁻¹]"]
        k_O_O2_M = 6.0e-34, [description = "O + O₂ + M → O₃ rate [cm⁶/s]"]
        k_NO_O3 = 1.8e-14, [description = "NO + O₃ → NO₂ rate [cm³/s]"]
    end

    # Input variables
    @variables begin
        NO(t), [description = "NO concentration [molecules/cm³]"]
        NO2(t), [description = "NO₂ concentration [molecules/cm³]"]
        O3(t), [description = "O₃ concentration [molecules/cm³]"]
        O2(t), [description = "O₂ concentration [molecules/cm³]"]
        M(t), [description = "Total air number density [molecules/cm³]"]
    end

    # Output variables
    @variables begin
        O(t), [description = "O atom concentration [molecules/cm³]"]
        O3_pss(t), [description = "Photostationary state O₃ [molecules/cm³]"]
        Φ(t), [description = "Photostationary state parameter"]
        P_O3(t), [description = "Net O₃ production rate [molecules/cm³/s]"]
    end

    # Equations
    eqs = [
        # Equation 6.5: O atom steady-state concentration
        O ~ j_NO2 * NO2 / (k_O_O2_M * O2 * M),

        # Equation 6.6: Photostationary state ozone concentration
        O3_pss ~ j_NO2 * NO2 / (k_NO_O3 * NO),

        # Equation 6.7: Photostationary state parameter
        Φ ~ j_NO2 * NO2 / (k_NO_O3 * NO * O3),

        # Equation 6.8: Net ozone production/loss rate
        P_O3 ~ j_NO2 * NO2 - k_NO_O3 * NO * O3,
    ]

    return ODESystem(eqs, t; name=name)
end

"""
    PhotostationaryState(; name)

Simplified model for analyzing photostationary state deviations.

This system calculates the deviation from photostationary state (Φ - 1),
which indicates the net effect of peroxy radical chemistry on the NO-NO₂-O₃ cycle.

When Φ > 1: Additional oxidants (HO₂, RO₂) are converting NO to NO₂
When Φ < 1: Additional reductants are present
When Φ = 1: Pure photostationary state (no net ozone production)
"""
function PhotostationaryState(; name=:PhotostationaryState)
    # Parameters
    @parameters begin
        j_NO2 = 8e-3, [description = "NO₂ photolysis rate [s⁻¹]"]
        k_NO_O3 = 1.8e-14, [description = "NO + O₃ → NO₂ rate [cm³/s]"]
    end

    # Input variables
    @variables begin
        NO(t), [description = "NO concentration [molecules/cm³]"]
        NO2(t), [description = "NO₂ concentration [molecules/cm³]"]
        O3(t), [description = "O₃ concentration [molecules/cm³]"]
    end

    # Output variables
    @variables begin
        Φ(t), [description = "Photostationary state parameter"]
        Φ_deviation(t), [description = "Deviation from photostationary state"]
    end

    eqs = [
        Φ ~ j_NO2 * NO2 / (k_NO_O3 * NO * O3),
        Φ_deviation ~ Φ - 1,
    ]

    return ODESystem(eqs, t; name=name)
end

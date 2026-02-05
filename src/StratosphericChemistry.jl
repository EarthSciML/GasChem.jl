# ============================================================================
# Stratospheric Chemistry
# ============================================================================
"""
Stratospheric ozone chemistry based on Chapter 5 of Seinfeld & Pandis (2006)
"Atmospheric Chemistry and Physics: From Air Pollution to Climate Change", 2nd Edition.

This file implements:
- Chapman mechanism for stratospheric ozone (Equations 5.1-5.17)
- NOx catalytic cycles (Equations 5.20-5.22)
- HOx catalytic cycles (Equations 5.23-5.28)
- ClOx catalytic cycles (Equations 5.29-5.30)
- BrOx catalytic cycles
- Heterogeneous chemistry on polar stratospheric clouds

Reference: Seinfeld, J.H. and Pandis, S.N. (2006), Chapter 5, pp. 138-203.
"""

# Export main systems
export ChapmanMechanism, NOxCycle, HOxCycle, ClOxCycle, BrOxCycle
export StratosphericOzoneSystem
export stratospheric_rate_coefficients

# ============================================================================
# Rate Coefficient Functions (Table B.1 and B.2 from Seinfeld & Pandis)
# ============================================================================

"""
    k_O_O2_M(T, M)

Rate coefficient for O + O2 + M → O3 + M (Reaction 2 in Chapman mechanism)
k2 = 6.0 × 10⁻³⁴ (T/300)⁻²·⁴ cm⁶ molecule⁻² s⁻¹

Reference: Table B.2, Seinfeld & Pandis (2006)
"""
function k_O_O2_M(T)
    return 6.0e-34 * (T / 300.0)^(-2.4)
end

"""
    k_O_O3(T)

Rate coefficient for O + O3 → O2 + O2 (Reaction 4 in Chapman mechanism)
k4 = 8.0 × 10⁻¹² exp(-2060/T) cm³ molecule⁻¹ s⁻¹

Reference: Table B.1, Seinfeld & Pandis (2006)
"""
function k_O_O3(T)
    return 8.0e-12 * exp(-2060.0 / T)
end

"""
    k_O1D_M(T, M_type)

Rate coefficient for O(¹D) + M → O + M (quenching)
For M = O2: k = 3.2 × 10⁻¹¹ exp(70/T) cm³ molecule⁻¹ s⁻¹
For M = N2: k = 1.8 × 10⁻¹¹ exp(110/T) cm³ molecule⁻¹ s⁻¹

Reference: Page 143, Seinfeld & Pandis (2006)
"""
function k_O1D_M(T, M_type::Symbol)
    if M_type == :O2
        return 3.2e-11 * exp(70.0 / T)
    elseif M_type == :N2
        return 1.8e-11 * exp(110.0 / T)
    else
        # Weighted average for air (0.21 O2 + 0.79 N2)
        k_O2 = 3.2e-11 * exp(70.0 / T)
        k_N2 = 1.8e-11 * exp(110.0 / T)
        return 0.21 * k_O2 + 0.79 * k_N2
    end
end

"""
    k_NO_O3(T)

Rate coefficient for NO + O3 → NO2 + O2
k = 3.0 × 10⁻¹² exp(-1500/T) cm³ molecule⁻¹ s⁻¹

Reference: Page 154, Seinfeld & Pandis (2006)
"""
function k_NO_O3(T)
    return 3.0e-12 * exp(-1500.0 / T)
end

"""
    k_NO2_O(T)

Rate coefficient for NO2 + O → NO + O2
k = 5.6 × 10⁻¹² exp(180/T) cm³ molecule⁻¹ s⁻¹

Reference: Page 154, Seinfeld & Pandis (2006)
"""
function k_NO2_O(T)
    return 5.6e-12 * exp(180.0 / T)
end

"""
    k_OH_O3(T)

Rate coefficient for OH + O3 → HO2 + O2
k = 1.7 × 10⁻¹² exp(-940/T) cm³ molecule⁻¹ s⁻¹

Reference: Page 161, Seinfeld & Pandis (2006)
"""
function k_OH_O3(T)
    return 1.7e-12 * exp(-940.0 / T)
end

"""
    k_HO2_O3(T)

Rate coefficient for HO2 + O3 → OH + O2 + O2
k ≈ 10⁻¹⁵ cm³ molecule⁻¹ s⁻¹

Reference: Page 160, Seinfeld & Pandis (2006)
"""
function k_HO2_O3(T)
    return 1.0e-15
end

"""
    k_HO2_O(T)

Rate coefficient for HO2 + O → OH + O2
k ≈ 3.0 × 10⁻¹¹ exp(200/T) cm³ molecule⁻¹ s⁻¹

Reference: Page 161, Seinfeld & Pandis (2006)
"""
function k_HO2_O(T)
    return 3.0e-11 * exp(200.0 / T)
end

"""
    k_HO2_NO(T)

Rate coefficient for HO2 + NO → NO2 + OH
k = 3.5 × 10⁻¹² exp(250/T) cm³ molecule⁻¹ s⁻¹

Reference: Page 158, Seinfeld & Pandis (2006)
"""
function k_HO2_NO(T)
    return 3.5e-12 * exp(250.0 / T)
end

"""
    k_Cl_O3(T)

Rate coefficient for Cl + O3 → ClO + O2
k = 2.3 × 10⁻¹¹ exp(-200/T) cm³ molecule⁻¹ s⁻¹

Reference: Page 162, Seinfeld & Pandis (2006)
"""
function k_Cl_O3(T)
    return 2.3e-11 * exp(-200.0 / T)
end

"""
    k_ClO_O(T)

Rate coefficient for ClO + O → Cl + O2
k = 3.0 × 10⁻¹¹ exp(70/T) cm³ molecule⁻¹ s⁻¹

Reference: Page 162, Seinfeld & Pandis (2006)
"""
function k_ClO_O(T)
    return 3.0e-11 * exp(70.0 / T)
end

"""
    k_ClO_NO(T)

Rate coefficient for ClO + NO → Cl + NO2
k = 6.4 × 10⁻¹² exp(290/T) cm³ molecule⁻¹ s⁻¹

Reference: Page 163, Seinfeld & Pandis (2006)
"""
function k_ClO_NO(T)
    return 6.4e-12 * exp(290.0 / T)
end

"""
    k_Cl_CH4(T)

Rate coefficient for Cl + CH4 → HCl + CH3
k ≈ 10⁻¹⁴ cm³ molecule⁻¹ s⁻¹ at 260 K

Reference: Page 169 footnote, Seinfeld & Pandis (2006)
"""
function k_Cl_CH4(T)
    return 1.0e-14  # Approximate value at stratospheric temperatures
end

"""
    k_OH_HCl(T)

Rate coefficient for OH + HCl → H2O + Cl
k = 2.6 × 10⁻¹² exp(-350/T) cm³ molecule⁻¹ s⁻¹

Reference: Page 168, Seinfeld & Pandis (2006)
"""
function k_OH_HCl(T)
    return 2.6e-12 * exp(-350.0 / T)
end

"""
    k_Br_O3(T)

Rate coefficient for Br + O3 → BrO + O2
k ≈ 7 × 10⁻¹³ cm³ molecule⁻¹ s⁻¹ at 250 K

Reference: Page 198 (Problem 5.11), Seinfeld & Pandis (2006)
"""
function k_Br_O3(T)
    return 7.0e-13
end

"""
    k_BrO_ClO(T)

Rate coefficient for BrO + ClO → products
Multiple channels with different products

Reference: Page 166, Seinfeld & Pandis (2006)
"""
function k_BrO_ClO_BrCl(T)
    return 2.0e-12  # → BrCl + O2
end

function k_BrO_ClO_ClOO(T)
    return 8.0e-12  # → ClOO + Br
end

"""
    k_ClO_ClO_M(T, M)

Rate coefficient for ClO + ClO + M → Cl2O2 + M
Termolecular reaction forming the ClO dimer

Reference: Page 172, Seinfeld & Pandis (2006)
"""
function k_ClO_ClO_M(T, M)
    # Troe formula approximation
    k0 = 1.6e-32 * (T / 300.0)^(-4.5)
    kinf = 2.0e-12 * (T / 300.0)^(-2.4)
    return k0 * M / (1.0 + k0 * M / kinf)
end

"""
    k_N2O5_H2O_het(gamma, T, Ap)

First-order rate coefficient for heterogeneous hydrolysis of N2O5
N2O5 + H2O(s) → 2 HNO3

k = (γ/4) × (8kT/πm)^(1/2) × Ap

Reference: Equation 5.31, Seinfeld & Pandis (2006)
"""
function k_N2O5_H2O_het(gamma, T, Ap)
    # Molecular mass of N2O5 in kg
    m_N2O5 = 108.0 * 1.66054e-27  # kg per molecule
    k_B = 1.38065e-23  # J/K
    mean_speed = sqrt(8.0 * k_B * T / (π * m_N2O5))  # m/s
    mean_speed_cm = mean_speed * 100.0  # cm/s
    return (gamma / 4.0) * mean_speed_cm * Ap
end

# ============================================================================
# Collect all rate coefficients into a dictionary
# ============================================================================

"""
    stratospheric_rate_coefficients(T, M)

Return a dictionary of all stratospheric rate coefficients at temperature T
and total number density M.
"""
function stratospheric_rate_coefficients(T, M)
    return Dict(
        :k_O_O2_M => k_O_O2_M(T),
        :k_O_O3 => k_O_O3(T),
        :k_O1D_M => k_O1D_M(T, :air),
        :k_NO_O3 => k_NO_O3(T),
        :k_NO2_O => k_NO2_O(T),
        :k_OH_O3 => k_OH_O3(T),
        :k_HO2_O3 => k_HO2_O3(T),
        :k_HO2_O => k_HO2_O(T),
        :k_HO2_NO => k_HO2_NO(T),
        :k_Cl_O3 => k_Cl_O3(T),
        :k_ClO_O => k_ClO_O(T),
        :k_ClO_NO => k_ClO_NO(T),
        :k_Cl_CH4 => k_Cl_CH4(T),
        :k_OH_HCl => k_OH_HCl(T),
        :k_Br_O3 => k_Br_O3(T),
        :k_BrO_ClO_BrCl => k_BrO_ClO_BrCl(T),
        :k_BrO_ClO_ClOO => k_BrO_ClO_ClOO(T),
        :k_ClO_ClO_M => k_ClO_ClO_M(T, M)
    )
end

# ============================================================================
# ModelingToolkit Systems
# ============================================================================

"""
    ChapmanMechanism(; name=:ChapmanMechanism)

Create a ModelingToolkit System for the Chapman mechanism.

The Chapman mechanism describes the basic production and destruction of ozone
in the stratosphere through photolysis of O2 and subsequent reactions.

## Reactions (Section 5.2, Seinfeld & Pandis 2006):

 1. O2 + hν → O + O                    (j_O2)
 2. O + O2 + M → O3 + M                (k2)
 3. O3 + hν → O + O2                   (j_O3)
 4. O + O3 → O2 + O2                   (k4)

## Rate Equations (Equations 5.1-5.2):

d[O]/dt = 2j_O2[O2] - k2[O][O2][M] + j_O3[O3] - k4[O][O3]
d[O3]/dt = k2[O][O2][M] - j_O3[O3] - k4[O][O3]

## Steady-State Ozone (Equation 5.13):

[O3]_ss = 0.21 × (k2 × j_O2 / (k4 × j_O3))^(1/2) × [M]^(3/2)
"""
@component function ChapmanMechanism(; name = :ChapmanMechanism)
    @parameters begin
        j_O2 = 1e-10, [unit = u"s^-1", description = "O2 photolysis rate"]
        j_O3 = 1e-3, [unit = u"s^-1", description = "O3 photolysis rate"]
        k2 = 6e-46,
        [unit = u"m^6/s",
            description = "O + O2 + M rate coefficient (Eq. 5.2, 6e-34 cm^6/molec^2/s)"]
        k4 = 8e-18,
        [
            unit = u"m^3/s", description = "O + O3 rate coefficient (Eq. 5.4, 8e-12 cm^3/molec/s)"]
        M = 3e23, [unit = u"m^-3", description = "Air number density (3e17 molec/cm^3)"]
        O2_mix = 0.21, [unit = u"1", description = "O2 mixing ratio (dimensionless)"]
    end

    @variables begin
        O(t) = 1e13,
        [unit = u"m^-3", description = "Atomic oxygen concentration (1e7 molec/cm^3)"]
        O3(t) = 3e18,
        [unit = u"m^-3", description = "Ozone concentration (3e12 molec/cm^3)"]
        Ox(t), [unit = u"m^-3", description = "Odd oxygen = O + O3"]
    end

    O2 = O2_mix * M

    eqs = [
        D(O) ~ 2 * j_O2 * O2 - k2 * O * O2 * M + j_O3 * O3 - k4 * O * O3,  # Eq. 5.1
        D(O3) ~ k2 * O * O2 * M - j_O3 * O3 - k4 * O * O3,  # Eq. 5.2
        Ox ~ O + O3  # Odd oxygen family
    ]

    return System(eqs, t; name)
end

"""
    NOxCycle(; name=:NOxCycle)

Create a ModelingToolkit System for the NOx catalytic ozone destruction cycle.

## NOx Cycle 1 (Page 154):

NO + O3 → NO2 + O2     (k1)
NO2 + O → NO + O2      (k2)
Net: O3 + O → O2 + O2

## NOx Source from N2O (Reaction 2a, Page 151):

N2O + O(¹D) → NO + NO  (k = 6.7 × 10⁻¹¹ cm³ molecule⁻¹ s⁻¹)

## Rate of Odd Oxygen Destruction (Equation 5.22):

d[Ox]/dt = -2 k2[NO2][O]
"""
@component function NOxCycle(; name = :NOxCycle)
    @parameters begin
        k_NO_O3 = 3e-18,
        [unit = u"m^3/s", description = "NO + O3 rate coefficient (3e-12 cm^3/molec/s)"]
        k_NO2_O = 5.6e-18,
        [unit = u"m^3/s", description = "NO2 + O rate coefficient (5.6e-12 cm^3/molec/s)"]
        j_NO2 = 1e-2, [unit = u"s^-1", description = "NO2 photolysis rate"]
        O = 1e13, [unit = u"m^-3", description = "Atomic oxygen (1e7 molec/cm^3)"]
        O3 = 3e18, [unit = u"m^-3", description = "Ozone (3e12 molec/cm^3)"]
    end

    @variables begin
        NO(t) = 1e15, [unit = u"m^-3", description = "NO concentration (1e9 molec/cm^3)"]
        NO2(t) = 1e15, [unit = u"m^-3", description = "NO2 concentration (1e9 molec/cm^3)"]
        NOx(t), [unit = u"m^-3", description = "NOx = NO + NO2"]
    end

    eqs = [
        D(NO) ~ -k_NO_O3 * NO * O3 + k_NO2_O * NO2 * O + j_NO2 * NO2,
        D(NO2) ~ k_NO_O3 * NO * O3 - k_NO2_O * NO2 * O - j_NO2 * NO2,
        NOx ~ NO + NO2  # NOx family
    ]

    return System(eqs, t; name)
end

"""
    HOxCycle(; name=:HOxCycle)

Create a ModelingToolkit System for the HOx catalytic ozone destruction cycle.

## HOx Cycle 1 (Page 159):

OH + O3 → HO2 + O2
HO2 + O → OH + O2
Net: O3 + O → O2 + O2

## HOx Cycle 2 (Page 159):

OH + O3 → HO2 + O2
HO2 + O3 → OH + O2 + O2
Net: O3 + O3 → O2 + O2 + O2

## Steady-State Ratio (Equation 5.28):

[HO2]/[OH] = k_OH+O3[O3] / (k_HO2+NO[NO])
"""
@component function HOxCycle(; name = :HOxCycle)
    @parameters begin
        k_OH_O3 = 1.7e-18,
        [unit = u"m^3/s", description = "OH + O3 rate coefficient (1.7e-12 cm^3/molec/s)"]
        k_HO2_O3 = 1e-21,
        [unit = u"m^3/s", description = "HO2 + O3 rate coefficient (1e-15 cm^3/molec/s)"]
        k_HO2_O = 3e-17,
        [unit = u"m^3/s", description = "HO2 + O rate coefficient (3e-11 cm^3/molec/s)"]
        k_HO2_NO = 3.5e-18,
        [unit = u"m^3/s", description = "HO2 + NO rate coefficient (3.5e-12 cm^3/molec/s)"]
        O = 1e13, [unit = u"m^-3", description = "Atomic oxygen (1e7 molec/cm^3)"]
        O3 = 3e18, [unit = u"m^-3", description = "Ozone (3e12 molec/cm^3)"]
        NO = 1e15, [unit = u"m^-3", description = "NO (1e9 molec/cm^3)"]
    end

    @variables begin
        OH(t) = 1e12, [unit = u"m^-3", description = "OH concentration (1e6 molec/cm^3)"]
        HO2(t) = 1e13, [unit = u"m^-3", description = "HO2 concentration (1e7 molec/cm^3)"]
        HOx(t), [unit = u"m^-3", description = "HOx = OH + HO2"]
    end

    eqs = [
        D(OH) ~
        -k_OH_O3 * OH * O3 + k_HO2_O * HO2 * O + k_HO2_O3 * HO2 * O3 + k_HO2_NO * HO2 * NO,
        D(HO2) ~
        k_OH_O3 * OH * O3 - k_HO2_O * HO2 * O - k_HO2_O3 * HO2 * O3 - k_HO2_NO * HO2 * NO,
        HOx ~ OH + HO2  # HOx family
    ]

    return System(eqs, t; name)
end

"""
    ClOxCycle(; name=:ClOxCycle)

Create a ModelingToolkit System for the ClOx catalytic ozone destruction cycle.

## ClOx Cycle 1 (Page 162):

Cl + O3 → ClO + O2    (k1 = 2.3 × 10⁻¹¹ exp(-200/T))
ClO + O → Cl + O2     (k2 = 3.0 × 10⁻¹¹ exp(70/T))
Net: O3 + O → O2 + O2

## Rate of Odd Oxygen Destruction (Equation 5.29):

d[Ox]/dt = -2 k2[ClO][O]

## Steady-State [Cl]/[ClO] Ratio (Equation 5.30):

[Cl]/[ClO] = (k_ClO+O[O] + k_ClO+NO[NO]) / (k_Cl+O3[O3])

## Reservoir Species:

  - HCl formed by Cl + CH4 → HCl + CH3
  - ClONO2 formed by ClO + NO2 + M → ClONO2 + M
"""
@component function ClOxCycle(; name = :ClOxCycle)
    @parameters begin
        k_Cl_O3 = 2.3e-17,
        [unit = u"m^3/s", description = "Cl + O3 rate coefficient (2.3e-11 cm^3/molec/s)"]
        k_ClO_O = 3e-17,
        [unit = u"m^3/s", description = "ClO + O rate coefficient (3e-11 cm^3/molec/s)"]
        k_ClO_NO = 6.4e-18,
        [unit = u"m^3/s", description = "ClO + NO rate coefficient (6.4e-12 cm^3/molec/s)"]
        k_Cl_CH4 = 1e-20,
        [unit = u"m^3/s", description = "Cl + CH4 rate coefficient (1e-14 cm^3/molec/s)"]
        k_OH_HCl = 2.6e-18,
        [unit = u"m^3/s", description = "OH + HCl rate coefficient (2.6e-12 cm^3/molec/s)"]
        k_ClO_NO2_M = 1.8e-43,
        [unit = u"m^6/s", description = "ClO + NO2 + M → ClONO2 (1.8e-31 cm^6/molec^2/s)"]
        j_ClONO2 = 1e-4, [unit = u"s^-1", description = "ClONO2 photolysis rate"]
        CH4 = 1e19,
        [unit = u"m^-3", description = "Methane concentration (1e13 molec/cm^3)"]
        M = 3.1e23,
        [unit = u"m^-3", description = "Air number density at 30 km (3.1e17 molec/cm^3)"]
        O = 1e13, [unit = u"m^-3", description = "Atomic oxygen (1e7 molec/cm^3)"]
        O3 = 3e18, [unit = u"m^-3", description = "Ozone (3e12 molec/cm^3)"]
        NO = 1e15, [unit = u"m^-3", description = "NO (1e9 molec/cm^3)"]
        NO2 = 1e15, [unit = u"m^-3", description = "NO2 (1e9 molec/cm^3)"]
        OH = 1e12, [unit = u"m^-3", description = "OH (1e6 molec/cm^3)"]
    end

    @variables begin
        Cl(t) = 1e10,
        [unit = u"m^-3", description = "Cl atom concentration (1e4 molec/cm^3)"]
        ClO(t) = 1e13, [unit = u"m^-3", description = "ClO concentration (1e7 molec/cm^3)"]
        HCl(t) = 1e15, [unit = u"m^-3", description = "HCl reservoir (1e9 molec/cm^3)"]
        ClONO2(t) = 1e15,
        [unit = u"m^-3", description = "ClONO2 reservoir (1e9 molec/cm^3)"]
        ClOx(t), [unit = u"m^-3", description = "ClOx = Cl + ClO"]
        Cly(t), [unit = u"m^-3", description = "Cly = Cl + ClO + HCl + ClONO2"]
    end

    eqs = [
        D(Cl) ~
        -k_Cl_O3 * Cl * O3 + k_ClO_O * ClO * O + k_ClO_NO * ClO * NO -
        k_Cl_CH4 * Cl * CH4 + k_OH_HCl * OH * HCl + j_ClONO2 * ClONO2,
        D(ClO) ~ k_Cl_O3 * Cl * O3 - k_ClO_O * ClO * O - k_ClO_NO * ClO * NO -
                 k_ClO_NO2_M * ClO * NO2 * M,  # ClO + NO2 + M → ClONO2 (Page 167)
        D(HCl) ~ k_Cl_CH4 * Cl * CH4 - k_OH_HCl * OH * HCl,
        D(ClONO2) ~ k_ClO_NO2_M * ClO * NO2 * M -  # ClO + NO2 + M → ClONO2 (formation)
                    j_ClONO2 * ClONO2,               # ClONO2 + hν → Cl + NO3 (photolysis)
        ClOx ~ Cl + ClO,  # ClOx family
        Cly ~ Cl + ClO + HCl + ClONO2  # Cly family (total inorganic chlorine)
    ]

    return System(eqs, t; name)
end

"""
    BrOxCycle(; name=:BrOxCycle)

Create a ModelingToolkit System for the BrOx catalytic ozone destruction cycle.

Bromine is approximately 50 times more effective than chlorine in destroying ozone
on an atom-for-atom basis (Page 169).

## Key Reactions:

Br + O3 → BrO + O2
BrO + ClO → Br + Cl + O2  (or BrCl + O2)
BrO + HO2 → HOBr + O2

## Note on Reservoir Species:

Unlike chlorine, bromine does not form a stable HBr reservoir because
Br + CH4 is endothermic and extremely slow (Page 169).
"""
@component function BrOxCycle(; name = :BrOxCycle)
    @parameters begin
        k_Br_O3 = 7e-19,
        [unit = u"m^3/s", description = "Br + O3 rate coefficient (7e-13 cm^3/molec/s)"]
        k_BrO_O = 5e-17,
        [unit = u"m^3/s", description = "BrO + O rate coefficient (5e-11 cm^3/molec/s)"]
        k_BrO_ClO = 2e-18,
        [unit = u"m^3/s", description = "BrO + ClO rate coefficient (2e-12 cm^3/molec/s)"]
        k_BrO_HO2 = 4e-17,
        [unit = u"m^3/s", description = "BrO + HO2 rate coefficient (4e-11 cm^3/molec/s)"]
        j_HOBr = 1e-3, [unit = u"s^-1", description = "HOBr photolysis rate"]
        O = 1e13, [unit = u"m^-3", description = "Atomic oxygen (1e7 molec/cm^3)"]
        O3 = 3e18, [unit = u"m^-3", description = "Ozone (3e12 molec/cm^3)"]
        ClO = 1e13, [unit = u"m^-3", description = "ClO (1e7 molec/cm^3)"]
        HO2 = 1e13, [unit = u"m^-3", description = "HO2 (1e7 molec/cm^3)"]
    end

    @variables begin
        Br(t) = 1e11,
        [unit = u"m^-3", description = "Br atom concentration (1e5 molec/cm^3)"]
        BrO(t) = 1e12, [unit = u"m^-3", description = "BrO concentration (1e6 molec/cm^3)"]
        HOBr(t) = 1e12,
        [unit = u"m^-3", description = "HOBr concentration (1e6 molec/cm^3)"]
        BrOx(t), [unit = u"m^-3", description = "BrOx = Br + BrO"]
        Bry(t), [unit = u"m^-3", description = "Bry = Br + BrO + HOBr"]
    end

    eqs = [
        D(Br) ~
        -k_Br_O3 * Br * O3 + k_BrO_O * BrO * O + k_BrO_ClO * BrO * ClO + j_HOBr * HOBr,
        D(BrO) ~
        k_Br_O3 * Br * O3 - k_BrO_O * BrO * O - k_BrO_ClO * BrO * ClO -
        k_BrO_HO2 * BrO * HO2,
        D(HOBr) ~ k_BrO_HO2 * BrO * HO2 - j_HOBr * HOBr,
        BrOx ~ Br + BrO,  # BrOx family
        Bry ~ Br + BrO + HOBr  # Bry family
    ]

    return System(eqs, t; name)
end

"""
    StratosphericOzoneSystem(; name=:StratosphericOzoneSystem)

Create a comprehensive ModelingToolkit System combining all stratospheric
ozone chemistry cycles.

This system includes:

  - Chapman mechanism (O, O3 production and loss)
  - O(¹D) photochemistry
  - NOx cycle (catalytic O3 destruction)
  - HOx cycle (catalytic O3 destruction)
  - ClOx cycle (catalytic O3 destruction)
  - BrOx cycle (catalytic O3 destruction)

All rate coefficients are temperature-dependent, computed from Arrhenius
parameters defined as `@constants`.

## Key Equations from Seinfeld & Pandis Chapter 5:

### Odd Oxygen Balance (Equation 5.9):

d[Ox]/dt = 2j_O2[O2] - 2k4[O][O3]

### Steady-State O3 (Equation 5.13):

[O3]_ss = 0.21 × (k2 × j_O2 / (k4 × j_O3))^(1/2) × [M]^(3/2)

### [O]/[O3] Ratio (Equation 5.7):

[O]/[O3] = j_O3 / (k2[O2][M])

### Time to Steady State (Equation 5.17):

τ_O3^ss = (1/4) × (k2[M] / (k4 × j_O2 × j_O3))^(1/2)

"""
@component function StratosphericOzoneSystem(; name = :StratosphericOzoneSystem)
    # =========================================================================
    # Physical and kinetic constants (Arrhenius parameters)
    # Rate expressions: k = A * exp(C / T) or k = A * (T_ref / T)^n
    # =========================================================================
    @constants begin
        T_ref = 300.0, [unit = u"K", description = "Reference temperature"]

        # Chapman mechanism (Table B.1/B.2, Seinfeld & Pandis 2006)
        # CGS→SI: cm^6/molec^2/s × 1e-12 → m^6/s; cm^3/molec/s × 1e-6 → m^3/s
        k2_A = 6.0e-46,
        [
            unit = u"m^6/s", description = "Pre-factor: O + O2 + M → O3 + M (6e-34 cm^6/molec^2/s)"]
        k4_A = 8.0e-18,
        [unit = u"m^3/s", description = "Pre-factor: O + O3 → 2O2 (8e-12 cm^3/molec/s)"]
        C_k4 = -2060.0, [unit = u"K", description = "exp(C/T) factor: O + O3"]

        # O(1D) quenching (Page 143)
        k_O1D_O2_A = 3.2e-17,
        [unit = u"m^3/s", description = "Pre-factor: O(1D) + O2 (3.2e-11 cm^3/molec/s)"]
        C_O1D_O2 = 70.0, [unit = u"K", description = "exp(C/T) factor: O(1D) + O2"]
        k_O1D_N2_A = 1.8e-17,
        [unit = u"m^3/s", description = "Pre-factor: O(1D) + N2 (1.8e-11 cm^3/molec/s)"]
        C_O1D_N2 = 110.0, [unit = u"K", description = "exp(C/T) factor: O(1D) + N2"]
        k_O1D_H2O_c = 2.2e-16,
        [unit = u"m^3/s", description = "O(1D) + H2O → 2OH (2.2e-10 cm^3/molec/s)"]

        # NOx (Page 154)
        k_NO_O3_A = 3.0e-18,
        [unit = u"m^3/s", description = "Pre-factor: NO + O3 (3e-12 cm^3/molec/s)"]
        C_NO_O3 = -1500.0, [unit = u"K", description = "exp(C/T) factor: NO + O3"]
        k_NO2_O_A = 5.6e-18,
        [unit = u"m^3/s", description = "Pre-factor: NO2 + O (5.6e-12 cm^3/molec/s)"]
        C_NO2_O = 180.0, [unit = u"K", description = "exp(C/T) factor: NO2 + O"]

        # HOx (Pages 159-161)
        k_OH_O3_A = 1.7e-18,
        [unit = u"m^3/s", description = "Pre-factor: OH + O3 (1.7e-12 cm^3/molec/s)"]
        C_OH_O3 = -940.0, [unit = u"K", description = "exp(C/T) factor: OH + O3"]
        k_HO2_O3_c = 1.0e-21,
        [unit = u"m^3/s", description = "HO2 + O3 rate (1e-15 cm^3/molec/s)"]
        k_HO2_O_A = 3.0e-17,
        [unit = u"m^3/s", description = "Pre-factor: HO2 + O (3e-11 cm^3/molec/s)"]
        C_HO2_O = 200.0, [unit = u"K", description = "exp(C/T) factor: HO2 + O"]
        k_HO2_NO_A = 3.5e-18,
        [unit = u"m^3/s", description = "Pre-factor: HO2 + NO (3.5e-12 cm^3/molec/s)"]
        C_HO2_NO = 250.0, [unit = u"K", description = "exp(C/T) factor: HO2 + NO"]

        # ClOx (Pages 162-169)
        k_Cl_O3_A = 2.3e-17,
        [unit = u"m^3/s", description = "Pre-factor: Cl + O3 (2.3e-11 cm^3/molec/s)"]
        C_Cl_O3 = -200.0, [unit = u"K", description = "exp(C/T) factor: Cl + O3"]
        k_ClO_O_A = 3.0e-17,
        [unit = u"m^3/s", description = "Pre-factor: ClO + O (3e-11 cm^3/molec/s)"]
        C_ClO_O = 70.0, [unit = u"K", description = "exp(C/T) factor: ClO + O"]
        k_ClO_NO_A = 6.4e-18,
        [unit = u"m^3/s", description = "Pre-factor: ClO + NO (6.4e-12 cm^3/molec/s)"]
        C_ClO_NO = 290.0, [unit = u"K", description = "exp(C/T) factor: ClO + NO"]
        k_Cl_CH4_c = 1.0e-20,
        [unit = u"m^3/s", description = "Cl + CH4 rate (1e-14 cm^3/molec/s)"]
        k_OH_HCl_A = 2.6e-18,
        [unit = u"m^3/s", description = "Pre-factor: OH + HCl (2.6e-12 cm^3/molec/s)"]
        C_OH_HCl = -350.0, [unit = u"K", description = "exp(C/T) factor: OH + HCl"]

        # ClO + NO2 + M → ClONO2 + M (Page 165, termolecular)
        # Approximate effective bimolecular rate at stratospheric conditions
        # k ≈ 1.8e-31 (T/300)^(-3.4) cm^6/molec^2/s (JPL recommendation)
        # At 30 km (M=3.1e17): k_eff ≈ 1.8e-31 * 3.1e17 ≈ 5.6e-14 cm^3/molec/s
        k_ClO_NO2_M_A = 1.8e-43,
        [unit = u"m^6/s", description = "Pre-factor: ClO + NO2 + M → ClONO2 (1.8e-31 cm^6/molec^2/s)"]

        # BrOx (Pages 166-169)
        k_Br_O3_c = 7.0e-19,
        [unit = u"m^3/s", description = "Br + O3 rate (7e-13 cm^3/molec/s)"]
        k_BrO_O_c = 5.0e-17,
        [unit = u"m^3/s", description = "BrO + O rate (5e-11 cm^3/molec/s)"]
        k_BrO_ClO_c = 2.0e-18,
        [unit = u"m^3/s", description = "BrO + ClO rate (2e-12 cm^3/molec/s)"]
        k_BrO_HO2_c = 4.0e-17,
        [unit = u"m^3/s", description = "BrO + HO2 rate (4e-11 cm^3/molec/s)"]
    end

    @parameters begin
        # Photolysis rates
        j_O2 = 1e-10, [unit = u"s^-1", description = "O2 photolysis rate"]
        j_O3 = 1e-3,
        [unit = u"s^-1", description = "O3 photolysis rate (total, both channels)"]
        j_O3_O1D = 5e-4, [unit = u"s^-1", description = "O3 → O(1D) photolysis rate"]
        j_NO2 = 1e-2, [unit = u"s^-1", description = "NO2 photolysis rate"]
        j_ClONO2 = 1e-4, [unit = u"s^-1", description = "ClONO2 photolysis rate"]
        j_HOBr = 1e-3, [unit = u"s^-1", description = "HOBr photolysis rate"]

        # Environmental parameters
        T = 227.0, [unit = u"K", description = "Temperature"]
        M = 3.1e23,
        [unit = u"m^-3", description = "Air number density at 30 km (3.1e17 molec/cm^3)"]
        O2_mix = 0.21, [unit = u"1", description = "O2 mixing ratio (dimensionless)"]
        N2_mix = 0.79, [unit = u"1", description = "N2 mixing ratio (dimensionless)"]
        CH4_mix = 1.6e-6, [unit = u"1", description = "CH4 mixing ratio (dimensionless)"]
        H2O_mix = 5e-6, [unit = u"1", description = "H2O mixing ratio (dimensionless)"]
    end

    # Derived concentrations
    O2_conc = O2_mix * M
    N2_conc = N2_mix * M
    CH4_conc = CH4_mix * M
    H2O_conc = H2O_mix * M

    @variables begin
        # Odd oxygen family
        O(t) = 1e13, [unit = u"m^-3", description = "Atomic oxygen O(3P) (1e7 molec/cm^3)"]
        O1D(t) = 5e7, [unit = u"m^-3", description = "Excited oxygen O(1D) (50 molec/cm^3)"]
        O3(t) = 3e18, [unit = u"m^-3", description = "Ozone (3e12 molec/cm^3)"]

        # Nitrogen oxides
        NO(t) = 1e15, [unit = u"m^-3", description = "Nitric oxide (1e9 molec/cm^3)"]
        NO2(t) = 1e15, [unit = u"m^-3", description = "Nitrogen dioxide (1e9 molec/cm^3)"]

        # Hydrogen radicals
        OH(t) = 1e12, [unit = u"m^-3", description = "Hydroxyl radical (1e6 molec/cm^3)"]
        HO2(t) = 1e13,
        [unit = u"m^-3", description = "Hydroperoxyl radical (1e7 molec/cm^3)"]

        # Chlorine species
        Cl(t) = 1e10, [unit = u"m^-3", description = "Chlorine atom (1e4 molec/cm^3)"]
        ClO(t) = 1e13, [unit = u"m^-3", description = "Chlorine monoxide (1e7 molec/cm^3)"]
        HCl(t) = 1e15, [unit = u"m^-3", description = "Hydrogen chloride (1e9 molec/cm^3)"]
        ClONO2(t) = 1e15,
        [unit = u"m^-3", description = "Chlorine nitrate (1e9 molec/cm^3)"]

        # Bromine species
        Br(t) = 1e11, [unit = u"m^-3", description = "Bromine atom (1e5 molec/cm^3)"]
        BrO(t) = 1e12, [unit = u"m^-3", description = "Bromine monoxide (1e6 molec/cm^3)"]
        HOBr(t) = 1e12, [unit = u"m^-3", description = "Hypobromous acid (1e6 molec/cm^3)"]

        # Chemical families (algebraic)
        Ox(t), [unit = u"m^-3", description = "Odd oxygen = O + O3"]
        NOx(t), [unit = u"m^-3", description = "NOx = NO + NO2"]
        HOx(t), [unit = u"m^-3", description = "HOx = OH + HO2"]
        ClOx(t), [unit = u"m^-3", description = "ClOx = Cl + ClO"]
        BrOx(t), [unit = u"m^-3", description = "BrOx = Br + BrO"]
    end

    # =========================================================================
    # Rate coefficients (temperature-dependent symbolic expressions)
    # Convention: k = A * exp(C / T) where C includes the sign
    # =========================================================================

    # Chapman mechanism
    k2 = k2_A * (T_ref / T)^2.4              # O + O2 + M → O3 + M
    k4 = k4_A * exp(C_k4 / T)               # O + O3 → 2O2

    # O(1D) quenching (weighted average for air)
    k_O1D_M = O2_mix * k_O1D_O2_A * exp(C_O1D_O2 / T) +
              N2_mix * k_O1D_N2_A * exp(C_O1D_N2 / T)
    k_O1D_H2O = k_O1D_H2O_c

    # NOx
    k_NO_O3 = k_NO_O3_A * exp(C_NO_O3 / T)  # NO + O3
    k_NO2_O = k_NO2_O_A * exp(C_NO2_O / T)  # NO2 + O

    # HOx
    k_OH_O3 = k_OH_O3_A * exp(C_OH_O3 / T)  # OH + O3
    k_HO2_O3 = k_HO2_O3_c                   # HO2 + O3
    k_HO2_O = k_HO2_O_A * exp(C_HO2_O / T)  # HO2 + O
    k_HO2_NO = k_HO2_NO_A * exp(C_HO2_NO / T)  # HO2 + NO

    # ClOx
    k_Cl_O3 = k_Cl_O3_A * exp(C_Cl_O3 / T)  # Cl + O3
    k_ClO_O = k_ClO_O_A * exp(C_ClO_O / T)  # ClO + O
    k_ClO_NO = k_ClO_NO_A * exp(C_ClO_NO / T)  # ClO + NO
    k_Cl_CH4 = k_Cl_CH4_c                   # Cl + CH4
    k_OH_HCl = k_OH_HCl_A * exp(C_OH_HCl / T)  # OH + HCl

    # ClO + NO2 + M → ClONO2 (termolecular, approximate as T-dependent)
    k_ClO_NO2_M = k_ClO_NO2_M_A * (T_ref / T)^3.4  # ClO + NO2 + M → ClONO2 + M

    # BrOx
    k_Br_O3 = k_Br_O3_c                     # Br + O3
    k_BrO_O = k_BrO_O_c                     # BrO + O
    k_BrO_ClO = k_BrO_ClO_c                 # BrO + ClO
    k_BrO_HO2 = k_BrO_HO2_c                 # BrO + HO2

    eqs = [
        # =================================================================
        # Odd Oxygen Species
        # =================================================================

        # Atomic oxygen O(3P) — Eq. 5.1
        # Note: (j_O3 - j_O3_O1D) gives only the O(3P) channel of O3 photolysis
        D(O) ~
        2 * j_O2 * O2_conc +              # O2 + hν → 2O (source)
        (j_O3 - j_O3_O1D) * O3 +           # O3 + hν → O(3P) + O2 (source)
        k_O1D_M * O1D * M -                # O(1D) + M → O(3P) + M (source)
        k2 * O * O2_conc * M -             # O + O2 + M → O3 (sink)
        k4 * O * O3 -                      # O + O3 → 2O2 (sink)
        k_NO2_O * NO2 * O -                # NO2 + O → NO + O2 (sink)
        k_ClO_O * ClO * O -                # ClO + O → Cl + O2 (sink)
        k_BrO_O * BrO * O -                # BrO + O → Br + O2 (sink)
        k_HO2_O * HO2 * O,                 # HO2 + O → OH + O2 (sink)

        # Excited oxygen O(1D)
        D(O1D) ~ j_O3_O1D * O3 -                  # O3 + hν → O(1D) + O2 (source)
                 k_O1D_M * O1D * M -              # O(1D) + M → O(3P) (sink)
                 k_O1D_H2O * O1D * H2O_conc,      # O(1D) + H2O → 2OH (sink)

        # Ozone — Eq. 5.2
        D(O3) ~
        k2 * O * O2_conc * M -           # O + O2 + M → O3 (source)
        j_O3 * O3 -                       # O3 photolysis (sink, total both channels)
        k4 * O * O3 -                     # O + O3 → 2O2 (sink)
        k_NO_O3 * NO * O3 -               # NO + O3 → NO2 + O2 (sink)
        k_OH_O3 * OH * O3 -               # OH + O3 → HO2 + O2 (sink)
        k_HO2_O3 * HO2 * O3 -             # HO2 + O3 → OH + 2O2 (sink)
        k_Cl_O3 * Cl * O3 -               # Cl + O3 → ClO + O2 (sink)
        k_Br_O3 * Br * O3,                # Br + O3 → BrO + O2 (sink)

        # =================================================================
        # Nitrogen Oxide Species
        # =================================================================

        # NO — all terms verified for correct sign
        D(NO) ~
        j_NO2 * NO2 +                    # NO2 + hν → NO + O (produces NO)
        k_NO2_O * NO2 * O -               # NO2 + O → NO + O2 (produces NO)
        k_NO_O3 * NO * O3 -               # NO + O3 → NO2 + O2 (consumes NO)
        k_ClO_NO * ClO * NO -             # ClO + NO → Cl + NO2 (consumes NO)
        k_HO2_NO * HO2 * NO,              # HO2 + NO → NO2 + OH (consumes NO)

        # NO2
        D(NO2) ~
        k_NO_O3 * NO * O3 +             # NO + O3 → NO2 + O2 (produces NO2)
        k_HO2_NO * HO2 * NO +            # HO2 + NO → NO2 + OH (produces NO2)
        k_ClO_NO * ClO * NO -             # ClO + NO → Cl + NO2 (produces NO2)
        j_NO2 * NO2 -                    # NO2 + hν (consumes NO2)
        k_NO2_O * NO2 * O -              # NO2 + O → NO + O2 (consumes NO2)
        k_ClO_NO2_M * ClO * NO2 * M,     # ClO + NO2 + M → ClONO2 + M (consumes NO2, Page 167)

        # =================================================================
        # HOx Species
        # =================================================================

        # OH
        D(OH) ~
        2 * k_O1D_H2O * O1D * H2O_conc + # O(1D) + H2O → 2OH (source)
        k_HO2_O * HO2 * O +               # HO2 + O → OH + O2
        k_HO2_O3 * HO2 * O3 +             # HO2 + O3 → OH + 2O2
        k_HO2_NO * HO2 * NO +             # HO2 + NO → NO2 + OH
        j_HOBr * HOBr -                   # HOBr + hν → OH + Br (produces OH)
        k_OH_O3 * OH * O3 -               # OH + O3 → HO2 + O2
        k_OH_HCl * OH * HCl,              # OH + HCl → H2O + Cl

        # HO2
        D(HO2) ~
        k_OH_O3 * OH * O3 -             # OH + O3 → HO2 + O2
        k_HO2_O * HO2 * O -              # HO2 + O → OH + O2
        k_HO2_O3 * HO2 * O3 -            # HO2 + O3 → OH + 2O2
        k_HO2_NO * HO2 * NO -            # HO2 + NO → NO2 + OH
        k_BrO_HO2 * BrO * HO2,           # BrO + HO2 → HOBr + O2

        # =================================================================
        # Chlorine Species
        # =================================================================

        # Cl
        D(Cl) ~
        k_ClO_O * ClO * O +              # ClO + O → Cl + O2
        k_ClO_NO * ClO * NO +             # ClO + NO → Cl + NO2
        k_OH_HCl * OH * HCl +             # OH + HCl → H2O + Cl
        j_ClONO2 * ClONO2 +               # ClONO2 + hν → Cl + ...
        k_BrO_ClO * BrO * ClO -           # BrO + ClO → Br + Cl + O2
        k_Cl_O3 * Cl * O3 -               # Cl + O3 → ClO + O2
        k_Cl_CH4 * Cl * CH4_conc,         # Cl + CH4 → HCl + CH3

        # ClO
        D(ClO) ~
        k_Cl_O3 * Cl * O3 -             # Cl + O3 → ClO + O2
        k_ClO_O * ClO * O -              # ClO + O → Cl + O2
        k_ClO_NO * ClO * NO -            # ClO + NO → Cl + NO2
        k_BrO_ClO * BrO * ClO -          # BrO + ClO → products
        k_ClO_NO2_M * ClO * NO2 * M,     # ClO + NO2 + M → ClONO2 + M (Page 165)

        # HCl reservoir
        D(HCl) ~ k_Cl_CH4 * Cl * CH4_conc -      # Cl + CH4 → HCl + CH3
                 k_OH_HCl * OH * HCl,             # OH + HCl → H2O + Cl

        # ClONO2 reservoir (Page 167-168)
        D(ClONO2) ~ k_ClO_NO2_M * ClO * NO2 * M -  # ClO + NO2 + M → ClONO2 + M (formation)
                    j_ClONO2 * ClONO2,               # ClONO2 + hν → Cl + NO3 (photolysis)

        # =================================================================
        # Bromine Species
        # =================================================================

        # Br
        D(Br) ~
        k_BrO_O * BrO * O +              # BrO + O → Br + O2
        k_BrO_ClO * BrO * ClO +           # BrO + ClO → Br + Cl + O2
        j_HOBr * HOBr -                   # HOBr + hν → OH + Br
        k_Br_O3 * Br * O3,                # Br + O3 → BrO + O2

        # BrO
        D(BrO) ~
        k_Br_O3 * Br * O3 -             # Br + O3 → BrO + O2
        k_BrO_O * BrO * O -              # BrO + O → Br + O2
        k_BrO_ClO * BrO * ClO -          # BrO + ClO → products
        k_BrO_HO2 * BrO * HO2,           # BrO + HO2 → HOBr + O2

        # HOBr
        D(HOBr) ~ k_BrO_HO2 * BrO * HO2 -        # BrO + HO2 → HOBr + O2
                  j_HOBr * HOBr,                  # HOBr + hν → OH + Br

        # =================================================================
        # Chemical Families (Algebraic)
        # =================================================================

        Ox ~ O + O3,
        NOx ~ NO + NO2,
        HOx ~ OH + HO2,
        ClOx ~ Cl + ClO,
        BrOx ~ Br + BrO
    ]

    return System(eqs, t; name)
end

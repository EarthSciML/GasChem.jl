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
    ChapmanMechanism(; name=:chapman)

Create a ModelingToolkit ODESystem for the Chapman mechanism.

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
function ChapmanMechanism(; name=:chapman)
    @parameters begin
        t
        j_O2 = 1e-10, [description = "O2 photolysis rate (s^-1)"]
        j_O3 = 1e-3, [description = "O3 photolysis rate (s^-1)"]
        k2 = 6e-34, [description = "O + O2 + M rate coefficient (cm^6/molec^2/s)"]
        k4 = 8e-12, [description = "O + O3 rate coefficient (cm^3/molec/s)"]
        M = 3e17, [description = "Air number density (molec/cm^3)"]
        O2_conc = 0.21, [description = "O2 mixing ratio"]
    end

    @variables begin
        O(t) = 1e7, [description = "Atomic oxygen concentration (molec/cm^3)"]
        O3(t) = 3e12, [description = "Ozone concentration (molec/cm^3)"]
        Ox(t), [description = "Odd oxygen = O + O3 (molec/cm^3)"]
    end

    # O2 concentration from mixing ratio and M
    O2 = O2_conc * M

    eqs = [
        # Equation 5.1: d[O]/dt
        D(O) ~ 2 * j_O2 * O2 - k2 * O * O2 * M + j_O3 * O3 - k4 * O * O3,

        # Equation 5.2: d[O3]/dt
        D(O3) ~ k2 * O * O2 * M - j_O3 * O3 - k4 * O * O3,

        # Odd oxygen family
        Ox ~ O + O3
    ]

    return ODESystem(eqs, t; name=name)
end

"""
    NOxCycle(; name=:nox)

Create a ModelingToolkit ODESystem for the NOx catalytic ozone destruction cycle.

## NOx Cycle 1 (Page 154):
NO + O3 → NO2 + O2     (k1)
NO2 + O → NO + O2      (k2)
Net: O3 + O → O2 + O2

## NOx Source from N2O (Reaction 2a, Page 151):
N2O + O(¹D) → NO + NO  (k = 6.7 × 10⁻¹¹ cm³ molecule⁻¹ s⁻¹)

## Rate of Odd Oxygen Destruction (Equation 5.22):
d[Ox]/dt = -2 k2[NO2][O]
"""
function NOxCycle(; name=:nox)
    @parameters begin
        t
        k_NO_O3 = 3e-12, [description = "NO + O3 rate (cm^3/molec/s)"]
        k_NO2_O = 5.6e-12, [description = "NO2 + O rate (cm^3/molec/s)"]
        j_NO2 = 1e-2, [description = "NO2 photolysis rate (s^-1)"]
    end

    @variables begin
        NO(t) = 1e9, [description = "NO concentration (molec/cm^3)"]
        NO2(t) = 1e9, [description = "NO2 concentration (molec/cm^3)"]
        O(t) = 1e7, [description = "Atomic oxygen (input) (molec/cm^3)"]
        O3(t) = 3e12, [description = "Ozone (input) (molec/cm^3)"]
        NOx(t), [description = "NOx = NO + NO2 (molec/cm^3)"]
    end

    eqs = [
        # NO rate equation
        D(NO) ~ -k_NO_O3 * NO * O3 + k_NO2_O * NO2 * O + j_NO2 * NO2,

        # NO2 rate equation
        D(NO2) ~ k_NO_O3 * NO * O3 - k_NO2_O * NO2 * O - j_NO2 * NO2,

        # NOx family
        NOx ~ NO + NO2
    ]

    return ODESystem(eqs, t; name=name)
end

"""
    HOxCycle(; name=:hox)

Create a ModelingToolkit ODESystem for the HOx catalytic ozone destruction cycle.

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
function HOxCycle(; name=:hox)
    @parameters begin
        t
        k_OH_O3 = 1.7e-12, [description = "OH + O3 rate (cm^3/molec/s)"]
        k_HO2_O3 = 1e-15, [description = "HO2 + O3 rate (cm^3/molec/s)"]
        k_HO2_O = 3e-11, [description = "HO2 + O rate (cm^3/molec/s)"]
        k_HO2_NO = 3.5e-12, [description = "HO2 + NO rate (cm^3/molec/s)"]
    end

    @variables begin
        OH(t) = 1e6, [description = "OH concentration (molec/cm^3)"]
        HO2(t) = 1e7, [description = "HO2 concentration (molec/cm^3)"]
        O(t) = 1e7, [description = "Atomic oxygen (input) (molec/cm^3)"]
        O3(t) = 3e12, [description = "Ozone (input) (molec/cm^3)"]
        NO(t) = 1e9, [description = "NO (input) (molec/cm^3)"]
        HOx(t), [description = "HOx = OH + HO2 (molec/cm^3)"]
    end

    eqs = [
        # OH rate equation
        D(OH) ~ -k_OH_O3 * OH * O3 + k_HO2_O * HO2 * O + k_HO2_O3 * HO2 * O3 + k_HO2_NO * HO2 * NO,

        # HO2 rate equation
        D(HO2) ~ k_OH_O3 * OH * O3 - k_HO2_O * HO2 * O - k_HO2_O3 * HO2 * O3 - k_HO2_NO * HO2 * NO,

        # HOx family
        HOx ~ OH + HO2
    ]

    return ODESystem(eqs, t; name=name)
end

"""
    ClOxCycle(; name=:clox)

Create a ModelingToolkit ODESystem for the ClOx catalytic ozone destruction cycle.

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
function ClOxCycle(; name=:clox)
    @parameters begin
        t
        k_Cl_O3 = 2.3e-11, [description = "Cl + O3 rate (cm^3/molec/s)"]
        k_ClO_O = 3e-11, [description = "ClO + O rate (cm^3/molec/s)"]
        k_ClO_NO = 6.4e-12, [description = "ClO + NO rate (cm^3/molec/s)"]
        k_Cl_CH4 = 1e-14, [description = "Cl + CH4 rate (cm^3/molec/s)"]
        k_OH_HCl = 2.6e-12, [description = "OH + HCl rate (cm^3/molec/s)"]
        j_ClONO2 = 1e-4, [description = "ClONO2 photolysis rate (s^-1)"]
        CH4 = 1e13, [description = "Methane concentration (molec/cm^3)"]
    end

    @variables begin
        Cl(t) = 1e4, [description = "Cl atom concentration (molec/cm^3)"]
        ClO(t) = 1e7, [description = "ClO concentration (molec/cm^3)"]
        HCl(t) = 1e9, [description = "HCl reservoir (molec/cm^3)"]
        ClONO2(t) = 1e9, [description = "ClONO2 reservoir (molec/cm^3)"]
        O(t) = 1e7, [description = "Atomic oxygen (input) (molec/cm^3)"]
        O3(t) = 3e12, [description = "Ozone (input) (molec/cm^3)"]
        NO(t) = 1e9, [description = "NO (input) (molec/cm^3)"]
        OH(t) = 1e6, [description = "OH (input) (molec/cm^3)"]
        ClOx(t), [description = "ClOx = Cl + ClO (molec/cm^3)"]
        Cly(t), [description = "Cly = Cl + ClO + HCl + ClONO2 (molec/cm^3)"]
    end

    eqs = [
        # Cl rate equation
        D(Cl) ~ -k_Cl_O3 * Cl * O3 + k_ClO_O * ClO * O + k_ClO_NO * ClO * NO -
                k_Cl_CH4 * Cl * CH4 + k_OH_HCl * OH * HCl + j_ClONO2 * ClONO2,

        # ClO rate equation
        D(ClO) ~ k_Cl_O3 * Cl * O3 - k_ClO_O * ClO * O - k_ClO_NO * ClO * NO,

        # HCl reservoir
        D(HCl) ~ k_Cl_CH4 * Cl * CH4 - k_OH_HCl * OH * HCl,

        # ClONO2 reservoir (simplified - full treatment requires NO2)
        D(ClONO2) ~ -j_ClONO2 * ClONO2,

        # ClOx family
        ClOx ~ Cl + ClO,

        # Cly family (total inorganic chlorine)
        Cly ~ Cl + ClO + HCl + ClONO2
    ]

    return ODESystem(eqs, t; name=name)
end

"""
    BrOxCycle(; name=:brox)

Create a ModelingToolkit ODESystem for the BrOx catalytic ozone destruction cycle.

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
function BrOxCycle(; name=:brox)
    @parameters begin
        t
        k_Br_O3 = 7e-13, [description = "Br + O3 rate (cm^3/molec/s)"]
        k_BrO_O = 5e-11, [description = "BrO + O rate (cm^3/molec/s)"]
        k_BrO_ClO = 2e-12, [description = "BrO + ClO rate (cm^3/molec/s)"]
        k_BrO_HO2 = 4e-11, [description = "BrO + HO2 rate (cm^3/molec/s)"]
        j_HOBr = 1e-3, [description = "HOBr photolysis rate (s^-1)"]
    end

    @variables begin
        Br(t) = 1e5, [description = "Br atom concentration (molec/cm^3)"]
        BrO(t) = 1e6, [description = "BrO concentration (molec/cm^3)"]
        HOBr(t) = 1e6, [description = "HOBr concentration (molec/cm^3)"]
        O(t) = 1e7, [description = "Atomic oxygen (input) (molec/cm^3)"]
        O3(t) = 3e12, [description = "Ozone (input) (molec/cm^3)"]
        ClO(t) = 1e7, [description = "ClO (input) (molec/cm^3)"]
        HO2(t) = 1e7, [description = "HO2 (input) (molec/cm^3)"]
        BrOx(t), [description = "BrOx = Br + BrO (molec/cm^3)"]
        Bry(t), [description = "Bry = Br + BrO + HOBr (molec/cm^3)"]
    end

    eqs = [
        # Br rate equation
        D(Br) ~ -k_Br_O3 * Br * O3 + k_BrO_O * BrO * O + k_BrO_ClO * BrO * ClO + j_HOBr * HOBr,

        # BrO rate equation
        D(BrO) ~ k_Br_O3 * Br * O3 - k_BrO_O * BrO * O - k_BrO_ClO * BrO * ClO - k_BrO_HO2 * BrO * HO2,

        # HOBr
        D(HOBr) ~ k_BrO_HO2 * BrO * HO2 - j_HOBr * HOBr,

        # BrOx family
        BrOx ~ Br + BrO,

        # Bry family
        Bry ~ Br + BrO + HOBr
    ]

    return ODESystem(eqs, t; name=name)
end

"""
    StratosphericOzoneSystem(; name=:strat_ozone)

Create a comprehensive ModelingToolkit ODESystem combining all stratospheric
ozone chemistry cycles.

This system includes:
- Chapman mechanism (O, O3 production and loss)
- NOx cycle (catalytic O3 destruction)
- HOx cycle (catalytic O3 destruction)
- ClOx cycle (catalytic O3 destruction)
- BrOx cycle (catalytic O3 destruction)

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
function StratosphericOzoneSystem(; name=:strat_ozone)
    @parameters begin
        t
        # Photolysis rates
        j_O2 = 1e-10, [description = "O2 photolysis rate (s^-1)"]
        j_O3 = 1e-3, [description = "O3 photolysis rate (total) (s^-1)"]
        j_O3_O1D = 5e-4, [description = "O3 -> O(1D) photolysis rate (s^-1)"]
        j_NO2 = 1e-2, [description = "NO2 photolysis rate (s^-1)"]
        j_ClONO2 = 1e-4, [description = "ClONO2 photolysis rate (s^-1)"]
        j_HOBr = 1e-3, [description = "HOBr photolysis rate (s^-1)"]

        # Environmental parameters
        T = 227.0, [description = "Temperature (K)"]
        M = 3.1e17, [description = "Air number density at 30 km (molec/cm^3)"]
        O2_mix = 0.21, [description = "O2 mixing ratio"]
        N2_mix = 0.79, [description = "N2 mixing ratio"]
        CH4_mix = 1.6e-6, [description = "CH4 mixing ratio"]
        H2O_mix = 5e-6, [description = "H2O mixing ratio"]
    end

    # Derived concentrations
    O2_conc = O2_mix * M
    N2_conc = N2_mix * M
    CH4_conc = CH4_mix * M
    H2O_conc = H2O_mix * M

    @variables begin
        # Odd oxygen family
        O(t) = 1e7, [description = "Atomic oxygen O(3P) (molec/cm^3)"]
        O1D(t) = 50.0, [description = "Excited oxygen O(1D) (molec/cm^3)"]
        O3(t) = 3e12, [description = "Ozone (molec/cm^3)"]

        # Nitrogen oxides
        NO(t) = 1e9, [description = "Nitric oxide (molec/cm^3)"]
        NO2(t) = 1e9, [description = "Nitrogen dioxide (molec/cm^3)"]
        NO3(t) = 1e6, [description = "Nitrate radical (molec/cm^3)"]
        N2O5(t) = 1e8, [description = "Dinitrogen pentoxide (molec/cm^3)"]
        HNO3(t) = 1e10, [description = "Nitric acid (molec/cm^3)"]

        # Hydrogen radicals
        OH(t) = 1e6, [description = "Hydroxyl radical (molec/cm^3)"]
        HO2(t) = 1e7, [description = "Hydroperoxyl radical (molec/cm^3)"]

        # Chlorine species
        Cl(t) = 1e4, [description = "Chlorine atom (molec/cm^3)"]
        ClO(t) = 1e7, [description = "Chlorine monoxide (molec/cm^3)"]
        HCl(t) = 1e9, [description = "Hydrogen chloride (molec/cm^3)"]
        ClONO2(t) = 1e9, [description = "Chlorine nitrate (molec/cm^3)"]

        # Bromine species
        Br(t) = 1e5, [description = "Bromine atom (molec/cm^3)"]
        BrO(t) = 1e6, [description = "Bromine monoxide (molec/cm^3)"]
        HOBr(t) = 1e6, [description = "Hypobromous acid (molec/cm^3)"]

        # Chemical families (algebraic)
        Ox(t), [description = "Odd oxygen = O + O3 (molec/cm^3)"]
        NOx(t), [description = "NOx = NO + NO2 (molec/cm^3)"]
        HOx(t), [description = "HOx = OH + HO2 (molec/cm^3)"]
        ClOx(t), [description = "ClOx = Cl + ClO (molec/cm^3)"]
        BrOx(t), [description = "BrOx = Br + BrO (molec/cm^3)"]
    end

    # Rate coefficients (temperature-dependent)
    k2 = 6.0e-34 * (T / 300.0)^(-2.4)  # O + O2 + M
    k4 = 8.0e-12 * exp(-2060.0 / T)     # O + O3
    k_O1D_M = 0.21 * 3.2e-11 * exp(70.0 / T) + 0.79 * 1.8e-11 * exp(110.0 / T)  # O(1D) quenching
    k_O1D_H2O = 2.2e-10  # O(1D) + H2O → 2 OH
    k_NO_O3 = 3.0e-12 * exp(-1500.0 / T)
    k_NO2_O = 5.6e-12 * exp(180.0 / T)
    k_OH_O3 = 1.7e-12 * exp(-940.0 / T)
    k_HO2_O3 = 1.0e-15
    k_HO2_O = 3.0e-11 * exp(200.0 / T)
    k_HO2_NO = 3.5e-12 * exp(250.0 / T)
    k_Cl_O3 = 2.3e-11 * exp(-200.0 / T)
    k_ClO_O = 3.0e-11 * exp(70.0 / T)
    k_ClO_NO = 6.4e-12 * exp(290.0 / T)
    k_Cl_CH4 = 1.0e-14
    k_OH_HCl = 2.6e-12 * exp(-350.0 / T)
    k_Br_O3 = 7.0e-13
    k_BrO_O = 5.0e-11
    k_BrO_ClO = 2.0e-12
    k_BrO_HO2 = 4.0e-11

    eqs = [
        # =====================================================================
        # Odd Oxygen Species
        # =====================================================================

        # Atomic oxygen O(3P) - Equations 5.1
        D(O) ~ 2 * j_O2 * O2_conc +              # O2 photolysis (source)
               j_O3 * O3 +                        # O3 photolysis (source)
               k_O1D_M * O1D * M -                # O(1D) quenching (source)
               k2 * O * O2_conc * M -             # O + O2 + M → O3 (sink)
               k4 * O * O3 -                      # O + O3 → 2 O2 (sink)
               k_NO2_O * NO2 * O -                # NO2 + O → NO + O2 (sink)
               k_ClO_O * ClO * O -                # ClO + O → Cl + O2 (sink)
               k_BrO_O * BrO * O -                # BrO + O → Br + O2 (sink)
               k_HO2_O * HO2 * O,                 # HO2 + O → OH + O2 (sink)

        # Excited oxygen O(1D)
        D(O1D) ~ j_O3_O1D * O3 -                  # O3 + hν → O(1D) + O2
                 k_O1D_M * O1D * M -              # Quenching to O(3P)
                 k_O1D_H2O * O1D * H2O_conc,      # Reaction with H2O

        # Ozone - Equation 5.2
        D(O3) ~ k2 * O * O2_conc * M -           # O + O2 + M → O3 (source)
                j_O3 * O3 -                       # O3 photolysis (sink)
                k4 * O * O3 -                     # O + O3 → 2 O2 (sink)
                k_NO_O3 * NO * O3 -               # NO + O3 → NO2 + O2 (sink)
                k_OH_O3 * OH * O3 -               # OH + O3 → HO2 + O2 (sink)
                k_HO2_O3 * HO2 * O3 -             # HO2 + O3 → OH + 2 O2 (sink)
                k_Cl_O3 * Cl * O3 -               # Cl + O3 → ClO + O2 (sink)
                k_Br_O3 * Br * O3,                # Br + O3 → BrO + O2 (sink)

        # =====================================================================
        # Nitrogen Oxide Species
        # =====================================================================

        # NO
        D(NO) ~ j_NO2 * NO2 +                    # NO2 photolysis
                k_NO2_O * NO2 * O -               # NO2 + O → NO + O2
                k_NO_O3 * NO * O3 +               # NO + O3 → NO2 + O2
                k_ClO_NO * ClO * NO +             # ClO + NO → Cl + NO2 (regenerates NO indirectly)
                k_HO2_NO * HO2 * NO,              # HO2 + NO → NO2 + OH (actually consumes NO)

        # NO2
        D(NO2) ~ k_NO_O3 * NO * O3 +             # NO + O3 → NO2 + O2
                 k_HO2_NO * HO2 * NO -            # HO2 + NO → NO2 + OH
                 j_NO2 * NO2 -                    # NO2 photolysis
                 k_NO2_O * NO2 * O,               # NO2 + O → NO + O2

        # =====================================================================
        # HOx Species
        # =====================================================================

        # OH - production from O(1D) + H2O
        D(OH) ~ 2 * k_O1D_H2O * O1D * H2O_conc + # O(1D) + H2O → 2 OH (source)
                k_HO2_O * HO2 * O +               # HO2 + O → OH + O2
                k_HO2_O3 * HO2 * O3 +             # HO2 + O3 → OH + 2 O2
                k_HO2_NO * HO2 * NO -             # HO2 + NO → NO2 + OH
                k_OH_O3 * OH * O3 -               # OH + O3 → HO2 + O2
                k_OH_HCl * OH * HCl,              # OH + HCl → H2O + Cl

        # HO2
        D(HO2) ~ k_OH_O3 * OH * O3 -             # OH + O3 → HO2 + O2
                 k_HO2_O * HO2 * O -              # HO2 + O → OH + O2
                 k_HO2_O3 * HO2 * O3 -            # HO2 + O3 → OH + 2 O2
                 k_HO2_NO * HO2 * NO -            # HO2 + NO → NO2 + OH
                 k_BrO_HO2 * BrO * HO2,           # BrO + HO2 → HOBr + O2

        # =====================================================================
        # Chlorine Species
        # =====================================================================

        # Cl
        D(Cl) ~ k_ClO_O * ClO * O +              # ClO + O → Cl + O2
                k_ClO_NO * ClO * NO +             # ClO + NO → Cl + NO2
                k_OH_HCl * OH * HCl +             # OH + HCl → H2O + Cl
                j_ClONO2 * ClONO2 +               # ClONO2 + hν → Cl + NO3
                k_BrO_ClO * BrO * ClO -           # BrO + ClO → Br + Cl + O2
                k_Cl_O3 * Cl * O3 -               # Cl + O3 → ClO + O2
                k_Cl_CH4 * Cl * CH4_conc,         # Cl + CH4 → HCl + CH3

        # ClO
        D(ClO) ~ k_Cl_O3 * Cl * O3 -             # Cl + O3 → ClO + O2
                 k_ClO_O * ClO * O -              # ClO + O → Cl + O2
                 k_ClO_NO * ClO * NO -            # ClO + NO → Cl + NO2
                 k_BrO_ClO * BrO * ClO,           # BrO + ClO → products

        # HCl reservoir
        D(HCl) ~ k_Cl_CH4 * Cl * CH4_conc -      # Cl + CH4 → HCl + CH3
                 k_OH_HCl * OH * HCl,             # OH + HCl → H2O + Cl

        # ClONO2 reservoir (simplified)
        D(ClONO2) ~ -j_ClONO2 * ClONO2,          # ClONO2 photolysis

        # =====================================================================
        # Bromine Species
        # =====================================================================

        # Br
        D(Br) ~ k_BrO_O * BrO * O +              # BrO + O → Br + O2
                k_BrO_ClO * BrO * ClO +           # BrO + ClO → Br + Cl + O2
                j_HOBr * HOBr -                   # HOBr + hν → OH + Br
                k_Br_O3 * Br * O3,                # Br + O3 → BrO + O2

        # BrO
        D(BrO) ~ k_Br_O3 * Br * O3 -             # Br + O3 → BrO + O2
                 k_BrO_O * BrO * O -              # BrO + O → Br + O2
                 k_BrO_ClO * BrO * ClO -          # BrO + ClO → products
                 k_BrO_HO2 * BrO * HO2,           # BrO + HO2 → HOBr + O2

        # HOBr
        D(HOBr) ~ k_BrO_HO2 * BrO * HO2 -        # BrO + HO2 → HOBr + O2
                  j_HOBr * HOBr,                  # HOBr photolysis

        # =====================================================================
        # Chemical Families (Algebraic Equations)
        # =====================================================================

        Ox ~ O + O3,
        NOx ~ NO + NO2,
        HOx ~ OH + HO2,
        ClOx ~ Cl + ClO,
        BrOx ~ Br + BrO
    ]

    return ODESystem(eqs, t; name=name)
end

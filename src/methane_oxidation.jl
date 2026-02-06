"""
    Methane Oxidation Mechanism

Implements Table 6.1 from Seinfeld & Pandis Chapter 6, Section 6.4.

Methane (CH₄) oxidation is the archetypal VOC oxidation mechanism in the
troposphere. The complete mechanism involves 17 reactions and proceeds through:

    CH₄ → CH₃O₂ → CH₃O → HCHO → HCO → CO → CO₂

Key intermediates:
- CH₃O₂: Methylperoxy radical
- CH₃O: Methoxy radical
- HCHO: Formaldehyde
- HCO: Formyl radical
- CO: Carbon monoxide

Each step can produce ozone when NOx is present through peroxy + NO reactions.

Reference: Seinfeld & Pandis (2006), Section 6.4, Table 6.1, pp. 219-227
"""

using ModelingToolkit: t, D

# ============================================================================
# Table 6.1: Complete Methane Oxidation Mechanism
# ============================================================================
# Reaction                                          | k at 298 K (cm³/molec/s)
# ============================================================================
# 1.  CH₄ + OH → CH₃ + H₂O                         | 6.3 × 10⁻¹⁵
# 2.  CH₃ + O₂ + M → CH₃O₂ + M                     | 1.0 × 10⁻³⁰ [M] (cm⁶)
# 3.  CH₃O₂ + NO → CH₃O + NO₂                      | 7.7 × 10⁻¹²
# 4.  CH₃O₂ + HO₂ → CH₃OOH + O₂                    | 5.2 × 10⁻¹²
# 5.  CH₃O₂ + CH₃O₂ → products                     | 3.5 × 10⁻¹³
# 6.  CH₃O + O₂ → HCHO + HO₂                       | 1.9 × 10⁻¹⁵
# 7.  CH₃OOH + OH → CH₃O₂ + H₂O                    | 3.8 × 10⁻¹²
# 8.  CH₃OOH + OH → HCHO + OH + H₂O                | 1.9 × 10⁻¹²
# 9.  CH₃OOH + hν → CH₃O + OH                      | j ≈ 5 × 10⁻⁶ s⁻¹
# 10. HCHO + OH → HCO + H₂O                        | 8.5 × 10⁻¹²
# 11. HCHO + hν → HCO + H                          | j ≈ 3 × 10⁻⁵ s⁻¹
# 12. HCHO + hν → H₂ + CO                          | j ≈ 5 × 10⁻⁵ s⁻¹
# 13. HCO + O₂ → CO + HO₂                          | 5.2 × 10⁻¹²
# 14. H + O₂ + M → HO₂ + M                         | 5.7 × 10⁻³² [M] (cm⁶)
# 15. HO₂ + NO → OH + NO₂                          | 8.1 × 10⁻¹²
# 16. NO₂ + hν → NO + O                            | j ≈ 8 × 10⁻³ s⁻¹
# 17. O + O₂ + M → O₃ + M                          | 6.0 × 10⁻³⁴ [M] (cm⁶)

"""
    MethaneOxidation(; name)

ModelingToolkit System implementing reaction rate diagnostics for the methane
oxidation mechanism from Table 6.1 of Seinfeld & Pandis Chapter 6.

This system computes individual reaction rates and diagnostic production/loss terms.

# Species (Input Variables)

  - CH4, CH3, CH3O2, CH3O, CH3OOH: Methane chain species [m⁻³]
  - HCHO, HCO, CO: Formaldehyde and products [m⁻³]
  - OH, HO2, H: HOx species [m⁻³]
  - NO, NO2, O, O2, O3: NOx and oxygen species [m⁻³]
  - M: Total air density [m⁻³]

# Diagnostics (Output Variables)

  - R1-R17: Individual reaction rates [m⁻³ s⁻¹]
  - P_O3_net: Net O₃ production [m⁻³ s⁻¹]
  - P_HCHO: HCHO production rate [m⁻³ s⁻¹]
  - L_CH4: CH₄ loss rate [m⁻³ s⁻¹]

# Rate Constants

All rate constants from Table 6.1 at 298 K are implemented as parameters.
Bimolecular rate constants converted from cm³/molec/s to m³/s (×10⁻⁶).
Termolecular rate constants converted from cm⁶/molec²/s to m⁶/s (×10⁻¹²).    # Parameters - Rate constants at 298 K from Table 6.1 (converted to SI)
"""
@component function MethaneOxidation(; name = :MethaneOxidation)
    # Parameters - Rate constants at 298 K from Table 6.1 (converted to SI)
    @parameters begin
        # Bimolecular reactions (m³ s⁻¹)
        k1 = 6.3e-15 * 1e-6,
        [description = "CH₄ + OH rate (6.3e-15 cm³/molec/s)", unit = u"m^3/s"]
        k3 = 7.7e-12 * 1e-6,
        [description = "CH₃O₂ + NO rate (7.7e-12 cm³/molec/s)", unit = u"m^3/s"]
        k4 = 5.2e-12 * 1e-6,
        [description = "CH₃O₂ + HO₂ rate (5.2e-12 cm³/molec/s)", unit = u"m^3/s"]
        k5 = 3.5e-13 * 1e-6,
        [description = "CH₃O₂ + CH₃O₂ rate (3.5e-13 cm³/molec/s)", unit = u"m^3/s"]
        k6 = 1.9e-15 * 1e-6,
        [description = "CH₃O + O₂ rate (1.9e-15 cm³/molec/s)", unit = u"m^3/s"]
        k7 = 3.8e-12 * 1e-6,
        [description = "CH₃OOH + OH → CH₃O₂ rate (3.8e-12 cm³/molec/s)", unit = u"m^3/s"]
        k8 = 1.9e-12 * 1e-6,
        [description = "CH₃OOH + OH → HCHO rate (1.9e-12 cm³/molec/s)", unit = u"m^3/s"]
        k10 = 8.5e-12 * 1e-6,
        [description = "HCHO + OH rate (8.5e-12 cm³/molec/s)", unit = u"m^3/s"]
        k13 = 5.2e-12 * 1e-6,
        [description = "HCO + O₂ rate (5.2e-12 cm³/molec/s)", unit = u"m^3/s"]
        k15 = 8.1e-12 * 1e-6,
        [description = "HO₂ + NO rate (8.1e-12 cm³/molec/s)", unit = u"m^3/s"]

        # Termolecular reactions (m⁶ s⁻¹)
        k2_0 = 1.0e-30 * 1e-12,
        [description = "CH₃ + O₂ + M rate (1.0e-30 cm⁶/molec²/s)", unit = u"m^6/s"]
        k14_0 = 5.7e-32 * 1e-12,
        [description = "H + O₂ + M rate (5.7e-32 cm⁶/molec²/s)", unit = u"m^6/s"]
        k17_0 = 6.0e-34 * 1e-12,
        [description = "O + O₂ + M rate (6.0e-34 cm⁶/molec²/s)", unit = u"m^6/s"]

        # Photolysis rates (s⁻¹)
        j9 = 5e-6, [description = "CH₃OOH photolysis rate", unit = u"s^-1"]
        j11 = 3e-5, [description = "HCHO → HCO + H photolysis rate", unit = u"s^-1"]
        j12 = 5e-5, [description = "HCHO → H₂ + CO photolysis rate", unit = u"s^-1"]
        j16 = 8e-3, [description = "NO₂ photolysis rate", unit = u"s^-1"]
    end

    # Species concentrations (input, in SI: m⁻³)
    @variables begin
        CH4(t), [description = "Methane", unit = u"m^-3"]
        CH3(t), [description = "Methyl radical", unit = u"m^-3"]
        CH3O2(t), [description = "Methylperoxy radical", unit = u"m^-3"]
        CH3O(t), [description = "Methoxy radical", unit = u"m^-3"]
        CH3OOH(t), [description = "Methyl hydroperoxide", unit = u"m^-3"]
        HCHO(t), [description = "Formaldehyde", unit = u"m^-3"]
        HCO(t), [description = "Formyl radical", unit = u"m^-3"]
        CO(t), [description = "Carbon monoxide", unit = u"m^-3"]
        OH(t), [description = "Hydroxyl radical", unit = u"m^-3"]
        HO2(t), [description = "Hydroperoxy radical", unit = u"m^-3"]
        H(t), [description = "Hydrogen atom", unit = u"m^-3"]
        NO(t), [description = "Nitric oxide", unit = u"m^-3"]
        NO2(t), [description = "Nitrogen dioxide", unit = u"m^-3"]
        O(t), [description = "Oxygen atom", unit = u"m^-3"]
        O2(t), [description = "Molecular oxygen", unit = u"m^-3"]
        O3(t), [description = "Ozone", unit = u"m^-3"]
        M(t), [description = "Total air density", unit = u"m^-3"]
    end

    # Reaction rates (output, in SI: m⁻³ s⁻¹)
    @variables begin
        R1(t), [description = "CH₄ + OH rate", unit = u"m^-3*s^-1"]
        R2(t), [description = "CH₃ + O₂ rate", unit = u"m^-3*s^-1"]
        R3(t), [description = "CH₃O₂ + NO rate", unit = u"m^-3*s^-1"]
        R4(t), [description = "CH₃O₂ + HO₂ rate", unit = u"m^-3*s^-1"]
        R5(t), [description = "CH₃O₂ + CH₃O₂ rate", unit = u"m^-3*s^-1"]
        R6(t), [description = "CH₃O + O₂ rate", unit = u"m^-3*s^-1"]
        R7(t), [description = "CH₃OOH + OH → CH₃O₂ rate", unit = u"m^-3*s^-1"]
        R8(t), [description = "CH₃OOH + OH → HCHO rate", unit = u"m^-3*s^-1"]
        R9(t), [description = "CH₃OOH photolysis rate", unit = u"m^-3*s^-1"]
        R10(t), [description = "HCHO + OH rate", unit = u"m^-3*s^-1"]
        R11(t), [description = "HCHO → HCO + H rate", unit = u"m^-3*s^-1"]
        R12(t), [description = "HCHO → H₂ + CO rate", unit = u"m^-3*s^-1"]
        R13(t), [description = "HCO + O₂ rate", unit = u"m^-3*s^-1"]
        R14(t), [description = "H + O₂ rate", unit = u"m^-3*s^-1"]
        R15(t), [description = "HO₂ + NO rate", unit = u"m^-3*s^-1"]
        R16(t), [description = "NO₂ photolysis rate", unit = u"m^-3*s^-1"]
        R17(t), [description = "O + O₂ rate", unit = u"m^-3*s^-1"]
        P_O3_net(t), [description = "Net O₃ production", unit = u"m^-3*s^-1"]
        P_HCHO(t), [description = "HCHO production", unit = u"m^-3*s^-1"]
        L_CH4(t), [description = "CH₄ loss rate", unit = u"m^-3*s^-1"]
    end

    # Equations
    eqs = [
        # Individual reaction rates (Table 6.1)
        R1 ~ k1 * CH4 * OH,                    # CH₄ + OH → CH₃ + H₂O
        R2 ~ k2_0 * M * CH3 * O2,              # CH₃ + O₂ + M → CH₃O₂ + M
        R3 ~ k3 * CH3O2 * NO,                  # CH₃O₂ + NO → CH₃O + NO₂
        R4 ~ k4 * CH3O2 * HO2,                 # CH₃O₂ + HO₂ → CH₃OOH + O₂
        R5 ~ k5 * CH3O2 * CH3O2,               # CH₃O₂ + CH₃O₂ → products
        R6 ~ k6 * CH3O * O2,                   # CH₃O + O₂ → HCHO + HO₂
        R7 ~ k7 * CH3OOH * OH,                 # CH₃OOH + OH → CH₃O₂ + H₂O
        R8 ~ k8 * CH3OOH * OH,                 # CH₃OOH + OH → HCHO + OH + H₂O
        R9 ~ j9 * CH3OOH,                      # CH₃OOH + hν → CH₃O + OH
        R10 ~ k10 * HCHO * OH,                 # HCHO + OH → HCO + H₂O
        R11 ~ j11 * HCHO,                      # HCHO + hν → HCO + H
        R12 ~ j12 * HCHO,                      # HCHO + hν → H₂ + CO
        R13 ~ k13 * HCO * O2,                  # HCO + O₂ → CO + HO₂
        R14 ~ k14_0 * M * H * O2,              # H + O₂ + M → HO₂ + M
        R15 ~ k15 * HO2 * NO,                  # HO₂ + NO → OH + NO₂
        R16 ~ j16 * NO2,                       # NO₂ + hν → NO + O
        R17 ~ k17_0 * M * O * O2,              # O + O₂ + M → O₃ + M

        # Diagnostic variables
        L_CH4 ~ R1,                            # CH₄ loss = R1
        P_HCHO ~ R6 + R8,                      # HCHO production
        P_O3_net ~ R15 + R3                   # Net O₃ production = HO₂+NO + CH₃O₂+NO (Eq. 6.9 analog)
    ]

    return System(eqs, t; name)
end

"""
    MethaneOxidationODE(; name)

Full ODE system for methane oxidation with species time derivatives.

This uses Catalyst.jl's reaction network DSL to define the 17 reactions from
Table 6.1 plus auxiliary reactions (CO+OH, OH+NO₂, HO₂+HO₂, NO+O₃) and an
external OH source. The reaction network is converted to an ODE system.

Note: This is a stiff system due to the wide range of timescales (radicals
have lifetimes of seconds, while CH₄ has a lifetime of years).

O₂ and M (total air density) are treated as parameters with default values
for termolecular reactions (reactions 2, 6, 13, 14, 17).
"""
@component function MethaneOxidationODE(; name = :MethaneOxidationODE)
    rn = @network_component MethaneOxidationRxns begin
        @ivs t [unit = u"s"]

        @parameters begin
            # Bimolecular reactions (m³ s⁻¹)
            k1 = 6.3e-15 * 1e-6,
            [description = "CH₄ + OH rate (6.3e-15 cm³/molec/s)", unit = u"m^3/s"]
            k3 = 7.7e-12 * 1e-6,
            [description = "CH₃O₂ + NO rate (7.7e-12 cm³/molec/s)", unit = u"m^3/s"]
            k4 = 5.2e-12 * 1e-6,
            [description = "CH₃O₂ + HO₂ rate (5.2e-12 cm³/molec/s)", unit = u"m^3/s"]
            k5 = 3.5e-13 * 1e-6,
            [description = "CH₃O₂ + CH₃O₂ rate (3.5e-13 cm³/molec/s)", unit = u"m^3/s"]
            k7 = 3.8e-12 * 1e-6,
            [
                description = "CH₃OOH + OH → CH₃O₂ rate (3.8e-12 cm³/molec/s)", unit = u"m^3/s"]
            k8 = 1.9e-12 * 1e-6,
            [description = "CH₃OOH + OH → HCHO rate (1.9e-12 cm³/molec/s)", unit = u"m^3/s"]
            k10 = 8.5e-12 * 1e-6,
            [description = "HCHO + OH rate (8.5e-12 cm³/molec/s)", unit = u"m^3/s"]
            k15 = 8.1e-12 * 1e-6,
            [description = "HO₂ + NO rate (8.1e-12 cm³/molec/s)", unit = u"m^3/s"]
            k_CO_OH = 2.4e-13 * 1e-6,
            [description = "CO + OH rate (2.4e-13 cm³/molec/s)", unit = u"m^3/s"]
            k_OH_NO2 = 1.0e-11 * 1e-6,
            [description = "OH + NO₂ rate (1.0e-11 cm³/molec/s)", unit = u"m^3/s"]
            k_HO2_HO2 = 2.9e-12 * 1e-6,
            [description = "HO₂ + HO₂ rate (2.9e-12 cm³/molec/s)", unit = u"m^3/s"]
            k_NO_O3 = 1.8e-14 * 1e-6,
            [description = "NO + O₃ rate (1.8e-14 cm³/molec/s)", unit = u"m^3/s"]

            # Termolecular reactions (m⁶ s⁻¹)
            k2_0 = 1.0e-30 * 1e-12,
            [description = "CH₃ + O₂ + M rate (1.0e-30 cm⁶/molec²/s)", unit = u"m^6/s"]
            k14_0 = 5.7e-32 * 1e-12,
            [description = "H + O₂ + M rate (5.7e-32 cm⁶/molec²/s)", unit = u"m^6/s"]
            k17_0 = 6.0e-34 * 1e-12,
            [description = "O + O₂ + M rate (6.0e-34 cm⁶/molec²/s)", unit = u"m^6/s"]

            # Effective bimolecular rates with O₂ folded in (m³ s⁻¹ × m⁻³ = s⁻¹)
            k6_eff = 1.9e-15 * 1e-6 * 5.25e18 * 1e6,
            [description = "CH₃O + O₂ effective rate (k6*[O₂])", unit = u"s^-1"]
            k13_eff = 5.2e-12 * 1e-6 * 5.25e18 * 1e6,
            [description = "HCO + O₂ effective rate (k13*[O₂])", unit = u"s^-1"]

            # Fixed concentrations as parameters (m⁻³)
            M_fixed = 2.5e19 * 1e6, [description = "Total air density", unit = u"m^-3"]
            O2_fixed = 5.25e18 * 1e6, [description = "O₂ concentration", unit = u"m^-3"]

            # Photolysis rates (s⁻¹)
            j9 = 5e-6, [description = "CH₃OOH photolysis rate", unit = u"s^-1"]
            j11 = 3e-5, [description = "HCHO → HCO + H rate", unit = u"s^-1"]
            j12 = 5e-5, [description = "HCHO → H₂ + CO rate", unit = u"s^-1"]
            j16 = 8e-3, [description = "NO₂ photolysis rate", unit = u"s^-1"]

            # External source for OH production (e.g., from O₃ photolysis)
            P_OH_ext = 1e6 * 1e6,
            [description = "External OH production", unit = u"m^-3*s^-1"]
        end

        @species begin
            CH4(t), [unit = u"m^-3", description = "Methane"]
            CH3(t), [unit = u"m^-3", description = "Methyl radical"]
            CH3O2(t), [unit = u"m^-3", description = "Methylperoxy radical"]
            CH3O(t), [unit = u"m^-3", description = "Methoxy radical"]
            CH3OOH(t), [unit = u"m^-3", description = "Methyl hydroperoxide"]
            HCHO(t), [unit = u"m^-3", description = "Formaldehyde"]
            HCO(t), [unit = u"m^-3", description = "Formyl radical"]
            CO(t), [unit = u"m^-3", description = "Carbon monoxide"]
            H2(t), [unit = u"m^-3", description = "Molecular hydrogen"]
            OH(t), [unit = u"m^-3", description = "Hydroxyl radical"]
            HO2(t), [unit = u"m^-3", description = "Hydroperoxy radical"]
            H(t), [unit = u"m^-3", description = "Hydrogen atom"]
            NO(t), [unit = u"m^-3", description = "Nitric oxide"]
            NO2(t), [unit = u"m^-3", description = "Nitrogen dioxide"]
            O(t), [unit = u"m^-3", description = "Oxygen atom"]
            O3(t), [unit = u"m^-3", description = "Ozone"]
            HNO3(t), [unit = u"m^-3", description = "Nitric acid"]
            H2O2(t), [unit = u"m^-3", description = "Hydrogen peroxide"]
        end

        # =====================================================================
        # Table 6.1: Methane Oxidation Mechanism (17 reactions)
        # =====================================================================

        # R1:  CH₄ + OH → CH₃ + H₂O
        k1, CH4 + OH --> CH3

        # R2:  CH₃ + O₂ + M → CH₃O₂ + M  (effective rate = k2_0 * [M] * [O₂])
        k2_0 * M_fixed * O2_fixed, CH3 --> CH3O2

        # R3:  CH₃O₂ + NO → CH₃O + NO₂
        k3, CH3O2 + NO --> CH3O + NO2

        # R4:  CH₃O₂ + HO₂ → CH₃OOH + O₂
        k4, CH3O2 + HO2 --> CH3OOH

        # R5:  CH₃O₂ + CH₃O₂ → products  (products leave system)
        k5, 2CH3O2 --> 0

        # R6:  CH₃O + O₂ → HCHO + HO₂  (effective rate = k6 * [O₂])
        k6_eff, CH3O --> HCHO + HO2

        # R7:  CH₃OOH + OH → CH₃O₂ + H₂O
        k7, CH3OOH + OH --> CH3O2

        # R8:  CH₃OOH + OH → HCHO + H₂O  (OH consumed, not regenerated per ODE)
        k8, CH3OOH + OH --> HCHO

        # R9:  CH₃OOH + hν → CH₃O + OH
        j9, CH3OOH --> CH3O + OH

        # R10: HCHO + OH → HCO + H₂O
        k10, HCHO + OH --> HCO

        # R11: HCHO + hν → HCO + H
        j11, HCHO --> HCO + H

        # R12: HCHO + hν → H₂ + CO
        j12, HCHO --> H2 + CO

        # R13: HCO + O₂ → CO + HO₂  (effective rate = k13 * [O₂])
        k13_eff, HCO --> CO + HO2

        # R14: H + O₂ + M → HO₂ + M  (effective rate = k14_0 * [M] * [O₂])
        k14_0 * M_fixed * O2_fixed, H --> HO2

        # R15: HO₂ + NO → OH + NO₂
        k15, HO2 + NO --> OH + NO2

        # R16: NO₂ + hν → NO + O
        j16, NO2 --> NO + O

        # R17: O + O₂ + M → O₃ + M  (effective rate = k17_0 * [M] * [O₂])
        k17_0 * M_fixed * O2_fixed, O --> O3

        # =====================================================================
        # Auxiliary reactions (beyond Table 6.1)
        # =====================================================================

        # R18: CO + OH → HO₂  (+ CO₂, not tracked)
        k_CO_OH, CO + OH --> HO2

        # R19: OH + NO₂ → HNO₃
        k_OH_NO2, OH + NO2 --> HNO3

        # R20: HO₂ + HO₂ → H₂O₂  (+ O₂, not tracked)
        k_HO2_HO2, 2HO2 --> H2O2

        # R21: NO + O₃ → NO₂  (+ O₂, not tracked)
        k_NO_O3, NO + O3 --> NO2

        # R22: External OH production (zero-order source, e.g., from O₃ photolysis)
        P_OH_ext, 0 --> OH
    end

    # Convert the reaction network to an ODE system.
    # combinatoric_ratelaws=false: use macroscopic rate laws (rate = k*[A]*[B])
    # rather than microscopic (rate = k*[A]*[B]/2 for A+A reactions).
    convert(Catalyst.ReactionRateSystem, complete(rn); combinatoric_ratelaws = false, name = name)
end

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

using ModelingToolkit
using Unitful
using ModelingToolkit: t_nounits as t, D_nounits as D

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
- CH4, CH3, CH3O2, CH3O, CH3OOH: Methane chain species
- HCHO, HCO, CO: Formaldehyde and products
- OH, HO2, H: HOx species
- NO, NO2, O, O2, O3: NOx and oxygen species
- M: Total air density [molecules/cm³]

# Diagnostics (Output Variables)
- R1-R17: Individual reaction rates [molecules/cm³/s]
- P_O3_net: Net O₃ production [molecules/cm³/s]
- P_HCHO: HCHO production rate [molecules/cm³/s]
- L_CH4: CH₄ loss rate [molecules/cm³/s]

# Rate Constants
All rate constants from Table 6.1 at 298 K are implemented as parameters.
"""
function MethaneOxidation(; name=:MethaneOxidation)
    # Parameters - Rate constants at 298 K from Table 6.1
    @parameters begin
        # Bimolecular reactions (cm³ molecule⁻¹ s⁻¹)
        k1 = 6.3e-15, [description = "CH₄ + OH rate [cm³/s]"]
        k3 = 7.7e-12, [description = "CH₃O₂ + NO rate [cm³/s]"]
        k4 = 5.2e-12, [description = "CH₃O₂ + HO₂ rate [cm³/s]"]
        k5 = 3.5e-13, [description = "CH₃O₂ + CH₃O₂ rate [cm³/s]"]
        k6 = 1.9e-15, [description = "CH₃O + O₂ rate [cm³/s]"]
        k7 = 3.8e-12, [description = "CH₃OOH + OH → CH₃O₂ rate [cm³/s]"]
        k8 = 1.9e-12, [description = "CH₃OOH + OH → HCHO rate [cm³/s]"]
        k10 = 8.5e-12, [description = "HCHO + OH rate [cm³/s]"]
        k13 = 5.2e-12, [description = "HCO + O₂ rate [cm³/s]"]
        k15 = 8.1e-12, [description = "HO₂ + NO rate [cm³/s]"]

        # Termolecular reactions (cm⁶ molecule⁻² s⁻¹, multiply by [M])
        k2_0 = 1.0e-30, [description = "CH₃ + O₂ + M rate [cm⁶/s]"]
        k14_0 = 5.7e-32, [description = "H + O₂ + M rate [cm⁶/s]"]
        k17_0 = 6.0e-34, [description = "O + O₂ + M rate [cm⁶/s]"]

        # Photolysis rates (s⁻¹)
        j9 = 5e-6, [description = "CH₃OOH photolysis rate [s⁻¹]"]
        j11 = 3e-5, [description = "HCHO → HCO + H photolysis rate [s⁻¹]"]
        j12 = 5e-5, [description = "HCHO → H₂ + CO photolysis rate [s⁻¹]"]
        j16 = 8e-3, [description = "NO₂ photolysis rate [s⁻¹]"]
    end

    # Species concentrations (input)
    @variables begin
        CH4(t), [description = "Methane [molecules/cm³]"]
        CH3(t), [description = "Methyl radical [molecules/cm³]"]
        CH3O2(t), [description = "Methylperoxy radical [molecules/cm³]"]
        CH3O(t), [description = "Methoxy radical [molecules/cm³]"]
        CH3OOH(t), [description = "Methyl hydroperoxide [molecules/cm³]"]
        HCHO(t), [description = "Formaldehyde [molecules/cm³]"]
        HCO(t), [description = "Formyl radical [molecules/cm³]"]
        CO(t), [description = "Carbon monoxide [molecules/cm³]"]
        OH(t), [description = "Hydroxyl radical [molecules/cm³]"]
        HO2(t), [description = "Hydroperoxy radical [molecules/cm³]"]
        H(t), [description = "Hydrogen atom [molecules/cm³]"]
        NO(t), [description = "Nitric oxide [molecules/cm³]"]
        NO2(t), [description = "Nitrogen dioxide [molecules/cm³]"]
        O(t), [description = "Oxygen atom [molecules/cm³]"]
        O2(t), [description = "Molecular oxygen [molecules/cm³]"]
        O3(t), [description = "Ozone [molecules/cm³]"]
        M(t), [description = "Total air density [molecules/cm³]"]
    end

    # Reaction rates (output)
    @variables begin
        R1(t), [description = "CH₄ + OH rate [molecules/cm³/s]"]
        R2(t), [description = "CH₃ + O₂ rate [molecules/cm³/s]"]
        R3(t), [description = "CH₃O₂ + NO rate [molecules/cm³/s]"]
        R4(t), [description = "CH₃O₂ + HO₂ rate [molecules/cm³/s]"]
        R5(t), [description = "CH₃O₂ + CH₃O₂ rate [molecules/cm³/s]"]
        R6(t), [description = "CH₃O + O₂ rate [molecules/cm³/s]"]
        R7(t), [description = "CH₃OOH + OH → CH₃O₂ rate [molecules/cm³/s]"]
        R8(t), [description = "CH₃OOH + OH → HCHO rate [molecules/cm³/s]"]
        R9(t), [description = "CH₃OOH photolysis rate [molecules/cm³/s]"]
        R10(t), [description = "HCHO + OH rate [molecules/cm³/s]"]
        R11(t), [description = "HCHO → HCO + H rate [molecules/cm³/s]"]
        R12(t), [description = "HCHO → H₂ + CO rate [molecules/cm³/s]"]
        R13(t), [description = "HCO + O₂ rate [molecules/cm³/s]"]
        R14(t), [description = "H + O₂ rate [molecules/cm³/s]"]
        R15(t), [description = "HO₂ + NO rate [molecules/cm³/s]"]
        R16(t), [description = "NO₂ photolysis rate [molecules/cm³/s]"]
        R17(t), [description = "O + O₂ rate [molecules/cm³/s]"]
        P_O3_net(t), [description = "Net O₃ production [molecules/cm³/s]"]
        P_HCHO(t), [description = "HCHO production [molecules/cm³/s]"]
        L_CH4(t), [description = "CH₄ loss rate [molecules/cm³/s]"]
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
        P_O3_net ~ R17 - R3,                   # Simplified net O₃ (production - titration)
    ]

    return ODESystem(eqs, t; name=name)
end

"""
    MethaneOxidationODE(; name)

Full ODE system for methane oxidation with species time derivatives.

This extends MethaneOxidation to include the actual differential equations
for time evolution of all species. This can be used to simulate the time
evolution of the methane oxidation mechanism.

Note: This is a stiff system due to the wide range of timescales (radicals
have lifetimes of seconds, while CH₄ has a lifetime of years).
"""
function MethaneOxidationODE(; name=:MethaneOxidationODE)
    # Parameters - Rate constants at 298 K from Table 6.1
    @parameters begin
        k1 = 6.3e-15, [description = "CH₄ + OH rate [cm³/s]"]
        k3 = 7.7e-12, [description = "CH₃O₂ + NO rate [cm³/s]"]
        k4 = 5.2e-12, [description = "CH₃O₂ + HO₂ rate [cm³/s]"]
        k5 = 3.5e-13, [description = "CH₃O₂ + CH₃O₂ rate [cm³/s]"]
        k6 = 1.9e-15, [description = "CH₃O + O₂ rate [cm³/s]"]
        k7 = 3.8e-12, [description = "CH₃OOH + OH → CH₃O₂ rate [cm³/s]"]
        k8 = 1.9e-12, [description = "CH₃OOH + OH → HCHO rate [cm³/s]"]
        k10 = 8.5e-12, [description = "HCHO + OH rate [cm³/s]"]
        k13 = 5.2e-12, [description = "HCO + O₂ rate [cm³/s]"]
        k15 = 8.1e-12, [description = "HO₂ + NO rate [cm³/s]"]
        k_CO_OH = 2.4e-13, [description = "CO + OH rate [cm³/s]"]
        k_OH_NO2 = 1.0e-11, [description = "OH + NO₂ rate [cm³/s]"]
        k_HO2_HO2 = 2.9e-12, [description = "HO₂ + HO₂ rate [cm³/s]"]
        k_NO_O3 = 1.8e-14, [description = "NO + O₃ rate [cm³/s]"]

        k2_0 = 1.0e-30, [description = "CH₃ + O₂ + M rate [cm⁶/s]"]
        k14_0 = 5.7e-32, [description = "H + O₂ + M rate [cm⁶/s]"]
        k17_0 = 6.0e-34, [description = "O + O₂ + M rate [cm⁶/s]"]

        j9 = 5e-6, [description = "CH₃OOH photolysis rate [s⁻¹]"]
        j11 = 3e-5, [description = "HCHO → HCO + H rate [s⁻¹]"]
        j12 = 5e-5, [description = "HCHO → H₂ + CO rate [s⁻¹]"]
        j16 = 8e-3, [description = "NO₂ photolysis rate [s⁻¹]"]

        # External source for OH production (e.g., from O₃ photolysis)
        P_OH_ext = 1e6, [description = "External OH production [molecules/cm³/s]"]

        # Fixed concentrations (effectively constant on CH4 oxidation timescale)
        M_fixed = 2.5e19, [description = "Total air density [molecules/cm³]"]
        O2_fixed = 5.25e18, [description = "O₂ concentration [molecules/cm³]"]
    end

    # State variables (species that evolve in time)
    @variables begin
        CH4(t), [description = "Methane [molecules/cm³]"]
        CH3(t), [description = "Methyl radical [molecules/cm³]"]
        CH3O2(t), [description = "Methylperoxy radical [molecules/cm³]"]
        CH3O(t), [description = "Methoxy radical [molecules/cm³]"]
        CH3OOH(t), [description = "Methyl hydroperoxide [molecules/cm³]"]
        HCHO(t), [description = "Formaldehyde [molecules/cm³]"]
        HCO(t), [description = "Formyl radical [molecules/cm³]"]
        CO(t), [description = "Carbon monoxide [molecules/cm³]"]
        H2(t), [description = "Molecular hydrogen [molecules/cm³]"]
        OH(t), [description = "Hydroxyl radical [molecules/cm³]"]
        HO2(t), [description = "Hydroperoxy radical [molecules/cm³]"]
        H(t), [description = "Hydrogen atom [molecules/cm³]"]
        NO(t), [description = "Nitric oxide [molecules/cm³]"]
        NO2(t), [description = "Nitrogen dioxide [molecules/cm³]"]
        O(t), [description = "Oxygen atom [molecules/cm³]"]
        O3(t), [description = "Ozone [molecules/cm³]"]
        HNO3(t), [description = "Nitric acid [molecules/cm³]"]
        H2O2(t), [description = "Hydrogen peroxide [molecules/cm³]"]
    end

    # Build the ODE system
    eqs = [
        # d[CH₄]/dt = -k1[CH₄][OH]
        D(CH4) ~ -k1 * CH4 * OH,

        # d[CH₃]/dt = k1[CH₄][OH] - k2[CH₃][O₂][M]
        D(CH3) ~ k1 * CH4 * OH - k2_0 * M_fixed * CH3 * O2_fixed,

        # d[CH₃O₂]/dt = k2[CH₃][O₂] - k3[CH₃O₂][NO] - k4[CH₃O₂][HO₂]
        #               - 2k5[CH₃O₂]² + k7[CH₃OOH][OH]
        D(CH3O2) ~ k2_0 * M_fixed * CH3 * O2_fixed - k3 * CH3O2 * NO - k4 * CH3O2 * HO2 - 2 * k5 * CH3O2^2 + k7 * CH3OOH * OH,

        # d[CH₃O]/dt = k3[CH₃O₂][NO] + j9[CH₃OOH] - k6[CH₃O][O₂]
        D(CH3O) ~ k3 * CH3O2 * NO + j9 * CH3OOH - k6 * CH3O * O2_fixed,

        # d[CH₃OOH]/dt = k4[CH₃O₂][HO₂] - k7[CH₃OOH][OH] - k8[CH₃OOH][OH] - j9[CH₃OOH]
        D(CH3OOH) ~ k4 * CH3O2 * HO2 - k7 * CH3OOH * OH - k8 * CH3OOH * OH - j9 * CH3OOH,

        # d[HCHO]/dt = k6[CH₃O][O₂] + k8[CH₃OOH][OH] - k10[HCHO][OH] - j11[HCHO] - j12[HCHO]
        D(HCHO) ~ k6 * CH3O * O2_fixed + k8 * CH3OOH * OH - k10 * HCHO * OH - j11 * HCHO - j12 * HCHO,

        # d[HCO]/dt = k10[HCHO][OH] + j11[HCHO] - k13[HCO][O₂]
        D(HCO) ~ k10 * HCHO * OH + j11 * HCHO - k13 * HCO * O2_fixed,

        # d[CO]/dt = k13[HCO][O₂] + j12[HCHO] - k_CO_OH[CO][OH]
        D(CO) ~ k13 * HCO * O2_fixed + j12 * HCHO - k_CO_OH * CO * OH,

        # d[H₂]/dt = j12[HCHO]
        D(H2) ~ j12 * HCHO,

        # d[H]/dt = j11[HCHO] - k14[H][O₂]
        D(H) ~ j11 * HCHO - k14_0 * M_fixed * H * O2_fixed,

        # d[OH]/dt = P_OH_ext + k15[HO₂][NO] + j9[CH₃OOH]
        #            - k1[CH₄][OH] - k7[CH₃OOH][OH] - k10[HCHO][OH]
        #            - k_CO_OH[CO][OH] - k_OH_NO2[OH][NO₂]
        D(OH) ~ P_OH_ext + k15 * HO2 * NO + j9 * CH3OOH - k1 * CH4 * OH - k7 * CH3OOH * OH - k8 * CH3OOH * OH - k10 * HCHO * OH - k_CO_OH * CO * OH - k_OH_NO2 * OH * NO2,

        # d[HO₂]/dt = k6[CH₃O][O₂] + k13[HCO][O₂] + k14[H][O₂] + k_CO_OH[CO][OH]
        #             - k4[CH₃O₂][HO₂] - k15[HO₂][NO] - 2k_HO2_HO2[HO₂]²
        D(HO2) ~ k6 * CH3O * O2_fixed + k13 * HCO * O2_fixed + k14_0 * M_fixed * H * O2_fixed + k_CO_OH * CO * OH - k4 * CH3O2 * HO2 - k15 * HO2 * NO - 2 * k_HO2_HO2 * HO2^2,

        # d[NO]/dt = j16[NO₂] - k3[CH₃O₂][NO] - k15[HO₂][NO] - k_NO_O3[NO][O₃]
        D(NO) ~ j16 * NO2 - k3 * CH3O2 * NO - k15 * HO2 * NO - k_NO_O3 * NO * O3,

        # d[NO₂]/dt = k3[CH₃O₂][NO] + k15[HO₂][NO] + k_NO_O3[NO][O₃] - j16[NO₂] - k_OH_NO2[OH][NO₂]
        D(NO2) ~ k3 * CH3O2 * NO + k15 * HO2 * NO + k_NO_O3 * NO * O3 - j16 * NO2 - k_OH_NO2 * OH * NO2,

        # d[O]/dt = j16[NO₂] - k17[O][O₂]
        D(O) ~ j16 * NO2 - k17_0 * M_fixed * O * O2_fixed,

        # d[O₃]/dt = k17[O][O₂] - k_NO_O3[NO][O₃]
        D(O3) ~ k17_0 * M_fixed * O * O2_fixed - k_NO_O3 * NO * O3,

        # d[HNO₃]/dt = k_OH_NO2[OH][NO₂]
        D(HNO3) ~ k_OH_NO2 * OH * NO2,

        # d[H₂O₂]/dt = k_HO2_HO2[HO₂]²
        D(H2O2) ~ k_HO2_HO2 * HO2^2,
    ]

    return ODESystem(eqs, t; name=name)
end

"""
    Combined Tropospheric Chemistry System

Integrates all the individual chemistry systems from Chapter 6 into a
comprehensive tropospheric chemistry model.

This combined system includes:
- OH Production from O₃ photolysis (Section 6.1)
- NOx Photochemical Cycle (Section 6.2)
- CO Oxidation Chemistry (Section 6.3)
- Methane Oxidation Mechanism (Section 6.4)

The systems are coupled through shared species (OH, HO₂, NO, NO₂, O₃).

Reference: Seinfeld & Pandis (2006), Chapter 6
"""

using ModelingToolkit
using Unitful
using ModelingToolkit: t_nounits as t, D_nounits as D

"""
    TroposphericChemistrySystem(; name)

Combined ModelingToolkit System for tropospheric chemistry diagnostics.

This system couples the OH production, NOx cycling, and CO/CH₄ oxidation
mechanisms to create a comprehensive diagnostic model of tropospheric photochemistry.

# Key coupling points:
- OH is produced from O₃ photolysis and HO₂ + NO reaction
- OH is consumed by CO, CH₄, and other species
- HO₂ is produced from CO/VOC oxidation
- HO₂ is lost to NO (producing OH) and self-reaction
- NOx cycles between NO and NO₂ through O₃ and peroxy radicals
- O₃ is produced when peroxy radicals oxidize NO to NO₂

# Input Variables
All species concentrations must be provided as inputs.

# Diagnostic outputs:
- P_O3_net: Net ozone production rate
- OPE: Ozone production efficiency
- HOx: Total HOx (OH + HO₂)
- chain_length: HOx chain length
"""
function TroposphericChemistrySystem(; name=:TroposphericChemistrySystem)
    # =========================================================================
    # Parameters
    # =========================================================================
    @parameters begin
        # OH Production (from O₃ photolysis)
        j_O3 = 1e-5, [description = "O₃ photolysis rate [s⁻¹]"]
        k_O1D_N2 = 2.6e-11, [description = "O(¹D) + N₂ rate [cm³/s]"]
        k_O1D_O2 = 4.0e-11, [description = "O(¹D) + O₂ rate [cm³/s]"]
        k_O1D_H2O = 2.2e-10, [description = "O(¹D) + H₂O rate [cm³/s]"]

        # NOx photochemistry
        j_NO2 = 8e-3, [description = "NO₂ photolysis rate [s⁻¹]"]
        k_O_O2_M = 6.0e-34, [description = "O + O₂ + M rate [cm⁶/s]"]
        k_NO_O3 = 1.8e-14, [description = "NO + O₃ rate [cm³/s]"]

        # CO oxidation
        k_CO_OH = 2.4e-13, [description = "CO + OH rate [cm³/s]"]
        k_HO2_NO = 8.1e-12, [description = "HO₂ + NO rate [cm³/s]"]
        k_HO2_HO2 = 2.9e-12, [description = "HO₂ + HO₂ rate [cm³/s]"]
        k_OH_NO2 = 1.0e-11, [description = "OH + NO₂ rate [cm³/s]"]
        k_HO2_O3 = 2.0e-15, [description = "HO₂ + O₃ rate [cm³/s]"]
        k_OH_O3 = 7.3e-14, [description = "OH + O₃ rate [cm³/s]"]

        # CH₄ oxidation
        k_CH4_OH = 6.3e-15, [description = "CH₄ + OH rate [cm³/s]"]
        k_CH3O2_NO = 7.7e-12, [description = "CH₃O₂ + NO rate [cm³/s]"]

        # Air composition
        f_N2 = 0.78, [description = "N₂ fraction"]
        f_O2 = 0.21, [description = "O₂ fraction"]
    end

    # =========================================================================
    # Input Variables (species concentrations)
    # =========================================================================
    @variables begin
        # Major species
        O3(t), [description = "Ozone [molecules/cm³]"]
        NO(t), [description = "Nitric oxide [molecules/cm³]"]
        NO2(t), [description = "Nitrogen dioxide [molecules/cm³]"]
        OH(t), [description = "Hydroxyl radical [molecules/cm³]"]
        HO2(t), [description = "Hydroperoxy radical [molecules/cm³]"]
        CO(t), [description = "Carbon monoxide [molecules/cm³]"]
        CH4(t), [description = "Methane [molecules/cm³]"]
        CH3O2(t), [description = "Methylperoxy radical [molecules/cm³]"]
        H2O(t), [description = "Water vapor [molecules/cm³]"]
        M(t), [description = "Total air density [molecules/cm³]"]
        O2(t), [description = "Molecular oxygen [molecules/cm³]"]
    end

    # Effective O(¹D) quenching rate for air
    k_O1D_M = f_N2 * k_O1D_N2 + f_O2 * k_O1D_O2

    # =========================================================================
    # Output Variables (diagnostics)
    # =========================================================================
    @variables begin
        # Intermediate species (steady-state)
        O1D(t), [description = "O(¹D) [molecules/cm³]"]
        O(t), [description = "O(³P) [molecules/cm³]"]

        # Derived quantities from Section 6.1 (OH Production)
        P_OH_O3(t), [description = "OH production from O₃ [molecules/cm³/s]"]
        ε_OH(t), [description = "OH yield from O(¹D)"]

        # Derived quantities from Section 6.2 (NOx)
        O3_pss(t), [description = "Photostationary O₃ [molecules/cm³]"]
        Φ(t), [description = "Photostationary state parameter"]

        # Derived quantities from Section 6.3 (CO oxidation)
        P_O3_CO(t), [description = "O₃ production from CO oxidation [molecules/cm³/s]"]
        L_HOx(t), [description = "HOx loss rate [molecules/cm³/s]"]

        # Combined diagnostics
        NOx(t), [description = "Total NOx [molecules/cm³]"]
        HOx(t), [description = "Total HOx [molecules/cm³]"]
        RO2(t), [description = "Total organic peroxy radicals [molecules/cm³]"]
        P_O3_total(t), [description = "Total O₃ production [molecules/cm³/s]"]
        L_O3_total(t), [description = "Total O₃ loss [molecules/cm³/s]"]
        P_O3_net(t), [description = "Net O₃ tendency [molecules/cm³/s]"]
        OPE(t), [description = "Ozone production efficiency"]
        chain_length(t), [description = "HOx chain length"]
        L_NOx(t), [description = "NOx loss rate [molecules/cm³/s]"]
    end

    # =========================================================================
    # Equations
    # =========================================================================
    eqs = [
        # --- Section 6.1: OH Production ---

        # Eq. 6.1: O(¹D) steady-state
        O1D ~ j_O3 * O3 / (k_O1D_M * M + k_O1D_H2O * H2O),

        # Eq. 6.4: OH yield
        ε_OH ~ k_O1D_H2O * H2O / (k_O1D_M * M + k_O1D_H2O * H2O),

        # Eq. 6.3: OH production from O₃ photolysis
        P_OH_O3 ~ 2 * j_O3 * O3 * ε_OH,

        # --- Section 6.2: NOx Photochemistry ---

        # Eq. 6.5: O atom steady-state
        O ~ j_NO2 * NO2 / (k_O_O2_M * O2 * M),

        # Eq. 6.6: Photostationary state O₃
        O3_pss ~ j_NO2 * NO2 / (k_NO_O3 * NO),

        # Eq. 6.7: Photostationary state parameter
        Φ ~ j_NO2 * NO2 / (k_NO_O3 * NO * O3),

        # --- Section 6.3: CO/HOx Chemistry ---

        # O₃ production from HO₂ + NO (main pathway)
        P_O3_CO ~ k_HO2_NO * HO2 * NO,

        # HOx termination (OH + NO₂ and HO₂ + HO₂)
        L_HOx ~ k_OH_NO2 * OH * NO2 + 2 * k_HO2_HO2 * HO2^2,

        # HOx chain length
        chain_length ~ P_O3_CO / L_HOx,

        # --- Combined Diagnostics ---

        # Total NOx and HOx
        NOx ~ NO + NO2,
        HOx ~ OH + HO2,

        # RO₂ (simplified, just CH₃O₂ for now)
        RO2 ~ CH3O2,

        # Total O₃ production (from all peroxy + NO reactions)
        P_O3_total ~ k_HO2_NO * HO2 * NO + k_CH3O2_NO * CH3O2 * NO,

        # Total O₃ loss
        L_O3_total ~ k_NO_O3 * NO * O3 + k_OH_O3 * OH * O3 + k_HO2_O3 * HO2 * O3,

        # Net O₃ production
        P_O3_net ~ P_O3_total - L_O3_total,

        # NOx loss (mainly HNO₃ formation)
        L_NOx ~ k_OH_NO2 * OH * NO2,

        # Ozone Production Efficiency
        OPE ~ P_O3_total / L_NOx,
    ]

    return ODESystem(eqs, t; name=name)
end

"""
    get_typical_conditions()

Returns a dictionary of typical lower troposphere conditions for use
with the TroposphericChemistrySystem.

Based on values from Seinfeld & Pandis Chapter 6.
"""
function get_typical_conditions()
    return Dict(
        # Concentrations in molecules cm⁻³
        :M => 2.5e19,      # Total air at STP
        :O2 => 5.25e18,    # 21% of M
        :H2O => 4e17,      # ~1% relative humidity equivalent
        :O3 => 1e12,       # ~40 ppb
        :NO => 2.5e9,      # ~0.1 ppb
        :NO2 => 2.5e10,    # ~1 ppb
        :CO => 2.5e12,     # ~100 ppb
        :CH4 => 4.5e13,    # ~1800 ppb
        :OH => 1e6,        # typical daytime
        :HO2 => 1e8,       # typical daytime
        :CH3O2 => 1e8,     # typical daytime
    )
end

"""
    get_urban_conditions()

Returns a dictionary of typical urban conditions with elevated NOx.
"""
function get_urban_conditions()
    return Dict(
        :M => 2.5e19,
        :O2 => 5.25e18,
        :H2O => 4e17,
        :O3 => 2e12,       # ~80 ppb (can be high in urban areas)
        :NO => 2.5e11,     # ~10 ppb
        :NO2 => 7.5e11,    # ~30 ppb
        :CO => 5e13,       # ~2 ppm
        :CH4 => 4.5e13,    # ~1800 ppb
        :OH => 5e5,        # reduced due to high NOx
        :HO2 => 5e7,       # reduced due to high NOx
        :CH3O2 => 5e7,
    )
end

"""
    get_remote_conditions()

Returns a dictionary of typical remote/background conditions with low NOx.
"""
function get_remote_conditions()
    return Dict(
        :M => 2.5e19,
        :O2 => 5.25e18,
        :H2O => 4e17,
        :O3 => 7.5e11,     # ~30 ppb
        :NO => 2.5e8,      # ~10 ppt
        :NO2 => 5e8,       # ~20 ppt
        :CO => 2e12,       # ~80 ppb
        :CH4 => 4.5e13,    # ~1800 ppb
        :OH => 1e6,
        :HO2 => 2e8,       # higher due to low NOx
        :CH3O2 => 2e8,
    )
end

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

using ModelingToolkit: t, D

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

# Subsystem Composition:

  - `oh`: OHProduction subsystem (Section 6.1)
  - `nox`: NOxPhotochemistry subsystem (Section 6.2)
  - `co`: COOxidation subsystem (Section 6.3)

# Input Variables

All species concentrations must be provided as inputs [m⁻³].

# Diagnostic outputs:

  - P_O3_net: Net ozone production rate [m⁻³ s⁻¹]
  - OPE: Ozone production efficiency [dimensionless]
  - HOx: Total HOx (OH + HO₂) [m⁻³]    # =========================================================================
  - chain_length: HOx chain length [dimensionless]    # Subsystems (composed from individual component functions)
"""
@component function TroposphericChemistrySystem(; name = :TroposphericChemistrySystem)
    # =========================================================================
    # Subsystems (composed from individual component functions)
    # =========================================================================
    oh_sys = OHProduction(; name = :oh)
    nox_sys = NOxPhotochemistry(; name = :nox)
    co_sys = COOxidation(; name = :co)

    # =========================================================================
    # Parameters (additional parameters for combined diagnostics, in SI)
    # =========================================================================
    @parameters begin
        # CH₄ oxidation
        k_CH4_OH = 6.3e-15 * 1e-6,
        [description = "CH₄ + OH rate constant (6.3e-15 cm³/molec/s)", unit = u"m^3/s"]
        k_CH3O2_NO = 7.7e-12 * 1e-6,
        [description = "CH₃O₂ + NO rate constant (7.7e-12 cm³/molec/s)", unit = u"m^3/s"]
    end

    # =========================================================================
    # Input Variables (species concentrations in SI: m⁻³)
    # =========================================================================
    @variables begin
        # Major species
        O3(t), [description = "Ozone", unit = u"m^-3"]
        NO(t), [description = "Nitric oxide", unit = u"m^-3"]
        NO2(t), [description = "Nitrogen dioxide", unit = u"m^-3"]
        OH(t), [description = "Hydroxyl radical", unit = u"m^-3"]
        HO2(t), [description = "Hydroperoxy radical", unit = u"m^-3"]
        CO(t), [description = "Carbon monoxide", unit = u"m^-3"]
        CH4(t), [description = "Methane", unit = u"m^-3"]
        CH3O2(t), [description = "Methylperoxy radical", unit = u"m^-3"]
        H2O(t), [description = "Water vapor", unit = u"m^-3"]
        M(t), [description = "Total air density", unit = u"m^-3"]
        O2(t), [description = "Molecular oxygen", unit = u"m^-3"]
    end

    # =========================================================================
    # Output Variables (combined diagnostics)
    # =========================================================================
    @variables begin
        # Combined diagnostics
        NOx(t), [description = "Total NOx", unit = u"m^-3"]
        HOx(t), [description = "Total HOx", unit = u"m^-3"]
        RO2(t), [description = "Total organic peroxy radicals", unit = u"m^-3"]
        P_O3_total(t), [description = "Total O₃ production", unit = u"m^-3*s^-1"]
        L_O3_total(t), [description = "Total O₃ loss", unit = u"m^-3*s^-1"]
        P_O3_net(t), [description = "Net O₃ tendency", unit = u"m^-3*s^-1"]
        OPE(t), [description = "Ozone production efficiency (dimensionless)", unit = u"1"]
        chain_length(t), [description = "HOx chain length (dimensionless)", unit = u"1"]
        L_NOx(t), [description = "NOx loss rate", unit = u"m^-3*s^-1"]
    end

    # =========================================================================
    # Coupling Equations
    # =========================================================================
    # Connect shared species from the combined level to subsystem inputs.
    # Subsystem variables are accessed via dot notation (e.g., oh_sys.O3).
    eqs = [
        # --- Coupling: Feed shared species into OH Production subsystem ---
        oh_sys.O3 ~ O3,
        oh_sys.H2O ~ H2O,
        oh_sys.M ~ M,

        # --- Coupling: Feed shared species into NOx Photochemistry subsystem ---
        nox_sys.NO ~ NO,
        nox_sys.NO2 ~ NO2,
        nox_sys.O3 ~ O3,
        nox_sys.O2 ~ O2,
        nox_sys.M ~ M,

        # --- Coupling: Feed shared species into CO Oxidation subsystem ---
        co_sys.CO ~ CO,
        co_sys.OH ~ OH,
        co_sys.HO2 ~ HO2,
        co_sys.NO ~ NO,
        co_sys.NO2 ~ NO2,
        co_sys.O3 ~ O3,

        # --- Combined Diagnostics ---

        # Total NOx and HOx
        NOx ~ NO + NO2,
        HOx ~ OH + HO2,

        # RO₂ (simplified, just CH₃O₂ for now)
        RO2 ~ CH3O2,

        # Total O₃ production (from all peroxy + NO reactions)
        # Uses co_sys rate constants for HO₂ + NO and adds CH₃O₂ + NO contribution
        P_O3_total ~ co_sys.k_HO2_NO * HO2 * NO + k_CH3O2_NO * CH3O2 * NO,

        # Total O₃ loss
        L_O3_total ~
        nox_sys.k_NO_O3 * NO * O3 + co_sys.k_OH_O3 * OH * O3 + co_sys.k_HO2_O3 * HO2 * O3,

        # Net O₃ production
        P_O3_net ~ P_O3_total - L_O3_total,

        # NOx loss (mainly HNO₃ formation)
        L_NOx ~ co_sys.k_OH_NO2 * OH * NO2,

        # Ozone Production Efficiency
        OPE ~ P_O3_total / L_NOx,

        # HOx chain length (from CO oxidation subsystem diagnostic, re-expressed at combined level)
        chain_length ~ co_sys.chain_length
    ]

    return System(eqs, t; systems = [oh_sys, nox_sys, co_sys], name)
end

"""
    TypicalConditions(; name)

ModelingToolkit System representing typical lower troposphere conditions.

All concentrations are parameters in SI units (m⁻³).

Based on values from Seinfeld & Pandis Chapter 6.
"""
@component function TypicalConditions(; name = :TypicalConditions)
    @parameters begin
        M = 2.5e25, [description = "Total air at STP", unit = u"m^-3"]
        O2 = 5.25e24, [description = "O₂, 21% of M", unit = u"m^-3"]
        H2O = 4e23, [description = "Water vapor, ~50% RH at 298 K", unit = u"m^-3"]
        O3 = 1e18, [description = "Ozone, ~40 ppb", unit = u"m^-3"]
        NO = 2.5e15, [description = "NO, ~0.1 ppb", unit = u"m^-3"]
        NO2 = 2.5e16, [description = "NO₂, ~1 ppb", unit = u"m^-3"]
        CO = 2.5e18, [description = "CO, ~100 ppb", unit = u"m^-3"]
        CH4 = 4.5e19, [description = "CH₄, ~1800 ppb", unit = u"m^-3"]
        OH = 1e12, [description = "OH, typical daytime", unit = u"m^-3"]
        HO2 = 1e14, [description = "HO₂, typical daytime", unit = u"m^-3"]
        CH3O2 = 1e14, [description = "CH₃O₂, typical daytime", unit = u"m^-3"]
    end
    @variables begin
        M_out(t), [description = "Total air number density", unit = u"m^-3"]
        O2_out(t), [description = "O₂ concentration", unit = u"m^-3"]
        H2O_out(t), [description = "Water vapor concentration", unit = u"m^-3"]
        O3_out(t), [description = "Ozone concentration", unit = u"m^-3"]
        NO_out(t), [description = "NO concentration", unit = u"m^-3"]
        NO2_out(t), [description = "NO₂ concentration", unit = u"m^-3"]
        CO_out(t), [description = "CO concentration", unit = u"m^-3"]
        CH4_out(t), [description = "CH₄ concentration", unit = u"m^-3"]
        OH_out(t), [description = "OH concentration", unit = u"m^-3"]
        HO2_out(t), [description = "HO₂ concentration", unit = u"m^-3"]
        CH3O2_out(t), [description = "CH₃O₂ concentration", unit = u"m^-3"]
    end
    eqs = [
        M_out ~ M, O2_out ~ O2, H2O_out ~ H2O, O3_out ~ O3,
        NO_out ~ NO, NO2_out ~ NO2, CO_out ~ CO, CH4_out ~ CH4,
        OH_out ~ OH, HO2_out ~ HO2, CH3O2_out ~ CH3O2,
    ]
    return System(eqs, t; name)
end

"""
    UrbanConditions(; name)

ModelingToolkit System representing typical urban conditions with elevated NOx.

All concentrations are parameters in SI units (m⁻³).
"""
@component function UrbanConditions(; name = :UrbanConditions)
    @parameters begin
        M = 2.5e25, [description = "Total air at STP", unit = u"m^-3"]
        O2 = 5.25e24, [description = "O₂, 21% of M", unit = u"m^-3"]
        H2O = 4e23, [description = "Water vapor", unit = u"m^-3"]
        O3 = 2e18, [description = "Ozone, ~80 ppb", unit = u"m^-3"]
        NO = 2.5e17, [description = "NO, ~10 ppb", unit = u"m^-3"]
        NO2 = 7.5e17, [description = "NO₂, ~30 ppb", unit = u"m^-3"]
        CO = 5e19, [description = "CO, ~2 ppm", unit = u"m^-3"]
        CH4 = 4.5e19, [description = "CH₄, ~1800 ppb", unit = u"m^-3"]
        OH = 5e11, [description = "OH, reduced due to high NOx", unit = u"m^-3"]
        HO2 = 5e13, [description = "HO₂, reduced due to high NOx", unit = u"m^-3"]
        CH3O2 = 5e13, [description = "CH₃O₂", unit = u"m^-3"]
    end
    @variables begin
        M_out(t), [description = "Total air number density", unit = u"m^-3"]
        O2_out(t), [description = "O₂ concentration", unit = u"m^-3"]
        H2O_out(t), [description = "Water vapor concentration", unit = u"m^-3"]
        O3_out(t), [description = "Ozone concentration", unit = u"m^-3"]
        NO_out(t), [description = "NO concentration", unit = u"m^-3"]
        NO2_out(t), [description = "NO₂ concentration", unit = u"m^-3"]
        CO_out(t), [description = "CO concentration", unit = u"m^-3"]
        CH4_out(t), [description = "CH₄ concentration", unit = u"m^-3"]
        OH_out(t), [description = "OH concentration", unit = u"m^-3"]
        HO2_out(t), [description = "HO₂ concentration", unit = u"m^-3"]
        CH3O2_out(t), [description = "CH₃O₂ concentration", unit = u"m^-3"]
    end
    eqs = [
        M_out ~ M, O2_out ~ O2, H2O_out ~ H2O, O3_out ~ O3,
        NO_out ~ NO, NO2_out ~ NO2, CO_out ~ CO, CH4_out ~ CH4,
        OH_out ~ OH, HO2_out ~ HO2, CH3O2_out ~ CH3O2,
    ]
    return System(eqs, t; name)
end

"""
    RemoteConditions(; name)

ModelingToolkit System representing typical remote/background conditions with low NOx.

All concentrations are parameters in SI units (m⁻³).
"""
@component function RemoteConditions(; name = :RemoteConditions)
    @parameters begin
        M = 2.5e25, [description = "Total air at STP", unit = u"m^-3"]
        O2 = 5.25e24, [description = "O₂, 21% of M", unit = u"m^-3"]
        H2O = 4e23, [description = "Water vapor", unit = u"m^-3"]
        O3 = 7.5e17, [description = "Ozone, ~30 ppb", unit = u"m^-3"]
        NO = 2.5e14, [description = "NO, ~10 ppt", unit = u"m^-3"]
        NO2 = 5e14, [description = "NO₂, ~20 ppt", unit = u"m^-3"]
        CO = 2e18, [description = "CO, ~80 ppb", unit = u"m^-3"]
        CH4 = 4.5e19, [description = "CH₄, ~1800 ppb", unit = u"m^-3"]
        OH = 1e12, [description = "OH, typical daytime", unit = u"m^-3"]
        HO2 = 2e14, [description = "HO₂, higher due to low NOx", unit = u"m^-3"]
        CH3O2 = 2e14, [description = "CH₃O₂", unit = u"m^-3"]
    end
    @variables begin
        M_out(t), [description = "Total air number density", unit = u"m^-3"]
        O2_out(t), [description = "O₂ concentration", unit = u"m^-3"]
        H2O_out(t), [description = "Water vapor concentration", unit = u"m^-3"]
        O3_out(t), [description = "Ozone concentration", unit = u"m^-3"]
        NO_out(t), [description = "NO concentration", unit = u"m^-3"]
        NO2_out(t), [description = "NO₂ concentration", unit = u"m^-3"]
        CO_out(t), [description = "CO concentration", unit = u"m^-3"]
        CH4_out(t), [description = "CH₄ concentration", unit = u"m^-3"]
        OH_out(t), [description = "OH concentration", unit = u"m^-3"]
        HO2_out(t), [description = "HO₂ concentration", unit = u"m^-3"]
        CH3O2_out(t), [description = "CH₃O₂ concentration", unit = u"m^-3"]
    end
    eqs = [
        M_out ~ M, O2_out ~ O2, H2O_out ~ H2O, O3_out ~ O3,
        NO_out ~ NO, NO2_out ~ NO2, CO_out ~ CO, CH4_out ~ CH4,
        OH_out ~ OH, HO2_out ~ HO2, CH3O2_out ~ CH3O2,
    ]
    return System(eqs, t; name)
end

"""
    get_conditions_dict(sys)

Helper to extract default parameter values from a conditions System as a Dict{Symbol,Float64}.
This is a convenience function for use in tests and analysis code.
"""
function get_conditions_dict(sys)
    params = parameters(sys)
    result = Dict{Symbol,Float64}()
    for p in params
        if ModelingToolkit.hasdefault(p)
            result[Symbol(p)] = Float64(ModelingToolkit.getdefault(p))
        end
    end
    return result
end

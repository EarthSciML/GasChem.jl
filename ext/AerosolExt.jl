module AerosolExt
using GasChem, Aerosol, ModelingToolkit, EarthSciMLBase, DynamicQuantities

# Register custom units if not already registered by GasChem.
# (Extension precompilation may run before GasChem's unit registration takes effect.)
if !(:molec in DynamicQuantities.Units.UNIT_SYMBOLS)
    @register_unit molec 1
end
if !(:mol_air in DynamicQuantities.Units.UNIT_SYMBOLS)
    @register_unit mol_air 1u"mol"
end
if !(:ppb in DynamicQuantities.Units.UNIT_SYMBOLS)
    @register_unit ppb 1u"mol/mol_air"
end

"""
    couple2(c::GEOSChemGasPhaseCoupler, a::AerosolDistributionCoupler)

Couple the GEOS-Chem gas-phase mechanism to an [`Aerosol.AerosolDistribution`](@ref)
so the heterogeneous uptake rate constants `k_het_*` are driven by the aerosol
surface-area concentration rather than left at their inert default of 0.

Each first-order uptake frequency follows the standard reactive-uptake form
(GEOS-Chem `Ars_L1k`, kinetic limit):

    k_het_X = S_t * c_bar_X * gamma_X / 4

where

  - `S_t`      = total aerosol surface-area concentration [m^-1] from the Aerosol
                 size distribution (Seinfeld & Pandis 2006, Eq. 8.5);
  - `c_bar_X`  = mean molecular speed of species X [m s^-1] from kinetic theory
                 (S&P Eq. 12.24), evaluated at the gas-phase temperature `c.T`;
  - `gamma_X`  = reactive uptake coefficient of X (representative constants below;
                 GEOS-Chem uses RH/composition-dependent parameterizations).

Restores the heterogeneous NOx/HOx sinks absent from the gas-phase-only port:
N2O5 + H2O(aq) -> 2 HNO3 (dominant nighttime NOx sink), HO2 uptake (HOx sink),
and NO2 / NO3 uptake.
"""
function EarthSciMLBase.couple2(
        c::GasChem.GEOSChemGasPhaseCoupler,
        a::Aerosol.AerosolDistributionCoupler
    )
    c, a = c.sys, a.sys

    @constants(
        k_B = 1.380649e-23,
        [unit = u"J/K", description = "Boltzmann constant"],
        N_A = 6.02214076e23,
        [unit = u"mol^-1", description = "Avogadro constant"],
        MW_N2O5 = 0.10801,
        [unit = u"kg/mol", description = "N2O5 molar mass"],
        MW_HO2 = 0.03301,
        [unit = u"kg/mol", description = "HO2 molar mass"],
        MW_NO2 = 0.046006,
        [unit = u"kg/mol", description = "NO2 molar mass"],
        MW_NO3 = 0.062004,
        [unit = u"kg/mol", description = "NO3 molar mass"],
        gamma_N2O5 = 0.02,
        [unit = u"1", description = "N2O5 reactive uptake coefficient (representative)"],
        gamma_HO2 = 0.2,
        [unit = u"1", description = "HO2 reactive uptake coefficient (representative, Mao2013)"],
        gamma_NO2 = 1.0e-4,
        [unit = u"1", description = "NO2 reactive uptake coefficient (representative)"],
        gamma_NO3 = 0.01,
        [unit = u"1", description = "NO3 reactive uptake coefficient (representative)"],
    )

    # Mean molecular speed of species with molar mass MW (kinetic theory, S&P Eq. 12.24)
    # evaluated at the gas-phase temperature.
    c_bar(MW) = sqrt(8 * k_B * c.T * N_A / (π * MW))

    c = param_to_var(c, :k_het_N2O5, :k_het_HO2, :k_het_NO2, :k_het_NO3)
    return ConnectorSystem(
        [
            c.k_het_N2O5 ~ a.S_t * c_bar(MW_N2O5) * gamma_N2O5 / 4,
            c.k_het_HO2 ~ a.S_t * c_bar(MW_HO2) * gamma_HO2 / 4,
            c.k_het_NO2 ~ a.S_t * c_bar(MW_NO2) * gamma_NO2 / 4,
            c.k_het_NO3 ~ a.S_t * c_bar(MW_NO3) * gamma_NO3 / 4,
        ],
        c, a
    )
end

"""
    couple2(c::GEOSChemGasPhaseCoupler, cc::CloudChemistryFixedpHCoupler)

Couple the GEOS-Chem gas-phase mechanism to [`Aerosol.CloudChemistryFixedpH`](@ref)
so the in-cloud sulfate rate constant `k_cld1` (SO2 + H2O2 -> SO4) is driven by the
aqueous-phase S(IV)+H2O2 oxidation rate rather than left at its inert default of 0.

The aqueous oxidation rate `R_H2O2` [mol m^-3 s^-1] is converted to a gas-phase loss
frequency (S&P Eq. 7.75) and divided by the gas concentrations to recover the
second-order rate constant the ported KPP reaction expects:

    k_cld1 = FC * (1e9 * w_L * R_gas * T * R_H2O2 / P) / ([SO2] * [H2O2])

The H2O2 pathway is pH-insensitive, so a prescribed cloud pH is used (no per-cell
electroneutrality solve). Cloud fraction `FC` and in-cloud LWC `L` are prescribed
representative constants here; in a CTM they should be supplied from the GEOSFP
A3cld cloud fields. `k_cld2` (O3 path) and `k_cld3` (Fe/Mn) are left at 0 pending the
dynamic-pH coupling. Other CloudChemistry inputs (Fe/Mn/NH3/CO2/HNO3) are prescribed
representative values; only the H2O2 path is consumed.
"""
function EarthSciMLBase.couple2(
        c::GasChem.GEOSChemGasPhaseCoupler,
        cc::Aerosol.CloudChemistryFixedpHCoupler
    )
    c, cc = c.sys, cc.sys

    @constants(
        inv_ppb = 1.0e-9,
        [unit = u"1", description = "ppb number -> mole fraction"],
        ppb_per1 = 1.0e9,
        [unit = u"ppb", description = "mole fraction -> ppb"],
        R_gas = 8.31446,
        [unit = u"m^3*Pa/mol/K", description = "Gas constant"],
        rho_w_inv = 1.0e-6,
        [unit = u"m^3/g", description = "Inverse water density (LWC conversion)"],
        FC = 1.0,
        [unit = u"1", description = "Cloud fraction (prescribed; use GEOSFP A3cld in a CTM)"],
        L_cld = 0.3,
        [unit = u"g/m^3", description = "In-cloud liquid water content (prescribed)"],
        pH_cld = 4.5,
        [unit = u"1", description = "Cloud droplet pH (prescribed; H2O2 path is pH-insensitive)"],
        Fe_cld = 1.0e-9,
        [unit = u"mol/m^3", description = "Fe(III) (prescribed ~0; H2O2 path unaffected)"],
        Mn_cld = 1.0e-9,
        [unit = u"mol/m^3", description = "Mn(II) (prescribed ~0)"],
        co2_frac = 4.0e-4,
        [unit = u"1", description = "CO2 mole fraction (~400 ppm)"],
        zero_pa = 0.0,
        [unit = u"Pa", description = "Zero partial pressure (NH3 not in mechanism)"],
        floor_c = 1.0e-3,
        [unit = u"ppb", description = "Denominator guard"],
    )

    # in-cloud liquid water volume mixing ratio (dimensionless)
    w_L = rho_w_inv * L_cld
    # aqueous oxidation rate (mol/m^3/s) -> gas-phase loss rate (ppb/s), per S&P Eq. 7.75
    Lpp(R) = ppb_per1 * w_L * R_gas * c.T * R / c.P

    c = param_to_var(c, :k_cld1)
    return ConnectorSystem(
        [
            cc.T ~ c.T,
            cc.pH_input ~ pH_cld,
            cc.L ~ L_cld,
            cc.xi_SO2 ~ c.SO2 * inv_ppb,
            cc.Fe_III ~ Fe_cld,
            cc.Mn_II ~ Mn_cld,
            cc.p_SO2 ~ c.SO2 * inv_ppb * c.P,
            cc.p_H2O2 ~ c.H2O2 * inv_ppb * c.P,
            cc.p_O3 ~ c.O3 * inv_ppb * c.P,
            cc.p_HNO3 ~ c.HNO3 * inv_ppb * c.P,
            cc.p_CO2 ~ co2_frac * c.P,
            cc.p_NH3 ~ zero_pa,
            c.k_cld1 ~ FC * Lpp(cc.R_H2O2) / ((c.SO2 + floor_c) * (c.H2O2 + floor_c)),
        ],
        c, cc
    )
end

end

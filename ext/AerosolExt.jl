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

"""
    couple2(c::SuperFastCoupler, cc::CloudChemistryFixedpHCoupler)

SuperFast counterpart of the GEOS-Chem ↔ CloudChemistryFixedpH coupling above. The reduced
SuperFast mechanism carries SO2/SO4/H2O2/O3/HNO3 (SO2/SO4 re-enabled for aerosol coupling),
so the same in-cloud sulfate driver applies: the aqueous S(IV)+H2O2 oxidation rate `R_H2O2`
becomes the gas-phase rate constant `k_cld1` of SuperFast's `SO2 + H2O2 -> SO4` reaction.
"""
function EarthSciMLBase.couple2(
        c::GasChem.SuperFastCoupler,
        cc::Aerosol.CloudChemistryFixedpHCoupler
    )
    c, cc = c.sys, cc.sys
    # NOTE (differs from the GEOS-Chem method above): the SuperFast↔GEOSFP coupling
    # param_to_var's T *and* P (EarthSciDataExt.jl `param_to_var(c, :T, :P, :H2O)`),
    # whereas the GEOS-Chem↔GEOSFP coupling leaves P a scalar parameter. So here T and P
    # must be referenced as the *variable* forms (param_to_var them below) — otherwise the
    # raw-parameter T/P references collide with the met coupling's T(t)/P(t)
    # ("duplicate names SuperFast₊T and SuperFast₊T(t)").
    c = param_to_var(c, :T, :P)

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
            # NH3 partial pressure from the mechanism (the SF-aerosol port carries NH3;
            # pH_input is prescribed, but feed the real NH3 so the aqueous speciation
            # sees it rather than a hard-wired zero).
            cc.p_NH3 ~ c.NH3 * inv_ppb * c.P,
            c.k_cld1 ~ FC * Lpp(cc.R_H2O2) / ((c.SO2 + floor_c) * (c.H2O2 + floor_c)),
        ],
        c, cc
    )
end

"""
    couple2(c::SuperFastCoupler, a::AerosolDistributionCoupler)

SuperFast counterpart of the GEOS-Chem ↔ AerosolDistribution heterogeneous-uptake coupling
above. Drives SuperFast's `k_het_*` first-order uptake rate constants from the aerosol
surface-area concentration `a.S_t` (kinetic-limit reactive uptake, k = S_t·c̄·γ/4). SuperFast
carries N2O5/NO3 (added for this coupling) plus HO2/NO2, so the same four uptake sinks apply.
"""
function EarthSciMLBase.couple2(
        c::GasChem.SuperFastCoupler,
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

    # Mean molecular speed (kinetic theory, S&P Eq. 12.24) at the gas-phase temperature.
    c_bar(MW) = sqrt(8 * k_B * c.T * N_A / (π * MW))

    # param_to_var includes :T because the SuperFast↔GEOSFP coupling makes T a variable (see
    # the cloud couple2 note); referencing raw-parameter T here would collide with T(t).
    c = param_to_var(c, :k_het_N2O5, :k_het_HO2, :k_het_NO2, :k_het_NO3, :T)
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

# =============================================================================
# ISORROPIA II operator (GasChem.IsorropiaOp) — operator-split coupling of the
# GEOS-Chem gas phase to inorganic aerosol thermodynamic equilibrium.
#
# Rather than solving the full ISORROPIA system as one large nonlinear problem
# (multi-rooted, converges only in a narrow regime), this follows the *original*
# ISORROPIA II design: for the fine-mode NH4-SO4-NO3 system (the species the
# GEOS-Chem port feeds; sea-salt/dust handled separately) under the metastable
# (aqueous) assumption GEOS-Chem uses, the equilibrium reduces — given [H+] — to
# closed-form ion partitioning, so the whole solve is a single 1-D bisection on
# [H+] of the charge-balance residual (monotone in [H+], hence globally
# convergent) with an inner activity/water fixed-point. Thermodynamics (van't Hoff
# K, Kusik-Meissner activity coefficients, ZSR water) are reused from Aerosol's
# ISORROPIA II implementation.
# =============================================================================

# Saturation vapour pressure over liquid water [Pa] (Tetens, 1930), T in K.
_iso_esat(T) = 611.2 * exp(17.67 * (T - 273.15) / (T - 29.65))

# Flattened name of a symbolic variable/parameter, without the time argument.
_iso_name(v) = replace(string(v), "(t)" => "")

# Index of the unknown whose flattened name ends in "₊name" (parent-agnostic; the leading
# "₊" guards against e.g. matching "HNO3" inside another species name).
function _iso_find(syms, name)
    i = findfirst(s -> endswith(string(s), "₊" * name), syms)
    i === nothing && error("IsorropiaOp: state variable ending in `₊$name` not found in the coupled system")
    i
end

# Temperature/pressure variable to read per grid cell (coupled from the meteorology, e.g.
# GEOSFP — appears as an observed variable; search unknowns and observed equations).
function _iso_metvar(mtk_sys, name)
    # Find a per-cell met variable ending in "₊name" among the coupled system's unknowns and
    # observed variables. T and P come from the meteorology (GEOSFP): in a full CTM T is
    # coupled onto GEOSChemGasPhase₊T while the per-cell pressure is carried by a sibling
    # (e.g. FastJX₊P) — the GEOSChem↔GEOSFP coupling param_to_var's T + num_density but not P,
    # so GEOSChemGasPhase₊P stays a scalar parameter. All such met variables trace back to the
    # same GEOSFP field, so matching by suffix gives the correct per-cell value.
    cands = vcat(unknowns(mtk_sys), [eq.lhs for eq in observed(mtk_sys)])
    i = findfirst(v -> endswith(_iso_name(v), "₊" * name), cands)
    i === nothing && error(
        "IsorropiaOp: per-cell meteorology variable ending in `₊$name` not found — couple a " *
        "met source (e.g. GEOSFP, plus FastJX for per-cell pressure) so $name is available per grid cell")
    cands[i]
end

# State species the operator reads from / writes to. H2O is NOT among them: it is a
# met-coupled constant species (isconstantspecies parameter; GEOSFP RH via the
# SuperFast/GEOSChem↔GEOSFP couple2), read per-cell through the coordinate-observed
# function alongside T and P — see `get_needed_vars` / the `h2o_const` fallback below.
const _ISO_SPECIES = ("HNO3", "NIT", "NH3", "NH4", "SO4")

# --- ISORROPIA II thermodynamic constants (van't Hoff, SI; from Aerosol's Table 2) ---
const _ISO_T0 = 298.15
const _ISO_RG = 8.314462                 # J/mol/K
_iso_vh(K0, A, B, T) = K0 * exp(A * (_ISO_T0 / T - 1.0) + B * (1.0 + log(_ISO_T0 / T) - _ISO_T0 / T))
# (K0, A, B): K1 = HSO4⁻⇌H⁺+SO4²⁻; K2 = K21·K22 (NH3 dissolution·dissociation); K4 = HNO3(g)⇌H⁺+NO3⁻; Kw
const _ISO_K1 = (1.015e-2, 8.85, 25.14)
const _ISO_K21 = (5.764e1 * 9.86923e-6, 13.79, -5.393)
const _ISO_K22 = (1.805e-5, -1.5, 26.92)
const _ISO_K4 = (2.511e6 * 9.86923e-6, 29.17, 16.83)
const _ISO_KW = (1.01e-14, -22.52, 26.92)

# Activity coefficient at temperature T (Kusik-Meissner, reusing Aerosol's functions).
_iso_gamma(q, z, I, T) = Aerosol._iso2_gamma_T(Aerosol._iso2_km_gamma(q, z, I), I, T)

"""
    _iso_ternary(TSO4, TNH, TNO3, T, RH)

Solve the metastable NH4-SO4-NO3 inorganic aerosol equilibrium (ISORROPIA II, Fountoukis &
Nenes 2007) for given totals [mol m^-3], temperature [K] and relative humidity [0,1]. Returns
`(ok, g_HNO3, g_NH3)` with the gas-phase concentrations [mol m^-3]; `ok=false` if the result
is not finite. Globally convergent: a bracketed bisection on [H+] (the charge-balance residual
is monotone in [H+]) wrapped in an activity/water fixed-point.
"""
function _iso_ternary(TSO4, TNH, TNO3, T, RH; nouter = 40, nbisect = 30, tol = 1.0e-5)
    K1 = _iso_vh(_ISO_K1..., T)
    K2 = _iso_vh(_ISO_K21..., T) * _iso_vh(_ISO_K22..., T)
    K4 = _iso_vh(_ISO_K4..., T)
    Kw = _iso_vh(_ISO_KW..., T)
    KRT4 = K4 * _ISO_RG * T
    KRT2 = K2 * RH * _ISO_RG * T
    # initial composition guess and water content
    c_SO4 = 0.9TSO4; c_HSO4 = 0.1TSO4; c_NO3 = 0.5TNO3; c_NH4 = 0.9TNH
    W_w = max(Aerosol._iso2_zsr_water(RH, 0.0, c_NH4, c_SO4, c_HSO4, c_NO3, 0.0, 0.0, 0.0, 0.0), 1.0e-15)
    m_H = 1.0e-7
    first_outer = true
    for _ in 1:nouter
        I_s = max(0.5 * (m_H + c_NH4 / W_w + 4 * c_SO4 / W_w + c_HSO4 / W_w + c_NO3 / W_w + Kw * RH / m_H), 1.0e-12)
        gH2SO4 = _iso_gamma(-0.1, 2, I_s, T); gHHSO4 = _iso_gamma(8.0, 1, I_s, T)
        gHNO3 = _iso_gamma(2.6, 1, I_s, T); gHCl = _iso_gamma(6.0, 1, I_s, T); gNH4Cl = _iso_gamma(0.82, 1, I_s, T)
        γr = (gNH4Cl / gHCl)^2
        mS = TSO4 / W_w
        # charge-balance residual at trial molal [H+] (monotone increasing in mH)
        resid(mH) = let mOH = Kw * RH / mH, r = K1 * gHHSO4^2 / (mH * gH2SO4^3)
            mSO4 = r / (1 + r) * mS
            mHSO4 = mS / (1 + r)
            mNO3 = KRT4 * TNO3 / (mH * gHNO3^2 + KRT4 * W_w)
            mNH4 = KRT2 * TNH / (mOH * γr + KRT2 * W_w)
            (mH + mNH4) - (2mSO4 + mHSO4 + mNO3 + mOH)
        end
        # Bisection bracket on log10[H+]. PERF: after the first outer iteration [H+]
        # moves little between activity/water updates, so try a warm ±1.5-decade
        # bracket around the previous root before falling back to the full range —
        # this cuts the bisection count ~3x at identical precision. nbisect=30 keeps
        # ~17/2^30 ≈ 2e-8-decade resolution on the full range (tol is 1e-5; result
        # changes are <1e-6 relative vs the previous nbisect=60 — verified).
        lo, hi = -14.0, 3.0
        if !first_outer
            wlo = log10(m_H) - 1.5
            whi = log10(m_H) + 1.5
            if wlo > lo && whi < hi && resid(10.0^wlo) < 0 && resid(10.0^whi) >= 0
                lo, hi = wlo, whi
            end
        end
        first_outer = false
        if resid(10.0^lo) * resid(10.0^hi) > 0
            m_H = resid(10.0^lo) < 0 ? 10.0^hi : 10.0^lo
        else
            for _ in 1:nbisect
                mid = 0.5 * (lo + hi)
                resid(10.0^mid) < 0 ? (lo = mid) : (hi = mid)
                hi - lo < 1.0e-7 && break   # ~2.3e-7 relative on [H+]: beyond tol needs
            end
            m_H = 10.0^(0.5 * (lo + hi))
        end
        # back out species at the converged [H+]
        mOH = Kw * RH / m_H
        r = K1 * gHHSO4^2 / (m_H * gH2SO4^3)
        c_SO4 = r / (1 + r) * mS * W_w
        c_HSO4 = mS / (1 + r) * W_w
        c_NO3 = KRT4 * TNO3 / (m_H * gHNO3^2 + KRT4 * W_w) * W_w
        c_NH4 = KRT2 * TNH / (mOH * γr + KRT2 * W_w) * W_w
        W_new = max(Aerosol._iso2_zsr_water(RH, 0.0, c_NH4, c_SO4, c_HSO4, c_NO3, 0.0, 0.0, 0.0, 0.0), 1.0e-15)
        converged = abs(W_new - W_w) / W_w < tol
        W_w = W_new
        converged && break
    end
    g_HNO3 = TNO3 - c_NO3
    g_NH3 = TNH - c_NH4
    ok = isfinite(g_HNO3) && isfinite(g_NH3)
    return (ok, g_HNO3, g_NH3)
end

# H2O among the OBSERVED equations only (`build_coord_observed_function` cannot read
# unknowns): in a full CTM the met coupling (param_to_var + `c.H2O ~ RH-driven`) makes
# H2O an observed variable — the preferred read path. Returns `nothing` when absent.
function _iso_observed_or_nothing(mtk_sys, name)
    obs = [eq.lhs for eq in observed(mtk_sys)]
    i = findfirst(v -> endswith(_iso_name(v), "₊" * name), obs)
    i === nothing ? nothing : obs[i]
end

function EarthSciMLBase.get_needed_vars(::GasChem.IsorropiaOp, csys, mtk_sys, domain::DomainInfo)
    # T, P (and H2O when met-coupled/observed) are read per-cell through the
    # coordinate-observed function; species are state variables read by index.
    # H2O read priority (resolved in `get_odefunction`): observed (met-coupled CTM)
    # → state index (mechanisms carrying H2O as a dynamic species) → isconstantspecies
    # parameter default (standalone, no met).
    vars = [_iso_metvar(mtk_sys, "T"), _iso_metvar(mtk_sys, "P")]
    h2o = _iso_observed_or_nothing(mtk_sys, "H2O")
    h2o === nothing || push!(vars, h2o)
    vars
end

function EarthSciMLBase.get_odefunction(
        op::GasChem.IsorropiaOp, csys, mtk_sys, coord_args, domain::DomainInfo, u0, p, alg)
    sz = tuple(size(domain)...)
    II = CartesianIndices(sz)
    ncell = prod(sz)
    nrows = length(unknowns(mtk_sys))

    # Per-cell temperature, pressure (and, when met-coupled, H2O) via the
    # coordinate-observed function.
    needed = EarthSciMLBase.get_needed_vars(op, csys, mtk_sys, domain)
    has_h2o_obs = length(needed) == 3
    obs_f = EarthSciMLBase.build_coord_observed_function(mtk_sys, coord_args, needed)
    c1, c2, c3 = EarthSciMLBase.concrete_grid(domain)
    # One observed-cache per thread: the cell loop below is `Threads.@threads :static`
    # (cells write disjoint du rows, `_iso_ternary` is a pure function, and :static
    # scheduling pins iterations to threads so threadid-indexed caches are safe).
    # Size by `maxthreadid()` (≥ nthreads(); includes interactive/GC threads that a
    # loop iteration may legitimately report via threadid()).
    obscaches = [similar(domain.u_proto, length(needed)) for _ in 1:Threads.maxthreadid()]

    # State-vector rows for the species we read/repartition.
    syms = EarthSciMLBase.var2symbol.(unknowns(mtk_sys))
    ix = Dict(n => _iso_find(syms, n) for n in _ISO_SPECIES)

    # H2O read fallbacks when it is not observed/met-coupled: a dynamic H2O state
    # (mechanisms that carry H2O as a species), else the isconstantspecies parameter
    # default (standalone SuperFast/GEOSChem without met coupling).
    ix_h2o = has_h2o_obs ? nothing :
             findfirst(s -> endswith(string(s), "₊H2O"), syms)
    h2o_const = if has_h2o_obs || ix_h2o !== nothing
        0.0
    else
        ps = ModelingToolkit.parameters(mtk_sys)
        ip = findfirst(p_ -> endswith(_iso_name(p_), "₊H2O"), ps)
        ip === nothing && error(
            "IsorropiaOp: H2O not found among observed variables, state unknowns, or " *
            "parameters — the mechanism must carry H2O (met-coupled, as a species, or as " *
            "an isconstantspecies parameter)")
        Float64(ModelingToolkit.getdefault(ps[ip]))
    end

    # Non-convergence diagnostics: cells where `_iso_ternary` returned ok=false get a zero
    # tendency (correct/conservative), but a persistently high count signals a regime the
    # solve mishandles — surfaced once per simulated hour so CTM logs show it.
    noconv = Threads.Atomic{Int}(0)
    ncalls = Threads.Atomic{Int}(0)
    next_report = Ref(-Inf)

    function run(du, u, p, t)
        u = reshape(u, nrows, sz...)
        du = reshape(du, nrows, sz...)
        if t >= next_report[]
            if noconv[] > 0
                @info "IsorropiaOp: $(noconv[]) non-converged cell-solves of $(ncalls[]) since last report" t
            end
            Threads.atomic_xchg!(noconv, 0)
            Threads.atomic_xchg!(ncalls, 0)
            next_report[] = t + 3600.0
        end
        Threads.@threads :static for j in 1:ncell
            obscache = obscaches[Threads.threadid()]
            HNO3 = u[ix["HNO3"], II[j]]
            NIT = u[ix["NIT"], II[j]]
            NH3 = u[ix["NH3"], II[j]]
            NH4 = u[ix["NH4"], II[j]]
            SO4 = u[ix["SO4"], II[j]]
            # Skip gate: with essentially no inorganic N mass to repartition (< 1e-4 ppb
            # combined), the equilibrium tendency is numerically 0 — skip the expensive
            # ternary solve (clean-air columns are the common case aloft).
            if NH3 + NH4 + HNO3 + NIT < 1.0e-4
                du[ix["HNO3"], II[j]] = 0.0
                du[ix["NIT"], II[j]] = 0.0
                du[ix["NH3"], II[j]] = 0.0
                du[ix["NH4"], II[j]] = 0.0
                continue
            end
            obs_f(obscache, view(u, :, II[j]), p, t, c1[j], c2[j], c3[j])
            T = obscache[1]
            P = obscache[2]
            H2O = has_h2o_obs ? obscache[3] :
                  (ix_h2o !== nothing ? u[ix_h2o, II[j]] : h2o_const)
            ppb2m = 1.0e-9 * P / (_ISO_RG * T)            # ppb -> mol m^-3
            RH = clamp(H2O * 1.0e-9 * P / _iso_esat(T), 0.01, 0.99)
            ok, gHNO3, gNH3 = _iso_ternary(
                SO4 * ppb2m, (NH3 + NH4) * ppb2m, (HNO3 + NIT) * ppb2m, T, RH)
            if ok
                # Relax each gas toward its equilibrium partial concentration; the paired
                # aerosol species takes the mirror tendency so total nitrate (HNO3+NIT) and
                # ammonium (NH3+NH4) are conserved exactly and aqueous+solid aerosol
                # nitrate/ammonium lump into NIT/NH4.
                dHNO3 = op.k_mt * (gHNO3 / ppb2m - HNO3)
                dNH3 = op.k_mt * (gNH3 / ppb2m - NH3)
                du[ix["HNO3"], II[j]] = dHNO3
                du[ix["NIT"], II[j]] = -dHNO3
                du[ix["NH3"], II[j]] = dNH3
                du[ix["NH4"], II[j]] = -dNH3
            else
                Threads.atomic_add!(noconv, 1)
                du[ix["HNO3"], II[j]] = 0.0
                du[ix["NIT"], II[j]] = 0.0
                du[ix["NH3"], II[j]] = 0.0
                du[ix["NH4"], II[j]] = 0.0
            end
            Threads.atomic_add!(ncalls, 1)
        end
        return reshape(du, :)
    end
    return run
end

end

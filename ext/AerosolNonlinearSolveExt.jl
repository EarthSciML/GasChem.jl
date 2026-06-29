module AerosolNonlinearSolveExt

# Operator-split coupling of the GEOS-Chem gas phase to ISORROPIA II inorganic aerosol
# thermodynamic equilibrium. Provides the `EarthSciMLBase.Operator` methods for
# `GasChem.IsorropiaOp` (the struct lives in the main package).
#
# Rather than solving the full ISORROPIA system as one large nonlinear problem (which is
# multi-rooted and converges only in a narrow regime), this follows the *original*
# ISORROPIA II design: for the fine-mode NH4-SO4-NO3 system (the species the GEOS-Chem port
# feeds; sea-salt/dust handled separately) under the metastable (aqueous) assumption GEOS-Chem
# uses, the equilibrium reduces — given [H+] — to closed-form ion partitioning, so the whole
# solve is a single 1-D bisection on [H+] of the charge-balance residual (monotone in [H+],
# hence globally convergent) with an inner activity/water fixed-point. The thermodynamics
# (van't Hoff K, Kusik-Meissner activity coefficients, ZSR water) are reused from
# `Aerosol`'s ISORROPIA II implementation. (NB: this extension is keyed on the [Aerosol,
# NonlinearSolve] weak-dep pair for now but no longer needs NonlinearSolve; the solver could
# move into Aerosol as an exported API and the extension be rekeyed on Aerosol alone.)

using GasChem, Aerosol, ModelingToolkit, EarthSciMLBase

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

# State species the operator reads from / writes to.
const _ISO_SPECIES = ("HNO3", "NIT", "NH3", "NH4", "SO4", "H2O")

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
function _iso_ternary(TSO4, TNH, TNO3, T, RH; nouter = 40, nbisect = 60, tol = 1.0e-5)
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
        lo, hi = -14.0, 3.0
        if resid(10.0^lo) * resid(10.0^hi) > 0
            m_H = resid(10.0^lo) < 0 ? 10.0^hi : 10.0^lo
        else
            for _ in 1:nbisect
                mid = 0.5 * (lo + hi)
                resid(10.0^mid) < 0 ? (lo = mid) : (hi = mid)
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

function EarthSciMLBase.get_needed_vars(::GasChem.IsorropiaOp, csys, mtk_sys, domain::DomainInfo)
    # T and P are read per-cell through the coordinate-observed function; species are state
    # variables read directly from the state vector by index.
    [_iso_metvar(mtk_sys, "T"), _iso_metvar(mtk_sys, "P")]
end

function EarthSciMLBase.get_odefunction(
        op::GasChem.IsorropiaOp, csys, mtk_sys, coord_args, domain::DomainInfo, u0, p, alg)
    sz = tuple(size(domain)...)
    II = CartesianIndices(sz)
    ncell = prod(sz)
    nrows = length(unknowns(mtk_sys))

    # Per-cell temperature and pressure via the coordinate-observed function.
    obs_f = EarthSciMLBase.build_coord_observed_function(
        mtk_sys, coord_args, EarthSciMLBase.get_needed_vars(op, csys, mtk_sys, domain))
    c1, c2, c3 = EarthSciMLBase.concrete_grid(domain)
    obscache = similar(domain.u_proto, 2)

    # State-vector rows for the species we read/repartition.
    syms = EarthSciMLBase.var2symbol.(unknowns(mtk_sys))
    ix = Dict(n => _iso_find(syms, n) for n in _ISO_SPECIES)

    function run(du, u, p, t)
        u = reshape(u, nrows, sz...)
        du = reshape(du, nrows, sz...)
        for j in 1:ncell
            obs_f(obscache, view(u, :, II[j]), p, t, c1[j], c2[j], c3[j])
            T, P = obscache
            HNO3 = u[ix["HNO3"], II[j]]
            NIT = u[ix["NIT"], II[j]]
            NH3 = u[ix["NH3"], II[j]]
            NH4 = u[ix["NH4"], II[j]]
            SO4 = u[ix["SO4"], II[j]]
            H2O = u[ix["H2O"], II[j]]
            ppb2m = 1.0e-9 * P / (_ISO_RG * T)            # ppb -> mol m^-3
            RH = clamp(H2O * 1.0e-9 * P / _iso_esat(T), 0.01, 0.99)
            ok, gHNO3, gNH3 = _iso_ternary(
                (HNO3 + NIT) * ppb2m, (NH3 + NH4) * ppb2m, SO4 * ppb2m, T, RH)
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
                du[ix["HNO3"], II[j]] = 0.0
                du[ix["NIT"], II[j]] = 0.0
                du[ix["NH3"], II[j]] = 0.0
                du[ix["NH4"], II[j]] = 0.0
            end
        end
        return reshape(du, :)
    end
    return run
end

end

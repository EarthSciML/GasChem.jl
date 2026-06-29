module AerosolNonlinearSolveExt

# Operator-split coupling of the GEOS-Chem gas phase to the ISORROPIA II inorganic
# aerosol thermodynamic-equilibrium model (`Aerosol.IsorropiaEquilibrium`). Provides the
# `EarthSciMLBase.Operator` methods for `GasChem.IsorropiaOp` (the struct lives in the main
# package). Requires both Aerosol and NonlinearSolve, hence a dedicated extension.

using GasChem, Aerosol, NonlinearSolve, ModelingToolkit, EarthSciMLBase
using SymbolicIndexingInterface: setp, setu, getu
using NonlinearSolve: ReturnCode   # re-exported from SciMLBase

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

# Temperature/pressure variable to read per grid cell. In a CTM these are coupled from the
# meteorology (e.g. GEOSFP) and appear as observed variables; search the unknowns and the
# observed equations so the per-cell read works wherever the met provider puts them.
function _iso_metvar(mtk_sys, name)
    cands = vcat(unknowns(mtk_sys), [eq.lhs for eq in observed(mtk_sys)])
    i = findfirst(v -> endswith(_iso_name(v), "₊" * name), cands)
    i === nothing && error(
        "IsorropiaOp: meteorology variable ending in `₊$name` not found — couple a met " *
        "source (e.g. GEOSFP) so $name is available per grid cell")
    cands[i]
end

# State species the operator reads from / writes to (all GEOS-Chem state variables).
const _ISO_SPECIES = ("HNO3", "NIT", "NH3", "NH4", "SO4", "H2O")

function EarthSciMLBase.get_needed_vars(::GasChem.IsorropiaOp, csys, mtk_sys, domain::DomainInfo)
    # T and P are read per-cell through the coordinate-observed function; the species are
    # state variables read directly from the state vector by index.
    [_iso_metvar(mtk_sys, "T"), _iso_metvar(mtk_sys, "P")]
end

function EarthSciMLBase.get_odefunction(
        op::GasChem.IsorropiaOp, csys, mtk_sys, coord_args, domain::DomainInfo, u0, p, alg)
    sz = tuple(size(domain)...)
    II = CartesianIndices(sz)
    ncell = prod(sz)
    nrows = length(unknowns(mtk_sys))
    R_gas = 8.314

    # Per-cell temperature and pressure via the coordinate-observed function.
    obs_f = EarthSciMLBase.build_coord_observed_function(
        mtk_sys, coord_args, EarthSciMLBase.get_needed_vars(op, csys, mtk_sys, domain))
    c1, c2, c3 = EarthSciMLBase.concrete_grid(domain)
    obscache = similar(domain.u_proto, 2)

    # State-vector rows for the species we read/repartition.
    syms = EarthSciMLBase.var2symbol.(unknowns(mtk_sys))
    ix = Dict(n => _iso_find(syms, n) for n in _ISO_SPECIES)

    # ISORROPIA II equilibrium problem, built once and warm-started per cell.
    iso = mtkcompile(Aerosol.IsorropiaEquilibrium())
    set_tot! = setp(iso, [iso.W_NO3_total, iso.W_NH3_total, iso.W_SO4_total, iso.RH, iso.T_env])
    # Seed a few key aqueous ions from the totals each solve. The equilibrium is nonlinear
    # with multiple roots; a default guess lets Newton converge to a non-physical partition.
    # This minimal seed (only the ions below; the rest keep their package defaults) is the one
    # validated to converge for the ammonia-rich / humid regime — over-specifying the seed
    # (e.g. also pinning c_NH4/c_HSO4) was found to push the solve into a spurious basin.
    set_guess! = setu(iso, [iso.c_SO4, iso.c_NO3, iso.c_H, iso.c_Cl, iso.I_s])
    get_gas = getu(iso, [iso.g_HNO3, iso.g_NH3])
    init = Pair{Any, Float64}[]
    for u in unknowns(iso)
        g = ModelingToolkit.getguess(u)
        g !== nothing && push!(init, u => Float64(g))
    end
    baseprob = NonlinearProblem(iso, init)
    warm = [copy(baseprob.u0) for _ in 1:ncell]   # per-cell warm-start cache

    # Minimal initial guess (mol m^-3): most sulfate as SO4^2-, a little aerosol nitrate, a
    # weakly acidic solution, modest ionic strength. Strictly positive (the activity model
    # takes logs of concentrations, so exact zeros diverge).
    iso_guess(n_no3, n_nhx, n_so4) =
        [max(1.0e-20, 0.95 * n_so4), max(1.0e-20, 0.6 * n_no3), 1.0e-11, 1.0e-20, 10.0]

    # Solve the equilibrium for one cell. Returns `(valid, g_HNO3, g_NH3)`. The ISORROPIA II
    # nonlinear system does not converge for every composition/RH (it is most robust for the
    # ammonia-rich, humid regime), and even a "converged" solve can land on a non-physical
    # root. We therefore accept a solve only if it both converged and is physically valid
    # (each gas concentration within [0, total]); otherwise the caller leaves the cell
    # unchanged (du = 0) rather than injecting a spurious partition.
    function solve_cell!(j, n_no3, n_nhx, n_so4, RH, T)
        prob = remake(baseprob; u0 = warm[j])
        set_tot!(prob, [n_no3, n_nhx, n_so4, RH, T])
        set_guess!(prob, iso_guess(n_no3, n_nhx, n_so4))
        s = solve(prob, NewtonRaphson(); maxiters = 100)
        if s.retcode != ReturnCode.Success
            prob = remake(baseprob; u0 = baseprob.u0)
            set_tot!(prob, [n_no3, n_nhx, n_so4, RH, T])
            set_guess!(prob, iso_guess(n_no3, n_nhx, n_so4))
            s = solve(prob, RobustMultiNewton(); maxiters = 10000)
        end
        gHNO3, gNH3 = get_gas(s)   # (g_HNO3, g_NH3) [mol m^-3]
        # Accept only a converged, physically-sane root: each gas within [0, total], and — for
        # a sulfate-rich (acidic) aerosol, which binds ammonia strongly — gas-phase ammonia
        # must be a minority of total NHx (a spurious "all-gas NH3" root can otherwise satisfy
        # the simple [0, total] bound). Anything else ⇒ treat as non-converged (du = 0).
        valid = s.retcode == ReturnCode.Success &&
                -1.0e-20 <= gHNO3 <= n_no3 * (1 + 1.0e-6) &&
                -1.0e-20 <= gNH3 <= n_nhx * (1 + 1.0e-6) &&
                !(n_nhx < 2 * n_so4 && gNH3 > 0.5 * n_nhx)
        valid && (warm[j] = copy(s.u))
        return (valid, gHNO3, gNH3)
    end

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
            ppb2m = 1.0e-9 * P / (R_gas * T)            # ppb -> mol m^-3
            RH = clamp(H2O * 1.0e-9 * P / _iso_esat(T), 0.01, 0.99)
            valid, gHNO3, gNH3 = solve_cell!(
                j, (HNO3 + NIT) * ppb2m, (NH3 + NH4) * ppb2m, SO4 * ppb2m, RH, T)
            if valid
                # Relax each gas toward its equilibrium partial concentration; the paired
                # aerosol species takes the mirror tendency so the total (HNO3+NIT, NH3+NH4)
                # is conserved exactly and aqueous+solid aerosol nitrate/ammonium lump into
                # NIT/NH4.
                dHNO3 = op.k_mt * (gHNO3 / ppb2m - HNO3)
                dNH3 = op.k_mt * (gNH3 / ppb2m - NH3)
                du[ix["HNO3"], II[j]] = dHNO3
                du[ix["NIT"], II[j]] = -dHNO3
                du[ix["NH3"], II[j]] = dNH3
                du[ix["NH4"], II[j]] = -dNH3
            else
                # Equilibrium did not converge to a physical partition: leave the cell
                # unchanged this step rather than injecting a spurious tendency.
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

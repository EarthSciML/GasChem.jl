module AerosolNonlinearSolveExt

# Operator-split coupling of the GEOS-Chem gas phase to the ISORROPIA II inorganic
# aerosol thermodynamic-equilibrium model (`Aerosol.IsorropiaEquilibrium`). Provides the
# `EarthSciMLBase.Operator` methods for `GasChem.IsorropiaOp` (the struct lives in the main
# package). Requires both Aerosol and NonlinearSolve, hence a dedicated extension.

using GasChem, Aerosol, NonlinearSolve, ModelingToolkit, EarthSciMLBase
using SymbolicIndexingInterface: setp, getu, getp
using NonlinearSolve: ReturnCode   # re-exported from SciMLBase

# Saturation vapour pressure over liquid water [Pa] (Tetens, 1930), T in K.
_iso_esat(T) = 611.2 * exp(17.67 * (T - 273.15) / (T - 29.65))

# Index of the unknown/parameter whose flattened name ends in "₊name" (parent-agnostic;
# the leading "₊" guards against e.g. matching "HNO3" inside another species name).
function _iso_find(syms, name)
    i = findfirst(s -> endswith(string(s), "₊" * name), syms)
    i === nothing && error("IsorropiaOp: variable/parameter ending in `₊$name` not found in the coupled system")
    i
end

# Species the operator reads from / writes to (all GEOS-Chem state variables).
const _ISO_SPECIES = ("HNO3", "NIT", "NH3", "NH4", "SO4", "H2O")

function EarthSciMLBase.get_needed_vars(::GasChem.IsorropiaOp, csys, mtk_sys, domain::DomainInfo)
    us = unknowns(mtk_sys)
    syms = EarthSciMLBase.var2symbol.(us)
    [us[_iso_find(syms, n)] for n in _ISO_SPECIES]
end

function EarthSciMLBase.get_odefunction(
        op::GasChem.IsorropiaOp, csys, mtk_sys, coord_args, domain::DomainInfo, u0, p, alg)
    sz = tuple(size(domain)...)
    II = CartesianIndices(sz)
    ncell = prod(sz)
    nrows = length(unknowns(mtk_sys))
    R_gas = 8.314

    # State-vector rows for the species we read/repartition.
    syms = EarthSciMLBase.var2symbol.(unknowns(mtk_sys))
    ix = Dict(n => _iso_find(syms, n) for n in _ISO_SPECIES)

    # Temperature and pressure are read from the gas-phase parameters (scalar in a box; in a
    # CTM they are supplied per-cell from GEOSFP — extend to a coordinate-observed read there).
    ps = parameters(mtk_sys)
    getT = getp(mtk_sys, ps[findfirst(s -> endswith(string(s), "₊T"), ps)])
    getP = getp(mtk_sys, ps[findfirst(s -> endswith(string(s), "₊P"), ps)])

    # ISORROPIA II equilibrium problem, built once and warm-started per cell.
    iso = mtkcompile(Aerosol.IsorropiaEquilibrium())
    set_tot! = setp(iso, [iso.W_NO3_total, iso.W_NH3_total, iso.W_SO4_total, iso.RH, iso.T_env])
    get_gas = getu(iso, [iso.g_HNO3, iso.g_NH3])
    init = Pair{Any, Float64}[]
    for u in unknowns(iso)
        g = ModelingToolkit.getguess(u)
        g !== nothing && push!(init, u => Float64(g))
    end
    baseprob = NonlinearProblem(iso, init)
    warm = [copy(baseprob.u0) for _ in 1:ncell]   # per-cell warm-start cache

    # Solve the equilibrium for one cell; warm-started Newton with a robust fallback.
    function solve_cell!(j, n_no3, n_nhx, n_so4, RH, T)
        prob = remake(baseprob; u0 = warm[j])
        set_tot!(prob, [n_no3, n_nhx, n_so4, RH, T])
        s = solve(prob, NewtonRaphson(); maxiters = 100)
        if s.retcode != ReturnCode.Success
            prob = remake(baseprob; u0 = baseprob.u0)
            set_tot!(prob, [n_no3, n_nhx, n_so4, RH, T])
            s = solve(prob, RobustMultiNewton(); maxiters = 10000)
        end
        s.retcode == ReturnCode.Success && (warm[j] = copy(s.u))
        get_gas(s)   # (g_HNO3, g_NH3) [mol m^-3]
    end

    function run(du, u, p, t)
        u = reshape(u, nrows, sz...)
        du = reshape(du, nrows, sz...)
        T = getT(p)
        P = getP(p)
        ppb2m = 1.0e-9 * P / (R_gas * T)            # ppb -> mol m^-3
        for j in 1:ncell
            HNO3 = u[ix["HNO3"], II[j]]
            NIT = u[ix["NIT"], II[j]]
            NH3 = u[ix["NH3"], II[j]]
            NH4 = u[ix["NH4"], II[j]]
            SO4 = u[ix["SO4"], II[j]]
            H2O = u[ix["H2O"], II[j]]
            RH = clamp(H2O * 1.0e-9 * P / _iso_esat(T), 0.01, 0.99)
            gHNO3, gNH3 = solve_cell!(
                j, (HNO3 + NIT) * ppb2m, (NH3 + NH4) * ppb2m, SO4 * ppb2m, RH, T)
            # Relax each gas toward its equilibrium partial concentration; the paired aerosol
            # species takes the mirror tendency so the total (HNO3+NIT, NH3+NH4) is conserved
            # exactly and aqueous+solid aerosol nitrate/ammonium are lumped into NIT/NH4.
            dHNO3 = op.k_mt * (gHNO3 / ppb2m - HNO3)
            dNH3 = op.k_mt * (gNH3 / ppb2m - NH3)
            du[ix["HNO3"], II[j]] = dHNO3
            du[ix["NIT"], II[j]] = -dHNO3
            du[ix["NH3"], II[j]] = dNH3
            du[ix["NH4"], II[j]] = -dNH3
        end
        return reshape(du, :)
    end
    return run
end

end

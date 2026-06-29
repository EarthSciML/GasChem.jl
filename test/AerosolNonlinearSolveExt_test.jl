# Tests for the GasChem↔Aerosol ISORROPIA operator (AerosolNonlinearSolveExt):
#   - operator-split coupling of GEOS-Chem to Aerosol.IsorropiaEquilibrium (ISORROPIA II)
#   - repartitions total nitrate (HNO3 ⇌ NIT) and ammonium (NH3 ⇌ NH4) toward equilibrium
#
# These require Aerosol (IsorropiaEquilibrium) + NonlinearSolve, so they exercise the
# AerosolNonlinearSolveExt weak-dependency extension; run locally against the dev'd Aerosol
# fork until an Aerosol release with IsorropiaEquilibrium is available.

@testsnippet IsorropiaOpSetup begin
    using GasChem, Aerosol, NonlinearSolve, EarthSciMLBase, ModelingToolkit, DomainSets
    using ModelingToolkit: t_nounits as t, D_nounits as D

    # Minimal inorganic-gas host exposing the species the operator reads. The parent name is
    # irrelevant: the operator matches species by their `₊<name>` suffix, so it works with any
    # gas mechanism (here a stand-in; in production it is the full GEOSChemGasPhase).
    function InorgGasHost(; name = :gas)
        @variables(HNO3(t)=1.0, NIT(t)=0.3, NH3(t)=1.0, NH4(t)=0.2, SO4(t)=2.0, H2O(t)=1.84e7)
        @parameters T=298.0 P=101325.0 lon=0.0 lat=0.0 lev=1.0
        System(
            [D(HNO3) ~ 0.0, D(NIT) ~ 0.0, D(NH3) ~ 0.0, D(NH4) ~ 0.0, D(SO4) ~ 0.0,
                D(H2O) ~ 1.0e-30 * (lon + lat + lev + T + P)],   # keep coord+T/P params alive
            t; name)
    end
    function small_domain()
        @parameters lon lat lev
        indep = t ∈ Interval(0.0, 3600.0)
        pdoms = [lon ∈ Interval(-1.0, 1.0), lat ∈ Interval(-1.0, 1.0), lev ∈ Interval(1.0, 1.0)]
        DomainInfo(constIC(0.0, indep), constBC(0.0, pdoms...); grid_spacing = [1.0, 1.0, 1.0])
    end
end

@testitem "ISORROPIA op: extension loads and registers operator methods" setup = [IsorropiaOpSetup] begin
    @test Base.get_extension(GasChem, :AerosolNonlinearSolveExt) !== nothing
    @test hasmethod(
        EarthSciMLBase.get_odefunction,
        Tuple{IsorropiaOp, Any, Any, Any, DomainInfo, Any, Any, Any})
    @test IsorropiaOp() isa EarthSciMLBase.Operator
    @test IsorropiaOp(k_mt = 0.5).k_mt == 0.5
end

@testitem "ISORROPIA op: GEOS-Chem exposes the required species" begin
    using GasChem, ModelingToolkit
    us = string.(unknowns(convert(System, GEOSChemGasPhase())))
    for n in ("HNO3", "NIT", "NH3", "NH4", "SO4", "H2O")
        @test "$n(t)" in us
    end
end

@testitem "ISORROPIA op: runs, conserves totals, partitions correctly" setup = [IsorropiaOpSetup] begin
    domain = small_domain()
    csys = couple(InorgGasHost(), IsorropiaOp(), domain)
    mtk = convert(System, csys)
    sc, ca = EarthSciMLBase._prepare_coord_sys(mtk, domain)
    pp = EarthSciMLBase.default_params(sc)
    u = EarthSciMLBase.init_u(sc, domain)

    opf = EarthSciMLBase.nonstiff_ops(csys, sc, ca, domain, reshape(u, :), pp, MapBroadcast())
    du = zero(reshape(u, :))
    opf(du, reshape(u, :), pp, 0.0)
    duM = reshape(du, length(unknowns(sc)), size(domain)...)

    sy = EarthSciMLBase.var2symbol.(unknowns(sc))
    row(n) = findfirst(s -> endswith(string(s), "₊" * n), sy)
    dHNO3 = duM[row("HNO3"), 1, 1, 1]
    dNIT = duM[row("NIT"), 1, 1, 1]
    dNH3 = duM[row("NH3"), 1, 1, 1]
    dNH4 = duM[row("NH4"), 1, 1, 1]

    # total nitrate (HNO3+NIT) and total ammonium (NH3+NH4) are conserved exactly
    @test dHNO3 + dNIT ≈ 0 atol = 1.0e-12
    @test dNH3 + dNH4 ≈ 0 atol = 1.0e-12
    # the operator is doing something (nonzero mass transfer)
    @test abs(dNH3) > 0
    # sulfate-rich with ammonia present: ammonia condenses onto the acidic aerosol
    @test dNH3 < 0
    @test dNH4 > 0
    # sulfate is non-volatile: no ISORROPIA tendency on SO4
    @test duM[row("SO4"), 1, 1, 1] == 0
end

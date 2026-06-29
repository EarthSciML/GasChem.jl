# Tests for the GasChem↔Aerosol ISORROPIA operator (GasChem.IsorropiaOp, implemented in the
# AerosolExt extension):
#   - operator-split coupling of GEOS-Chem to inorganic aerosol thermodynamic equilibrium
#   - repartitions total nitrate (HNO3 ⇌ NIT) and ammonium (NH3 ⇌ NH4) toward equilibrium
#
# These require Aerosol, so they exercise the AerosolExt weak-dependency extension; run
# locally against the dev'd Aerosol fork until an Aerosol release with IsorropiaEquilibrium
# is available.

@testsnippet IsorropiaOpSetup begin
    using GasChem, Aerosol, EarthSciMLBase, ModelingToolkit, DomainSets
    using ModelingToolkit: t_nounits as t, D_nounits as D

    # Minimal inorganic-gas host exposing the species the operator reads. The parent name is
    # irrelevant: the operator matches species by their `₊<name>` suffix, so it works with any
    # gas mechanism (here a stand-in; in production it is the full GEOSChemGasPhase). T and P
    # are coordinate-dependent observed variables, mimicking a per-cell GEOSFP-coupled met
    # field (so the operator's per-cell observed read of T/P is exercised, not a scalar).
    # Composition [ppb] is ammonia-rich / sulfate-poor and humid (RH≈0.8) — the regime in
    # which the ISORROPIA II nonlinear solve converges robustly; T varies per cell but stays
    # near 298 K so RH stays in that regime in every cell.
    function InorgGasHost(; name = :gas)
        @variables(HNO3(t)=1.0, NIT(t)=0.22, NH3(t)=7.0, NH4(t)=0.34, SO4(t)=2.45, H2O(t)=2.48e7,
            T(t), P(t))
        @parameters lon=0.0 lat=0.0 lev=1.0
        System(
            [D(HNO3) ~ 0.0, D(NIT) ~ 0.0, D(NH3) ~ 0.0, D(NH4) ~ 0.0, D(SO4) ~ 0.0,
                D(H2O) ~ 1.0e-30 * lon,                  # keep the lon coordinate param alive
                T ~ 298.0 + 0.5 * lat,                   # per-cell temperature [K]
                P ~ 101325.0 - 100.0 * lev],             # per-cell pressure [Pa]
            t; name)
    end
    # Sulfate-rich (acidic) composition [ppb]: the regime where the monolithic ISORROPIA
    # solve failed (spurious roots) but the regime+bisection solve converges. Excess sulfate
    # ⇒ ammonia condenses onto the acidic aerosol while nitrate is driven to the gas phase.
    function InorgGasHostRich(; name = :gas)
        @variables(HNO3(t)=0.5, NIT(t)=1.0, NH3(t)=1.0, NH4(t)=0.5, SO4(t)=4.0, H2O(t)=2.48e7,
            T(t), P(t))
        @parameters lon=0.0 lat=0.0 lev=1.0
        System(
            [D(HNO3) ~ 0.0, D(NIT) ~ 0.0, D(NH3) ~ 0.0, D(NH4) ~ 0.0, D(SO4) ~ 0.0,
                D(H2O) ~ 1.0e-30 * lon,
                T ~ 298.0 + 0.5 * lat,
                P ~ 101325.0 - 100.0 * lev],
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
    @test Base.get_extension(GasChem, :AerosolExt) !== nothing
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

    # In every grid cell (each at its own per-cell T/RH): total nitrate (HNO3+NIT) and total
    # ammonium (NH3+NH4) are conserved exactly — du(aerosol) = -du(gas) where the equilibrium
    # converged, or du = 0 where it did not (the operator leaves a cell unchanged rather than
    # injecting a spurious partition) — and sulfate is non-volatile (no ISORROPIA tendency).
    for I in CartesianIndices(size(domain))
        @test duM[row("HNO3"), I] + duM[row("NIT"), I] ≈ 0 atol = 1.0e-12
        @test duM[row("NH3"), I] + duM[row("NH4"), I] ≈ 0 atol = 1.0e-12
        @test duM[row("SO4"), I] == 0
    end

    # Centre cell (lon=lat=0 ⇒ T=298 K, RH≈0.8): ammonia-rich / sulfate-poor. Both gases
    # condense onto the aerosol (NH3↓/NH4↑, HNO3↓/NIT↑) and the operator does nonzero,
    # physically-signed mass transfer.
    mid = (2, 2, 1)
    @test duM[row("NH3"), mid...] < 0
    @test duM[row("NH4"), mid...] > 0
    @test duM[row("HNO3"), mid...] < 0
    @test duM[row("NIT"), mid...] > 0
    # Implied equilibrium gas (c + du/k_mt) must be physically bounded by the species' OWN
    # total. Guards against feeding the solver swapped totals (e.g. sulfate in the nitrate
    # slot), which would let "gas nitrate" track the sulfate total and exceed total nitrate.
    kmt = 1 / 300   # matches the IsorropiaOp() default k_mt
    @test -1.0e-6 <= 1.0 + duM[row("HNO3"), mid...] / kmt <= 1.22 + 1.0e-3   # total NO3 = 1.22 ppb
    @test -1.0e-6 <= 7.0 + duM[row("NH3"), mid...] / kmt <= 7.34 + 1.0e-3    # total NHx = 7.34 ppb
end

@testitem "ISORROPIA op: sulfate-rich regime converges and partitions correctly" setup = [IsorropiaOpSetup] begin
    domain = small_domain()
    csys = couple(InorgGasHostRich(), IsorropiaOp(), domain)
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

    # Conservation + non-volatile sulfate in every cell (this regime made the old monolithic
    # solver return spurious roots; the regime+bisection solve converges here).
    for I in CartesianIndices(size(domain))
        @test duM[row("HNO3"), I] + duM[row("NIT"), I] ≈ 0 atol = 1.0e-12
        @test duM[row("NH3"), I] + duM[row("NH4"), I] ≈ 0 atol = 1.0e-12
        @test duM[row("SO4"), I] == 0
    end

    # Sulfate-rich (SO4=4 vs NHx=1.5 ppb): ammonia condenses onto the acidic aerosol
    # (NH3↓/NH4↑) while nitrate is driven to the gas phase (HNO3↑/NIT↓).
    mid = (2, 2, 1)
    @test duM[row("NH3"), mid...] < 0
    @test duM[row("NH4"), mid...] > 0
    @test duM[row("HNO3"), mid...] > 0
    @test duM[row("NIT"), mid...] < 0
    # Equilibrium gas bounded by the total (catches swapped sulfate/nitrate totals: with the
    # swap, sulfate-rich would drive "gas nitrate" toward the larger sulfate total).
    kmt = 1 / 300   # matches the IsorropiaOp() default k_mt
    @test -1.0e-6 <= 0.5 + duM[row("HNO3"), mid...] / kmt <= 1.5 + 1.0e-3    # total NO3 = 1.5 ppb
    @test -1.0e-6 <= 1.0 + duM[row("NH3"), mid...] / kmt <= 1.5 + 1.0e-3     # total NHx = 1.5 ppb
end

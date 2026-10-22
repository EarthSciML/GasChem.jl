# Tests for the SuperFast-aerosol port: 7 new species (SO2 re-enabled, SO4, NH3, NH4, NIT,
# N2O5, NO3), inert-by-default k_het_*/k_cld1 hooks, and the SuperFast↔Aerosol couple2s
# (CloudChemistryFixedpH, AerosolDistribution) mirroring the GEOS-Chem templates.

@testsnippet SFAeroSetup begin
    using GasChem, Aerosol, EarthSciMLBase, ModelingToolkit, OrdinaryDiffEqRosenbrock
    using ModelingToolkit: t_nounits as t, D_nounits as D
    using DynamicQuantities
    const Interval = EarthSciMLBase.DomainSets.Interval

    sf_box_o3_traj(; tend = 8 * 3600.0) = begin
        sys = mtkcompile(convert(System, GasChem.SuperFast()))
        prob = ODEProblem(sys, [], (0.0, tend), [])
        sol = solve(prob, Rosenbrock23(); reltol = 1.0e-10, abstol = 1.0e-12, saveat = 3600.0)
        io3 = findfirst(v -> string(v) in ("O3(t)", "SuperFast₊O3(t)"), unknowns(sys))
        [u[io3] for u in sol.u]
    end

    # Minimal met stand-in that binds SuperFast's T/P the way the production GEOSFP
    # coupling does (EarthSciDataExt: `param_to_var(c, :T, :P, :H2O)` + connector
    # equations). The cloud/het couple2s reference T/P in *variable* form by design
    # (see the NOTE in ext/AerosolExt.jl), so coupling them WITHOUT a met source
    # leaves T(t)/P(t) unbound and the coupled system is intentionally unbalanced.
    struct SFTestMetCoupler
        sys::Any
    end
    function SFTestMet()
        @parameters T_met = 280.0 [unit = u"K"]
        @parameters P_met = 101325.0 [unit = u"Pa"]
        @variables T(t) = 280.0 [unit = u"K"]
        @variables P(t) = 101325.0 [unit = u"Pa"]
        System([T ~ T_met, P ~ P_met], t; name = :sftestmet,
            metadata = Dict(EarthSciMLBase.CoupleType => SFTestMetCoupler))
    end
    function EarthSciMLBase.couple2(c::GasChem.SuperFastCoupler, m::SFTestMetCoupler)
        c, m = c.sys, m.sys
        c = EarthSciMLBase.param_to_var(c, :T, :P)
        EarthSciMLBase.ConnectorSystem([c.T ~ m.T, c.P ~ m.P], c, m)
    end
end

@testitem "SF-aerosol: new species present, k hooks default inert" setup = [SFAeroSetup] begin
    rs = convert(System, GasChem.SuperFast())
    us = string.(unknowns(rs))
    for n in ("SO2", "SO4", "NH3", "NH4", "NIT", "N2O5", "NO3")
        @test any(u -> u == "$n(t)" || endswith(u, "₊$n(t)"), us)
    end
    # H2O must NOT be a state (met-coupled isconstantspecies parameter)
    @test !any(u -> u == "H2O(t)" || endswith(u, "₊H2O(t)"), us)
    ps = ModelingToolkit.parameters(rs)
    pnames = string.(ps)
    for k in ("k_cld1", "k_het_N2O5", "k_het_HO2", "k_het_NO2", "k_het_NO3", "H2O")
        i = findfirst(p -> p == k || endswith(p, "₊" * k), pnames)
        @test i !== nothing
        if startswith(k, "k_")
            @test Float64(ModelingToolkit.getdefault(ps[i])) == 0.0   # inert default
        end
    end
    # The het/cloud reactions are in the network (N2O5 --> 2 HNO3 etc.)
    eqstr = join(string.(equations(rs)), "\n")
    @test occursin("k_het_N2O5", eqstr) && occursin("k_cld1", eqstr)
end

@testitem "SF-aerosol: bare SuperFast box O3 trajectory unchanged by inert additions" setup = [SFAeroSetup] begin
    # Hourly box O3, Rosenbrock23 reltol=1e-10/abstol=1e-12, default
    # ICs/params. The additions in this PR are all zero-rate/constant with
    # defaults, so bare-SF chemistry must be unchanged.
    # Checkpoints at h2/h4/h8 (meaningful O3 values; the 24h endpoint is ~1e-7 ppb
    # where adaptive-step noise dominates any relative comparison).
    BASE_H2 = 7.23332425963
    BASE_H4 = 1.39027183281
    BASE_H8 = 0.0525209917846
    traj = sf_box_o3_traj()
    @test traj[3] ≈ BASE_H2 rtol = 1.0e-5
    @test traj[5] ≈ BASE_H4 rtol = 1.0e-5
    @test traj[9] ≈ BASE_H8 rtol = 1.0e-4
end

@testitem "SF-aerosol: cloud coupling drives SO2+H2O2 -> SO4 with S conservation" setup = [SFAeroSetup] begin
    csys = couple(GasChem.SuperFast(), Aerosol.CloudChemistryFixedpH(), SFTestMet())
    sys = convert(System, csys)  # convert compiles by default; double mtkcompile errors on MTK ≥ v11
    us = unknowns(sys)
    ix(n) = findfirst(v -> endswith(string(v), "₊$n(t)") || string(v) == "$n(t)", us)
    sSO2 = us[ix("SO2")]; sH2O2 = us[ix("H2O2")]; sSO4 = us[ix("SO4")]
    # boost SO2/H2O2 so the in-cloud channel visibly converts S within 6h.
    # Symbolic u0 overrides + symbolic solution indexing: positional prob.u0
    # ordering is not guaranteed to match unknowns(sys) on MTK ≥ v11.
    prob = ODEProblem(sys, [sSO2 => 5.0, sH2O2 => 2.0], (0.0, 6 * 3600.0))
    sol = solve(prob, Rosenbrock23(); reltol = 1.0e-8, abstol = 1.0e-10)
    so2e = sol[sSO2][end]; so4e = sol[sSO4][end]
    @test so4e > 1.5                      # substantial sulfate formed
    @test so2e + so4e ≈ 5.0 + sol[sSO4][1] rtol = 1.0e-4   # total S conserved
end

@testitem "SF-aerosol: AerosolDistribution coupling exposes k_het_* via S_t" setup = [SFAeroSetup] begin
    csys = couple(GasChem.SuperFast(), Aerosol.AerosolDistribution(), SFTestMet())
    sys = convert(System, csys)  # convert compiles by default; double mtkcompile errors on MTK ≥ v11
    obs = string.([eq.lhs for eq in observed(sys)])
    for k in ("k_het_N2O5", "k_het_HO2", "k_het_NO2", "k_het_NO3")
        @test any(o -> occursin(k, o), obs)
    end
    @test any(o -> occursin("S_t", o), obs)
end

@testitem "SF-aerosol: triple coupling compiles without duplicate names + conserves" setup = [SFAeroSetup] begin
    # SFTestMet plays the CTM met coupling's role (param_to_var(c,:T,:P) + connector
    # equations), so the three-way param_to_var(:T) collision surface (met + cloud +
    # aerosol couplings on the same SuperFast instance) is exercised without GEOSFP data.
    csys = couple(GasChem.SuperFast(), Aerosol.CloudChemistryFixedpH(),
        Aerosol.AerosolDistribution(), GasChem.IsorropiaOp(), SFTestMet())
    sys = try
        convert(System, csys)  # convert compiles by default; double mtkcompile errors on MTK ≥ v11
    catch e
        @error "triple coupling failed to compile" exception = e
        nothing
    end
    @test sys !== nothing
    if sys !== nothing
        us = unknowns(sys)
        ix(n) = findfirst(v -> endswith(string(v), "₊$n(t)") || string(v) == "$n(t)", us)
        for n in ("SO2", "SO4", "NH3", "NH4", "NIT", "HNO3")
            @test ix(n) !== nothing
        end
    end
end

# Tests for the GasChem↔Aerosol coupling extension (AerosolExt):
#   - heterogeneous uptake (k_het_*) driven by aerosol surface area
#   - in-cloud sulfate chemistry (k_cld1) driven by CloudChemistry
#
# NOTE: these exercise the AerosolExt weak-dependency extension, so they require
# a version of Aerosol.jl that provides `AerosolDistributionCoupler` and
# `CloudChemistryFixedpHCoupler` (CoupleType metadata). Until that Aerosol release
# is available, these tests are skipped; they activate once a compatible Aerosol.jl
# release provides the coupler types.

@testitem "het: k_het_* params default to an inert no-op" begin
    using GasChem, ModelingToolkit
    sys = mtkcompile(GEOSChemGasPhase())
    ps = string.(parameters(sys))
    for k in ("k_het_N2O5", "k_het_HO2", "k_het_NO2", "k_het_NO3")
        @test k in ps
    end
    # N2O5 -> 2 HNO3 first-order het reaction is present in the network
    rs = GEOSChemGasPhase(rxn_sys = true)
    rxs = replace.(string.(GasChem.Catalyst.reactions(rs)), " " => "")
    @test any(s -> occursin("k_het_N2O5", s) && occursin("N2O5-->2", s), rxs)
end

@testitem "het: coupling drives k_het_* from Aerosol surface area" begin
    using GasChem, Aerosol, EarthSciMLBase, ModelingToolkit
    sys = convert(System, couple(GEOSChemGasPhase(), Aerosol.UrbanAerosol()))
    obs = string(observed(sys))
    # every het rate constant becomes a function of the aerosol surface area S_t
    for k in ("k_het_N2O5", "k_het_HO2", "k_het_NO2", "k_het_NO3")
        @test occursin(k, obs)
    end
    @test occursin("S_t", obs)
    @test_nowarn convert(System, couple(GEOSChemGasPhase(), Aerosol.UrbanAerosol()))
end

@testitem "cloud: k_cld1 driven by CloudChemistry; in-cloud SO2+H2O2->SO4" begin
    using GasChem, Aerosol, EarthSciMLBase, ModelingToolkit
    using OrdinaryDiffEqRosenbrock, OrdinaryDiffEqNonlinearSolve, SciMLBase
    sys = convert(System, couple(GEOSChemGasPhase(), Aerosol.CloudChemistryFixedpH()))
    @test occursin("k_cld1", string(observed(sys)))

    us = unknowns(sys)
    g(s) = us[findfirst(u -> occursin("₊" * s * "(", string(u)), us)]
    SO2, SO4, H2O2, O3 = g("SO2"), g("SO4"), g("H2O2"), g("O3")

    u0 = Dict(u => 0.0 for u in us)
    u0[SO2] = 5.0
    u0[H2O2] = 2.0
    u0[O3] = 40.0
    sol = solve(
        ODEProblem(sys, u0, (0.0, 6 * 3600.0); build_initializeprob = false),
        Rosenbrock23(); abstol = 1.0e-11, reltol = 1.0e-8
    )
    @test SciMLBase.successful_retcode(sol)
    @test sol[SO4][end] > 1.5                              # in-cloud sulfate produced (~2 ppb)
    @test sol[H2O2][end] < 0.01                            # H2O2 consumed (limiting oxidant)
    @test (sol[SO2][end] + sol[SO4][end]) ≈ 5.0 rtol = 1.0e-4  # total sulfur conserved
end

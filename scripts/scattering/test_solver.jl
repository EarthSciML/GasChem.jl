# Validation tests for the ported Cloud-J scattering solver.
# Run: julia --project=. scripts/scattering/test_solver.jl

include(joinpath(@__DIR__, "cloudj_solver.jl"))
using Test

@testset "energy conservation: conservative atmosphere + mirror surface" begin
    # Pure Rayleigh scattering (ss-albedo = 1) with albedo = 1: nothing absorbs,
    # so ALL incident solar flux (= u0, normalized) must exit the top as
    # diffuse flux FJTOP (plus any surviving direct beam is impossible - the
    # only exits are TOA diffuse). Solver discretization tolerance ~1%.
    nlay = 40
    ND = 2 * nlay + 1
    for (u0, tau_tot) in ((0.9, 1.0), (0.5, 2.0), (0.3, 0.5))
        # uniform grid in cumulative tau, TOA at LZ=1
        ztau = collect(range(0.0, tau_tot, length = ND))
        fz = [u0 > 0 ? exp(-t / u0) : 0.0 for t in ztau]
        pom = zeros(M2_, ND)
        pom[1, :] .= 1.0      # w0 = 1
        pom[3, :] .= 0.5      # Rayleigh P2 moment
        fsbot = fz[end] * u0
        rfl = fill(1.0, 5)    # mirror surface
        FJ, FJTOP, FJBOT, FIBOT = miesct(pom, fz, ztau, fsbot, rfl, u0, ND)
        @test isapprox(FJTOP, u0; rtol = 0.02)
    end
end

@testset "no scattering: diffuse field vanishes, J = direct beam" begin
    nlay = 30
    ND = 2 * nlay + 1
    ztau = collect(range(0.0, 3.0, length = ND))
    u0 = 0.7
    fz = [exp(-t / u0) for t in ztau]
    pom = zeros(M2_, ND)              # w0 = 0: pure absorption
    rfl = fill(0.0, 5)
    FJ, FJTOP, FJBOT, _ = miesct(pom, fz, ztau, fz[end] * u0, rfl, u0, ND)
    @test maximum(abs.(FJ)) < 1.0e-12
    @test abs(FJTOP) < 1.0e-12
end

@testset "optically thin limit: J -> direct beam ~ 1" begin
    J, _ = solve_band(0.9, 18)         # band 18: tiny Rayleigh + O3 OD
    @test all(isfinite, J)
    @test all(J .>= 0.0)
    @test 0.95 < J[end] < 1.3          # TOA: direct 1 + small upwelling
end

@testset "real column: diffuse enhancement bounded and positive" begin
    for band in (11, 12, 15, 17, 18), u0 in (0.9, 0.5, 0.1)
        J, d = solve_band(u0, band)
        @test all(isfinite, J)
        @test all(J .>= -1.0e-12)
        # J can exceed the unattenuated beam (backscatter+albedo) but not wildly
        @test maximum(J) < 3.0
        @test d.FJTOP >= -1.0e-12
    end
end

@testset "direct component matches calc_direct_flux exactly" begin
    # solve_band's FZ reuses sphere2J; at the surface edge, FZ (normalized)
    # must equal calc_direct_flux / top_flux at the surface pressure.
    for band in (12, 17), u0 in (0.9, 0.42255961917649837)
        ftau = solar_beam_edges(u0, band)
        direct = GasChem.calc_direct_flux(u0, GasChem.P_levels[1], band) /
            GasChem.top_flux[band]
        @test isapprox(ftau[1], direct; rtol = 1.0e-10)
    end
end

@testset "twilight (u0 < 0): finite, non-negative, dark surface UV" begin
    for band in (12, 17)
        J, d = solve_band(-0.1, band)
        @test all(isfinite, J)
        @test all(J .>= -1.0e-12)
        @test d.FSBOT == 0.0
    end
end

println("ALL SOLVER TESTS COMPLETE")

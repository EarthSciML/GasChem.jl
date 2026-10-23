# Regenerate src/tropospheric_interpolation_data.bson with multiple scattering
# from the ported Cloud-J Feautrier solver (see cloudj_solver.jl / README.md).
#
# Axes are taken from the EXISTING table so the result is a drop-in
# replacement: 18 bands x (tropospheric_P x cosSZA_vals), absolute actinic
# flux (top_flux-scaled). Writes to this directory; copy to src/ after review.
#
# Run: julia --project=. scripts/scattering/precompute_table.jl

include(joinpath(@__DIR__, "cloudj_solver.jl"))
using BSON
using Interpolations: interpolate, Gridded, Linear, extrapolate, Flat

const ALBEDO = 0.1

# --- axes from the existing table ------------------------------------------
BSON.@load joinpath(@__DIR__, "..", "..", "src", "tropospheric_interpolation_data.bson") Z_all tropospheric_P cosSZA_vals
const Z_old = Z_all
const P_nodes = collect(tropospheric_P)      # Pa, ascending 1000..100000
const U0_nodes = collect(cosSZA_vals)        # -0.2 .. 1.0

# Edge pressures of the solver column (Pa): edge 1 = surface .. 74 = TOA
# (the 74th edge is the ZZHT top-of-atmosphere pad; assign it P = 0).
const P_edges = vcat(collect(GasChem.P_levels), 0.0)

println(
    "Precomputing $(length(U0_nodes)) SZA nodes x 18 bands x ",
    length(P_nodes), " pressure nodes (albedo = $ALBEDO, AOD = 0) ..."
)

# Two products:
#   Z_diffuse - the DIFFUSE-ONLY field 4*FJ (what a box model cannot compute
#     locally; read at runtime by the online FastJX and added to its own
#     online direct-beam calculation).
#   Z_new     - the TOTAL flux table for FastJX_interpolation*: the offline
#     precompute of the online FastJX, i.e. calc_direct_flux evaluated at
#     each grid node plus the diffuse field at that node. This keeps
#     "interpolated = precompute of the online model" exact at the nodes.
Z_diffuse = [zeros(length(P_nodes), length(U0_nodes)) for _ in 1:18]
Z_new = [zeros(length(P_nodes), length(U0_nodes)) for _ in 1:18]

t0 = time()
for (ju, u0) in enumerate(U0_nodes)
    ftau_all = solar_beam_edges_allbands(u0)
    for band in 1:18
        ftau = ftau_all[:, band]
        J_edges, _ = solve_band(u0, band; albedo = ALBEDO, ftau = ftau)
        # diffuse-only component at the solver's edge pressures (absolute)
        Dabs = (J_edges .- ftau) .* GasChem.top_flux[band]
        # interpolate onto the table's pressure axis (Gridded needs
        # ascending knots; P_edges is descending in index)
        itp = extrapolate(
            interpolate((reverse(P_edges),), reverse(Dabs), Gridded(Linear())),
            Flat()
        )
        for (jp, P) in enumerate(P_nodes)
            D = max(itp(P), 0.0)
            Z_diffuse[band][jp, ju] = D
            Z_new[band][jp, ju] = GasChem.calc_direct_flux(u0, P, band) + D
        end
    end
end
println("done in ", round(time() - t0, digits = 1), " s")

# --- summary vs the table currently in src/ --------------------------------
println("\nnew/current table | band | surface (u0=1) | 500 hPa (u0=0.5) | max ratio")
for band in 1:18
    iP_s = length(P_nodes)               # 100000 Pa
    iP_m = findfirst(==(50500.0), P_nodes)
    iU_1 = length(U0_nodes)              # u0 = 1.0
    iU_h = findfirst(==(0.5), U0_nodes)
    olds = Z_old[band][iP_s, iU_1]
    news = Z_new[band][iP_s, iU_1]
    oldm = Z_old[band][iP_m, iU_h]
    newm = Z_new[band][iP_m, iU_h]
    r = filter(isfinite, vec(Z_new[band] ./ max.(Z_old[band], 1.0e-30)))
    println(
        rpad(band, 5), "| ",
        rpad(olds > 0 ? round(news / olds, digits = 3) : "-", 19), "| ",
        rpad(oldm > 0 ? round(newm / oldm, digits = 3) : "-", 21), "| ",
        round(
            maximum(
                r[vec(Z_old[band]) .> 1.0e-30 * maximum(Z_old[band])];
                init = 0.0
            ), digits = 2
        )
    )
end

# --- write -----------------------------------------------------------------
provenance = Dict(
    "generator" => "scripts/scattering/precompute_table.jl",
    "solver" => "Cloud-J OPMIE/MIESCT/BLKSLV port (geoschem/cloud-j@16d18b07, GEOS-Chem 14.7.1 pin)",
    "atmosphere" => "fixed climatological column: Rayleigh + O2/O3 absorption, 313.4 DU O3, T_profile clim.",
    "albedo" => ALBEDO,
    "aerosol_od" => 0.0,
    "solf" => "not included (applied at runtime by solar_flux_factor)",
)

# diffuse-only field (consumed online by FastJX)
Z_all = Z_diffuse
out_d = joinpath(@__DIR__, "diffuse_flux_data.bson")
BSON.@save out_d Z_all tropospheric_P cosSZA_vals provenance
println("\nwrote ", out_d)

# total-flux table (consumed by FastJX_interpolation*): precompute of the
# online FastJX (calc_direct_flux + diffuse)
Z_all = Z_new
out = joinpath(@__DIR__, "tropospheric_interpolation_data_scattering.bson")
BSON.@save out Z_all tropospheric_P cosSZA_vals provenance
println("wrote ", out)

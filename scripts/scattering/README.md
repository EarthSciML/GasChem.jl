# Clear-sky scattering for Fast-JX: first-principles precompute

Adds Rayleigh multiple scattering + surface reflection to this package's
Fast-JX photolysis from first principles: a Julia port of the Cloud-J
Feautrier radiative-transfer solver, run once over the package's fixed
climatological column. It replaces the empirical GEOS-Chem-fitted table of
EarthSciML/GasChem.jl#240/#243 - no GEOS-Chem output data is used, every
atmospheric assumption is explicit, and everything regenerates from this repo
with one script (rerun it after any change to `top_flux`, `σ_Raylay` or the
O2/O3 attenuation tables).

## Architecture

The diffuse (scattered) flux is a whole-column boundary-value problem that a
box model cannot compute locally - but for a fixed column it is a pure
function of (P, cosSZA), so it can be precomputed. Two data files, both
written by `precompute_table.jl`:

| File | Contents | Consumed by |
|---|---|---|
| `src/diffuse_flux_data.bson` | diffuse-only field `4*FJ(P, cosSZA)` per band | `FastJX` (default `fluxes = :scattering`): **computes** its spherical Beer-Lambert direct beam per evaluation (any pressure) and **adds** this field. `fluxes = :direct` drops it (clear-beam baseline). |
| `src/tropospheric_interpolation_data.bson` | total flux = `calc_direct_flux` + diffuse at each grid node | `FastJX_interpolation*`: compute no radiation, read the offline precompute of the online `FastJX`. |

`diffuse + direct = total` holds exactly at the nodes, so
`FastJX` and `FastJX_interpolation_troposphere` agree exactly at the nodes and
to within ~0.5% between them (e.g. 500 hPa) on all shared J's - the residual is
interpolation error of the direct part, the defining property of an offline
surrogate. Diffuse field valid
10-1000 hPa (held flat outside); axes 23 P-nodes x 61 cosSZA-nodes
(-0.2..1.0, twilight included). SOLFX is excluded from both files and
applied at runtime by `solar_flux_factor` (#224 convention).

## Regenerating

```
julia --project=. scripts/scattering/precompute_table.jl   # ~2 s
cp scripts/scattering/diffuse_flux_data.bson src/
cp scripts/scattering/tropospheric_interpolation_data_scattering.bson \
   src/tropospheric_interpolation_data.bson
```

`include_dependency` on both files makes the swap take effect on the next
precompile. Generated copies in this directory are gitignored.

## Fixed atmosphere (all explicit)

| Parameter | Value |
|---|---|
| Pressure grid | GEOS 72-layer hybrid at Psurf = 101325 Pa (`P_levels`) |
| Temperature | climatological `T_profile` (attenuation only; species sigma-phi use local T at runtime) |
| O3 profile | `mean_o3_profile`: **313.4 DU** total, peak 8.14 ppm @ ~6 hPa, surface 24 ppb (verified midlatitude-mean reasonable) |
| Rayleigh | `sigma_Raylay` per band (GEOS-Chem v2024-09 FJX_spec.dat grid, bins 17-18 at 380/574 nm); ssalb = 1; phase moments (1, 0, 0.5, 0, ...) |
| O2/O3 absorption | `sigma_O2_interp` / `sigma_O3_interp` at climatological T |
| Aerosol | **AOD = 0** (deliberate aerosol-free baseline) |
| Surface albedo | **0.10** Lambertian, spectrally flat |

Per-layer single-scattering albedo: `omega0 = OD_ray / (OD_ray + OD_abs)`;
`POMEGA = omega0 * [1, 0, 0.5, 0, ...]` (pure Rayleigh).

## Solver provenance

`cloudj_solver.jl` is a port of `OPMIE`/`MIESCT`/`LEGND0`/`GEN_ID`/`BLKSLV`
(Feautrier method, 4 Gauss points x 8-term Legendre) from
`cloud-j/src/Core/cldj_fjx_sub_mod.F90` at geoschem/cloud-j commit
`16d18b07` - the Cloud-J pinned by GEOS-Chem Classic 14.7.1 (science
codebase `b9f570e2`); reference clones under `~/geoschem-refs/`.
Clear-sky simplifications (documented in-file): JXTRA = 0 (cloud sub-layer
insertion unneeded), AMG = 1 (plane-parallel scattering); the solar beam
keeps full spherical geometry via the package's `sphere2J`, so the direct
component matches `calc_direct_flux` by construction.

## Validation

- `test_solver.jl` (78 tests): energy conservation - a conservative column
  over a mirror surface returns all incident flux out the top (rtol 2%);
  zero-scattering limit exactly reproduces the direct beam; optically-thin
  limit; twilight finite and non-negative; direct component matches
  `calc_direct_flux` to 1e-10.
- vs GEOS-Chem J (reconstructed from #243's fitted table with the
  cross-sections it was fitted with; #240 reports 0.95-0.99 held-out
  agreement), ratio ours/GEOS-Chem for the SuperFast channels:

  | | NO2 | H2O2 | H2CO b | CH3OOH | O3->O1D |
  |---|---|---|---|---|---|
  | surface, u0 = 1 | 1.07 | 1.06 | 1.09 | 1.06 | 0.89 |
  | surface, u0 = 0.5 | 1.07 | 1.14 | 1.13 | 1.12 | 1.17 |
  | surface, u0 = 0.2 | 1.14 | 1.23 | 1.23 | 1.19 | 1.34 |
  | 500 hPa, u0 = 1 | 1.04 | 0.97 | 1.02 | 0.98 | 0.77 |

  The low-sun excess is the expected AOD = 0 signature (the fit embeds one
  April day of CONUS aerosol).
- End-to-end (surface, u0 = 1, 298 K): J_NO2 **1.10e-2 s^-1**; enhancement
  over the beam-only baseline 2.1x (J_NO2), 2.6-3.0x (UV channels: H2O2,
  H2CO, CH3OOH, O3->O1D).
- Pinned 24-h box-model test solutions (compose_fastjx_superfast_test.jl,
  pollu_test.jl) are insensitive at their rtol = 1e-4 and pass unchanged.

Possible extensions: node-level check against Cloud-J standalone compiled
from the pinned clone; albedo (0.05/0.10/0.8-snow) and O3-column (+/-20%)
sensitivity tables.

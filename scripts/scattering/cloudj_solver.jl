# Julia port of the Cloud-J multiple-scattering solver (Feautrier method,
# 4 Gauss points x 8-term Legendre phase expansion) for the fixed clear-sky
# column of GasChem's Fast-JX port.
#
# Source: geoschem/cloud-j @ 16d18b07 (the Cloud-J pinned by GEOS-Chem Classic
# 14.7.1), src/Core/cldj_fjx_sub_mod.F90: OPMIE / MIESCT / LEGND0 / GEN_ID /
# BLKSLV. See README.md in this directory for the port plan and provenance.
#
# Simplifications relative to the Fortran (all valid for this use):
#   * JXTRA = 0 everywhere - the sub-layer insertion (EXTRAL1) exists for
#     optically thick cloud layers; our clear-sky Rayleigh column max layer
#     OD is << 1.
#   * AMG = 1 (plane-parallel scattering geometry, as Cloud-J uses when the
#     'geom' model is off). The *solar beam* attenuation FZ retains full
#     spherical geometry via GasChem's existing sphere2J machinery, so the
#     direct component reproduces GasChem.calc_direct_flux by construction.
#   * One wavelength band per solve (we loop bands in plain Julia).
#
# The solver works in Fast-JX normalized units: FZ = 1 at top of atmosphere;
# mean actinic flux at a point = FZ + 4*FJ (direct + diffuse). Multiply by
# GasChem.top_flux[band] for absolute flux.

using GasChem
using LinearAlgebra
using StaticArrays

# ---- Gauss quadrature constants (cldj_cmn_mod.F90) -------------------------
const M_ = 4
const M2_ = 8
const EMU = SVector{4, Float64}(
    0.06943184420297, 0.33000947820757, 0.66999052179243, 0.93056815579703
)
const WT = SVector{4, Float64}(
    0.17392742256873, 0.32607257743127, 0.32607257743127, 0.17392742256873
)

# ---- LEGND0: ordinary Legendre polynomials P_0..P_{n-1}(x) -----------------
function legnd0(x::Float64, n::Int)
    pl = zeros(n)
    pl[1] = 1.0
    pl[2] = x
    for i in 3:n
        den = i - 1
        pl[i] = pl[i - 1] * x * (2.0 - 1.0 / den) - pl[i - 2] * (1.0 - 1.0 / den)
    end
    return pl
end

# ---- GEN_ID: build the block tri-diagonal Feautrier system -----------------
#   A(L)*X(L-1) + B(L)*X(L) + C(L)*X(L+1) = H(L)
# pomega: (8, ND) phase fn x ss-albedo on the doubled grid
# fz:     (ND,) attenuated solar beam
# ztau:   (ND,) cumulative optical depth (0 at top, LZ=1)
# zflux:  direct solar flux onto surface (FSBOT)
# rfl:    (5,) Lambertian surface albedo for the 4 Gauss angles + direct beam
function gen_id(pomega, fz, ztau, zflux, rfl, PM, PM0, ND)
    B = zeros(M_, M_, ND); AA = zeros(M_, M_, ND); CC = zeros(M_, M_, ND)
    A = zeros(M_, ND); C = zeros(M_, ND); H = zeros(M_, ND)
    S = zeros(M_, M_); T = zeros(M_, M_); U = zeros(M_, M_)
    V = zeros(M_, M_); W = zeros(M_, M_)

    odd = (1, 3, 5, 7); even = (2, 4, 6, 8)
    sumP(l, terms, i) = sum(pomega[m, l] * PM[i, m] * PM0[m] for m in terms)
    sumPP(l, terms, i, j) = sum(pomega[m, l] * PM[i, m] * PM[j, m] for m in terms)

    # -- upper boundary (L1=1, L2=2), 2nd order ------------------------------
    for (L1, L2, lower) in ((1, 2, false), (ND, ND - 1, true))
        for i in 1:M_
            sum0 = sumP(L1, odd, i); sum2 = sumP(L2, odd, i)
            sum1 = sumP(L1, even, i); sum3 = sumP(L2, even, i)
            H[i, L1] = 0.5 * (sum0 * fz[L1] + sum2 * fz[L2])
            A[i, L1] = 0.5 * (sum1 * fz[L1] + sum3 * fz[L2])
        end
        for i in 1:M_, j in 1:i
            sum0 = sumPP(L1, odd, i, j); sum2 = sumPP(L2, odd, i, j)
            sum1 = sumPP(L1, even, i, j); sum3 = sumPP(L2, even, i, j)
            S[i, j] = -sum2 * WT[j]; S[j, i] = -sum2 * WT[i]
            T[i, j] = -sum1 * WT[j]; T[j, i] = -sum1 * WT[i]
            V[i, j] = -sum3 * WT[j]; V[j, i] = -sum3 * WT[i]
            B[i, j, L1] = -0.5 * (sum0 + sum2) * WT[j]
            B[j, i, L1] = -0.5 * (sum0 + sum2) * WT[i]
        end
        for i in 1:M_
            S[i, i] += 1.0; T[i, i] += 1.0; V[i, i] += 1.0; B[i, i, L1] += 1.0
            C[i, L1] = sum(S[i, k] * A[k, L1] / EMU[k] for k in 1:M_)
        end
        for i in 1:M_, j in 1:M_
            W[j, i] = sum(S[j, k] * T[k, i] / EMU[k] for k in 1:M_)
            U[j, i] = sum(S[j, k] * V[k, i] / EMU[k] for k in 1:M_)
        end
        if !lower
            deltau = ztau[L2] - ztau[L1]
            d2 = 0.25 * deltau
            for i in 1:M_, j in 1:M_
                B[i, j, L1] += d2 * W[i, j]
                CC[i, j, L1] = d2 * U[i, j]
            end
            for i in 1:M_
                H[i, L1] += 2.0 * d2 * C[i, L1]
                A[i, L1] = 0.0
            end
            for i in 1:M_
                d1 = EMU[i] / deltau
                B[i, i, L1] += d1
                CC[i, i, L1] -= d1
            end
        else
            # -- lower boundary with Lambertian albedo (v7.6 form) -----------
            deltau = ztau[L1] - ztau[L2]
            d2 = 0.25 * deltau
            sumrfl = sum(rfl[j] * EMU[j] * WT[j] for j in 1:M_)  # = avg(RFL)/2
            surfac = 4.0 / (1.0 + 2.0 * sumrfl)
            for i in 1:M_
                d1 = EMU[i] / deltau
                sum0 = d1 + d2 * (W[i, 1] + W[i, 2] + W[i, 3] + W[i, 4])
                for j in 1:M_
                    AA[i, j, L1] = -d2 * U[i, j]
                    B[i, j, L1] += d2 * W[i, j] - sum0 * surfac * rfl[j] * EMU[j] * WT[j]
                end
                H[i, L1] += -2.0 * d2 * C[i, L1] + sum0 * surfac * 0.25 * rfl[5] * zflux
            end
            for i in 1:M_
                d1 = EMU[i] / deltau
                AA[i, i, L1] += d1
                B[i, i, L1] += d1
                C[i, L1] = 0.0
            end
        end
    end

    # -- interior points: even LZ ('h', Legendre 2,4,6,8), odd LZ ('j', 1,3,5,7)
    for LL in 2:(ND - 1)
        terms = iseven(LL) ? even : odd
        deltau = ztau[LL + 1] - ztau[LL - 1]
        for i in 1:M_
            A[i, LL] = EMU[i] / deltau
            C[i, LL] = -A[i, LL]
            H[i, LL] = fz[LL] * sumP(LL, terms, i)
        end
        for i in 1:M_, j in 1:i
            sum0 = sumPP(LL, terms, i, j)
            B[i, j, LL] = -sum0 * WT[j]
            B[j, i, LL] = -sum0 * WT[i]
        end
        for i in 1:M_
            B[i, i, LL] += 1.0
        end
    end
    return B, CC, AA, A, H, C
end

# ---- BLKSLV: block tri-diagonal Thomas solve -------------------------------
# Returns FJ (mean diffuse intensity on the doubled grid), FJTOP, FJBOT, FIBOT.
function blkslv(pomega, fz, ztau, fsbot, rfl, PM, PM0, ND)
    B, CC, AA, A, H, C = gen_id(pomega, fz, ztau, fsbot, rfl, PM, PM0, ND)

    DD = zeros(M_, M_, ND)
    RR = zeros(M_, ND)

    # L = 1 (uses full CC block)
    E = inv(SMatrix{4, 4}(view(B, :, :, 1)))
    for j in 1:M_, i in 1:M_
        DD[i, j, 1] = -sum(E[i, k] * CC[k, j, 1] for k in 1:M_)
    end
    for j in 1:M_
        RR[j, 1] = sum(E[j, k] * H[k, 1] for k in 1:M_)
    end

    # L = 2 .. ND-1 (diagonal A and C)
    for L in 2:(ND - 1)
        for j in 1:M_
            for i in 1:M_
                B[i, j, L] += A[i, L] * DD[i, j, L - 1]
            end
            H[j, L] -= A[j, L] * RR[j, L - 1]
        end
        E = inv(SMatrix{4, 4}(view(B, :, :, L)))
        for j in 1:M_, i in 1:M_
            DD[i, j, L] = -E[i, j] * C[j, L]
        end
        for j in 1:M_
            RR[j, L] = sum(E[j, k] * H[k, L] for k in 1:M_)
        end
    end

    # L = ND (uses full AA block)
    L = ND
    for j in 1:M_
        for i in 1:M_
            B[i, j, L] += sum(AA[i, k, L] * DD[k, j, L - 1] for k in 1:M_)
        end
        H[j, L] -= sum(AA[j, k, L] * RR[k, L - 1] for k in 1:M_)
    end
    E = inv(SMatrix{4, 4}(view(B, :, :, L)))
    for j in 1:M_
        RR[j, L] = sum(E[j, k] * H[k, L] for k in 1:M_)
    end

    # back substitution
    for L in (ND - 1):-1:1
        for j in 1:M_
            RR[j, L] += sum(DD[j, k, L] * RR[k, L + 1] for k in 1:M_)
        end
    end

    # mean J (odd points) and H (even points)
    FJ = zeros(ND)
    for L in 1:2:ND
        FJ[L] = sum(RR[k, L] * WT[k] for k in 1:M_)
    end
    for L in 2:2:ND
        FJ[L] = sum(RR[k, L] * WT[k] * EMU[k] for k in 1:M_)
    end

    FJTOP = 4.0 * sum(RR[k, 1] * WT[k] * EMU[k] for k in 1:M_)
    sumb = sum(RR[k, ND] * WT[k] * EMU[k] for k in 1:M_)
    sumbr = sum(RR[k, ND] * WT[k] * EMU[k] * rfl[k] for k in 1:M_)
    sumrf = sum(WT[k] * EMU[k] * rfl[k] for k in 1:M_)
    sumbx = (4.0 * sumbr + fsbot * rfl[5]) / (1.0 + 2.0 * sumrf)
    FJBOT = 4.0 * sumb - sumbx
    FIBOT = zeros(5)
    FIBOT[5] = sumbx
    for j in 1:4
        FIBOT[j] = 2.0 * RR[j, ND] - sumbx
    end
    return FJ, FJTOP, FJBOT, FIBOT
end

# ---- MIESCT: Legendre setup + solve (single band) --------------------------
function miesct(pomega, fz, ztau, fsbot, rfl, u0, ND)
    PM = zeros(M_, M2_)
    for i in 1:M_
        PM[i, :] = legnd0(EMU[i], M2_)
    end
    PM0 = 0.25 .* legnd0(-u0, M2_)
    return blkslv(pomega, fz, ztau, fsbot, rfl, PM, PM0, ND)
end

# ---- Column assembly (OPMIE grid, JXTRA = 0) -------------------------------
# The doubled ("Feautrier") grid: ND = 2*L1 + 1 points, LZ = 1 at TAU = 0
# (top), LZ = ND at the surface. Odd LZ are layer edges ('j' points), even LZ
# mid-layer ('h') points. Mapping: edge L (1 = surface .. L1+1 = TOA edge)
# sits at LZ = ND + 1 - 2*(L - 1)... following the F90: LZ = ND + 2 - 2*L.
#
# Per-band inputs from GasChem's fixed column (all module constants):
#   layer OD:       GasChem.OD_ray_profile (18 x 73) + GasChem.OD_abs_profile
#   solar beam:     spherical-geometry attenuation via GasChem.sphere2J
const NLAY = 73                     # GEOS 72-layer + top pad, as in direct_flux.jl
const NEDGE = NLAY + 1
const ND_CLR = 2 * NLAY + 1          # 147

# Attenuated solar beam at each layer edge (normalized, spherical geometry).
# Edge 1 = surface .. edge 74 = TOA. Reuses the package's fine-grid air-mass
# machinery so the direct beam is identical to calc_direct_flux.
function solar_beam_edges(u0, band)
    return solar_beam_edges_allbands(u0)[:, band]
end

# All 18 bands at once - the air-mass factors depend only on (u0, edge), so
# compute each edge's AMF once and reuse across bands.
function solar_beam_edges_allbands(u0)
    ftau = zeros(NEDGE, 18)
    zhl = GasChem.z_profile
    ng = 2 * NLAY + 1
    for L in 1:NEDGE
        fine_index = 2L - 1
        amf = GasChem.sphere2J(u0, zhl, fine_index)
        amf[fine_index] > 0.0 || continue
        for band in 1:18
            dtau = view(GasChem.OD_total, :, band)
            tau = 0.5 * sum(dtau[div(i + 1, 2)] * amf[i] for i in 1:ng)
            ftau[L, band] = tau < 76.0 ? exp(-tau) : 0.0
        end
    end
    return ftau
end

# Rayleigh phase-function moments (P0 = 1, P2 = 1/2) times ss-albedo.
function pomega_layers(band)
    pomx = zeros(M2_, NLAY)
    for L in 1:NLAY
        od_ray = GasChem.OD_ray_profile[band, L]
        od_tot = GasChem.OD_total[L, band]
        w0 = od_tot > 0 ? od_ray / od_tot : 0.0
        pomx[1, L] = w0
        pomx[3, L] = 0.5 * w0
    end
    return pomx
end

"""
    solve_band(u0, band; albedo = 0.10)

Solve the clear-sky Rayleigh + O2/O3-absorption column for one wavelength
band at cosine-SZA `u0`. Returns `(J_edges, diag)` where `J_edges[L]` is the
normalized mean actinic flux (direct + diffuse; TOA incident = 1) at layer
edge L (1 = surface .. 74 = TOA), and `diag` holds (FJTOP, FJBOT, FSBOT).
Multiply by `GasChem.top_flux[band]` for absolute flux.
"""
function solve_band(u0, band; albedo = 0.1, ftau = nothing)
    ND = ND_CLR
    ftau = ftau === nothing ? solar_beam_edges(u0, band) : ftau
    pomx = pomega_layers(band)

    # cumulative column OD (TTAU): 0 at TOA edge, increasing downward
    ttau_edge = zeros(NEDGE)                       # index 1 = surface
    for L in NLAY:-1:1
        ttau_edge[L] = ttau_edge[L + 1] + GasChem.OD_total[L, band]
    end

    # phase fn interpolated to edges (F90 POMEGA1)
    pom_edge = zeros(M2_, NEDGE)
    pom_edge[:, 1] .= pomx[:, 1]
    pom_edge[:, NEDGE] .= pomx[:, NLAY]
    for L in 2:NLAY
        d0 = GasChem.OD_total[L, band]; d1 = GasChem.OD_total[L - 1, band]
        pom_edge[:, L] .= (pomx[:, L] .* d0 .+ pomx[:, L - 1] .* d1) ./ (d0 + d1)
    end

    # move onto the reversed doubled grid (odd points)
    ztau = zeros(ND); fz = zeros(ND); pom = zeros(M2_, ND)
    for L in 1:NEDGE
        LZ = ND + 2 - 2L
        ztau[LZ] = ttau_edge[L]
        fz[LZ] = ftau[L]
        pom[:, LZ] .= pom_edge[:, L]
    end
    # fill even 'h' points
    for LZ in 2:2:(ND - 1)
        ztau[LZ] = 0.5 * (ztau[LZ - 1] + ztau[LZ + 1])
        fz[LZ] = sqrt(fz[LZ - 1] * fz[LZ + 1])
        pom[:, LZ] .= 0.5 .* (pom[:, LZ - 1] .+ pom[:, LZ + 1])
    end

    # direct flux onto surface (F90: FTAU(1)/AMF(1,1); plane-parallel ~ u0)
    fsbot = u0 > 0.0 ? ftau[1] * u0 : 0.0
    rfl = fill(albedo, 5)

    FJ, FJTOP, FJBOT, FIBOT = miesct(pom, fz, ztau, fsbot, rfl, u0, ND)

    # mean actinic flux at layer edges: J = FZ + 4*FJ at the odd points
    J_edges = zeros(NEDGE)
    for L in 1:NEDGE
        LZ = ND + 2 - 2L
        J_edges[L] = fz[LZ] + 4.0 * FJ[LZ]
    end
    return J_edges, (FJTOP = FJTOP, FJBOT = FJBOT, FSBOT = fsbot)
end

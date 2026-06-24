export FastJX_interpolation_troposphere

BSON.@load joinpath(@__DIR__, "tropospheric_interpolation_data.bson") Z_all tropospheric_P cosSZA_vals
# Z_all is a vector of 18 matrices, each of which represents the actinic flux at different CSZA and Pressure.

interpolations_18_troposphere = []
for i in 1:18
    itp = interpolate(Z_all[i], BSpline(Linear()), OnGrid())
    f_in = Interpolations.scale(itp, tropospheric_P, cosSZA_vals)
    f_ext = extrapolate(f_in, Flat())
    push!(interpolations_18_troposphere, f_ext)
end

const interpolations_18_const = tuple(interpolations_18_troposphere...)

# Create symbolic wrapper functions for each interpolation
flux_interp_1(P, csa) = interpolations_18_const[1](ustrip(P), ustrip(csa))
flux_interp_2(P, csa) = interpolations_18_const[2](ustrip(P), ustrip(csa))
flux_interp_3(P, csa) = interpolations_18_const[3](ustrip(P), ustrip(csa))
flux_interp_4(P, csa) = interpolations_18_const[4](ustrip(P), ustrip(csa))
flux_interp_5(P, csa) = interpolations_18_const[5](ustrip(P), ustrip(csa))
flux_interp_6(P, csa) = interpolations_18_const[6](ustrip(P), ustrip(csa))
flux_interp_7(P, csa) = interpolations_18_const[7](ustrip(P), ustrip(csa))
flux_interp_8(P, csa) = interpolations_18_const[8](ustrip(P), ustrip(csa))
flux_interp_9(P, csa) = interpolations_18_const[9](ustrip(P), ustrip(csa))
flux_interp_10(P, csa) = interpolations_18_const[10](ustrip(P), ustrip(csa))
flux_interp_11(P, csa) = interpolations_18_const[11](ustrip(P), ustrip(csa))
flux_interp_12(P, csa) = interpolations_18_const[12](ustrip(P), ustrip(csa))
flux_interp_13(P, csa) = interpolations_18_const[13](ustrip(P), ustrip(csa))
flux_interp_14(P, csa) = interpolations_18_const[14](ustrip(P), ustrip(csa))
flux_interp_15(P, csa) = interpolations_18_const[15](ustrip(P), ustrip(csa))
flux_interp_16(P, csa) = interpolations_18_const[16](ustrip(P), ustrip(csa))
flux_interp_17(P, csa) = interpolations_18_const[17](ustrip(P), ustrip(csa))
flux_interp_18(P, csa) = interpolations_18_const[18](ustrip(P), ustrip(csa))

@register_symbolic flux_interp_1(P, csa)
@register_symbolic flux_interp_2(P, csa)
@register_symbolic flux_interp_3(P, csa)
@register_symbolic flux_interp_4(P, csa)
@register_symbolic flux_interp_5(P, csa)
@register_symbolic flux_interp_6(P, csa)
@register_symbolic flux_interp_7(P, csa)
@register_symbolic flux_interp_8(P, csa)
@register_symbolic flux_interp_9(P, csa)
@register_symbolic flux_interp_10(P, csa)
@register_symbolic flux_interp_11(P, csa)
@register_symbolic flux_interp_12(P, csa)
@register_symbolic flux_interp_13(P, csa)
@register_symbolic flux_interp_14(P, csa)
@register_symbolic flux_interp_15(P, csa)
@register_symbolic flux_interp_16(P, csa)
@register_symbolic flux_interp_17(P, csa)
@register_symbolic flux_interp_18(P, csa)

# Symbolic equations for actinic flux
function flux_eqs_interpolation(csa, P)
    flux_vals = []
    flux_vars = []
    @constants c_flux = 1.0 [
        unit = u"s^-1", description = "Constant actinic flux (for unit conversion)",
    ]

    interpolation_funcs = [
        flux_interp_1, flux_interp_2, flux_interp_3,
        flux_interp_4, flux_interp_5, flux_interp_6,
        flux_interp_7, flux_interp_8, flux_interp_9, flux_interp_10, flux_interp_11, flux_interp_12,
        flux_interp_13, flux_interp_14, flux_interp_15, flux_interp_16, flux_interp_17, flux_interp_18,
    ]

    for i in 1:18
        f = interpolation_funcs[i](P, csa)
        wl = WL[i]
        n1 = Symbol("F_", Int(round(wl)))
        v1 = @variables $n1(t) [unit = u"s^-1", description = "Actinic flux at $wl nm"]
        push!(flux_vars, only(v1))
        push!(flux_vals, f)
    end

    return flux_vars, (flux_vars .~ collect(flux_vals) .* c_flux), c_flux # TODO(CT): remove "collect" when https://github.com/SciML/ModelingToolkit.jl/issues/3888 is fixed.
end

# Photolysis species exposed by the interpolated Fast-JX, as (name, cross-section
# function) pairs. `:all` is the complete set (couples to the full
# `GEOSChemGasPhase` mechanism); `:superfast` is the reduced set used by
# `SuperFast` / `Pollu`. `j_o32OH` (= j_O31D adjusted for the O(1D) + H2O branch)
# is always added separately, so `:O31D` must be present in either list.
const _FJX_INTERP_J_FULL = [
    (:ActAld, j_mean_ActAld),
    (:PAN, j_mean_PAN),
    (:O3, j_mean_O3),
    (:NO3b, j_mean_NO3b),
    (:NO3a, j_mean_NO3a),
    (:N2O5, j_mean_N2O5),
    (:H2O2, j_mean_H2O2),
    (:H2COa, j_mean_H2COa),
    (:H2COb, j_mean_H2COb),
    (:O31D, j_mean_O31D),
    (:CH3OOH, j_mean_CH3OOH),
    (:NO2, j_mean_NO2),
    (:HOCl, j_mean_HOCl),
    (:MeAcr, j_mean_MeAcr),
    (:H1301, j_mean_H1301),
    (:CFCl3, j_mean_CFCl3),
    (:NO, j_mean_NO),
    (:Glyxlc, j_mean_Glyxlc),
    (:F114, j_mean_F114),
    (:CH3NO3, j_mean_CH3NO3),
    (:CHBr3, j_mean_CHBr3),
    (:F123, j_mean_F123),
    (:CHF2Cl, j_mean_CHF2Cl),
    (:OClO, j_mean_OClO),
    (:H1211, j_mean_H1211),
    (:BrO, j_mean_BrO),
    (:CH3Cl, j_mean_CH3Cl),
    (:MEKeto, j_mean_MEKeto),
    (:H2402, j_mean_H2402),
    (:PrAld, j_mean_PrAld),
    (:MeVKa, j_mean_MeVKa),
    (:MeVKb, j_mean_MeVKb),
    (:MeVKc, j_mean_MeVKc),
    (:ClNO3b, j_mean_ClNO3b),
    (:F113, j_mean_F113),
    (:HNO4, j_mean_HNO4),
    (:ClO, j_mean_ClO),
    (:CH2Br2, j_mean_CH2Br2),
    (:OCS, j_mean_OCS),
    (:F142b, j_mean_F142b),
    (:F115, j_mean_F115),
    (:CF3I, j_mean_CF3I),
    (:Glyxla, j_mean_Glyxla),
    (:CCl4, j_mean_CCl4),
    (:Cl2, j_mean_Cl2),
    (:CH3I, j_mean_CH3I),
    (:HNO2, j_mean_HNO2),
    (:Aceta, j_mean_Aceta),
    (:MeCCl3, j_mean_MeCCl3),
    (:Cl2O2, j_mean_Cl2O2),
    (:CH3Br, j_mean_CH3Br),
    (:HNO3, j_mean_HNO3),
    (:CF2Cl2, j_mean_CF2Cl2),
    (:Glyxlb, j_mean_Glyxlb),
    (:F141b, j_mean_F141b),
    (:ClNO3a, j_mean_ClNO3a),
    (:CH2Cl2, j_mean_CH2Cl2),
    (:O2, j_mean_O2),
    (:BrNO3, j_mean_BrNO3),
    (:GlyAld, j_mean_GlyAld),
    (:MGlyxl, j_mean_MGlyxl),
    (:HOBr, j_mean_HOBr),
    (:Acetb, j_mean_Acetb),
    (:BrCl, j_mean_BrCl),
]

const _FJX_INTERP_J_SUPERFAST = [
    (:H2O2, j_mean_H2O2),
    (:H2COa, j_mean_H2COa),
    (:H2COb, j_mean_H2COb),
    (:O31D, j_mean_O31D),
    (:CH3OOH, j_mean_CH3OOH),
    (:NO2, j_mean_NO2),
    (:ActAld, j_mean_ActAld),
    (:PAN, j_mean_PAN),
    (:NO3b, j_mean_NO3b),
    (:NO3a, j_mean_NO3a),
    (:N2O5, j_mean_N2O5),
    (:O3, j_mean_O3),
]

"""
    FastJX_interpolation_troposphere(t_ref; name=:FastJX, domaininfo=nothing, mech=:all)

Fast-JX photolysis using **interpolated** actinic fluxes from a precomputed
`(pressure, cosSZA)` table (`flux_eqs_interpolation`) rather than the online
direct-beam radiative-transfer integral of [`FastJX`](@ref). The
temperature-dependent cross sections and quantum yields are applied identically
to [`FastJX`](@ref), so the only approximation relative to the online scheme is
the flux interpolation. The shipped table spans ~10-1000 hPa (flux held flat
above the table top), so this targets tropospheric / lower-stratospheric columns.

`mech` selects which photolysis rates are exposed:

  - `:all` (default) -- the complete set, so this couples to the full
    `GEOSChemGasPhase` mechanism in addition to `SuperFast` / `Pollu`.
  - `:superfast` -- the reduced set needed by `SuperFast` / `Pollu`.

Both share the same flux table and cross sections, so any species present in both
produces identical j-values. Passing a `DomainInfo` attaches it as `SysDomainInfo`
metadata, mirroring [`FastJX`](@ref).

# Example

```julia
fj = FastJX_interpolation_troposphere(DateTime(2000, 1, 1))                    # full set
fj = FastJX_interpolation_troposphere(DateTime(2000, 1, 1); mech = :superfast)
```
"""
function FastJX_interpolation_troposphere(t_ref::AbstractFloat; name = :FastJX,
        domaininfo = nothing, mech::Symbol = :all)
    jlist = mech === :all ? _FJX_INTERP_J_FULL :
            mech === :superfast ? _FJX_INTERP_J_SUPERFAST :
            throw(ArgumentError("`mech` must be :all or :superfast, got :$(mech)"))

    consts = @constants begin
        T_unit = 1.0, [unit = u"K", description = "Unit temperature (for unit conversion)"]
        P_unit = 1.0, [unit = u"Pa", description = "Unit pressure"]
    end
    params = @parameters begin
        T = 298.0, [unit = u"K", description = "Temperature"]
        lat = 40.0, [description = "Latitude (Degrees)"]
        long = -97.0, [description = "Longitude (Degrees)"]
        P = 101325, [unit = u"Pa", description = "Pressure"]
        H2O = 450, [unit = u"ppb"]
        t_ref = t_ref, [unit = u"s", description = "Reference Unix time"]
    end
    @variables cosSZA(t) [description = "Cosine of the solar zenith angle"]
    @variables j_o32OH(t) [unit = u"s^-1"]

    flux_vars, fluxeqs, c_flux = flux_eqs_interpolation(cosSZA, P / P_unit)
    j_o31D_adj = adjust_j_o31D(ParentScope(T), ParentScope(P), ParentScope(H2O))

    # Build the requested j-variables and their cross-section equations.
    jvars = Num[]
    jeqs = Equation[]
    for (sp, meanf) in jlist
        nm = Symbol(:j_, sp)
        v = only(@variables $nm(t) [unit = u"s^-1"])
        push!(jvars, v)
        push!(jeqs, v ~ meanf(T / T_unit, flux_vars))
    end
    iO31D = findfirst(p -> p[1] === :O31D, jlist)
    @assert iO31D !== nothing "species list must include :O31D (needed for j_o32OH)"

    eqs = [
        cosSZA ~ cos_solar_zenith_angle(t + t_ref, lat, long);
        fluxeqs;
        j_o32OH ~ jvars[iO31D] * j_o31D_adj.j_O31D_adj;
        jeqs...
    ]

    fjx = System(
        eqs,
        t,
        [cosSZA; j_o32OH; jvars; flux_vars],
        [params; consts; c_flux];
        name = name,
        metadata = isnothing(domaininfo) ? Dict(CoupleType => FastJXCoupler) :
                   Dict(CoupleType => FastJXCoupler, SysDomainInfo => domaininfo),
        systems = [j_o31D_adj]
    )
    return flatten(fjx) # Need to do flatten because otherwise coupling doesn't work correctly
end
function FastJX_interpolation_troposphere(t_ref::DateTime; kwargs...)
    return FastJX_interpolation_troposphere(datetime2unix(t_ref); kwargs...)
end
FastJX_interpolation_troposphere(domain::DomainInfo; kwargs...) =
    FastJX_interpolation_troposphere(get_tref(domain); domaininfo = domain, kwargs...)

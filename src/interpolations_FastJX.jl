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
    @constants c_flux = 1.0 [unit = u"s^-1", description = "Constant actinic flux (for unit conversion)"]

    interpolation_funcs = [flux_interp_1, flux_interp_2, flux_interp_3, flux_interp_4, flux_interp_5, flux_interp_6,
                          flux_interp_7, flux_interp_8, flux_interp_9, flux_interp_10, flux_interp_11, flux_interp_12,
                          flux_interp_13, flux_interp_14, flux_interp_15, flux_interp_16, flux_interp_17, flux_interp_18]

    for i in 1:18
        f = interpolation_funcs[i](P, csa)
        wl = WL[i]
        n1 = Symbol("F_", Int(round(wl)))
        v1 = @variables $n1(t) [unit = u"s^-1", description = "Actinic flux at $wl nm"]
        push!(flux_vars, only(v1))
        push!(flux_vals, f)
    end

    flux_vars, (flux_vars .~ collect(flux_vals) .* c_flux), c_flux # TODO(CT): remove "collect" when https://github.com/SciML/ModelingToolkit.jl/issues/3888 is fixed.
end

"""
Description: This is a box model used to calculate the photolysis reaction rate constant using the Fast-JX scheme
(Neu, J. L., Prather, M. J., and Penner, J. E. (2007), Global atmospheric chemistry: Integrating over fractional cloud cover, J. Geophys. Res., 112, D11306, doi:10.1029/2006JD008007.)

Argument:

  - `t_ref`: Reference time for the model, can be a `DateTime` or a Unix timestamp (in seconds).

# Example

Build Fast-JX model:

```julia
fj = FastJX(DateTime(2000, 1, 1))
```
"""
function FastJX_interpolation_troposphere(t_ref::AbstractFloat; name = :FastJX)
    @constants T_unit = 1.0 [
        unit = u"K",
        description = "Unit temperature (for unit conversion)"
    ]
    @parameters T = 298.0 [unit = u"K", description = "Temperature"]
    @parameters lat = 40.0 [description = "Latitude (Degrees)"]
    @parameters long = -97.0 [description = "Longitude (Degrees)"]
    @parameters P = 101325 [unit = u"Pa", description = "Pressure"]
    @constants P_unit = 1.0 [unit = u"Pa", description = "Unit pressure"]
    @parameters H2O = 450 [unit = u"ppb"]
    @parameters t_ref = t_ref [unit = u"s", description = "Reference Unix time"]

    @variables j_H2O2(t) [unit = u"s^-1"]
    @variables j_H2COa(t) [unit = u"s^-1"]
    @variables j_H2COb(t) [unit = u"s^-1"]
    @variables j_O31D(t) [unit = u"s^-1"]
    @variables j_o32OH(t) [unit = u"s^-1"]
    @variables j_CH3OOH(t) [unit = u"s^-1"]
    @variables j_NO2(t) [unit = u"s^-1"]
    @variables cosSZA(t) [description = "Cosine of the solar zenith angle"]

    flux_vars, fluxeqs, c_flux = flux_eqs_interpolation(cosSZA, P/P_unit)
    j_o31D_adj = adjust_j_o31D(ParentScope(T), ParentScope(P), ParentScope(H2O))


    eqs = [cosSZA ~ cos_solar_zenith_angle(t + t_ref, lat, long);
           fluxeqs;
           j_H2O2 ~ j_mean_H2O2(T/T_unit, flux_vars);
           j_H2COa ~ j_mean_H2COa(T/T_unit, flux_vars);
           j_H2COb ~ j_mean_H2COb(T/T_unit, flux_vars);
           j_O31D ~ j_mean_O31D(T/T_unit, flux_vars);
           j_o32OH ~ j_O31D * j_o31D_adj.j_O31D_adj;
           j_CH3OOH ~ j_mean_CH3OOH(T/T_unit, flux_vars); 
           j_NO2 ~ j_mean_NO2(T/T_unit, flux_vars)]

    fjx = System(
        eqs,
        t,
        [j_H2O2, j_H2COa, j_H2COb, j_o32OH, j_O31D, j_CH3OOH, j_NO2, cosSZA, flux_vars...],
        [lat, long, T, P, H2O, t_ref, c_flux, T_unit, P_unit];
        name = name,
        metadata = Dict(CoupleType => FastJXCoupler),
        systems = [j_o31D_adj],
    )
    return flatten(fjx) # Need to do flatten because otherwise coupling doesn't work correctly
end
FastJX_interpolation_troposphere(t_ref::DateTime; kwargs...) = FastJX_interpolation_troposphere(datetime2unix(t_ref); kwargs...)

export FastJX_interpolation_troposphere

BSON.@load joinpath(@__DIR__, "interpolations_troposphere.bson") interpolations_18_troposphere
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

    flux_vars, (flux_vars .~ flux_vals .* c_flux)
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

    @variables j_h2o2(t) [unit = u"s^-1"]
    @variables j_CH2Oa(t) [unit = u"s^-1"]
    @variables j_CH2Ob(t) [unit = u"s^-1"]
    @variables j_o31D(t) [unit = u"s^-1"]
    @variables j_o32OH(t) [unit = u"s^-1"]
    @variables j_CH3OOH(t) [unit = u"s^-1"]
    @variables j_NO2(t) [unit = u"s^-1"]
    @variables cosSZA(t) [description = "Cosine of the solar zenith angle"]

    flux_vars, fluxeqs = flux_eqs_interpolation(cosSZA, P/P_unit)

    eqs = [cosSZA ~ cos_solar_zenith_angle(t + t_ref, lat, long);
           fluxeqs;
           j_h2o2 ~ j_mean_H2O2(T/T_unit, flux_vars)*0.0557; #0.0557 is a parameter to adjust the calculated H2O2 photolysis to appropriate magnitudes.
           j_CH2Oa ~ j_mean_CH2Oa(T/T_unit, flux_vars)*0.945; #0.945 is a parameter to adjust the calculated CH2Oa photolysis to appropriate magnitudes.
           j_CH2Ob ~ j_mean_CH2Ob(T/T_unit, flux_vars)*0.813; #0.813 is a parameter to adjust the calculated CH2Ob photolysis to appropriate magnitudes.
           j_o31D ~ j_mean_o31D(T/T_unit, flux_vars)*2.33e-21; #2.33e-21 is a parameter to adjust the calculated O(^3)1D photolysis to appropriate magnitudes.
           j_o32OH ~ j_o31D*adjust_j_o31D(T, P, H2O);
           j_CH3OOH ~ j_mean_CH3OOH(T/T_unit, flux_vars)*0.0931; #0.0931 is a parameter to adjust the calculated CH3OOH photolysis to appropriate magnitudes.
           j_NO2 ~ j_mean_NO2(T/T_unit, flux_vars)*0.444]

    ODESystem(
        eqs,
        t,
        [j_h2o2, j_CH2Oa, j_CH2Ob, j_o32OH, j_o31D, j_CH3OOH, j_NO2, cosSZA, flux_vars...],
        [lat, long, T, P, H2O, t_ref];
        name = name,
        metadata = Dict(:coupletype => FastJXCoupler)
    )
end
FastJX_interpolation_troposphere(t_ref::DateTime; kwargs...) = FastJX_interpolation_troposphere(datetime2unix(t_ref); kwargs...)

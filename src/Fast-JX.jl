export FastJX

# Effective wavelength in 18 bins covering 177–850 nm
const WL = SA_F32[187, 191, 193, 196, 202, 208, 211, 214, 261, 267, 277, 295, 303, 310, 316, 333, 380, 574]
# const WL_edges = SVector{19, Float32}([177.0; [(WL[i+1] + WL[i])/2 for i ∈ 1:length(WL)-1]; 850])
const ΔWL = SVector{17, Float32}([WL[i+1] - WL[i] for i ∈ 1:length(WL)-1])

"""
Create a vector of interpolators to interpolate the cross sections σ (TODO: What are the units?) for different wavelengths (in nm) and temperatures (in K).

We use use linear interpolation with flat extrapolation.
"""
function create_fjx_interp(temperatures::Vector{Float32}, cross_sections::Vector{SVector{18, Float32}})
    [linear_interpolation(temperatures, [x[i] for x ∈ cross_sections], extrapolation_bc=Flat()) for i ∈ 1:length(WL)]
end

#   Cross sections and quantum yield from GEOS-CHEM "FJX_spec.dat" for photo-reactive species mentioned in Superfast:

# Cross sections σ for O3(1D) for different wavelengths(18 bins), temperatures
# O3 -> O2 + O(1D) (superfast: O3 -> 2OH  including O(1D) + H2O -> 2OH)
const ϕ_o31D_jx = 1.0f0
σ_o31D_interp = create_fjx_interp([200.0f0, 260.0f0, 320.0f0], [
    SA_F32[4.842, 4.922, 5.071, 5.228, 6.040, 6.803, 7.190, 7.549, 9.000, 8.989, 8.929, 9.000, 8.901, 4.130, 8.985 * 0.1, 6.782 * 0.1, 0.000, 0.000] * 0.1f0,
    SA_F32[4.843, 4.922, 5.072, 5.229, 6.040, 6.802, 7.189, 7.549, 9.000, 8.989, 8.929, 9.000, 8.916, 4.656, 1.417, 6.995 * 0.1, 0.000, 0.000] * 0.1f0,
    SA_F32[4.843, 4.922, 5.072, 5.229, 6.040, 6.805, 7.189, 7.550, 9.000, 8.989, 8.929, 9.000, 8.967, 5.852, 2.919, 7.943, 0.000, 0.000] * 0.1f0,
])

# Cross sections σ for H2O2 for different wavelengths(18 bins), temperatures
# H2O2 -> OH + OH (same with superfast)
const ϕ_H2O2_jx = 1.0f0
σ_H2O2_interp = create_fjx_interp([200.0f0, 300.0f0], [
    SA_F32[2.325, 4.629, 5.394, 5.429, 4.447, 3.755, 3.457, 3.197, 5.346 * 0.1, 4.855 * 0.1, 3.423 * 0.1, 8.407 * 0.01, 5.029 * 0.01, 3.308 * 0.01, 2.221 * 0.01, 8.598 * 0.001, 1.807 * 0.0001, 0] * 10.0f0^-19.0f0,
    SA_F32[2.325, 4.629, 5.394, 5.429, 4.447, 3.755, 3.457, 3.197, 5.465 * 0.1, 4.966 * 0.1, 3.524 * 0.1, 9.354 * 0.01, 5.763 * 0.01, 3.911 * 0.01, 2.718 * 0.01, 1.138 * 0.01, 2.419 * 0.0001, 0] * 10.0f0^-19.0f0,
])

# Cross sections σ for CH2Oa for different wavelengths(18 bins), temperatures
# CH2O -> H + HO2 + CO (superfast: CH2O -> CO + 2HO2 including H + O2 -> HO2)
const ϕ_CH2Oa_jx = 1.0f0
σ_CH2Oa_interp = create_fjx_interp([223.0f0, 298.0f0], [
    SA_F32[0.0, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 3.143 * 0.1, 1.021, 1.269, 2.323, 2.498, 1.133, 2.183, 4.746 * 0.1, 0.000, 0.000] * 10.0f0^-20.0f0,
    SA_F32[0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 3.147 * 0.1, 1.018, 1.266, 2.315, 2.497, 1.131, 2.189, 4.751 * 0.1, 0.000, 0.000] * 10.0f0^-20.0f0,
])


# Cross sections σ for CH2Ob for different wavelengths(18 bins), temperatures
# CH2O -> CO + H2 （same with superfast）
const ϕ_CH2Ob_jx = 1.0f0
σ_CH2Ob_interp = create_fjx_interp([223.0f0, 298.0f0], [
    SA_F32[0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 3.642 * 10^-21, 5.787 * 10^-21, 5.316 * 10^-21, 8.181 * 10^-21, 7.917 * 10^-21, 4.011 * 10^-21, 1.081 * 10^-20, 1.082^-20, 2.088^-22, 0.000],
    SA_F32[0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 3.649 * 10^-21, 5.768 * 10^-21, 5.305 * 10^-21, 8.154 * 10^-21, 7.914 * 10^-21, 4.002 * 10^-21, 1.085 * 10^-20, 1.085 * 10^-20, 2.081^-22, 0.000],
])

# Cross sections σ for CH3OOH for different wavelengths(18 bins), temperatures
# CH3OOH -> OH + HO2 + CH2O (same with superfast)
const ϕ_CH3OOH_jx = 1.0f0
const σ_CH3OOH = SA_F32[0.000, 0.000, 0.000, 0.000, 0.000, 3.120, 2.882, 2.250, 2.716 * 0.1, 2.740 * 0.1, 2.143 * 0.1, 5.624 * 0.01, 3.520 * 0.01, 2.403 * 0.01, 1.697 * 0.01, 7.230 * 0.001, 6.973 * 0.0001, 0.000] * 10.f0^-19f0 # 298K
σ_CH3OOH_interp = [(T) -> σ_CH3OOH[i] for i in 1:18]

# Cross sections σ for NO2 for different wavelengths(18 bins), temperatures
# NO2 -> NO + O (superfast: NO2 -> NO + O3  including O + O2 -> O3)
const ϕ_NO2_jx = 1.0f0
σ_NO2_interp = create_fjx_interp([200.0f0, 294.0f0], [
    SA_F32[0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.835 * 0.1, 4.693 * 0.1, 7.705 * 0.1, 1.078, 1.470, 1.832, 2.181, 3.138, 4.321, 1.386 * 0.001] * 10.0f0^-19.0f0,
    SA_F32[0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.313 * 0.1, 4.694 * 0.1, 7.553 * 0.1, 1.063, 1.477, 1.869, 2.295, 3.448, 4.643, 4.345 * 0.001] * 10.0f0^-19.f0,
])


const actinic_flux = SA_F32[1.391E+12, 1.627E+12, 1.664E+12, 9.278E+11, 7.842E+12, 4.680E+12, 9.918E+12, 1.219E+13, 6.364E+14, 4.049E+14, 3.150E+14, 5.889E+14, 7.678E+14, 5.045E+14, 8.902E+14, 3.853E+15, 1.547E+16, 2.131E+17]

"""
    cos_solar_zenith_angle(lat, LST, DOY)
This function is to compute the cosine of the solar zenith angle, given the day of the year and local solar hour
The input variables: lat=latitude(°), LST=Local Solar Time(hour), DOY=Day of Year
    the cosine of the solar zenith angle (SZA) is given by:                                                                            .
           cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR)
                                                                            
           where LAT = the latitude angle,
                 DEC = the solar declination angle,
                 AHR = the hour angle, all in radians.  All in radians
"""
function cos_solar_zenith_angle(lat, time)
    dt = Dates.unix2datetime(time)
    LST = Dates.hour(dt) + Dates.minute(dt) / 60 + Dates.second(dt) / 3600.0 # Hour of day
    DOY = dayofyear(dt) # Day of year

    LAT = abs(lat * pi / 180.0) #lat>0, northern hemisphere; lat<0, southern hemisphere
    DEC = -23.45 * pi / 180.0 * cos(360.0 / 365.0 * (DOY + 10.0))
    AHR = 15.0 * pi / 180.0 * (LST - 12.0)
    CSZA = sin(LAT) * sin(DEC) + cos(LAT) * cos(DEC) * cos(AHR)
    return CSZA
end

"""
calculate actinic flux at the given cosine of the solar zenith angle `csa`
maximium actinic flux `max_actinic_flux`
"""
function calc_flux(csa, max_actinic_flux)
    max(zero(max_actinic_flux), max_actinic_flux * csa)
end

"""   
Get mean photolysis rates at different times
"""
function j_mean(σ_interp, ϕ, time, lat, Temperature)
    csa = cos_solar_zenith_angle(lat, time)
    j = zero(Temperature)
    for i in 1:17
        j1 = calc_flux(csa, actinic_flux[i]) * σ_interp[i](Temperature) * ϕ / WL[i]
        j2 = calc_flux(csa, actinic_flux[i+1]) * σ_interp[i+1](Temperature) * ϕ / WL[i+1]
        j += (j1 + j2)/2 * ΔWL[i]
    end
    j
end

j_mean_H2O2(t, lat, T) = j_mean(σ_H2O2_interp, ϕ_H2O2_jx, t, lat, T)
@register j_mean_H2O2(t, lat, T)
j_mean_CH2Oa(t, lat, T) = j_mean(σ_CH2Oa_interp, ϕ_CH2Oa_jx, t, lat, T)
@register j_mean_CH2Oa(t, lat, T)
j_mean_CH3OOH(t, lat, T) = j_mean(σ_CH3OOH_interp, ϕ_CH3OOH_jx, t, lat, T)
@register j_mean_CH3OOH(t, lat, T)
j_mean_NO2(t, lat, T) = j_mean(σ_NO2_interp, ϕ_NO2_jx,  t, lat, T)
@register j_mean_NO2(t, lat, T)
j_mean_o31D(t, lat, T) = j_mean(σ_o31D_interp, ϕ_o31D_jx,  t, lat, T)
@register j_mean_o31D(t, lat, T)

"""
Description: This is a box model used to calculate the photolysis reaction rate constant using the Fast-JX scheme 
(Neu, J. L., Prather, M. J., and Penner, J. E. (2007), Global atmospheric chemistry: Integrating over fractional cloud cover, J. Geophys. Res., 112, D11306, doi:10.1029/2006JD008007.)

Build Fast-JX model
# Example
``` julia
    @parameters t
    fj = FastJX(t)
```
"""
struct FastJX <: EarthSciMLODESystem
    sys::ODESystem
    function FastJX(t)
        @parameters T = 298
        @parameters lat = 30
        @parameters j_unit = 1 [unit = u"s^-1"]

        @variables j_h2o2(t) = 1.0097 * 10.0^-5 [unit = u"s^-1"]
        @variables j_CH2Oa(t) = 0.00014 [unit = u"s^-1"]
        @variables j_o31D(t) = 4.0 * 10.0^-3 [unit = u"s^-1"]
        @variables j_CH3OOH(t) = 8.9573 * 10.0^-6 [unit = u"s^-1"]
        @variables j_NO2(t) = 0.0149 [unit = u"s^-1"]
        # TODO(JL): What's difference between the two photolysis reactions of CH2O, do we really need both? 
        # (@variables j_CH2Ob(t) = 0.00014 [unit = u"s^-1"]) (j_CH2Ob ~ mean_J_CH2Ob(t,lat,T)*j_unit)

        eqs = [
            j_h2o2 ~  j_mean_H2O2(t, lat, T) * j_unit
            j_CH2Oa ~ j_mean_CH2Oa(t, lat, T) * j_unit
            j_o31D ~ j_mean_o31D(t, lat, T) * j_unit
            j_CH3OOH ~ j_mean_CH3OOH(t, lat, T) * j_unit
            j_NO2 ~ j_mean_NO2(t, lat, T) * j_unit
        ]

        new(ODESystem(eqs, t, [j_h2o2, j_CH2Oa, j_o31D, j_CH3OOH, j_NO2], [lat, j_unit, T]; name=:fastjx))
    end
end
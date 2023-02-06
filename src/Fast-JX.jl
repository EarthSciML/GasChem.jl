export FastJX

# Basic Input information: 

# Effective wavelength in 18 bins covering 177–850 nm
const WL = [187, 191, 193, 196, 202, 208, 211, 214, 261, 267, 277, 295, 303, 310, 316, 333, 380, 574]

#   Cross sections and quantum yield from GEOS-CHEM "FJX_spec.dat" for photo-reactive species mentioned in Superfast:

# Cross sections σ for O3(1D) for different wavelengths(18 bins), temperatures
# O3 -> O2 + O(1D) (superfast: O3 -> 2OH  including O(1D) + H2O -> 2OH)
const ϕ_o31D_jx = 1.0
const table_σ_o31D_jx = [
    # 200K
    [4.842, 4.922, 5.071, 5.228, 6.040, 6.803, 7.190, 7.549, 9.000, 8.989, 8.929, 9.000, 8.901, 4.130, 8.985 * 0.1, 6.782 * 0.1, 0.000, 0.000] * 0.1,
    # 260K
    [4.843, 4.922, 5.072, 5.229, 6.040, 6.802, 7.189, 7.549, 9.000, 8.989, 8.929, 9.000, 8.916, 4.656, 1.417, 6.995 * 0.1, 0.000, 0.000] * 0.1,
    # 320K
    [4.843, 4.922, 5.072, 5.229, 6.040, 6.805, 7.189, 7.550, 9.000, 8.989, 8.929, 9.000, 8.967, 5.852, 2.919, 7.943, 0.000, 0.000] * 0.1,
]
# Cross sections σ for H2O2 for different wavelengths(18 bins), temperatures
# H2O2 -> OH + OH (same with superfast)
const ϕ_H2O2_jx = 1.0
const table_σ_H2O2_jx = [
    # 200K
    [2.325, 4.629, 5.394, 5.429, 4.447, 3.755, 3.457, 3.197, 5.346 * 0.1, 4.855 * 0.1, 3.423 * 0.1, 8.407 * 0.01, 5.029 * 0.01, 3.308 * 0.01, 2.221 * 0.01, 8.598 * 0.001, 1.807 * 0.0001, 0] * 10^-19.0,
    # 300K
    [2.325, 4.629, 5.394, 5.429, 4.447, 3.755, 3.457, 3.197, 5.465 * 0.1, 4.966 * 0.1, 3.524 * 0.1, 9.354 * 0.01, 5.763 * 0.01, 3.911 * 0.01, 2.718 * 0.01, 1.138 * 0.01, 2.419 * 0.0001, 0] * 10^-19.0,
]

# Cross sections σ for CH2Oa for different wavelengths(18 bins), temperatures
# CH2O -> H + HO2 + CO (superfast: CH2O -> CO + 2HO2 including H + O2 -> HO2)
const ϕ_CH2Oa_jx = 1.0
const table_σ_CH2Oa_jx = [
    # 223K
    [0.0, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 3.143 * 0.1, 1.021, 1.269, 2.323, 2.498, 1.133, 2.183, 4.746 * 0.1, 0.000, 0.000] * 10^-20.0,
    # 298K
    [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 3.147 * 0.1, 1.018, 1.266, 2.315, 2.497, 1.131, 2.189, 4.751 * 0.1, 0.000, 0.000] * 10^-20.0,
]

# Cross sections σ for CH2Ob for different wavelengths(18 bins), temperatures
# CH2O -> CO + H2 （same with superfast）
const ϕ_CH2Ob_jx = 1.0
const table_σ_CH2Ob_jx = [
    # 223K
    [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 3.642 * 10^-21, 5.787 * 10^-21, 5.316 * 10^-21, 8.181 * 10^-21, 7.917 * 10^-21, 4.011 * 10^-21, 1.081 * 10^-20, 1.082^-20, 2.088^-22, 0.000],
    # 298K
    [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 3.649 * 10^-21, 5.768 * 10^-21, 5.305 * 10^-21, 8.154 * 10^-21, 7.914 * 10^-21, 4.002 * 10^-21, 1.085 * 10^-20, 1.085 * 10^-20, 2.081^-22, 0.000],
]

# Cross sections σ for CH3OOH for different wavelengths(18 bins), temperatures
# CH3OOH -> OH + HO2 + CH2O (same with superfast)
const ϕ_CH3OOH_jx = 1.0
const table_σ_CH3OOH_jx =
    [0.000, 0.000, 0.000, 0.000, 0.000, 3.120, 2.882, 2.250, 2.716 * 0.1, 2.740 * 0.1, 2.143 * 0.1, 5.624 * 0.01, 3.520 * 0.01, 2.403 * 0.01, 1.697 * 0.01, 7.230 * 0.001, 6.973 * 0.0001, 0.000] * 10^-19 # 298K

# Cross sections σ for NO2 for different wavelengths(18 bins), temperatures
# NO2 -> NO + O (superfast: NO2 -> NO + O3  including O + O2 -> O3)
const ϕ_NO2_jx = 1.0
const table_σ_NO2_jx = [
    # 200K
    [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.835 * 0.1, 4.693 * 0.1, 7.705 * 0.1, 1.078, 1.470, 1.832, 2.181, 3.138, 4.321, 1.386 * 0.001] * 10^-19.0,
    # 294K
    [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.313 * 0.1, 4.694 * 0.1, 7.553 * 0.1, 1.063, 1.477, 1.869, 2.295, 3.448, 4.643, 4.345 * 0.001] * 10^-19,
]

actinic_flux = [1.391E+12, 1.627E+12, 1.664E+12, 9.278E+11, 7.842E+12, 4.680E+12, 9.918E+12, 1.219E+13, 6.364E+14, 4.049E+14, 3.150E+14, 5.889E+14, 7.678E+14, 5.045E+14, 8.902E+14, 3.853E+15, 1.547E+16, 2.131E+17]

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
function cos_solar_zenith_angle(lat, LST, DOY)
    LAT = abs(lat * pi / 180) #lat>0, northern hemisphere; lat<0, southern hemisphere
    DEC = -23.45 * pi / 180 * cos(360 / 365 * (DOY + 10))
    AHR = 15 * pi / 180 * (LST - 12)
    CSZA = sin(LAT) * sin(DEC) + cos(LAT) * cos(DEC) * cos(AHR)
    return CSZA
end

"""
Functions to interpolate cross sections & quantum yields to calculate J values:
"""

function calc_J_H2O2(actinic_flux, T)
    if T <= 200
        σ = table_σ_H2O2_jx[1]
    elseif T >= 300
        σ = table_σ_H2O2_jx[2]
    elseif 200 < T < 300
        σ =
            (table_σ_H2O2_jx[2] .- table_σ_H2O2_jx[1]) * (T - 200) / (300 - 200) .+
            table_σ_H2O2_jx[1]
    end
    r = actinic_flux .* σ * ϕ_H2O2_jx
    return r
end

@register calc_J_H2O2(actinic_flux, T)

"""
    calc_J_o31D(actinic_flux,T)
This function is to calculate the J values of O3(1D) reaction with given actinic flux(unit: quanta*cm^-2*s^-1) and temperature(unit:K)
"""
function calc_J_o31D(actinic_flux, T)
    if T <= 200
        σ = table_σ_o31D_jx[1]
    elseif T > 320
        σ = table_σ_o31D_jx[3]
    elseif 200 < T <= 260
        σ =
            (table_σ_o31D_jx[2] .- table_σ_o31D_jx[1]) * (T - 200) / (260 - 200) .+
            table_σ_o31D_jx[1]
    elseif 260 < T <= 320
        σ =
            (table_σ_o31D_jx[3] .- table_σ_o31D_jx[2]) * (T - 260) / (320 - 260) .+
            table_σ_o31D_jx[2]
    end
    r = actinic_flux .* σ * ϕ_o31D_jx
    return r
end

@register calc_J_o31D(actinic_flux, T)

"""
    calc_J_CH2Oa(actinic_flux, T)
This function is to calculate the J values of CH2O(a) reaction with given actinic flux(unit: quanta*cm^-2*s^-1) and temperature(unit:K)
"""
function calc_J_CH2Oa(actinic_flux, T)
    if T <= 223
        r = actinic_flux .* table_σ_CH2Oa_jx[1] * ϕ_CH2Oa_jx
    elseif T >= 298
        r = actinic_flux .* table_σ_CH2Oa_jx[2] * ϕ_CH2Oa_jx
    elseif 223 < T < 298
        k = (T - 223) / (298 - 223)
        σ = (table_σ_CH2Oa_jx[2] - table_σ_CH2Oa_jx[1]) * k + table_σ_CH2Oa_jx[1]
        r = actinic_flux .* σ * ϕ_CH2Oa_jx
    end
    return r
end

@register calc_J_CH2Oa(actinic_flux, T)

"""
    calc_J_CH2Ob(actinic_flux, T)
This function is to calculate the J values of CH2O(a) reaction with given actinic flux(unit: quanta*cm^-2*s^-1) and temperature(unit:K)
"""
function calc_J_CH2Ob(actinic_flux, T)
    if T <= 223
        r = actinic_flux .* table_σ_CH2Ob_jx[1] * ϕ_CH2Ob_jx
    elseif T >= 298
        r = actinic_flux .* table_σ_CH2Ob_jx[2] * ϕ_CH2Ob_jx
    elseif 223 < T < 298
        k = (T - 223) / (298 - 223)
        σ = (table_σ_CH2Ob_jx[2] - table_σ_CH2Ob_jx[1]) * k + table_σ_CH2Ob_jx[1]
        r = actinic_flux .* σ * ϕ_CH2Ob_jx
    end
    return r
end

@register calc_J_CH2Ob(actinic_flux, T)

"""
    calc_J_CH3OOH(actinic_flux,T)
This function is to calculate the J values of CH3OOH reaction with given actinic flux(unit: quanta*cm^-2*s^-1) and temperature(unit:K)
"""
function calc_J_CH3OOH(actinic_flux, T)
    r = actinic_flux .* table_σ_CH3OOH_jx * ϕ_CH3OOH_jx
    return r
end

@register calc_J_CH3OOH(actinic_flux, T)

"""
    calc_J_NO2(actinic_flux, t)
This function is to calculate the J values of NO2 reaction with given actinic flux(unit: quanta*cm^-2*s^-1) and temperature(unit:K)
"""
function calc_J_NO2(actinic_flux, T)
    if T <= 200
        r = actinic_flux .* table_σ_NO2_jx[1] * ϕ_NO2_jx
    elseif T >= 294
        r = actinic_flux .* table_σ_NO2_jx[2] * ϕ_NO2_jx
    elseif 200 < T < 294
        k = (T - 200) / (294 - 200)
        σ = (table_σ_NO2_jx[2] - table_σ_NO2_jx[1]) * k + table_σ_NO2_jx[1]
        r = actinic_flux .* σ * ϕ_NO2_jx
    end
    return r
end

@register calc_J_NO2(actinic_flux, T)

"""
Get Doy of Year
"""
function Get_DOY(t)
    dayofyear(Dates.unix2datetime(t))
end

@register Get_DOY(t)

"""    
Get Hour of Day
"""
function Get_LST(t)
    LST =
        Dates.hour(Dates.unix2datetime(t)) +
        Dates.minute(Dates.unix2datetime(t)) / 60 +
        Dates.second(Dates.unix2datetime(t)) / 3600
end

@register Get_LST(t)

"""
calculate actinic flux at different time
"""
function calc_flux(t, lat)
    LST = Get_LST(t)
    DOY = Get_DOY(t)
    CSA = cos_solar_zenith_angle(lat, LST, DOY)
    r = zeros(18)
    r = max.(0, actinic_flux * CSA)
    return r
end

@register calc_flux(t, lat)

"""   
Get mean photolysis rates at different time
"""

function mean_J_H2O2(t, lat, T)
    flux = calc_flux(t, lat)
    j = calc_J_H2O2(flux, T) ./ WL
    r = 0
    for i = 1:17
        r += 1 / 2 * (j[i] + j[i+1]) * (WL[i+1] - WL[i])
    end
    return r
end

@register mean_J_H2O2(t, lat, T)

function mean_J_o31D(t, lat, T)
    flux = calc_flux(t, lat)
    j = calc_J_o31D(flux, T) ./ WL
    r = 0
    for i = 1:17
        r += 1 / 2 * (j[i] + j[i+1]) * (WL[i+1] - WL[i])
    end
    return r
end

@register mean_J_o31D(t, lat, T)

function mean_J_CH2Oa(t, lat, T)
    flux = calc_flux(t, lat)
    j = calc_J_CH2Oa(flux, T) ./ WL
    r = 0
    for i = 1:17
        r += 1 / 2 * (j[i] + j[i+1]) * (WL[i+1] - WL[i])
    end
    return r
end

@register mean_J_CH2Oa(t, lat, T)

function mean_J_CH2Ob(t, lat, T)
    flux = calc_flux(t, lat)
    j = calc_J_CH2Ob(flux, T) ./ WL
    r = 0
    for i = 1:17
        r += 1 / 2 * (j[i] + j[i+1]) * (WL[i+1] - WL[i])
    end
    return r
end

@register mean_J_CH2Ob(t, lat, T)

function mean_J_CH3OOH(t, lat, T)
    flux = calc_flux(t, lat)
    j = calc_J_CH3OOH(flux, T) ./ WL
    r = 0
    for i = 1:17
        r += 1 / 2 * (j[i] + j[i+1]) * (WL[i+1] - WL[i])
    end
    return r
end

@register mean_J_CH3OOH(t, lat, T)

function mean_J_NO2(t, lat, T)
    flux = calc_flux(t, lat)
    j = calc_J_NO2(flux, T) ./ WL
    r = 0
    for i = 1:17
        r += 1 / 2 * (j[i] + j[i+1]) * (WL[i+1] - WL[i])
    end
    return r
end

@register mean_J_NO2(t, lat, T)

"""
Description: This is a box model used to calculate the photolysis reaction rate constant using the Fast-JX scheme 
(Neu, J. L., Prather, M. J., and Penner, J. E. (2007), Global atmospheric chemistry: Integrating over fractional cloud cover, J. Geophys. Res., 112, D11306, doi:10.1029/2006JD008007.)

Build Fast-JX model
# Example
```
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
            j_h2o2 ~ mean_J_H2O2(t, lat, T) * j_unit
            j_CH2Oa ~ mean_J_CH2Oa(t, lat, T) * j_unit
            j_o31D ~ mean_J_o31D(t, lat, T) * j_unit
            j_CH3OOH ~ mean_J_CH3OOH(t, lat, T) * j_unit
            j_NO2 ~ mean_J_NO2(t, lat, T) * j_unit
        ]

        new(ODESystem(eqs, t, [j_h2o2, j_CH2Oa, j_o31D, j_CH3OOH, j_NO2], [lat, j_unit, T]; name=:fastjx))
    end
end
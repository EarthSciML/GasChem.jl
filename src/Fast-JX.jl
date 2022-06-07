export calc_J_o31D, calc_J_H2O2, calc_J_CH2Oa, calc_J_CH2Ob, calc_J_CH3OOH, calc_J_NO2
using StaticArrays
#   Description: This is a box model used to calculate the photolysis reaction rate constant using the Fast-JX scheme (Neu, J. L., Prather, M. J., and Penner, J. E. (2007), Global atmospheric chemistry: Integrating over fractional cloud cover, J. Geophys. Res., 112, D11306, doi:10.1029/2006JD008007.)

#   Basic Input infromation: 

    # Effective wavelength in 18 bins
    const WL = SA_F32[187,191,193,196,202,208,211,214,261,267,277,295,303,310,316,333,380,574] 

    #   Cross sections and quantum yield from GEOS-CHEM "FJX_spec.dat" for photo-reactive species mentioned in Superfast:

    # Cross sections σ for O3(1D) for different wavelengths(18 bins), temperatures
    # O3 -> O2 + O(1D) (superfast: O3 -> 2OH  including O(1D) + H2O -> 2OH)
    const ϕ_o31D_jx = 1.0
    const table_σ_o31D_jx = [[4.842, 4.922, 5.071, 5.228, 6.040, 6.803, 7.190, 7.549, 9.000, 8.989, 8.929, 9.000, 8.901, 4.130, 8.985*0.1, 6.782*0.1, 0.000, 0.000]*0.1 [4.843, 4.922, 5.072, 5.229, 6.040, 6.802, 7.189, 7.549, 9.000, 8.989, 8.929, 9.000, 8.916, 4.656, 1.417, 6.995*0.1, 0.000, 0.000]*0.1 [4.843, 4.922, 5.072, 5.229, 6.040, 6.805, 7.189, 7.550, 9.000, 8.989, 8.929, 9.000, 8.967, 5.852, 2.919, 7.943, 0.000, 0.000]*0.1]
    # first column: 200K
    # second column: 260K
    # third column: 320K

    # Cross sections σ for H2O2 for different wavelengths(18 bins), temperatures
    # H2O2 -> OH + OH (same with superfast)
    const ϕ_H2O2_jx = 1.0
    const table_σ_H2O2_jx = [[2.325, 4.629, 5.394, 5.429, 4.447, 3.755, 3.457, 3.197, 5.346*0.1, 4.855*0.1, 3.423*0.1, 8.407*0.01, 5.029*0.01, 3.308*0.01, 2.221*0.01, 8.598*0.001, 1.807*0.0001,0]*10^-19.0 [2.325, 4.629, 5.394, 5.429, 4.447, 3.755, 3.457, 3.197, 5.465*0.1, 4.966*0.1, 3.524*0.1, 9.354*0.01, 5.763*0.01, 3.911*0.01, 2.718*0.01, 1.138*0.01, 2.419*0.0001, 0]*10^-19.0]
    # first column : 200K
    # second column: 300K

    # Cross sections σ for CH2Oa for different wavelengths(18 bins), temperatures
    # CH2O -> H + HO2 + CO (superfast: CH2O -> CO + 2HO2)
    const ϕ_CH2Oa_jx = 1.0
    const table_σ_CH2Oa_jx = [[0.0, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 3.143*0.1, 1.021, 1.269, 2.323, 2.498, 1.133, 2.183, 4.746*0.1, 0.000, 0.000]*10^-20.0 [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 3.147*0.1, 1.018, 1.266, 2.315, 2.497, 1.131, 2.189, 4.751*0.1, 0.000, 0.000]*10^-20.0]
    # first column : 223K
    # second column: 298K

    # Cross sections σ for CH2Ob for different wavelengths(18 bins), temperatures
    # CH2O -> CO + H2 （same with superfast）
    const ϕ_CH2Ob_jx = 1.0
    const table_σ_CH2Ob_jx = [[0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 3.642*10^-21, 5.787*10^-21, 5.316*10^-21, 8.181*10^-21, 7.917*10^-21, 4.011*10^-21, 1.081*10^-20, 1.082^-20, 2.088^-22, 0.000] [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 3.649*10^-21, 5.768*10^-21, 5.305*10^-21, 8.154*10^-21, 7.914*10^-21, 4.002*10^-21, 1.085*10^-20, 1.085*10^-20, 2.081^-22, 0.000]]
    # first column: 223K
    # second column: 298K
    

    # Cross sections σ for CH3OOH for different wavelengths(18 bins), temperatures
    # CH3OOH -> OH + HO2 + CH2O (same with superfast)
    const ϕ_CH3OOH_jx = 1.0
    const table_σ_CH3OOH_jx = [0.000, 0.000, 0.000, 0.000, 0.000, 3.120, 2.882, 2.250, 2.716*0.1, 2.740*0.1, 2.143*0.1, 5.624*0.01, 3.520*0.01, 2.403*0.01, 1.697*0.01, 7.230*0.001, 6.973*0.0001, 0.000]*10^-19
    # 298K

    # Cross sections σ for NO2 for different wavelengths(18 bins), temperatures
    # NO2 -> NO + O (superfast: NO2 -> NO + O3  including O + O2 -> O3)
    const ϕ_NO2_jx = 1.0
    const table_σ_NO2_jx = [[0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.835*0.1, 4.693*0.1, 7.705*0.1, 1.078, 1.470, 1.832, 2.181, 3.138, 4.321, 1.386*0.001]*10^-19.0 [0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 2.313*0.1, 4.694*0.1, 7.553*0.1, 1.063, 1.477, 1.869, 2.295, 3.448, 4.643, 4.345*0.001]*10^-19] 
    # first column: 200K
    # second column:294K

#   Functions to calculate actinic flux
     
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
        LAT = abs(lat*pi/180) #lat>0, northern hemisphere; lat<0, southern hemisphere
        DEC = -23.45*pi/180*cos(360/365*(DOY+10))
        AHR = 15*pi/180*(LST-12)
        CSZA=sin(LAT)*sin(DEC)+cos(LAT)*cos(DEC)*cos(AHR)
        return CSZA
    end



#   Functions to interpolate cross sections & quantum yields to calculate J values:
     
"""
    calc_J_o31D(actinic_flux,T)
This function is to calculate the J values of O3(1D) reaction with given actinic flux(unit: quanta*cm^-2*s^-1) and temperature(unit:K)
"""
    function calc_J_o31D(actinic_flux,T)
        r = zeros(18)
        for i in 1:18
            if T <= 200
                r[i] = actinic_flux * table_σ_o31D_jx[i,1] * ϕ_o31D_jx
            elseif T > 320
                r[i] = actinic_flux * table_σ_o31D_jx[i,3] * ϕ_o31D_jx
            elseif 200 < T <= 260
                k = (T-200)/(260-200)
                σ = (table_σ_o31D_jx[i,2]-table_σ_o31D_jx[i,1])*k + table_σ_o31D_jx[i,1]
                r[i] = actinic_flux * σ * ϕ_o31D_jx
            elseif 260 < T <= 320
                k = (T-260)/(320-260)
                σ = (table_σ_o31D_jx[i,3]-table_σ_o31D_jx[i,2])*k + table_σ_o31D_jx[i,2]
                r[i] = actinic_flux * σ * ϕ_o31D_jx
            end
        end
        return r
    end

"""
    calc_J_H2O2(actinic_flux,T)
This function is to calculate the J values of H2O2 reaction with given actinic flux(unit: quanta*cm^-2*s^-1) and temperature(unit:K)
"""
    function calc_J_H2O2(actinic_flux,T)
        r = zeros(18)
        for i in 1:18
            if T <= 200
                r[i] = actinic_flux * table_σ_H2O2_jx[i,1] * ϕ_H2O2_jx
            elseif T >= 300
                r[i] = actinic_flux * table_σ_H2O2_jx[i,2] * ϕ_H2O2_jx
            elseif 200 < T < 300
                k = (T-200)/(300-200)
                σ = (table_σ_H2O2_jx[i,2]-table_σ_H2O2_jx[i,1])*k + table_σ_H2O2_jx[i,1]
                r[i] = actinic_flux * σ * ϕ_H2O2_jx
            end
        end
        return r
    end

"""
    calc_J_CH2Oa(actinic_flux, T)
This function is to calculate the J values of CH2O(a) reaction with given actinic flux(unit: quanta*cm^-2*s^-1) and temperature(unit:K)
"""
    function calc_J_CH2Oa(actinic_flux, T)
        r = zeros(18)
        for i in 1:18
            if T <= 223
                r[i] = actinic_flux * table_σ_CH2Oa_jx[i,1] * ϕ_CH2Oa_jx
            elseif T >= 298
                r[i] = actinic_flux * table_σ_CH2Oa_jx[i,2] * ϕ_CH2Oa_jx
            elseif 223 < T < 298
                k = (T-223)/(298-223)
                σ = (table_σ_CH2Oa_jx[i,2]-table_σ_CH2Oa_jx[i,1])*k + table_σ_CH2Oa_jx[i,1]
                r[i] = actinic_flux * σ * ϕ_CH2Oa_jx
            end
        end
        return r
    end

"""
    calc_J_CH2Ob(actinic_flux, T)
This function is to calculate the J values of CH2O(a) reaction with given actinic flux(unit: quanta*cm^-2*s^-1) and temperature(unit:K)
"""
    function calc_J_CH2Ob(actinic_flux, T)
        r = zeros(18)
        for i in 1:18
            if T <= 223
                r[i] = actinic_flux * table_σ_CH2Ob_jx[i,1] * ϕ_CH2Ob_jx
            elseif T >= 298
                r[i] = actinic_flux * table_σ_CH2Ob_jx[i,2] * ϕ_CH2Ob_jx
            elseif 223 < T < 298
                k = (T-223)/(298-223)
                σ = (table_σ_CH2Ob_jx[i,2]-table_σ_CH2Ob_jx[i,1])*k + table_σ_CH2Ob_jx[i,1]
                r[i] = actinic_flux * σ * ϕ_CH2Ob_jx
            end
        end
        return r
    end

"""
    calc_J_CH3OOH(actinic_flux,T)
This function is to calculate the J values of CH3OOH reaction with given actinic flux(unit: quanta*cm^-2*s^-1) and temperature(unit:K)
"""
    function calc_J_CH3OOH(actinic_flux,T)
        r = zeros(18)
        for i in 1:18
            r[i,1] = actinic_flux * table_σ_CH3OOH_jx[i] * ϕ_CH3OOH_jx
        end
        return r
    end

"""
    calc_J_NO2(actinic_flux, t)
This function is to calculate the J values of NO2 reaction with given actinic flux(unit: quanta*cm^-2*s^-1) and temperature(unit:K)
"""
    function calc_J_NO2(actinic_flux, T)
        r = zeros(18)
        for i in 1:18
            if T <= 200
                r[i] = actinic_flux * table_σ_NO2_jx[i,1] * ϕ_NO2_jx
            elseif T >= 294
                r[i] = actinic_flux * table_σ_NO2_jx[i,2] * ϕ_NO2_jx
            elseif 200 < T < 294
                k = (T-200)/(294-200)
                σ = (table_σ_NO2_jx[i,2]-table_σ_NO2_jx[i,1])*k + table_σ_NO2_jx[i,1]
                r[i] = actinic_flux * σ * ϕ_NO2_jx
            end
        end
        return r
    end





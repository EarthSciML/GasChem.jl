module EarthSciDataExt
using GasChem, EarthSciData, ModelingToolkit, EarthSciMLBase, DynamicQuantities
@register_unit molec 1
@register_unit mol_air 1u"mol"
@register_unit ppb 1u"mol/mol_air"

function EarthSciMLBase.couple2(
        c::GasChem.SuperFastCoupler,
        e::EarthSciData.NEI2016MonthlyEmisCoupler
)
    c, e = c.sys, e.sys

    @constants(MW_NO2=46.0055e-3,
        [unit=u"kg/mol", description="NO2 molar mass"],
        MW_NO=30.01e-3,
        [unit=u"kg/mol", description="NO molar mass"],
        MW_FORM=30.0260e-3,
        [unit=u"kg/mol", description="Formaldehyde molar mass"],
        #MW_CH4 = 16.0425e-3, [unit = u"kg/mol", description="Methane molar mass"], # CH4 is currently a constant in SuperFast.
        MW_CO=28.0101e-3,
        [unit=u"kg/mol", description="Carbon monoxide molar mass"],
        MW_SO2=64.0638e-3,
        [unit=u"kg/mol", description="Sulfur dioxide molar mass"],
        MW_ISOP=68.12e-3,
        [unit=u"kg/mol", description="Isoprene molar mass"],
        nmolpermol=1e9,
        [unit=u"ppb", description="nmol/mol, Conversion factor from mol to nmol"],
        R=8.31446261815324,
        [unit=u"m^3*Pa/mol/K", description="Ideal gas constant"],)

    # Emissions are in units of "kg/m3/s" and need to be converted to "ppb/s" or "nmol/mol/s".
    # To do this we need to convert kg of emissions to nmol of emissions,
    # and we need to convert m3 of air to mol of air.
    # nmol_emissions = kg_emissions * gperkg / MW_emission * nmolpermol = kg / kg/mol * nmol/mol = nmol
    # mol_air = m3_air / R / T * P = m3 / (m3*Pa/mol/K) / K * Pa = mol
    # So, the overall conversion is:
    # nmol_emissions / mol_air = (kg_emissions / MW_emission * nmolpermol) / (m3_air / R / T * P)
    uconv = nmolpermol * R * c.T / c.P # Conversion factor with MW factored out.
    operator_compose(
        convert(ODESystem, c),
        e,
        Dict(
            c.NO2 => e.NO2 => uconv / MW_NO2,
            c.NO => e.NO => uconv / MW_NO,
            c.CH2O => e.FORM => uconv / MW_FORM,
            #c.CH4 => e.CH4 => uconv / MW_CH4, # CH4 is currently a constant in SuperFast.
            c.CO => e.CO => uconv / MW_CO,
            c.ISOP => e.ISOP => uconv / MW_ISOP
        )
    )
end

function EarthSciMLBase.couple2(
        c::GasChem.GEOSChemGasPhaseCoupler,
        e::EarthSciData.NEI2016MonthlyEmisCoupler
)
    c, e = c.sys, e.sys

    @constants(MW_ACET=58.09*10^-3,
        [unit=u"kg/mol", description="Acetone molar mass"],
        MW_ALD2=44.06*10^-3,
        [unit=u"kg/mol", description="Acetaldehyde molar mass"],
        MW_BENZ=78.12*10^-3,
        [unit=u"kg/mol", description="Benzene molar mass"],
        MW_NO2=46.0055*10^-3,
        [unit=u"kg/mol", description="NO2 molar mass"],
        MW_NO=30.01*10^-3,
        [unit=u"kg/mol", description="NO molar mass"],
        MW_FORM=30.0260*10^-3,
        [unit=u"kg/mol", description="Formaldehyde molar mass"],
        MW_CH4=16.0425*10^-3,
        [unit=u"kg/mol", description="Methane molar mass"],
        MW_CO=28.0101*10^-3,
        [unit=u"kg/mol", description="Carbon monoxide molar mass"],
        MW_SO2=64.0638*10^-3,
        [unit=u"kg/mol", description="Sulfur dioxide molar mass"],
        MW_ISOP=68.12*10^-3,
        [unit=u"kg/mol", description="Isoprene molar mass"],
        nmolpermol=1e9,
        [unit=u"ppb", description="nmol/mol, Conversion factor from mol to nmol"],
        R=8.31446261815324,
        [unit=u"m^3*Pa/mol/K", description="Ideal gas constant"],)

    # Emissions are in units of "kg/m3/s" and need to be converted to "ppb/s" or "nmol/mol/s".
    # To do this we need to convert kg of emissions to nmol of emissions,
    # and we need to convert m3 of air to mol of air.
    # nmol_emissions = kg_emissions * gperkg / MW_emission * nmolpermol = kg / kg/mol * nmol/mol = nmol
    # mol_air = m3_air / R / T * P = m3 / (m3*Pa/mol/K) / K * Pa = mol
    # So, the overall conversion is:
    # nmol_emissions / mol_air = (kg_emissions / MW_emission * nmolpermol) / (m3_air / R / T * P)
    @constants(P_conv=1/6.02214076e+23,
        [
            unit=us"mol/cm^2",
            description="Number density conversion (handles lack of units in GEOS-Chem pressure)"
        ],
        mw_air=28.97e-3,
        [unit=u"kg/mol", description="Molar mass of air"],
        g=9.80665,
        [unit=u"m/s^2", description="Acceleration due to gravity"],)
    P = c.num_density * P_conv * mw_air * g
    uconv = nmolpermol * R * c.T / P # Conversion factor with MW factored out.
    #TODO(CT): Add missing couplings.
    operator_compose(
        convert(ODESystem, c),
        e,
        Dict(
            c.ACET => e.ACET => uconv / MW_ACET,
            c.ALD2 => e.ALD2 => uconv / MW_ALD2,
            c.BENZ => e.BENZ => uconv / MW_BENZ,
            c.NO2 => e.NO2 => uconv / MW_NO2,
            c.NO => e.NO => uconv / MW_NO,
            c.CH2O => e.FORM => uconv / MW_FORM,
            c.CH4 => e.CH4 => uconv / MW_CH4,
            c.CO => e.CO => uconv / MW_CO,
            c.SO2 => e.SO2 => uconv / MW_SO2,
            c.SO2 => e.SULF => uconv / MW_SO2,
            c.ISOP => e.ISOP => uconv / MW_ISOP
        )
    )
end

function EarthSciMLBase.couple2(c::GasChem.SuperFastCoupler, g::EarthSciData.GEOSFPCoupler)
    c, g = c.sys, g.sys
    @constants(T_inv=1,
        [unit=u"K^-1", description="Inverse of temperature"],
        P_inv=1,
        [unit=u"Pa^-1", description="Pressure"],
        ppb_unit=1,
        [unit=u"ppb"],)
    function water_concentration_ppb(RH, p, T)
        Tc = T - 273.15
        es_hPa = 6.112 * exp((17.62 * Tc) / (Tc + 243.12))
        es = es_hPa * 100.0          # Convert hPa to Pa
        e = RH * es                  # Actual water vapor partial pressure
        return (e / p) * 1e9       # ppb
    end

    c = param_to_var(c, :T, :P, :H2O)
    ConnectorSystem(
        [
            c.T ~ g.I3₊T,
            c.P ~ g.P,
            c.H2O ~ water_concentration_ppb(g.A3dyn₊RH, g.P*P_inv, g.I3₊T*T_inv)*ppb_unit
        ],
        c,
        g
    )
end

function EarthSciMLBase.couple2(
        c::GasChem.GEOSChemGasPhaseCoupler,
        gfp::EarthSciData.GEOSFPCoupler
)
    c, gfp = c.sys, gfp.sys

    #TODO(CT): Add missing couplings.
    c = param_to_var(c, :T, :num_density)
    @constants(P_conv=1/6.02214076e+23,
        [
            unit=us"mol/cm^2",
            description="Number density conversion (handles lack of units in GEOS-Chem pressure)"
        ],
        mw_air=28.97e-3,
        [unit=u"kg/mol", description="Molar mass of air"],
        g=9.80665,
        [unit=u"m/s^2", description="Acceleration due to gravity"],)
    ConnectorSystem([c.T ~ gfp.I3₊T, c.num_density * P_conv * mw_air * g ~ gfp.P], c, gfp)
end

function EarthSciMLBase.couple2(f::GasChem.FastJXCoupler, g::EarthSciData.GEOSFPCoupler)
    f, g = f.sys, g.sys

    @constants(T_inv=1,
        [unit=u"K^-1", description="Inverse of temperature"],
        P_inv=1,
        [unit=u"Pa^-1", description="Pressure"],
        ppb_unit=1,
        [unit=u"ppb"],)
    function water_concentration_ppb(RH, p, T)
        Tc = T - 273.15
        es_hPa = 6.112 * exp((17.62 * Tc) / (Tc + 243.12))
        es = es_hPa * 100.0          # Convert hPa to Pa
        e = RH * es                  # Actual water vapor partial pressure
        return (e / p) * 1e9       # ppb
    end

    f = param_to_var(f, :T, :lat, :long, :P, :H2O)
    ConnectorSystem(
        [
            f.T ~ g.I3₊T,
            f.P ~ g.P,
            f.lat ~ rad2deg(g.lat),
            f.long ~ rad2deg(g.lon),
            f.H2O ~ water_concentration_ppb(g.A3dyn₊RH, g.P*P_inv, g.I3₊T*T_inv)*ppb_unit
        ],
        f,
        g
    )
end

end

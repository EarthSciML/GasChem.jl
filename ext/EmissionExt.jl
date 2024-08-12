module EmissionExt
using GasChem, EarthSciData, ModelingToolkit, EarthSciMLBase, Unitful

function EarthSciMLBase.couple2(c::GasChem.SuperFastCoupler, e::EarthSciData.NEI2016MonthlyEmisCoupler)
    c, e = c.sys, e.sys

    @constants(
        MW_NO2 = 46.0055, [unit = u"g/mol", description="NO2 molar mass"],
        MW_NO = 30.01, [unit = u"g/mol", description="NO molar mass"],
        MW_FORM = 30.0260, [unit = u"g/mol", description="Formaldehyde molar mass"],
        MW_CH4 = 16.0425, [unit = u"g/mol", description="Methane molar mass"],
        MW_CO = 28.0101, [unit = u"g/mol", description="Carbon monoxide molar mass"],
        MW_SO2 = 64.0638, [unit = u"g/mol", description="Sulfur dioxide molar mass"],
        MW_ISOP = 68.12, [unit = u"g/mol", description="Isoprene molar mass"],

        gperkg = 1e3, [unit = u"g/kg", description="Conversion factor from kg to g"],
        nmolpermol = 1e9, [unit = u"nmol/mol", description="Conversion factor from mol to nmol"],
        R = 8.31446261815324, [unit = u"m^3*Pa/mol/K", description="Ideal gas constant"],
    )

    # Emissions are in units of "kg/m3/s" and need to be converted to "ppb/s" or "nmol/mol/s".
    # To do this we need to convert kg of emissions to nmol of emissions,
    # and we need to convert m3 of air to mol of air.
    # nmol_emissions = kg_emissions * gperkg / MW_emission * nmolpermol = kg * g/kg / g/mol * nmol/mol = nmol
    # mol_air = m3_air / R / T * P = m3 / (m3*Pa/mol/K) / K * Pa = mol
    # So, the overall conversion is:
    # nmol_emissions / mol_air = (kg_emissions * gperkg / MW_emission * nmolpermol) / (m3_air / R / T * P)
    uconv = gperkg * nmolpermol * R * c.T / c.P # Conversion factor with MW factored out.
    operator_compose(convert(ODESystem, c), e, Dict(
        c.NO2 => e.NO2 => uconv / MW_NO2,
        c.NO => e.NO => uconv / MW_NO,
        c.CH2O => e.FORM => uconv / MW_FORM,
        c.CH4 => e.CH4 => uconv / MW_CH4,
        c.CO => e.CO => uconv / MW_CO,
        c.SO2 => e.SO2 => uconv / MW_SO2,
        c.ISOP => e.ISOP => uconv / MW_ISOP,
    ))
end

end
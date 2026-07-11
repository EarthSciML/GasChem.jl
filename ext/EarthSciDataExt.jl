module EarthSciDataExt
using GasChem, EarthSciData, ModelingToolkit, EarthSciMLBase, DynamicQuantities

# Register custom units if not already registered by GasChem.
# This is needed because extension precompilation may occur before GasChem's
# unit registration takes effect.
if !(:molec in DynamicQuantities.Units.UNIT_SYMBOLS)
    @register_unit molec 1
end
if !(:mol_air in DynamicQuantities.Units.UNIT_SYMBOLS)
    @register_unit mol_air 1u"mol"
end
if !(:ppb in DynamicQuantities.Units.UNIT_SYMBOLS)
    @register_unit ppb 1u"mol/mol_air"
end

function EarthSciMLBase.couple2(
        c::GasChem.SuperFastCoupler,
        e::EarthSciData.NEI2016MonthlyEmisCoupler
    )
    c, e = c.sys, e.sys

    @constants(
        MW_NO2 = 46.0055e-3,
        [unit = u"kg/mol", description = "NO2 molar mass"],
        MW_NO = 30.01e-3,
        [unit = u"kg/mol", description = "NO molar mass"],
        MW_FORM = 30.026e-3,
        [unit = u"kg/mol", description = "Formaldehyde molar mass"],
        #MW_CH4 = 16.0425e-3, [unit = u"kg/mol", description="Methane molar mass"], # CH4 is currently a constant in SuperFast.
        MW_CO = 28.0101e-3,
        [unit = u"kg/mol", description = "Carbon monoxide molar mass"],
        MW_SO2 = 64.0638e-3,
        [unit = u"kg/mol", description = "Sulfur dioxide molar mass"],
        MW_SULF = 98.079e-3,
        [unit = u"kg/mol", description = "Sulfuric acid (NEI SULF sulfate emission) molar mass"],
        MW_NH3 = 17.031e-3,
        [unit = u"kg/mol", description = "Ammonia molar mass"],
        MW_ISOP = 68.12e-3,
        [unit = u"kg/mol", description = "Isoprene molar mass"],
        MW_Air = 28.97e-3,
        [unit = u"kg/mol", description = "Molar mass of air"],
        nmolpermol = 1.0e9,
        [unit = u"ppb", description = "nmol/mol, Conversion factor from mol to nmol"],
    )

    # Emissions are in units of "kg/kg air/s" and need to be converted to "ppb/s" or "nmol/mol/s".
    # NOTE: the translate table is a VECTOR of pairs, not a Dict — `c.NH3` appears twice
    # (NH3 + NH3_FERT sectors) and a Dict would silently keep only the last entry
    # (`operator_compose`'s `normalize_translate(::AbstractVector)` + `findall` handle
    # duplicate left-hand species natively).
    uconv = nmolpermol * MW_Air # Conversion factor with MW factored out.
    return operator_compose(
        convert(System, c),
        e,
        [
            c.NO2 => e.NO2 => uconv / MW_NO2,
            c.NO => e.NO => uconv / MW_NO,
            c.CH2O => e.FORM => uconv / MW_FORM,
            #c.CH4 => e.CH4 => uconv / MW_CH4, # CH4 is currently a constant in SuperFast.
            c.CO => e.CO => uconv / MW_CO,
            c.ISOP => e.ISOP => uconv / MW_ISOP,
            # SF-aerosol port: SO2 (gas) + SULF (direct sulfate) + NH3 (livestock/other +
            # fertilizer sectors) feed the het/cloud/ISORROPIA components.
            c.SO2 => e.SO2 => uconv / MW_SO2,
            c.SO4 => e.SULF => uconv / MW_SULF,
            c.NH3 => e.NH3 => uconv / MW_NH3,
            c.NH3 => e.NH3_FERT => uconv / MW_NH3
        ]
    )
end

function EarthSciMLBase.couple2(
        c::GasChem.GEOSChemGasPhaseCoupler,
        e::EarthSciData.NEI2016MonthlyEmisCoupler
    )
    c, e = c.sys, e.sys

    @constants(
        MW_ACET = 58.09e-3,
        [unit = u"kg/mol", description = "Acetone molar mass"],
        MW_ALD2 = 44.06e-3,
        [unit = u"kg/mol", description = "Acetaldehyde molar mass"],
        MW_BENZ = 78.12e-3,
        [unit = u"kg/mol", description = "Benzene molar mass"],
        MW_NO2 = 46.0055e-3,
        [unit = u"kg/mol", description = "NO2 molar mass"],
        MW_NO = 30.01e-3,
        [unit = u"kg/mol", description = "NO molar mass"],
        MW_FORM = 30.026e-3,
        [unit = u"kg/mol", description = "Formaldehyde molar mass"],
        MW_CH4 = 16.0425e-3,
        [unit = u"kg/mol", description = "Methane molar mass"],
        MW_CO = 28.0101e-3,
        [unit = u"kg/mol", description = "Carbon monoxide molar mass"],
        MW_SO2 = 64.0638e-3,
        [unit = u"kg/mol", description = "Sulfur dioxide molar mass"],
        MW_ISOP = 68.12e-3,
        [unit = u"kg/mol", description = "Isoprene molar mass"],
        MW_NH3 = 17.031e-3,
        [unit = u"kg/mol", description = "Ammonia molar mass"],
        MW_Air = 28.97e-3,
        [unit = u"kg/mol", description = "Molar mass of air"],
        nmolpermol = 1.0e9,
        [unit = u"ppb", description = "nmol/mol, Conversion factor from mol to nmol"],
    )

    # Emissions are in units of "kg/kg_air/s" and need to be converted to "ppb/s" or "nmol/mol/s".
    uconv = nmolpermol * MW_Air
    #TODO(CT): Add missing couplings.
    # NOTE: the translate table is a VECTOR of pairs, not a Dict. With a Dict, the duplicate
    # left-hand species below (c.SO2 ×2, c.NH3 ×2) silently collapsed last-wins — dropping
    # the e.SO2 gas emission (~98% of NEI sulfur) and the livestock e.NH3 sector entirely.
    # `operator_compose` handles duplicate left-hand entries natively for vectors
    # (`normalize_translate(::AbstractVector)` + `findall`).
    return operator_compose(
        convert(System, c),
        e,
        [
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
            c.ISOP => e.ISOP => uconv / MW_ISOP,
            # NEI ammonia: the merged grid keeps anthropogenic (livestock/other, `NH3`) and
            # fertilizer (`NH3_FERT`) sectors separate, so sum both into total NH3 — same
            # two-source pattern as SO2+SULF above. NH3 feeds the ISORROPIA aerosol partition
            # (IsorropiaOp). (If a NEI product reports NH3 inclusive of fertilizer, drop the
            # NH3_FERT term to avoid double-counting.)
            c.NH3 => e.NH3 => uconv / MW_NH3,
            c.NH3 => e.NH3_FERT => uconv / MW_NH3
        ]
    )
end

function EarthSciMLBase.couple2(c::GasChem.SuperFastCoupler, g::EarthSciData.GEOSFPCoupler)
    c, g = c.sys, g.sys
    @constants(
        T_inv = 1,
        [unit = u"K^-1", description = "Inverse of temperature"],
        P_inv = 1,
        [unit = u"Pa^-1", description = "Pressure"],
        ppb_unit = 1,
        [unit = u"ppb"],
    )
    function water_concentration_ppb(RH, p, T)
        Tc = T - 273.15
        es_hPa = 6.112 * exp((17.62 * Tc) / (Tc + 243.12))
        es = es_hPa * 100.0          # Convert hPa to Pa
        e = RH * es                  # Actual water vapor partial pressure
        return (e / p) * 1.0e9       # ppb
    end

    # H2O is met-coupled: an isconstantspecies parameter param_to_var'd here
    # and driven from GEOSFP RH. IsorropiaOp reads it per-cell through its
    # coordinate-observed function (AerosolExt.jl), so it does NOT need H2O in the state
    # vector — the observed `c.H2O ~ ...` equation is exactly what it consumes.
    c = param_to_var(c, :T, :P, :H2O)
    return ConnectorSystem(
        [
            c.T ~ g.I3₊T,
            c.P ~ g.P,
            c.H2O ~ water_concentration_ppb(g.A3dyn₊RH, g.P * P_inv, g.I3₊T * T_inv) * ppb_unit,
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
    @constants(
        R = 8.31446261815324,
        [unit = u"m^3*Pa/mol/K", description = "Ideal gas constant"],
    )
    return ConnectorSystem(
        [
            c.T ~ gfp.I3₊T,
            c.num_density ~ (gfp.P / R / gfp.I3₊T),
        ], c, gfp
    )
end

function EarthSciMLBase.couple2(f::GasChem.FastJXCoupler, g::EarthSciData.GEOSFPCoupler)
    f, g = f.sys, g.sys

    @constants(
        T_inv = 1,
        [unit = u"K^-1", description = "Inverse of temperature"],
        P_inv = 1,
        [unit = u"Pa^-1", description = "Pressure"],
        ppb_unit = 1,
        [unit = u"ppb"],
    )
    function water_concentration_ppb(RH, p, T)
        Tc = T - 273.15
        es_hPa = 6.112 * exp((17.62 * Tc) / (Tc + 243.12))
        es = es_hPa * 100.0          # Convert hPa to Pa
        e = RH * es                  # Actual water vapor partial pressure
        return (e / p) * 1.0e9       # ppb
    end

    f = param_to_var(f, :T, :lat, :long, :P, :H2O)
    return ConnectorSystem(
        [
            f.T ~ g.I3₊T,
            f.P ~ g.P,
            f.lat ~ rad2deg(g.lat),
            f.long ~ rad2deg(g.lon),
            f.H2O ~ water_concentration_ppb(g.A3dyn₊RH, g.P * P_inv, g.I3₊T * T_inv) * ppb_unit,
        ],
        f,
        g
    )
end

function EarthSciMLBase.couple2(
        c::GasChem.PolluCoupler,
        e::EarthSciData.NEI2016MonthlyEmisCoupler
    )
    c, e = c.sys, e.sys

    @constants(
        MW_NO2 = 46.0055e-3,
        [unit = u"kg/mol", description = "NO2 molar mass"],
        MW_NO = 30.01e-3,
        [unit = u"kg/mol", description = "NO molar mass"],
        MW_FORM = 30.026e-3,
        [unit = u"kg/mol", description = "Formaldehyde molar mass"],
        #MW_CH4 = 16.0425e-3, [unit = u"kg/mol", description="Methane molar mass"], # CH4 is currently a constant in SuperFast.
        MW_CO = 28.0101e-3,
        [unit = u"kg/mol", description = "Carbon monoxide molar mass"],
        # MW_SO2=64.0638e-3,
        # [unit=u"kg/mol", description="Sulfur dioxide molar mass"],
        MW_ALD2 = 44.052e-3,
        [unit = u"kg/mol", description = "Aldehyde molar mass"],
        # MW_ALDX=58.08e-3,
        # [unit=u"kg/mol", description="Aldehyde molar mass"],
        MW_Air = 28.97e-3,
        [unit = u"kg/mol", description = "Molar mass of air"],
        nmolpermol = 1.0e9,
        [unit = u"ppb", description = "nmol/mol, Conversion factor from mol to nmol"],
    )

    # Emissions are in units of "kg/kg air/s" and need to be converted to "ppb/s" or "nmol/mol/s".
    uconv = nmolpermol * MW_Air # Conversion factor with MW factored out.
    return operator_compose(
        convert(System, c),
        e,
        Dict(
            c.NO2 => e.NO2 => uconv / MW_NO2,
            c.NO => e.NO => uconv / MW_NO,
            c.CH2O => e.FORM => uconv / MW_FORM,
            c.CO => e.CO => uconv / MW_CO,
            c.ALD => e.ALD2 => uconv / MW_ALD2,
            # c.ALD => e.ALDX => uconv / MW_ALDX,
        )
    )
end

end

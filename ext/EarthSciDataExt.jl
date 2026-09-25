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
        MW_ISOP = 68.12e-3,
        [unit = u"kg/mol", description = "Isoprene molar mass"],
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
        MW_SULF = 98.079e-3, [unit = u"kg/mol", description = "H2SO4; Sulfuric acid"],
        MW_C2H4 = 28.054e-3, [unit = u"kg/mol", description = "C2H4; Ethylene"],
        MW_C2H6 = 30.070e-3, [unit = u"kg/mol", description = "C2H6; Ethane"],
        MW_C2H2 = 26.038e-3, [unit = u"kg/mol", description = "C2H2; Acetylene"],
        MW_C3H8 = 44.096e-3, [unit = u"kg/mol", description = "C3H8; Propane"],
        MW_TOLU = 92.141e-3, [unit = u"kg/mol", description = "C7H8; Toluene"],
        MW_MOH = 32.042e-3, [unit = u"kg/mol", description = "CH3OH; Methanol"],
        MW_EOH = 46.069e-3, [unit = u"kg/mol", description = "C2H5OH; Ethanol"],
        MW_NAP = 128.174e-3, [unit = u"kg/mol", description = "C10H8; Naphthalene"],
        MW_HNO2 = 47.014e-3, [unit = u"kg/mol", description = "HONO; Nitrous acid"],
        MW_HCl = 36.461e-3, [unit = u"kg/mol", description = "HCl; Hydrochloric acid"],
        MW_Cl2 = 70.906e-3, [unit = u"kg/mol", description = "Cl2; Molecular chlorine"],
        MW_CO2 = 44.010e-3, [unit = u"kg/mol", description = "CO2; Carbon dioxide"],
        MW_N2O = 44.013e-3, [unit = u"kg/mol", description = "N2O; Nitrous oxide"],
        MW_XYLE = 106.167e-3, [unit = u"kg/mol", description = "C8H10; Xylene (lumped)"],
        MW_MTPA = 136.234e-3, [unit = u"kg/mol", description = "C10H16; Monoterpenes"],
        MW_RCHO = 58.080e-3, [unit = u"kg/mol", description = "C2H5CHO; Propionaldehyde"],
        MW_PRPE = 42.081e-3, [unit = u"kg/mol", description = "C3H6; Propylene"],
        MW_ALK4 = 58.12e-3, [unit = u"kg/mol", description = "Lumped C4+C5 alkanes (GC species_database MW_g)"],
        MW_HNO4 = 79.01e-3, [unit = u"kg/mol", description = "HNO4; Peroxynitric acid"],
        MW_MEK = 72.107e-3, [unit = u"kg/mol", description = "C4H8O; Methyl ethyl ketone"],
        MW_Air = 28.97e-3,
        [unit = u"kg/mol", description = "Molar mass of air"],
        nmolpermol = 1.0e9,
        [unit = u"ppb", description = "nmol/mol, Conversion factor from mol to nmol"],
    )

    # Emissions are in units of "kg/kg_air/s" and need to be converted to "ppb/s" or "nmol/mol/s".
    uconv = nmolpermol * MW_Air
    #TODO(CT): Add missing couplings.
    # NOTE: the translate table is a VECTOR of pairs, not a Dict. Upstream's Dict literal
    # keyed c.SO2 twice (=> e.SULF and => e.SO2), silently collapsing last-wins — dropping
    # the e.SO2 gas emission (~98% of NEI sulfur). The Vector form preserves every entry,
    # and `operator_compose` handles duplicate left-hand entries natively for vectors
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
            c.SO4 => e.SULF => uconv / MW_SULF,   # NEI SULF is direct sulfate emission -> SO4
            c.ISOP => e.ISOP => uconv / MW_ISOP,
            # --- the 18 remaining NEI species not yet in the base coupling ---
            c.C2H4 => e.ETH => uconv / MW_C2H4,
            c.C2H6 => e.ETHA => uconv / MW_C2H6,
            c.C2H2 => e.ETHY => uconv / MW_C2H2,
            c.C3H8 => e.PRPA => uconv / MW_C3H8,
            c.TOLU => e.TOL => uconv / MW_TOLU,
            c.MOH => e.MEOH => uconv / MW_MOH,
            c.EOH => e.ETOH => uconv / MW_EOH,
            c.NAP => e.NAPH => uconv / MW_NAP,
            c.HNO2 => e.HONO => uconv / MW_HNO2,
            c.HCl => e.HCL => uconv / MW_HCl,
            c.Cl2 => e.CL2 => uconv / MW_Cl2,
            c.CO2 => e.CO2_INV => uconv / MW_CO2,
            c.N2O => e.N2O_INV => uconv / MW_N2O,
            c.XYLE => e.XYLMN => uconv / MW_XYLE,
            c.MTPA => e.TERP => uconv / MW_MTPA,
            c.RCHO => e.ALDX => uconv / MW_RCHO,
            c.PRPE => e.OLE => uconv / MW_PRPE,
            # NEI variables GEOS-Chem maps but this table did not. HEMCO applies only
            # temporal/spatial scale factors to each (26 = time-of-day, 212/213 = NEI99
            # day-of-week, 254 = VOC year scale, 1007 = CONUS mask) — no mass or carbon
            # conversion — so these are plain mass -> mole conversions like the rest.
            c.ALK4 => e.PAR => uconv / MW_ALK4,    # HEMCO EPA16_ALK4__*PAR; ALK4 had ZERO anthropogenic emission
            c.PRPE => e.IOLE => uconv / MW_PRPE,   # HEMCO maps both OLE and IOLE to PRPE
            c.HNO4 => e.PNA => uconv / MW_HNO4,    # HEMCO EPA16_HNO4__*PNA
            c.MEK => e.KET => uconv / MW_MEK,
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

    # Note on P vs SuperFast:
    # P is defined as @parameter but unused in any reaction rate (num_density
    # is used instead), so the compiler strips it. Use GEOSFP₊P directly.
    #
    # H2O is a constant species (@parameter with isconstantspecies=true, the
    # GEOS-Chem met-driven fixed-species treatment) coupled from GEOSFP meteorology via param_to_var,
    # same as SuperFast. This drives the explicit O1D + H2O --> 2OH OH source with
    # real per-cell humidity instead of a frozen boundary-layer constant.
    @constants(
        R = 8.31446261815324,
        [unit = u"m^3*Pa/mol/K", description = "Ideal gas constant"],
        T_inv = 1,
        [unit = u"K^-1", description = "Inverse of temperature"],
        P_inv = 1,
        [unit = u"Pa^-1", description = "Inverse of pressure"],
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

    c = param_to_var(c, :T, :num_density, :H2O)
    return ConnectorSystem(
        [
            c.T ~ gfp.I3₊T,
            c.num_density ~ (gfp.P / R / gfp.I3₊T),
            c.H2O ~ water_concentration_ppb(gfp.A3dyn₊RH, gfp.P * P_inv, gfp.I3₊T * T_inv) * ppb_unit,
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

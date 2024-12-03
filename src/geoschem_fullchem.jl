export GEOSChemGasPhase

struct GEOSChemGasPhaseCoupler
    sys
end

"""
GEOS-Chem full-chem mechanism adapted from GEOS-Chem version 14.1.1
* Adapted from file https://github.com/geoschem/geos-chem/blob/4722f288e90291ba904222f4bbe4fc216d17c34a/KPP/fullchem/fullchem.eqn
* The GEOS-Chem license applies: https://github.com/geoschem/geos-chem/blob/main/LICENSE.txt

If the keyword argument `rxn_sys` is set to `true`, the function will return a reaction system instead of an ODE system.

===============================================================================
REFERENCES (alphabetical order)
===============================================================================
* Atkinson2003: Atkinson and Arey, Chem. Rev., doi:10.1021/cr0206420, 2003.
* Bates2014:    Bates et al., J. Phys. Chem A, 118, doi:10.1021/jp4107958, 2014.
* Bates2019:    Bates and Jacob, Atmos. Chem. Phys., doi:10.5194/acp-19-9613-2019, 2019.
* Bates2021a:   Bates et al, JGR, https://doi.org/10.1029/2020JD033439, 2021.
* Bates2021b:   Bates et al, ACP, https://doi.org/10.5194/acp-2021-605, 2021.
* Browne2011:   Browne et al., Atmos. Chem. Phys., doi:10.5194/acp-11-4209-2011, 2011.
* Browne2014:   Browne et al., Atmos. Chem. Phys., doi:10.5194/acp-14-1225-2014, 2014.
* Chen2017:     Chen et al., Geophys. Res. Lett., doi:10.1002/2017GL073812, 2017.
* Crounse2012:  Crounse et al., J. Phys. Chem. A, doi:10.1021/jp211560u, 2012.
* Eastham2014:  Eastham et al., Atmos. Env., doi:10.1016/j.atmosenv.2014.02.001, 2014.
* Fischer2014:  Fischer et al., Atmos. Chem. Phys., doi:10.5194/acp-14-2679-2014, 2014.
* Fisher2016:   Fisher et al., Atmos. Chem. Phys., doi:10.5194/acp-16-5969-2016, 2016.
* Fisher2018:   Fisher et al., J. Geophys. Res., doi:10.1029/2018JD029046, 2018.
* Fry2014:      Fry et al. Environ. Sci. Technol., doi:10.1021/es502204x, 2014.
* Gill2002:     Gill and Hites, J. Phys. Chem. A, doi:10.1021/jp013532, 2002.
* Goliff2013:   Goliff et al., Atmos. Environ., doi:10.1016/j.atmosenv.2012.11.038, 2013.
* Jacobs2014:   Jacobs et al., Atmos. Chem. Phys., doi:10.5194/acp-14-8933-2014, 2014.
* Jenkin2015:   Jenkin et al., Atmos. Chem. Phys., doi:10.5194/acp-15-11433-2015, 2015.
* Kasibhatla2018: Kasibhatla et al., Atmos. Chem. Phys., doi:10.5194/acp-18-11185-2018, 2018
* IUPAC ROO_19: https://iupac-aeris.ipsl.fr/htdocs/datasheets/pdf/ROO_19_CH3O2_NO3.pdf
* JPL 10-6:     JPL Publication 10-6, https://jpldataeval.jpl.nasa.gov/previous_evaluations.html, 2011.
* JPL 15-10:    JPL Publication 15-10, https://jpldataeval.jpl.nasa.gov, 2015.
* Kwon2020:     Kwon et al, Elementa, https://doi.org/10.1525/elementa.2021.00109, 2020.
* Lee2014:      Lee et al., J. Phys. Chem. A, doi:10.1021/jp4107603, 2014.
* Marais2016:   Marais et al., Atmos. Chem. Phys, doi:10.5194/acp-16-1603-2016, 2016.
* Miller2017:   Miller et al., Atmos. Chem. Phys. Discuss., doi:10.5194/acp-2016-1042, 2017.
* Millet2015:   Millet et al., Atmos. Chem. Phys., doi:10.5194/acp-15-6283-2015, 2015.
Moch et al, JGR, https, * Moch2020 # //doi.org/10.1029/2020JD032706, 2020.
* Muller2014:   Muller et al., Atmos. Chem. Phys., doi:10.5194/acp-14-2497-2014, 2014.
* Parrella2012: Parrella et al. Atmos. Chem. Phys, doi:10.5194/acp-12-6723-2012, 2012.
* Paulot2009:   Paulot et al., Atmos. Chem. Phys., doi:10.5194/acp-9-1479-2009, 2009a and
                Paulot et al., Science, doi:10.1126/science.1172910, 2009b.
* Peeters2010:  Peeters and Muller, Phys. Chem. Chem. Phys., doi:10.1039/C0CP00811G, 2010.
* Peeters2014:  Peeters et al., J. Phys. Chem. A, doi:10.1021/jp5033146, 2014.
* Pye2010:      Pye et al., Atmos. Chem. Phys., doi:10.5194/acp-10-11261-2010, 2010.
* Roberts1992:  Roberts and Bertman, Int. J. Chem. Kinet., doi:10.1002/kin.550240307, 1992.
* Sherwen2016b: Sherwen et al., Atmos. Chem. Phys., doi:10.5194/acp-16-12239-2016, 2016b.
* Sherwen2017:  Sherwen et al., Faraday Discuss., doi:10.1039/C7FD00026J, 2017.
* StClair2016:  St. Clair et al., J. Phys. Chem. A, doi:10.1021/acs.jpca.5b065322016, 2016.
* Travis2016:   Travis et al., Atmos. Chem. Phys., doi:10.5194/acp-16-13561-2016, 2016.
* Wolfe2012:    Wolfe et al., Phys. Chem. Chem. Phys., doi: 10.1039/C2CP40388A, 2012.
* Xie2013:      Xie et al., Atmos. Chem. Phys., doi:10.5194/acp-13-8439-2013, 2013.
"""
function GEOSChemGasPhase(; name=:GEOSChemGasPhase, rxn_sys=false)

    # Create reaction rate constant system constructors
    rate_systems = []
    i = 1
    function arr(t, T, a0, b0, c0)
        sys = arrhenius(t, T, a0, b0, c0, name=Symbol(:arrhenius_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function arrnodim(t, T, a0, b0, c0)
        sys = arrhenius(t, T, a0, b0, c0, unit=u"ppb*s^-1" / u"ppb", name=Symbol(:arrhenius_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function arrsq(t, T, a0, b0, c0)
        sys = arrhenius(t, T, a0, b0, c0, unit=u"ppb^-2*s^-1", name=Symbol(:arrhenius_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function c(c)
        sys = constant_k(t, c, name=Symbol(:constant_k_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function regressT(t, T, a_0, b_0, T_0)
        sys = regress_T(t, T, a_0, b_0, T_0; name=Symbol(:acet_oh_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function arr3(t, T, num_density, a1, b1, c1, a2, b2, c2, fv)
        sys = arr_3rdbody(t, T, num_density, a1, b1, c1, a2, b2, c2, fv;
            name=Symbol(:arr_3rdbody_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function arr3nodim(t, T, num_density, a1, b1, c1, a2, b2, c2, fv)
        sys = arr_3rdbody(t, T, num_density, a1, b1, c1, a2, b2, c2, fv;
            unit=u"ppb*s^-1" / u"ppb", name=Symbol(:arr_3rdbody_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rHO2HO2(t, T, num_density, H2O, a0, c0, a1, c1)
        sys = rate_HO2HO2(t, T, num_density, H2O, a0, c0, a1, c1,
            name=Symbol(:rate_HO2HO2_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rOHCO(t, T, num_density)
        sys = rate_OHCO(t, T, num_density;
            name=Symbol(:rate_OHCO_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rRO2NOa1(t, T, a0, c0)
        sys = rate_RO2NO_a1(t, T, a0, c0; name=Symbol(:rate_RO2NO_a1_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rRO2NOb1(t, T, a0, c0)
        sys = rate_RO2NO_b1(t, T, a0, c0; name=Symbol(:rate_RO2NO_b1_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rRO2NOa2(t, T, num_density, a0, c0, a1)
        sys = rate_RO2NO_a2(t, T, num_density, a0, c0, a1; name=Symbol(:rate_RO2NO_a2_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rRO2NOb2(t, T, num_density, a0, c0, a1)
        sys = rate_RO2NO_b2(t, T, num_density, a0, c0, a1; name=Symbol(:rate_RO2NO_b2_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function tbrnch(t, T, a0, b0, c0, a1, b1, c1)
        sys = tbranch(t, T, a0, b0, c0, a1, b1, c1; name=Symbol(:tbranch_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rOHHNO3(t, T, num_density, a0, c0, a1, c1, a2, c2)
        sys = rate_OHHNO3(t, T, num_density, a0, c0, a1, c1, a2, c2;
            name=Symbol(:rate_OHHNO3_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function eqconst(t, T, num_density, a0, c0, a1, b1, a2, b2, fv)
        sys = eq_const(t, T, num_density, a0, c0, a1, b1, a2, b2, fv,
            name=Symbol(:eq_const_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function eqconstnodim(t, T, num_density, a0, c0, a1, b1, a2, b2, fv)
        sys = eq_const(t, T, num_density, a0, c0, a1, b1, a2, b2, fv;
            unit=u"ppb*s^-1" / u"ppb", name=Symbol(:eq_const_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rRO2HO2(t, T, a0, c0, a1)
        sys = rate_RO2HO2(t, T, a0, c0, a1; name=Symbol(:rate_RO2HO2_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rGLYCOH_a(t, T, a0)
        sys = rate_GLYCOH_a(t, T, a0; name=Symbol(:rate_GLYCOH_a_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rGLYCOH_b(t, T, a0)
        sys = rate_GLYCOH_b(t, T, a0; name=Symbol(:rate_GLYCOH_b_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rHACOH_a(t, T, a0, c0)
        sys = rate_HACOH_a(t, T, a0, c0; name=Symbol(:rate_HACOH_a_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rHACOH_b(t, T, a0, c0)
        sys = rate_HACOH_b(t, T, a0, c0; name=Symbol(:rate_HACOH_b_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rDMSOH(t, T, num_density, a0, c0, a1, c1)
        sys = rate_DMSOH(t, T, num_density, a0, c0, a1, c1; name=Symbol(:rate_DMSOH_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rGLYXNO3(t, T, num_density, a0, c0)
        sys = rate_GLYXNO3(t, T, num_density, a0, c0; name=Symbol(:rate_GLYXNO3_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function arrpls(t, T, a0, b0, c0, d0, e0)
        sys = arrplus(t, T, a0, b0, c0, d0, e0; name=Symbol(:arrplus_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function arrplsnodim(t, T, a0, b0, c0, d0, e0)
        sys = arrplus(t, T, a0, b0, c0, d0, e0; unit=u"ppb*s^-1" / u"ppb",
            name=Symbol(:arrplus_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function tunpls(t, T, a0, b0, c0, d0, e0)
        sys = tunplus(t, T, a0, b0, c0, d0, e0; name=Symbol(:tunplus_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function tunplsnodim(t, T, a0, b0, c0, d0, e0)
        sys = tunplus(t, T, a0, b0, c0, d0, e0;
            unit=u"ppb*s^-1" / u"ppb", name=Symbol(:tunplus_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rISO1(t, T, a0, b0, c0, d0, e0, f0, g0)
        sys = rate_ISO1(t, T, a0, b0, c0, d0, e0, f0, g0; name=Symbol(:rate_ISO1_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rISO2(t, T, a0, b0, c0, d0, e0, f0, g0)
        sys = rate_ISO2(t, T, a0, b0, c0, d0, e0, f0, g0; name=Symbol(:rate_ISO2_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rEPO(t, T, num_density, a1, e1, m1)
        sys = rate_EPO(t, T, num_density, a1, e1, m1; name=Symbol(:rate_EPO_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rPAN_abab(t, T, num_density, a0, b0, a1, b1, cf)
        sys = rate_PAN_abab(t, T, num_density, a0, b0, a1, b1, cf; name=Symbol(:rate_PAN_abab_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rPAN_ababnodim(t, T, num_density, a0, b0, a1, b1, cf)
        sys = rate_PAN_abab(t, T, num_density, a0, b0, a1, b1, cf;
            unit=u"ppb*s^-1" / u"ppb", name=Symbol(:rate_PAN_abab_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rPAN_acac(t, T, num_density, a0, c0, a1, c1, cf)
        sys = rate_PAN_acac(t, T, num_density, a0, c0, a1, c1, cf;
            name=Symbol(:rate_PAN_acac_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rNIT(t, T, num_density, a0, b0, c0, n, x0, y0)
        sys = rate_NIT(t, T, num_density, a0, b0, c0, n, x0, y0;
            name=Symbol(:rate_NIT_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end
    function rALK(t, T, num_density, a0, b0, c0, n, x0, y0)
        sys = rate_ALK(t, T, num_density, a0, b0, c0, n, x0, y0; name=Symbol(:rate_ALK_, i))
        i += 1
        push!(rate_systems, sys)
        return sys.k
    end

    rx_sys = @reaction_network GEOSChemGasPhase begin
        # Comment format is:
        # Species   - Molecular formula; full name
        # Equations - Date modified; Reference; Developer initials

        @ivs t [unit = u"s", description = "Time"]
        @species(
            A3O2(t) = 59.469, [unit = u"ppb", description = "CH3CH2CH2OO; Primary RO2 from C3H8"],
            ACET(t) = 87.473, [unit = u"ppb", description = "CH3C(O)CH3; Acetone"],
            ACTA(t) = 70.187, [unit = u"ppb", description = "CH3C(O)OH; Acetic acid"],
            #AERI(t) = 30.013, [unit = u"ppb", description = "I; Dissolved iodine"],
            ALD2(t) = 9.4976, [unit = u"ppb", description = "CH3CHO; Acetaldehyde"],
            ALK4(t) = 8.4827, [unit = u"ppb", description = ">= C4 alkanes"],
            #AONITA(t) = 77.604, [unit = u"ppb", description = "Aerosol-phase organic nitrate from aromatic precursors"],
            AROMRO2(t) = 77.132, [unit = u"ppb", description = "generic peroxy radical from aromatic oxidation"],
            AROMP4(t) = 52.322, [unit = u"ppb", description = "Generic C4 product from aromatic oxidation"],
            AROMP5(t) = 4.8975, [unit = u"ppb", description = "Generic C5 product from aromatic oxidation"],
            ATO2(t) = 88.056, [unit = u"ppb", description = "CH3C(O)CH2O2; RO2 from acetone"],
            ATOOH(t) = 95.672, [unit = u"ppb", description = "CH3C(O)CH2OOH; ATO2 peroxide"],
            B3O2(t) = 44.692, [unit = u"ppb", description = "CH3CH(OO)CH3; Secondary RO2 from C3H8"],
            BALD(t) = 35.859, [unit = u"ppb", description = "benzaldehyde and tolualdehyde"],
            BENZ(t) = 18.166, [unit = u"ppb", description = "C6H6; Benzene"],
            BENZO(t) = 92.094, [unit = u"ppb", description = "C6H5O radical"],
            BENZO2(t) = 45.770, [unit = u"ppb", description = "C6H5O2 radical"],
            BENZP(t) = 45.646, [unit = u"ppb", description = "hydroperoxide from BENZO2"],
            Br(t) = 56.283, [unit = u"ppb", description = "Br; Atomic bromine"],
            Br2(t) = 98.265, [unit = u"ppb", description = "Br2; Molecular bromine"],
            BrCl(t) = 70.858, [unit = u"ppb", description = "BrCl; Bromine chloride"],
            BrNO2(t) = 15.765, [unit = u"ppb", description = "BrNO2; Nitryl bromide"],
            BrNO3(t) = 8.8808, [unit = u"ppb", description = "BrNO3; Bromine nitrate"],
            BrO(t) = 90.871, [unit = u"ppb", description = "BrO; Bromine monoxide"],
            BRO2(t) = 90.134, [unit = u"ppb", description = "C6H5O2 ; Peroxy radical from BENZ oxidation"],
            #BrSALA(t) = 5.6082, [unit = u"ppb", description = "Br; Fine sea salt bromine"],
            #BrSALC(t) = 23.602, [unit = u"ppb", description = "Br; Course sea salt bromine"],
            BZCO3(t) = 34.301, [unit = u"ppb", description = "benzoylperoxy radical"],
            BZCO3H(t) = 18.280, [unit = u"ppb", description = "perbenzoic acid"],
            BZPAN(t) = 25.224, [unit = u"ppb", description = "peroxybenzoyl nitrate"],
            C2H2(t) = 64.023, [unit = u"ppb", description = "C2H2; Ethyne"],
            C2H4(t) = 48.832, [unit = u"ppb", description = "Ethylene"],
            C2H6(t) = 95.500, [unit = u"ppb", description = "C2H6; Ethane"],
            C3H8(t) = 93.457, [unit = u"ppb", description = "C3H8; Propane"],
            C4HVP1(t) = 90.362, [unit = u"ppb", description = "C4 hydroxy-vinyl-peroxy radicals from HPALDs"],
            C4HVP2(t) = 60.972, [unit = u"ppb", description = "C4 hydroxy-vinyl-peroxy radicals from HPALDs"],
            CCl4(t) = 38.249, [unit = u"ppb", description = "CCl4; Carbon tetrachloride"],
            CFC11(t) = 16.012, [unit = u"ppb", description = "CCl3F ; CFC-11, R-11, Freon 11"],
            CFC12(t) = 60.683, [unit = u"ppb", description = "CCl2F2; CFC-12, R-12, Freon 12"],
            CFC113(t) = 99.375, [unit = u"ppb", description = "C2Cl3F3; CFC-113, Freon 113"],
            CFC114(t) = 82.883, [unit = u"ppb", description = "C2Cl2F4; CFC-114, Freon 114"],
            CFC115(t) = 63.617, [unit = u"ppb", description = "C2ClF5; CFC-115, Freon 115"],
            CH2Br2(t) = 53.658, [unit = u"ppb", description = "CH3Br2; Dibromomethane"],
            CH2Cl2(t) = 6.5485, [unit = u"ppb", description = "CH2Cl2; Dichloromethane"],
            CH2I2(t) = 62.463, [unit = u"ppb", description = "CH2I2; Diiodomethane"],
            CH2IBr(t) = 98.965, [unit = u"ppb", description = "CH2IBr; Bromoiodomethane"],
            CH2ICl(t) = 33.656, [unit = u"ppb", description = "CH2ICl; Chloroiodomethane"],
            CH2O(t) = 34.136, [unit = u"ppb", description = "CH2O; Formaldehyde"],
            CH2OO(t) = 72.666, [unit = u"ppb", description = "CH2OO; Criegee intermediate"],
            CH3Br(t) = 76.768, [unit = u"ppb", description = "CH3Br; Methyl bromide"],
            CH3CCl3(t) = 99.842, [unit = u"ppb", description = "CH3CCl3; Methyl chloroform"],
            CH3CHOO(t) = 27.830, [unit = u"ppb", description = "CH3CHOO; Criegee intermediate"],
            CH3Cl(t) = 14.055, [unit = u"ppb", description = "CH3Cl; Chloromethane"],
            CH3I(t) = 19.671, [unit = u"ppb", description = "CH3I; Methyl iodide"],
            CH4(t) = 6.9330, [unit = u"ppb", description = "CH4; Methane"],
            CHBr3(t) = 84.204, [unit = u"ppb", description = "CHBr3; Tribromethane"],
            CHCl3(t) = 79.435, [unit = u"ppb", description = "CHCl3; Chloroform"],
            Cl(t) = 47.484, [unit = u"ppb", description = "Cl; Atomic chlorine"],
            Cl2(t) = 71.292, [unit = u"ppb", description = "Cl2; Molecular chlorine"],
            Cl2O2(t) = 96.248, [unit = u"ppb", description = "Cl2O2; Dichlorine dioxide"],
            ClNO2(t) = 36.482, [unit = u"ppb", description = "ClNO2; Nitryl chloride"],
            ClNO3(t) = 55.020, [unit = u"ppb", description = "ClONO2; Chlorine nitrate"],
            ClO(t) = 0.9863, [unit = u"ppb", description = "ClO; Chlorine monoxide"],
            ClOO(t) = 53.802, [unit = u"ppb", description = "ClOO; Chlorine dioxide"],
            CO(t) = 32.978, [unit = u"ppb", description = "CO; Carbon monoxide"],
            CO2(t) = 55.611, [unit = u"ppb", description = "CO2; Carbon dioxide"],
            CSL(t) = 23.332, [unit = u"ppb", description = "cresols and xylols"],
            DMS(t) = 5.8615, [unit = u"ppb", description = "(CH3)2S; Dimethylsulfide"],
            EOH(t) = 18.988, [unit = u"ppb", description = "C2H5OH; Ethanol"],
            ETHLN(t) = 69.851, [unit = u"ppb", description = "CHOCH2ONO2; Ethanal nitrate"],
            ETHN(t) = 24.767, [unit = u"ppb", description = "stable hydroxy-nitrooxy-ethane"],
            ETHP(t) = 59.843, [unit = u"ppb", description = "stable hydroxy-hydroperoxy-ethane"],
            ETNO3(t) = 11.558, [unit = u"ppb", description = "C2H5ONO2; Ethyl nitrate"],
            ETO(t) = 85.116, [unit = u"ppb", description = "hydroxy-alkoxy-ethane radical"],
            ETOO(t) = 24.842, [unit = u"ppb", description = "hydroxy-peroxy-ethane radical, formed from ethene + OH"],
            ETO2(t) = 29.300, [unit = u"ppb", description = "CH3CH2OO; Ethylperoxy radical"],
            ETP(t) = 81.950, [unit = u"ppb", description = "CH3CH2OOH; Ethylhydroperoxide"],
            GLYC(t) = 20.407, [unit = u"ppb", description = "HOCH2CHO; Glycoaldehyde"],
            GLYX(t) = 81.571, [unit = u"ppb", description = "CHOCHO; Glyoxal"],
            H(t) = 7.6410, [unit = u"ppb", description = "H; Atomic hydrogen"],
            H1211(t) = 46.027, [unit = u"ppb", description = "CBrClF2; H-1211"],
            H1301(t) = 36.584, [unit = u"ppb", description = "CBrF3; H-1301"],
            H2402(t) = 52.639, [unit = u"ppb", description = "C2Br2F4; H-2402"],
            H2O(t) = 56.623, [unit = u"ppb", description = "H2O; Water vapor"],
            H2O2(t) = 44.325, [unit = u"ppb", description = "H2O2; Hydrogen peroxide"],
            HAC(t) = 88.335, [unit = u"ppb", description = "HOCH2C(O)CH3; Hydroxyacetone"],
            HBr(t) = 15.715, [unit = u"ppb", description = "HBr; Hypobromic acid"],
            HC5A(t) = 44.638, [unit = u"ppb", description = "C5H8O2; Isoprene-4,1-hydroxyaldehyde"],
            HCFC123(t) = 22.296, [unit = u"ppb", description = "C2HCl2F3; HCFC-123, R-123, Freon 123"],
            HCFC141b(t) = 93.091, [unit = u"ppb", description = "C(CH3)Cl2F; HCFC-141b, R-141b, Freon 141b"],
            HCFC142b(t) = 81.415, [unit = u"ppb", description = "C(CH3)ClF2; HCFC-142b, R-142b, Freon 142b"],
            HCFC22(t) = 16.333, [unit = u"ppb", description = "CHClF2 ; HCFC-22, R-22, Freon 22"],
            HCl(t) = 87.644, [unit = u"ppb", description = "HCl; Hydrochloric acid"],
            HCOOH(t) = 99.489, [unit = u"ppb", description = "HCOOH; Formic acid"],
            HI(t) = 68.377, [unit = u"ppb", description = "HI; Hydrogen iodide"],
            HMHP(t) = 16.437, [unit = u"ppb", description = "HOCH2OOH; Hydroxymethyl hydroperoxide"],
            HMML(t) = 33.654, [unit = u"ppb", description = "C4H6O3; Hydroxymethyl-methyl-a-lactone"],
            HMS(t) = 54.099, [unit = u"ppb", description = "HOCH2SO3-; hydroxymethanesulfonate"],
            HNO2(t) = 18.174, [unit = u"ppb", description = "HONO; Nitrous acid"],
            HNO3(t) = 62.170, [unit = u"ppb", description = "HNO3; Nitric acid"],
            HNO4(t) = 39.302, [unit = u"ppb", description = "HNO4; Pernitric acid"],
            HO2(t) = 58.306, [unit = u"ppb", description = "HO2; Hydroperoxyl radical"],
            HOBr(t) = 59.529, [unit = u"ppb", description = "HOBr; Hypobromous acid"],
            HOCl(t) = 90.397, [unit = u"ppb", description = "HOCl; Hypochlorous acid"],
            HOI(t) = 35.827, [unit = u"ppb", description = "HOI; Hypoiodous acid"],
            HONIT(t) = 58.760, [unit = u"ppb", description = "2nd gen monoterpene organic nitrate"],
            HPALD1(t) = 73.164, [unit = u"ppb", description = "O=CHC(CH3)=CHCH2OOH; d-4,1-C5-hydroperoxyaldehyde"],
            HPALD1OO(t) = 34.759, [unit = u"ppb", description = "peroxy radicals from HPALD1"],
            HPALD2(t) = 1.0291, [unit = u"ppb", description = "HOOCH2C(CH3)=CHCH=O; d-1,4-C5-hydroperoxyaldehyde"],
            HPALD2OO(t) = 28.332, [unit = u"ppb", description = "peroxy radicals from HPALD2"],
            HPALD3(t) = 91.404, [unit = u"ppb", description = "O=CHC(CH3)OOHCH=CH2; b-2,1-C5-hydroperoxyaldehyde"],
            HPALD4(t) = 14.949, [unit = u"ppb", description = "CH2=C(CH3)CHOOHCH=O; b-3,4-C5-hydroperoxyaldehyde"],
            HPETHNL(t) = 90.330, [unit = u"ppb", description = "CHOCH2OOH; hydroperoxyethanal"],
            I(t) = 85.092, [unit = u"ppb", description = "I; Atmoic iodine"],
            I2(t) = 49.639, [unit = u"ppb", description = "I2; Molecular iodine"],
            I2O2(t) = 70.687, [unit = u"ppb", description = "I2O2; Diiodine dioxide"],
            I2O3(t) = 85.950, [unit = u"ppb", description = "I2O3; Diiodine trioxide"],
            I2O4(t) = 36.304, [unit = u"ppb", description = "I2O4; Diiodine tetraoxide"],
            IBr(t) = 46.529, [unit = u"ppb", description = "IBr; Iodine monobromide"],
            ICHE(t) = 17.739, [unit = u"ppb", description = "C5H8O3; Isoprene hydroxy-carbonyl-epoxides"],
            ICHOO(t) = 75.448, [unit = u"ppb", description = "peroxy radical from IEPOXD"],
            ICl(t) = 32.124, [unit = u"ppb", description = "ICl; Iodine monochloride"],
            ICN(t) = 22.128, [unit = u"ppb", description = "C5H7NO4; Lumped isoprene carbonyl nitrates"],
            ICNOO(t) = 12.649, [unit = u"ppb", description = "peroxy radicals from ICN"],
            ICPDH(t) = 71.736, [unit = u"ppb", description = "C5H10O5; Isoprene dihydroxy hydroperoxycarbonyl"],
            IDC(t) = 28.840, [unit = u"ppb", description = "C5H6O2; Lumped isoprene dicarbonyls"],
            IDCHP(t) = 2.3035, [unit = u"ppb", description = "C5H8O5; Isoprene dicarbonyl hydroxy dihydroperoxide"],
            IDHDP(t) = 21.649, [unit = u"ppb", description = "C5H12O6; Isoprene dihydroxy dihydroperoxide"],
            IDHNBOO(t) = 25.086, [unit = u"ppb", description = "peroxy radicals from INPB"],
            IDHNDOO1(t) = 15.276, [unit = u"ppb", description = "peroxy radicals from INPD"],
            IDHNDOO2(t) = 61.961, [unit = u"ppb", description = "peroxy radicals from INPD"],
            IDHPE(t) = 94.216, [unit = u"ppb", description = "C5H10O5; Isoprene dihydroxy hydroperoxy epoxide"],
            IDN(t) = 13.328, [unit = u"ppb", description = "C5H8N2O6; Lumped isoprene dinitrates"],
            IDNOO(t) = 44.188, [unit = u"ppb", description = "peroxy radicals from IDN"],
            IEPOXA(t) = 60.554, [unit = u"ppb", description = "C5H10O3; trans-Beta isoprene epoxydiol"],
            IEPOXAOO(t) = 27.948, [unit = u"ppb", description = "peroxy radical from trans-Beta isoprene epoxydiol"],
            IEPOXB(t) = 38.609, [unit = u"ppb", description = "C5H10O3; cis-Beta isoprene epoxydiol"],
            IEPOXBOO(t) = 87.094, [unit = u"ppb", description = "peroxy radical from cis-Beta isoprene epoxydiol"],
            IEPOXD(t) = 8.0346, [unit = u"ppb", description = "C5H10O3; Delta isoprene epoxydiol"],
            IHN1(t) = 11.028, [unit = u"ppb", description = "C5H9NO4; Isoprene-d-4-hydroxy-1-nitrate"],
            IHN2(t) = 71.286, [unit = u"ppb", description = "C5H9NO4; Isoprene-b-1-hydroxy-2-nitrate"],
            IHN3(t) = 82.443, [unit = u"ppb", description = "C5H9NO4; Isoprene-b-4-hydroxy-3-nitrate"],
            IHN4(t) = 91.975, [unit = u"ppb", description = "C5H9NO4; Isoprene-d-1-hydroxy-4-nitrate"],
            IHOO1(t) = 27.517, [unit = u"ppb", description = "peroxy radical from OH addition to isoprene at C1"],
            IHOO4(t) = 51.624, [unit = u"ppb", description = "peroxy radical from OH addition to isoprene at C4"],
            IHPNBOO(t) = 66.363, [unit = u"ppb", description = "peroxy radicals from INPB"],
            IHPNDOO(t) = 63.447, [unit = u"ppb", description = "peroxy radicals from INPD"],
            IHPOO1(t) = 79.476, [unit = u"ppb", description = "peroxy radical from ISOPOOH"],
            IHPOO2(t) = 58.328, [unit = u"ppb", description = "peroxy radical from ISOPOOH"],
            IHPOO3(t) = 5.9515, [unit = u"ppb", description = "peroxy radical from ISOPOOH"],
            INA(t) = 26.490, [unit = u"ppb", description = "alkoxy radical from INO2D"],
            INDIOL(t) = 52.424, [unit = u"ppb", description = "Generic aerosol phase organonitrate hydrolysis product"],
            INO(t) = 73.106, [unit = u"ppb", description = "INO; Nitrosyl iodide"],
            INO2B(t) = 25.723, [unit = u"ppb", description = "beta-peroxy radicals from isoprene + NO3"],
            INO2D(t) = 30.524, [unit = u"ppb", description = "delta-peroxy radicals from isoprene + NO3"],
            INPB(t) = 96.224, [unit = u"ppb", description = "C5H9NO5; Lumped isoprene beta-hydroperoxy nitrates"],
            INPD(t) = 80.193, [unit = u"ppb", description = "C5H9NO5; Lumped isoprene delta-hydroperoxy nitrates"],
            IO(t) = 26.615, [unit = u"ppb", description = "IO; Iodine monoxide"],
            IONITA(t) = 25.539, [unit = u"ppb", description = "Aerosol-phase organic nitrate from isoprene precursors"],
            IONO(t) = 56.756, [unit = u"ppb", description = "IONO; Nitryl iodide"],
            IONO2(t) = 53.992, [unit = u"ppb", description = "IONO2; Iodine nitrate"],
            IPRNO3(t) = 22.638, [unit = u"ppb", description = "C3H8ONO2; Isopropyl nitrate"],
            #ISALA(t) = 79.006, [unit = u"ppb", description = "I; Fine sea salt iodine"],
            #ISALC(t) = 13.184, [unit = u"ppb", description = "I; Coarse sea salt iodine"],
            ISOP(t) = 83.240, [unit = u"ppb", description = "CH2=C(CH3)CH=CH2; Isoprene"],
            ISOPNOO1(t) = 60.043, [unit = u"ppb", description = "peroxy radicals from IHN2"],
            ISOPNOO2(t) = 71.160, [unit = u"ppb", description = "peroxy radicals from IHN3"],
            ITCN(t) = 26.206, [unit = u"ppb", description = "C5H9NO7; Lumped tetrafunctional isoprene carbonyl-nitrates"],
            ITHN(t) = 11.572, [unit = u"ppb", description = "C5H11NO7; Lumped tetrafunctional isoprene hydroxynitrates"],
            KO2(t) = 84.868, [unit = u"ppb", description = "RO2 from >3 ketones"],
            LBRO2H(t) = 61.633, [unit = u"ppb", description = "Dummy spc to track oxidation of BRO2 by HO2"],
            LBRO2N(t) = 11.255, [unit = u"ppb", description = "Dummy spc to track oxidation of BRO2 by NO"],
            LIMO(t) = 53.446, [unit = u"ppb", description = "C10H16; Limonene"],
            LIMO2(t) = 92.794, [unit = u"ppb", description = "RO2 from LIMO"],
            LISOPOH(t) = 92.199, [unit = u"ppb", description = "Dummy spc to track oxidation of ISOP by OH"],
            LISOPNO3(t) = 29.918, [unit = u"ppb", description = "Dummy spc to track oxidation of ISOP by NO3"],
            LNRO2H(t) = 29.939, [unit = u"ppb", description = "Dummy spc to track oxidation of NRO2 by HO2"],
            LNRO2N(t) = 57.661, [unit = u"ppb", description = "Dummy spc to track oxidation of NRO2 by NO"],
            LTRO2H(t) = 46.441, [unit = u"ppb", description = "Dummy spc to track oxidation of TRO2 by HO2"],
            LTRO2N(t) = 28.826, [unit = u"ppb", description = "Dummy spc to track oxidation of TRO2 by NO"],
            LVOC(t) = 47.079, [unit = u"ppb", description = "C5H14O5; Gas-phase low-volatility non-IEPOX product of ISOPOOH (RIP) oxidation"],
            #LVOCOA(t) = 97.788, [unit = u"ppb", description = "C5H14O5; Aerosol-phase low-volatility non-IEPOX product of ISOPOOH (RIP) oxidation"],
            LXRO2H(t) = 62.462, [unit = u"ppb", description = "Dummy spc to track oxidation of XRO2 by HO2"],
            LXRO2N(t) = 39.337, [unit = u"ppb", description = "Dummy spc to track oxidation of XRO2 by NO"],
            MACR(t) = 17.935, [unit = u"ppb", description = "CH2=C(CH3)CHO; Methacrolein"],
            MACR1OO(t) = 91.114, [unit = u"ppb", description = "peroxyacyl radical from MACR + OH"],
            MACR1OOH(t) = 61.547, [unit = u"ppb", description = "CH2=C(CH3)C(O)OOH; Peracid from MACR"],
            MACRNO2(t) = 97.585, [unit = u"ppb", description = "Product of MCRHN + OH"],
            MAP(t) = 61.421, [unit = u"ppb", description = "CH3C(O)OOH; Peroxyacetic acid"],
            MCO3(t) = 0.6220, [unit = u"ppb", description = "CH3C(O)OO; Peroxyacetyl radical"],
            MCRDH(t) = 23.756, [unit = u"ppb", description = "C4H8O3; Dihydroxy-MACR"],
            MCRENOL(t) = 43.151, [unit = u"ppb", description = "C4H6O2; Lumped enols from MVK/MACR"],
            MCRHN(t) = 5.5067, [unit = u"ppb", description = "HOCH2C(ONO2)(CH3)CHO; Hydroxynitrate from MACR"],
            MCRHNB(t) = 2.8277, [unit = u"ppb", description = "O2NOCH2C(OH)(CH3)CHO; Hydroxynitrate from MACR"],
            MCRHP(t) = 0.2048, [unit = u"ppb", description = "HOCH2C(OOH)(CH3)CHO; Hydroxy-hydroperoxy-MACR"],
            MCROHOO(t) = 49.422, [unit = u"ppb", description = "peroxy radical from MACR + OH"],
            MCT(t) = 77.753, [unit = u"ppb", description = "methylcatechols"],
            MEK(t) = 61.732, [unit = u"ppb", description = "RC(O)R; Methyl ethyl ketone"],
            MENO3(t) = 89.773, [unit = u"ppb", description = "CH3ONO2; methyl nitrate"],
            MGLY(t) = 93.209, [unit = u"ppb", description = "CH3COCHO; Methylglyoxal"],
            MO2(t) = 73.881, [unit = u"ppb", description = "CH3O2; Methylperoxy radical"],
            MOH(t) = 55.585, [unit = u"ppb", description = "CH3OH; Methanol"],
            MONITA(t) = 61.052, [unit = u"ppb", description = "Aerosol-phase organic nitrate from monoterpene precursors"],
            MONITS(t) = 71.179, [unit = u"ppb", description = "Saturated 1st gen monoterpene organic nitrate"],
            MONITU(t) = 23.294, [unit = u"ppb", description = "Unsaturated 1st gen monoterpene organic nitrate"],
            MP(t) = 25.386, [unit = u"ppb", description = "CH3OOH; Methylhydroperoxide"],
            MPAN(t) = 89.575, [unit = u"ppb", description = "CH2=C(CH3)C(O)OONO2; Peroxymethacroyl nitrate (PMN)"],
            MPN(t) = 59.399, [unit = u"ppb", description = "CH3O2NO2; Methyl peroxy nitrate"],
            MSA(t) = 82.711, [unit = u"ppb", description = "CH4SO3; Methanesulfonic acid"],
            MTPA(t) = 25.413, [unit = u"ppb", description = "Lumped monoterpenes: a-pinene, b-pinene, sabinene, carene"],
            MTPO(t) = 8.9141, [unit = u"ppb", description = "Other monoterpenes: Terpinene, terpinolene, myrcene, ocimene, other monoterpenes"],
            MVK(t) = 62.109, [unit = u"ppb", description = "CH2=CHC(=O)CH3; Methyl vinyl ketone"],
            MVKDH(t) = 97.635, [unit = u"ppb", description = "HOCH2CH2OHC(O)CH3; Dihydroxy-MVK"],
            MVKHC(t) = 98.159, [unit = u"ppb", description = "C4H6O3; MVK hydroxy-carbonyl"],
            MVKHCB(t) = 50.621, [unit = u"ppb", description = "C4H6O3; MVK hydroxy-carbonyl"],
            MVKHP(t) = 35.334, [unit = u"ppb", description = "C4H8O4; MVK hydroxy-hydroperoxide"],
            MVKN(t) = 48.570, [unit = u"ppb", description = "HOCH2CH(ONO2)C(=O)CH3; Hydroxynitrate from MVK"],
            MVKOHOO(t) = 41.490, [unit = u"ppb", description = "peroxy radical from MVK + OH"],
            MVKPC(t) = 8.2590, [unit = u"ppb", description = "OCHCH(OOH)C(O)CH3; MVK hydroperoxy-carbonyl"],
            N(t) = 56.891, [unit = u"ppb", description = "N; Atomic nitrogen"],
            N2O(t) = 50.751, [unit = u"ppb", description = "N2O; Nitrous oxide"],
            N2O5(t) = 19.118, [unit = u"ppb", description = "N2O5; Dinitrogen pentoxide"],
            NAP(t) = 33.270, [unit = u"ppb", description = "C10H8; Naphthalene; IVOC surrogate"],
            NIT(t) = 77.088, [unit = u"ppb", description = "NIT; Fine mode inorganic nitrate"],
            NITs(t) = 56.341, [unit = u"ppb", description = "NITs; Coarse mode inorganic nitrate"],
            NO(t) = 54.454, [unit = u"ppb", description = "NO; Nitric oxide"],
            NO2(t) = 5.5605, [unit = u"ppb", description = "NO2; Nitrogen dioxide"],
            NO3(t) = 31.739, [unit = u"ppb", description = "NO3; Nitrate radical"],
            NPHEN(t) = 87.002, [unit = u"ppb", description = "nitrophenols"],
            NPRNO3(t) = 58.488, [unit = u"ppb", description = "C3H8ONO2; n-propyl nitrate"],
            NRO2(t) = 83.618, [unit = u"ppb", description = "Peroxy radical from NAP oxidation"],
            O(t) = 14.164, [unit = u"ppb", description = "O(3P); Ground state atomic oxygen"],
            O1D(t) = 17.648, [unit = u"ppb", description = "O(1D); Excited atomic oxygen"],
            O3(t) = 2.1326, [unit = u"ppb", description = "O3; Ozone"],
            #O3A(t) = 28.223, [unit = u"ppb", description = "O3; Ozone in accum seasalt"],
            #O3C(t) = 49.728, [unit = u"ppb", description = "O3; Ozone in coarse seasalt"],
            OClO(t) = 60.037, [unit = u"ppb", description = "OClO; Chlorine dioxide"],
            OCS(t) = 20.867, [unit = u"ppb", description = "COS; Carbonyl sulfide"],
            OH(t) = 19.293, [unit = u"ppb", description = "OH; Hydroxyl radical"],
            OIO(t) = 51.406, [unit = u"ppb", description = "OIO; Iodine dioxide"],
            OLND(t) = 63.247, [unit = u"ppb", description = "Monoterpene-derived NO3-alkene adduct"],
            OLNN(t) = 95.667, [unit = u"ppb", description = "Monoterpene-derived NO3 adduct"],
            OTHRO2(t) = 74.722, [unit = u"ppb", description = "Other C2 RO2 not from C2H6 oxidation"],
            PAN(t) = 59.945, [unit = u"ppb", description = "CH3C(O)OONO2; Peroxyacetylnitrate"],
            PHEN(t) = 59.463, [unit = u"ppb", description = "phenol"],
            PIO2(t) = 94.953, [unit = u"ppb", description = "RO2 from MTPA"],
            PIP(t) = 96.809, [unit = u"ppb", description = "Peroxides from MTPA"],
            PO2(t) = 25.827, [unit = u"ppb", description = "HOCH2CH(OO)CH3; RO2 from propene"],
            PP(t) = 61.804, [unit = u"ppb", description = "HOCH2CH(OOH)CH3; Peroxide from PO2"],
            PPN(t) = 94.436, [unit = u"ppb", description = "CH3CH2C(O)OONO2; Peroxypropionylnitrate"],
            PRN1(t) = 87.618, [unit = u"ppb", description = "O2NOCH2CH(OO)CH3; RO2 from propene + NO3"],
            PROPNN(t) = 33.842, [unit = u"ppb", description = "CH3C(=O)CH2ONO2; Propanone nitrate"],
            PRPE(t) = 63.131, [unit = u"ppb", description = "C3H6; >= C3 alkenes"],
            PRPN(t) = 43.083, [unit = u"ppb", description = "O2NOCH2CH(OOH)CH3; Peroxide from PRN1"],
            PYAC(t) = 79.104, [unit = u"ppb", description = "CH3COCOOH; Pyruvic acid"],
            R4N1(t) = 37.438, [unit = u"ppb", description = "RO2 from R4N2"],
            R4N2(t) = 46.153, [unit = u"ppb", description = "RO2NO; >= C4 alkylnitrates"],
            R4O2(t) = 36.162, [unit = u"ppb", description = "RO2 from ALK4"],
            R4P(t) = 40.804, [unit = u"ppb", description = "CH3CH2CH2CH2OOH; Peroxide from R4O2"],
            RA3P(t) = 57.025, [unit = u"ppb", description = "CH3CH2CH2OOH; Peroxide from A3O2"],
            RB3P(t) = 72.421, [unit = u"ppb", description = "CH3CH(OOH)CH3; Peroxide from B3O2"],
            RCHO(t) = 37.680, [unit = u"ppb", description = "CH3CH2CHO; >= C3 aldehydes"],
            RCO3(t) = 63.108, [unit = u"ppb", description = "CH3CH2C(O)OO; Peroxypropionyl radical"],
            RIPA(t) = 42.881, [unit = u"ppb", description = "HOCH2C(OOH)(CH3)CH=CH2; 1,2-ISOPOOH"],
            RIPB(t) = 24.737, [unit = u"ppb", description = "HOCH2C(OOH)(CH3)CH=CH2; 4,3-ISOPOOH"],
            RIPC(t) = 57.374, [unit = u"ppb", description = "C5H10O3; d(1,4)-ISOPOOH"],
            RIPD(t) = 5.9941, [unit = u"ppb", description = "C5H10O3; d(4,1)-ISOPOOH"],
            ROH(t) = 67.750, [unit = u"ppb", description = "C3H7OH; > C2 alcohols"],
            RP(t) = 38.598, [unit = u"ppb", description = "CH3CH2C(O)OOH; Peroxide from RCO3"],
            SALAAL(t) = 37.928, [unit = u"ppb", description = "Accumulation mode seasalt aerosol alkalinity"],
            SALCAL(t) = 13.651, [unit = u"ppb", description = "Coarse mode seasalt aerosol alkalinity"],
            SALACL(t) = 41.613, [unit = u"ppb", description = "Cl; Fine chloride"],
            SALCCL(t) = 41.813, [unit = u"ppb", description = "Cl; Coarse chloride"],
            #SALASO2(t) = 46.971, [unit = u"ppb", description = "SO2; Fine seasalt"],
            #SALCSO2(t) = 98.028, [unit = u"ppb", description = "SO2; Coarse seasalt"],
            #SALASO3(t) = 18.107, [unit = u"ppb", description = "SO3--; Fine seasalt"],
            #SALCSO3(t) = 51.008, [unit = u"ppb", description = "SO3--; Coarse chloride"],
            SO2(t) = 6.2180, [unit = u"ppb", description = "SO2; Sulfur dioxide"],
            SO4(t) = 97.494, [unit = u"ppb", description = "SO4; Sulfate"],
            SO4s(t) = 74.244, [unit = u"ppb", description = "SO4 on sea-salt; Sulfate"],
            #SOAGX(t) = 82.469, [unit = u"ppb", description = "CHOCHO; Aerosol-phase glyoxal"],
            #SOAIE(t) = 81.289, [unit = u"ppb", description = "C5H10O3; Aerosol-phase IEPOX"],
            TOLU(t) = 58.265, [unit = u"ppb", description = "C7H8; Toluene"],
            TRO2(t) = 10.659, [unit = u"ppb", description = "Peroxy radical from TOLU oxidation"],
            XYLE(t) = 24.563, [unit = u"ppb", description = "C8H10; Xylene"],
            XRO2(t) = 25.058, [unit = u"ppb", description = "Peroxy radical from XYLE oxidation"],
            H2(t) = 97.924, [unit = u"ppb", description = "H2; Molecular hydrogen"],
            N2(t) = 36.998, [unit = u"ppb", description = "N2; Molecular nitrogen"],
            O2(t) = 51.201, [unit = u"ppb", description = "O2; Molecular oxygen"],
            RCOOH(t) = 12.892, [unit = u"ppb", description = "C2H5C(O)OH; > C2 organic acids "],
        )

        @parameters(
            T = 298.15, [unit = u"K", description = "Temperature"],
            num_density = 2.7e19, [description = "Number density of air (The units should be molecules/cm^3 but the equations here treat it as unitless)."],
            k_mt1 = 0, [unit = u"ppb^-1*s^-1", description = "Atmospheric molecular density (i.e. pressure)"],
            k_mt2 = 0, [unit = u"ppb^-1*s^-1", description = "Seasalt rate constant"],
            k_mt3 = 0, [unit = u"ppb^-1*s^-1", description = "Seasalt rate constant"],
            k_mt4 = 0, [unit = u"ppb^-1*s^-1", description = "Seasalt rate constant"],
            k_mt5 = 0, [unit = u"ppb^-1*s^-1", description = "Seasalt rate constant"],
            k_mt6 = 0, [unit = u"ppb^-1*s^-1", description = "Seasalt rate constant"],
            k_cld1 = 0, [unit = u"ppb^-1*s^-1", description = "Cloud rate constant"],
            k_cld2 = 0, [unit = u"ppb^-1*s^-1", description = "Cloud rate constant"],
            k_cld3 = 0, [unit = u"s^-1", description = "Cloud rate constant"],
            k_cld4 = 0, [unit = u"ppb^-1*s^-1", description = "HMS rate constant"],
            k_cld5 = 0, [unit = u"s^-1", description = "HMS rate constant"],
            k_cld6 = 0, [unit = u"ppb^-1*s^-1", description = "HMS rate constant"],
            j_2 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_3 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_1 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_11 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_9 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_10 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_7 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_8 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_16 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_15 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_17 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_18 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_12 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_13 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_14 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_61 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_62 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_59 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_70 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_76 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_77 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_69 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_68 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_72 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_73 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_74 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_71 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_63 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_64 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_65 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_66 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_75 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_79 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_80 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_81 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_82 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_83 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_84 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_85 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_98 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_99 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_23 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_28 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_32 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_29 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_30 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_31 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_56 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_55 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_50 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_43 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_45 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_33 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_22 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_27 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_25 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_26 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_21 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_19 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_20 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_24 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_44 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_42 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_37 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_38 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_39 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_40 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_41 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_47 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_48 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_49 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_46 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_53 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_51 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_54 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_101 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_114 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_115 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_116 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_117 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_118 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_119 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_120 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_121 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_122 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_123 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_124 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_125 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_126 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_127 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_128 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_129 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_103 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_104 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_97 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_36 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_34 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_100 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_6 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_105 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_107 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_111 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_112 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_113 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_130 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_131 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_132 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_133 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_134 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_135 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_136 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_137 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_86 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_87 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_88 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_89 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_90 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_91 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_92 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_93 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_94 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_95 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_96 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_108 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_109 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_110 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_138 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_139 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_140 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_141 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_142 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_143 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_144 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_145 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_146 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_147 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_148 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_149 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_150 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_151 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_152 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_106 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_78 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_153 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_154 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_155 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_156 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_157 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_158 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_159 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_160 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_161 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_162 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_163 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_164 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_165 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
            j_166 = 0, [unit = u"s^-1", description = "Photolysis rate constant"],
        )

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # %%%%% Reactions extracted from sulfate_mod.F90 (MSL, BMY)             %%%%%
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #==
        NOTE(CT): KPP uses minus signs to denote species that are consumed in a reaction but do not
        effect the reaction rate, e.g.:
            SO2  + SALAAL + O3  = SO4 - SALAAL : K_MT(1);
        (as described here: https://kpp.readthedocs.io/en/stable/using_kpp/04_input_for_kpp.html?highlight=minus#equations)
        To get the same effect here, we divide the reaction rate constant by the species that was
        subtracted in KPP, e.g.:
            k_mt1 / SALAAL, SO2  + SALAAL + O3 --> SO4
        ==#
        #
        # Seasalt
        k_mt1 / SALAAL, SO2 + SALAAL + O3 --> SO4
        k_mt2, HCl + SALAAL --> SALACL
        k_mt3, HNO3 + SALAAL --> NIT
        k_mt4 / SALCAL, SO2 + SALCAL + O3 --> SO4s
        k_mt5, HCl + SALCAL --> SALCCL
        k_mt6, HNO3 + SALCAL --> NITs
        #
        # Cloud
        # S(IV) --> S(VI)
        k_cld1, SO2 + H2O2 --> SO4
        k_cld2, SO2 + O3 --> SO4
        k_cld3, SO2 --> SO4 #==Mn & Fe catalysis + HET_DROP_CHEM()==#
        #
        # HMS
        k_cld4, CH2O + SO2 --> HMS #==Sep 2021; Moch2020; MSL==#
        k_cld5, HMS --> SO2 + CH2O #==Sep 2021; Moch2020; MSL==#
        k_cld6 / SO2, HMS + OH + SO2 --> 2SO4 + CH2O #==Sep 2021; Moch2020; MSL==#
        #
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # %%%%% Gas-phase chemistry reactions                                   %%%%%
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        arr(t, T, 3.00e-12, 0.0, -1500.0e0), O3 + NO --> NO2 + O2
        arr(t, T, 1.70e-12, 0.0, -940.0e0), O3 + OH --> HO2 + O2
        arr(t, T, 1.00e-14, 0.0, -490.0e0), O3 + HO2 --> OH + O2 + O2
        arr(t, T, 1.20e-13, 0.0, -2450.0e0), O3 + NO2 --> O2 + NO3
        arr(t, T, 2.90e-16, 0.0, -1000.0e0), O3 + MO2 --> CH2O + HO2 + O2 #==2014/02/03; Eastham2014; SDE==#
        c(1.80e-12), OH + OH --> H2O + O #==2014/02/03; Eastham2014; SDE==#
        arr3(t, T, num_density, 6.90e-31, 1.0e+00, 0.0, 2.6e-11, 0.0, 0.0, 0.6e0), OH + OH --> H2O2 #==+M==#
        arr(t, T, 4.80e-11, 0.0, 250.0e0), OH + HO2 --> H2O + O2
        c(1.80e-12), OH + H2O2 --> H2O + HO2
        arr(t, T, 3.30e-12, 0.0, 270.0e0), HO2 + NO --> OH + NO2 #==2013/02/12; JPL 10-6; BHH,JMAO,EAM==#
        rHO2HO2(t, T, num_density, H2O, 3.00e-13, 460.0e0, 2.1e-33, 920.0e0), HO2 + HO2 --> H2O2 + O2 #==2014/02/03; Eastham2014; SDE==#
        rOHCO(t, T, num_density), OH + CO --> HO2 + CO2 #==2017/02/22; JPL 15-10; BHH,MJE==#
        arr(t, T, 2.45e-12, 0.0, -1775.0e0), OH + CH4 --> MO2 + H2O
        rRO2NOb1(t, T, 2.80e-12, 300.0e0), MO2 + NO --> CH2O + HO2 + NO2 #==2019/05/10; Fisher2018; JAF==#
        rRO2NOa1(t, T, 2.80e-12, 300.0e0), MO2 + NO --> MENO3 #==2019/05/10; Fisher2018; JAF==#
        arr(t, T, 4.10e-13, 0.0e0, 750.0e0), MO2 + HO2 --> MP + O2
        tbrnch(t, T, 9.50e-14, 0, 390.0e0, 2.62e1, 0, -1130.0e0), MO2 + MO2 --> MOH + CH2O + O2
        tbrnch(t, T, 9.50e-14, 0, 390.0e0, 4.0e-2, 0, 1130.0e0), MO2 + MO2 --> 2.000CH2O + 2.000HO2
        c(1.60e-10), MO2 + OH --> 0.13MOH + 0.87CH2O + 1.74HO2 #==2021/09/22; Bates2021a; KHB,MSL==#
        arr(t, T, 2.66e-12, 0.0, 200.0e0), MP + OH --> MO2 + H2O
        arr(t, T, 1.14e-12, 0.0, 200.0e0), MP + OH --> CH2O + OH + H2O
        arr(t, T, 2.66e-12, 0.0, 200.0e0), ATOOH + OH --> ATO2 + H2O #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        arr(t, T, 1.14e-12, 0.0, 200.0e0), ATOOH + OH --> MGLY + OH + H2O #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        arr(t, T, 5.50e-12, 0.0, 125.0e0), CH2O + OH --> CO + HO2 + H2O
        arr3(t, T, num_density, 1.80e-30, 3.0e+00, 0.0, 2.8e-11, 0.0, 0.0, 0.6e0), NO2 + OH --> HNO3 #==+M==#
        rOHHNO3(t, T, num_density, 2.41e-14, 460.0e0, 2.69e-17, 2199.0e0, 6.51e-34, 1335.0e0), HNO3 + OH --> H2O + NO3
        arr3(t, T, num_density, 7.00e-31, 2.6e+00, 0.0, 3.60e-11, 0.1e0, 0.0, 0.6e0), NO + OH --> HNO2 #==+M==#
        arr(t, T, 1.80e-11, 0.0, -390.0e0), HNO2 + OH --> H2O + NO2
        arr3(t, T, num_density, 1.90e-31, 3.4e+00, 0.0, 4.0e-12, 0.3e0, 0.0, 0.6e0), HO2 + NO2 --> HNO4 #==2017/02/22; JPL 15-10; BHH,MJE==#
        arr3nodim(t, T, num_density, 9.05e-05, 3.4e0, -10900.0e0, 1.90e15, 0.3e0, -10900.0e0, 0.6e0), HNO4 --> HO2 + NO2 #==2017/02/22; JPL 15-10; BHH,MJE==#
        arr(t, T, 1.30e-12, 0.0, 380.0e0), HNO4 + OH --> H2O + NO2 + O2
        c(3.50e-12), HO2 + NO3 --> OH + NO2 + O2
        arr(t, T, 1.50e-11, 0.0, 170.0e0), NO + NO3 --> 2.000NO2
        c(2.20e-11), OH + NO3 --> HO2 + NO2
        arr3(t, T, num_density, 2.40e-30, 3.0e+00, 0.0, 1.6e-12, -0.1e0, 0.0, 0.6e0), NO2 + NO3 --> N2O5 #==2017/02/22; JPL 15-10; BHH,MJE==#
        arr3nodim(t, T, num_density, 4.14e-04, 3.0e0, -10840.0e0, 2.76e14, -0.1e0, -10840.0e0, 0.6e0), N2O5 --> NO2 + NO3 #==2017/02/22; JPL 15-10; BHH,MJE==#
        c(4.00e-13), HCOOH + OH --> H2O + CO2 + HO2 #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        arr(t, T, 2.90e-12, 0.0, -345.0e0), MOH + OH --> HO2 + CH2O
        arr(t, T, 4.50e-14, 0.0, -1260.0e0), NO2 + NO3 --> NO + NO2 + O2
        c(5.80e-16), NO3 + CH2O --> HNO3 + HO2 + CO
        arr(t, T, 4.63e-12, 0.0, 350.0e0), ALD2 + OH --> 0.950MCO3 + 0.050CH2O + 0.050CO + 0.050HO2 + H2O #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 1.40e-12, 0.0, -1900.0e0), ALD2 + NO3 --> HNO3 + MCO3
        arr3(t, T, num_density, 9.70e-29, 5.6e+00, 0.0, 9.3e-12, 1.5e0, 0.0, 0.6e0), MCO3 + NO2 --> PAN #==JPL Eval 17==#
        eqconstnodim(t, T, num_density, 9.30e-29, 14000.0e0, 9.7e-29, 5.6e0, 9.3e-12, 1.5e0, 0.6e0), PAN --> MCO3 + NO2
        arr(t, T, 8.10e-12, 0.0, 270.0e0), MCO3 + NO --> MO2 + NO2 + CO2
        arr(t, T, 7.66e-12, 0.0, -1020.0e0), C2H6 + OH --> ETO2 + H2O #==2013/02/12; JPL 10-6; BHH,JMAO,EAM==#
        rRO2NOb2(t, T, num_density, 2.60e-12, 365.0e0, 2.0e0), ETO2 + NO --> ALD2 + NO2 + HO2 #==2019/05/10; Fisher2018; JAF==#
        rRO2NOa2(t, T, num_density, 2.60e-12, 365.0e0, 2.0e0), ETO2 + NO --> ETNO3 #==2019/05/10; Fisher2018; JAF==#
        arr(t, T, 2.60e-12, 0.0, 365.0e0), OTHRO2 + NO --> ALD2 + NO2 + HO2 #==2019/05/10; Fisher2018; JAF==#
        tbrnch(t, T, 7.60e-12, 0, -585.0e0, 5.87e0, 0.64e0, -816.0e0), C3H8 + OH --> B3O2
        tbrnch(t, T, 7.60e-12, 0, -585.0e0, 1.7e-1, -0.64e0, 816.0e0), C3H8 + OH --> A3O2
        rRO2NOb2(t, T, num_density, 2.90e-12, 350.0e0, 3.0e0), A3O2 + NO --> NO2 + HO2 + RCHO #==2019/05/10; Fisher2018; JAF==#
        rRO2NOa2(t, T, num_density, 2.90e-12, 350.0e0, 3.0e0), A3O2 + NO --> NPRNO3 #==2019/05/10; Fisher2018; JAF==#
        arr(t, T, 2.70e-12, 0.0, 350.0e0), PO2 + NO --> NO2 + HO2 + CH2O + ALD2
        arr(t, T, 9.10e-12, 0.0, -405.0e0), ALK4 + OH --> R4O2
        rRO2NOb2(t, T, num_density, 2.70e-12, 350.0e0, 4.5e0), R4O2 + NO --> NO2 + 0.320ACET + 0.190MEK + 0.190MO2 + 0.270HO2 + 0.320ALD2 + 0.140RCHO + 0.050A3O2 + 0.180B3O2 + 0.320OTHRO2 #==2017/02/23; ALK4 lumping fix; BHH==#
        rRO2NOa2(t, T, num_density, 2.70e-12, 350.0e0, 4.5e0), R4O2 + NO --> R4N2
        arr(t, T, 2.70e-12, 0.0, 350.0e0), R4N1 + NO --> 2.000NO2 + 0.570RCHO + 0.860ALD2 + 0.570CH2O #==2017/07/27; Fix C creation; SAS,BHH,MJE==#
        arr(t, T, 2.80e-12, 0.0, 300.0e0), ATO2 + NO --> NO2 + CH2O + MCO3 #==2017/07/27; Fix C creation; SAS,BHH,MJE==#
        arr(t, T, 2.70e-12, 0.0, 350.0e0), KO2 + NO --> 0.930NO2 + 0.930ALD2 + 0.930MCO3 + 0.070R4N2
        rRO2NOa2(t, T, num_density, 2.70e-12, 360.0e0, 3.0e0), B3O2 + NO --> NO2 + HO2 + ACET #==2019/05/10; Fisher2018; JAF==#
        rRO2NOa2(t, T, num_density, 2.70e-12, 360.0e0, 3.0e0), B3O2 + NO --> IPRNO3 #==2019/05/10; Fisher2018; JAF==#
        arr(t, T, 2.70e-12, 0.0, 350.0e0), PRN1 + NO --> 2.000NO2 + CH2O + ALD2
        arr(t, T, 2.80e-12, 0.0, -3280.0e0), ALK4 + NO3 --> HNO3 + R4O2
        c(1.60e-12), R4N2 + OH --> R4N1 + H2O
        arr(t, T, 3.15e-14, 0.0, 920.0e0), ACTA + OH --> MO2 + CO2 + H2O #==2013/02/12; JPL 10-6; BHH,JMAO,EAM==#
        arr(t, T, 6.00e-12, 0.0, 410.0e0), OH + RCHO --> RCO3 + H2O
        arr3(t, T, num_density, 9.00e-28, 8.9e0, 0.0, 7.7e-12, 0.2e0, 0.0, 0.6e0), RCO3 + NO2 --> PPN #==JPL Eval 17==#
        eqconstnodim(t, T, num_density, 9.00e-29, 14000.0e0, 9.00e-28, 8.9e0, 7.7e-12, 0.2e0, 0.6e0), PPN --> RCO3 + NO2
        arr(t, T, 6.70e-12, 0.0, 340.0e0), RCO3 + NO --> NO2 + 0.500OTHRO2 + 0.070A3O2 + 0.270B3O2 #==2019/05/10; Fisher2018; JAF==#
        c(6.50e-15), RCHO + NO3 --> HNO3 + RCO3
        regressT(t, T, 1.33e-13, 3.82e-11, -2000.0e0), ACET + OH --> ATO2 + H2O #==JPL Eval 17, p1-62-D31; EVF==#
        c(5.92e-13), A3O2 + MO2 --> HO2 + 0.750CH2O + 0.750RCHO + 0.250MOH + 0.250ROH
        c(5.92e-13), PO2 + MO2 --> HO2 + 0.500ALD2 + 1.250CH2O + 0.160HAC + 0.090RCHO + 0.250MOH + 0.250ROH
        arr(t, T, 7.40e-13, 0.0, 700.0e0), R4O2 + HO2 --> R4P
        arr(t, T, 7.40e-13, 0.0, 700.0e0), R4N1 + HO2 --> R4N2
        arr(t, T, 8.60e-13, 0.0, 700.0e0), ATO2 + HO2 --> 0.150MCO3 + 0.150OH + 0.150CH2O + 0.850ATOOH #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        rRO2HO2(t, T, 2.91e-13, 1300.0e0, 4.0e0), KO2 + HO2 --> 0.150OH + 0.150ALD2 + 0.150MCO3 + 0.850ATOOH #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        rRO2HO2(t, T, 2.91e-13, 1300.0e0, 3.0e0), B3O2 + HO2 --> RB3P #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        rRO2HO2(t, T, 2.91e-13, 1300.0e0, 3.0e0), PRN1 + HO2 --> PRPN #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        arr(t, T, 1.30e-12, 0.0, -25.0e0), MEK + OH --> KO2 + H2O
        c(3.00e-13), MO2 + ETO2 --> 0.750CH2O + 0.750ALD2 + HO2 + 0.250MOH + 0.250EOH
        c(3.00e-13), MO2 + OTHRO2 --> 0.750CH2O + 0.750ALD2 + HO2 + 0.250MOH + 0.250EOH #==2019/05/10; Fisher2018; JAF==#
        c(8.00e-16), MEK + NO3 --> HNO3 + KO2
        c(8.37e-14), R4O2 + MO2 --> 0.160ACET + 0.100MEK + 0.090MO2 + 0.140HO2 + 0.160ALD2 + 0.070RCHO + 0.030A3O2 + 0.090B3O2 + 0.160OTHRO2 + 0.250MEK + 0.750CH2O + 0.250MOH + 0.250ROH + 0.500HO2
        c(8.37e-14), R4N1 + MO2 --> NO2 + 0.200CH2O + 0.380ALD2 + 0.290RCHO + 0.150R4O2 + 0.250RCHO + 0.750CH2O + 0.250MOH + 0.250ROH + 0.500HO2
        arr(t, T, 7.50e-13, 0.0, 500.0e0), ATO2 + MO2 --> 0.300HO2 + 0.300CH2O + 0.300MCO3 + 0.200HAC + 0.200CH2O + 0.500MGLY + 0.500MOH
        c(8.37e-14), KO2 + MO2 --> 0.500ALD2 + 0.500MCO3 + 0.250MEK + 0.750CH2O + 0.250MOH + 0.250ROH + 0.500HO2
        c(8.37e-14), B3O2 + MO2 --> 0.500HO2 + 0.500ACET + 0.250ACET + 0.750CH2O + 0.250MOH + 0.250ROH + 0.500HO2
        c(8.37e-14), PRN1 + MO2 --> NO2 + 0.500CH2O + 0.500ALD2 + 0.250RCHO + 0.750CH2O + 0.250MOH + 0.250ROH + 0.500HO2
        c(3.35e-12), EOH + OH --> HO2 + ALD2 #==2013/02/12; JPL 10-6; BHH,JMAO,EAM==#
        arr(t, T, 4.60e-12, 0.0, 70.0e0), ROH + OH --> HO2 + RCHO
        c(4.10e-14), ETO2 + ETO2 --> 2.000ALD2 + 2.000HO2
        c(4.10e-14), OTHRO2 + OTHRO2 --> 2.000ALD2 + 2.000HO2 #==2019/05/10; Fisher2018; JAF==#
        c(2.70e-14), ETO2 + ETO2 --> EOH + ALD2
        c(2.70e-14), OTHRO2 + OTHRO2 --> EOH + ALD2 #==2019/05/10; Fisher2018; JAF==#
        arr(t, T, 7.40e-13, 0.0, 700.0e0), HO2 + ETO2 --> ETP
        arr(t, T, 7.40e-13, 0.0, 700.0e0), HO2 + OTHRO2 --> ETP #==2019/05/10; Fisher2018; JAF==#
        rRO2HO2(t, T, 2.91e-13, 1300.0e0, 3.0e0), A3O2 + HO2 --> RA3P #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        rRO2HO2(t, T, 2.91e-13, 1300.0e0, 3.0e0), PO2 + HO2 --> PP #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        arr(t, T, 4.30e-13, 0.0, 1040.0e0), RCO3 + HO2 --> 0.410RP + 0.150RCOOH + 0.150O3 + 0.440OH + 0.220OTHRO2 + 0.030A3O2 + 0.120B3O2 #==2019/05/10; Fisher2018; JAF==#
        arr3(t, T, num_density, 4.60e-27, 4.0e0, 0.0, 2.6e-11, 1.3e0, 0.0, 0.5e0), PRPE + OH --> PO2 #==2017/02/22; JPL 15-10; BHH,MJE==#
        arr(t, T, 5.50e-15, 0.0, -1880.0e0), PRPE + O3 --> 0.500ALD2 + 0.500CH2O + 0.120CH3CHOO + 0.100CH4 + 0.120CH2OO + 0.280MO2 + 0.560CO + 0.280HO2 + 0.360OH #==2015/09/25; Millet2015; DBM,EAM==#
        rGLYCOH_a(t, T, 8.00e-12), GLYC + OH --> 0.732CH2O + 0.361CO2 + 0.505CO + 0.227OH + 0.773HO2 + 0.134GLYX + 0.134HCOOH #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        rGLYCOH_b(t, T, 8.00e-12), GLYC + OH --> HCOOH + OH + CO #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        arr(t, T, 4.59e-13, 0.0, -1156.0e0), PRPE + NO3 --> PRN1
        arr(t, T, 3.10e-12, 0.0, 340.0e0), GLYX + OH --> HO2 + 2.000CO #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        c(1.50e-11), MGLY + OH --> MCO3 + CO
        rGLYXNO3(t, T, num_density, 1.40e-12, -1860.0e0), GLYX + NO3 --> HNO3 + HO2 + 2.000CO
        arr(t, T, 3.36e-12, 0.0, -1860.0e0), MGLY + NO3 --> HNO3 + CO + MCO3 #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        rHACOH_a(t, T, 2.15e-12, 305.0e0), HAC + OH --> MGLY + HO2 #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        rHACOH_b(t, T, 2.15e-12, 305.0e0), HAC + OH --> 0.500HCOOH + OH + 0.500ACTA + 0.500CO2 + 0.500CO + 0.500MO2 #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        arr(t, T, 1.68e-12, 0.0, 500.0e0), MCO3 + A3O2 --> MO2 + RCHO + HO2
        arr(t, T, 1.68e-12, 0.0, 500.0e0), MCO3 + PO2 --> MO2 + ALD2 + CH2O + HO2
        arr(t, T, 1.87e-13, 0.0, 500.0e0), MCO3 + A3O2 --> ACTA + RCHO
        arr(t, T, 1.87e-13, 0.0, 500.0e0), MCO3 + PO2 --> ACTA + 0.350RCHO + 0.650HAC
        arr(t, T, 1.68e-12, 0.0, 500.0e0), RCO3 + MO2 --> CH2O + HO2 + 0.500OTHRO2 + 0.070A3O2 + 0.270B3O2 #==2019/05/10; Fisher2018; JAF==#
        arr(t, T, 1.87e-13, 0.0, 500.0e0), RCO3 + MO2 --> RCOOH + CH2O
        arr(t, T, 8.78e-12, 0.0, 200.0e0), PRPN + OH --> 0.209PRN1 + 0.791OH + 0.791PROPNN #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        arr(t, T, 5.18e-12, 0.0, 200.0e0), ETP + OH --> 0.640OH + 0.360OTHRO2 + 0.640ALD2 #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        arr(t, T, 5.18e-12, 0.0, 200.0e0), RA3P + OH --> 0.640OH + 0.360A3O2 + 0.640RCHO #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        arr(t, T, 8.78e-12, 0.0, 200.0e0), RB3P + OH --> 0.791OH + 0.209B3O2 + 0.791ACET #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        arr(t, T, 8.78e-12, 0.0, 200.0e0), R4P + OH --> 0.791OH + 0.209R4O2 + 0.791RCHO #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        arr(t, T, 6.13e-13, 0.0, 200.0e0), RP + OH --> RCO3 #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        arr(t, T, 8.78e-12, 0.0, 200.0e0), PP + OH --> 0.791OH + 0.209PO2 + 0.791HAC #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        arr(t, T, 4.82e-11, 0.0, -400.0e0), LVOC + OH --> OH #==2017/06/14; Marais2016; EAM==#
        arr(t, T, 6.13e-13, 0.0, 200.0e0), OH + MAP --> MCO3 #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        c(1.40e-18), C2H6 + NO3 --> ETO2 + HNO3 #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        arr(t, T, 2.50e-12, 0.0, 500.0e0), MCO3 + MCO3 --> 2.000MO2
        arr(t, T, 1.80e-12, 0.0, 500.0e0), MCO3 + MO2 --> CH2O + MO2 + HO2
        arr(t, T, 2.00e-13, 0.0, 500.0e0), MCO3 + MO2 --> ACTA + CH2O
        arr(t, T, 1.68e-12, 0.0, 500.0e0), R4O2 + MCO3 --> MO2 + 0.320ACET + 0.190MEK + 0.270HO2 + 0.320ALD2 + 0.130RCHO + 0.050A3O2 + 0.180B3O2 + 0.320OTHRO2 #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        arr(t, T, 1.68e-12, 0.0, 500.0e0), ATO2 + MCO3 --> MO2 + MCO3 + CH2O #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        arr(t, T, 1.68e-12, 0.0, 500.0e0), KO2 + MCO3 --> MO2 + ALD2 + MCO3
        arr(t, T, 1.68e-12, 0.0, 500.0e0), B3O2 + MCO3 --> MO2 + HO2 + ACET
        arr(t, T, 1.68e-12, 0.0, 500.0e0), R4N1 + MCO3 --> MO2 + NO2 + 0.390CH2O + 0.750ALD2 + 0.570RCHO + 0.300R4O2
        arr(t, T, 1.68e-12, 0.0, 500.0e0), PRN1 + MCO3 --> MO2 + NO2 + CH2O + ALD2
        arr(t, T, 1.87e-13, 0.0, 500.0e0), R4O2 + MCO3 --> MEK + ACTA
        arr(t, T, 1.87e-13, 0.0, 500.0e0), ATO2 + MCO3 --> MGLY + ACTA #==2017/07/27; Fix C creation; SAS,BHH,MJE==#
        arr(t, T, 1.87e-13, 0.0, 500.0e0), KO2 + MCO3 --> MEK + ACTA
        arr(t, T, 1.87e-13, 0.0, 500.0e0), R4N1 + MCO3 --> RCHO + ACTA + NO2
        arr(t, T, 1.87e-13, 0.0, 500.0e0), PRN1 + MCO3 --> RCHO + ACTA + NO2
        arr(t, T, 1.87e-13, 0.0, 500.0e0), B3O2 + MCO3 --> ACET + ACTA
        arr(t, T, 1.68e-12, 0.0, 500.0e0), MCO3 + ETO2 --> MO2 + ALD2 + HO2
        arr(t, T, 1.68e-12, 0.0, 500.0e0), MCO3 + OTHRO2 --> MO2 + ALD2 + HO2 #==2019/05/10; Fisher2018; JAF==#
        arr(t, T, 1.87e-13, 0.0, 500.0e0), MCO3 + ETO2 --> ACTA + ALD2
        arr(t, T, 1.87e-13, 0.0, 500.0e0), MCO3 + OTHRO2 --> ACTA + ALD2 #==2019/05/10; Fisher2018; JAF==#
        arr(t, T, 2.50e-12, 0.0, 500.0e0), RCO3 + MCO3 --> MO2 + 0.500OTHRO2 + 0.070A3O2 + 0.270B3O2 #==2019/05/10; Fisher2018; JAF==#
        arr(t, T, 8.50e-13, 0.0, -2450.0e0), NO3 + NO3 --> 2.000NO2 + O2
        arr3(t, T, num_density, 1.00e-30, 4.8e+00, 0.0, 7.2e-12, 2.1e0, 0.0, 0.6e0), MO2 + NO2 --> MPN #==2012/02/12; Browne2011; ECB==#
        arr3nodim(t, T, num_density, 1.05e-02, 4.8e+00, -11234.0e0, 7.58e16, 2.1e0, -11234.0e0, 0.6e0), MPN --> MO2 + NO2 #==2012/02/12; Browne2011; ECB==#
        arr(t, T, 1.20e-11, 0.0, -280.0e0), DMS + OH --> SO2 + MO2 + CH2O
        rDMSOH(t, T, num_density, 8.20e-39, 5376.0e0, 1.05e-5, 3644.0e0), DMS + OH --> 0.750SO2 + 0.250MSA + MO2
        arr(t, T, 1.90e-13, 0.0, 530.0e0), DMS + NO3 --> SO2 + HNO3 + MO2 + CH2O
        arr3(t, T, num_density, 3.30e-31, 4.3e+00, 0.0, 1.6e-12, 0.0, 0.0, 0.6e0), SO2 + OH --> SO4 + HO2 #==+M==#
        arr(t, T, 1.60e-11, 0.0, -780.0e0), Br + O3 --> BrO + O2 #==2012/06/07; Parrella2012; JPP==#
        arr(t, T, 4.50e-12, 0.0, 460.0e0), BrO + HO2 --> HOBr + O2 #==2012/06/07; Parrella2012; JPP==#
        arr(t, T, 4.80e-12, 0.0, -310.0e0), Br + HO2 --> HBr + O2 #==2012/06/07; Parrella2012; JPP==#
        arr(t, T, 5.50e-12, 0.0, 200.0e0), HBr + OH --> Br + H2O #==2012/06/07; Parrella2012; JPP==#
        arr(t, T, 2.40e-12, 0.0, 40.0e0), BrO + BrO --> 2.000Br + O2 #==2012/06/07; Parrella2012; JPP==#
        arr(t, T, 2.80e-14, 0.0, 860.0e0), BrO + BrO --> Br2 + O2 #==2012/06/07; Parrella2012; JPP==#
        arr(t, T, 8.80e-12, 0.0, 260.0e0), BrO + NO --> Br + NO2 #==2012/06/07; Parrella2012; JPP==#
        c(4.90e-11), Br + BrNO3 --> Br2 + NO3 #==2012/06/07; Parrella2012; JPP==#
        arr(t, T, 2.10e-11, 0.0, 240.0e0), Br2 + OH --> HOBr + Br #==2012/06/07; Parrella2012; JPP==#
        arr(t, T, 1.20e-10, 0.0, -430.0e0), HOBr + O --> OH + BrO #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 5.80e-12, 0.0, -1500.0e0), HBr + O --> OH + Br #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 1.70e-11, 0.0, 250.0e0), BrO + OH --> Br + HO2 #==2012/06/07; Parrella2012; JPP==#
        c(1.60e-11), Br + NO3 --> BrO + NO2 #==2012/06/07; Parrella2012; JPP==#
        arr(t, T, 1.70e-11, 0.0, -800.0e0), Br + CH2O --> HBr + HO2 + CO #==2012/06/07; Parrella2012; JPP==#
        arr(t, T, 1.80e-11, 0.0, -460.0e0), Br + ALD2 --> HBr + MCO3 #==2017/07/27; Parrella2012,Fix C creation; SAS,BHH,MJE==#
        arr(t, T, 1.66e-10, 0.0, -7000.0e0), Br + ACET --> HBr + ATO2 #==2017/07/27; Parrella2012,Fix C creation; SAS,BHH,MJE==#
        arr(t, T, 2.36e-10, 0.0, -6411.0e0), Br + C2H6 --> HBr + ETO2 #==2017/07/27; Parrella2012,Fix C creation; SAS,BHH,MJE==#
        arr(t, T, 8.77e-11, 0.0, -4330.0e0), Br + C3H8 --> HBr + A3O2 #==2017/07/27; Parrella2012,Fix C creation; SAS,BHH,MJE==#
        arr3(t, T, num_density, 4.20e-31, 2.4e0, 0.0, 2.7e-11, 0.0, 0.0, 0.6e0), Br + NO2 --> BrNO2 #==2012/06/07; Parrella2012; JPP==#
        arr3(t, T, num_density, 5.40e-31, 3.1e0, 0.0, 6.5e-12, 2.9e0, 0.0, 0.6e0), BrO + NO2 --> BrNO3 #==2017/02/22; JPL 15-10; BHH,MJE==#
        arr(t, T, 9.00e-13, 0.0, -360.0e0), CHBr3 + OH --> 3.000Br #==2017/02/22; JPL 15-10; BHH,MJE==#
        arr(t, T, 2.00e-12, 0.0, -840.0e0), CH2Br2 + OH --> 2.000Br #==2012/06/07; Parrella2012; JPP==#
        arr(t, T, 1.42e-12, 0.0, -1150.0e0), CH3Br + OH --> Br + H2O + HO2 #==2017/03/08; JPL 15-10; TS,BHH,MJE==#
        arr(t, T, 1.63e-10, 0.0, 60.0e0), O1D + H2O --> 2.000OH #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 2.15e-11, 0.0, 110.0e0), O1D + N2 --> O + N2 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 3.30e-11, 0.0, 55.0e0), O1D + O2 --> O + O2 #==2014/02/03; Eastham2014; SDE==#
        c(1.20e-10), O1D + H2 --> H + OH #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 4.63e-11, 0.0, 20.0e0), O1D + N2O --> N2 + O2 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 7.25e-11, 0.0, 20.0e0), O1D + N2O --> 2.000NO #==2014/02/03; Eastham2014; SDE==#
        c(1.31e-10), O1D + CH4 --> MO2 + OH #==2014/02/03; Eastham2014; SDE==#
        c(0.09e-10), O1D + CH4 --> CH2O + H2 #==2014/02/03; Eastham2014; SDE==#
        c(0.35e-10), O1D + CH4 --> CH2O + H + HO2 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 6.00e-34, 2.4e0, 0) * num_density, O + O2 --> O3 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 8.00e-12, 0.0, -2060.0e0), O + O3 --> 2.000O2 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 2.80e-12, 0.0, -1800.0e0), OH + H2 --> H2O + H #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 1.80e-11, 0.0, 180.0e0), O + OH --> O2 + H #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 3.00e-11, 0.0, 200.0e0), HO2 + O --> OH + O2 #==2014/02/03; Eastham2014; SDE==#
        c(1.20e-10), O1D + O3 --> 2.000O2 #==2014/02/03; Eastham2014; SDE==#
        c(1.20e-10), O1D + O3 --> 2.000O + O2 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 2.10e-11, 0.0, -2200.0e0), OCS + O --> CO + SO2 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 1.10e-13, 0.0, -1200.0e0), OCS + OH --> CO2 + SO2 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 5.10e-12, 0.0, 210.0e0), NO2 + O --> NO + O2 #==2014/02/03; Eastham2014; SDE==#
        c(1.00e-11), NO3 + O --> NO2 + O2 #==2014/02/03; Eastham2014; SDE==#
        arr3(t, T, num_density, 9.00e-32, 1.5e+00, 0.0, 3.0e-11, 0.0, 0.0, 0.6e0), NO + O --> NO2 #==2014/02/03; Eastham2014; SDE==#
        arr3(t, T, num_density, 2.50e-31, 1.8e+00, 0.0, 2.2e-11, 0.7e0, 0.0, 0.6e0), NO2 + O --> NO3 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 1.40e-12, 0.0, -2000.0e0), H2O2 + O --> OH + HO2 #==2014/02/03; Eastham2014; SDE==#
        arr3(t, T, num_density, 4.40e-32, 1.3e+00, 0.0, 7.5e-11, -0.2e0, 0.0, 0.6e0), H + O2 --> HO2 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 1.40e-10, 0.0, -470.0e0), H + O3 --> OH + O2 #==2014/02/03; Eastham2014; SDE==#
        c(7.20e-11), H + HO2 --> 2.000OH #==2014/02/03; Eastham2014; SDE==#
        c(1.60e-12), H + HO2 --> O + H2O #==2014/02/03; Eastham2014; SDE==#
        c(6.90e-12), H + HO2 --> H2 + O2 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 1.50e-11, 0.0, -3600.0e0), N + O2 --> NO + O #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 2.10e-11, 0.0, 100.0e0), N + NO --> N2 + O #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 5.80e-12, 0.0, 220.0e0), N + NO2 --> N2O + O #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 1.90e-11, 0.0, 230.0e0), BrO + O --> Br + O2 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 3.40e-11, 0.0, -1600.0e0), CH2O + O --> CO + HO2 + OH #==2014/02/03; Eastham2014; SDE==#
        c(1.50e-10), O1D + HCl --> 0.090O + 0.090HCl + 0.240H + 0.670Cl + 0.240ClO + 0.670OH #==2014/02/03; Eastham2014; SDE==#
        c(1.50e-10), O1D + HBr --> 0.200O + 0.200HBr + 0.150BrO + 0.650OH + 0.150H + 0.650Br #==2014/02/03; Eastham2014; SDE==#
        c(2.70e-10), O1D + Cl2 --> 0.250O + 0.250Cl2 + 0.750Cl + 0.750ClO #==2014/02/03; Eastham2014; SDE==#
        c(3.30e-10), O1D + CCl4 --> 0.140O + 0.140CCl4 + 0.860ClO + 2.580Cl #==2014/02/03; Eastham2014; SDE==#
        c(1.80e-10), O1D + CH3Br --> 0.440BrO + MO2 + 0.560Br #==2014/02/03; Eastham2014; SDE==#
        c(2.70e-10), O1D + CH2Br2 --> 0.050O + 0.050CH2Br2 + 0.950BrO + 0.950Br #==2014/02/03; Eastham2014; SDE==#
        c(6.60e-10), O1D + CHBr3 --> 0.320O + 0.320CHBr3 + 0.680BrO + 1.360Br #==2014/02/03; Eastham2014; SDE==#
        c(1.02e-10), O1D + HCFC22 --> 0.280O + 0.280HCFC22 + 0.550ClO + 0.170Cl #==2017/02/22; JPL 15-10; BHH,MJE==#
        c(2.30e-10), O1D + CFC11 --> 0.120O + 0.120CFC11 + 0.880ClO + 1.760Cl #==2014/02/03; Eastham2014; SDE==#
        c(1.40e-10), O1D + CFC12 --> 0.140O + 0.140CFC12 + 0.860ClO + 0.860Cl #==2014/02/03; Eastham2014; SDE==#
        c(1.50e-10), O1D + H1211 --> 0.360O + 0.360H1211 + 0.310BrO + 0.310Cl + 0.330Br + 0.330ClO #==2014/02/03; Eastham2014; SDE==#
        c(1.00e-10), O1D + H1301 --> 0.590O + 0.590H1301 + 0.410BrO #==2014/02/03; Eastham2014; SDE==#
        c(2.60e-10), O1D + HCFC141b --> 0.310O + 0.310HCFC141b + 0.690ClO + 0.690Cl #==2014/02/03; Eastham2014; SDE==#
        c(2.00e-10), O1D + HCFC142b --> 0.260O + 0.260HCFC142b + 0.740ClO #==2017/02/22; JPL 15-10; BHH,MJE==#
        c(2.00e-10), O1D + HCFC123 --> 0.210O + 0.210HCFC123 + 0.790Cl + 0.790ClO #==2014/02/03; Eastham2014; SDE==#
        c(2.32e-10), O1D + CFC113 --> 0.250O + 0.250CFC113 + 1.500Cl + 0.750ClO #==2017/02/22; JPL 15-10; BHH,MJE==#
        arr(t, T, 1.30e-10, 0.0, -25.0e0), O1D + CFC114 --> 0.250O + 0.250CFC114 + 0.750Cl + 0.750ClO #==2017/02/22; JPL 15-10; BHH,MJE==#
        arr(t, T, 5.40e-11, 0.0, -30.0e0), O1D + CFC115 --> 0.700O + 0.700CFC115 + 0.300ClO #==2017/02/22; JPL 15-10; BHH,MJE==#
        arr(t, T, 1.60e-10, 0.0, 0.0e0), O1D + H2402 --> 0.250O + 0.250H2402 + 0.750Br + 0.750BrO #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 2.60e-12, 0.0, -1100.0e0), OH + Cl2 --> HOCl + Cl #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 1.80e-11, 0.0, -600.0e0), MO2 + ClO --> ClOO + HO2 + CH2O #==2017/03/20; JPL 15-10; TS,BHH,MJE==#
        arr(t, T, 7.40e-12, 0.0, 270.0e0), OH + ClO --> HO2 + Cl #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 6.00e-13, 0.0, 230.0e0), OH + ClO --> HCl + O2 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 1.40e-12, 0.0, 600.0e0), OH + OClO --> HOCl + O2 #==2017/02/22; JPL 15-10; BHH,MJE==#
        arr(t, T, 6.00e-13, 0.0, 670.0e0), OH + Cl2O2 --> HOCl + ClOO #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 1.80e-12, 0.0, -250.0e0), OH + HCl --> H2O + Cl #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 3.00e-12, 0.0, -500.0e0), OH + HOCl --> H2O + ClO #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 2.40e-12, 0.0, -1250.0e0), OH + ClNO2 --> HOCl + NO2 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 1.20e-12, 0.0, -330.0e0), OH + ClNO3 --> HOCl + NO3 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 1.96e-12, 0.0, -1200.0e0), OH + CH3Cl --> Cl + HO2 + H2O #==2017/02/22; JPL 15-10; BHH,MJE==#
        arr(t, T, 2.61e-12, 0.0, -944.0e0), OH + CH2Cl2 --> 2.000Cl + HO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 4.69e-12, 0.0, -1134.0e0), OH + CHCl3 --> 3.000Cl + HO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 1.64e-12, 0.0, -1520.0e0), OH + CH3CCl3 --> 3.000Cl + H2O #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 9.20e-13, 0.0, -1560.0e0), OH + HCFC22 --> Cl + H2O #==2017/02/22; JPL 15-10; BHH,MJE==#
        arr(t, T, 1.25e-12, 0.0, -1600.0e0), OH + HCFC141b --> 2.000Cl + H2O #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 1.30e-12, 0.0, -1770.0e0), OH + HCFC142b --> Cl + H2O #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 7.40e-13, 0.0, -900.0e0), OH + HCFC123 --> 2.000Cl + H2O #==2017/02/22; JPL 15-10; BHH,MJE==#
        arr(t, T, 7.10e-12, 0.0, -1270.0e0), CH4 + Cl --> HCl + MO2 #==2017/03/08; JPL 15-10; TS,BHH,MJE==#
        arr(t, T, 7.32e-11, 0.0, -30.0e0), CH2O + Cl --> CO + HCl + HO2 #==2017/09/22; Sherwen2016b; TS,JAS,SDE==#
        arr(t, T, 2.30e-11, 0.0, -200.0e0), Cl + O3 --> ClO + O2 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 3.05e-11, 0.0, -2270.0e0), Cl + H2 --> H + HCl #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 1.10e-11, 0.0, -980.0e0), Cl + H2O2 --> HO2 + HCl #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 1.40e-11, 0.0, 270.0e0), Cl + HO2 --> O2 + HCl #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 3.60e-11, 0.0, -375.0e0), Cl + HO2 --> OH + ClO #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 2.80e-11, 0.0, 85.0e0), ClO + O --> Cl + O2 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 2.60e-12, 0.0, 290.0e0), ClO + HO2 --> O2 + HOCl #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 6.40e-12, 0.0, 290.0e0), ClO + NO --> Cl + NO2 #==2014/02/03; Eastham2014; SDE==#
        arr3(t, T, num_density, 1.80e-31, 3.4e+00, 0.0, 1.50e-11, 1.9e0, 0.0, 0.6e0), ClO + NO2 --> ClNO3 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 1.00e-12, 0.0, -1590.0e0), ClO + ClO --> Cl2 + O2 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 3.00e-11, 0.0, -2450.0e0), ClO + ClO --> Cl + ClOO #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 3.50e-13, 0.0, -1370.0e0), ClO + ClO --> OClO + Cl #==2014/02/03; Eastham2014; SDE==#
        arr3(t, T, num_density, 2.20e-33, 3.1e+00, 0.0, 1.8e-10, 0.0, 0.0, 0.6e0), Cl + O2 --> ClOO #==2014/02/03; Eastham2014; SDE==#
        eqconstnodim(t, T, num_density, 6.60e-25, 2502.0e0, 2.20e-33, 3.1e+00, 1.8e-10, 0.0e0, 0.6e0), ClOO --> Cl + O2 #==JPL 15-10; XW==#
        arr3(t, T, num_density, 1.90e-32, 3.6e+00, 0.0, 3.7e-12, 1.6e0, 0.0, 0.6e0), ClO + ClO --> Cl2O2 #==2017/02/22; JPL 15-10; BHH,MJE==#
        eqconstnodim(t, T, num_density, 2.16e-27, 8537.0e0, 1.90e-32, 3.6e+00, 3.7e-12, 1.6e0, 0.6e0), Cl2O2 --> 2.000ClO #==JPL 15-10; XW==#
        c(2.30e-10), ClOO + Cl --> Cl2 + O2 #==2014/02/03; Eastham2014; SDE==#
        c(1.20e-11), ClOO + Cl --> 2.000ClO #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 9.50e-13, 0.0, 550.0e0), ClO + BrO --> Br + OClO #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 2.30e-12, 0.0, 260.0e0), ClO + BrO --> Br + ClOO #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 4.10e-13, 0.0, 290.0e0), ClO + BrO --> BrCl + O2 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 3.60e-12, 0.0, -840.0e0), ClNO3 + O --> ClO + NO3 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 6.50e-12, 0.0, 135.0e0), ClNO3 + Cl --> Cl2 + NO3 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 2.17e-11, 0.0, -1130.0e0), CH3Cl + Cl --> CO + 2.000HCl + HO2 #==2014/02/03; Eastham2014; SDE==#
        arr(t, T, 1.24e-12, 0.0, -1070.0e0), CH2Cl2 + Cl --> CO + HCl + 2.000Cl + HO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 3.77e-12, 0.0, -1011.0e0), CHCl3 + Cl --> CO + HCl + 3.000Cl + HO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        c(2.00e-13), Cl + HCOOH --> HCl + CO2 + H2O #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        c(1.60e-10), Cl + MO2 --> ClO + CH2O + HO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        c(5.7e-11), Cl + MP --> HCl + MO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 7.2e-11, 0.0, -70.0e0), Cl + C2H6 --> HCl + ETO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        c(7.4e-11), Cl + ETO2 --> ClO + HO2 + ALD2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        c(7.4e-11), Cl + OTHRO2 --> ClO + HO2 + ALD2 #==2019/05/10; Fisher2018; JAF==#
        c(5.5e-11), Cl + MOH --> HCl + CH2O + HO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        c(9.6e-11), Cl + EOH --> HCl + ALD2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        c(2.8e-14), Cl + ACTA --> HCl + MO2 + CO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 6.54e-11, 0.0, 60.0e0), Cl + C3H8 --> HCl + B3O2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 8.12e-11, 0.0, -90.0e0), Cl + C3H8 --> HCl + A3O2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 7.70e-11, 0.0, -1000.0e0), Cl + ACET --> HCl + ATO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 7.60e-11, 0.0, 500.0e0), Cl + ISOP --> HCl + 0.5IHOO1 + 0.5IHOO4 #==2019/11/06; Sherwen2016b;KHB,TS,JAS,SDE==#
        c(2.05e-10), Cl + ALK4 --> HCl + R4O2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr3(t, T, num_density, 4.00e-28, 0.0, 0.0, 2.8e-10, 0.0, 0.0, 0.6e0), Cl + PRPE --> HCl + PO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        c(3.60e-12), Br + PRPE --> HBr + PO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr3(t, T, num_density, 1.80e-32, 1.0e0, 0.0, 1.77e-11, 0.0, 0.0, 0.6e0), I + NO --> INO #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 8.40e-11, 0.0, -2620.0e0), INO + INO --> I2 + 2.000NO #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr3(t, T, num_density, 3.00e-31, 1.0e0, 0.0, 6.6e-11, 0.0, 0.0, 0.63e0), I + NO2 --> IONO #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arrnodim(t, T, 9.94e+17, 0.0, -11859.0e0), IONO --> I + NO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 2.90e-11, 0.0, -2600.0e0), IONO + IONO --> I2 + 2.000NO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        c(1.50e-12), I2 + NO3 --> I + IONO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr3(t, T, num_density, 7.50e-31, 3.5e0, 0.0, 7.6e-12, 1.5e0, 0.0, 0.6e0), IO + NO2 --> IONO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arrnodim(t, T, 2.10e+15, 0.0, -13670.0e0), IONO2 --> IO + NO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 9.10e-11, 0.0, -146.0e0), IONO2 + I --> I2 + NO3 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        c(1.20e-11), I + BrO --> IO + Br #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 3.00e-12, 0.0, 510.0e0), IO + BrO --> Br + I + O2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 1.20e-11, 0.0, 510.0e0), IO + BrO --> Br + OIO #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        c(1.00e-10), IO + OIO --> I2O3 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        c(1.50e-10), OIO + OIO --> I2O4 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        c(3.80e-02), I2O4 --> 2.000OIO #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 1.10e-12, 0.0, 542.0e0), OIO + NO --> IO + NO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 5.10e-12, 0.0, 280.0e0), IO + ClO --> I + OClO #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 2.81e-12, 0.0, 280.0e0), IO + ClO --> I + Cl + O2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 1.02e-12, 0.0, 280.0e0), IO + ClO --> ICl + O2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 2.30e-11, 0.0, -870.0e0), I + O3 --> IO + O2 #==2017/09/22; Sherwen2017;TS,JAS,SDE==#
        arr(t, T, 1.50e-11, 0.0, -1090.0e0), I + HO2 --> HI + O2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        c(1.80e-10), I2 + OH --> HOI + I #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        c(3.00e-11), HI + OH --> I + H2O #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        c(5.00e-12), HOI + OH --> IO + H2O #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 1.30e-11, 0.0, 570.0e0), IO + HO2 --> HOI + O2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 9.10e-12, 0.0, 240.0e0), IO + NO --> I + NO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 6.00e-12, 0.0, 500.0e0), IO + IO --> I + OIO #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 9.00e-12, 0.0, 500.0e0), IO + IO --> I2O2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arrnodim(t, T, 1.00e+12, 0.0, -9770.0e0), I2O2 --> 2.000IO #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arrnodim(t, T, 2.50e+14, 0.0, -9770.0e0), I2O2 --> OIO + I #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        arr(t, T, 2.90e-12, 0.0, -1100.0e0), CH3I + OH --> H2O + I + MO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        c(2.40e-12), ETHLN + OH --> CH2O + CO2 + NO2 #==2017/06/15, Marais2016, EAM==#
        c(6.70e-13), PROPNN + OH --> NO2 + MGLY #==2017/07/14; MCMv3.3; KRT,JAF,CCM,EAM,KHB,RHS==#
        c(1.20e-15), CH2OO + CO --> CH2O #==2015/09/25; Millet2015; DBM,EAM==#
        c(1.00e-14), CH2OO + NO --> CH2O + NO2 #==2015/09/25; Millet2015; DBM,EAM==#
        c(1.00e-15), CH2OO + NO2 --> CH2O + NO3 #==2015/09/25; Millet2015; DBM,EAM==#
        c(1.70e-15), CH2OO + H2O --> 0.730HMHP + 0.210HCOOH + 0.060CH2O + 0.060H2O2 #==2019/11/06; Bates2019; KHB==#
        arrsq(t, T, 2.88e-35, 0.0, 1391.0e0), CH2OO + H2O + H2O --> 0.400HMHP + 0.540HCOOH + 0.060CH2O + 0.060H2O2 #==2019/11/06; Bates2019; KHB==#
        c(1.40e-12), CH2OO + O3 --> CH2O #==2019/11/06; Bates2019; KHB==#
        c(3.70e-11), CH2OO + SO2 --> CH2O + SO4 #==2019/11/06; Bates2019; KHB==#
        c(1.20e-15), CH3CHOO + CO --> ALD2 #==2015/09/25; Millet2015; DBM,EAM==#
        c(1.00e-14), CH3CHOO + NO --> ALD2 + NO2 #==2015/09/25; Millet2015; DBM,EAM==#
        c(1.00e-15), CH3CHOO + NO2 --> ALD2 + NO3 #==2015/09/25; Millet2015; DBM,EAM==#
        c(7.00e-14), CH3CHOO + SO2 --> ALD2 + SO4 #==2015/09/25; Millet2015; DBM,EAM==#
        c(6.00e-18), CH3CHOO + H2O --> ALD2 + H2O2 #==2015/09/25; Millet2015; DBM,EAM==#
        c(1.00e-17), CH3CHOO + H2O --> ACTA #==2015/09/25; Millet2015; DBM,EAM==#
        arr(t, T, 1.21e-11, 0.0, 440.0e0), MTPA + OH --> PIO2 #==2017/03/23; IUPAC2010; EVF==#
        arr(t, T, 1.21e-11, 0.0, 440.0e0), MTPO + OH --> PIO2 #==2017/03/23; IUPAC2010; EVF==#
        c(4.00e-12), PIO2 + NO --> 0.820HO2 + 0.820NO2 + 0.230CH2O + 0.430RCHO + 0.110ACET + 0.440MEK + 0.070HCOOH + 0.120MONITS + 0.060MONITU #==2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS==#
        c(1.50e-11), PIO2 + HO2 --> PIP #==2017/03/23; Roberts1992; EVF==#
        arr(t, T, 3.56e-14, 0.0, 708.0e0), PIO2 + MO2 --> HO2 + 0.750CH2O + 0.250MOH + 0.250ROH + 0.750RCHO + 0.750MEK #==2017/07/14; Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 7.40e-13, 0.0, 765.0e0), PIO2 + MCO3 --> 0.500HO2 + 0.500MO2 + RCHO + MEK + RCOOH #==2017/03/23; Roberts1992; EVF==#
        c(1.20e-12), PIO2 + NO3 --> HO2 + NO2 + RCHO + MEK #==2017/03/23; Roberts1992; EVF==#
        arr(t, T, 5.00e-16, 0.0, -530.0e0), MTPA + O3 --> 0.850OH + 0.100HO2 + 0.620KO2 + 0.140CO + 0.020H2O2 + 0.650RCHO + 0.530MEK #==2017/07/14; Atkinson2003; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 5.00e-16, 0.0, -530.0e0), MTPO + O3 --> 0.850OH + 0.100HO2 + 0.620KO2 + 0.140CO + 0.020H2O2 + 0.650RCHO + 0.530MEK #==2017/07/14; Atkinson2003; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 8.33e-13, 0.0, 490.0e0), MTPA + NO3 --> 0.100OLNN + 0.900OLND #==2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 8.33e-13, 0.0, 490.0e0), MTPO + NO3 --> 0.100OLNN + 0.900OLND #==2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 4.20e-11, 0.0, 401.0e0), LIMO + OH --> LIMO2 #==2017/07/14; Gill2002; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 2.95e-15, 0.0, -783.0e0), LIMO + O3 --> 0.850OH + 0.100HO2 + 0.160OTHRO2 + 0.420KO2 + 0.020H2O2 + 0.140CO + 0.460PRPE + 0.040CH2O + 0.790MACR + 0.010HCOOH + 0.070RCOOH #==2017/07/14; Atkinson2003; KRT,JAF,CCM,EAM,KHB,RHS==#
        c(1.22e-11), LIMO + NO3 --> 0.500OLNN + 0.500OLND #==2017/07/14; Fry2014,Atkinson2003; KRT,JAF,CCM,EAM,KHB,RHS==#
        c(4.00e-12), LIMO2 + NO --> 0.686HO2 + 0.780NO2 + 0.220MONITU + 0.289PRPE + 0.231CH2O + 0.491RCHO + 0.058HAC + 0.289MEK #==2017/07/14; Browne2014,Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS==#
        c(1.50e-11), LIMO2 + HO2 --> PIP #==2017/07/14; Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 3.56e-14, 0.0, 708.0e0), LIMO2 + MO2 --> HO2 + 0.192PRPE + 1.040CH2O + 0.308MACR + 0.250MOH + 0.250ROH #==2017/07/14; Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 7.40e-13, 0.0, 765.0e0), LIMO2 + MCO3 --> 0.500HO2 + 0.500MO2 + 0.192PRPE + 0.385CH2O + 0.308MACR + 0.500RCOOH #==2017/07/14; Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS==#
        c(1.20e-12), LIMO2 + NO3 --> HO2 + NO2 + 0.385PRPE + 0.385CH2O + 0.615MACR #==2017/07/14; Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 3.40e-12, 0.0, 190.0e0), PIP + OH --> 0.490OH + 0.440R4O2 + 0.080RCHO + 0.410MEK #==2017/07/14; Goliff2013; KRT,JAF,CCM,EAM,KHB,RHS==#
        c(4.00e-12), OLNN + NO --> HO2 + NO2 + MONITS #==2017/07/14; Browne2014,Goliff2013; KRT,JAF,CCM,EAM,KHB,RHS==#
        c(4.00e-12), OLND + NO --> 2.000NO2 + 0.287CH2O + 1.240RCHO + 0.464MEK #==2017/07/14; Goliff2013; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 1.66e-13, 0.0, 1300.0e0), OLNN + HO2 --> 0.700MONITS + 0.300MONITU #==2017/07/14; Browne2014,Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 1.66e-13, 0.0, 1300.0e0), OLND + HO2 --> 0.700MONITS + 0.300MONITU #==2017/07/14; Browne2014,Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 1.60e-13, 0.0, 708.0e0), OLNN + MO2 --> 2.000HO2 + CH2O + 0.700MONITS + 0.300MONITU #==2017/07/14; Browne2014,Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 9.68e-14, 0.0, 708.0e0), OLND + MO2 --> 0.500HO2 + 0.500NO2 + 0.965CH2O + 0.930RCHO + 0.348MEK + 0.250MOH + 0.250ROH + 0.350MONITS + 0.150MONITU #==2017/07/14; Browne2014,Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 8.85e-13, 0.0, 765.0e0), OLNN + MCO3 --> HO2 + MO2 + 0.700MONITS + 0.300MONITU #==2017/07/14; Browne2014,Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 5.37e-13, 0.0, 765.0e0), OLND + MCO3 --> 0.500MO2 + NO2 + 0.287CH2O + 1.240RCHO + 0.464MEK + 0.500RCOOH #==2017/07/14; Browne2014,Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS==#
        c(1.20e-12), OLNN + NO3 --> HO2 + NO2 + 0.700MONITS + 0.300MONITU #==2017/07/14; Browne2014,Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS==#
        c(1.20e-12), OLND + NO3 --> 2.000NO2 + 0.287CH2O + 1.240RCHO + 0.464MEK #==2017/07/14; Browne2014,Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 7.00e-14, 0.0, 1000.0e0), OLNN + OLNN --> HO2 + 1.400MONITS + 0.600MONITU #==2017/07/14; Browne2014,Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 4.25e-14, 0.0, 1000.0e0), OLNN + OLND --> 0.500HO2 + 0.500NO2 + 0.202CH2O + 0.640RCHO + 0.149MEK + 1.050MONITS + 0.450MONITU #==2017/07/14; Browne2014,Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 2.96e-14, 0.0, 1000.0e0), OLND + OLND --> NO2 + 0.504CH2O + 1.210RCHO + 0.285MEK + 0.700MONITS + 0.300MONITU #==2017/07/14; Browne2014,Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS==#
        c(4.80e-12), MONITS + OH --> HONIT #==2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS==#
        c(7.29e-11), MONITU + OH --> HONIT #==2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS==#
        c(1.67e-16), MONITU + O3 --> HONIT #==2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 3.15e-13, 0.0, -448.0e0), MONITU + NO3 --> HONIT #==2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 3.15e-13, 0.0, -448.0e0), MONITS + NO3 --> HONIT #==2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS==#
        c(2.78e-04), IONITA --> INDIOL + HNO3 #==2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS==#
        c(2.78e-04), MONITA --> INDIOL + HNO3 #==2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS==#
        rOHHNO3(t, T, num_density, 2.41e-14, 460.0e0, 2.69e-17, 2199.0e0, 6.51e-34, 1335.0e0), HONIT + OH --> NO3 + HAC #==2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS==#
        arr(t, T, 8.00e-13, 0.0, -1000.0e0), MENO3 + OH --> CH2O + NO2 #==2019/05/16; JPL 15-10,Fisher2018; JAF==#
        arr(t, T, 1.00e-12, 0.0, -490.0e0), ETNO3 + OH --> ALD2 + NO2 #==2019/05/16; JPL 15-10,Fisher2018; JAF==#
        arr(t, T, 1.20e-12, 0.0, -320.0e0), IPRNO3 + OH --> ACET + NO2 #==2019/05/16; JPL 15-10,Fisher2018; JAF==#
        c(7.10e-13), NPRNO3 + OH --> RCHO + NO2 #==2019/05/16; JPL 15-10,Fisher2018; JAF==#
        c(1.3e-17), ISOP + O3 --> 0.416MACR + 0.177MVK + 0.28OH + 0.407CO2 + 0.407CO + 0.407MO2 + 0.16HO2 + 0.58CH2OO + 0.827CH2O + 0.013H2O2 #==2019/11/06; Bates2019; KHB==#
        rISO1(t, T, 1.7e-11, 3.90e2, 9.33e-2, 5.05e15, -1.22e4, 1.79e14, -8.830e3), ISOP + OH --> LISOPOH + IHOO1 #==2019/11/06; Bates2019; KHB==#
        rISO1(t, T, 1.0e-11, 3.90e2, 2.26e-1, 2.22e9, -7.160e3, 1.75e14, -9.054e3), ISOP + OH --> LISOPOH + IHOO4 #==2019/11/06; Bates2019; KHB==#
        rISO2(t, T, 1.7e-11, 3.90e2, 9.33e-2, 5.05e15, -1.22e4, 1.79e14, -8.830e3), ISOP + OH --> 0.3MCO3 + 0.3MGLY + 0.3CH2O + 0.15HPALD3 + 0.25HPALD1 + 0.4HO2 + 0.6CO + 1.5OH + 0.3HPETHNL + LISOPOH #==2019/11/06; Bates2019; KHB==#
        rISO2(t, T, 1.0e-11, 3.90e2, 2.26e-1, 2.22e9, -7.160e3, 1.75e14, -9.054e3), ISOP + OH --> 0.3CH2O + 0.15HPALD4 + 0.25HPALD2 + 1.5OH + 0.9CO + 0.7HO2 + 0.3MGLY + 0.3ATOOH + LISOPOH #==2019/11/06; Bates2019; KHB==#
        arrpls(t, T, 2.12e-13, -1300e0, 0, 1.1644e0, -7.0485e-4), IHOO1 + HO2 --> 0.063MVK + 0.063OH + 0.063HO2 + 0.063CH2O + 0.937RIPA #==2019/11/06; Bates2019; KHB==#
        arrpls(t, T, 2.12e-13, -1300e0, 0, -0.1644e0, 7.0485e-4), IHOO1 + HO2 --> RIPC #==2019/11/06; Bates2019; KHB==#
        arrpls(t, T, 2.12e-13, -1300e0, 0, 1.2038e0, -9.0435e-4), IHOO4 + HO2 --> 0.063MACR + 0.063OH + 0.063HO2 + 0.063CH2O + 0.937RIPB #==2019/11/06; Bates2019; KHB==#
        arrpls(t, T, 2.12e-13, -1300e0, 0, -0.2038e0, 9.0435e-4), IHOO4 + HO2 --> RIPD #==2019/11/06; Bates2019; KHB==#
        arrplsnodim(t, T, 1.04e11, 9.746e3, 0, 1.1644e0, -7.0485e-4), IHOO1 --> CH2O + OH + MVK #==2019/11/06; Bates2019; KHB==#
        tunplsnodim(t, T, 5.05e15, -1.22e4, 1.0e8, -0.0128e0, 5.1242e-5), IHOO1 --> 0.15HPALD3 + 0.25HPALD1 + 0.4HO2 + 0.6CO + 1.5OH + 0.3CH2O + 0.3MGLY + 0.3HPETHNL + 0.3MCO3 #==2019/11/06; Bates2019; KHB==#
        arrplsnodim(t, T, 1.88e11, 9.752e3, 0, 1.2038e0, -9.0435e-4), IHOO4 --> MACR + OH + CH2O #==2019/11/06; Bates2019; KHB==#
        tunplsnodim(t, T, 2.22e9, -7.160e3, 1.0e8, -0.0306e0, 1.1346e-4), IHOO4 --> 0.15HPALD4 + 0.25HPALD2 + 1.5OH + 0.3CH2O + 0.9CO + 0.7HO2 + 0.3MGLY + 0.3ATOOH #==2019/11/06; Bates2019; KHB==#
        arrpls(t, T, 6.92e-14, 0, 0, 1.1644e0, -7.0485e-4), IHOO1 + IHOO1 --> 2MVK + 2HO2 + 2CH2O #==2019/11/06; Bates2019; KHB==#
        arrpls(t, T, 5.74e-12, 0, 0, 1.2038e0, -9.0435e-4), IHOO4 + IHOO4 --> 2MACR + 2HO2 + 2CH2O #==2019/11/06; Bates2019; KHB==#
        arrpls(t, T, 1.54e-12, 0, 0, 2.3682e0, -1.6092e-3), IHOO1 + IHOO4 --> MACR + MVK + 2HO2 + 2CH2O #==2019/11/06; Bates2019; KHB==#
        arrpls(t, T, 2.49e-12, 0, 0, -0.1644e0, 7.0485e-4), IHOO1 + IHOO1 --> HO2 + HC5A + CO + OH + MVKHP #==2019/11/06; Bates2019; KHB==#
        arrpls(t, T, 3.94e-12, 0, 0, -0.2038e0, 9.0435e-4), IHOO4 + IHOO4 --> HO2 + HC5A + CO + OH + MCRHP #==2019/11/06; Bates2019; KHB==#
        arrpls(t, T, 1.54e-12, 0, 0, -0.3682e0, 1.6092e-3), IHOO1 + IHOO4 --> HO2 + HC5A + CO + OH + 0.5MVKHP + 0.5MCRHP #==2019/11/06; Bates2019; KHB==#
        arrpls(t, T, 2.0e-12, 0, 0, 1.1644e0, -7.0485e-4), IHOO1 + MO2 --> MVK + 2HO2 + 2CH2O #==2019/11/06; Bates2019; KHB==#
        arrpls(t, T, 2.0e-12, 0, 0, -0.1644e0, 7.0485e-4), IHOO1 + MO2 --> CH2O + 0.5HC5A + 1.5HO2 + 0.5MVKHP + 0.5CO + 0.5OH #==2019/11/06; Bates2019; KHB==#
        arrpls(t, T, 2.0e-12, 0, 0, 1.2038e0, -9.0435e-4), IHOO4 + MO2 --> MACR + 2HO2 + 2CH2O #==2019/11/06; Bates2019; KHB==#
        arrpls(t, T, 2.0e-12, 0, 0, -0.2038e0, 9.0435e-4), IHOO4 + MO2 --> CH2O + 0.5HC5A + 1.5HO2 + 0.5MCRHP + 0.5CO + 0.5OH #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 3.50e2, 1.19e0, 6.0e0, 1.1644e0, 7.05e-4), IHOO1 + NO --> IHN2 #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 3.50e2, 1.19e0, 6.0e0, 1.1644e0, 7.05e-4), IHOO1 + NO --> NO2 + MVK + HO2 + CH2O #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 3.50e2, 1.421e0, 6.0e0, -0.1644e0, -7.05e-4), IHOO1 + NO --> IHN4 #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 3.50e2, 1.421e0, 6.0e0, -0.1644e0, -7.05e-4), IHOO1 + NO --> NO2 + 0.45HC5A + 0.45HO2 + 0.55MVKHP + 0.55CO + 0.55OH #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 3.50e2, 1.297e0, 6.0e0, 1.2038e0, 9.04e-4), IHOO4 + NO --> IHN3 #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 3.50e2, 1.297e0, 6.0e0, 1.2038e0, 9.04e-4), IHOO4 + NO --> NO2 + MACR + HO2 + CH2O #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 3.50e2, 1.421e0, 6.0e0, -0.2038e0, -9.04e-4), IHOO4 + NO --> IHN1 #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 3.50e2, 1.421e0, 6.0e0, -0.2038e0, -9.04e-4), IHOO4 + NO --> NO2 + 0.45HO2 + 0.45HC5A + 0.55MCRHP + 0.55CO + 0.55OH #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 1.17e-11, 0.0, 450.0e0), HPALD1 + OH --> 0.035MVK + 0.315HPALD1OO + 0.15IDC + 0.33MVKHP + 0.085HO2 + 0.085CH2O + 0.085MGLY + 0.085ICHE + 1.085OH + 0.45CO #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 1.17e-11, 0.0, 450.0e0), HPALD2 + OH --> 0.035MACR + 0.315HPALD2OO + 0.15IDC + 0.17MCRHP + 0.165HO2 + 0.165CH2O + 0.165MGLY + 0.165ICHE + 1.165OH + 0.37CO #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.20e-11, 0.0, 390.0e0), HPALD3 + OH --> OH + 0.230MVK + 0.420CO + 0.190MVKHP + 0.580ICHE #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 3.50e-11, 0.0, 390.0e0), HPALD4 + OH --> OH + 0.770ICHE + 0.230CO + 0.090MCRHP + 0.140MACR #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 4.64e-12, 0.0, 650.0e0), HC5A + OH --> 1.065OH + 0.355CO2 + 0.638CO + 0.355MGLY + 0.283HO2 + 0.294IEPOXAOO + 0.125MVKHP + 0.158MCRHP + 0.068IEPOXBOO #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 9.85e-12, 0.0, 410.0e0), ICHE + OH --> OH + 1.5CO + 0.5CH2O + 0.5MGLY + 0.5HAC #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 3.00e-12, 0.0, 650.0e0), IDC + OH --> CO + HO2 + MVKPC #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.47e-12, 0.0, 390.0e0), RIPA + OH --> 0.655IHPOO3 + 0.345IHPOO1 + 0.005LVOC #==2019/11/06; Bates2019; KHB==#
        rEPO(t, T, num_density, 1.62e-11, 3.90e2, 4.77e-21), RIPA + OH --> 0.67IEPOXA + 0.33IEPOXB + OH + 0.005LVOC #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 4.35e-12, 0.0, 390.0e0), RIPB + OH --> 0.655IHPOO3 + 0.345IHPOO2 + 0.005LVOC #==2019/11/06; Bates2019; KHB==#
        rEPO(t, T, num_density, 2.85e-11, 390.0e0, 4.77e-21), RIPB + OH --> 0.68IEPOXA + 0.32IEPOXB + OH + 0.005LVOC #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 6.10e-12, 0.0, 200.0e0), RIPA + OH --> 0.75IHOO1 + 0.125MVK + 0.25CO + 0.125MVKHP + 0.25HO2 + 0.005LVOC #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 4.10e-12, 0.0, 200.0e0), RIPB + OH --> 0.51IHOO4 + 0.16ICHOO + 0.33CO + 0.33HO2 + 0.165MACR + 0.165MCRHP + 0.005LVOC #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 3.53e-11, 0.0, 390.0e0), RIPC + OH --> 0.595IHPOO1 + 0.03IHOO1 + 0.06HC5A + 0.024HO2 + 0.009HPALD3 + 0.015HPALD1 + 0.405OH + 0.036CO + 0.018CH2O + 0.018MGLY + 0.018HPETHNL + 0.018MCO3 + 0.255IEPOXD + 0.005LVOC #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 3.53e-11, 0.0, 390.0e0), RIPD + OH --> 0.255IHPOO2 + 0.03IHOO4 + 0.745OH + 0.06HC5A + 0.009HPALD4 + 0.015HPALD2 + 0.042HO2 + 0.018CH2O + 0.054CO + 0.018MGLY + 0.018ATOOH + 0.595IEPOXD + 0.005LVOC #==2019/11/06; Bates2019; KHB==#
        arrnodim(t, T, 1.59e+13, 0.0, -10000.0e0), IHPOO1 --> 0.176ICPDH + 0.824IDHPE + OH #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 3.50e2, 2.1e0, 9.0e0, 1.0e0, 0.0e0), IHPOO1 + NO --> 0.716MCRHP + 0.716CH2O + 0.284HPETHNL + 0.284HAC + NO2 + HO2 #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 3.50e2, 2.1e0, 9.0e0, 1.0e0, 0.0e0), IHPOO1 + NO --> ITHN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.47e-13, 0.0, 1300.0e0), IHPOO1 + HO2 --> 0.725IDHDP + 0.14MCRHP + 0.14CH2O + 0.135HPETHNL + 0.135HAC + 0.275OH + 0.275HO2 #==2019/11/06; Bates2019; KHB==#
        arrnodim(t, T, 2.91e+13, 0.0, -10000.0e0), IHPOO2 --> 0.548ICPDH + 0.452IDHPE + OH #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 3.50e2, 2.315e0, 9.0e0, 1.0e0, 0.0e0), IHPOO2 + NO --> 0.706MVKHP + 0.706CH2O + 0.294GLYC + 0.294ATOOH + NO2 + HO2 #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 3.50e2, 2.315e0, 9.0e0, 1.0e0, 0.0e0), IHPOO2 + NO --> ITHN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.47e-13, 0.0, 1300.0e0), IHPOO2 + HO2 --> 0.725IDHDP + 0.14MVKHP + 0.14CH2O + 0.135GLYC + 0.135ATOOH + 0.275OH + 0.275HO2 #==2019/11/06; Bates2019; KHB==#
        arrnodim(t, T, 1.875e+13, 0.0, -10000.0e0), IHPOO3 --> IDHPE #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 3.50e2, 3.079e0, 9.0e0, 1.0e0, 0.0e0), IHPOO3 + NO --> GLYC + HAC + NO2 + OH #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 3.50e2, 3.079e0, 9.0e0, 1.0e0, 0.0e0), IHPOO3 + NO --> ITHN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.47e-13, 0.0, 1300.0e0), IHPOO3 + HO2 --> 0.35IDHDP + 0.65GLYC + 0.65HAC + 1.3OH #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 3.22e-11, 0.0, -400.0e0), IEPOXD + OH --> 0.75ICHE + 0.75HO2 + 0.25ICHOO #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 1.05e-11, 0.0, -400.0e0), IEPOXA + OH --> ICHE + HO2 #==2019/11/06; Bates2019; KHB==#
        rEPO(t, T, num_density, 5.82e-11, -4.00e2, 1.14e-20), IEPOXA + OH --> 0.67IEPOXAOO + 0.33IEPOXBOO #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 8.25e-12, 0.0, -400.0e0), IEPOXB + OH --> ICHE + HO2 #==2019/11/06; Bates2019; KHB==#
        rEPO(t, T, num_density, 3.75e-11, -4.00e2, 8.91e-21), IEPOXB + OH --> 0.81IEPOXAOO + 0.19IEPOXBOO #==2019/11/06; Bates2019; KHB==#
        arrnodim(t, T, 1.875e+13, 0.0, -10000.0e0), IEPOXAOO --> IDCHP + HO2 #==2019/11/06; Bates2019; KHB==#
        arrnodim(t, T, 1.0e+7, 0.0, -5000.0e0), IEPOXAOO --> OH + CO + MVKDH #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.38e-13, 0.0, 1300.0e0), IEPOXAOO + HO2 --> 0.13CO + 0.65OH + 0.65HO2 + 0.13MVKDH + 0.52GLYC + 0.52MGLY + 0.35ICPDH #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 3.50e2, 13.098e0, 8.0e0, 1.0e0, 0.0e0), IEPOXAOO + NO --> 0.2MVKDH + HO2 + NO2 + 0.2CO + 0.8GLYC + 0.8MGLY #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 3.50e2, 13.098e0, 8.0e0, 1.0e0, 0.0e0), IEPOXAOO + NO --> ITCN #==2019/11/06; Bates2019; KHB==#
        arrnodim(t, T, 1.875e+13, 0.0, -10000.0e0), IEPOXBOO --> IDCHP + HO2 #==2019/11/06; Bates2019; KHB==#
        arrnodim(t, T, 1.0e+7, 0.0, -5000.0e0), IEPOXBOO --> CO + OH + MCRDH #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 3.50e2, 16.463e0, 8.0e0, 1.0e0, 0.0e0), IEPOXBOO + NO --> NO2 + HO2 + 0.8GLYX + 0.8HAC + 0.2CO + 0.2MCRDH #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 3.50e2, 16.463e0, 8.0e0, 1.0e0, 0.0e0), IEPOXBOO + NO --> ITCN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.38e-13, 0.0, 1300.0e0), IEPOXBOO + HO2 --> 0.13CO + 0.65OH + 0.65HO2 + 0.13MCRDH + 0.52HAC + 0.52GLYX + 0.35ICPDH #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.38e-13, 0.0, 1300.0e0), ICHOO + HO2 --> 0.35ICPDH + 0.65OH + 0.52CO + 0.13MVKHC + 0.65CH2O + 0.65HO2 + 0.52HAC #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 3.50e2, 13.098e0, 8.0e0, 1.0e0, 0.0e0), ICHOO + NO --> ITCN #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 3.50e2, 13.098e0, 8.0e0, 1.0e0, 0.0e0), ICHOO + NO --> NO2 + 0.8HAC + 0.8CO + CH2O + HO2 + 0.2MVKHC #==2019/11/06; Bates2019; KHB==#
        arrnodim(t, T, 1.875e+13, 0.0, -10000.0e0), ICHOO --> HO2 + 2.000CO + HAC + OH #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.70e-12, 0.0, 350.0e0), HPALD1OO + NO --> NO2 + OH + CO2 + MVK #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.38e-13, 0.0, 1300.0e0), HPALD1OO + HO2 --> OH + OH + CO2 + MVK #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.70e-12, 0.0, 350.0e0), HPALD2OO + NO --> NO2 + OH + CO2 + MACR #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.38e-13, 0.0, 1300.0e0), HPALD2OO + HO2 --> OH + OH + CO2 + MACR #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 7.14e-12, 0.0, 390.0e0), IHN2 + OH --> ISOPNOO1 #==2019/11/06; Bates2019; KHB==#
        rEPO(t, T, num_density, 6.30e-12, 390.0e0, 1.62e-19), IHN2 + OH --> 0.67IEPOXA + 0.33IEPOXB + NO2 #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 1.02e-11, 0.0, 390.0e0), IHN3 + OH --> ISOPNOO2 #==2019/11/06; Bates2019; KHB==#
        rEPO(t, T, num_density, 1.05e-11, 390.0e0, 2.49e-19), IHN3 + OH --> 0.67IEPOXA + 0.33IEPOXB + NO2 #==2019/11/06; Bates2019; KHB==#
        rEPO(t, T, num_density, 1.55e-11, 390.0e0, 2.715e-19), IHN1 + OH --> IEPOXD + NO2 #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.04e-11, 0.0, 390.0e0), IHN1 + OH --> IDHNDOO1 #==2019/11/06; Bates2019; KHB==#
        rEPO(t, T, num_density, 9.52e-12, 390.0e0, 2.715e-19), IHN4 + OH --> IEPOXD + NO2 #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.95e-11, 0.0, 390.0e0), IHN4 + OH --> IDHNDOO2 #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 7.5e-12, 0.0, 20.0e0), IHN1 + OH --> 0.6OH + 0.6CO + 0.6MCRHNB + 0.4HO2 + 0.4ICN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 7.5e-12, 0.0, 20.0e0), IHN4 + OH --> 0.6OH + 0.6CO + 0.6MVKN + 0.4HO2 + 0.4ICN #==2019/11/06; Bates2019; KHB==#
        arrnodim(t, T, 1.875e+13, 0.0, -10000.0e0), ISOPNOO1 --> ITCN + HO2 #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.60e-13, 0.0, 1300.0e0), ISOPNOO1 + HO2 --> 0.482ITHN + 0.059MCRHN + 0.059CH2O + 0.459GLYC + 0.459HAC + 0.059HO2 + 0.459NO2 + 0.518OH #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 350.0e0, 6.32e0, 11.0e0, 1.0e0, 0.0e0), ISOPNOO1 + NO --> 0.272MCRHN + 0.272CH2O + 0.728GLYC + 0.728HAC + 0.272HO2 + 1.728NO2 #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 350.0e0, 6.32e0, 11.0e0, 1.0e0, 0.0e0), ISOPNOO1 + NO --> IDN #==2019/11/06; Bates2019; KHB==#
        arrnodim(t, T, 1.875e+13, 0.0, -10000.0e0), ISOPNOO2 --> ITCN + HO2 #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.60e-13, 0.0, 1300.0e0), ISOPNOO2 + HO2 --> 0.401ITHN + 0.599MVKN + 0.599CH2O + 0.599HO2 + 0.599OH #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 350.0e0, 7.941e0, 11.0e0, 1.0e0, 0.0e0), ISOPNOO2 + NO --> MVKN + CH2O + HO2 + NO2 #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 350.0e0, 7.941e0, 11.0e0, 1.0e0, 0.0e0), ISOPNOO2 + NO --> IDN #==2019/11/06; Bates2019; KHB==#
        arrnodim(t, T, 1.256e+13, 0.0, -10000.0e0), IDHNDOO1 --> ITCN + HO2 #==2019/11/06; Bates2019; KHB==#
        arrnodim(t, T, 5.092e+12, 0.0, -10000.0e0), IDHNDOO2 --> ITCN + HO2 #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.60e-13, 0.0, 1300.0e0), IDHNDOO1 + HO2 --> 0.418ITHN + 0.551PROPNN + 0.551GLYC + 0.031MCRHNB + 0.031CH2O + 0.582HO2 + 0.582OH #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 350.0e0, 4.712e0, 11.0e0, 1.0e0, 0.0e0), IDHNDOO1 + NO --> 0.935PROPNN + 0.935GLYC + 0.065MCRHNB + 0.065CH2O + HO2 + NO2 #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 350.0e0, 4.712e0, 11.0e0, 1.0e0, 0.0e0), IDHNDOO1 + NO --> IDN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.60e-13, 0.0, 1300.0e0), IDHNDOO2 + HO2 --> 0.494ITHN + 0.441HAC + 0.441ETHLN + 0.065MVKN + 0.065CH2O + 0.506OH + 0.506HO2 #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 350.0e0, 2.258e0, 11.0e0, 1.0e0, 0.0e0), IDHNDOO2 + NO --> 0.858HAC + 0.858ETHLN + 0.142MVKN + 0.142CH2O + HO2 + NO2 #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 350.0e0, 2.258e0, 11.0e0, 1.0e0, 0.0e0), IDHNDOO2 + NO --> IDN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.60e-13, 0.0, 1300.0e0), IDHNBOO + HO2 --> 0.379HO2 + 0.379OH + 0.621ITHN + 0.094MCRHNB + 0.242GLYC + 0.242PROPNN + 0.010MVKN + 0.033HAC + 0.033ETHLN + 0.104CH2O #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 350.0e0, 1.851e0, 11.0e0, 1.0e0, 0.0e0), IDHNBOO + NO --> 0.355MCRHNB + 0.546PROPNN + 0.546GLYC + 0.028MVKN + 0.071ETHLN + 0.071HAC + HO2 + NO2 + 0.383CH2O #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 350.0e0, 1.851e0, 11.0e0, 1.0e0, 0.0e0), IDHNBOO + NO --> IDN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.95e-12, 0.0, 450.0e0), ISOP + NO3 --> 0.465INO2B + 0.535INO2D + LISOPNO3 #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.47e-13, 0.0, 1300.0e0), INO2B + HO2 --> 0.473INPB + 0.048MACR + 0.479MVK + 0.527OH + 0.527CH2O + 0.527NO2 #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.47e-13, 0.0, 1300.0e0), INO2D + HO2 --> INPD #==2019/11/06; Bates2019; KHB==#
        c(1.61e-12), INO2B + INO2B --> 1.737MVK + 0.123MACR + 1.860CH2O + 1.860NO2 + 0.070INPB + 0.070ICN #==2019/11/06; Bates2019; KHB==#
        c(2.56e-12), INO2B + INO2D --> 0.399INPB + 0.544MVK + 0.532ICN + 0.563NO2 + 0.474INA + 0.089HO2 + 0.019MACR + 0.563CH2O + 0.032IHN1 #==2019/11/06; Bates2019; KHB==#
        c(3.71e-12), INO2D + INO2D --> 0.064HO2 + 0.340INA + 0.861ICN + 0.671IHN1 + 0.127IHN4 #==2019/11/06; Bates2019; KHB==#
        c(1.18e-12), INO2D + MO2 --> 0.298IHN1 + 0.057IHN4 + 0.244INA + 0.401ICN + 0.355MOH + 0.336HO2 + 0.645CH2O #==2019/11/06; Bates2019; KHB==#
        c(2.80e-13), INO2B + MO2 --> 0.355INPB + 0.583MVK + 0.028MACR + 0.034ICN + 0.611HO2 + 1.577CH2O + 0.611NO2 + 0.034MOH #==2019/11/06; Bates2019; KHB==#
        c(1.92e-12), INO2B + MCO3 --> CH2O + NO2 + MO2 + 0.903MVK + 0.097MACR #==2019/11/06; Bates2019; KHB==#
        c(7.71e-12), INO2D + MCO3 --> MO2 + 0.841INA + 0.159HO2 + 0.159ICN #==2019/11/06; Bates2019; KHB==#
        c(2.3e-12), INO2B + NO3 --> CH2O + 2NO2 + 0.903MVK + 0.097MACR #==2019/11/06; Bates2019; KHB==#
        c(2.3e-12), INO2D + NO3 --> NO2 + 0.841INA + 0.159HO2 + 0.159ICN #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 350.0e0, 12.915e0, 9.0e0, 1.0e0, 0.0e0), INO2B + NO --> 2NO2 + CH2O + 0.096MACR + 0.904MVK #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 350.0e0, 12.915e0, 9.0e0, 1.0e0, 0.0e0), INO2B + NO --> IDN #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 350.0e0, 1.412e0, 9.0e0, 1.0e0, 0.0e0), INO2D + NO --> NO2 + 0.159HO2 + 0.159ICN + 0.841INA #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 350.0e0, 1.412e0, 9.0e0, 1.0e0, 0.0e0), INO2D + NO --> IDN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.50e-14, 0.0, -300.0e0), INA + O2 --> ICN + HO2 #==2019/11/06; Bates2019; KHB==#
        arrnodim(t, T, 1.00e+20, 0.0, -10000.0e0), INA --> IDHNBOO #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 5.88e-12, 0.0, 390.0e0), INPB + OH --> 0.670IHPNBOO + 0.33IDHNBOO #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 1.61e-11, 0.0, 390.0e0), INPD + OH --> IHPNDOO #==2019/11/06; Bates2019; KHB==#
        rEPO(t, T, num_density, 4.471e-12, 390.0e0, 2.28e-20), INPB + OH --> OH + ITHN #==2019/11/06; Bates2019; KHB==#
        rEPO(t, T, num_density, 8.77e-12, 390.0e0, 2.185e-20), INPD + OH --> OH + ITHN #==2019/11/06; Bates2019; KHB==#
        rEPO(t, T, num_density, 1.493e-11, 390.0e0, 2.715e-19), INPD + OH --> NO2 + ICHE #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.278e-12, 0.0, 200.0e0), INPB + OH --> INO2B #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 3.40e-12, 0.0, 200.0e0), INPD + OH --> INO2D #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 7.50e-12, 0.0, 20.0e0), INPD + OH --> ICN + OH #==2019/11/06; Bates2019; KHB==#
        arrnodim(t, T, 6.55e+12, 0.0, -10000.0e0), IHPNDOO --> OH + ITCN #==2019/11/06; Bates2019; KHB==#
        arrnodim(t, T, 8.72e+12, 0.0, -10000.0e0), IHPNBOO --> OH + 0.5ITCN + 0.5ITHN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.64e-13, 0.0, 1300.0e0), IHPNBOO + HO2 --> 0.234ITHN + 0.060MCRHNB + 0.340GLYC + 0.249HPETHNL + 0.004MCRHP + 0.008MVKN + 0.009ATOOH + 0.054MVKHP + 0.042HAC + 1.147OH + 0.326HO2 + 0.058NO2 + 0.126CH2O + 0.589PROPNN + 0.051ETHLN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.64e-13, 0.0, 1300.0e0), IHPNDOO + HO2 --> 0.387ITHN + 0.073MCRHNB + 0.471HPETHNL + 0.015MVKN + 0.054ATOOH + 0.646OH + 0.580HO2 + 0.088CH2O + 0.471PROPNN + 0.054ETHLN #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 350.0e0, 6.092e0, 12.0e0, 1.0e0, 0.0e0), IHPNBOO + NO --> 0.384GLYC + 0.170MCRHNB + 0.303HPETHNL + 0.014MVKN + 0.051HAC + 0.013ATOOH + 0.059MVKHP + 0.006MCRHP + 0.687PROPNN + 0.064ETHLN + 0.249CH2O + 1.065NO2 + 0.500HO2 + 0.435OH #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 350.0e0, 6.092e0, 12.0e0, 1.0e0, 0.0e0), IHPNBOO + NO --> IDN #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 350.0e0, 4.383e0, 12.0e0, 1.0e0, 0.0e0), IHPNDOO + NO --> 0.291MCRHNB + 0.590HPETHNL + 0.070ATOOH + 0.049MVKN + 0.590PROPNN + 0.070ETHLN + 0.340CH2O + 1.000NO2 + 0.904HO2 + 0.096OH #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 350.0e0, 4.383e0, 12.0e0, 1.0e0, 0.0e0), IHPNDOO + NO --> IDN #==2019/11/06; Bates2019; KHB==#
        rEPO(t, T, num_density, 2.97e-12, 390.0e0, 2.715e-19), ICN + OH --> NO2 + ICHE #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 9.35e-12, 0.0, 390.0e0), ICN + OH --> 0.244OH + 0.539CO + 0.295HO2 + 0.378MCRHNB + 0.461ICNOO + 0.161MVKN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.70e-12, 0.0, 350.0e0), ICNOO + NO --> 0.67ICNOO + 0.33CO2 + 0.33CO + 0.33HO2 + 0.231PROPNN + NO2 + 0.099ETHLN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.54e-13, 0.0, 1300.0e0), ICNOO + HO2 --> 0.67ICNOO + 0.33CO2 + 0.33CO + 0.33HO2 + 0.231PROPNN + OH + 0.099ETHLN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 1.00e-11, 0.0, 0.0e0), IDN + OH --> 0.565NO2 + 0.565ITHN + 0.435IDNOO #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.70e-12, 0.0, 350.0e0), IDNOO + NO --> PROPNN + 1.11NO2 + 0.11GLYC + 0.89ETHLN + 0.89HO2 #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.71e-13, 0.0, 1300.0e0), IDNOO + HO2 --> 0.18IDN + 0.09NO2 + 0.09GLYC + 0.82OH + 0.73HO2 + 0.82PROPNN + 0.73ETHLN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.60e-12, 0.0, 610.0e0), MVK + OH --> MVKOHOO #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 8.50e-16, 0.0, -1520.0e0), MVK + O3 --> 0.545MGLY + 0.500CH2OO + 0.600CH2O + 0.380MCO3 + 0.100HO2 + 0.080OH + 0.180CO + 0.075PYAC + 0.045H2O2 #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 4.40e-12, 0.0, 380.0e0), MACR + OH --> 0.036ATOOH + 0.036CO + 0.036HO2 + 0.964MCROHOO #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.70e-12, 0.0, 470.0e0), MACR + OH --> MACR1OO #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 1.40e-15, 0.0, -2100.0e0), MACR + O3 --> 0.880MGLY + 0.880CH2OO + 0.120CH2O + 0.120OH + 0.120CO + 0.120MCO3 #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 1.80e-13, 0.0, -1190.0e0), MACR + NO3 --> 0.320HNO3 + 0.320MACR1OO + 0.680OH + 0.680CO + 0.680PROPNN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 1.24e-12, 0.0, 380.0e0), MVKN + OH --> 0.241CH2O + 0.690NO3 + 0.020OH + 0.449MGLY + 0.449HCOOH + 0.241PYAC + 0.290MVKHCB + 0.310NO2 + 0.040MCO3 #==2019/11/06; Bates2019; KHB==#
        c(5.77e-11), MVKHP + OH --> 0.53MVKHC + 0.47MVKHCB + OH #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.70e-12, 0.0, 470.0e0), MCRHP + OH --> 0.77CO + OH + 0.77HAC + 0.23ATOOH + 0.23CO2 #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 1.39e-11, 0.0, 380.0e0), MCRHN + OH --> MACRNO2 #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.70e-12, 0.0, 470.0e0), MCRHNB + OH --> 0.250CO + OH + PROPNN + 0.750CO2 #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.70e-12, 0.0, 350.0e0), C4HVP1 + NO --> NO2 + MVKOHOO #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 1.93e-13, 0.0, 1300.0e0), C4HVP1 + HO2 --> OH + MVKOHOO #==2019/11/06; Bates2019; KHB==#
        c(9.00e-12), C4HVP1 + NO2 --> MVKN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.70e-12, 0.0, 350.0e0), C4HVP2 + NO --> NO2 + MCROHOO #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 1.93e-13, 0.0, 1300.0e0), C4HVP2 + HO2 --> OH + MCROHOO #==2019/11/06; Bates2019; KHB==#
        c(9.00e-12), C4HVP2 + NO2 --> MCRHN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 3.71e-12, 0.0, 983.0e0), MCRENOL + OH --> 0.75CO + 0.285OH + 0.715HO2 + 0.653PYAC + 0.097CO2 + 0.097MCO3 + 0.063MVKHCB + 0.187MGLY + 0.187HCOOH #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 5.00e-12, 0.0, 470.0e0), MVKPC + OH --> OH + CO + MGLY #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 8.70e-12, 0.0, 70.0e0), MVKDH + OH --> 0.4MVKHCB + 0.6MVKHC + HO2 #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 5.00e-12, 0.0, 470.0e0), MVKHCB + OH --> OH + MGLY #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.00e-12, 0.0, 70.0e0), MVKHC + OH --> 2CO + HO2 + MCO3 #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.4e-11, 0.0, 70.0e0), MCRDH + OH --> 0.16MVKHCB + HO2 + 0.84HAC + 0.84CO #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.12e-13, 0.0, 1300.0e0), MVKOHOO + HO2 --> 0.360MCO3 + 0.360GLYC + 0.665OH + 0.305HO2 + 0.255MVKHC + 0.335MVKHP + 0.050MGLY + 0.050CH2O #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 350.0e0, 4.573e0, 6.0e0, 1.0e0, 0.0e0), MVKOHOO + NO --> 0.758MCO3 + 0.758GLYC + 0.242MGLY + 0.242CH2O + 0.242HO2 + NO2 #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 350.0e0, 4.573e0, 6.0e0, 1.0e0, 0.0e0), MVKOHOO + NO --> 0.438MVKN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.12e-13, 0.0, 1300.0e0), MCROHOO + HO2 --> 0.41MCRHP + 0.507HAC + 0.507CO + 0.507HO2 + 0.59OH + 0.59O2 + 0.083MGLY + 0.083CH2O #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 3.14e-12, 0.0, 580.0e0), MACR1OO + HO2 --> 0.5MACR1OOH + 0.5CH2O + 0.325CO + 0.325MO2 + 0.175MCO3 + 0.5CO2 + 0.5OH + 0.13O3 #==2019/11/06; Bates2019; KHB==#
        c(1.66e-11), MACR1OOH + OH --> 0.165MACR1OO + 0.585OH + 0.488HAC + 0.488CO + 0.098HMML + 0.415CO2 + 0.25CH2O + 0.087MCO3 + 0.162MO2 #==2019/11/06; Bates2019; KHB==#
        arrnodim(t, T, 2.90e+7, 0.0, -5297.0e0), MCROHOO --> HAC + CO + OH #==2019/11/06; Bates2019; KHB==#
        rALK(t, T, num_density, 2.7e-12, 350.0e0, 2.985e0, 6.0e0, 1.0e0, 0.0e0), MCROHOO + NO --> 0.86HAC + 0.86CO + 0.86HO2 + NO2 + 0.14MGLY + 0.14CH2O #==2019/11/06; Bates2019; KHB==#
        rNIT(t, T, num_density, 2.7e-12, 350.0e0, 2.985e0, 6.0e0, 1.0e0, 0.0e0), MCROHOO + NO --> MCRHN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 8.7e-12, 0.0, 290.0e0), MACR1OO + NO --> 0.35MCO3 + 0.65MO2 + 0.65CO + CH2O + CO2 + NO2 #==2019/11/06; Bates2019; KHB==#
        rPAN_acac(t, T, num_density, 2.591e-28, -6.87e0, 1.125e-11, -1.105e0, 0.3e0), MACR1OO + NO2 --> MPAN #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 3.14e-12, 0.0, 580.0e0), MACRNO2 + HO2 --> 0.5HAC + 0.5OH + 0.5CO2 + 0.5NO2 + 0.13O3 + 0.37MCRHN + 0.13MCRHNB #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 7.50e-12, 0.0, 290.0e0), MACRNO2 + NO --> HAC + 2NO2 + CO2 #==2019/11/06; Bates2019; KHB==#
        rPAN_acac(t, T, num_density, 2.591e-28, -6.87e0, 1.125e-11, -1.105e0, 0.3e0), MACRNO2 + NO2 --> MPAN + NO2 #==2019/11/06; Bates2019; KHB==#
        c(4.00e-12), MACRNO2 + NO3 --> HAC + 2NO2 + CO2 #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 2.9e-12, 0.0, 500.0e0), MACRNO2 + MO2 --> 0.7HAC + 0.7CO2 + 0.7NO2 + 0.7HO2 + CH2O + 0.3MCRHNB #==2019/11/06; Bates2019; KHB==#
        arrnodim(t, T, 1.58e+16, 0.0, -13500.0e0), MPAN --> MACR1OO + NO2 #==2019/11/06; Bates2019; KHB==#
        c(2.90e-11), MPAN + OH --> 0.75HMML + NO3 + 0.25HAC + 0.25CO #==2019/11/06; Bates2019; KHB==#
        c(4.33e-12), HMML + OH --> 0.700MGLY + 0.700OH + 0.300MCO3 + 0.300HCOOH #==2019/11/06; Bates2019; KHB==#
        c(1.00e-11), ICPDH + OH --> CO + 0.5HO2 + 0.5OH + 0.5MCRHP + 0.35MVKDH + 0.15MCRDH #==2019/11/06; Bates2019; KHB==#
        c(2.25e-11), IDCHP + OH --> 0.888CO + 0.444OH + 0.444HO2 + 0.318MVKHC + 0.08IEPOXAOO + 0.126MVKHCB + 0.444MVKPC + 0.032IEPOXBOO #==2019/11/06; Bates2019; KHB==#
        c(3.00e-12), IDHDP + OH --> OH + 0.333ICPDH + 0.667IDHPE #==2019/11/06; Bates2019; KHB==#
        c(3.00e-12), IDHPE + OH --> OH + CO2 + 0.571MCRHP + 0.429MVKHP #==2019/11/06; Bates2019; KHB==#
        c(1.00e-11), ITCN + OH --> CO + NO2 + 0.75MVKHP + 0.25MCRHP #==2019/11/06; Bates2019; KHB==#
        c(3.00e-12), ITHN + OH --> 0.300OH + 0.620HO2 + 0.920ITCN + 0.037IDHNBOO + 0.041ICNOO + 0.022MCRENOL + 0.022NO2 + 0.022CH2O #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 1.40e-12, 0.0, -1860.0e0), ETHLN + NO3 --> HNO3 + NO2 + MCO3 #==2019/11/06; Bates2019; KHB==#
        c(8.00e-13), PYAC + OH --> MCO3 + CO2 #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 1.30e-12, 0.0, 500.0e0), HMHP + OH --> 0.5CH2O + 0.5HO2 + 0.5HCOOH + 0.5OH #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 3.14e-12, 0.0, 580.0e0), MCO3 + HO2 --> 0.13O3 + 0.13ACTA + 0.37MAP + 0.5MO2 + 0.5CO2 + 0.5OH #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 1.55e-12, 0.0, 340.0e0), HPETHNL + OH --> CO + OH + CH2O #==2019/11/06; Bates2019; KHB==#
        c(2.91e-11), HPETHNL + OH --> GLYX + OH #==2019/11/06; Bates2019; KHB==#
        arr(t, T, 1.56e-11, 0.0, 117.0e0), NAP + OH --> NRO2 + OH #==2013/08/12; Pye2010; HOTP==#
        arr(t, T, 1.40e-12, 0.0, 700.0e0), NRO2 + HO2 --> LNRO2H + HO2 #==2013/08/12; Pye2010; HOTP==#
        arr(t, T, 2.60e-12, 0.0, 350.0e0), NRO2 + NO --> LNRO2N + NO #==2013/08/12; Pye2010; HOTP==#

        # --- C2H2 & C2H4 chemistry (per KHB)
        arr(t, T, 9.10e-15, 0.0e0, -2580.0e0), C2H4 + O3 --> CH2O + CH2OO #==2021/09/22; Kwon2020; KHB,MSL==#
        arr3(t, T, num_density, 1.10e-28, 3.5e+00, 0.0, 8.4e-12, 1.75e0, 0.0, 0.5e0), C2H4 + OH --> ETOO #==2021/09/22; Kwon2020; KHB,MSL==#
        arr3(t, T, num_density, 5.50e-30, 0.0e0, 0.0, 8.3e-13, -2.0e0, 0.0, 0.5e0), C2H2 + OH --> 0.636GLYX + 0.636OH + 0.364CO + 0.364HO2 + 0.364HCOOH #==2021/09/22; Kwon2020; KHB,MSL==#
        arr(t, T, 1.53e-13, 0.0e0, 1300.0e0), ETOO + HO2 --> ETHP #==2021/09/22; Kwon2020; KHB,MSL==#
        arr(t, T, 2.7e-12, 0.0e+00, 360.0e0), ETOO + NO --> 0.995ETO + 0.995NO2 + 0.005ETHN #==2021/09/22; Kwon2020; KHB,MSL==#
        c(2.3e-12), ETOO + NO3 --> ETO + NO2 #==2021/09/22; Kwon2020; KHB,MSL==#
        c(6.00e-13), ETOO + MO2 --> 0.6ETO + 0.6HO2 + 0.8CH2O + 0.2MOH + 0.2ETHP + 0.2GLYC #==2021/09/22; Kwon2020; KHB,MSL==#
        arrnodim(t, T, 9.5e+13, 0.0e0, -5988.0e0), ETO --> HO2 + 2.000CH2O #==2021/09/22; Kwon2020; KHB,MSL==#
        arr(t, T, 2.5e-14, 0.0e0, -300.0e0), ETO + O2 --> GLYC + HO2 #==2021/09/22; Kwon2020; KHB,MSL==#
        c(8.40e-13), ETHN + OH --> GLYC + NO2 #==2021/09/22; Kwon2020; KHB,MSL==#
        arr(t, T, 1.90e-12, 0.0e+00, 190.0e0), ETHP + OH --> ETOO #==2021/09/22; Kwon2020; KHB,MSL==#
        c(1.38e-11), ETHP + OH --> OH + GLYC #==2021/09/22; Kwon2020; KHB,MSL==#
        #
        # --- Aromatic Chemistry (per KHB)
        arr(t, T, 2.3e-12, 0.0e0, -193.0e0), BENZ + OH --> BRO2 + 0.54PHEN + 0.54HO2 + 0.46AROMRO2 + 0.18GLYX + 0.2CO + 0.56AROMP4 #==2021/09/29; Bates2021b; KHB,MSL==#
        arr(t, T, 1.8e-12, 0.0e0, 340.0e0), TOLU + OH --> TRO2 + 0.19CSL + 0.19HO2 + 0.81AROMRO2 + 0.06BALD + 0.12GLYX + 0.12MGLY + 0.27CO + 0.04MVK + 0.3AROMP5 + 0.68AROMP4 #==2021/09/29; Bates2021b; KHB,MSL==#
        c(1.7e-11), XYLE + OH --> XRO2 + 0.15CSL + 0.15HO2 + 0.85AROMRO2 + 0.06BALD + 0.1GLYX + 0.2MGLY + 0.3CO + 0.04MVK + 0.56AROMP5 + 0.28AROMP4 + 0.45RCOOH #==2021/09/29; Bates2021b; KHB,MSL==#
        regressT(t, T, 0, 2.91e-13 * 0.82e0, 1300.0e0), AROMRO2 + HO2 --> OH + HO2 #==2021/09/29; Bates2021b; KHB,MSL==#
        arr(t, T, 2.60e-12, 0.0e+00, 365.0e0), AROMRO2 + NO --> NO2 + HO2 #==2021/09/29; Bates2021b; KHB,MSL==#
        c(2.30e-12), AROMRO2 + NO3 --> NO2 + HO2 #==2021/09/29; Bates2021b; KHB,MSL==#
        arr(t, T, 1.70e-14, 0.0e0, 220.0e0), AROMRO2 + MO2 --> CH2O + HO2 + HO2 #==2021/09/29; Bates2021b; KHB,MSL==#
        arr(t, T, 4.20e-14, 0.0e0, 220.0e0), AROMRO2 + MCO3 --> MO2 + HO2 #==2021/09/29; Bates2021b; KHB,MSL==#
        arr(t, T, 4.70e-13, 0.0e0, 1220.0e0), PHEN + OH --> 0.06BENZO + 0.06GLYX + 0.18AROMP4 + 0.14AROMRO2 + 0.8MCT + 0.8HO2 #==2021/09/29; Bates2021b; KHB,MSL==#
        c(3.8e-12), PHEN + NO3 --> 0.258NPHEN + 0.742HNO3 + 0.742BENZO #==2021/09/29; Bates2021b; KHB,MSL==#
        c(4.7e-11), CSL + OH --> 0.727MCT + 0.727HO2 + 0.2AROMRO2 + 0.073BENZO + 0.44AROMP5 #==2021/09/29; Bates2021b; KHB,MSL==#
        c(1.4e-11), CSL + NO3 --> 0.5NPHEN + 0.2AROMRO2 + 0.5HNO3 + 0.3BENZO + 0.44AROMP5 #==2021/09/29; Bates2021b; KHB,MSL==#
        c(2.0e-11), MCT + OH --> 0.3BENZO + 0.7AROMRO2 + 1.05AROMP4 #==2021/09/29; Bates2021b; KHB,MSL==#
        c(9.2e-18), MCT + O3 --> GLYC + HO2 + OH + AROMP4 #==2021/09/29; Bates2021b; KHB,MSL==#
        c(9.9e-11), MCT + NO3 --> 0.5NPHEN + 0.5HNO3 + 0.3BENZO + 0.2AROMRO2 + 0.3AROMP4 #==2021/09/29; Bates2021b; KHB,MSL==#
        arr(t, T, 5.90e-12, 0.0e0, 225.0e0), BALD + OH --> BZCO3 #==2021/09/29; Bates2021b; KHB,MSL==#
        c(2.4e-15), BALD + NO3 --> BZCO3 + HNO3 #==2021/09/29; Bates2021b; KHB,MSL==#
        arr(t, T, 1.10e-11, 0.0e0, 340.0e0), BZCO3 + HO2 --> 0.35CO2 + 0.2BENZO2 + 0.15O3 + 0.2OH + 0.15BENZP + 0.65BZCO3H #==2021/09/29; Bates2021b; KHB,MSL==#
        arr(t, T, 7.50e-12, 0.0e0, 290.0e0), BZCO3 + NO --> NO2 + CO2 + BENZO2 #==2021/09/29; Bates2021b; KHB,MSL==#
        rPAN_acac(t, T, num_density, 3.28e-28, -6.87e0, 1.125e-11, -1.105e0, 0.3e0), BZCO3 + NO2 --> BZPAN #==2021/09/29; Bates2021b; KHB,MSL==#
        c(4.66e-12), BZCO3H + OH --> BZCO3 #==2021/09/29; Bates2021b; KHB,MSL==#
        rPAN_ababnodim(t, T, num_density, 1.10e-5, -10100.0e0, 1.90e+17, -14100.0e0, 0.3e0) * 0.67e0, BZPAN --> BZCO3 + NO2 #==2021/09/29; Bates2021b; KHB,MSL==#
        c(1.06e-12), BZPAN + OH --> BENZP + CO2 + NO2 #==2021/09/29; Bates2021b; KHB,MSL==#
        c(7.00e-12), BENZO2 + NO2 --> BENZO + NO3 #==2021/09/29; Bates2021b; KHB,MSL==#
        arr(t, T, 2.670e-12, 0.0e0, 365.0e0), BENZO2 + NO --> BENZO + NO2 #==2021/09/29; Bates2021b; KHB,MSL==#
        c(2.30e-12), BENZO2 + NO3 --> BENZO + NO2 #==2021/09/29; Bates2021b; KHB,MSL==#
        arr(t, T, 2.24e-13, 0.0e0, 1300.0e0), BENZO2 + HO2 --> BENZP #==2021/09/29; Bates2021b; KHB,MSL==#
        c(3.60e-12), BENZP + OH --> BENZO2 #==2021/09/29; Bates2021b; KHB,MSL==#
        c(2.86e-13), BENZO + O3 --> BENZO2 #==2021/09/29; Bates2021b; KHB,MSL==#
        c(2.08e-12), BENZO + NO2 --> NPHEN #==2021/09/29; Bates2021b; KHB,MSL==#
        c(3.47e-12), NPHEN + OH --> 0.5R4N1 + AROMP4 + 0.5NO2 #==2021/09/29; Bates2021b; KHB,MSL==#
        c(2.60e-12), NPHEN + NO3 --> 0.5HNO3 + NO2 + 0.5R4N1 + AROMP4 #==2021/09/29; Bates2021b; KHB,MSL==#
        arr(t, T, 2.670e-13, 0.0e0, 365.0e0), BENZO2 + MO2 --> BENZO + HO2 + CH2O #==2021/09/29; Bates2021b; KHB,MSL==#
        arr(t, T, 2.670e-12, 0.0e0, 365.0e0), BZCO3 + MO2 --> BENZO2 + CO2 + HO2 + CH2O #==2021/09/29; Bates2021b; KHB,MSL==#
        c(5.0e-11), AROMP4 + OH --> 0.6GLYX + 0.25CO + 0.25HCOOH + 0.25OH + 0.33HO2 + 0.33RCO3 + 0.45RCOOH #==2021/09/29; Bates2021b; KHB,MSL==#
        c(8.0e-16), AROMP4 + O3 --> 0.5HCOOH + 0.5CO + 0.6GLYX + 0.9GLYC + 0.1HO2 + 0.1OH #==2021/09/29; Bates2021b; KHB,MSL==#
        c(1.5e-3), AROMP4 --> 0.2HO2 + 0.2GLYX + 1.2RCHO #==2021/09/29; Bates2021b; KHB,MSL==#
        c(5.0e-11), AROMP5 + OH --> 0.6MGLY + 0.15ACTA + 0.1HCOOH + 0.25OH + 0.33HO2 + 0.33RCO3 + 0.25CO + 0.52RCOOH #==2021/09/29; Bates2021b; KHB,MSL==#
        c(8.0e-16), AROMP5 + O3 --> 0.6MGLY + 0.3ACTA + 0.2HCOOH + 0.5CO + 0.95GLYC + 0.1HO2 + 0.1OH #==2021/09/29; Bates2021b; KHB,MSL==#
        c(1.5e-3), AROMP5 --> 0.2HO2 + 0.2R4O2 + 0.2MGLY + 1.2RCHO #==2021/09/29; Bates2021b; KHB,MSL==#
        #
        # KHB -- "we still need to include the dummy species for aromatic oxidation
        #         to make the complex SOA code work. Hopefully this will be changed
        #         very soon when Jared Brewer updates the aromatic SOA, but I think it's
        #         still necessary, in which case, we need to add the following reactions too.
        #         (If I'm wrong, we can delete XRO2, TRO2, BRO2, LXRO2N, LXRO2H,
        #         LTRO2N, LTRO2H, LBRO2N, and LBRO2H from the species list, delete
        #         XRO2, TRO2, and BRO2 as products from the BENZ + OH, TOLU + OH,
        #         and XYLE + OH reactions above, and not include the following reactions)"
        #
        arr(t, T, 1.40e-12, 0.0e0, 700.0e0), BRO2 + HO2 --> HO2 + LBRO2H #==2021/09/29; Bates2021b; KHB,MSL==#
        arr(t, T, 2.60e-12, 0.0e0, 350.0e0), BRO2 + NO --> NO + LBRO2N #==2021/09/29; Bates2021b; KHB,MSL==#
        arr(t, T, 1.40e-12, 0.0e0, 700.0e0), TRO2 + HO2 --> HO2 + LTRO2H #==2021/09/29; Bates2021b; KHB,MSL==#
        arr(t, T, 2.60e-12, 0.0e0, 350.0e0), TRO2 + NO --> NO + LTRO2N #==2021/09/29; Bates2021b; KHB,MSL==#
        arr(t, T, 1.40e-12, 0.0e0, 700.0e0), XRO2 + HO2 --> HO2 + LXRO2H #==2021/09/29; Bates2021b; KHB,MSL==#
        arr(t, T, 2.60e-12, 0.0e0, 350.0e0), XRO2 + NO --> NO + LXRO2N #==2021/09/29; Bates2021b; KHB,MSL==#
        c(1.20e-12), MO2 + NO3 --> NO2 + CH2O + HO2 #==2022/10/18: IUPAC ROO_19; KHB,BMY==#
        #
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # %%%%% Heterogeneous chemistry reactions                               %%%%%
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #
        # HO2uptk1stOrd( State_Het ), HO2 --> H2O #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        # NO2uptk1stOrdAndCloud( State_Het ), NO2 --> 0.500HNO3 + 0.500HNO2
        # NO3uptk1stOrdAndCloud( State_Het ), NO3 --> HNO3
        # NO3hypsisClonSALA( State_Het ), NO3 --> NIT #==2018/03/16; XW==#
        # NO3hypsisClonSALC( State_Het ), NO3 --> NITs #==2018/03/16; XW==#
        # N2O5uptkByH2O( State_Het ), N2O5 + H2O --> 2.000HNO3
        # N2O5uptkByStratHCl( State_Het ), N2O5 + HCl --> ClNO2 + HNO3 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # N2O5uptkByCloud( State_Het ), N2O5 --> 2.000HNO3 #==2018/10/17; Cloud uptake, CDH==#
        # N2O5uptkBySALACl( State_Het ), N2O5 + SALACL --> ClNO2 + HNO3 #==2018/01/19; Sherwen2017;TS,JAS,SDE,XW==#
        # N2O5uptkBySALCCl( State_Het ), N2O5 + SALCCL --> ClNO2 + HNO3 #==2018/01/19; Sherwen2017;TS,JAS,SDE,XW==#
        # OHuptkBySALACl( State_Het ), OH + SALACL --> 0.500Cl2 #==2018/03/12; XW==#
        # OHuptkBySALCCl( State_Het ), OH + SALCCL --> 0.500Cl2 #==2018/03/12; XW==#
        # BrNO3uptkByH2O( State_Het ), BrNO3 + H2O --> HOBr + HNO3 #==2014/02/03; Eastham2014; SDE==#
        # BrNO3uptkByHCl( State_Het ), BrNO3 + HCl --> BrCl + HNO3 #==2014/02/03; Eastham2014; SDE==#
        # ClNO3uptkByH2O( State_Het ), ClNO3 + H2O --> HOCl + HNO3 #==2014/02/03; Eastham2014; SDE==#
        # ClNO3uptkByHCl( State_Het ), ClNO3 + HCl --> Cl2 + HNO3 #==2014/02/03; Eastham2014; SDE==#
        # ClNO3uptkByHBr( State_Het ), ClNO3 + HBr --> BrCl + HNO3 #==2014/02/03; Eastham2014; SDE==#
        # ClNO3uptkByBrSALA( State_Het ), ClNO3 + BrSALA --> BrCl + HNO3 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # ClNO3uptkByBrSALC( State_Het ), ClNO3 + BrSALC --> BrCl + HNO3 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # ClNO3uptkBySALACL( State_Het ), ClNO3 + SALACL --> Cl2 + HNO3 #==2018/01/22; XW==#
        # ClNO3uptkBySALCCL( State_Het ), ClNO3 + SALCCL --> Cl2 + HNO3 #==2018/01/22; XW==#
        # ClNO2uptkBySALACL( State_Het ), ClNO2 + SALACL --> Cl2 + HNO2 #==2018/01/22; XW==#
        # ClNO2uptkBySALCCL( State_Het ), ClNO2 + SALCCL --> Cl2 + HNO2 #==2018/01/22; XW==#
        # ClNO2uptkByHCl( State_Het ), ClNO2 + HCl --> Cl2 + HNO2 #==2018/01/22; XW==#
        # ClNO2uptkByBrSALA( State_Het ), ClNO2 + BrSALA --> BrCl + HNO2 #==2018/01/22; XW==#
        # ClNO2uptkByBrSALC( State_Het ), ClNO2 + BrSALC --> BrCl + HNO2 #==2018/01/22; XW==#
        # ClNO2uptkByHBr( State_Het ), ClNO2 + HBr --> BrCl + HNO2 #==2018/01/22; XW==#
        # HOClUptkByHCl( State_Het ), HOCl + HCl --> Cl2 + H2O #==2014/02/03; Eastham2014; SDE==#
        # HOClUptkByHBr( State_Het ), HOCl + HBr --> BrCl + H2O #==2014/02/03; Eastham2014; SDE==#
        # HOClUptkBySALACL( State_Het ), HOCl + SALACL --> Cl2 + H2O #==2018/01/22; XW==#
        # HOClUptkBySALCCL( State_Het ), HOCl + SALCCL --> Cl2 + H2O #==2018/01/22; XW==#
        # HOClUptkByHSO3m( State_Het ) + HOClUptkBySO3mm( State_Het ), HOCl + SO2 --> SO4 + HCl #==2018/11/08; XW; June 6, 2021, MSL==#
        # HOBrUptkByHBr( State_Het ), HOBr + HBr --> Br2 + H2O #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # HOBrUptkByHCl( State_Het ), HOBr + HCl --> BrCl + H2O #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # HOBrUptkBySALACL( State_Het ), HOBr + SALACL --> BrCl + H2O #==2018/01/22; Sherwen2017;TS,JAS,SDE;XW==#
        # HOBrUptkBySALCCL( State_Het ), HOBr + SALCCL --> BrCl + H2O #==2018/01/22; Sherwen2017;TS,JAS,SDE,XW==#
        # HOBrUptkByBrSALA( State_Het ), HOBr + BrSALA --> Br2 #==2017/09/22; Sherwen2017;TS,JAS,SDE==#
        # HOBrUptkByBrSALC( State_Het ), HOBr + BrSALC --> Br2 #==2017/09/22; Sherwen2017;TS,JAS,SDE==#
        # HOBrUptkByHSO3m( State_Het ) + HOBrUptkBySO3mm( State_Het ), HOBr + SO2 --> SO4 + HBr #==2017/11/15; Chen2017; QJC; June 6, 2021, MSL==#
        # O3uptkByHBr( State_Het ), O3 + HBr --> HOBr #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # O3uptkByBrSALA( State_Het ), O3 + BrSALA --> HOBr #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # O3uptkByBrSALC( State_Het ), O3 + BrSALC --> HOBr #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # HBrUptkBySALA( State_Het ), HBr --> BrSALA #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # HBrUptkBySALC( State_Het ), HBr --> BrSALC #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # IuptkBySulf1stOrd( SR_MW(ind_HI), 0.10, State_Het ), HI --> AERI #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # IuptkBySALA1stOrd( SR_MW(ind_HI), 0.10, State_Het ), HI --> ISALA #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # IuptkBySALC1stOrd( SR_MW(ind_HI), 0.10, State_Het ), HI --> ISALC #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # IuptkBySulf1stOrd( SR_MW(ind_I2O2), 0.02, State_Het ), I2O2 --> 2.000AERI #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # IuptkBySALA1stOrd( SR_MW(ind_I2O2), 0.02, State_Het ), I2O2 --> 2.000ISALA #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # IuptkBySALC1stOrd( SR_MW(ind_I2O2), 0.02, State_Het ), I2O2 --> 2.000ISALC #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # IuptkBySulf1stOrd( SR_MW(ind_I2O3), 0.02, State_Het ), I2O3 --> 2.000AERI #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # IuptkBySALA1stOrd( SR_MW(ind_I2O3), 0.02, State_Het ), I2O3 --> 2.000ISALA #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # IuptkBySALC1stOrd( SR_MW(ind_I2O3), 0.02, State_Het ), I2O3 --> 2.000ISALC #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # IuptkBySulf1stOrd( SR_MW(ind_I2O4), 0.02, State_Het ), I2O4 --> 2.000AERI #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # IuptkBySALA1stOrd( SR_MW(ind_I2O4), 0.02, State_Het ), I2O4 --> 2.000ISALA #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # IuptkBySALC1stOrd( SR_MW(ind_I2O4), 0.02, State_Het ), I2O4 --> 2.000ISALC #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        # IONO2uptkByH2O( State_Het ), IONO2 + H2O --> HOI + HNO3 #==2021/09/16 XW, TSherwen==#
        # IbrkdnByAcidBrSALA( SR_MW(ind_IONO), C(ind_IONO), 0.02, State_Het ), IONO + BrSALA --> IBr + HNO2 #==2017/09/22; Sherwen2017;TS,JAS,SDE,XW==#
        # IbrkdnByAcidBrSALC( SR_MW(ind_IONO), C(ind_IONO), 0.02, State_Het ), IONO + BrSALC --> IBr + HNO2 #==2017/09/22; Sherwen2017;TS,JAS,SDE,XW==#
        # IbrkdnByAcidSALACl( SR_MW(ind_IONO), C(ind_IONO), 0.02, State_Het ), IONO + SALACL --> ICl + HNO2 #==2017/09/22; Sherwen2017;TS,JAS,SDE,XW==#
        # IbrkdnByAcidSALCCl( SR_MW(ind_IONO), C(ind_IONO), 0.02, State_Het ), IONO + SALCCL --> ICl + HNO2 #==2017/09/22; Sherwen2017;TS,JAS,SDE,XW==#
        # IbrkdnByAcidBrSALA( SR_MW(ind_IONO2), C(ind_IONO2), 0.01, State_Het ), IONO2 + BrSALA --> IBr + HNO3 #==2017/09/22; Sherwen2017;TS,JAS,SDE,XW==#
        # IbrkdnByAcidBrSALC( SR_MW(ind_IONO2), C(ind_IONO2), 0.01, State_Het ), IONO2 + BrSALC --> IBr + HNO3 #==2017/09/22; Sherwen2017;TS,JAS,SDE,XW==#
        # IbrkdnByAcidSALACl( SR_MW(ind_IONO2), C(ind_IONO2), 0.01, State_Het ), IONO2 + SALACL --> ICl + HNO3 #==2017/09/22; Sherwen2017;TS,JAS,SDE,XW==#
        # IbrkdnByAcidSALCCl( SR_MW(ind_IONO2), C(ind_IONO2), 0.01, State_Het ), IONO2 + SALCCL --> ICl + HNO3 #==2017/09/22; Sherwen2017;TS,JAS,SDE,XW==#
        # IbrkdnByAcidBrSALA( SR_MW(ind_HOI), C(ind_HOI), 0.01, State_Het ), HOI + BrSALA --> IBr #==2017/09/22; Sherwen2017;TS,JAS,SDE,XW==#
        # IbrkdnByAcidBrSALC( SR_MW(ind_HOI), C(ind_HOI), 0.01, State_Het ), HOI + BrSALC --> IBr #==2017/09/22; Sherwen2017;TS,JAS,SDE,XW==#
        # IbrkdnByAcidSALACl( SR_MW(ind_HOI), C(ind_HOI), 0.01, State_Het ), HOI + SALACL --> ICl #==2017/09/22; Sherwen2017;TS,JAS,SDE,XW==#
        # IbrkdnByAcidSALCCl( SR_MW(ind_HOI), C(ind_HOI), 0.01, State_Het ), HOI + SALCCL --> ICl #==2017/09/22; Sherwen2017;TS,JAS,SDE,XW==#
        # GLYXuptk1stOrd( SR_MW(ind_GLYX), State_Het), GLYX --> SOAGX #==2017/06/15; Marais2016, EAM==#
        # MGLYuptk1stOrd( SR_MW(ind_MGLY), State_Het), MGLY --> SOAGX #==2017/06/15; Marais2016, EAM==#
        # IEPOXuptk1stOrd( SR_MW(ind_IEPOXA), false, State_Het ), IEPOXA --> SOAIE #==2017/06/15; Marais2016, EAM==#
        # IEPOXuptk1stOrd( SR_MW(ind_IEPOXB), false, State_Het ), IEPOXB --> SOAIE #==2017/06/15; Marais2016, EAM==#
        # IEPOXuptk1stOrd( SR_MW(ind_IEPOXD), false, State_Het ), IEPOXD --> SOAIE #==2017/06/15; Marais2016, EAM==#
        # VOCuptk1stOrd( SR_MW(ind_LVOC), 1.0, State_Het ), LVOC --> LVOCOA #==2017/06/15; Marais2016, EAM==#
        # VOCuptk1stOrd( SR_MW(ind_MVKN), 5.0E-3, State_Het ), MVKN --> IONITA #==2017/06/15; Marais2016, EAM==#
        # VOCuptk1stOrd( SR_MW(ind_R4N2), 5.0E-3, State_Het ), R4N2 --> IONITA #==2017/06/15; Marais2016, EAM==#
        # VOCuptk1stOrd( SR_MW(ind_MONITS), 1.0E-2, State_Het ), MONITS --> MONITA #==2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS==#
        # VOCuptk1stOrd( SR_MW(ind_MONITU), 1.0E-2, State_Het ), MONITU --> MONITA #==2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS==#
        # VOCuptk1stOrd( SR_MW(ind_HONIT), 1.0E-2, State_Het ), HONIT --> MONITA #==2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS==#
        # MGLYuptk1stOrd( SR_MW(ind_PYAC), State_Het ), PYAC --> SOAGX #==2019/11/06; Bates2019; KHB==#
        # IEPOXuptk1stOrd( SR_MW(ind_HMML), true, State_Het), HMML --> SOAIE #==2019/11/06; Bates2019; KHB==#
        # VOCuptk1stOrd( SR_MW(ind_IHN1), 5.0E-3, State_Het ), IHN1 --> IONITA #==2019/11/06; Bates2019; KHB==#
        # VOCuptk1stOrd( SR_MW(ind_IHN2), 5.0E-2, State_Het ), IHN2 --> IONITA #==2019/11/06; Bates2019; KHB==#
        # VOCuptk1stOrd( SR_MW(ind_IHN3), 5.0E-3, State_Het ), IHN3 --> IONITA #==2019/11/06; Bates2019; KHB==#
        # VOCuptk1stOrd( SR_MW(ind_IHN4), 5.0E-3, State_Het ), IHN4 --> IONITA #==2019/11/06; Bates2019; KHB==#
        # IEPOXuptk1stOrd( SR_MW(ind_ICHE), false, State_Het ), ICHE --> SOAIE #==2019/11/06; Bates2019; KHB==#
        # VOCuptk1stOrd( SR_MW(ind_INPD), 5.0E-3, State_Het ), INPD --> IONITA #==2019/11/06; Bates2019; KHB==#
        # VOCuptk1stOrd( SR_MW(ind_INPB), 5.0E-3, State_Het ), INPB --> IONITA #==2019/11/06; Bates2019; KHB==#
        # VOCuptk1stOrd( SR_MW(ind_IDN), 5.0E-3, State_Het ), IDN --> IONITA #==2019/11/06; Bates2019; KHB==#
        # VOCuptk1stOrd( SR_MW(ind_ITCN), 5.0E-3, State_Het ), ITCN --> IONITA #==2019/11/06; Bates2019; KHB==#
        # VOCuptk1stOrd( SR_MW(ind_ITHN), 5.0E-3, State_Het ), ITHN --> IONITA #==2019/11/06; Bates2019; KHB==#
        # VOCuptk1stOrd( SR_MW(ind_MCRHNB), 5.0E-3, State_Het ), MCRHNB --> IONITA #==2019/11/06; Bates2019; KHB==#
        # VOCuptk1stOrd( SR_MW(ind_MCRHN), 5.0E-3, State_Het ), MCRHN --> IONITA #==2019/11/06; Bates2019; KHB==#
        # VOCuptk1stOrd( SR_MW(ind_NPHEN), 1.0E-2, State_Het ), NPHEN --> AONITA #==2021/09/29; Bates2021b; KHB,MSL==#
        #
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # %%%%% Photolysis reactions                                            %%%%%
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        #
        # NOTE(CT): I am removing hv from the photolysis reactions (e.g. j_2, O3 + hv --> O + O2 is changed to j_2, O3 --> O + O2)
        #   because hv is not a chemical species.
        j_2, O3 --> O + O2 #==2014/02/03; Eastham2014; SDE==#
        j_3, O3 --> O1D + O2 #==2014/02/03; Eastham2014; SDE==#
        j_1, O2 --> 2.000O #==2014/02/03; Eastham2014; SDE==#
        j_11, NO2 --> NO + O #==2014/02/03; Eastham2014; SDE==#
        j_9, H2O2 --> OH + OH
        j_10, MP --> CH2O + HO2 + OH
        j_7, CH2O --> HO2 + H + CO #==2014/02/03; Eastham2014; SDE==#
        j_8, CH2O --> H2 + CO
        j_16, HNO3 --> OH + NO2
        j_15, HNO2 --> OH + NO
        j_17, HNO4 --> OH + NO3
        j_18, HNO4 --> HO2 + NO2
        j_12, NO3 --> NO2 + O #==2014/02/03; Eastham2014; SDE==#
        j_13, NO3 --> NO + O2
        j_14, N2O5 --> NO3 + NO2
        j_61, ALD2 --> 0.880MO2 + HO2 + 0.880CO + 0.120MCO3 #==2014/12/19; FAST-JX v7.0 fix; JMAO==#
        j_62, ALD2 --> CH4 + CO
        j_59, PAN --> 0.700MCO3 + 0.700NO2 + 0.300MO2 + 0.300NO3 #==2014/05/23; Eastham2014; JMAO,SDE==#
        j_70, RCHO --> 0.500OTHRO2 + HO2 + CO + 0.070A3O2 + 0.270B3O2 #==2019/05/10; Fisher2018; JAF==#
        j_76, ACET --> MCO3 + MO2
        j_77, ACET --> 2.000MO2 + CO
        j_69, MEK --> 0.850MCO3 + 0.425OTHRO2 + 0.150MO2 + 0.150RCO3 + 0.060A3O2 + 0.230B3O2 #==2019/05/10; Fisher2018; JAF==#
        j_68, GLYC --> 0.900CH2O + 1.730HO2 + CO + 0.070OH + 0.100MOH #==2014/05/23; Eastham2014; JMAO,SDE==#
        j_72, GLYX --> 2.000HO2 + 2.000CO
        j_73, GLYX --> H2 + 2.000CO
        j_74, GLYX --> CH2O + CO
        j_71, MGLY --> MCO3 + CO + HO2
        j_63, MVK --> PRPE + CO
        j_64, MVK --> MCO3 + CH2O + CO + HO2
        j_65, MVK --> MO2 + RCO3 #==2014/05/23; Eastham2014; JMAO,SDE==#
        j_66, MACR --> CO + HO2 + CH2O + MCO3 #==2014/05/23; Eastham2014; JMAO,SDE==#
        j_75, HAC --> MCO3 + CH2O + HO2
        j_79, PRPN --> OH + HO2 + RCHO + NO2
        j_80, ETP --> OH + HO2 + ALD2
        j_81, RA3P --> OH + HO2 + RCHO
        j_82, RB3P --> OH + HO2 + ACET
        j_83, R4P --> OH + HO2 + RCHO
        j_84, PP --> OH + HO2 + ALD2 + CH2O
        j_85, RP --> OH + HO2 + ALD2
        j_98, R4N2 --> NO2 + 0.320ACET + 0.190MEK + 0.180MO2 + 0.270HO2 + 0.320ALD2 + 0.130RCHO + 0.050A3O2 + 0.180B3O2 + 0.320OTHRO2
        j_99, MAP --> OH + MO2
        j_23, Br2 --> 2.000Br #==2012/06/07; Parrella2012; JPP==#
        j_28, BrO --> Br + O #==2014/02/03; Eastham2014; SDE==#
        j_32, HOBr --> Br + OH #==2012/06/07; Parrella2012; JPP==#
        j_29, BrNO3 --> Br + NO3 #==2012/06/07; Parrella2012; JPP==#
        j_30, BrNO3 --> BrO + NO2 #==2012/06/07; Parrella2012; JPP==#
        j_31, BrNO2 --> Br + NO2 #==2012/06/07; Parrella2012; JPP==#
        j_56, CHBr3 --> 3.000Br #==2012/06/07; Parrella2012; JPP==#
        j_55, CH2Br2 --> 2.000Br #==2014/02/03; Eastham2014; SDE==#
        j_50, CH3Br --> MO2 + Br #==2014/02/03; Eastham2014; SDE==#
        j_43, CH3Cl --> MO2 + Cl #==2014/02/03; Eastham2014; SDE==#
        j_45, CH2Cl2 --> 2.000Cl #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        j_33, BrCl --> Br + Cl #==2014/02/03; Eastham2014; SDE==#
        j_22, Cl2 --> 2.000Cl #==2014/02/03; Eastham2014; SDE==#
        j_27, ClO --> Cl + O #==2014/02/03; Eastham2014; SDE==#
        j_25, OClO --> ClO + O #==2014/02/03; Eastham2014; SDE==#
        j_26, Cl2O2 --> Cl + ClOO #==2014/02/03; Eastham2014; SDE==#
        j_21, ClNO2 --> Cl + NO2 #==2014/02/03; Eastham2014; SDE==#
        j_19, ClNO3 --> Cl + NO3 #==2014/02/03; Eastham2014; SDE==#
        j_20, ClNO3 --> ClO + NO2 #==2014/02/03; Eastham2014; SDE==#
        j_24, HOCl --> Cl + OH #==2014/02/03; Eastham2014; SDE==#
        j_44, CH3CCl3 --> 3.000Cl #==2014/02/03; Eastham2014; SDE==#
        j_42, CCl4 --> 4.000Cl #==2014/02/03; Eastham2014; SDE==#
        j_37, CFC11 --> 3.000Cl #==2014/02/03; Eastham2014; SDE==#
        j_38, CFC12 --> 2.000Cl #==2014/02/03; Eastham2014; SDE==#
        j_39, CFC113 --> 3.000Cl #==2014/02/03; Eastham2014; SDE==#
        j_40, CFC114 --> 2.000Cl #==2014/02/03; Eastham2014; SDE==#
        j_41, CFC115 --> Cl #==2014/02/03; Eastham2014; SDE==#
        j_47, HCFC123 --> 2.000Cl #==2014/02/03; Eastham2014; SDE==#
        j_48, HCFC141b --> 2.000Cl #==2014/02/03; Eastham2014; SDE==#
        j_49, HCFC142b --> Cl #==2014/02/03; Eastham2014; SDE==#
        j_46, HCFC22 --> Cl #==2014/02/03; Eastham2014; SDE==#
        j_53, H1301 --> Br #==2014/02/03; Eastham2014; SDE==#
        j_51, H1211 --> Cl + Br #==2014/02/03; Eastham2014; SDE==#
        j_54, H2402 --> 2.000Br #==2014/02/03; Eastham2014; SDE==#
        j_101, ClOO --> Cl + O2 #==2014/02/03; Eastham2014; SDE==#
        j_114, I2 --> 2.000I #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        j_115, HOI --> I + OH #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        j_116, IO --> I + O #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        j_117, OIO --> I + O2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        j_118, INO --> I + NO #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        j_119, IONO --> I + NO2 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        j_120, IONO2 --> I + NO3 #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        j_121, I2O2 --> I + OIO #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        j_122, CH3I --> I #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        j_123, CH2I2 --> 2.000I #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        j_124, CH2ICl --> I + Cl #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        j_125, CH2IBr --> I + Br #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        j_126, I2O4 --> 2.000OIO #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        j_127, I2O3 --> OIO + IO #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        j_128, IBr --> I + Br #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        j_129, ICl --> I + Cl #==2017/09/22; Sherwen2016b;TS,JAS,SDE==#
        j_103, MPN --> CH2O + NO3 + HO2 #==2012/02/12; Browne2011; ECB==#
        j_104, MPN --> MO2 + NO2 #==2012/02/12; Browne2011; ECB==#
        j_97, ATOOH --> OH + CH2O + MCO3 #==2013/03/22; Paulot2009; FP,EAM,JMAO,MJE==#
        j_36, N2O --> N2 + O1D #==2014/02/03; Eastham2014; SDE==#
        j_34, OCS --> SO2 + CO #==2014/02/03; Eastham2014; SDE==#
        j_100, SO4 --> SO2 + 2.000OH #==2014/02/03; Eastham2014; SDE==#
        j_6, NO --> O + N #==2014/02/03; Eastham2014; SDE==#
        j_105, PIP --> RCHO + OH + HO2 #==2017/03/23; Fischer2014; EVF==#
        j_107, ETHLN --> NO2 + CH2O + CO + HO2 #==2017/06/15; Marais2016; EAM==#
        j_111, MONITS --> MEK + NO2 #==2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS==#
        j_112, MONITU --> RCHO + NO2 #==2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS==#
        j_113, HONIT --> HAC + NO2 #==2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS==#
        j_130, NITs --> HNO2 #==2018/07/19; Kasibhatla2018; PK, TMS==#
        j_131, NITs --> NO2 #==2018/07/19; Kasibhatla2018; PK, TMS==#
        j_132, NIT --> HNO2 #==2018/07/19; Kasibhatla2018; PK, TMS==#
        j_133, NIT --> NO2 #==2018/07/19; Kasibhatla2018; PK, TMS==#
        j_134, MENO3 --> NO2 + HO2 + CH2O #==2019/07/11; Fisher2018; JAF==#
        j_135, ETNO3 --> NO2 + HO2 + ALD2 #==2019/07/11; Fisher2018; JAF==#
        j_136, IPRNO3 --> NO2 + HO2 + ACET #==2019/07/11; Fisher2018; JAF==#
        j_137, NPRNO3 --> NO2 + HO2 + RCHO #==2019/07/11; Fisher2018; JAF==#
        j_86, HMHP --> 2OH + CH2O #==2019/11/06; Bates2019; KHB==#
        j_87, HPETHNL --> OH + CO + HO2 + CH2O #==2019/11/06; Bates2019; KHB==#
        j_88, PYAC --> MCO3 + CO2 + HO2 #==2019/11/06; Bates2019; KHB==#
        j_89, PROPNN --> NO2 + CH2O + MCO3 #==2019/11/06; Bates2019; KHB==#
        j_90, MVKHC --> CO + HO2 + CH2O + MCO3 #==2019/11/06; Bates2019; KHB==#
        j_91, MVKHCB --> 0.5GLYX + 1.5HO2 + 0.5MCO3 + 0.5CO + 0.5MGLY #==2019/11/06; Bates2019; KHB==#
        j_92, MVKHP --> 0.53MCO3 + 0.53GLYC + OH + 0.47HO2 + 0.47CH2O + 0.47MGLY #==2019/11/06; Bates2019; KHB==#
        j_93, MVKPC --> OH + 0.571CO + 0.571MGLY + 0.571HO2 + 0.429GLYX + 0.429MCO3 #==2019/11/06; Bates2019; KHB==#
        j_94, MCRENOL --> 0.875CO + 0.75PYAC + 1.75OH + 0.125MGLY + 0.125HO2 + 0.125MCO3 + 0.125GLYX #==2019/11/06; Bates2019; KHB==#
        j_95, MCRHP --> OH + 0.77CO + HO2 + 0.77HAC + 0.23MGLY + 0.23CH2O #==2019/11/06; Bates2019; KHB==#
        j_96, MACR1OOH --> 0.75OH + 1.238CO2 + 0.488MO2 + 0.75CH2O + 0.262MCO3 + 0.25MACR1OOH #==2019/11/06; Bates2019; KHB==#
        j_108, MVKN --> 0.290HO2 + 0.010OH + 0.700NO2 + 1.010MCO3 + 0.690GLYC + 0.300ETHLN #==2019/11/06; Bates2019; KHB==#
        j_109, MCRHN --> HAC + CO + HO2 + NO2 #==2019/11/06; Bates2019; KHB==#
        j_110, MCRHNB --> PROPNN + OH + CO + HO2 #==2019/11/06; Bates2019; KHB==#
        j_138, RIPA --> MVK + CH2O + HO2 + OH #==2019/11/06; Bates2019; KHB==#
        j_139, RIPB --> MACR + CH2O + HO2 + OH #==2019/11/06; Bates2019; KHB==#
        j_140, RIPC --> OH + HO2 + HC5A #==2019/11/06; Bates2019; KHB==#
        j_141, RIPD --> OH + HO2 + HC5A #==2019/11/06; Bates2019; KHB==#
        j_142, HPALD1 --> 0.888CO + 1.662OH + 0.112HO2 + 0.112IDC + 0.112MVKPC + 0.552MCRENOL + 0.224C4HVP1 #==2019/11/06; Bates2019; KHB==#
        j_143, HPALD2 --> 0.818CO + 1.637OH + 0.182HO2 + 0.182IDC + 0.182MVKPC + 0.455MCRENOL + 0.182C4HVP2 #==2019/11/06; Bates2019; KHB==#
        j_144, HPALD3 --> CO + OH + HO2 + MVK #==2019/11/06; Bates2019; KHB==#
        j_145, HPALD4 --> CO + OH + HO2 + MACR #==2019/11/06; Bates2019; KHB==#
        j_146, IHN1 --> NO2 + 0.45HC5A + 0.45HO2 + 0.55MVKHP + 0.55CO + 0.55OH #==2019/11/06; Bates2019; KHB==#
        j_147, IHN2 --> NO2 + MVK + HO2 + CH2O #==2019/11/06; Bates2019; KHB==#
        j_148, IHN3 --> NO2 + MACR + HO2 + CH2O #==2019/11/06; Bates2019; KHB==#
        j_149, IHN4 --> NO2 + 0.45HC5A + 0.45HO2 + 0.55MCRHP + 0.55CO + 0.55OH #==2019/11/06; Bates2019; KHB==#
        j_150, INPB --> NO2 + CH2O + 0.097MACR + 0.903MVK + 0.67OH + 0.33HO2 #==2019/11/06; Bates2019; KHB==#
        j_151, INPD --> OH + 0.159HO2 + 0.159ICN + 0.841INA #==2019/11/06; Bates2019; KHB==#
        j_152, INPD --> NO2 + 0.841IHOO1 + 0.159IHOO4 #==2019/11/06; Bates2019; KHB==#
        j_106, ICN --> NO2 + 0.839CO + 0.645OH + 0.161HO2 + 0.161IDC + 0.162MVKPC + 0.481MCRENOL + 0.128C4HVP2 + 0.068C4HVP1 #==2019/11/06; Bates2019; KHB==#
        j_78, IDN --> 1.555NO2 + 0.5GLYC + 0.5HAC + 0.05MVK + 0.005MACR + 0.055CH2O + 0.227INA + 0.228ICN + 0.228HO2 #==2019/11/06; Bates2019; KHB==#
        j_153, ICPDH --> CO + 1.5HO2 + 0.5OH + 0.5MCRHP + 0.35MVKDH + 0.15MCRDH #==2019/11/06; Bates2019; KHB==#
        j_154, ICPDH --> OH + HO2 + 0.122CO + 0.1CH2O + 0.1MVKHCB + 0.438HAC + 0.438GLYX + 0.088GLYC + 0.088MGLY + 0.122MCRDH #==2019/11/06; Bates2019; KHB==#
        j_155, IDHDP --> 1.25OH + 0.25GLYC + 0.25HAC + 0.75ICPDH + 0.75HO2 #==2019/11/06; Bates2019; KHB==#
        j_156, IDHPE --> OH + HO2 + 0.429MGLY + 0.429GLYC + 0.571GLYX + 0.571HAC #==2019/11/06; Bates2019; KHB==#
        j_157, IDCHP --> 0.546OH + CO + 1.454HO2 + 0.391MVKHC + 0.155MVKHCB + 0.454MVKPC #==2019/11/06; Bates2019; KHB==#
        j_158, ITHN --> OH + 0.7HO2 + 0.55CH2O + 0.5MCRHN + 0.3GLYC + 0.45HAC + 0.3NO2 + 0.15ETHLN + 0.05MVKN #==2019/11/06; Bates2019; KHB==#
        j_159, ITHN --> NO2 + 0.8HAC + 0.7HO2 + 0.5HPETHNL + 0.35GLYC + 0.15CH2O + 0.15MCRHP + 0.05ATOOH + 0.3OH #==2019/11/06; Bates2019; KHB==#
        j_160, ITCN --> MGLY + OH + NO2 + GLYC #==2019/11/06; Bates2019; KHB==#
        j_161, ITCN --> 0.5MVKHP + 0.5MCRHP + CO + NO2 + HO2 #==2019/11/06; Bates2019; KHB==#
        j_162, ETHP --> ETO + OH #==2021/09/22; Bates2021a; KHB,MSL==#
        j_163, BALD --> BENZO2 + CO + HO2 #==2021/09/29; Bates2021b; KHB,MSL==#
        j_164, BZCO3H --> BENZO2 + OH + CO2 #==2021/09/29; Bates2021b; KHB,MSL==#
        j_165, BENZP --> BENZO #==2021/09/29; Bates2021b; KHB,MSL==#
        j_166, NPHEN --> HNO2 + CO + CO2 + AROMP4 + HO2 #==2021/09/29; Bates2021b; KHB,MSL==#
    end
    rxns = compose(rx_sys, rate_systems)
    if rxn_sys
        return rxns
    end
    sys = convert(ODESystem, rxns, combinatoric_ratelaws=false, name=name,
        metadata=Dict(:coupletype => GEOSChemGasPhaseCoupler))
end

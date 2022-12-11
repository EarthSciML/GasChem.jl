"""
This julia file define the model containing Gas-phase part of full chemistry model.
The rates fuctions are in "lowfunctions.jl", "ratesfunctions.jl", and "RateLawUtilFuncs.jl"
The original model is form:
https://github.com/geoschem

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
* JPL 10-6:     JPL Publication 10-6, https://jpldataeval.jpl.nasa.gov/previous_evaluations.html, 2011.
* JPL 15-10:    JPL Publication 15-10, https://jpldataeval.jpl.nasa.gov, 2015.
* Kwon2020:     Kwon et al, Elementa, https://doi.org/10.1525/elementa.2021.00109, 2020.
* Lee2014:      Lee et al., J. Phys. Chem. A, doi:10.1021/jp4107603, 2014.
* Marais2016:   Marais et al., Atmos. Chem. Phys, doi:10.5194/acp-16-1603-2016, 2016.
* Miller2017:   Miller et al., Atmos. Chem. Phys. Discuss., doi:10.5194/acp-2016-1042, 2017.
* Millet2015:   Millet et al., Atmos. Chem. Phys., doi:10.5194/acp-15-6283-2015, 2015.
* Moch2020:     Moch et al, JGR, https;//doi.org/10.1029/2020JD032706, 2020.
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



using Catalyst
using Unitful
using OrdinaryDiffEq
using DifferentialEquations
using PyPlot

include("Lowfunctions.jl")
include("Ratesfunctions.jl")
include("RateLawUtilFuncs.jl")

export fullchem

# Add unit "ppb" to Unitful 
module MyUnits
using Unitful
@unit ppb "ppb" Number 1 / 1000000000 false
end
Unitful.register(MyUnits)




begin
    function fullchem()
    
        @parameters t [unit = u"s"]
        
        @variables A3O2(t)        =59.469  [unit = u"ppb"]   #CH3CH2CH2OO; Primary RO2 from C3H8
        @variables ACET(t)        =87.473  [unit = u"ppb"]   #CH3C(O)CH3; Acetone
        @variables ACTA(t)        =70.187  [unit = u"ppb"]   #CH3C(O)OH; Acetic acid
        @variables AERI(t)        =30.013  [unit = u"ppb"]   #I; Dissolved iodine
        @variables ALD2(t)        =9.4976  [unit = u"ppb"]   #CH3CHO; Acetaldehyde
        @variables ALK4(t)        =8.4827  [unit = u"ppb"]   #>= C4 alkanes
        @variables AONITA(t)      =77.604  [unit = u"ppb"]   #Aerosol-phase organic nitrate from aromatic precursors
        @variables AROMRO2(t)     =77.132  [unit = u"ppb"]   #generic peroxy radical from aromatic oxidation
        @variables AROMP4(t)      =52.322  [unit = u"ppb"]   #Generic C4 product from aromatic oxidation
        @variables AROMP5(t)      =4.8975  [unit = u"ppb"]   #Generic C5 product from aromatic oxidation
        @variables ATO2(t)        =88.056  [unit = u"ppb"]   #CH3C(O)CH2O2; RO2 from acetone
        @variables ATOOH(t)       =95.672  [unit = u"ppb"]   #CH3C(O)CH2OOH; ATO2 peroxide
        @variables B3O2(t)        =44.692  [unit = u"ppb"]   #CH3CH(OO)CH3; Secondary RO2 from C3H8
        @variables BALD(t)        =35.859  [unit = u"ppb"]   #benzaldehyde and tolualdehyde
        @variables BENZ(t)        =18.166  [unit = u"ppb"]   #C6H6; Benzene
        @variables BENZO(t)       =92.094  [unit = u"ppb"]   #C6H5O radical
        @variables BENZO2(t)      =45.770  [unit = u"ppb"]   #C6H5O2 radical
        @variables BENZP(t)       =45.646  [unit = u"ppb"]   #hydroperoxide from BENZO2
        @variables Br(t)          =56.283  [unit = u"ppb"]   #Br; Atomic bromine
        @variables Br2(t)         =98.265  [unit = u"ppb"]   #Br2; Molecular bromine
        @variables BrCl(t)        =70.858  [unit = u"ppb"]   #BrCl; Bromine chloride
        @variables BrNO2(t)       =15.765  [unit = u"ppb"]   #BrNO2; Nitryl bromide
        @variables BrNO3(t)       =8.8808  [unit = u"ppb"]   #BrNO3; Bromine nitrate
        @variables BrO(t)         =90.871  [unit = u"ppb"]   #BrO; Bromine monoxide
        @variables BRO2(t)        =90.134  [unit = u"ppb"]   #C6H5O2 ; Peroxy radical from BENZ oxidation
        @variables BrSALA(t)      =5.6082  [unit = u"ppb"]   #Br; Fine sea salt bromine
        @variables BrSALC(t)      =23.602  [unit = u"ppb"]   #Br; Course sea salt bromine
        @variables BZCO3(t)       =34.301  [unit = u"ppb"]   #benzoylperoxy radical
        @variables BZCO3H(t)      =18.280  [unit = u"ppb"]   #perbenzoic acid
        @variables BZPAN(t)       =25.224  [unit = u"ppb"]   #peroxybenzoyl nitrate
        @variables C2H2(t)        =64.023  [unit = u"ppb"]   #C2H2; Ethyne
        @variables C2H4(t)        =48.832  [unit = u"ppb"]   #Ethylene
        @variables C2H6(t)        =95.500  [unit = u"ppb"]   #C2H6; Ethane
        @variables C3H8(t)        =93.457  [unit = u"ppb"]   #C3H8; Propane
        @variables C4HVP1(t)      =90.362  [unit = u"ppb"]   #C4 hydroxy-vinyl-peroxy radicals from HPALDs
        @variables C4HVP2(t)      =60.972  [unit = u"ppb"]   #C4 hydroxy-vinyl-peroxy radicals from HPALDs
        @variables CCl4(t)        =38.249  [unit = u"ppb"]   #CCl4; Carbon tetrachloride
        @variables CFC11(t)       =16.012  [unit = u"ppb"]   #CCl3F ; CFC-11, R-11, Freon 11
        @variables CFC12(t)       =60.683  [unit = u"ppb"]   #CCl2F2; CFC-12, R-12, Freon 12
        @variables CFC113(t)      =99.375  [unit = u"ppb"]   #C2Cl3F3; CFC-113, Freon 113
        @variables CFC114(t)      =82.883  [unit = u"ppb"]   #C2Cl2F4; CFC-114, Freon 114
        @variables CFC115(t)      =63.617  [unit = u"ppb"]   #C2ClF5; CFC-115, Freon 115
        @variables CH2Br2(t)      =53.658  [unit = u"ppb"]   #CH3Br2; Dibromomethane
        @variables CH2Cl2(t)      =6.5485  [unit = u"ppb"]   #CH2Cl2; Dichloromethane
        @variables CH2I2(t)       =62.463  [unit = u"ppb"]   #CH2I2; Diiodomethane
        @variables CH2IBr(t)      =98.965  [unit = u"ppb"]   #CH2IBr; Bromoiodomethane
        @variables CH2ICl(t)      =33.656  [unit = u"ppb"]   #CH2ICl; Chloroiodomethane
        @variables CH2O(t)        =34.136  [unit = u"ppb"]   #CH2O; Formaldehyde
        @variables CH2OO(t)       =72.666  [unit = u"ppb"]   #CH2OO; Criegee intermediate
        @variables CH3Br(t)       =76.768  [unit = u"ppb"]   #CH3Br; Methyl bromide
        @variables CH3CCl3(t)     =99.842  [unit = u"ppb"]   #CH3CCl3; Methyl chloroform
        @variables CH3CHOO(t)     =27.830  [unit = u"ppb"]   #CH3CHOO; Criegee intermediate
        @variables CH3Cl(t)       =14.055  [unit = u"ppb"]   #CH3Cl; Chloromethane
        @variables CH3I(t)        =19.671  [unit = u"ppb"]   #CH3I; Methyl iodide
        @variables CH4(t)         =6.9330  [unit = u"ppb"]   #CH4; Methane
        @variables CHBr3(t)       =84.204  [unit = u"ppb"]   #CHBr3; Tribromethane
        @variables CHCl3(t)       =79.435  [unit = u"ppb"]   #CHCl3; Chloroform
        @variables Cl(t)          =47.484  [unit = u"ppb"]   #Cl; Atomic chlorine
        @variables Cl2(t)         =71.292  [unit = u"ppb"]   #Cl2; Molecular chlorine
        @variables Cl2O2(t)       =96.248  [unit = u"ppb"]   #Cl2O2; Dichlorine dioxide
        @variables ClNO2(t)       =36.482  [unit = u"ppb"]   #ClNO2; Nitryl chloride
        @variables ClNO3(t)       =55.020  [unit = u"ppb"]   #ClONO2; Chlorine nitrate
        @variables ClO(t)         =0.9863  [unit = u"ppb"]   #ClO; Chlorine monoxide
        @variables ClOO(t)        =53.802  [unit = u"ppb"]   #ClOO; Chlorine dioxide
        @variables CO(t)          =32.978  [unit = u"ppb"]   #CO; Carbon monoxide
        @variables CO2(t)         =55.611  [unit = u"ppb"]   #CO2; Carbon dioxide
        @variables CSL(t)         =23.332  [unit = u"ppb"]   #cresols and xylols
        @variables DMS(t)         =5.8615  [unit = u"ppb"]   #(CH3)2S; Dimethylsulfide
        @variables EOH(t)         =18.988  [unit = u"ppb"]   #C2H5OH; Ethanol
        @variables ETHLN(t)       =69.851  [unit = u"ppb"]   #CHOCH2ONO2; Ethanal nitrate
        @variables ETHN(t)        =24.767  [unit = u"ppb"]   #stable hydroxy-nitrooxy-ethane
        @variables ETHP(t)        =59.843  [unit = u"ppb"]   #stable hydroxy-hydroperoxy-ethane
        @variables ETNO3(t)       =11.558  [unit = u"ppb"]   #C2H5ONO2; Ethyl nitrate
        @variables ETO(t)         =85.116  [unit = u"ppb"]   #hydroxy-alkoxy-ethane radical
        @variables ETOO(t)        =24.842  [unit = u"ppb"]   #hydroxy-peroxy-ethane radical, formed from ethene + OH
        @variables ETO2(t)        =29.300  [unit = u"ppb"]   #CH3CH2OO; Ethylperoxy radical
        @variables ETP(t)         =81.950  [unit = u"ppb"]   #CH3CH2OOH; Ethylhydroperoxide
        @variables GLYC(t)        =20.407  [unit = u"ppb"]   #HOCH2CHO; Glycoaldehyde
        @variables GLYX(t)        =81.571  [unit = u"ppb"]   #CHOCHO; Glyoxal
        @variables H(t)           =7.6410  [unit = u"ppb"]   #H; Atomic hydrogen
        @variables H1211(t)       =46.027  [unit = u"ppb"]   #CBrClF2; H-1211
        @variables H1301(t)       =36.584  [unit = u"ppb"]   #CBrF3; H-1301
        @variables H2402(t)       =52.639  [unit = u"ppb"]   #C2Br2F4; H-2402
        @variables H2O(t)         =56.623  [unit = u"ppb"]   #H2O; Water vapor
        @variables H2O2(t)        =44.325  [unit = u"ppb"]   #H2O2; Hydrogen peroxide
        @variables HAC(t)         =88.335  [unit = u"ppb"]   #HOCH2C(O)CH3; Hydroxyacetone
        @variables HBr(t)         =15.715  [unit = u"ppb"]   #HBr; Hypobromic acid
        @variables HC5A(t)        =44.638  [unit = u"ppb"]   #C5H8O2; Isoprene-4,1-hydroxyaldehyde
        @variables HCFC123(t)     =22.296  [unit = u"ppb"]   #C2HCl2F3; HCFC-123, R-123, Freon 123
        @variables HCFC141b(t)    =93.091  [unit = u"ppb"]   #C(CH3)Cl2F; HCFC-141b, R-141b, Freon 141b
        @variables HCFC142b(t)    =81.415  [unit = u"ppb"]   #C(CH3)ClF2; HCFC-142b, R-142b, Freon 142b
        @variables HCFC22(t)      =16.333  [unit = u"ppb"]   #CHClF2 ; HCFC-22, R-22, Freon 22
        @variables HCl(t)         =87.644  [unit = u"ppb"]   #HCl; Hydrochloric acid
        @variables HCOOH(t)       =99.489  [unit = u"ppb"]   #HCOOH; Formic acid
        @variables HI(t)          =68.377  [unit = u"ppb"]   #HI; Hydrogen iodide
        @variables HMHP(t)        =16.437  [unit = u"ppb"]   #HOCH2OOH; Hydroxymethyl hydroperoxide
        @variables HMML(t)        =33.654  [unit = u"ppb"]   #C4H6O3; Hydroxymethyl-methyl-a-lactone
        @variables HMS(t)         =54.099  [unit = u"ppb"]   #HOCH2SO3-; hydroxymethanesulfonate
        @variables HNO2(t)        =18.174  [unit = u"ppb"]   #HONO; Nitrous acid
        @variables HNO3(t)        =62.170  [unit = u"ppb"]   #HNO3; Nitric acid
        @variables HNO4(t)        =39.302  [unit = u"ppb"]   #HNO4; Pernitric acid
        @variables HO2(t)         =58.306  [unit = u"ppb"]   #HO2; Hydroperoxyl radical
        @variables HOBr(t)        =59.529  [unit = u"ppb"]   #HOBr; Hypobromous acid
        @variables HOCl(t)        =90.397  [unit = u"ppb"]   #HOCl; Hypochlorous acid
        @variables HOI(t)         =35.827  [unit = u"ppb"]   #HOI; Hypoiodous acid
        @variables HONIT(t)       =58.760  [unit = u"ppb"]   #2nd gen monoterpene organic nitrate
        @variables HPALD1(t)      =73.164  [unit = u"ppb"]   #O=CHC(CH3)=CHCH2OOH; d-4,1-C5-hydroperoxyaldehyde
        @variables HPALD1OO(t)    =34.759  [unit = u"ppb"]   #peroxy radicals from HPALD1
        @variables HPALD2(t)      =1.0291  [unit = u"ppb"]   #HOOCH2C(CH3)=CHCH=O; d-1,4-C5-hydroperoxyaldehyde
        @variables HPALD2OO(t)    =28.332  [unit = u"ppb"]   #peroxy radicals from HPALD2
        @variables HPALD3(t)      =91.404  [unit = u"ppb"]   #O=CHC(CH3)OOHCH=CH2; b-2,1-C5-hydroperoxyaldehyde
        @variables HPALD4(t)      =14.949  [unit = u"ppb"]   #CH2=C(CH3)CHOOHCH=O; b-3,4-C5-hydroperoxyaldehyde
        @variables HPETHNL(t)     =90.330  [unit = u"ppb"]   #CHOCH2OOH; hydroperoxyethanal
        @variables I(t)           =85.092  [unit = u"ppb"]   #I; Atmoic iodine
        @variables I2(t)          =49.639  [unit = u"ppb"]   #I2; Molecular iodine
        @variables I2O2(t)        =70.687  [unit = u"ppb"]   #I2O2; Diiodine dioxide
        @variables I2O3(t)        =85.950  [unit = u"ppb"]   #I2O3; Diiodine trioxide
        @variables I2O4(t)        =36.304  [unit = u"ppb"]   #I2O4; Diiodine tetraoxide
        @variables IBr(t)         =46.529  [unit = u"ppb"]   #IBr; Iodine monobromide
        @variables ICHE(t)        =17.739  [unit = u"ppb"]   #C5H8O3; Isoprene hydroxy-carbonyl-epoxides
        @variables ICHOO(t)       =75.448  [unit = u"ppb"]   #peroxy radical from IEPOXD
        @variables ICl(t)         =32.124  [unit = u"ppb"]   #ICl; Iodine monochloride
        @variables ICN(t)         =22.128  [unit = u"ppb"]   #C5H7NO4; Lumped isoprene carbonyl nitrates
        @variables ICNOO(t)       =12.649  [unit = u"ppb"]   #peroxy radicals from ICN
        @variables ICPDH(t)       =71.736  [unit = u"ppb"]   #C5H10O5; Isoprene dihydroxy hydroperoxycarbonyl
        @variables IDC(t)         =28.840  [unit = u"ppb"]   #C5H6O2; Lumped isoprene dicarbonyls
        @variables IDCHP(t)       =2.3035  [unit = u"ppb"]   #C5H8O5; Isoprene dicarbonyl hydroxy dihydroperoxide
        @variables IDHDP(t)       =21.649  [unit = u"ppb"]   #C5H12O6; Isoprene dihydroxy dihydroperoxide
        @variables IDHNBOO(t)     =25.086  [unit = u"ppb"]   #peroxy radicals from INPB
        @variables IDHNDOO1(t)    =15.276  [unit = u"ppb"]   #peroxy radicals from INPD
        @variables IDHNDOO2(t)    =61.961  [unit = u"ppb"]   #peroxy radicals from INPD
        @variables IDHPE(t)       =94.216  [unit = u"ppb"]   #C5H10O5; Isoprene dihydroxy hydroperoxy epoxide
        @variables IDN(t)         =13.328  [unit = u"ppb"]   #C5H8N2O6; Lumped isoprene dinitrates
        @variables IDNOO(t)       =44.188  [unit = u"ppb"]   #peroxy radicals from IDN
        @variables IEPOXA(t)      =60.554  [unit = u"ppb"]   #C5H10O3; trans-Beta isoprene epoxydiol
        @variables IEPOXAOO(t)    =27.948  [unit = u"ppb"]   #peroxy radical from trans-Beta isoprene epoxydiol
        @variables IEPOXB(t)      =38.609  [unit = u"ppb"]   #C5H10O3; cis-Beta isoprene epoxydiol
        @variables IEPOXBOO(t)    =87.094  [unit = u"ppb"]   #peroxy radical from cis-Beta isoprene epoxydiol
        @variables IEPOXD(t)      =8.0346  [unit = u"ppb"]   #C5H10O3; Delta isoprene epoxydiol
        @variables IHN1(t)        =11.028  [unit = u"ppb"]   #C5H9NO4; Isoprene-d-4-hydroxy-1-nitrate
        @variables IHN2(t)        =71.286  [unit = u"ppb"]   #C5H9NO4; Isoprene-b-1-hydroxy-2-nitrate
        @variables IHN3(t)        =82.443  [unit = u"ppb"]   #C5H9NO4; Isoprene-b-4-hydroxy-3-nitrate
        @variables IHN4(t)        =91.975  [unit = u"ppb"]   #C5H9NO4; Isoprene-d-1-hydroxy-4-nitrate
        @variables IHOO1(t)       =27.517  [unit = u"ppb"]   #peroxy radical from OH addition to isoprene at C1
        @variables IHOO4(t)       =51.624  [unit = u"ppb"]   #peroxy radical from OH addition to isoprene at C4
        @variables IHPNBOO(t)     =66.363  [unit = u"ppb"]   #peroxy radicals from INPB
        @variables IHPNDOO(t)     =63.447  [unit = u"ppb"]   #peroxy radicals from INPD
        @variables IHPOO1(t)      =79.476  [unit = u"ppb"]   #peroxy radical from ISOPOOH
        @variables IHPOO2(t)      =58.328  [unit = u"ppb"]   #peroxy radical from ISOPOOH
        @variables IHPOO3(t)      =5.9515  [unit = u"ppb"]   #peroxy radical from ISOPOOH
        @variables INA(t)         =26.490  [unit = u"ppb"]   #alkoxy radical from INO2D
        @variables INDIOL(t)      =52.424  [unit = u"ppb"]   #Generic aerosol phase organonitrate hydrolysis product
        @variables INO(t)         =73.106  [unit = u"ppb"]   #INO; Nitrosyl iodide
        @variables INO2B(t)       =25.723  [unit = u"ppb"]   #beta-peroxy radicals from isoprene + NO3
        @variables INO2D(t)       =30.524  [unit = u"ppb"]   #delta-peroxy radicals from isoprene + NO3
        @variables INPB(t)        =96.224  [unit = u"ppb"]   #C5H9NO5; Lumped isoprene beta-hydroperoxy nitrates
        @variables INPD(t)        =80.193  [unit = u"ppb"]   #C5H9NO5; Lumped isoprene delta-hydroperoxy nitrates
        @variables IO(t)          =26.615  [unit = u"ppb"]   #IO; Iodine monoxide
        @variables IONITA(t)      =25.539  [unit = u"ppb"]   #Aerosol-phase organic nitrate from isoprene precursors
        @variables IONO(t)        =56.756  [unit = u"ppb"]   #IONO; Nitryl iodide
        @variables IONO2(t)       =53.992  [unit = u"ppb"]   #IONO2; Iodine nitrate
        @variables IPRNO3(t)      =22.638  [unit = u"ppb"]   #C3H8ONO2; Isopropyl nitrate
        @variables ISALA(t)       =79.006  [unit = u"ppb"]   #I; Fine sea salt iodine
        @variables ISALC(t)       =13.184  [unit = u"ppb"]   #I; Coarse sea salt iodine
        @variables ISOP(t)        =83.240  [unit = u"ppb"]   #CH2=C(CH3)CH=CH2; Isoprene
        @variables ISOPNOO1(t)    =60.043  [unit = u"ppb"]   #peroxy radicals from IHN2
        @variables ISOPNOO2(t)    =71.160  [unit = u"ppb"]   #peroxy radicals from IHN3
        @variables ITCN(t)        =26.206  [unit = u"ppb"]   #C5H9NO7; Lumped tetrafunctional isoprene carbonyl-nitrates
        @variables ITHN(t)        =11.572  [unit = u"ppb"]   #C5H11NO7; Lumped tetrafunctional isoprene hydroxynitrates
        @variables KO2(t)         =84.868  [unit = u"ppb"]   #RO2 from >3 ketones
        @variables LBRO2H(t)      =61.633  [unit = u"ppb"]   #Dummy spc to track oxidation of BRO2 by HO2
        @variables LBRO2N(t)      =11.255  [unit = u"ppb"]   #Dummy spc to track oxidation of BRO2 by NO
        @variables LIMO(t)        =53.446  [unit = u"ppb"]   #C10H16; Limonene
        @variables LIMO2(t)       =92.794  [unit = u"ppb"]   #RO2 from LIMO
        @variables LISOPOH(t)     =92.199  [unit = u"ppb"]   #Dummy spc to track oxidation of ISOP by OH
        @variables LISOPNO3(t)    =29.918  [unit = u"ppb"]   #Dummy spc to track oxidation of ISOP by NO3
        @variables LNRO2H(t)      =29.939  [unit = u"ppb"]   #Dummy spc to track oxidation of NRO2 by HO2
        @variables LNRO2N(t)      =57.661  [unit = u"ppb"]   #Dummy spc to track oxidation of NRO2 by NO
        @variables LTRO2H(t)      =46.441  [unit = u"ppb"]   #Dummy spc to track oxidation of TRO2 by HO2
        @variables LTRO2N(t)      =28.826  [unit = u"ppb"]   #Dummy spc to track oxidation of TRO2 by NO
        @variables LVOC(t)        =47.079  [unit = u"ppb"]   #C5H14O5; Gas-phase low-volatility non-IEPOX product of ISOPOOH (RIP) oxidation
        @variables LVOCOA(t)      =97.788  [unit = u"ppb"]   #C5H14O5; Aerosol-phase low-volatility non-IEPOX product of ISOPOOH (RIP) oxidation
        @variables LXRO2H(t)      =62.462  [unit = u"ppb"]   #Dummy spc to track oxidation of XRO2 by HO2
        @variables LXRO2N(t)      =39.337  [unit = u"ppb"]   #Dummy spc to track oxidation of XRO2 by NO
        @variables MACR(t)        =17.935  [unit = u"ppb"]   #CH2=C(CH3)CHO; Methacrolein
        @variables MACR1OO(t)     =91.114  [unit = u"ppb"]   #peroxyacyl radical from MACR + OH
        @variables MACR1OOH(t)    =61.547  [unit = u"ppb"]   #CH2=C(CH3)C(O)OOH; Peracid from MACR
        @variables MACRNO2(t)     =97.585  [unit = u"ppb"]   #Product of MCRHN + OH
        @variables MAP(t)         =61.421  [unit = u"ppb"]   #CH3C(O)OOH; Peroxyacetic acid
        @variables MCO3(t)        =0.6220  [unit = u"ppb"]   #CH3C(O)OO; Peroxyacetyl radical
        @variables MCRDH(t)       =23.756  [unit = u"ppb"]   #C4H8O3; Dihydroxy-MACR
        @variables MCRENOL(t)     =43.151  [unit = u"ppb"]   #C4H6O2; Lumped enols from MVK/MACR
        @variables MCRHN(t)       =5.5067  [unit = u"ppb"]   #HOCH2C(ONO2)(CH3)CHO; Hydroxynitrate from MACR
        @variables MCRHNB(t)      =2.8277  [unit = u"ppb"]   #O2NOCH2C(OH)(CH3)CHO; Hydroxynitrate from MACR
        @variables MCRHP(t)       =0.2048  [unit = u"ppb"]   #HOCH2C(OOH)(CH3)CHO; Hydroxy-hydroperoxy-MACR
        @variables MCROHOO(t)     =49.422  [unit = u"ppb"]   #peroxy radical from MACR + OH
        @variables MCT(t)         =77.753  [unit = u"ppb"]   #methylcatechols
        @variables MEK(t)         =61.732  [unit = u"ppb"]   #RC(O)R; Methyl ethyl ketone
        @variables MENO3(t)       =89.773  [unit = u"ppb"]   #CH3ONO2; methyl nitrate
        @variables MGLY(t)        =93.209  [unit = u"ppb"]   #CH3COCHO; Methylglyoxal
        @variables MO2(t)         =73.881  [unit = u"ppb"]   #CH3O2; Methylperoxy radical
        @variables MOH(t)         =55.585  [unit = u"ppb"]   #CH3OH; Methanol
        @variables MONITA(t)      =61.052  [unit = u"ppb"]   #Aerosol-phase organic nitrate from monoterpene precursors
        @variables MONITS(t)      =71.179  [unit = u"ppb"]   #Saturated 1st gen monoterpene organic nitrate
        @variables MONITU(t)      =23.294  [unit = u"ppb"]   #Unsaturated 1st gen monoterpene organic nitrate
        @variables MP(t)          =25.386  [unit = u"ppb"]   #CH3OOH; Methylhydroperoxide
        @variables MPAN(t)        =89.575  [unit = u"ppb"]   #CH2=C(CH3)C(O)OONO2; Peroxymethacroyl nitrate (PMN)
        @variables MPN(t)         =59.399  [unit = u"ppb"]   #CH3O2NO2; Methyl peroxy nitrate
        @variables MSA(t)         =82.711  [unit = u"ppb"]   #CH4SO3; Methanesulfonic acid
        @variables MTPA(t)        =25.413  [unit = u"ppb"]   #Lumped monoterpenes: a-pinene, b-pinene, sabinene, carene
        @variables MTPO(t)        =8.9141  [unit = u"ppb"]   #Other monoterpenes: Terpinene, terpinolene, myrcene, ocimene, other monoterpenes
        @variables MVK(t)         =62.109  [unit = u"ppb"]   #CH2=CHC(=O)CH3; Methyl vinyl ketone
        @variables MVKDH(t)       =97.635  [unit = u"ppb"]   #HOCH2CH2OHC(O)CH3; Dihydroxy-MVK
        @variables MVKHC(t)       =98.159  [unit = u"ppb"]   #C4H6O3; MVK hydroxy-carbonyl
        @variables MVKHCB(t)      =50.621  [unit = u"ppb"]   #C4H6O3; MVK hydroxy-carbonyl
        @variables MVKHP(t)       =35.334  [unit = u"ppb"]   #C4H8O4; MVK hydroxy-hydroperoxide
        @variables MVKN(t)        =48.570  [unit = u"ppb"]   #HOCH2CH(ONO2)C(=O)CH3; Hydroxynitrate from MVK
        @variables MVKOHOO(t)     =41.490  [unit = u"ppb"]   #peroxy radical from MVK + OH
        @variables MVKPC(t)       =8.2590  [unit = u"ppb"]   #OCHCH(OOH)C(O)CH3; MVK hydroperoxy-carbonyl
        @variables N(t)           =56.891  [unit = u"ppb"]   #N; Atomic nitrogen
        @variables N2O(t)         =50.751  [unit = u"ppb"]   #N2O; Nitrous oxide
        @variables N2O5(t)        =19.118  [unit = u"ppb"]   #N2O5; Dinitrogen pentoxide
        @variables NAP(t)         =33.270  [unit = u"ppb"]   #C10H8; Naphthalene; IVOC surrogate
        @variables NIT(t)         =77.088  [unit = u"ppb"]   #NIT; Fine mode inorganic nitrate
        @variables NITs(t)        =56.341  [unit = u"ppb"]   #NITs; Coarse mode inorganic nitrate
        @variables NO(t)          =54.454  [unit = u"ppb"]   #NO; Nitric oxide
        @variables NO2(t)         =5.5605  [unit = u"ppb"]   #NO2; Nitrogen dioxide
        @variables NO3(t)         =31.739  [unit = u"ppb"]   #NO3; Nitrate radical
        @variables NPHEN(t)       =87.002  [unit = u"ppb"]   #nitrophenols
        @variables NPRNO3(t)      =58.488  [unit = u"ppb"]   #C3H8ONO2; n-propyl nitrate
        @variables NRO2(t)        =83.618  [unit = u"ppb"]   #Peroxy radical from NAP oxidation
        @variables O(t)           =14.164  [unit = u"ppb"]   #O(3P); Ground state atomic oxygen
        @variables O1D(t)         =17.648  [unit = u"ppb"]   #O(1D); Excited atomic oxygen
        @variables O3(t)          =2.1326  [unit = u"ppb"]   #O3; Ozone
        @variables O3A(t)         =28.223  [unit = u"ppb"]   #O3; Ozone in accum seasalt
        @variables O3C(t)         =49.728  [unit = u"ppb"]   #O3; Ozone in coarse seasalt
        @variables OClO(t)        =60.037  [unit = u"ppb"]   #OClO; Chlorine dioxide
        @variables OCS(t)         =20.867  [unit = u"ppb"]   #COS; Carbonyl sulfide
        @variables OH(t)          =19.293  [unit = u"ppb"]   #OH; Hydroxyl radical
        @variables OIO(t)         =51.406  [unit = u"ppb"]   #OIO; Iodine dioxide
        @variables OLND(t)        =63.247  [unit = u"ppb"]   #Monoterpene-derived NO3-alkene adduct
        @variables OLNN(t)        =95.667  [unit = u"ppb"]   #Monoterpene-derived NO3 adduct
        @variables OTHRO2(t)      =74.722  [unit = u"ppb"]   #Other C2 RO2 not from C2H6 oxidation
        @variables PAN(t)         =59.945  [unit = u"ppb"]   #CH3C(O)OONO2; Peroxyacetylnitrate
        @variables PHEN(t)        =59.463  [unit = u"ppb"]   #phenol
        @variables PIO2(t)        =94.953  [unit = u"ppb"]   #RO2 from MTPA
        @variables PIP(t)         =96.809  [unit = u"ppb"]   #Peroxides from MTPA
        @variables PO2(t)         =25.827  [unit = u"ppb"]   #HOCH2CH(OO)CH3; RO2 from propene
        @variables PP(t)          =61.804  [unit = u"ppb"]   #HOCH2CH(OOH)CH3; Peroxide from PO2
        @variables PPN(t)         =94.436  [unit = u"ppb"]   #CH3CH2C(O)OONO2; Peroxypropionylnitrate
        @variables PRN1(t)        =87.618  [unit = u"ppb"]   #O2NOCH2CH(OO)CH3; RO2 from propene + NO3
        @variables PROPNN(t)      =33.842  [unit = u"ppb"]   #CH3C(=O)CH2ONO2; Propanone nitrate
        @variables PRPE(t)        =63.131  [unit = u"ppb"]   #C3H6; >= C3 alkenes
        @variables PRPN(t)        =43.083  [unit = u"ppb"]   #O2NOCH2CH(OOH)CH3; Peroxide from PRN1
        @variables PYAC(t)        =79.104  [unit = u"ppb"]   #CH3COCOOH; Pyruvic acid
        @variables R4N1(t)        =37.438  [unit = u"ppb"]   #RO2 from R4N2
        @variables R4N2(t)        =46.153  [unit = u"ppb"]   #RO2NO; >= C4 alkylnitrates
        @variables R4O2(t)        =36.162  [unit = u"ppb"]   #RO2 from ALK4
        @variables R4P(t)         =40.804  [unit = u"ppb"]   #CH3CH2CH2CH2OOH; Peroxide from R4O2
        @variables RA3P(t)        =57.025  [unit = u"ppb"]   #CH3CH2CH2OOH; Peroxide from A3O2
        @variables RB3P(t)        =72.421  [unit = u"ppb"]   #CH3CH(OOH)CH3; Peroxide from B3O2
        @variables RCHO(t)        =37.680  [unit = u"ppb"]   #CH3CH2CHO; >= C3 aldehydes
        @variables RCO3(t)        =63.108  [unit = u"ppb"]   #CH3CH2C(O)OO; Peroxypropionyl radical
        @variables RIPA(t)        =42.881  [unit = u"ppb"]   #HOCH2C(OOH)(CH3)CH=CH2; 1,2-ISOPOOH
        @variables RIPB(t)        =24.737  [unit = u"ppb"]   #HOCH2C(OOH)(CH3)CH=CH2; 4,3-ISOPOOH
        @variables RIPC(t)        =57.374  [unit = u"ppb"]   #C5H10O3; d(1,4)-ISOPOOH
        @variables RIPD(t)        =5.9941  [unit = u"ppb"]   #C5H10O3; d(4,1)-ISOPOOH
        @variables ROH(t)         =67.750  [unit = u"ppb"]   #C3H7OH; > C2 alcohols
        @variables RP(t)          =38.598  [unit = u"ppb"]   #CH3CH2C(O)OOH; Peroxide from RCO3
        @variables SALAAL(t)      =37.928  [unit = u"ppb"]   #Accumulation mode seasalt aerosol alkalinity
        @variables SALCAL(t)      =13.651  [unit = u"ppb"]   #Coarse mode seasalt aerosol alkalinity
        @variables SALACL(t)      =41.613  [unit = u"ppb"]   #Cl; Fine chloride
        @variables SALCCL(t)      =41.813  [unit = u"ppb"]   #Cl; Coarse chloride
        @variables SALASO2(t)     =46.971  [unit = u"ppb"]   #SO2; Fine seasalt
        @variables SALCSO2(t)     =98.028  [unit = u"ppb"]   #SO2; Coarse seasalt
        @variables SALASO3(t)     =18.107  [unit = u"ppb"]   #SO3--; Fine seasalt
        @variables SALCSO3(t)     =51.008  [unit = u"ppb"]   #SO3--; Coarse chloride
        @variables SO2(t)         =6.2180  [unit = u"ppb"]   #SO2; Sulfur dioxide
        @variables SO4(t)         =97.494  [unit = u"ppb"]   #SO4; Sulfate
        @variables SO4s(t)        =74.244  [unit = u"ppb"]   #SO4 on sea-salt; Sulfate
        @variables SOAGX(t)       =82.469  [unit = u"ppb"]   #CHOCHO; Aerosol-phase glyoxal
        @variables SOAIE(t)       =81.289  [unit = u"ppb"]   #C5H10O3; Aerosol-phase IEPOX
        @variables TOLU(t)        =58.265  [unit = u"ppb"]   #C7H8; Toluene
        @variables TRO2(t)        =10.659  [unit = u"ppb"]   #Peroxy radical from TOLU oxidation
        @variables XYLE(t)        =24.563  [unit = u"ppb"]   #C8H10; Xylene
        @variables XRO2(t)        =25.058  [unit = u"ppb"]   #Peroxy radical from XYLE oxidation
        @variables H2(t)          =97.924  [unit = u"ppb"]   #H2; Molecular hydrogen
        @variables N2(t)          =36.998  [unit = u"ppb"]   #N2; Molecular nitrogen
        @variables O2(t)          =51.201  [unit = u"ppb"]   #O2; Molecular oxygen
        @variables RCOOH(t)       =12.892  [unit = u"ppb"]   #C2H5C(O)OH; > C2 organic acids 
        
        

        
        fc = [

            
            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            # %%%%% Gas-phase chemistry reactions                                   %%%%%
            # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
            # NOTES:
    
            # To avoid useless CPU cycles, we have introduced new rate law functions
            # that skip computing Arrhenius terms (and other terms) that would
            # evaluate to 1.  The Arrhenius terms that are passed to the function
            # are in most cases now noted in the function name (e.g. GCARR_abc takes
            # Arrhenius A, B, C parameters but GCARR_ac only passes A and C
            # parameters because B=0 and the (300/T)*B would evaluate to 1).
            # This should be much more computationally efficient, as these functions
            # are called (sometimes multiple times) for each grid box where we
            # perform chemistry.
            # -- Bob Yantosca (25 Jan 2020)
            

            Reaction(GCARR_ac(3.00e-12, -1500.0e0), [O3, NO], [NO2, O2], [1, 1], [1, 1])  #  O3 + NO = NO2 + O2   
            Reaction(GCARR_ac(1.70e-12, -940.0e0), [O3, OH], [HO2, O2], [1, 1], [1, 1])  #  O3 + OH = HO2 + O2   
            Reaction(GCARR_ac(1.00e-14, -490.0e0), [O3, HO2], [OH, O2, O2], [1, 1], [1, 1, 1])  #  O3 + HO2 = OH + O2 + O2   
            Reaction(GCARR_ac(1.20e-13, -2450.0e0), [O3, NO2], [O2, NO3], [1, 1], [1, 1])  #  O3 + NO2 = O2 + NO3   
            Reaction(GCARR_ac(2.90e-16, -1000.0e0), [O3, MO2], [CH2O, HO2, O2], [1, 1], [1, 1, 1])   #2014/02/03; Eastham2014; SDE  #  O3 + MO2 = CH2O + HO2 + O2   
            Reaction(1.80e-12, [OH, OH], [H2O, O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + OH = H2O + O   
            Reaction(GCJPLPR_aba(6.90e-31, 1.0e+00, 2.6e-11, 0.6e0), [OH, OH], [H2O2], [1, 1], [1])  #  OH + OH = H2O2   
            Reaction(GCARR_ac(4.80e-11, 250.0e0), [OH, HO2], [H2O, O2], [1, 1], [1, 1])  #  OH + HO2 = H2O + O2   
            Reaction(1.80e-12, [OH, H2O2], [H2O, HO2], [1, 1], [1, 1])  #  OH + H2O2 = H2O + HO2   
            Reaction(GCARR_ac(3.30e-12, 270.0e0), [HO2, NO], [OH, NO2], [1, 1], [1, 1])   #2013/02/12; JPL 10-6; BHH,JMAO,EAM  #  HO2 + NO = OH + NO2   
            Reaction(GC_HO2HO2_acac(3.00e-13, 460.0e0, 2.1e-33, 920.0e0), [HO2, HO2], [H2O2, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  HO2 + HO2 = H2O2 + O2   
            Reaction(GC_OHCO_a(1.50e-13), [OH, CO], [HO2, CO2], [1, 1], [1, 1])   #2017/02/22; JPL 15-10; BHH,MJE  #  OH + CO = HO2 + CO2   
            Reaction(GCARR_ac(2.45e-12, -1775.0e0), [OH, CH4], [MO2, H2O], [1, 1], [1, 1])  #  OH + CH4 = MO2 + H2O   
            Reaction(GC_RO2NO_B1_ac(2.80e-12, 300.0e0), [MO2, NO], [CH2O, HO2, NO2], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF  #  MO2 + NO = CH2O + HO2 + NO2   
            Reaction(GC_RO2NO_A1_ac(2.80e-12, 300.0e0), [MO2, NO], [MENO3], [1, 1], [1])   #2019/05/10; Fisher2018; JAF  #  MO2 + NO = MENO3   
            Reaction(GCARR_abc(4.10e-13, 0.0e0, 750.0e0), [MO2, HO2], [MP, O2], [1, 1], [1, 1])  #  MO2 + HO2 = MP + O2   
            Reaction(GC_TBRANCH_1_acac(9.50e-14, 390.0e0, 2.62e1, -1130.0e0), [MO2, MO2], [MOH, CH2O, O2], [1, 1], [1, 1, 1])  #  MO2 + MO2 = MOH + CH2O + O2   
            Reaction(GC_TBRANCH_1_acac(9.50e-14, 390.0e0, 4.0e-2, 1130.0e0), [MO2, MO2], [CH2O, HO2], [1, 1], [2.000, 2.000])  #  MO2 + MO2 = 2.000CH2O + 2.000HO2   
            Reaction(1.60e-10 , [MO2, OH], [MOH, CH2O, HO2], [1, 1], [0.13, 0.87, 1.74])   #2021/09/22; Bates2021a; KHB,MSL  #  MO2 + OH = 0.13MOH + 0.87CH2O + 1.74HO2   
            Reaction(GCARR_ac(2.66e-12, 200.0e0), [MP, OH], [MO2, H2O], [1, 1], [1, 1])  #  MP + OH = MO2 + H2O   
            Reaction(GCARR_ac(1.14e-12, 200.0e0), [MP, OH], [CH2O, OH, H2O], [1, 1], [1, 1, 1])  #  MP + OH = CH2O + OH + H2O   
            Reaction(GCARR_ac(2.66e-12, 200.0e0), [ATOOH, OH], [ATO2, H2O], [1, 1], [1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  ATOOH + OH = ATO2 + H2O   
            Reaction(GCARR_ac(1.14e-12, 200.0e0), [ATOOH, OH], [MGLY, OH, H2O], [1, 1], [1, 1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  ATOOH + OH = MGLY + OH + H2O   
            Reaction(GCARR_ac(5.50e-12, 125.0e0), [CH2O, OH], [CO, HO2, H2O], [1, 1], [1, 1, 1])  #  CH2O + OH = CO + HO2 + H2O   
            Reaction(GCJPLPR_aba(1.80e-30, 3.0e+00, 2.8e-11, 0.6e0), [NO2, OH], [HNO3], [1, 1], [1])  #  NO2 + OH = HNO3   
            Reaction(GC_OHHNO3_acacac(2.41e-14, 460.0e0, 2.69e-17, 2199.0e0, 6.51e-34, 1335.0e0), [HNO3, OH], [H2O, NO3], [1, 1], [1, 1])  #  HNO3 + OH = H2O + NO3   
            Reaction(GCJPLPR_abab(7.00e-31, 2.6e+00, 3.60e-11, 0.1e0, 0.6e0), [NO, OH], [HNO2], [1, 1], [1])  #  NO + OH = HNO2   
            Reaction(GCARR_ac(1.80e-11, -390.0e0), [HNO2, OH], [H2O, NO2], [1, 1], [1, 1])  #  HNO2 + OH = H2O + NO2   
            Reaction(GCJPLPR_abab(1.90e-31, 3.4e+00, 4.0e-12, 0.3e0, 0.6e0), [HO2, NO2], [HNO4], [1, 1], [1])   #2017/02/22; JPL 15-10; BHH,MJE  #  HO2 + NO2 = HNO4   
            Reaction(GCJPLPR_abcabc(9.05e-05, 3.4e0, -10900.0e0, 1.90e15, 0.3e0, -10900.0e0, 0.6e0), [HNO4], [HO2, NO2], [1], [1, 1])   #2017/02/22; JPL 15-10; BHH,MJE  #  HNO4 = HO2 + NO2   
            Reaction(GCARR_ac(1.30e-12, 380.0e0), [HNO4, OH], [H2O, NO2, O2], [1, 1], [1, 1, 1])  #  HNO4 + OH = H2O + NO2 + O2   
            Reaction(3.50e-12, [HO2, NO3], [OH, NO2, O2], [1, 1], [1, 1, 1])  #  HO2 + NO3 = OH + NO2 + O2   
            Reaction(GCARR_ac(1.50e-11, 170.0e0), [NO, NO3], [NO2], [1, 1], [2.000])  #  NO + NO3 = 2.000NO2   
            Reaction(2.20e-11, [OH, NO3], [HO2, NO2], [1, 1], [1, 1])  #  OH + NO3 = HO2 + NO2   
            Reaction(GCJPLPR_abab(2.40e-30, 3.0e+00, 1.6e-12, -0.1e0, 0.6e0), [NO2, NO3], [N2O5], [1, 1], [1])   #2017/02/22; JPL 15-10; BHH,MJE  #  NO2 + NO3 = N2O5   
            Reaction(GCJPLPR_abcabc(4.14e-04, 3.0e0, -10840.0e0, 2.76e14, -0.1e0, -10840.0e0, 0.6e0), [N2O5], [NO2, NO3], [1], [1, 1])   #2017/02/22; JPL 15-10; BHH,MJE  #  N2O5 = NO2 + NO3   
            Reaction(4.00e-13, [HCOOH, OH], [H2O, CO2, HO2], [1, 1], [1, 1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  HCOOH + OH = H2O + CO2 + HO2   
            Reaction(GCARR_ac(2.90e-12, -345.0e0), [MOH, OH], [HO2, CH2O], [1, 1], [1, 1])  #  MOH + OH = HO2 + CH2O   
            Reaction(GCARR_ac(4.50e-14, -1260.0e0), [NO2, NO3], [NO, NO2, O2], [1, 1], [1, 1, 1])  #  NO2 + NO3 = NO + NO2 + O2   
            Reaction(5.80e-16, [NO3, CH2O], [HNO3, HO2, CO], [1, 1], [1, 1, 1])  #  NO3 + CH2O = HNO3 + HO2 + CO   
            Reaction(GCARR_ac(1.40e-12, -1900.0e0), [ALD2, NO3], [HNO3, MCO3], [1, 1], [1, 1])  #  ALD2 + NO3 = HNO3 + MCO3   
            Reaction(GCJPLPR_abab(9.70e-29, 5.6e+00, 9.3e-12, 1.5e0, 0.6e0), [MCO3, NO2], [PAN], [1, 1], [1])   #JPL Eval 17  #  MCO3 + NO2 = PAN   
            Reaction(GCJPLEQ_acabab(9.30e-29, 14000.0e0, 9.7e-29, 5.6e0, 9.3e-12, 1.5e0, 0.6e0), [PAN], [MCO3, NO2], [1], [1, 1])  #  PAN = MCO3 + NO2   
            Reaction(GCARR_ac(8.10e-12, 270.0e0), [MCO3, NO], [MO2, NO2, CO2], [1, 1], [1, 1, 1])  #  MCO3 + NO = MO2 + NO2 + CO2   
            Reaction(GCARR_ac(7.66e-12, -1020.0e0), [C2H6, OH], [ETO2, H2O], [1, 1], [1, 1])   #2013/02/12; JPL 10-6; BHH,JMAO,EAM  #  C2H6 + OH = ETO2 + H2O   
            Reaction(GC_RO2NO_B2_aca(2.60e-12, 365.0e0, 2.0e0), [ETO2, NO], [ALD2, NO2, HO2], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF  #  ETO2 + NO = ALD2 + NO2 + HO2   
            Reaction(GC_RO2NO_A2_aca(2.60e-12, 365.0e0, 2.0e0), [ETO2, NO], [ETNO3], [1, 1], [1])   #2019/05/10; Fisher2018; JAF  #  ETO2 + NO = ETNO3   
            Reaction(GCARR_ac(2.60e-12, 365.0e0), [OTHRO2, NO], [ALD2, NO2, HO2], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF  #  OTHRO2 + NO = ALD2 + NO2 + HO2   
            Reaction(GC_TBRANCH_2_acabc(7.60e-12, -585.0e0, 5.87e0, 0.64e0, -816.0e0), [C3H8, OH], [B3O2], [1, 1], [1])  #  C3H8 + OH = B3O2   
            Reaction(GC_TBRANCH_2_acabc(7.60e-12, -585.0e0, 1.7e-1, -0.64e0, 816.0e0), [C3H8, OH], [A3O2], [1, 1], [1])  #  C3H8 + OH = A3O2   
            Reaction(GC_RO2NO_B2_aca(2.90e-12, 350.0e0, 3.0e0), [A3O2, NO], [NO2, HO2, RCHO], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF  #  A3O2 + NO = NO2 + HO2 + RCHO   
            Reaction(GC_RO2NO_A2_aca(2.90e-12, 350.0e0, 3.0e0), [A3O2, NO], [NPRNO3], [1, 1], [1])   #2019/05/10; Fisher2018; JAF  #  A3O2 + NO = NPRNO3   
            Reaction(GCARR_ac(2.70e-12, 350.0e0), [PO2, NO], [NO2, HO2, CH2O, ALD2], [1, 1], [1, 1, 1, 1])  #  PO2 + NO = NO2 + HO2 + CH2O + ALD2   
            Reaction(GCARR_ac(9.10e-12, -405.0e0), [ALK4, OH], [R4O2], [1, 1], [1])  #  ALK4 + OH = R4O2   
            Reaction(GC_RO2NO_A2_aca(2.70e-12, 350.0e0, 4.5e0), [R4O2, NO], [R4N2], [1, 1], [1])  #  R4O2 + NO = R4N2   
            Reaction(GCARR_ac(2.80e-12, 300.0e0), [ATO2, NO], [NO2, CH2O, MCO3], [1, 1], [1, 1, 1])   #2017/07/27; Fix C creation; SAS,BHH,MJE  #  ATO2 + NO = NO2 + CH2O + MCO3   
            Reaction(GC_RO2NO_B2_aca(2.70e-12, 360.0e0, 3.0e0), [B3O2, NO], [NO2, HO2, ACET], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF  #  B3O2 + NO = NO2 + HO2 + ACET   
            Reaction(GC_RO2NO_A2_aca(2.70e-12, 360.0e0, 3.0e0), [B3O2, NO], [IPRNO3], [1, 1], [1])   #2019/05/10; Fisher2018; JAF  #  B3O2 + NO = IPRNO3   
            Reaction(GCARR_ac(2.70e-12, 350.0e0), [PRN1, NO], [NO2, CH2O, ALD2], [1, 1], [2.000, 1, 1])  #  PRN1 + NO = 2.000NO2 + CH2O + ALD2   
            Reaction(GCARR_ac(2.80e-12, -3280.0e0), [ALK4, NO3], [HNO3, R4O2], [1, 1], [1, 1])  #  ALK4 + NO3 = HNO3 + R4O2   
            Reaction(1.60e-12, [R4N2, OH], [R4N1, H2O], [1, 1], [1, 1])  #  R4N2 + OH = R4N1 + H2O   
            Reaction(GCARR_ac(3.15e-14, 920.0e0), [ACTA, OH], [MO2, CO2, H2O], [1, 1], [1, 1, 1])   #2013/02/12; JPL 10-6; BHH,JMAO,EAM  #  ACTA + OH = MO2 + CO2 + H2O   
            Reaction(GCARR_ac(6.00e-12, 410.0e0), [OH, RCHO], [RCO3, H2O], [1, 1], [1, 1])  #  OH + RCHO = RCO3 + H2O   
            Reaction(GCJPLPR_abab(9.00e-28, 8.9e0, 7.7e-12, 0.2e0, 0.6e0), [RCO3, NO2], [PPN], [1, 1], [1])   #JPL Eval 17  #  RCO3 + NO2 = PPN   
            Reaction(GCJPLEQ_acabab(9.00e-29, 14000.0e0, 9.00e-28, 8.9e0, 7.7e-12, 0.2e0, 0.6e0), [PPN], [RCO3, NO2], [1], [1, 1])  #  PPN = RCO3 + NO2   
            Reaction(6.50e-15, [RCHO, NO3], [HNO3, RCO3], [1, 1], [1, 1])  #  RCHO + NO3 = HNO3 + RCO3   
            Reaction(1.33e-13 + 3.82e-11*exp(-2000.0e0/288.15), [ACET, OH], [ATO2, H2O], [1, 1], [1, 1])   #JPL Eval 17, p1-62-D31; EVF  #  ACET + OH = ATO2 + H2O   
            Reaction(GCARR_ac(7.40e-13, 700.0e0), [R4O2, HO2], [R4P], [1, 1], [1])  #  R4O2 + HO2 = R4P   
            Reaction(GCARR_ac(7.40e-13, 700.0e0), [R4N1, HO2], [R4N2], [1, 1], [1])  #  R4N1 + HO2 = R4N2   
            Reaction(GC_RO2HO2_aca(2.91e-13, 1300.0e0, 3.0e0), [B3O2, HO2], [RB3P], [1, 1], [1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  B3O2 + HO2 = RB3P   
            Reaction(GC_RO2HO2_aca(2.91e-13, 1300.0e0, 3.0e0), [PRN1, HO2], [PRPN], [1, 1], [1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  PRN1 + HO2 = PRPN   
            Reaction(GCARR_ac(1.30e-12, -25.0e0), [MEK, OH], [KO2, H2O], [1, 1], [1, 1])  #  MEK + OH = KO2 + H2O   
            Reaction(8.00e-16, [MEK, NO3], [HNO3, KO2], [1, 1], [1, 1])  #  MEK + NO3 = HNO3 + KO2   
            Reaction(3.35e-12, [EOH, OH], [HO2, ALD2], [1, 1], [1, 1])   #2013/02/12; JPL 10-6; BHH,JMAO,EAM  #  EOH + OH = HO2 + ALD2   
            Reaction(GCARR_ac(4.60e-12, 70.0e0), [ROH, OH], [HO2, RCHO], [1, 1], [1, 1])  #  ROH + OH = HO2 + RCHO   
            Reaction(4.10e-14, [ETO2, ETO2], [ALD2, HO2], [1, 1], [2.000, 2.000])  #  ETO2 + ETO2 = 2.000ALD2 + 2.000HO2   
            Reaction(4.10e-14, [OTHRO2, OTHRO2], [ALD2, HO2], [1, 1], [2.000, 2.000])   #2019/05/10; Fisher2018; JAF  #  OTHRO2 + OTHRO2 = 2.000ALD2 + 2.000HO2   
            Reaction(2.70e-14, [ETO2, ETO2], [EOH, ALD2], [1, 1], [1, 1])  #  ETO2 + ETO2 = EOH + ALD2   
            Reaction(2.70e-14, [OTHRO2, OTHRO2], [EOH, ALD2], [1, 1], [1, 1])   #2019/05/10; Fisher2018; JAF  #  OTHRO2 + OTHRO2 = EOH + ALD2   
            Reaction(GCARR_ac(7.40e-13, 700.0e0), [HO2, ETO2], [ETP], [1, 1], [1])  #  HO2 + ETO2 = ETP   
            Reaction(GCARR_ac(7.40e-13, 700.0e0), [HO2, OTHRO2], [ETP], [1, 1], [1])   #2019/05/10; Fisher2018; JAF  #  HO2 + OTHRO2 = ETP   
            Reaction(GC_RO2HO2_aca(2.91e-13, 1300.0e0, 3.0e0), [A3O2, HO2], [RA3P], [1, 1], [1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  A3O2 + HO2 = RA3P   
            Reaction(GC_RO2HO2_aca(2.91e-13, 1300.0e0, 3.0e0), [PO2, HO2], [PP], [1, 1], [1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  PO2 + HO2 = PP   
            Reaction(GCJPLPR_abab(4.60e-27, 4.0e0, 2.6e-11, 1.3e0, 0.5e0), [PRPE, OH], [PO2], [1, 1], [1])   #2017/02/22; JPL 15-10; BHH,MJE  #  PRPE + OH = PO2   
            Reaction(GC_GLYCOH_B_a(8.00e-12), [GLYC, OH], [HCOOH, OH, CO], [1, 1], [1, 1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  GLYC + OH = HCOOH + OH + CO   
            Reaction(GCARR_ac(4.59e-13, -1156.0e0), [PRPE, NO3], [PRN1], [1, 1], [1])  #  PRPE + NO3 = PRN1   
            Reaction(GCARR_ac(3.10e-12, 340.0e0), [GLYX, OH], [HO2, CO], [1, 1], [1, 2.000])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  GLYX + OH = HO2 + 2.000CO   
            Reaction(1.50e-11, [MGLY, OH], [MCO3, CO], [1, 1], [1, 1])  #  MGLY + OH = MCO3 + CO   
            Reaction(GC_GLYXNO3_ac(1.40e-12, -1860.0e0), [GLYX, NO3], [HNO3, HO2, CO], [1, 1], [1, 1, 2.000])  #  GLYX + NO3 = HNO3 + HO2 + 2.000CO   
            Reaction(GCARR_ac(3.36e-12, -1860.0e0), [MGLY, NO3], [HNO3, CO, MCO3], [1, 1], [1, 1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  MGLY + NO3 = HNO3 + CO + MCO3   
            Reaction(GC_HACOH_A_ac(2.15e-12, 305.0e0), [HAC, OH], [MGLY, HO2], [1, 1], [1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  HAC + OH = MGLY + HO2   
            Reaction(GCARR_ac(1.68e-12, 500.0e0), [MCO3, A3O2], [MO2, RCHO, HO2], [1, 1], [1, 1, 1])  #  MCO3 + A3O2 = MO2 + RCHO + HO2   
            Reaction(GCARR_ac(1.68e-12, 500.0e0), [MCO3, PO2], [MO2, ALD2, CH2O, HO2], [1, 1], [1, 1, 1, 1])  #  MCO3 + PO2 = MO2 + ALD2 + CH2O + HO2   
            Reaction(GCARR_ac(1.87e-13, 500.0e0), [MCO3, A3O2], [ACTA, RCHO], [1, 1], [1, 1])  #  MCO3 + A3O2 = ACTA + RCHO   
            Reaction(GCARR_ac(1.87e-13, 500.0e0), [MCO3, PO2], [ACTA, RCHO, HAC], [1, 1], [1, 0.350, 0.650])  #  MCO3 + PO2 = ACTA + 0.350RCHO + 0.650HAC   
            Reaction(GCARR_ac(1.87e-13, 500.0e0), [RCO3, MO2], [RCOOH, CH2O], [1, 1], [1, 1])  #  RCO3 + MO2 = RCOOH + CH2O   
            Reaction(GCARR_ac(8.78e-12, 200.0e0), [R4P, OH], [OH, R4O2, RCHO], [1, 1], [0.791, 0.209, 0.791])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  R4P + OH = 0.791OH + 0.209R4O2 + 0.791RCHO   
            Reaction(GCARR_ac(6.13e-13, 200.0e0), [RP, OH], [RCO3], [1, 1], [1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  RP + OH = RCO3   
            Reaction(GCARR_ac(8.78e-12, 200.0e0), [PP, OH], [OH, PO2, HAC], [1, 1], [0.791, 0.209, 0.791])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  PP + OH = 0.791OH + 0.209PO2 + 0.791HAC   
            Reaction(GCARR_ac(4.82e-11, -400.0e0), [LVOC, OH], [OH], [1, 1], [1])   #2017/06/14; Marais2016; EAM  #  LVOC + OH = OH   
            Reaction(GCARR_ac(6.13e-13, 200.0e0), [OH, MAP], [MCO3], [1, 1], [1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  OH + MAP = MCO3   
            Reaction(1.40e-18, [C2H6, NO3], [ETO2, HNO3], [1, 1], [1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  C2H6 + NO3 = ETO2 + HNO3   
            Reaction(GCARR_ac(2.50e-12, 500.0e0), [MCO3, MCO3], [MO2], [1, 1], [2.000])  #  MCO3 + MCO3 = 2.000MO2   
            Reaction(GCARR_ac(1.80e-12, 500.0e0), [MCO3, MO2], [CH2O, MO2, HO2], [1, 1], [1, 1, 1])  #  MCO3 + MO2 = CH2O + MO2 + HO2   
            Reaction(GCARR_ac(2.00e-13, 500.0e0), [MCO3, MO2], [ACTA, CH2O], [1, 1], [1, 1])  #  MCO3 + MO2 = ACTA + CH2O   
            Reaction(GCARR_ac(1.68e-12, 500.0e0), [ATO2, MCO3], [MO2, MCO3, CH2O], [1, 1], [1, 1, 1])   #2013/03/22; Paulot2009; FP,EAM,JMAO,MJE  #  ATO2 + MCO3 = MO2 + MCO3 + CH2O   
            Reaction(GCARR_ac(1.68e-12, 500.0e0), [KO2, MCO3], [MO2, ALD2, MCO3], [1, 1], [1, 1, 1])  #  KO2 + MCO3 = MO2 + ALD2 + MCO3   
            Reaction(GCARR_ac(1.68e-12, 500.0e0), [B3O2, MCO3], [MO2, HO2, ACET], [1, 1], [1, 1, 1])  #  B3O2 + MCO3 = MO2 + HO2 + ACET   
            Reaction(GCARR_ac(1.68e-12, 500.0e0), [PRN1, MCO3], [MO2, NO2, CH2O, ALD2], [1, 1], [1, 1, 1, 1])  #  PRN1 + MCO3 = MO2 + NO2 + CH2O + ALD2   
            Reaction(GCARR_ac(1.87e-13, 500.0e0), [R4O2, MCO3], [MEK, ACTA], [1, 1], [1, 1])  #  R4O2 + MCO3 = MEK + ACTA   
            Reaction(GCARR_ac(1.87e-13, 500.0e0), [ATO2, MCO3], [MGLY, ACTA], [1, 1], [1, 1])   #2017/07/27; Fix C creation; SAS,BHH,MJE  #  ATO2 + MCO3 = MGLY + ACTA   
            Reaction(GCARR_ac(1.87e-13, 500.0e0), [KO2, MCO3], [MEK, ACTA], [1, 1], [1, 1])  #  KO2 + MCO3 = MEK + ACTA   
            Reaction(GCARR_ac(1.87e-13, 500.0e0), [R4N1, MCO3], [RCHO, ACTA, NO2], [1, 1], [1, 1, 1])  #  R4N1 + MCO3 = RCHO + ACTA + NO2   
            Reaction(GCARR_ac(1.87e-13, 500.0e0), [PRN1, MCO3], [RCHO, ACTA, NO2], [1, 1], [1, 1, 1])  #  PRN1 + MCO3 = RCHO + ACTA + NO2   
            Reaction(GCARR_ac(1.87e-13, 500.0e0), [B3O2, MCO3], [ACET, ACTA], [1, 1], [1, 1])  #  B3O2 + MCO3 = ACET + ACTA   
            Reaction(GCARR_ac(1.68e-12, 500.0e0), [MCO3, ETO2], [MO2, ALD2, HO2], [1, 1], [1, 1, 1])  #  MCO3 + ETO2 = MO2 + ALD2 + HO2   
            Reaction(GCARR_ac(1.68e-12, 500.0e0), [MCO3, OTHRO2], [MO2, ALD2, HO2], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF  #  MCO3 + OTHRO2 = MO2 + ALD2 + HO2   
            Reaction(GCARR_ac(1.87e-13, 500.0e0), [MCO3, ETO2], [ACTA, ALD2], [1, 1], [1, 1])  #  MCO3 + ETO2 = ACTA + ALD2   
            Reaction(GCARR_ac(1.87e-13, 500.0e0), [MCO3, OTHRO2], [ACTA, ALD2], [1, 1], [1, 1])   #2019/05/10; Fisher2018; JAF  #  MCO3 + OTHRO2 = ACTA + ALD2   
            Reaction(GCARR_ac(8.50e-13, -2450.0e0), [NO3, NO3], [NO2, O2], [1, 1], [2.000, 1])  #  NO3 + NO3 = 2.000NO2 + O2   
            Reaction(GCJPLPR_abab(1.00e-30, 4.8e+00, 7.2e-12, 2.1e0, 0.6e0), [MO2, NO2], [MPN], [1, 1], [1])   #2012/02/12; Browne2011; ECB  #  MO2 + NO2 = MPN   
            Reaction(GCJPLPR_abcabc(1.05e-02, 4.8e+00, -11234.0e0, 7.58e16, 2.1e0, -11234.0e0, 0.6e0), [MPN], [MO2, NO2], [1], [1, 1])   #2012/02/12; Browne2011; ECB  #  MPN = MO2 + NO2   
            Reaction(GCARR_ac(1.20e-11, -280.0e0), [DMS, OH], [SO2, MO2, CH2O], [1, 1], [1, 1, 1])  #  DMS + OH = SO2 + MO2 + CH2O   
            Reaction(GC_DMSOH_acac(8.20e-39, 5376.0e0, 1.05e-5, 3644.0e0), [DMS, OH], [SO2, MSA, MO2], [1, 1], [0.750, 0.250, 1])  #  DMS + OH = 0.750SO2 + 0.250MSA + MO2   
            Reaction(GCARR_ac(1.90e-13, 530.0e0), [DMS, NO3], [SO2, HNO3, MO2, CH2O], [1, 1], [1, 1, 1, 1])  #  DMS + NO3 = SO2 + HNO3 + MO2 + CH2O   
            Reaction(GCJPLPR_aba(3.30e-31, 4.3e+00, 1.6e-12, 0.6e0), [SO2, OH], [SO4, HO2], [1, 1], [1, 1])  #  SO2 + OH = SO4 + HO2   
            Reaction(GCARR_ac(1.60e-11, -780.0e0), [Br, O3], [BrO, O2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  Br + O3 = BrO + O2   
            Reaction(GCARR_ac(4.50e-12, 460.0e0), [BrO, HO2], [HOBr, O2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  BrO + HO2 = HOBr + O2   
            Reaction(GCARR_ac(4.80e-12, -310.0e0), [Br, HO2], [HBr, O2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  Br + HO2 = HBr + O2   
            Reaction(GCARR_ac(5.50e-12, 200.0e0), [HBr, OH], [Br, H2O], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  HBr + OH = Br + H2O   
            Reaction(GCARR_ac(2.40e-12,  40.0e0), [BrO, BrO], [Br, O2], [1, 1], [2.000, 1])   #2012/06/07; Parrella2012; JPP  #  BrO + BrO = 2.000Br + O2   
            Reaction(GCARR_ac(2.80e-14, 860.0e0), [BrO, BrO], [Br2, O2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  BrO + BrO = Br2 + O2   
            Reaction(GCARR_ac(8.80e-12, 260.0e0), [BrO, NO], [Br, NO2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  BrO + NO = Br + NO2   
            Reaction(4.90e-11, [Br, BrNO3], [Br2, NO3], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  Br + BrNO3 = Br2 + NO3   
            Reaction(GCARR_ac(2.10e-11, 240.0e0), [Br2, OH], [HOBr, Br], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  Br2 + OH = HOBr + Br   
            Reaction(GCARR_ac(1.20e-10, -430.0e0), [HOBr, O], [OH, BrO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  HOBr + O = OH + BrO   
            Reaction(GCARR_ac(5.80e-12, -1500.0e0), [HBr, O], [OH, Br], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  HBr + O = OH + Br   
            Reaction(GCARR_ac(1.70e-11, 250.0e0), [BrO, OH], [Br, HO2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  BrO + OH = Br + HO2   
            Reaction(1.60e-11, [Br, NO3], [BrO, NO2], [1, 1], [1, 1])   #2012/06/07; Parrella2012; JPP  #  Br + NO3 = BrO + NO2   
            Reaction(GCARR_ac(1.70e-11, -800.0e0), [Br, CH2O], [HBr, HO2, CO], [1, 1], [1, 1, 1])   #2012/06/07; Parrella2012; JPP  #  Br + CH2O = HBr + HO2 + CO   
            Reaction(GCARR_ac(1.80e-11, -460.0e0), [Br, ALD2], [HBr, MCO3], [1, 1], [1, 1])   #2017/07/27; Parrella2012,Fix C creation; SAS,BHH,MJE  #  Br + ALD2 = HBr + MCO3   
            Reaction(GCARR_ac(1.66e-10, -7000.0e0), [Br, ACET], [HBr, ATO2], [1, 1], [1, 1])   #2017/07/27; Parrella2012,Fix C creation; SAS,BHH,MJE  #  Br + ACET = HBr + ATO2   
            Reaction(GCARR_ac(2.36e-10, -6411.0e0), [Br, C2H6], [HBr, ETO2], [1, 1], [1, 1])   #2017/07/27; Parrella2012,Fix C creation; SAS,BHH,MJE  #  Br + C2H6 = HBr + ETO2   
            Reaction(GCARR_ac(8.77e-11, -4330.0e0), [Br, C3H8], [HBr, A3O2], [1, 1], [1, 1])   #2017/07/27; Parrella2012,Fix C creation; SAS,BHH,MJE  #  Br + C3H8 = HBr + A3O2   
            Reaction(GCJPLPR_aba(4.20e-31, 2.4e0, 2.7e-11, 0.6e0), [Br, NO2], [BrNO2], [1, 1], [1])   #2012/06/07; Parrella2012; JPP  #  Br + NO2 = BrNO2   
            Reaction(GCJPLPR_abab(5.40e-31, 3.1e0, 6.5e-12, 2.9e0, 0.6e0), [BrO, NO2], [BrNO3], [1, 1], [1])   #2017/02/22; JPL 15-10; BHH,MJE  #  BrO + NO2 = BrNO3   
            Reaction(GCARR_ac(9.00e-13, -360.0e0), [CHBr3, OH], [Br], [1, 1], [3.000])   #2017/02/22; JPL 15-10; BHH,MJE  #  CHBr3 + OH = 3.000Br   
            Reaction(GCARR_ac(2.00e-12, -840.0e0), [CH2Br2, OH], [Br], [1, 1], [2.000])   #2012/06/07; Parrella2012; JPP  #  CH2Br2 + OH = 2.000Br   
            Reaction(GCARR_ac(1.42e-12, -1150.0e0), [CH3Br, OH], [Br, H2O, HO2], [1, 1], [1, 1, 1])   #2017/03/08; JPL 15-10; TS,BHH,MJE  #  CH3Br + OH = Br + H2O + HO2   
            Reaction(GCARR_ac(1.63e-10, 60.0e0), [O1D, H2O], [OH], [1, 1], [2.000])   #2014/02/03; Eastham2014; SDE  #  O1D + H2O = 2.000OH   
            Reaction(GCARR_ac(2.15e-11, 110.0e0), [O1D, N2], [O, N2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  O1D + N2 = O + N2   
            Reaction(GCARR_ac(3.30e-11, 55.0e0), [O1D, O2], [O, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  O1D + O2 = O + O2   
            Reaction(1.20e-10, [O1D, H2], [H, OH], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  O1D + H2 = H + OH   
            Reaction(GCARR_ac(4.63e-11, 20.0e0), [O1D, N2O], [N2, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  O1D + N2O = N2 + O2   
            Reaction(GCARR_ac(7.25e-11, 20.0e0), [O1D, N2O], [NO], [1, 1], [2.000])   #2014/02/03; Eastham2014; SDE  #  O1D + N2O = 2.000NO   
            Reaction(1.31e-10, [O1D, CH4], [MO2, OH], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  O1D + CH4 = MO2 + OH   
            Reaction(0.09e-10, [O1D, CH4], [CH2O, H2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  O1D + CH4 = CH2O + H2   
            Reaction(0.35e-10, [O1D, CH4], [CH2O, H, HO2], [1, 1], [1, 1, 1])   #2014/02/03; Eastham2014; SDE  #  O1D + CH4 = CH2O + H + HO2   
            Reaction(GCARR_ab(6.00e-34, 2.4e0)*NUMDEN, [O, O2], [O3], [1, 1], [1])   #2014/02/03; Eastham2014; SDE  #  O + O2 = O3   
            Reaction(GCARR_ac(8.00e-12, -2060.0e0), [O, O3], [O2], [1, 1], [2.000])   #2014/02/03; Eastham2014; SDE  #  O + O3 = 2.000O2   
            Reaction(GCARR_ac(2.80e-12, -1800.0e0), [OH, H2], [H2O, H], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + H2 = H2O + H   
            Reaction(GCARR_ac(1.80e-11, 180.0e0), [O, OH], [O2, H], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  O + OH = O2 + H   
            Reaction(GCARR_ac(3.00e-11, 200.0e0), [HO2, O], [OH, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  HO2 + O = OH + O2   
            Reaction(1.20e-10, [O1D, O3], [O2], [1, 1], [2.000])   #2014/02/03; Eastham2014; SDE  #  O1D + O3 = 2.000O2   
            Reaction(1.20e-10, [O1D, O3], [O, O2], [1, 1], [2.000, 1])   #2014/02/03; Eastham2014; SDE  #  O1D + O3 = 2.000O + O2   
            Reaction(GCARR_ac(2.10e-11, -2200.0e0), [OCS, O], [CO, SO2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OCS + O = CO + SO2   
            Reaction(GCARR_ac(1.10e-13, -1200.0e0), [OCS, OH], [CO2, SO2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OCS + OH = CO2 + SO2   
            Reaction(GCARR_ac(5.10e-12, 210.0e0), [NO2, O], [NO, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  NO2 + O = NO + O2   
            Reaction(1.00e-11, [NO3, O], [NO2, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  NO3 + O = NO2 + O2   
            Reaction(GCJPLPR_aba(9.00e-32, 1.5e+00, 3.0e-11, 0.6e0), [NO, O], [NO2], [1, 1], [1])   #2014/02/03; Eastham2014; SDE  #  NO + O = NO2   
            Reaction(GCJPLPR_abab(2.50e-31, 1.8e+00, 2.2e-11, 0.7e0, 0.6e0), [NO2, O], [NO3], [1, 1], [1])   #2014/02/03; Eastham2014; SDE  #  NO2 + O = NO3   
            Reaction(GCARR_ac(1.40e-12, -2000.0e0), [H2O2, O], [OH, HO2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  H2O2 + O = OH + HO2   
            Reaction(GCJPLPR_abab(4.40e-32, 1.3e+00, 7.5e-11, -0.2e0, 0.6e0), [H, O2], [HO2], [1, 1], [1])   #2014/02/03; Eastham2014; SDE  #  H + O2 = HO2   
            Reaction(GCARR_ac(1.40e-10, -470.0e0), [H, O3], [OH, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  H + O3 = OH + O2   
            Reaction(7.20e-11, [H, HO2], [OH], [1, 1], [2.000])   #2014/02/03; Eastham2014; SDE  #  H + HO2 = 2.000OH   
            Reaction(1.60e-12, [H, HO2], [O, H2O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  H + HO2 = O + H2O   
            Reaction(6.90e-12, [H, HO2], [H2, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  H + HO2 = H2 + O2   
            Reaction(GCARR_ac(1.50e-11, -3600.0e0), [N, O2], [NO, O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  N + O2 = NO + O   
            Reaction(GCARR_ac(2.10e-11, 100.0e0), [N, NO], [N2, O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  N + NO = N2 + O   
            Reaction(GCARR_ac(5.80e-12, 220.0e0), [N, NO2], [N2O, O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  N + NO2 = N2O + O   
            Reaction(GCARR_ac(1.90e-11, 230.0e0), [BrO, O], [Br, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  BrO + O = Br + O2   
            Reaction(GCARR_ac(3.40e-11, -1600.0e0), [CH2O, O], [CO, HO2, OH], [1, 1], [1, 1, 1])   #2014/02/03; Eastham2014; SDE  #  CH2O + O = CO + HO2 + OH   
            Reaction(1.80e-10, [O1D, CH3Br], [BrO, MO2, Br], [1, 1], [0.440, 1, 0.560])   #2014/02/03; Eastham2014; SDE  #  O1D + CH3Br = 0.440BrO + MO2 + 0.560Br   
            Reaction(GCARR_ac(2.60e-12, -1100.0e0), [OH, Cl2], [HOCl, Cl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + Cl2 = HOCl + Cl   
            Reaction(GCARR_ac(1.80e-11, -600.0e0), [MO2, ClO], [ClOO, HO2, CH2O], [1, 1], [1, 1, 1])   #2017/03/20; JPL 15-10; TS,BHH,MJE  #  MO2 + ClO = ClOO + HO2 + CH2O   
            Reaction(GCARR_ac(7.40e-12, 270.0e0), [OH, ClO], [HO2, Cl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + ClO = HO2 + Cl   
            Reaction(GCARR_ac(6.00e-13, 230.0e0), [OH, ClO], [HCl, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + ClO = HCl + O2   
            Reaction(GCARR_ac(1.40e-12, 600.0e0), [OH, OClO], [HOCl, O2], [1, 1], [1, 1])   #2017/02/22; JPL 15-10; BHH,MJE  #  OH + OClO = HOCl + O2   
            Reaction(GCARR_ac(6.00e-13, 670.0e0), [OH, Cl2O2], [HOCl, ClOO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + Cl2O2 = HOCl + ClOO   
            Reaction(GCARR_ac(1.80e-12, -250.0e0), [OH, HCl], [H2O, Cl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + HCl = H2O + Cl   
            Reaction(GCARR_ac(3.00e-12, -500.0e0), [OH, HOCl], [H2O, ClO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + HOCl = H2O + ClO   
            Reaction(GCARR_ac(2.40e-12, -1250.0e0), [OH, ClNO2], [HOCl, NO2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + ClNO2 = HOCl + NO2   
            Reaction(GCARR_ac(1.20e-12, -330.0e0), [OH, ClNO3], [HOCl, NO3], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + ClNO3 = HOCl + NO3   
            Reaction(GCARR_ac(1.96e-12, -1200.0e0), [OH, CH3Cl], [Cl, HO2, H2O], [1, 1], [1, 1, 1])   #2017/02/22; JPL 15-10; BHH,MJE  #  OH + CH3Cl = Cl + HO2 + H2O   
            Reaction(GCARR_ac(2.61e-12, -944.0e0), [OH, CH2Cl2], [Cl, HO2], [1, 1], [2.000, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  OH + CH2Cl2 = 2.000Cl + HO2   
            Reaction(GCARR_ac(4.69e-12, -1134.0e0), [OH, CHCl3], [Cl, HO2], [1, 1], [3.000, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  OH + CHCl3 = 3.000Cl + HO2   
            Reaction(GCARR_ac(1.64e-12, -1520.0e0), [OH, CH3CCl3], [Cl, H2O], [1, 1], [3.000, 1])   #2014/02/03; Eastham2014; SDE  #  OH + CH3CCl3 = 3.000Cl + H2O   
            Reaction(GCARR_ac(9.20e-13, -1560.0e0), [OH, HCFC22], [Cl, H2O], [1, 1], [1, 1])   #2017/02/22; JPL 15-10; BHH,MJE  #  OH + HCFC22 = Cl + H2O   
            Reaction(GCARR_ac(1.25e-12, -1600.0e0), [OH, HCFC141b], [Cl, H2O], [1, 1], [2.000, 1])   #2014/02/03; Eastham2014; SDE  #  OH + HCFC141b = 2.000Cl + H2O   
            Reaction(GCARR_ac(1.30e-12, -1770.0e0), [OH, HCFC142b], [Cl, H2O], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  OH + HCFC142b = Cl + H2O   
            Reaction(GCARR_ac(7.40e-13, -900.0e0), [OH, HCFC123], [Cl, H2O], [1, 1], [2.000, 1])   #2017/02/22; JPL 15-10; BHH,MJE  #  OH + HCFC123 = 2.000Cl + H2O   
            Reaction(GCARR_ac(7.10e-12, -1270.0e0), [CH4, Cl], [HCl, MO2], [1, 1], [1, 1])   #2017/03/08; JPL 15-10; TS,BHH,MJE  #  CH4 + Cl = HCl + MO2   
            Reaction(GCARR_ac(7.32e-11, -30.0e0), [CH2O, Cl], [CO, HCl, HO2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b; TS,JAS,SDE  #  CH2O + Cl = CO + HCl + HO2   
            Reaction(GCARR_ac(2.30e-11, -200.0e0), [Cl, O3], [ClO, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  Cl + O3 = ClO + O2   
            Reaction(GCARR_ac(3.05e-11, -2270.0e0), [Cl, H2], [H, HCl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  Cl + H2 = H + HCl   
            Reaction(GCARR_ac(1.10e-11, -980.0e0), [Cl, H2O2], [HO2, HCl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  Cl + H2O2 = HO2 + HCl   
            Reaction(GCARR_ac(1.40e-11, 270.0e0), [Cl, HO2], [O2, HCl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  Cl + HO2 = O2 + HCl   
            Reaction(GCARR_ac(3.60e-11, -375.0e0), [Cl, HO2], [OH, ClO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  Cl + HO2 = OH + ClO   
            Reaction(GCARR_ac(2.80e-11, 85.0e0), [ClO, O], [Cl, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClO + O = Cl + O2   
            Reaction(GCARR_ac(2.60e-12, 290.0e0), [ClO, HO2], [O2, HOCl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClO + HO2 = O2 + HOCl   
            Reaction(GCARR_ac(6.40e-12, 290.0e0), [ClO, NO], [Cl, NO2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClO + NO = Cl + NO2   
            Reaction(GCJPLPR_abab(1.80e-31, 3.4e+00, 1.50e-11, 1.9e0, 0.6e0), [ClO, NO2], [ClNO3], [1, 1], [1])   #2014/02/03; Eastham2014; SDE  #  ClO + NO2 = ClNO3   
            Reaction(GCARR_ac(1.00e-12, -1590.0e0), [ClO, ClO], [Cl2, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClO + ClO = Cl2 + O2   
            Reaction(GCARR_ac(3.00e-11, -2450.0e0), [ClO, ClO], [Cl, ClOO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClO + ClO = Cl + ClOO   
            Reaction(GCARR_ac(3.50e-13, -1370.0e0), [ClO, ClO], [OClO, Cl], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClO + ClO = OClO + Cl   
            Reaction(GCJPLPR_aba(2.20e-33, 3.1e+00, 1.8e-10, 0.6e0), [Cl, O2], [ClOO], [1, 1], [1])   #2014/02/03; Eastham2014; SDE  #  Cl + O2 = ClOO   
            Reaction(GCJPLEQ_acabab(6.60e-25, 2502.0e0, 2.20e-33, 3.1e+00, 1.8e-10, 0.0e0, 0.6e0), [ClOO], [Cl, O2], [1], [1, 1])   #JPL 15-10; XW  #  ClOO = Cl + O2   
            Reaction(GCJPLPR_abab(1.90e-32, 3.6e+00, 3.7e-12, 1.6e0, 0.6e0), [ClO, ClO], [Cl2O2], [1, 1], [1])   #2017/02/22; JPL 15-10; BHH,MJE  #  ClO + ClO = Cl2O2   
            Reaction(GCJPLEQ_acabab(2.16e-27, 8537.0e0, 1.90e-32, 3.6e+00, 3.7e-12, 1.6e0, 0.6e0), [Cl2O2], [ClO], [1], [2.000])   #JPL 15-10; XW  #  Cl2O2 = 2.000ClO   
            Reaction(2.30e-10, [ClOO, Cl], [Cl2, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClOO + Cl = Cl2 + O2   
            Reaction(1.20e-11, [ClOO, Cl], [ClO], [1, 1], [2.000])   #2014/02/03; Eastham2014; SDE  #  ClOO + Cl = 2.000ClO   
            Reaction(GCARR_ac(9.50e-13, 550.0e0), [ClO, BrO], [Br, OClO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClO + BrO = Br + OClO   
            Reaction(GCARR_ac(2.30e-12, 260.0e0), [ClO, BrO], [Br, ClOO], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClO + BrO = Br + ClOO   
            Reaction(GCARR_ac(4.10e-13, 290.0e0), [ClO, BrO], [BrCl, O2], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClO + BrO = BrCl + O2   
            Reaction(GCARR_ac(3.60e-12, -840.0e0), [ClNO3, O], [ClO, NO3], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClNO3 + O = ClO + NO3   
            Reaction(GCARR_ac(6.50e-12, 135.0e0), [ClNO3, Cl], [Cl2, NO3], [1, 1], [1, 1])   #2014/02/03; Eastham2014; SDE  #  ClNO3 + Cl = Cl2 + NO3   
            Reaction(GCARR_ac(2.17e-11, -1130.0e0), [CH3Cl, Cl], [CO, HCl, HO2], [1, 1], [1, 2.000, 1])   #2014/02/03; Eastham2014; SDE  #  CH3Cl + Cl = CO + 2.000HCl + HO2   
            Reaction(GCARR_ac(1.24e-12, -1070.0e0), [CH2Cl2, Cl], [CO, HCl, Cl, HO2], [1, 1], [1, 1, 2.000, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  CH2Cl2 + Cl = CO + HCl + 2.000Cl + HO2   
            Reaction(GCARR_ac(3.77e-12, -1011.0e0), [CHCl3, Cl], [CO, HCl, Cl, HO2], [1, 1], [1, 1, 3.000, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  CHCl3 + Cl = CO + HCl + 3.000Cl + HO2   
            Reaction(2.00e-13, [Cl, HCOOH], [HCl, CO2, H2O], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + HCOOH = HCl + CO2 + H2O   
            Reaction(1.60e-10, [Cl, MO2], [ClO, CH2O, HO2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + MO2 = ClO + CH2O + HO2   
            Reaction(5.7e-11, [Cl, MP], [HCl, MO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + MP = HCl + MO2   
            Reaction(GCARR_ac(7.2e-11, -70.0e0), [Cl, C2H6], [HCl, ETO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + C2H6 = HCl + ETO2   
            Reaction(7.4e-11, [Cl, ETO2], [ClO, HO2, ALD2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + ETO2 = ClO + HO2 + ALD2   
            Reaction(7.4e-11, [Cl, OTHRO2], [ClO, HO2, ALD2], [1, 1], [1, 1, 1])   #2019/05/10; Fisher2018; JAF  #  Cl + OTHRO2 = ClO + HO2 + ALD2   
            Reaction(5.5e-11, [Cl, MOH], [HCl, CH2O, HO2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + MOH = HCl + CH2O + HO2   
            Reaction(9.6e-11, [Cl, EOH], [HCl, ALD2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + EOH = HCl + ALD2   
            Reaction(2.8e-14, [Cl, ACTA], [HCl, MO2, CO2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + ACTA = HCl + MO2 + CO2   
            Reaction(GCARR_ac(6.54e-11, 60.0e0), [Cl, C3H8], [HCl, B3O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + C3H8 = HCl + B3O2   
            Reaction(GCARR_ac(8.12e-11, -90.0e0), [Cl, C3H8], [HCl, A3O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + C3H8 = HCl + A3O2   
            Reaction(GCARR_ac(7.70e-11, -1000.0e0), [Cl, ACET], [HCl, ATO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + ACET = HCl + ATO2   
            Reaction(GCARR_ac(7.60e-11, 500.0e0), [Cl, ISOP], [HCl, IHOO1, IHOO4], [1, 1], [1, 0.5, 0.5])   #2019/11/06; Sherwen2016b;KHB,TS,JAS,SDE  #  Cl + ISOP = HCl + 0.5IHOO1 + 0.5IHOO4   
            Reaction(2.05e-10, [Cl, ALK4], [HCl, R4O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + ALK4 = HCl + R4O2   
            Reaction(GCJPLPR_aa(4.00e-28, 2.8e-10, 0.6e0), [Cl, PRPE], [HCl, PO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Cl + PRPE = HCl + PO2   
            Reaction(3.60e-12, [Br, PRPE], [HBr, PO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  Br + PRPE = HBr + PO2   
            Reaction(GCJPLPR_aba(1.80e-32, 1.0e0, 1.77e-11, 0.6e0), [I, NO], [INO], [1, 1], [1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I + NO = INO   
            Reaction(GCARR_ac(8.40e-11, -2620.0e0), [INO, INO], [I2, NO], [1, 1], [1, 2.000])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  INO + INO = I2 + 2.000NO   
            Reaction(GCJPLPR_aba(3.00e-31, 1.0e0, 6.6e-11, 0.63e0), [I, NO2], [IONO], [1, 1], [1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I + NO2 = IONO   
            Reaction(GCARR_ac(9.94e+17, -11859.0e0), [IONO], [I, NO2], [1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IONO = I + NO2   
            Reaction(GCARR_ac(2.90e-11, -2600.0e0), [IONO, IONO], [I2, NO2], [1, 1], [1, 2.000])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IONO + IONO = I2 + 2.000NO2   
            Reaction(1.50e-12, [I2, NO3], [I, IONO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I2 + NO3 = I + IONO2   
            Reaction(GCJPLPR_abab(7.50e-31, 3.5e0, 7.6e-12, 1.5e0, 0.6e0), [IO, NO2], [IONO2], [1, 1], [1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IO + NO2 = IONO2   
            Reaction(GCARR_ac(2.10e+15, -13670.0e0), [IONO2], [IO, NO2], [1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IONO2 = IO + NO2   
            Reaction(GCARR_ac(9.10e-11, -146.0e0), [IONO2, I], [I2, NO3], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IONO2 + I = I2 + NO3   
            Reaction(1.20e-11, [I, BrO], [IO, Br], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I + BrO = IO + Br   
            Reaction(GCARR_ac(3.00e-12, 510.0e0), [IO, BrO], [Br, I, O2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IO + BrO = Br + I + O2   
            Reaction(GCARR_ac(1.20e-11, 510.0e0), [IO, BrO], [Br, OIO], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IO + BrO = Br + OIO   
            Reaction(1.00e-10, [IO, OIO], [I2O3], [1, 1], [1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IO + OIO = I2O3   
            Reaction(1.50e-10, [OIO, OIO], [I2O4], [1, 1], [1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  OIO + OIO = I2O4   
            Reaction(3.80e-02, [I2O4], [OIO], [1], [2.000])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I2O4 = 2.000OIO   
            Reaction(GCARR_ac(1.10e-12, 542.0e0), [OIO, NO], [IO, NO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  OIO + NO = IO + NO2   
            Reaction(GCARR_ac(5.10e-12, 280.0e0), [IO, ClO], [I, OClO], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IO + ClO = I + OClO   
            Reaction(GCARR_ac(2.81e-12, 280.0e0), [IO, ClO], [I, Cl, O2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IO + ClO = I + Cl + O2   
            Reaction(GCARR_ac(1.02e-12, 280.0e0), [IO, ClO], [ICl, O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IO + ClO = ICl + O2   
            Reaction(GCARR_ac(2.30e-11, -870.0e0), [I, O3], [IO, O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2017;TS,JAS,SDE  #  I + O3 = IO + O2   
            Reaction(GCARR_ac(1.50e-11, -1090.0e0), [I, HO2], [HI, O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I + HO2 = HI + O2   
            Reaction(1.80e-10, [I2, OH], [HOI, I], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I2 + OH = HOI + I   
            Reaction(3.00e-11, [HI, OH], [I, H2O], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  HI + OH = I + H2O   
            Reaction(5.00e-12, [HOI, OH], [IO, H2O], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  HOI + OH = IO + H2O   
            Reaction(GCARR_ac(1.30e-11, 570.0e0), [IO, HO2], [HOI, O2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IO + HO2 = HOI + O2   
            Reaction(GCARR_ac(9.10e-12, 240.0e0), [IO, NO], [I, NO2], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IO + NO = I + NO2   
            Reaction(GCARR_ac(6.00e-12, 500.0e0), [IO, IO], [I, OIO], [1, 1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IO + IO = I + OIO   
            Reaction(GCARR_ac(9.00e-12, 500.0e0), [IO, IO], [I2O2], [1, 1], [1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  IO + IO = I2O2   
            Reaction(GCARR_ac(1.00e+12, -9770.0e0), [I2O2], [IO], [1], [2.000])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I2O2 = 2.000IO   
            Reaction(GCARR_ac(2.50e+14, -9770.0e0), [I2O2], [OIO, I], [1], [1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  I2O2 = OIO + I   
            Reaction(GCARR_ac(2.90e-12, -1100.0e0), [CH3I, OH], [H2O, I, MO2], [1, 1], [1, 1, 1])   #2017/09/22; Sherwen2016b;TS,JAS,SDE  #  CH3I + OH = H2O + I + MO2   
            Reaction(2.40e-12, [ETHLN, OH], [CH2O, CO2, NO2], [1, 1], [1, 1, 1])   #2017/06/15, Marais2016, EAM  #  ETHLN + OH = CH2O + CO2 + NO2   
            Reaction(6.70e-13, [PROPNN, OH], [NO2, MGLY], [1, 1], [1, 1])   #2017/07/14; MCMv3.3; KRT,JAF,CCM,EAM,KHB,RHS  #  PROPNN + OH = NO2 + MGLY   
            Reaction(1.20e-15, [CH2OO, CO], [CH2O], [1, 1], [1])   #2015/09/25; Millet2015; DBM,EAM  #  CH2OO + CO = CH2O   
            Reaction(1.00e-14, [CH2OO, NO], [CH2O, NO2], [1, 1], [1, 1])   #2015/09/25; Millet2015; DBM,EAM  #  CH2OO + NO = CH2O + NO2   
            Reaction(1.00e-15, [CH2OO, NO2], [CH2O, NO3], [1, 1], [1, 1])   #2015/09/25; Millet2015; DBM,EAM  #  CH2OO + NO2 = CH2O + NO3   
            Reaction(1.40e-12, [CH2OO, O3], [CH2O], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  CH2OO + O3 = CH2O   
            Reaction(3.70e-11, [CH2OO, SO2], [CH2O, SO4], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  CH2OO + SO2 = CH2O + SO4   
            Reaction(1.20e-15, [CH3CHOO, CO], [ALD2], [1, 1], [1])   #2015/09/25; Millet2015; DBM,EAM  #  CH3CHOO + CO = ALD2   
            Reaction(1.00e-14, [CH3CHOO, NO], [ALD2, NO2], [1, 1], [1, 1])   #2015/09/25; Millet2015; DBM,EAM  #  CH3CHOO + NO = ALD2 + NO2   
            Reaction(1.00e-15, [CH3CHOO, NO2], [ALD2, NO3], [1, 1], [1, 1])   #2015/09/25; Millet2015; DBM,EAM  #  CH3CHOO + NO2 = ALD2 + NO3   
            Reaction(7.00e-14, [CH3CHOO, SO2], [ALD2, SO4], [1, 1], [1, 1])   #2015/09/25; Millet2015; DBM,EAM  #  CH3CHOO + SO2 = ALD2 + SO4   
            Reaction(6.00e-18, [CH3CHOO, H2O], [ALD2, H2O2], [1, 1], [1, 1])   #2015/09/25; Millet2015; DBM,EAM  #  CH3CHOO + H2O = ALD2 + H2O2   
            Reaction(1.00e-17, [CH3CHOO, H2O], [ACTA], [1, 1], [1])   #2015/09/25; Millet2015; DBM,EAM  #  CH3CHOO + H2O = ACTA   
            Reaction(GCARR_ac(1.21e-11, 440.0e0), [MTPA, OH], [PIO2], [1, 1], [1])   #2017/03/23; IUPAC2010; EVF  #  MTPA + OH = PIO2   
            Reaction(GCARR_ac(1.21e-11, 440.0e0), [MTPO, OH], [PIO2], [1, 1], [1])   #2017/03/23; IUPAC2010; EVF  #  MTPO + OH = PIO2   
            Reaction(1.50e-11, [PIO2, HO2], [PIP], [1, 1], [1])   #2017/03/23; Roberts1992; EVF  #  PIO2 + HO2 = PIP   
            Reaction(1.20e-12, [PIO2, NO3], [HO2, NO2, RCHO, MEK], [1, 1], [1, 1, 1, 1])   #2017/03/23; Roberts1992; EVF  #  PIO2 + NO3 = HO2 + NO2 + RCHO + MEK   
            Reaction(GCARR_ac(8.33e-13, 490.0e0), [MTPA, NO3], [OLNN, OLND], [1, 1], [0.100, 0.900])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS  #  MTPA + NO3 = 0.100OLNN + 0.900OLND   
            Reaction(GCARR_ac(8.33e-13, 490.0e0), [MTPO, NO3], [OLNN, OLND], [1, 1], [0.100, 0.900])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS  #  MTPO + NO3 = 0.100OLNN + 0.900OLND   
            Reaction(GCARR_ac(4.20e-11, 401.0e0), [LIMO, OH], [LIMO2], [1, 1], [1])   #2017/07/14; Gill2002; KRT,JAF,CCM,EAM,KHB,RHS  #  LIMO + OH = LIMO2   
            Reaction(1.22e-11, [LIMO, NO3], [OLNN, OLND], [1, 1], [0.500, 0.500])   #2017/07/14; Fry2014,Atkinson2003; KRT,JAF,CCM,EAM,KHB,RHS  #  LIMO + NO3 = 0.500OLNN + 0.500OLND   
            Reaction(1.50e-11, [LIMO2, HO2], [PIP], [1, 1], [1])   #2017/07/14; Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS  #  LIMO2 + HO2 = PIP   
            Reaction(4.00e-12, [OLNN, NO], [HO2, NO2, MONITS], [1, 1], [1, 1, 1])   #2017/07/14; Browne2014,Goliff2013; KRT,JAF,CCM,EAM,KHB,RHS  #  OLNN + NO = HO2 + NO2 + MONITS   
            Reaction(GCARR_ac(1.66e-13, 1300.0e0), [OLNN, HO2], [MONITS, MONITU], [1, 1], [0.700, 0.300])   #2017/07/14; Browne2014,Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS  #  OLNN + HO2 = 0.700MONITS + 0.300MONITU   
            Reaction(GCARR_ac(1.66e-13, 1300.0e0), [OLND, HO2], [MONITS, MONITU], [1, 1], [0.700, 0.300])   #2017/07/14; Browne2014,Roberts1992; KRT,JAF,CCM,EAM,KHB,RHS  #  OLND + HO2 = 0.700MONITS + 0.300MONITU   
            Reaction(4.80e-12, [MONITS, OH], [HONIT], [1, 1], [1])   #2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS  #  MONITS + OH = HONIT   
            Reaction(7.29e-11, [MONITU, OH], [HONIT], [1, 1], [1])   #2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS  #  MONITU + OH = HONIT   
            Reaction(1.67e-16, [MONITU, O3], [HONIT], [1, 1], [1])   #2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS  #  MONITU + O3 = HONIT   
            Reaction(GCARR_ac(3.15e-13, -448.0e0), [MONITU, NO3], [HONIT], [1, 1], [1])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS  #  MONITU + NO3 = HONIT   
            Reaction(GCARR_ac(3.15e-13, -448.0e0), [MONITS, NO3], [HONIT], [1, 1], [1])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS  #  MONITS + NO3 = HONIT   
            Reaction(2.78e-04, [IONITA], [INDIOL, HNO3], [1], [1, 1])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS  #  IONITA = INDIOL + HNO3   
            Reaction(2.78e-04, [MONITA], [INDIOL, HNO3], [1], [1, 1])   #2017/07/14; Fisher2016; KRT,JAF,CCM,EAM,KHB,RHS  #  MONITA = INDIOL + HNO3   
            Reaction(GC_OHHNO3_acacac(2.41e-14, 460.0e0, 2.69e-17, 2199.0e0, 6.51e-34, 1335.0e0), [HONIT, OH], [NO3, HAC], [1, 1], [1, 1])   #2017/07/14; Browne2014; KRT,JAF,CCM,EAM,KHB,RHS  #  HONIT + OH = NO3 + HAC   
            Reaction(GCARR_ac(8.00e-13, -1000.0e0), [MENO3, OH], [CH2O, NO2], [1, 1], [1, 1])   #2019/05/16; JPL 15-10,Fisher2018; JAF  #  MENO3 + OH = CH2O + NO2   
            Reaction(GCARR_ac(1.00e-12, -490.0e0), [ETNO3, OH], [ALD2, NO2], [1, 1], [1, 1])   #2019/05/16; JPL 15-10,Fisher2018; JAF  #  ETNO3 + OH = ALD2 + NO2   
            Reaction(GCARR_ac(1.20e-12, -320.0e0), [IPRNO3, OH], [ACET, NO2], [1, 1], [1, 1])   #2019/05/16; JPL 15-10,Fisher2018; JAF  #  IPRNO3 + OH = ACET + NO2   
            Reaction(7.10e-13, [NPRNO3, OH], [RCHO, NO2], [1, 1], [1, 1])   #2019/05/16; JPL 15-10,Fisher2018; JAF  #  NPRNO3 + OH = RCHO + NO2   
            Reaction(GC_ISO1(1.7e-11, 3.90e2, 9.33e-2, 5.05e15, -1.22e4, 1.79e14, -8.830e3), [ISOP, OH], [LISOPOH, IHOO1], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  ISOP + OH = LISOPOH + IHOO1   
            Reaction(GC_ISO1(1.0e-11, 3.90e2, 2.26e-1, 2.22e9, -7.160e3, 1.75e14, -9.054e3), [ISOP, OH], [LISOPOH, IHOO4], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  ISOP + OH = LISOPOH + IHOO4   
            Reaction(ARRPLUS_abde(2.12e-13, -1300e0, -0.1644e0, 7.0485e-4), [IHOO1, HO2], [RIPC], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHOO1 + HO2 = RIPC   
            Reaction(ARRPLUS_abde(2.12e-13, -1300e0, -0.2038e0, 9.0435e-4), [IHOO4, HO2], [RIPD], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHOO4 + HO2 = RIPD   
            Reaction(ARRPLUS_abde(1.04e11, 9.746e3,  1.1644e0, -7.0485e-4), [IHOO1], [CH2O, OH, MVK], [1], [1, 1, 1])   #2019/11/06; Bates2019; KHB  #  IHOO1 = CH2O + OH + MVK   
            Reaction(ARRPLUS_abde(1.88e11, 9.752e3, 1.2038e0, -9.0435e-4), [IHOO4], [MACR, OH, CH2O], [1], [1, 1, 1])   #2019/11/06; Bates2019; KHB  #  IHOO4 = MACR + OH + CH2O   
            Reaction(ARRPLUS_ade(6.92e-14, 1.1644e0, -7.0485e-4), [IHOO1, IHOO1], [MVK, HO2, CH2O], [1, 1], [2, 2, 2])   #2019/11/06; Bates2019; KHB  #  IHOO1 + IHOO1 = 2MVK + 2HO2 + 2CH2O   
            Reaction(ARRPLUS_ade(5.74e-12, 1.2038e0, -9.0435e-4), [IHOO4, IHOO4], [MACR, HO2, CH2O], [1, 1], [2, 2, 2])   #2019/11/06; Bates2019; KHB  #  IHOO4 + IHOO4 = 2MACR + 2HO2 + 2CH2O   
            Reaction(ARRPLUS_ade(1.54e-12, 2.3682e0, -1.6092e-3), [IHOO1, IHOO4], [MACR, MVK, HO2, CH2O], [1, 1], [1, 1, 2, 2])   #2019/11/06; Bates2019; KHB  #  IHOO1 + IHOO4 = MACR + MVK + 2HO2 + 2CH2O   
            Reaction(ARRPLUS_ade(2.0e-12, 1.1644e0, -7.0485e-4), [IHOO1, MO2], [MVK, HO2, CH2O], [1, 1], [1, 2, 2])   #2019/11/06; Bates2019; KHB  #  IHOO1 + MO2 = MVK + 2HO2 + 2CH2O   
            Reaction(ARRPLUS_ade(2.0e-12, 1.2038e0, -9.0435e-4), [IHOO4, MO2], [MACR, HO2, CH2O], [1, 1], [1, 2, 2])   #2019/11/06; Bates2019; KHB  #  IHOO4 + MO2 = MACR + 2HO2 + 2CH2O   
            Reaction(GC_NIT(2.7e-12, 3.50e2, 1.19e0,  6.0e0, 1.1644e0, 7.05e-4), [IHOO1, NO], [IHN2], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHOO1 + NO = IHN2   
            Reaction(GC_ALK(2.7e-12, 3.50e2, 1.19e0,  6.0e0, 1.1644e0, 7.05e-4), [IHOO1, NO], [NO2, MVK, HO2, CH2O], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  IHOO1 + NO = NO2 + MVK + HO2 + CH2O   
            Reaction(GC_NIT(2.7e-12, 3.50e2, 1.421e0, 6.0e0, -0.1644e0, -7.05e-4), [IHOO1, NO], [IHN4], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHOO1 + NO = IHN4   
            Reaction(GC_NIT(2.7e-12, 3.50e2, 1.297e0, 6.0e0, 1.2038e0, 9.04e-4), [IHOO4, NO], [IHN3], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHOO4 + NO = IHN3   
            Reaction(GC_ALK(2.7e-12, 3.50e2, 1.297e0, 6.0e0, 1.2038e0, 9.04e-4), [IHOO4, NO], [NO2, MACR, HO2, CH2O], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  IHOO4 + NO = NO2 + MACR + HO2 + CH2O   
            Reaction(GC_NIT(2.7e-12, 3.50e2, 1.421e0, 6.0e0, -0.2038e0, -9.04e-4), [IHOO4, NO], [IHN1], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHOO4 + NO = IHN1   
            Reaction(GCARR_ac(3.00e-12, 650.0e0), [IDC, OH], [CO, HO2, MVKPC], [1, 1], [1, 1, 1])   #2019/11/06; Bates2019; KHB  #  IDC + OH = CO + HO2 + MVKPC   
            Reaction(GCARR_ac(1.59e+13, -10000.0e0), [IHPOO1], [ICPDH, IDHPE, OH], [1], [0.176, 0.824, 1])   #2019/11/06; Bates2019; KHB  #  IHPOO1 = 0.176ICPDH + 0.824IDHPE + OH   
            Reaction(GC_NIT(2.7e-12, 3.50e2, 2.1e0, 9.0e0, 1.0e0, 0.0e0), [IHPOO1, NO], [ITHN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHPOO1 + NO = ITHN   
            Reaction(GCARR_ac(2.91e+13, -10000.0e0), [IHPOO2], [ICPDH, IDHPE, OH], [1], [0.548, 0.452, 1])   #2019/11/06; Bates2019; KHB  #  IHPOO2 = 0.548ICPDH + 0.452IDHPE + OH   
            Reaction(GC_NIT(2.7e-12, 3.50e2, 2.315e0, 9.0e0, 1.0e0, 0.0e0), [IHPOO2, NO], [ITHN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHPOO2 + NO = ITHN   
            Reaction(GCARR_ac(1.875e+13, -10000.0e0), [IHPOO3], [IDHPE], [1], [1])   #2019/11/06; Bates2019; KHB  #  IHPOO3 = IDHPE   
            Reaction(GC_ALK(2.7e-12, 3.50e2, 3.079e0, 9.0e0, 1.0e0, 0.0e0), [IHPOO3, NO], [GLYC, HAC, NO2, OH], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  IHPOO3 + NO = GLYC + HAC + NO2 + OH   
            Reaction(GC_NIT(2.7e-12, 3.50e2, 3.079e0, 9.0e0, 1.0e0, 0.0e0), [IHPOO3, NO], [ITHN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHPOO3 + NO = ITHN   
            Reaction(GCARR_ac(1.05e-11, -400.0e0), [IEPOXA, OH], [ICHE, HO2], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  IEPOXA + OH = ICHE + HO2   
            Reaction(GC_EPO_a(5.82e-11, -4.00e2, 1.14e-20), [IEPOXA, OH], [IEPOXAOO, IEPOXBOO], [1, 1], [0.67, 0.33])   #2019/11/06; Bates2019; KHB  #  IEPOXA + OH = 0.67IEPOXAOO + 0.33IEPOXBOO   
            Reaction(GCARR_ac(8.25e-12, -400.0e0), [IEPOXB, OH], [ICHE, HO2], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  IEPOXB + OH = ICHE + HO2   
            Reaction(GC_EPO_a(3.75e-11, -4.00e2, 8.91e-21), [IEPOXB, OH], [IEPOXAOO, IEPOXBOO], [1, 1], [0.81, 0.19])   #2019/11/06; Bates2019; KHB  #  IEPOXB + OH = 0.81IEPOXAOO + 0.19IEPOXBOO   
            Reaction(GCARR_ac(1.875e+13, -10000.0e0), [IEPOXAOO], [IDCHP, HO2], [1], [1, 1])   #2019/11/06; Bates2019; KHB  #  IEPOXAOO = IDCHP + HO2   
            Reaction(GCARR_ac(1.0e+7, -5000.0e0), [IEPOXAOO], [OH, CO, MVKDH], [1], [1, 1, 1])   #2019/11/06; Bates2019; KHB  #  IEPOXAOO = OH + CO + MVKDH   
            Reaction(GC_NIT(2.7e-12, 3.50e2, 13.098e0, 8.0e0, 1.0e0, 0.0e0), [IEPOXAOO, NO], [ITCN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IEPOXAOO + NO = ITCN   
            Reaction(GCARR_ac(1.875e+13, -10000.0e0), [IEPOXBOO], [IDCHP, HO2], [1], [1, 1])   #2019/11/06; Bates2019; KHB  #  IEPOXBOO = IDCHP + HO2   
            Reaction(GCARR_ac(1.0e+7, -5000.0e0), [IEPOXBOO], [CO, OH, MCRDH], [1], [1, 1, 1])   #2019/11/06; Bates2019; KHB  #  IEPOXBOO = CO + OH + MCRDH   
            Reaction(GC_NIT(2.7e-12, 3.50e2, 16.463e0, 8.0e0, 1.0e0, 0.0e0), [IEPOXBOO, NO], [ITCN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IEPOXBOO + NO = ITCN   
            Reaction(GC_NIT(2.7e-12, 3.50e2, 13.098e0, 8.0e0, 1.0e0, 0.0e0), [ICHOO, NO], [ITCN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  ICHOO + NO = ITCN   
            Reaction(GCARR_ac(1.875e+13, -10000.0e0), [ICHOO], [HO2, CO, HAC, OH], [1], [1, 2.000, 1, 1])   #2019/11/06; Bates2019; KHB  #  ICHOO = HO2 + 2.000CO + HAC + OH   
            Reaction(GCARR_ac(2.70e-12, 350.0e0), [HPALD1OO, NO], [NO2, OH, CO2, MVK], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  HPALD1OO + NO = NO2 + OH + CO2 + MVK   
            Reaction(GCARR_ac(2.38e-13, 1300.0e0), [HPALD1OO, HO2], [OH, OH, CO2, MVK], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  HPALD1OO + HO2 = OH + OH + CO2 + MVK   
            Reaction(GCARR_ac(2.70e-12, 350.0e0), [HPALD2OO, NO], [NO2, OH, CO2, MACR], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  HPALD2OO + NO = NO2 + OH + CO2 + MACR   
            Reaction(GCARR_ac(2.38e-13, 1300.0e0), [HPALD2OO, HO2], [OH, OH, CO2, MACR], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  HPALD2OO + HO2 = OH + OH + CO2 + MACR   
            Reaction(GCARR_ac(7.14e-12, 390.0e0), [IHN2, OH], [ISOPNOO1], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHN2 + OH = ISOPNOO1   
            Reaction(GC_EPO_a(6.30e-12, 390.0e0, 1.62e-19), [IHN2, OH], [IEPOXA, IEPOXB, NO2], [1, 1], [0.67, 0.33, 1])   #2019/11/06; Bates2019; KHB  #  IHN2 + OH = 0.67IEPOXA + 0.33IEPOXB + NO2   
            Reaction(GCARR_ac(1.02e-11, 390.0e0), [IHN3, OH], [ISOPNOO2], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHN3 + OH = ISOPNOO2   
            Reaction(GC_EPO_a(1.05e-11, 390.0e0, 2.49e-19), [IHN3, OH], [IEPOXA, IEPOXB, NO2], [1, 1], [0.67, 0.33, 1])   #2019/11/06; Bates2019; KHB  #  IHN3 + OH = 0.67IEPOXA + 0.33IEPOXB + NO2   
            Reaction(GC_EPO_a(1.55e-11, 390.0e0, 2.715e-19), [IHN1, OH], [IEPOXD, NO2], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  IHN1 + OH = IEPOXD + NO2   
            Reaction(GCARR_ac(2.04e-11, 390.0e0), [IHN1, OH], [IDHNDOO1], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHN1 + OH = IDHNDOO1   
            Reaction(GC_EPO_a(9.52e-12, 390.0e0, 2.715e-19), [IHN4, OH], [IEPOXD, NO2], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  IHN4 + OH = IEPOXD + NO2   
            Reaction(GCARR_ac(2.95e-11, 390.0e0), [IHN4, OH], [IDHNDOO2], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHN4 + OH = IDHNDOO2   
            Reaction(GCARR_ac(1.875e+13, -10000.0e0), [ISOPNOO1], [ITCN, HO2], [1], [1, 1])   #2019/11/06; Bates2019; KHB  #  ISOPNOO1 = ITCN + HO2   
            Reaction(GC_NIT(2.7e-12, 350.0e0, 6.32e0, 11.0e0, 1.0e0, 0.0e0), [ISOPNOO1, NO], [IDN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  ISOPNOO1 + NO = IDN   
            Reaction(GCARR_ac(1.875e+13, -10000.0e0), [ISOPNOO2], [ITCN, HO2], [1], [1, 1])   #2019/11/06; Bates2019; KHB  #  ISOPNOO2 = ITCN + HO2   
            Reaction(GC_ALK(2.7e-12, 350.0e0, 7.941e0, 11.0e0, 1.0e0, 0.0e0), [ISOPNOO2, NO], [MVKN, CH2O, HO2, NO2], [1, 1], [1, 1, 1, 1])   #2019/11/06; Bates2019; KHB  #  ISOPNOO2 + NO = MVKN + CH2O + HO2 + NO2   
            Reaction(GC_NIT(2.7e-12, 350.0e0, 7.941e0, 11.0e0, 1.0e0, 0.0e0), [ISOPNOO2, NO], [IDN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  ISOPNOO2 + NO = IDN   
            Reaction(GCARR_ac(1.256e+13, -10000.0e0), [IDHNDOO1], [ITCN, HO2], [1], [1, 1])   #2019/11/06; Bates2019; KHB  #  IDHNDOO1 = ITCN + HO2   
            Reaction(GCARR_ac(5.092e+12, -10000.0e0), [IDHNDOO2], [ITCN, HO2], [1], [1, 1])   #2019/11/06; Bates2019; KHB  #  IDHNDOO2 = ITCN + HO2   
            Reaction(GC_NIT(2.7e-12, 350.0e0, 4.712e0, 11.0e0, 1.0e0, 0.0e0), [IDHNDOO1, NO], [IDN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IDHNDOO1 + NO = IDN   
            Reaction(GC_NIT(2.7e-12, 350.0e0, 2.258e0, 11.0e0, 1.0e0, 0.0e0), [IDHNDOO2, NO], [IDN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IDHNDOO2 + NO = IDN   
            Reaction(GC_NIT(2.7e-12, 350.0e0, 1.851e0, 11.0e0, 1.0e0, 0.0e0), [IDHNBOO, NO], [IDN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IDHNBOO + NO = IDN   
            Reaction(GCARR_ac(2.47e-13, 1300.0e0), [INO2D, HO2], [INPD], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  INO2D + HO2 = INPD   
            Reaction(GC_NIT(2.7e-12, 350.0e0, 12.915e0, 9.0e0, 1.0e0, 0.0e0), [INO2B, NO], [IDN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  INO2B + NO = IDN   
            Reaction(GC_NIT(2.7e-12, 350.0e0, 1.412e0, 9.0e0, 1.0e0, 0.0e0), [INO2D, NO], [IDN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  INO2D + NO = IDN   
            Reaction(GCARR_ac(2.50e-14, -300.0e0), [INA, O2], [ICN, HO2], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  INA + O2 = ICN + HO2   
            Reaction(GCARR_ac(1.00e+20, -10000.0e0), [INA], [IDHNBOO], [1], [1])   #2019/11/06; Bates2019; KHB  #  INA = IDHNBOO   
            Reaction(GCARR_ac(5.88e-12, 390.0e0), [INPB, OH], [IHPNBOO, IDHNBOO], [1, 1], [0.670, 0.33])   #2019/11/06; Bates2019; KHB  #  INPB + OH = 0.670IHPNBOO + 0.33IDHNBOO   
            Reaction(GCARR_ac(1.61e-11, 390.0e0), [INPD, OH], [IHPNDOO], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  INPD + OH = IHPNDOO   
            Reaction(GC_EPO_a(4.471e-12, 390.0e0, 2.28e-20), [INPB, OH], [OH, ITHN], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  INPB + OH = OH + ITHN   
            Reaction(GC_EPO_a(8.77e-12,  390.0e0, 2.185e-20), [INPD, OH], [OH, ITHN], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  INPD + OH = OH + ITHN   
            Reaction(GC_EPO_a(1.493e-11, 390.0e0, 2.715e-19), [INPD, OH], [NO2, ICHE], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  INPD + OH = NO2 + ICHE   
            Reaction(GCARR_ac(2.278e-12, 200.0e0), [INPB, OH], [INO2B], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  INPB + OH = INO2B   
            Reaction(GCARR_ac(3.40e-12, 200.0e0), [INPD, OH], [INO2D], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  INPD + OH = INO2D   
            Reaction(GCARR_ac(7.50e-12, 20.0e0), [INPD, OH], [ICN, OH], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  INPD + OH = ICN + OH   
            Reaction(GCARR_ac(6.55e+12, -10000.0e0), [IHPNDOO], [OH, ITCN], [1], [1, 1])   #2019/11/06; Bates2019; KHB  #  IHPNDOO = OH + ITCN   
            Reaction(GCARR_ac(8.72e+12, -10000.0e0), [IHPNBOO], [OH, ITCN, ITHN], [1], [1, 0.5, 0.5])   #2019/11/06; Bates2019; KHB  #  IHPNBOO = OH + 0.5ITCN + 0.5ITHN   
            Reaction(GC_NIT(2.7e-12, 350.0e0, 6.092e0, 12.0e0, 1.0e0, 0.0e0), [IHPNBOO, NO], [IDN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHPNBOO + NO = IDN   
            Reaction(GC_NIT(2.7e-12, 350.0e0, 4.383e0, 12.0e0, 1.0e0, 0.0e0), [IHPNDOO, NO], [IDN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  IHPNDOO + NO = IDN   
            Reaction(GC_EPO_a(2.97e-12, 390.0e0, 2.715e-19), [ICN, OH], [NO2, ICHE], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  ICN + OH = NO2 + ICHE   
            Reaction(GCARR_ac(2.60e-12, 610.0e0), [MVK, OH], [MVKOHOO], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  MVK + OH = MVKOHOO   
            Reaction(GCARR_ac(2.70e-12, 470.0e0), [MACR, OH], [MACR1OO], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  MACR + OH = MACR1OO   
            Reaction(5.77e-11, [MVKHP, OH], [MVKHC, MVKHCB, OH], [1, 1], [0.53, 0.47, 1])   #2019/11/06; Bates2019; KHB  #  MVKHP + OH = 0.53MVKHC + 0.47MVKHCB + OH   
            Reaction(GCARR_ac(1.39e-11, 380.0e0), [MCRHN, OH], [MACRNO2], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  MCRHN + OH = MACRNO2   
            Reaction(GCARR_ac(2.70e-12, 350.0e0), [C4HVP1, NO], [NO2, MVKOHOO], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  C4HVP1 + NO = NO2 + MVKOHOO   
            Reaction(GCARR_ac(1.93e-13, 1300.0e0), [C4HVP1, HO2], [OH, MVKOHOO], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  C4HVP1 + HO2 = OH + MVKOHOO   
            Reaction(9.00e-12, [C4HVP1, NO2], [MVKN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  C4HVP1 + NO2 = MVKN   
            Reaction(GCARR_ac(2.70e-12, 350.0e0), [C4HVP2, NO], [NO2, MCROHOO], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  C4HVP2 + NO = NO2 + MCROHOO   
            Reaction(GCARR_ac(1.93e-13, 1300.0e0), [C4HVP2, HO2], [OH, MCROHOO], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  C4HVP2 + HO2 = OH + MCROHOO   
            Reaction(9.00e-12, [C4HVP2, NO2], [MCRHN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  C4HVP2 + NO2 = MCRHN   
            Reaction(GCARR_ac(5.00e-12, 470.0e0), [MVKPC, OH], [OH, CO, MGLY], [1, 1], [1, 1, 1])   #2019/11/06; Bates2019; KHB  #  MVKPC + OH = OH + CO + MGLY   
            Reaction(GCARR_ac(8.70e-12, 70.0e0), [MVKDH, OH], [MVKHCB, MVKHC, HO2], [1, 1], [0.4, 0.6, 1])   #2019/11/06; Bates2019; KHB  #  MVKDH + OH = 0.4MVKHCB + 0.6MVKHC + HO2   
            Reaction(GCARR_ac(5.00e-12, 470.0e0), [MVKHCB, OH], [OH, MGLY], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  MVKHCB + OH = OH + MGLY   
            Reaction(GCARR_ac(2.00e-12, 70.0e0), [MVKHC, OH], [CO, HO2, MCO3], [1, 1], [2, 1, 1])   #2019/11/06; Bates2019; KHB  #  MVKHC + OH = 2CO + HO2 + MCO3   
            Reaction(GC_NIT(2.7e-12, 350.0e0, 4.573e0, 6.0e0, 1.0e0, 0.0e0), [MVKOHOO, NO], [MVKN], [1, 1], [0.438])   #2019/11/06; Bates2019; KHB  #  MVKOHOO + NO = 0.438MVKN   
            Reaction(GCARR_ac(2.90e+7, -5297.0e0), [MCROHOO], [HAC, CO, OH], [1], [1, 1, 1])   #2019/11/06; Bates2019; KHB  #  MCROHOO = HAC + CO + OH   
            Reaction(GC_NIT(2.7e-12, 350.0e0, 2.985e0, 6.0e0, 1.0e0, 0.0e0), [MCROHOO, NO], [MCRHN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  MCROHOO + NO = MCRHN   
            Reaction(GC_PAN_acac(2.591e-28, -6.87e0, 1.125e-11, -1.105e0, 0.3e0), [MACR1OO, NO2], [MPAN], [1, 1], [1])   #2019/11/06; Bates2019; KHB  #  MACR1OO + NO2 = MPAN   
            Reaction(GCARR_ac(7.50e-12, 290.0e0), [MACRNO2, NO], [HAC, NO2, CO2], [1, 1], [1, 2, 1])   #2019/11/06; Bates2019; KHB  #  MACRNO2 + NO = HAC + 2NO2 + CO2   
            Reaction(GC_PAN_acac(2.591e-28, -6.87e0, 1.125e-11, -1.105e0, 0.3e0), [MACRNO2, NO2], [MPAN, NO2], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  MACRNO2 + NO2 = MPAN + NO2   
            Reaction(4.00e-12, [MACRNO2, NO3], [HAC, NO2, CO2], [1, 1], [1, 2, 1])   #2019/11/06; Bates2019; KHB  #  MACRNO2 + NO3 = HAC + 2NO2 + CO2   
            Reaction(GCARR_ac(1.58e+16, -13500.0e0), [MPAN], [MACR1OO, NO2], [1], [1, 1])   #2019/11/06; Bates2019; KHB  #  MPAN = MACR1OO + NO2   
            Reaction(3.00e-12, [IDHDP, OH], [OH, ICPDH, IDHPE], [1, 1], [1, 0.333, 0.667])   #2019/11/06; Bates2019; KHB  #  IDHDP + OH = OH + 0.333ICPDH + 0.667IDHPE   
            Reaction(GCARR_ac(1.40e-12, -1860.0e0), [ETHLN, NO3], [HNO3, NO2, MCO3], [1, 1], [1, 1, 1])   #2019/11/06; Bates2019; KHB  #  ETHLN + NO3 = HNO3 + NO2 + MCO3   
            Reaction(8.00e-13, [PYAC, OH], [MCO3, CO2], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  PYAC + OH = MCO3 + CO2   
            Reaction(GCARR_ac(1.55e-12, 340.0e0), [HPETHNL, OH], [CO, OH, CH2O], [1, 1], [1, 1, 1])   #2019/11/06; Bates2019; KHB  #  HPETHNL + OH = CO + OH + CH2O   
            Reaction(2.91e-11, [HPETHNL, OH], [GLYX, OH], [1, 1], [1, 1])   #2019/11/06; Bates2019; KHB  #  HPETHNL + OH = GLYX + OH   
            Reaction(GCARR_ac(1.56e-11, 117.0e0), [NAP, OH], [NRO2, OH], [1, 1], [1, 1])   #2013/08/12; Pye2010; HOTP  #  NAP + OH = NRO2 + OH   
            Reaction(GCARR_ac(1.40e-12, 700.0e0), [NRO2, HO2], [LNRO2H, HO2], [1, 1], [1, 1])   #2013/08/12; Pye2010; HOTP  #  NRO2 + HO2 = LNRO2H + HO2   
            Reaction(GCARR_ac(2.60e-12, 350.0e0), [NRO2, NO], [LNRO2N, NO], [1, 1], [1, 1])   #2013/08/12; Pye2010; HOTP  #  NRO2 + NO = LNRO2N + NO   
            Reaction(GCARR_abc(9.10e-15, 0.0e0, -2580.0e0), [C2H4, O3], [CH2O, CH2OO], [1, 1], [1, 1])   #2021/09/22; Kwon2020; KHB,MSL  #  C2H4 + O3 = CH2O + CH2OO   
            Reaction(GCJPLPR_abab(1.10e-28, 3.5e+00, 8.4e-12, 1.75e0, 0.5e0), [C2H4, OH], [ETOO], [1, 1], [1])   #2021/09/22; Kwon2020; KHB,MSL  #  C2H4 + OH = ETOO   
            Reaction(GCARR_abc(1.53e-13, 0.0e0, 1300.0e0), [ETOO, HO2], [ETHP], [1, 1], [1])   #2021/09/22; Kwon2020; KHB,MSL  #  ETOO + HO2 = ETHP   
            Reaction(2.3e-12, [ETOO, NO3], [ETO, NO2], [1, 1], [1, 1])   #2021/09/22; Kwon2020; KHB,MSL  #  ETOO + NO3 = ETO + NO2   
            Reaction(GCARR_abc(9.5e-13, 0.0e0, -5988.0e0), [ETO], [HO2, CH2O], [1], [1, 2.000])   #2021/09/22; Kwon2020; KHB,MSL  #  ETO = HO2 + 2.000CH2O   
            Reaction(GCARR_abc(2.5e-14, 0.0e0, -300.0e0), [ETO, O2], [GLYC, HO2], [1, 1], [1, 1])   #2021/09/22; Kwon2020; KHB,MSL  #  ETO + O2 = GLYC + HO2   
            Reaction(8.40e-13, [ETHN, OH], [GLYC, NO2], [1, 1], [1, 1])   #2021/09/22; Kwon2020; KHB,MSL  #  ETHN + OH = GLYC + NO2   
            Reaction(GCARR_abc(1.90e-12, 0.0e+00, 190.0e0), [ETHP, OH], [ETOO], [1, 1], [1])   #2021/09/22; Kwon2020; KHB,MSL  #  ETHP + OH = ETOO   
            Reaction(1.38e-11, [ETHP, OH], [OH, GLYC], [1, 1], [1, 1])   #2021/09/22; Kwon2020; KHB,MSL  #  ETHP + OH = OH + GLYC   
            Reaction(2.91e-13 * exp( 1300.0/ 288.15) * 0.82, [AROMRO2, HO2], [OH, HO2], [1, 1], [1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  AROMRO2 + HO2 = OH + HO2   
            Reaction(GCARR_abc(2.60e-12, 0.0e+00, 365.0e0), [AROMRO2, NO], [NO2, HO2], [1, 1], [1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  AROMRO2 + NO = NO2 + HO2   
            Reaction(2.30e-12, [AROMRO2, NO3], [NO2, HO2], [1, 1], [1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  AROMRO2 + NO3 = NO2 + HO2   
            Reaction(GCARR_abc(1.70e-14, 0.0e0, 220.0e0), [AROMRO2, MO2], [CH2O, HO2, HO2], [1, 1], [1, 1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  AROMRO2 + MO2 = CH2O + HO2 + HO2   
            Reaction(GCARR_abc(4.20e-14, 0.0e0, 220.0e0), [AROMRO2, MCO3], [MO2, HO2], [1, 1], [1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  AROMRO2 + MCO3 = MO2 + HO2   
            Reaction(9.2e-18, [MCT, O3], [GLYC, HO2, OH, AROMP4], [1, 1], [1, 1, 1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  MCT + O3 = GLYC + HO2 + OH + AROMP4   
            Reaction(GCARR_abc(5.90e-12, 0.0e0, 225.0e0), [BALD, OH], [BZCO3], [1, 1], [1])   #2021/09/29; Bates2021b; KHB,MSL  #  BALD + OH = BZCO3   
            Reaction(2.4e-15, [BALD, NO3], [BZCO3, HNO3], [1, 1], [1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  BALD + NO3 = BZCO3 + HNO3   
            Reaction(GCARR_abc(7.50e-12, 0.0e0, 290.0e0), [BZCO3, NO], [NO2, CO2, BENZO2], [1, 1], [1, 1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  BZCO3 + NO = NO2 + CO2 + BENZO2   
            Reaction(GC_PAN_acac(3.28e-28, -6.87e0, 1.125e-11, -1.105e0, 0.3e0), [BZCO3, NO2], [BZPAN], [1, 1], [1])   #2021/09/29; Bates2021b; KHB,MSL  #  BZCO3 + NO2 = BZPAN   
            Reaction(4.66e-12, [BZCO3H, OH], [BZCO3], [1, 1], [1])   #2021/09/29; Bates2021b; KHB,MSL  #  BZCO3H + OH = BZCO3   
            Reaction(GC_PAN_abab(1.10e-5, -10100.0e0, 1.90e+17, -14100.0e0, 0.3e0)*0.67e0, [BZPAN], [BZCO3, NO2], [1], [1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  BZPAN = BZCO3 + NO2   
            Reaction(1.06e-12, [BZPAN, OH], [BENZP, CO2, NO2], [1, 1], [1, 1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  BZPAN + OH = BENZP + CO2 + NO2   
            Reaction(7.00e-12, [BENZO2, NO2], [BENZO, NO3], [1, 1], [1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  BENZO2 + NO2 = BENZO + NO3   
            Reaction(GCARR_abc(2.670e-12, 0.0e0, 365.0e0), [BENZO2, NO], [BENZO, NO2], [1, 1], [1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  BENZO2 + NO = BENZO + NO2   
            Reaction(2.30e-12, [BENZO2, NO3], [BENZO, NO2], [1, 1], [1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  BENZO2 + NO3 = BENZO + NO2   
            Reaction(GCARR_abc(2.24e-13, 0.0e0, 1300.0e0), [BENZO2, HO2], [BENZP], [1, 1], [1])   #2021/09/29; Bates2021b; KHB,MSL  #  BENZO2 + HO2 = BENZP   
            Reaction(3.60e-12, [BENZP, OH], [BENZO2], [1, 1], [1])   #2021/09/29; Bates2021b; KHB,MSL  #  BENZP + OH = BENZO2   
            Reaction(2.86e-13, [BENZO, O3], [BENZO2], [1, 1], [1])   #2021/09/29; Bates2021b; KHB,MSL  #  BENZO + O3 = BENZO2   
            Reaction(2.08e-12, [BENZO, NO2], [NPHEN], [1, 1], [1])   #2021/09/29; Bates2021b; KHB,MSL  #  BENZO + NO2 = NPHEN   
            Reaction(3.47e-12, [NPHEN, OH], [R4N1, AROMP4, NO2], [1, 1], [0.5, 1, 0.5])   #2021/09/29; Bates2021b; KHB,MSL  #  NPHEN + OH = 0.5R4N1 + AROMP4 + 0.5NO2   
            Reaction(GCARR_abc(2.670e-13, 0.0e0, 365.0e0), [BENZO2, MO2], [BENZO, HO2, CH2O], [1, 1], [1, 1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  BENZO2 + MO2 = BENZO + HO2 + CH2O   
            Reaction(GCARR_abc(2.670e-12, 0.0e0, 365.0e0), [BZCO3, MO2], [BENZO2, CO2, HO2, CH2O], [1, 1], [1, 1, 1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  BZCO3 + MO2 = BENZO2 + CO2 + HO2 + CH2O   
            Reaction(1.5e-3, [AROMP4], [HO2, GLYX, RCHO], [1], [0.2, 0.2, 1.2])   #2021/09/29; Bates2021b; KHB,MSL  #  AROMP4 = 0.2HO2 + 0.2GLYX + 1.2RCHO   
            Reaction(GCARR_abc(1.40e-12, 0.0e0, 700.0e0), [BRO2, HO2], [HO2, LBRO2H], [1, 1], [1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  BRO2 + HO2 = HO2 + LBRO2H   
            Reaction(GCARR_abc(2.60e-12, 0.0e0, 350.0e0), [BRO2, NO], [NO, LBRO2N], [1, 1], [1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  BRO2 + NO  = NO + LBRO2N    
            Reaction(GCARR_abc(1.40e-12, 0.0e0, 700.0e0), [TRO2, HO2], [HO2, LTRO2H], [1, 1], [1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  TRO2 + HO2 = HO2 + LTRO2H   
            Reaction(GCARR_abc(2.60e-12, 0.0e0, 350.0e0), [TRO2, NO], [NO, LTRO2N], [1, 1], [1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  TRO2 + NO  = NO + LTRO2N    
            Reaction(GCARR_abc(1.40e-12, 0.0e0, 700.0e0), [XRO2, HO2], [HO2, LXRO2H], [1, 1], [1, 1])   #2021/09/29; Bates2021b; KHB,MSL  #  XRO2 + HO2 = HO2 + LXRO2H   
            Reaction(GCARR_abc(2.60e-12, 0.0e0, 350.0e0), [XRO2, NO], [NO, LXRO2N], [1, 1], [1, 1])  #  XRO2 + NO  = NO + LXRO2N    

    
        ]
    
        @named fullchem = ReactionSystem(fc,t)
    end
end

rs = fullchem()


#Test, plot the changing of O3, NO2 and ISOP
@unpack O3, NO2, ISOP = rs
rs_solved = solve(ODEProblem(rs, [], (0,1), [], combinatoric_ratelaws=false),Tsit5())
plt.plot(rs_solved.t, rs_solved[O3])
plt.plot(rs_solved.t, rs_solved[NO2])
plt.plot(rs_solved.t, rs_solved[ISOP])
plt.show()
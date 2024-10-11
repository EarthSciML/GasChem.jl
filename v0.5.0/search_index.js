var documenterSearchIndex = {"docs":
[{"location":"superfast/#SuperFast-Gas-Phase-Atmospheric-Chemical-Mechanism","page":"SuperFast","title":"SuperFast Gas-Phase Atmospheric Chemical Mechanism","text":"","category":"section"},{"location":"superfast/","page":"SuperFast","title":"SuperFast","text":"The Super Fast Chemical Mechanism is one of the simplest representations of atmospheric chemistry. It can efficiently simulate background tropheric ozone chemistry and perform well for those species included in the mechanism. The chemical equations used is included in the supporting table S2 of the paper, \"Evaluating simplified chemical mechanisms within present-day simulations of the Community Earth System Model version 1.2 with CAM4 (CESM1.2 CAM-chem): MOZART-4 vs. Reduced Hydrocarbon vs. Super-Fast chemistry\" (2018), Benjamin Brown-Steiner, Noelle E. Selin, Ronald G. Prinn, Simone Tilmes, Louisa Emmons, Jean-François Lamarque, and Philip Cameron-Smith.","category":"page"},{"location":"superfast/#Illustrative-Example","page":"SuperFast","title":"Illustrative Example","text":"","category":"section"},{"location":"superfast/","page":"SuperFast","title":"SuperFast","text":"Here is a simple example of initializing the SuperFast model and running a simulation. First, we can look at the reaction equations:","category":"page"},{"location":"superfast/","page":"SuperFast","title":"SuperFast","text":"using GasChem, EarthSciMLBase, ModelingToolkit, Unitful, DifferentialEquations\nusing Catalyst\n\n@parameters t [unit = u\"s\", description = \"Time\"]\nmodel = SuperFast(t)\n\nmodel.rxn_sys","category":"page"},{"location":"superfast/","page":"SuperFast","title":"SuperFast","text":"We can also look at them as a graph:","category":"page"},{"location":"superfast/","page":"SuperFast","title":"SuperFast","text":"Graph(model.rxn_sys)","category":"page"},{"location":"superfast/","page":"SuperFast","title":"SuperFast","text":"Before running any simulations with the model, we need to convert it into a system of differential equations. We can solve it using the default values for variables and parameters. However, by using the @unpack command, we can assign new values to specific variables and parameters, allowing for simulations under varied conditions.","category":"page"},{"location":"superfast/","page":"SuperFast","title":"SuperFast","text":"We can visualize the differential equation version of the system as follows:","category":"page"},{"location":"superfast/","page":"SuperFast","title":"SuperFast","text":"sys = structural_simplify(get_mtk(model))\n@unpack O3, T = sys\ntspan = (0.0, 3600*24)\nu0 = [O3 => 15] # Change the initial concentration of O₃ to 15 ppb\np0 = [T => 293] # temperature = 293K\nprob = ODEProblem(sys, u0, tspan, p0)\n\nequations(sys)","category":"page"},{"location":"superfast/","page":"SuperFast","title":"SuperFast","text":"We can finally solve the system and plot the result as","category":"page"},{"location":"superfast/","page":"SuperFast","title":"SuperFast","text":"sol = solve(prob,AutoTsit5(Rosenbrock23()), saveat=10.0)\n\nusing Plots\nplot(sol, ylim = (0,50), xlabel = \"Time\", ylabel = \"Concentration (ppb)\", legend=:outerright)","category":"page"},{"location":"superfast/#Variables-and-parameters","page":"SuperFast","title":"Variables and parameters","text":"","category":"section"},{"location":"superfast/","page":"SuperFast","title":"SuperFast","text":"The species included in the superfast model are: O₃, OH, HO₂, O₂, NO, NO₂, CH₄, CH₃O₂, H₂O, CH₂O, CO, CH₃OOH, CH₃O, DMS, SO₂, ISOP, O₁d, H₂O₂.","category":"page"},{"location":"superfast/","page":"SuperFast","title":"SuperFast","text":"The parameters in the model that are not constant are the photolysis reaction rates jO31D, j2OH, jH2O2, jNO2, jCH2Oa, jCH3OOH and temperature T","category":"page"},{"location":"superfast/","page":"SuperFast","title":"SuperFast","text":"states(sys) # Give you the variables in the system\nparameters(sys) # Give you the parameters in the system","category":"page"},{"location":"superfast/","page":"SuperFast","title":"SuperFast","text":"Let's run some simulation with different values for parameter T.","category":"page"},{"location":"superfast/","page":"SuperFast","title":"SuperFast","text":"p1 = [T => 273]\np2 = [T => 298]\nsol1 = solve(ODEProblem(sys, [], tspan, p1),AutoTsit5(Rosenbrock23()), saveat=10.0)\nsol2 = solve(ODEProblem(sys, [], tspan, p2),AutoTsit5(Rosenbrock23()), saveat=10.0)\n\nplot([sol1[O3],sol2[O3]], label = [\"T=273K\" \"T=298K\"], title = \"Change of O3 concentration at different temperatures\", xlabel=\"Time (second)\", ylabel=\"concentration (ppb)\")","category":"page"},{"location":"api/#Index","page":"API","title":"Index","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"","category":"page"},{"location":"api/#API-Documentation","page":"API","title":"API Documentation","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [GasChem]","category":"page"},{"location":"api/#GasChem.FastJX","page":"API","title":"GasChem.FastJX","text":"Description: This is a box model used to calculate the photolysis reaction rate constant using the Fast-JX scheme  (Neu, J. L., Prather, M. J., and Penner, J. E. (2007), Global atmospheric chemistry: Integrating over fractional cloud cover, J. Geophys. Res., 112, D11306, doi:10.1029/2006JD008007.)\n\nBuild Fast-JX model\n\nExample\n\n    @parameters t\n    fj = FastJX(t)\n\n\n\n\n\n","category":"type"},{"location":"api/#GasChem.GEOSChemGasPhase","page":"API","title":"GasChem.GEOSChemGasPhase","text":"GEOS-Chem full-chem mechanism adapted from GEOS-Chem version 14.1.1\n\nAdapted from file https://github.com/geoschem/geos-chem/blob/4722f288e90291ba904222f4bbe4fc216d17c34a/KPP/fullchem/fullchem.eqn\nThe GEOS-Chem license applies: https://github.com/geoschem/geos-chem/blob/main/LICENSE.txt\n\n=============================================================================== REFERENCES (alphabetical order) ===============================================================================\n\nAtkinson2003: Atkinson and Arey, Chem. Rev., doi:10.1021/cr0206420, 2003.\nBates2014:    Bates et al., J. Phys. Chem A, 118, doi:10.1021/jp4107958, 2014.\nBates2019:    Bates and Jacob, Atmos. Chem. Phys., doi:10.5194/acp-19-9613-2019, 2019.\nBates2021a:   Bates et al, JGR, https://doi.org/10.1029/2020JD033439, 2021.\nBates2021b:   Bates et al, ACP, https://doi.org/10.5194/acp-2021-605, 2021.\nBrowne2011:   Browne et al., Atmos. Chem. Phys., doi:10.5194/acp-11-4209-2011, 2011.\nBrowne2014:   Browne et al., Atmos. Chem. Phys., doi:10.5194/acp-14-1225-2014, 2014.\nChen2017:     Chen et al., Geophys. Res. Lett., doi:10.1002/2017GL073812, 2017.\nCrounse2012:  Crounse et al., J. Phys. Chem. A, doi:10.1021/jp211560u, 2012.\nEastham2014:  Eastham et al., Atmos. Env., doi:10.1016/j.atmosenv.2014.02.001, 2014.\nFischer2014:  Fischer et al., Atmos. Chem. Phys., doi:10.5194/acp-14-2679-2014, 2014.\nFisher2016:   Fisher et al., Atmos. Chem. Phys., doi:10.5194/acp-16-5969-2016, 2016.\nFisher2018:   Fisher et al., J. Geophys. Res., doi:10.1029/2018JD029046, 2018.\nFry2014:      Fry et al. Environ. Sci. Technol., doi:10.1021/es502204x, 2014.\nGill2002:     Gill and Hites, J. Phys. Chem. A, doi:10.1021/jp013532, 2002.\nGoliff2013:   Goliff et al., Atmos. Environ., doi:10.1016/j.atmosenv.2012.11.038, 2013.\nJacobs2014:   Jacobs et al., Atmos. Chem. Phys., doi:10.5194/acp-14-8933-2014, 2014.\nJenkin2015:   Jenkin et al., Atmos. Chem. Phys., doi:10.5194/acp-15-11433-2015, 2015.\nKasibhatla2018: Kasibhatla et al., Atmos. Chem. Phys., doi:10.5194/acp-18-11185-2018, 2018\nIUPAC ROO19: https://iupac-aeris.ipsl.fr/htdocs/datasheets/pdf/ROO19CH3O2NO3.pdf\nJPL 10-6:     JPL Publication 10-6, https://jpldataeval.jpl.nasa.gov/previous_evaluations.html, 2011.\nJPL 15-10:    JPL Publication 15-10, https://jpldataeval.jpl.nasa.gov, 2015.\nKwon2020:     Kwon et al, Elementa, https://doi.org/10.1525/elementa.2021.00109, 2020.\nLee2014:      Lee et al., J. Phys. Chem. A, doi:10.1021/jp4107603, 2014.\nMarais2016:   Marais et al., Atmos. Chem. Phys, doi:10.5194/acp-16-1603-2016, 2016.\nMiller2017:   Miller et al., Atmos. Chem. Phys. Discuss., doi:10.5194/acp-2016-1042, 2017.\nMillet2015:   Millet et al., Atmos. Chem. Phys., doi:10.5194/acp-15-6283-2015, 2015.\n\nMoch et al, JGR, https, * Moch2020 # //doi.org/10.1029/2020JD032706, 2020.\n\nMuller2014:   Muller et al., Atmos. Chem. Phys., doi:10.5194/acp-14-2497-2014, 2014.\nParrella2012: Parrella et al. Atmos. Chem. Phys, doi:10.5194/acp-12-6723-2012, 2012.\nPaulot2009:   Paulot et al., Atmos. Chem. Phys., doi:10.5194/acp-9-1479-2009, 2009a and               Paulot et al., Science, doi:10.1126/science.1172910, 2009b.\nPeeters2010:  Peeters and Muller, Phys. Chem. Chem. Phys., doi:10.1039/C0CP00811G, 2010.\nPeeters2014:  Peeters et al., J. Phys. Chem. A, doi:10.1021/jp5033146, 2014.\nPye2010:      Pye et al., Atmos. Chem. Phys., doi:10.5194/acp-10-11261-2010, 2010.\nRoberts1992:  Roberts and Bertman, Int. J. Chem. Kinet., doi:10.1002/kin.550240307, 1992.\nSherwen2016b: Sherwen et al., Atmos. Chem. Phys., doi:10.5194/acp-16-12239-2016, 2016b.\nSherwen2017:  Sherwen et al., Faraday Discuss., doi:10.1039/C7FD00026J, 2017.\nStClair2016:  St. Clair et al., J. Phys. Chem. A, doi:10.1021/acs.jpca.5b065322016, 2016.\nTravis2016:   Travis et al., Atmos. Chem. Phys., doi:10.5194/acp-16-13561-2016, 2016.\nWolfe2012:    Wolfe et al., Phys. Chem. Chem. Phys., doi: 10.1039/C2CP40388A, 2012.\nXie2013:      Xie et al., Atmos. Chem. Phys., doi:10.5194/acp-13-8439-2013, 2013.\n\n\n\n\n\n","category":"type"},{"location":"api/#GasChem.SuperFast","page":"API","title":"GasChem.SuperFast","text":"SuperFast(t)\n\nThis atmospheric chemical system model is built based on the Super Fast Chemical Mechanism, which is one of the simplest representations of atmospheric chemistry. It can efficiently simulate background tropheric ozone chemistry and perform well for those species included in the mechanism. The chemical equations we used is included in the supporting table S2 of the paper:\n\n\"Evaluating simplified chemical mechanisms within present-day simulations of the Community Earth System Model version 1.2 with CAM4 (CESM1.2 CAM-chem): MOZART-4 vs. Reduced Hydrocarbon vs. Super-Fast chemistry\" (2018), Benjamin Brown-Steiner, Noelle E. Selin, Ronald G. Prinn, Simone Tilmes, Louisa Emmons, Jean-François Lamarque, and Philip Cameron-Smith.\n\nThe input of the function is Temperature, concentrations of all chemicals, and reaction rates of photolysis reactions \n\nExample\n\nusing GasChem, EarthSciMLBase, OrdinaryDiffEq, Plots\n@variables t [unit = u\"s\"]\nrs = SuperFast(t)\nsol = solve(ODEProblem(get_mtk(rs), [], (0,360), [], combinatoric_ratelaws=false), Tsit5())\nplot(sol)\n\nWe set combinatoric_ratelaws=false because we are modeling macroscopic rather than microscopic behavior.  See here  and here.\n\n\n\n\n\n","category":"type"},{"location":"api/#Base.:+-Tuple{GEOSChemGasPhase, FastJX}","page":"API","title":"Base.:+","text":"Compose GEOSChemGasPhase and FastJX models together.\n\nExample\n\nusing GasChem, EarthSciMLBase, DifferentialEquations\n\n@parameters t \ngf = GEOSChemGasPhase(t) + FastJX(t)\ntspan = (0.0, 3600*24*2)\nsys = structural_simplify(get_mtk(gf))\nsol = solve(ODEProblem(sys, [], tspan, []),Rosenbrock23(), saveat=10.0)\nusing Plots\nplot(sol,ylims=(0,20),xlabel=\"Time (second)\", ylabel=\"concentration (ppb)\",legend=:outertopright)\n\n\n\n\n\n","category":"method"},{"location":"api/#Base.:+-Tuple{SuperFast, FastJX}","page":"API","title":"Base.:+","text":"Compose superfast and fast-jx models together.\n\nExample\n\nusing GasChem, EarthSciMLBase\n\n@parameters t \nsf = SuperFast(t) + FastJX(t)\ntspan = (0.0, 3600*24*2)\nsol = solve(ODEProblem(structural_simplify(get_mtk(sf)), [], tspan, []),Tsit5(), saveat=10.0)\nusing Plots\nplot(sol,ylims=(0,20),xlabel=\"Time (second)\", ylabel=\"concentration (ppb)\",legend=:outertopright)\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.arr_3rdbody-NTuple{10, Any}","page":"API","title":"GasChem.arr_3rdbody","text":"Third body effect for pressure dependence of rate coefficients. a1, b1, c1 are the Arrhenius parameters for the lower-limit rate. a2, b2, c2 are the Arrhenius parameters for the upper-limit rate. fv         is the falloff curve paramter, (see ATKINSON ET. AL (1992)            J. Phys. Chem. Ref. Data 21, P. 1145). Usually fv = 0.6.\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.arrhenius-NTuple{5, Any}","page":"API","title":"GasChem.arrhenius","text":"Arrhenius equation:\n\n    k = a0 * exp( c0  T ) * (T300)^b0\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.arrplus-NTuple{7, Any}","page":"API","title":"GasChem.arrplus","text":"Modified Arrhenius law.\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.calc_flux-Union{Tuple{T}, Tuple{T, T}} where T","page":"API","title":"GasChem.calc_flux","text":"calculate actinic flux at the given cosine of the solar zenith angle csa maximium actinic flux max_actinic_flux\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.constant_k-Tuple{Any, Any}","page":"API","title":"GasChem.constant_k","text":"Function to create a constant rate coefficient\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.cos_solar_zenith_angle-Union{Tuple{T2}, Tuple{T}, Tuple{T, T2, T}} where {T, T2}","page":"API","title":"GasChem.cos_solar_zenith_angle","text":"cos_solar_zenith_angle(lat, t, long)\n\nThis function is to compute the cosine of the solar zenith angle, given the unixtime, latitude and longitude The input variables: lat=latitude(°), long=longitude(°), t=unixtime(s)     the cosine of the solar zenith angle (SZA) is given by:                                                                            .            cos(SZA) = sin(LAT)sin(DEC) + cos(LAT)cos(DEC)*cos(AHR)\n\n       where LAT = the latitude angle,\n             DEC = the solar declination angle,\n             AHR = the hour angle, all in radians.  All in radians\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.create_fjx_interp-Tuple{Vector{Float32}, Vector{StaticArraysCore.SVector{18, Float32}}}","page":"API","title":"GasChem.create_fjx_interp","text":"Create a vector of interpolators to interpolate the cross sections σ (TODO: What are the units?) for different wavelengths (in nm) and temperatures (in K).\n\nWe use use linear interpolation with flat extrapolation.\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.eq_const-NTuple{10, Any}","page":"API","title":"GasChem.eq_const","text":"Calculates the equilibrium constant Find the backwards reaction by K=kforward/kbackwards Calculates the rate constant of the forward reaction\n\nUsed to compute the rate for these reactions:    PPN        = RCO3 + NO2    PAN        = MCO3 + NO2    ClOO  {+M} = Cl   + O2 {+M}    Cl2O2 {+M} = 2ClO      {+M}\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.j_mean-Union{Tuple{T2}, Tuple{T}, Tuple{Any, Float32, T2, T, T, T}} where {T, T2}","page":"API","title":"GasChem.j_mean","text":"Get mean photolysis rates at different times\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.rate_ALK-NTuple{9, Any}","page":"API","title":"GasChem.rate_ALK","text":"Used to compute the rate for these reactions:     IHOO1    + NO =      NO2 + ...     IHOO4    + NO =      NO2 + ...     IHP001   + NO =      NO2 + ...     IHP002   + NO =      NO2 + ...     IHP003   + NO =      NO2 + ...     IEPOXAOO + NO =      NO2 + ...     IEPOXBOO + NO =      NO2 + ...     ICHOO    + NO =      NO2 + ...     ISOPNOO1 + NO = 1.728NO2 + ...     ISOPNOO2 + NO =      NO2 + ...     IDHNDOO1 + NO =      NO2 + ...     IDHNDOO2 + NO =      NO2 + ...     IDHNBOO  + NO =      NO2 + ...     IDHNDOO  + NO =      NO2 + ...     INO2B    + NO = 2.000NO2 + ...     INO2D    + NO =      NO2 + ...     IHPNBOO  + NO = 1.065NO2 + ...     IHPNDOO  + NO =      NO2 + ...     MVKOHOO  + NO =      NO2 + ...     MCROHOO  + NO =      NO2 + ...\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.rate_DMSOH-NTuple{7, Any}","page":"API","title":"GasChem.rate_DMSOH","text":"Reaction rate for:     DMS + OH = 0.750SO2 + 0.250MSA + MO2\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.rate_EPO-NTuple{6, Any}","page":"API","title":"GasChem.rate_EPO","text":"Used to compute the rate for these reactions:     RIPA   + OH = 0.67IEPOXA   + 0.33IEPOXB   + OH + 0.005LVOC     RIPB   + OH = 0.68IEPOXA   + 0.321IEPOB   + OH + 0.005LVOC     IEPOXA + OH = 0.67IEPOXA00 + 0.33IEPOXB00     IEPOXB + OH = 0.81IEPOXA00 + 0.19IEPOXB00     IHN2   + OH = 0.67IEPOXA   + 0.33IEPOXB   + NO2     IHN3   + OH = 0.67IEPOXA   + 0.33IEPOXB   + NO2     IHN1   + OH = IEPOXD       + NO2     IHN4   + OH = IEPOXD       + NO2     INPB   + OH = OH           + ITHN     INPD   + OH = OH           + ITHN     INPD   + OH = NO2          + ICHE     ICN    + OH = NO2          + ICHE\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.rate_GLYCOH_a-Tuple{Any, Any, Any}","page":"API","title":"GasChem.rate_GLYCOH_a","text":"Used to compute the rate for this reaction:     GLYC + OH = 0.732CH2O + 0.361CO2  + 0.505CO    + 0.227OH               + 0.773HO2  + 0.134GLYX + 0.134HCOOH  which is the \"A\" branch of GLYC + OH.\n\nFor this reaction, these Arrhenius law terms evaluate to 1:     (300/T)^b0 * exp(c0/T)  Because b0 = c0 = 0.\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.rate_GLYCOH_b-Tuple{Any, Any, Any}","page":"API","title":"GasChem.rate_GLYCOH_b","text":"Used to compute the rate for this reaction:     GLYC + OH = HCOOH + OH + CO  which is the \"B\" branch of GLYC + OH.\n\nFor this reaction, these Arrhenius law terms evaluate to 1:     (300/T)^b0 * exp(c0/T)  Because b0 = c0 = 0.\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.rate_GLYXNO3-NTuple{5, Any}","page":"API","title":"GasChem.rate_GLYXNO3","text":"Reaction rate for:     GLYX + NO3 = HNO3 + HO2 + 2CO     i.e. the HO2 + 2*CO branch\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.rate_HACOH_a-NTuple{4, Any}","page":"API","title":"GasChem.rate_HACOH_a","text":"Used to compute the rate for this reaction:     HAC + OH = MGLY + HO2  which is the \"A\" branch of HAC + OH.\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.rate_HACOH_b-NTuple{4, Any}","page":"API","title":"GasChem.rate_HACOH_b","text":"Used to compute the rate for this reaction:     HAC + OH = 0.5HCOOH + OH + 0.5ACTA + 0.5CO2 + 0.5CO + 0.5MO2  which is the \"B\" branch of HAC + OH.\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.rate_HO2HO2-NTuple{8, Any}","page":"API","title":"GasChem.rate_HO2HO2","text":"Used to compute the rate for this reactions:     HO2 + HO2 = H2O2 + O2\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.rate_ISO1-NTuple{9, Any}","page":"API","title":"GasChem.rate_ISO1","text":"Used to compute the rate for these reactions:     ISOP + OH = LISOPOH + IHOO1     ISOP + OH = LISOPOH + IHOO4\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.rate_ISO2-NTuple{9, Any}","page":"API","title":"GasChem.rate_ISO2","text":"Used to compute the rate for these reactions:     ISOP + OH = 0.3MCO3 + 0.3MGLY + 0.3CH2O               + 0.15HPALD3 + 0.25HPALD1 + 0.4HO2               + 0.6CO + 1.5OH + 0.3HPETHNL + LISOPOH     ISOP + OH = 0.3CH2O + 0.15HPALD4 + 0.25HPALD2               + 1.5OH + 0.9CO + 0.7HO2 + 0.3MGLY               + 0.3ATOOH + LISOPOH\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.rate_NIT-NTuple{9, Any}","page":"API","title":"GasChem.rate_NIT","text":"Used to compute the rate for these reactions:     IHOO1    + NO = IHN2     IHOO4    + NO = IHN4     IHPOO1   + NO = IHTN     IHPOO2   + NO = IHTN     IHPOO2   + NO = IHTN     IEPOXAOO + NO = IHTN     IEPOXBOO + NO = IHTN     IHCOO    + NO = IHTN     ISOPNOO1 + NO = IDN     ISOPNOO2 + NO = IDN     IDHNDOO1 + NO = IDN     IDHNDOO2 + NO = IDN     INO2B    + NO = IDN     INO2D    + NO = IDN     IHPNBOO  + NO = IDN     IHPNDOO  + NO = IDN     MVK0HOO  + NO = 0.438MVKN     MCROHOO  + NO = MCRHN\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.rate_OHCO-Tuple{Any, Any, Any}","page":"API","title":"GasChem.rate_OHCO","text":"Reaction rate for:    OH + CO = HO2 + CO2 (cf. JPL 15-10)\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.rate_OHHNO3-NTuple{9, Any}","page":"API","title":"GasChem.rate_OHHNO3","text":"Used to compute the rate for these reactions:    HNO3  + OH = H2O + NO3    HONIT + OH = NO3 + HAC\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.rate_RO2HO2-NTuple{5, Any}","page":"API","title":"GasChem.rate_RO2HO2","text":"Carbon Dependence of RO2+HO2, used in these reactions:    A3O2 + HO2 = RA3P    PO2  + HO2 = PP    KO2  + HO2 = 0.150OH + 0.150ALD2 + 0.150MCO3 + 0.850ATOOH    B3O2 + HO2 = RB3P    PRN1 + HO2 = PRPN\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.rate_RO2NO_a1-NTuple{4, Any}","page":"API","title":"GasChem.rate_RO2NO_a1","text":"Reaction rate for the \"A\" branch of these RO2 + NO reactions:    MO2  + NO = MENO3 in which the \"a1\" parameter equals exactly 1.\n\nFor these reactions, these Arrhenius law terms evaluate to 1:    (300/T)^b0    (300/T)^b1 * exp(c1/T) because b0 = b1 = c1 = 0.\n\nSpecial treatment for methyl nitrate based on observations as Carter and Atkinson formulation does not apply to C1. Value based on upper limit of Flocke et al. 1998 as applied in Fisher et al. 2018\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.rate_RO2NO_a2-NTuple{6, Any}","page":"API","title":"GasChem.rate_RO2NO_a2","text":"\" Reaction rate for the \"A\" branch of these RO2 + NO reactions,     ETO2 + NO = ETNO3     A3O2 + NO = NPRNO3     R4O2 + NO = R4N2     B3O2 + NO = IPRNO3  in which the \"a1\" parameter is greater than 1.0.\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.rate_RO2NO_b1-NTuple{4, Any}","page":"API","title":"GasChem.rate_RO2NO_b1","text":"Reaction rate for the \"B\" branch of these RO2 + NO reactions:     MO2 + NO = CH2O + NO2 + HO2  in which the \"a1\" parameter equals exactly 1.\n\nFor these reactions, these Arrhenius law terms evaluate to 1:     (300/T)^b0     (300/T)^b1 * exp(c1/T)  because b0 = c0 = c1 = 0.\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.rate_RO2NO_b2-NTuple{6, Any}","page":"API","title":"GasChem.rate_RO2NO_b2","text":"Reaction rate for the \"B\" branch of these RO2 + NO reactions:     ETO2 + NO = NO2 +     HO2 + ...     A3O2 + NO = NO2 +     HO2 + ...     R4O2 + NO = NO2 + 0.27HO2 + ...     B3O2 + NO = NO2 +     HO2 + ...  in which the \"a1\" parameter is greater than 1.0.\n\nUse this function when a1 input argument is greater than 1.0.\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.tbranch-NTuple{8, Any}","page":"API","title":"GasChem.tbranch","text":"Temperature Dependent Branching Ratio\n\n\n\n\n\n","category":"method"},{"location":"api/#GasChem.tunplus-NTuple{7, Any}","page":"API","title":"GasChem.tunplus","text":"Used to compute the rate for these reactions:     IHOO1 = 1.5OH + ...     IHOO4 = 1.5OH + ...\n\n\n\n\n\n","category":"method"},{"location":"composing_models/","page":"Composing models","title":"Composing models","text":"CurrentModule = GasChem","category":"page"},{"location":"composing_models/#Composing-models","page":"Composing models","title":"Composing models","text":"","category":"section"},{"location":"composing_models/#Illustrative-Example","page":"Composing models","title":"Illustrative Example","text":"","category":"section"},{"location":"composing_models/","page":"Composing models","title":"Composing models","text":"Here is the complete example of composing, visualizing and solving the SuperFast model and the Fast-JX model, with explanation to follow:","category":"page"},{"location":"composing_models/","page":"Composing models","title":"Composing models","text":"using EarthSciMLBase, GasChem, ModelingToolkit, OrdinaryDiffEq, Dates, Unitful, DifferentialEquations\n\n\n@parameters t\ncomposed_ode = SuperFast(t) + FastJX(t) # Compose two models simply use the \"+\" operator\n\nstart = Dates.datetime2unix(Dates.DateTime(2024, 2, 29))\ntspan = (start, start+3600*24*3)\nsys = structural_simplify(get_mtk(composed_ode)) # Define the coupled system  \nsol = solve(ODEProblem(sys, [], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0) # Solve the coupled system","category":"page"},{"location":"composing_models/","page":"Composing models","title":"Composing models","text":"In the composed system, the variable name for O₃ is not O3 but superfast₊O3(t). So we need some preparation of the result before visualizing. ","category":"page"},{"location":"composing_models/","page":"Composing models","title":"Composing models","text":"vars = states(sys)  # Get the variables in the composed system\nvar_dict = Dict(string(var) => var for var in vars)\npols = [\"O3\", \"OH\", \"NO\", \"NO2\", \"CH4\", \"CH3O2\", \"CO\",\"CH3OOH\", \"CH3O\", \"DMS\", \"SO2\", \"ISOP\"]\nvar_names_p = [\"superfast₊$(v)(t)\" for v in pols]\n\nx_t = unix2datetime.(sol[t]) # Convert from unixtime to date time for visualizing ","category":"page"},{"location":"composing_models/","page":"Composing models","title":"Composing models","text":"Then, we could plot the results as:","category":"page"},{"location":"composing_models/","page":"Composing models","title":"Composing models","text":"using EarthSciMLBase, GasChem, ModelingToolkit, OrdinaryDiffEq, Dates, Unitful, DifferentialEquations\n\n\n@parameters t\ncomposed_ode = SuperFast(t) + FastJX(t) # Compose two models simply use the \"+\" operator\n\nstart = Dates.datetime2unix(Dates.DateTime(2024, 2, 29))\ntspan = (start, start+3600*24*3)\nsys = structural_simplify(get_mtk(composed_ode)) # Define the coupled system  \nsol = solve(ODEProblem(sys, [], tspan, []),AutoTsit5(Rosenbrock23()), saveat=10.0) # Solve the coupled system\n\nvars = states(sys)  # Get the variables in the composed system\nvar_dict = Dict(string(var) => var for var in vars)\npols = [\"O3\", \"OH\", \"NO\", \"NO2\", \"CH4\", \"CH3O2\", \"CO\",\"CH3OOH\", \"CH3O\", \"DMS\", \"SO2\", \"ISOP\"]\nvar_names_p = [\"superfast₊$(v)(t)\" for v in pols]\n\nx_t = unix2datetime.(sol[t]) # Convert from unixtime to date time for visualizing ","category":"page"},{"location":"composing_models/","page":"Composing models","title":"Composing models","text":"using Plots\npp = []\nfor (i, v) in enumerate(var_names_p)\n    name = pols[i]\n    push!(pp, Plots.plot(x_t,sol[var_dict[v]],label = \"$name\", size = (1000, 600), xrotation=45))\nend\nPlots.plot(pp..., layout=(3, 4))","category":"page"},{"location":"geoschem/overview/#GEOS-Chem-\"fullchem\"-Gas-Phase-Mechanism","page":"Overview","title":"GEOS-Chem \"fullchem\" Gas-Phase Mechanism","text":"","category":"section"},{"location":"geoschem/overview/","page":"Overview","title":"Overview","text":"This is an implementation of the GEOS-Chem \"fullchem\" gas-phase mechanism in Julia.  The original version is used in the GEOS-Chem global 3-D atmospheric chemistry transport model. The version here is adapted from GEOS-Chem version 14.1.1. This mechanism is the result of many journal articles which are cited in API documentation for the GEOSChemGasPhase type.","category":"page"},{"location":"geoschem/overview/","page":"Overview","title":"Overview","text":"warning: Warning\nThis implementation is a work in progress. In particular, it does not yet include heterogeneous chemistry.","category":"page"},{"location":"geoschem/overview/#System-overview","page":"Overview","title":"System overview","text":"","category":"section"},{"location":"geoschem/overview/","page":"Overview","title":"Overview","text":"First, let's initialize the model and inspect the first few reactions in the mechanism:","category":"page"},{"location":"geoschem/overview/","page":"Overview","title":"Overview","text":"using GasChem, EarthSciMLBase\nusing DifferentialEquations, ModelingToolkit\nusing Unitful, Plots\n\ntspan = (0.0, 360.0)\n@variables t [unit = u\"s\", description = \"Time\"]\ngc = GEOSChemGasPhase(t)\n\nequations(gc.rxn_sys)[1:5]","category":"page"},{"location":"geoschem/overview/","page":"Overview","title":"Overview","text":"We can also look at the first few equations after converting the reaction network to a system of ODEs:","category":"page"},{"location":"geoschem/overview/","page":"Overview","title":"Overview","text":"equations(gc.sys)[1:5]","category":"page"},{"location":"geoschem/overview/","page":"Overview","title":"Overview","text":"You can see that each reaction has a rate constant; rate constants are specified at the end of the list of equations:","category":"page"},{"location":"geoschem/overview/","page":"Overview","title":"Overview","text":"equations(gc.sys)[end-3:end]","category":"page"},{"location":"geoschem/overview/#Simulation","page":"Overview","title":"Simulation","text":"","category":"section"},{"location":"geoschem/overview/","page":"Overview","title":"Overview","text":"Now, let's run a simulation and plot the results:","category":"page"},{"location":"geoschem/overview/","page":"Overview","title":"Overview","text":"sys = structural_simplify(get_mtk(gc))\nvals = ModelingToolkit.get_defaults(sys)\nfor k in setdiff(states(sys),keys(vals))\n    vals[k] = 0 # Set variables with no default to zero.\nend\nprob = ODEProblem(sys, vals, tspan, vals)\nsol = solve(prob, AutoTsit5(Rosenbrock23()))\nplot(sol, legend = :outertopright, xlabel = \"Time (s)\", \n        ylabel = \"Concentration (nmol/mol)\")","category":"page"},{"location":"geoschem/states/#Chemical-Species","page":"State Variables","title":"Chemical Species","text":"","category":"section"},{"location":"geoschem/states/","page":"State Variables","title":"State Variables","text":"Here is a list of the chemical species in the mechanism:","category":"page"},{"location":"geoschem/states/","page":"State Variables","title":"State Variables","text":"using GasChem, DataFrames, EarthSciMLBase, ModelingToolkit, Unitful\n@variables t [unit = u\"s\", description = \"Time\"]\ngc = structural_simplify(get_mtk(GEOSChemGasPhase(t)))\nvars = states(gc)\nDataFrame(\n        :Name => [string(Symbolics.tosymbol(v, escape=false)) for v ∈ vars],\n        :Units => [ModelingToolkit.get_unit(v) for v ∈ vars],\n        :Description => [ModelingToolkit.getdescription(v) for v ∈ vars],\n        :Default => [ModelingToolkit.getdefault(v) for v ∈ vars])","category":"page"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = GasChem","category":"page"},{"location":"#GasChem:-Gas-Phase-Atmospheric-Chemical-Mechanisms","page":"Home","title":"GasChem: Gas-Phase Atmospheric Chemical Mechanisms","text":"","category":"section"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"To install GasChem.jl, use the Julia package manager:","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"GasChem\")","category":"page"},{"location":"#Features","page":"Home","title":"Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Currently, we have implemented versions of the SuperFast and GEOS-Chem chemical mechanisms, which can be optionally coupled with the Fast-JX photolysis model, which is also implemented here.","category":"page"},{"location":"#Contributing","page":"Home","title":"Contributing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"...coming soon","category":"page"},{"location":"#Reproducibility","page":"Home","title":"Reproducibility","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"<details><summary>The documentation of this EarthSciML package was built using these direct dependencies,</summary>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg # hide\nPkg.status() # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"</details>","category":"page"},{"location":"","page":"Home","title":"Home","text":"<details><summary>and using this machine and Julia version.</summary>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using InteractiveUtils # hide\nversioninfo() # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"</details>","category":"page"},{"location":"","page":"Home","title":"Home","text":"<details><summary>A more complete overview of all dependencies and their versions is also provided.</summary>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Pkg # hide\nPkg.status(; mode = PKGMODE_MANIFEST) # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"</details>","category":"page"},{"location":"","page":"Home","title":"Home","text":"using TOML\nusing Markdown\nversion = TOML.parse(read(\"../../Project.toml\", String))[\"version\"]\nname = TOML.parse(read(\"../../Project.toml\", String))[\"name\"]\nlink_manifest = \"https://github.com/EarthSciML/\" * name * \".jl/tree/gh-pages/v\" * version *\n                \"/assets/Manifest.toml\"\nlink_project = \"https://github.com/EarthSciML/\" * name * \".jl/tree/gh-pages/v\" * version *\n               \"/assets/Project.toml\"\nMarkdown.parse(\"\"\"You can also download the\n[manifest]($link_manifest)\nfile and the\n[project]($link_project)\nfile.\n\"\"\")","category":"page"},{"location":"geoschem/params/#Model-Parameters","page":"Parameters","title":"Model Parameters","text":"","category":"section"},{"location":"geoschem/params/","page":"Parameters","title":"Parameters","text":"The two main parameters are temperature (T (K)) and pressure (num_density (molecules/cm³)). We can explore what happens when we change them:","category":"page"},{"location":"geoschem/params/","page":"Parameters","title":"Parameters","text":"warning: Warning\nThis demonstration does not current work.","category":"page"},{"location":"geoschem/params/","page":"Parameters","title":"Parameters","text":"using GasChem, EarthSciMLBase\nusing DifferentialEquations, ModelingToolkit\nusing Unitful, Plots\n\ntspan = (0.0, 60.0*60*24*4) # 4 day simulation\n@variables t [unit = u\"s\", description = \"Time\"]\nsys = get_mtk(GEOSChemGasPhase(t))\n\n# Convert parameters to variables so we can change them over time.\nsys = param_to_var(sys, :T, :num_density)\n\n# Vary temperature and pressure over time.\n@constants T_0 = 300 [unit=u\"K\"]\n@constants t_0 = 1 [unit=u\"s\"]\neqs = [\n    sys.T ~ T_0 + T_0 / 150 * sin(2π*t/t_0/(60*60*24)),\n    sys.num_density ~ 2.7e19 - 1e19*t/t_0/(60*60*24*4),\n]\nsys = extend(sys,ODESystem(eqs, t; name=:var_T))\n\n# Run the simulation.\nsys = structural_simplify(sys)\nvals = ModelingToolkit.get_defaults(sys)\nfor k in setdiff(states(sys),keys(vals))\n    vals[k] = 0 # Set variables with no default to zero.\nend\nprob = ODEProblem(sys, vals, tspan, vals)\nsol = solve(prob, AutoTsit5(Rosenbrock23()))\np1 = plot(sol, legend = :outertopright, xticks=:none, \n        ylabel = \"Concentration (nmol/mol)\")\n\np2 = plot(sol.t, sol[sys.T], label = \"T\", xticks=:none)\np3 = plot(sol.t, sol[sys.numden], label = \"numden\", xlabel = \"Time (s)\")\n\nplot(p1, p2, p3, layout=grid(3, 1, heights=[0.7, 0.15, 0.15]))","category":"page"},{"location":"geoschem/params/","page":"Parameters","title":"Parameters","text":"Here is a list of all of the model parameters:","category":"page"},{"location":"geoschem/params/","page":"Parameters","title":"Parameters","text":"using GasChem, DataFrames, EarthSciMLBase, ModelingToolkit, Unitful\n@variables t [unit = u\"s\", description = \"Time\"]\ngc = structural_simplify(get_mtk(GEOSChemGasPhase(t)))\nvars = parameters(gc)\nDataFrame(\n        :Name => [string(Symbolics.tosymbol(v, escape=false)) for v ∈ vars],\n        :Units => [ModelingToolkit.get_unit(v) for v ∈ vars],\n        :Description => [ModelingToolkit.getdescription(v) for v ∈ vars],\n        :Default => [ModelingToolkit.getdefault(v) for v ∈ vars])","category":"page"}]
}
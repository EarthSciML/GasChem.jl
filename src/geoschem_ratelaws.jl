# GEOS-Chem rate laws
# Adapted from https://github.com/geoschem/geos-chem/tree/4722f288e90291ba904222f4bbe4fc216d17c34a/KPP/fullchem
# The GEOS-Chem license applies: https://github.com/geoschem/geos-chem/blob/main/LICENSE.txt

""" 
Arrhenius equation:
``` math
    k = a0 * exp( c0 / T ) * (T/300)^b0
```
"""
function arrhenius(T, a0, b0, c0)
    K_300 = 300 # [unit=u"K"]
    k = a0 * exp(c0 / T) * (K_300 / T)^b0
end

"""
Third body effect for pressure dependence of rate coefficients.
a1, b1, c1 are the Arrhenius parameters for the lower-limit rate.
a2, b2, c2 are the Arrhenius parameters for the upper-limit rate.
fv         is the falloff curve paramter, (see ATKINSON ET. AL (1992)
           J. Phys. Chem. Ref. Data 21, P. 1145). Usually fv = 0.6.
"""
function arr_3rdbody(T, num_density, a1, b1, c1, a2, b2, c2, fv)
    rlow = arrhenius(T, a1, b1, c1) * num_density
    rhigh = arrhenius(T, a2, b2, c2)
    xyrat = rlow / rhigh
    blog = log10(xyrat)
    fexp = 1.0 / (1.0 + (blog * blog))
    k = rlow * (fv^fexp) / (1.0 + xyrat)
end


"""
 Used to compute the rate for this reactions:
    HO2 + HO2 = H2O2 + O2
"""
function rate_HO2HO2(T, num_density, H2O, a0, c0, a1, c1)
    k0 = arrhenius(T, a0, 0.0, c0)
    k1 = arrhenius(T, a1, 0.0, c1)
    k = (k0 + k1 * num_density) * (1.0 + 1.4E-21 * H2O * exp(2200.0 / T))
end

"""
Reaction rate for:
   OH + CO = HO2 + CO2 (cf. JPL 15-10)
"""
function rate_OHCO(T, num_density)
    klo1 = arrhenius(T, 5.9E-33, 1, 0)
    khi1 = arrhenius(T, 1.1E-12, -1.3, 0)
    xyrat1 = klo1 * num_density / khi1
    blog1 = log10(xyrat1)
    fexp1 = 1.0 / (1.0 + blog1 * blog1)
    kco1 = klo1 * num_density * 0.6^fexp1 / (1.0 + xyrat1)
    klo2 = arrhenius(T, 1.5E-13, 0, 0)
    khi2 = arrhenius(T, 2.1E+09, -6.1, 0)
    xyrat2 = klo2 * num_density / khi2
    blog2 = log10(xyrat2)
    fexp2 = 1.0 / (1.0 + blog2 * blog2)
    kco2 = klo2 * 0.6^fexp2 / (1.0 + xyrat2)
    k = kco1 + kco2
end

"""
Reaction rate for the "A" branch of these RO2 + NO reactions:
   MO2  + NO = MENO3
in which the "a1" parameter equals exactly 1.

For these reactions, these Arrhenius law terms evaluate to 1:
   (300/T)^b0
   (300/T)^b1 * exp(c1/T)
because b0 = b1 = c1 = 0.

Special treatment for methyl nitrate based on observations
as Carter and Atkinson formulation does not apply to C1.
Value based on upper limit of Flocke et al. 1998 as applied
in Fisher et al. 2018
"""
function rate_RO2NO_a1(T, a0, c0)
    k = a0 * exp(c0 / T) * 3.0e-4
end

"""
Reaction rate for the "B" branch of these RO2 + NO reactions:
    MO2 + NO = CH2O + NO2 + HO2
 in which the "a1" parameter equals exactly 1.
 
 For these reactions, these Arrhenius law terms evaluate to 1:
    (300/T)^b0
    (300/T)^b1 * exp(c1/T)
 because b0 = c0 = c1 = 0.
 """
function rate_RO2NO_b1(T, a0, c0)
    one_minus_fyrno3 = 1.0 - 3.0e-4
    k = a0 * exp(c0 / T) * one_minus_fyrno3
end

"""
Reaction rate for the "A" branch of these RO2 + NO reactions,
    ETO2 + NO = ETNO3
    A3O2 + NO = NPRNO3
    R4O2 + NO = R4N2
    B3O2 + NO = IPRNO3
 in which the "a1" parameter is greater than 1.0.
"""
function rate_RO2NO_a2(T, num_density, a0, c0, a1)
    k0 = arrhenius(T, a0, 0, c0)
    xxyn = 1.94e-22 * exp(0.97 * a1) * num_density
    yyyn = arrhenius(T, 0.826, 8.1, 0.0)
    aaa = log10(xxyn / yyyn)
    zzyn = (1.0 / (1.0 + (aaa * aaa)))
    rarb = (xxyn / (1.0 + (xxyn / yyyn))) * (0.411^zzyn)
    fyrno3 = (rarb / (1.0 + rarb))
    k = k0 * fyrno3
end

"""
Reaction rate for the "B" branch of these RO2 + NO reactions:
    ETO2 + NO = NO2 +     HO2 + ...
    A3O2 + NO = NO2 +     HO2 + ...
    R4O2 + NO = NO2 + 0.27HO2 + ...
    B3O2 + NO = NO2 +     HO2 + ...
 in which the "a1" parameter is greater than 1.0.
 
 Use this function when a1 input argument is greater than 1.0.
 """
function rate_RO2NO_b2(T, num_density, a0, c0, a1)
    k0 = arrhenius(T, a0, 0.0, c0)
    xxyn = 1.94e-22 * exp(0.97 * a1) * num_density
    yyyn = arrhenius(T, 0.826, 8.1, 0.0)
    aaa = log10(xxyn / yyyn)
    zzyn = (1.0 / (1.0 + (aaa * aaa)))
    rarb = (xxyn / (1.0 + (xxyn / yyyn))) * (0.411^zzyn)
    fyrno3 = (rarb / (1.0 + rarb))
    k = k0 * (1.0 - fyrno3)
end

"""
Temperature Dependent Branching Ratio
"""
function tbranch(T, a0, b0, c0, a1, b1, c1)
    k0 = arrhenius(T, a0, b0, c0)
    k1 = arrhenius(T, a1, b1, c1)
    k = k0 / (1.0 + k1)
end

"""
Used to compute the rate for these reactions:
   HNO3  + OH = H2O + NO3
   HONIT + OH = NO3 + HAC
"""
function rate_OHHNO3(T, num_density, a0, c0, a1, c1, a2, c2)
    # ---  OH + HNO3:   K = K0 + K3[M] / (1 + K3[M]/K2)  ------
    k0 = arrhenius(T, a0, 0, c0)
    k1 = arrhenius(T, a1, 0, c1)
    k2 = num_density * arrhenius(T, a2, 0, c2)
    k = k0 + k2 / (1.0 + k2 / k1)
end

"""
Calculates the equilibrium constant
Find the backwards reaction by K=kforward/kbackwards
Calculates the rate constant of the forward reaction

Used to compute the rate for these reactions:
   PPN        = RCO3 + NO2
   PAN        = MCO3 + NO2
   ClOO  {+M} = Cl   + O2 {+M}
   Cl2O2 {+M} = 2ClO      {+M}
"""
function eq_const(T, num_density, a0, c0, a1, b1, a2, b2, fv)
    k0 = arrhenius(T, a0, 0, c0)               # backwards rxn rate
    k1 = arr_3rdbody(T, num_density, a1, b1, 0, a2, b2, 0, fv)  # forwards rxn rate
    k = k1 / k0
end

"""
Carbon Dependence of RO2+HO2, used in these reactions:
   A3O2 + HO2 = RA3P
   PO2  + HO2 = PP
   KO2  + HO2 = 0.150OH + 0.150ALD2 + 0.150MCO3 + 0.850ATOOH
   B3O2 + HO2 = RB3P
   PRN1 + HO2 = PRPN
"""
function rate_RO2HO2(T, a0, c0, a1)
    k = arrhenius(T, a0, 0, c0)
    k * (1.0 - exp(-0.245 * a1))
end

"""
Used to compute the rate for this reaction:
    GLYC + OH = 0.732CH2O + 0.361CO2  + 0.505CO    + 0.227OH
              + 0.773HO2  + 0.134GLYX + 0.134HCOOH
 which is the "A" branch of GLYC + OH.
 
 For this reaction, these Arrhenius law terms evaluate to 1:
    (300/T)^b0 * exp(c0/T)
 Because b0 = c0 = 0.
"""
function rate_GLYCOH_a(T, a0)
    exp_arg = -1.0 / 73.0
    glyc_frac = 1.0 - 11.0729 * exp(exp_arg * T)
    glyc_frac = max(glyc_frac, 0.0)
    k = a0 * glyc_frac
end

"""
Used to compute the rate for this reaction:
    GLYC + OH = HCOOH + OH + CO
 which is the "B" branch of GLYC + OH.
 
 For this reaction, these Arrhenius law terms evaluate to 1:
    (300/T)^b0 * exp(c0/T)
 Because b0 = c0 = 0.
"""
function rate_GLYCOH_b(T, a0)
    exp_arg = -1.0 / 73.0
    glyc_frac = 1.0 - 11.0729 * exp(exp_arg * T)
    glyc_frac = max(glyc_frac, 0.0)
    k = a0 * (1.0 - glyc_frac)
end

"""
Used to compute the rate for this reaction:
    HAC + OH = MGLY + HO2
 which is the "A" branch of HAC + OH.
 """
function rate_HACOH_a(T, a0, c0)
    exp_arg = -1.0 / 60.0
    k0 = arrhenius(T, a0, 0, c0)
    hac_frac = 1.0 - 23.7 * exp(exp_arg * T)
    hac_frac = max(hac_frac, 0.0)
    k = k0 * hac_frac
end

"""
Used to compute the rate for this reaction:
    HAC + OH = 0.5HCOOH + OH + 0.5ACTA + 0.5CO2 + 0.5CO + 0.5MO2
 which is the "B" branch of HAC + OH.
"""
function rate_HACOH_b(T, a0, c0)
    exp_arg = -1.0 / 60.0
    k0 = arrhenius(T, a0, 0, c0)
    hac_frac = 1.0 - 23.7 * exp(exp_arg * T)
    hac_frac = max(hac_frac, 0.0)
    k = k0 * (1.0 - hac_frac)
end

"""
Reaction rate for:
    DMS + OH = 0.750SO2 + 0.250MSA + MO2
"""
function rate_DMSOH(T, num_density, a0, c0, a1, c1)
    k0 = arrhenius(T, a0, 0, c0)
    k1 = arrhenius(T, a1, 0, c1)
    k = (k0 * num_density * 0.2095e0) / (1.0 + k1 * 0.2095e0)
end

"""
Reaction rate for:
    GLYX + NO3 = HNO3 + HO2 + 2CO
    i.e. the HO2 + 2*CO branch
"""
function rate_GLYXNO3(T, num_density, a0, c0)
    # ---  K = K1*([O2]+3.5D18)/(2*[O2]+3.5D18)
    O2 = num_density * 0.2095
    k = arrhenius(T, a0, 0, c0)
    k = k * (O2 + 3.5E+18) / (2.0 * O2 + 3.5E+18)
end

"""
Modified Arrhenius law.
"""
function arrplus(T, a0, b0, c0, d0, e0)
    K_300 = 300 # [unit=u"K"]
    k = a0 * (d0 + (T * e0)) * exp(-b0 / T) * (T / K_300)^c0
    k = max(k, 0.0)
end

"""
Used to compute the rate for these reactions:
    IHOO1 = 1.5OH + ...
    IHOO4 = 1.5OH + ...
"""
function tunplus(T, a0, b0, c0, d0, e0)
    k = a0 * (d0 + (T * e0))
    k = k * exp(b0 / T) * exp(c0 / T^3)
    k = max(k, 0.0)
end

"""
Used to compute the rate for these reactions:
    ISOP + OH = LISOPOH + IHOO1
    ISOP + OH = LISOPOH + IHOO4
"""
function rate_ISO1(T, a0, b0, c0, d0, e0, f0, g0)
    k0 = d0 * exp(e0 / T) * exp(1.0E8 / T^3)
    k1 = f0 * exp(g0 / T)
    k2 = c0 * k0 / (k0 + k1)
    k = a0 * exp(b0 / T) * (1.0 - k2)
end

"""
Used to compute the rate for these reactions:
    ISOP + OH = 0.3MCO3 + 0.3MGLY + 0.3CH2O
              + 0.15HPALD3 + 0.25HPALD1 + 0.4HO2
              + 0.6CO + 1.5OH + 0.3HPETHNL + LISOPOH
    ISOP + OH = 0.3CH2O + 0.15HPALD4 + 0.25HPALD2
              + 1.5OH + 0.9CO + 0.7HO2 + 0.3MGLY
              + 0.3ATOOH + LISOPOH
"""
function rate_ISO2(T, a0, b0, c0, d0, e0, f0, g0)
    k0 = d0 * exp(e0 / T) * exp(1.0E8 / T^3)
    k1 = f0 * exp(g0 / T)
    k2 = c0 * k0 / (k0 + k1)
    k = a0 * exp(b0 / T) * k2
end

"""
Used to compute the rate for these reactions:
    RIPA   + OH = 0.67IEPOXA   + 0.33IEPOXB   + OH + 0.005LVOC
    RIPB   + OH = 0.68IEPOXA   + 0.321IEPOB   + OH + 0.005LVOC
    IEPOXA + OH = 0.67IEPOXA00 + 0.33IEPOXB00
    IEPOXB + OH = 0.81IEPOXA00 + 0.19IEPOXB00
    IHN2   + OH = 0.67IEPOXA   + 0.33IEPOXB   + NO2
    IHN3   + OH = 0.67IEPOXA   + 0.33IEPOXB   + NO2
    IHN1   + OH = IEPOXD       + NO2
    IHN4   + OH = IEPOXD       + NO2
    INPB   + OH = OH           + ITHN
    INPD   + OH = OH           + ITHN
    INPD   + OH = NO2          + ICHE
    ICN    + OH = NO2          + ICHE
"""
function rate_EPO(T, num_density, a1, e1, m1)
    k1 = 1.0 / (m1 * num_density + 1.0)
    k = a1 * exp(e1 / T) * k1
end

function rate_PAN_abab(T, num_density, a0, b0, a1, b1, cf)
    k0 = a0 * exp(b0 / T)
    k1 = a1 * exp(b1 / T)
    k0 = k0 * num_density
    kr = k0 / k1
    nc = 0.75 - 1.27 * (log10(cf))
    f = 10.0^(log10(cf) / (1.0 + (log10(kr) / nc)^2))
    k = k0 * k1 * f / (k0 + k1)
end

function rate_PAN_acac(T, num_density, a0, c0, a1, c1, cf)
    K_300 = 300 # [unit=u"K"]
    k0 = a0 * (T / K_300)^c0
    k1 = a1 * (T / K_300)^c1
    k0 = k0 * num_density
    kr = k0 / k1
    nc = 0.75 - 1.27 * (log10(cf))
    f = 10.0^(log10(cf) / (1.0 + (log10(kr) / nc)^2))
    k = k0 * k1 * f / (k0 + k1)
end

"""
Used to compute the rate for these reactions:
    IHOO1    + NO = IHN2
    IHOO4    + NO = IHN4
    IHPOO1   + NO = IHTN
    IHPOO2   + NO = IHTN
    IHPOO2   + NO = IHTN
    IEPOXAOO + NO = IHTN
    IEPOXBOO + NO = IHTN
    IHCOO    + NO = IHTN
    ISOPNOO1 + NO = IDN
    ISOPNOO2 + NO = IDN
    IDHNDOO1 + NO = IDN
    IDHNDOO2 + NO = IDN
    INO2B    + NO = IDN
    INO2D    + NO = IDN
    IHPNBOO  + NO = IDN
    IHPNDOO  + NO = IDN
    MVK0HOO  + NO = 0.438MVKN
    MCROHOO  + NO = MCRHN
"""
function rate_NIT(T, num_density, a0, b0, c0, n, x0, y0)
    k0 = 2.0E-22 * exp(n)
    k1 = 4.3E-1 * (T / 298.0)^(-8)
    k0 = k0 * num_density
    k1 = k0 / k1
    k2 = (k0 / (1.0 + k1)) * 4.1E-1^(1.0 / (1.0 + (log10(k1))^2))
    k3 = k2 / (k2 + c0)
    k4 = a0 * (x0 - T * y0)
    k = k4 * exp(b0 / T) * k3
    k = max(k, 0.0)
end

"""
Used to compute the rate for these reactions:
    IHOO1    + NO =      NO2 + ...
    IHOO4    + NO =      NO2 + ...
    IHP001   + NO =      NO2 + ...
    IHP002   + NO =      NO2 + ...
    IHP003   + NO =      NO2 + ...
    IEPOXAOO + NO =      NO2 + ...
    IEPOXBOO + NO =      NO2 + ...
    ICHOO    + NO =      NO2 + ...
    ISOPNOO1 + NO = 1.728NO2 + ...
    ISOPNOO2 + NO =      NO2 + ...
    IDHNDOO1 + NO =      NO2 + ...
    IDHNDOO2 + NO =      NO2 + ...
    IDHNBOO  + NO =      NO2 + ...
    IDHNDOO  + NO =      NO2 + ...
    INO2B    + NO = 2.000NO2 + ...
    INO2D    + NO =      NO2 + ...
    IHPNBOO  + NO = 1.065NO2 + ...
    IHPNDOO  + NO =      NO2 + ...
    MVKOHOO  + NO =      NO2 + ...
    MCROHOO  + NO =      NO2 + ...
"""
function rate_ALK(T, num_density, a0, b0, c0, n, x0, y0)
    k0 = 2.0E-22 * exp(n)
    k1 = 4.3E-1 * (T / 298.0)^(-8)
    k0 = k0 * num_density
    k1 = k0 / k1
    k2 = (k0 / (1.0 + k1)) * 4.1E-1^(1.0 / (1.0 + (log10(k1))^2))
    k3 = c0 / (k2 + c0)
    k4 = a0 * (x0 - T * y0)
    k = k4 * exp(b0 / T) * k3
    k = max(k, 0.0)
end
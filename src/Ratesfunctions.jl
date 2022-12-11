"""
The original model is from:
https://github.com/geoschem

"""


include("Lowfunctions.jl")
include("RateLawUtilFuncs.jl")


DU1 = 1  # Dust (Reff = 0.151 um)
DU2 = 2  # Dust (Reff = 0.253 um)
DU3 = 3  # Dust (Reff = 0.402 um)
DU4 = 4  # Dust (Reff = 0.818 um)
DU5 = 5  # Dust (Reff = 1.491 um)
DU6 = 6  # Dust (Reff = 2.417 um)
DU7 = 7  # Dust (Reff = 3.721 um)
SUL = 8  # Tropospheric Sulfate
BKC = 9  # Black Carbon
ORC = 10 # Organic Carbon
SSA = 11 # Accum-mode sea salt
SSC = 12 # Coarse-mode sea salt
SLA = 13 # Strat sulfate liq aer
IIC = 14 # Irregular ice cloud

SS_FINE = 1
SS_COARSE = 2

#Indices for the KHETI_SLA array
H2O = 1
N2O5_plus_H2O = 1
N2O5_plus_HCl = 2   # KHETI_SLA(2) = 0
ClNO3_plus_H2O = 3
ClNO3_plus_HCl = 4
ClNO3_plus_HBr = 5
BrNO3_plus_H2O = 6
BrNO3_plus_HCl = 7
HOCl_plus_HCl = 8
HOCl_plus_HBr = 9   # KHETI_SLA(9) = 0
HOBr_plus_HCl = 10
HOBr_plus_HBr = 11  # KHETI_SLA(11)= 0

# Critical RH [%] for uptake of GLYX, MGLYX, and GLYC:
# REAL(dp), PRIVATE, PARAMETER :: 
CRITRH = 35.0

# Conversion factor from atm to bar
# REAL(dp), PRIVATE, PARAMETER :: 
CON_ATM_BAR = 1.0 / 1.01325

# Reference temperature used in Henry's law
# REAL(dp), PRIVATE, PARAMETER :: 
INV_T298 = 1.0 / 298.15


##########################################################################
######          RATE-LAW FUNCTIONS FOR GAS-PHASE REACTIONS           #####
######   Some common functions are defined in rateLawUtilFuncs.F90   #####
##########################################################################


function ARRPLUS_ade(a0, d0, e0)
    # Modified Arrhenius law, skipping computation of exp( -b0/T )
    # and ( 300/T )^c0 terms, which evaluate to 1 when b0 = c0 = 0.
    # This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # Used to compute the rate for these reactions:
    #    IHOO1 + IHOO1 = 2MVK  + 2HO2 + 2CH2O
    #    IHOO4 + IHOO4 = 2MACR + 2HO2 + 2CH2O
    #    IHOO1 + IHOO4 = MACR + MVK + 2HO2 + 2CH2O
    #    IHOO1 + IHOO1 = HO2 + HC5A + CO + OH +  MVKHP
    #    IHOO4 + IHOO4 = HO2 + HC5A + CO + OH +  MCRHP
    #    IHOO1 + IHOO4 = HO2 + HC5A + CO + OH +  0.5MVKHP + 0.5MCRHP
    #    IHOO1 + MO2   = MVK + 2HO2 + 2CH2O :
    #    IHOO1 + MO2   = CH2O + 0.5HC5A + 1.5HO2 + 0.5MVKHP + 0.5CO + 0.5OH
    #    IHOO4 + MO2   = MACR + 2HO2 + 2CH2O
    #    IHOO4 + MO2   = CH2O + 0.5HC5A + 1.5HO2 +  0.5MCRHP + 0.5CO + 0.5OH
    #
    #   REAL(dp), INTENT(IN) :: a0, d0, e0
    #   REAL(dp)             :: k
    #
    k = a0 * (d0 + (temp * e0))
    k = max(k, 0.0)
    return k
end

function ARRPLUS_abde(a0, b0, d0, e0)
    # Modified Arrhenius law, skipping computation of ( T/300 )^c0,
    # which evaluates to 1 when c0=0.  This avoids excess CPU cycles.
    # (bmy, 12/18/20)
    #
    # Used to compute the rate for these reactions:
    #    IHOO1 + HO2 = 0.063MVK + 0.063OH + 0.063HO2 + 0.063CH2O + 0.937RIPA
    #    IHOO1 + HO2 = RIPC
    #    IHOO4 + HO2 = 0.063MACR + 0.063OH + 0.063HO2 + 0.063CH2O + 0.937RIPB
    #    IHOO4 + HO2 = RIPD
    #    IHOO1       = CH2O + OH + MVK
    #    IHOO4       = MACR + OH + CH2O
    #
    # REAL(dp), INTENT(IN) :: a0, b0, d0, e0
    # REAL(dp)             :: k
    #
    k = a0 * (d0 + (temp * e0)) * exp(-b0 / temp)
    k = max(k, 0.0)
    return k
end

function TUNPLUS_abcde(a0, b0, c0, d0, e0)
    # Used to compute the rate for these reactions:
    #    IHOO1 = 1.5OH + ...
    #    IHOO4 = 1.5OH + ...
    #
    # REAL(dp), INTENT(IN) :: a0, b0, c0, d0, e0
    # REAL(dp)             :: k
    #
    k = a0 * (d0 + (temp * e0))
    k = k * exp(b0 / temp) * exp(c0 / temp^3)
    k = max(k, 0.0)
    return k
end

function GC_ISO1(a0, b0, c0, d0, e0, f0, g0)
    # Used to compute the rate for these reactions:
    #    ISOP + OH = LISOPOH + IHOO1
    #    ISOP + OH = LISOPOH + IHOO4
    #
    # REAL(dp), INTENT(IN) :: a0, b0, c0, d0, e0, f0, g0
    # REAL(dp)             :: k0, k1, k2, k
    #
    k0 = d0 * exp(e0 / temp) * exp(1.0E8 / temp^3)
    k1 = f0 * exp(g0 / temp)
    k2 = c0 * k0 / (k0 + k1)
    k = a0 * exp(b0 / temp) * (1.0 - k2)
    return k
end

function GC_ISO2(a0, b0, c0, d0, e0, f0, g0)
    # Used to compute the rate for these reactions:
    #    ISOP + OH = 0.3MCO3 + 0.3MGLY + 0.3CH2O
    #              + 0.15HPALD3 + 0.25HPALD1 + 0.4HO2
    #              + 0.6CO + 1.5OH + 0.3HPETHNL + LISOPOH
    #    ISOP + OH = 0.3CH2O + 0.15HPALD4 + 0.25HPALD2
    #              + 1.5OH + 0.9CO + 0.7HO2 + 0.3MGLY
    #              + 0.3ATOOH + LISOPOH
    #
    # REAL(dp), INTENT(IN) :: a0, b0, c0, d0, e0, f0, g0
    # REAL(dp)             :: k0, k1, k2, k
    #
    k0 = d0 * exp(e0 / temp) * exp(1.0E8 / temp^3)
    k1 = f0 * exp(g0 / temp)
    k2 = c0 * k0 / (k0 + k1)
    k = a0 * exp(b0 / temp) * k2
    return k
end

function GC_EPO_a(a1, e1, m1)
    # Used to compute the rate for these reactions:
    #    RIPA   + OH = 0.67IEPOXA   + 0.33IEPOXB   + OH + 0.005LVOC
    #    RIPB   + OH = 0.68IEPOXA   + 0.321IEPOB   + OH + 0.005LVOC
    #    IEPOXA + OH = 0.67IEPOXA00 + 0.33IEPOXB00
    #    IEPOXB + OH = 0.81IEPOXA00 + 0.19IEPOXB00
    #    IHN2   + OH = 0.67IEPOXA   + 0.33IEPOXB   + NO2
    #    IHN3   + OH = 0.67IEPOXA   + 0.33IEPOXB   + NO2
    #    IHN1   + OH = IEPOXD       + NO2
    #    IHN4   + OH = IEPOXD       + NO2
    #    INPB   + OH = OH           + ITHN
    #    INPD   + OH = OH           + ITHN
    #    INPD   + OH = NO2          + ICHE
    #    ICN    + OH = NO2          + ICHE
    #
    # REAL(dp), INTENT(IN) :: a1, e1, m1
    # REAL(dp)             :: k1, k
    #
    k1 = 1.0 / (m1 * NUMDEN + 1.0)
    k = a1 * exp(e1 / temp) * k1
    return k
end

function GC_PAN_abab(a0, b0, a1, b1, cf)
    # Used to compute the rate for these reactions:
    #    MACR1OO + NO2 = MPAN
    #    MACRNO2 + NO2 = MPAN + NO2
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    exp(b0/T)
    #    exp(b1/T)
    # because b0 = b1 = 0.  Therefore we can skip computing these
    # terms.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # Sept 27 2012: Added GC_PAN_abab per Kelvin Bates' requirements
    #               for aromatic chem.
    # REAL(dp), INTENT(IN) :: a0, b0, a1, b1, cf
    # REAL(dp)             :: k0, k1, kr, nc, f,  k
    #
    k0 = a0 * exp(b0 / temp)
    k1 = a1 * exp(b1 / temp)
    k0 = k0 * NUMDEN
    kr = k0 / k1
    nc = 0.75 - 1.27 * (log10(cf))
    f = 10.0^(log10(cf) / (1.0 + (log10(kr) / nc)^2))
    k = k0 * k1 * f / (k0 + k1)
    return k
end

function GC_PAN_acac(a0, c0, a1, c1, cf)
    # Used to compute the rate for these reactions:
    #    MACR1OO + NO2 = MPAN
    #    MACRNO2 + NO2 = MPAN + NO2
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    exp(b0/T)
    #    exp(b1/T)
    # because b0 = b1 = 0.  Therefore we can skip computing these
    # terms.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # REAL(dp), INTENT(IN) :: a0, c0, a1, c1, cf
    # REAL(dp)             :: k0, k1, kr, nc, f,  k
    #
    k0 = a0 * temp_over_K300^c0
    k1 = a1 * temp_over_K300^c1
    k0 = k0 * NUMDEN
    kr = k0 / k1
    nc = 0.75 - 1.27 * (log10(cf))
    f = 10.0^(log10(cf) / (1.0 + (log10(kr) / nc)^2))
    k = k0 * k1 * f / (k0 + k1)
    return k
end

function GC_NIT(a0, b0, c0, n, x0, y0)
    # Used to compute the rate for these reactions:
    #    IHOO1    + NO = IHN2
    #    IHOO4    + NO = IHN4
    #    IHPOO1   + NO = IHTN
    #    IHPOO2   + NO = IHTN
    #    IHPOO2   + NO = IHTN
    #    IEPOXAOO + NO = IHTN
    #    IEPOXBOO + NO = IHTN
    #    IHCOO    + NO = IHTN
    #    ISOPNOO1 + NO = IDN
    #    ISOPNOO2 + NO = IDN
    #    IDHNDOO1 + NO = IDN
    #    IDHNDOO2 + NO = IDN
    #    INO2B    + NO = IDN
    #    INO2D    + NO = IDN
    #    IHPNBOO  + NO = IDN
    #    IHPNDOO  + NO = IDN
    #    MVk0HOO  + NO = 0.438MVKN
    #    MCROHOO  + NO = MCRHN
    #
    # REAL(dp), INTENT(IN) :: a0, b0, c0, n,  x0, y0
    # REAL(dp)             :: k0, k1, k2, k3, k4, k
    #
    k0 = 2.0e-22 * exp(n)
    k1 = 4.3e-1 * (temp / 298.0)^(-8)
    k0 = k0 * NUMDEN
    k1 = k0 / k1
    k2 = (k0 / (1.0 + k1)) * 4.1e-1^(1.0 / (1.0 + (log10(k1))^2))
    k3 = k2 / (k2 + c0)
    k4 = a0 * (x0 - temp * y0)
    k = k4 * exp(b0 / temp) * k3
    k = max(k, 0.0)
    return k
end

function GC_ALK(a0, b0, c0, n, x0, y0)
    # Used to compute the rate for these reactions:
    #   IHOO1    + NO =      NO2 + ...
    #   IHOO4    + NO =      NO2 + ...
    #   IHP001   + NO =      NO2 + ...
    #   IHP002   + NO =      NO2 + ...
    #   IHP003   + NO =      NO2 + ...
    #   IEPOXAOO + NO =      NO2 + ...
    #   IEPOXBOO + NO =      NO2 + ...
    #   ICHOO    + NO =      NO2 + ...
    #   ISOPNOO1 + NO = 1.728NO2 + ...
    #   ISOPNOO2 + NO =      NO2 + ...
    #   IDHNDOO1 + NO =      NO2 + ...
    #   IDHNDOO2 + NO =      NO2 + ...
    #   IDHNBOO  + NO =      NO2 + ...
    #   IDHNDOO  + NO =      NO2 + ...
    #   INO2B    + NO = 2.000NO2 + ...
    #   INO2D    + NO =      NO2 + ...
    #   IHPNBOO  + NO = 1.065NO2 + ...
    #   IHPNDOO  + NO =      NO2 + ...
    #   MVKOHOO  + NO =      NO2 + ...
    #   MCROHOO  + NO =      NO2 + ...
    #
    # REAL(dp), INTENT(IN) :: a0, b0, c0, n,  x0, y0
    # REAL(dp)             :: k0, k1, k2, k3, k4, k
    #
    k0 = 2.0e-22 * exp(n)
    k1 = 4.3e-1 * (temp / 298.0)^(-8)
    k0 = k0 * NUMDEN
    k1 = k0 / k1
    k2 = (k0 / (1.0 + k1)) * 4.1e-1^(1.0 / (1.0 + (log10(k1))^2))
    k3 = c0 / (k2 + c0)
    k4 = a0 * (x0 - temp * y0)
    k = k4 * exp(b0 / temp) * k3
    k = max(k, 0.0)
    return k
end

function GC_HO2HO2_acac(a0, c0, a1, c1)
    # Used to compute the rate for these reactions:
    #    HO2 + HO2 = H2O2 + O2
    #
    # For this reaction, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0
    #    (300/T)^b1
    # because b0 = b1 = 0.  Therefore we can skip computing these
    # terms.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # REAL(dp), INTENT(IN) :: a0, c0, a1, c1
    # REAL(dp)             :: k0, k1, k
    #
    k0 = a0 * exp(c0 / temp)
    k1 = a1 * exp(c1 / temp)
    k = (k0 + k1 * NUMDEN) * (1.0 + 1.4E-21 * H2O * exp(2200.0 / temp))
    return k
end

function GC_TBRANCH_1_acac(a0, c0, a1, c1)
    # temperature Dependent Branching Ratio, used for reactions:
    #    MO2 + MO2 = CH2O  + MOH + O2
    #    MO2 + MO2 = 2CH2O + 2HO2
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0
    #    (300/T)^b1
    # because b0 = b1 = 0.  Therefore we can skip computing these
    # terms.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # REAL(dp), INTENT(IN) :: a0, c0, a1, c1
    # REAL(dp)             :: k0, k1, k
    #
    k0 = a0 * exp(c0 / temp)
    k1 = a1 * exp(c1 / temp)
    k = k0 / (1.0 + k1)
    return k
end

function GC_TBRANCH_2_acabc(a0, c0, a1, b1, c1)
    # temperature Dependent Branching Ratio, used for reactions:
    #    C3H8 + OH = B3O2
    #    C3H8 + OH = A3O2
    #
    # For these reactions, this Arrhenius law term evaluates to 1:
    #    (300/T)^b0
    # because b0 = 0.  Therefore we can skip computing this
    # term.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # REAL(dp), INTENT(IN) :: a0, c0, a1, b1, c1
    # REAL(dp)             :: k0, k1, k
    #
    k0 = a0 * exp(c0 / temp)
    k1 = a1 * exp(c1 / temp) * k300_over_temp^b1
    k = k0 / (1.0 + k1)
    return k
end

function GC_RO2HO2_aca(a0, c0, a1)
    # Carbon Dependence of RO2+HO2, used in these reactions:
    #    A3O2 + HO2 = RA3P
    #    PO2  + HO2 = PP
    #    KO2  + HO2 = 0.150OH + 0.150ALD2 + 0.150MCO3 + 0.850ATOOH
    #    B3O2 + HO2 = RB3P
    #    PRN1 + HO2 = PRPN
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0
    #    (300/T)^b1 * exp(c1/T)
    # Because b0 = b1 = c1 = 0.  Therefore we can skip computing these
    # terms.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # REAL(dp), INTENT(IN) :: a0, c0, a1
    # REAL(dp)             :: k
    #
    k = a0 * exp(c0 / temp)
    k = k * (1.0 - exp(-0.245 * a1))
    return k
end

function GC_DMSOH_acac(a0, c0, a1, c1)
    # Reaction rate for:
    #    DMS + OH = 0.750SO2 + 0.250MSA + MO2
    #
    # For this reaction, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0
    #    (300/T)^b1
    # Because b0 = b1 = 0.  Therefore we can skip computing these
    # terms.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # REAL(dp), INTENT(IN) :: a0, c0, a1, c1
    # REAL(dp)             :: k0, k1, k
    #
    k0 = a0 * exp(c0 / temp)
    k1 = a1 * exp(c1 / temp)
    k = (k0 * NUMDEN * 0.2095e0) / (1.0 + k1 * 0.2095e0)
    return k
end

function GC_GLYXNO3_ac(a0, c0)
    # Reaction rate for:
    #    GLYX + NO3 = HNO3 + HO2 + 2CO
    #    i.e. the HO2 + 2*CO branch
    #
    # For this reaction, this Arrhenius term evaluates to 1:
    #    (300/T)^b0
    # because b0 = 0.  Therefore we can skip computing this
    # term.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # REAL(dp), INTENT(IN) :: a0, c0
    # REAL(dp)             :: O2, k
    #
    # ---  K = k1*([O2]+3.5D18)/(2*[O2]+3.5D18)
    O2 = NUMDEN * 0.2095
    k = a0 * exp(c0 / temp)
    k = k * (O2 + 3.5e+18) / (2.0 * O2 + 3.5e+18)
    return k
end

function GC_OHHNO3_acacac(a0, c0, a1, c1, a2, c2)
    # Used to compute the rate for these reactions:
    #    HNO3  + OH = H2O + NO3
    #    HONIT + OH = NO3 + HAC
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0
    #    (300/T)^b1
    #    (300/T)^b2
    # Because b0 = b1 = b2 = 0.  Therefore we can skip computing
    # these terms.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # REAL(dp), INTENT(IN) :: a0, c0, a1, c1, a2, c2
    # REAL(dp)             :: k0, k1, k2, k
    #
    # ---  OH + HNO3:   K = k0 + K3[M] / (1 + K3[M]/K2)  ------
    k0 = a0 * exp(c0 / temp)
    k1 = a1 * exp(c1 / temp)
    k2 = NUMDEN * (a2 * exp(c2 / temp))
    k = k0 + k2 / (1.0 + k2 / k1)
    return k
end

function GC_GLYCOH_A_a(a0)
    # Used to compute the rate for this reaction:
    #    GLYC + OH = 0.732CH2O + 0.361CO2  + 0.505CO    + 0.227OH
    #              + 0.773HO2  + 0.134GLYX + 0.134HCOOH
    # which is the "A" branch of GLYC + OH.
    #
    # For this reaction, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0 * exp(c0/T)
    # Because b0 = c0 = 0.  Therefore we can skip computing these
    # terms.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # REAL(dp), INTENT(IN) :: a0
    # REAL(dp)             :: glyc_frac, k
    # REAL(dp), PARAMETER  :: 
    exp_arg = -1.0 / 73.0
    # #
    glyc_frac = 1.0 - 11.0729 * exp(exp_arg * temp)
    glyc_frac = max(glyc_frac, 0.0)
    k = a0 * glyc_frac
    return k
end

function GC_GLYCOH_B_a(a0)
    # Used to compute the rate for this reaction:
    #    GLYC + OH = HCOOH + OH + CO
    # which is the "B" branch of GLYC + OH.
    #
    # For this reaction, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0 * exp(c0/T)
    # Because b0 = c0 = 0.  Therefore we can skip computing these
    # terms.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # REAL(dp), INTENT(IN) :: a0
    # REAL(dp)             :: glyc_frac, k
    # REAL(dp), PARAMETER  :: 
    exp_arg = -1.0 / 73.0
    #
    glyc_frac = 1.0 - 11.0729 * exp(exp_arg * temp)
    glyc_frac = max(glyc_frac, 0.0)
    k = a0 * (1.0 - glyc_frac)
    return k
end

function GC_HACOH_A_ac(a0, c0)
    # Used to compute the rate for this reaction:
    #    HAC + OH = MGLY + HO2
    # which is the "A" branch of HAC + OH.
    #
    # For this reaction, this Arrhenius law term evaluates to 1:
    #    (300/T)^b0
    # because b0 = 0.  Therefore we can skip computing this
    # term.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # REAL(dp), INTENT(IN) :: a0, c0
    # REAL(dp)             :: k0, hac_frac, k
    # REAL(dp), PARAMETER  :: 
    exp_arg = -1.0 / 60.0
    #
    k0 = a0 * exp(c0 / temp)
    hac_frac = 1.0 - 23.7 * exp(exp_arg * temp)
    hac_frac = max(hac_frac, 0.0)
    k = k0 * hac_frac
    return k
end

function GC_HACOH_B_ac(a0, c0)
    # Used to compute the rate for this reaction:
    #    HAC + OH = 0.5HCOOH + OH + 0.5ACTA + 0.5CO2 + 0.5CO + 0.5MO2
    # which is the "B" branch of HAC + OH.
    #
    # For this reaction, this Arrhenius law term evaluates to 1:
    #    (300/T)^b0}
    # because b0 = 0.  Therefore we can skip computing this
    # term.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # REAL(dp), INTENT(IN) :: a0, c0
    # REAL(dp)             :: k0, hac_frac, k
    # REAL(dp), PARAMETER  :: 
    exp_arg = -1.0 / 60.0
    #
    k0 = a0 * exp(c0 / temp)
    hac_frac = 1.0 - 23.7 * exp(exp_arg * temp)
    hac_frac = max(hac_frac, 0.0)
    k = k0 * (1.0 - hac_frac)
    return k
end

function GC_OHCO_a(a0)
    # Reaction rate for:
    #    OH + CO = HO2 + CO2 (cf. JPL 15-10)
    #
    # For this reaction, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0 * exp(c0/T)
    # because b0 = c0 = 0.  Therefore we can skip computing these
    # terms.  This avoids excess CPU cycles. (bmy, 12/18/20)
    #
    # REAL(dp), INTENT(IN) :: a0
    # #
    # REAL(dp)             :: klo1,   klo2,   khi1,  khi2
    # REAL(dp)             :: xyrat1, xyrat2, blog1, blog2,   fexp1
    # REAL(dp)             :: fexp2,  kco1,   kco2,  temp300, k
    #
    klo1 = 5.9E-33 * k300_over_temp
    khi1 = 1.1E-12 * k300_over_temp^(-1.3)
    xyrat1 = klo1 * NUMDEN / khi1
    blog1 = log10(xyrat1)
    fexp1 = 1.0 / (1.0 + blog1 * blog1)
    kco1 = klo1 * NUMDEN * 0.6^fexp1 / (1.0 + xyrat1)
    klo2 = 1.5E-13
    khi2 = 2.1E+09 * k300_over_temp^(-6.1)
    xyrat2 = klo2 * NUMDEN / khi2
    blog2 = log10(xyrat2)
    fexp2 = 1.0 / (1.0 + blog2 * blog2)
    kco2 = klo2 * 0.6^fexp2 / (1.0 + xyrat2)
    k = kco1 + kco2
    return k
end

function GC_RO2NO_A1_ac(a0, c0)
    # Reaction rate for the "A" branch of these RO2 + NO reactions:
    #    MO2  + NO = MENO3
    # in which the "a1" parameter equals exactly 1.
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0
    #    (300/T)^b1 * exp(c1/T)
    # because b0 = b1 = c1 = 0.  Therefore we can skip computing
    # these terms.  This avoids excess CPU cycles. (bmy, 1/4/20)
    #
    # Special treatment for methyl nitrate based on observations
    # as Carter and Atkinson formulation does not apply to C1.
    # Value based on upper limit of Flocke et al. 1998 as applied
    # in Fisher et al. 2018
    #
    # REAL(dp), INTENT(IN) :: a0, c0
    # REAL(dp)             :: k
    #
    k = a0 * exp(c0 / temp) * 3.0e-4
    return k
end

function GC_RO2NO_B1_ac(a0, c0)
    # Reaction rate for the "B" branch of these RO2 + NO reactions:
    #    MO2 + NO = CH2O + NO2 + HO2
    # in which the "a1" parameter equals exactly 1.
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0
    #    (300/T)^b1 * exp(c1/T)
    # because b0 = c0 = c1 = 0.  Therefore we can skip computing
    # these terms.  This avoids excess CPU cycles. (bmy, 1/4/20)
    #
    # REAL(dp), INTENT(IN) :: a0, c0
    # REAL(dp), PARAMETER  :: 
    one_minus_fyrno3 = 1.0 - 3.0e-4
    # REAL(dp)             :: k
    #
    k = a0 * exp(c0 / temp) * one_minus_fyrno3
    return k
end

function GC_RO2NO_A2_aca(a0, c0, a1)
    # Reaction rate for the "A" branch of these RO2 + NO reactions,
    #    ETO2 + NO = ETNO3
    #    A3O2 + NO = NPRNO3
    #    R4O2 + NO = R4N2
    #    B3O2 + NO = IPRNO3
    # in which the "a1" parameter is greater than 1.0.
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0
    #    (300/T)^b1 * exp(c1/T)
    # because b0 = b1 = c1 = 0.  Therefore we can skip computing
    # these terms.  This avoids excess CPU cycles. (bmy, 1/4/20)
    #
    # REAL(dp), INTENT(IN) :: a0,  c0,   a1
    # REAL(dp)             :: k0,  k, yyyn, xxyn
    # REAL(dp)             :: aaa, rarb, zzyn, fyrno3
    #
    k0 = a0 * exp(c0 / temp)
    xxyn = 1.94e-22 * exp(0.97 * a1) * NUMDEN
    yyyn = 0.826 * ((300.0 / temp)^8.1)
    aaa = log10(xxyn / yyyn)
    zzyn = (1.0 / (1.0 + (aaa * aaa)))
    rarb = (xxyn / (1.0 + (xxyn / yyyn))) * (0.411^zzyn)
    fyrno3 = (rarb / (1.0 + rarb))
    k = k0 * fyrno3
    return k
end

function GC_RO2NO_B2_aca(a0, c0, a1)
    # Reaction rate for the "B" branch of these RO2 + NO reactions:
    #    ETO2 + NO = NO2 +     HO2 + ...
    #    A3O2 + NO = NO2 +     HO2 + ...
    #    R4O2 + NO = NO2 + 0.27HO2 + ...
    #    B3O2 + NO = NO2 +     HO2 + ...
    # in which the "a1" parameter is greater than 1.0.
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0
    #    (300/T)^b1 * exp(c1/T)
    # because b0 = c0 = c1 = 0.  Therefore we can skip computing
    # these terms.  This avoids excess CPU cycles. (bmy, 1/4/20)
    #
    # Use this function when a1 input argument is greater than 1.0.
    # This avoids IF statements, which saves CPU cycles (bmy, 1/4/20)
    #
    # REAL(dp), INTENT(IN) :: a0,  c0,   a1
    # REAL(dp)             :: k0,  k, yyyn, xxyn
    # REAL(dp)             :: aaa, rarb, zzyn, fyrno3
    #
    k0 = a0 * exp(c0 / temp)
    xxyn = 1.94e-22 * exp(0.97 * a1) * NUMDEN
    yyyn = 0.826 * (k300_over_temp^8.1)
    aaa = log10(xxyn / yyyn)
    zzyn = (1.0 / (1.0 + (aaa * aaa)))
    rarb = (xxyn / (1.0 + (xxyn / yyyn))) * (0.411^zzyn)
    fyrno3 = (rarb / (1.0 + rarb))
    k = k0 * (1.0 - fyrno3)
    return k
end

function GCJPLEQ_acabab(a0, c0, a1, b1, a2, b2, fv)
    # Calculates the equilibrium constant
    # Find the backwards reaction by K=kforward/kbackwards
    # Calculates the rate constant of the forward reaction
    #
    # Used to compute the rate for these reactions:
    #    PPN        = RCO3 + NO2
    #    PAN        = MCO3 + NO2
    #    ClOO  {+M} = Cl   + O2 {+M}
    #    Cl2O2 {+M} = 2ClO      {+M}
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b0
    #    exp(c1/T)
    #    exp(c2/T)
    # because b0 = c1 = c2 = 0.  Therefore we can skip computing these terms.
    # Also, fct1 = fct2 = 0, so we will skip those terms as well.  This is
    # more computationally efficient. (bmy, 1/25/20)
    #
    # REAL(dp), INTENT(IN) :: a0, c0, a1, b1, a2, b2, fv
    # REAL(dp)             :: k0, k1, k
    #
    k0 = a0 * exp(c0 / temp)               # backwards rxn rate
    k1 = GCJPLPR_abab(a1, b1, a2, b2, fv)  # forwards rxn rate
    k = k1 / k0
    return k
end

function GCJPLPR_aa(a1, a2, fv)
    # Third body effect for pressure dependence of rate coefficients.
    # a1 is Arrhenius parameters for the lower-limit rate.
    # a2 is Arrhenius parameters for the upper-limit rate.
    # fv is the falloff curve paramter, (see ATKINSON ET. AL (1992)
    # J. Phys. Chem. Ref. Data 21, P. 1145). Usually fv = 0.6.
    #
    # Used to compute the rate for this reaction:
    #    Cl + PRPE {+M} = HCl + PO2 {+M}
    #
    # For this reactions, these Arrhenius law terms evaluate to 1:
    #    (300/T)^b1 * exp(c1/T)
    #    (300/T)^b2 * exp(c2/T)
    # because b1 = b2 = c1 = c2 = 0.  Therefore we can skip computing
    # these terms.  Also, fct1 = fct2 = 0, so we will skip computing
    # these terms as well.  This is more computationally efficient.
    # (bmy, 1/25/20)
    #
    # REAL(dp), INTENT(IN) :: a1,   a2,    fv
    # REAL(dp)             :: rlow, xyrat, blog, fexp, k
    #
    rlow = a1 * NUMDEN
    xyrat = rlow / a2         # rhigh = a2
    blog = log10(xyrat)
    fexp = 1.0 / (1.0 + (blog * blog))
    k = rlow * (fv^fexp) / (1.0 + xyrat)
    return k
end

function GCJPLPR_aba(a1, b1, a2, fv)
    # Third body effect for pressure dependence of rate coefficients.
    # a1, b1 are the Arrhenius parameters for the lower-limit rate.
    # a2     is  the Arrhenius parameters for the upper-limit rate.
    # fv     is the falloff curve paramter, (see ATKINSON ET. AL (1992)
    #        J. Phys. Chem. Ref. Data 21, P. 1145). Usually fv = 0.6.
    #
    # Used to compute the rate for these reactions:
    #    OH  + OH  {+M} = H2O2
    #    NO2 + OH  {+M} = HNO3       {+M}
    #    Cl  + O2  {+M} = ClOO       {+M}
    #    SO2 + OH  {+M} = SO4  + HO2
    #    Br  + NO2 {+M} = BrNO2      {+M}
    #    NO  + O   {+M} = NO2        {+M}
    #    I   + NO2 {+M} = IONO       {+M}
    #    I   + NO  {+M} = INO        {+M}
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    exp(c1/T)
    #    (300/T)^b2 * exp(c2/T)
    # because b2 = c1 = c2 = 0.  Therefore we can skip computing these
    # terms.  Also, fct1 = fct2 = 0, so we will skip computing these
    # terms as well.  This is more computationally efficient.
    # (bmy, 1/25/20)
    #
    # REAL(dp), INTENT(IN) :: a1,   b1,    a2,   fv
    # REAL(dp)             :: rlow, xyrat, blog, fexp, k
    #
    rlow = a1 * (k300_over_temp^1) * NUMDEN
    xyrat = rlow / a2   #rhigh = a2
    blog = log10(xyrat)
    fexp = 1.0 / (1.0 + (blog * blog))
    k = rlow * (fv^fexp) / (1.0 + xyrat)
    return k
end

function GCJPLPR_abab(a1, b1, a2, b2, fv)
    # Third body effect for pressure dependence of rate coefficients.
    # a1, b1 are the Arrhenius parameters for the lower-limit rate.
    # a2, b2 are the Arrhenius parameters for the upper-limit rate.
    # fv     is the falloff curve paramter, (see ATKINSON ET. AL (1992)
    #        J. Phys. Chem. Ref. Data 21, P. 1145). Usually fv = 0.6.
    #
    # Used to compute the rate for these reactions:
    #    NO   + OH  {+M} = HNO2  {+M}
    #    HO2  + NO2 {+M} = HNO4
    #    NO2  + NO3 {+M} = N2O5
    #    ClO  + NO2 {+M} = ClNO3 {+M}
    #    MCO3 + NO2 {+M} = PAN
    #    RCO3 + NO2 {+M} = PPN
    #    PRPE + OH  {+M} = PO2
    #    MO2  + NO2 {+M} = MPN   {+M}
    #    BrO  + NO2 {+M} = BrNO3 {+M}
    #    NO2  + O   {+M} = NO3   {+M}
    #    H    + O2  {+M} = HO2   {+M}
    #    IO   + NO2 {+M} = IONO2 {+M}
    #
    # For these reactions, these Arrhenius law terms evaluate to 1:
    #    exp(c1/T)
    #    exp(c2/T)
    # because c1 = c2 = 0.  Therefore we can skip computing these
    # terms.  Also, fct1 = fct2 = 0, so we will skip computing these
    # terms as well.  This is more computationally efficient.
    # (bmy, 1/25/20)
    #
    # REAL(dp), INTENT(IN) :: a1,   b1,    a2,    b2,   fv
    # REAL(dp)             :: rlow, rhigh, xyrat, blog, fexp, k
    #
    rlow = a1 * (k300_over_temp^b1) * NUMDEN
    rhigh = a2 * (k300_over_temp^b2)
    xyrat = rlow / rhigh
    blog = log10(xyrat)
    fexp = 1.0 / (1.0 + (blog * blog))
    k = rlow * (fv^fexp) / (1.0 + xyrat)
    return k
end

function GCJPLPR_abcabc(a1, b1, c1, a2, b2, c2, fv)
    # Third body effect for pressure dependence of rate coefficients.
    # a1, b1, c1 are the Arrhenius parameters for the lower-limit rate.
    # a2, b2, c2 are the Arrhenius parameters for the upper-limit rate.
    # fv         is the falloff curve paramter, (see ATKINSON ET. AL (1992)
    #           J. Phys. Chem. Ref. Data 21, P. 1145). Usually fv = 0.6.
    #
    # Used to compute the rate for these reactions:
    #    HNO4 {+M} = HO2 + NO2
    #    N2O5 {+M} = NO2 + NO3
    #    MPN  {+M} = MO2 + NO2
    #
    # REAL(dp), INTENT(IN) :: a1,   b1,    c1,    a2,   b2,   c2,  fv
    # REAL(dp)             :: rlow, rhigh, xyrat, blog, fexp, k
    #
    rlow = a1 * (k300_over_temp^b1) * exp(c1 / temp) * NUMDEN
    rhigh = a2 * (k300_over_temp^b2) * exp(c2 / temp)
    xyrat = rlow / rhigh
    blog = log10(xyrat)
    fexp = 1.0 / (1.0 + (blog * blog))
    k = rlow * (fv^fexp) / (1.0 + xyrat)
    return k
end

##########################################################################
######        RATE-LAW FUNCTIONS FOR HETEROGENEOUS REACTIONS         #####
######   Some common functions are defined in rateLawUtilFuncs.F90   #####
##########################################################################

#   =========================================================================
#    Hetchem rate-law functions for BrNO3
#   =========================================================================

function BrNO3uptkByH2O(H)
    #
    # Computes the uptake rate [1/s] for BrNO3 + H2O  (cf. Johan Schmidt)
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    # REAL(dp) :: gamma, gamLiq, gamIce, srMW      # local vars
    #
    k = 0.0
    gamLiq = 0.0021 * temp - 0.561          # Rxn prob, liq (Deiber 2004)
    gamIce = 5.3e-4 * exp(1100.0 / temp) # Rxn prob on ice
    srMw = SR_MW(ind_BrNO3)
    #
    # BrNO3 + H2O on sulfate and sea salt (clear sky)
    gamma = gamLiq
    k = k + Ars_L1K(H % ClearFr * H % xArea(SUL), H % xRadi(SUL), gamma, srMw)
    k = k + Ars_L1K(H % ClearFr * H % xArea(SSA), H % xRadi(SSA), gamma, srMw)
    k = k + Ars_L1K(H % ClearFr * H % xArea(SSC), H % xRadi(SSC), gamma, srMw)
    k = k + H % xArea(SLA) * H % KHETI_SLA(BrNO3_plus_H2O)
    #
    # BrNO3 + H2O uptake on irregular ice cloud (clear sky)
    gamma = 0.3                               # rxn prob, ice [1]
    if H % NatSurface
        gamma = 0.001         # rxn prob, NAT [1]
    end
    k = k + Ars_L1K(H % ClearFr * H % xArea(IIC), H % xRadi(IIC), gamma, srMw)
    #
    # BrNO3 + H2O in tropospheric cloud
    k = k + CloudHet(H, srMw, gamLiq, gamIce, 1.0, 1.0)
    #
    # Assume BrNO3 is limiting, so update the removal rate accordingly
    k = kIIR1Ltd(C(ind_BrNO3), C(ind_H2O), k)
    return k
end

function BrNO3uptkByHCl(H)
    #
    # Computes uptake rate for BrNO3(g) + HCl(l,s)
    # in polar stratospheric clouds and on tropospheric sulfate.
    #
    # TYPE(HetState), INTENT(IN) :: H               # Hetchem State
    # REAL(dp)                   :: k               # rxn prob[1], rxn rate [1/s]
    # REAL(dp)                   :: srMw            # local vars
    #
    k = 0.0
    srMw = SR_MW(ind_BrNO3)
    #
    # Apply BrNO3 uptake in stratosphere
    # NOTE: NAT and ICE both use the same gamma = 0.3

    # % ：： .
    if H % stratBox
        k = k + Ars_L1K(H % xArea(SUL), H % xRadi(SUL), 0.9, srMw)
        k = k + H % xArea(SLA) * H % KHETI_SLA(BrNO3_plus_HCl)
        k = k + Ars_L1K(H % xArea(IIC), H % xRadi(IIC), 0.3, srMw)
    end

    # Assume BrNO3 is limiting, so update the removal rate accordingly
    k = kIIR1Ltd(C(ind_BrNO3), C(ind_HCl), k)
    return k
end

# =========================================================================
# Hetchem rate-law functions for ClNO2
# =========================================================================

function Gam_ClNO2(H, radius, pH, C_Cl, C_Br, gamma, branchCl, branchBr)
    #
    # Calculates reactive uptake coefficient [1] for
    # ClNO2 + Cl- and ClNO2 + Br-.
    #
    # TYPE(HetState), INTENT(IN)   :: H            # Hetchem State
    # REAL(dp),       INTENT(IN)   :: Radius       # Radius [cm]
    # REAL(dp),       INTENT(IN)   :: C_Cl         # Cl- conc [mol/L]
    # REAL(dp),       INTENT(IN)   :: C_Br         # Br- conc [mol/L]
    # REAL(dp),       INTENT(IN)   :: pH           # H+ conc
    # REAL(dp),       INTENT(OUT)  :: gamma        # Rxn prob [1]
    # REAL(dp),       INTENT(OUT)  :: branchCl     # Branching ratio, Cl path
    # REAL(dp),       INTENT(OUT)  :: branchBr     # Branching ratio, Br path
    #
    # REAL(dp), PARAMETER :: INV_AB = 1.0  / 0.01  # 1/mass accom coeff
    # REAL(dp), PARAMETER :: D_l    = 1.0e-5         # Liq phase diffusion coef
    #
    # REAL(dp) :: cavg, gb_tot, H_X, k_Cl
    # REAL(dp) :: k_Br, k_tot,  l_r, M_X
    #
    # Thermal velocity (cm/s)
    M_X = MW(ind_ClNO2) * 1.0e-3
    cavg = sqrt(EIGHT_RSTARG_T / (H % Pi * M_X)) * 100.0
    #
    # Henry's law [M/bar]
    H_X = 4.5e-2 * CON_ATM_BAR
    #
    # Reaction rates (Cl path, Br path, total)
    k_Cl = 1.0e+7 * C_Cl
    if (pH >= 2.0)
        k_Cl = 0.0
    end
    k_Br = (1.01e-1 / (H_X * H_X * D_l)) * C_Br
    k_tot = k_Cl + k_Br
    #
    # Uptake coefficient [1] and branching ratios [1] for Cl, Br paths
    # Prevent div by zero
    gamma = 0.0
    branchCl = 0.0
    branchBr = 0.0
    if (k_tot > 0.0)
        l_r = sqrt(D_l / k_tot)
        gb_tot = FOUR_R_T * H_X * l_r * k_tot / cavg
        gb_tot = gb_tot * ReactoDiff_Corr(radius, l_r)
        gamma = 1.0 / (INV_AB + 1.0 / gb_tot)
        branchCl = k_Cl / k_tot
        branchBr = k_Br / k_tot
    end
end

"""
function ClNO2uptkByBrSALA( H )  
    #
    # Computes the uptake rate [1/s] of ClNO2 + BrSALA.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # Rxn rate [1/s]
    #
    # REAL(dp) :: area,     gamma, branch
    # REAL(dp) :: branchBr, dummy, srMw
    #
    k    = 0.0 
    srMw = SR_MW(ind_ClNO2)
    #
    # ClNO2 + BrSALA uptake rate [1/s] in tropospheric cloud
    if  H%stratBox 
        boolifelse = true
    else
       Gam_ClNO2(H, H%rLiq, H%phCloud, H%Cl_conc_Cld, H%Br_conc_Cld, gamma, dummy, branchBr)
       branch = branchBr * H%frac_Br_CldA
       k      = k + CloudHet( H, srMw, gamma, 0.0 , branch, 0.0  )
    end
    #
    # ClNO2 + BrSALA uptake rate [1/s] on fine sea salt aerosol, in clear sky
    CALL Gam_ClNO2(H,H%aClRadi,H%phSSA(1), H%Cl_conc_SSA,H%Br_Conc_SSA, gamma,dummy,branchBr)
    area = H%ClearFr * H%aClArea
    k    = k + Ars_L1K( area, H%aClRadi, gamma, srMw ) * branchBr

    # Assume ClNO2 is limiting, so recompute reaction rate accordingly
    k = kIIR1Ltd( C(ind_ClNO2), C(ind_BrSALA), k )
    return k 
end

function ClNO2uptkByBrSALC( H )  
    #
    # Computes the uptake rate [1/s] of ClNO2 + BrSALC.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # Rxn rate [1/s]
    #
    # REAL(dp) :: area,     gamma, branch
    # REAL(dp) :: branchBr, dummy, srMw
    #
    k    = 0.0 
    srMw = SR_MW(ind_ClNO2)
    #
    # ClNO2 + BrSALA uptake rate [1/s] in tropospheric cloud
    if  H%stratBox 
        boolifelse = true
    else
       Gam_ClNO2(H,H%rLiq, H%phCloud, H%Cl_conc_Cld,H%Br_conc_Cld,gamma,dummy,branchBr)
       branch = branchBr * H%frac_Br_CldC
       k      = k + CloudHet( H, srMw, gamma, 0.0 , branch, 0.0  )
    end
    #
    # ClNO2 + BrSALA uptake rate [1/s] on fine sea salt aerosol, in clear sky
    Gam_ClNO2(H,H%xRadi(SSC), H%phSSA(2), H%Cl_conc_SSC,H%Br_Conc_SSC, gamma,dummy,branchBr)
    area = H%ClearFr * H%xArea(SSC)
    k    = k + Ars_L1K( area, H%xRadi(SSC), gamma, srMw ) * branchBr

    # Assume ClNO2 is limiting, so recompute reaction rate accordingly
    k = kIIR1Ltd( C(ind_ClNO2), C(ind_BrSALC), k )
    return k 
end

function ClNO2uptkByHBr( H )  
    # Computes the uptake rate [1/s] of ClNO2 + HCl.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # Rxn rate [1/s]
    #
    # REAL(dp) :: gamma, branch, branchBr, dummy, srMw
    #
    k    = 0.0 
    srMw = SR_MW(ind_ClNO2)
    #
    # ClNO2 + HCl uptake rate [1/s] in tropospheric cloud
    if  H%stratBox 
        boolifelse = true
    else
       Gam_ClNO2(H,H%rLiq, H%phCloud, H%Cl_conc_Cld,H%Br_conc_Cld, gamma,dummy,branchBr)
       branch = branchBr * H%frac_Br_CldG
       k      = k + CloudHet( H, srMw, gamma, 0.0 , branch, 0.0  )
    end
    #
    # Assume ClNO2 is limiting, so recompute reaction rate accordingly
    k = kIIR1Ltd( C(ind_ClNO2), C(ind_HBr), k )
    return k 
end

function ClNO2uptkBySALACL( H )  
    #
    # Computes the uptake rate [1/s] of ClNO2 + SALACL.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # Rxn rate [1/s]
    #
    # REAL(dp) :: area,     gamma, branch
    # REAL(dp) :: branchCl, dummy, srMw
    #
    k    = 0.0 
    srMw = SR_MW(ind_ClNO2)
    #
    # ClNO2 + SALACl uptake rate [1/s] in tropospheric cloud
    if  H%stratBox 
        boolifelse = true
    else
       Gam_ClNO2( H,H%rLiq, H%phCloud, H%Cl_conc_Cld,H%Br_conc_Cld, gamma,branchCl,dummy)
       branch = branchCl * H%frac_Cl_CldA
       k      = k + CloudHet( H, srMw, gamma, 0.0 , branch, 0.0  )
    end
    #
    # ClNO2 + SALACL uptake rate [1/s] on fine sea salt aerosol, in clear sky
    Gam_ClNO2(H,H%aClRadi,H%phSSA(1),H%Cl_conc_SSA,H%Br_Conc_SSA, gamma,branchCl,dummy)
    area = H%ClearFr * H%aClArea
    k    = k + Ars_L1K( area, H%aClRadi, gamma, srMw ) * branchCl

    # Assume ClNO2 is limiting, so recompute reaction rate accordingly
    k = kIIR1Ltd( C(ind_ClNO2), C(ind_SALACL), k )
    return k 
end

function ClNO2uptkBySALCCL( H ) 
    # 
    #  Computes the uptake rate [1/s] of ClNO2 + SALACL.
    # 
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # Rxn rate [1/s]
    # 
    # REAL(dp) :: area,     gamma, branch
    # REAL(dp) :: branchCl, dummy, srMw
    # 
    k    = 0.0 
    srMw = SR_MW(ind_ClNO2)
    # #
    # # ClNO2 + SALCCL uptake rate [1/s] in tropospheric cloud
    if H%stratBox 
        ifcontinuelog = 1 #just move
    else
        Gam_ClNO2(H,H%rLiq,H%phCloud,H%Cl_conc_Cld,H%Br_conc_Cld,gamma,branchCl,dummy)
        branch = branchCl * H%frac_Cl_CldC
        k      = k + CloudHet( H, srMw, gamma, 0.0 , branch, 0.0  )
    end
    # #
    # # Assume ClNO2 is limiting, so recompute reaction rate accordingly
    k = kIIR1Ltd( C(ind_ClNO2), C(ind_SALCCL), k )
    return k
end

function ClNO2uptkByHCl( H )  
    # #
    # # Computes the uptake rate [1/s] of ClNO2 + HCl.
    # #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # Rxn rate [1/s]
    # #
    # REAL(dp) :: gamma, branch, branchCl, dummy, srMw
    # #
    k    = 0.0 
    srMw = SR_MW(ind_ClNO2)
    #
    # ClNO2 + HCl uptake rate [1/s] in tropospheric cloud
    if H%stratBox 
        ifcontinuelog = 1 #just move
    else
        Gam_ClNO2(H,H%rLiq, H%phCloud, H%Cl_conc_Cld,H%Br_conc_Cld, gamma, branchCl,  dummy)
        branch = branchCl * H%frac_Cl_CldG
        k      = k + CloudHet( H, srMw, gamma, 0.0 , branch, 0.0  )
    end
    #
    # Assume ClNO2 is limiting, so recompute reaction rate accordingly
    k = kIIR1Ltd( C(ind_ClNO2), C(ind_HCl), k )
    return k
end

function  Gam_ClNO3_Aer( H, C_Br, gamma, branchBr )
    #
    # Calculates reactive uptake coefficients [1] for ClNO3 + Br-.
    #
    #
    # TYPE(HetState), INTENT(IN) :: H
    # REAL(dp), INTENT(IN)       :: C_Br     # Br concentration (mol/L)
    # REAL(dp), INTENT(OUT)      :: gamma    # Rxn prob [1]
    # REAL(dp), INTENT(OUT)      :: branchBr # ClNO3 + HBr- branch ratio [1]
    #
    # REAL(dp), PARAMETER :: INV_AB = 1.0  / 0.108    # 1 / mass accum coeff
    # REAL(dp), PARAMETER :: K_0    = 1.2e+5  ^ 2.0  # H2k0
    # REAL(dp), PARAMETER :: D_l    = 5.0e-6            # Deiber et al 2004
    #
    # REAL(dp) :: M_X, cavg, k_Br, k_tot, gb0, gb2, gb_tot #gbr
    #
    # thermal velocity (cm/s)
    M_X      = MW(ind_ClNO3) * 1.0e-3 
    cavg     = sqrt( EIGHT_RSTARG_T / ( H%PI * M_X ) ) * 100.0 
    #
    # H2k2br cm2 s-1.
    k_Br     = 1.0e+12  * C_Br
    #
    # Calculate gb1 for ClNO3 + Cl-
    # Following [Deiber et al., 2004], gamma is not significantly different
    # from ClNO3 + H2O (gamma = 0.0244) independent of Cl- concentration,
    # but Cl2 rather than HOCl formed. gb2 can be calculated reversely from
    # gb1 = gb0 hydrolysis
    gb0      = FOUR_R_T * 1.2e+5  *  sqrt( D_l ) / cavg
    k_tot    = K_0 + k_Br                                    #H2(k0+k2Br)
    gb_tot   = FOUR_R_T *  sqrt( k_tot * D_l ) / cavg
    #
    # Reaction probability for ClNO3 + Br- [1]
    gamma    = 1.0  / ( INV_AB + 1.0  / gb_tot )
    #
    # Branching ratio for ClNO3 + HBr-
    # BOTE: ClNO3 + Cl- branch ratio = 1.0 - branchBr
    branchBr = k_Br / k_tot
end 

function Gam_ClNO3_Ice( H, gamma, brHCl, brHBr, brH2O )
    #
    # Computes the reactive uptake probability and branching ratio
    # for ClNO3 + H2O, ClNO3 + HCl, and ClNO3 + HBr in ice clouds
    #
    # TYPE(HetState), INTENT(IN)  :: H             # Hetchem State
    # REAL(dp),       INTENT(OUT) :: gamma         # Uptake prob [1]
    # REAL(dp),       INTENT(OUT) :: brHCl         # ClNO3 + HCl branch ratio
    # REAL(dp),       INTENT(OUT) :: brHBr         # ClNO3 + HBr branch ratio
    # REAL(dp),       INTENT(OUT) :: brH2O         # ClNO3 + H2O branch ratio
    #
    # REAL(dp), PARAMETER :: twenty = 1.0  / 0.5 
    #
    # REAL(dp) :: cavg, g1, g2, g3, H2Os, kks, M_X
    #
    # ClNO3 + HCl uptake probability [1]
    g1    = 0.24  * H%HCl_theta
    #
    # ClNO3 + HBr uptake probability [1]
    g2    = 0.56  * H%HBr_theta
    #
    # ClNO3 + H2O uptake probability [1]
    M_X   = MW(ind_ClNO3) * 1.0e-3 
    cavg  = sqrt( EIGHT_RSTARG_T / ( H%PI * M_X ) ) * 100.0 
    H2Os  = 1e+15  - ( 3.0  * 2.7e+14  * H%HNO3_theta )
    kks   = 4.0  * 5.2e-17  * exp( 2032.0  / temp )
    g3    = 1.0  / ( twenty + cavg / ( kks * H2Os ) )   # 1.0/0.5 = 20
    #
    # Total reaction probability
    gamma = g1 + g2 + g3
    #
    # Branching ratios for each path (HCl, HBr, H2O)
    brHCl = g1 / gamma
    brHBr = g2 / gamma
    brH2O = g3 / gamma
end


function ClNO3uptkByH2O( H )  
    # Computes the hydrolysis reaction rate [1/s] of ClNO3 + H2O.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # Rxn rate [1/s]
    # #
    # REAL(dp) :: area,       branchBr, branchLiq
    # REAL(dp) :: branchIce,  dum1,     dum2
    # REAL(dp) :: gamma,      gammaIce, srMw
    #
    k    = 0.0 
    srMw = SR_MW(ind_ClNO3)
    #
    # Rxn rate of ClNO3 + H2O on fine sea salt in clear sky
    Gam_ClNO3_Aer( H, H%Br_conc_SSA, gamma, branchBr )
    area      = H%ClearFr * H%aClArea
    branchLiq = ( 1.0  - branchBr ) * ( 1.0 - H%frac_SALACL )
    k         = k + Ars_L1K( area, H%aClRadi, gamma, srMw ) * branchLiq
    #
    # Rate of ClNO3 + H2O on stratospheric liquid aerosol
    k = k + H%xArea(SLA) * H%KHETI_SLA(ClNO3_plus_H2O)
    #
    # Rate of ClNO3 + H2O on irregular ice cloud
    gamma = 0.3                                # Rxn prob, ice [1]
    if H%NatSurface 
        gamma = 0.004          # Rxn prob, NAT [1]
    end
    k = k + Ars_L1K( H%xArea(IIC), H%xRadi(IIC), gamma, srMw )
    #
    if H%stratBox 
        ifcontinuelog = 1 #just move
    else
        #
        # ClNO3 + H2O uptake prob [1] in liquid tropospheric cloud
        Gam_ClNO3_Aer( H, H%Br_conc_Cld, gamma, branchBr )
        branchLiq = 1.0  - branchBr
        #
        # ClNO3 + H2O uptake prob [1] in tropospheric ice cloud
        Gam_ClNO3_Ice( H, gammaIce, dum1, dum2, branchIce )
        #
        # ClNO3 + H2O rxn rate in cloudy tropopsheric grid box
        k = k + CloudHet( H, srMw, gamma, gammaIce, branchLiq, branchIce )
    end
    #
    # Assume ClNO3 is limiting, so recompute reaction rate accordingly
    k = kIIR1Ltd( C(ind_ClNO3), C(ind_H2O), k )
    return k
end

function ClNO3uptkByHCl( H )  
    #
    # Computes the rate [1/s] of ClNO3(g) + HCl(l,s).
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # Rxn rate [1/s]
    # #
    # REAL(dp) :: branchIce, dum1,     dum2
    # REAL(dp) :: gamma,     gammaIce, srMw
    #
    k    = 0.0 
    srMw = SR_MW(ind_ClNO3)
    #
    if H%stratBox
    # IF ( H%stratBox ) THEN
        #
        # Rxn rate of ClNO3 + HCl on tropospheric sulfate in stratosphere
        gamma = 0.1e-4 
        k = k + Ars_L1K( H%xArea(SUL), H%xRadi(SUL), gamma, srMw )
        #
        # Rate of ClNO3 + HCl on stratospheric liquid aerosol
        k = k + H%xArea(SLA) * H%KHETI_SLA(ClNO3_plus_HCl)
        #
        # Rate of ClNO3 + HCl on irregular ice cloud
        gamma = 0.3                                # Rxn prob, ice [1]
        if H%NatSurface
            gamma = 0.2            # Rxn prob, NAT [1]
        end
        k = k + Ars_L1K( H%xArea(IIC), H%xRadi(IIC), gamma, srMw )
    else
       #
       # NOTE: No ClNO3 + HCl uptake in tropospheric liquid cloud
       #
       # ClNO3 + HCl uptake rate in tropospheric ice cloud
       CALL Gam_ClNO3_Ice( H, gammaIce, branchIce, dum1, dum2 )
       k = k + CloudHet( H, srMw, 0.0 , gammaIce, 0.0 , branchIce )
    end
    #
    # Assume ClNO3 is limiting, so recompute reaction rate accordingly
    k = kIIR1Ltd( C(ind_ClNO3), C(ind_HCl), k )
    return k
end

function ClNO3uptkByHBr( H ) 
    # Computes the reaction rate [1/s] of ClNO3(g) + HBr-.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # Rxn rate [1/s]
    # #
    # REAL(dp) :: branchBr, branchLiq, branchIce, dum1
    # REAL(dp) :: dum2,     gamma,     gammaIce,  srMw
    #
    k    = 0.0 
    srMw = SR_MW(ind_ClNO3)
    #
    if H%stratBox
       # ClNO3 + HBr uptake rate on stratospheric liquid aerosol
       k = k + H%xArea(SLA) * H%KHETI_SLA(ClNO3_plus_HBr)
       #
       # ClNO3 + HBr uptake rate  on irregular ice cloud
       gamma = 0.3              # Rxn prob, ice and NAT [1]
       k = k + Ars_L1K( H%xArea(IIC), H%xRadi(IIC), gamma, srMw )
    else
       #
       # ClNO3 + HBr uptake rate in tropospheric liquid cloud
       Gam_ClNO3_Aer( H, H%Br_conc_Cld, gamma, branchBr )
       branchLiq = branchBr * H%frac_Br_CldG
       #
       # ClNO3 + HBr uptake rate in tropospheric ice cloud
       Gam_ClNO3_Ice( H, gammaIce, dum1, branchIce, dum2 )
       #
       # ClNO3 + HBr overall uptake rate, accounting for cloud fraction
       k = CloudHet( H, srMw, gamma, gammaIce, branchLiq, branchIce )
    end
    #
    # Assume ClNO3 is limiting, so recompute reaction rate accordingly
    k = kIIR1Ltd( C(ind_ClNO3), C(ind_HBr), k )
    return k
end

function FUNCTION ClNO3uptkByBrSALA( H ) 
    #
    # Computes rxn rate [1/s] of ClNO3 + BrSALA.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # Rxn rate [1/s]
    #
    # REAL(dp) :: area, branch, branchBr, gamma, srMw
    #
    k    = 0.0 
    srMw = SR_MW(ind_ClNO3)
    #
    # First compute uptake of ClNO3 + BrSALA in tropospheric cloud
    if H%stratBox
        ifcontinuelog = 1 #continue
    else
        #
        # Compute ClNO3 + BrSALA uptake rate & branching ratio
        Gam_ClNO3_Aer( H, H%Br_conc_Cld, gamma, branchBr )
        branch = branchBr * H%frac_Br_CldA
        #
        # Compute ClNO3 + BrSALA uptake rate accounting for cloud fraction
        k = k + CloudHet( H, srMw, gamma, 0.0 , branch, 0.0  )
    end
    #
    # Compute uptake rate of ClNO3 + BrSALA in clear sky
    Gam_ClNO3_Aer( H, H%Br_conc_SSA, gamma, branchBr )
    area = H%ClearFr * H%aClArea
    k    = k + Ars_L1K( area, H%aClRadi, gamma, srMw ) * branchBr
    #
    # Assume ClNO3 is limiting, so recompute reaction rate accordingly
    k = kIIR1Ltd( C(ind_ClNO3), C(ind_BrSALA), k )
    return k
end

function ClNO3uptkByBrSALC( H )  
    #
    # Computes rxn rate [1/s] of ClNO3 + BrSALC.
    #
    # TYPE(HetState), INTENT(IN) :: H
    # REAL(dp)                   :: k
    # #
    # REAL(dp) :: area, branch, branchBr, gamma, srMw
    #
    k    = 0.0 
    srMw = SR_MW(ind_ClNO3)
    #
    # First compute uptake of ClNO3 + BrSALA in tropospheric cloud
    if H%stratBox
        ifcontinuelog = 1 #continue
    else
       #
       # Compute ClNO3 + BrSALA uptake rate & branching ratio
       Gam_ClNO3_Aer( H, H%Br_conc_Cld, gamma, branchBr )
       branch = branchBr * H%frac_Br_CldC
       #
       # Compute ClNO3 + BrSALA uptake rate accounting for cloud fraction
       k = k + CloudHet( H, srMw, gamma, 0.0 , branch, 0.0  )
    end
    #
    # Compute uptake rate of ClNO3 + BrSALA in clear sky
    Gam_ClNO3_Aer( H, H%Br_conc_SSC, gamma, branchBr )
    area = H%ClearFr * H%xArea(SSC)
    k    = k + Ars_L1K( area, H%xRadi(SSC), gamma, srMw ) * branchBr
    #
    # Assume ClNO3 is limiting, so recompute reaction rate accordingly
    k = kIIR1Ltd( C(ind_ClNO3), C(ind_BrSALC), k )
    return k
end

function ClNO3uptkBySALACL( H )  
     #
     # Computes rxn rate [1/s] of ClNO3 + SALACL.
     #
    # TYPE(HetState), INTENT(IN) :: H
    # REAL(dp)                   :: k
    #  #
    # REAL(dp) :: area, branch, branchBr, gamma, srMw
     #
    k    = 0.0 
    srMw = SR_MW(ind_ClNO3)
     #
     # Compute uptake rate of ClNO3 + BrSALA in clear sky
    Gam_ClNO3_Aer( H, H%Br_conc_SSA, gamma, branchBr )
    area   = H%ClearFr * H%aClArea
    branch = ( 1.0  - branchBr ) * H%frac_SALACL
    k      = k + Ars_L1K( area, H%aClRadi, gamma, srMw )* branch
     #
     # Assume ClNO3 is limiting, so recompute reaction rate accordingly
    k = kIIR1Ltd( C(ind_ClNO3), C(ind_SALACL), k )
    return k
end

function ClNO3uptkBySALCCL( H ) RESULT( k )
    #
    # Computes rxn rate [1/s] of ClNO3 + SALCCL.
    #
    TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    REAL(dp)                   :: k              # Rxn rate [1/s]
    #
    REAL(dp) :: area, branch, branchBr, gamma, srMw
    #
    k    = 0.0 
    srMw = SR_MW(ind_ClNO3)
    #
    # Compute uptake rate of ClNO3 + BrSALA in clear sky
    CALL Gam_ClNO3_Aer( H, H%Br_conc_SSC, gamma, branchBr )
    area   = H%ClearFr * H%xArea(SSC)
    branch = 1.0 - branchBr
    k      = k + Ars_L1K( area, H%xRadi(SSC), gamma, srMw )* branch
    #
    # Assume ClNO3 is limiting, so recompute reaction rate accordingly
    k = kIIR1Ltd( C(ind_ClNO3), C(ind_SALCCL), k )
    return k
end

# #=========================================================================
# # Hetchem rate-law functions for HBr
# #=========================================================================

function HBrUptkBySALA( H )  
    #
    # Computes uptake rate of HBr on fine sea salt (in clear-sky).
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    # REAL(dp)                   :: area, gamma    # local vars
    #
    area  = H%ClearFr * H%aClArea
    gamma = 1.3e-8  * exp( 4290.0  / temp )
    k     = Ars_L1K( area, H%aClRadi, gamma, SR_MW(ind_HBr) )
    return k
end

function HBrUptkBySALC( H ) 
    #
    # Computes uptake rate of HBr on coarse sea salt (in clear-sky).
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    # REAL(dp)                   :: area, gamma    # local vars
    #
    area  = H%ClearFr * H%xArea(SSC)
    gamma = 1.3e-8  * exp( 4290.0  / temp )
    k     = Ars_L1K( area, H%xRadi(SSC), gamma, SR_MW(ind_HBr) )
    return k
end


# =========================================================================
#  Hetchem rate-law functions for HO2
# =========================================================================

function HO2uptk1stOrd( H )  
  #
  # Computes the reaction rate [1/s] for 1st order uptake of HO2.
  #
#   TYPE(HetState), INTENT(IN) :: H              # HetChem State
#   REAL(dp)                   :: srMw, k        # sqrt(mol wt), rxn rate [1/s]
  #
    k    = 0.0 
    srMw = SR_MW(ind_HO2)
    #
    # Uptake by various aerosol types
    k = k + Ars_L1k( H%xArea(DU1), H%xRadi(DU1), H%gamma_HO2, srMw )
    k = k + Ars_L1k( H%xArea(DU2), H%xRadi(DU2), H%gamma_HO2, srMw )
    k = k + Ars_L1k( H%xArea(DU3), H%xRadi(DU3), H%gamma_HO2, srMw )
    k = k + Ars_L1k( H%xArea(DU4), H%xRadi(DU4), H%gamma_HO2, srMw )
    k = k + Ars_L1k( H%xArea(DU5), H%xRadi(DU5), H%gamma_HO2, srMw )
    k = k + Ars_L1k( H%xArea(DU6), H%xRadi(DU6), H%gamma_HO2, srMw )
    k = k + Ars_L1k( H%xArea(DU7), H%xRadi(DU7), H%gamma_HO2, srMw )
    k = k + Ars_L1k( H%xArea(SUL), H%xRadi(SUL), H%gamma_HO2, srMw )
    k = k + Ars_L1k( H%xArea(BKC), H%xRadi(BKC), H%gamma_HO2, srMw )
    k = k + Ars_L1k( H%xArea(ORC), H%xRadi(ORC), H%gamma_HO2, srMw )
    k = k + Ars_L1k( H%xArea(SSA), H%xRadi(SSA), H%gamma_HO2, srMw )
    k = k + Ars_L1k( H%xArea(SSC), H%xRadi(SSC), H%gamma_HO2, srMw )
    return k
end

# =========================================================================
#  Hetchem rate-law functions for HOBr
# =========================================================================
function Br2_Yield( Br_over_Cl )  
    #
    # Returns yield [1] of Br2 from the Br- / Cl- ratio.
    #
    # REAL(dp), INTENT(IN) :: Br_over_Cl           # Br- / Cl- ratio
    # REAL(dp)             :: Y_Br2                # local vars

    # Yield of Br2
    Y_Br2 = 0.0 
    if  Br_over_Cl > 0.0
       Y_Br2 = 0.41  * log10( Br_over_Cl ) + 2.25 
       Y_Br2 = max( min( Y_Br2, 0.9  ), 0.0  )
    end
    return Y_Br2
end

function Gam_HOBr_Aer( H, radius, C_Hp, C_Clm, C_Brm, gamma )
    #
    # Returns uptake probability [1] for HOBr on aerosols.
    #
    # TYPE(HetState), INTENT(IN)  :: H             # Hetchem State
    # REAL(dp),       INTENT(IN)  :: radius        # Aerosol radius
    # REAL(dp),       INTENT(IN)  :: C_Hp          # H+ concentration
    # REAL(dp),       INTENT(IN)  :: C_Clm         # Cl- concentration
    # REAL(dp),       INTENT(IN)  :: C_Brm         # Br- concentration
    # REAL(dp),       INTENT(OUT) :: gamma         # Uptake probability [1/s]
    #
    # REAL(dp) :: M_X,   cavg,      H_X
    # REAL(dp) :: l_r,   C_Hp1,     C_Hp2
    # REAL(dp) :: k_tot, k_HOBr_Cl, k_HOBr_Br, gb_tot

    #
    # REAL(dp), PARAMETER :: INV_AB = 1.0  / 0.6  # Inv. mass accum coef
    # REAL(dp), PARAMETER :: D_l    = 1.4e-5        # Amman et al, ACP, 2013
    #
    # Henry's law
    H_X = ( HENRY_k0(ind_HOBr) * CON_ATM_BAR )* exp( HENRY_CR(ind_HOBr) * ( 1.0 /temp - INV_T298 ) )
    #
    # Thermal velocity [cm/s]
    M_X       = MW(ind_HOBr) * 1.0e-3 
    cavg      = sqrt( EIGHT_RSTARG_T / ( H%PI * M_X ) ) * 100.0 
    #
    # Follow Roberts et al, (2014)
    C_Hp1     = max( min( C_Hp, 1.0e-6  ), 1.0e-9  )
    C_Hp2     = max( min( C_Hp, 1.0e-2  ), 1.0e-6  )
    #
    # Rates for each HOBr + {Cl-, Br-} rxn
    k_HOBr_Cl = 2.3e+10  * C_Clm     * C_Hp1  # Liu & Margerum, EST, 2001
    k_HOBr_Br = 1.6e+10  * C_Brm     * C_Hp2  # ??
    k_tot     = k_HOBr_Cl  + k_HOBr_Br
    #
    # Compute reactive uptake coefficient [unitless], prevent div by zero
    # l_r is diffusive length scale [cm];
    # gb is Bulk reaction coefficient [unitless]
    gamma     = 0.0 
    if  k_tot > 0.0 
    l_r    = sqrt( D_l / k_tot )
    gb_tot = FOUR_R_T * H_X * l_r * k_tot / cavg
    gb_tot = gb_tot * ReactoDiff_Corr( radius, l_r )
    gamma  = 1.0  / ( INV_AB + 1.0 / gb_tot )
    end
end


function Gam_HOBr_Cld( H,gamma,k_tot,k_HOBr_Cl, k_HOBr_Br, k_HOBr_HSO3, k_HOBr_HSO3_2 )
    #
    # Returns uptake probability for HOBr in clouds.
    #
    # TYPE(HetState), INTENT(IN)  :: H        # Hetchem species metadata
    # REAL(dp),       INTENT(OUT) :: gamma
    # REAL(dp),       INTENT(OUT) :: k_tot
    # REAL(dp),       INTENT(OUT) :: k_HOBr_Cl
    # REAL(dp),       INTENT(OUT) :: k_HOBr_Br
    # REAL(dp),       INTENT(OUT) :: k_HOBr_HSO3
    # REAL(dp),       INTENT(OUT) :: k_HOBr_HSO3_2
    #
    # REAL(dp) :: gd,   M_X,  cavg,  H_X,   gb_tot
    # REAL(dp) :: ybr2, l_r,  C_Hp1, C_Hp2, Br_over_Cl
    #
    INV_AB = 1.0  / 0.6  # Inv. mass accum coef
    D_l    = 1.4e-5        # Amman et al, ACP, 2013
    #
    # Henry's law
    H_X    = ( HENRY_k0(ind_HOBr) * CON_ATM_BAR )* exp( HENRY_CR(ind_HOBr) * ( 1.0 /temp - INV_T298 ) )
    #
    # Thermal velocity [cm/s]
    M_X    = MW(ind_HOBr) * 1.0e-3 
    cavg   = sqrt( EIGHT_RSTARG_T / ( H%PI * M_X ) ) * 100.0 
    #
    # Follow Roberts et al, (2014)
    C_Hp1  = min( H%H_conc_lCl, 1.0e-6  )
    C_Hp2  = min( H%H_conc_lCl, 1.0e-2  )
    C_Hp1  = max( C_Hp1,        1.0e-9  )
    C_Hp2  = max( C_Hp2,        1.0e-6  )
    #
    # Rates for each HOBr + {Cl-, Br-, HSO3-, HSO3--} rxn
    k_HOBr_Cl     = 2.3e+10  * H%Cl_conc_Cld * C_Hp1  # Liu & Margerum, EST, 2001
    k_HOBr_Br     = 1.6e+10  * H%Br_conc_Cld * C_Hp2  # ??
    k_HOBr_HSO3   = 2.6e+7   * H%HSO3_aq # Liu and Abbatt, GRL, 2020
    k_HOBr_HSO3_2 = 5.0e+9   * H%SO3_aq  # Troy & Margerum, Inorg. Chem., 1991
    #
    # Total rate
    k_tot  = k_HOBr_Cl + k_HOBr_Br + k_HOBr_HSO3 + k_HOBr_HSO3_2
    #
    # Compue reactive uptake coefficient [unitless], prevent div by zero
    # l_r is diffusive length scale [cm];
    # gb is Bulk reaction coefficient [unitless]
    gamma  = 0.0 
    IF ( k_tot > 0.0  ) THEN
    l_r    = sqrt( D_l / k_tot )
    gb_tot = FOUR_R_T * H_X * l_r * k_tot / cavg
    gb_tot = gb_tot * ReactoDiff_Corr( H%rLiq, l_r )
    gamma  = 1.0  / ( INV_AB + 1.0 / gb_tot )
    end
end

function Gam_HOBr_Ice( H, gamma, branch_HCl, branch_HBr )
    #
    # Calculates total reactive uptake coefficient for
    # HOBr + HCl and HOBr + HBr in ice clouds.
    #
    #   TYPE(HetState), INTENT(IN)  :: H              # Hetchem State
    #   REAL(dp),       INTENT(OUT) :: gamma          # Total rxn prob [1]
    #   REAL(dp),       INTENT(OUT) :: branch_HCl     # HCl branch ratio
    #   REAL(dp),       INTENT(OUT) :: branch_HBr     # HBr branch ratio
    #
    #   REAL(dp) :: gamma_HCl, gamma_HBr              # local vars
    #
    # Overall uptake prob. of  HOBr+HCl and HOBr+HBr together
    gamma_HCl = H%HCl_theta * 0.25 
    gamma_HBr = H%HBr_theta * 4.8e-4  * exp( 1240.0  / temp )
    gamma     = gamma_HCl + gamma_HBr
    #
    # Branching ratios HCl/total and HBr/total
    branch_HCl = 0.0 
    branch_HBr = 0.0 
    if gamma > 0.0 
        branch_HCl = gamma_HCl / gamma
        branch_HBr = gamma_HBr / gamma
    end
end


function HOBrUptkByHBr( H )  
    #
    # Computes the uptake rate [1/s] for the HOBr + HBr reaction
    # in the stratosphere and in tropospheric clouds.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    # REAL(dp) :: branch,       branch_0,     brIce
    # REAL(dp) :: brLiq,        dummy,        gammaLiq
    # REAL(dp) :: gammaIce,     k_HOBr_Cl,    k_HOBr_Br
    # REAL(dp) :: k_HOBr_HSO3m, k_HOBr_SO3mm, k_tot
    # REAL(dp) :: srMw
    #
    k        = 0.0 
    brIce    = 0.0 
    brLiq    = 0.0 
    gammaIce = 0.0 
    gammaLiq = 0.0 
    srMw     = SR_MW(ind_HOBr)
    #
    if H%stratBox  
    #
    # Uptake on tropospheric (origin) sulfate in stratosphere
        gammaLiq = 0.25 
        k = k + Ars_L1k( H%xArea(SUL), H%xRadi(SUL), gammaLiq, srMw )
        #
        # Uptake on strat sulfate liquid aerosol
        k = k + H%xArea(SLA) * H%KHETI_SLA(HOBr_plus_HBr)
        #
        # Uptake on irregular ice cloud
        gammaIce = 0.3 
        if H%natSurface 
            gammaIce = 0.001 
        end
        k = k + Ars_L1k( H%xArea(IIC), H%xRadi(IIC), gammaIce, srMw )

    else
            #
            # HOBr + HBr rxn probability in tropospheric liquid cloud
                Gam_HOBr_CLD(H,gammaLiq,k_tot,k_HOBr_Cl, k_HOBr_Br, k_HOBr_HSO3m, k_HOBr_SO3mm)
            #
            # Branching ratio for liquid path of HOBr + HBr
            branch_0 = ( k_HOBr_Cl + k_HOBr_Br ) / k_tot
            branch   = branch_0 * 0.9 
            if H%Br_over_Cl_Cld <=  5.0e-4 
                branch = branch_0 * Br2_Yield( H%Br_over_Cl_Cld )
            end
            brLiq = branch * H%frac_Br_CldG
            #
            # Overall probability of HOBr uptake and
            # ice-path branching ratio for HOBr + HBr
            CALL Gam_HOBr_Ice( H, gammaIce, dummy, brIce )
            #
            # Compute overall HOBr removal rate in cloud
            k = k + CloudHet( H, srMw, gammaLiq, gammaIce, brLiq, brIce )
            #
    end

    # Assume HOBr is limiting, so update the removal rate accordingly
    k = kIIR1Ltd( C(ind_HOBr), C(ind_HBr), k )
    return k 
end


function HOBrUptkByHCl( H ) 
    # Computes the uptake rate [1/s] for the HOBr + HCl reaction
    # which only occurs in the stratosphere.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    # REAL(dp) :: branch,       branch_0,     brIce
    # REAL(dp) :: brLiq,        dummy,        gammaLiq
    # REAL(dp) :: gammaIce,     k_HOBr_Cl,    k_HOBr_Br
    # REAL(dp) :: k_HOBr_HSO3m, k_HOBr_SO3mm, k_tot
    # REAL(dp) :: srMw
    #
    k        = 0.0 
    brIce    = 0.0 
    brLiq    = 0.0 
    gammaIce = 0.0 
    gammaLiq = 0.0 
    srMw     = SR_MW(ind_HOBr)
    #
    if H%stratBox 
        
        #
        # Uptake on tropospheric (origin) sulfate in stratosphere
        gammaLiq = 0.2 
        k = k + Ars_L1k( H%xArea(SUL), H%xRadi(SUL), gammaLiq, srMw )
        #
        # Uptake on strat sulfate liquid aerosol
        k = k + H%xArea(SLA) * H%KHETI_SLA(HOBr_plus_HCl)
        #
        # Uptake on irregular ice cloud
        gammaIce = 0.3 
        if H%natSurface 
            gammaIce = 0.1 
        end
        k = k + Ars_L1k( H%xArea(IIC), H%xRadi(IIC), gammaIce, srMw )
        #
    else
        #
        # HOBr + HBr rxn probability in tropospheric liquid cloud
        Gam_HOBr_CLD(H,gammaLiq,k_tot,k_HOBr_Cl, k_HOBr_Br, k_HOBr_HSO3m, k_HOBr_SO3mm)
        #
        # Branching ratio for liquid path of HOBr + HCl
        branch_0 = ( k_HOBr_Cl + k_HOBr_Br ) / k_tot
        branch   = branch_0 * 0.1 
        if H%Br_over_Cl_Cld <= 5.0e-4 
            branch = branch_0 * ( 1.0  - Br2_Yield( H%Br_over_Cl_Cld ) )
        end 
 
        brLiq = branch * H%frac_Cl_CldG
        #
        # Overall probability of HOBr uptake and
        # ice-path branching ratio for HOBr + HCl
        Gam_HOBr_Ice( H, gammaIce, brIce, dummy )
        #
        # Compute overall HOBr removal rate in cloud
        k = k + CloudHet( H, srMw, gammaLiq, gammaIce, brLiq, brIce )
        #
    end

    # Assume HOBr is limiting, so update the removal rate accordingly
    k = kIIR1Ltd( C(ind_HOBr), C(ind_HCl), k )
    return k 
end

function HOBrUptkByBrSALA( H )  
    #
    # Computes the uptake rate [1/s] for the HOBr + BrSALA reaction.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    # REAL(dp) :: area,         branch,    branch_0
    # REAL(dp) :: brLiq,        gammaAer,  gammaLiq
    # REAL(dp) :: k_HOBr_Cl,    k_HOBr_Br, k_HOBr_HSO3m
    # REAL(dp) :: k_HOBr_SO3mm, k_tot,     srMw
    #
    k        = 0.0 
    brLiq    = 0.0 
    gammaAer = 0.0 
    gammaLiq = 0.0 
    srMw     = SR_MW(ind_HOBr)
    #
    if H%stratBox
        ifcontinuelog = 1 #continue
    else  
        #
        # HOBr + HBr rxn probability in tropospheric liquid cloud
        Gam_HOBr_CLD(H,gammaLiq,k_tot,k_HOBr_Cl, k_HOBr_Br, k_HOBr_HSO3m, k_HOBr_SO3mm)
        #
        # Branching ratio for liquid path of HOBr + BrSALA in cloud
        branch_0 = ( k_HOBr_Cl + k_HOBr_Br ) / k_tot
        branch   = branch_0 * 0.9 
        if H%Br_over_Cl_Cld <= 5.0e-4   
            branch = branch_0 * Br2_Yield( H%Br_over_Cl_Cld )
        end
        brLiq = branch * H%frac_Br_CldA
        #
        # Compute overall HOBr removal rate in cloud
        k = k + CloudHet( H, srMw, gammaLiq, 0.0 , brLiq, 0.0  )
        #
    end
    #
    # Now consider HOBr uptake by acidic BrSALA in clear-sky
    if H%SSA_is_Acid  
       #
       # Uptake probability [1]
        Gam_HOBr_Aer( H,H%aClRadi,H%H_conc_SSA,H%Cl_conc_SSA, H%Br_conc_SSA, gammaAer)
        #
        # Branching ratio (depends on Br- / Cl- ratio)
        branch = 0.9 
        if H%Br_over_Cl_SSA <= 5.0e-4 
            branch = Br2_Yield( H%Br_over_Cl_SSA )
        end
        #
        # Uptake rate [1/s]
        area = H%ClearFr * H%aClArea
        k    = k + Ars_L1K( area, H%aClRadi, gammaAer, srMw ) * branch
    end

    # Assume HOBr is limiting, so update the removal rate accordingly
    k = kIIR1Ltd( C(ind_HOBr), C(ind_BrSALA), k )
    return k 
end

function HOBrUptkByBrSALC( H )  
    #
    # Computes the uptake rate [1/s] for the HOBr + BrSALC reaction.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    # REAL(dp) :: area,         branch,    branch_0
    # REAL(dp) :: brLiq,        gammaAer,  gammaLiq
    # REAL(dp) :: k_HOBr_Cl,    k_HOBr_Br, k_HOBr_HSO3m
    # REAL(dp) :: k_HOBr_SO3mm, k_tot,     srMw
    #
    k        = 0.0 
    brLiq    = 0.0 
    gammaAer = 0.0 
    gammaLiq = 0.0 
    srMw     = SR_MW(ind_HOBr)
    #
    if H%stratBox
        ifcontinuelog = 1 #continue
    else
        # HOBr + HBr rxn probability in tropospheric liquid cloud
        Gam_HOBr_CLD(H,gammaLiq,k_tot,k_HOBr_Cl, k_HOBr_Br, k_HOBr_HSO3m, k_HOBr_SO3mm )
        #
        # Branching ratio for liquid path of HOBr + BrSALC in cloud
        branch_0 = ( k_HOBr_Cl + k_HOBr_Br ) / k_tot
        branch   = branch_0 * 0.9 
        if H%Br_over_Cl_Cld <= 5.0e-4 
            branch = branch_0 * Br2_Yield( H%Br_over_Cl_Cld )
        end
        brLiq = branch * H%frac_Br_CldC
        #
        # Compute overall HOBr removal rate in cloud
        k = k + CloudHet( H, srMw, gammaLiq, 0.0 , brLiq, 0.0  )
        #
    end
    #
    # Now consider HOBr uptake by acidic BrSALC in clear-sky
    if H%SSC_is_Acid  
        #
        # Uptake probability [1]
        Gam_HOBr_Aer( H,H%xRadi(SSC),H%H_conc_SSC,H%Cl_conc_SSC, H%Br_conc_SSC, gammaAer)
        #
        # Branching ratio (depends on Br- / Cl- ratio)
        branch = 0.9 
        if H%Br_over_Cl_SSC <= 5.0e-4 
            branch = Br2_Yield( H%Br_over_Cl_SSC )
        end
        #
        # Uptake rate [1/s]
        area = H%ClearFr * H%xArea(SSC)
        k    = k + Ars_L1K( area, H%xRadi(SSC), gammaAer, srMw ) * branch
    end

    # Assume HOBr is limiting, so update the removal rate accordingly
    k = kIIR1Ltd( C(ind_HOBr), C(ind_BrSALC), k )
    return k 
end

function HOBrUptkBySALACL( H ) 
    #
    # Computes the uptake rate [1/s] for the HOBr + SALACL reaction.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    # REAL(dp) :: area,         branch,    branch_0
    # REAL(dp) :: brLiq,        gammaAer,  gammaLiq
    # REAL(dp) :: k_HOBr_Cl,    k_HOBr_Br, k_HOBr_HSO3m
    # REAL(dp) :: k_HOBr_SO3mm, k_tot,     srMw
    #
    k        = 0.0 
    brLiq    = 0.0 
    gammaAer = 0.0 
    gammaLiq = 0.0 
    srMw     = SR_MW(ind_HOBr)
    #
    if H%stratBox 
        ifcontinuelog = 1 #continue
    else

        #
        # HOBr + HBr rxn probability in tropospheric liquid cloud
        Gam_HOBr_CLD(H,gammaLiq,k_tot,k_HOBr_Cl, k_HOBr_Br, k_HOBr_HSO3m, k_HOBr_SO3mm )
        #
        # Branching ratio for liquid path of HOBr + SALACL in cloud
        branch_0 = ( k_HOBr_Cl + k_HOBr_Br ) / k_tot
        branch   = branch_0 * 0.1 
        if H%Br_over_Cl_Cld <= 5.0e-4   
            branch = branch_0 * ( 1.0  - Br2_Yield( H%Br_over_Cl_Cld )  )
        end
        brLiq = branch * H%frac_Cl_CldA
        #
        # Compute overall HOBr removal rate in cloud
        k = k + CloudHet( H, srMw, gammaLiq, 0.0 , brLiq, 0.0  )
        #
    end
    #
    # Now consider HOBr uptake by acidic SALACL in clear-sky
    if H%SSA_is_Acid  
        #
        # Uptake probability [1]
        Gam_HOBr_Aer( H,H%aClRadi, H%H_conc_SSA, H%Cl_conc_SSA, H%Br_conc_SSA, gammaAer )
        #
        # Branching ratio (depends on Br- / Cl- ratio)
        branch = 0.1 
        if H%Br_over_Cl_SSA <= 5.0e-4 
            branch = 1.0  - Br2_Yield( H%Br_over_Cl_SSA )
        end
        #
        # Uptake rate [1/s]
        area = H%ClearFr * H%aClArea
        k    = k + Ars_L1K( area, H%aClRadi, gammaAer, srMw ) * branch
    end

    # Assume HOBr is limiting, so update the removal rate accordingly
    k = kIIR1Ltd( C(ind_HOBr), C(ind_SALACL), k )
    return k 
end

function HOBrUptkBySALCCL( H )  
    #
    # Computes the uptake rate [1/s] for the HOBr + SALCCL reaction.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    # REAL(dp) :: area,         branch,    branch_0
    # REAL(dp) :: brLiq,        gammaAer,  gammaLiq
    # REAL(dp) :: k_HOBr_Cl,    k_HOBr_Br, k_HOBr_HSO3m
    # REAL(dp) :: k_HOBr_SO3mm, k_tot,     srMw
    #
    k        = 0.0 
    brLiq    = 0.0 
    gammaAer = 0.0 
    gammaLiq = 0.0 
    srMw     = SR_MW(ind_HOBr)
 
    if H%stratBox
        ifcontinuelog = 1 #continue
    else

        #
        # HOBr + HBr rxn probability in tropospheric liquid cloud
        Gam_HOBr_CLD(H,gammaLiq,k_tot,k_HOBr_Cl, k_HOBr_Br, k_HOBr_HSO3m, k_HOBr_SO3mm )
        #
        # Branching ratio for liquid path of HOBr + SALACL in cloud
        branch_0 = ( k_HOBr_Cl + k_HOBr_Br ) / k_tot
        branch   = branch_0 * 0.1 
        if H%Br_over_Cl_Cld <= 5.0e-4  
            branch = branch_0 * ( 1.0  -  Br2_Yield( H%Br_over_Cl_Cld ) )
        end
        brLiq = branch * H%frac_Cl_CldC
        #
        # Compute overall HOBr removal rate in cloud
        k = k + CloudHet( H, srMw, gammaLiq, 0.0 , brLiq, 0.0  )
       #
    end
    #
    # Now consider HOBr uptake by acidic SALCCL in clear-sky
    if H%SSC_is_Acid 
        #
        # Uptake probability [1]
        Gam_HOBr_Aer( H,H%xRadi(SSC),H%H_conc_SSC, H%Cl_conc_SSC, H%Br_conc_SSC, gammaAer            )
        #
        # Branching ratio (depends on Br- / Cl- ratio)
        branch = 0.1 
        if H%Br_over_Cl_SSC <= 5.0e-4  
            branch = 1.0  - Br2_Yield( H%Br_over_Cl_SSC )
        end
        #
        # Uptake rate [1/s]
        area = H%ClearFr * H%xArea(SSC)
        k    = k + Ars_L1K( area, H%xRadi(SSC), gammaAer, srMw ) * branch
    end

        # Assume HOBr is limiting, so update the removal rate accordingly
        k = kIIR1Ltd( C(ind_HOBr), C(ind_SALCCL), k )
    return k 
end

function HOBrUptkByHSO3m( H )  
    #
    # Computes the uptake rate [1/s] for the HOBr + HSO3(-) reaction.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    # REAL(dp) :: brLiq,     gammaLiq,     k_HOBr_Cl
    # REAL(dp) :: k_HOBr_Br, k_HOBr_HSO3m, k_HOBr_SO3mm
    # REAL(dp) :: k_tot,     srMw
    #
    k        = 0.0 
    brLiq    = 0.0 
    gammaLiq = 0.0 
    srMw     = SR_MW(ind_HOBr)
    #
    if H%stratBox
        ifcontinuelog = 1 #continue
    else
        #
        # HOBr + HBr rxn probability in tropospheric liquid cloud
        Gam_HOBr_CLD(                                                    &
            H,         gammaLiq,  k_tot,                                     &
            k_HOBr_Cl, k_HOBr_Br, k_HOBr_HSO3m, k_HOBr_SO3mm )
        #
        # Branching ratio for liquid path of HOBr + HSO3- in cloud
        brLiq = k_HOBr_HSO3m / k_tot
        #
        # Compute overall HOBr removal rate in cloud
        k = k + CloudHet( H, srMw, gammaLiq, 0.0 , brLiq, 0.0  )
        #
    end

    # Assume HOBr is limiting, so update the removal rate accordingly
    # Convert SO2 to HSO3- with the HSO3-/SO2 ratio
    k = kIIR1Ltd( C(ind_HOBr), C(ind_SO2), k ) * H%HSO3m
    return k 
end

function HOBrUptkBySO3mm( H )  
    #
    # Computes the uptake rate [1/s] for the HOBr + HSO3(-) reaction.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    # REAL(dp) :: brLiq,     gammaLiq,     k_HOBr_Cl
    # REAL(dp) :: k_HOBr_Br, k_HOBr_HSO3m, k_HOBr_SO3mm
    # REAL(dp) :: k_tot,     srMw
    #
    k    = 0.0 
    srMw = SR_MW(ind_HOBr)
    #
    if H%stratBox
        ifcontinuelog = 1 #continue
    else
       #
       # HOBr + HBr rxn probability in tropospheric liquid cloud
       Gam_HOBr_CLD( H,gammaLiq,k_tot,k_HOBr_Cl, k_HOBr_Br, k_HOBr_HSO3m, k_HOBr_SO3mm )
       #
       # Branching ratio for liquid path of HOBr + SO3-- in cloud
       brLiq = k_HOBr_SO3mm / k_tot
       #
       # Compute overall HOBr removal rate in cloud
       k = k + CloudHet( H, srMw, gammaLiq, 0.0 , brLiq, 0.0  )
       #
    end

    # Assume HOBr is limiting, so update the removal rate accordingly
    # Convert SO2 to SO3-- with the SO3--/SO2 ratio
    k = kIIR1Ltd( C(ind_HOBr), C(ind_SO2), k ) * H%SO3mm

    return k 
end

# =========================================================================
#  Hetchem rate-law functions for HOCl
# =========================================================================

function Gam_HOCl_Cld( H, gamma, branchCl, branchSO3 )
    #
    # Computes the uptake coefficient [1] of HOCl + HSO3 and HOCl + SO3.
    # Returns reaction rates of Cl and SO3 paths to compute branching ratios.
    #
    # TYPE(HetState), INTENT(IN)  :: H             # Hetchem State
    # REAL(dp),       INTENT(OUT) :: gamma         # Rxn prob [1]
    # REAL(dp),       INTENT(OUT) :: branchCl      # Branch ratio, Cl path [1]
    # REAL(dp),       INTENT(OUT) :: branchSO3     # Branch ratio, SO3 path [1]
    #
    INV_AB = 1.0  / 0.8   # 1/mass accom coeff
    D_l    = 2.0e-5         # Liq diff phase coeff
    #
    # REAL(dp) :: cavg,  gb_tot, H_X, k_Cl
    # REAL(dp) :: k_SO3, k_tot,  l_r, M_X
    #
    # Reaction rates, Cl and SO3 paths [1/s]
    k_Cl      = 1.5e+4  * H%H_Conc_LCL  * H%Cl_conc_Cld
    k_SO3     = 2.8e+5  * H%TSO3_aq
    k_tot     = k_Cl + k_SO3
    #
    # Compute reactive uptake coefficient [1] and branching ratio [1], Cl path
    # but avoid division by zero
    gamma     = 0.0 
    branchCl  = 0.0 
    branchSO3 = 0.0 
    if k_tot > 0.0 
       #
       # thermal velocity (cm/s)
       M_X       = MW(ind_HOCl) * 1.0e-3 
       cavg      = sqrt( EIGHT_RSTARG_T / ( H%Pi * M_X ) ) * 100.0 
       #
       # Henry's law
       H_X       = ( HENRY_k0(ind_HOCl) * CON_ATM_BAR )  * exp( HENRY_CR(ind_HOCl) * ( INV_temp - INV_T298 ) )
       #
       l_r       = sqrt( D_l / k_tot )
       gb_tot    = FOUR_R_T * H_X * l_r * k_tot / cavg
       gb_tot    = gb_tot * ReactoDiff_Corr( H%rLiq, l_r )
       #
       gamma     = 1.0  / ( INV_AB + 1.0  / gb_tot )
       branchCl  = k_Cl   / k_tot
       branchSO3 = k_SO3  / k_tot
    end
end 

function Gam_HOCl_AER( H, radius, C_Hp, C_Cl, gamma )
    #
    # Calculates reactive uptake coefficients [1] for the reactions
    # HOCl + SALACL and HOCl + SALCCL.
    #
    # TYPE(HetState), INTENT(IN)  :: H             # Hetchem State
    # REAL(dp),       INTENT(IN)  :: radius        # Radius [cm]
    # REAL(dp),       INTENT(IN)  :: C_Hp          # H+ conc [mol/L]
    # REAL(dp),       INTENT(IN)  :: C_Cl          # Cl- conc [mol/L]
    # REAL(dp),       INTENT(OUT) :: gamma
    #
    INV_AB = 1.0  / 0.8  # 1/mass accum coeff
    D_l    = 2.0e-5        # Liq phase diffusion coeff
    K_TER  = 1.5e+4        # Units: M-1 s-1
    #
    # REAL(dp) :: cavg, gb, H_X, l_r, M_X
    #
    gamma = 0.0 
    #
    # If C_Cl is zero, gamma is zero#
    
    if C_Cl > 0.0
        #
        # Thermal velocity (cm/s)
        M_X   = MW(ind_HOCl) * 1.0e-3 
        cavg  = sqrt( EIGHT_RSTARG_T / ( H%PI * M_X ) ) * 100.0 
        #
        # Henry's law
        H_X   = ( HENRY_k0(ind_HOCl) * CON_ATM_BAR ) * exp( HENRY_CR(ind_HOCl) * ( INV_temp - INV_T298 ) )
        l_r   = sqrt( D_l / ( K_TER * C_Hp * C_Cl ) )
        gb    = FOUR_R_T * H_X * l_r * K_TER * C_Hp * C_Cl / cavg
        gb    = gb * ReactoDiff_Corr( radius, l_r )
        #
        # Reactive uptake coefficient [1]
        gamma = 1.0  / ( INV_AB  +  1.0  / gb )
    end
end

function HOClUptkByHCl( H ) RESULT( k )
    #
    # Computes the uptake rate [1/s] for the HOCl + HBr reaction.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    # REAL(dp) :: branch, branchCl, brIce, dummy
    # REAL(dp) :: gamma,  gammaIce, srMw
    #
    k    = 0.0 
    srMw = SR_MW(ind_HOCl)
    #
    if  H%stratBox
        #
        # HOCl + HBr on tropospheric sulfate in stratosphere
        gamma = 0.8 
        k = k + Ars_L1K( H%xArea(SUL), H%xRadi(SUL), gamma, srMw )
        #
        # HOCl + HBr on stratopsheric liquid aerosol
        k = k + H%xArea(SLA) * H%KHETI_SLA(HOCl_plus_HCl)
        #
        # HOCl + HBr on irregular ice cloud
        gamma = 0.2                           # Rxn prob, ice
        if H%natSurface 
            gamma = 0.1       # Rxn prob, NAT
        end
        k = k + Ars_L1K( H%xArea(IIC), H%xRadi(IIC), gamma, srMw )
    else
        #
        # HOCl + HCl uptake coeff [1] & branch ratio [1] in trop liquid cloud
        CALL Gam_HOCl_Cld( H, gamma, branchCl, dummy )
        branch = branchCl * H%frac_Cl_CldG
        #
        # HOCl + HCl uptake coeff [1] & branch ratio [1] in trop ice cloud
        gammaIce = 0.22  * H%HCl_theta
        brIce    = 1.0 
        #
        # Compute overall HOCl + HCl uptake rate accounting for cloud fraction
        k = k + CloudHet( H, srMw, gamma, gammaIce, branch, brIce )
    end
    #
    # Assume HOCl is limiting, so recompute reaction rate accordingly
    k = kIIR1Ltd( C(ind_HOCl), C(ind_HCl), k )
    return k 
end

function HOClUptkByHBr( H )  
    #
    # Computes the uptake rate [1/s] for the HOCl + HBr reaction.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    # REAL(dp) :: gamma, srMw
    #
    k    = 0.0 
    srMw = SR_MW(ind_HOCl)
    #
    if H%stratBox  
        #
        # HOCl + HBr on tropospheric sulfate in stratosphere
        gamma = 0.8 
        k = k + Ars_L1K( H%xArea(SUL), H%xRadi(SUL), gamma, srMw )
        #
        # HOCl + HBr on stratopsheric liquid aerosol
        k = k + H%xArea(SLA) * H%KHETI_SLA(HOCl_plus_HBr)
        #
        # HOCl + HBr on irregular ice cloud (ice and NAT surface)
        gamma = 0.3 
        k = k + Ars_L1K( H%xArea(IIC), H%xRadi(IIC), gamma, srMw )
    end
    #
    # Assume HOCl is limiting, so recompute reaction rate accordingly
    k = kIIR1Ltd( C(ind_HOCl), C(ind_HBr), k )
    return k 
end

function HOClUptkBySALACL( H ) 
    #
    # Computes the uptake rate [1/s] for the HOCl + SALACL reaction.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    # REAL(dp) :: area,  branch, branchCl
    # REAL(dp) :: dummy, gamma,  srMw
    #
    k    = 0.0 
    srMw = SR_MW(ind_HOCl)
    #
    # Compute HOCl + SALACL uptake in tropospheric liquid cloud
    if H%stratBox
        ifcontinuelog = 1 #continue
    else
        #
        # HOCl + SALACL uptake coeff [1] & branch ratio [1], liquid path
        Gam_HOCl_Cld( H, gamma, branchCl, dummy )
        branch = branchCl * H%frac_Cl_CldA
        #
        # HOCl + HCl uptake rate [1/s] accounting for cloud fraction
        k = k + CloudHet( H, srMw, gamma, 0.0 , branch, 0.0  )
    end
    #
    # Compute HOCl + SALACL uptake rate [1/s] on acidic aerosols in clear-sky

    if H%SSA_is_Acid  
        CALL Gam_HOCl_Aer( H, H%aClRadi, H%H_conc_SSA, H%Cl_conc_SSA, gamma )
        area = H%ClearFr * H%aClArea
        k    = k + Ars_L1k( area, H%aClRadi, gamma, srMw )
    end
    #
    # Assume HOCl is limiting, so recompute reaction rate accordingly
    k = kIIR1Ltd( C(ind_HOCl), C(ind_SALACL), k )
    return k 
end

function HOClUptkBySALCCL( H )  
    #
    # Computes the uptake rate [1/s] for the HOCl + SALCCL reaction.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    # REAL(dp) :: area,  branch, branchCl
    # REAL(dp) :: dummy, gamma,  srMw
    #
    k    = 0.0 
    srMw = SR_MW(ind_HOCl)
    #
    # Compute HOCl + SALACL uptake in tropospheric liquid cloud
    if H%stratBox
        ifcontinuelog = 1 #continue
    else
       #
       # HOCl + SALACL uptake coeff [1] & branch ratio [1], liquid path
       Gam_HOCl_Cld( H, gamma, branchCl, dummy )
       branch = branchCl * H%frac_Cl_CldC
       #
       # HOCl + HCl uptake rate [1/s] accounting for cloud fraction
       k = k + CloudHet( H, srMw, gamma, 0.0 , branch, 0.0  )
    end 
    #
    # Compute HOCl + SALCCL uptake rate [1/s] on acidic aerosols in clear-sky
    if  H%SSC_is_Acid  
       Gam_HOCl_Aer( H, H%xRadi(SSC), H%H_conc_SSC, H%Cl_conc_SSC, gamma )
       area = H%ClearFr * H%xArea(SSC)
       k    = k + Ars_L1k( area, H%xRadi(SSC), gamma, srMw )
    end
    #
    # Assume HOCl is limiting, so recompute reaction rate accordingly
    k = kIIR1Ltd( C(ind_HOCl), C(ind_SALCCL), k )
    return k 
end

function HOClUptkByHSO3m( H )  
    #
    # Computes the uptake rate [1/s] for the HOCl + HSO3- reaction.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    # REAL(dp) :: area,  branch, branchSO3
    # REAL(dp) :: dummy, gamma,  srMw
    #
    k    = 0.0 
    srMw = SR_MW(ind_HOCl)
    #
    # Compute HOCl + HSO3- uptake in tropospheric liquid cloud
    if H%stratBox
        ifcontinuelog = 1 #continue
    else
       # HOCl + HSO3- uptake coeff [1] & branch ratio [1], liquid path
       Gam_HOCl_Cld( H, gamma, dummy, branchSO3 )
       branch = branchSO3 * H%frac_HSO3_aq
       #
       # HOCl + HSO3- uptake rate [1/s] accounting for cloud fraction
       k = k + CloudHet( H, srMw, gamma, 0.0 , branch, 0.0  )
    end
    #
    # Assume HOCl is limiting, so recompute reaction rate accordingly
    # Convert SO2 to HSO3- with the HSO3-/SO2 ratio
    k = kIIR1Ltd( C(ind_HOCl), C(ind_SO2), k ) * H%HSO3m
    return k 
end

function HOClUptkBySO3mm( H ) 
    # Computes the uptake rate [1/s] for the HOCl + SO3-- reaction.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    # REAL(dp) :: area,  branch, branchSO3
    # REAL(dp) :: dummy, gamma,  srMw
    #
    k    = 0.0 
    srMw = SR_MW(ind_HOCl)
    #
    # Compute HOCl + SO3-- uptake in tropospheric liquid cloud
    if H%stratBox
        ifcontinuelog = 1 #continue
    else
        # HOCl + SO3-- uptake coeff [1] & branch ratio [1], liquid path
        CALL Gam_HOCl_Cld( H, gamma, dummy, branchSO3 )
        branch = branchSO3 * H%frac_SO3_aq
        #
        # HOCl + SO3-- uptake rate [1/s] accounting for cloud fraction
        k = k + CloudHet( H, srMw, gamma, 0.0 , branch, 0.0  )
    end
    #
    # Assume HOCl is limiting, so recompute reaction rate accordingly
    # Convert SO2 to SO3-- with the SO3--/SO2 ratio
    k = kIIR1Ltd( C(ind_HOCl), C(ind_SO2), k ) * H%SO3mm
    return k 
end


# =========================================================================
# Hetchem rate-law functions for iodine species
# (HI, HOI, I2O2, I2O3, I2O4, IONO2, IONO3)
# =========================================================================

function IuptkBySulf1stOrd( srMw, gamma, H ) 
    #
    # Computes the reaction rate [1/s] for uptake of iodine species
    # by sulfate (aerosol types #8 and #13).
    #
    #   REAL(dp),       INTENT(IN) :: srMw, gamma    # sqrt( mol wt ), rxn prob
    #   TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    #   REAL(dp)                   :: k              # rxn rate [1/s]
    #
    # Uptake rate of iodine by tropospheric sulfate
    k = Ars_L1k( H%xArea(SUL), H%xRadi(SUL), gamma, srMw )

    # For UCX-based mechanisms also allow reaction on stratospheric
    # sulfate liq aerosol if tropospheric sulfate is requested
    k = k + Ars_L1k( H%xArea(SLA), H%xRadi(SLA), gamma, srMw )
    return k 
end

function IuptkBySALA1stOrd( srMw, gamma, H )
    #
    # Computes the reaction rate [1/s] for uptake of iodine species
    # by accumulation-mode (aka fine) sea-salt aerosol.
    #
    # REAL(dp),       INTENT(IN) :: srMw, gamma    # sqrt( mol wt ) rxn prob
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    k = Ars_L1k( H%xArea(SSA), H%xRadi(SSA), gamma, srMw )
    return k
end

function IuptkByAlkSALA1stOrd( srMw, gamma, H )  
    #
    # Computes the reaction rate [1/s] for uptake of iodine species
    # by alkaline accumulation-mode (aka fine) sea-salt aerosol.
    #
    # REAL(dp),       INTENT(IN) :: srMw, gamma    # sqrt( mol wt ) rxn prob
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    k = 0.0 
    if H%SSA_is_Alk  
       k = IuptkBySALA1stOrd( srMw, gamma, H )
    end 

    return k
end

function IuptkBySALC1stOrd( srMw, gamma, H )  
    #
    # Computes the reaction rate [1/s] for uptake of iodine species
    # by coarse-mode sea-salt aerosol.
    #
    # REAL(dp),       INTENT(IN) :: srMw, gamma    # sqrt( mol wt ), rxn prob
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    k = Ars_L1k( H%xArea(SSC), H%xRadi(SSC), gamma, srMw )
    return k
end

function IuptkByAlkSALC1stOrd( srMw, gamma, H )  
    #
    # Computes the reaction rate [1/s] for uptake of iodine species
    # by alkaline coarse-mode sea-salt aerosol.
    #
    # REAL(dp),       INTENT(IN) :: srMw, gamma    # sqrt( mol wt ), rxn prob
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    k = 0.0 
    if H%SSC_is_Alk 
        k = IuptkBySALC1stOrd( srMw, gamma, H )
    end
    return k
end

function IbrkdnByAcidBrSALA( srMw, conc, gamma, H )  
    #
    # Breakdown of iodine species on acidic sea-salt (accumulation mode)
    # Assume a ratio of IBr:ICl = 0.15:0.85
    #
    # REAL(dp),       INTENT(IN) :: srMw, conc, gamma
    # TYPE(HetState), INTENT(IN) :: H
    # REAL(dp)                   :: k
    #
    k = 0.0 
    if H%SSA_is_Acid 
       k = 0.15  * IuptkBySALA1stOrd( srMw, gamma, H )
       k = kIIR1Ltd( conc, C(ind_BrSALA), k ) # conc is limiting, so update k
    end
    return k
end

function IbrkdnByAcidBrSALC( srMw, conc, gamma, H ) 
    #
    # Breakdown of iodine species on acidic sea-salt (accumulation mode)
    # Assume a ratio of IBr:ICl = 0.15:0.85
    #
    # REAL(dp),       INTENT(IN) :: srMw, conc, gamma
    # TYPE(HetState), INTENT(IN) :: H
    # REAL(dp)                   :: k
    #
    k = 0.0 
    if  H%SSC_is_Acid  
       k = 0.15  * IuptkBySALC1stOrd( srMw, gamma, H )
       k = kIIR1Ltd( conc, C(ind_BrSALC), k ) # conc is limiting, so update k
    end
    return k
end

function IbrkdnByAcidSALACl( srMw, conc, gamma, H ) RESULT( k )
    #
    # Breakdown of iodine species on acidic sea-salt (accumulation mode)
    # Assume a ratio of IBr:ICl = 0.15:0.85
    #
    # REAL(dp),       INTENT(IN) :: srMw, conc, gamma
    # TYPE(HetState), INTENT(IN) :: H
    # REAL(dp)                   :: k
    #
    k = 0.0 
    if H%SSA_is_Acid  
       k = 0.85  * IuptkBySALA1stOrd( srMw, gamma, H )
       k = kIIR1Ltd( conc, C(ind_SALACl), k ) # conc is limiting, so update k
    end
    return k
end

function IbrkdnByAcidSALCCl( srMw, conc, gamma, H )  
    #
    # Breakdown of iodine species on acidic sea-salt (accumulation mode)
    # Assume a ratio of IBr:ICl = 0.15:0.85
    #
    # REAL(dp),       INTENT(IN) :: srMw, conc, gamma
    # TYPE(HetState), INTENT(IN) :: H
    # REAL(dp)                   :: k
    #
    k = 0.0 
    if   H%SSC_is_Acid  
        k = 0.85  * IuptkBySALC1stOrd( srMw, gamma, H )
        k = kIIR1Ltd( conc, C(ind_SALCCl), k ) # conc is limiting, so update k
    end
    return k
end

function IONO2uptkByH2O( H ) 
    #
    # Computes the reaction rate [1/s] for IONO2 + H2O = HOI + HNO3
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # Rxn rate [1/s]
    #
    # REAL(dp) :: area, conc, gamma, srMw
    #
    k    = 0.0 
    srMw = SR_MW(ind_IONO2)
    #
    # Tropopsheric sulfate (use T-dependent gamma, cf Deiber et al 2004)
    area  = H%ClearFr * H%xArea(SUL)
    gamma = max( ( 0.0021  * temp - 0.561  ), 0.0  )
    k     = k + Ars_L1K( area, H%xRadi(SUL), gamma, srMw )
    #
    # Alkaline fine sea salt (use gamma from Sherwen et al 2016)
    if H%SSA_is_Alk  
       area  = H%ClearFr * H%xArea(SSA)
       gamma = 0.01 
       k = k + Ars_L1K( area, H%xRadi(SSA), gamma, srMw )
    end
    #
    # Alkaline coarse sea salt (use gamma from Sherwen et al 2016)
    if H%SSC_is_Alk 
       area  = H%ClearFr * H%xArea(SSC)
       gamma = 0.01 
       k = k + Ars_L1K( area, H%xRadi(SSC), gamma, srMw )
    end
    #
    # Stratospheric liquid aerosol
    k = k + H%xArea(SLA) * H%KHETI_SLA(BrNO3_plus_H2O)
    #
    # Irregular ice cloud
    area  = H%ClearFr * H%xArea(IIC)
    gamma = 0.3 
    if H%natSurface
        gamma = 0.001
    end
    k = k + Ars_L1K( area, H%xRadi(IIC), gamma, srMw )
    #
    # Also account for cloudy grid box
    # Use gamma(liquid) = gamma(ice) = 0.01, to make the uptake coefficient
    # consistent with hydrolysis in aerosols (T. Sherwen, 28 Sep 2021)
    k = k + CloudHet( H, srMw, 0.01 , 0.01 , 1.0 , 1.0  )
    #
    # Assume IONO2 is limiting, so update the reaction rate accordingly
    k = kIIR1Ltd( C(ind_IONO2), C(ind_H2O), k )
    return k
end


# =========================================================================
# Hetchem rate-law functions for N2O5
# =========================================================================

function N2O5uptkByH2O( H )  
    #
    # Set heterogenous chemistry rate for N2O5.
    #
    # TYPE(HetState), INTENT(INOUT) :: H
    # REAL(dp)                      :: k
    #
    # REAL(dp) :: Y_ClNO2, Rp, SA, SA_sum, area, gamma, srMw, ktmp
    #
    k    = 0.0 
    srMw = SR_MW(ind_N2O5)
    #
    # Uptake on mineral dust
    gamma = 0.02 
    k = k + Ars_L1K( H%ClearFr * H%xArea(DU1), H%xRadi(DU1), gamma, srMw )
    k = k + Ars_L1K( H%ClearFr * H%xArea(DU2), H%xRadi(DU2), gamma, srMw )
    k = k + Ars_L1K( H%ClearFr * H%xArea(DU3), H%xRadi(DU3), gamma, srMw )
    k = k + Ars_L1K( H%ClearFr * H%xArea(DU4), H%xRadi(DU4), gamma, srMw )
    k = k + Ars_L1K( H%ClearFr * H%xArea(DU5), H%xRadi(DU5), gamma, srMw )
    k = k + Ars_L1K( H%ClearFr * H%xArea(DU6), H%xRadi(DU6), gamma, srMw )
    k = k + Ars_L1K( H%ClearFr * H%xArea(DU7), H%xRadi(DU7), gamma, srMw )
    #
    # Uptake on tropospheric sulfate
    # Reduce the rate of the HNO3 pathway in accordinace with
    # the ClNO2 yield on SNA + ORG aerosol
    # Reduce ClNO2 yield by 75% (cf McDuffie et al, JGR, 2018)
    N2O5_InorgOrg(H,H%AClVol, H%xVol(ORC), H%xH2O(SUL),  H%xH2O(ORC), H%AClRadi, C(ind_NIT), C(ind_SALACL), gamma,Y_ClNO2, Rp, SA   )
    #
    ktmp = Ars_L1K( H%ClearFr * SA, Rp, gamma, srMw )
    k = k + ktmp - ( ktmp * Y_ClNO2 * 0.25  )
    #
    # Uptake on black carbon
    gamma = 0.005 
    k = k + Ars_L1K( H%ClearFr * H%xArea(BKC), H%xRadi(BKC), gamma, srMw )
    #
    # Uptake on coarse sea salt (aerosol type #12)
    # Reduce the rate of this HNO3 pathway in accordance with the yield
    CALL N2O5_InorgOrg(H, H%xVol(SSC),0.0 , H%xH2O(SSC), 0.0 , H%xRadi(SSC), C(ind_NITs), C(ind_SALCCL), gamma,  Y_ClNO2, Rp,  SA )
    #
    ktmp = Ars_L1k( H%ClearFr * SA, Rp, gamma, srMw )
    k    = k + ktmp - ( ktmp * Y_ClNO2 )
    #
    # Uptake on stratopsheric liquid aerosol
    k = k + H%xArea(SLA) * H%KHETI_SLA(N2O5_plus_H2O)
    #
    # Uptake on irregular ice cloud
    gamma = 0.02                            # Ice
    if H%natSurface 
        gamma = 4.0e-4     # NAT
    end
    k = k + Ars_L1K( H%ClearFr * H%xArea(IIC), H%xRadi(IIC), gamma, srMw )

    # Assume N2O5 is limiting, so update the removal rate accordingly
    k = kIIR1Ltd( C(ind_N2O5), C(ind_H2O), k )
    return k 
end

function N2O5uptkBySALACl( H ) RESULT( k )
    #
    # Computes uptake rate of N2O5 on Cl- in fine sea salt.
    # This reaction follows the N2O5 + Cl- channel.
    #
    TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    REAL(dp)                   :: k              # Rxn rate [1/s]
    REAL(dp) :: gamma, Y_ClNO2, Rp, SA           # local vars
    #
    # Initialize
    k = 0.0 
    #
    # Properties of inorganic (SNA) sea salt coated with organics
    N2O5_InorgOrg(H, H%AClVol,  H%xVol(ORC), H%xH2O(SUL), H%xH2O(ORC), H%aClRadi, C(ind_NIT),  C(ind_SALACL),gamma,Y_ClNO2, Rp,SA )

    # Total loss rate of N2O5 (kN2O5) on SNA+ORG+SSA aerosol.
    # Reduce ClNO2 production yield on fine inorganic+organic
    # aerosol by 75% (cf. McDuffie et al, JGR, 2018).
    k = Ars_L1K( H%ClearFr * SA, Rp, gamma, SR_MW(ind_N2O5)                 )
    k = k * Y_ClNO2 * 0.25 

    # Assume N2O5 is limiting, so update the removal rate accordingly
    k = kIIR1Ltd( C(ind_N2O5), C(ind_SALACL), k )
    return k 
end

function N2O5uptkBySALCCl( H )  
    #
    # Computes uptake rate of N2O5 on Cl- in coarse sea salt.
    # This reaction follows the N2O5 + Cl- channel.
    #
    TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    REAL(dp)                   :: k              # Rxn rate [1/s]
    REAL(dp) :: gamma, Y_ClNO2, Rp, SA           # local vars
    #
    # Initialize
    k = 0.0 
    #
    # Properties of inorganic (SNA) sea salt coated with organics
    N2O5_InorgOrg( H,H%xVol(SSC),0.0 , H%xH2O(SSC), 0.0 , H%xRadi(SSC), C(ind_NITs), C(ind_SALCCL),gamma,Y_ClNO2,Rp,SA)

    # Total loss rate of N2O5 (kN2O5) on SNA+ORG+SSA aerosol
    k = Ars_L1k( H%ClearFr * SA, Rp, gamma, SR_MW(ind_N2O5) )
    k = k * Y_ClNO2

    # Assume N2O5 is limiting, so update the removal rate accordingly
    k = kIIR1Ltd( C(ind_N2O5), C(ind_SALCCL), k )
    return k 
end

function N2O5_InorgOrg( H,      volInorg, volOrg, H2Oinorg,              &
    H2Oorg, Rcore,    NIT,    Cl,                    &
    gamma,  Y_ClNO2,  rp,     areaTotal             )
    #
    # Computes the GAMMA reaction probability for N2O5 loss in inorganic
    # (sulfate-nitrate-ammonium-sea salt) aerosols with organic coatings,
    # based on the recommendation of McDuffie (2018) JGR.
    # The inorganic core is based on Bertram and Thornton ACP (2009).
    #
    # TYPE(HetState), INTENT(IN) :: H    # HetState Object
    # REAL(dp), INTENT(IN)  :: volInorg  # vol of wet inorg aerosol core  [cm3/cm3]
    # REAL(dp), INTENT(IN)  :: volOrg    # vol of wet org aerosol coating [cm3/cm3]
    # REAL(dp), INTENT(IN)  :: H2Oinorg  # vol of H2O in inorg core [cm3/cm3]
    # REAL(dp), INTENT(IN)  :: H2Oorg    # vol of H2O in org coating [cm3/cm3]
    # REAL(dp), INTENT(IN)  :: Rcore     # radius of inorg core [cm]
    # REAL(dp), INTENT(IN)  :: NIT       # aer nitrate conc [molec/cm3]
    # REAL(dp), INTENT(IN)  :: Cl        # aer chloride conc [molecule/cm3]
    # REAL(dp), INTENT(OUT) :: gamma     # [1]
    # REAL(dp), INTENT(OUT) :: Y_ClNO2   # [1]
    # REAL(dp), INTENT(OUT) :: rp        # [cm]
    # REAL(dp), INTENT(OUT) :: areaTotal # [cm2/cm3]
    #
    # Parameters from Bertram and Thornton (2009) ACP and McDuffie 2018 JGR
    KH    = 5.1e+1    #unitless
    k3k2b = 4.0e-2    #unitless
    beta  = 1.15e+6   #[s-1]
    delta = 1.3e-1    #[M-1]

    # Organic Parameters from Antilla et al., 2006 and Riemer et al., 2009
    Haq  = 5e+3  # Henry coef [mol/m3/atm], Antilla
    Daq  = 1e-9  # Aq diff coef [m2/s], Riemer
    ONE_THIRD = 1.0  / 3.0 
    #
    # REAL(dp) :: k2f,      A,        speed,       gamma_core, gamma_coat
    # REAL(dp) :: volTotal, H2Ototal, volRatioDry, l,          eps
    # REAL(dp) :: OCratio,  M_H2O,    M_NIT,       M_Cl,       M_N2O5

    #------------------------------------------------------------------------
    # Concentrations, thickness, etc.
    #------------------------------------------------------------------------

    # Total volume (organic + inorganic), cm3(aerosol)/cm3(air)
    volTotal = volInorg + volOrg

    # Total H2O (organic + inorganic), cm3(H2O)/cm3(air)
    H2Ototal = H2Oinorg + H2Oorg

    # Ratio of inorganic to total (organic+inorganic) volumes when dry, unitless
    volRatioDry = SafeDiv( max( volInorg - H2Oinorg, 0.0  ),               &
    max( volTotal - H2Ototal, 0.0  ), 0.0        )

    # Particle radius, cm
    # Derived from spherical geometry
    # [note: The radius and surface area of a wet particle are
    # properly calculated from the wet volume volume ratio (including water).
    # We use the dry volume ratio here because McDuffie et al. (2018) fitted
    # the N2O5 gamma parameters to field data in a model using the
    # dry ratio. cdholmes 7/22/2019]
    Rp = SafeDiv( Rcore, volRatioDry^ONE_THIRD, Rcore )

    # Coating thickness, cm
    l = Rp - Rcore

    # mean molecular speed [m s-1]
    # sqrt( 8RT / (pi M) )
    M_N2O5  = MW(ind_N2O5) * 1.0e-3 
    speed = sqrt( EIGHT_RSTARG_T / ( H%PI * M_N2O5 ) )

    # Concentrations [mol/L]
    M_H2O = H2Ototal / 18e+0  / volTotal * 1000.0    # H2O
    M_NIT = NIT / volTotal / H%AVO * 1000.0              # Nitrate
    M_Cl  = Cl  / volTotal / H%AVO * 1000.0              # Chloride

    #------------------------------------------------------------------------
    # Gamma for the organic shell (cf McDuffie (2018) JGR)
    #------------------------------------------------------------------------

    #O:C ratio from Eq. 10 of Canagaratna et al., 2015 (ACP)
    # Take average OM/OC ratio from /GeosCore/aerosol_mod.F90
    OCratio = ((( H%OMOC_POA + H%OMOC_OPOA ) / 2.0  ) - 1.17  ) / 1.29 

    # organic scaling factor (eps(Haq*Daq) = Horg*Dorg)
    # from McDuffie (2018) JGR
    eps = 1.5e-1  * OCratio + 1.6e-3  * RELHUM

    # Gamma for coating
    # [Rcore, Rp, and l converted cm -> m here]
    if l <= 0.0e+0 
        gamma_coat = 0.0 
    else
        gamma_coat = ( FOUR_RGASLATM_T * 1.0e-3  * eps * Haq * Daq * Rcore /100.0  )/( speed * l/100.0  * Rp/100.0                                 )
    end

    # Total particle surface area, cm2/cm3
    areaTotal = 3.0  * volTotal / Rp

    #------------------------------------------------------------------------
    # Gamma for the inorganic core
    # Implements recommendations by McDuffie (2018) JGR,
    # following the general approach from Bertram and Thornton ACP (2009).
    #------------------------------------------------------------------------

    # Select dry or deliquesed aerosol based on molar concentration of H2O
    if M_H2O < 0.1 

    # When H2O is nearly zero, use dry aerosol value
        gamma_core = 0.005 

    else

       # mean molecular speed [cm/s]
        speed = speed * 1e+2 

        # A factor from Bertram and Thornton (2009), s
        # Their paper suggested an approximated value of A = 3.2D-8
        A = ( ( 4.0  * volTotal ) / ( speed * areaTotal ) ) * KH

        # Cap A at 3.2D-8 to prevent high gamma values when vol/area is high.
        # See Github issue: https://github.com/geoschem/geos-chem/issues/907
        A = min( A, 3.2e-8  )

        # k2f - reaction rate constant of N2O5 with H2O
        # From McDuffie (2018): k2f = 2.14D5 * H2O
        # This linear water dependence is not accurate at large
        # (>20 M) aerosol water concentrations. Therefore, k2f is
        # calculated following Bertram and Thornton ACP (2009).
        # Eq 11 from Bertram and Thronton (2009):
        # Modified to avoid underflow when exp(-delta*H2O) ~1
        if delta * M_H2O < 1e-2 
            k2f = beta * ( delta * M_H2O )
        else
            k2f = beta * ( 1e+0  - exp( -delta * M_H2O ) )
        end

        # Eq 12 from Bertram and Thornton (2009)
        # Use safe_div to avoid overflow when NIT ~ 0
        gamma_core = A * k2f * (1.0  - 1.0  / (1.0  + SafeDiv(k3k2b*M_H2O, M_NIT, 1.0e+30 )))
    end

    #------------------------------------------------------------------------
    # Gamma for overall uptake
    #------------------------------------------------------------------------
    if gamma_coat <= 0.0 
        gamma = gamma_core
    elseif   gamma_core <= 0.0   
        gamma = 0.0 
    else
        gamma = 1.0  / ( ( 1.0 /gamma_core ) + ( 1.0 /gamma_coat )      )
    end

    #------------------------------------------------------------------------
    # ClNO2 yield
    #------------------------------------------------------------------------

    # Calculate the ClNO2 yield following Bertram and Thornton 2009 ACP
    Y_ClNO2 = ClNO2_BT( M_Cl, M_H2O )
end

function ClNO2_BT( Cl, H2O )  
    #
    # Computes the PHI production yield of ClNO2 from N2O5 loss
    # in sulfate-nitrate-ammonium (SNA) aerosols based on the
    # recommendation of Bertram and Thornton (2009) ACP.
    #
    # REAL(dp), INTENT(IN) :: Cl, H2O  # [mol/L]
    # REAL(dp)             :: PHI
    # REAL(dp), PARAMETER  :: k2k3  = 1.0  / 4.5e+2   # BT 2009
    #
    # When H2O is nearly zero, assign phi accordingly and exit
    if H2O < 0.1   
        phi = 0.0 
        if Cl > 1e-3  
            phi = 1.0 
            return phi
        end
    end
    # Eq from Bertram and Thronton (2009); avoid overflow
    phi = 1.0  / ( 1.0  + k2k3 * SafeDiv( H2O, Cl, 1.0e+30  )          )
    return phi 
end

function N2O5uptkByCloud( H )  
    #
    # Computes uptake of N2O5 on liquid water cloud.
    #
    # TYPE(HetState), INTENT(IN) :: H
    # REAL(dp)                   :: gamma, k
    const = 0.03  / 0.019 

    # Rxn probability is 0.03 at 298 K (JPL, Burkholder et al., 2015).
    # For temperature dependence, JPL recommends the same as sulfuric acid
    # aerosol at zero percent H2SO4, which is 0.019 at 298 K.
    # Then apply constant scale factor (0.03/0.019)
    gamma = const * exp( -25.5265  + 9283.76 /temp - 851801.0 /temp^2)
    #
    # Removal rate of N2O5 in liquid water cloud
    k = CloudHet( H, SR_MW(ind_N2O5), gamma, 0.02 , 1.0 , 1.0         )
    return k 
end

function N2O5uptkByStratHCl( H )  
    #
    # Sets heterogenous chemistry rate for N2O5(g) + HCl(l,s)
    # in polar stratospheric clouds and on tropospheric sulfate aerosol.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: gamma, k       # Rxn prob [1], Rxn rate [1/s]
    #
    k = 0.0 
    #
    if H%stratBox 
        #
        # Uptake on stratospheric liquid aerosol
        k = k + ( H%xArea(SLA) * H%KHETI_SLA(N2O5_plus_HCl) )
        #
        # Uptake on irregular ice cloud
        gamma = 0.03                          # Ice
        if H%natSurface 
            gamma = 0.003    # NAT
        end
        k = k + Ars_L1K( H%xArea(IIC), H%xRadi(IIC), gamma, SR_MW(ind_N2O5)  )
    end

    # Assume N2O5 is limiting, so update the removal rate accordingly
    k = kIIR1Ltd( C(ind_N2O5), C(ind_HCl), k )
    return k 
end

# #=========================================================================
# # Rate-law functions for NO2 and NO3
# #=========================================================================

function NO2uptk1stOrdAndCloud( H )  
    #
    # Computes the reaction rate [1/s] for 1st-order uptake of NO2.
    #
    # TYPE(HetState), INTENT(IN) :: H              # HetChem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    # REAL(dp)                   :: gamma, srMw    # local vars
    #
    k    = 0.0 
    srMw = SR_MW(ind_NO2)
    #
    # Uptake by mineral dust (aerosol types 1-7)
    gamma = 1.0e-8 
    k = k + Ars_L1k( H%xArea(DU1), H%xRadi(DU1), gamma, srMw )
    k = k + Ars_L1k( H%xArea(DU2), H%xRadi(DU2), gamma, srMw )
    k = k + Ars_L1k( H%xArea(DU3), H%xRadi(DU3), gamma, srMw )
    k = k + Ars_L1k( H%xArea(DU4), H%xRadi(DU4), gamma, srMw )
    k = k + Ars_L1k( H%xArea(DU5), H%xRadi(DU5), gamma, srMw )
    k = k + Ars_L1k( H%xArea(DU6), H%xRadi(DU6), gamma, srMw )
    k = k + Ars_L1k( H%xArea(DU7), H%xRadi(DU7), gamma, srMw )
    #
    # Uptake by tropospheric sulfate (aerosol type 8)
    gamma = 5e-6 
    k = k + Ars_L1k( H%xArea(SUL), H%xRadi(SUL), gamma, srMw )
    #
    # Uptake by black carbon (aerosol type 9)
    gamma = 1e-4 
    k = k + Ars_L1k( H%xArea(BKC), H%xRadi(BKC), gamma, srMw )
    #
    # Uptake by organic carbon (aerosol type 10)
    gamma = 1e-6 
    k = k + Ars_L1k( H%xArea(ORC), H%xRadi(ORC), gamma, srMw )
    #
    # Uptake by fine & coarse sea salt (aerosol types 11-12)
    if relhum < 40.0   
       gamma = 1.0e-8 
    elseif relhum > 70.0   
       gamma = 1.0e-4 
    else
       gamma = 1.0e-8  + (1e-4  - 1e-8 ) * (relhum - 40.0 )/30.0 
    end
    k = k + Ars_L1k( H%xArea(SSA), H%xRadi(SSA), gamma, srMw )
    k = k + Ars_L1k( H%xArea(SSC), H%xRadi(SSC), gamma, srMw )
    #
    # Uptake by stratospheric sulfate (aerosol type 13)
    # and by irregular ice cloud (aerosol type 14)
    gamma = 1.0e-4 
    k = k + H%xArea(SLA) * gamma
    k = k + Ars_L1k( H%xArea(IIC), H%xRadi(IIC), gamma, srMw )

    # Uptake of NO2 in cloud (liquid branch only)
    k = k + CloudHet( H, SR_MW(ind_NO2), 1.0e-8 , 0.0 , 1.0 , 0.0  )
    return k 
end

function Gam_NO3( aArea, aRadi, aWater, C_X, H ) 
    #
    # Calculates reactive uptake coef. for NO3 on salts and water
    #
    # REAL(dp),       INTENT(IN) :: aArea, aRadi, aWater, C_X
    # TYPE(Hetstate), INTENT(IN) :: H
    #
    # REAL(dp)            :: gamma
    # REAL(dp)            :: M_X, k_tot, H_X,    cavg
    # REAL(dp)            :: gb,  l_r,   WaterC, Vol, corr
    INV_AB = 1.0  / 1.3e-2 
    #
    Vol      = aArea * aRadi * 1.0e-3  / 3.0        # L/cm3 air
    WaterC   = aWater / 18.0e+12  / Vol               # mol/L aerosol
    #
    # Thermal velocity [cm/s]
    M_X      = MW(ind_NO3) * 1.0e-3                   # NO3 mol wt kg/mol
    cavg     = sqrt( EIGHT_RSTARG_T / ( H%Pi * M_X ) ) * 1.0e2 
    #
    k_tot    = ( 2.76e+6  * C_X ) + ( 23.0  * WaterC )
    #
    # Compute reactive uptake coefficient [1], but prevent div by zero
    gamma   = 0.0 
    if k_tot > 0.0   
       H_X   = 0.6  * CON_ATM_BAR            # M/bar
       l_r   = sqrt( 1.0e-5  / k_tot )       # diff const = 1e-5 for NO3
       #
       gb    = FOUR_R_T * H_X * l_r * k_tot / cavg
       corr  = Reactodiff_Corr( aRadi, l_r )
       gb    = gb * corr
       #
       gamma = 1.0  / ( INV_AB + 1.0  / gb )
    end
    return gamma
end

function NO3uptk1stOrdAndCloud( H )  
    # Computes reaction rate [1/s] for 1st-order uptake of NO3
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    # REAL(dp)                   :: gamma, srMw    # local vars
    #
    k    = 0.0 
    srMw = SR_MW(ind_NO3)
    #
    # Uptake by mineral dust bins 1-7
    gamma = 0.01 
    k = k + Ars_L1k( H%xArea(DU1), H%xRadi(DU1), gamma, srMw )
    k = k + Ars_L1k( H%xArea(DU2), H%xRadi(DU2), gamma, srMw )
    k = k + Ars_L1k( H%xArea(DU3), H%xRadi(DU3), gamma, srMw )
    k = k + Ars_L1k( H%xArea(DU4), H%xRadi(DU4), gamma, srMw )
    k = k + Ars_L1k( H%xArea(DU5), H%xRadi(DU5), gamma, srMw )
    k = k + Ars_L1k( H%xArea(DU6), H%xRadi(DU6), gamma, srMw )
    k = k + Ars_L1k( H%xArea(DU7), H%xRadi(DU7), gamma, srMw )
    #
    # Uptake by black carbon
    if relhum < 50.0 
       gamma = 2.0e-4 
    else
       gamma = 1.0e-3 
    end
    k = k + Ars_L1k( H%xArea(BKC), H%xRadi(BKC), gamma, srMw )
    #
    # Uptake by organic carbon
    gamma = 0.005 
    k = k + Ars_L1k( H%xArea(ORC), H%xRadi(ORC), gamma, srMw )
    #
    # Uptake by stratospheric sulfate liquid aerosol
    # and by irregular ice cloud
    gamma = 0.1 
    k = k + H%xArea(SLA) * gamma
    k = k + Ars_L1k( H%xArea(IIC), H%xRadi(IIC), gamma, srMw )

    # Uptake of NO3 in cloud (liquid and ice branches)
    k = k + CloudHet( H, SR_MW(ind_NO3), 0.002 , 0.001 , 1.0 , 1.0  )
    return k
end

function NO3hypsisClonSALA( H )  
    #
    # Computes the NO3(g) hypsis rate [1/s]  for Cl-
    # reacting on surface of fine sea-salt aerosol (SALA).
    #
    TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    REAL(dp)                   :: k              # rxn rate [1/s]
    REAL(dp) :: area, conc, gamma, radi, water   # local vars
    #
    # Compute reactive uptake coefficient [1]
    area  = H%aClArea
    radi  = H%aClRadi
    water = H%aWater(SS_FINE)
    conc  = H%Cl_conc_SSA
    gamma = Gam_NO3( area, radi, water, conc, H ) * 0.01 
    #
    # Reaction rate for surface of aerosol [1/s]
    area  = H%ClearFr * area
    k     = Ars_L1k( area, radi, gamma, SR_MW(ind_NO3) )
    return k 
end

function NO3hypsisClonSALC( H ) 
    #
    # Computes the NO3(g) hypsis rate [1/s]  for Cl-
    # reacting on surface of coarse sea-salt aerosol (SALC).
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    # REAL(dp) :: area, conc, gamma, radi, water   # local vars
    #
    # Compute reactive uptake coefficient [1]
    area  = H%xArea(SSC)
    radi  = H%xRadi(SSC)
    water = H%aWater(SS_COARSE)
    conc  = H%Cl_conc_SSC
    gamma = Gam_NO3( area, radi, water, conc, H ) * 0.01 
    #
    # Reaction rate for surface of aerosol [1/s]
    area  = H%ClearFr * area
    k     = Ars_L1k( area, radi, gamma, SR_MW(ind_NO3) )
    return k 
end

# =========================================================================
#  Hetchem rate-law functions for O3
# =========================================================================

function O3uptkByBrInTropCloud( H, Br_branch )  
    #
    # Computes the Sets the O3 uptake rate in tropospheric cloud
    # by Br- in either fine sea salt, coarse sea salt, or gas phase.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp),       INTENT(IN) :: Br_branch      # Branching ratio (A,C,G)
    # REAL(dp)                   :: k              # rxn rate [1/s]
    # REAL(dp)                   :: gamma          # local vars
    #
    k = 0.0 
    #
    # Exit if we are not in the troposphere
    if H%StratBox 
        return k 
    end
    #
    # Compute uptake of O3 by Br- in cloud
    gamma  = Gamma_O3_Br( H, H%rLiq, H%Br_conc_Cld )
    k      = CloudHet( H, SR_MW(ind_O3), gamma, 0.0 , Br_branch, 0.0  )
    return k 
end

function O3uptkByHBr( H )  
    #
    # Computes the O3 + HBr uptake rate in tropospheric cloud.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    # O3 + HBr uptake rate (gas-phase path), in trop cloud
    k = O3uptkByBrInTropCloud( H, H%frac_Br_CldG )
    #
    # Assume OH is limiting, so update the removal rate accordingly
    k = kIIR1Ltd( C(ind_O3), C(ind_HBr), k )
    return k 
end

function O3uptkByBrSALA( H )  
    #
    # Computes the uptake rate of O3 + Br- (in tropospheric
    # cloud) and on acidic fine sea salt (in clear sky).
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    # REAL(dp)                   :: area, gamma    # local vars
    #
    k = 0.0 
    #
    # O3 + Br- uptake by acidic fine sea salt, in trop cloud
    k = k + O3uptkByBrInTropCloud( H, H%frac_Br_CldA )
    #
    # O3 + Br- uptake on acidic fine sea-salt, clear sky
    if H%SSA_is_Acid 
       area  = H%ClearFr * H%aClArea
       gamma = Gamma_O3_Br( H, H%aClRadi, H%Br_conc_SSA )
       k     = k + Ars_L1K( area, H%aClRadi, gamma, SR_MW(ind_O3) )
    end

    # Assume OH is limiting, so update the removal rate accordingly
    k = kIIR1Ltd( C(ind_O3), C(ind_BrSALA), k )
    return k 
end

function O3uptkByBrSALC( H )  
    #
    # Computes the uptake rate of O3 + Br- in tropospheric
    # cloud and on acidic coarse sea salt.
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    # REAL(dp)                   :: area, gamma    # local vars
    #
    k = 0.0 
    #
    # O3 + Br- uptake by acidic coarse sea salt, in trop cloud
    k = k + O3uptkByBrInTropCloud( H, H%frac_Br_CldC )
    #
    # O3 + Br- uptake on acidic coarse sea salt, clear sky
    if H%SSC_is_Acid 
       area  = H%ClearFr * H%xArea(SSC)
       gamma = Gamma_O3_Br( H, H%xRadi(SSC), H%Br_conc_SSC )
       k     = k + Ars_L1K( area, H%xRadi(SSC), gamma, SR_MW(ind_O3) )
    end

    # Assume OH is limiting, so update the removal rate accordingly
    k = kIIR1Ltd( C(ind_O3), C(ind_BrSALC), k )
    return k 
end

function Gamma_O3_Br( H, Radius, C_Br )  
    #
    # Computes reactive uptake coefficient for Br- oxidation by O3.
    #
    TYPE(HetState), INTENT(IN) :: H             # Hetchem State
    REAL(dp),       INTENT(IN) :: radius        # Radius in cm
    REAL(dp),       INTENT(IN) :: C_Br          # Br- concentration
    REAL(dp)                   :: gamma         # rxn prob [1]
    #
    REAL(dp), PARAMETER  :: k0_O3 = 1.1e-2  * CON_ATM_BAR # Henry k0(O3)
    #
    REAL(dp) :: ab,  gb,     gd,  gs,       cavg,  H_X
    REAL(dp) :: M_X, KLangC, k_s, C_Br_surf, Nmax,  k_b,   D_l, l_r
    #
    gamma     = 0.0 
    #
    # If C_Br is zero, gamma is zero
    if C_Br <= 0.0  
        return gamma 
    end
    #
    # Henry's law for O3 (use constants for numerical stability)
    H_X       = k0_O3 * exp( 2300.0  * ( INV_temp - INV_T298 ) )
    #
    # Thermal velocity (cm/s)
    M_X       = MW(ind_O3) * 1.0e-3 
    cavg      = sqrt( EIGHT_RSTARG_T / ( H%Pi * M_X ) ) * 100.0 
    #
    Nmax      = 3.0e+14   # #/cm2
    KLangC    = 1.0e-13  #cm3
    k_s       = 1.0e-16  #cm2s-1, from ks*Nmax=0.03s-1
    #
    # [Br-(surf)] = 3.41E14 cm-2/M * [Br-(bulk)], but not gt Nmax.
    C_Br_surf = min( 3.41e+14  * C_Br, Nmax )
    gs        = ( 4.0  * k_s * C_Br_surf * KLangC * Nmax ) / ( cavg * ( 1.0  + KLangC * C(ind_O3)  ) )
    #
    k_b       = 6.3e+8  * exp(-4.45e+3  / temp )       # M-1 s-1
    D_l       = 8.9e-6                                   # cm2 s-1.
    #
    l_r    = sqrt( D_l / ( k_b * C_Br ) ) # cm
    gb     = FOUR_R_T * H_X * l_r * k_b * C_Br / cavg
    gb     = gb * ReactoDiff_Corr( Radius, l_r )
    #
    # Reactive uptake coefficient [1]
    gamma     = gb + gs
    return gamma
end

# =========================================================================
#  Hetchem rate-law functions for OH
# =========================================================================

function OHuptkBySALACl( H )  
    #
    # Computes uptake rate of OH + Cl on accumulation-mode sea-salt
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: gamma, k       # rxn prob [1], rxn rate [1/s]
    #
    # Compute uptake; gamma is from cf Knipping & Dabdub, 2002
    gamma = 0.04  * H%Cl_conc_SSA
    k = Ars_L1k( H%aClArea, H%aClRadi, gamma, SR_MW(ind_OH) )
    #
    # Assume OH is limiting, so update the removal rate accordingly
    k = kIIR1Ltd( C(ind_OH), C(ind_SALACL), k )
    return k
end

function OHuptkBySALCCl( H ) 
    #
    # Computes uptake rate of OH + Cl on coarse-mode sea-salt
    #
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: gamma, k       # rxn prob [1], rxn rate [1/s]
    #
    # Compute uptake; gamma is from cf Knipping & Dabdub, 2002
    gamma = 0.04  * H%Cl_conc_SSC
    k = Ars_L1k( H%xArea(SSC), H%xRadi(SSC), gamma, SR_MW(ind_OH) )
    #
    # Assume OH is limiting, so update the removal rate accordingly
    k = kIIR1Ltd( C(ind_OH), C(ind_SALCCL), k )
    return k
end

# =========================================================================
# Hetchem rate-law functions for VOC species
# =========================================================================

function GLYXuptk1stOrd( srMw, H ) 
    # Computes the reaction rate [1/s] for 1st-order uptake of GLYX.
    # Only consider inorganic aqueous aerosols with RH > 35%
    # and use diffe rent uptake for day & night.
    #
    # REAL(dp),       INTENT(IN) :: srMw   # sqrt( mol wt )
    # TYPE(HetState), INTENT(IN) :: H      # Hetchem State
    # REAL(dp)                   :: k      # rxn rate [1/s]
    # REAL(dp)                   :: gamma  # local vars
    #
    k     = 0.0 
    gamma = 0.0 
    #
    # Uptake by tropospheric sulfate
    if RELHUM >= CRITRH 
        if SUNCOS > 0.0   
            gamma = 4.4e-3    # cf Liggio et al 2005
        else
            gamma = 8.0e-6    # F. McNeill, to E. Marais (2015)
        end
        k = k + Ars_L1k( H%xArea(SUL), H%xRadi(SUL), gamma, srMw )
    end
    return k
end

function  EpoxUptkGamma( srMw, H )  
    #
    # Gomputes the GAMMA uptake probability for EPOXUPTK hydrolysis to
    # form 2-methyltetrols (AITET). (eam, 2014).
    #
    # Calculation is only done for inorganic aqueous phase aerosols.
    # This calculation uses the parameterization of Gaston et al., EST, 2014.
    # Redistribution of products (e.g. AITET) to yield organosulfates and
    # organonitrates is done in SOA_CHEMISTRY in carbon_mod.F.
    # This is only done for IEPOX and HMML if it's an SOA simulation
    #
    # REAL(dp),       INTENT(IN) :: srMw           # sqrt(mol wt)
    # TYPE(HetState), INTENT(IN) :: H              # HetChem State
    # REAL(dp)                   :: gamma          # Uptake prob [1]
    # REAL(dp) :: aervol, kpart, xmms, val1, val2, val3, valtmp # local vars
    #
    # Gas-phase diffusion constant [cm2/s]:
    DIFF_N2O5_STD = 1.0e-1 
    #
    # Mass accommodation coefficient [unitless]:
    MACOEFF = 1.0e-1 
    K_HPLUS = 3.6e-2 
    K_NUC   = 2.0e-4 
    K_HSO4  = 7.3e-4 
    K_HYDRO = 0.0e+0 
    #
    # Effective Henry's Law constant of IEPOX for reactive uptake to aqueous
    # aerosols (M/atm).  Eloise Marais (2015/07) reset this to the value from
    # [Nguyen et al., 2014] in order to accomodate reduction in yields of RIP
    # (which is the precursor of IEPOX).
    HSTAR_EPOX  = 1.7e+7 
    #
    # Initialize
    gamma  = 0.0 
    valTmp = 0.0 
    #
    # Calculate aerosol volume (use formula in aerosol_mod.F):
    aerVol = ( H%xArea(SUL) *  H%xRadi(SUL) ) / 3.0 
    #
    # Calculate mean molecular speed [cm/s]:
    xmms = sqrt( ( 2.117e+8  * temp ) / ( srMw * srMw ) )
    #
    # Calculate first-order particle-phase reaction rate:
    # (assume [H+] = proton activity)
    # KHYDRO is only important for alkylnitrates (not currently used).
    kPart = ( K_HPLUS * H%H_PLUS)+( K_NUC* H%H_PLUS*( H%NO3_molal + H%SO4_molal )) + ( K_HSO4  * H%HSO4_molal)+(K_HYDRO)
    #
    # Calculate the first uptake parameterization term:
    val1 = ( H%xRadi(SUL) * xmms ) / ( 4.0  * DIFF_N2O5_STD )
    #
    # Calculate the second uptake parameterization term:
    val2 = ( 1.0  / MACOEFF )
    #
    # Calculate the third uptake parameterization term:
    if  H%xArea(SUL) > 0.0 && XMMS > 0.0 
        valTmp = ( FOUR_RGASLATM_T * aerVol * HSTAR_EPOX * kPart ) / ( H%xArea(SUL) * xmms)
    end
    #
    val3 = 0.0 
    if valTmp > 0.0
        val3 = 1.0  / valtmp
    end
    #
    # Account for small reaction rates:
    gamma = 0.0 
    if kPart >= 1.e-8
        gamma = 1.0  / ( val1 + val2 + val3 )
    end
    if gamma <  0.0   
        gamma = 0.0 
    end
    return gamma 
end

function IEPOXuptk1stOrd( srMw, doScale, H ) 
    #
    # Sets the heterogenous chemistry rate for first-order
    # uptake of ICHE, IEPOXA, IEPOXB, and IEPOXD.
    #
    # REAL(dp),       INTENT(IN) :: srMw           # sqrt( mol wt )
    # LOGICAL,        INTENT(IN) :: doScale        # =T for HMML, else F
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    # REAL(dp)                   :: gamma          # local vars
    #
    k     = 0.0 
    gamma = 0.0 
    #
    # Only consider inorganic aqueous aerosols with RH > 35%.
    if RELHUM >= CRITRH  
       #
       # Get GAMMA for IEPOX hydrolysis
       gamma = EpoxUptkGamma( srMw, H )
       #
       # Scale down gamma if [H+] > 8d-5 (cf Riedel et al, 2015)
       if doScale && ( H%H_PLUS > 8.0e-5  ) )  
          gamma = gamma / 30.0 
       end
       #
       # Uptake by tropospheric sulfate
       k = k + Ars_L1k( H%xArea(SUL), H%xRadi(SUL), gamma, srMw )
    end
    return k 
end

function MGLYuptk1stOrd( srMw, H ) 
    #
    # Computes the reaction rate [1/s] for 1st order uptake of MGLY and PYAC.
    #
    # REAL(dp),       INTENT(IN) :: srMw           # sqrt( mol wt )
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: gamma, k       # rxn prob [1], rxn rate [1/s]
    #
    k     = 0.0 
    gamma = 0.0 
    #
    # Only consider inorganic aqueous aerosols with RH > 35%.
    if RELHUM >= CRITRH 
       #
       # Define gamma for MGLY: Obtained by scaling gamma GLYX by the
       # ratio of effective Henry's law constants for GLYX (3d7) and
       # MGLY (3.7d3) (eam, 02/2015):
       gamma = 3.6e-7 
       #
       # Uptake by tropospheric sulfate
       k = k + Ars_L1k( H%xArea(SUL), H%xRadi(SUL), gamma, srMw )
    end
    return k 
end

function VOCuptk1stOrd( srMw, gamma, H ) 
    #
    # Computes the rxn rate [1/s]for 1st-order uptake of several VOC species.
    #
    # REAL(dp),       INTENT(IN) :: srMw, gamma    # sqrt( mol wt ), rxn prob
    # TYPE(HetState), INTENT(IN) :: H              # Hetchem State
    # REAL(dp)                   :: k              # rxn rate [1/s]
    #
    # Initialize
    k  = 0.0_dp
    #
    # Only consider inorganic aqueous aerosols with RH > 35%.
    if RELHUM >= CRITRH  
       k = k + Ars_L1k( H%xArea(SUL), H%xRadi(SUL), gamma, srMw )
       k = k + Ars_L1k( H%xArea(BKC), H%xRadi(BKC), gamma, srMw )
       k = k + Ars_L1k( H%xArea(ORC), H%xRadi(ORC), gamma, srMw )
       k = k + Ars_L1k( H%xArea(SSA), H%xRadi(SSA), gamma, srMw )
       k = k + Ars_L1k( H%xArea(SSC), H%xRadi(SSC), gamma, srMw )
       k = k + H%xArea(SLA) * gamma
       k = k + Ars_L1k( H%xArea(IIC), H%xRadi(IIC), gamma, srMw )
    end
    return k 
end

"""
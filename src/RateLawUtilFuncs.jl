"""
The original model is from:
https://github.com/geoschem


define based on Set_Kpp_GridBox_Values
https://github.com/geoschem/geos-chem/blob/0f85431c4beedb0af8325b54d789f910cd3d6357/GeosCore/fullchem_mod.F90

-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
                 geos - chem global chemical transport model                 
-  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 

#IROUTINE: rateLawUtilFuncs

#DESCRIPTION: Provides common functions for computing reaction rates.
"""

temp = 288.15
INV_TEMP = 1.0 / temp
temp_over_K300 = temp / 300.0
k300_over_temp = 300.0 / temp
SR_TEMP = sqrt(temp)
NUMDEN = 1
NSPEC = 1
NREACT = 1


##########################################################################
######                   ARRHENIUS FUNCTIONS                         #####
##########################################################################

function GCARR_ab(a0, b0)
  # Arrhenius function, skipping computation of EXP( c0 / T ),
  # which evaluates to 1 when c0 = 0.0  This avoids excess CPU
  # cycles. (bmy, 12 / 18 / 20)

  k = a0 * k300_over_temp^b0
  return k
end

function GCARR_ac(a0, c0)
  # Arrhenius function, skipping computation of ( 300 / T )^b0,
  # which evaluates to 1 when b0 = 0.0  This avoids excess CPU
  # cycles (bmy, 12 / 18 / 20)
  #


  k = a0 * exp(c0 / temp)
  return k
end

function GCARR_abc(a0, b0, c0)
  # Arrhenius function, using all 3 terms.
  # Use this when a0, b0, c0 are all nonzero.
  #

  k = a0 * exp(c0 / temp) * k300_over_temp^b0
  return k
end

##########################################################################
######         COMMON FUNCTIONS FOR COMPUTING UPTAKE RATES           #####
##########################################################################

function ars_l1k(area, radius, Γ, srmw)
  #
  # Calculates the 1st - order loss rate of species on wet aerosol surface.
  #

  # If Γ or radius is very small, set rate to zero and return
  if (Γ < 1.0e-30 || radius < 1.0e-30)
    k = 0.0
    return k
  end
  #
  # DFKG = Gas phase diffusion coeff [cm2 / s] (order of 0.1)
  dfkg = (9.45e + 17 / numden) * sr_temp * sqrt(3.472e - 2 + 1.0 / (srmw * srmw))
  #
  # Compute ArsL1k according to the formula listed above
  k = area / ((radius / dfkg) + 2.749064e - 4 * srmw / (Γ * sr_temp))
  return k
end

function kiir1ltd(concgas, conceduct, kisource)
  #
  # Determine removal rates for both species in an uptake reaction.
  # - Assume that the 1st reactant (concGas) is limiting.
  # - Assume that the 2nd reactant (concEduct) is "abundant".
  # - Calculate the overall rate (kII) based only on the uptake
  #   rate of the first reactant.
  # NOTE: Rewritten for computational efficiency (bmy, 5 / 13 / 21)
  #

  kigas = 0.0
  kieduct = 0.0
  kii = 0.0
  #
  # Prevent div by zero.  Now use 1.0 as the error trap for concEduct.
  # 100 and 1e-8 (the previous criteria) were too large and too small,
  # respectively.  See https: / / github.com / geoschem / geos - chem / issues / 1115.0
  # -  - Seb Eastham, Bob Yantosca (09 Feb 2022)
  if (conceduct < 1.0)
    return kii
    if (!is_safediv(concgas * kisource, conceduct))
      return kii
      #
      # Compute rates
      kigas = kisource
      kieduct = kigas * concgas / conceduct
      kii = kigas / conceduct
      #
      # Enforce a minimum lifetime?
      if (kigas > 0.0)
        #
        # Calculate lifetime of each reactant against removal
        lifea = safediv(1.0, kigas, 0.0)
      end
      lifeb = safediv(1.0, kieduct, 0.0)
    end
    #
    # Check if either lifetime is "too short"
    if ((lifea < lifeb) && (lifea < het_min_life))
      kii = safediv(het_min_rate, conceduct, 0.0)
    elseif (lifeb < het_min_life)
      kii = safediv(het_min_rate, concgas, 0.0)
    end
  end
end

##########################################################################
######         COMMON FUNCTIONS FOR COMPUTING UPTAKE RATES           #####
##########################################################################

function coth(x)
  result(f_x)
  #
  # Hyperbolic cotangent = [1 + exp( - 2x)] / [1 - exp( - 2x)]

  y = exp(-2.0 * x)
  f_x = (1.0 + y) / (1.0 - y)
end

function reactodiff_corr(radius, l)
  result(corr)
  #
  # For x = radius / l, correction =  COTH( x ) - ( 1 / x )
  # Correction approaches 1 as x becomes large, corr(x>1000)~1
  # Correction approaches x / 3 as x goes towards 0

  x = radius / l
  if (x > 1000.0)
    corr = 1.0
    return
  end
  if (x < 0.1)
    corr = x / 3.0
    return
  end
  corr = coth(x) - (1.0 / x)
end

##########################################################################
######   COMMON FUNCTIONS FOR ENFORCING SAFE NUMERICAL OPERATIONS    #####
##########################################################################

function safediv(num, denom, alt)
  result(quot)
  #
  # Performs "safe division", that is to prevent overflow, underlow,
  # NaN, or infinity errors.  An alternate value is returned if the
  # division cannot be performed.

  # Exponent difference (base 2)
  # For REAL * 8, max exponent = 1024 and min = - 1021
  ediff = exponent(num) - exponent(denom)
  #  end
  #
  if (ediff > 1023 || denom == 0.0)

    quot = alt

  elseif (ediff < -1020)
    quot = 0.0
  else
    quot = num / denom
  end
end

function is_safediv(num, denom)
  #
  # Returns TRUE if a division can be performed safely.

  # Exponent difference (base 2)
  # For REAL * 8, max exponent = 1024 and min = - 1021
  safe = true
  ediff = exponent(num) - exponent(denom)

  #
  if (ediff < -1020 || ediff > 1023 || denom == 0.0)

    safe = false

  end
  return safe

  function issafeexp(x)
    result(safe)
    #
    # Returns TRUE if an exponential can be performed safely

    # Note EXP( 708 ) = 8.2e + 307 and EXP( - 708 ) = 3.3e-308, which are
    # very close to the maximum representable values at double precision.
    safe = (abs(x) < 709.0)
  end

  function safeexp(x, alt)
    result(y)
    #
    # Performs a "safe exponential", that is to prevent overflow, underflow,
    # underlow, NaN, or infinity errors when taking the value EXP( x ).  An
    # alternate value is returned if the exponential cannot be performed.
    #
    # real(dp), intent(in) :: x, alt
    # real(dp)             :: y
    #
    y = alt
    if (abs(x) < 709.0)
      y = exp(x)
    end
  end
  #EOC
end
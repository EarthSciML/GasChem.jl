"""
The original model is from:
https://github.com/geoschem

"""

include("RateLawUtilFuncs.jl")

GCARR_ab(a0, b0) = a0 * K300_OVER_temp^b0

GCARR_ab(a0, c0) = a0 * exp(c0 / TEMP)

GCARR_ab(a0, b0, c0) = a0 * exp(c0 / TEMP) * K300_OVER_temp^b0


function Ars_L1k(area, radius, gamma, srMw)
    # !
    # ! Calculates the 1st-order loss rate of species on wet aerosol surface.
    # !
    # ! If gamma or radius is very small, set rate to zero and return
    if (gamma < 1.0e-30 || radius < 1.0e-30)
        k = 0.0
    else
        dfkg = (9.45e+17 / NUMDEN) * SR_TEMP * SQRT(3.472E-2 + 1.0 / (srMw * srMw))
        k = area / ((radius / dfkg) + 2.749064E-4 * srMw / (gamma * SR_TEMP))
        return K
    end
end

function kIIR1Ltd(concGas, concEduct, kISource)
    # !
    # ! Determine removal rates for both species in an uptake reaction.
    # ! - Assume that the 1st reactant (concGas) is limiting.
    # ! - Assume that the 2nd reactant (concEduct) is "abundant".
    # ! - Calculate the overall rate (kII) based only on the uptake
    # !   rate of the first reactant.
    # ! NOTE: Rewritten for computational efficiency (bmy, 5/13/21)
    # !
    KiGas = 0
    KIEduct = 0
    KII = 0

    # !
    # ! Prevent div by zero.  Now use 1.0 as the error trap for concEduct.
    # ! 100 and 1e-8 (the previous criteria) were too large and too small,
    # ! respectively.  See https://github.com/geoschem/geos-chem/issues/1115.
    # !  -- Seb Eastham, Bob Yantosca (09 Feb 2022)

    if (concEduct < 1.0)
        return KII
    end
    if Is.SafeDiv(concGas * KISource, concEduct)
        continuetemp = 1
    else
        return KII
    end

    KiGas = kISource
    KIEduct = KIGas * concGas / concEduct
    KII = KIGas / concEduct
    # ! Enforce a minimum lifetime?
    if KIGas > 0.0
        lifeA = SafeDiv(1.0, kIGas, 0.0)
        lifeB = SafeDiv(1.0, kIEduct, 0.0)
        if (lifeA < lifeB) && lifeA > HET_MIN_LIFE
            kII = SafeDiv(HET_MIN_RATE, concEduct, 0.0)
        elseif lifeB < HET_MIN_LIFE
            kII = SafeDiv(HET_MIN_RATE, concGas, 0.0)
        end

    end
    return KII
end

function CloudHet(H, srMw, gamLiq, gamIce, brLiq, brIce)

    # ! Function CloudHet calculates the loss frequency (1/s) of gas species
    # ! due to heterogeneous chemistry on clouds in a partially cloudy grid
    # ! cell. The function uses the "entrainment limited uptake" equations of
    # ! Holmes et al. (2019). Both liquid and ice water clouds are treated.
    # !
    # ! For gasses that are that are consumed in multiple aqueous reactions
    # ! with different products, CloudHet can provide the loss frequency for
    # ! each reaction branch using the branch ratios (branchLiq, branchIce).
    # !
    # ! Reference:
    # ! Holmes, C.D., Bertram, T. H., Confer, K. L., Ronan, A. C., Wirks,
    # !   C. K., Graham, K. A., Shah, V. (2019) The role of clouds in the
    # !   tropospheric NOx cycle: a new modeling approach for cloud chemistry
    # !   and its global implications, Geophys. Res. Lett. 46, 4980-4990,
    # !   https://doi.org/10.1029/2019GL081990
    # TYPE(HetState), INTENT(IN) :: H              ! Hetchem State object
    # REAL(dp),       INTENT(IN) :: srMw           ! SQRT( mol wt [g/mole] )
    # REAL(dp),       INTENT(IN) :: gamLiq         ! Rxn prob, liquid [1]
    # REAL(dp),       INTENT(IN) :: gamIce         ! Rxn prob, ice [1]
    # REAL(dp),       INTENT(IN) :: brLiq          ! Frac of reactant consumed
    # REAL(dp),       INTENT(IN) :: brIce          !  in liq & ice branches [0-1]
    # REAL(dp)                   :: kHet           ! Grid-avg loss frequency [1/s]

    # ! If cloud fraction < 0.0001 (0.01%) or there is zero cloud surface
    # ! area, then return zero uptake
    if (H % CldFr < 0.0001) || (H % aLiq + H % aIce <= 0.0)
        kHet = 0.0
        return kHet
    end

    kI = 0.0
    kIb = 0.0
    ktmp = 0.0
    kHet = 0.0

    # !-----------------------------------------------------------------------
    # ! Liquid branch (skip if the liquid branching ratio is zero)
    # !-----------------------------------------------------------------------
    if brLiq > 0.0
        #    ! Convert grid-average cloud condensate surface area density
        #    ! to in-cloud surface area density
        area = SafeDiv(H % aLiq, H % CldFr, 0.0)

        #    ! Skip if no area
        if area > 0.0

            #   ! In-cloud loss frequency [1/s]
            ktmp = Ars_L1K(area, H % rLiq, gamLiq, srMw)
            kI = kI + ktmp

            #   ! In-cloud loss frequency for liquid rxn branch [1/s]
            kIb = kIb + (ktmp * brLiq)
        end
    end

    # !------------------------------------------------------------------
    # ! Ice branch (skip if the ice branching ratio is zero)
    # !------------------------------------------------------------------
    if brIce > 0.0

        #    ! Convert grid-average cloud condensate surface area density
        #    ! to in-cloud surface area density
        area = SafeDiv(H % aIce, H % CldFr, 0.0)

        #    ! Skip if no area
        if area > 0.0

            #   ! In-cloud loss frequency [1/s]
            ktmp = Ars_L1K(area, H % rIce, gamIce, srMw)
            kI = kI + ktmp

            #   ! In-continue loud loss frequency for ice rxn branch [1/s]
            kIb = kIb + (ktmp * brIce)
        end
    end

    # !------------------------------------------------------------------
    # ! Mean branch ratio for reaction of interest in cloud
    # ! (averaged over ice and liquid)
    # !
    # ! If the division can't be done, set kHet = 0 and return
    # !------------------------------------------------------------------
    # branch = SafeDiv( kiB, kI, 0.0 )  #mark KIB  typo
    branch = SafeDiv(kIb, kI, 0.0)  #mark KIB  typo
    if !branch > 0.0
        kHet = 0.0
        return kHet
    end


    # !------------------------------------------------------------------------
    # ! Grid-average loss frequency
    # !
    # ! EXACT expression for entrainment-limited uptake
    # !------------------------------------------------------------------------

    # ! Ratio (in cloud) of heterogeneous loss to detrainment, s/s
    kk = kI * tauc

    # ! Ratio of volume inside to outside cloud
    # ! ff has a range [0,+inf], so cap it at 1e30
    ff = SafeDiv(H % CldFr, H % ClearFr, 1.0e+30)
    ff = min(ff, 1.0e+30)

    # ! Ratio of mass inside to outside cloud
    # ! xx has range [0,+inf], but ff is capped at 1e30, so shouldn't overflow.
    xx = (ff - kk - 1.0) / 2.0 + sqrt(1.0 + ff * ff + kk * kk + 2.0 * ff + 2.0 * kk - 2.0 * ff * kk) / 2.0

    # ! Do not let xx go negative, as this can cause numerical instability.
    # ! See https://github.com/geoschem/geos-chem/issues/1205
    xx = max(xx, 0.0)

    # ! Overall heterogeneous loss rate, grid average, 1/s
    # ! kHet = kI * xx / ( 1d0 + xx )
    # !  Since the expression ( xx / (1+xx) ) may behave badly when xx>>1,
    # !  use the equivalent 1 / (1 + 1/x) with an upper bound on 1/x
    kHet = kI / (1.0 + SafeDiv(1.0, xx, 1.0e+30))

    # ! Overall loss rate in a particular reaction branch, 1/s
    kHet = kHet * branch
    return kHet
end

function Cld_Params(AD, CLDF, FRLAND, FROCEAN, QI, QL, T, H)

    # !
    # ! Returns ice and liquid cloud parameters (based on State_Met)
    # ! for cloud particles.
    # !
    # ! References:
    # !  Heymsfield, A. J., Winker, D., Avery, M., et al. (2014). Relationships
    # !   between ice water content and volume extinction coefficient from in
    # !   situ observations for temperatures from 0° to –86°C: implications
    # !   for spaceborne lidar retrievals. Journal of Applied Meteorology and
    # !   Climatology, 53(2), 479–505. https://doi.org/10.1175/JAMC-D-13-087.1
    # !
    # !  Schmitt, C. G., & Heymsfield, A. J. (2005). Total Surface Area Estimates
    # !   for Individual Ice Particles and Particle Populations. Journal of
    # !   Applied Meteorology, 44(4), 467–474. https://doi.org/10.1175/JAM2209.1
    # !
    # REAL(dp),       INTENT(IN)    :: AD          ! Air mass [kg]
    # REAL(dp),       INTENT(IN)    :: CLDF        ! Cloud fraction [1]
    # REAL(dp),       INTENT(IN)    :: FRLAND      ! Land fraction [1]
    # REAL(dp),       INTENT(IN)    :: FROCEAN     ! Ocean fraction [1]
    # REAL(dp),       INTENT(IN)    :: QI          ! Ice mixing ratio [kg/kg]
    # REAL(dp),       INTENT(IN)    :: QL          ! Liquid mixing ratio [kg/kg]
    # REAL(dp),       INTENT(IN)    :: T           ! Temperature [K]

    # output OUTPUT PARAMETERS:


    # !
    # TYPE(HetState), INTENT(INOUT) :: H          # ! Hetchem State object
    #https://docs.julialang.org/en/v1/manual/types/#Mutable-Composite-Types

    # !
    # ! !REMARKS:
    # !EOP
    # !------------------------------------------------------------------------------
    # !BOC
    # !
    # ! !DEFINED PARAMETERS:
    # !
    # ! Cloud droplet radius in continental warm clouds [cm]
    CLDR_CONT = 6.0e-4

    # ! Cloud droplet radius in marine warm clouds [cm]
    CLDR_MARI = 10.0e-4

    # ! Ice cloud droplet radius [cm]
    CLDR_ICE = 38.5e-4

    # ! Density of H2O liquid [kg/cm3]
    DENS_LIQ = 0.001

    # ! Density of H2O ice [kg/cm3]
    DENS_ICE = 0.91e-3
    # !
    # ! !LOCAL VARIABLES:
    # !
    # REAL(dp) :: alpha, beta

    # !=======================================================================
    # ! CLD_PARAMS begins here!
    # !=======================================================================

    # ! Exit if there is no cloud
    if (QL + QI <= 0.0) || (CLDF <= 0.0)
        H % rLiq = CLDR_CONT
        H % rIce = CLDR_ICE
        H % ALiq = 0.0
        H % VLiq = 0.0
        H % AIce = 0.0
        H % VIce = 0.0
        return
    end

    # !-----------------------------------------------------------------------
    # ! In GC 12.0 and earlier, the liquid water volume was set to zero at
    # ! temperatures colder than 258K and over land ice (Antarctica &
    # ! Greenland). That was likely legacy code from GEOS-4, which provided
    # ! no information on cloud phase. As of GC 12.0, all met data sources
    # ! provide explicit liquid and ice condensate amounts, so we use those
    # ! as provided. (C.D. Holmes)
    # !
    # ! Liquid water clouds
    # !
    # ! Droplets are spheres, so
    # ! Surface area = 3 * Volume / Radius
    # !
    # ! Surface area density = Surface area / Grid volume
    # !-----------------------------------------------------------------------
    if (FRLAND > FROCEAN)
        THEN
        H % rLiq = CLDR_CONT    #  ! Continental cloud droplet radius [cm]
    else
        H % rLiq = CLDR_MARI    #  ! Marine cloud droplet radius [cm]

    end

    # ! get the volume of cloud condensate [cm3(condensate)/cm3(air)]
    # ! QL is [g/g]
    H % VLiq = QL * AD / DENS_LIQ / H % vAir
    H % VIce = QI * AD / DENS_ICE / H % vAir
    H % ALiq = 3.0_dp * H % vLiq / H % rLiq #??

    # !-----------------------------------------------------------------------
    # ! Ice water clouds
    # !
    # ! Surface area calculation requires information about ice crystal size
    # ! and shape, which is a function of temperature. Use Heymsfield (2014)
    # ! empirical relationships between temperature, effective radius,
    # ! surface area and ice water content.
    # !
    # ! Schmitt and Heymsfield (2005) found that ice surface area is about
    # ! 9 times its cross-sectional area.
    # !
    # ! For any shape,
    # !   Cross section area = pi * (Effective Radius)^2, so
    # !   Cross section area = 3 * Volume / ( 4 * Effective Radius ).
    # !
    # ! Thus, for ice
    # !   Surface area = 9 * Cross section area
    # !                = 2.25 * 3 * Volume / Effective Radius
    # ! (C.D. Holmes)
    # !-----------------------------------------------------------------------

    # ! Heymsfield (2014) ice size parameters
    if T < 202.0        # ! -71 C
        alpha = 83.3
        beta = 0.0184
    elseif T < 217.0    # ! -56 C
        alpha = 9.1744e+4_dp
        beta = 0.117_dp
    else
        alpha = 308.4_dp
        beta = 0.0152_dp
    end

    # ! Effective radius, cm
    H % rIce = 0.5 * alpha * exp(beta * (T - 273.15)) / 1e+4

    # ! Ice surface area density, cm2/cm3
    H % aIce = 3.0 * H % vIce / H % rIce * 2.25

end

#########################################################################
#####         COMMON FUNCTIONS FOR COMPUTING UPTAKE RATES           #####
#########################################################################

function coth(x)

    # Hyperbolic cotangent = [1 + exp(-2x)] / [1 - exp(-2x)]


    y = exp(-2.0 * x)
    f_x = (1.0 + y) / (1.0 - y)

    return f_x
end

function ReactoDiff_Corr(radius, l)
    # !
    # ! For x = radius / l, correction =  COTH( x ) - ( 1/x )
    # ! Correction approaches 1 as x becomes large, corr(x>1000)~1
    # ! Correction approaches x/3 as x goes towards 0
    # !
    # REAL(dp), INTENT(IN)  :: l, radius           ! [cm] and [cm]
    # REAL(dp)              :: x, corr
    # !
    x = radius / l
    if x > 1000.0
        corr = 1.0
        return corr
    end

    if x < 0.1
        corr = x / 3.0
        return corr
    end
    corr = coth(x) - (1.0 / x)
    return corr
end


function SafeDiv(num, denom, alt)
    # !
    # ! Performs "safe division", that is to prevent overflow, underlow,
    # ! NaN, or infinity errors.  An alternate value is returned if the
    # ! division cannot be performed.
    # REAL(dp), INTENT(IN) :: num, denom, alt
    # REAL(dp)             :: ediff, quot
    # !
    # ! Exponent difference (base 2)
    # ! For REAL*8, max exponent = 1024 and min = -1021
    #ediff = EXPONENT( num ) - EXPONENT( denom )

    ediff = exponent(num) - exponent(denom)

    if ediff > 1023 || denom == 0.0
        quot = alt
    elseif (ediff < -1020)
        THEN
        quot = 0.0
    else
        quot = num / denom
    end
    return quot
end

function Is_SafeDiv(num, denom)
    # !
    # ! Returns TRUE if a division can be performed safely.
    # REAL(dp), INTENT(IN) :: num, denom
    # LOGICAL              :: safe
    # REAL(dp)             :: ediff
    # !
    # ! Exponent difference (base 2)
    # ! For REAL*8, max exponent = 1024 and min = -1021
    safe = true
    ediff = exponent(num) - exponent(denom)

    if ediff < -1020 || ediff > 1023 || denom == 0.0
        safe = false
    end
    return safe
end

function IsSafeExp(x)
    # !
    # ! Returns TRUE if an exponential can be performed safely
    # !
    # REAL(dp), INTENT(IN) :: x
    # LOGICAL              :: safe
    # !
    # ! Note EXP( 708 ) = 8.2e+307 and EXP( -708 ) = 3.3e-308, which are
    # ! very close to the maximum representable values at double precision.
    safe = (abs(x) < 709.0)
    return safe
end


function SafeExp(x, alt)
    # !
    # ! Performs a "safe exponential", that is to prevent overflow, underflow,
    # ! underlow, NaN, or infinity errors when taking the value EXP( x ).  An
    # ! alternate value is returned if the exponential cannot be performed.
    # !

    y = alt
    if abs(X) < 709.0
        y = exp(x)
    end
    return y
end

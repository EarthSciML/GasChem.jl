export FastJX

# Effective wavelength in 18 bins covering 177–850 nm
const WL = SA_F32[
    187,
    191,
    193,
    196,
    202,
    208,
    211,
    214,
    261,
    267,
    277,
    295,
    303,
    310,
    316,
    333,
    380,
    574
]

# Top of the atmosphere solar flux in 18 bins
const top_flux = SA_F32[
    1.391E+12,
    1.627E+12,
    1.664E+12,
    9.278E+11,
    7.842E+12,
    4.680E+12,
    9.918E+12,
    1.219E+13,
    6.364E+14,
    4.049E+14,
    3.150E+14,
    5.889E+14,
    7.678E+14,
    5.045E+14,
    8.902E+14,
    3.853E+15,
    1.547E+16,
    2.131E+17
]

#   Cross sections and quantum yield from GEOS-CHEM "FJX_spec.dat" for photo-reactive species included in SuperFast:




# HOCl=>OH+Cl      JPL10
const ϕ_HOCl_jx = 1.0f0
const σ_HOCl = SA_F32[0, 3.695e-21, 1.571e-20, 2.435e-20, 5.887e-20, 5.424e-20, 5.798e-20, 6.694e-20, 1.019e-19, 6.541e-20, 5.43e-20, 5.57e-20, 6.067e-20, 5.955e-20, 5.376e-20, 3.12e-20, 2.197e-21, 0]
const σ_HOCl_interp = [(T) -> σ_HOCl[i] for i in 1:18]

# H2CO=>H2+CO      JPL10
const ϕ_H2COb_jx = 1.0f0
const σ_H2COb_interp = create_fjx_interp(
    [223.0f0, 298.0f0],
    [
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 3.642e-21, 5.787e-21, 5.316e-21, 8.181e-21, 7.917e-21, 4.011e-21, 1.081e-20, 1.082e-20, 6.842e-23, 0],
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 3.649e-21, 5.768e-21, 5.305e-21, 8.154e-21, 7.914e-21, 4.002e-21, 1.085e-20, 1.085e-20, 6.819e-23, 0],
    ]
)
const ϕ_CH2Ob_jx = ϕ_H2COb_jx
const σ_CH2Ob_interp = σ_H2COb_interp

# CH2C(CH3)CHO >   Methacrolein => CH2=C(CH3)+HCO (q=0.003) JPL10
const ϕ_MeAcr_jx = 1.0f0
const σ_MeAcr = SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 9.45e-24, 2.852e-23, 5.252e-23, 7.912e-23, 1.182e-22, 1.5e-22, 1.779e-22, 1.959e-22, 1.305e-23, 0]
const σ_MeAcr_interp = [(T) -> σ_MeAcr[i] for i in 1:18]

# N2O5=>NO2+NO3    JPL10
const ϕ_N2O5_jx = 1.0f0
const σ_N2O5_interp = create_fjx_interp(
    [233.0f0, 300.0f0],
    [
        SA_F32[0, 0, 8.922e-19, 1.183e-18, 5.868e-18, 4.682e-18, 3.395e-18, 2.613e-18, 2.138e-19, 2.155e-19, 1.988e-19, 3.772e-20, 2.182e-20, 1.334e-20, 8.419e-21, 2.621e-21, 4.355e-23, 0],
        SA_F32[0, 0, 1.078e-18, 1.429e-18, 7.088e-18, 5.655e-18, 4.101e-18, 3.156e-18, 2.606e-19, 2.645e-19, 2.454e-19, 5.154e-20, 3.21e-20, 2.135e-20, 1.468e-20, 5.902e-21, 2.025e-22, 0],
    ]
)


# CF3Br=>          CF3Br = halon 1301         JPL10
const ϕ_H1301_jx = 1.0f0
const σ_H1301_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[6.268e-20, 8.611e-20, 9.979e-20, 1.158e-19, 1.4e-19, 1.404e-19, 1.282e-19, 1.096e-19, 6.655e-22, 5.389e-21, 6.89e-21, 1.311e-25, 1.58e-26, 0, 0, 0, 0, 0],
        SA_F32[5.426e-20, 7.273e-20, 8.354e-20, 9.598e-20, 1.173e-19, 1.227e-19, 1.165e-19, 1.047e-19, 1.413e-21, 7.024e-21, 7.262e-21, 4.028e-25, 4.856e-26, 0, 0, 0, 0, 0],
    ]
)


# CFCl3=>CFCl2+Cl  CFC-11    JPL10
const ϕ_CFCl3_jx = 1.0f0
const σ_CFCl3_interp = create_fjx_interp(
    [220.0f0, 298.0f0],
    [
        SA_F32[2.177e-18, 1.608e-18, 1.302e-18, 9.616e-19, 4.538e-19, 1.638e-19, 9.562e-20, 7.113e-20, 4.112e-23, 7.026e-22, 1.892e-21, 0, 0, 0, 0, 0, 0, 0],
        SA_F32[2.237e-18, 1.652e-18, 1.339e-18, 9.882e-19, 4.884e-19, 2.026e-19, 1.266e-19, 9.361e-20, 9.732e-23, 1.289e-21, 2.955e-21, 0, 0, 0, 0, 0, 0, 0],
    ]
)


# NO=N+O          delta(1,0)&(0,0) scaled by 0.6 (WACCM in PhotoComp2008)
const ϕ_NO_jx = 1.0f0
const σ_NO = SA_F32[6.054e-19, 6.936e-19, 3.475e-19, 0, 0, 0, 0, 8.892e-21, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
const σ_NO_interp = [(T) -> σ_NO[i] for i in 1:18]

# CHOCHO=>HCO+HCO  Glyoxal[c] = CHOCHO => HCHO+CO      JPL10
const ϕ_Glyxlc_jx = 1.0f0
const σ_Glyxlc_interp = create_fjx_interp(
    [177.0f0, 999.0f0],
    [
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 2.968e-21, 1.178e-20, 1.5e-20, 1.704e-20, 1.636e-20, 1.462e-20, 1.148e-20, 4.054e-21, 1.431e-22, 0],
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 2.968e-21, 1.178e-20, 1.5e-20, 1.704e-20, 1.636e-20, 1.462e-20, 1.148e-20, 4.054e-21, 9.496e-23, 0],
    ]
)


# CF3CFCl2=>       CF3CFCl2 = CFC-114   JPL10
const ϕ_F114_jx = 1.0f0
const σ_F114_interp = create_fjx_interp(
    [210.0f0, 300.0f0],
    [
        SA_F32[1.291e-19, 5.656e-20, 3.403e-20, 1.585e-20, 4.491e-21, 9.317e-22, 4.48e-22, 6.972e-22, 9.733e-26, 2.472e-24, 6.812e-24, 0, 0, 0, 0, 0, 0, 0],
        SA_F32[1.203e-19, 5.566e-20, 3.52e-20, 1.825e-20, 6.078e-21, 1.642e-21, 8.743e-22, 9.578e-22, 2.385e-25, 6.055e-24, 1.608e-23, 0, 0, 0, 0, 0, 0, 0],
    ]
)


# CH3ONO2>CH3O+NO2 Methyl nitrate JPL10
const ϕ_CH3NO3_jx = 1.0f0
const σ_CH3NO3_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[4.459e-18, 8.856e-18, 1.02e-17, 1.021e-17, 6.757e-18, 3.48e-18, 2.248e-18, 1.458e-18, 2.364e-20, 4.404e-20, 6.296e-20, 3.955e-21, 1.626e-21, 6.61e-22, 2.584e-22, 2.194e-23, 8.431e-28, 0],
        SA_F32[6.234e-18, 1.238e-17, 1.426e-17, 1.427e-17, 9.445e-18, 4.866e-18, 3.143e-18, 2.038e-18, 3.219e-20, 6.182e-20, 8.858e-20, 6.137e-21, 2.74e-21, 1.236e-21, 5.475e-22, 6.702e-23, 4.217e-27, 0],
    ]
)


# CHBr3=>          CHBr3 = bromoform   JPL10
const ϕ_CHBr3_jx = 1.0f0
const σ_CHBr3_interp = create_fjx_interp(
    [210.0f0, 300.0f0],
    [
        SA_F32[5.928e-18, 4.527e-18, 4.157e-18, 3.943e-18, 4.498e-18, 4.775e-18, 5.072e-18, 5.503e-18, 7.881e-19, 9.16e-19, 5.449e-19, 7.716e-21, 1.996e-21, 6.036e-22, 2.461e-22, 3.24e-23, 1.119e-25, 0],
        SA_F32[5.63e-18, 4.299e-18, 3.948e-18, 3.744e-18, 4.272e-18, 4.535e-18, 4.817e-18, 5.227e-18, 8.969e-19, 9.049e-19, 5.338e-19, 1.676e-20, 5.066e-21, 1.783e-21, 7.401e-22, 9.748e-23, 3.365e-25, 0],
    ]
)


# CF3CCHCl2=>      CF3CCHCl2 = HCFC-123  JPL10
const ϕ_F123_jx = 1.0f0
const σ_F123_interp = create_fjx_interp(
    [210.0f0, 295.0f0],
    [
        SA_F32[9.371e-19, 5.088e-19, 3.401e-19, 1.837e-19, 6.049e-20, 1.614e-20, 8.774e-21, 9.599e-21, 4.437e-24, 6.935e-23, 1.699e-22, 0, 0, 0, 0, 0, 0, 0],
        SA_F32[1.037e-18, 6.011e-19, 4.213e-19, 2.461e-19, 8.86e-20, 2.646e-20, 1.495e-20, 1.438e-20, 1.191e-23, 1.409e-22, 3.152e-22, 0, 0, 0, 0, 0, 0, 0],
    ]
)


# CHF2Cl=>CHF2+Cl  CHF2Cl = HCFC-22    JPL10
const ϕ_CHF2Cl_jx = 1.0f0
const σ_CHF2Cl_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[7.852e-21, 2.816e-21, 1.551e-21, 6.716e-22, 1.891e-22, 4.359e-23, 2.467e-23, 3.483e-23, 0, 2.454e-27, 4.684e-25, 0, 0, 0, 0, 0, 0, 0],
        SA_F32[7.528e-21, 2.775e-21, 1.578e-21, 7.386e-22, 2.375e-22, 6.207e-23, 3.512e-23, 4.089e-23, 0, 3.494e-27, 6.67e-25, 0, 0, 0, 0, 0, 0, 0],
    ]
)


# OClO=>O+ClO      JPL10
const ϕ_OClO_jx = 1.0f0
const σ_OClO = SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 4.061e-19, 5.616e-19, 8.776e-19, 1.211e-18, 1.857e-18, 2.606e-18, 2.522e-18, 4.227e-18, 1.25e-18, 0]
const σ_OClO_interp = [(T) -> σ_OClO[i] for i in 1:18]

# CF2ClBr=>        CF2ClBr = halon 1211       JPL10
const ϕ_H1211_jx = 1.0f0
const σ_H1211_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[5.81e-19, 6.607e-19, 7.844e-19, 9.493e-19, 1.29e-18, 1.385e-18, 1.313e-18, 1.181e-18, 2.933e-20, 8.981e-20, 8.335e-20, 2.038e-23, 3.698e-24, 1.109e-24, 3.929e-25, 5.01e-27, 0, 0],
        SA_F32[5.063e-19, 5.757e-19, 6.828e-19, 8.271e-19, 1.105e-18, 1.166e-18, 1.105e-18, 9.995e-19, 3.541e-20, 8.319e-20, 7.215e-20, 1.068e-22, 2.727e-23, 8.654e-24, 3.067e-24, 3.91e-26, 0, 0],
    ]
)


# BrO=>Br+O        JPL10
const ϕ_BrO_jx = 1.0f0
const σ_BrO = SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 0, 5.62e-19, 1.202e-18, 2.008e-18, 3.239e-18, 4.52e-18, 5.064e-18, 5.809e-18, 2.408e-19, 0]
const σ_BrO_interp = [(T) -> σ_BrO[i] for i in 1:18]

# CH3Cl=>CH3+Cl    CH3Cl = Methyl chloride   JPL10
const ϕ_CH3Cl_jx = 1.0f0
const σ_CH3Cl_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[2.838e-19, 1.31e-19, 8.192e-20, 3.94e-20, 1.098e-20, 2.102e-21, 9.71e-22, 1.654e-21, 1.353e-25, 4.702e-24, 1.424e-23, 0, 0, 0, 0, 0, 0, 0],
        SA_F32[2.842e-19, 1.308e-19, 8.238e-20, 4.125e-20, 1.284e-20, 3.135e-21, 1.619e-21, 2.046e-21, 2.691e-25, 9.348e-24, 2.828e-23, 0, 0, 0, 0, 0, 0, 0],
    ]
)


# CH3COC2H5 >      Methylethyl Ketone => C2H5+CH3CO[0.85] CH3+C2H5CO[0.15]  X6
const ϕ_MEKeto_jx = 1.0f0
const σ_MEKeto_interp = create_fjx_interp(
    [177.0f0, 999.0f0],
    [
        SA_F32[0, 0, 4.387e-22, 0, 4.393e-21, 1.291e-21, 1.231e-21, 1.356e-21, 3.133e-20, 3.51e-20, 3.397e-20, 3.077e-20, 1.981e-20, 1.158e-20, 5.697e-21, 6.008e-22, 9.031e-26, 0],
        SA_F32[0, 0, 1.982e-22, 0, 1.984e-21, 5.829e-22, 5.559e-22, 6.123e-22, 1.415e-20, 1.586e-20, 1.534e-20, 1.39e-20, 8.947e-21, 5.232e-21, 2.574e-21, 2.714e-22, 4.076e-26, 0],
    ]
)


# CH3C(O)COONO2=>  PeroxyAcetylNitrate =CH3C(O)O2+NO2[0.70] =CH3C(O)O+NO3[0.30
const ϕ_PAN_jx = 1.0f0
const σ_PAN_interp = create_fjx_interp(
    [250.0f0, 298.0f0],
    [
        SA_F32[7.378e-19, 1.68e-18, 2.179e-18, 3.774e-18, 3.022e-18, 2.163e-18, 1.735e-18, 1.343e-18, 6.568e-20, 9.393e-20, 8.602e-20, 2.421e-21, 9.352e-22, 4.32e-22, 2.291e-22, 5.393e-23, 2.381e-25, 0],
        SA_F32[8.121e-19, 1.847e-18, 2.379e-18, 4.129e-18, 3.193e-18, 2.239e-18, 1.803e-18, 1.407e-18, 7.914e-20, 1.055e-19, 9.304e-20, 3.531e-21, 1.414e-21, 6.696e-22, 3.632e-22, 9.109e-23, 4.339e-25, 0],
    ]
)


# CF2BrCF2Br=>     CF2BrCF2Br = halon 2402    JPL10
const ϕ_H2402_jx = 1.0f0
const σ_H2402_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[1.038e-18, 1.137e-18, 1.167e-18, 1.193e-18, 1.145e-18, 1.013e-18, 9.185e-19, 8.183e-19, 1.999e-20, 6.373e-20, 5.698e-20, 2.284e-24, 6.437e-25, 2.441e-25, 1.156e-25, 1.879e-27, 0, 0],
        SA_F32[1.023e-18, 1.133e-18, 1.175e-18, 1.221e-18, 1.213e-18, 1.101e-18, 1.008e-18, 9.029e-19, 2.581e-20, 7.236e-20, 6.36e-20, 4.568e-23, 1.287e-23, 4.881e-24, 2.311e-24, 3.757e-26, 0, 0],
    ]
)


# C2H5CHO >C2H5+   Propionaldehyde(propanal) => C2H5+HCO  JPL10
const ϕ_PrAld_jx = 1.0f0
const σ_PrAld = SA_F32[0, 0, 1.797e-23, 0, 3.43e-22, 5.337e-22, 6.041e-22, 6.842e-22, 2.988e-20, 4.52e-20, 5.122e-20, 5.545e-20, 4.628e-20, 3.576e-20, 2.436e-20, 5.836e-21, 4.083e-24, 0]
const σ_PrAld_interp = [(T) -> σ_PrAld[i] for i in 1:18]

# CH3C(O)CH=CH2 >  MethylVinyl Ketone = CH3C(O)CH=CH2    JPL10
const ϕ_MeVK_jx = 1.0f0
const σ_MeVK_interp = create_fjx_interp(
    [177.0f0, 566.0f0, 999.0f0],
    [
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 1.8e-22, 4.76e-21, 5.468e-21, 5.834e-21, 5.243e-21, 4.484e-21, 3.669e-21, 1.71e-21, 2.97e-23, 0],
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 1.011e-22, 2.673e-21, 3.07e-21, 3.275e-21, 2.943e-21, 2.517e-21, 2.061e-21, 9.6e-22, 1.664e-23, 0],
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 6.933e-23, 1.833e-21, 2.106e-21, 2.247e-21, 2.019e-21, 1.725e-21, 1.414e-21, 6.587e-22, 1.139e-23, 0],
    ]
)


# ClNO3=ClO+NO2    JPL10
const ϕ_ClNO3b_jx = 1.0f0
const σ_ClNO3b_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[2.147e-19, 4.899e-19, 6.457e-19, 1.109e-18, 1.076e-18, 1.183e-18, 1.281e-18, 1.339e-18, 1.299e-19, 1.867e-19, 1.364e-19, 1.362e-20, 7.256e-21, 4.187e-21, 2.372e-21, 7.063e-22, 8.654e-24, 0],
        SA_F32[2.351e-19, 5.361e-19, 7.034e-19, 1.21e-18, 1.136e-18, 1.203e-18, 1.285e-18, 1.334e-18, 1.546e-19, 2.022e-19, 1.427e-19, 1.858e-20, 1.043e-20, 6.303e-21, 3.68e-21, 1.052e-21, 1.091e-23, 0],
    ]
)


# CF3CCl3=>        CF3CCl3 = CFC-113    JPL10
const ϕ_F113_jx = 1.0f0
const σ_F113_interp = create_fjx_interp(
    [210.0f0, 300.0f0],
    [
        SA_F32[1.061e-18, 6.852e-19, 4.983e-19, 2.81e-19, 8.627e-20, 2.02e-20, 1.043e-20, 1.272e-20, 4.158e-24, 6.808e-23, 1.8e-22, 0, 0, 0, 0, 0, 0, 0],
        SA_F32[9.495e-19, 5.932e-19, 4.317e-19, 2.618e-19, 9.3e-20, 2.661e-20, 1.453e-20, 1.395e-20, 6.689e-24, 1.046e-22, 2.685e-22, 0, 0, 0, 0, 0, 0, 0],
    ]
)


# HO2NO2=>HO2+NO2  JPL10 & scaled 1-micron J added to last bin
const ϕ_HNO4_jx = 1.0f0
const σ_HNO4 = SA_F32[4.631e-18, 7.576e-18, 8.407e-18, 7.479e-18, 4.744e-18, 2.854e-18, 2.217e-18, 1.843e-18, 2.799e-19, 2.21e-19, 1.539e-19, 2.685e-20, 1.171e-20, 5.683e-21, 3.102e-21, 8.269e-22, 1.309e-23, 6.723e-23]
const σ_HNO4_interp = [(T) -> σ_HNO4[i] for i in 1:18]

# ClO=>Cl+O        JPL10
const ϕ_ClO_jx = 1.0f0
const σ_ClO = SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 4.737e-18, 3.262e-18, 1.858e-18, 1.428e-18, 5.982e-19, 3.548e-19, 1.427e-19, 0, 0, 0]
const σ_ClO_interp = [(T) -> σ_ClO[i] for i in 1:18]

# H2O2=>OH+OH      JPL10
const ϕ_H2O2_jx = 1.0f0
const σ_H2O2_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[2.325e-19, 4.629e-19, 5.394e-19, 5.429e-19, 4.447e-19, 3.755e-19, 3.457e-19, 3.197e-19, 5.346e-20, 4.855e-20, 3.423e-20, 8.407e-21, 5.029e-21, 3.308e-21, 2.221e-21, 8.598e-22, 5.921e-24, 0],
        SA_F32[2.325e-19, 4.629e-19, 5.394e-19, 5.429e-19, 4.447e-19, 3.755e-19, 3.457e-19, 3.197e-19, 5.465e-20, 4.966e-20, 3.524e-20, 9.354e-21, 5.763e-21, 3.911e-21, 2.718e-21, 1.138e-21, 7.927e-24, 0],
    ]
)
const ϕ_h2o2_jx = ϕ_H2O2_jx
const σ_h2o2_interp = σ_H2O2_interp

# CH2Br2=>         CH2Br2 = dibromo methane    JPL10
const ϕ_CH2Br2_jx = 1.0f0
const σ_CH2Br2_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[2.106e-18, 2.1e-18, 2.205e-18, 2.384e-18, 2.343e-18, 2.413e-18, 2.691e-18, 3.182e-18, 1.489e-19, 4.916e-19, 3.187e-19, 5.186e-23, 8.043e-24, 0, 0, 0, 0, 0],
        SA_F32[1.983e-18, 1.977e-18, 2.076e-18, 2.244e-18, 2.206e-18, 2.272e-18, 2.404e-18, 2.561e-18, 1.724e-19, 3.708e-19, 2.349e-19, 2.67e-22, 4.141e-23, 0, 0, 0, 0, 0],
    ]
)


# OCS=>CO+S        OCS = carbonyl sulfide   JPL10
const ϕ_OCS_jx = 1.0f0
const σ_OCS_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[5.275e-20, 3.879e-20, 2.92e-20, 2.343e-20, 6.049e-20, 1.232e-19, 1.691e-19, 2.181e-19, 9.974e-21, 4.126e-20, 2.417e-20, 9.292e-25, 0, 0, 0, 0, 0, 0],
        SA_F32[7.547e-20, 5.314e-20, 3.637e-20, 2.436e-20, 6.009e-20, 1.222e-19, 1.689e-19, 2.186e-19, 1.382e-20, 4.616e-20, 2.47e-20, 3.294e-24, 7.435e-26, 0, 0, 0, 0, 0],
    ]
)


# CH3CF2Cl=>       CH3CF2Cl = HCFC-142b    JPL10
const ϕ_F142b_jx = 1.0f0
const σ_F142b_interp = create_fjx_interp(
    [210.0f0, 298.0f0],
    [
        SA_F32[2.253e-20, 8.483e-21, 4.768e-21, 2.113e-21, 6.068e-22, 1.276e-22, 6.211e-23, 9.375e-23, 0, 1.766e-25, 8.252e-25, 0, 0, 0, 0, 0, 0, 0],
        SA_F32[2.842e-20, 1.172e-20, 6.981e-21, 3.37e-21, 1.059e-21, 2.543e-22, 1.303e-22, 1.602e-22, 0, 4.755e-25, 1.964e-24, 0, 0, 0, 0, 0, 0, 0],
    ]
)


# CF3CF2Cl=>       CF3CF2Cl = CFC-115   JPL10
const ϕ_F115_jx = 1.0f0
const σ_F115 = SA_F32[6.876e-21, 2.986e-21, 1.891e-21, 1.018e-21, 3.531e-22, 1.124e-22, 6.079e-23, 6.176e-23, 5.728e-27, 6.317e-25, 1.345e-24, 0, 0, 0, 0, 0, 0, 0]
const σ_F115_interp = [(T) -> σ_F115[i] for i in 1:18]

# Qyld O3=O(1D)+O2 JPL10 3/2013   (7/2020 interp fix ~1e-4)
const ϕ_O31D_jx = 1.0f0
const σ_O31D_interp = create_fjx_interp(
    [218.0f0, 258.0f0, 298.0f0],
    [
        SA_F32[0.4896, 0.5071, 0.5266, 0.5646, 0.655, 0.7351, 0.7763, 0.8162, 0.9, 0.9, 0.8984, 0.9, 0.8948, 0.4755, 0.1023, 0.07773, 0, 0],
        SA_F32[0.4896, 0.5071, 0.5267, 0.5646, 0.6549, 0.735, 0.7763, 0.8162, 0.9, 0.9, 0.8984, 0.9, 0.8954, 0.5095, 0.1444, 0.08415, 0, 0],
        SA_F32[0.4896, 0.5072, 0.5268, 0.5647, 0.6549, 0.735, 0.7763, 0.8162, 0.9, 0.9, 0.8985, 0.9, 0.8968, 0.5709, 0.2309, 0.09735, 0, 0],
    ]
)
const ϕ_o31D_jx = ϕ_O31D_jx
const σ_o31D_interp = σ_O31D_interp

# CF3I=>CF3+I      CF3I = tri-fluoro iodo methane   JPL10
const ϕ_CF3I_jx = 1.0f0
const σ_CF3I_interp = create_fjx_interp(
    [243.0f0, 300.0f0],
    [
        SA_F32[8.229e-21, 3.885e-21, 2.404e-21, 1.653e-21, 1.772e-21, 2.856e-21, 4.178e-21, 6.205e-21, 5.54e-19, 3.605e-19, 2.108e-19, 1.33e-19, 4.992e-20, 2.036e-20, 9.495e-21, 1.567e-21, 6.121e-24, 0],
        SA_F32[8.229e-21, 3.885e-21, 2.404e-21, 1.653e-21, 1.772e-21, 2.856e-21, 4.178e-21, 6.205e-21, 5.164e-19, 3.545e-19, 2.275e-19, 1.544e-19, 6.552e-20, 2.972e-20, 1.471e-20, 2.62e-21, 1.292e-23, 0],
    ]
)


# CHOCHO=>HCO+HCO  Glyoxal[a] = CHOCHO => HCO+HCO      JPL10
const ϕ_Glyxla_jx = 1.0f0
const σ_Glyxla_interp = create_fjx_interp(
    [177.0f0, 999.0f0],
    [
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 1.173e-21, 4.354e-21, 5.353e-21, 6.331e-21, 7.09e-21, 7.781e-21, 7.452e-21, 4.214e-21, 1.9e-21, 0],
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 1.173e-21, 4.354e-21, 5.353e-21, 6.331e-21, 7.09e-21, 7.781e-21, 7.452e-21, 4.214e-21, 1.025e-21, 0],
    ]
)


# CCl4=>           CCl4 = Carbon-tetrachloride   JPL10
const ϕ_CCl4_jx = 1.0f0
const σ_CCl4_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[3.61e-18, 1.67e-18, 1.1e-18, 7.155e-19, 6.122e-19, 4.744e-19, 3.774e-19, 2.864e-19, 6.486e-22, 7.133e-21, 1.323e-20, 0, 0, 0, 0, 0, 0, 0],
        SA_F32[3.707e-18, 1.712e-18, 1.124e-18, 7.237e-19, 6.267e-19, 5.208e-19, 4.392e-19, 3.544e-19, 1.614e-21, 1.234e-20, 1.854e-20, 0, 0, 0, 0, 0, 0, 0],
    ]
)


# Cl2=>Cl+Cl       JPL10
const ϕ_Cl2_jx = 1.0f0
const σ_Cl2_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 5.482e-21, 2.419e-20, 4.98e-20, 8.413e-20, 1.393e-19, 1.878e-19, 2.265e-19, 2.524e-19, 2.388e-20, 0],
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 6.002e-21, 2.708e-20, 5.409e-20, 8.855e-20, 1.404e-19, 1.844e-19, 2.185e-19, 2.411e-19, 2.455e-20, 0],
    ]
)


# CH3I=>CH3+I      CH3I = Methyl iodide    JPL10
const ϕ_CH3I_jx = 1.0f0
const σ_CH3I_interp = create_fjx_interp(
    [243.0f0, 300.0f0],
    [
        SA_F32[0, 0, 0, 0, 0, 0, 3.958e-20, 4.769e-20, 8.257e-19, 2.154e-19, 7.349e-20, 3.351e-20, 1.034e-20, 4.435e-21, 2.272e-21, 4.393e-22, 1.328e-24, 0],
        SA_F32[0, 0, 0, 0, 0, 0, 3.958e-20, 4.769e-20, 8.376e-19, 2.438e-19, 9.355e-20, 4.68e-20, 1.529e-20, 6.532e-21, 3.328e-21, 6.783e-22, 3.401e-24, 0],
    ]
)


# HONO=>OH+NO      JPL10
const ϕ_HNO2_jx = 1.0f0
const σ_HNO2 = SA_F32[9.693e-19, 1.387e-18, 1.574e-18, 1.808e-18, 2.154e-18, 2.188e-18, 2.079e-18, 1.911e-18, 1.171e-19, 1.849e-19, 1.489e-19, 5.487e-21, 9.353e-21, 1.895e-20, 3.195e-20, 9.008e-20, 2.261e-20, 0]
const σ_HNO2_interp = [(T) -> σ_HNO2[i] for i in 1:18]

# Acetn=CH3CO+CH3  Acetone=CH3C(O)CH3          JPL10
const ϕ_Aceta_jx = 1.0f0
const σ_Aceta_interp = create_fjx_interp(
    [177.0f0, 566.0f0, 999.0f0],
    [
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 0, 2.312e-20, 2.824e-20, 1.98e-20, 5.927e-21, 6e-22, 5.868e-23, 5.934e-25, 0, 0],
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 0, 1.486e-20, 1.723e-20, 1.24e-20, 4.464e-21, 7.146e-22, 1.171e-22, 2.202e-24, 0, 0],
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 0, 1.053e-20, 1.237e-20, 9.213e-21, 3.702e-21, 7.1e-22, 1.357e-22, 3.115e-24, 0, 0],
    ]
)


# N2O=>N2+O        JPL10
const ϕ_N2O_jx = 1.0f0
const σ_N2O_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[1.076e-19, 8.49e-20, 6.974e-20, 5.2e-20, 2.171e-20, 6.405e-21, 3.097e-21, 2.473e-21, 3.308e-25, 1.048e-23, 4.114e-23, 0, 0, 0, 0, 0, 0, 0],
        SA_F32[1.225e-19, 1.006e-19, 8.493e-20, 6.658e-20, 3.135e-20, 1.15e-20, 6.282e-21, 4.614e-21, 1.501e-24, 4.07e-23, 1.173e-22, 0, 0, 0, 0, 0, 0, 0],
    ]
)


# CH3CCl3=>CH3CCl2 CH3CCl3 = Methyl chloroform   JPL10
const ϕ_MeCCl3_jx = 1.0f0
const σ_MeCCl3_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[2.72e-18, 2.025e-18, 1.688e-18, 1.354e-18, 7.434e-19, 3.348e-19, 2.026e-19, 1.374e-19, 5.478e-23, 1.44e-21, 4.133e-21, 0, 0, 0, 0, 0, 0, 0],
        SA_F32[2.722e-18, 2.02e-18, 1.68e-18, 1.338e-18, 7.42e-19, 3.515e-19, 2.232e-19, 1.574e-19, 2.197e-22, 2.618e-21, 5.521e-21, 0, 0, 0, 0, 0, 0, 0],
    ]
)


# ClOOCl=>ClO+ClO  JPL10
const ϕ_Cl2O2_jx = 1.0f0
const σ_Cl2O2 = SA_F32[0, 2.177e-19, 9.321e-19, 1.435e-18, 3.578e-18, 3.009e-18, 2.648e-18, 2.357e-18, 4.063e-18, 1.971e-18, 1.236e-18, 8.681e-19, 5.719e-19, 4.163e-19, 3.234e-19, 1.921e-19, 1.541e-20, 0]
const σ_Cl2O2_interp = [(T) -> σ_Cl2O2[i] for i in 1:18]

# CH3Br=>CH3+Br    CH3Br = Methyl bromide     JPL10
const ϕ_CH3Br_jx = 1.0f0
const σ_CH3Br_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[1.064e-18, 4.918e-19, 5.748e-19, 6.814e-19, 7.562e-19, 6.868e-19, 6.096e-19, 5.216e-19, 6.038e-21, 3.133e-20, 3.378e-20, 9.878e-26, 0, 0, 0, 0, 0, 0],
        SA_F32[1.091e-18, 5.041e-19, 5.891e-19, 6.986e-19, 7.717e-19, 6.999e-19, 6.239e-19, 5.378e-19, 7.764e-21, 3.424e-20, 3.54e-20, 7.365e-25, 0, 0, 0, 0, 0, 0],
    ]
)


# HONO2=>OH+NO2    JPL10
const ϕ_HNO3_jx = 1.0f0
const σ_HNO3_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[4.169e-18, 8.113e-18, 8.882e-18, 8.211e-18, 3.777e-18, 1.288e-18, 6.882e-19, 4.998e-19, 1.611e-20, 1.794e-20, 1.848e-20, 3.383e-21, 1.374e-21, 5.442e-22, 2.087e-22, 2.118e-23, 2.685e-26, 0],
        SA_F32[4.284e-18, 8.419e-18, 9.368e-18, 9.112e-18, 4.462e-18, 1.551e-18, 8.426e-19, 5.852e-19, 1.855e-20, 2.162e-20, 2.299e-20, 4.371e-21, 1.92e-21, 8.312e-22, 3.573e-22, 4.706e-23, 8.32e-26, 0],
    ]
)


# CF2Cl2=>CF2Cl+   CFC-12   JPL10
const ϕ_CF2Cl2_jx = 1.0f0
const σ_CF2Cl2_interp = create_fjx_interp(
    [220.0f0, 300.0f0],
    [
        SA_F32[7.001e-19, 4.044e-19, 2.78e-19, 1.507e-19, 4.142e-20, 6.627e-21, 2.876e-21, 5.395e-21, 2.698e-25, 8.767e-24, 3.338e-23, 0, 0, 0, 0, 0, 0, 0],
        SA_F32[1.012e-18, 5.844e-19, 4.02e-19, 2.178e-19, 6.291e-20, 1.258e-20, 5.959e-21, 8.779e-21, 8.279e-25, 2.466e-23, 8.261e-23, 0, 0, 0, 0, 0, 0, 0],
    ]
)


# CHOCHO=>HCO+HCO  Glyoxal[b] = CHOCHO => H2+CO+CO     JPL10
const ϕ_Glyxlb_jx = 1.0f0
const σ_Glyxlb_interp = create_fjx_interp(
    [177.0f0, 999.0f0],
    [
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 3.378e-21, 1.177e-20, 1.388e-20, 1.453e-20, 1.224e-20, 9.558e-21, 6.758e-21, 2.443e-21, 2.613e-23, 0],
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 3.378e-21, 1.177e-20, 1.388e-20, 1.453e-20, 1.224e-20, 9.558e-21, 6.758e-21, 2.443e-21, 2.613e-23, 0],
    ]
)


# CH3CFCl2=>       CH3CFCl2 = HCFC-141b    JPL10
const ϕ_F141b_jx = 1.0f0
const σ_F141b_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[1.165e-18, 7.326e-19, 5.257e-19, 3.034e-19, 9.839e-20, 2.384e-20, 1.224e-20, 1.398e-20, 2.005e-24, 6.588e-23, 2.015e-22, 0, 0, 0, 0, 0, 0, 0],
        SA_F32[1.126e-18, 7.466e-19, 5.559e-19, 3.448e-19, 1.221e-19, 3.303e-20, 1.757e-20, 1.737e-20, 5.164e-24, 1.146e-22, 3.133e-22, 0, 0, 0, 0, 0, 0, 0],
    ]
)


# O3-total         JPl10+IUPAC 2/2013 & rev 6/2017 &7/2020 (<1e-3)
const ϕ_O3_jx = 1.0f0
const σ_O3_interp = create_fjx_interp(
    [218.0f0, 258.0f0, 298.0f0],
    [
        SA_F32[5.992e-19, 4.866e-19, 4.31e-19, 3.657e-19, 3.403e-19, 4.829e-19, 6.515e-19, 9.272e-19, 8.764e-18, 3.528e-18, 1.511e-18, 7.955e-19, 2.468e-19, 8.939e-20, 3.675e-20, 4.56e-21, 2.125e-22, 2.325e-21],
        SA_F32[5.993e-19, 4.869e-19, 4.317e-19, 3.67e-19, 3.414e-19, 4.824e-19, 6.5e-19, 9.251e-19, 8.837e-18, 3.582e-18, 1.551e-18, 8.304e-19, 2.636e-19, 9.81e-20, 4.173e-20, 5.578e-21, 2.125e-22, 2.325e-21],
        SA_F32[5.994e-19, 4.873e-19, 4.323e-19, 3.681e-19, 3.425e-19, 4.82e-19, 6.485e-19, 9.231e-19, 8.904e-18, 3.632e-18, 1.589e-18, 8.626e-19, 2.791e-19, 1.062e-19, 4.634e-20, 6.52e-21, 2.125e-22, 2.325e-21],
    ]
)


# ClNO3=Cl+NO3     JPL10
const ϕ_ClNO3a_jx = 1.0f0
const σ_ClNO3a_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[3.222e-19, 7.351e-19, 9.688e-19, 1.664e-18, 1.614e-18, 1.774e-18, 1.921e-18, 2.009e-18, 1.948e-19, 2.801e-19, 2.046e-19, 2.043e-20, 1.088e-20, 6.655e-21, 4.564e-21, 2.324e-21, 3.331e-22, 0],
        SA_F32[3.527e-19, 8.041e-19, 1.055e-18, 1.815e-18, 1.705e-18, 1.805e-18, 1.928e-18, 2.001e-18, 2.319e-19, 3.033e-19, 2.14e-19, 2.787e-20, 1.564e-20, 1.002e-20, 7.086e-21, 3.404e-21, 4.302e-22, 0],
    ]
)


# Acetaldhyde      CH3CHO=CH4+CO/CH3+CHO/H+CH3CO q=0.0/0.88/0.12 IUPAC14 P2 29
const ϕ_ActAld_jx = 1.0f0
const σ_ActAld_interp = create_fjx_interp(
    [177.0f0, 566.0f0, 999.0f0],
    [
        SA_F32[0, 0, 1.989e-23, 0, 3.699e-22, 4.938e-22, 4.737e-22, 4.659e-22, 2.45e-20, 3.409e-20, 3.82e-20, 3.732e-20, 2.707e-20, 1.579e-20, 6.566e-21, 3.883e-22, 1.862e-26, 0],
        SA_F32[0, 0, 1.903e-23, 0, 3.539e-22, 4.725e-22, 4.533e-22, 4.458e-22, 2.27e-20, 2.985e-20, 3.199e-20, 2.987e-20, 1.923e-20, 9.497e-21, 3.45e-21, 1.914e-22, 1.233e-26, 0],
        SA_F32[0, 0, 1.822e-23, 0, 3.389e-22, 4.525e-22, 4.34e-22, 4.269e-22, 2.112e-20, 2.647e-20, 2.74e-20, 2.479e-20, 1.485e-20, 6.739e-21, 2.319e-21, 1.258e-22, 9.142e-27, 0],
    ]
)


# CH2Cl2=>CH2Cl+Cl CH2Cl2 = di-chloromethane   JPL10
const ϕ_CH2Cl2_jx = 1.0f0
const σ_CH2Cl2_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[1.091e-18, 6.168e-19, 4.117e-19, 2.063e-19, 5.999e-20, 1.421e-20, 7.37e-21, 9.789e-21, 4.003e-24, 5.376e-23, 1.269e-22, 0, 0, 0, 0, 0, 0, 0],
        SA_F32[9.995e-19, 5.79e-19, 3.97e-19, 2.144e-19, 7.138e-20, 2.124e-20, 1.203e-20, 1.237e-20, 7.64e-24, 1.026e-22, 2.385e-22, 0, 0, 0, 0, 0, 0, 0],
    ]
)


# O2=O+O           JPL10
const ϕ_O2_jx = 1.0f0
const σ_O2_interp = create_fjx_interp(
    [180.0f0, 260.0f0, 300.0f0],
    [
        SA_F32[1.727e-21, 1.989e-22, 3.004e-23, 9.833e-24, 7.306e-24, 6.827e-24, 6.238e-24, 5.748e-24, 1.153e-25, 5.03e-25, 4.15e-25, 0, 0, 0, 0, 0, 0, 0],
        SA_F32[2.273e-21, 3.07e-22, 4.943e-23, 1.408e-23, 7.688e-24, 6.827e-24, 6.238e-24, 5.888e-24, 1.153e-25, 5.03e-25, 4.15e-25, 0, 0, 0, 0, 0, 0, 0],
        SA_F32[2.763e-21, 4.269e-22, 7.478e-23, 2.1e-23, 8.35e-24, 6.827e-24, 6.238e-24, 5.994e-24, 1.153e-25, 5.03e-25, 4.15e-25, 0, 0, 0, 0, 0, 0, 0],
    ]
)


# BrONO2=>BrO+NO2  JPL10
const ϕ_BrNO3_jx = 1.0f0
const σ_BrNO3_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[0, 0, 5.484e-19, 7.245e-19, 3.702e-18, 3.475e-18, 3.182e-18, 2.978e-18, 5.304e-19, 6.086e-19, 4.489e-19, 1.963e-19, 1.584e-19, 1.307e-19, 1.11e-19, 8.033e-20, 1.68e-20, 0],
        SA_F32[0, 0, 8.026e-19, 1.071e-18, 5.166e-18, 4.19e-18, 3.467e-18, 3.039e-18, 5.567e-19, 5.989e-19, 4.528e-19, 2.098e-19, 1.705e-19, 1.425e-19, 1.207e-19, 8.648e-20, 1.87e-20, 0],
    ]
)


# NO2=>NO+O        JPL10 (1/2016) log-extrap 200-300K (rev.integ. over wavel)
const ϕ_NO2_jx = 1.0f0
const σ_NO2_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 1.834e-20, 4.696e-20, 7.707e-20, 1.078e-19, 1.47e-19, 1.832e-19, 2.181e-19, 3.138e-19, 1.422e-19, 0],
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 2.354e-20, 4.697e-20, 7.546e-20, 1.063e-19, 1.477e-19, 1.872e-19, 2.303e-19, 3.469e-19, 1.546e-19, 0],
    ]
)


# CH3OOH=>CH3O+OH 
const ϕ_CH3OOH_jx = 1.0f0
const σ_CH3OOH = SA_F32[0, 0, 0, 0, 0, 3.12e-19, 2.882e-19, 2.25e-19, 2.716e-20, 2.74e-20, 2.143e-20, 5.624e-21, 3.52e-21, 2.403e-21, 1.697e-21, 7.23e-22, 2.285e-23, 0]
const σ_CH3OOH_interp = [(T) -> σ_CH3OOH[i] for i in 1:18]

# HOCH2CHO >       Glycol Aldehyde => CH2OH+HCO[0.83] CH3OH+CO[0.10] OH+CH2CHO
const ϕ_GlyAld_jx = 1.0f0
const σ_GlyAld = SA_F32[0, 0, 0, 0, 7.413e-20, 1.879e-19, 1.038e-19, 4.946e-20, 5.168e-20, 5.719e-20, 5.581e-20, 4.81e-20, 2.912e-20, 1.567e-20, 7.019e-21, 7.099e-22, 0, 0]
const σ_GlyAld_interp = [(T) -> σ_GlyAld[i] for i in 1:18]

# H2CO=>H+HCO      JPL10
const ϕ_H2COa_jx = 1.0f0
const σ_H2COa_interp = create_fjx_interp(
    [223.0f0, 298.0f0],
    [
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 3.143e-21, 1.021e-20, 1.269e-20, 2.323e-20, 2.498e-20, 1.133e-20, 2.183e-20, 4.746e-21, 0, 0],
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 3.147e-21, 1.018e-20, 1.266e-20, 2.315e-20, 2.497e-20, 1.131e-20, 2.189e-20, 4.751e-21, 0, 0],
    ]
)
const ϕ_CH2Oa_jx = ϕ_H2COa_jx
const σ_CH2Oa_interp = σ_H2COa_interp

# CH3COC(O)H >HCO+ Methyl glyoxal = CH3COC(O)H => CH3CO + HCO   JPL10
const ϕ_MGlyxl_jx = 1.0f0
const σ_MGlyxl_interp = create_fjx_interp(
    [177.0f0, 566.0f0, 999.0f0],
    [
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 2.085e-20, 3.955e-20, 4.398e-20, 4.414e-20, 3.504e-20, 2.355e-20, 1.811e-20, 6.006e-21, 7.416e-21, 0],
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 2.085e-20, 3.955e-20, 4.398e-20, 4.414e-20, 3.504e-20, 2.355e-20, 1.811e-20, 6.006e-21, 3.855e-21, 0],
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 2.085e-20, 3.955e-20, 4.398e-20, 4.414e-20, 3.504e-20, 2.355e-20, 1.811e-20, 6.006e-21, 2.687e-21, 0],
    ]
)


# HOBr=>OH+Br      JPL10
const ϕ_HOBr_jx = 1.0f0
const σ_HOBr =  SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 1.324e-19, 2.011e-19, 2.202e-19, 2.196e-19, 1.726e-19, 1.367e-19, 1.157e-19, 1.125e-19, 3.274e-20, 0]
const σ_HOBr_interp = [(T) -> σ_HOBr[i] for i in 1:18]

# NO3=NO2+O/NO+O2  JPL10, calculated q(NO2+O)=0.882, q(NO+O2)=0.112
const ϕ_NO3_jx = 1.0f0
const σ_NO3_interp = create_fjx_interp(
    [190.0f0, 298.0f0],
    [
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3.622e-19, 1.48e-18],
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2.597e-19, 1.178e-18],
    ]
)


# Acetn=CH3+CH3+CO Acetone=CH3C(O)CH3           JPL10
const ϕ_Acetb_jx = 1.0f0
const σ_Acetb_interp = create_fjx_interp(
    [235.0f0, 260.0f0, 298.0f0],
    [
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 0, 3.678e-22, 2.462e-22, 1.158e-22, 2.648e-23, 6.014e-24, 1.502e-24, 4.211e-26, 0, 0],
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 0, 1.249e-21, 1.02e-21, 5.664e-22, 1.681e-22, 4.919e-23, 1.477e-23, 5.602e-25, 0, 0],
        SA_F32[0, 0, 0, 0, 0, 0, 0, 0, 0, 4.209e-21, 4.23e-21, 2.804e-21, 1.092e-21, 4.079e-22, 1.496e-22, 7.707e-24, 0, 0],
    ]
)


# BrCl=>Br+Cl      JPL10
const ϕ_BrCl_jx = 1.0f0
const σ_BrCl_interp = create_fjx_interp(
    [200.0f0, 300.0f0],
    [
        SA_F32[0, 0, 3.138e-21, 3.852e-21, 2.798e-20, 4.258e-20, 4.915e-20, 5.482e-20, 2.293e-20, 1.531e-20, 7.62e-21, 1.983e-21, 4.115e-21, 9.571e-21, 2.126e-20, 1.102e-19, 1.966e-19, 0],
        SA_F32[0, 0, 3.356e-21, 4.152e-21, 2.906e-20, 4.234e-20, 4.785e-20, 5.244e-20, 2.43e-20, 1.596e-20, 8.718e-21, 3.854e-21, 7.582e-21, 1.55e-20, 3.038e-20, 1.223e-19, 1.966e-19, 0],
    ]
)





"""
    cos_solar_zenith_angle(lat, t, long)

This function is to compute the cosine of the solar zenith angle, given the unixtime, latitude and longitude
The input variables: lat=latitude(°), long=longitude(°), t=unixtime(s)
the cosine of the solar zenith angle (SZA) is given by:                                                                            .
cos(SZA) = sin(LAT)*sin(DEC) + cos(LAT)*cos(DEC)*cos(AHR)

           where LAT = the latitude angle,
                 DEC = the solar declination angle,
                 AHR = the hour angle, all in radians.  All in radians
"""
function cos_solar_zenith_angle(t, lat, long)
    ut = Dates.unix2datetime(t)
    DOY = dayofyear(ut)
    hours = Dates.hour(ut) + Dates.minute(ut) / 60 + Dates.second(ut) / 3600
    y = Dates.year(ut)
    yearday_temp = 2 * pi * (DOY - 1 + (hours - 12) / 24)
    γ = ifelse(
        mod(y, 4) == 0, # the fraction year in radians
        yearday_temp / 366,
        yearday_temp / 365
    )
    Eot = 229.18 * (
        0.000075 + 0.001868 * cos(γ) - 0.032077 * sin(γ) - 0.014615 * cos(γ * 2) -
        0.040849 * sin(γ * 2)
    )

    timezone = floor(long / 15)
    time_offset = Eot + 4 * (long - 15 * timezone) # in minutes
    dt = floor(long / 15) # in hours
    t_local = t + dt * 3600 # in seconds
    ut_local = Dates.unix2datetime(t_local)
    tst = Dates.hour(ut_local) +
          Dates.minute(ut_local) / 60 +
          Dates.second(ut_local) / 3600 +
          time_offset / 60
    AHR = deg2rad(15) * (tst - 12) # in radians

    LAT = abs(lat * pi / 180) #lat>0, northern hemisphere; lat<0, southern hemisphere
    DEC = asin(
        sin(deg2rad(-23.44)) * cos(
        deg2rad(
        360 / 365.24 * (DOY + 10) +
        360 / pi * 0.0167 * sin(deg2rad(360 / 365.24 * (DOY - 2))),
    ),
    ),
    )
    CSZA = sin(LAT) * sin(DEC) + cos(LAT) * cos(DEC) * cos(AHR)
    return CSZA
end

@register_symbolic cos_solar_zenith_angle(t, lat, long)

# Dummy function for unit validation. Basically ModelingToolkit
# will call the function with a DynamicQuantities.Quantity or an integer to
# get information about the type and units of the output.
cos_solar_zenith_angle(t::DynamicQuantities.Quantity, lat, long) = 1.0

function calc_direct_flux(CSZA, P, i::Int)
    # calculate direct flux attenuation factor
    single_direct_flux_factor = direct_solar_beam_box_singlewavelength(
        OD_total, CSZA, z_profile, P, i)
    fluxes = top_flux[i] * single_direct_flux_factor

    return fluxes
end
@register_symbolic calc_direct_flux(CSZA, P, i::Int)

#Dummy function for unit validation.
function calc_direct_flux(
        CSZA::DynamicQuantities.Quantity,
        P::DynamicQuantities.Quantity,
        i::DynamicQuantities.Quantity
)
    1.0
end

function calc_direct_fluxes(CSZA, P)
    # calculate direct flux attenuation factor
    direct_flux_factor = direct_solar_beam_box(OD_total, CSZA, z_profile, P)
    fluxes = top_flux .* direct_flux_factor
    return fluxes
end
@register_symbolic calc_direct_fluxes(CSZA, P)

# Symbolic equations for actinic flux
function flux_eqs(csa, P)
    flux_vals = []
    flux_vars = []
    @constants c_flux = 1.0 [
        unit = u"s^-1",
        description = "Constant actinic flux (for unit conversion)"
    ]
    for i in 1:18
        wl = WL[i]
        n = Symbol("F_", Int(round(wl)))
        v = @variables $n(t) [unit = u"s^-1", description = "Actinic flux at $wl nm"]
        push!(flux_vars, only(v))
        push!(flux_vals, calc_direct_flux(csa, P, i))
    end
    flux_vars, (flux_vars .~ flux_vals .* c_flux)
end

"""
Get mean photolysis rates at different times
"""
function j_mean(σ_interp, ϕ, Temperature, fluxes)
    j = zero(Temperature)
    for i in 1:18
        j += fluxes[i] * σ_interp[i](Temperature) * ϕ
    end
    j
end

j_mean_CH2Oa(T, fluxes) = j_mean(σ_CH2Oa_interp, ϕ_CH2Oa_jx, T, fluxes)
j_mean_CH2Ob(T, fluxes) = j_mean(σ_CH2Ob_interp, ϕ_CH2Ob_jx, T, fluxes)
j_mean_o31D(T, fluxes) = j_mean(σ_o31D_interp, ϕ_o31D_jx, T, fluxes)
j_mean_HOCl(T, fluxes) = j_mean(σ_HOCl_interp, ϕ_HOCl_jx, T, fluxes)
j_mean_H2COb(T, fluxes) = j_mean(σ_H2COb_interp, ϕ_H2COb_jx, T, fluxes)
j_mean_MeAcr(T, fluxes) = j_mean(σ_MeAcr_interp, ϕ_MeAcr_jx, T, fluxes)
j_mean_N2O5(T, fluxes) = j_mean(σ_N2O5_interp, ϕ_N2O5_jx, T, fluxes)
j_mean_H1301(T, fluxes) = j_mean(σ_H1301_interp, ϕ_H1301_jx, T, fluxes)
j_mean_CFCl3(T, fluxes) = j_mean(σ_CFCl3_interp, ϕ_CFCl3_jx, T, fluxes)
j_mean_NO(T, fluxes) = j_mean(σ_NO_interp, ϕ_NO_jx, T, fluxes)
j_mean_Glyxlc(T, fluxes) = j_mean(σ_Glyxlc_interp, ϕ_Glyxlc_jx, T, fluxes)
j_mean_F114(T, fluxes) = j_mean(σ_F114_interp, ϕ_F114_jx, T, fluxes)
j_mean_CH3NO3(T, fluxes) = j_mean(σ_CH3NO3_interp, ϕ_CH3NO3_jx, T, fluxes)
j_mean_CHBr3(T, fluxes) = j_mean(σ_CHBr3_interp, ϕ_CHBr3_jx, T, fluxes)
j_mean_F123(T, fluxes) = j_mean(σ_F123_interp, ϕ_F123_jx, T, fluxes)
j_mean_CHF2Cl(T, fluxes) = j_mean(σ_CHF2Cl_interp, ϕ_CHF2Cl_jx, T, fluxes)
j_mean_OClO(T, fluxes) = j_mean(σ_OClO_interp, ϕ_OClO_jx, T, fluxes)
j_mean_H1211(T, fluxes) = j_mean(σ_H1211_interp, ϕ_H1211_jx, T, fluxes)
j_mean_BrO(T, fluxes) = j_mean(σ_BrO_interp, ϕ_BrO_jx, T, fluxes)
j_mean_CH3Cl(T, fluxes) = j_mean(σ_CH3Cl_interp, ϕ_CH3Cl_jx, T, fluxes)
j_mean_MEKeto(T, fluxes) = j_mean(σ_MEKeto_interp, ϕ_MEKeto_jx, T, fluxes)
j_mean_PAN(T, fluxes) = j_mean(σ_PAN_interp, ϕ_PAN_jx, T, fluxes)
j_mean_H2402(T, fluxes) = j_mean(σ_H2402_interp, ϕ_H2402_jx, T, fluxes)
j_mean_PrAld(T, fluxes) = j_mean(σ_PrAld_interp, ϕ_PrAld_jx, T, fluxes)
j_mean_MeVKa(T, fluxes) = j_mean(σ_MeVK_interp, ϕ_MeVK_jx, T, fluxes).*0.6
j_mean_MeVKb(T, fluxes) = j_mean(σ_MeVK_interp, ϕ_MeVK_jx, T, fluxes).*0.2
j_mean_MeVKc(T, fluxes) = j_mean(σ_MeVK_interp, ϕ_MeVK_jx, T, fluxes).*0.2
j_mean_ClNO3b(T, fluxes) = j_mean(σ_ClNO3b_interp, ϕ_ClNO3b_jx, T, fluxes)
j_mean_F113(T, fluxes) = j_mean(σ_F113_interp, ϕ_F113_jx, T, fluxes)
j_mean_HNO4(T, fluxes) = j_mean(σ_HNO4_interp, ϕ_HNO4_jx, T, fluxes)
j_mean_ClO(T, fluxes) = j_mean(σ_ClO_interp, ϕ_ClO_jx, T, fluxes)
j_mean_H2O2(T, fluxes) = j_mean(σ_H2O2_interp, ϕ_H2O2_jx, T, fluxes)
j_mean_CH2Br2(T, fluxes) = j_mean(σ_CH2Br2_interp, ϕ_CH2Br2_jx, T, fluxes)
j_mean_OCS(T, fluxes) = j_mean(σ_OCS_interp, ϕ_OCS_jx, T, fluxes)
j_mean_F142b(T, fluxes) = j_mean(σ_F142b_interp, ϕ_F142b_jx, T, fluxes)
j_mean_F115(T, fluxes) = j_mean(σ_F115_interp, ϕ_F115_jx, T, fluxes)
j_mean_O31D(T, fluxes) = j_mean(σ_O31D_interp, ϕ_O31D_jx, T, fluxes)
j_mean_CF3I(T, fluxes) = j_mean(σ_CF3I_interp, ϕ_CF3I_jx, T, fluxes)
j_mean_Glyxla(T, fluxes) = j_mean(σ_Glyxla_interp, ϕ_Glyxla_jx, T, fluxes)
j_mean_CCl4(T, fluxes) = j_mean(σ_CCl4_interp, ϕ_CCl4_jx, T, fluxes)
j_mean_Cl2(T, fluxes) = j_mean(σ_Cl2_interp, ϕ_Cl2_jx, T, fluxes)
j_mean_CH3I(T, fluxes) = j_mean(σ_CH3I_interp, ϕ_CH3I_jx, T, fluxes)
j_mean_HNO2(T, fluxes) = j_mean(σ_HNO2_interp, ϕ_HNO2_jx, T, fluxes)
j_mean_Aceta(T, fluxes) = j_mean(σ_Aceta_interp, ϕ_Aceta_jx, T, fluxes)
j_mean_N2O(T, fluxes) = j_mean(σ_N2O_interp, ϕ_N2O_jx, T, fluxes)
j_mean_MeCCl3(T, fluxes) = j_mean(σ_MeCCl3_interp, ϕ_MeCCl3_jx, T, fluxes)
j_mean_Cl2O2(T, fluxes) = j_mean(σ_Cl2O2_interp, ϕ_Cl2O2_jx, T, fluxes)
j_mean_CH3Br(T, fluxes) = j_mean(σ_CH3Br_interp, ϕ_CH3Br_jx, T, fluxes)
j_mean_HNO3(T, fluxes) = j_mean(σ_HNO3_interp, ϕ_HNO3_jx, T, fluxes)
j_mean_CF2Cl2(T, fluxes) = j_mean(σ_CF2Cl2_interp, ϕ_CF2Cl2_jx, T, fluxes)
j_mean_Glyxlb(T, fluxes) = j_mean(σ_Glyxlb_interp, ϕ_Glyxlb_jx, T, fluxes)
j_mean_F141b(T, fluxes) = j_mean(σ_F141b_interp, ϕ_F141b_jx, T, fluxes)
j_mean_O3(T, fluxes) = j_mean(σ_O3_interp, ϕ_O3_jx, T, fluxes)
j_mean_ClNO3a(T, fluxes) = j_mean(σ_ClNO3a_interp, ϕ_ClNO3a_jx, T, fluxes)
j_mean_ActAld(T, fluxes) = j_mean(σ_ActAld_interp, ϕ_ActAld_jx, T, fluxes)
j_mean_CH2Cl2(T, fluxes) = j_mean(σ_CH2Cl2_interp, ϕ_CH2Cl2_jx, T, fluxes)
j_mean_O2(T, fluxes) = j_mean(σ_O2_interp, ϕ_O2_jx, T, fluxes)
j_mean_BrNO3(T, fluxes) = j_mean(σ_BrNO3_interp, ϕ_BrNO3_jx, T, fluxes)
j_mean_NO2(T, fluxes) = j_mean(σ_NO2_interp, ϕ_NO2_jx, T, fluxes)
j_mean_CH3OOH(T, fluxes) = j_mean(σ_CH3OOH_interp, ϕ_CH3OOH_jx, T, fluxes)
j_mean_GlyAld(T, fluxes) = j_mean(σ_GlyAld_interp, ϕ_GlyAld_jx, T, fluxes)
j_mean_H2COa(T, fluxes) = j_mean(σ_H2COa_interp, ϕ_H2COa_jx, T, fluxes)
j_mean_MGlyxl(T, fluxes) = j_mean(σ_MGlyxl_interp, ϕ_MGlyxl_jx, T, fluxes)
j_mean_HOBr(T, fluxes) = j_mean(σ_HOBr_interp, ϕ_HOBr_jx, T, fluxes)
j_mean_NO3a(T, fluxes) = j_mean(σ_NO3_interp, ϕ_NO3_jx, T, fluxes) .* 0.886
j_mean_NO3b(T, fluxes) = j_mean(σ_NO3_interp, ϕ_NO3_jx, T, fluxes) .* 0.114
j_mean_Acetb(T, fluxes) = j_mean(σ_Acetb_interp, ϕ_Acetb_jx, T, fluxes)
j_mean_BrCl(T, fluxes) = j_mean(σ_BrCl_interp, ϕ_BrCl_jx, T, fluxes)


"""
    adjust_j_o31D(T, P, H2O)

Adjust the photolysis rate of O3 -> O2 + O(1D) to represent the effective rate for O3 -> 2OH.
This adjustment is based on the fraction of O(1D) that reacts with H2O to produce 2 OH.
"""
function adjust_j_o31D(T, P, H2O)
    @constants(T_unit=1,
        [unit=u"K", description="unit of Temperature"],
        A=6.02e23,
        [unit=u"molec/mol", description="Avogadro's number"],
        R=8.314e6,
        [unit=u"(Pa*cm^3)/(K*mol)", description="universal gas constant"],
        ppb_unit=1e-9,
        [unit=u"ppb", description="Convert from mol/mol_air to ppb"],
        num_density_unit_inv=1,
        [
            unit=u"cm^3/molec",
            description="multiply by num_density to obtain the unitless value of num_density"
        ],
        ppb_inv=1,
        [unit=u"ppb^-1"],)
    num_density_unitless = A*P/(R*T)*num_density_unit_inv

    # Define species concentrations value in unit of molec/cm3
    C_H2O = H2O*ppb_inv*1e-9*num_density_unitless # convert value of H2O concentration in unit of ppb to unit of molec/cm3, but here is unitless
    C_O2 = 0.2095 * num_density_unitless
    C_N2 = 0.7808 * num_density_unitless
    C_H2 = 0.5e-6 * num_density_unitless

    # Define rate constants for reactions involving O(1D)
    RO1DplH2O = 1.63e-10 * exp(60.0*T_unit / T) * C_H2O
    RO1DplH2 = 1.2e-10 * C_H2
    RO1DplN2 = 2.15e-11 * exp(110.0*T_unit / T) * C_N2
    RO1DplO2 = 3.30e-11 * exp(55.0*T_unit / T) * C_O2

    # Total rate constant for O(1D)
    RO1D = RO1DplH2O + RO1DplH2 + RO1DplN2 + RO1DplO2

    # Prevent division by zero
    return RO1DplH2O / RO1D
end

struct FastJXCoupler
    sys::Any
end

"""
Description: This is a box model used to calculate the photolysis reaction rate constant using the Fast-JX scheme
(Neu, J. L., Prather, M. J., and Penner, J. E. (2007), Global atmospheric chemistry: Integrating over fractional cloud cover, J. Geophys. Res., 112, D11306, doi:10.1029/2006JD008007.)

Argument:

  - `t_ref`: Reference time for the model, can be a `DateTime` or a Unix timestamp (in seconds).

# Example

Build Fast-JX model:

```julia
fj = FastJX(DateTime(2000, 1, 1))
```
"""
function FastJX(t_ref::AbstractFloat; name = :FastJX)
    @constants T_unit = 1.0 [
        unit = u"K",
        description = "Unit temperature (for unit conversion)"
    ]
    @parameters T = 298.0 [unit = u"K", description = "Temperature"]
    @parameters lat = 40.0 [description = "Latitude (Degrees)"]
    @parameters long = -97.0 [description = "Longitude (Degrees)"]
    @parameters P = 101325 [unit = u"Pa", description = "Pressure"]
    @constants P_unit = 1.0 [unit = u"Pa", description = "Unit pressure"]
    @parameters H2O = 450 [unit = u"ppb"]
    @parameters t_ref = t_ref [unit = u"s", description = "Reference Unix time"]

    @variables j_h2o2(t) [unit = u"s^-1"]
    @variables j_CH2Oa(t) [unit = u"s^-1"]
    @variables j_CH2Ob(t) [unit = u"s^-1"]
    @variables j_o31D(t) [unit = u"s^-1"]
    @variables j_o32OH(t) [unit = u"s^-1"]
    @variables j_NO2(t) [unit = u"s^-1"]
    @variables cosSZA(t) [description = "Cosine of the solar zenith angle"]
    @variables j_HOCl(t) [unit = u"s^-1"]
    @variables j_H2COb(t) [unit = u"s^-1"]
    @variables j_MeAcr(t) [unit = u"s^-1"]
    @variables j_N2O5(t) [unit = u"s^-1"]
    @variables j_H1301(t) [unit = u"s^-1"]
    @variables j_CFCl3(t) [unit = u"s^-1"]
    @variables j_NO(t) [unit = u"s^-1"]
    @variables j_Glyxlc(t) [unit = u"s^-1"]
    @variables j_F114(t) [unit = u"s^-1"]
    @variables j_CH3NO3(t) [unit = u"s^-1"]
    @variables j_CHBr3(t) [unit = u"s^-1"]
    @variables j_F123(t) [unit = u"s^-1"]
    @variables j_CHF2Cl(t) [unit = u"s^-1"]
    @variables j_OClO(t) [unit = u"s^-1"]
    @variables j_H1211(t) [unit = u"s^-1"]
    @variables j_BrO(t) [unit = u"s^-1"]
    @variables j_CH3Cl(t) [unit = u"s^-1"]
    @variables j_MEKeto(t) [unit = u"s^-1"]
    @variables j_PAN(t) [unit = u"s^-1"]
    @variables j_H2402(t) [unit = u"s^-1"]
    @variables j_PrAld(t) [unit = u"s^-1"]
    @variables j_MeVKa(t) [unit = u"s^-1"]
    @variables j_MeVKb(t) [unit = u"s^-1"]
    @variables j_MeVKc(t) [unit = u"s^-1"]
    @variables j_ClNO3b(t) [unit = u"s^-1"]
    @variables j_F113(t) [unit = u"s^-1"]
    @variables j_HNO4(t) [unit = u"s^-1"]
    @variables j_ClO(t) [unit = u"s^-1"]
    @variables j_H2O2(t) [unit = u"s^-1"]
    @variables j_CH2Br2(t) [unit = u"s^-1"]
    @variables j_OCS(t) [unit = u"s^-1"]
    @variables j_F142b(t) [unit = u"s^-1"]
    @variables j_F115(t) [unit = u"s^-1"]
    @variables j_O31D(t) [unit = u"s^-1"]
    @variables j_CF3I(t) [unit = u"s^-1"]
    @variables j_Glyxla(t) [unit = u"s^-1"]
    @variables j_CCl4(t) [unit = u"s^-1"]
    @variables j_Cl2(t) [unit = u"s^-1"]
    @variables j_CH3I(t) [unit = u"s^-1"]
    @variables j_HNO2(t) [unit = u"s^-1"]
    @variables j_Aceta(t) [unit = u"s^-1"]
    @variables j_N2O(t) [unit = u"s^-1"]
    @variables j_MeCCl3(t) [unit = u"s^-1"]
    @variables j_Cl2O2(t) [unit = u"s^-1"]
    @variables j_CH3Br(t) [unit = u"s^-1"]
    @variables j_HNO3(t) [unit = u"s^-1"]
    @variables j_CF2Cl2(t) [unit = u"s^-1"]
    @variables j_Glyxlb(t) [unit = u"s^-1"]
    @variables j_F141b(t) [unit = u"s^-1"]
    @variables j_O3(t) [unit = u"s^-1"]
    @variables j_ClNO3a(t) [unit = u"s^-1"]
    @variables j_ActAld(t) [unit = u"s^-1"]
    @variables j_CH2Cl2(t) [unit = u"s^-1"]
    @variables j_O2(t) [unit = u"s^-1"]
    @variables j_BrNO3(t) [unit = u"s^-1"]
    @variables j_CH3OOH(t) [unit = u"s^-1"]
    @variables j_GlyAld(t) [unit = u"s^-1"]
    @variables j_H2COa(t) [unit = u"s^-1"]
    @variables j_MGlyxl(t) [unit = u"s^-1"]
    @variables j_HOBr(t) [unit = u"s^-1"]
    @variables j_NO3a(t) [unit = u"s^-1"]
    @variables j_NO3b(t) [unit = u"s^-1"]
    @variables j_Acetb(t) [unit = u"s^-1"]
    @variables j_BrCl(t) [unit = u"s^-1"]


    flux_vars, fluxeqs = flux_eqs(cosSZA, P/P_unit)

    eqs = [cosSZA ~ cos_solar_zenith_angle(t + t_ref, lat, long);
           fluxeqs;
           j_h2o2 ~ j_mean_H2O2(T/T_unit, flux_vars)*0.0557; #0.0557 is a parameter to adjust the calculated H2O2 photolysis to appropriate magnitudes.
           j_CH2Oa ~ j_mean_CH2Oa(T/T_unit, flux_vars)*0.945; #0.945 is a parameter to adjust the calculated CH2Oa photolysis to appropriate magnitudes.
           j_CH2Ob ~ j_mean_CH2Ob(T/T_unit, flux_vars)*0.813; #0.813 is a parameter to adjust the calculated CH2Ob photolysis to appropriate magnitudes.
           j_o31D ~ j_mean_o31D(T/T_unit, flux_vars)*2.33e-21; #2.33e-21 is a parameter to adjust the calculated O(^3)1D photolysis to appropriate magnitudes.
           j_o32OH ~ j_o31D*adjust_j_o31D(T, P, H2O);
           j_CH3OOH ~ j_mean_CH3OOH(T/T_unit, flux_vars)*0.0931; #0.0931 is a parameter to adjust the calculated CH3OOH photolysis to appropriate magnitudes.
           j_NO2 ~ j_mean_NO2(T/T_unit, flux_vars)*0.444
           j_HOCl ~ j_mean_HOCl(T/T_unit, flux_vars);
           j_H2COb ~ j_mean_H2COb(T/T_unit, flux_vars);
           j_MeAcr ~ j_mean_MeAcr(T/T_unit, flux_vars);
           j_N2O5 ~ j_mean_N2O5(T/T_unit, flux_vars);
           j_H1301 ~ j_mean_H1301(T/T_unit, flux_vars);
           j_CFCl3 ~ j_mean_CFCl3(T/T_unit, flux_vars);
           j_NO ~ j_mean_NO(T/T_unit, flux_vars);
           j_Glyxlc ~ j_mean_Glyxlc(T/T_unit, flux_vars);
           j_F114 ~ j_mean_F114(T/T_unit, flux_vars);
           j_CH3NO3 ~ j_mean_CH3NO3(T/T_unit, flux_vars);
           j_CHBr3 ~ j_mean_CHBr3(T/T_unit, flux_vars);
           j_F123 ~ j_mean_F123(T/T_unit, flux_vars);
           j_CHF2Cl ~ j_mean_CHF2Cl(T/T_unit, flux_vars);
           j_OClO ~ j_mean_OClO(T/T_unit, flux_vars);
           j_H1211 ~ j_mean_H1211(T/T_unit, flux_vars);
           j_BrO ~ j_mean_BrO(T/T_unit, flux_vars);
           j_CH3Cl ~ j_mean_CH3Cl(T/T_unit, flux_vars);
           j_MEKeto ~ j_mean_MEKeto(T/T_unit, flux_vars);
           j_PAN ~ j_mean_PAN(T/T_unit, flux_vars);
           j_H2402 ~ j_mean_H2402(T/T_unit, flux_vars);
           j_PrAld ~ j_mean_PrAld(T/T_unit, flux_vars);
           j_MeVKa ~ j_mean_MeVKa(T/T_unit, flux_vars);
           j_MeVKb ~ j_mean_MeVKb(T/T_unit, flux_vars);
           j_MeVKc ~ j_mean_MeVKc(T/T_unit, flux_vars);
           j_ClNO3b ~ j_mean_ClNO3b(T/T_unit, flux_vars);
           j_F113 ~ j_mean_F113(T/T_unit, flux_vars);
           j_HNO4 ~ j_mean_HNO4(T/T_unit, flux_vars);
           j_ClO ~ j_mean_ClO(T/T_unit, flux_vars);
           j_H2O2 ~ j_mean_H2O2(T/T_unit, flux_vars);
           j_CH2Br2 ~ j_mean_CH2Br2(T/T_unit, flux_vars);
           j_OCS ~ j_mean_OCS(T/T_unit, flux_vars);
           j_F142b ~ j_mean_F142b(T/T_unit, flux_vars);
           j_F115 ~ j_mean_F115(T/T_unit, flux_vars);
           j_O31D ~ j_mean_O31D(T/T_unit, flux_vars);
           j_CF3I ~ j_mean_CF3I(T/T_unit, flux_vars);
           j_Glyxla ~ j_mean_Glyxla(T/T_unit, flux_vars);
           j_CCl4 ~ j_mean_CCl4(T/T_unit, flux_vars);
           j_Cl2 ~ j_mean_Cl2(T/T_unit, flux_vars);
           j_CH3I ~ j_mean_CH3I(T/T_unit, flux_vars);
           j_HNO2 ~ j_mean_HNO2(T/T_unit, flux_vars);
           j_Aceta ~ j_mean_Aceta(T/T_unit, flux_vars);
           j_MeCCl3 ~ j_mean_MeCCl3(T/T_unit, flux_vars);
           j_Cl2O2 ~ j_mean_Cl2O2(T/T_unit, flux_vars);
           j_CH3Br ~ j_mean_CH3Br(T/T_unit, flux_vars);
           j_HNO3 ~ j_mean_HNO3(T/T_unit, flux_vars);
           j_CF2Cl2 ~ j_mean_CF2Cl2(T/T_unit, flux_vars);
           j_Glyxlb ~ j_mean_Glyxlb(T/T_unit, flux_vars);
           j_F141b ~ j_mean_F141b(T/T_unit, flux_vars);
           j_O3 ~ j_mean_O3(T/T_unit, flux_vars);
           j_ClNO3a ~ j_mean_ClNO3a(T/T_unit, flux_vars);
           j_ActAld ~ j_mean_ActAld(T/T_unit, flux_vars);
           j_CH2Cl2 ~ j_mean_CH2Cl2(T/T_unit, flux_vars);
           j_O2 ~ j_mean_O2(T/T_unit, flux_vars);
           j_BrNO3 ~ j_mean_BrNO3(T/T_unit, flux_vars);
           j_GlyAld ~ j_mean_GlyAld(T/T_unit, flux_vars);
           j_H2COa ~ j_mean_H2COa(T/T_unit, flux_vars);
           j_MGlyxl ~ j_mean_MGlyxl(T/T_unit, flux_vars);
           j_HOBr ~ j_mean_HOBr(T/T_unit, flux_vars);
           j_NO3a ~ j_mean_NO3a(T/T_unit, flux_vars);
           j_NO3b ~ j_mean_NO3b(T/T_unit, flux_vars);
           j_Acetb ~ j_mean_Acetb(T/T_unit, flux_vars);
           j_BrCl ~ j_mean_BrCl(T/T_unit, flux_vars);

           ]

    ODESystem(
        eqs,
        t,
        [j_h2o2, j_CH2Oa, j_CH2Ob, j_o31D, j_o32OH, j_CH3OOH, j_NO2, 
        j_HOCl, j_H2COb, j_MeAcr, j_N2O5, j_H1301, j_CFCl3, j_NO, j_Glyxlc, j_F114, 
        j_CH3NO3, j_CHBr3, j_F123, j_CHF2Cl, j_OClO, j_H1211, j_BrO, j_CH3Cl, j_MEKeto, 
        j_PAN, j_H2402, j_PrAld, j_MeVKa, j_MeVKb, j_MeVKc, j_ClNO3b, j_F113, j_HNO4, j_ClO, j_H2O2, j_CH2Br2, 
        j_OCS, j_F142b, j_F115, j_O31D, j_CF3I, j_Glyxla, j_CCl4, j_Cl2, j_CH3I, j_HNO2, 
        j_Aceta, j_MeCCl3, j_Cl2O2, j_CH3Br, j_HNO3, j_CF2Cl2, j_Glyxlb, 
        j_F141b, j_O3, j_ClNO3a, j_ActAld, j_CH2Cl2, j_O2, j_BrNO3, 
        j_GlyAld, j_H2COa, j_MGlyxl, j_HOBr, j_NO3a, j_NO3b, j_Acetb, j_BrCl, 
        cosSZA, flux_vars...],
        [lat, long, T, P, H2O, t_ref];
        name = name,
        metadata = Dict(:coupletype => FastJXCoupler)
    )
end
FastJX(t_ref::DateTime; kwargs...) = FastJX(datetime2unix(t_ref); kwargs...)

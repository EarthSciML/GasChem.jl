# Hybrid grid parameters from https://wiki.seas.harvard.edu/geos-chem/index.php/GEOS-Chem_vertical_grids
const Ap = SVector{73}([
    0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01, 1.961311e+01, 2.609201e+01,
    3.257081e+01, 3.898201e+01, 4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01,
    7.062198e+01, 7.883422e+01, 8.909992e+01, 9.936521e+01, 1.091817e+02, 1.189586e+02,
    1.286959e+02, 1.429100e+02, 1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02,
    2.032590e+02, 2.121500e+02, 2.187760e+02, 2.238980e+02, 2.243630e+02, 2.168650e+02,
    2.011920e+02, 1.769300e+02, 1.503930e+02, 1.278370e+02, 1.086630e+02, 9.236572e+01,
    7.851231e+01, 6.660341e+01, 5.638791e+01, 4.764391e+01, 4.017541e+01, 3.381001e+01,
    2.836781e+01, 2.373041e+01, 1.979160e+01, 1.645710e+01, 1.364340e+01, 1.127690e+01,
    9.292942e+00, 7.619842e+00, 6.216801e+00, 5.046801e+00, 4.076571e+00, 3.276431e+00,
    2.620211e+00, 2.084970e+00, 1.650790e+00, 1.300510e+00, 1.019440e+00, 7.951341e-01,
    6.167791e-01, 4.758061e-01, 3.650411e-01, 2.785261e-01, 2.113490e-01, 1.594950e-01,
    1.197030e-01, 8.934502e-02, 6.600001e-02, 4.758501e-02, 3.270000e-02, 2.000000e-02,
    1.000000e-02] .* 100) # Pa

const Bp = @SVector [
    1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01, 9.203870e-01, 8.989080e-01,
    8.774290e-01, 8.560180e-01, 8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01,
    7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01, 6.158184e-01, 5.810415e-01,
    5.463042e-01, 4.945902e-01, 4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01,
    2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01, 6.372006e-02, 2.801004e-02,
    6.960025e-03, 8.175413e-09, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
    0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
    0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
    0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
    0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
    0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
    0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
    0.000000e+00]

const P_levels = Ap + Bp .* 101325 # Pa
const P_levels_mid = (P_levels[1:end-1] + P_levels[2:end]) / 2 # Pa
const T_profile = [279.57358, 279.2307, 278.61282, 277.94315, 277.26285, 276.6288, 276.03024, 275.4809, 274.95572, 274.52994, 274.11856, 273.64145, 272.94287, 271.98056, 270.7896, 269.47098, 268.0806, 266.63275, 264.8483, 262.46146, 259.6963, 256.6932, 253.46152, 250.04655, 246.3538, 242.32751, 238.16446, 233.11511, 227.47636, 222.1113, 217.82246, 215.36874, 212.70734, 210.53001, 208.99005, 207.88849, 207.85678, 208.95131, 210.57327, 212.48679, 213.93152, 215.19081, 216.08836, 216.68036, 217.47514, 218.05592, 218.91273, 220.38663, 221.94107, 223.34564, 224.70229, 227.0882, 229.8271, 233.01927, 237.65991, 244.24698, 250.52005, 255.36069, 257.55826, 257.14658, 253.79321, 248.75366, 243.75272, 239.24557, 236.23053, 234.71906, 233.12137, 230.1704, 225.42912, 220.58194, 216.58136, 214.11511]
const T_profile_top = [279.57358, 279.2307, 278.61282, 277.94315, 277.26285, 276.6288, 276.03024, 275.4809, 274.95572, 274.52994, 274.11856, 273.64145, 272.94287, 271.98056, 270.7896, 269.47098, 268.0806, 266.63275, 264.8483, 262.46146, 259.6963, 256.6932, 253.46152, 250.04655, 246.3538, 242.32751, 238.16446, 233.11511, 227.47636, 222.1113, 217.82246, 215.36874, 212.70734, 210.53001, 208.99005, 207.88849, 207.85678, 208.95131, 210.57327, 212.48679, 213.93152, 215.19081, 216.08836, 216.68036, 217.47514, 218.05592, 218.91273, 220.38663, 221.94107, 223.34564, 224.70229, 227.0882, 229.8271, 233.01927, 237.65991, 244.24698, 250.52005, 255.36069, 257.55826, 257.14658, 253.79321, 248.75366, 243.75272, 239.24557, 236.23053, 234.71906, 233.12137, 230.1704, 225.42912, 220.58194, 216.58136, 214.11511, 214.11511]

function path_density(P)
    A = 6.022140857e23# description = "Avogadro's number", particles/mol
    Air_M = 28.97 # molecular weight of air, g/mol
    g0 = 9.80665 # gravity, m/s²

    MASFAC = A / (Air_M * g0 * 10) #the factor 10 accounting for converting from (kg/m²) to (g/cm²) (i.e. combining the 1000 g/kg and 10⁴ cm²/m² factors)
    path_density = zeros(length(P))
    for i in 1:length(P)-1
        path_density[i] = MASFAC * (P[i] - P[i+1])
    end
    path_density[end] = MASFAC * P[end]
    return path_density #in unit of molecules/cm^2
end

# Cross sections σ for Rayleigh scattering for different wavelengths(18 bins), temperatures
const σ_Raylay = SA_F32[5.073, 4.479, 4.196, 3.906, 3.355, 2.929, 2.736, 2.581, 1.049, 9.492*0.1, 8.103*0.1, 6.131*0.1, 5.422*0.1, 4.923*0.1, 4.514*0.1, 3.643*0.1, 2.087*0.1, 3.848*0.01] * 10.0f0^-25.0f0

function Rayleigh_OD(N)
    OD = σ_Raylay * N
    return OD
end

# Cross sections σ for O2 for different wavelengths(18 bins), temperatures
const σ_O2_interp = create_fjx_interp([180.0f0, 260.0f0, 300.0f0], [
    SA_F32[1.727, 1.989*0.1, 3.004*0.01, 9.833*0.001, 7.306*0.001, 6.827*0.001, 6.238*0.001, 5.748*0.001, 1.153*0.0001, 5.030*0.0001, 4.150*0.0001, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000] * 10.0f0^-21.0f0,
    SA_F32[1.727, 1.989*0.1, 3.004*0.01, 9.833*0.001, 7.306*0.001, 6.827*0.001, 6.238*0.001, 5.748*0.001, 1.153*0.0001, 5.030*0.0001, 4.150*0.0001, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000] * 10.0f0^-21.0f0,
    SA_F32[2.763, 4.269*0.1, 7.478*0.01, 2.100*0.01, 8.350*0.001, 6.827*0.001, 6.238*0.001, 5.994*0.001, 1.153*0.0001, 5.030*0.0001, 4.150*0.0001, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000] * 10.0f0^-21.0f0,
])

function OD_abs(T, N)
    XQO2 = [σ_O2_interp[i](T) for i in 1:18]
    OD_O2 = XQO2 * N * 0.20948
    OD_abs = OD_O2 # TODO not included O3 absorption optical depth right now
    return OD_abs
end

function ZHL(P, T)
    A = 6.022140857e23 # description = "Avogadro's number", [particles/mol]
    Air_M = 28.97 # molecular weight of air, [g/mol]
    g0 = 9.80665 # gravity, [m/s²]
    BOLTZ = 1.38064852e-23 #Boltzmann's constant [J/K]
    ZZHT = 5e5 #scale height [cm]

    MASFAC = A / (Air_M * g0 * 10) #the factor 10 accounting for converting from (kg/m²) to (g/cm²) (i.e. combining the 1000 g/kg and 10⁴ cm²/m² factors)
    Z_CLIM = zeros(Float64, length(P) + 1)
    for i in 1:length(T)
        SCALEH = BOLTZ * 1e6 * MASFAC * T[i] # 1e6 is a conversion factor, converts Boltzmann's constant from SI units (Pa·m³/K) into units of (Pa·cm³/K)
        Z_CLIM[i+1] = Z_CLIM[i] - (log(P[i+1] / P[i]) * SCALEH)
    end
    Z_CLIM[end] = Z_CLIM[end-1] + ZZHT #An extra altitude increment (ZZHT) is added at the end to define the very top of the profile.
    return Z_CLIM
end

# using default P_levels, T_profile
const N_profile = SVector{73}(path_density(P_levels)) #molecules/cm^2 in each layer
const z_profile = SVector{74}(ZHL(P_levels, T_profile))

# calculate optical depth
const OD_ray_profile = hcat(Rayleigh_OD.(N_profile)...) # Rayleigh OD in each layer
const OD_abs_profile = SMatrix{18,73}(hcat(OD_abs.(T_profile_top, N_profile)...)) # O2 absorption OD in each layer
const OD_total = SMatrix{74,18}(
    vcat((OD_abs_profile + OD_ray_profile)', zeros(1, 18))
    # Build the extended optical depth grid by appending a zero row (top-of-atmosphere).
)
# TODO: add O3 absorption OD & aerosols and cloud OD


# Convert CTM heights to absolute radii
function calcRZ(ZHL, i)
    # Earth's radius in cm.
    RAD = 6375e5
    RAD + ZHL[i]
end
# Build the fine vertical grid RZ2 with 2*L1U+1 points:
# Odd indices correspond to CTM edges; even indices are midpoints.
function calcRZ2(ZHL, i)
    L1U = 73 # L1U: The CTM grid dimension (levels+1).
    if i == 2 * L1U + 1
        return calcRZ(ZHL, L1U + 1)
    elseif i % 2 == 0 # midpoint
        j = i ÷ 2
        return 0.5 * (calcRZ(ZHL, j) + calcRZ(ZHL, j + 1))
    else # edge
        return calcRZ(ZHL, (i + 1) ÷ 2)
    end
end

# --- Ascending (Upward) Calculation ---
function ascending_mass_factor(XMU1, ZHL, I)
    RZ2I, RZ2Ip1, RQ2 = sphere2info(ZHL, I)
    XMU2 = sqrt(1.0 - RQ2 * (1.0 - XMU1^2))
    diff = RZ2Ip1 - RZ2I
    AMF2I = (RZ2Ip1 * XMU2 - RZ2I * XMU1) / diff
    XMU1 = XMU2
    return AMF2I, XMU1
end

# --- Descending (Twilight) Calculation (only for U0 < 0) ---
function descending_mass_factor(XMU1, ZHL, II)
    RZ2II, RZ2IIp1, RQ2 = sphere2info(ZHL, II)
    diff = RZ2IIp1 - RZ2II
    DIFF = RZ2IIp1 * sqrt(1.0 - XMU1^2) - RZ2II
    done = false
    if II == 1
        DIFF = max(DIFF, 0.0)
    end
    if DIFF < 0.0
        XMU2 = sqrt(1.0 - (1.0 - XMU1^2) / RQ2)
        XL = abs(RZ2IIp1 * XMU1 - RZ2II * XMU2)
        AMF2II = 2.0 * XL / diff
        XMU1 = XMU2
    else
        AMF2II = 2.0 * RZ2IIp1 * XMU1 / diff
        done = true
    end
    return AMF2II, XMU1, done
end

function sphere2info(ZHL, i)
    # Build the fine vertical grid RZ2 with 2*L1U+1 points:
    # Odd indices correspond to CTM edges; even indices are midpoints.
    RZ2i = calcRZ2(ZHL, i)
    RZ2ip1 = calcRZ2(ZHL, i+1)

    # Pre-calculate squared ratios for each sublayer.
    RQ2 = (RZ2i / RZ2ip1)^2

    return RZ2i, RZ2ip1, RQ2
end

# Compute shadow height for sun below the horizon.
shadht(U0, RZ21) = ifelse(U0 < 0.0, RZ21 / sqrt(1.0 - U0^2), 0.0)

# Calculate row J in the air mass factor matrix using spherical geometry
function sphere2J(U0, ZHL, J)
    # U0 is cosine of the solar zenith angle
    # J is row in matrix to calculate for.
    L2 = 73 * 2
    # Define the size of the output grid. (Original code sets n = 2*LJX1U+1.)
    n = 2 * 73 + 1
    AMF2 = MVector{n,Float64}(undef)
    AMF2 .= 0.0

    if calcRZ2(ZHL, J) < shadht(U0, calcRZ2(ZHL, 1))
        return AMF2
    end

    # --- Ascending (Upward) Calculation ---
    XMU1 = abs(U0)
    for I in J:L2
        AMF2[I], XMU1 = ascending_mass_factor(XMU1, ZHL, I)
    end
    AMF2[L2+1] = 1.0

    # --- Descending (Twilight) Calculation (only for U0 < 0) ---
    if U0 >= 0.0
        return AMF2
    end

    XMU1 = abs(U0)
    for II in (J-1):-1:1
        AMF2[II], XMU1, done = descending_mass_factor(XMU1, ZHL, II)
        if done
            break
        end
    end
    return AMF2
end

function find_closest_pressure_index(pressure, P_levels)
    n = length(P_levels)
    # For each index i in 1:n, we add 1 if pressure <= P_levels[i].
    # Then the interpolation index is 1 plus that sum.
    index = ifelse(pressure > P_levels[1], 1, sum(ifelse(pressure - P_levels[i] <= 0, 1, 0) for i in 1:n))

    return index
end

function direct_solar_beam_box(DTAU, U0, ZHL, P; threshold::Float64=76.0)
    # Calculate the direct attenuated solar beam (i.e. the unscattered or "direct" beam)by computing an effective optical depth along the slanted light path and applying the Beer–Lambert law.
    # DTAU: Optical depth per CTM layer (dimensions: n_layers+1 × n_wave)
    # U0: Cosine of the solar zenith angle.
    # ZHL: Vector (length L1U+1) containing the heights (in cm) of the bottom edge of each CTM level and the top-of-atmosphere.
    # P: Pressure (unit: Pa)
    # threshold: A cutoff value for the effective optical depth (default 76.0)

    closest_index = find_closest_pressure_index(P, P_levels)
    fine_index = 2 * closest_index - 1

    # Determine dimensions from DTAU.
    n_layers, n_wave = size(DTAU)
    n_layers -= 1
    # n_edge is assumed to be the number of CTM grid edges.

    AMF2J = sphere2J(U0, ZHL, fine_index)

    # Define the number of fine-grid points.
    ng = 2 * n_layers + 1
    FTAU2 = zeros(Float64, n_wave)
    if P < 1.0
        FTAU2[:] .= 1.0  # Top-of-atmosphere: exp(0) = 1.
    else
        # For each wavelength and fine-grid point (except the top level),
        # compute the effective optical depth using the AMF2 weights via a comprehension.
        for k in 1:n_wave
            if AMF2J[fine_index] > 0.0
                tau = 0.5 * sum(DTAU[div(i + 1, 2), k] * AMF2J[i] for i in 1:ng)
                FTAU2[k] = tau < threshold ? exp(-tau) : 0.0
            end
        end
    end
    return FTAU2
end

function direct_solar_beam_box_singlewavelength(DTAU, U0, ZHL, P, k; threshold::Float64=76.0)
    # Calculate the direct attenuated solar beam (i.e. the unscattered or "direct" beam)
    # by computing an effective optical depth along the slanted light path and applying
    # the Beer–Lambert law.
    # DTAU: Optical depth per CTM layer (dimensions: n_layers+1 × n_wave)
    # U0: Cosine of the solar zenith angle.
    # ZHL: Vector (length L1U+1) containing the heights (in cm) of the bottom edge of each CTM level and the top-of-atmosphere.
    # P: Pressure (unit: Pa)
    # threshold: A cutoff value for the effective optical depth (default 76.0)

    closest_index = find_closest_pressure_index(P, P_levels)
    fine_index = 2 * closest_index - 1

    # Determine dimensions from DTAU.
    n_layers = size(DTAU, 1) - 1
    # n_edge is assumed to be the number of CTM grid edges.

    AMF2J = sphere2J(U0, ZHL, fine_index)

    # Define the number of fine-grid points.
    ng = 2 * n_layers + 1
    FTAU2 = 0
    if P < 1.0
        FTAU2 = 1.0  # Top-of-atmosphere: exp(0) = 1.
    else
        # For certain single wavelength and fine-grid point (except the top level),
        # compute the effective optical depth using the AMF2 weights via a comprehension.
        if AMF2J[fine_index] > 0.0
            tau = 0.5 * sum(DTAU[div(i + 1, 2), Int(k)] * AMF2J[i] for i in 1:ng)
            FTAU2 = tau < threshold ? exp(-tau) : 0.0
        end
    end
    return FTAU2
end

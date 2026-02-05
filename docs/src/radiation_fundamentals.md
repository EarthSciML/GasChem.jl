# Radiation Fundamentals

## Overview

This module implements the fundamental atmospheric radiation equations from Chapter 4 "Atmospheric Radiation and Photochemistry" of Seinfeld and Pandis (2006). These equations form the foundation for understanding radiative transfer in the atmosphere, including solar radiation, thermal emission, and planetary energy balance.

The equations describe:

  - **Photon energy relationships** (Eq. 4.1): The quantum nature of electromagnetic radiation
  - **Planck's blackbody radiation law** (Eq. 4.2): Spectral distribution of thermal radiation
  - **Wien's displacement law** (Eq. 4.3): Peak wavelength of blackbody emission
  - **Stefan-Boltzmann law** (Eq. 4.4): Total thermal radiation power
  - **Planetary energy balance** (Eqs. 4.5-4.7): Earth's equilibrium temperature
  - **Climate sensitivity** (Eqs. 4.8-4.10): Temperature response to radiative forcing
  - **Top-of-atmosphere radiative forcing** (Eq. 4.11): Net energy flux at TOA

**Reference**: Seinfeld, J.H. and Pandis, S.N. (2006). *Atmospheric Chemistry and Physics: From Air Pollution to Climate Change*, 2nd Edition. John Wiley & Sons, Inc., Hoboken, New Jersey. Chapter 4, pp. 98-105.

```@docs
PhotonEnergy
BlackbodyRadiation
WienDisplacement
StefanBoltzmann
PlanetaryEnergyBalance
ClimateSensitivity
TOARadiativeForcing
RadiationFundamentals
```

### Physical Constants

The following fundamental constants are used throughout:

| Constant                     | Symbol | Value          | Units       |
|:---------------------------- |:------ |:-------------- |:----------- |
| Planck's constant            | h      | 6.626 x 10^-34 | J s         |
| Speed of light               | c      | 2.9979 x 10^8  | m/s         |
| Boltzmann constant           | k      | 1.381 x 10^-23 | J/K         |
| Stefan-Boltzmann constant    | sigma  | 5.671 x 10^-8  | W m^-2 K^-4 |
| Wien's displacement constant | b      | 2.897 x 10^-3  | m K         |

* * *

## Implementation

The radiation fundamentals are implemented as a set of ModelingToolkit.jl component systems. Each component encapsulates a specific physical relationship and can be used independently or combined into larger systems.

### Component Systems

The module provides seven individual component systems plus one composed system that combines them all:

 1. **PhotonEnergy** - Photon energy-frequency-wavelength relations (Eq. 4.1)
 2. **BlackbodyRadiation** - Planck's law for spectral radiance (Eq. 4.2)
 3. **WienDisplacement** - Peak emission wavelength (Eq. 4.3)
 4. **StefanBoltzmann** - Total emissive power (Eq. 4.4)
 5. **PlanetaryEnergyBalance** - Earth's energy balance (Eqs. 4.5-4.7)
 6. **ClimateSensitivity** - Temperature response to forcing (Eqs. 4.8-4.10)
 7. **TOARadiativeForcing** - Net flux at top of atmosphere (Eq. 4.11)
 8. **RadiationFundamentals** - Composed system combining all components

### PhotonEnergy System

Implements the photon energy equation (Eq. 4.1):

```math
\Delta\varepsilon = h\nu = \frac{hc}{\lambda}
```

where h is Planck's constant, nu is frequency, c is the speed of light, and lambda is wavelength.

```@example radiation
using GasChem
using DataFrames, ModelingToolkit, Symbolics, DynamicQuantities, NonlinearSolve
using ModelingToolkit: t
using GasChem: PhotonEnergy, BlackbodyRadiation, WienDisplacement, StefanBoltzmann
using GasChem: PlanetaryEnergyBalance, ClimateSensitivity, TOARadiativeForcing,
               RadiationFundamentals

@named photon = PhotonEnergy()
sys = mtkcompile(photon)

vars = unknowns(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars],
    :Description => [ModelingToolkit.getdescription(v) for v in vars]
)
```

**Parameters:**

```@example radiation
params = parameters(sys)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params],
    :Default => [ModelingToolkit.getdefault(p) for p in params],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params],
    :Description => [ModelingToolkit.getdescription(p) for p in params]
)
```

**Equations:**

```@example radiation
equations(photon)
```

### BlackbodyRadiation System

Implements Planck's blackbody radiation law (Eq. 4.2):

```math
F_B(\lambda) = \frac{2\pi c^2 h \lambda^{-5}}{\exp\left(\frac{ch}{k\lambda T}\right) - 1}
```

This equation gives the monochromatic emissive power of a blackbody at temperature T for a given wavelength lambda.

```@example radiation
@named blackbody = BlackbodyRadiation()
sys_bb = mtkcompile(blackbody)

vars_bb = unknowns(sys_bb)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars_bb],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars_bb],
    :Description => [ModelingToolkit.getdescription(v) for v in vars_bb]
)
```

**Parameters:**

```@example radiation
params_bb = parameters(sys_bb)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params_bb],
    :Default => [ModelingToolkit.getdefault(p) for p in params_bb],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params_bb],
    :Description => [ModelingToolkit.getdescription(p) for p in params_bb]
)
```

**Equations:**

```@example radiation
equations(blackbody)
```

### WienDisplacement System

Implements Wien's displacement law (Eq. 4.3):

```math
\lambda_{max} = \frac{2.897 \times 10^{-3}}{T}
```

This gives the wavelength at which the blackbody emission spectrum peaks for a given temperature T.

```@example radiation
@named wien = WienDisplacement()
sys_wien = mtkcompile(wien)

vars_wien = unknowns(sys_wien)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars_wien],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars_wien],
    :Description => [ModelingToolkit.getdescription(v) for v in vars_wien]
)
```

**Parameters:**

```@example radiation
params_wien = parameters(sys_wien)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params_wien],
    :Default => [ModelingToolkit.getdefault(p) for p in params_wien],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params_wien],
    :Description => [ModelingToolkit.getdescription(p) for p in params_wien]
)
```

**Equations:**

```@example radiation
equations(wien)
```

### StefanBoltzmann System

Implements the Stefan-Boltzmann law (Eq. 4.4):

```math
F_B = \sigma T^4
```

The total emissive power of a blackbody integrated over all wavelengths.

```@example radiation
@named sb = StefanBoltzmann()
sys_sb = mtkcompile(sb)

vars_sb = unknowns(sys_sb)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars_sb],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars_sb],
    :Description => [ModelingToolkit.getdescription(v) for v in vars_sb]
)
```

**Parameters:**

```@example radiation
params_sb = parameters(sys_sb)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params_sb],
    :Default => [ModelingToolkit.getdefault(p) for p in params_sb],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params_sb],
    :Description => [ModelingToolkit.getdescription(p) for p in params_sb]
)
```

**Equations:**

```@example radiation
equations(sb)
```

### PlanetaryEnergyBalance System

Implements the planetary energy balance equations (Eqs. 4.5-4.7):

```math
F_S = \frac{S_0}{4}(1 - R_p) \quad \text{(Eq. 4.5)}
```

```math
F_L = \sigma T_e^4 \quad \text{(Eq. 4.6)}
```

```math
T_e = \left[\frac{(1-R_p)S_0}{4\sigma}\right]^{1/4} \quad \text{(Eq. 4.7)}
```

At equilibrium, the absorbed solar flux equals the emitted longwave flux: F_S = F_L.

```@example radiation
@named balance = PlanetaryEnergyBalance()
sys_balance = mtkcompile(balance)

vars_balance = unknowns(sys_balance)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars_balance],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars_balance],
    :Description => [ModelingToolkit.getdescription(v) for v in vars_balance]
)
```

**Parameters:**

```@example radiation
params_balance = parameters(sys_balance)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params_balance],
    :Default => [ModelingToolkit.getdefault(p) for p in params_balance],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params_balance],
    :Description => [ModelingToolkit.getdescription(p) for p in params_balance]
)
```

**Equations:**

```@example radiation
equations(balance)
```

### ClimateSensitivity System

Implements the climate sensitivity equations (Eqs. 4.8-4.10):

```math
\Delta F_{net} = \Delta F_S - \Delta F_L \quad \text{(Eq. 4.8)}
```

```math
\Delta T_e = \lambda_0 \Delta F_{net} \quad \text{(Eq. 4.9)}
```

```math
\lambda_0 = \frac{1}{4\sigma T_e^3} = \frac{T_e}{4F_L} \quad \text{(Eq. 4.10)}
```

```@example radiation
@named sensitivity = ClimateSensitivity()
sys_sensitivity = mtkcompile(sensitivity)

vars_sensitivity = unknowns(sys_sensitivity)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars_sensitivity],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars_sensitivity],
    :Description => [ModelingToolkit.getdescription(v) for v in vars_sensitivity]
)
```

**Parameters:**

```@example radiation
params_sensitivity = parameters(sys_sensitivity)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params_sensitivity],
    :Default => [ModelingToolkit.getdefault(p) for p in params_sensitivity],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params_sensitivity],
    :Description => [ModelingToolkit.getdescription(p) for p in params_sensitivity]
)
```

**Equations:**

```@example radiation
equations(sensitivity)
```

### TOARadiativeForcing System

Implements the top of atmosphere radiative forcing based on Eq. 4.11 from Seinfeld & Pandis:

```math
F_{net} = \frac{S_0}{4}(1 - R_p) - F_L
```

Note: Seinfeld & Pandis write Eq. 4.11 as ``-F_{net} = S_0/4(1-R_p) - F_L``, using the convention that ``-F_{net}`` represents net downward flux. Here we define ``F_{net}`` as net incoming flux directly, so positive ``F_{net}`` indicates the planet is gaining energy (warming).

```@example radiation
@named toa = TOARadiativeForcing()
sys_toa = mtkcompile(toa)

vars_toa = unknowns(sys_toa)
DataFrame(
    :Name => [string(Symbolics.tosymbol(v, escape = false)) for v in vars_toa],
    :Units => [dimension(ModelingToolkit.get_unit(v)) for v in vars_toa],
    :Description => [ModelingToolkit.getdescription(v) for v in vars_toa]
)
```

**Parameters:**

```@example radiation
params_toa = parameters(sys_toa)
DataFrame(
    :Name => [string(Symbolics.tosymbol(p, escape = false)) for p in params_toa],
    :Default => [ModelingToolkit.getdefault(p) for p in params_toa],
    :Units => [dimension(ModelingToolkit.get_unit(p)) for p in params_toa],
    :Description => [ModelingToolkit.getdescription(p) for p in params_toa]
)
```

**Equations:**

```@example radiation
equations(toa)
```

### RadiationFundamentals Composed System

The composed system combines all individual components:

```@example radiation
@named radiation = RadiationFundamentals()

subsystems = ModelingToolkit.get_systems(radiation)
DataFrame(
    :Subsystem => [string(nameof(s)) for s in subsystems],
    :Equations => ["Eq. 4.1", "Eq. 4.2", "Eq. 4.3", "Eq. 4.4",
        "Eqs. 4.5-4.7", "Eqs. 4.8-4.10", "Eq. 4.11"]
)
```

* * *

## Analysis

This section demonstrates the physics captured by the radiation fundamentals equations through numerical examples and visualizations that reproduce key results from Seinfeld and Pandis (2006).

### Blackbody Radiation Spectra (cf. Figures 4.2-4.3)

The Planck function describes the spectral distribution of radiation emitted by a blackbody. The following plot shows blackbody spectra at different temperatures, demonstrating how the peak wavelength shifts according to Wien's law and total power increases with temperature according to the Stefan-Boltzmann law.

```@example radiation
using Plots

# Use the BlackbodyRadiation MTK component to compute spectra
@named bb_analysis = BlackbodyRadiation()
compiled_bb = mtkcompile(bb_analysis)

# Also use WienDisplacement for peak wavelengths
@named wien_analysis = WienDisplacement()
compiled_wien_analysis = mtkcompile(wien_analysis)

# Wavelength range (in meters)
lambda_range = 10 .^ range(-7, -4, length = 500)  # 100 nm to 100 um

# Temperatures to plot
temperatures = [5800, 1000, 500, 300]
colors = [:orange, :red, :darkred, :blue]
labels = ["Sun (5800 K)", "1000 K", "500 K", "Earth (300 K)"]

p1 = plot(xlabel = "Wavelength (m)", ylabel = "Spectral Radiance (W/m³)",
    title = "Blackbody Radiation Spectra (Planck's Law, Eq. 4.2)",
    xscale = :log10, yscale = :log10,
    xlims = (1e-7, 1e-4), ylims = (1e0, 1e14),
    legend = :topright, size = (700, 500))

bb_prob = NonlinearProblem(compiled_bb, Dict(); build_initializeprob = false)
wien_prob = NonlinearProblem(compiled_wien_analysis, Dict(); build_initializeprob = false)

for (T, col, lab) in zip(temperatures, colors, labels)
    # Compute spectrum using BlackbodyRadiation component
    F = Float64[]
    for λ in lambda_range
        prob_i = remake(bb_prob, p = [compiled_bb.T => Float64(T), compiled_bb.λ => λ])
        sol_i = solve(prob_i)
        push!(F, sol_i[compiled_bb.F_B_λ])
    end
    plot!(p1, lambda_range, F, label = lab, color = col, linewidth = 2)

    # Mark peak wavelength using WienDisplacement component
    wien_prob_i = remake(wien_prob, p = [compiled_wien_analysis.T => Float64(T)])
    wien_sol = solve(wien_prob_i)
    lambda_max = wien_sol[compiled_wien_analysis.λ_max]

    bb_peak_prob = remake(bb_prob, p = [
        compiled_bb.T => Float64(T), compiled_bb.λ => lambda_max])
    bb_peak_sol = solve(bb_peak_prob)
    F_max = bb_peak_sol[compiled_bb.F_B_λ]
    scatter!(p1, [lambda_max], [F_max], color = col, markersize = 6, label = "")
end

# Add Wien's displacement law annotation
annotate!(p1, 5e-7, 1e12, text("Wien: lambda_max = 2897/T", 8, :left))

savefig(p1, "blackbody_spectra.svg")
p1
```

**Figure 1**: Blackbody radiation spectra at different temperatures (cf. S&P Figures 4.2-4.3). The dots mark the peak wavelength given by Wien's displacement law (Eq. 4.3). The Sun at 5800 K peaks in the visible range (~500 nm), while Earth at ~300 K peaks in the thermal infrared (~10 um). Note: These are ideal blackbody curves; Figure 4.2 in S&P also shows the observed solar spectrum with Fraunhofer absorption lines.

### Wien's Displacement Law Demonstration (Eq. 4.3)

```@example radiation
# Demonstrate Wien's displacement law
@named wien_sys = WienDisplacement()
compiled_wien = mtkcompile(wien_sys)

temperatures_K = [5800, 3000, 1000, 500, 300, 255]
peak_wavelengths_nm = Float64[]

for T in temperatures_K
    prob = NonlinearProblem(compiled_wien, Dict(); build_initializeprob = false)
    prob = remake(prob, p = [compiled_wien.T => T])
    sol = solve(prob)
    push!(peak_wavelengths_nm, sol[compiled_wien.λ_max] * 1e9)  # Convert to nm
end

DataFrame(
    :Temperature_K => temperatures_K,
    :Peak_Wavelength_nm => round.(peak_wavelengths_nm, digits = 1),
    :Spectral_Region => [
        "Visible (yellow)", "Near-IR", "Near-IR", "Mid-IR", "Thermal IR", "Thermal IR"]
)
```

### Planetary Energy Balance (Eqs. 4.5-4.7)

The equilibrium temperature of Earth can be calculated from the balance between absorbed solar radiation and emitted thermal radiation.

```@example radiation
@named energy_sys = PlanetaryEnergyBalance()
compiled_energy = mtkcompile(energy_sys)

# Solve with default parameters (S_0 = 1370 W/m², R_p = 0.3)
prob = NonlinearProblem(compiled_energy, Dict(compiled_energy.T_e => 250.0); build_initializeprob = false)
sol = solve(prob)

println("Planetary Energy Balance Results:")
println("================================")
println("Solar constant (S_0): 1370 W/m²")
println("Planetary albedo (R_p): 0.3")
println("Absorbed solar flux (F_S): ", round(sol[compiled_energy.F_S], digits = 2), " W/m²")
println("Emitted longwave flux (F_L): ", round(sol[compiled_energy.F_L], digits = 2), " W/m²")
println("Equilibrium temperature (T_e): ", round(sol[compiled_energy.T_e], digits = 1), " K")
println("Equilibrium temperature (T_e): ", round(sol[compiled_energy.T_e] - 273.15, digits = 1), " °C")
```

The calculated equilibrium temperature of ~255 K (-18 C) is significantly colder than Earth's actual average surface temperature of ~288 K (15 C). This 33 K difference is due to the greenhouse effect, which is not included in this simple energy balance model.

### Effect of Albedo on Equilibrium Temperature

```@example radiation
albedos = 0.0:0.05:0.8
T_equilibrium = Float64[]

for R_p in albedos
    prob = NonlinearProblem(compiled_energy, Dict(compiled_energy.T_e => 250.0); build_initializeprob = false)
    prob = remake(prob, p = [compiled_energy.R_p => R_p])
    sol = solve(prob)
    push!(T_equilibrium, sol[compiled_energy.T_e])
end

p2 = plot(albedos, T_equilibrium,
    xlabel = "Planetary Albedo (R_p)",
    ylabel = "Equilibrium Temperature (K)",
    title = "Effect of Albedo on Planetary Temperature (Eq. 4.7)",
    label = "T_e = [(1-R_p)S_0/(4σ)]^(1/4)",
    linewidth = 2, color = :blue,
    size = (600, 400))

# Mark Earth's approximate values
scatter!(p2, [0.3], [255.0], color = :red, markersize = 8, label = "Earth (R_p ≈ 0.3)")

# Add horizontal line for freezing point
hline!(p2, [273.15], linestyle = :dash, color = :gray, label = "Freezing (273 K)")

savefig(p2, "albedo_temperature.svg")
p2
```

**Figure 2**: Equilibrium temperature as a function of planetary albedo. Higher albedo reflects more solar radiation, resulting in a colder equilibrium temperature. Earth with an albedo of ~0.3 has an equilibrium temperature of ~255 K.

### Climate Sensitivity Analysis (Eqs. 4.8-4.10)

The climate sensitivity parameter lambda_0 describes how the equilibrium temperature changes in response to radiative forcing.

```@example radiation
@named climate_sys = ClimateSensitivity()
compiled_climate = mtkcompile(climate_sys)

# Calculate climate sensitivity at Earth's equilibrium temperature
prob = NonlinearProblem(compiled_climate, Dict(); build_initializeprob = false)
sol = solve(prob)

lambda_0 = sol[compiled_climate.λ_0]
println("Climate Sensitivity Analysis (no-feedback case):")
println("================================================")
println("Reference temperature (T_e): 255 K")
println("Climate sensitivity (λ_0): ", round(lambda_0, digits = 4), " K/(W/m²)")
println("")
println("For a radiative forcing of 4 W/m² (typical CO2 doubling):")
println("Temperature change (ΔT_e): ", round(sol[compiled_climate.ΔT_e], digits = 2), " K")
```

### Temperature Response to Radiative Forcing

```@example radiation
# Calculate temperature response for different forcing levels
forcings = 0:0.5:10  # W/m²
delta_T = Float64[]

for ΔF in forcings
    prob = NonlinearProblem(compiled_climate, Dict(); build_initializeprob = false)
    prob = remake(prob, p = [compiled_climate.ΔF_S => ΔF])
    sol = solve(prob)
    push!(delta_T, sol[compiled_climate.ΔT_e])
end

p3 = plot(forcings, delta_T,
    xlabel = "Radiative Forcing (W/m²)",
    ylabel = "Temperature Change (K)",
    title = "Climate Sensitivity (Eq. 4.9): ΔT = λ₀ × ΔF",
    label = "No-feedback response",
    linewidth = 2, color = :red,
    size = (600, 400))

# Mark typical CO2 doubling forcing
vline!(p3, [4.0], linestyle = :dash, color = :gray, label = "CO₂ doubling (~4 W/m²)")
scatter!(p3, [4.0], [4.0 * lambda_0], color = :blue, markersize = 8, label = "ΔT ≈ 1.1 K")

savefig(p3, "climate_sensitivity.svg")
p3
```

**Figure 3**: Temperature change as a function of radiative forcing for the no-feedback case. The climate sensitivity lambda_0 = 1/(4*sigma*T_e^3) = 0.266 K/(W/m^2) at T_e = 255 K (S&P approximate this as ~0.3). For a CO2 doubling forcing of ~4.6 W/m^2 (S&P p. 105), this gives ΔT_e ≈ 1.2 K with the exact lambda_0 (S&P report ~1.4 K using the rounded value). Note: Real climate sensitivity is higher due to positive feedbacks (water vapor, ice-albedo, etc.).

### Stefan-Boltzmann Law: Total Emissive Power

```@example radiation
@named sb_sys = StefanBoltzmann()
compiled_sb = mtkcompile(sb_sys)

temperatures = 200:10:400  # K
emissive_power = Float64[]

for T in temperatures
    prob = NonlinearProblem(compiled_sb, Dict(); build_initializeprob = false)
    prob = remake(prob, p = [compiled_sb.T => T])
    sol = solve(prob)
    push!(emissive_power, sol[compiled_sb.F_B])
end

p4 = plot(temperatures, emissive_power,
    xlabel = "Temperature (K)",
    ylabel = "Total Emissive Power (W/m²)",
    title = "Stefan-Boltzmann Law (Eq. 4.4): F_B = σT⁴",
    label = "Blackbody emission",
    linewidth = 2, color = :orange,
    size = (600, 400))

# Mark key temperatures using the StefanBoltzmann component
sb_prob = NonlinearProblem(compiled_sb, Dict(); build_initializeprob = false)

sb_eq_prob = remake(sb_prob, p = [compiled_sb.T => 255.0])
sb_eq_sol = solve(sb_eq_prob)
F_eq = sb_eq_sol[compiled_sb.F_B]

sb_surf_prob = remake(sb_prob, p = [compiled_sb.T => 288.0])
sb_surf_sol = solve(sb_surf_prob)
F_surf = sb_surf_sol[compiled_sb.F_B]

scatter!(p4, [255], [F_eq], color = :blue, markersize = 8, label = "Earth eq. (255 K)")
scatter!(p4, [288], [F_surf], color = :red, markersize = 8, label = "Earth surface (288 K)")

savefig(p4, "stefan_boltzmann.svg")
p4
```

**Figure 4**: Total blackbody emissive power as a function of temperature. The T^4 dependence means that small temperature increases lead to significantly higher emission, which acts as a stabilizing feedback in the climate system.

### Photon Energy Across the Electromagnetic Spectrum

```@example radiation
@named photon_sys = PhotonEnergy()
compiled_photon = mtkcompile(photon_sys)

# Wavelengths spanning UV to thermal IR
wavelengths_m = [1e-8, 1e-7, 4e-7, 5e-7, 7e-7, 1e-6, 1e-5, 1e-4]
wavelength_labels = ["X-ray (10 nm)", "UV (100 nm)", "Violet (400 nm)", "Green (500 nm)",
    "Red (700 nm)", "Near-IR (1 μm)", "Thermal IR (10 μm)", "Far-IR (100 μm)"]

energies_eV = Float64[]
frequencies_Hz = Float64[]

for λ in wavelengths_m
    prob = NonlinearProblem(compiled_photon, Dict(); build_initializeprob = false)
    prob = remake(prob, p = [compiled_photon.λ => λ])
    sol = solve(prob)
    push!(energies_eV, sol[compiled_photon.Δε] / 1.602e-19)  # Convert J to eV
    push!(frequencies_Hz, sol[compiled_photon.ν])
end

DataFrame(
    :Wavelength => wavelength_labels,
    :Wavelength_m => wavelengths_m,
    :Frequency_Hz => frequencies_Hz,
    :Energy_eV => round.(energies_eV, digits = 4)
)
```

This table demonstrates the inverse relationship between wavelength and photon energy (Eq. 4.1). UV and visible photons have sufficient energy to break chemical bonds and drive photochemistry, while infrared photons primarily contribute to thermal processes.

* * *

## Summary

The radiation fundamentals equations from Seinfeld and Pandis Chapter 4 provide the physical basis for understanding:

 1. **Quantum nature of light**: Photon energy is proportional to frequency (Eq. 4.1)
 2. **Thermal emission spectra**: Planck's law describes spectral distribution (Eq. 4.2)
 3. **Peak emission wavelength**: Wien's law relates peak wavelength to temperature (Eq. 4.3)
 4. **Total thermal power**: Stefan-Boltzmann law gives T^4 dependence (Eq. 4.4)
 5. **Planetary temperature**: Energy balance determines equilibrium temperature (Eqs. 4.5-4.7)
 6. **Climate response**: Sensitivity relates temperature change to forcing (Eqs. 4.8-4.10)
 7. **Radiative imbalance**: TOA flux determines warming or cooling tendency (Eq. 4.11)

The ModelingToolkit.jl implementation allows these equations to be easily combined, modified, and solved for various atmospheric and climate applications.

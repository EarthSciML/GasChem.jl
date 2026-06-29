export IsorropiaOp

"""
    IsorropiaOp(; k_mt = 1.0e-2)

An [`EarthSciMLBase.Operator`](@ref) that couples the GEOS-Chem gas phase to the
ISORROPIA II inorganic aerosol thermodynamic-equilibrium model
([`Aerosol.IsorropiaEquilibrium`](@ref)) as an operator-split (Strang) step,
repartitioning total ammonium (`NH3` ⇌ `NH4`) and total nitrate (`HNO3` ⇌ `NIT`)
toward thermodynamic equilibrium given the sulfate load, temperature, and humidity.

This is the EarthSciML idiom for an instantaneous-equilibrium process: rather than
adding the ~30 nonlinear equilibrium constraints to the gas-phase ODE as a stiff DAE,
the operator solves the equilibrium per grid cell once per Strang step (with a
warm-started Newton solve) and returns a mass-transfer tendency

    du(X) = k_mt * (X_eq - X)

that relaxes each species toward its equilibrium partition while conserving the
totals (`d[HNO3]/dt + d[NIT]/dt = 0`, `d[NH3]/dt + d[NH4]/dt = 0`). This mirrors
GEOS-Chem's dynamic gas–aerosol mass transfer. The relaxation rate `k_mt` [s⁻¹]
should be fast relative to transport but slow enough to remain non-stiff in the
operator slot.

The actual implementation lives in the `AerosolExt` package extension (it requires
`Aerosol`); load both `GasChem` and `Aerosol` to use it:

```julia
using GasChem, Aerosol
model = couple(GEOSChemGasPhase(), IsorropiaOp(), domain)
```

The operator reads `T`, `P`, `H2O`, `HNO3`, `NIT`, `NH3`, `NH4`, `SO4` from the
GEOS-Chem state; relative humidity is derived from `H2O`/`T`/`P` (Tetens saturation
vapour pressure). Sea-salt (`Na`/`Cl`) and crustal (`Ca`/`K`/`Mg`) totals are taken
as zero (the port does not yet speciate them), so the coupled subsystem is the
dominant fine-mode NH₄–NO₃–SO₄(–H₂O) aerosol.
"""
struct IsorropiaOp <: EarthSciMLBase.Operator
    k_mt::Float64
end
IsorropiaOp(; k_mt = 1.0e-2) = IsorropiaOp(k_mt)

# [Methane Oxidation](@id methane_oxidation)

Implementation of the methane oxidation mechanism, based on Seinfeld & Pandis
Chapter 6, Section 6.4 (pp. 219-227) and Table 6.1.

## Background

Methane (CH₄) is the most abundant hydrocarbon in the atmosphere (~1.8 ppm).
Its oxidation is the archetypal VOC oxidation mechanism and proceeds through
several intermediates:

```
CH₄ → CH₃O₂ → CH₃O → HCHO → HCO → CO → CO₂
```

Each step can produce ozone when NOx is present through peroxy + NO reactions.

## Table 6.1: Complete Mechanism

| # | Reaction | k at 298 K | Units |
|---|----------|------------|-------|
| 1 | CH₄ + OH → CH₃ + H₂O | 6.3 × 10⁻¹⁵ | cm³/s |
| 2 | CH₃ + O₂ + M → CH₃O₂ + M | 1.0 × 10⁻³⁰ [M] | cm⁶/s |
| 3 | CH₃O₂ + NO → CH₃O + NO₂ | 7.7 × 10⁻¹² | cm³/s |
| 4 | CH₃O₂ + HO₂ → CH₃OOH + O₂ | 5.2 × 10⁻¹² | cm³/s |
| 5 | CH₃O₂ + CH₃O₂ → products | 3.5 × 10⁻¹³ | cm³/s |
| 6 | CH₃O + O₂ → HCHO + HO₂ | 1.9 × 10⁻¹⁵ | cm³/s |
| 7 | CH₃OOH + OH → CH₃O₂ + H₂O | 3.8 × 10⁻¹² | cm³/s |
| 8 | CH₃OOH + OH → HCHO + OH + H₂O | 1.9 × 10⁻¹² | cm³/s |
| 9 | CH₃OOH + hν → CH₃O + OH | ~5 × 10⁻⁶ | s⁻¹ |
| 10 | HCHO + OH → HCO + H₂O | 8.5 × 10⁻¹² | cm³/s |
| 11 | HCHO + hν → HCO + H | ~3 × 10⁻⁵ | s⁻¹ |
| 12 | HCHO + hν → H₂ + CO | ~5 × 10⁻⁵ | s⁻¹ |
| 13 | HCO + O₂ → CO + HO₂ | 5.2 × 10⁻¹² | cm³/s |
| 14 | H + O₂ + M → HO₂ + M | 5.7 × 10⁻³² [M] | cm⁶/s |
| 15 | HO₂ + NO → OH + NO₂ | 8.1 × 10⁻¹² | cm³/s |
| 16 | NO₂ + hν → NO + O | ~8 × 10⁻³ | s⁻¹ |
| 17 | O + O₂ + M → O₃ + M | 6.0 × 10⁻³⁴ [M] | cm⁶/s |

## Key Features

### Methylperoxy Radical (CH₃O₂) Fate

CH₃O₂ has three loss pathways with different implications:

1. **CH₃O₂ + NO → CH₃O + NO₂** (dominant at high NOx)
   - Produces O₃ (via NO₂ photolysis)
   - Regenerates OH chain

2. **CH₃O₂ + HO₂ → CH₃OOH** (dominant at low NOx)
   - Terminates radical chain
   - CH₃OOH can photolyze to regenerate radicals

3. **CH₃O₂ + CH₃O₂ → products**
   - Minor pathway

### Formaldehyde (HCHO) as Key Intermediate

HCHO is produced from methoxy radical oxidation and has three loss pathways:

1. **HCHO + OH → HCO + H₂O**: Continues oxidation
2. **HCHO + hν → HCO + H** (radical channel): Produces HOx
3. **HCHO + hν → H₂ + CO** (molecular channel): Non-radical

The photolysis branching ratio affects the HOx budget.

### O₃ Yield from CH₄ Oxidation

At high NOx, complete CH₄ oxidation can produce 3-5 O₃ molecules:
- CH₃O₂ + NO → NO₂ → O₃ (1st)
- HO₂ + NO → NO₂ → O₃ (from CH₃O + O₂)
- HO₂ + NO → NO₂ → O₃ (from HCO + O₂)
- Additional O₃ from HCHO photolysis products

At low NOx, O₃ yield is much lower due to peroxy + HO₂ termination.

## Oxidation Chain Diagram

```
CH₄
 ↓ + OH
CH₃
 ↓ + O₂ + M
CH₃O₂ ──────────┬──────────┐
 ↓ + NO         │ + HO₂    │ + CH₃O₂
CH₃O + NO₂      │          │
 ↓ + O₂     CH₃OOH      products
HCHO + HO₂      │
 ↓              ↓ hν
HCO + H₂O   CH₃O + OH
 ↓ + O₂
CO + HO₂
 ↓ + OH
CO₂
```

## Usage

```julia
using TroposphericChemistry
using ModelingToolkit

# Create the algebraic system (steady-state diagnostics)
@named ch4_sys = MethaneOxidation()

# Or create the full ODE system for time evolution
@named ch4_ode = MethaneOxidationODE()
```

## API

```@docs
MethaneOxidation
MethaneOxidationODE
```

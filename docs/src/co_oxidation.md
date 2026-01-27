# [CO Oxidation](@id co_oxidation)

Implementation of CO oxidation chemistry, based on Seinfeld & Pandis
Chapter 6, Section 6.3 (pp. 212-219).

## Background

CO oxidation is the simplest hydrocarbon oxidation mechanism and illustrates
the fundamental HOx cycling that drives tropospheric ozone production:

```
CO + OH → CO₂ + H         (R6.9)
H + O₂ + M → HO₂ + M      (R6.10)
HO₂ + NO → OH + NO₂       (R6.11)
NO₂ + hν → NO + O         (R6.12)
O + O₂ + M → O₃ + M       (R6.13)
─────────────────────────────────
Net: CO + 2O₂ + hν → CO₂ + O₃
```

The key feature is that OH is regenerated through the HO₂ + NO reaction,
allowing **catalytic** ozone production.

## Equations

### Equation 6.14: HO₂ Steady-State (High NOx)

```math
[HO_2] = \frac{k_9[CO][OH]}{k_{11}[NO]}
```

At high NOx, HO₂ is removed primarily by reaction with NO, and the HO₂
concentration is inversely proportional to NO.

### HOx Termination Reactions

The HOx cycle is terminated by:
- **High NOx**: OH + NO₂ + M → HNO₃ + M
- **Low NOx**: HO₂ + HO₂ → H₂O₂ + O₂

The relative importance determines the chain length.

### Chain Length

```math
\text{Chain length} = \frac{k_{HO_2+NO}[HO_2][NO]}{L_{HOx}}
```

where L_HOx is the total HOx loss rate.

Typical values:
- High NOx (urban): 1-5 cycles
- Low NOx (remote): 100+ cycles

### Equations 6.21-6.24: Ozone Production Efficiency

```math
OPE = \frac{\Delta[O_3]}{\Delta[NO_x]} = \frac{P(O_3)}{L(NO_x)}
```

where:
- P(O₃) = gross ozone production rate (from peroxy + NO reactions)
- L(NOx) = NOx loss rate (mainly OH + NO₂ → HNO₃)

| Regime | NOx Level | OPE |
|--------|-----------|-----|
| VOC-limited | High | 1-3 |
| Transition | Medium | 3-10 |
| NOx-limited | Low | 10-30 |

## Rate Constants at 298 K

| Reaction | Rate Constant | Units |
|----------|---------------|-------|
| CO + OH | 2.4 × 10⁻¹³ | cm³ molecule⁻¹ s⁻¹ |
| HO₂ + NO | 8.1 × 10⁻¹² | cm³ molecule⁻¹ s⁻¹ |
| HO₂ + HO₂ | 2.9 × 10⁻¹² | cm³ molecule⁻¹ s⁻¹ |
| OH + NO₂ + M | 1.0 × 10⁻¹¹ | cm³ molecule⁻¹ s⁻¹ |
| HO₂ + O₃ | 2.0 × 10⁻¹⁵ | cm³ molecule⁻¹ s⁻¹ |

## HOx Cycling Diagram

```
                    O₃
                     ↓ (photolysis)
                   O(¹D)
                     ↓ (+ H₂O)
     ┌───────────→ OH ←────────────┐
     │              ↓               │
     │            + CO              │
     │              ↓               │
     │             HO₂              │
     │              ↓               │
     │            + NO              │
     │              ↓               │
     └───────────→ OH + NO₂ ───────┘
                        ↓
                      + hν
                        ↓
                    O₃ produced
```

## Usage

```julia
using TroposphericChemistry
using ModelingToolkit

# Create the CO oxidation system
@named co_sys = COOxidation()

# Or use the OPE diagnostic system
@named ope_sys = OzoneProductionEfficiency()
```

## API

```@docs
COOxidation
OzoneProductionEfficiency
```

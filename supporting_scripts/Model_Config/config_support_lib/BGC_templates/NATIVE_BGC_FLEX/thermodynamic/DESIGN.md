# Thermodynamic-Kinetic Templates for BGC_FLEX

## Overview

These templates implement thermodynamic control of biogeochemical reaction rates using the thermodynamic potential factor (F_T) approach from [Jin & Bethke (2003, 2007)](https://aem.asm.org/content/69/4/2340).

**Key Concept**: Reaction rates are modified by F_T, which ranges from 0 to 1:
- F_T = 1: Reaction proceeds at maximum rate (thermodynamically favorable)
- F_T = 0: Reaction stops (thermodynamic equilibrium reached)

## Available Templates

### Individual Processes (in `templates/` folder)

| Template | Reaction | ΔG° (kJ/mol e-) | F_T at 20°C |
|----------|----------|-----------------|-------------|
| `aerobic_respiration.json` | O2 → H2O | -125 | ≈ 1.0 |
| `nitrification.json` | NH4 → NO2 → NO3 | -275 / -74 | ≈ 1.0 |
| `denitrification_thermodynamic.json` | NO3 → N2 | -119 | ≈ 1.0 |
| `anammox.json` | NH4 + NO2 → N2 | -357 | ≈ 1.0 |
| `manganese_reduction.json` | Mn(IV) → Mn(II) | -95 | ≈ 1.0 |
| `iron_reduction.json` | Fe(III) → Fe(II) | -50 | ≈ 0.64 |
| `sulfate_reduction.json` | SO4 → H2S | -25 | < 0 * |
| `methanogenesis.json` | CO2 → CH4 | -20 | < 0 * |

\* These processes require concentration-dependent ΔG adjustment to proceed

### Integrated Models

| Template | Description |
|----------|-------------|
| `redox_cascade_thermodynamic.json` | O2, NO3, SO4 competition |
| `full_sediment_diagenesis.json` | Complete redox ladder with 13 reactions |

### Documentation

| File | Description |
|------|-------------|
| `COMPARISON_kinetic_vs_thermodynamic.md` | Side-by-side comparison of approaches |

## The Redox Ladder

Electron acceptors are used in order of energy yield:

```
1. O2      ΔG° = -125 kJ/mol e-  (always favorable)
2. NO3     ΔG° = -119 kJ/mol e-  (always favorable)
3. Mn(IV)  ΔG° = -95  kJ/mol e-  (always favorable)
4. Fe(III) ΔG° = -50  kJ/mol e-  (moderately limited)
5. SO4     ΔG° = -25  kJ/mol e-  (strongly limited)
6. CO2     ΔG° = -20  kJ/mol e-  (most limited)
```

F_T automatically enforces this ladder.

---

## Technical Details

## Current Approach (Pure Kinetic)

Currently, BGC_FLEX uses empirical rate expressions:

```
dC/dt = k * C * θ^(T-20)
```

Where `k` is an empirical rate constant that must be calibrated. The reaction proceeds regardless of thermodynamic feasibility.

**Problem**: A reaction can proceed even when it's thermodynamically unfavorable (ΔG > 0).

## Thermodynamic-Kinetic Approach

### The Thermodynamic Potential Factor (F_T)

The key innovation is adding a **thermodynamic potential factor** F_T to the rate law:

```
Rate = k * [X] * F_D * F_A * F_T
```

Where:
- `k` = maximum rate constant
- `[X]` = microbial biomass
- `F_D` = donor (substrate) limitation factor (0-1)
- `F_A` = acceptor (electron acceptor) limitation factor (0-1)
- **`F_T` = thermodynamic potential factor (0-1)**

### F_T Definition

```
F_T = 1 - exp(-(ΔG_available - ΔG_required) / (χ * R * T))
```

Where:
- `ΔG_available` = Gibbs free energy available from the reaction (kJ/mol)
- `ΔG_required` = minimum energy needed for ATP synthesis (~45-50 kJ/mol)
- `χ` = average stoichiometric number (typically 1-4)
- `R` = gas constant (8.314 J/mol·K)
- `T` = temperature (K)

### Key Behaviors of F_T

| Condition | F_T Value | Interpretation |
|-----------|-----------|----------------|
| ΔG_avail >> ΔG_req | → 1 | Thermodynamics not limiting |
| ΔG_avail ≈ ΔG_req | → 0.5 | Partial thermodynamic limitation |
| ΔG_avail < ΔG_req | → 0 | Reaction thermodynamically blocked |

## Implementation in BGC_FLEX

### Option 1: Explicit in Kinetics Expression (Current Parser)

Users can already implement this using the exprtk parser:

```json
{
  "DENITRIFICATION": {
    "CONSUMED": "NO3_N",
    "PRODUCED": "N2",
    "KINETICS": "k_max * X_denit * (NO3_N/(Km_NO3+NO3_N)) * (DOC/(Km_DOC+DOC)) * (1 - exp(-(deltaG_denit - 45)/(2*8.314*T_K/1000)))",
    "PARAMETERS": {
      "k_max": {"VALUE": 10.0, "UNITS": "1/day"},
      "Km_NO3": {"VALUE": 0.5, "UNITS": "mg/L"},
      "Km_DOC": {"VALUE": 2.0, "UNITS": "mg/L"},
      "X_denit": {"VALUE": 1.0, "UNITS": "mg/L", "DESCRIPTION": "Denitrifier biomass"},
      "deltaG_denit": {"VALUE": -110, "UNITS": "kJ/mol", "DESCRIPTION": "ΔG of denitrification"}
    }
  }
}
```

**Limitation**: ΔG is hardcoded, not calculated from actual concentrations.

### Option 2: New THERMODYNAMIC Block (Recommended Enhancement)

Add a new optional block that auto-calculates ΔG from species concentrations:

```json
{
  "DENITRIFICATION": {
    "CONSUMED": "NO3_N",
    "PRODUCED": "N2",

    "THERMODYNAMIC": {
      "REACTION": "NO3- + 1.25 CH2O + H+ -> 0.5 N2 + 1.25 CO2 + 1.75 H2O",
      "DELTA_G0": -110.0,
      "CHI": 2,
      "DELTA_G_MIN": 45.0,
      "SPECIES_MAPPING": {
        "NO3-": "NO3_N",
        "CH2O": "DOC",
        "H+": "pH",
        "N2": "ATMOSPHERIC",
        "CO2": "DIC"
      }
    },

    "KINETICS": "k_max * X * F_D * F_A * F_T",
    "F_D": "NO3_N / (Km_NO3 + NO3_N)",
    "F_A": "DOC / (Km_DOC + DOC)",

    "PARAMETERS": {
      "k_max": {"VALUE": 10.0, "UNITS": "1/day"},
      "Km_NO3": {"VALUE": 0.5, "UNITS": "mg/L"},
      "Km_DOC": {"VALUE": 2.0, "UNITS": "mg/L"},
      "X": {"VALUE": 1.0, "UNITS": "mg/L"}
    }
  }
}
```

The code would automatically:
1. Calculate ΔG from actual concentrations using: `ΔG = ΔG° + RT·ln(Q)`
2. Compute F_T from ΔG
3. Apply F_T to the rate expression

## Example: Denitrification with Thermodynamic Control

### Standard Gibbs Energies (ΔG°) for Common Reactions

| Reaction | ΔG° (kJ/mol e-) |
|----------|-----------------|
| Aerobic respiration (O2) | -125 |
| Denitrification (NO3→N2) | -110 |
| Mn(IV) reduction | -95 |
| Fe(III) reduction | -50 |
| Sulfate reduction | -25 |
| Methanogenesis | -20 |

### Why This Matters

In an anoxic sediment with both NO3 and SO4:
- **Pure kinetic model**: Both denitrification and sulfate reduction proceed simultaneously based on k values
- **Thermodynamic model**: Denitrification dominates because ΔG is more negative; sulfate reduction is suppressed until NO3 is depleted

This correctly captures the **redox ladder** observed in nature.

## Calculating ΔG from Concentrations

For a reaction: `aA + bB → cC + dD`

```
ΔG = ΔG° + R·T·ln(Q)

where Q = ([C]^c · [D]^d) / ([A]^a · [B]^b)
```

### Unit Conversions
- Concentrations should be in mol/L (not mg/L)
- Gases use partial pressures (atm)
- R = 8.314 J/mol·K = 0.008314 kJ/mol·K

## Implementation Roadmap

### Phase 1: Documentation & Templates (Current)
- Create thermodynamic templates showing manual F_T implementation
- Document ΔG° values for common reactions

### Phase 2: Helper Functions
- Add unit conversion utilities (mg/L → mol/L)
- Add ΔG calculation functions
- Add F_T calculation functions

### Phase 3: Native THERMODYNAMIC Block
- Extend JSON parser to recognize THERMODYNAMIC block
- Auto-calculate ΔG from species concentrations at each timestep
- Compute and apply F_T automatically

## Comparison: Kinetic vs Thermodynamic-Kinetic

### Simple Kinetic (Current)
```json
"KINETICS": "k_denitrif * NO3_N * 1.07^(T - 20)"
```
- **Pros**: Simple, fewer parameters
- **Cons**: Rate proceeds regardless of thermodynamic feasibility

### Thermodynamic-Kinetic (Proposed)
```json
"KINETICS": "k_max * (NO3_N/(Km+NO3_N)) * (DOC/(Km_DOC+DOC)) * F_T"
```
With automatic F_T calculation from ΔG.

- **Pros**: More mechanistic, captures redox competition, fewer calibration parameters (ΔG° is known)
- **Cons**: More complex, requires concentration→activity conversion

## References

1. Jin Q, Bethke CM (2003). A new rate law describing microbial respiration. Applied and Environmental Microbiology 69:2340-2348.

2. Jin Q, Bethke CM (2007). The thermodynamics and kinetics of microbial metabolism. American Journal of Science 307:643-677.

3. LaRowe DE, Van Cappellen P (2011). Degradation of natural organic matter: A thermodynamic analysis. Geochimica et Cosmochimica Acta 75:2030-2042.

4. Regnier P et al. (2011). Modeling reactive transport in sediments subject to bioturbation and compaction. Geochimica et Cosmochimica Acta 75:4894-4916.

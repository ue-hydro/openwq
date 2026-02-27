# Comparison: Pure Kinetic vs Thermodynamic-Kinetic

## The Same Reaction - Two Approaches

### Denitrification: NO3 + organic C → N2 + CO2

---

## Approach 1: Pure Kinetic (Current BGC_FLEX)

```json
"kinetics": ["k_denit * NO3_N * 1.07^(T - 20)", "mg-N/L/day"]
```

**Rate equation:**
```
Rate = k × [NO3] × θ^(T-20)
```

**Characteristics:**
- Simple first-order kinetics
- Rate depends only on NO3 concentration and temperature
- Reaction proceeds as long as NO3 > 0
- No consideration of:
  - Electron donor (DOC) availability
  - Product accumulation (N2)
  - Thermodynamic feasibility
  - Microbial biomass

**Parameters to calibrate:** 1 (k_denit)

---

## Approach 2: Michaelis-Menten with Inhibition (Common Enhancement)

```json
"kinetics": ["k_max * (NO3_N/(Km_NO3+NO3_N)) * (DOC/(Km_DOC+DOC)) * (Ki_O2/(Ki_O2+DO))", "mg-N/L/day"]
```

**Rate equation:**
```
Rate = k_max × [NO3]/(Km+[NO3]) × [DOC]/(Km_DOC+[DOC]) × Ki_O2/(Ki_O2+[O2])
```

**Characteristics:**
- Saturating kinetics (Michaelis-Menten)
- Considers both electron donor and acceptor
- O2 inhibition included
- Still no thermodynamic constraint
- Reaction can proceed even at equilibrium

**Parameters to calibrate:** 4 (k_max, Km_NO3, Km_DOC, Ki_O2)

---

## Approach 3: Thermodynamic-Kinetic (Proposed)

```json
"kinetics": ["k_max * X * (NO3_N/(Km_NO3+NO3_N)) * (DOC/(Km_DOC+DOC)) * (1 - exp(-((119-45)/(2*0.008314*(T+273)))))", "mg-N/L/day"]
```

**Rate equation:**
```
Rate = k_max × [X] × F_D × F_A × F_T

where:
  F_D = [NO3]/(Km + [NO3])           # Donor limitation
  F_A = [DOC]/(Km_DOC + [DOC])       # Acceptor limitation
  F_T = 1 - exp(-(ΔG - ΔG_min)/(χRT)) # Thermodynamic factor
```

**Characteristics:**
- Includes microbial biomass [X]
- Saturating kinetics for substrates
- **Thermodynamic constraint (F_T):**
  - F_T → 1 when ΔG << ΔG_min (plenty of energy)
  - F_T → 0 when ΔG → ΔG_min (approaching equilibrium)
- Rate naturally decreases as products accumulate
- No need for arbitrary inhibition constants

**Parameters:** 4 kinetic (k_max, Km_NO3, Km_DOC, X) + 2 thermodynamic (ΔG°, χ)

**But:** ΔG° is known from thermodynamic tables, not calibrated!

---

## Numerical Comparison

### Scenario: Anoxic sediment, T = 20°C

| Parameter | Value |
|-----------|-------|
| NO3 | 2.0 mg-N/L |
| DOC | 5.0 mg-C/L |
| O2 | 0.1 mg/L |
| k_max | 15 /day |
| Km_NO3 | 0.5 mg-N/L |
| Km_DOC | 2.0 mg-C/L |

### Pure Kinetic (k = 0.1 /day)
```
Rate = 0.1 × 2.0 × 1.07^0 = 0.2 mg-N/L/day
```

### Michaelis-Menten
```
F_D = 2.0/(0.5+2.0) = 0.80
F_A = 5.0/(2.0+5.0) = 0.71
Rate = 15 × 0.80 × 0.71 = 8.5 mg-N/L/day
```

### Thermodynamic-Kinetic (assuming X = 1 mg/L)
```
F_D = 0.80
F_A = 0.71
F_T = 1 - exp(-(119-45)/(2×0.008314×293))
    = 1 - exp(-74/4.87)
    = 1 - exp(-15.2)
    = 1 - 2.5×10^-7
    ≈ 1.0

Rate = 15 × 1 × 0.80 × 0.71 × 1.0 = 8.5 mg-N/L/day
```

**Under standard conditions, F_T ≈ 1, so rates are similar!**

---

## When Does F_T Matter?

### Case 1: Low substrate, high products

If DOC is limiting and N2 accumulates (closed system):

```
ΔG = ΔG° + RT×ln(Q)
   = -119 + 2.5×ln([N2][CO2]/[NO3][DOC])
```

As products accumulate, ΔG becomes less negative, F_T decreases.

### Case 2: Competition between electron acceptors

**With O2 present (ΔG° = -125 kJ/mol):**
- Aerobic: F_T_aer = 1 - exp(-(125-45)/(2×0.008314×293)) ≈ 1.0
- Denitrif: F_T_den = suppressed by O2 inhibition anyway

**Without O2, with NO3 present:**
- Denitrif: F_T = 1.0 (proceeds)
- Sulfate red: F_T = 1 - exp(-(25-45)/(2×0.008314×293))
            = 1 - exp(+20/4.87) = 1 - exp(4.1) ≈ **negative → 0**

**Sulfate reduction is thermodynamically blocked when denitrification is active!**

This is the key insight: F_T automatically enforces the redox ladder.

---

## Summary: Is It Worth It?

### When Pure Kinetic is Sufficient:
- Surface waters (well-mixed, oxic)
- Simple screening models
- When calibration data is limited
- Far from thermodynamic equilibrium

### When Thermodynamic-Kinetic Adds Value:
- Sediment diagenesis models
- Redox transition zones
- Anoxic/hypoxic systems
- Groundwater contamination
- When you want mechanistic predictions, not just curve fitting

### Implementation Effort:
- **Current BGC_FLEX can already do it** - just use the full expression
- No code changes needed
- Just need thermodynamic templates as examples
- Future: native ΔG calculation from concentrations would be cleaner

---

## Recommended Thermodynamic Constants

| Reaction | ΔG° (kJ/mol e-) | χ | ΔG_min |
|----------|-----------------|---|--------|
| O2 → H2O | -125 | 2 | 45 |
| NO3 → N2 | -119 | 2 | 45 |
| NO3 → NO2 | -92 | 2 | 45 |
| Mn(IV) → Mn(II) | -95 | 2 | 45 |
| Fe(III) → Fe(II) | -50 | 2 | 45 |
| SO4 → H2S | -25 | 2 | 45 |
| CO2 → CH4 | -20 | 2 | 45 |

**Source:** Jin & Bethke (2007), Thullner et al. (2007)

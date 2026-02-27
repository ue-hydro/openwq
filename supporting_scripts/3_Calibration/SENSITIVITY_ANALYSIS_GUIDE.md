# OpenWQ Parameter Sensitivity Analysis Guide

## Table of Contents
- [Overview](#overview)
- [Parameter Summary](#parameter-summary)
- [Sensitivity Analysis Methods](#sensitivity-analysis-methods)
- [Parameter Categories and Sensitivity Rankings](#parameter-categories-and-sensitivity-rankings)
  - [BGC Parameters (NATIVE_BGC_FLEX)](#bgc-parameters-native_bgc_flex)
  - [PHREEQC Parameters](#phreeqc-parameters)
  - [Sorption Parameters](#sorption-parameters)
  - [Sediment Transport Parameters](#sediment-transport-parameters)
  - [Transport Dissolved Parameters](#transport-dissolved-parameters)
  - [Lateral Exchange Parameters](#lateral-exchange-parameters)
  - [Source/Sink Parameters](#sourcesink-parameters)
- [Recommended Calibration Strategy](#recommended-calibration-strategy)
- [Running Sensitivity Analysis](#running-sensitivity-analysis)

---

## Overview

This document provides a comprehensive sensitivity analysis reference for all calibratable parameters in OpenWQ. Understanding parameter sensitivity is essential for:

1. **Efficient Calibration**: Focus calibration effort on influential parameters
2. **Model Understanding**: Identify which processes drive model behavior
3. **Uncertainty Quantification**: Prioritize parameters for uncertainty analysis
4. **Computational Efficiency**: Reduce the parameter space for optimization

**Total Calibratable Parameters: 106+**

---

## Parameter Summary

| Category | Count | Primary Influence |
|----------|-------|-------------------|
| BGC (NATIVE_BGC_FLEX) | 26 | Nutrient cycling, DO dynamics |
| PHREEQC | 22 | Geochemical speciation, mineral equilibria |
| Sorption | 16 | Nutrient retention, P dynamics |
| Sediment Transport | 14 | Sediment-bound pollutants, erosion |
| Transport Dissolved | 5 | Mixing, dispersion |
| Lateral Exchange | 4 | Compartment interactions |
| Source/Sink | 21 | External loads, land-use contributions |

---

## Sensitivity Analysis Methods

### Morris Screening (Elementary Effects)

**Purpose**: Identify non-influential parameters for elimination from calibration

**Characteristics**:
- Computationally efficient: O(r × (k+1)) model runs (r=10-20 trajectories, k=parameters)
- Provides two metrics:
  - **μ*** (mu-star): Mean absolute elementary effect → Overall influence
  - **σ** (sigma): Standard deviation → Non-linearity/interactions

**Interpretation**:
- High μ*, Low σ: Linear influence → Include in calibration
- High μ*, High σ: Non-linear/interactive → Include, may need careful handling
- Low μ*, Low σ: Non-influential → Fix at default value
- Low μ*, High σ: Interactive but weak → Investigate further

**Recommended for**: Initial screening of 20+ parameters

### Sobol Sensitivity Analysis

**Purpose**: Quantify parameter importance with interaction effects

**Characteristics**:
- More computationally expensive: O(N × (2k+2)) model runs (N=1000-5000)
- Provides:
  - **S1** (First-order index): Direct parameter effect
  - **ST** (Total-order index): Total effect including interactions

**Interpretation**:
- S1 ≈ ST: Parameter acts independently
- ST >> S1: Strong interactions with other parameters
- Sum(S1) < 1: Significant parameter interactions exist

**Recommended for**: Detailed analysis of top 10-15 parameters after Morris screening

---

## Parameter Categories and Sensitivity Rankings

### BGC Parameters (NATIVE_BGC_FLEX)

**Module**: Biogeochemical reaction rates for nitrogen, phosphorus, and oxygen cycling

| Parameter | Units | Range | Default | Expected Sensitivity | Priority |
|-----------|-------|-------|---------|---------------------|----------|
| **k_nitrification** | 1/day | 0.001 - 0.5 | 0.03 | 🔴 **HIGH** | 1 |
| **k_denitrification** | 1/day | 1e-8 - 0.01 | 1e-6 | 🔴 **HIGH** | 1 |
| **k_mineralization_active** | 1/day | 0.01 - 1.0 | 0.1 | 🟡 MEDIUM | 2 |
| **k_mineralization_stable** | 1/day | 0.0001 - 0.01 | 0.001 | 🟢 LOW | 3 |
| **k_immobilization** | 1/day | 0.001 - 0.5 | 0.05 | 🟡 MEDIUM | 2 |
| **k_plant_uptake_N** | 1/day | 0.0 - 0.5 | 0.1 | 🟡 MEDIUM (seasonal) | 2 |
| **k_volatilization** | 1/day | 0.0 - 0.1 | 0.01 | 🟢 LOW | 3 |
| **k_fixation** | mg/L/day | 0.0 - 0.5 | 0.01 | 🟢 LOW | 3 |
| **k_P_mineralization** | 1/day | 0.001 - 0.2 | 0.02 | 🟡 MEDIUM | 2 |
| **k_P_immobilization** | 1/day | 0.001 - 0.1 | 0.01 | 🟢 LOW | 3 |
| **k_P_adsorption** | 1/day | 0.01 - 5.0 | 0.5 | 🔴 **HIGH** | 1 |
| **k_P_desorption** | 1/day | 0.001 - 0.5 | 0.05 | 🟡 MEDIUM | 2 |
| **k_BOD_decay** | 1/day | 0.01 - 1.0 | 0.2 | 🟡 MEDIUM | 2 |
| **k_reaeration** | 1/day | 0.1 - 20.0 | 2.0 | 🔴 **HIGH** (for DO) | 1 |
| **k_SOD** | g/m²/day | 0.1 - 10.0 | 1.5 | 🟡 MEDIUM | 2 |
| **theta_nitrification** | - | 1.02 - 1.10 | 1.08 | 🟡 MEDIUM | 2 |
| **theta_denitrification** | - | 1.02 - 1.12 | 1.07 | 🟡 MEDIUM | 2 |
| **theta_BOD** | - | 1.02 - 1.10 | 1.047 | 🟢 LOW | 3 |
| **K_NH4** | mg/L | 0.01 - 2.0 | 0.5 | 🟡 MEDIUM | 2 |
| **K_NO3** | mg/L | 0.01 - 5.0 | 1.0 | 🟡 MEDIUM | 2 |
| **K_O2_nitrification** | mg/L | 0.1 - 2.0 | 0.5 | 🟡 MEDIUM | 2 |
| **K_O2_denitrification** | mg/L | 0.1 - 1.0 | 0.2 | 🟡 MEDIUM | 2 |

**Key Findings**:
- **Highest Sensitivity**: k_nitrification, k_denitrification, k_P_adsorption, k_reaeration
- These control the primary N and P transformation pathways
- Temperature coefficients (theta) become important for long-term simulations with seasonal temperature variation

**Calibration Notes**:
- Use **log transform** for rate constants (k_*) due to wide ranges
- k_nitrification and k_denitrification often interact - calibrate together
- k_reaeration is site-specific; derive from O'Connor-Dobbins or measured values when possible

---

### PHREEQC Parameters

**Module**: Geochemical speciation and equilibrium (when using PHREEQC BGC module)

| Parameter | Units | Range | Default | Expected Sensitivity | Priority |
|-----------|-------|-------|---------|---------------------|----------|
| **Ca** | mg/L | 5 - 500 | 40 | 🟡 MEDIUM | 2 |
| **Mg** | mg/L | 1 - 200 | 10 | 🟢 LOW | 3 |
| **Na** | mg/L | 1 - 10000 | 20 | 🟢 LOW | 3 |
| **K** | mg/L | 0.1 - 100 | 5 | 🟢 LOW | 3 |
| **Cl** | mg/L | 1 - 20000 | 30 | 🟢 LOW | 3 |
| **S_6_as_SO4** | mg/L | 1 - 2000 | 25 | 🟢 LOW | 3 |
| **C_4_as_HCO3** | mg/L | 10 - 500 | 100 | 🟡 MEDIUM | 2 |
| **N_5_as_NO3** | mg/L-N | 0.01 - 50 | 2 | 🔴 **HIGH** | 1 |
| **N_minus3_as_NH4** | mg/L-N | 0.001 - 10 | 0.1 | 🔴 **HIGH** | 1 |
| **P_as_PO4** | mg/L-P | 0.001 - 5 | 0.05 | 🔴 **HIGH** | 1 |
| **Fe_2** | mg/L | 0.001 - 50 | 0.1 | 🟡 MEDIUM (redox) | 2 |
| **Fe_3** | mg/L | 0.001 - 10 | 0.01 | 🟡 MEDIUM | 2 |
| **Mn_2** | mg/L | 0.001 - 10 | 0.05 | 🟢 LOW | 3 |
| **pCO2_log_atm** | log10(atm) | -4.5 to -1.5 | -3.5 | 🔴 **HIGH** | 1 |
| **pO2_log_atm** | log10(atm) | -6.0 to 0.0 | -0.7 | 🔴 **HIGH** | 1 |
| **Calcite_SI** | - | -3.0 to 1.0 | 0.0 | 🟡 MEDIUM | 2 |
| **Gypsum_SI** | - | -3.0 to 0.0 | -1.0 | 🟢 LOW | 3 |
| **kinetic_rate_denitrif** | mol/L/s | 1e-12 to 1e-6 | 1e-9 | 🔴 **HIGH** | 1 |
| **kinetic_rate_nitrif** | mol/L/s | 1e-10 to 1e-5 | 1e-7 | 🔴 **HIGH** | 1 |
| **kinetic_rate_organic_decay** | 1/s | 1e-9 to 1e-5 | 1e-7 | 🟡 MEDIUM | 2 |
| **Hfo_sites** | mol/mol Fe | 0.01 - 0.5 | 0.2 | 🟡 MEDIUM (P binding) | 2 |
| **Hfo_specific_area** | m²/g | 100 - 800 | 600 | 🟡 MEDIUM | 2 |

**Key Findings**:
- **Initial N/P concentrations** (NO3, NH4, PO4) are highly sensitive - they set the baseline for all reactions
- **Gas partial pressures** (pCO2, pO2) control redox conditions and pH
- **Kinetic rate constants** for nitrification/denitrification dominate nitrogen dynamics
- Major ions (Ca, Mg, Na, K, Cl) have low direct sensitivity but affect ionic strength

**Calibration Notes**:
- Use measured values for major ions when available
- pCO2 is often the most uncertain parameter - varies 10-100× from atmosphere in organic-rich waters
- Iron/manganese concentrations control P solubility through surface complexation

---

### Sorption Parameters

**Module**: Equilibrium and kinetic sorption (FREUNDLICH or LANGMUIR isotherms)

| Parameter | Units | Range | Default | Expected Sensitivity | Priority |
|-----------|-------|-------|---------|---------------------|----------|
| **Kfr_NH4** | L/kg | 0.01 - 100 | 1.0 | 🔴 **HIGH** | 1 |
| **Nfr_NH4** | - | 0.3 - 1.0 | 0.7 | 🟡 MEDIUM | 2 |
| **Kfr_PO4** | L/kg | 0.1 - 500 | 10 | 🔴 **HIGH** | 1 |
| **Nfr_PO4** | - | 0.2 - 0.8 | 0.4 | 🟡 MEDIUM | 2 |
| **Kfr_organic** | L/kg | 0.001 - 1000 | 1.0 | 🟡 MEDIUM | 2 |
| **Nfr_organic** | - | 0.5 - 1.0 | 0.9 | 🟢 LOW | 3 |
| **qmax_NH4** | mg/kg | 10 - 2000 | 200 | 🔴 **HIGH** | 1 |
| **KL_NH4** | L/mg | 0.001 - 1.0 | 0.05 | 🟡 MEDIUM | 2 |
| **qmax_PO4** | mg/kg | 50 - 5000 | 500 | 🔴 **HIGH** | 1 |
| **KL_PO4** | L/mg | 0.01 - 5.0 | 0.5 | 🟡 MEDIUM | 2 |
| **Kadsdes_NH4** | 1/s | 1e-6 - 0.1 | 0.001 | 🟡 MEDIUM | 2 |
| **Kadsdes_PO4** | 1/s | 1e-7 - 0.01 | 0.0001 | 🟡 MEDIUM | 2 |
| **bulk_density** | kg/m³ | 1000 - 2000 | 1500 | 🔴 **HIGH** | 1 |
| **layer_thickness** | m | 0.01 - 0.5 | 0.1 | 🔴 **HIGH** | 1 |

**Key Findings**:
- **Kfr (Freundlich K)** and **qmax (Langmuir capacity)** are the most influential sorption parameters
- These directly control how much nutrient is retained vs. transported
- **bulk_density** and **layer_thickness** multiply the sorption effect - very sensitive
- **Freundlich exponent Nfr** matters most at high concentrations

**Calibration Notes**:
- Measure or estimate from soil type: clay >> silt >> sand
- P sorption parameters are typically more sensitive than NH4 (P is strongly sorbed)
- If sediment P is the focus, prioritize Kfr_PO4 and qmax_PO4

---

### Sediment Transport Parameters

**Module**: HYPE_HBVSED or HYPE_MMF erosion/transport models

| Parameter | Units | Range | Default | Expected Sensitivity | Priority |
|-----------|-------|-------|---------|---------------------|----------|
| **erosion_index** | - | 0.1 - 2.0 | 0.5 | 🔴 **HIGH** | 1 |
| **slope** | m/m | 0.001 - 0.5 | 0.05 | 🔴 **HIGH** | 1 |
| **soil_erosion_factor_land** | - | 0.1 - 2.0 | 0.5 | 🟡 MEDIUM | 2 |
| **soil_erosion_factor_soil** | - | 0.1 - 2.0 | 0.5 | 🟡 MEDIUM | 2 |
| **slope_erosion_exponent** | - | 0.5 - 3.0 | 1.5 | 🔴 **HIGH** | 1 |
| **precip_erosion_exponent** | - | 0.5 - 3.0 | 1.5 | 🔴 **HIGH** | 1 |
| **param_scaling_erosion** | - | 0.1 - 2.0 | 0.5 | 🔴 **HIGH** | 1 |
| **monthly_erosion_jan** | - | -0.5 - 2.0 | 0.0 | 🟡 MEDIUM | 2 |
| **monthly_erosion_jul** | - | -0.5 - 2.0 | 0.5 | 🟡 MEDIUM | 2 |
| **cohesion** (MMF) | kPa | 1 - 30 | 7.5 | 🟡 MEDIUM | 2 |
| **erodibility** (MMF) | g/J | 0.5 - 10 | 2.0 | 🔴 **HIGH** | 1 |
| **sreroexp** (MMF) | - | 0.5 - 2.5 | 1.2 | 🟡 MEDIUM | 2 |
| **cropcover** | fraction | 0 - 1 | 0.5 | 🔴 **HIGH** | 1 |
| **transport_factor_1** | - | 0.05 - 0.5 | 0.2 | 🟡 MEDIUM | 2 |

**Key Findings**:
- **erosion_index**, **slope**, and **exponents** are consistently the most sensitive
- Erosion is highly non-linear with slope and precipitation (exponents > 1)
- **cropcover** matters most during growing season - strong temporal interaction
- param_scaling_erosion acts as a global multiplier - easy to calibrate

**Calibration Notes**:
- Calibrate erosion_index and param_scaling_erosion first as they scale all erosion
- Slope should be measured/derived from DEM, not calibrated
- Monthly factors introduce seasonality - important for snowmelt-driven erosion

---

### Transport Dissolved Parameters

**Module**: OPENWQ_NATIVE_TD_ADVDISP advection-dispersion

| Parameter | Units | Range | Default | Expected Sensitivity | Priority |
|-----------|-------|-------|---------|---------------------|----------|
| **dispersion_x** | m²/s | 0.01 - 100 | 1.0 | 🔴 **HIGH** | 1 |
| **dispersion_y** | m²/s | 0.001 - 10 | 0.1 | 🟡 MEDIUM | 2 |
| **dispersion_z** | m²/s | 1e-5 - 0.1 | 0.001 | 🟢 LOW | 3 |
| **characteristic_length** | m | 10 - 10000 | 100 | 🟡 MEDIUM | 2 |
| **molecular_diffusion** | m²/s | 1e-10 - 1e-8 | 1e-9 | 🟢 LOW | 3 |

**Key Findings**:
- **Longitudinal dispersion (dispersion_x)** dominates in river systems
- dispersion_z only matters in stratified water bodies
- molecular_diffusion is negligible compared to turbulent dispersion

**Calibration Notes**:
- Dispersion is scale-dependent - calibrate per reach if possible
- For rivers: D_x ≈ 0.011 × U × W (Fischer formula)
- For well-mixed systems, high dispersion approximates complete mixing

---

### Lateral Exchange Parameters

**Module**: NATIVE_LE_BOUNDMIX boundary exchange

| Parameter | Units | Range | Default | Expected Sensitivity | Priority |
|-----------|-------|-------|---------|---------------------|----------|
| **K_val_surface_soil** | 1/s | 1e-12 - 1e-6 | 1e-9 | 🟡 MEDIUM | 2 |
| **K_val_soil_groundwater** | 1/s | 1e-14 - 1e-8 | 1e-10 | 🟢 LOW | 3 |
| **K_val_surface_groundwater** | 1/s | 1e-12 - 1e-7 | 1e-9 | 🟡 MEDIUM | 2 |
| **K_val_lake_sediment** | 1/s | 1e-10 - 1e-6 | 1e-8 | 🟡 MEDIUM (lakes) | 2 |

**Key Findings**:
- Sensitivity depends strongly on system configuration
- In river-dominated systems, lateral exchange is less sensitive
- In lake/wetland systems, lake_sediment exchange becomes important

**Calibration Notes**:
- These are "tuning" parameters - no direct measurement possible
- Calibrate last, after primary BGC parameters are established
- Use log transform due to extremely wide ranges

---

### Source/Sink Parameters

**Module**: External loading (CSV, Copernicus LULC, ML model)

| Parameter | Units | Range | Default | Expected Sensitivity | Priority |
|-----------|-------|-------|---------|---------------------|----------|
| **ss_csv_scale_TN** | - | 0.1 - 5.0 | 1.0 | 🔴 **HIGH** | 1 |
| **ss_csv_scale_TP** | - | 0.1 - 5.0 | 1.0 | 🔴 **HIGH** | 1 |
| **ss_csv_scale_all** | - | 0.1 - 5.0 | 1.0 | 🔴 **HIGH** | 1 |
| **cropland_TN_coeff** | kg/ha/yr | 5 - 80 | 25 | 🔴 **HIGH** | 1 |
| **cropland_TP_coeff** | kg/ha/yr | 0.5 - 10 | 2 | 🔴 **HIGH** | 1 |
| **forest_TN_coeff** | kg/ha/yr | 0.5 - 10 | 2 | 🟡 MEDIUM | 2 |
| **forest_TP_coeff** | kg/ha/yr | 0.05 - 1.0 | 0.2 | 🟡 MEDIUM | 2 |
| **grassland_TN_coeff** | kg/ha/yr | 2 - 30 | 10 | 🟡 MEDIUM | 2 |
| **grassland_TP_coeff** | kg/ha/yr | 0.2 - 3.0 | 0.8 | 🟡 MEDIUM | 2 |
| **urban_TN_coeff** | kg/ha/yr | 5 - 50 | 15 | 🔴 **HIGH** | 1 |
| **urban_TP_coeff** | kg/ha/yr | 0.5 - 5.0 | 1.5 | 🔴 **HIGH** | 1 |
| **wetland_TN_coeff** | kg/ha/yr | 0.1 - 5.0 | 1.0 | 🟢 LOW | 3 |
| **wetland_TP_coeff** | kg/ha/yr | 0.01 - 0.5 | 0.1 | 🟢 LOW | 3 |
| **precip_scaling_power** | - | 0.5 - 2.5 | 1.0 | 🔴 **HIGH** | 1 |
| **Q10_biological** | - | 1.5 - 4.0 | 2.0 | 🟡 MEDIUM | 2 |
| **T_reference** | °C | 10 - 25 | 15 | 🟡 MEDIUM | 2 |
| **ml_n_estimators** | count | 50 - 500 | 200 | 🟢 LOW | 3 |
| **ml_max_depth** | count | 3 - 20 | 8 | 🟡 MEDIUM | 2 |
| **ml_min_samples_split** | count | 2 - 50 | 10 | 🟢 LOW | 3 |
| **ml_learning_rate** | - | 0.01 - 0.3 | 0.1 | 🟡 MEDIUM | 2 |

**Key Findings**:
- **Source/sink parameters are typically THE MOST SENSITIVE** parameters in OpenWQ
- External loading often dominates internal cycling, especially in impacted watersheds
- **Cropland and urban export coefficients** strongly affect downstream concentrations
- **precip_scaling_power** controls how loads respond to wet vs. dry years

**Calibration Notes**:
- **ALWAYS calibrate SS parameters first** - they set the mass balance
- Use literature values (export coefficients) as starting points
- Scale factors (ss_csv_scale_*) are ideal for simple calibration
- For Copernicus methods: calibrate dominant land use class first (often cropland)

---

## Recommended Calibration Strategy

### Hierarchical Workflow (for 20+ parameters)

```
┌─────────────────────────────────────────────────────────────┐
│ Step 1: MORRIS SCREENING (~200 model evaluations)          │
│   - Test all ~100 parameters                                │
│   - Identify non-influential (μ* < threshold)               │
│   - FIX low-sensitivity params at defaults                  │
│   - Typically reduces to 30-40 parameters                   │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│ Step 2: SOBOL SENSITIVITY (~1000 evaluations)               │
│   - Analyze remaining parameters                            │
│   - Compute S1 (direct) and ST (total) indices              │
│   - Rank by ST to identify top 10-15 parameters             │
│   - Identify parameter interactions                         │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│ Step 3: DDS CALIBRATION (~300-500 evaluations)              │
│   - Optimize only top 10-15 influential parameters          │
│   - Fix all others at default/screening values              │
│   - Use checkpointing for restart capability                │
└─────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────┐
│ Step 4: VALIDATION                                          │
│   - Test calibrated parameters on independent period        │
│   - Compute performance metrics (RMSE, NSE, KGE)            │
│   - Assess parameter uncertainty                            │
└─────────────────────────────────────────────────────────────┘
```

### Priority-Based Quick Calibration (for limited resources)

If computational resources are limited, calibrate these parameters first:

**Tier 1 (Must Calibrate)** - 5-8 parameters:
1. `ss_csv_scale_all` or dominant land-use export coefficient
2. `k_nitrification` (for N dynamics)
3. `k_denitrification` (for N dynamics)
4. `Kfr_PO4` or `qmax_PO4` (for P dynamics)
5. `dispersion_x` (for transport)
6. `erosion_index` (if sediment transport is used)

**Tier 2 (Important)** - Next 5-8 parameters:
- Temperature coefficients (theta_*)
- Additional BGC rates (k_mineralization, k_P_adsorption)
- Secondary sorption parameters
- Monthly erosion factors

**Tier 3 (Refinement)** - Remaining parameters:
- Lateral exchange coefficients
- Half-saturation constants
- Minor BGC processes

---

## Running Sensitivity Analysis

### Using the Calibration Framework

```python
from calibration_lib.sensitivity.morris_screening import MorrisScreening
from calibration_lib.sensitivity.sobol_analysis import SobolAnalysis
from calibration_lib.parameter_defaults import (
    BGC_PARAMETERS,
    SOURCE_SINK_PARAMETERS,
    create_parameter_entry
)

# Define parameters to analyze
parameters = [
    create_parameter_entry(
        name="k_nitrif",
        category="bgc",
        param_name="k_nitrification",
        file_type="bgc_json",
        path={"reaction": "nitrification", "param": "k"}
    ),
    # ... add more parameters
]

# Morris Screening
morris = MorrisScreening(
    model_runner=model_runner,
    parameters=parameters,
    n_trajectories=20,
    output_dir="./sensitivity_results"
)
morris_results = morris.run()
morris.plot_results()

# Filter to influential parameters
influential_params = [p for p, r in zip(parameters, morris_results)
                      if r.mu_star > 0.1]

# Sobol Analysis on reduced set
sobol = SobolAnalysis(
    model_runner=model_runner,
    parameters=influential_params,
    n_samples=2048,
    output_dir="./sensitivity_results"
)
sobol_results = sobol.run()
sobol.plot_results()
```

### Interpreting Results

```python
# Print sensitivity rankings
print("=" * 60)
print("Parameter Sensitivity Ranking (by Total-Order Index)")
print("=" * 60)

for param, result in sorted(zip(parameters, sobol_results),
                             key=lambda x: x[1].ST, reverse=True):
    interaction = "Yes" if result.ST - result.S1 > 0.1 else "No"
    print(f"{param['name']:25s}  S1={result.S1:.3f}  ST={result.ST:.3f}  "
          f"Interactions: {interaction}")
```

---

## References

1. Chapra, S.C. (2008). *Surface Water-Quality Modeling*. Waveland Press.
2. Thomann, R.V. & Mueller, J.A. (1987). *Principles of Surface Water Quality Modeling*. Harper & Row.
3. Saltelli, A. et al. (2008). *Global Sensitivity Analysis: The Primer*. Wiley.
4. Morris, M.D. (1991). Factorial sampling plans for preliminary computational experiments. *Technometrics* 33(2), 161-174.
5. Sobol, I.M. (2001). Global sensitivity indices for nonlinear mathematical models. *Mathematics and Computers in Simulation* 55, 271-280.
6. HYPE Model Documentation: https://hypeweb.smhi.se/
7. SWAT Theoretical Documentation: https://swat.tamu.edu/

---

*Last updated: 2026*
*OpenWQ Calibration Framework v1.0*

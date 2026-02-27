# OpenWQ BGC_FLEX Reaction Kinetics Improvements

## Overview

This document provides recommendations for improving reaction kinetics in OpenWQ's **NATIVE_BGC_FLEX** module, along with example JSON configurations for common biogeochemical cycles.

**Note**: BGC_FLEX and PHREEQC are mutually exclusive options. This document focuses exclusively on BGC_FLEX.

---

## Current State

OpenWQ's NATIVE_BGC_FLEX module currently supports:

**Dependency variables from host model** (via `dependVar`):
- **Temperature** (K or °C)
- **Shortwave Radiation** (W/m²)
- **Soil Moisture** (fraction) - available from SUMMA

**Internal chemical species** (tracked by OpenWQ):
- All defined species can be used in reaction expressions (e.g., NH4_N, NO3_N, DO, etc.)

These can be used in reaction rate expressions through the exprtk parser.

---

## Proposed Kinetics Improvements

Based on analysis of established water quality models ([WASP](https://www.epa.gov/hydrowq/water-quality-analysis-simulation-program-wasp), [SWAT](https://swat.tamu.edu/media/99192/swat2009-theory.pdf), [HYPE](https://hypeweb.smhi.se/model-water/about-hype-code/), [DNDC](https://www.dndc.sr.unh.edu/papers/DNDC_Scientific_Basis_and_Processes.pdf)), the following improvements would enhance BGC_FLEX accuracy.

### 1. Dissolved Oxygen (DO) as Chemical Species

**Scientific Basis**: Nitrification requires O₂ (aerobic), while denitrification is inhibited by O₂ (anaerobic).

**Implementation**: Define DO as an OpenWQ species and model its dynamics:
- Reaeration from atmosphere
- Consumption by BOD oxidation
- Consumption by nitrification
- Use DO in other reaction expressions via Michaelis-Menten kinetics

### 2. Michaelis-Menten Kinetics

**Scientific Basis**: Most BGC reactions saturate at high substrate concentrations.

**Current**: Linear kinetics (rate ∝ concentration)
**Improved**: `rate = k * S / (Km + S)`

### 3. Oxygen Limitation/Inhibition

```
// Aerobic processes (nitrification, BOD decay)
f_O2 = DO / (K_O2 + DO)

// Anaerobic processes (denitrification)
f_O2_inhib = K_O2_inhib / (K_O2_inhib + DO)
```

### 4. Soil Moisture Effects

**Source**: From host model (SUMMA provides soil moisture fraction)

```
// Optimal moisture for nitrification (~60% WFPS)
f_moist = (SM / SM_opt) * exp(1 - SM / SM_opt)

// Threshold for denitrification (>70% WFPS)
f_moist_anaer = max(0, (SM - SM_thresh) / (1 - SM_thresh))
```

### 5. Temperature Corrections

**Arrhenius/van't Hoff**:
```
f_temp = theta^(T - T_ref)
```

**Q10 formulation**:
```
f_temp = Q10^((T - T_ref) / 10)
```

---

## Example BGC_FLEX JSON Configurations

### Example 1: Simple Nitrogen Cycle (N-only)

Basic nitrogen transformations without oxygen limitation.

```json
{
  "BIOGEOCHEMISTRY_CONFIGURATION": {
    "CHEMICAL_SPECIES": {
      "NH4_N": {
        "DESCRIPTION": "Ammonium nitrogen",
        "UNITS": "mg/L"
      },
      "NO3_N": {
        "DESCRIPTION": "Nitrate nitrogen",
        "UNITS": "mg/L"
      },
      "ORG_N": {
        "DESCRIPTION": "Organic nitrogen",
        "UNITS": "mg/L"
      }
    },

    "CYCLING_FRAMEWORKS": {
      "SIMPLE_N_CYCLE": {
        "MINERALIZATION": {
          "CONSUMED": "ORG_N",
          "PRODUCED": "NH4_N",
          "KINETICS": "k_min * ORG_N * 1.08^(T - 20)",
          "PARAMETERS": {
            "k_min": {"VALUE": 0.05, "UNITS": "1/day", "DESCRIPTION": "Mineralization rate"}
          }
        },
        "NITRIFICATION": {
          "CONSUMED": "NH4_N",
          "PRODUCED": "NO3_N",
          "KINETICS": "k_nitrif * NH4_N * 1.08^(T - 20)",
          "PARAMETERS": {
            "k_nitrif": {"VALUE": 0.1, "UNITS": "1/day", "DESCRIPTION": "Nitrification rate"}
          }
        },
        "DENITRIFICATION": {
          "CONSUMED": "NO3_N",
          "PRODUCED": "N2_loss",
          "KINETICS": "k_denitrif * NO3_N * 1.07^(T - 20)",
          "PARAMETERS": {
            "k_denitrif": {"VALUE": 0.02, "UNITS": "1/day", "DESCRIPTION": "Denitrification rate"}
          }
        }
      }
    },

    "SOIL": {
      "CYCLING_FRAMEWORK": ["SIMPLE_N_CYCLE"]
    }
  }
}
```

---

### Example 2: Nitrogen Cycle with Oxygen Limitation

Full nitrogen cycle with DO dynamics and oxygen-dependent kinetics.

```json
{
  "BIOGEOCHEMISTRY_CONFIGURATION": {
    "CHEMICAL_SPECIES": {
      "DO": {
        "DESCRIPTION": "Dissolved oxygen",
        "UNITS": "mg/L",
        "INITIAL_CONDITION": 8.0
      },
      "BOD": {
        "DESCRIPTION": "Biochemical oxygen demand",
        "UNITS": "mg/L"
      },
      "NH4_N": {
        "DESCRIPTION": "Ammonium nitrogen",
        "UNITS": "mg/L"
      },
      "NO3_N": {
        "DESCRIPTION": "Nitrate nitrogen",
        "UNITS": "mg/L"
      },
      "ORG_N": {
        "DESCRIPTION": "Organic nitrogen",
        "UNITS": "mg/L"
      }
    },

    "CYCLING_FRAMEWORKS": {
      "DO_DYNAMICS": {
        "REAERATION": {
          "CONSUMED": "NONE",
          "PRODUCED": "DO",
          "KINETICS": "k_reaer * (DO_sat - DO) * 1.024^(T - 20)",
          "PARAMETERS": {
            "k_reaer": {"VALUE": 2.0, "UNITS": "1/day", "DESCRIPTION": "Reaeration coefficient"},
            "DO_sat": {"VALUE": 9.0, "UNITS": "mg/L", "DESCRIPTION": "Saturation DO at 20C"}
          },
          "NOTES": "DO_sat should ideally be computed from temperature"
        },
        "BOD_DECAY": {
          "CONSUMED": "DO",
          "PRODUCED": "NONE",
          "KINETICS": "k_BOD * BOD * (DO / (K_O2_BOD + DO)) * 1.047^(T - 20)",
          "PARAMETERS": {
            "k_BOD": {"VALUE": 0.2, "UNITS": "1/day", "DESCRIPTION": "BOD decay rate"},
            "K_O2_BOD": {"VALUE": 0.5, "UNITS": "mg/L", "DESCRIPTION": "Half-sat O2 for BOD"}
          }
        }
      },

      "N_CYCLE_WITH_O2": {
        "MINERALIZATION": {
          "CONSUMED": "ORG_N",
          "PRODUCED": "NH4_N",
          "KINETICS": "k_min * ORG_N * 1.08^(T - 20)",
          "PARAMETERS": {
            "k_min": {"VALUE": 0.05, "UNITS": "1/day"}
          }
        },
        "NITRIFICATION": {
          "CONSUMED": "NH4_N",
          "PRODUCED": "NO3_N",
          "KINETICS": "k_nitrif * NH4_N * (DO / (K_O2_nitrif + DO)) * 1.08^(T - 20)",
          "PARAMETERS": {
            "k_nitrif": {"VALUE": 0.15, "UNITS": "1/day", "DESCRIPTION": "Max nitrification rate"},
            "K_O2_nitrif": {"VALUE": 1.0, "UNITS": "mg/L", "DESCRIPTION": "Half-sat O2 for nitrification"}
          },
          "NOTES": "Oxygen-limited via Michaelis-Menten"
        },
        "NITRIFICATION_O2_DEMAND": {
          "CONSUMED": "DO",
          "PRODUCED": "NONE",
          "KINETICS": "4.57 * k_nitrif * NH4_N * (DO / (K_O2_nitrif + DO)) * 1.08^(T - 20)",
          "PARAMETERS": {
            "k_nitrif": {"VALUE": 0.15, "UNITS": "1/day"},
            "K_O2_nitrif": {"VALUE": 1.0, "UNITS": "mg/L"}
          },
          "NOTES": "4.57 mg O2 consumed per mg NH4-N oxidized"
        },
        "DENITRIFICATION": {
          "CONSUMED": "NO3_N",
          "PRODUCED": "N2_loss",
          "KINETICS": "k_denitrif * NO3_N * (K_O2_inhib / (K_O2_inhib + DO)) * 1.07^(T - 20)",
          "PARAMETERS": {
            "k_denitrif": {"VALUE": 0.1, "UNITS": "1/day", "DESCRIPTION": "Max denitrification rate"},
            "K_O2_inhib": {"VALUE": 0.5, "UNITS": "mg/L", "DESCRIPTION": "O2 inhibition constant"}
          },
          "NOTES": "Inhibited by oxygen - only active when DO is low"
        }
      }
    },

    "RIVER": {
      "CYCLING_FRAMEWORK": ["DO_DYNAMICS", "N_CYCLE_WITH_O2"]
    }
  }
}
```

---

### Example 3: Nitrogen + Phosphorus Cycle

Combined N and P cycling with sorption effects.

```json
{
  "BIOGEOCHEMISTRY_CONFIGURATION": {
    "CHEMICAL_SPECIES": {
      "NH4_N": {"DESCRIPTION": "Ammonium nitrogen", "UNITS": "mg/L"},
      "NO3_N": {"DESCRIPTION": "Nitrate nitrogen", "UNITS": "mg/L"},
      "ORG_N": {"DESCRIPTION": "Organic nitrogen", "UNITS": "mg/L"},
      "PO4_P": {"DESCRIPTION": "Dissolved phosphate", "UNITS": "mg/L"},
      "ORG_P": {"DESCRIPTION": "Organic phosphorus", "UNITS": "mg/L"},
      "PART_P": {"DESCRIPTION": "Particulate/sorbed phosphorus", "UNITS": "mg/L"}
    },

    "CYCLING_FRAMEWORKS": {
      "N_CYCLE": {
        "MINERALIZATION_N": {
          "CONSUMED": "ORG_N",
          "PRODUCED": "NH4_N",
          "KINETICS": "k_min_N * ORG_N * 1.08^(T - 20)",
          "PARAMETERS": {"k_min_N": {"VALUE": 0.05, "UNITS": "1/day"}}
        },
        "NITRIFICATION": {
          "CONSUMED": "NH4_N",
          "PRODUCED": "NO3_N",
          "KINETICS": "k_nitrif * NH4_N * 1.08^(T - 20)",
          "PARAMETERS": {"k_nitrif": {"VALUE": 0.1, "UNITS": "1/day"}}
        },
        "DENITRIFICATION": {
          "CONSUMED": "NO3_N",
          "PRODUCED": "N2_loss",
          "KINETICS": "k_denitrif * NO3_N * 1.07^(T - 20)",
          "PARAMETERS": {"k_denitrif": {"VALUE": 0.02, "UNITS": "1/day"}}
        },
        "PLANT_UPTAKE_NH4": {
          "CONSUMED": "NH4_N",
          "PRODUCED": "NONE",
          "KINETICS": "k_uptake_N * NH4_N * f_season",
          "PARAMETERS": {
            "k_uptake_N": {"VALUE": 0.1, "UNITS": "1/day"},
            "f_season": {"VALUE": 1.0, "DESCRIPTION": "Seasonal factor 0-1"}
          }
        },
        "PLANT_UPTAKE_NO3": {
          "CONSUMED": "NO3_N",
          "PRODUCED": "NONE",
          "KINETICS": "k_uptake_N * NO3_N * f_season",
          "PARAMETERS": {
            "k_uptake_N": {"VALUE": 0.1, "UNITS": "1/day"},
            "f_season": {"VALUE": 1.0}
          }
        }
      },

      "P_CYCLE": {
        "MINERALIZATION_P": {
          "CONSUMED": "ORG_P",
          "PRODUCED": "PO4_P",
          "KINETICS": "k_min_P * ORG_P * 1.06^(T - 20)",
          "PARAMETERS": {"k_min_P": {"VALUE": 0.02, "UNITS": "1/day"}}
        },
        "P_SORPTION": {
          "CONSUMED": "PO4_P",
          "PRODUCED": "PART_P",
          "KINETICS": "k_ads * PO4_P",
          "PARAMETERS": {"k_ads": {"VALUE": 0.5, "UNITS": "1/day", "DESCRIPTION": "Adsorption rate"}}
        },
        "P_DESORPTION": {
          "CONSUMED": "PART_P",
          "PRODUCED": "PO4_P",
          "KINETICS": "k_des * PART_P",
          "PARAMETERS": {"k_des": {"VALUE": 0.05, "UNITS": "1/day", "DESCRIPTION": "Desorption rate"}}
        },
        "PLANT_UPTAKE_P": {
          "CONSUMED": "PO4_P",
          "PRODUCED": "NONE",
          "KINETICS": "k_uptake_P * PO4_P * f_season",
          "PARAMETERS": {
            "k_uptake_P": {"VALUE": 0.05, "UNITS": "1/day"},
            "f_season": {"VALUE": 1.0}
          }
        }
      }
    },

    "SOIL": {
      "CYCLING_FRAMEWORK": ["N_CYCLE", "P_CYCLE"]
    }
  }
}
```

---

### Example 4: Full Streeter-Phelps DO-BOD Model

Classic dissolved oxygen sag curve model for rivers.

```json
{
  "BIOGEOCHEMISTRY_CONFIGURATION": {
    "CHEMICAL_SPECIES": {
      "DO": {
        "DESCRIPTION": "Dissolved oxygen",
        "UNITS": "mg/L",
        "INITIAL_CONDITION": 8.0
      },
      "BOD_fast": {
        "DESCRIPTION": "Fast-decaying BOD (labile)",
        "UNITS": "mg/L"
      },
      "BOD_slow": {
        "DESCRIPTION": "Slow-decaying BOD (refractory)",
        "UNITS": "mg/L"
      }
    },

    "CYCLING_FRAMEWORKS": {
      "STREETER_PHELPS": {
        "REAERATION": {
          "CONSUMED": "NONE",
          "PRODUCED": "DO",
          "KINETICS": "k_reaer * (DO_sat - DO) * 1.024^(T - 20)",
          "PARAMETERS": {
            "k_reaer": {"VALUE": 3.0, "UNITS": "1/day", "DESCRIPTION": "Reaeration rate (depth/velocity dependent)"},
            "DO_sat": {"VALUE": 9.0, "UNITS": "mg/L", "DESCRIPTION": "Saturation DO"}
          }
        },
        "BOD_FAST_DECAY": {
          "CONSUMED": "BOD_fast",
          "PRODUCED": "NONE",
          "KINETICS": "k_BOD_fast * BOD_fast * 1.047^(T - 20)",
          "PARAMETERS": {
            "k_BOD_fast": {"VALUE": 0.3, "UNITS": "1/day", "DESCRIPTION": "Fast BOD decay rate"}
          }
        },
        "BOD_SLOW_DECAY": {
          "CONSUMED": "BOD_slow",
          "PRODUCED": "NONE",
          "KINETICS": "k_BOD_slow * BOD_slow * 1.047^(T - 20)",
          "PARAMETERS": {
            "k_BOD_slow": {"VALUE": 0.05, "UNITS": "1/day", "DESCRIPTION": "Slow BOD decay rate"}
          }
        },
        "O2_CONSUMPTION_FAST": {
          "CONSUMED": "DO",
          "PRODUCED": "NONE",
          "KINETICS": "k_BOD_fast * BOD_fast * 1.047^(T - 20)",
          "PARAMETERS": {
            "k_BOD_fast": {"VALUE": 0.3, "UNITS": "1/day"}
          },
          "NOTES": "1:1 stoichiometry assumed (1 mg BOD consumes 1 mg O2)"
        },
        "O2_CONSUMPTION_SLOW": {
          "CONSUMED": "DO",
          "PRODUCED": "NONE",
          "KINETICS": "k_BOD_slow * BOD_slow * 1.047^(T - 20)",
          "PARAMETERS": {
            "k_BOD_slow": {"VALUE": 0.05, "UNITS": "1/day"}
          }
        },
        "SEDIMENT_O2_DEMAND": {
          "CONSUMED": "DO",
          "PRODUCED": "NONE",
          "KINETICS": "SOD / depth * 1.065^(T - 20)",
          "PARAMETERS": {
            "SOD": {"VALUE": 1.0, "UNITS": "g/m2/day", "DESCRIPTION": "Sediment oxygen demand"},
            "depth": {"VALUE": 2.0, "UNITS": "m", "DESCRIPTION": "Water depth"}
          }
        }
      }
    },

    "RIVER": {
      "CYCLING_FRAMEWORK": ["STREETER_PHELPS"]
    }
  }
}
```

---

### Example 5: Carbon-Nitrogen Coupled Cycle

Organic matter decomposition with C:N stoichiometry effects.

```json
{
  "BIOGEOCHEMISTRY_CONFIGURATION": {
    "CHEMICAL_SPECIES": {
      "DOC": {"DESCRIPTION": "Dissolved organic carbon", "UNITS": "mg/L"},
      "POC": {"DESCRIPTION": "Particulate organic carbon", "UNITS": "mg/L"},
      "DON": {"DESCRIPTION": "Dissolved organic nitrogen", "UNITS": "mg/L"},
      "PON": {"DESCRIPTION": "Particulate organic nitrogen", "UNITS": "mg/L"},
      "NH4_N": {"DESCRIPTION": "Ammonium", "UNITS": "mg/L"},
      "NO3_N": {"DESCRIPTION": "Nitrate", "UNITS": "mg/L"},
      "CO2": {"DESCRIPTION": "Dissolved CO2 (respiration product)", "UNITS": "mg/L"}
    },

    "CYCLING_FRAMEWORKS": {
      "C_N_COUPLED": {
        "POC_HYDROLYSIS": {
          "CONSUMED": "POC",
          "PRODUCED": "DOC",
          "KINETICS": "k_hydrol * POC * 1.05^(T - 20)",
          "PARAMETERS": {
            "k_hydrol": {"VALUE": 0.1, "UNITS": "1/day", "DESCRIPTION": "Hydrolysis rate"}
          }
        },
        "PON_HYDROLYSIS": {
          "CONSUMED": "PON",
          "PRODUCED": "DON",
          "KINETICS": "k_hydrol * PON * 1.05^(T - 20)",
          "PARAMETERS": {"k_hydrol": {"VALUE": 0.1, "UNITS": "1/day"}}
        },
        "DOC_RESPIRATION": {
          "CONSUMED": "DOC",
          "PRODUCED": "CO2",
          "KINETICS": "k_resp * DOC * 1.07^(T - 20)",
          "PARAMETERS": {
            "k_resp": {"VALUE": 0.15, "UNITS": "1/day", "DESCRIPTION": "Respiration rate"}
          }
        },
        "DON_MINERALIZATION": {
          "CONSUMED": "DON",
          "PRODUCED": "NH4_N",
          "KINETICS": "k_min * DON * 1.08^(T - 20)",
          "PARAMETERS": {
            "k_min": {"VALUE": 0.08, "UNITS": "1/day", "DESCRIPTION": "N mineralization rate"}
          },
          "NOTES": "Releases NH4 when C:N < critical ratio (~20)"
        },
        "N_IMMOBILIZATION": {
          "CONSUMED": "NH4_N",
          "PRODUCED": "DON",
          "KINETICS": "k_immob * NH4_N * (DOC / (DOC + K_DOC)) * 1.05^(T - 20)",
          "PARAMETERS": {
            "k_immob": {"VALUE": 0.05, "UNITS": "1/day", "DESCRIPTION": "Immobilization rate"},
            "K_DOC": {"VALUE": 5.0, "UNITS": "mg/L", "DESCRIPTION": "Half-sat for DOC"}
          },
          "NOTES": "Consumes NH4 when C:N > critical ratio"
        },
        "NITRIFICATION": {
          "CONSUMED": "NH4_N",
          "PRODUCED": "NO3_N",
          "KINETICS": "k_nitrif * NH4_N * 1.08^(T - 20)",
          "PARAMETERS": {"k_nitrif": {"VALUE": 0.1, "UNITS": "1/day"}}
        },
        "DENITRIFICATION": {
          "CONSUMED": "NO3_N",
          "PRODUCED": "N2_loss",
          "KINETICS": "k_denitrif * NO3_N * (DOC / (DOC + K_DOC_denit)) * 1.07^(T - 20)",
          "PARAMETERS": {
            "k_denitrif": {"VALUE": 0.08, "UNITS": "1/day"},
            "K_DOC_denit": {"VALUE": 2.0, "UNITS": "mg/L", "DESCRIPTION": "DOC limitation for denitrification"}
          },
          "NOTES": "Denitrification requires organic carbon as electron donor"
        }
      }
    },

    "SOIL": {
      "CYCLING_FRAMEWORK": ["C_N_COUPLED"]
    }
  }
}
```

---

### Example 6: Soil Moisture-Dependent N Cycle (for SUMMA coupling)

Nitrogen transformations controlled by soil moisture from the host model.

```json
{
  "BIOGEOCHEMISTRY_CONFIGURATION": {
    "CHEMICAL_SPECIES": {
      "NH4_N": {"DESCRIPTION": "Ammonium", "UNITS": "mg/L"},
      "NO3_N": {"DESCRIPTION": "Nitrate", "UNITS": "mg/L"},
      "ORG_N": {"DESCRIPTION": "Organic nitrogen", "UNITS": "mg/L"},
      "N2O": {"DESCRIPTION": "Nitrous oxide emissions", "UNITS": "mg/L"}
    },

    "CYCLING_FRAMEWORKS": {
      "MOISTURE_DEPENDENT_N": {
        "MINERALIZATION": {
          "CONSUMED": "ORG_N",
          "PRODUCED": "NH4_N",
          "KINETICS": "k_min * ORG_N * f_moist_aerobic * 1.08^(T - 20)",
          "PARAMETERS": {
            "k_min": {"VALUE": 0.05, "UNITS": "1/day"},
            "f_moist_aerobic": {"VALUE": 1.0, "DESCRIPTION": "Moisture factor for aerobic (from SM dependency)"}
          },
          "NOTES": "f_moist_aerobic = (SM/0.6) * exp(1 - SM/0.6) for optimal at 60% saturation"
        },
        "NITRIFICATION": {
          "CONSUMED": "NH4_N",
          "PRODUCED": "NO3_N",
          "KINETICS": "k_nitrif * NH4_N * f_moist_aerobic * 1.08^(T - 20)",
          "PARAMETERS": {
            "k_nitrif": {"VALUE": 0.15, "UNITS": "1/day"},
            "f_moist_aerobic": {"VALUE": 1.0}
          },
          "NOTES": "Nitrification optimal at ~60% WFPS, inhibited when saturated"
        },
        "NITRIFICATION_N2O": {
          "CONSUMED": "NH4_N",
          "PRODUCED": "N2O",
          "KINETICS": "k_nitrif * NH4_N * f_N2O_nitrif * f_moist_aerobic * 1.08^(T - 20)",
          "PARAMETERS": {
            "k_nitrif": {"VALUE": 0.15, "UNITS": "1/day"},
            "f_N2O_nitrif": {"VALUE": 0.01, "DESCRIPTION": "Fraction of nitrification producing N2O"},
            "f_moist_aerobic": {"VALUE": 1.0}
          }
        },
        "DENITRIFICATION": {
          "CONSUMED": "NO3_N",
          "PRODUCED": "N2_loss",
          "KINETICS": "k_denitrif * NO3_N * f_moist_anaerobic * 1.07^(T - 20)",
          "PARAMETERS": {
            "k_denitrif": {"VALUE": 0.2, "UNITS": "1/day"},
            "f_moist_anaerobic": {"VALUE": 0.0, "DESCRIPTION": "Moisture factor for anaerobic"}
          },
          "NOTES": "f_moist_anaerobic = max(0, (SM - 0.7) / 0.3) for onset at 70% saturation"
        },
        "DENITRIFICATION_N2O": {
          "CONSUMED": "NO3_N",
          "PRODUCED": "N2O",
          "KINETICS": "k_denitrif * NO3_N * f_N2O_denitrif * f_moist_anaerobic * 1.07^(T - 20)",
          "PARAMETERS": {
            "k_denitrif": {"VALUE": 0.2, "UNITS": "1/day"},
            "f_N2O_denitrif": {"VALUE": 0.1, "DESCRIPTION": "Fraction producing N2O vs N2"},
            "f_moist_anaerobic": {"VALUE": 0.0}
          },
          "NOTES": "N2O:N2 ratio increases at intermediate moisture"
        }
      }
    },

    "SOIL": {
      "CYCLING_FRAMEWORK": ["MOISTURE_DEPENDENT_N"]
    }
  }
}
```

---

## Summary: Recommended Implementation Order

| Priority | Feature | Impact | Implementation |
|----------|---------|--------|----------------|
| 🔴 1 | DO as chemical species | Very High | Define DO, model reaeration/consumption |
| 🔴 1 | Michaelis-Menten kinetics | High | Use `S/(Km+S)` in expressions |
| 🔴 1 | Oxygen limitation/inhibition | Very High | Use DO in N-cycle expressions |
| 🟡 2 | Soil moisture from host model | High | Use SM dependency in expressions |
| 🟡 2 | Temperature corrections | Medium | Already available (theta^(T-20)) |
| 🟢 3 | C:N stoichiometry | Medium | Track DOC/DON, use ratios |
| 🟢 3 | Light attenuation | Low | Use SW_RAD dependency |

---

## Parameter Reference

### Typical Rate Constants

| Process | k (1/day) | θ (temp coeff) | Reference |
|---------|-----------|----------------|-----------|
| Nitrification | 0.05 - 0.3 | 1.08 | Chapra 2008 |
| Denitrification | 0.02 - 0.2 | 1.07 | SWAT |
| Mineralization | 0.02 - 0.1 | 1.08 | HYPE |
| BOD decay | 0.1 - 0.5 | 1.047 | Streeter-Phelps |
| Reaeration | 0.5 - 10 | 1.024 | O'Connor-Dobbins |
| P sorption | 0.1 - 1.0 | 1.04 | Literature |

### Half-Saturation Constants

| Substrate | Km (mg/L) | Process |
|-----------|-----------|---------|
| O2 for nitrification | 0.5 - 2.0 | Nitrification |
| O2 for BOD decay | 0.5 - 1.0 | Heterotrophic respiration |
| O2 inhibition (denitrif) | 0.1 - 0.5 | Denitrification |
| NH4 for nitrification | 0.5 - 2.0 | Nitrification |
| NO3 for denitrification | 0.5 - 5.0 | Denitrification |
| DOC for denitrification | 2.0 - 10.0 | Heterotrophic denitrification |

---

## References

1. [WASP 8 Model Documentation - EPA](https://www.epa.gov/hydrowq/water-quality-analysis-simulation-program-wasp)
2. [SWAT Theoretical Documentation](https://swat.tamu.edu/media/99192/swat2009-theory.pdf)
3. [HYPE Model - SMHI](https://hypeweb.smhi.se/model-water/about-hype-code/)
4. [DNDC Scientific Basis](https://www.dndc.sr.unh.edu/papers/DNDC_Scientific_Basis_and_Processes.pdf)
5. Chapra, S.C. (2008). Surface Water-Quality Modeling. Waveland Press.
6. Thomann, R.V. & Mueller, J.A. (1987). Principles of Surface Water Quality Modeling.

---

*Last updated: 2026*
*OpenWQ Development Team*

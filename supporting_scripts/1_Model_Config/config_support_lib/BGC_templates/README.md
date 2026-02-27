# OpenWQ BGC Framework Examples

This directory contains ready-to-use JSON templates for biogeochemical cycling in OpenWQ's **NATIVE_BGC_FLEX** module.

## Directory Structure

```
examples_BGC_framework/
├── README.md                           # This file
│
├── individual_chemicals/               # Single-chemical templates
│   ├── nitrogen_simple.json           # Basic N cycle
│   ├── nitrogen_full.json             # Complete N cycle with all forms
│   ├── phosphorus_simple.json         # Basic P cycle
│   ├── phosphorus_full.json           # Complete P cycle with sorption
│   ├── carbon_organic.json            # DOC/POC cycling
│   ├── dissolved_oxygen.json          # DO dynamics (Streeter-Phelps)
│   ├── sediment_oxygen_demand.json    # SOD modeling
│   ├── silica.json                    # Dissolved silica
│   ├── iron_manganese.json            # Fe/Mn redox cycling
│   ├── sulfur.json                    # S cycling and sulfide
│   ├── pathogens.json                 # E. coli, fecal coliforms
│   ├── pesticides.json                # Generic pesticide decay
│   ├── heavy_metals.json              # Cu, Zn, Pb, Cd, Hg
│   ├── pharmaceuticals.json           # Emerging contaminants
│   └── microplastics.json             # Microplastic transport
│
├── model_frameworks/                   # Complete model-mimicking frameworks
│   ├── SWAT_full_nutrients.json       # SWAT-style N/P/C/sediment
│   ├── QUAL2E_stream.json             # QUAL2E stream water quality
│   ├── QUAL2K_enhanced.json           # QUAL2K with sediment diagenesis
│   ├── HYPE_nutrients.json            # HYPE N/P cycling
│   ├── WASP8_eutrophication.json      # WASP8 eutrophication module
│   ├── CE_QUAL_W2_reservoir.json      # CE-QUAL-W2 reservoir model
│   ├── INCA_nitrogen.json             # INCA-N catchment model
│   ├── SPARROW_regional.json          # SPARROW regional loading
│   └── DNDC_soil_nitrogen.json        # DNDC soil N/C cycling
│
└── phreeqc_templates/                  # PHREEQC geochemical templates
    ├── major_ions_equilibrium.pqi     # Ca, Mg, Na, K, Cl, SO4, HCO3
    ├── carbonate_system.pqi           # CO2-HCO3-CO3 equilibrium
    ├── nitrogen_speciation.pqi        # N species + redox
    ├── phosphorus_minerals.pqi        # P with Fe/Al minerals
    ├── iron_manganese_redox.pqi       # Fe/Mn redox cycling
    ├── sulfur_redox.pqi               # S(-II) to S(VI) cycling
    ├── trace_metals.pqi               # Heavy metal speciation
    ├── acid_mine_drainage.pqi         # AMD chemistry
    └── seawater_mixing.pqi            # Freshwater-seawater mixing
```

## Usage

### For BGC_FLEX (JSON templates)

Copy the desired template to your OpenWQ configuration directory and modify as needed:

```bash
cp individual_chemicals/nitrogen_full.json /path/to/your/openwq/config/bgc_config.json
```

Then reference it in your master OpenWQ configuration.

### For PHREEQC (.pqi templates)

Copy the template and modify the SOLUTION, EQUILIBRIUM_PHASES, and KINETICS blocks:

```bash
cp phreeqc_templates/nitrogen_speciation.pqi /path/to/your/openwq/config/phreeqc_input.pqi
```

## Species Naming Convention

Templates use the **"as element"** convention for nutrient species, following
standard water quality practice:

| Template Name | Meaning | Units |
|---------------|---------|-------|
| `NO3-N` | Nitrate as nitrogen | mg-N/L |
| `NH4-N` | Ammonium as nitrogen | mg-N/L |
| `NO2-N` | Nitrite as nitrogen | mg-N/L |
| `PO4-P` | Phosphate as phosphorus | mg-P/L |
| `SO4-S` | Sulfate as sulfur | mg-S/L |

This convention separates the ion from its element using a hyphen (e.g. `NO3-N`).
The GRQA observation database reports the full ion (e.g. `NO3` in mg-NO3/L).
To compare model output with GRQA observations, a stoichiometric conversion is
applied:

```
model_value (mg-N/L) = GRQA_value (mg-NO3/L) × MW(N) / MW(NO3)
```

| GRQA → Model | Factor | Formula |
|-------------|--------|---------|
| NO3 → NO3-N | 0.2259 | 14.007 / 62.004 |
| NH4 → NH4-N | 0.7764 | 14.007 / 18.039 |
| NO2 → NO2-N | 0.3045 | 14.007 / 46.006 |
| PO4 → PO4-P | 0.3261 | 30.974 / 94.971 |
| SO4 → SO4-S | 0.3338 | 32.065 / 96.06 |

Note: Some templates (e.g. WASP8) use full-ion names (`NO3`, `NH4`, `PO4`)
without the element suffix. These correspond directly to GRQA codes and need
no conversion.

## Parameter Conventions

All templates use consistent parameter naming:
- `k_*` - First-order rate constants (1/day)
- `K_*` - Half-saturation constants (mg/L)
- `theta_*` - Temperature correction factors (dimensionless)
- `f_*` - Fraction or factor (dimensionless)

Temperature correction follows: `rate = rate_20 * theta^(T - 20)`

## References

1. SWAT: https://swat.tamu.edu/media/99192/swat2009-theory.pdf
2. QUAL2E: https://www.epa.gov/ceam/qual2e-windows-versions
3. QUAL2K: https://www.epa.gov/ceam/qual2k
4. HYPE: https://hypeweb.smhi.se/model-water/
5. WASP8: https://www.epa.gov/ceam/water-quality-analysis-simulation-program-wasp
6. CE-QUAL-W2: http://www.ce.pdx.edu/w2/
7. PHREEQC: https://www.usgs.gov/software/phreeqc-version-3

---
*OpenWQ Development Team - 2026*

# OpenWQ PHREEQC Module Setup Guide

## Overview

OpenWQ supports the PHREEQC geochemical engine for advanced water quality simulations. PHREEQC allows modeling of:
- Aqueous speciation and saturation indices
- Equilibrium with mineral phases and gases
- Ion exchange and surface complexation
- Kinetic reactions (dissolution, precipitation, redox)
- Full thermodynamic consistency across chemical species

This guide covers how to set up and configure PHREEQC with OpenWQ.

## Prerequisites

1. **OpenWQ compiled with PHREEQC support** (PhreeqcRM library)
2. **PHREEQC thermodynamic database** (e.g., `phreeqc.dat`)
3. **PHREEQC input file** (`.pqi` format)

## Quick Start

### Step 1: Obtain PHREEQC Database

Download the official PHREEQC database from USGS:

```bash
# Download phreeqc.dat (main thermodynamic database)
curl -L -o phreeqc.dat https://raw.githubusercontent.com/usgs-coupled/phreeqc3/main/database/phreeqc.dat
```

Alternative databases available:
- `llnl.dat` - Lawrence Livermore National Lab database (more minerals)
- `wateq4f.dat` - WATEQ4F database (trace metals)
- `minteq.dat` - MINTEQA2 database
- `pitzer.dat` - High ionic strength solutions

### Step 2: Create PHREEQC Input File (.pqi)

Create a `.pqi` file defining your initial solution chemistry. Example for river water:

```
TITLE River water chemistry - OpenWQ PHREEQC example

# ============================================================
# SOLUTION 1: Define initial water composition
# This solution is replicated to all grid cells by PhreeqcRM.
# ============================================================
SOLUTION 1 River water
    temp      15.0           # Temperature in Celsius
    pH        7.5
    pe        8.0
    redox     pe
    units     mg/kgw         # milligrams per kilogram water
    density   1.0

    # Nitrogen species
    N(5)      2.0 as N       # Nitrate as mg-N/kgw
    N(-3)     0.5 as N       # Ammonium as mg-N/kgw

    # Major ions
    Ca        40.0
    Mg        10.0
    Na        15.0
    K         3.0
    Cl        20.0
    S(6)      15.0           # Sulfate
    Alkalinity 100.0 as HCO3

    # Dissolved oxygen
    O(0)      8.0 as O

# ============================================================
# EQUILIBRIUM_PHASES: Gas-water equilibrium (optional)
# ============================================================
EQUILIBRIUM_PHASES 1
    CO2(g)    -3.5   10.0    # log10 partial pressure, initial moles
    O2(g)     -0.7    1.0

SAVE SOLUTION 1
END
```

### Step 3: Configure OpenWQ Master JSON

In your `openWQ_master.json`, set the biogeochemistry module to PHREEQC:

```json
{
    "MODULES": {
        "BIOGEOCHEMISTRY": {
            "MODULE_NAME": "PHREEQC",
            "MODULE_CONFIG_FILEPATH": "openwq_in/openWQ_MODULE_PHREEQC.json"
        }
    }
}
```

### Step 4: Create PHREEQC Module Configuration

Create `openwq_in/openWQ_MODULE_PHREEQC.json`:

```json
{
    "MODULE_NAME": "PHREEQC",
    "FILEPATH": "openwq_in/phreeqc_river.pqi",
    "DATABASE": "openwq_in/phreeqc.dat",
    "COMPONENT_H2O": true,
    "BGC_GENERAL_MOBILE_SPECIES": ["H", "O", "Charge", "Ca", "Mg", "Na"],
    "TEMPERATURE": {
        "RIVER_NETWORK_REACHES": "Treach_K"
    }
}
```

#### Configuration Options

| Key | Type | Description |
|-----|------|-------------|
| `MODULE_NAME` | string | Must be `"PHREEQC"` |
| `FILEPATH` | string | Path to `.pqi` input file (relative to master JSON) |
| `DATABASE` | string | Path to thermodynamic database (relative to master JSON) |
| `COMPONENT_H2O` | bool | Include H2O as a trackable component |
| `BGC_GENERAL_MOBILE_SPECIES` | array | Species names to transport (case-insensitive) |
| `TEMPERATURE` | object | Map compartments to temperature dependency variables |
| `PRESSURE` | object | (Optional) Map compartments to pressure dependency variables |

### Step 5: Configure Initial Conditions

In `openwq_in/openWQ_config.json`, set initial concentrations for mobile species:

```json
{
    "BIOGEOCHEMISTRY_CONFIGURATION": {
        "RIVER_NETWORK_REACHES": {
            "INITIAL_CONDITIONS": {
                "Data_Format": "JSON",
                "Data": {
                    "H": {
                        "1": ["all", "all", "all", 0, "g/l"]
                    },
                    "O": {
                        "1": ["all", "all", "all", 0, "g/l"]
                    },
                    "Charge": {
                        "1": ["all", "all", "all", 0, "g/l"]
                    },
                    "Ca": {
                        "1": ["all", "all", "all", 0.04, "g/l"]
                    },
                    "Mg": {
                        "1": ["all", "all", "all", 0.01, "g/l"]
                    },
                    "Na": {
                        "1": ["all", "all", "all", 0.015, "g/l"]
                    }
                }
            }
        }
    }
}
```

**Note:** Initial conditions in the config JSON override values in the `.pqi` file.

### Step 6: Configure Output

In `openWQ_master.json`, specify which species to output:

```json
{
    "OPENWQ_OUTPUT": {
        "CHEMICAL_SPECIES": ["H", "O", "Charge", "Ca", "Mg", "Na"],
        "UNITS": "G/L"
    }
}
```

## Unit Conversion (Automatic)

OpenWQ automatically handles concentration unit conversion between OpenWQ and PHREEQC:

| System | Internal Units | Notes |
|--------|----------------|-------|
| **OpenWQ** | mg/L (or g/L) | Mass stored in grams, concentration = mass/volume |
| **PHREEQC** | mol/kgw | Moles per kilogram water (PhreeqcRM default) |

### Conversion Process

OpenWQ retrieves **Gram Formula Weights (GFW)** from PHREEQC for each chemical component and applies automatic conversion:

**Before PHREEQC reactions (OpenWQ → PHREEQC):**
```
conc_mol_kgw = conc_g_L / GFW
```

**After PHREEQC reactions (PHREEQC → OpenWQ):**
```
conc_g_L = conc_mol_kgw × GFW
mass_g = conc_g_L × volume_L
```

### Example GFW Values

| Component | GFW (g/mol) |
|-----------|-------------|
| Ca | 40.08 |
| Mg | 24.312 |
| Na | 22.9898 |
| K | 39.102 |
| Cl | 35.453 |
| S (as SO4) | 32.064 |
| C (as CO3) | 12.0111 |
| N (as NO3) | 14.0067 |

The GFW values are automatically retrieved from PHREEQC during initialization and logged to the console.

### Initial Conditions

When setting initial conditions in `openWQ_config.json`, use **g/L** units:

```json
{
    "Ca": {
        "1": ["all", "all", "all", 0.04, "g/l"]
    }
}
```

This 0.04 g/L of Ca is automatically converted to:
- `0.04 / 40.08 = 0.000998 mol/kgw` for PHREEQC

## Temperature Handling

OpenWQ automatically handles temperature unit conversion:

- **Host model provides Kelvin** (MIZUROUTE's `Treach_K`, SUMMA's state variables)
- **PHREEQC expects Celsius**
- **OpenWQ auto-detects** values > 100 as Kelvin and converts to Celsius

### Temperature Dependency Variables

| Host Model | Compartment | Variable Name | Units |
|------------|-------------|---------------|-------|
| MIZUROUTE | RIVER_NETWORK_REACHES | `Treach_K` | Kelvin |
| MIZUROUTE | LAKES | `Tlake_K` | Kelvin |
| SUMMA | Various | `Tair_K`, `Tsoil_K` | Kelvin |

Example configuration:
```json
{
    "TEMPERATURE": {
        "RIVER_NETWORK_REACHES": "Treach_K",
        "LAKES": "Tlake_K"
    }
}
```

## Complete PHREEQC Species Reference

PHREEQC discovers chemical components automatically from the input file and database.
Species names in `BGC_GENERAL_MOBILE_SPECIES` are **case-insensitive** (both "Ca" and "CA" work).

### Master Species (Elements)

Use these element names in `BGC_GENERAL_MOBILE_SPECIES`:

| Element | Description | Default Species | GFW (g/mol) |
|---------|-------------|-----------------|-------------|
| `H` | Hydrogen | H+ | 1.008 |
| `O` | Oxygen | H2O | 16.0 |
| `Ca` | Calcium | Ca+2 | 40.08 |
| `Mg` | Magnesium | Mg+2 | 24.312 |
| `Na` | Sodium | Na+ | 22.9898 |
| `K` | Potassium | K+ | 39.102 |
| `Fe` | Iron | Fe+2 | 55.847 |
| `Mn` | Manganese | Mn+2 | 54.938 |
| `Al` | Aluminum | Al+3 | 26.9815 |
| `Ba` | Barium | Ba+2 | 137.34 |
| `Sr` | Strontium | Sr+2 | 87.62 |
| `Si` | Silicon (silica) | H4SiO4 | 28.0843 |
| `Cl` | Chloride | Cl- | 35.453 |
| `C` | Carbon (inorganic) | CO3-2 | 12.0111 |
| `S` | Sulfur | SO4-2 | 32.064 |
| `N` | Nitrogen | NO3- | 14.0067 |
| `B` | Boron | H3BO3 | 10.81 |
| `P` | Phosphorus | PO4-3 | 30.9738 |
| `F` | Fluoride | F- | 18.9984 |
| `Li` | Lithium | Li+ | 6.939 |
| `Br` | Bromide | Br- | 79.904 |
| `Zn` | Zinc | Zn+2 | 65.37 |
| `Cd` | Cadmium | Cd+2 | 112.4 |
| `Pb` | Lead | Pb+2 | 207.19 |
| `Cu` | Copper | Cu+2 | 63.546 |
| `Charge` | Electrical charge balance | -- | -- |
| `Alkalinity` | Alkalinity (as CaCO3) | CO3-2 | 50.05 |

### Redox States

For redox-sensitive elements, PHREEQC tracks multiple oxidation states automatically:

| Element | Redox States |
|---------|--------------|
| `H` | H(0) = H2, H(1) = H+ |
| `O` | O(0) = O2, O(-2) = H2O |
| `Fe` | Fe(+2) = Fe+2, Fe(+3) = Fe+3 |
| `Mn` | Mn(+2) = Mn+2, Mn(+3) = Mn+3 |
| `Cu` | Cu(+1) = Cu+, Cu(+2) = Cu+2 |
| `C` | C(+4) = CO3-2/HCO3-, C(-4) = CH4 |
| `S` | S(6) = SO4-2, S(-2) = HS-/H2S |
| `N` | N(+5) = NO3-, N(+3) = NO2-, N(0) = N2, N(-3) = NH4+ |

### Mineral Phases (for EQUILIBRIUM_PHASES)

**Carbonate minerals:**
Calcite, Aragonite, Dolomite, Siderite, Rhodochrosite, Strontianite, Witherite, Cerussite, Smithsonite, Otavite

**Sulfate minerals:**
Gypsum, Anhydrite, Celestite, Barite, Anglesite, Arcanite, Mirabilite, Thenardite, Epsomite, Melanterite

**Silicate minerals:**
SiO2(a), Chalcedony, Quartz, Gibbsite, Kaolinite, Albite, Anorthite, K-feldspar, K-mica, Illite, Ca-Montmorillonite, Chlorite(14A), Talc, Chrysotile, Sepiolite

**Iron/Manganese minerals:**
Hematite, Goethite, Fe(OH)3(a), Pyrite, FeS(ppt), Mackinawite, Hausmannite, Manganite, Pyrochroite

**Other minerals:**
Hydroxyapatite, Vivianite, Fluorite, Halite, Sylvite, Sulfur, Sphalerite, Willemite, Alunite, Jarosite-K

### Gas Phases

Available for equilibrium with aqueous solutions:
- CO2(g), O2(g), H2(g), N2(g), H2S(g), CH4(g), NH3(g), H2O(g)

Redox-uncoupled gas species:
- Hdg (H2), Oxg (O2), Mtg (CH4), Sg (H2S), Ntg (N2)

### Recommended Species by Application

| Application | Recommended Mobile Species |
|-------------|---------------------------|
| Major ion chemistry | Ca, Mg, Na, K, Cl, S, C |
| Nutrient modeling | N, P, C |
| Trace metal transport | Fe, Mn, Zn, Cu, Pb, Cd |
| Complete water quality | Ca, Mg, Na, K, Cl, S, C, N, P, Si, Fe, Mn |

## Using the Python Configuration Generator

The `supporting_scripts/Model_Config/` folder contains a Python helper:

```python
from config_support_lib.Gen_PHREEQCmodule_file import create_phreeqc_module_json

create_phreeqc_module_json(
    bgc_config_filepath="openwq_in/openWQ_MODULE_PHREEQC.json",
    json_header_comment=["# PHREEQC module configuration"],
    phreeqc_input_filepath="/path/to/phreeqc_river.pqi",
    phreeqc_database_filepath="/path/to/phreeqc.dat",
    phreeqc_mobile_species=["Ca", "Mg", "Na", "K", "Cl", "N", "S"],
    phreeqc_component_h2o=True,
    phreeqc_temperature_mapping={
        "RIVER_NETWORK_REACHES": "Treach_K"
    }
)
```

## Calibration Support

The OpenWQ Calibration Framework supports PHREEQC parameters:

```python
# In calibration_config.py
parameters = [
    # Initial solution concentrations
    {
        "name": "initial_NO3",
        "file_type": "phreeqc_pqi",
        "path": {"block": "SOLUTION", "species": "N(5)"},
        "initial": 2.0,
        "bounds": (0.1, 10.0)
    },

    # Equilibrium phase partial pressures
    {
        "name": "pCO2_log",
        "file_type": "phreeqc_pqi",
        "path": {"block": "EQUILIBRIUM_PHASES", "phase": "CO2(g)", "field": "si"},
        "initial": -3.5,
        "bounds": (-4.5, -2.5)
    },

    # Kinetic rate constants (if KINETICS block present)
    {
        "name": "denitrification_rate",
        "file_type": "phreeqc_pqi",
        "path": {"block": "KINETICS", "reaction": "Denitrification", "param": "-parms", "index": 0},
        "initial": 1e-6,
        "bounds": (1e-8, 1e-4),
        "transform": "log"
    }
]
```

## Troubleshooting

### "PHREEQC found 0 components"

**Cause:** Invalid `.pqi` file or database path.

**Solution:**
1. Verify the `.pqi` file contains a valid `SOLUTION` block
2. Verify the database file exists and is complete
3. Check file paths are relative to `openWQ_master.json`

### "Species not found in PHREEQC component list"

**Cause:** Mobile species name doesn't match PHREEQC components.

**Solution:**
1. Run once to see the discovered components in the log:
   ```
   <OpenWQ> PHREEQC: Found 12 components: [1]=H2O [2]=H [3]=O [4]=Charge ...
   ```
2. Use exact component names (case-insensitive)
3. Species like "NO3" should be specified as "N" (the element)

### PHREEQC Convergence Errors

**Cause:** Extreme chemical conditions or incompatible phases.

**Solution:**
1. Convergence warnings are often acceptable for early timesteps
2. Check that initial conditions are chemically reasonable
3. Reduce timestep if errors persist
4. Simplify the chemical system (fewer equilibrium phases)

### Temperature Not Updating

**Cause:** Missing or incorrect temperature mapping.

**Solution:**
1. Verify `TEMPERATURE` is defined in the PHREEQC module JSON
2. Ensure the dependency variable name matches exactly
3. Check that the host model exports the temperature variable

## Example Test Case

A complete PHREEQC test case is available at:
```
test_case_phreeqc/
├── openWQ_master.json
├── openwq_in/
│   ├── openWQ_config.json
│   ├── openWQ_MODULE_PHREEQC.json
│   ├── phreeqc_river.pqi
│   ├── phreeqc.dat
│   └── ...
└── openwq_out/
```

Run with:
```bash
./mizuroute_lakes_cslm_openwq_fast -g 1 1 -m file_manager_openwq.txt
```

## References

- [PHREEQC Version 3 (USGS)](https://www.usgs.gov/software/phreeqc-version-3)
- [PhreeqcRM Documentation](https://water.usgs.gov/water-resources/software/PHREEQC/documentation/phreeqcrm-4.0.0.pdf)
- [OpenWQ Documentation](https://openwq.readthedocs.io)

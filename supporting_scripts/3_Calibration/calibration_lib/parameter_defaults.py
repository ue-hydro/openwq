# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.

# This program, openWQ, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Parameter Defaults and Reasonable Ranges
========================================

This module provides scientifically-grounded min/max bounds for all calibratable
parameters in OpenWQ. These ranges are derived from peer-reviewed literature
and best practices in water quality modeling.

Usage:
    from calibration_lib.parameter_defaults import (
        BGC_PARAMETERS,
        PHREEQC_PARAMETERS,
        SORPTION_PARAMETERS,
        SEDIMENT_PARAMETERS,
        TRANSPORT_PARAMETERS,
        LATERAL_EXCHANGE_PARAMETERS,
        SOURCE_SINK_PARAMETERS
    )

Each parameter dictionary contains:
    - description: What the parameter represents
    - units: Physical units
    - min_val: Minimum reasonable value
    - max_val: Maximum reasonable value
    - default: Typical/default value
    - transform: "linear" or "log" (for optimization)
    - scientific_basis: Literature reference or rationale

References:
    - SWAT theoretical documentation
    - HYPE model documentation
    - PHREEQC database documentation
    - Chapra, S.C. (2008). Surface Water-Quality Modeling
    - Thomann, R.V. & Mueller, J.A. (1987). Principles of Surface Water Quality Modeling
"""

from typing import Dict, Any, List


# =============================================================================
# NATIVE_BGC_FLEX PARAMETERS (Biogeochemical reaction rates)
# =============================================================================

BGC_PARAMETERS: Dict[str, Dict[str, Any]] = {
    # Nitrogen cycle rates
    "k_nitrification": {
        "description": "First-order nitrification rate (NH4 -> NO3)",
        "units": "1/day",
        "min_val": 0.001,
        "max_val": 0.5,
        "default": 0.03,
        "transform": "log",
        "scientific_basis": "Chapra (2008): 0.01-0.2/day typical; can reach 0.5/day in warm water"
    },
    "k_denitrification": {
        "description": "First-order denitrification rate (NO3 -> N2)",
        "units": "1/day",
        "min_val": 1e-8,
        "max_val": 0.01,
        "default": 1e-6,
        "transform": "log",
        "scientific_basis": "Highly variable; 0.001-0.01/day in anoxic sediments (SWAT)"
    },
    "k_mineralization_active": {
        "description": "Active organic N mineralization rate (org_N -> NH4)",
        "units": "1/day",
        "min_val": 0.01,
        "max_val": 1.0,
        "default": 0.1,
        "transform": "log",
        "scientific_basis": "Fast pool: 0.05-0.5/day (Thomann & Mueller 1987)"
    },
    "k_mineralization_stable": {
        "description": "Stable organic N mineralization rate",
        "units": "1/day",
        "min_val": 0.0001,
        "max_val": 0.01,
        "default": 0.001,
        "transform": "log",
        "scientific_basis": "Refractory pool: 10-100x slower than active pool"
    },
    "k_immobilization": {
        "description": "NH4/NO3 immobilization to organic N rate",
        "units": "1/day",
        "min_val": 0.001,
        "max_val": 0.5,
        "default": 0.05,
        "transform": "log",
        "scientific_basis": "Microbial uptake: similar magnitude to mineralization"
    },
    "k_plant_uptake_N": {
        "description": "Plant nitrogen uptake rate",
        "units": "1/day",
        "min_val": 0.0,
        "max_val": 0.5,
        "default": 0.1,
        "transform": "linear",
        "scientific_basis": "Seasonal; 0 in winter, 0.1-0.3/day in growing season"
    },
    "k_volatilization": {
        "description": "Ammonia volatilization rate",
        "units": "1/day",
        "min_val": 0.0,
        "max_val": 0.1,
        "default": 0.01,
        "transform": "linear",
        "scientific_basis": "pH and temp dependent; 0.01-0.05/day typical"
    },
    "k_fixation": {
        "description": "Nitrogen fixation rate (N2 -> org_N)",
        "units": "mg/L/day",
        "min_val": 0.0,
        "max_val": 0.5,
        "default": 0.01,
        "transform": "linear",
        "scientific_basis": "0.01-0.1 mg/L/day in eutrophic waters"
    },

    # Phosphorus cycle rates
    "k_P_mineralization": {
        "description": "Organic P mineralization rate",
        "units": "1/day",
        "min_val": 0.001,
        "max_val": 0.2,
        "default": 0.02,
        "transform": "log",
        "scientific_basis": "Similar to N mineralization but slower (Chapra 2008)"
    },
    "k_P_immobilization": {
        "description": "Phosphate immobilization rate",
        "units": "1/day",
        "min_val": 0.001,
        "max_val": 0.1,
        "default": 0.01,
        "transform": "log",
        "scientific_basis": "Microbial P uptake"
    },
    "k_P_adsorption": {
        "description": "Phosphate adsorption rate",
        "units": "1/day",
        "min_val": 0.01,
        "max_val": 5.0,
        "default": 0.5,
        "transform": "log",
        "scientific_basis": "Fast process; 0.1-2.0/day typical"
    },
    "k_P_desorption": {
        "description": "Phosphate desorption rate",
        "units": "1/day",
        "min_val": 0.001,
        "max_val": 0.5,
        "default": 0.05,
        "transform": "log",
        "scientific_basis": "Slower than adsorption; 0.01-0.2/day typical"
    },

    # Carbon/Oxygen cycle rates
    "k_BOD_decay": {
        "description": "BOD decay rate (first-order)",
        "units": "1/day",
        "min_val": 0.01,
        "max_val": 1.0,
        "default": 0.2,
        "transform": "log",
        "scientific_basis": "0.1-0.5/day for labile BOD (Chapra 2008)"
    },
    "k_reaeration": {
        "description": "Reaeration coefficient",
        "units": "1/day",
        "min_val": 0.1,
        "max_val": 20.0,
        "default": 2.0,
        "transform": "log",
        "scientific_basis": "Depth/velocity dependent; 0.5-10/day typical"
    },
    "k_SOD": {
        "description": "Sediment oxygen demand rate",
        "units": "g/m2/day",
        "min_val": 0.1,
        "max_val": 10.0,
        "default": 1.5,
        "transform": "linear",
        "scientific_basis": "0.5-5.0 g/m2/day for organic-rich sediments"
    },

    # Temperature coefficients
    "theta_nitrification": {
        "description": "Temperature correction factor for nitrification",
        "units": "dimensionless",
        "min_val": 1.02,
        "max_val": 1.10,
        "default": 1.08,
        "transform": "linear",
        "scientific_basis": "Arrhenius factor; 1.04-1.08 typical"
    },
    "theta_denitrification": {
        "description": "Temperature correction factor for denitrification",
        "units": "dimensionless",
        "min_val": 1.02,
        "max_val": 1.12,
        "default": 1.07,
        "transform": "linear",
        "scientific_basis": "Arrhenius factor; 1.04-1.10 typical"
    },
    "theta_BOD": {
        "description": "Temperature correction factor for BOD decay",
        "units": "dimensionless",
        "min_val": 1.02,
        "max_val": 1.10,
        "default": 1.047,
        "transform": "linear",
        "scientific_basis": "Standard value 1.047 (Streeter-Phelps)"
    },

    # Half-saturation constants
    "K_NH4": {
        "description": "Half-saturation for NH4 nitrification",
        "units": "mg/L",
        "min_val": 0.01,
        "max_val": 2.0,
        "default": 0.5,
        "transform": "log",
        "scientific_basis": "0.1-1.0 mg/L typical (Michaelis-Menten)"
    },
    "K_NO3": {
        "description": "Half-saturation for NO3 denitrification",
        "units": "mg/L",
        "min_val": 0.01,
        "max_val": 5.0,
        "default": 1.0,
        "transform": "log",
        "scientific_basis": "0.5-2.0 mg/L typical"
    },
    "K_O2_nitrification": {
        "description": "Half-saturation of O2 for nitrification",
        "units": "mg/L",
        "min_val": 0.1,
        "max_val": 2.0,
        "default": 0.5,
        "transform": "linear",
        "scientific_basis": "0.3-1.0 mg/L; aerobic process"
    },
    "K_O2_denitrification": {
        "description": "Half-saturation of O2 for denitrification inhibition",
        "units": "mg/L",
        "min_val": 0.1,
        "max_val": 1.0,
        "default": 0.2,
        "transform": "linear",
        "scientific_basis": "Inhibited when O2 > ~0.2 mg/L"
    },
}


# =============================================================================
# PHREEQC PARAMETERS (Geochemical model)
# =============================================================================

PHREEQC_PARAMETERS: Dict[str, Dict[str, Any]] = {
    # Initial solution concentrations (mmol/kgw in PHREEQC, mg/L internally)
    "Ca": {
        "description": "Calcium concentration",
        "units": "mg/L",
        "min_val": 5.0,
        "max_val": 500.0,
        "default": 40.0,
        "transform": "log",
        "scientific_basis": "Natural waters: 10-200 mg/L typical"
    },
    "Mg": {
        "description": "Magnesium concentration",
        "units": "mg/L",
        "min_val": 1.0,
        "max_val": 200.0,
        "default": 10.0,
        "transform": "log",
        "scientific_basis": "Natural waters: 5-50 mg/L typical"
    },
    "Na": {
        "description": "Sodium concentration",
        "units": "mg/L",
        "min_val": 1.0,
        "max_val": 10000.0,
        "default": 20.0,
        "transform": "log",
        "scientific_basis": "Freshwater 1-50 mg/L, brackish up to 10000 mg/L"
    },
    "K": {
        "description": "Potassium concentration",
        "units": "mg/L",
        "min_val": 0.1,
        "max_val": 100.0,
        "default": 5.0,
        "transform": "log",
        "scientific_basis": "Natural waters: 1-10 mg/L typical"
    },
    "Cl": {
        "description": "Chloride concentration",
        "units": "mg/L",
        "min_val": 1.0,
        "max_val": 20000.0,
        "default": 30.0,
        "transform": "log",
        "scientific_basis": "Freshwater 5-100 mg/L, saline much higher"
    },
    "S_6_as_SO4": {
        "description": "Sulfate concentration",
        "units": "mg/L",
        "min_val": 1.0,
        "max_val": 2000.0,
        "default": 25.0,
        "transform": "log",
        "scientific_basis": "Natural waters: 5-200 mg/L typical"
    },
    "C_4_as_HCO3": {
        "description": "Bicarbonate alkalinity",
        "units": "mg/L",
        "min_val": 10.0,
        "max_val": 500.0,
        "default": 100.0,
        "transform": "log",
        "scientific_basis": "Natural waters: 50-300 mg/L typical"
    },
    "N_5_as_NO3": {
        "description": "Nitrate concentration",
        "units": "mg/L as N",
        "min_val": 0.01,
        "max_val": 50.0,
        "default": 2.0,
        "transform": "log",
        "scientific_basis": "Pristine <0.5 mg/L, agricultural up to 50 mg/L"
    },
    "N_minus3_as_NH4": {
        "description": "Ammonium concentration",
        "units": "mg/L as N",
        "min_val": 0.001,
        "max_val": 10.0,
        "default": 0.1,
        "transform": "log",
        "scientific_basis": "Surface water typically <0.5 mg/L"
    },
    "P_as_PO4": {
        "description": "Phosphate concentration",
        "units": "mg/L as P",
        "min_val": 0.001,
        "max_val": 5.0,
        "default": 0.05,
        "transform": "log",
        "scientific_basis": "Oligotrophic <0.01, eutrophic 0.1-1.0 mg/L"
    },
    "Fe_2": {
        "description": "Ferrous iron concentration",
        "units": "mg/L",
        "min_val": 0.001,
        "max_val": 50.0,
        "default": 0.1,
        "transform": "log",
        "scientific_basis": "Anoxic waters: 1-50 mg/L; oxic <0.1 mg/L"
    },
    "Fe_3": {
        "description": "Ferric iron concentration",
        "units": "mg/L",
        "min_val": 0.001,
        "max_val": 10.0,
        "default": 0.01,
        "transform": "log",
        "scientific_basis": "Low solubility in oxic neutral waters"
    },
    "Mn_2": {
        "description": "Manganese concentration",
        "units": "mg/L",
        "min_val": 0.001,
        "max_val": 10.0,
        "default": 0.05,
        "transform": "log",
        "scientific_basis": "Anoxic waters: 0.5-5 mg/L; oxic <0.1 mg/L"
    },

    # Equilibrium phase SI values (saturation indices)
    "pCO2_log_atm": {
        "description": "CO2 partial pressure (log10 of atm)",
        "units": "log10(atm)",
        "min_val": -4.5,
        "max_val": -1.5,
        "default": -3.5,
        "transform": "linear",
        "scientific_basis": "Atmospheric: -3.5; organic-rich water up to -1.5"
    },
    "pO2_log_atm": {
        "description": "O2 partial pressure (log10 of atm)",
        "units": "log10(atm)",
        "min_val": -6.0,
        "max_val": 0.0,
        "default": -0.7,
        "transform": "linear",
        "scientific_basis": "Atmospheric: -0.7 (21%); anoxic approaches -6"
    },
    "Calcite_SI": {
        "description": "Calcite saturation index",
        "units": "dimensionless",
        "min_val": -3.0,
        "max_val": 1.0,
        "default": 0.0,
        "transform": "linear",
        "scientific_basis": "Equilibrium: 0; undersaturated <0; supersaturated >0"
    },
    "Gypsum_SI": {
        "description": "Gypsum saturation index",
        "units": "dimensionless",
        "min_val": -3.0,
        "max_val": 0.0,
        "default": -1.0,
        "transform": "linear",
        "scientific_basis": "Usually undersaturated in freshwater"
    },

    # Kinetic rate parameters (if KINETICS block used)
    "kinetic_rate_denitrif": {
        "description": "Kinetic denitrification rate constant",
        "units": "mol/L/s",
        "min_val": 1e-12,
        "max_val": 1e-6,
        "default": 1e-9,
        "transform": "log",
        "scientific_basis": "Microbially-mediated; highly variable"
    },
    "kinetic_rate_nitrif": {
        "description": "Kinetic nitrification rate constant",
        "units": "mol/L/s",
        "min_val": 1e-10,
        "max_val": 1e-5,
        "default": 1e-7,
        "transform": "log",
        "scientific_basis": "Microbially-mediated; highly variable"
    },
    "kinetic_rate_organic_decay": {
        "description": "Organic matter decay rate constant",
        "units": "1/s",
        "min_val": 1e-9,
        "max_val": 1e-5,
        "default": 1e-7,
        "transform": "log",
        "scientific_basis": "First-order decay of particulate/dissolved organic matter"
    },

    # Surface complexation parameters
    "Hfo_sites": {
        "description": "Hydrous ferric oxide surface site density",
        "units": "mol/mol Fe",
        "min_val": 0.01,
        "max_val": 0.5,
        "default": 0.2,
        "transform": "log",
        "scientific_basis": "0.1-0.3 mol/mol typical for HFO"
    },
    "Hfo_specific_area": {
        "description": "HFO specific surface area",
        "units": "m2/g",
        "min_val": 100.0,
        "max_val": 800.0,
        "default": 600.0,
        "transform": "linear",
        "scientific_basis": "Fresh HFO: 600 m2/g; aged: 100-300 m2/g"
    },
}


# =============================================================================
# SORPTION ISOTHERM PARAMETERS
# =============================================================================

SORPTION_PARAMETERS: Dict[str, Dict[str, Any]] = {
    # Freundlich isotherm: q = Kfr * C^Nfr
    "Kfr_NH4": {
        "description": "Freundlich K for ammonium",
        "units": "L/kg",
        "min_val": 0.01,
        "max_val": 100.0,
        "default": 1.0,
        "transform": "log",
        "scientific_basis": "Soil-dependent; sandy 0.1-1, clay 5-50 L/kg"
    },
    "Nfr_NH4": {
        "description": "Freundlich exponent for ammonium",
        "units": "dimensionless",
        "min_val": 0.3,
        "max_val": 1.0,
        "default": 0.7,
        "transform": "linear",
        "scientific_basis": "0.5-1.0 typical; 1.0 = linear"
    },
    "Kfr_PO4": {
        "description": "Freundlich K for phosphate",
        "units": "L/kg",
        "min_val": 0.1,
        "max_val": 500.0,
        "default": 10.0,
        "transform": "log",
        "scientific_basis": "Highly variable; depends on Fe/Al content"
    },
    "Nfr_PO4": {
        "description": "Freundlich exponent for phosphate",
        "units": "dimensionless",
        "min_val": 0.2,
        "max_val": 0.8,
        "default": 0.4,
        "transform": "linear",
        "scientific_basis": "Often 0.3-0.5 for P sorption"
    },
    "Kfr_organic": {
        "description": "Freundlich K for organic compounds",
        "units": "L/kg",
        "min_val": 0.001,
        "max_val": 1000.0,
        "default": 1.0,
        "transform": "log",
        "scientific_basis": "Depends on Kow; pesticides 1-1000 L/kg"
    },
    "Nfr_organic": {
        "description": "Freundlich exponent for organics",
        "units": "dimensionless",
        "min_val": 0.5,
        "max_val": 1.0,
        "default": 0.9,
        "transform": "linear",
        "scientific_basis": "Often near-linear for organics"
    },

    # Langmuir isotherm: q = qmax * KL * C / (1 + KL * C)
    "qmax_NH4": {
        "description": "Maximum sorption capacity for NH4",
        "units": "mg/kg",
        "min_val": 10.0,
        "max_val": 2000.0,
        "default": 200.0,
        "transform": "log",
        "scientific_basis": "Sandy 50-200, clay 500-2000 mg/kg"
    },
    "KL_NH4": {
        "description": "Langmuir affinity for NH4",
        "units": "L/mg",
        "min_val": 0.001,
        "max_val": 1.0,
        "default": 0.05,
        "transform": "log",
        "scientific_basis": "0.01-0.2 L/mg typical"
    },
    "qmax_PO4": {
        "description": "Maximum sorption capacity for PO4",
        "units": "mg/kg",
        "min_val": 50.0,
        "max_val": 5000.0,
        "default": 500.0,
        "transform": "log",
        "scientific_basis": "Depends on Fe/Al oxides content"
    },
    "KL_PO4": {
        "description": "Langmuir affinity for PO4",
        "units": "L/mg",
        "min_val": 0.01,
        "max_val": 5.0,
        "default": 0.5,
        "transform": "log",
        "scientific_basis": "Strong P-mineral affinity: 0.1-2.0 L/mg"
    },

    # Kinetic sorption rates
    "Kadsdes_NH4": {
        "description": "NH4 adsorption/desorption rate",
        "units": "1/s",
        "min_val": 1e-6,
        "max_val": 0.1,
        "default": 0.001,
        "transform": "log",
        "scientific_basis": "Equilibrium in minutes to hours"
    },
    "Kadsdes_PO4": {
        "description": "PO4 adsorption/desorption rate",
        "units": "1/s",
        "min_val": 1e-7,
        "max_val": 0.01,
        "default": 0.0001,
        "transform": "log",
        "scientific_basis": "Slower than NH4; hours to days"
    },

    # Soil properties
    "bulk_density": {
        "description": "Soil bulk density",
        "units": "kg/m3",
        "min_val": 1000.0,
        "max_val": 2000.0,
        "default": 1500.0,
        "transform": "linear",
        "scientific_basis": "Sandy 1600-1800, clay 1200-1400 kg/m3"
    },
    "layer_thickness": {
        "description": "Active sediment layer thickness",
        "units": "m",
        "min_val": 0.01,
        "max_val": 0.5,
        "default": 0.1,
        "transform": "linear",
        "scientific_basis": "Typically 0.05-0.2 m active layer"
    },
}


# =============================================================================
# SEDIMENT TRANSPORT PARAMETERS
# =============================================================================

SEDIMENT_PARAMETERS: Dict[str, Dict[str, Any]] = {
    # HYPE_HBVSED module
    "erosion_index": {
        "description": "Soil erosion index",
        "units": "dimensionless",
        "min_val": 0.1,
        "max_val": 2.0,
        "default": 0.5,
        "transform": "linear",
        "scientific_basis": "Calibration coefficient; HYPE 0.2-1.0 typical"
    },
    "slope": {
        "description": "Average slope",
        "units": "m/m",
        "min_val": 0.001,
        "max_val": 0.5,
        "default": 0.05,
        "transform": "linear",
        "scientific_basis": "Land slope; 0.01-0.30 for most terrain"
    },
    "soil_erosion_factor_land": {
        "description": "Land-use dependent erosion factor",
        "units": "dimensionless",
        "min_val": 0.1,
        "max_val": 2.0,
        "default": 0.5,
        "transform": "linear",
        "scientific_basis": "HYPE: 0.1 (forest) to 1.0 (bare soil)"
    },
    "soil_erosion_factor_soil": {
        "description": "Soil-type dependent erosion factor",
        "units": "dimensionless",
        "min_val": 0.1,
        "max_val": 2.0,
        "default": 0.5,
        "transform": "linear",
        "scientific_basis": "HYPE: 0.1 (sand) to 1.0 (silt)"
    },
    "slope_erosion_exponent": {
        "description": "Slope erosion exponent",
        "units": "dimensionless",
        "min_val": 0.5,
        "max_val": 3.0,
        "default": 1.5,
        "transform": "linear",
        "scientific_basis": "USLE uses ~1.4; HYPE 1.0-2.0"
    },
    "precip_erosion_exponent": {
        "description": "Precipitation erosion exponent",
        "units": "dimensionless",
        "min_val": 0.5,
        "max_val": 3.0,
        "default": 1.5,
        "transform": "linear",
        "scientific_basis": "Erosivity exponent; ~1.5 typical"
    },
    "param_scaling_erosion": {
        "description": "Global erosion scaling parameter",
        "units": "dimensionless",
        "min_val": 0.1,
        "max_val": 2.0,
        "default": 0.5,
        "transform": "linear",
        "scientific_basis": "Calibration multiplier"
    },

    # Monthly erosion factors (seasonal adjustment)
    "monthly_erosion_jan": {
        "description": "January erosion factor",
        "units": "dimensionless",
        "min_val": -0.5,
        "max_val": 2.0,
        "default": 0.0,
        "transform": "linear",
        "scientific_basis": "Seasonal adjustment; frozen ground reduces erosion"
    },
    "monthly_erosion_jul": {
        "description": "July erosion factor",
        "units": "dimensionless",
        "min_val": -0.5,
        "max_val": 2.0,
        "default": 0.5,
        "transform": "linear",
        "scientific_basis": "Peak growing season; crop cover reduces erosion"
    },

    # HYPE_MMF module
    "cohesion": {
        "description": "Soil cohesion",
        "units": "kPa",
        "min_val": 1.0,
        "max_val": 30.0,
        "default": 7.5,
        "transform": "linear",
        "scientific_basis": "Sandy 2-5, loam 5-15, clay 15-30 kPa"
    },
    "erodibility": {
        "description": "Soil erodibility coefficient",
        "units": "g/J",
        "min_val": 0.5,
        "max_val": 10.0,
        "default": 2.0,
        "transform": "linear",
        "scientific_basis": "MMF model; 0.5-5.0 g/J typical"
    },
    "sreroexp": {
        "description": "Surface runoff erosion exponent",
        "units": "dimensionless",
        "min_val": 0.5,
        "max_val": 2.5,
        "default": 1.2,
        "transform": "linear",
        "scientific_basis": "MMF model calibration parameter"
    },
    "cropcover": {
        "description": "Crop canopy cover fraction",
        "units": "fraction",
        "min_val": 0.0,
        "max_val": 1.0,
        "default": 0.5,
        "transform": "linear",
        "scientific_basis": "0 = bare, 1 = full cover"
    },
    "transport_factor_1": {
        "description": "Transport capacity factor 1",
        "units": "dimensionless",
        "min_val": 0.05,
        "max_val": 0.5,
        "default": 0.2,
        "transform": "linear",
        "scientific_basis": "MMF transport coefficient"
    },
    "transport_factor_2": {
        "description": "Transport capacity factor 2",
        "units": "dimensionless",
        "min_val": 0.1,
        "max_val": 1.0,
        "default": 0.5,
        "transform": "linear",
        "scientific_basis": "MMF transport coefficient"
    },
}


# =============================================================================
# TRANSPORT DISSOLVED PARAMETERS
# =============================================================================

TRANSPORT_PARAMETERS: Dict[str, Dict[str, Any]] = {
    "dispersion_x": {
        "description": "Longitudinal dispersion coefficient",
        "units": "m2/s",
        "min_val": 0.01,
        "max_val": 100.0,
        "default": 1.0,
        "transform": "log",
        "scientific_basis": "Rivers: 1-100 m2/s; scale-dependent"
    },
    "dispersion_y": {
        "description": "Transverse dispersion coefficient",
        "units": "m2/s",
        "min_val": 0.001,
        "max_val": 10.0,
        "default": 0.1,
        "transform": "log",
        "scientific_basis": "Typically 10x smaller than longitudinal"
    },
    "dispersion_z": {
        "description": "Vertical dispersion coefficient",
        "units": "m2/s",
        "min_val": 1e-5,
        "max_val": 0.1,
        "default": 0.001,
        "transform": "log",
        "scientific_basis": "Stratified: 1e-5; well-mixed: 0.01-0.1 m2/s"
    },
    "characteristic_length": {
        "description": "Characteristic mixing length",
        "units": "m",
        "min_val": 10.0,
        "max_val": 10000.0,
        "default": 100.0,
        "transform": "log",
        "scientific_basis": "Reach length scale; 10-1000 m typical"
    },
    "molecular_diffusion": {
        "description": "Molecular diffusion coefficient",
        "units": "m2/s",
        "min_val": 1e-10,
        "max_val": 1e-8,
        "default": 1e-9,
        "transform": "log",
        "scientific_basis": "Physical constant ~1e-9 m2/s at 20C"
    },
}


# =============================================================================
# LATERAL EXCHANGE PARAMETERS
# =============================================================================

LATERAL_EXCHANGE_PARAMETERS: Dict[str, Dict[str, Any]] = {
    "K_val_surface_soil": {
        "description": "Exchange coefficient: surface water <-> soil",
        "units": "1/s",
        "min_val": 1e-12,
        "max_val": 1e-6,
        "default": 1e-9,
        "transform": "log",
        "scientific_basis": "Slow diffusive exchange; 1e-10 to 1e-8 1/s"
    },
    "K_val_soil_groundwater": {
        "description": "Exchange coefficient: soil <-> groundwater",
        "units": "1/s",
        "min_val": 1e-14,
        "max_val": 1e-8,
        "default": 1e-10,
        "transform": "log",
        "scientific_basis": "Very slow; 1e-12 to 1e-9 1/s"
    },
    "K_val_surface_groundwater": {
        "description": "Exchange coefficient: surface <-> groundwater",
        "units": "1/s",
        "min_val": 1e-12,
        "max_val": 1e-7,
        "default": 1e-9,
        "transform": "log",
        "scientific_basis": "Hyporheic exchange; highly variable"
    },
    "K_val_lake_sediment": {
        "description": "Exchange coefficient: lake <-> sediment",
        "units": "1/s",
        "min_val": 1e-10,
        "max_val": 1e-6,
        "default": 1e-8,
        "transform": "log",
        "scientific_basis": "Sediment-water exchange"
    },
}


# =============================================================================
# SOURCE/SINK PARAMETERS
# =============================================================================

SOURCE_SINK_PARAMETERS: Dict[str, Dict[str, Any]] = {
    # CSV load scaling
    "ss_csv_scale_TN": {
        "description": "Total nitrogen load scaling factor",
        "units": "dimensionless",
        "min_val": 0.1,
        "max_val": 5.0,
        "default": 1.0,
        "transform": "linear",
        "scientific_basis": "Calibration multiplier for CSV loads"
    },
    "ss_csv_scale_TP": {
        "description": "Total phosphorus load scaling factor",
        "units": "dimensionless",
        "min_val": 0.1,
        "max_val": 5.0,
        "default": 1.0,
        "transform": "linear",
        "scientific_basis": "Calibration multiplier for CSV loads"
    },
    "ss_csv_scale_all": {
        "description": "Global load scaling factor (all species)",
        "units": "dimensionless",
        "min_val": 0.1,
        "max_val": 5.0,
        "default": 1.0,
        "transform": "linear",
        "scientific_basis": "Calibration multiplier for all loads"
    },

    # Copernicus LULC export coefficients (static) - kg/ha/yr
    "cropland_TN_coeff": {
        "description": "Cropland TN export coefficient",
        "units": "kg/ha/yr",
        "min_val": 5.0,
        "max_val": 80.0,
        "default": 25.0,
        "transform": "linear",
        "scientific_basis": "Literature: 10-50 kg/ha/yr; fertilized up to 80"
    },
    "cropland_TP_coeff": {
        "description": "Cropland TP export coefficient",
        "units": "kg/ha/yr",
        "min_val": 0.5,
        "max_val": 10.0,
        "default": 2.0,
        "transform": "linear",
        "scientific_basis": "Literature: 0.5-5 kg/ha/yr; erosion-prone higher"
    },
    "forest_TN_coeff": {
        "description": "Forest TN export coefficient",
        "units": "kg/ha/yr",
        "min_val": 0.5,
        "max_val": 10.0,
        "default": 2.0,
        "transform": "linear",
        "scientific_basis": "Natural background: 1-5 kg/ha/yr"
    },
    "forest_TP_coeff": {
        "description": "Forest TP export coefficient",
        "units": "kg/ha/yr",
        "min_val": 0.05,
        "max_val": 1.0,
        "default": 0.2,
        "transform": "linear",
        "scientific_basis": "Natural background: 0.1-0.5 kg/ha/yr"
    },
    "grassland_TN_coeff": {
        "description": "Grassland TN export coefficient",
        "units": "kg/ha/yr",
        "min_val": 2.0,
        "max_val": 30.0,
        "default": 10.0,
        "transform": "linear",
        "scientific_basis": "5-20 kg/ha/yr; grazed higher"
    },
    "grassland_TP_coeff": {
        "description": "Grassland TP export coefficient",
        "units": "kg/ha/yr",
        "min_val": 0.2,
        "max_val": 3.0,
        "default": 0.8,
        "transform": "linear",
        "scientific_basis": "0.2-2 kg/ha/yr; grazed higher"
    },
    "urban_TN_coeff": {
        "description": "Urban TN export coefficient",
        "units": "kg/ha/yr",
        "min_val": 5.0,
        "max_val": 50.0,
        "default": 15.0,
        "transform": "linear",
        "scientific_basis": "Stormwater: 5-30 kg/ha/yr"
    },
    "urban_TP_coeff": {
        "description": "Urban TP export coefficient",
        "units": "kg/ha/yr",
        "min_val": 0.5,
        "max_val": 5.0,
        "default": 1.5,
        "transform": "linear",
        "scientific_basis": "Stormwater: 0.5-3 kg/ha/yr"
    },
    "wetland_TN_coeff": {
        "description": "Wetland TN export coefficient",
        "units": "kg/ha/yr",
        "min_val": 0.1,
        "max_val": 5.0,
        "default": 1.0,
        "transform": "linear",
        "scientific_basis": "Often net sink; low export 0.5-3 kg/ha/yr"
    },
    "wetland_TP_coeff": {
        "description": "Wetland TP export coefficient",
        "units": "kg/ha/yr",
        "min_val": 0.01,
        "max_val": 0.5,
        "default": 0.1,
        "transform": "linear",
        "scientific_basis": "Often net sink; low export"
    },

    # Climate response parameters (dynamic Copernicus)
    "precip_scaling_power": {
        "description": "Precipitation load scaling exponent",
        "units": "dimensionless",
        "min_val": 0.5,
        "max_val": 2.5,
        "default": 1.0,
        "transform": "linear",
        "scientific_basis": "Linear: 1.0; nonlinear erosion: 1.5-2.0"
    },
    "Q10_biological": {
        "description": "Temperature response Q10 factor",
        "units": "dimensionless",
        "min_val": 1.5,
        "max_val": 4.0,
        "default": 2.0,
        "transform": "linear",
        "scientific_basis": "Standard Q10: 2.0; microbial 1.5-3.0"
    },
    "T_reference": {
        "description": "Reference temperature for Q10",
        "units": "degC",
        "min_val": 10.0,
        "max_val": 25.0,
        "default": 15.0,
        "transform": "linear",
        "scientific_basis": "Regional baseline temperature"
    },

    # ML model hyperparameters
    "ml_n_estimators": {
        "description": "Number of trees in Random Forest",
        "units": "count",
        "min_val": 50,
        "max_val": 500,
        "default": 200,
        "transform": "linear",
        "dtype": "int",
        "scientific_basis": "Typically 100-300; diminishing returns >500"
    },
    "ml_max_depth": {
        "description": "Maximum tree depth",
        "units": "count",
        "min_val": 3,
        "max_val": 20,
        "default": 8,
        "transform": "linear",
        "dtype": "int",
        "scientific_basis": "Prevent overfitting; 5-15 typical"
    },
    "ml_min_samples_split": {
        "description": "Minimum samples to split internal node",
        "units": "count",
        "min_val": 2,
        "max_val": 50,
        "default": 10,
        "transform": "linear",
        "dtype": "int",
        "scientific_basis": "Regularization parameter"
    },
    "ml_learning_rate": {
        "description": "Gradient boosting learning rate",
        "units": "dimensionless",
        "min_val": 0.01,
        "max_val": 0.3,
        "default": 0.1,
        "transform": "log",
        "scientific_basis": "Lower = more trees needed but better generalization"
    },
}


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_parameter_info(category: str, param_name: str) -> Dict[str, Any]:
    """
    Get information about a specific parameter.

    Parameters
    ----------
    category : str
        One of: "bgc", "phreeqc", "sorption", "sediment",
                "transport", "lateral_exchange", "source_sink"
    param_name : str
        Parameter name

    Returns
    -------
    Dict
        Parameter information including bounds, default, etc.
    """
    category_map = {
        "bgc": BGC_PARAMETERS,
        "phreeqc": PHREEQC_PARAMETERS,
        "sorption": SORPTION_PARAMETERS,
        "sediment": SEDIMENT_PARAMETERS,
        "transport": TRANSPORT_PARAMETERS,
        "lateral_exchange": LATERAL_EXCHANGE_PARAMETERS,
        "source_sink": SOURCE_SINK_PARAMETERS,
    }

    if category not in category_map:
        raise ValueError(f"Unknown category: {category}")

    params = category_map[category]
    if param_name not in params:
        raise ValueError(f"Unknown parameter: {param_name} in category {category}")

    return params[param_name]


def get_bounds(category: str, param_name: str) -> tuple:
    """Get (min, max) bounds for a parameter."""
    info = get_parameter_info(category, param_name)
    return (info["min_val"], info["max_val"])


def get_default(category: str, param_name: str) -> float:
    """Get default value for a parameter."""
    return get_parameter_info(category, param_name)["default"]


def get_transform(category: str, param_name: str) -> str:
    """Get optimization transform for a parameter."""
    return get_parameter_info(category, param_name).get("transform", "linear")


def list_all_parameters() -> Dict[str, List[str]]:
    """List all available parameters by category."""
    return {
        "bgc": list(BGC_PARAMETERS.keys()),
        "phreeqc": list(PHREEQC_PARAMETERS.keys()),
        "sorption": list(SORPTION_PARAMETERS.keys()),
        "sediment": list(SEDIMENT_PARAMETERS.keys()),
        "transport": list(TRANSPORT_PARAMETERS.keys()),
        "lateral_exchange": list(LATERAL_EXCHANGE_PARAMETERS.keys()),
        "source_sink": list(SOURCE_SINK_PARAMETERS.keys()),
    }


def create_parameter_entry(
    name: str,
    category: str,
    param_name: str,
    file_type: str,
    path: Any,
    initial: float = None,
    bounds: tuple = None,
    transform: str = None
) -> Dict:
    """
    Create a calibration parameter entry with defaults from this module.

    Parameters
    ----------
    name : str
        Unique name for calibration
    category : str
        Parameter category (bgc, phreeqc, etc.)
    param_name : str
        Parameter name in this module
    file_type : str
        OpenWQ file type (bgc_json, phreeqc_pqi, etc.)
    path : Any
        Path specification for parameter_handler
    initial : float, optional
        Override default initial value
    bounds : tuple, optional
        Override (min, max) bounds
    transform : str, optional
        Override transform (linear/log)

    Returns
    -------
    Dict
        Complete parameter entry for calibration_parameters list
    """
    info = get_parameter_info(category, param_name)

    return {
        "name": name,
        "file_type": file_type,
        "path": path,
        "initial": initial if initial is not None else info["default"],
        "bounds": bounds if bounds is not None else (info["min_val"], info["max_val"]),
        "transform": transform if transform is not None else info.get("transform", "linear"),
    }


def validate_parameter_bounds(param: Dict) -> List[str]:
    """
    Validate that a parameter's bounds are reasonable.

    Returns list of warning messages (empty if valid).
    """
    warnings = []

    name = param.get("name", "unknown")
    bounds = param.get("bounds", (0, 1))
    initial = param.get("initial", (bounds[0] + bounds[1]) / 2)

    # Check initial is within bounds
    if not (bounds[0] <= initial <= bounds[1]):
        warnings.append(f"{name}: initial value {initial} outside bounds {bounds}")

    # Check bounds are valid
    if bounds[0] >= bounds[1]:
        warnings.append(f"{name}: invalid bounds {bounds} (min >= max)")

    # Check for very wide bounds (may cause optimization issues)
    if bounds[1] / bounds[0] > 1e6 and param.get("transform") != "log":
        warnings.append(f"{name}: very wide bounds {bounds} - consider log transform")

    return warnings


# Print summary when run directly
if __name__ == "__main__":
    print("=" * 70)
    print("OpenWQ Calibration Parameter Reference")
    print("=" * 70)

    all_params = list_all_parameters()
    total = 0

    for category, params in all_params.items():
        print(f"\n{category.upper()}: {len(params)} parameters")
        total += len(params)
        for p in params[:3]:
            info = get_parameter_info(category, p)
            print(f"  - {p}: [{info['min_val']}, {info['max_val']}] {info['units']}")
        if len(params) > 3:
            print(f"  ... and {len(params) - 3} more")

    print(f"\nTotal: {total} calibratable parameters")

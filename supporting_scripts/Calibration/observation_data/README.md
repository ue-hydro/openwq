# Observation Data Format

This directory contains observation data files for OpenWQ calibration.

## Required Format

Observation data must be in CSV format with the following columns:

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| datetime | datetime | Yes | Observation timestamp (YYYY-MM-DD HH:MM:SS) |
| reach_id | integer | Yes | River network reach ID matching model output |
| species | string | Yes | Chemical species name (must match model species) |
| value | float | Yes | Measured concentration value |
| units | string | Yes | Concentration units (e.g., mg/l, ug/l) |
| source | string | No | Data source identifier |
| uncertainty | float | No | Measurement uncertainty (same units as value) |
| quality_flag | string | No | Quality flag (GOOD, SUSPECT, BAD) |

## Example

```csv
datetime,reach_id,species,value,units,source,uncertainty,quality_flag
2010-06-01 00:00:00,1200014181,NO3-N,2.50,mg/l,USGS_station_A,0.25,GOOD
2010-06-01 00:00:00,1200014181,NH4-N,0.15,mg/l,USGS_station_A,0.02,GOOD
```

## Notes

1. **Species names** must match exactly with species defined in your model configuration
2. **Reach IDs** must match the reach identifiers in your model output (HDF5 files)
3. **Units** should be consistent with model output units (typically mg/l or MG/L)
4. **datetime** format must be parseable by pandas (ISO format recommended)

## Supported Species (examples)

Common species names used in OpenWQ:
- NO3-N (Nitrate-nitrogen)
- NH4-N (Ammonium-nitrogen)
- N_ORG_active (Active organic nitrogen)
- N_ORG_stable (Stable organic nitrogen)
- PO4-P (Phosphate-phosphorus)
- TSS (Total suspended solids)
- DOC (Dissolved organic carbon)

## Quality Control

The calibration framework will:
1. Filter out observations with `quality_flag = 'BAD'` by default
2. Use uncertainty values for weighted objective functions (if provided)
3. Match observations to model output based on reach_id, species, and nearest timestamp

# Observation Data for OpenWQ Calibration

This directory contains tools and data for preparing observation data for OpenWQ calibration.

## Data Sources

### Option 1: Manual CSV File (`observation_data_source = "csv"`)

Prepare observation data manually in CSV format (see format below).

### Option 2: GRQA Database Extraction (`observation_data_source = "grqa"`)

Use the included scripts to automatically extract data from the **Global River Water Quality Archive (GRQA)** database (https://zenodo.org/records/15335450).

**Scripts:**
- `grqa_extract_stations.py` - Downloads/loads and extracts GRQA data for your study area
- `prepare_calibration_observations.py` - Converts GRQA data to calibration format

**Data Source Options:**
- **Download from Zenodo**: Set `local_data_path = None` (default) - data will be downloaded automatically
- **Use local data**: Set `local_data_path` to your GRQA data folder if already downloaded

**Usage:**

```python
# In your calibration configuration:
observation_data_source = "grqa"

grqa_config = {
    # Use local data (if already downloaded) or None to download
    "local_data_path": "/data/GRQA",  # or None to download from Zenodo

    # River network shapefile (bounding box derived from shapefile extent)
    "river_network_shapefile": "/path/to/river_network.shp",
    "reach_id_column": "seg_id",
    "start_date": "2000-01-01",
    "end_date": "2020-12-31",
    "max_station_distance_m": 500,
    "species_mapping": {
        "NO3": "NO3-N",
        "NH4": "NH4-N",
        "TN": "TN",
        "PO4": "PO4-P",
    },
}
```

**Expected local data structure:**
```
local_data_path/
├── GRQA_data_v1.3/ or GRQA_data_v1.4/
│   ├── GRQA_data_v1.3_NO3.csv
│   ├── GRQA_data_v1.3_NH4.csv
│   └── ... (other parameter files)
└── GRQA_meta_v1.3.csv (optional)
```

```bash
# Extract GRQA data only
python my_calibration.py --extract-grqa-only

# Or run calibration (extraction is automatic)
python my_calibration.py
```

---

## Required CSV Format

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

## GRQA Species Mapping

When using GRQA extraction, map GRQA parameter names to your model species names:

| GRQA Parameter | Description | Typical Model Species |
|----------------|-------------|----------------------|
| NO3 | Nitrate | NO3-N |
| NH4 | Ammonium | NH4-N |
| TN | Total Nitrogen | TN |
| NO2 | Nitrite | NO2-N |
| DON | Dissolved Organic Nitrogen | DON |
| PO4 | Orthophosphate | PO4-P |
| TP | Total Phosphorus | TP |
| DP | Dissolved Phosphorus | DP |
| DO | Dissolved Oxygen | DO |
| BOD | Biochemical Oxygen Demand | BOD |
| COD | Chemical Oxygen Demand | COD |
| DOC | Dissolved Organic Carbon | DOC |
| TOC | Total Organic Carbon | TOC |
| TSS | Total Suspended Solids | TSS |
| Chl-a | Chlorophyll-a | Chla |
| EC | Electrical Conductivity | EC |
| pH | pH | pH |

## Notes

1. **Species names** must match exactly with species defined in your model configuration
2. **Reach IDs** must match the reach identifiers in your model output (HDF5 files)
3. **Units** should be consistent with model output units (typically mg/l or MG/L)
4. **datetime** format must be parseable by pandas (ISO format recommended)

## Common Species Names in OpenWQ

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

## Dependencies

### For GRQA Extraction

```bash
pip install geopandas requests shapely pandas numpy
```

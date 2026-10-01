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
# !/usr/bin/env python3
"""
Gen_Observation_Sources.py — multi-source water-quality observation orchestrator.

The user selects one or more observation datasets in the model config. This
module:
  * documents every available source (OBSERVATION_SOURCES registry),
  * flags REDUNDANT selections (e.g. WQP is already inside GRQA) and, in an
    interactive terminal, asks whether to remove the duplication,
  * extracts each selected source, clipped to the basin/river search area, into
    ONE harmonized schema, then merges (optionally de-duplicating rows).

Harmonized schema (every adapter emits exactly these columns):
    station_id, lat, lon, parameter, year, month, day, minute, value, units, source
`parameter` uses the model/BGC species names; `source` records provenance.

Adapters are best-effort for the large bulk downloads — their live column names
come from each dataset's published data descriptor and may need a one-line tweak
in the registry on first download. GRQA, WQP, the redundancy/dedup logic and the
harmonized schema are the validated core.
"""

import io
import os
import sys
import json
import shutil
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

HARMONIZED_COLUMNS = [
    "station_id", "lat", "lon", "parameter",
    "year", "month", "day", "minute", "value", "units", "source",
]

# Dedup key: two observations are "the same" if they agree on all of these.
DEDUP_KEY = ["lat", "lon", "parameter", "year", "month", "day", "value"]


# ─────────────────────────────────────────────────────────────────────────────
#  Source registry
#  access:   "auto"  = downloads itself (open API or static file)
#            "manual"= online but gated (Terms/registration) -> needs a local
#                      path the user downloaded once
#  overlaps: other source keys whose observations this source (partly) contains
#            or is contained by -> used to flag redundancy
#  fmt:      "long"  = one row per observation (param + value columns)
#            "wide"  = one column per parameter (melted on read)
# ─────────────────────────────────────────────────────────────────────────────
OBSERVATION_SOURCES: Dict[str, dict] = {
    "grqa": {
        "name": "GRQA — Global River Water Quality Archive",
        "description": ("Aggregated global archive: 17M+ obs, 43 parameters "
                        "(nutrients, carbon, oxygen, sediment, ions). Combines "
                        "WQP + GEMStat + GLORICH + Waterbase + CESI; frozen ~2020."),
        "coverage": "Global (US/EU/Canada-dense)",
        "variables": ["nutrients", "carbon", "oxygen", "sediment", "ions"],
        "access": "auto",
        "overlaps": ["wqp", "gemstat", "glorich", "waterbase", "cesi",
                     "aquasat", "camels_chem"],
        "adapter": "grqa",
    },
    "wqp": {
        "name": "WQP — Water Quality Portal (USGS/EPA)",
        "description": ("US in-situ, open REST API. Nutrients, ions, metals, "
                        "sediment, DO, pH, etc. Most current US data."),
        "coverage": "United States",
        "variables": ["nutrients", "ions", "metals", "sediment", "oxygen"],
        "access": "auto",
        "overlaps": ["grqa"],
        "adapter": "wqp",
    },
    "waterbase": {
        "name": "Waterbase — EEA Water Quality (Europe)",
        "description": ("European rivers/lakes disaggregated water-quality "
                        "observations (EEA WISE, 1900-2019): nutrients, oxygen, "
                        "carbon, ions. Site coordinates come from the EEA "
                        "Discodata spatial table; measurements from the "
                        "disaggregated CSV (first run downloads a 733 MB zip)."),
        "coverage": "Europe",
        "variables": ["nutrients", "ions", "metals", "oxygen"],
        "access": "auto",
        "overlaps": ["grqa"],
        "adapter": "waterbase",
        # Correct live link (the old sdi.eea.europa.eu path 404s): the WISE-6
        # disaggregated-data zip. Coordinates are joined from Discodata, not here.
        "url": "https://cmshare.eea.europa.eu/s/yyAYF45k6nFAsZB/download",
        "cache": "waterbase_disaggregated.zip",
        "member": "Waterbase_v2020_1_T_WISE6_DisaggregatedData.csv",
        "sep": ",",
        "fmt": "long",
    },
    "camels_chem": {
        "name": "CAMELS-Chem — US catchment stream chemistry",
        "description": ("516 US headwater catchments, 1980-2018, stream "
                        "chemistry (NO3, TN/TON/DON, DOC/TOC, major ions, DO, "
                        "pH, Si). Gauge lat/lon are resolved from the USGS NWIS "
                        "site service (gauge_id = USGS site number)."),
        "coverage": "United States (headwaters)",
        "variables": ["nutrients", "carbon", "ions", "oxygen"],
        "access": "auto",
        "overlaps": ["wqp", "grqa"],
        "adapter": "camels",
        "url": ("https://www.hydroshare.org/resource/841f5e85085c423f889ac809c1bed4ac/"
                "data/contents/Camels_Chem_%20Dataset/Camels_chem_1980_2018.csv"),
        "cache": "camels_chem.csv",
        "sep": ",",
        "fmt": "wide",   # constituent-per-column; coords via USGS (see adapter)
    },
    "gloria": {
        "name": "GLORIA — global optical water-quality (in-situ)",
        "description": ("Global hyperspectral reflectance with co-located "
                        "chlorophyll-a, TSS, CDOM and Secchi depth (450 water "
                        "bodies). Remote-sensing oriented."),
        "coverage": "Global",
        "variables": ["sediment", "chlorophyll", "cdom", "clarity"],
        "access": "auto",
        "overlaps": [],
        "adapter": "bulk",
        "url": "https://download.pangaea.de/dataset/948492/files/GLORIA-2022.zip",
        "cache": "gloria_2022.zip",
        "member": "GLORIA_2022/GLORIA_meta_and_lab.csv",
        "sep": ",",
        "fmt": "wide",
        "cols": {"station": "GLORIA_ID", "lat": "Latitude", "lon": "Longitude",
                 "date": "Date_Time_UTC"},
        "wide_map": {"TSS": "TSS", "Chla": "chl-a", "aCDOM440": "CDOM", "Secchi_depth": "Secchi"},
        "units": {"TSS": "mg/L", "chl-a": "ug/L", "CDOM": "1/m", "Secchi": "m"},
    },
    "aquasat": {
        "name": "AquaSat — US remote-sensing matchups",
        "description": ("600k+ US matchups (1984-2019) of in-situ TSS, "
                        "chlorophyll-a, DOC, CDOM, Secchi with Landsat. In-situ "
                        "part sourced from WQP."),
        "coverage": "United States",
        "variables": ["sediment", "chlorophyll", "carbon", "cdom", "clarity"],
        "access": "auto",
        "overlaps": ["wqp", "grqa"],
        "adapter": "bulk",
        "url": "https://ndownloader.figshare.com/files/18733733",  # sr_wq_rs_join.csv (291 MB)
        "cache": "aquasat_insitu.csv",
        "sep": ",",
        "fmt": "wide",
        "cols": {"station": "SiteID", "lat": "lat", "lon": "long", "date": "date"},
        "wide_map": {"tss": "TSS", "chl_a": "chl-a", "doc": "DOC",
                     "cdom": "CDOM", "secchi": "Secchi"},
        "units": {"TSS": "mg/L", "chl-a": "ug/L", "DOC": "mg/L",
                  "CDOM": "1/m", "Secchi": "m"},
    },
    "gemstat": {
        "name": "GEMStat — UNEP GEMS/Water (open subset)",
        "description": ("Global in-situ freshwater quality, open (CC-BY) subset "
                        "on Zenodo. Nutrients, ions, oxygen, carbon."),
        "coverage": "Global",
        "variables": ["nutrients", "ions", "oxygen", "carbon"],
        "access": "auto",
        "overlaps": ["grqa"],
        "adapter": "gemstat",
        "url": "https://zenodo.org/records/14230628/files/GFQA_v2.zip?download=1",
        "cache": "gemstat_GFQA_v2.zip",
    },
    # ── gated: online but behind a Terms/registration wall → manual download ──
    "gemstat_full": {
        "name": "GEMStat — full global database (gated)",
        "description": ("Complete GEMStat (all licences). Portal caps at 675 "
                        "stations; global data only via email request to "
                        "gwdc@bafg.de. Provide the downloaded file path."),
        "coverage": "Global",
        "variables": ["nutrients", "ions", "oxygen", "carbon", "metals"],
        "access": "manual",
        "overlaps": ["grqa", "gemstat"],
        "adapter": "manual",
        "url": "https://gemstat.org/data-gemstat/data-portal/custom-data-request/",
        "fmt": "long",
        "cols": {"station": "GEMS_station_number", "lat": "Latitude", "lon": "Longitude",
                 "date": "Sample_Date", "param": "Parameter_Code",
                 "value": "Value", "units": "Unit"},
    },
    "grdc": {
        "name": "GRDC — Global Runoff Data Centre (discharge, gated)",
        "description": ("River DISCHARGE (not concentration) — needed for "
                        "load-based calibration. No API: download via portal "
                        "form + Terms, link emailed in 24h. Provide the path."),
        "coverage": "Global",
        "variables": ["discharge"],
        "access": "manual",
        "overlaps": [],
        "adapter": "manual",
        "url": "https://grdc.bafg.de/data/data_portal/",
        "fmt": "long",
        "cols": {"station": "grdc_no", "lat": "lat", "lon": "long",
                 "date": "date", "value": "value", "units": "units", "param": "parameter"},
    },
}

# Back-compat / friendly aliases accepted from the config.
SOURCE_ALIASES = {
    "user_csv": "user_csv", "skip": "skip", "none": "skip",
    "gems": "gemstat", "gemstat_open": "gemstat",
    "eea": "waterbase", "camelschem": "camels_chem", "camels-chem": "camels_chem",
}

# ── WQP CharacteristicName → model species. Keys are the EXACT WQP
#    controlled-vocabulary names (used verbatim in the query — do NOT re-case
#    them, WQP rejects unknown names with HTTP 400). Response matching is
#    case-insensitive via _WQP_CHAR_LOWER. ──
WQP_PARAM_MAP = {
    "Nitrate": "NO3-N",
    "Nitrate as N": "NO3-N",
    "Nitrite": "NO2-N",
    "Ammonia and ammonium": "NH4-N",
    "Ammonia-nitrogen as N": "NH4-N",
    "Orthophosphate": "PO4-P",
    "Phosphate-phosphorus": "PO4-P",
    "Phosphorus": "TP",
    "Total suspended solids": "TSS",
    "Suspended sediment concentration (SSC)": "TSS",
    "Dissolved oxygen (DO)": "DO",
    "pH": "pH",
    "Organic carbon": "DOC",
    "Chlorophyll a": "chl-a",
    "Temperature, water": "WTEMP",
}
_WQP_CHAR_LOWER = {k.lower(): v for k, v in WQP_PARAM_MAP.items()}

# NOTE: Waterbase, CAMELS-Chem and GEMStat use dedicated adapters with their own
# parameter maps (WATERBASE_LABEL_MAP, CAMELS_WIDE_MAP, GEMSTAT_FILE_CODES). The
# generic long-format `param_map` path in _adapter_bulk stays available for any
# future long bulk source that sets `param_map` in its registry entry.


# ─────────────────────────────────────────────────────────────────────────────
#  Selection normalization + redundancy resolution
# ─────────────────────────────────────────────────────────────────────────────
def normalize_selection(selected) -> List[str]:
    """Accept a single string (back-compat) or a list; return validated keys.

    'user_csv' and 'skip' are passed through untouched (handled by the caller).
    Unknown names are dropped with a warning.
    """
    if selected is None:
        return ["grqa"]
    if isinstance(selected, str):
        selected = [selected]
    out: List[str] = []
    for s in selected:
        key = str(s).strip().lower()
        if not key:
            continue
        key = SOURCE_ALIASES.get(key, key)
        if key in ("user_csv", "skip"):
            out.append(key)
        elif key in OBSERVATION_SOURCES:
            out.append(key)
        else:
            print(f"  WARNING: unknown observation source '{s}' — ignored. "
                  f"Valid: {', '.join(list(OBSERVATION_SOURCES) + ['user_csv', 'skip'])}")
    # de-duplicate while preserving order
    seen = set()
    return [x for x in out if not (x in seen or seen.add(x))]


# Aggregation level — higher = larger/aggregate dataset. When two sources
# overlap, the lower-rank one is (mostly) contained in the higher-rank one, which
# is what makes the "X is inside Y" flag point the right way.
_SOURCE_RANK = {
    "grqa": 3, "gemstat_full": 3,
    "wqp": 2, "waterbase": 2, "gemstat": 2,
    "aquasat": 1, "camels_chem": 1, "gloria": 1, "grdc": 1,
}


def detect_overlaps(selected: List[str]) -> List[Tuple[str, str]]:
    """Return (contained, container) pairs among the selected data sources.

    Two sources overlap if either lists the other in its `overlaps` set; the
    lower-`_SOURCE_RANK` source is reported as contained in the higher-rank one.
    Each contained source is reported once, against its largest container.
    """
    sel = [s for s in selected if s in OBSERVATION_SOURCES]
    seen = set()
    best: Dict[str, str] = {}
    for a in sel:
        for b in sel:
            if a == b:
                continue
            related = (b in OBSERVATION_SOURCES[a].get("overlaps", [])
                       or a in OBSERVATION_SOURCES[b].get("overlaps", []))
            if not related:
                continue
            fkey = frozenset((a, b))
            if fkey in seen:
                continue
            seen.add(fkey)
            ra, rb = _SOURCE_RANK.get(a, 1), _SOURCE_RANK.get(b, 1)
            inside, container = (a, b) if ra <= rb else (b, a)
            # keep the largest container for each contained source
            if (inside not in best
                    or _SOURCE_RANK.get(container, 1) > _SOURCE_RANK.get(best[inside], 1)):
                best[inside] = container
    return [(k, v) for k, v in best.items()]


def resolve_sources(selected, interactive: bool = True,
                    dedup_mode: str = "ask") -> Tuple[List[str], bool]:
    """Normalize the selection, flag redundant overlaps, and decide dedup.

    Returns (data_source_keys, do_dedup). 'user_csv'/'skip' are kept in the list
    for the caller to handle but never trigger the dedup prompt on their own.
    dedup_mode: 'ask' (prompt if interactive TTY), 'true', or 'false'.
    """
    norm = normalize_selection(selected)
    data_sources = [s for s in norm if s in OBSERVATION_SOURCES]
    overlaps = detect_overlaps(data_sources)

    if not overlaps:
        return norm, False

    print("\n  ⚠  Redundant/overlapping observation sources selected:")
    for a, b in overlaps:
        print(f"       • {OBSERVATION_SOURCES[a]['name']}  is already inside  "
              f"{OBSERVATION_SOURCES[b]['name']}")
    print("     Keeping both may duplicate the same observations.")

    mode = str(dedup_mode).strip().lower()
    if mode in ("true", "yes", "1", "remove"):
        print("     -> observation_dedup_overlaps=true: duplicate rows will be removed.")
        return norm, True
    if mode in ("false", "no", "0", "keep"):
        print("     -> observation_dedup_overlaps=false: keeping all rows (may duplicate).")
        return norm, False

    # mode == "ask"
    if interactive and sys.stdin and sys.stdin.isatty():
        try:
            ans = input("     Remove duplicate observations across overlapping "
                        "sources? [Y/n]: ").strip().lower()
        except (EOFError, KeyboardInterrupt):
            ans = "y"
        do_dedup = ans in ("", "y", "yes")
        print(f"     -> {'removing duplicates' if do_dedup else 'keeping all rows'}.")
        return norm, do_dedup

    # non-interactive default: dedup (safer than silent duplication)
    print("     (non-interactive: removing duplicate rows by default; set "
          "observation_dedup_overlaps='false' to keep all.)")
    return norm, True


# ─────────────────────────────────────────────────────────────────────────────
#  Harmonization helpers
# ─────────────────────────────────────────────────────────────────────────────
def _empty_harmonized() -> pd.DataFrame:
    return pd.DataFrame(columns=HARMONIZED_COLUMNS)


def _split_datetime(series) -> pd.DataFrame:
    # Cast to string FIRST. Some sources store the sampling date as a compact
    # integer YYYYMMDD (e.g. Waterbase '20061004'); read as int64, pandas'
    # to_datetime treats those as NANOSECONDS-since-epoch → every row collapses
    # to 1970-01-01. As strings, ISO ('2008-04-07'), US ('1/31/1980'), ISO-T and
    # YYYYMMDD all parse correctly.
    s = pd.Series(series).astype("string").str.strip()
    dt = pd.to_datetime(s, errors="coerce")
    # Explicit fallback for any 8-digit YYYYMMDD the generic parser missed.
    _miss = dt.isna() & s.str.fullmatch(r"\d{8}").fillna(False)
    if _miss.any():
        dt = dt.fillna(pd.to_datetime(s.where(_miss), format="%Y%m%d",
                                      errors="coerce"))
    return pd.DataFrame({
        "year": dt.dt.year, "month": dt.dt.month, "day": dt.dt.day,
        "minute": dt.dt.hour.fillna(0) * 60 + dt.dt.minute.fillna(0),
    })


def _finalize(df: pd.DataFrame, source: str) -> pd.DataFrame:
    """Coerce a partly-built frame to the harmonized schema."""
    if df is None or len(df) == 0:
        return _empty_harmonized()
    for c in HARMONIZED_COLUMNS:
        if c not in df.columns:
            df[c] = None
    df["source"] = source
    df = df[HARMONIZED_COLUMNS].copy()
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    df = df.dropna(subset=["lat", "lon", "parameter", "value"])
    for c in ("year", "month", "day", "minute"):
        df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
    df["station_id"] = df["station_id"].astype(str)
    df["parameter"] = df["parameter"].astype(str)
    return df


def _clip_bbox(df, lat_col, lon_col, bbox):
    """Keep rows inside bbox=(lonW, latS, lonE, latN)."""
    if bbox is None:
        return df
    w, s, e, n = bbox
    la = pd.to_numeric(df[lat_col], errors="coerce")
    lo = pd.to_numeric(df[lon_col], errors="coerce")
    return df[(la >= s) & (la <= n) & (lo >= w) & (lo <= e)].copy()


# ─────────────────────────────────────────────────────────────────────────────
#  Adapters — each returns a harmonized DataFrame
# ─────────────────────────────────────────────────────────────────────────────
def _download_cached(url: str, dest: Path) -> Optional[Path]:
    """Download `url` to `dest` (cached). Returns None on failure."""
    if dest.exists() and dest.stat().st_size > 0:
        print(f"  Using cached file: {dest}")
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    print(f"  Downloading {url}\n    -> {dest}")
    try:
        # Browser-like UA + streaming: some hosts (e.g. Zenodo) 403 the default
        # Python agent, and streaming avoids loading big files fully into RAM.
        req = urllib.request.Request(url, headers={
            "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                          "AppleWebKit/537.36"})
        with urllib.request.urlopen(req, timeout=300) as r, open(dest, "wb") as fh:  # noqa: S310
            shutil.copyfileobj(r, fh)
        return dest
    except Exception as exc:
        print(f"  WARNING: download failed ({exc}). Skipping this source.")
        try:
            dest.unlink()
        except OSError:
            pass
        return None


def _read_table(path: Path, sep: str, member: str = None,
                encoding: str = None) -> Optional[pd.DataFrame]:
    """Read a CSV, or a specific `member` (else the first CSV) inside a zip."""
    try:
        if str(path).lower().endswith(".zip"):
            import zipfile
            with zipfile.ZipFile(path) as zf:
                name = member or next((n for n in zf.namelist()
                             if n.lower().endswith((".csv", ".txt", ".tab"))), None)
                if not name or name not in zf.namelist():
                    print(f"  WARNING: member '{member}' not in {path}.")
                    return None
                with zf.open(name) as fh:
                    return pd.read_csv(fh, sep=sep, low_memory=False, encoding=encoding)
        return pd.read_csv(path, sep=sep, low_memory=False, encoding=encoding)
    except Exception as exc:
        print(f"  WARNING: could not parse {path} ({exc}).")
        return None


def _adapter_grqa(key, ctx) -> pd.DataFrame:
    """Delegate to the existing GRQA extractor, then harmonize its output."""
    try:
        from Gen_GRQA_Extract import GRQACalibrationExtractor, SpeciesMapper
    except Exception as exc:
        print(f"  WARNING: GRQA extractor unavailable ({exc}).")
        return _empty_harmonized()

    grqa_mapping = ctx.get("grqa_species_mapping")     # {grqa_param: model_name}
    if not grqa_mapping:
        print("  WARNING: GRQA selected but no GRQA species mapping provided; skipping.")
        return _empty_harmonized()
    mapper = SpeciesMapper(mapping=grqa_mapping)
    out_dir = os.path.join(ctx["output_dir"], "obs_cache_grqa")
    extractor = GRQACalibrationExtractor(
        output_dir=out_dir, species_mapper=mapper,
        local_data_path=ctx.get("grqa_local_data_path"),
        buffer_distance_m=ctx["buffer_m"])
    stations, obs = extractor.extract_stations_and_observations(ctx["search_area_gdf"])
    if obs is None or len(obs) == 0:
        return _empty_harmonized()
    cm = extractor.column_map or {}
    dt = _split_datetime(obs[cm.get("obs_date", "obs_date")])
    out = pd.DataFrame({
        "station_id": obs[cm.get("site_id", "site_id")].astype(str),
        "lat": obs[cm.get("lat", "lat_wgs84")],
        "lon": obs[cm.get("lon", "lon_wgs84")],
        "parameter": obs["model_species"],
        "value": obs[cm.get("obs_value", "obs_value")],
        "units": obs[cm.get("unit", "unit")] if cm.get("unit") in obs.columns else "",
    })
    out = pd.concat([out.reset_index(drop=True), dt.reset_index(drop=True)], axis=1)
    return _finalize(out, "GRQA")


def _adapter_wqp(key, ctx) -> pd.DataFrame:
    """Water Quality Portal live REST API (US), clipped to the bbox."""
    bbox = ctx["bbox"]
    if bbox is None:
        print("  WARNING: WQP needs a bounding box; skipping.")
        return _empty_harmonized()
    w, s, e, n = bbox
    # only request characteristics that map back to a target model species
    target_species = set(ctx.get("species") or [])
    chars = sorted({c for c in WQP_PARAM_MAP if WQP_PARAM_MAP[c] in target_species})
    if not chars:
        print("  WQP: none of the target species are WQP-mappable; skipping.")
        return _empty_harmonized()
    params = [("bBox", f"{w},{s},{e},{n}"), ("mimeType", "csv"), ("zip", "no"),
              ("dataProfile", "resultPhysChem"), ("providers", "NWIS"),
              ("providers", "STORET")]
    for c in chars:
        params.append(("characteristicName", c))   # exact WQP name, verbatim
    yrs = ctx.get("years")
    if yrs:
        params += [("startDateLo", f"01-01-{yrs[0]}"), ("startDateHi", f"12-31-{yrs[1]}")]
    from urllib.parse import urlencode
    url = "https://www.waterqualitydata.us/data/Result/search?" + urlencode(params)
    print(f"  Querying WQP API (bbox={w:.2f},{s:.2f},{e:.2f},{n:.2f})...")
    try:
        with urllib.request.urlopen(url, timeout=120) as resp:  # noqa: S310
            raw = resp.read().decode("utf-8", "replace")
        df = pd.read_csv(io.StringIO(raw), low_memory=False)
    except Exception as exc:
        print(f"  WARNING: WQP query failed ({exc}). Skipping.")
        ctx["_status"] = "download_failed"
        return _empty_harmonized()
    if df is None or len(df) == 0:
        print("  WQP: no records in area.")
        return _empty_harmonized()
    lat_c = _first_present(df, ["ActivityLocation/LatitudeMeasure",
                                "LatitudeMeasure", "lat"])
    lon_c = _first_present(df, ["ActivityLocation/LongitudeMeasure",
                                "LongitudeMeasure", "lon"])
    if lat_c is None or lon_c is None:
        print("  WQP: response lacks coordinates; skipping.")
        return _empty_harmonized()
    dt = _split_datetime(df.get("ActivityStartDate"))
    char = df.get("CharacteristicName", pd.Series([""] * len(df))).astype(str).str.lower()
    out = pd.DataFrame({
        "station_id": df.get("MonitoringLocationIdentifier", "").astype(str),
        "lat": df[lat_c], "lon": df[lon_c],
        "parameter": char.map(_WQP_CHAR_LOWER),
        "value": df.get("ResultMeasureValue"),
        "units": df.get("ResultMeasure/MeasureUnitCode", ""),
    })
    out = pd.concat([out.reset_index(drop=True), dt.reset_index(drop=True)], axis=1)
    out = out.dropna(subset=["parameter"])
    return _finalize(out, "WQP")


def _first_present(df, names):
    for c in names:
        if c in df.columns:
            return c
    return None


def _adapter_bulk(key, ctx) -> pd.DataFrame:
    """Generic 'download big file → clip by bbox → map params → harmonize'."""
    cfg = OBSERVATION_SOURCES[key]
    cache = Path(ctx["cache_dir"]) / cfg.get("cache", f"{key}.csv")
    path = _download_cached(cfg["url"], cache)
    if path is None:
        ctx["_status"] = "download_failed"
        return _empty_harmonized()
    df = _read_table(path, cfg.get("sep", ","), cfg.get("member"), cfg.get("encoding"))
    if df is None or len(df) == 0:
        ctx["_status"] = "download_failed"
        return _empty_harmonized()

    cols = cfg["cols"]
    missing = [cols[k] for k in ("lat", "lon") if cols.get(k) not in df.columns]
    if missing:
        print(f"  WARNING: {key}: expected columns {missing} not found "
              f"(available: {list(df.columns)[:12]}...). Adjust registry 'cols'. Skipping.")
        ctx["_status"] = "unavailable"
        return _empty_harmonized()

    df = _clip_bbox(df, cols["lat"], cols["lon"], ctx["bbox"])
    if len(df) == 0:
        print(f"  {key}: no observations in area.")
        return _empty_harmonized()

    dt = _split_datetime(df[cols["date"]]) if cols.get("date") in df.columns \
        else pd.DataFrame({"year": 0, "month": 0, "day": 0, "minute": 0}, index=df.index)

    if cfg.get("fmt") == "wide":
        # one column per parameter -> melt to long
        frames = []
        for src_col, species in cfg.get("wide_map", {}).items():
            if src_col not in df.columns:
                continue
            sub = pd.DataFrame({
                "station_id": df[cols["station"]].astype(str),
                "lat": df[cols["lat"]], "lon": df[cols["lon"]],
                "parameter": species, "value": df[src_col],
                "units": cfg.get("units", {}).get(species, ""),
            })
            sub = pd.concat([sub.reset_index(drop=True), dt.reset_index(drop=True)], axis=1)
            frames.append(sub)
        out = pd.concat(frames, ignore_index=True) if frames else _empty_harmonized()
    else:
        # Map native parameter names → BGC/model species. With a param_map,
        # unmapped names become NaN and are dropped by _finalize, so the output
        # only carries BGC-species names (what calibration matches against).
        pmap = cfg.get("param_map")
        if cols.get("param") in df.columns:
            raw = df[cols["param"]].astype(str).str.strip()
            param = (raw.str.lower().map({k.lower(): v for k, v in pmap.items()})
                     if pmap else raw)
        else:
            param = key
        out = pd.DataFrame({
            "station_id": df[cols["station"]].astype(str),
            "lat": df[cols["lat"]], "lon": df[cols["lon"]],
            "parameter": param,
            "value": df[cols["value"]] if cols.get("value") in df.columns else None,
            "units": df[cols["units"]] if cols.get("units") in df.columns else "",
        })
        out = pd.concat([out.reset_index(drop=True), dt.reset_index(drop=True)], axis=1)

    return _finalize(out, cfg["name"].split(" —")[0])


def _adapter_manual(key, ctx) -> pd.DataFrame:
    """Gated source: ingest a user-downloaded file, or explain how to get it."""
    cfg = OBSERVATION_SOURCES[key]
    path = (ctx.get("manual_paths") or {}).get(key, "")
    if not path or not os.path.isfile(path):
        print(f"\n  ℹ  '{key}' ({cfg['name']}) is a gated dataset with no open API.")
        print(f"     Download once from: {cfg['url']}")
        print(f"     Then set observation_source_manual_paths['{key}'] = "
              f"'/path/to/file.csv'. Skipping for now.")
        ctx["_status"] = "needs_manual_path"
        return _empty_harmonized()
    df = _read_table(Path(path), cfg.get("sep", ","))
    if df is None or len(df) == 0:
        ctx["_status"] = "download_failed"
        return _empty_harmonized()
    cols = cfg["cols"]
    if cols.get("lat") not in df.columns or cols.get("lon") not in df.columns:
        print(f"  WARNING: {key}: '{path}' missing lat/lon columns "
              f"({cols.get('lat')}, {cols.get('lon')}). Skipping.")
        ctx["_status"] = "unavailable"
        return _empty_harmonized()
    df = _clip_bbox(df, cols["lat"], cols["lon"], ctx["bbox"])
    dt = _split_datetime(df[cols["date"]]) if cols.get("date") in df.columns \
        else pd.DataFrame({"year": 0, "month": 0, "day": 0, "minute": 0}, index=df.index)
    out = pd.DataFrame({
        "station_id": df[cols["station"]].astype(str) if cols.get("station") in df.columns else "NA",
        "lat": df[cols["lat"]], "lon": df[cols["lon"]],
        "parameter": df[cols["param"]] if cols.get("param") in df.columns else key,
        "value": df[cols["value"]] if cols.get("value") in df.columns else None,
        "units": df[cols["units"]] if cols.get("units") in df.columns else "",
    })
    out = pd.concat([out.reset_index(drop=True), dt.reset_index(drop=True)], axis=1)
    return _finalize(out, cfg["name"].split(" —")[0])


# GEMStat GFQA_v2.zip layout: one semicolon-delimited (latin-1) CSV per parameter
# group, keyed by GEMS.Station.Number; lat/lon live in GEMStat_station_metadata.csv
# with COMMA decimals. Map: {file: {Parameter.Code: BGC species}}.
GEMSTAT_FILE_CODES = {
    "Oxidized_Nitrogen.csv": {"NO3N": "NO3-N", "NO2N": "NO2-N"},
    "Other_Nitrogen.csv":    {"NH4N": "NH4-N", "NH3N": "NH4-N", "TN": "TN", "TDN": "TN"},
    "Phosphorus.csv":        {"DRP": "PO4-P", "DIP": "PO4-P", "TP": "TP", "TDP": "TP"},
    "pH.csv":                {"pH": "pH"},
    "Dissolved_Gas.csv":     {"O2-Dis": "DO"},
}


def _adapter_gemstat(key, ctx) -> pd.DataFrame:
    """GEMStat (UNEP GEMS/Water open subset) — GFQA_v2.zip: per-parameter CSVs
    joined to station metadata for lat/lon. Reads only the parameter files that
    contain a target BGC species."""
    import zipfile
    cfg = OBSERVATION_SOURCES[key]
    cache = Path(ctx["cache_dir"]) / cfg.get("cache", "GFQA_v2.zip")
    path = _download_cached(cfg["url"], cache)
    if path is None:
        ctx["_status"] = "download_failed"
        return _empty_harmonized()
    try:
        z = zipfile.ZipFile(path)
        with z.open("GEMStat_station_metadata.csv") as f:
            st = pd.read_csv(f, sep=";", encoding="latin-1", low_memory=False)
    except Exception as exc:
        print(f"  WARNING: GEMStat archive unreadable ({exc}).")
        ctx["_status"] = "download_failed"
        return _empty_harmonized()

    def _num(series):
        return pd.to_numeric(series.astype(str).str.replace(",", ".", regex=False),
                             errors="coerce")
    st["_lat"] = _num(st["Latitude"])
    st["_lon"] = _num(st["Longitude"])
    stmap = st.dropna(subset=["_lat", "_lon"]).set_index(
        "GEMS Station Number")[["_lat", "_lon"]]

    target = set(ctx.get("species") or [])
    frames = []
    for fname, code_map in GEMSTAT_FILE_CODES.items():
        wanted = {c: sp for c, sp in code_map.items() if sp in target}
        if not wanted:
            continue
        try:
            with z.open(fname) as f:
                d = pd.read_csv(f, sep=";", encoding="latin-1", low_memory=False)
        except Exception:
            continue
        d = d[d["Parameter.Code"].isin(wanted.keys())]
        if len(d) == 0:
            continue
        d = d.join(stmap, on="GEMS.Station.Number")
        d = _clip_bbox(d, "_lat", "_lon", ctx["bbox"])
        if len(d) == 0:
            continue
        dt = _split_datetime(d["Sample.Date"])
        out = pd.DataFrame({
            "station_id": d["GEMS.Station.Number"].astype(str),
            "lat": d["_lat"], "lon": d["_lon"],
            "parameter": d["Parameter.Code"].map(wanted),
            "value": _num(d["Value"]),
            "units": d["Unit"],
        })
        out = pd.concat([out.reset_index(drop=True), dt.reset_index(drop=True)], axis=1)
        frames.append(out)
    if not frames:
        return _empty_harmonized()
    return _finalize(pd.concat(frames, ignore_index=True), "GEMStat")


# ─────────────────────────────────────────────────────────────────────────────
#  CAMELS-Chem (US) — wide constituent columns keyed by USGS gauge_id. The file
#  carries NO coordinates, so lat/lon are resolved from the USGS NWIS site
#  service (gauge_id = USGS site number with leading zeros dropped) and cached.
# ─────────────────────────────────────────────────────────────────────────────
_BROWSER_UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
               "AppleWebKit/537.36")

# CAMELS-Chem column -> model/BGC species. `o` is dissolved oxygen (values
# ~0-15 mg/L, NOT the δ18O isotope — verified against the data distribution).
CAMELS_WIDE_MAP = {
    "no3": "NO3-N", "tdn": "TN", "ton": "TON", "don": "DON", "doc": "DOC",
    "toc": "TOC", "ph": "pH", "o": "DO", "so4": "SO4", "cl": "Cl",
    "si": "Si", "na": "Na", "ca": "Ca", "k": "K", "mg": "Mg",
}
CAMELS_UNITS = {"pH": "", "DO": "mg/L", "DOC": "mg/L", "TOC": "mg/L",
                "NO3-N": "mg/L", "TN": "mg/L", "TON": "mg/L", "DON": "mg/L",
                "SO4": "mg/L", "Cl": "mg/L", "Si": "mg/L", "Na": "mg/L",
                "Ca": "mg/L", "K": "mg/L", "Mg": "mg/L"}


def _usgs_gauge_coords(gauge_ids, cache_path: Path) -> Dict[str, Tuple[float, float]]:
    """Resolve USGS site coordinates for CAMELS gauge_ids (8-digit, zero-padded)
    via the NWIS site service, cached to `cache_path`.

    Returns {padded_site_no: (lat, lon)}.
    """
    coords: Dict[str, Tuple[float, float]] = {}
    if cache_path.exists():
        try:
            c = pd.read_csv(cache_path, dtype={"site_no": str})
            coords = {str(r.site_no): (float(r.lat), float(r.lon))
                      for r in c.itertuples()}
        except Exception:
            coords = {}
    need = sorted({g for g in gauge_ids if g not in coords})
    for i in range(0, len(need), 150):          # NWIS accepts many sites/call
        batch = need[i:i + 150]
        url = ("https://waterservices.usgs.gov/nwis/site/?sites="
               + ",".join(batch) + "&format=rdb")
        if i == 0:
            print(f"  Resolving {len(need)} USGS gauge coordinates...")
        try:
            req = urllib.request.Request(url, headers={"User-Agent": _BROWSER_UA})
            with urllib.request.urlopen(req, timeout=90) as r:   # noqa: S310
                txt = r.read().decode("utf-8", "replace")
        except Exception as exc:
            print(f"  WARNING: USGS site lookup failed for a batch ({exc}).")
            continue
        rows = [ln for ln in txt.splitlines() if ln and not ln.startswith("#")]
        if len(rows) < 3:                        # header + dtype row + data
            continue
        hdr = rows[0].split("\t")
        try:
            si, la, lo = (hdr.index("site_no"), hdr.index("dec_lat_va"),
                          hdr.index("dec_long_va"))
        except ValueError:
            continue
        for ln in rows[2:]:                       # row[1] is the RDB dtype row
            f = ln.split("\t")
            try:
                coords[f[si]] = (float(f[la]), float(f[lo]))
            except (ValueError, IndexError):
                pass
    if need and coords:
        try:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            pd.DataFrame([{"site_no": k, "lat": v[0], "lon": v[1]}
                          for k, v in coords.items()]).to_csv(cache_path, index=False)
        except Exception:
            pass
    return coords


def _adapter_camels(key, ctx) -> pd.DataFrame:
    """CAMELS-Chem: download the wide CSV, resolve gauge lat/lon from USGS,
    clip to the bbox, melt constituent columns to species."""
    cfg = OBSERVATION_SOURCES[key]
    cache_dir = Path(ctx["cache_dir"])
    csv_path = _download_cached(cfg["url"],
                                cache_dir / cfg.get("cache", "camels_chem.csv"))
    if csv_path is None:
        ctx["_status"] = "download_failed"
        return _empty_harmonized()
    df = _read_table(csv_path, cfg.get("sep", ","))
    if df is None or len(df) == 0 or "gauge_id" not in df.columns:
        print("  WARNING: CAMELS-Chem file unreadable or missing 'gauge_id'.")
        ctx["_status"] = "download_failed"
        return _empty_harmonized()
    # USGS site numbers are 8-digit; CAMELS gauge_id drops leading zeros.
    df["_site"] = df["gauge_id"].apply(
        lambda g: str(g).strip().split(".")[0].zfill(8))
    coords = _usgs_gauge_coords(set(df["_site"]),
                                cache_dir / "camels_gauge_coords.csv")
    if not coords:
        print("  WARNING: could not resolve CAMELS gauge coordinates; skipping.")
        ctx["_status"] = "download_failed"
        return _empty_harmonized()
    df["lat"] = df["_site"].map(lambda s: coords.get(s, (None, None))[0])
    df["lon"] = df["_site"].map(lambda s: coords.get(s, (None, None))[1])
    df = df.dropna(subset=["lat", "lon"])
    df = _clip_bbox(df, "lat", "lon", ctx["bbox"])
    if len(df) == 0:
        print(f"  {key}: no observations in area.")
        return _empty_harmonized()
    dt = _split_datetime(df["sample_start_dt"]) if "sample_start_dt" in df.columns \
        else pd.DataFrame({"year": 0, "month": 0, "day": 0, "minute": 0},
                          index=df.index)
    frames = []
    for src_col, species in CAMELS_WIDE_MAP.items():
        if src_col not in df.columns:
            continue
        sub = pd.DataFrame({
            "station_id": df["gauge_id"].astype(str),
            "lat": df["lat"], "lon": df["lon"],
            "parameter": species, "value": df[src_col],
            "units": CAMELS_UNITS.get(species, ""),
        })
        sub = pd.concat([sub.reset_index(drop=True), dt.reset_index(drop=True)],
                        axis=1)
        frames.append(sub)
    if not frames:
        return _empty_harmonized()
    return _finalize(pd.concat(frames, ignore_index=True), "CAMELS-Chem")


# ─────────────────────────────────────────────────────────────────────────────
#  Waterbase (EEA WISE, Europe). The disaggregated measurement table has NO
#  coordinates, and the EEA Discodata SQL endpoint times out on measurement
#  pulls. So: fetch region site coordinates from Discodata's small SPATIAL table
#  (reliable, bbox-filtered), then read the measurements from the disaggregated
#  CSV (cached 733 MB zip → 4.3 GB) IN CHUNKS, keeping only rows at region sites
#  with a target determinand. Memory-safe; the first run downloads the zip.
# ─────────────────────────────────────────────────────────────────────────────
_DISCODATA_SQL = "https://discodata.eea.europa.eu/sql"
_WATERBASE_SPATIAL = ("[WISE_SOE].[latest].[Waterbase_S_WISE_SpatialObject_"
                      "DerivedData]")

# WISE determinand label -> model/BGC species. Units are carried through as
# reported (WISE gives e.g. Nitrate in mg{NO3}/L, Ammonium in mg{NH4}/L).
WATERBASE_LABEL_MAP = {
    "Nitrate": "NO3-N", "Nitrite": "NO2-N", "Ammonium": "NH4-N",
    "Total nitrogen": "TN", "Total oxidised nitrogen": "NO3-N",
    "Phosphate": "PO4-P", "Orthophosphate": "PO4-P", "Total phosphorus": "TP",
    "Dissolved oxygen": "DO", "pH": "pH",
    "Dissolved organic carbon": "DOC", "Total organic carbon": "TOC",
}


def _waterbase_site_coords(bbox) -> Dict[str, Tuple[float, float]]:
    """Region monitoring-site coordinates from the Waterbase spatial table
    (Discodata). Returns {monitoringSiteIdentifier: (lat, lon)}; {} on failure
    (caller treats empty as 'no sites / service down')."""
    import urllib.parse
    if bbox is None:
        return {}
    w, s, e, n = bbox
    coords: Dict[str, Tuple[float, float]] = {}
    for page in range(1, 26):                    # up to 25*2000 = 50k sites
        sql = (f"SELECT monitoringSiteIdentifier AS sid, lat, lon FROM "
               f"{_WATERBASE_SPATIAL} WHERE lat BETWEEN {s} AND {n} "
               f"AND lon BETWEEN {w} AND {e}")
        url = _DISCODATA_SQL + "?" + urllib.parse.urlencode(
            {"query": sql, "p": page, "nrOfHits": 2000})
        try:
            req = urllib.request.Request(url, headers={"User-Agent": _BROWSER_UA})
            with urllib.request.urlopen(req, timeout=60) as r:   # noqa: S310
                j = json.loads(r.read().decode("utf-8", "replace"))
        except Exception as exc:
            print(f"  WARNING: Waterbase site lookup failed ({exc}).")
            break
        if isinstance(j, dict) and j.get("errors"):
            print(f"  WARNING: Waterbase spatial query error "
                  f"({j['errors'][0].get('error')}).")
            break
        rows = j.get("results", []) if isinstance(j, dict) else []
        for row in rows:
            try:
                coords[str(row["sid"])] = (float(row["lat"]), float(row["lon"]))
            except (TypeError, ValueError, KeyError):
                pass
        if len(rows) < 2000:
            break
    return coords


def _adapter_waterbase(key, ctx) -> pd.DataFrame:
    """Waterbase (Europe): Discodata spatial coords + chunked read of the
    disaggregated CSV, filtered to region sites + target determinands."""
    cfg = OBSERVATION_SOURCES[key]
    bbox = ctx["bbox"]
    if bbox is None:
        print("  WARNING: Waterbase needs a bounding box; skipping.")
        return _empty_harmonized()
    coords = _waterbase_site_coords(bbox)
    if not coords:
        print("  Waterbase: no monitoring sites in area "
              "(or coordinate service unreachable).")
        return _empty_harmonized()
    region_sites = set(coords)

    cache = Path(ctx["cache_dir"]) / cfg.get("cache", "waterbase_disaggregated.zip")
    path = _download_cached(cfg["url"], cache)
    if path is None:
        ctx["_status"] = "download_failed"
        return _empty_harmonized()

    usecols = ["monitoringSiteIdentifier", "observedPropertyDeterminandLabel",
               "resultObservedValue", "resultUom", "phenomenonTimeSamplingDate"]
    frames, kept, cap = [], 0, 200000
    try:
        import zipfile
        with zipfile.ZipFile(path) as zf:
            member = cfg.get("member") or next(
                (nm for nm in zf.namelist() if nm.lower().endswith(".csv")), None)
            if member is None or member not in zf.namelist():
                print("  WARNING: Waterbase archive has no CSV member; skipping.")
                ctx["_status"] = "download_failed"
                return _empty_harmonized()
            print(f"  Scanning {member} for {len(region_sites)} region site(s)...")
            with zf.open(member) as fh:
                for chunk in pd.read_csv(fh, usecols=usecols,
                                         sep=cfg.get("sep", ","),
                                         chunksize=500000, low_memory=False):
                    sub = chunk[
                        chunk["monitoringSiteIdentifier"].astype(str).isin(region_sites)
                        & chunk["observedPropertyDeterminandLabel"].isin(WATERBASE_LABEL_MAP)]
                    if len(sub):
                        frames.append(sub)
                        kept += len(sub)
                        if kept >= cap:
                            print(f"  Waterbase: reached {cap}-row cap; "
                                  f"stopping scan early.")
                            break
    except Exception as exc:
        print(f"  WARNING: Waterbase CSV read failed ({exc}); skipping.")
        ctx["_status"] = "download_failed"
        return _empty_harmonized()

    if not frames:
        print("  Waterbase: no target determinands at region sites.")
        return _empty_harmonized()
    raw = pd.concat(frames, ignore_index=True)
    sid = raw["monitoringSiteIdentifier"].astype(str)
    dt = _split_datetime(raw["phenomenonTimeSamplingDate"])
    out = pd.DataFrame({
        "station_id": sid,
        "lat": sid.map(lambda s: coords.get(s, (None, None))[0]),
        "lon": sid.map(lambda s: coords.get(s, (None, None))[1]),
        "parameter": raw["observedPropertyDeterminandLabel"].map(WATERBASE_LABEL_MAP),
        "value": raw["resultObservedValue"],
        "units": raw["resultUom"],
    })
    out = pd.concat([out.reset_index(drop=True), dt.reset_index(drop=True)], axis=1)
    return _finalize(out, "Waterbase")


_ADAPTERS = {
    "grqa": _adapter_grqa, "wqp": _adapter_wqp, "gemstat": _adapter_gemstat,
    "camels": _adapter_camels, "waterbase": _adapter_waterbase,
    "bulk": _adapter_bulk, "manual": _adapter_manual,
}


def extract_one(key: str, ctx: dict):
    """Run a single source's adapter; never raises. Returns (df, status).

    status ∈ {'ok', 'no_data', 'download_failed', 'needs_manual_path',
    'unavailable'} — so the report can tell "the source has no data here" from
    "the source could not be fetched".
    """
    cfg = OBSERVATION_SOURCES.get(key)
    if cfg is None:
        return _empty_harmonized(), "unavailable"
    fn = _ADAPTERS.get(cfg.get("adapter", "bulk"))
    print(f"\n  → {cfg['name']}")
    ctx["_status"] = "ok"   # adapters overwrite this at their failure points
    try:
        df = fn(key, ctx)
    except Exception as exc:
        print(f"  WARNING: source '{key}' failed ({exc}). Skipping.")
        return _empty_harmonized(), "unavailable"
    if df is None or len(df) == 0:
        st = ctx.get("_status", "ok")
        return _empty_harmonized(), (st if st != "ok" else "no_data")
    return df, "ok"


# ─────────────────────────────────────────────────────────────────────────────
#  Merge + dedup
# ─────────────────────────────────────────────────────────────────────────────
def merge_and_dedup(frames: List[pd.DataFrame], dedup: bool) -> Tuple[pd.DataFrame, int]:
    """Concatenate harmonized frames.

    If dedup, remove CROSS-SOURCE duplicates only: when the same observation
    (rounded lat/lon + parameter + date + value) is reported by MORE THAN ONE
    source, keep the copy from the first-listed source and drop the copies from
    the others. A source's OWN rows are never removed — so adding a source can
    only ADD or keep the total the same, never reduce it (a single source's
    internal duplicates are that source's data, not a cross-source artifact).

    Returns (merged_df, n_cross_source_duplicates_removed).
    """
    frames = [f for f in frames if f is not None and len(f) > 0]
    if not frames:
        return _empty_harmonized(), 0
    merged = pd.concat(frames, ignore_index=True)
    # Nothing to reconcile unless dedup is on AND ≥2 sources actually contributed.
    if not dedup or merged["source"].nunique() < 2:
        return merged, 0

    before = len(merged)
    k = merged.copy()
    k["_lat"] = pd.to_numeric(k["lat"], errors="coerce").round(4)
    k["_lon"] = pd.to_numeric(k["lon"], errors="coerce").round(4)
    key_cols = ["_lat", "_lon", "parameter", "year", "month", "day", "value"]
    # Source priority = order of first appearance (i.e. the selection order).
    src_rank = {s: i for i, s in enumerate(dict.fromkeys(merged["source"]))}
    k["_rank"] = merged["source"].map(src_rank)
    # For each observation key, the best (lowest) source rank present. A row is
    # kept iff it belongs to that winning source — so ALL of the winning
    # source's rows for the key survive (within-source dupes included), and only
    # the redundant copies from other sources are dropped.
    win = k.groupby(key_cols, dropna=False)["_rank"].transform("min")
    keep = (k["_rank"] == win).to_numpy()
    result = merged[keep].reset_index(drop=True)
    return result, before - len(result)


# ─────────────────────────────────────────────────────────────────────────────
#  Outputs: polygon refine, clipped CSVs, stations GeoJSON, stats
# ─────────────────────────────────────────────────────────────────────────────
def _refine_to_polygon(merged: pd.DataFrame, search_area_gdf) -> pd.DataFrame:
    """Best-effort: keep only points inside the (non-rectangular) search polygon.

    Adapters clip by bounding box; this tightens to the actual basin/buffer shape
    when shapely is available. Silently returns the input on any failure.
    """
    if search_area_gdf is None or len(merged) == 0:
        return merged
    try:
        from shapely.geometry import Point
        poly = search_area_gdf.geometry.unary_union
        mask = [poly.contains(Point(xy)) for xy in zip(
            pd.to_numeric(merged["lon"], errors="coerce"),
            pd.to_numeric(merged["lat"], errors="coerce"))]
        return merged[mask].reset_index(drop=True)
    except Exception:
        return merged


# ---------------------------------------------------------------------------
# Where the basin-clipped observations live — the ONE place that knows the
# folder/file names. Everything downstream (config report snippet, calibration
# setup + results reports, the results plot driver) resolves them through
# find_clipped_obs_dir() / load_clipped_observations(), so a layout change
# here never silently breaks a consumer again.
# ---------------------------------------------------------------------------
CLIPPED_OBS_DIRNAME = "obs_clipped_data"            # multi-source layout (current)
LEGACY_CLIPPED_OBS_DIRNAME = "grqa_clipped_data"    # single-GRQA layout (older runs)
MERGED_OBS_FILENAME = "observations_all_sources.csv"
LEGACY_OBS_FILENAME = "grqa_clipped_observations.csv"
LEGACY_STN_FILENAME = "grqa_clipped_stations.csv"


def _has_clipped_files(d: str) -> bool:
    return bool(d) and os.path.isdir(d) and (
        os.path.isfile(os.path.join(d, MERGED_OBS_FILENAME))
        or os.path.isfile(os.path.join(d, LEGACY_OBS_FILENAME)))


def find_clipped_obs_dir(base) -> Optional[str]:
    """Folder holding a run's basin-clipped observations, or None.

    ``base`` may be the run folder (``dir2save_input_files``), its
    ``openwq_in/``, or a clipped folder itself (even one named with the other
    layout's name — the sibling is tried). Current layout first, legacy second.
    """
    if not base:
        return None
    base = os.path.normpath(str(base))
    names = (CLIPPED_OBS_DIRNAME, LEGACY_CLIPPED_OBS_DIRNAME)
    cands = []
    if os.path.basename(base) in names:
        parent = os.path.dirname(base)
        cands += [base] + [os.path.join(parent, n) for n in names]
    for root in (base, os.path.join(base, "openwq_in")):
        cands += [os.path.join(root, n) for n in names]
    for c in cands:
        if _has_clipped_files(c):
            return c
    return None


def load_clipped_observations(base) -> Optional[pd.DataFrame]:
    """The clipped observations of a run as ONE harmonized frame
    (``station_id, lat, lon, parameter, year, month, day, minute, value, units,
    source``), whichever layout is on disk; None when there is none."""
    d = find_clipped_obs_dir(base)
    if not d:
        return None
    merged = os.path.join(d, MERGED_OBS_FILENAME)
    if os.path.isfile(merged):
        df = pd.read_csv(merged)
        df.columns = [str(c).strip() for c in df.columns]
        for c in HARMONIZED_COLUMNS:
            if c not in df.columns:
                df[c] = None
        df = df[HARMONIZED_COLUMNS].copy()
        df["station_id"] = df["station_id"].astype(str)
        df["parameter"] = df["parameter"].astype(str)
        df["value"] = pd.to_numeric(df["value"], errors="coerce")
        for c in ("year", "month", "day", "minute"):
            df[c] = pd.to_numeric(df[c], errors="coerce").fillna(0).astype(int)
        return df
    # legacy single-GRQA layout: station-based CSV (+ optional stations CSV)
    obs = pd.read_csv(os.path.join(d, LEGACY_OBS_FILENAME))
    cols = {str(c).lower(): c for c in obs.columns}
    sid = cols.get("site_id") or cols.get("station_id")
    par = cols.get("model_species") or cols.get("parameter")
    val = cols.get("obs_value") or cols.get("value")
    dat = cols.get("obs_date") or cols.get("datetime") or cols.get("date")
    if not (sid and par and val and dat):
        return None
    out = pd.DataFrame({
        "station_id": obs[sid].astype(str),
        "parameter": obs[par].astype(str),
        "value": pd.to_numeric(obs[val], errors="coerce"),
        "units": obs[cols["unit"]] if "unit" in cols else (obs[cols["units"]] if "units" in cols else None),
    })
    latc = cols.get("lat_wgs84") or cols.get("lat")
    lonc = cols.get("lon_wgs84") or cols.get("lon")
    if latc and lonc:
        out["lat"] = pd.to_numeric(obs[latc], errors="coerce")
        out["lon"] = pd.to_numeric(obs[lonc], errors="coerce")
    else:
        out["lat"] = None
        out["lon"] = None
        stn_p = os.path.join(d, LEGACY_STN_FILENAME)
        if os.path.isfile(stn_p):
            stn = pd.read_csv(stn_p)
            sc = {str(c).lower(): c for c in stn.columns}
            ssid = sc.get("site_id") or sc.get("station_id")
            slat = sc.get("lat_wgs84") or sc.get("lat")
            slon = sc.get("lon_wgs84") or sc.get("lon")
            if ssid and slat and slon:
                lut = stn.drop_duplicates(subset=[ssid]).set_index(stn[ssid].astype(str))
                out["lat"] = out["station_id"].map(pd.to_numeric(lut[slat], errors="coerce"))
                out["lon"] = out["station_id"].map(pd.to_numeric(lut[slon], errors="coerce"))
    parts = _split_datetime(obs[dat])
    for c in parts.columns:
        out[c] = parts[c].values
    return _finalize(out, "grqa")


def write_clipped_csvs(merged: pd.DataFrame, output_dir: str) -> str:
    """Write merged observations (all sources) + per-species files.

    Returns the directory written. Mirrors the harmonized schema so the
    calibration side can ingest it exactly like the user_csv path.
    """
    obs_dir = os.path.join(output_dir, CLIPPED_OBS_DIRNAME)
    os.makedirs(obs_dir, exist_ok=True)
    merged.to_csv(os.path.join(obs_dir, MERGED_OBS_FILENAME), index=False)
    for sp in merged["parameter"].unique():
        sub = merged[merged["parameter"] == sp]
        safe = str(sp).replace("/", "_").replace(" ", "_")
        sub.to_csv(os.path.join(obs_dir, f"{safe}_observations.csv"), index=False)
    print(f"  Wrote {len(merged)} merged observations to {obs_dir}")
    return obs_dir


def build_stations_geojson(merged: pd.DataFrame) -> dict:
    """One GeoJSON point per unique station, with per-species presence + source."""
    features = []
    grp = merged.groupby(["station_id", "source"], dropna=False)
    for (sid, src), sub in grp:
        lat = pd.to_numeric(sub["lat"], errors="coerce").median()
        lon = pd.to_numeric(sub["lon"], errors="coerce").median()
        if pd.isna(lat) or pd.isna(lon):
            continue
        props = {"station_id": str(sid), "source": str(src),
                 "n_observations": int(len(sub)),
                 "parameters": ", ".join(sorted(sub["parameter"].unique()))}
        for sp in sub["parameter"].unique():
            props["has_" + str(sp)] = True
        features.append({"type": "Feature",
                         "geometry": {"type": "Point", "coordinates": [float(lon), float(lat)]},
                         "properties": props})
    return {"type": "FeatureCollection", "features": features}


def build_plot_data(df, *, param_col="parameter", value_col="value",
                    source_col="source", default_source="Observations",
                    date_col=None, year_col="year", month_col="month",
                    day_col="day", cap=4000):
    """Build the per-species conc-vs-time chart data used by the config report:
    ``{species: {source: {"x": [ISO dates], "y": [values]}}}``.

    Works for ANY observation table (multi-source merged, a user CSV, or the
    pre-extracted fallback): give a `date_col`, or year/month/day columns. Rows
    with no date/value are dropped; each species is capped/sampled to `cap`
    points so the HTML stays light. If `source_col` is absent, everything is
    attributed to `default_source`.
    """
    if df is None or len(df) == 0:
        return {}
    if date_col and date_col in df.columns:
        iso = pd.to_datetime(df[date_col], errors="coerce").dt.strftime("%Y-%m-%d")
    else:
        yy = pd.to_numeric(df.get(year_col), errors="coerce")
        mm = pd.to_numeric(df.get(month_col), errors="coerce").fillna(1).clip(1, 12)
        dd = pd.to_numeric(df.get(day_col), errors="coerce").fillna(1).clip(1, 28)
        iso = (yy.round().astype("Int64").astype(str).str.zfill(4) + "-"
               + mm.round().astype(int).astype(str).str.zfill(2) + "-"
               + dd.round().astype(int).astype(str).str.zfill(2))
        iso = iso.where(yy > 0)
    d = pd.DataFrame({
        "_p": df[param_col].astype(str).values,
        "_v": pd.to_numeric(df[value_col], errors="coerce").values,
        "_s": (df[source_col].astype(str).values
               if source_col in df.columns else default_source),
        "_d": iso.values,
    }).dropna(subset=["_v", "_d"])
    out = {}
    for sp in sorted(d["_p"].unique()):
        sp_obs = d[d["_p"] == sp]
        if len(sp_obs) > cap:
            sp_obs = sp_obs.sample(cap, random_state=0)
        by_src = {}
        for src, g in sp_obs.groupby("_s"):
            g = g.sort_values("_d")
            by_src[str(src)] = {"x": g["_d"].tolist(),
                                "y": [round(float(v), 6) for v in g["_v"]]}
        out[str(sp)] = by_src
    return out


def compute_stats(merged, species, n_deduped, output_dir,
                  buffer_km=None, unmapped_species=None, sources_attempted=None):
    """Assemble the stats dict consumed by the setup report.

    `species` is the list of target model-species names; source provenance is
    read straight off the merged frame's `source` column. `sources_attempted`
    is the ordered [(name, extracted_count)] list of EVERY selected source
    (including those that returned 0 rows), so the report can list them all.
    """
    species_stats = []
    found = set()
    for sp in sorted(merged["parameter"].unique()):
        sp_obs = merged[merged["parameter"] == sp]
        yrs = sp_obs["year"].replace(0, pd.NA).dropna()
        species_stats.append({
            "species": sp,
            "n_stations": sp_obs["station_id"].nunique(),
            "n_observations": int(len(sp_obs)),
            "year_start": int(yrs.min()) if len(yrs) else None,
            "year_end": int(yrs.max()) if len(yrs) else None,
            "sources": ", ".join(sorted(sp_obs["source"].unique())),
        })
        found.add(sp)
    model_species = list(species.values()) if isinstance(species, dict) else list(species or [])
    all_years = merged["year"].replace(0, pd.NA).dropna()
    plot_data = build_plot_data(merged)   # per-species conc-vs-time, by source

    return {
        "n_stations": merged["station_id"].nunique(),
        "n_observations": int(len(merged)),
        "year_start": int(all_years.min()) if len(all_years) else None,
        "year_end": int(all_years.max()) if len(all_years) else None,
        "species_stats": species_stats,
        "no_data_species": [s for s in model_species if s not in found],
        "unmapped_species": unmapped_species or [],
        "sources_used": sorted(merged["source"].dropna().astype(str).unique().tolist()),
        "per_source_counts": {str(k): int(v) for k, v
                              in merged.groupby("source").size().items()},
        "sources_attempted": [[str(n), int(c), str(s)]
                              for n, c, s in (sources_attempted or [])],
        "plot_data": plot_data,
        "n_deduplicated": int(n_deduped),
        "buffer_km": buffer_km,
        "output_dir": output_dir,
    }


def load_user_csv_harmonized(path: str) -> pd.DataFrame:
    """Read a user-provided observation CSV (already the 10-column schema) into
    the harmonized frame, tagged source='User CSV'. Case-insensitive headers."""
    if not path or not os.path.isfile(path):
        return _empty_harmonized()
    try:
        df = pd.read_csv(path)
    except Exception as exc:
        print(f"  WARNING: could not read user CSV '{path}' ({exc}).")
        return _empty_harmonized()
    df.columns = [str(c).strip().lower() for c in df.columns]
    return _finalize(df, "User CSV")


def extract_observations(selected, *, search_area_gdf, bbox, species,
                         output_dir, grqa_species_mapping=None, years=None,
                         cache_dir=None, manual_paths=None, extra_frames=None,
                         interactive=True, dedup_mode="ask", buffer_m=100000,
                         grqa_local_data_path=None, buffer_km=None,
                         unmapped_species=None):
    """Resolve the selected sources, extract + merge them, and return the
    (stations_geojson_str_or_None, stats_dict) the setup report consumes.

    `species` = target model-species names (list). `grqa_species_mapping` =
    {grqa_code: model_name} used only by the GRQA adapter. `extra_frames` =
    already-harmonized frames (e.g. a user CSV) folded into the merge + dedup.
    """
    resolved, do_dedup = resolve_sources(selected, interactive, dedup_mode)
    data_sources = [s for s in resolved if s in OBSERVATION_SOURCES]
    extra_frames = [f for f in (extra_frames or []) if f is not None and len(f) > 0]
    empty_stats = {
        "n_stations": 0, "n_observations": 0, "year_start": None, "year_end": None,
        "species_stats": [], "no_data_species": list(species or []),
        "unmapped_species": unmapped_species or [],
        "sources_used": [], "per_source_counts": {}, "sources_attempted": [],
        "plot_data": {}, "n_deduplicated": 0, "buffer_km": buffer_km,
        "output_dir": output_dir,
    }
    if not data_sources and not extra_frames:
        return None, empty_stats

    cache_dir = cache_dir or os.path.join(output_dir, "obs_cache")
    ctx = {
        "search_area_gdf": search_area_gdf, "bbox": bbox,
        "species": list(species or []), "grqa_species_mapping": grqa_species_mapping,
        "years": years, "cache_dir": cache_dir, "output_dir": output_dir,
        "manual_paths": manual_paths, "buffer_m": buffer_m,
        "grqa_local_data_path": grqa_local_data_path,
    }
    frames = list(extra_frames)
    attempted = []   # (display_name, extracted_row_count, status) in selection order
    if extra_frames:
        _uc = int(sum(len(f) for f in extra_frames))
        attempted.append(("User CSV", _uc, "ok" if _uc else "no_data"))
    for key in data_sources:
        df, status = extract_one(key, ctx)
        _nm = OBSERVATION_SOURCES[key]["name"].split("—")[0].strip()
        attempted.append((_nm, int(len(df)), status))
        frames.append(df)

    merged, n_dropped = merge_and_dedup(frames, do_dedup)
    merged = _refine_to_polygon(merged, search_area_gdf)
    if len(merged) == 0:
        print("\n  No observations found in the search area across the selected sources.")
        empty_stats["sources_attempted"] = attempted
        return None, empty_stats

    write_clipped_csvs(merged, output_dir)
    geojson = build_stations_geojson(merged)
    stats = compute_stats(merged, species, n_dropped, output_dir,
                          buffer_km, unmapped_species, attempted)
    print(f"\n  Observations: {stats['n_observations']} at {stats['n_stations']} "
          f"stations from {len(stats['sources_used'])} source(s)"
          + (f"; removed {n_dropped} duplicate rows" if n_dropped else "") + ".")
    return json.dumps(geojson), stats

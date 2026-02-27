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
OpenWQ ML-Based Source/Sink Generator

Trains an XGBoost or Random Forest model from water quality monitoring data,
then generates:
  1. An OpenWQ source/sink JSON configuration from predictions
  2. A serialized model file (XGBoost text format) for optional C++ runtime
     inference using FastForest

The ML approach is based on the EPA XGBest methodology, which demonstrated
that tree-based models outperform LOADEST and WRTDS for daily nutrient
load prediction when integrating hydrological and landscape data.

References:
    - EPA XGBest: https://cfpub.epa.gov/si/si_public_record_Report.cfm?dirEntryId=364132
    - FastForest (C++ XGBoost deployment): https://github.com/guitargeek/XGBoost-FastForest
"""

import json
import re
from pathlib import Path
from typing import List, Dict, Union, Optional

import numpy as np
import pandas as pd


def _check_ml_dependencies():
    """Check that ML dependencies are available."""
    missing = []
    try:
        import sklearn  # noqa: F401
    except ImportError:
        missing.append("scikit-learn")
    try:
        import xgboost  # noqa: F401
    except ImportError:
        missing.append("xgboost")

    if missing:
        raise ImportError(
            f"ML source/sink method requires: {', '.join(missing)}\n"
            f"Install with: pip install {' '.join(missing)}"
        )


def _prepare_features(df: pd.DataFrame,
                      feature_columns: Optional[List[str]],
                      target_columns: List[str]) -> tuple:
    """
    Prepare feature matrix and identify feature columns.

    If feature_columns is None, auto-detect all numeric columns
    that are not in target_columns or 'date'.

    Also engineers temporal features from the 'date' column:
    - month (1-12)
    - day_of_year (1-366)
    - sin_month / cos_month (cyclical encoding)

    Returns:
        (feature_columns, X DataFrame)
    """
    # Engineer temporal features if 'date' column exists
    if 'date' in df.columns:
        df['date'] = pd.to_datetime(df['date'])
        df['month'] = df['date'].dt.month
        df['day_of_year'] = df['date'].dt.dayofyear
        df['sin_month'] = np.sin(2 * np.pi * df['month'] / 12)
        df['cos_month'] = np.cos(2 * np.pi * df['month'] / 12)

    if feature_columns is None:
        # Auto-detect: all numeric columns except targets and date
        exclude = set(target_columns) | {'date'}
        feature_columns = [
            col for col in df.select_dtypes(include=[np.number]).columns
            if col not in exclude
        ]
        print(f"  Auto-detected feature columns: {feature_columns}")

    # Ensure all feature columns exist
    missing_feats = [f for f in feature_columns if f not in df.columns]
    if missing_feats:
        raise ValueError(
            f"Feature columns not found in CSV: {missing_feats}\n"
            f"Available columns: {list(df.columns)}"
        )

    return feature_columns, df[feature_columns]


def train_ml_model(
        training_data_csv: str,
        model_type: str = "xgboost",
        target_species: Optional[List[str]] = None,
        feature_columns: Optional[List[str]] = None,
        n_estimators: int = 200,
        max_depth: int = 6,
        output_dir: Optional[Path] = None
) -> tuple:
    """
    Train an ML model for nutrient concentration/load prediction.

    Parameters:
    -----------
    training_data_csv : str
        Path to CSV with monitoring data. Must contain:
        - 'date' column (datetime, for temporal feature engineering)
        - Feature columns (e.g., discharge_m3s, precip_mm, temp_c)
        - Target columns (nutrient concentrations, e.g., NO3-N, TP)
    model_type : str
        "xgboost" or "random_forest"
    target_species : list of str, optional
        Column names for target variables. If None, auto-detected as
        columns not in common feature names.
    feature_columns : list of str, optional
        Feature column names. If None, auto-detected.
    n_estimators : int
        Number of trees
    max_depth : int
        Maximum tree depth
    output_dir : Path, optional
        Directory to save model files and reports

    Returns:
    --------
    tuple : (models_dict, feature_columns, metrics_dict)
        - models_dict: {species_name: trained_model}
        - feature_columns: list of feature names used
        - metrics_dict: {species_name: {'r2': float, 'rmse': float, 'mae': float}}
    """
    _check_ml_dependencies()
    from sklearn.model_selection import train_test_split
    from sklearn.metrics import r2_score, mean_squared_error, mean_absolute_error

    print("\n" + "=" * 60)
    print("TRAINING ML MODEL FOR NUTRIENT PREDICTION")
    print("=" * 60)
    print(f"Model type: {model_type}")
    print(f"Training data: {training_data_csv}")

    # Load data
    df = pd.read_csv(training_data_csv)
    print(f"Loaded {len(df)} records with columns: {list(df.columns)}")

    # Auto-detect target species if not provided
    common_features = {
        'date', 'discharge_m3s', 'precip_mm', 'temp_c', 'streamflow',
        'precipitation', 'temperature', 'month', 'day_of_year',
        'sin_month', 'cos_month', 'area_km2', 'slope', 'elevation'
    }
    if target_species is None:
        target_species = [
            col for col in df.select_dtypes(include=[np.number]).columns
            if col.lower() not in {f.lower() for f in common_features}
        ]
        print(f"  Auto-detected target species: {target_species}")

    if not target_species:
        raise ValueError(
            "No target species found. Provide target_species parameter or ensure "
            "CSV contains nutrient concentration columns."
        )

    # Prepare features
    feature_columns, X = _prepare_features(df, feature_columns, target_species)

    # Drop rows with NaN in features or targets
    valid_mask = X.notna().all(axis=1)
    for sp in target_species:
        valid_mask &= df[sp].notna()

    X = X[valid_mask].copy()
    df_valid = df[valid_mask].copy()

    print(f"  Valid records after dropping NaN: {len(X)}")

    # Train/test split
    X_train, X_test, idx_train, idx_test = train_test_split(
        X, X.index, test_size=0.2, random_state=42
    )

    models = {}
    metrics = {}

    for species in target_species:
        print(f"\n  --- Training model for: {species} ---")

        y_train = df_valid.loc[idx_train, species].values
        y_test = df_valid.loc[idx_test, species].values

        # Remove negative values (concentrations should be >= 0)
        y_train = np.maximum(y_train, 0)
        y_test = np.maximum(y_test, 0)

        if model_type == "xgboost":
            import xgboost as xgb
            model = xgb.XGBRegressor(
                n_estimators=n_estimators,
                max_depth=max_depth,
                learning_rate=0.1,
                subsample=0.8,
                colsample_bytree=0.8,
                random_state=42,
                n_jobs=-1
            )
        elif model_type == "random_forest":
            from sklearn.ensemble import RandomForestRegressor
            model = RandomForestRegressor(
                n_estimators=n_estimators,
                max_depth=max_depth,
                random_state=42,
                n_jobs=-1
            )
        else:
            raise ValueError(
                f"Unknown model_type: '{model_type}'. Options: 'xgboost', 'random_forest'"
            )

        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        y_pred = np.maximum(y_pred, 0)  # Clamp to non-negative

        # Compute metrics
        r2 = r2_score(y_test, y_pred)
        rmse = np.sqrt(mean_squared_error(y_test, y_pred))
        mae = mean_absolute_error(y_test, y_pred)

        metrics[species] = {'r2': r2, 'rmse': rmse, 'mae': mae}
        models[species] = model

        print(f"    R2:   {r2:.4f}")
        print(f"    RMSE: {rmse:.4f}")
        print(f"    MAE:  {mae:.4f}")

        # Feature importance
        if hasattr(model, 'feature_importances_'):
            importances = model.feature_importances_
            sorted_idx = np.argsort(importances)[::-1]
            print(f"    Top features:")
            for i in range(min(5, len(sorted_idx))):
                fi = sorted_idx[i]
                print(f"      {feature_columns[fi]}: {importances[fi]:.4f}")

    # Save models if output_dir provided
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        for species, model in models.items():
            safe_name = species.replace("/", "_").replace("-", "_")

            if model_type == "xgboost":
                # Save in XGBoost text format (compatible with FastForest C++ library)
                model_path = output_dir / f"ml_model_{safe_name}.txt"
                model.get_booster().dump_model(str(model_path))
                print(f"\n  Saved XGBoost model (FastForest-compatible): {model_path.name}")

                # Also save in binary format for Python reuse
                bin_path = output_dir / f"ml_model_{safe_name}.json"
                model.save_model(str(bin_path))
                print(f"  Saved XGBoost binary model: {bin_path.name}")

            elif model_type == "random_forest":
                import joblib
                model_path = output_dir / f"ml_model_{safe_name}.joblib"
                joblib.dump(model, str(model_path))
                print(f"\n  Saved Random Forest model: {model_path.name}")

        # Save feature list
        features_path = output_dir / "ml_model_features.json"
        with open(features_path, 'w') as f:
            json.dump({
                'feature_columns': feature_columns,
                'target_species': target_species,
                'model_type': model_type,
                'n_estimators': n_estimators,
                'max_depth': max_depth,
                'metrics': {sp: {k: float(v) for k, v in m.items()}
                           for sp, m in metrics.items()}
            }, f, indent=2)
        print(f"  Saved feature config: {features_path.name}")

    return models, feature_columns, metrics


def generate_predictions_json(
        models: dict,
        training_data_csv: str,
        feature_columns: List[str],
        ss_config_filepath: str,
        json_header_comment: List[str],
        compartment_name: str,
        ss_metadata_comment: str = "ML-predicted nutrient loading",
        ss_metadata_source: str = "XGBoost/RF model trained on monitoring data"
) -> None:
    """
    Generate OpenWQ source/sink JSON from ML model predictions.

    Uses the trained models to predict concentrations for the entire
    training dataset, then converts to loads using discharge.

    Parameters:
    -----------
    models : dict
        {species_name: trained_model}
    training_data_csv : str
        Path to the original training data CSV
    feature_columns : list of str
        Feature column names
    ss_config_filepath : str
        Output JSON path
    json_header_comment : list of str
        Comment header lines
    compartment_name : str
        OpenWQ compartment name
    ss_metadata_comment : str
        Metadata comment
    ss_metadata_source : str
        Metadata source
    """
    print("\n" + "=" * 60)
    print("GENERATING SS JSON FROM ML PREDICTIONS")
    print("=" * 60)

    df = pd.read_csv(training_data_csv)

    # Engineer temporal features
    if 'date' in df.columns:
        df['date'] = pd.to_datetime(df['date'])
        df['month'] = df['date'].dt.month
        df['day_of_year'] = df['date'].dt.dayofyear
        df['sin_month'] = np.sin(2 * np.pi * df['month'] / 12)
        df['cos_month'] = np.cos(2 * np.pi * df['month'] / 12)

    X = df[feature_columns]

    config = {
        "METADATA": {
            "Comment": ss_metadata_comment,
            "Source": ss_metadata_source
        }
    }

    entry_idx = 1

    for species, model in models.items():
        print(f"\n  Predicting: {species}")

        predictions = model.predict(X)
        predictions = np.maximum(predictions, 0)

        # Convert concentration to load if discharge is available
        # Load (kg/s) = concentration (mg/L) * discharge (m3/s) * 1e-6
        if 'discharge_m3s' in df.columns:
            loads_kg_s = predictions * df['discharge_m3s'].values * 1e-6
            units = "kg"
            # Aggregate to monthly loads
            df_pred = df[['date']].copy()
            df_pred['load_kg_s'] = loads_kg_s
            df_pred['year'] = df_pred['date'].dt.year
            df_pred['month'] = df_pred['date'].dt.month

            # Sum daily loads to monthly (kg/s * 86400 s/day * days)
            monthly = df_pred.groupby(['year', 'month']).agg(
                load_sum=('load_kg_s', lambda x: x.sum() * 86400)
            ).reset_index()
        else:
            # No discharge: use predictions as concentration source
            units = "mg/l"
            df_pred = df[['date']].copy()
            df_pred['conc'] = predictions
            df_pred['year'] = df_pred['date'].dt.year
            df_pred['month'] = df_pred['date'].dt.month

            monthly = df_pred.groupby(['year', 'month']).agg(
                load_sum=('conc', 'mean')
            ).reset_index()

        # Build JSON data entries
        data_entries = {}
        sub_idx = 1

        for _, row in monthly.iterrows():
            data_entries[str(sub_idx)] = [
                int(row['year']), int(row['month']), 1, 1, "all", "all",
                "all", "all", "all",
                float(row['load_sum']),
                "continuous" if units == "mg/l" else "discrete"
            ]
            sub_idx += 1

        if data_entries:
            config[str(entry_idx)] = {
                "CHEMICAL_NAME": species,
                "COMPARTMENT_NAME": compartment_name,
                "COMMENT": f"ML-predicted {species} for {compartment_name}",
                "TYPE": "source",
                "UNITS": units,
                "DATA_FORMAT": "JSON",
                "DATA": data_entries
            }
            print(f"    Added {len(data_entries)} monthly prediction entries")
            entry_idx += 1

    # Write JSON
    json_string = json.dumps(config, indent=4)

    def compress_array(match):
        array_content = match.group(0)
        compressed = re.sub(r'\s+', ' ', array_content)
        compressed = re.sub(r'\[\s+', '[', compressed)
        compressed = re.sub(r'\s+\]', ']', compressed)
        compressed = re.sub(r'\s*,\s*', ', ', compressed)
        return compressed

    json_string = re.sub(r'\[[^\[\]]*\]', compress_array, json_string)

    output_path = Path(ss_config_filepath)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(ss_config_filepath, 'w') as f:
        for comment in json_header_comment:
            f.write(comment + "\n")
        f.write(json_string)
        f.write("\n")

    print(f"\n  SS JSON saved to: {ss_config_filepath}")


def train_and_generate_ss_json(
        ss_config_filepath: str,
        json_header_comment: List[str],
        training_data_csv: str,
        model_type: str = "xgboost",
        target_species: Optional[List[str]] = None,
        feature_columns: Optional[List[str]] = None,
        n_estimators: int = 200,
        max_depth: int = 6,
        ss_metadata_comment: str = "ML-predicted nutrient loading",
        ss_metadata_source: str = "XGBoost/RF model trained on monitoring data",
        compartment_name: str = "RIVER_NETWORK_REACHES"
) -> None:
    """
    Main entry point: train ML model and generate SS JSON.

    This is the function called by Gen_Input_Driver when ss_method = "ml_model".

    Steps:
    1. Train XGBoost/RF model from monitoring CSV
    2. Save trained model (XGBoost text format for FastForest C++ deployment)
    3. Generate SS JSON from model predictions
    4. Save model metrics and feature configuration

    Parameters:
    -----------
    ss_config_filepath : str
        Path for output SS JSON file
    json_header_comment : list of str
        Header comments for JSON file
    training_data_csv : str
        Path to monitoring data CSV
    model_type : str
        "xgboost" or "random_forest"
    target_species : list of str, optional
        Target nutrient column names in CSV
    feature_columns : list of str, optional
        Feature column names (auto-detected if None)
    n_estimators : int
        Number of trees
    max_depth : int
        Max tree depth
    ss_metadata_comment : str
        Metadata comment for JSON
    ss_metadata_source : str
        Metadata source for JSON
    compartment_name : str
        OpenWQ compartment name
    """
    # Output directory for model files
    output_dir = Path(ss_config_filepath).parent / 'ss_ml_model_files'
    output_dir.mkdir(parents=True, exist_ok=True)

    # Step 1: Train model
    models, feat_cols, metrics = train_ml_model(
        training_data_csv=training_data_csv,
        model_type=model_type,
        target_species=target_species,
        feature_columns=feature_columns,
        n_estimators=n_estimators,
        max_depth=max_depth,
        output_dir=output_dir
    )

    # Step 2: Generate SS JSON from predictions
    generate_predictions_json(
        models=models,
        training_data_csv=training_data_csv,
        feature_columns=feat_cols,
        ss_config_filepath=ss_config_filepath,
        json_header_comment=json_header_comment,
        compartment_name=compartment_name,
        ss_metadata_comment=ss_metadata_comment,
        ss_metadata_source=ss_metadata_source
    )

    # Print summary
    print("\n" + "=" * 60)
    print("ML SOURCE/SINK GENERATION COMPLETE!")
    print("=" * 60)
    print(f"SS JSON: {ss_config_filepath}")
    print(f"Model files: {output_dir}/")
    print(f"\nModel performance:")
    for species, m in metrics.items():
        print(f"  {species}: R2={m['r2']:.3f}, RMSE={m['rmse']:.3f}, MAE={m['mae']:.3f}")

    if model_type == "xgboost":
        print(f"\nFor C++ runtime inference:")
        print(f"  Model files saved in XGBoost text format (FastForest-compatible)")
        print(f"  See: https://github.com/guitargeek/XGBoost-FastForest")
        print(f"  Feature config: {output_dir / 'ml_model_features.json'}")

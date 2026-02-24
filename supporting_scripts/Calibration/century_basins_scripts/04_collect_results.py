#!/usr/bin/env python3
# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of OpenWQ model.
#
# This program, openWQ, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.

"""
04_collect_results.py
=====================

Aggregate calibration results across all basins and variants, and generate
cross-variant comparison analyses.

Outputs:
  1. Per-variant master CSV with basin ID, KGE_cal, KGE_eval, all params
  2. Cross-variant comparison:
     - KGE distributions (A vs B vs C) by basin scale
     - Complexity vs performance scatter
     - Equifinality analysis
     - Parameter identifiability per variant
  3. Spatial maps (KGE by basin for each variant)
  4. Summary recommendation report

USAGE:
  python 04_collect_results.py \\
      --prepared-dir /path/to/century_basins_prepared \\
      --output-dir /path/to/results \\
      [--variants A,B,C,D] \\
      [--plot]
"""

import sys
import os
import argparse
import json
import glob
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import pandas as pd
    import numpy as np
except ImportError:
    print("ERROR: pandas and numpy required.")
    sys.exit(1)

SCRIPT_DIR = Path(__file__).resolve().parent

# Variant metadata
VARIANT_META = {
    'A': {
        'label': 'Full N Cycle',
        'n_bgc_params': 11,
        'n_ss_params': 7,
        'max_evals': 500,
        'color': '#1f77b4',
    },
    'B': {
        'label': 'SWAT Full Nutrients',
        'n_bgc_params': 14,
        'n_ss_params': 7,
        'max_evals': 600,
        'color': '#ff7f0e',
    },
    'C': {
        'label': 'Thermodynamic N Cycle',
        'n_bgc_params': 9,
        'n_ss_params': 7,
        'max_evals': 450,
        'color': '#2ca02c',
    },
}


# =============================================================================
# RESULT COLLECTION
# =============================================================================

def collect_variant_results(
    prepared_dir: str,
    variant_key: str,
) -> pd.DataFrame:
    """
    Collect calibration results for one variant across all basins.

    Looks for:
      basin_*/calibration_workspace_{variant_key}/results/calibration_results.json
      basin_*/calibration_workspace_{variant_key}/results/best_parameters.json
      basin_*/calibration_workspace_{variant_key}/calibration_state.json (checkpoint)

    Returns DataFrame with columns:
      basin_id, country, scale, station_id,
      kge_best, kge_cal, kge_eval,
      n_evaluations, converged,
      param1, param2, ..., paramN
    """
    rows = []
    basin_dirs = sorted(glob.glob(os.path.join(prepared_dir, "basin_*")))

    for basin_dir in basin_dirs:
        basin_id = Path(basin_dir).name.replace('basin_', '')
        parts = basin_id.split('_')
        country = parts[0] if len(parts) >= 1 else 'UNKNOWN'
        station_id = parts[1] if len(parts) >= 2 else 'UNKNOWN'
        scale = parts[2] if len(parts) >= 3 else 'unknown'

        row = {
            'basin_id': basin_id,
            'country': country,
            'station_id': station_id,
            'scale': scale,
            'variant': variant_key,
            'status': 'missing',
        }

        # Try to find results
        workspace_dir = os.path.join(basin_dir, f"calibration_workspace_{variant_key}")
        results_dir = os.path.join(workspace_dir, "results")

        # Check calibration results JSON
        results_file = os.path.join(results_dir, "calibration_results.json")
        if os.path.exists(results_file):
            try:
                with open(results_file, 'r') as f:
                    results = json.load(f)

                row['status'] = 'completed'
                row['kge_best'] = results.get('best_objective', np.nan)
                row['n_evaluations'] = results.get('n_evaluations', 0)
                row['converged'] = results.get('converged', False)

                # Best parameters
                best_params = results.get('best_params', {})
                for pname, pval in best_params.items():
                    row[f'param_{pname}'] = pval

                # Per-species metrics
                species_metrics = results.get('species_metrics', {})
                for sp, metrics in species_metrics.items():
                    for metric_name, metric_val in metrics.items():
                        row[f'{sp}_{metric_name}'] = metric_val

            except Exception as e:
                row['status'] = 'error'
                row['error'] = str(e)
        else:
            # Check checkpoint for in-progress runs
            checkpoint_file = os.path.join(workspace_dir, "calibration_state.json")
            if os.path.exists(checkpoint_file):
                try:
                    with open(checkpoint_file, 'r') as f:
                        state = json.load(f)
                    row['status'] = 'in_progress'
                    row['n_evaluations'] = state.get('current_evaluation', 0)
                    row['kge_best'] = state.get('best_objective', np.nan)
                except Exception:
                    row['status'] = 'checkpoint_error'
            else:
                row['status'] = 'not_started'

        # Check sensitivity analysis results
        sa_file = os.path.join(results_dir, "sensitivity_results.json")
        if os.path.exists(sa_file):
            try:
                with open(sa_file, 'r') as f:
                    sa = json.load(f)
                row['sa_influential_params'] = len(sa.get('influential_params', []))
                row['sa_total_params'] = sa.get('total_params', 0)
            except Exception:
                pass

        rows.append(row)

    return pd.DataFrame(rows)


def compute_cross_variant_stats(
    results: Dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """
    Compute cross-variant comparison statistics.

    For each basin, compare KGE across variants.
    """
    # Merge all variants on basin_id
    all_variants = []
    for v_key, df in results.items():
        subset = df[['basin_id', 'kge_best', 'n_evaluations', 'status']].copy()
        subset = subset.rename(columns={
            'kge_best': f'kge_{v_key}',
            'n_evaluations': f'n_eval_{v_key}',
            'status': f'status_{v_key}',
        })
        all_variants.append(subset)

    if not all_variants:
        return pd.DataFrame()

    merged = all_variants[0]
    for df in all_variants[1:]:
        merged = merged.merge(df, on='basin_id', how='outer')

    # Compute comparison metrics
    kge_cols = [c for c in merged.columns if c.startswith('kge_')]
    if len(kge_cols) >= 2:
        merged['best_variant'] = merged[kge_cols].idxmax(axis=1).str.replace('kge_', '')
        merged['kge_range'] = merged[kge_cols].max(axis=1) - merged[kge_cols].min(axis=1)
        merged['kge_mean'] = merged[kge_cols].mean(axis=1)

    return merged


# =============================================================================
# REPORTING
# =============================================================================

def generate_summary_report(
    results: Dict[str, pd.DataFrame],
    cross_variant: pd.DataFrame,
    output_dir: str,
):
    """Generate comprehensive text summary report."""
    report_path = os.path.join(output_dir, "calibration_summary_report.txt")

    with open(report_path, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("CENTURY BASINS CALIBRATION — CROSS-VARIANT COMPARISON REPORT\n")
        f.write("=" * 80 + "\n\n")

        # Per-variant summary
        for v_key, df in results.items():
            meta = VARIANT_META[v_key]
            completed = df[df['status'] == 'completed']

            f.write(f"\n{'─' * 80}\n")
            f.write(f"VARIANT {v_key}: {meta['label']}\n")
            f.write(f"{'─' * 80}\n")
            f.write(f"  BGC parameters:    {meta['n_bgc_params']}\n")
            f.write(f"  SS parameters:     {meta['n_ss_params']}\n")
            f.write(f"  Total parameters:  {meta['n_bgc_params'] + meta['n_ss_params']}\n")
            f.write(f"  Max evaluations:   {meta['max_evals']}\n\n")

            f.write(f"  Status:\n")
            for status, count in df['status'].value_counts().items():
                f.write(f"    {status:20s}: {count:3d}\n")

            if not completed.empty and 'kge_best' in completed.columns:
                kge = completed['kge_best'].dropna()
                if not kge.empty:
                    f.write(f"\n  KGE Statistics (completed basins):\n")
                    f.write(f"    Mean:     {kge.mean():.3f}\n")
                    f.write(f"    Median:   {kge.median():.3f}\n")
                    f.write(f"    Std:      {kge.std():.3f}\n")
                    f.write(f"    Min:      {kge.min():.3f}\n")
                    f.write(f"    Max:      {kge.max():.3f}\n")
                    f.write(f"    >0.5:     {(kge > 0.5).sum()}/{len(kge)} "
                           f"({100*(kge > 0.5).mean():.0f}%)\n")
                    f.write(f"    >0.7:     {(kge > 0.7).sum()}/{len(kge)} "
                           f"({100*(kge > 0.7).mean():.0f}%)\n")

                    # By scale
                    f.write(f"\n  KGE by basin scale:\n")
                    for scale, group in completed.groupby('scale'):
                        scale_kge = group['kge_best'].dropna()
                        if not scale_kge.empty:
                            f.write(f"    {scale:12s}: {scale_kge.mean():.3f} "
                                   f"(n={len(scale_kge)}, "
                                   f"range=[{scale_kge.min():.3f}, {scale_kge.max():.3f}])\n")

        # Cross-variant comparison
        if not cross_variant.empty:
            f.write(f"\n{'═' * 80}\n")
            f.write(f"CROSS-VARIANT COMPARISON\n")
            f.write(f"{'═' * 80}\n")

            kge_cols = [c for c in cross_variant.columns if c.startswith('kge_') and len(c) == 5]
            completed_all = cross_variant.dropna(subset=kge_cols)

            if not completed_all.empty:
                f.write(f"\nBasins with all variants completed: {len(completed_all)}\n\n")

                # Best variant distribution
                if 'best_variant' in completed_all.columns:
                    f.write(f"Best variant distribution:\n")
                    for v, count in completed_all['best_variant'].value_counts().items():
                        pct = 100 * count / len(completed_all)
                        f.write(f"  Variant {v}: {count:3d} basins ({pct:.0f}%)\n")

                # KGE range (equifinality indicator)
                if 'kge_range' in completed_all.columns:
                    kge_range = completed_all['kge_range'].dropna()
                    f.write(f"\nEquifinality analysis (KGE range across variants):\n")
                    f.write(f"  Mean range:   {kge_range.mean():.3f}\n")
                    f.write(f"  Median range: {kge_range.median():.3f}\n")
                    f.write(f"  Range < 0.05: {(kge_range < 0.05).sum()} basins "
                           f"(variants nearly equivalent)\n")
                    f.write(f"  Range > 0.2:  {(kge_range > 0.2).sum()} basins "
                           f"(significant variant differences)\n")

                # Complexity vs performance
                f.write(f"\nComplexity vs Performance:\n")
                for v_key in sorted(VARIANT_META.keys()):
                    meta = VARIANT_META[v_key]
                    kge_col = f'kge_{v_key}'
                    if kge_col in completed_all.columns:
                        mean_kge = completed_all[kge_col].mean()
                        total_params = meta['n_bgc_params'] + meta['n_ss_params']
                        f.write(f"  [{v_key}] {meta['label']:25s}: "
                               f"params={total_params:2d}, "
                               f"mean KGE={mean_kge:.3f}, "
                               f"KGE/param={mean_kge/total_params:.4f}\n")

        f.write(f"\n{'═' * 80}\n")
        f.write(f"Report generated: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

    print(f"Summary report: {report_path}")
    return report_path


def generate_plots(
    results: Dict[str, pd.DataFrame],
    cross_variant: pd.DataFrame,
    output_dir: str,
):
    """Generate comparison plots (requires matplotlib)."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from matplotlib.gridspec import GridSpec
    except ImportError:
        print("WARNING: matplotlib not available, skipping plots.")
        return

    plots_dir = os.path.join(output_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)

    # 1. KGE Boxplots by variant
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)

    for idx, (v_key, df) in enumerate(sorted(results.items())):
        meta = VARIANT_META[v_key]
        completed = df[df['status'] == 'completed']
        if completed.empty or 'kge_best' not in completed.columns:
            continue

        kge = completed['kge_best'].dropna()

        # Overall boxplot
        bp = axes[idx].boxplot([kge.values], labels=[f'All\n(n={len(kge)})'],
                               patch_artist=True, widths=0.4)
        bp['boxes'][0].set_facecolor(meta['color'])

        # By scale
        scales_data = []
        scales_labels = []
        for scale in ['headwater', 'meso', 'macro']:
            scale_kge = completed[completed['scale'] == scale]['kge_best'].dropna()
            if not scale_kge.empty:
                scales_data.append(scale_kge.values)
                scales_labels.append(f'{scale}\n(n={len(scale_kge)})')

        if scales_data:
            positions = list(range(2, 2 + len(scales_data)))
            bp2 = axes[idx].boxplot(scales_data, labels=scales_labels,
                                     positions=positions,
                                     patch_artist=True, widths=0.4)
            for box in bp2['boxes']:
                box.set_facecolor(meta['color'])
                box.set_alpha(0.6)

        axes[idx].set_title(f'Variant {v_key}: {meta["label"]}')
        axes[idx].axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
        axes[idx].axhline(y=0.7, color='gray', linestyle=':', alpha=0.5)
        axes[idx].set_ylabel('KGE' if idx == 0 else '')
        axes[idx].set_ylim(-0.5, 1.0)

    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "kge_boxplots_by_variant.png"), dpi=150)
    plt.close()
    print(f"  Plot: kge_boxplots_by_variant.png")

    # 2. Cross-variant scatter (pairwise)
    kge_cols = [c for c in cross_variant.columns if c.startswith('kge_') and len(c) == 5]
    if len(kge_cols) >= 2:
        n_pairs = len(kge_cols) * (len(kge_cols) - 1) // 2
        fig, axes = plt.subplots(1, n_pairs, figsize=(5*n_pairs, 5))
        if n_pairs == 1:
            axes = [axes]

        pair_idx = 0
        for i in range(len(kge_cols)):
            for j in range(i+1, len(kge_cols)):
                ax = axes[pair_idx]
                valid = cross_variant.dropna(subset=[kge_cols[i], kge_cols[j]])
                if not valid.empty:
                    # Color by scale
                    for scale, marker in [('headwater', 'o'), ('meso', 's'), ('macro', '^')]:
                        parts = valid['basin_id'].str.split('_')
                        scale_mask = parts.str[-1] == scale if not parts.empty else pd.Series(False, index=valid.index)
                        subset = valid[scale_mask]
                        if not subset.empty:
                            ax.scatter(subset[kge_cols[i]], subset[kge_cols[j]],
                                      marker=marker, alpha=0.6, label=scale, s=30)

                    ax.plot([-.5, 1], [-.5, 1], 'k--', alpha=0.3)
                    ax.set_xlabel(f'Variant {kge_cols[i][-1]}')
                    ax.set_ylabel(f'Variant {kge_cols[j][-1]}')
                    ax.set_xlim(-0.5, 1)
                    ax.set_ylim(-0.5, 1)
                    ax.legend(fontsize=8)
                    ax.set_aspect('equal')
                pair_idx += 1

        plt.tight_layout()
        plt.savefig(os.path.join(plots_dir, "cross_variant_scatter.png"), dpi=150)
        plt.close()
        print(f"  Plot: cross_variant_scatter.png")

    # 3. Complexity vs Performance
    fig, ax = plt.subplots(figsize=(8, 5))
    for v_key in sorted(VARIANT_META.keys()):
        meta = VARIANT_META[v_key]
        if v_key not in results:
            continue
        completed = results[v_key][results[v_key]['status'] == 'completed']
        if completed.empty or 'kge_best' not in completed.columns:
            continue

        total_params = meta['n_bgc_params'] + meta['n_ss_params']
        kge_vals = completed['kge_best'].dropna()
        if not kge_vals.empty:
            ax.boxplot([kge_vals.values], positions=[total_params],
                      widths=1.5, patch_artist=True,
                      boxprops=dict(facecolor=meta['color'], alpha=0.7),
                      medianprops=dict(color='black'),
                      labels=[f'{v_key}\n({total_params})'])

    ax.set_xlabel('Number of Calibration Parameters')
    ax.set_ylabel('KGE')
    ax.set_title('Complexity vs Performance')
    ax.axhline(y=0.5, color='gray', linestyle='--', alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, "complexity_vs_performance.png"), dpi=150)
    plt.close()
    print(f"  Plot: complexity_vs_performance.png")

    # 4. Parameter identifiability (coefficient of variation across basins)
    for v_key, df in results.items():
        completed = df[df['status'] == 'completed']
        param_cols = [c for c in completed.columns if c.startswith('param_')]

        if not param_cols or completed.empty:
            continue

        cv_data = {}
        for col in param_cols:
            vals = completed[col].dropna()
            if len(vals) >= 5:
                mean = vals.mean()
                if abs(mean) > 1e-10:
                    cv_data[col.replace('param_', '')] = vals.std() / abs(mean)

        if cv_data:
            fig, ax = plt.subplots(figsize=(10, max(4, len(cv_data) * 0.3)))
            sorted_cv = sorted(cv_data.items(), key=lambda x: x[1], reverse=True)
            names = [x[0] for x in sorted_cv]
            values = [x[1] for x in sorted_cv]

            bars = ax.barh(range(len(names)), values,
                          color=VARIANT_META[v_key]['color'], alpha=0.7)
            ax.set_yticks(range(len(names)))
            ax.set_yticklabels(names, fontsize=8)
            ax.set_xlabel('Coefficient of Variation')
            ax.set_title(f'Parameter Identifiability — Variant {v_key}: '
                        f'{VARIANT_META[v_key]["label"]}')
            ax.axvline(x=0.5, color='red', linestyle='--', alpha=0.5,
                       label='High CV threshold')
            ax.legend()
            plt.tight_layout()
            plt.savefig(os.path.join(plots_dir, f"param_identifiability_{v_key}.png"), dpi=150)
            plt.close()
            print(f"  Plot: param_identifiability_{v_key}.png")


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Collect and compare calibration results across variants",
    )
    parser.add_argument('--prepared-dir', required=True,
                       help='Prepared basins directory')
    parser.add_argument('--output-dir', default=None,
                       help='Results output directory (default: prepared-dir/results)')
    parser.add_argument('--variants', default='A,B,C,D',
                       help='Comma-separated variants')
    parser.add_argument('--plot', action='store_true',
                       help='Generate comparison plots')

    args = parser.parse_args()

    variants = [v.strip().upper() for v in args.variants.split(',')]
    prepared_dir = args.prepared_dir
    output_dir = args.output_dir or os.path.join(prepared_dir, "results")
    os.makedirs(output_dir, exist_ok=True)

    print("\n" + "=" * 70)
    print("COLLECTING CALIBRATION RESULTS")
    print("=" * 70)

    # Collect results per variant
    all_results = {}
    for v_key in variants:
        print(f"\nCollecting Variant {v_key}...")
        df = collect_variant_results(prepared_dir, v_key)
        all_results[v_key] = df

        # Save per-variant CSV
        csv_path = os.path.join(output_dir, f"results_variant_{v_key}.csv")
        df.to_csv(csv_path, index=False)
        print(f"  Saved: {csv_path}")
        print(f"  Basins: {len(df)}")
        print(f"  Status: {df['status'].value_counts().to_dict()}")

    # Cross-variant comparison
    print("\nComputing cross-variant comparison...")
    cross_variant = compute_cross_variant_stats(all_results)

    cross_csv = os.path.join(output_dir, "cross_variant_comparison.csv")
    cross_variant.to_csv(cross_csv, index=False)
    print(f"  Saved: {cross_csv}")

    # Generate summary report
    print("\nGenerating summary report...")
    generate_summary_report(all_results, cross_variant, output_dir)

    # Generate plots
    if args.plot:
        print("\nGenerating plots...")
        generate_plots(all_results, cross_variant, output_dir)

    print(f"\nAll results saved to: {output_dir}")


if __name__ == "__main__":
    main()

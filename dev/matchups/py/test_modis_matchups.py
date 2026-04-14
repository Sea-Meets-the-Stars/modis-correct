"""
test_modis_matchups.py - Test the bow-tie match-up finder.

Usage:
    python test_modis_matchups.py --d_max 500
    python test_modis_matchups.py --d_max 1000 --n_bt 200

Outputs figures and results to dev/matchups/outputs/,
with filenames tagged by d_max.
"""

import argparse
import os
import time

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from modis_matchups import (
    load_latlon,
    reshape_to_bowties,
    find_matchups_bruteforce,
    build_candidate_table,
    find_matchups,
    validate,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Test MODIS bow-tie match-up finder.')
    parser.add_argument('--d_max', type=float, required=True,
                        help='Maximum separation distance in meters.')
    parser.add_argument('--n_bt', type=int, default=500,
                        help='Number of bow-ties to use (default: 500).')
    return parser.parse_args()


# ============================================================
# Plotting (writes to outdir with d_max in filename)
# ============================================================

def plot_separation_vs_column(matchups_df, d_max, outdir):
    """Mean separation vs column for each detector pair."""
    fig, ax = plt.subplots(figsize=(12, 6))

    for (di, dj), grp in matchups_df.groupby(['det_i', 'det_j']):
        mean_sep = grp.groupby('col_P')['sep_m'].mean()
        ax.plot(mean_sep.index, mean_sep.values, '.', ms=2,
                label=f'D{di}-D{dj}')

    ax.axvline(x=677, color='gray', ls='--', alpha=0.5, label='Nadir')
    ax.set_xlabel('Column P')
    ax.set_ylabel('Mean separation (m)')
    ax.set_title(f'Bow-tie match-up separation by column (d_max={d_max:.0f} m)')
    ax.legend(fontsize=8, ncol=3, loc='upper center')
    ax.set_ylim(0, d_max * 1.1)
    plt.tight_layout()

    path = os.path.join(outdir, f'separation_vs_column_dmax{d_max:.0f}.png')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f'  Saved {path}')


def plot_matchups_per_column(matchups_df, d_max, outdir):
    """Match-up count per column across the swath."""
    fig, ax = plt.subplots(figsize=(12, 6))

    counts = matchups_df.groupby('col_P').size()
    ax.bar(counts.index, counts.values, width=1, alpha=0.7, color='steelblue')
    ax.axvline(x=677, color='red', ls='--', alpha=0.5, label='Nadir')
    ax.set_xlabel('Column P')
    ax.set_ylabel('Number of match-ups')
    ax.set_title(f'Match-up count by column (d_max={d_max:.0f} m)')
    ax.legend()
    plt.tight_layout()

    path = os.path.join(outdir, f'matchups_per_column_dmax{d_max:.0f}.png')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f'  Saved {path}')


def plot_matchup_map(lat_bt, lon_bt, matchups_df, d_max, bt_idx, outdir):
    """Map of match-up pairs at one bow-tie boundary."""
    mu = matchups_df[matchups_df['bowtie_k'] == bt_idx]
    if len(mu) == 0:
        print(f'  No match-ups at bow-tie {bt_idx}, skipping map.')
        return

    fig, ax = plt.subplots(figsize=(14, 6))

    # Background: trailing detectors of B_k, leading of B_{k+1}
    for d in range(6, 10):
        ax.plot(lon_bt[bt_idx, d, :], lat_bt[bt_idx, d, :],
                '.', ms=0.3, color='blue', alpha=0.3)
    for d in range(0, 3):
        ax.plot(lon_bt[bt_idx + 1, d, :], lat_bt[bt_idx + 1, d, :],
                '.', ms=0.3, color='red', alpha=0.3)

    # Draw match-up connections
    for _, row in mu.iterrows():
        di0 = int(row.det_i) - 1
        dj0 = int(row.det_j) - 1
        c0 = int(row.col_P) - 1
        ax.plot(
            [lon_bt[bt_idx, di0, c0], lon_bt[bt_idx + 1, dj0, c0]],
            [lat_bt[bt_idx, di0, c0], lat_bt[bt_idx + 1, dj0, c0]],
            'g-', lw=0.5, alpha=0.5)

    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')
    ax.set_title(f'Match-ups at bow-tie boundary k={bt_idx} '
                 f'({len(mu)} pairs, d_max={d_max:.0f} m)')
    legend_el = [
        Line2D([0], [0], marker='.', color='blue', ls='', alpha=0.5,
               label='B_k  D7-10'),
        Line2D([0], [0], marker='.', color='red', ls='', alpha=0.5,
               label='B_{k+1} D1-3'),
        Line2D([0], [0], color='green', lw=1, label='Match-up'),
    ]
    ax.legend(handles=legend_el, fontsize=8)
    plt.tight_layout()

    path = os.path.join(outdir, f'matchup_map_bt{bt_idx}_dmax{d_max:.0f}.png')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f'  Saved {path}')


# ============================================================
# Main
# ============================================================

def main():
    args = parse_args()
    d_max = args.d_max
    n_bt_test = args.n_bt

    # Output directory
    outdir = os.path.join(os.path.dirname(__file__), '..', 'outputs')
    os.makedirs(outdir, exist_ok=True)

    # Data file
    data_file = os.path.join(os.path.dirname(__file__), '..', 'data', 'latlon.mat')

    # --- Load and reshape ---
    print(f'Loading lat/lon data from {data_file}...')
    lat, lon = load_latlon(data_file)
    print(f'  {lat.shape[0]} scan lines x {lat.shape[1]} pixels')

    lat_bt, lon_bt, n_skip = reshape_to_bowties(lat, lon)
    print(f'  Skipped {n_skip} lines (partial first bow-tie)')
    print(f'  {lat_bt.shape[0]} complete bow-ties')

    # Subset
    n_bt_use = min(n_bt_test, lat_bt.shape[0])
    lat_bt_sub = lat_bt[:n_bt_use]
    lon_bt_sub = lon_bt[:n_bt_use]
    print(f'  Using {n_bt_use} bow-ties')

    # --- Brute-force ---
    print(f'\nBrute-force search (d_max={d_max:.0f} m)...')
    t0 = time.time()
    mu_bf = find_matchups_bruteforce(lat_bt_sub, lon_bt_sub, d_max)
    t_bf = time.time() - t0
    print(f'  {len(mu_bf)} match-ups in {t_bf:.2f} s')

    if len(mu_bf) == 0:
        print('  No match-ups found. Try a larger d_max.')
        return

    # Summary by detector pair
    print('\n  Detector pair summary:')
    for (di, dj), grp in mu_bf.groupby(['det_i', 'det_j']):
        print(f'    D{di}-D{dj}: {len(grp):>7d} match-ups  '
              f'cols [{grp.col_P.min():>4d}, {grp.col_P.max():>4d}]  '
              f'sep [{grp.sep_m.min():>6.0f}, {grp.sep_m.max():>6.0f}] m')

    # --- Build candidate table and run optimized ---
    print('\nBuilding candidate table (margin=1.5)...')
    cand = build_candidate_table(mu_bf, margin_factor=1.5)
    for (di, dj), cols in sorted(cand.items()):
        print(f'  D{di}-D{dj}: cols {cols[0]:>4d}-{cols[-1]:>4d}')

    print(f'\nOptimized search (d_max={d_max:.0f} m)...')
    t0 = time.time()
    mu_opt = find_matchups(lat_bt_sub, lon_bt_sub, d_max, cand)
    t_opt = time.time() - t0
    print(f'  {len(mu_opt)} match-ups in {t_opt:.2f} s')
    print(f'  Speedup: {t_bf / t_opt:.1f}x')

    # --- Validate ---
    print('\nValidation:')
    validate(mu_bf, mu_opt)

    # --- Save results ---
    results_path = os.path.join(outdir, f'matchups_dmax{d_max:.0f}.csv')
    mu_bf.to_csv(results_path, index=False)
    print(f'\nResults saved to {results_path}')

    # --- Figures ---
    print('\nGenerating figures...')
    plot_separation_vs_column(mu_bf, d_max, outdir)
    plot_matchups_per_column(mu_bf, d_max, outdir)
    plot_matchup_map(lat_bt_sub, lon_bt_sub, mu_bf, d_max,
                     bt_idx=50, outdir=outdir)

    print('\nDone.')


if __name__ == '__main__':
    main()

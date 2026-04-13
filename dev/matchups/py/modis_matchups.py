"""
modis_matchups.py - Find bow-tie match-ups in MODIS satellite data.

A match-up is a pair of pixels from adjacent bow-ties whose
center-to-center separation distance is less than d_max meters.
Uses the Haversine formula on pixel-center lat/lon.

Implements:
  - Approach 3 (brute-force): all 12 detector pairs at every column
  - Approach 2 (optimized): candidate table restricts detector pairs per column
  - Validation: confirms optimized results match brute-force
"""

import numpy as np
import pandas as pd
import scipy.io
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import time
import os

# ============================================================
# Constants
# ============================================================

N_DETECTORS = 10
N_PIXELS = 1354
FIRST_DETECTOR = 6   # first scan line of each orbit is D=6
R_EARTH = 6_371_000.0  # Earth radius in meters

# Brute-force detector pairs: trailing D_i (7-10) vs leading D_j (1-3)
# Stored as 0-indexed detector indices within a bow-tie
BRUTE_FORCE_PAIRS = [
    (d_i, d_j)
    for d_i in range(6, 10)   # D=7,8,9,10
    for d_j in range(0, 3)    # D=1,2,3
]

# Figure output directory
FIG_DIR = '/home/xavier/Projects/overleaf/modis-correct'


# ============================================================
# Core functions
# ============================================================

def haversine(lat1, lon1, lat2, lon2):
    """Vectorized Haversine great-circle distance in meters.

    All inputs in degrees; arrays must be broadcastable.
    """
    lat1_r = np.radians(lat1)
    lat2_r = np.radians(lat2)
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)
    a = (np.sin(dlat / 2.0) ** 2
         + np.cos(lat1_r) * np.cos(lat2_r) * np.sin(dlon / 2.0) ** 2)
    return 2.0 * R_EARTH * np.arcsin(np.sqrt(a))


def load_latlon(filepath):
    """Load latitude and longitude from a .mat file.

    Returns arrays of shape (N_lines, 1354) in float32.
    """
    data = scipy.io.loadmat(filepath)
    # .mat file stores shape (1354, N_lines); transpose to (N_lines, 1354)
    lat = data['lat'].T.astype(np.float32)
    lon = data['lon'].T.astype(np.float32)
    return lat, lon


def reshape_to_bowties(lat, lon, first_detector=FIRST_DETECTOR):
    """Group scan lines into complete 10-detector bow-ties.

    The first scan line has detector D=first_detector (1-indexed).
    The leading partial bow-tie and any trailing partial are discarded.

    Returns
    -------
    lat_bt, lon_bt : ndarray, shape (N_bt, 10, 1354)
    n_skip : int
        Number of scan lines discarded from the start.
    """
    n_lines = lat.shape[0]
    # Partial first bow-tie: detectors D=first_detector through D=10
    n_skip = N_DETECTORS - (first_detector - 1)  # = 11 - first_detector
    n_remaining = n_lines - n_skip
    n_bt = n_remaining // N_DETECTORS
    end = n_skip + n_bt * N_DETECTORS

    lat_bt = lat[n_skip:end].reshape(n_bt, N_DETECTORS, N_PIXELS)
    lon_bt = lon[n_skip:end].reshape(n_bt, N_DETECTORS, N_PIXELS)
    return lat_bt, lon_bt, n_skip


def find_matchups_bruteforce(lat_bt, lon_bt, d_max):
    """Approach 3: brute-force search over all 12 detector pairs at every column.

    For each pair (D_i, D_j), computes Haversine across all bow-tie
    boundaries and all 1354 columns in one vectorized call.

    Returns a DataFrame of match-ups.
    """
    n_bt = lat_bt.shape[0]
    # Mirror sides alternate starting at M=1
    mirror = np.empty(n_bt, dtype=np.int8)
    mirror[0::2] = 1
    mirror[1::2] = 2

    results = []
    for d_i, d_j in BRUTE_FORCE_PAIRS:
        # shape (N_bt - 1, 1354)
        lat_i = lat_bt[:-1, d_i, :]
        lon_i = lon_bt[:-1, d_i, :]
        lat_j = lat_bt[1:, d_j, :]
        lon_j = lon_bt[1:, d_j, :]

        dist = haversine(lat_i, lon_i, lat_j, lon_j)
        bt_idx, col_idx = np.where(dist < d_max)

        if len(bt_idx) > 0:
            results.append(pd.DataFrame({
                'bowtie_k': bt_idx,
                'col_P': col_idx + 1,           # 1-indexed
                'det_i': d_i + 1,                # 1-indexed detector
                'det_j': d_j + 1,
                'sep_m': dist[bt_idx, col_idx],
                'mirror_k': mirror[bt_idx],
                'mirror_k1': mirror[bt_idx + 1],
            }))

    if not results:
        return pd.DataFrame(columns=[
            'bowtie_k', 'col_P', 'det_i', 'det_j',
            'sep_m', 'mirror_k', 'mirror_k1'])

    df = pd.concat(results, ignore_index=True)
    df.sort_values(['bowtie_k', 'col_P'], inplace=True)
    return df.reset_index(drop=True)


def build_candidate_table(matchups_df, margin_factor=1.5):
    """Build a candidate table from brute-force results.

    Returns a dict mapping (det_i, det_j) -> array of 1-indexed column
    indices where that pair should be checked.
    """
    table = {}
    for (di, dj), grp in matchups_df.groupby(['det_i', 'det_j']):
        cols = grp['col_P'].values
        cmin, cmax = cols.min(), cols.max()
        # Expand range by margin
        mid = (cmin + cmax) / 2.0
        half = (cmax - cmin) / 2.0
        exp_half = half * margin_factor
        new_min = max(1, int(mid - exp_half))
        new_max = min(N_PIXELS, int(mid + exp_half))
        table[(di, dj)] = np.arange(new_min, new_max + 1)
    return table


def find_matchups(lat_bt, lon_bt, d_max, candidate_table):
    """Approach 2: optimized search using candidate table.

    Only checks detector pairs at columns where the table says
    match-ups are possible.  Haversine is vectorized across all
    bow-tie boundaries and candidate columns in one call per pair.

    Returns a DataFrame of match-ups.
    """
    n_bt = lat_bt.shape[0]
    mirror = np.empty(n_bt, dtype=np.int8)
    mirror[0::2] = 1
    mirror[1::2] = 2

    results = []
    for (det_i, det_j), cols in candidate_table.items():
        d_i = det_i - 1  # 0-indexed
        d_j = det_j - 1
        ci = cols - 1     # 0-indexed column indices

        # shape (N_bt - 1, len(cols))
        lat_i = lat_bt[:-1, d_i, :][:, ci]
        lon_i = lon_bt[:-1, d_i, :][:, ci]
        lat_j = lat_bt[1:, d_j, :][:, ci]
        lon_j = lon_bt[1:, d_j, :][:, ci]

        dist = haversine(lat_i, lon_i, lat_j, lon_j)
        bt_idx, c_idx = np.where(dist < d_max)

        if len(bt_idx) > 0:
            results.append(pd.DataFrame({
                'bowtie_k': bt_idx,
                'col_P': cols[c_idx],       # already 1-indexed
                'det_i': det_i,
                'det_j': det_j,
                'sep_m': dist[bt_idx, c_idx],
                'mirror_k': mirror[bt_idx],
                'mirror_k1': mirror[bt_idx + 1],
            }))

    if not results:
        return pd.DataFrame(columns=[
            'bowtie_k', 'col_P', 'det_i', 'det_j',
            'sep_m', 'mirror_k', 'mirror_k1'])

    df = pd.concat(results, ignore_index=True)
    df.sort_values(['bowtie_k', 'col_P'], inplace=True)
    return df.reset_index(drop=True)


def validate(matchups_bf, matchups_opt):
    """Confirm optimized results contain every brute-force match-up."""
    keys = ['bowtie_k', 'col_P', 'det_i', 'det_j']
    merged = matchups_bf[keys].merge(
        matchups_opt[keys], on=keys, how='left', indicator=True)
    n_missing = (merged['_merge'] == 'left_only').sum()
    if n_missing == 0:
        print(f'  PASSED: all {len(matchups_bf)} brute-force '
              f'match-ups found in optimized results.')
        return True
    print(f'  FAILED: {n_missing}/{len(matchups_bf)} match-ups missing.')
    return False


# ============================================================
# Plotting
# ============================================================

def plot_separation_vs_column(matchups_df, d_max, outdir=FIG_DIR):
    """Scatter of mean separation vs column for each detector pair."""
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

    path = os.path.join(outdir, 'matchup_separation_vs_column.png')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f'  Saved {path}')


def plot_matchups_per_column(matchups_df, outdir=FIG_DIR):
    """Bar chart of total match-up count per column."""
    fig, ax = plt.subplots(figsize=(12, 6))

    counts = matchups_df.groupby('col_P').size()
    ax.bar(counts.index, counts.values, width=1, alpha=0.7, color='steelblue')
    ax.axvline(x=677, color='red', ls='--', alpha=0.5, label='Nadir')
    ax.set_xlabel('Column P')
    ax.set_ylabel('Number of match-ups')
    ax.set_title('Match-up count by column (all bow-tie boundaries)')
    ax.legend()
    plt.tight_layout()

    path = os.path.join(outdir, 'matchups_per_column.png')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f'  Saved {path}')


def plot_matchup_map(lat_bt, lon_bt, matchups_df, bt_idx=50,
                     outdir=FIG_DIR):
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
    ax.set_title(f'Match-ups at bow-tie boundary k={bt_idx} ({len(mu)} pairs)')
    legend_el = [
        Line2D([0], [0], marker='.', color='blue', ls='', alpha=0.5,
               label='B_k  D7-10'),
        Line2D([0], [0], marker='.', color='red', ls='', alpha=0.5,
               label='B_{k+1} D1-3'),
        Line2D([0], [0], color='green', lw=1, label='Match-up'),
    ]
    ax.legend(handles=legend_el, fontsize=8)
    plt.tight_layout()

    path = os.path.join(outdir, f'matchup_map_bt{bt_idx}.png')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f'  Saved {path}')


# ============================================================
# Main
# ============================================================

def main():
    data_file = os.path.join(
        os.path.dirname(__file__), '..', 'data', 'latlon.mat')
    d_max = 1000.0    # meters
    n_bt_test = 500   # bow-ties for development testing

    # --- Load ---
    print('Loading lat/lon data...')
    lat, lon = load_latlon(data_file)
    print(f'  Raw shape: {lat.shape}  ({lat.shape[0]} scan lines)')

    # --- Reshape ---
    lat_bt, lon_bt, n_skip = reshape_to_bowties(lat, lon)
    print(f'  Skipped {n_skip} lines (partial first bow-tie, D=6-10)')
    print(f'  Complete bow-ties: {lat_bt.shape[0]}')

    # --- Subset for testing ---
    lat_bt_sub = lat_bt[:n_bt_test]
    lon_bt_sub = lon_bt[:n_bt_test]
    print(f'  Using first {n_bt_test} bow-ties for testing')

    # --- Brute-force (Approach 3) ---
    print(f'\nBrute-force search (d_max={d_max:.0f} m)...')
    t0 = time.time()
    mu_bf = find_matchups_bruteforce(lat_bt_sub, lon_bt_sub, d_max)
    t_bf = time.time() - t0
    print(f'  Found {len(mu_bf)} match-ups in {t_bf:.2f} s')

    # Summary by detector pair
    print('\n  Detector pair summary:')
    for (di, dj), grp in mu_bf.groupby(['det_i', 'det_j']):
        print(f'    D{di}-D{dj}: {len(grp):>7d} match-ups  '
              f'cols [{grp.col_P.min():>4d}, {grp.col_P.max():>4d}]  '
              f'sep [{grp.sep_m.min():>6.0f}, {grp.sep_m.max():>6.0f}] m')

    # --- Build candidate table ---
    print('\nBuilding candidate table (margin=1.5)...')
    cand = build_candidate_table(mu_bf, margin_factor=1.5)
    for (di, dj), cols in sorted(cand.items()):
        print(f'  D{di}-D{dj}: cols {cols[0]:>4d}-{cols[-1]:>4d} '
              f'({len(cols)} cols)')

    # --- Optimized (Approach 2) ---
    print(f'\nOptimized search (d_max={d_max:.0f} m)...')
    t0 = time.time()
    mu_opt = find_matchups(lat_bt_sub, lon_bt_sub, d_max, cand)
    t_opt = time.time() - t0
    print(f'  Found {len(mu_opt)} match-ups in {t_opt:.2f} s')
    print(f'  Speedup: {t_bf / t_opt:.1f}x')

    # --- Validate ---
    print('\nValidation:')
    validate(mu_bf, mu_opt)

    # --- Figures ---
    print('\nGenerating figures...')
    os.makedirs(FIG_DIR, exist_ok=True)
    plot_separation_vs_column(mu_bf, d_max)
    plot_matchups_per_column(mu_bf)
    plot_matchup_map(lat_bt_sub, lon_bt_sub, mu_bf, bt_idx=50)

    print('\nDone.')


if __name__ == '__main__':
    main()

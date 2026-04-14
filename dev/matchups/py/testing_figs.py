"""
testing_figs.py - Generate diagnostic figures from match-up results.

Reads match-up CSV from dev/matchups/outputs/ and produces figures
for visual inspection.
"""

import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Directories
SCRIPT_DIR = os.path.dirname(__file__)
OUTPUT_DIR = os.path.join(SCRIPT_DIR, '..', 'outputs')
FIG_DIR = '/home/xavier/Projects/overleaf/modis-correct'


def plot_d10_d1_matchups(csv_path, outdir=FIG_DIR):
    """Scatter plot of D10-D1 match-ups: Column vs Row, colored by separation.

    Row is defined as (bowtie_k - 1) * 10 + detector number.
    Detector 10 is in B_k, detector 1 is in B_{k+1}.
    We plot the row of the D10 pixel.
    """
    df = pd.read_csv(csv_path)

    # Filter to D10-D1 pairs
    mask = (df['det_i'] == 10) & (df['det_j'] == 1)
    mu = df[mask].copy()
    print(f'D10-D1 match-ups: {len(mu)}')

    # Compute row: (bowtie_k) * 10 + detector number
    # D10 is detector 10 in bow-tie B_k (bowtie_k is 0-indexed)
    mu['row'] = mu['bowtie_k'] * 10 + 10

    # Extract d_max from filename for the title
    basename = os.path.basename(csv_path)
    d_max_str = basename.replace('matchups_dmax', '').replace('.csv', '')

    fig, ax = plt.subplots(figsize=(12, 8))
    sc = ax.scatter(mu['col_P'], mu['row'], c=mu['sep_m'],
                    cmap='viridis_r', s=1, alpha=0.7)
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label('Separation (m)')
    ax.set_xlabel('Column P')
    ax.set_ylabel('Row  [(bowtie_k) * 10 + detector]')
    ax.set_title(f'D10-D1 match-ups (d_max={d_max_str} m)')
    plt.tight_layout()

    path = os.path.join(outdir, f'd10_d1_matchups_dmax{d_max_str}.png')
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f'Saved {path}')


if __name__ == '__main__':
    csv_path = os.path.join(OUTPUT_DIR, 'matchups_dmax40.csv')
    os.makedirs(FIG_DIR, exist_ok=True)
    plot_d10_d1_matchups(csv_path)

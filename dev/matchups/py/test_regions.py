"""
test_regions.py - Diagnostic figures for the overlap-region algorithm.

Provides methods to visualize and validate the intersection-finding
step of the Approach B algorithm (polynomial curve-fit version).
"""

import numpy as np
import matplotlib.pyplot as plt
import os

from modis_matchups import (
    haversine, load_latlon, reshape_to_bowties, N_PIXELS
)
from find_regions import (
    find_intersection, find_regions_one_boundary,
    DETECTOR_PAIRS, NADIR_COL
)

FIG_DIR = '/home/xavier/Projects/overleaf/modis-correct'


def plot_intersection_diagnostic(lat_bt, lon_bt, k, det_i=10, det_j=1,
                                 side='left', poly_deg=8, outdir=FIG_DIR):
    """Generate a diagnostic figure for the polynomial curve-fit intersection.

    Shows three panels:
      1. Smooth polynomial fits to both scan lines' latitude, zoomed to the
         crossing region, with raw data points underneath.
      2. Smooth delta_lat (fitted meas - fitted ref) vs column, showing the
         single clean zero-crossing from the polynomial roots.
      3. Raw delta_lat for comparison, showing the many spurious sign changes
         from geolocation quantization that motivated the new approach.

    Parameters
    ----------
    lat_bt, lon_bt : (N_bt, 10, 1354) bow-tie arrays
    k : int, bow-tie boundary index (0-indexed)
    det_i, det_j : 1-indexed detector numbers (trailing, leading)
    side : 'left' or 'right'
    poly_deg : polynomial degree for curve fits
    outdir : directory for saving the figure
    """
    # Extract scan lines
    lat_meas = lat_bt[k, det_i - 1, :]
    lon_meas = lon_bt[k, det_i - 1, :]
    lat_ref = lat_bt[k + 1, det_j - 1, :]
    lon_ref = lon_bt[k + 1, det_j - 1, :]

    # Column range
    if side == 'left':
        cols = np.arange(0, NADIR_COL)
    else:
        cols = np.arange(NADIR_COL, N_PIXELS)

    x = cols.astype(np.float64)

    # --- Polynomial fits ---
    p_ref = np.polyfit(x, lat_ref[cols].astype(np.float64), poly_deg)
    p_meas = np.polyfit(x, lat_meas[cols].astype(np.float64), poly_deg)
    lat_ref_fit = np.polyval(p_ref, x)
    lat_meas_fit = np.polyval(p_meas, x)
    delta_fit = lat_meas_fit - lat_ref_fit

    # Run the actual algorithm
    cross_col = find_intersection(lat_ref, lon_ref, lat_meas, lon_meas,
                                  cols, poly_deg=poly_deg)

    # Raw delta_lat for comparison
    delta_raw = (lat_meas[cols] - lat_ref[cols]).astype(np.float64)
    raw_signs = np.sign(delta_raw)
    raw_sign_changes = np.where(np.diff(raw_signs) != 0)[0]

    # --- Figure: 3 panels ---
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Panel 1: Polynomial fits to latitude curves, zoomed to crossing
    ax = axes[0]
    # Determine zoom range around the crossing
    if cross_col is not None:
        zoom_center = int(round(cross_col))
    else:
        zoom_center = len(cols) // 2
    zoom_hw = 60
    z0 = max(0, zoom_center - cols[0] - zoom_hw)
    z1 = min(len(cols), zoom_center - cols[0] + zoom_hw + 1)
    zslice = slice(z0, z1)

    # Raw points
    ax.plot(cols[zslice] + 1, lat_ref[cols[zslice]], '.', ms=3,
            color='lightsalmon', alpha=0.5, label='D$_j$ raw')
    ax.plot(cols[zslice] + 1, lat_meas[cols[zslice]], '.', ms=3,
            color='lightblue', alpha=0.5, label='D$_i$ raw')
    # Fitted curves (thicker)
    ax.plot(cols[zslice] + 1, lat_ref_fit[zslice], '-', lw=2, color='red',
            label=f'D{det_j} B$_{{k+1}}$ fit (deg {poly_deg})')
    ax.plot(cols[zslice] + 1, lat_meas_fit[zslice], '-', lw=2, color='blue',
            label=f'D{det_i} B$_k$ fit (deg {poly_deg})')
    if cross_col is not None:
        # Mark crossing
        lat_at_cross = np.polyval(p_ref, cross_col)
        ax.plot(cross_col + 1, lat_at_cross, '*', ms=15, color='green',
                markeredgecolor='black', markeredgewidth=0.5, zorder=5,
                label=f'Crossing P={cross_col+1:.1f}')
    ax.set_xlabel('Column P')
    ax.set_ylabel('Latitude (deg)')
    ax.set_title('Smooth polynomial fits (zoomed)')
    ax.legend(fontsize=7, loc='best')

    # Panel 2: Smooth delta_lat from polynomial fits
    ax = axes[1]
    ax.plot(cols + 1, delta_fit * 111_000, '-', lw=1.2, color='steelblue',
            label='$\\Delta$lat fitted')
    ax.axhline(0, color='red', ls='--', lw=0.8)
    if cross_col is not None:
        ax.axvline(cross_col + 1, color='green', lw=1.5,
                   label=f'Crossing P={cross_col+1:.1f}')
        ax.plot(cross_col + 1, 0, '*', ms=15, color='green',
                markeredgecolor='black', markeredgewidth=0.5, zorder=5)
    ax.set_xlabel('Column P')
    ax.set_ylabel('$\\Delta$lat fitted (approx m)')
    ax.set_title('Smooth $\\Delta$lat: single clean crossing')
    ax.legend(fontsize=8)

    # Panel 3: Raw delta_lat showing the noise problem
    ax = axes[2]
    ax.plot(cols + 1, delta_raw * 111_000, '.', ms=1, color='steelblue',
            alpha=0.6, label='Raw $\\Delta$lat')
    ax.plot(cols + 1, delta_fit * 111_000, '-', lw=1, color='darkorange',
            alpha=0.8, label='Fitted')
    ax.axhline(0, color='red', ls='--', lw=0.8)
    # Mark all raw sign changes
    for sc in raw_sign_changes:
        ax.axvline(cols[sc] + 1, color='orange', alpha=0.15, lw=0.8)
    if len(raw_sign_changes) > 0:
        ax.axvline(cols[raw_sign_changes[0]] + 1, color='orange', alpha=0.15,
                   lw=0.8, label=f'{len(raw_sign_changes)} raw sign changes')
    if cross_col is not None:
        ax.axvline(cross_col + 1, color='green', lw=1.5,
                   label=f'Poly-fit crossing P={cross_col+1:.1f}')
    ax.set_xlabel('Column P')
    ax.set_ylabel('$\\Delta$lat (approx m)')
    ax.set_title('Raw vs fitted: noise eliminated')
    ax.legend(fontsize=7)

    plt.suptitle(f'Intersection diagnostic: D{det_i}-D{det_j} at boundary k={k} ({side})',
                 fontsize=13)
    plt.tight_layout()
    if outdir:
        path = os.path.join(outdir, 'region_intersection_diagnostic.png')
        fig.savefig(path, dpi=150)
        print(f'Saved {path}')
    plt.show()
    return fig


def plot_intersection_multi_boundary(lat_bt, lon_bt, boundaries,
                                     det_i=10, det_j=1, side='left',
                                     poly_deg=8, outdir=FIG_DIR):
    """Show the polynomial-fit intersection across multiple bow-tie boundaries.

    For each boundary, plots the smooth delta_lat curve with the crossing
    marked, confirming stability across different latitudes.
    """
    n = len(boundaries)
    fig, axes = plt.subplots(1, n, figsize=(4.5 * n, 4), squeeze=False)
    axes = axes[0]

    if side == 'left':
        cols = np.arange(0, NADIR_COL)
    else:
        cols = np.arange(NADIR_COL, N_PIXELS)

    x = cols.astype(np.float64)

    for i, k in enumerate(boundaries):
        lat_meas = lat_bt[k, det_i - 1, :]
        lon_meas = lon_bt[k, det_i - 1, :]
        lat_ref = lat_bt[k + 1, det_j - 1, :]
        lon_ref = lon_bt[k + 1, det_j - 1, :]

        # Polynomial fits
        p_ref = np.polyfit(x, lat_ref[cols].astype(np.float64), poly_deg)
        p_meas = np.polyfit(x, lat_meas[cols].astype(np.float64), poly_deg)
        delta_fit = np.polyval(p_meas, x) - np.polyval(p_ref, x)

        cross = find_intersection(lat_ref, lon_ref, lat_meas, lon_meas,
                                  cols, poly_deg=poly_deg)
        lat_mid = float(lat_bt[k, 5, 677])

        ax = axes[i]
        ax.plot(cols + 1, delta_fit * 111_000, '-', lw=1, color='steelblue')
        ax.axhline(0, color='red', ls='--', lw=0.8)
        if cross is not None:
            ax.axvline(cross + 1, color='green', lw=1.5,
                       label=f'P={cross+1:.0f}')
            ax.plot(cross + 1, 0, '*', ms=12, color='green',
                    markeredgecolor='black', markeredgewidth=0.5)
        ax.set_xlabel('Column P')
        if i == 0:
            ax.set_ylabel('$\\Delta$lat fitted (m)')
        ax.set_title(f'k={k}, lat={lat_mid:+.0f}$^\\circ$')
        ax.legend(fontsize=8)

    plt.suptitle(f'D{det_i}-D{det_j} {side}: polynomial-fit crossing across latitudes '
                 f'(deg {poly_deg})', fontsize=12)
    plt.tight_layout()
    if outdir:
        path = os.path.join(outdir, 'region_intersection_multi.png')
        fig.savefig(path, dpi=150)
        print(f'Saved {path}')
    plt.show()
    return fig


# ============================================================
# Main: run diagnostics when executed as a script
# ============================================================

if __name__ == '__main__':
    data_file = os.path.join(os.path.dirname(__file__), '..', 'data', 'latlon.mat')
    print('Loading data...')
    lat, lon = load_latlon(data_file)
    lat_bt, lon_bt, n_skip = reshape_to_bowties(lat, lon)
    print(f'{lat_bt.shape[0]} bow-ties')

    os.makedirs(FIG_DIR, exist_ok=True)

    # Detailed diagnostic for one boundary
    print('\n--- Single-boundary diagnostic ---')
    plot_intersection_diagnostic(lat_bt, lon_bt, k=900,
                                 det_i=10, det_j=1, side='left')

    # Multi-boundary stability check
    print('\n--- Multi-boundary stability ---')
    plot_intersection_multi_boundary(
        lat_bt, lon_bt,
        boundaries=[50, 500, 900, 2000, 3500],
        det_i=10, det_j=1, side='left')

    # Region counts with polynomial-fit intersection
    print('\n--- Region counts with polynomial-fit intersection ---')
    print(f'{"k":>5s} {"lat":>7s} {"regions":>8s}')
    for k in [1, 50, 100, 500, 900, 2000, 3000, 4000]:
        if k + 1 >= lat_bt.shape[0]:
            continue
        regions = find_regions_one_boundary(
            lat_bt[k], lon_bt[k], lat_bt[k + 1], lon_bt[k + 1],
            w=30, q_max=0.05)
        lat_mid = float(lat_bt[k, 5, 677])
        print(f'{k:5d} {lat_mid:+7.1f} {len(regions):8d}')

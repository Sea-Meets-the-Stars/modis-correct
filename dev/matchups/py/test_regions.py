"""
test_regions.py - Diagnostic figures for the overlap-region algorithm.

Provides methods to visualize and validate the intersection-finding
step of the Approach B algorithm (lat vs lon polynomial curve-fit version).
"""

import numpy as np
import matplotlib.pyplot as plt
import os

from modis_matchups import (
    haversine, load_latlon, reshape_to_bowties, N_PIXELS
)
from find_regions import (
    find_intersection, find_regions_one_boundary, _unwrap_lon,
    DETECTOR_PAIRS, NADIR_COL
)

FIG_DIR = '/home/xavier/Projects/overleaf/modis-correct'


def plot_intersection_diagnostic(lat_bt, lon_bt, k, det_i=10, det_j=1,
                                 side='left', poly_deg=8, outdir=FIG_DIR):
    """Generate a diagnostic figure for the lat-vs-lon curve-fit intersection.

    Shows three panels:
      1. Both scan lines in lat/lon space with polynomial fits, zoomed to
         the crossing region.  The fitted curves clearly cross each other.
      2. Smooth delta_lat as a function of longitude, showing the single
         clean zero-crossing from the polynomial roots.
      3. Raw lat/lon scatter overlaid with fitted curves, wider view to
         show the overall scan-line geometry.

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
    lat_meas_raw = lat_bt[k, det_i - 1, :]
    lon_meas_raw = lon_bt[k, det_i - 1, :]
    lat_ref_raw = lat_bt[k + 1, det_j - 1, :]
    lon_ref_raw = lon_bt[k + 1, det_j - 1, :]

    # Column range
    if side == 'left':
        cols = np.arange(0, NADIR_COL)
    else:
        cols = np.arange(NADIR_COL, N_PIXELS)

    lat_r = lat_ref_raw[cols].astype(np.float64)
    lon_r = _unwrap_lon(lon_ref_raw[cols].astype(np.float64))
    lat_m = lat_meas_raw[cols].astype(np.float64)
    lon_m = _unwrap_lon(lon_meas_raw[cols].astype(np.float64))

    # Center/scale longitude for conditioning (same as find_regions.py)
    lon_all = np.concatenate([lon_r, lon_m])
    lon_mu = np.mean(lon_all)
    lon_sd = np.std(lon_all)
    lon_r_n = (lon_r - lon_mu) / lon_sd
    lon_m_n = (lon_m - lon_mu) / lon_sd

    # Fit lat = f(lon_normalized) for each scan line
    p_ref = np.polyfit(lon_r_n, lat_r, poly_deg)
    p_meas = np.polyfit(lon_m_n, lat_m, poly_deg)

    # Evaluate fitted curves on a common dense longitude grid
    lon_lo = max(lon_r.min(), lon_m.min())
    lon_hi = min(lon_r.max(), lon_m.max())
    lon_dense = np.linspace(lon_lo, lon_hi, 2000)
    lon_dense_n = (lon_dense - lon_mu) / lon_sd
    lat_ref_fit = np.polyval(p_ref, lon_dense_n)
    lat_meas_fit = np.polyval(p_meas, lon_dense_n)
    delta_fit = lat_meas_fit - lat_ref_fit

    # Run the actual algorithm to get crossing column
    cross_col = find_intersection(lat_ref_raw, lon_ref_raw,
                                  lat_meas_raw, lon_meas_raw,
                                  cols, poly_deg=poly_deg)

    # Crossing longitude (for plotting)
    if cross_col is not None:
        cross_lon = lon_r[int(np.argmin(np.abs(cols - cross_col)))]
        cross_lat = np.polyval(p_ref, (cross_lon - lon_mu) / lon_sd)
    else:
        cross_lon = cross_lat = None

    # --- Figure: 3 panels ---
    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))

    # Panel 1: Zoomed lat/lon view of the two scan lines crossing
    ax = axes[0]
    if cross_col is not None:
        ci = int(np.argmin(np.abs(cols - cross_col)))
        zoom_hw = 50
        z0 = max(0, ci - zoom_hw)
        z1 = min(len(cols), ci + zoom_hw + 1)
    else:
        z0, z1 = 0, len(cols)
    zslice = slice(z0, z1)

    # Raw pixels
    ax.plot(lon_r[zslice], lat_r[zslice], '.', ms=3, color='lightsalmon',
            alpha=0.5)
    ax.plot(lon_m[zslice], lat_m[zslice], '.', ms=3, color='lightblue',
            alpha=0.5)
    # Fitted curves
    lon_z = lon_dense[(lon_dense >= lon_r[zslice].min()) &
                      (lon_dense <= lon_r[zslice].max())]
    lon_z_n = (lon_z - lon_mu) / lon_sd
    lat_ref_z = np.polyval(p_ref, lon_z_n)
    lat_meas_z = np.polyval(p_meas, lon_z_n)
    ax.plot(lon_z, lat_ref_z, '-', lw=2, color='red',
            label=f'D{det_j} B$_{{k+1}}$ fit')
    ax.plot(lon_z, lat_meas_z, '-', lw=2, color='blue',
            label=f'D{det_i} B$_k$ fit')
    if cross_lon is not None:
        ax.plot(cross_lon, cross_lat, '*', ms=15, color='green',
                markeredgecolor='black', markeredgewidth=0.5, zorder=5,
                label=f'Crossing P={cross_col + 1}')
    ax.set_xlabel('Longitude (deg)')
    ax.set_ylabel('Latitude (deg)')
    ax.set_title('Scan-line crossing in lat/lon (zoomed)')
    ax.legend(fontsize=7, loc='best')

    # Panel 2: Smooth delta_lat as function of longitude
    ax = axes[1]
    ax.plot(lon_dense, delta_fit * 111_000, '-', lw=1.2, color='steelblue',
            label='$\\Delta$lat(lon) fitted')
    ax.axhline(0, color='red', ls='--', lw=0.8)
    if cross_lon is not None:
        ax.axvline(cross_lon, color='green', lw=1.5,
                   label=f'Crossing lon={cross_lon:.2f}$^\\circ$')
        ax.plot(cross_lon, 0, '*', ms=15, color='green',
                markeredgecolor='black', markeredgewidth=0.5, zorder=5)
    ax.set_xlabel('Longitude (deg)')
    ax.set_ylabel('$\\Delta$lat fitted (approx m)')
    ax.set_title('Smooth $\\Delta$lat(lon): single clean crossing')
    ax.legend(fontsize=8)

    # Panel 3: Wider view — fitted curves over full side
    ax = axes[2]
    ax.plot(lon_r, lat_r, '.', ms=0.5, color='lightsalmon', alpha=0.3)
    ax.plot(lon_m, lat_m, '.', ms=0.5, color='lightblue', alpha=0.3)
    ax.plot(lon_dense, np.polyval(p_ref, lon_dense_n), '-', lw=1.2, color='red',
            label=f'D{det_j} fit')
    ax.plot(lon_dense, np.polyval(p_meas, lon_dense_n), '-', lw=1.2, color='blue',
            label=f'D{det_i} fit')
    if cross_lon is not None:
        ax.plot(cross_lon, cross_lat, '*', ms=12, color='green',
                markeredgecolor='black', markeredgewidth=0.5, zorder=5,
                label=f'Crossing P={cross_col + 1}')
    ax.set_xlabel('Longitude (deg)')
    ax.set_ylabel('Latitude (deg)')
    ax.set_title('Full-side view with fits')
    ax.legend(fontsize=7, loc='best')

    lat_mid = float(lat_bt[k, 5, 677])
    plt.suptitle(f'Intersection diagnostic: D{det_i}-D{det_j} at k={k} '
                 f'({side}, lat={lat_mid:+.0f}$^\\circ$, deg {poly_deg})',
                 fontsize=12)
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
    """Show the lat/lon curve-fit intersection across multiple boundaries.

    For each boundary, plots both scan lines in lat/lon space with fitted
    curves, confirming the crossing is found stably across latitudes.
    """
    n = len(boundaries)
    fig, axes = plt.subplots(1, n, figsize=(4.5 * n, 4), squeeze=False)
    axes = axes[0]

    for i, k in enumerate(boundaries):
        lat_meas_raw = lat_bt[k, det_i - 1, :]
        lon_meas_raw = lon_bt[k, det_i - 1, :]
        lat_ref_raw = lat_bt[k + 1, det_j - 1, :]
        lon_ref_raw = lon_bt[k + 1, det_j - 1, :]

        if side == 'left':
            cols = np.arange(0, NADIR_COL)
        else:
            cols = np.arange(NADIR_COL, N_PIXELS)

        lat_r = lat_ref_raw[cols].astype(np.float64)
        lon_r = _unwrap_lon(lon_ref_raw[cols].astype(np.float64))
        lat_m = lat_meas_raw[cols].astype(np.float64)
        lon_m = _unwrap_lon(lon_meas_raw[cols].astype(np.float64))

        lon_all = np.concatenate([lon_r, lon_m])
        lon_mu = np.mean(lon_all)
        lon_sd = np.std(lon_all)
        p_ref = np.polyfit((lon_r - lon_mu) / lon_sd, lat_r, poly_deg)
        p_meas = np.polyfit((lon_m - lon_mu) / lon_sd, lat_m, poly_deg)

        lon_lo = max(lon_r.min(), lon_m.min())
        lon_hi = min(lon_r.max(), lon_m.max())
        lon_d = np.linspace(lon_lo, lon_hi, 1000)
        lon_d_n = (lon_d - lon_mu) / lon_sd

        cross_col = find_intersection(lat_ref_raw, lon_ref_raw,
                                      lat_meas_raw, lon_meas_raw,
                                      cols, poly_deg=poly_deg)
        lat_mid = float(lat_bt[k, 5, 677])

        ax = axes[i]
        # Plot fitted curves (zoomed to crossing)
        delta = np.polyval(p_meas, lon_d_n) - np.polyval(p_ref, lon_d_n)
        ax.plot(lon_d, delta * 111_000, '-', lw=1, color='steelblue')
        ax.axhline(0, color='red', ls='--', lw=0.8)
        if cross_col is not None:
            cl = lon_r[int(np.argmin(np.abs(cols - cross_col)))]
            ax.axvline(cl, color='green', lw=1.5,
                       label=f'P={cross_col + 1}')
            ax.plot(cl, 0, '*', ms=12, color='green',
                    markeredgecolor='black', markeredgewidth=0.5)
        ax.set_xlabel('Longitude (deg)')
        if i == 0:
            ax.set_ylabel('$\\Delta$lat fitted (m)')
        ax.set_title(f'k={k}, lat={lat_mid:+.0f}$^\\circ$')
        ax.legend(fontsize=8)

    plt.suptitle(f'D{det_i}-D{det_j} {side}: lat/lon curve-fit crossing '
                 f'across latitudes (deg {poly_deg})', fontsize=12)
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

    # Region counts
    print('\n--- Region counts with lat/lon curve-fit intersection ---')
    print(f'{"k":>5s} {"lat":>7s} {"regions":>8s}')
    for k in [1, 50, 100, 500, 900, 2000, 3000, 4000]:
        if k + 1 >= lat_bt.shape[0]:
            continue
        regions = find_regions_one_boundary(
            lat_bt[k], lon_bt[k], lat_bt[k + 1], lon_bt[k + 1],
            w=30, q_max=0.05)
        lat_mid = float(lat_bt[k, 5, 677])
        print(f'{k:5d} {lat_mid:+7.1f} {len(regions):8d}')

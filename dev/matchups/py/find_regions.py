"""
find_regions.py - Find overlap regions between adjacent MODIS bow-ties.

Implements Approach B from find_regions_plan.tex:
  Intersection-centered line-fitting with perpendicular distance.

Algorithm for each detector pair on each side of nadir:
  1. Find the column where two scan lines cross (delta_lat sign change)
  2. Center a w=30 pixel window on the crossing point
  3. Fit a reference line through the leading scan line (OLS, equirectangular)
  4. Compute perpendicular distance from the trailing scan to the fitted line
  5. Normalize: q = |d_perp| / pixel_spacing
  6. Expand from crossing to define contiguous overlap region where q < q_max

Reuses load_latlon, reshape_to_bowties, haversine from modis_matchups.
"""

import numpy as np
from modis_matchups import (
    haversine, load_latlon, reshape_to_bowties,
    N_DETECTORS, N_PIXELS, R_EARTH
)


# ============================================================
# Constants
# ============================================================

# Detector pairs: (trailing D from B_k, leading D from B_{k+1}), 1-indexed
DETECTOR_PAIRS = [
    (di, dj)
    for di in range(7, 11)   # D=7,8,9,10 (trailing edge of B_k)
    for dj in range(1, 4)    # D=1,2,3   (leading edge of B_{k+1})
]

NADIR_COL = 677  # 1-indexed column at nadir


# ============================================================
# Projection
# ============================================================

def equirect_project(lat, lon, lat0, lon0):
    """Project lat/lon to local equirectangular (x_east, y_north) in meters."""
    x = np.radians(lon - lon0) * np.cos(np.radians(lat0)) * R_EARTH
    y = np.radians(lat - lat0) * R_EARTH
    return x, y


# ============================================================
# Step 1: Find intersection
# ============================================================

def find_intersection(lat_ref, lon_ref, lat_meas, lon_meas, cols,
                      poly_deg=8):
    """Find where two scan lines cross by fitting smooth polynomial curves.

    Fits degree-8 polynomials to the latitude of each scan line as a
    function of column index, then finds where the two fitted curves
    intersect (delta_lat_fit = 0).  This completely eliminates the
    ~111 m geolocation quantization zigzag that plagued earlier
    approaches, producing a single clean crossing per side of nadir.

    Parameters
    ----------
    lat_ref, lon_ref : 1D array (1354,), reference scan line (leading D)
    lat_meas, lon_meas : 1D array (1354,), measured scan line (trailing D)
    cols : 1D array of 0-indexed columns to search (one side of nadir)
    poly_deg : int, polynomial degree for the curve fits (default 8)

    Returns
    -------
    cross_col : float (0-indexed, interpolated) or None if no crossing
    """
    if len(cols) < poly_deg + 2:
        return None

    x = cols.astype(np.float64)

    # Fit smooth polynomials to lat(P) for each scan line
    p_ref = np.polyfit(x, lat_ref[cols].astype(np.float64), poly_deg)
    p_meas = np.polyfit(x, lat_meas[cols].astype(np.float64), poly_deg)

    # Delta polynomial: where meas_fit(P) - ref_fit(P) = 0
    p_delta = np.polysub(p_meas, p_ref)

    # Find roots of the delta polynomial
    roots = np.roots(p_delta)

    # Keep only real roots within the column range (with margin)
    margin = 10
    real_mask = np.abs(roots.imag) < 1e-6
    real_roots = roots[real_mask].real
    valid = real_roots[(real_roots >= cols[0] + margin) &
                       (real_roots <= cols[-1] - margin)]

    if len(valid) == 0:
        return None

    if len(valid) == 1:
        return float(valid[0])

    # Multiple crossings: pick the most interior one (furthest from edges)
    center = (cols[0] + cols[-1]) / 2.0
    return float(valid[np.argmin(np.abs(valid - center))])


# ============================================================
# Step 2: Fit reference line
# ============================================================

def fit_line(lat_ref, lon_ref, center_col, w=30):
    """Fit a straight line through w pixels of the reference scan line.

    Projects to a local equirectangular frame and fits y = slope*x + intercept.

    Parameters
    ----------
    lat_ref, lon_ref : 1D array (1354,)
    center_col : int, 0-indexed center of the window
    w : int, window width

    Returns
    -------
    dict with slope, intercept, lat0, lon0, cols, x_ref, y_ref
    """
    half = w // 2
    col_start = max(0, center_col - half)
    col_end = min(N_PIXELS - 1, center_col + half - 1)
    cols = np.arange(col_start, col_end + 1)

    # Projection center at window midpoint
    mid = len(cols) // 2
    lat0 = float(lat_ref[cols[mid]])
    lon0 = float(lon_ref[cols[mid]])
    x, y = equirect_project(lat_ref[cols], lon_ref[cols], lat0, lon0)

    # OLS line fit
    slope, intercept = np.polyfit(x, y, 1)

    return {
        'slope': slope,
        'intercept': intercept,
        'lat0': lat0,
        'lon0': lon0,
        'cols': cols,
        'x_ref': x,
        'y_ref': y,
    }


# ============================================================
# Step 3: Perpendicular distance
# ============================================================

def perp_distance(lat, lon, line):
    """Signed perpendicular distance (m) from points to the fitted line.

    Positive means above the line in the local projected frame.
    """
    x, y = equirect_project(lat[line['cols']], lon[line['cols']],
                            line['lat0'], line['lon0'])
    m, b = line['slope'], line['intercept']
    # Line equation: m*x - y + b = 0
    return (m * x - y + b) / np.sqrt(m**2 + 1)


# ============================================================
# Step 4: Pixel spacing
# ============================================================

def compute_pixel_spacing(lat, lon, cols):
    """Along-scan pixel spacing (m) at each column via Haversine.

    Centered differences in interior; one-sided at edges.
    """
    sp = haversine(lat[cols[:-1]], lon[cols[:-1]],
                   lat[cols[1:]], lon[cols[1:]])
    result = np.empty(len(cols))
    result[0] = sp[0]
    result[-1] = sp[-1]
    result[1:-1] = 0.5 * (sp[:-1] + sp[1:])
    return result


# ============================================================
# Steps 5-6: Quality metric and overlap region
# ============================================================

def find_region_one_side(lat_ref, lon_ref, lat_meas, lon_meas,
                         cols, w=30, q_max=0.01):
    """Find the overlap region on one side of nadir for one detector pair.

    Parameters
    ----------
    lat_ref, lon_ref : reference scan line (leading D from B_{k+1})
    lat_meas, lon_meas : measured scan line (trailing D from B_k)
    cols : 0-indexed columns for one side of nadir
    w : line-fit window width
    q_max : quality threshold (fraction of pixel spacing)

    Returns
    -------
    dict with region info, or None if no valid overlap found
    """
    # Step 1: intersection
    cross_col = find_intersection(lat_ref, lon_ref, lat_meas, lon_meas, cols)
    if cross_col is None:
        return None
    center = int(round(cross_col))

    # Step 2: fit reference line once
    line = fit_line(lat_ref, lon_ref, center, w=w)

    # Step 3: perpendicular distances
    d_perp = perp_distance(lat_meas, lon_meas, line)

    # Step 4: pixel spacing and quality metric
    win_cols = line['cols']
    ps = compute_pixel_spacing(lat_ref, lon_ref, win_cols)
    q = np.abs(d_perp) / ps

    # Step 5: contiguous region expanding from best-overlap column
    # Use the minimum-q column (best overlap) rather than the crossing
    # column, since delta_lat is only an approximate proxy for along-track
    # displacement and the crossing may not coincide with minimum q at
    # high latitudes or unusual orbit geometries.
    best_idx = int(np.argmin(q))
    mask = q < q_max

    if not mask[best_idx]:
        return None

    # Expand outward from the best-overlap column
    left = best_idx
    while left > 0 and mask[left - 1]:
        left -= 1
    right = best_idx
    while right < len(mask) - 1 and mask[right + 1]:
        right += 1

    overlap_cols = win_cols[left:right + 1]

    return {
        'cols': overlap_cols,                       # 0-indexed columns in overlap
        'cols_P': overlap_cols + 1,                 # 1-indexed
        'q': q[left:right + 1],                     # quality in overlap
        'd_perp': d_perp[left:right + 1],           # perpendicular distance (m)
        'pixel_spacing': ps[left:right + 1],        # pixel spacing (m)
        'cross_col': cross_col,                     # interpolated crossing (0-idx)
        'cross_col_P': cross_col + 1,               # 1-indexed
        'center_col': center,                       # integer center (0-idx)
        'line': line,                               # fitted line parameters
        # Full window diagnostics
        'q_full': q,
        'd_perp_full': d_perp,
        'ps_full': ps,
        'win_cols': win_cols,
    }


# ============================================================
# Main entry point
# ============================================================

def find_regions_one_boundary(lat_k, lon_k, lat_k1, lon_k1,
                              w=30, q_max=0.01):
    """Find all overlap regions at one bow-tie boundary.

    Parameters
    ----------
    lat_k, lon_k : (10, 1354) arrays for bow-tie B_k
    lat_k1, lon_k1 : (10, 1354) arrays for bow-tie B_{k+1}
    w : line-fit window width
    q_max : quality threshold

    Returns
    -------
    list of region dicts, each augmented with det_i, det_j, side
    """
    # 0-indexed column ranges for each side of nadir
    left_cols = np.arange(0, NADIR_COL)           # cols 0..676  (P=1..677)
    right_cols = np.arange(NADIR_COL, N_PIXELS)   # cols 677..1353 (P=678..1354)

    regions = []
    for det_i, det_j in DETECTOR_PAIRS:
        # Trailing D from B_k (measured against the fitted line)
        lat_meas = lat_k[det_i - 1]
        lon_meas = lon_k[det_i - 1]
        # Leading D from B_{k+1} (reference line is fit through this scan)
        lat_ref = lat_k1[det_j - 1]
        lon_ref = lon_k1[det_j - 1]

        for side, col_range in [('left', left_cols), ('right', right_cols)]:
            result = find_region_one_side(
                lat_ref, lon_ref, lat_meas, lon_meas,
                col_range, w=w, q_max=q_max)
            if result is not None:
                result['det_i'] = det_i   # trailing D from B_k
                result['det_j'] = det_j   # leading D from B_{k+1}
                result['side'] = side
                regions.append(result)

    return regions

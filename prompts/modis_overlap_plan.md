# MODIS Detector Overlap Matchup Analysis Prompt

---

## Scientific Objective

Use the natural crossovers of scan lines from adjacent 10-detector scan groups to characterize three key instrument parameters:

1. **Scan mirror reflectance characteristics** – Differences between the two sides (A and B) of the rotating scan mirror as they view the same Earth scene
2. **Detector-to-detector relative bias** – Systematic radiometric differences between detector pairs observing the same scene
3. **Detector-to-detector relative variance** – Random noise differences between detector pairs

The analysis leverages the bow-tie effect inherent to MODIS scanning geometry, where detectors viewing different parts of the scan are offset along the scan line but nearly perfectly aligned in the across-track direction. By comparing retrieved SST values from overlapping pixels, we can isolate and quantify instrument-related contributions to SST product uncertainty and detect mirror degradation over the mission lifetime.

---

## MODIS Dataset Geometry and Bow-Tie Effect

### Scan Configuration

MODIS acquires data with 10 detectors arranged along the cross-track direction. Detectors are organized into scan groups:
- **Scan group 1**: Detectors 1-10 (scan lines 1-10)
- **Scan group 2**: Detectors 1-10 (scan lines 11-20)
- **... continuing for ~4100 groups** (41,000 total scan lines per orbit)
- Each scan line is 1354 pixels in the along-scan direction

### Mirror Rotation

The MODIS scan mirror rotates continuously, alternating between two reflective surfaces:
- **Mirror side A**: Views and scans during rotations corresponding to groups 1, 3, 5, ... (odd-numbered groups)
- **Mirror side B**: Views and scans during rotations corresponding to groups 2, 4, 6, ... (even-numbered groups)

### Bow-Tie Effect and Scan Line Crossing Geometry

Due to the push-broom scanning geometry with a rotating mirror:
- Scan lines from consecutive 10-detector groups **cross one another** as they sweep across the Earth's surface
- At the crossing point, pixels from two different detector groups observe the **same geographic location** (within a few kilometers)
- **Across-track alignment**: Pixels are nearly perfectly co-registered in the across-track (nadir) direction (~km-scale precision)
- **Along-scan offset**: Due to the staggered detector arrangement and mirror geometry, pixels are offset by approximately **0.5 pixels in the along-scan direction** (~500 m at nadir)
- This offset is **consistent for each detector pair type** across the entire mission (e.g., detector 10 from group k always has a fixed offset relative to detector 1 from group k-1)

### Spatial and Temporal Structure

- **~2050 mirror rotations per orbit** (41,000 scan lines ÷ 20 lines per rotation)
- **~10-12 distinct detector pair crossovers per rotation** (pairs like 10-1, 10-2, 9-1, 9-2, 8-1, etc.)
- **~20,000-25,000 potential crossover matchups per orbit** (across all detector pairs)
- **100,000 orbits × 2.5 billion matchups = 2.5 billion matchup instances** available for analysis

---

## Matchup Definition and Quality Criteria

A valid **matchup** is defined as follows:

### Scan Selection
- **Scan i**: Any scan line whose detector number is in the range **1-5**
- **Scan j**: Any scan line from a **previous group** (different 10-detector group) such that the scan line j crosses the geographic location of scan line i

### Spatial Overlap Identification
For each scan i and potential scan j pair:
1. Identify the geographic crossing point where the two scan lines intersect (using latitude/longitude geolocation data)
2. Determine the range of pixels from scan i that participate in this intersection (typically ~20-50 pixels in the along-scan direction)
3. For each pixel k in scan i at position (lat_i(k), lon_i(k)):
   - Find the corresponding pixel in scan j at the same **distance from nadir** (same pixel index along the scan line)
   - Compute the distance from this scan j pixel to a **locally fitted straight line through scan i** at the same across-track position
   - **Keep pixel k in the overlap region** only if this distance is within a **user-specified tolerance** (e.g., 0.5%, 1%, or 2% of the pixel spacing)

### Quality Filtering
- **Exclude the entire matchup** if **any SST value** (from either scan i or scan j in the overlap region) is NaN or missing
- Only retain matchups with valid, non-NaN SST data in the complete overlap region

---

## Output Data Format

For each valid matchup, save the following information:

| Field | Description | Dimensions |
|-------|-------------|-----------|
| `scan_i` | Scan line number ranging from 1 to approximately 40000 | scalar |
| `scan_j` | Scan line number for detector group crossing | scalar |
| `detector_i` | Detector number in scan i (1-5) | scalar |
| `detector_j` | Detector number in scan j | scalar |
| `mirror_side_i` | Mirror side for scan i ('A' or 'B') | scalar |
| `mirror_side_j` | Mirror side for scan j ('A' or 'B') | scalar |
| `orbit_date` | Date of the orbit (YYYYMMDD or similar) | scalar |
| `orbit_number` | MODIS orbit identifier | scalar |
| `pixel_indices_from_nadir` | Array of pixel indices (distance from nadir) for pixels in overlap | vector, ~30-50 elements |
| `latitude` | Latitude of each overlap pixel | vector |
| `longitude` | Longitude of each overlap pixel | vector |
| `SST_i` | Retrieved SST values from scan i for overlap pixels | vector |
| `SST_j` | Retrieved SST values from scan j for overlap pixels | vector |
| `offset_pixels_along_scan` | The fixed along-scan pixel offset for this detector pair type | scalar |
| `overlap_quality` | Quality metric indicating how well pixels are aligned (distance to fit line / pixel spacing) | vector or scalar |

All data should be organized into a structured array or table, with each row representing one valid matchup instance.

---

## Code Request

Please write a MATLAB script that:

1. **Loads MODIS L2 data** (SST product) for a single orbit or multiple orbits, along with the corresponding geolocation (latitude, longitude)

2. **Pre-computes detector pair offsets** (one-time, from representative orbits):
   - For each detector pair type that crosses (10-1, 10-2, 9-1, etc.), determines the fixed along-scan pixel offset using geolocation data

3. **Extracts all valid matchups** by:
   - Iterating through all scan groups
   - For each scan i with detector ∈ [1,5], identifying all previous scans j that cross it
   - Identifying the pixel range where the scan lines intersect geographically
   - Filtering pixels based on distance to a locally fitted straight line (within user-specified tolerance)
   - Applying the NaN filtering criterion

4. **Saves matchup data** in a structured format (MATLAB table, array of structs, or HDF5) containing all fields listed above

5. **Provides progress reporting** showing:
   - Number of orbits processed
   - Number of matchups extracted per orbit
   - Number of matchups rejected due to NaN filtering
   - Estimated total matchup count

6. **Handles edge cases**:
   - Orbits with fewer than expected scan lines
   - Missing geolocation data
   - Scan lines at swath edges (where geometry may be different)

The code should be modular, well-commented, and structured for efficient processing of the full 100,000-orbit dataset (with options for testing on smaller subsets first).

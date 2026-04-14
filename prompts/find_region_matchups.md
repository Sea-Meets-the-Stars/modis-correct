# Find bow-tie match-ups in regions

Be sure to re-read the CLAUDE.md file to understand the nomenclature and goals
of the overall project.


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


## Output Data Format

For each valid matchup, save the following information:

| Field | Description | Dimensions |
|-------|-------------|-----------|
| `scan_i` | Scan line number ranging from 1 to approximately 40000 | scalar |
| `scan_j` | Scan line number for detector group crossing | scalar |
| `detector_i` | Detector number in scan i (1-5) | scalar |
| `detector_j` | Detector number in scan j | scalar |
| `bowtie_i` | Bow-tie number for scan i | scalar |
| `bowtie_j` | Bow-tie number for scan j | scalar |
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

## Data

### Development

There is a MatLab data file in modis-correct/dev/matchups/data/ named latlon.mat which contains the latitude and longitude of the center of each pixel in the orbit.

## Code

### Existing code

When possible, use any existing code in the modis-correct/dev/matchups/py directory.

### Writing

Here are guidelines for writing code:

- Add inline comments to explain the effort
- Reuse existing code when possible
- Use methods, not classes

## Overleaf

Place any Latex files in /home/xavier/Projects/overleaf/modis-correct/

The planning file for this activity is find_regions_plan.tex.  Make it a standlone LateX file.

Output PNG figures to the Overleaf project folder.

You may push to git as you work.  The access token is in my .bashrc profile with the name OVERLEAF.

## Requirements

## Development phase

1. Generate an algorithm to implement Approach B in our plan.  Use qmax=1% for now.
Put any new methods into the modis-correct/dev/matchups/py directory in a module named find_regions.py. 
For testing purposes, generate an IPython Notebook that operates on
one pair of bow-ties using the data in the modis-correct/dev/matchups/data/latlon.mat file.  
Name the Notebook Test_region.ipynb and put it in modis-correct/dev/matchups/nb.  Include figures that demonstrate the algorithm works.

2. Please make the modifications to the Notebook

- Use q_max = 5% instead of 1%
- For Step 1 in the Notebook(scan-line intersection), add a figure that shows the two scan lines intersecting lat/lon space.  Zoom-in as appropriate.
- For Step 2 in the Notebook, add a separate figure showing the fitted line on the scan line.  Zoom in as appropriate

3.  Generate a MATLAB file that performs the same functionality as provided by the Python code.  And create a test script that mirrors the IPython Notebook.  Put those files in 
modis-correct/dev/matchups/mat.

4. The algorithm you generated for finding the intersecton -- find_intersection() -- is not sufficiently accurate and too unstable to small fluctuations in the scanline lat/long values.  Please generate a new version.  Generate a new module called test_regions.py in py/ and add a method that shows one or more figures describing the intersection.

5. The intersection algorithm is still flawed.  It would be best to fit a smooth curve to each scan line and then find the intersection of the two curves.  Try this approach and update the find_intersection() method and diagnostic figures.

6. Ok, that is better.  But you will need to fit the scanlines in lat vs. lon, not lat vs. column P.  Then find where those 2 curves intersect and figure out the best P to the nearest integer. 

## Prompts

### Brainstorming

1. Read this file and brainstorm several approaches to generating an algorithm that meets algorithmic goals.  Put these in find_regions_plan.tex.  Push to Overleaf.
2. Make the following modifications to Approach B:

- Step 2 should find the interasection (if any) between the two scan lines on one side of nadir and then the other.  
- Use w=30 and center on the intersection point.
- Fit the line only once, not in a sliding window.

### Development

1. Re-read this file and generate the code to implement the first description under the Development phase above.
2. Re-read this file and perform item 2 of the Development phase.
3. Re-read this file and perform item 3 of the Development phase.
4. Re-read this file and perform item 4 of the Development phase.
5. Re-read this file and perform item 5 of the Development phase.
6. Re-read this file and perform item 6 of the Development phase.
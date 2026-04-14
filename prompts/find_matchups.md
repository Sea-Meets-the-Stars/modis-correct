# Finding the bow-tie match-ups

Be sure to re-read the CLAUDE.md file to understand the nomenclature and goals
of the overall project.

## Goals

Here are the goals of this task:

- Find all match-ups between bow-ties in a given orbit that have a separation distance in meters less than an input value
- Write a code that is efficient in compute and memory

## Data

### Development

There is a MatLab file in modis-correct/dev/matchups/data/ named latlon.mat which contains the latitude and longitude of the center of each pixel in the orbit.

## Code

## Advice

You will find that you will only need to compare adjacent bow-ties in the orbit for the match-ups.

Match-ups will only occur in the same column of the orbit, i.e. the same number of pixels from nadir. 

### Writing

Here are guidelines for writing code:

- Add inline comments to explain the effort
- Reuse existing code when possible
- Use methods, not classes

## Overleaf

Place any Latex files in /home/xavier/Projects/overleaf/modis-correct/

The planning file for this activity is find_matchups_plan.tex.  Make it a standlone LateX file.

Output PNG figures to the Overleaf project folder.

You may push to git as you work.  The access token is in my .bashrc profile with the name OVERLEAF.

## Requirements

These requirements implement Approach 2 (Column-Wise Vectorized Haversine) from our plan.

### Inputs

1. Path to a single MODIS L2 orbit file containing latitude and longitude arrays of shape (N_lines, 1354)
2. Maximum separation distance `d_max` in meters
3. Candidate table: a pre-computed structure mapping pixel index P to a list of detector pairs (D_i, D_j) to compare at that column

### Data Loading

4. Read the latitude and longitude arrays from the data file (float32)
5. Handle the partial first bow-tie: the first scan line of each orbit has D=6, so the first bow-tie contains only detectors D=6-10 (5 lines). Discard or flag this partial bow-tie — do not treat it as a complete 10-detector group.
6. Group the remaining scan lines into complete bow-ties of 10 consecutive lines, each corresponding to detectors D=1-10
7. Handle a possible partial bow-tie at the end of the orbit (if N_lines is not evenly divisible after accounting for the offset)
8. Reshape the latitude and longitude arrays into bow-tie arrays of shape (N_bt, 10, 1354)

### Match-Up Computation

9. Iterate over each candidate detector pair (D_i, D_j) from the candidate table
10. For each pair, extract the latitude/longitude of detector D_i in bow-tie B_k and detector D_j in bow-tie B_{k+1}, for all adjacent bow-tie boundaries and all candidate columns P simultaneously — these are 2-D arrays of shape (N_bt - 1, N_P)
11. Compute the Haversine great-circle distance between pixel centers using vectorized NumPy operations across both bow-tie boundaries and columns in a single call
12. Identify all pairs where d < d_max and record the match-ups

### Mirror Side Tracking

13. Assign mirror sides M=1 and M=2 to alternating bow-ties (the first complete bow-tie gets M=1, the next M=2, etc.)
14. Record the mirror side for both pixels in each match-up

### Output

15. For each match-up, record:
    - Orbit identifier (filename or orbit number)
    - Bow-tie index k (of the first pixel)
    - Pixel column P
    - Detector D_i (in B_k) and D_j (in B_{k+1})
    - Separation distance d in meters
    - Mirror sides M_k and M_{k+1}
16. Return results as a structured array or DataFrame, sorted by bow-tie index then column

### Performance

17. Use vectorized NumPy operations throughout — no Python-level loops over columns or bow-tie boundaries
18. The only Python loop should be over candidate detector pairs (~6-12 iterations)
19. Target processing time: < 10 seconds per orbit on a single core
20. Peak memory: < 1 GB per orbit (the lat/lon arrays are ~200 MB each in float32)
21. Code should be structured to allow future parallelization across orbits

### Candidate Table

22. The candidate table can be built empirically from Approach 3 (brute-force) on a representative set of orbits, or analytically from MODIS geometry
23. Use a generous margin (e.g. 1.5 x d_max) when building the table to avoid missing edge-case match-ups
24. The table maps each column P to the list of (D_i, D_j) pairs to check, where D_i is in {7,...,10} and D_j is in {1,...,3}
25. Exploit nadir symmetry: the candidate pairs at column P are the same as at column 1355 - P

### Validation

26. Provide a method to run Approach 3 (brute-force: all 12 detector pairs at every column) on a single orbit and compare results against Approach 2, confirming no match-ups are missed

## Development phase

### Notebooks

1. Generate a Python Notebook that does the following:

- Load the data from the data/latlon.mat file
- Given a bowtie number B, plot the lat/lon of that bowtie for each detector as one color and those of the next bowtie B+1 as a different color.
- Calculate the haversine distances, restricted to columns between the two bowties and generate a historgram of these.
- As possible, use existing code in modis_matchups.py
- Name the Notebook Explore_bowtie.ipynb and put it in modis-correct/dev/matchups/nb/


## Prompts

### Brainstorming

1. Read this file and brainstorm several approaches to generating an algorithm that meets the task goals.  Put these in find_matchups_plan.tex.  Push to Overleaf.
2. We have revised the goals of the task to focus on separation distance in meters.  Revise the plan to reflect the new goals.  Push to Overleaf.  We will use the latitude and longitude of the center of each pixel to calculate the separation distance in meters.
3. Consider the advice that I have added above and refine the plan.  Push to Overleaf.

### Requirements

1. Please generate a set of requirements using Approach 2 in our plan.  Add those to the Requirements section above.
2. Modify the requirements to reflect that the search for match-ups only need to occur within a given column and one must search over all columns.
3. Revert back to vectorize over columns.

### Development

1.  Re read this doc. Generate the first Notebook described in Development phase/Notebooks above.

### Coding

1. Reread the doc. Generate the code to meet the requirements.  Push to the dev/matchups/py directory. Call it modis_matchups.py.  Use a portion of the data in modis-correct/dev/matchups/data/ to test the code.  Generate figures to show the match-ups.
2. Add a new script -- py/test_modis_matchups.py -- to test the code.  It should input a d_max value and run the code and place figures and results in the dev/matchups/outputs directory.  Name the files according to d_max.
3. In the test outputs, I see many examples with 0m separation.  Check whether these are in error.  Such small values are unlikely unless you are comparing the same pixel.

### Testing

1. I am seeing many tens of matchups for bowtie 4021 in the matchups_dmax40.csv table in outputs/  This seems to be too many.  Please explore.
2. Generate a figure showing the Column and row of the matchups between detector 10 and 1.  Row is defined as (bowtie number - 1) * 10 + detector number.  Use the outputs in outputs/matchups_dmax40.csv.  Save the code to generate the figure in the dev/matchups/py/testing_figs.py.  Color the dots according to the matchup separation.
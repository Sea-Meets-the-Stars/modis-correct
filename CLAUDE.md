# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Interactions

- Do not make up data
- Talk to me directly
- Be concise and to the point
- Be critical of my requests and your own work

## Project Overview

**modis-correct** is a Python and MatLab package for analyzing overlapping bowties of MODIS satellite data, specifically sea surface temperature (SST).

## Primary context

I have 100,000 satellite MODIS L2 orbits, one orbit/file. The orbits are 40,000 scan lines long and 1354 pixels/line. MODIS is a 10 detector scanner hence suffers from the bow-tie effect.  As the satellite orbits, it creates a bow-tie of 10 scan lines using one side of a rotating mirror and then another 10 with the other side of the mirror.

I want to make use of the overlap to determine the noise associated with some of the detectors. Because of the bow-tie effect, at a given number of pixels from nadir along the line, the pixel of the 10th detector in one group of 10 pixels overlaps the same pixel in the 1st scan line of the next group to within a small number of meters.

Going farther from nadir, the pixel at that location on the 10th detector overlaps the corresponding pixel of the 2nd detector. At the swath edge the 10 detector overlaps the 7th. I would like you to write a Python script to first find the number of pixels from nadir for the 10-1 overlap with x m, the 10-2 overlap, the 10-3 overlap, the 10-4 overlap if there is one, the 9-1 overlap, 9-2 overlap, 9-3 overlap, 8-1, 8-2, 8-3, 7-1 and 7-2. 

The next step would be to generate a dataset of orbit number, scan group--there are about 4000 groups of 10 scan lines--and there are probably 12  overlaps, 6 on each side of nadir. For each orbit there will about 50,000 overlaps. I would want you to write the script and test it on some of my data, generating test scrips along the way. 

When we are convinced that it is working properly, I then like you to run the whole thing; i.e., on all 100,000 orbits. This this information I can get great stats on the noise in some of the detectors and how it is changing in time, to some extent along the scan line and in each orbit.

# Data

## Coordinates

The latitude and longitude of the center of each pixel is stored in the data files.

# Nomenclature

## Detectors

There are 10 detectors in the MODIS scanner.
We will label these D=1-10

## Bow-ties

A bow-tie is a group of 10 scan lines, one from each detector. 
These occur sequentially during an orbit.

For each orbit, the first scan line has D=6.

## Mirror

There are two sides to the mirror.  
We will label these M=1 and M=2.
Each bow-tie uses either M=1 or M=2.

## Scan-lines

Each scan line from a detector has 1354 pixels.  
We will refer to these as P=1-1354.

## Match-up

A match-up is two pixels from separate bow-ties that overlap spatially.  They are defined by their separation distance in meters.

## Overlap

The overlap is the area of the two pixels that overlap relative to the mean area of the two pixels.  Overlap=1 is perfect overlap.

# Code guidelines

- Reuse existing code when possible

## Python

Adhere to following:

- Generate methods when possible, not classes
- Include inline comments
- Use matplotlib for plotting

# Overleaf

Place any Latex files in /home/xavier/Projects/overleaf/modis-correct/

Output PNG figures to the Overleaf project folder.

You may push to git as you work.  The access token is in my .bashrc profile with the name OVERLEAF.

Embed the figures within the text.


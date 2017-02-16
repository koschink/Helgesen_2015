# Helgesen_2015
A script to measure distances of replication foci in E.coli using ImageJ and Python. 
For a more complete description, see our publication:
Dynamic Escherichia coli SeqA complexes organize the newly replicated DNA at a considerable distance from the replisome
Helgesen et al, 2015 
(https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4357733/)

This script perfomes automatic measurements of the distance between neighboring spots that are registered in two different fluorescence channels. The script utilizes the “Find maxima” function to detect spots and then measures the distance to the nearest neighbouring spot. The resulting output can give distances
that are below the resolution limit. This is rational since the neighbouring spots are located in different channels (3-5). 

The input to the script is a composite image with two fluorescence channels, a list of cell outlines (regions-of-interest) and noise-levels for the “Find Maxima” function applied by the script to detect centres-of-mass. The script then measures the distance from the peak of every spot detected in the first channel to the nearest peak of spots detected in the second channel, per cell outline. It also provides a per-cell count of number of spots, as well as their position relative to the cell centre.

The measured distances can also be used to estimate the level of colocalization (object-based olocalization). Objects (or foci) will be considered colocalized if the distance separating them is ithin a reference distance defined by the optical resolution of the system. The reference distance in our analyses was set to the resolution limit based on longest emission wavelength of the two fluorophores investigated. Cell length was estimated through Fit Ellipse. Fit Ellipse is a function that replaces a cell outline
with a best fit elliptical shape, where the major axis-length of this shape corresponds to cell
length.Detailed information about the script code can be found in the script file as comments in green
text, prefixed with hash tags (#)

For further information or comments, please contact the author (kay.oliver.schink@rr-research.no)

You are free to use this software in your research under the BSD3 license. If you find this software useful, please cite our publication

[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

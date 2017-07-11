# Single-Particle-Raster-Image-Analysis
-------------------------------------
This repository contains a Matlab program for Single Particle Raster Image Analysis as described in:

Longfils, M., Schuster, E., Lorén, N., Särkkä, A. & Rudemo, M. (2017) Single particle
raster image analysis of diffusion. Journal of Microscopy, 266, 3–14.

Longfils, M., Röding, M., Altskär, A., Schuster, E., Lorén, N., Särkkä, A. & Rudemo, M. (2017) Single particle
raster image analysis of diffusion for particle mixtures. Under revision.


Requirements
---------------------------------------
MATLAB (8.1)
Image Processing Toolbox (8.2)
Optimization Toolbox (6.3)

Introduction
---------------------------------------
These codes perform the Single Particle Raster Image Analysis as described
in the paper listed above. In particular, it estimates the diffusion coefficients
of single particles in a raster image collected with a CLSM (confocal laser scanning
microscope) and fit mixture models with up to 4 diffusing species.

The result will be printed in the .txt file called "run_SPRIA" and the postscript
"SPRIA_fit" and "SPRIA_scatterplot"

Installation
---------------------------------------
To install add the folder "Single Particle Raster Image Analysis" to the Matlab path. To run
this software load the desired data and type "extraction" at the Matlab prompt. Then type 
"SPRIA_analysis" at the Matlab prompt and open the file "run_SPRIA.txt" with a text editor.

File Support
---------------------------------------
Supports data saved in a .mat file containing a structure "A" with at least the following fields:
- Imgs : a MxMxN matrix containing N raster images with resolution MxM pixels
- Sx : pixel size in the x direction
- Sy : pixel size in the y direction
- Tl : line dwell time used to collect the raster images
- Tp : pixel dwell time used to collect the raster images


Algorithm
---------------------------------------
The centroid method is used to estimate the position of the particles. For
model selection and estimation in case of mixture models the maximum likelihood
method is used.

Contributors
---------------------------------------
Marco Longfils
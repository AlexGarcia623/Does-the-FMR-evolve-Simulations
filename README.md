# Does the FMR evolve with redshift? A Perspective From Cosmological Simulations

This repo contains the analysis scripts, figures, and data products from the paper Garcia et al. (2024b)

## Contents

### Data:

Contains minimal working example data for generation of the figures associated with this work for Illustris, IllustrisTNG, and EAGLE simulations from redshift 0-8 (technically 0-10, though none of 9 and 10 are used in this work).

All credit for simulation data goes to the respective collaborations: see [Illustris](https://www.illustris-project.org/), [IllustrisTNG](https://www.tng-project.org/), and [EAGLE](https://icc.dur.ac.uk/Eagle/).

### Data Reduction

!!! Not required to run the scripts to generate plots !!!
!!! Is required if you would like to significantly alter selection criteria !!!

Contains scripts used to generate data products from "raw" simulation data 

### Figures (pdfs)

Contains pdf versions of all figures (exception of Figure 1)

### Rest of files

- appendix_A.py -- Generate Figure A1
- appendix_B.py -- Generate Figure B1
- Figure1.tex -- Generate Figure 1
- Figure2.py -- Generate Figure 2
- Figure3.py -- Generate Figure 3

- getAlpha.py -- Contains most of the logic for generating Figures 1, 2, A1, and B1
- helpers.py -- Contains useful functions for getAlpha.py

## Questions Regarding scripts/paper

Please contact [alexgarcia@virginia.edu](mailto:alexgarcia@virginia.edu), Alex Garcia PhD Student University of Virginia
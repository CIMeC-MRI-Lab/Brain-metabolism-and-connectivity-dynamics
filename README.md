# Brain-metabolism-and-connectivity-dynamics

## Description
This repository contains the scripts used for the study _"Multiscale excitation-inhibition balance dynamics: integrating metabolite kinetics with time-varying executive networks"_.
It implements the analysis of novel interleaved dynamic functional MRI and functional edited-MRS to investigate temporal variation in neurometabolites and neurovascular compartments as a function of cognitive load increase. 

Sample code to reproduce results of the manuscript "X" by Saviola et al., bioRxiv 2024. https://biorxiv.org/cgi/content/short/2024.10.30.621153v1

The code generates Fig. 1 of the paper (together with Supplementary Figures) and performed the analysis for EIB kinetics detection (i.e. _"EIB_kinetics"_) and the statics on behavioural performance (i.e. _"Behav"_) and its respective interaction with EIB and FPN-CAPs (i.e. _"Behav_interactions"_). 
More comment and details are provided within the main scripts readme (i.e., "README.md"). 
Some sample data (i.e. mainly behaviour) are included in the repo to make the code standalone.  

Please note that the code for functional and structural image preprocessing can be downloaded from: [Lnifmri_prep](https://github.com/tambalostefano/lnifmri_prep)  

Please note that the code for the dynamic connectivity analysis (Fig. 2) can be downloaded from: [TbCAPs](https://github.com/MIPLabCH/TbCAPs)  

Please note that the code for Partial Least Square Correlation (Fig. 3) can be downloaded from: [myPLS](https://github.com/MIPLabCH/myPLS)


Code Author: Francesca Saviola, Stefano Tambalo, Barbara Cassone & Asia Ferrari

Any comments/queries can be sent to: francesca.saviola@epfl.ch

version 1.1 (November, 2024)

Please cite us! 

Francesca Saviola, Stefano Tambalo, Laura Beghini, Asia Ferrari, Barbara Cassone, Dimitri Van De Ville, Jorge Jovicich
"Multiscale excitation-inhibition balance dynamics: integrating metabolite kinetics with time-varying executive networks", bioRxiv 2024.10.30.621153
https://biorxiv.org/cgi/content/short/2024.10.30.621153v1


## Requirements
- MATLAB R2020b or later
- R 4.0.0 or later
- Bash shell (Unix-based system or Windows Subsystem for Linux)
- SPM12 toolbox (https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) for MATLAB
- Gannet (https://markmikkelsen.github.io/Gannet-docs/index.html)
- gramm (https://github.com/piermorel/gramm)
- R packages: dplyr, lme4, ggplot2, knitr, kableExtra, stats, tidyr, pacman, reshape2, ez, lmerTest, rmarkdown, lattice, DACF (install via `install.packages(c("dplyr", "lme4", "ggplot2", "knitr", "kableExtra", "stats", "tidyr", "pacman", "reshape2", "ez", "lmerTest", "rmarkdown", "lattice", "DACF"))`)

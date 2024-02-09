

The R package **calfsbm** provides functionality for the CALF-SBM model 
in network analysis:

1. Allows the user to generate simulated networks, with customizable size, 
density, and signal strength of covariates;
2. Run MCMC chains for both real-world and simulated networks with help from
the NIMBLE package (Note: the chains run slowly and may need to be run on a 
high-performance computer);
3. Post-processing of MCMC samples, accounting for label-switching

# Installation
One may install the current development version of this package with the 
following command in R:
```R
devtools::install_github('sydneylouit539/calfsbm')
```

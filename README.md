# MAC-PD Treatment Strategy Model

[![DOI](https://zenodo.org/badge/914634018.svg)](https://doi.org/10.5281/zenodo.17989200)

This repository provides the complete code to run the simulation and analysis in the following paper:

[Arbiv OA, Postill G, Im Y, Jhun BW, Pechlivanoglou P, Brode SK, Bayoumi AM, Marras TK. Antimicrobial Treatment of Mild Mycobacterium avium Complex–Pulmonary Disease Predicted to Increase Survival and Quality Adjusted Life Years: A Microsimulation Decision Analysis Model. CHEST. 2026 Mar 7.](https://journal.chestnet.org/article/S0012-3692(26)00280-1/)

The code is organized in the following way:

- `microsim.R` is a function containing the main simulation.
- `microsim_psa.R` is the script for the probabilistic analysis (primary outcome).
- `microsim_dsa.R` is a script that analyzes the base case scenario and sensitivity analyses.

Helper functions and baseline data are found in the `Functions` and `values` folders, respectively.

I welcome your contributions and comments.
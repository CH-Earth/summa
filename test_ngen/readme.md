# SUMMA case studies
This folder contains a case study to show how a typical SUMMA setup looks in the NextGen setup.The folder serves a double purpose as a way to track default versions of certain input files, such as the Noah-MP tables and spatially constant parameter files. These files are in /settings/meta.


## Settings
The meta folder contains the setting files that typically do not change for different model applications. Currently these include:
- `GENPARM.TBL`: lookup table for general parameters (legacy, currently unused)
- `MPTABLE.TBL`: lookup table for vegetation parameters
- `SOILPARM.TBL`: lookup table for soil parameters
- `VEGPARM.TBL`: lookup table for vegetation parameters

## Case studies 
For other (non-NextGen) examples the user is directed to clone the Laugh-Tests github repository at https://git.cs.usask.ca/numerical_simulations_lab/hydrology/laugh_tests
# SUMMA case studies
This folder contains a case study to show how a typical SUMMA setup looks. The folder serves a double purpose as a way to track default versions of certain input files, such as the Noah-MP tables and spatially constant parameter files. 

For more details about SUMMA inputs, see the [online documentation](https://summa.readthedocs.io/en/latest/input_output/SUMMA_input/).


## Base settings
This folder contains the setting files that typically do not change for different model applications. Currently these include:
- `GENPARM.TBL`: lookup table for general parameters (legacy, currently unused)
- `MPTABLE.TBL`: lookup table for vegetation parameters
- `SOILPARM.TBL`: lookup table for soil parameters
- `VEGPARM.TBL`: lookup table for vegetation parameters
- `basinParamInfo.txt`: default values for GRU(basin)-level parameters
- `localParamInfo.txt`: default values for HRU(local)-level parameters


## Case studies 
Instead of providing an example in the code, the user is directed to clone the Laugh-Tests github repository at https://git.cs.usask.ca/numerical_simulations_lab/laugh_tests.git.
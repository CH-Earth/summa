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


## Reynolds Mountain East
Contains input files needed to reproduce Figure 6a in Clark et al. (2015):
- Folder `forcing` contains forcing data, surprisingly
- `fileManager.txt`: main configuration file SUMMA expects as a command line argument
- `modelDecisions.txt`: specifies which modeling options (parametrizations, numerical method, etc.) to activate
- `forcingFileList.txt`: list of files in the forcing folder
- `outputControl.txt`: specifies which variables to output and, optionally, at which aggregation level
- `initialState.nc`: initial model states on a per-HRU basis
- `localAttributes.nc`: domain settings on a per-HRU basis
- `trialParams.nc`: experiment-specific parameter values on a per-HRU basis (can be empty if default parameters and lookup table values are used)


## References
 * Clark, M. P., B. Nijssen, J. D. Lundquist, D. Kavetski, D. E. Rupp, R. A. Woods, J. E. Freer, E. D. Gutmann, A. W. Wood, D. J. Gochis, R. M. Rasmussen, D. G. Tarboton, V. Mahat, G. N. Flerchinger, D. G. Marks, 2015: A unified approach for process-based hydrologic modeling: Part 2. Model implementation and case studies. _Water Resources Research_, [doi:10.1002/2015WR017200](http://dx.doi.org/10.1002/2015WR017200).
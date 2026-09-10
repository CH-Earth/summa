# SUMMA Output Files

<a id="outfile_file_formats"></a>
## Output file formats
All SUMMA output files are in [NetCDF format](SUMMA_input#infile_format_nc).

<a id="outfile_dimensions"></a>
## Output file dimensions
SUMMA output files can have the following dimensions (as defined in `build/source/netcdf/def_output.f90`). Dimensions may be present even in output files where they are not actually used. Most of these dimensions are pretty self-explanatory, except perhaps the `[mid|ifc][Snow|Soil|Toto]` dimensions, which are for depth information. The dimensions indicated by `ifc` are associated with variables that are specified at the interfaces between layers including the very top and bottom. For example, the flux into or out of a layer would be arranged along an `ifc` dimension. The dimensions indicated by `mid` are associated with variables that are specified at the mid-point of each layer (or layer-average). `Snow`, `Soil`, `Glce`, `Lake`, and `Toto` indicate snow layers, soil layers, glacier ice layers, lake layers and all layers, respectively.

| Dimension | long name | notes |
|-----------|-----------|-------|
| gru       | dimension for the GRUs | Variables and parameters that vary by GRU |
| hru       | dimension for the HRUs | Variables and parameters that vary by HRU |
| dom       | dimension for domain | Variables and parameters that vary by domain |
| glac      | dimension for the number of glaciers | Variables and parameters that vary by glacier
| depth     | dimension for soil depth | Variables and parameters that are defined for a fixed number of layers |
| scalarv   | dimension for scalar variables | Scalar variables and parameters (degenerate dimension) |
| spectral  | dimension for the number of spectral bands | Variables and parameters that vary for different spectral regimes |
| time      | dimension for the time step | Time-varying variables and parameters |
| tdh       | dimension for the time delay routing vectors | Variables and parameters that are held in memory as part of routing routines |
| midSnow   | dimension for midSnow | Variables and parameters at the mid-point of each snow layer |
| midSoil   | dimension for midSoil | Variables and parameters at the mid-point of each soil layer |
| midGlce   | dimension for midGlce | Variables and parameters at the mid-point of each glacier ice layer |
| midLake   | dimension for midLake | Variables and parameters at the mid-point of each lake layer (un-used currently)|
| midToto   | dimension for midToto | Variables and parameters at the mid-point of each layer in the combined layer profile |
| ifcSnow   | dimension for ifcSnow | Variables and parameters at the interfaces between snow layers (including top and bottom) |
| ifcSoil   | dimension for ifcSoil | Variables and parameters at the interfaces between soil layers (including top and bottom) |
| midGlce   | dimension for midGlce | Variables and parameters at the interfaces between glacier ice layers |
| midLake   | dimension for midLake | Variables and parameters at the interfaces between lake layers (un-used currently)|
| ifcToto   | dimension for ifcToto | Variables and parameters at the interfaces between all layers in the profile (including top and bottom) |
| grid      | dimension for the grid | Variables and parameters that vary by grid
| xgrid     | dimension for the x direction of the grid | Variables and parameters by x direction of grid
| ygrid     | dimension for the y direction of the grid | Variables and parameters by y direction of grid

<a id="outfile_restart"></a>
## Restart or state file
A SUMMA restart file is in [NetCDF forma](SUMMA_input#infile_format_nc) and is written by `build/source/netcdf/modelwrite.f90:writeRestart()`. This file is also an input file because it specifies the initial conditions at the start of a model simulation. It is described in more detail in the [SUMMA input](SUMMA_input#infile_initial_conditions) documentation. Note that when the file is written, the time for which it is valid is included as part of the model file name.

<a id="outfile_history"></a>
## Model history files
SUMMA history files are in [NetCDF format](SUMMA_input#infile_format_nc) and describe the time evolution of SUMMA variables and parameters. The files are written by the `writeParam`, `writeData`, and `writeTime` subroutines in `build/source/netcdf/modelwrite.f90`. SUMMA output is pretty flexible. You can output many time-varying model variables and parameters, including summary statistics. You can specify what you want to output in the [output control file](SUMMA_input#infile_output_control), which is one of SUMMA's required input files.

# SUMMA Input Files

SUMMA has a large number of input files that configure the model and provide the necessary initial conditions and time-varying boundary conditions to make a model simulation. This can at times be confusing. We encourage the user to look at the [SUMMA test cases](../installation/SUMMA_test_cases.md), which provide working SUMMA setups.

## Input file formats
<a id="infile_file_formats"></a>
SUMMA input files are either ASCII format or NetCDF. The general characteristics of these files are described in the next two subsections, while the contents of the individual input files are described after that.

### ASCII
<a id="infile_format_ASCII"></a>
ASCII or text files are in a format that can be modified using a text editor. Comments can be added to any SUMMA text file by starting the comments with  a `!`. Anything after the `!` will be discarded till the end of the line. You can include as many comments as you want, as they will be stripped as SUMMA processes the file.


### NetCDF
<a id="infile_format_nc"></a>
[NetCDF](https://www.unidata.ucar.edu/software/netcdf/) or Network Common Data Format is a file format that is widely used in geosciences to organize large data sets. The main advantages of using NetCDF files is that they are machine independent, they allow the user to include meta data directly in the data file, and they can be read by, visualized and analyzed using a large number of freely available software packages. The SUMMA documentation is not the place to learn about NetCDF. We assume that you know the difference between NetCDF dimensions, NetCDF variables, and NetCDF attributes (global and local). If you don't, then there are many tutorials available online. Note that the latter are different from the local SUMMA attributes that we are describing [below](#infile_local_attributes).


## Master configuration file
<a id="infile_master_configuration"></a>

The master configuration file is an [ASCII file](#infile_format_ASCII) and is provided to SUMMA at run-time as a command-line option. The path to this file needs to be supplied with the `-m` or `--master` command-line flag. The contents of this file orchestrate the remainder of the SUMMA run and are processed by the code in `build/source/hookup/summaFileManager.f90`. The file contents mostly consist of file paths that provide the actual information about the model configuration.

The following items must be provided in order in the master configuration file. Each item must be on its own line, but may be followed by a comment and you can add lines of comments between the items. Each entry must be enclosed in single quotes `'entry'`. In the following, I start each enumerated entry with the actual variable name that is used in the SUMMA source code to refer to each of the entries (in `summaFileManager.f90`) and its default value in case you are trying to trace a problem.

1. `summaFileManagerHeader`: Version of the file manager that should be used to process the master configuration file. At this time, this string should be equal to `'SUMMA_FILE_MANAGER_V1.0'`.

1. `SETNGS_PATH`: Base path for the configuration files. Most of the file paths in the remainder of the master configuration file are relative to this path (except `INPUT_PATH` and `OUTPUT_PATH`).

1. `INPUT_PATH`: Base path for the meteorological forcing files specified in the `FORCING_FILELIST`.

1. `OUTPUT_PATH`: Base path for the SUMMA output files.

1. `M_DECISIONS`: File path for the [model decisions file](#infile_model_decisions) (relative to `SETNGS_PATH`).

1. `META_TIME`: No longer used as of SUMMA 2.x - simply specify `'N/A'`.

1. `META_ATTR`: No longer used as of SUMMA 2.x - simply specify `'N/A'`.

1. `META_TYPE`: No longer used as of SUMMA 2.x - simply specify `'N/A'`.

1. `META_FORCE`: No longer used as of SUMMA 2.x - simply specify `'N/A'`.

1. `META_LOCALPARAM`: No longer used as of SUMMA 2.x - simply specify `'N/A'`.

1. `OUTPUT_CONTROL`: File path for the [output control file](#infile_output_control) (relative to `SETNGS_PATH`).

1. `META_LOCALINDEX`: No longer used as of SUMMA 2.x - simply specify `'N/A'`.

1. `META_BASINPARAM`: No longer used as of SUMMA 2.x - simply specify `'N/A'`.

1. `META_BASINMVAR`: No longer used as of SUMMA 2.x - simply specify `'N/A'`.

1. `LOCAL_ATTRIBUTES`: File path for the [local attributes file](#infile_local_attributes) (relative to `SETNGS_PATH`).

1. `LOCALPARAM_INFO`: File path for the [local parameters file](#infile_local_parameters) (relative to `SETNGS_PATH`).

1. `BASINPARAM_INFO`: File path for the [basin parameters file](#infile_basin_parameters) (relative to `SETNGS_PATH`).

1. `FORCING_FILELIST`: File path for the [list of forcing files file](#infile_forcing_list) (relative to `SETNGS_PATH`).

1. `MODEL_INITCOND`: File path for the [initial conditions file](#infile_initial_conditions) (relative to `SETNGS_PATH`).

1. `PARAMETER_TRIAL`: File path for the [trial parameters file](#infile_trial_parameters) (relative to `SETNGS_PATH`).

1. `OUTPUT_PREFIX`: Text string prepended to each output filename to identify a specific model setup. Note that the user can further modify the output file name at run-time by using the `-s|--suffix` command-line option.


## Model decisions file
<a id="infile_model_decisions"></a>
The model decisions file is an [ASCII file](#infile_format_ASCII) that indicates the model decisions with which SUMMA is configured. The model decisions file is parsed by `build/source/engine/mDecisions.f90`, which also serves as the file of record for all available options for the individual model decisions. The names for the model decisions are found in `build/source/dshare/get_ixname.f90:function get_ixdecisions(varName)`. Detailed information about the individual model decisions and their associated options can be found in the [configuration section](../configuration/SUMMA_model_decisions.md).

Model decisions can be specified in any order with one decision per line. The decisions take the form `<keyword> <value>`, where `<keyword>` is the decision to be made and `<value>` is the option that is selected for that decision. For example, the line `f_Richards mixdform` indicates that the mixed form (liquid/frozen) of the Richards's equation should be used in the simulation(`mixdform` option for the `f_Richards` decision). Another option for this model decision would be `moisture`, which would be the moisture-based form.

The model decisions file must also contain the start (`simulStart`) and end (`simulFinsh`) times of the simulation. These are specified as `'YYYY-MM-DD hh:mm'` and must be enclosed in single quotes. They are typically the first model decisions to be specified.

The model decisions and their options or values are listed in the following tables. Note that the decisions and their options are **case sensitive**. For details about each option see the [configuration section](../configuration/SUMMA_model_decisions.md).

| Decision  | option/value  | notes |
|---|---|---|
|simulStart | 'YYYY-MM-DD hh:mm' | ( 1) simulation start time
|simulFinsh | 'YYYY-MM-DD hh:mm' | ( 2) simulation end time
|soilCatTbl | STAS <br> STAS-RUC <br> ROSETTA | ( 3) soil-category dataset |
|vegeParTbl | USGS <nr> MODIFIED_IGBP_MODIS_NOAH | ( 4) vegetation category dataset
|soilStress | NoahType <br> CLM_Type <br> SiB_Type | ( 5) choice of function for the soil moisture control on stomatal resistance
|stomResist | BallBerry <br> Jarvis <br> simpleResistance <br> BallBerryFlex <br> BallBerryTest | ( 6) choice of function for stomatal resistance
|bbTempFunc | q10Func <br> Arrhenius | ( 7) Ball-Berry: leaf temperature controls on photosynthesis + stomatal resistance
|bbHumdFunc | humidLeafSurface <br> scaledHyperbolic | ( 8) Ball-Berry: humidity controls on stomatal resistance
|bbElecFunc | linear <br> linearJmax <br> quadraticJmax | ( 9) Ball-Berry: dependence of photosynthesis on PAR
|bbCO2point | origBWB <br> Leuning | (10) Ball-Berry: use of CO2 compensation point to calculate stomatal resistance
|bbNumerics | NoahMPsolution <br> newtonRaphson | (11) Ball-Berry: iterative numerical solution method
|bbAssimFnc | colimitation <br> minFunc | (12) Ball-Berry: controls on carbon assimilation
|bbCanIntg8 | constantScaling <br> laiScaling | (13) Ball-Berry: scaling of photosynthesis from the leaf to the canopy
|num_method | itertive <br> non_iter <br> itersurf | (14) choice of numerical method
|fDerivMeth | numericl <br> analytic | (15) choice of method to calculate flux derivatives
|LAI_method | monTable <br> specified | (16) choice of method to determine LAI and SAI
|cIntercept | sparseCanopy <br> storageFunc <br> notPopulatedYet | (17) choice of parameterization for canopy interception
|f_Richards | moisture <br> mixdform | (18) form of Richards' equation
|groundwatr | qTopmodl <br> bigBuckt <br> noXplict | (19) choice of groundwater parameterization
|hc_profile | constant <br> pow_prof | (20) choice of hydraulic conductivity profile
|bcUpprTdyn | presTemp <br> nrg_flux <br> zeroFlux | (21) type of upper boundary condition for thermodynamics
|bcLowrTdyn | presTemp <br> zeroFlux | (22) type of lower boundary condition for thermodynamics
|bcUpprSoiH | presHead <br> liq_flux | (23) type of upper boundary condition for soil hydrology
|bcLowrSoiH | presHead <br> bottmPsi <br> drainage <br> zeroFlux | (24) type of lower boundary condition for soil hydrology
|veg_traits | Raupach_BLM1994 <br> CM_QJRMS1998 <br> vegTypeTable | (25) choice of parameterization for vegetation roughness length and displacement height
|rootProfil | powerLaw <br> doubleExp | (26) choice of parameterization for the rooting profile
|canopyEmis | simplExp <br> difTrans | (27) choice of parameterization for canopy emissivity
|snowIncept | stickySnow <br> lightSnow | (28) choice of parameterization for snow interception
|windPrfile | exponential <br> logBelowCanopy | (29) choice of canopy wind profile
|astability | standard <br> louisinv <br> mahrtexp | (30) choice of stability function
|compaction | consettl <br> anderson | (31) choice of compaction routine
|snowLayers | jrdn1991 <br> CLM_2010 | (32) choice of method to combine and sub-divide snow layers
|thCondSnow | tyen1965 <br> melr1977 <br> jrdn1991 <br> smnv2000 | (33) choice of thermal conductivity representation for snow
|thCondSoil | funcSoilWet <br> mixConstit <br> hanssonVZJ | (34) choice of thermal conductivity representation for soil
|canopySrad | noah_mp <br> CLM_2stream <br> UEB_2stream <br> NL_scatter <br> BeersLaw | (35) choice of method for canopy shortwave radiation
|alb_method | conDecay <br> varDecay | (36) choice of albedo representation
|spatial_gw | localColumn <br> singleBasin | (37) choice of method for spatial representation of groundwater
|subRouting | timeDlay <br> qInstant | (38) choice of method for sub-grid routing
|snowDenNew | hedAndPom <br> anderson <br> pahaut_76 <br> constDens | (39) choice of method for new snow density

The model decisions for each simulation are included as global attributes in [SUMMA output files](SUMMA_output.md).

## Output control file
<a id="infile_output_control"></a>
The output control file is an [ASCII file](#infile_format_ASCII) that specifies which variables are retained in the [SUMMA output files](SUMMA_output.md). The output control file is parsed by `build/source/dshare/popMetadat.f90:read_output_file()`

SUMMA is pretty flexible in its output. There are many variables that you can output and for most of them you can also choose to record summary statistics. For example, you can configure the model to run with meteorological forcings that are defined every hour, but only save summary output with a daily time step. This flexibility comes at the small price that you need to be clear in specifying what output you want.

The output control file includes a listing of model variables that you would like to store, with one model variable per line. The variables that are available for output are the individual entries in the data structures specified in `build/source/dshare/var_lookup.f90`. Because there are many, there is not much point in repeating them here, but we direct the user to the model code. Any of the variables specified in the following structures can be specified in the output control file: `iLook_time`, `iLook_force`, `iLook_attr`, `iLook_type`, `iLook_param`, `iLook_index`, `iLook_prog`, `iLook_diag`, `iLook_flux`, `iLook_bpar`, `iLook_bvar`, `iLook_deriv`. SUMMA will print an error message if a specific variable cannot be output, so the faster way may be to select any variable in `build/source/dshare/var_lookup.f90` and remove it if it is not available for output. Some of these variables are time constant parameters, so it may not make sense to output them as part of your time-varying output. In addition, some of the variables may only be useful for debugging use, but that is up to the user.

At a minimum, each line in the output control file will contain two fields, separated by a `|`. The first field will be the variable name as specified in `build/source/dshare/var_lookup.f90` (case-sensitive). The second field will be the frequency of the model output specified as a multiple of the time resolution in the model forcing files. Thus, if you want to output data for every forcing time step, then this value should be equal to 1. If you want daily output and your forcing frequency is 3 hours, then this value should be equal to 8. Note that you can specify different output frequencies for separate variables, but at this time you can specify only a single output frequency for each variable. For example, you can store `scalarSenHeatTotal` with an output frequency of 1 and `scalarLatHeatTotal` with an output frequency of 8, but you cannot specify two different output frequencies for `scalarSenHeatTotal`.

For most variables you can also output a statistical summary if you output variables at a lower frequency than your forcing frequency. To do this, you extend the number of fields you specify in the output control file, with all fields separated by a `|`. For the fields after the first two, you specify a series of 0's and 1's, which indicate that a specific statistic should not (0) or should be stored (1). The available statistics are (in order) the instantaneous value, the sum over the interval, the mean, the variance, the minimum, the maximum and the mode. So, a complete line in the output control file would be
```
! varName          | outFreq | inst | sum | mean | var | min | max | mode
scalarSenHeatTotal | 24      | 0    | 1   | 1    | 0   | 1   | 1   | 0
```
In this example, the first line is a comment (starts with `!`) and then the sum, mean, min, max are calculated for `scalarSenHeatTotal` across 24 forcing time steps and written to the output file.

## List of forcing files file
<a id="infile_forcing_list"></a>

## Initial conditions file
<a id="infile_initial_conditions"></a>

## Attribute and parameter files
SUMMA uses a number of files to specify model attributes and parameters. Although SUMMA's distinction between attributes and parameters is somewhat arbitrary, attributes generally describe chracteristics of the model domain that are time-invariant during the simulation, such as GRU and HRU identifiers, spatial organization, an topography. The important part for understanding the organization of the SUMMA input files is that the values specified in the [local attributes file](#infile_local_attributes) do not overlap with those in the various parameter files. Thus, these values do not overwrite any attributes specified elsewhere. In contrast, the various parameter file are read in sequence (as explained in the next paragraph) and parameter values that are read in from the input files successively overwrite values that have been specified earlier.

The figure below shows the order in which SUMMA processes the various attribute and parameter files. First, the [local attributes file](#infile_local_attributes) is processed, which provides information about the organization of the GRUs and HRUs as well as some other information. Then, SUMMA parses the [local parameters file](#infile_local_parameters), which provides spatially constant values for all SUMMA parameters that need to be specified at the HRU level. SUMMA then parses the [basin parameters file](#infile_basin_parameters), which provides spatially constant values for all SUMMA parameters that need to be specified at the GRU level. In this case, it does not really matter which files is parsed first. The information in these two files does not overlap. At this point in SUMMA's initialization, all GRU and HRU parameters have been initialized to spatially constant values. SUMMA has inherited some routines from the NOAH land surface model and the next step is to parse the [NOAH parameter tables](#infile_noah_tables). The information in these tables is used to overwrite the spatially constant values that have already been initialized for each HRU. Finally, the [trial parameters file](#infile_trial_parameters) is parsed to provide additional GRU and HRU specific information. The values from this file will overwrite existing values. The number of variables specified in the [trial parameters file](#infile_trial_parameters) will vary with the amount of location-specific information that you have available for your simulation.

![Order in which SUMMA model attributes and parameters are specified and processed](../assets/img/SUMMA_parameters_spec_order.png)<a id="SUMMA_parameters_spec_order"></a>
*Order in which SUMMA model attributes and parameters are specified and processed.*

### Local attributes file
<a id="infile_local_attributes"></a>
The local attributes file is a [NetCDF file](#infile_format_nc) that specifies model element attributes for GRUs and individual HRUs. The local attributes file is parsed by `build/source/driver/multi_driver.f90` and `build/source/engine/read_attrb.f90`. As described above, the attributes specified in this file are separate from the values specified in the various parameter files.

The local attributes file contains a `gru` and an `hru` dimension as specified in the table below. All variables in the local attributes file must be specified.

| Variable | dimension | type | units | long name | notes |
|----------|-----------|------|-------|-----------|-------|
hruId | hru | int | - | Index of hydrological response unit (HRU) | Unique numeric ID for each HRU |
gruId | gru | int | - | Index of grouped response unit (GRU) | Unique numeric ID for each GRU |
hru2gruId | hru | int | - | Index of GRU to which the HRU belongs | gruId of the GRU to which the HRU belongs |
downHRUindex | hru | int | - | Index of downslope HRU (0 = basin outlet) | Downslope HRU must be within the same GRU. If the value is 0, then there is no exchange to a neighboring HRU. Setting this value to 0 for all HRUs emulates a series of independent columns |
longitude | hru | double | Decimal degree east | Longitude of HRU's centroid | West is negative or greater than 180 |
latitude | hru | double | Decimal degree north | Latitude of HRU's centroid | South is negative |
elevation | hru | double | m | Elevation of HRU's centroid | |
HRUarea | hru | double | m^2 | Area of HRU | |
tan_slope | hru | double | m m-1 | Average tangent slope of HRU | |
contourLength | hru | double | m | Contour length of HRU | Width of a hillslope (m) parallel to a stream. Used in `groundwatr.f90`|
slopeTypeIndex | hru | int | - |	Index defining slope | |
soilTypeIndex | hru | int | - |	Index defining soil type | |
vegTypeIndex | hru | int | - |	Index defining vegetation type | |
mHeight | hru | double | m | Measurement height above bare ground | |

Below is a sample layout of the local attributes file (the output of running `ncdump -h`). In this case,  both the gru and hru dimension are of size 1 (the example is taken from one of the [test cases](../installation/SUMMA_test_cases.md), most of which are point model simulations), but of course there can be many GRUs and HRUs.

```
netcdf sample_local_attributes_file_layout {
dimensions:
	hru = 1 ;
	gru = 1 ;
variables:
	int hruId(hru) ;
		hruId:long_name = "Index of hydrological response unit (HRU)" ;
		hruId:units = "-" ;
		hruId:v_type = "scalarv" ;
	int gruId(gru) ;
		gruId:long_name = "Index of grouped response unit (GRU)" ;
		gruId:units = "-" ;
		gruId:v_type = "scalarv" ;
	int hru2gruId(hru) ;
		hru2gruId:long_name = "Index of GRU to which the HRU belongs" ;
		hru2gruId:units = "-" ;
	int downHRUindex(hru) ;
		downHRUindex:long_name = "Index of downslope HRU (0 = basin outlet)" ;
		downHRUindex:units = "-" ;
	double longitude(hru) ;
		longitude:_FillValue = NaN ;
		longitude:long_name = "Longitude of HRU\'s centroid" ;
		longitude:units = "Decimal degree east" ;
	double latitude(hru) ;
		latitude:_FillValue = NaN ;
		latitude:long_name = "Latitude of HRU\'s centroid" ;
		latitude:units = "Decimal degree north" ;
	double elevation(hru) ;
		elevation:_FillValue = NaN ;
		elevation:long_name = "Elevation of HRU\'s centroid" ;
		elevation:units = "m" ;
	double HRUarea(hru) ;
		HRUarea:_FillValue = NaN ;
		HRUarea:long_name = "Area of HRU" ;
		HRUarea:units = "m^2" ;
	double tan_slope(hru) ;
		tan_slope:_FillValue = NaN ;
		tan_slope:long_name = "Average tangent slope of HRU" ;
		tan_slope:units = "m m-1" ;
	double contourLength(hru) ;
		contourLength:_FillValue = NaN ;
		contourLength:long_name = "Contour length of HRU" ;
		contourLength:units = "m" ;
	int slopeTypeIndex(hru) ;
		slopeTypeIndex:long_name = "Index defining slope" ;
		slopeTypeIndex:units = "-" ;
	int soilTypeIndex(hru) ;
		soilTypeIndex:long_name = "Index defining soil type" ;
		soilTypeIndex:units = "-" ;
	int vegTypeIndex(hru) ;
		vegTypeIndex:long_name = "Index defining vegetation type" ;
		vegTypeIndex:units = "-" ;
	double mHeight(hru) ;
		mHeight:_FillValue = NaN ;
		mHeight:long_name = "Measurement height above bare ground" ;
		mHeight:units = "m" ;
}
```

### Local parameters file
<a id="infile_local_parameters"></a>
The local parameters file is an [ASCII file](#infile_format_ASCII) that specifies spatially constant parameter values for which variables are retained in the [SUMMA output files](SUMMA_output.md). The output control file is parsed by `build/source/dshare/popMetadat.f90:read_output_file()`

### Basin parameters file
<a id="infile_basin_parameters"></a>

### NOAH tables
<a id="infile_noah_tables"></a>

### Trial parameters file
<a id="infile_trial_parameters"></a>

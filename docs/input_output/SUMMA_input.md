# SUMMA Input Files

SUMMA has a large number of input files that configure the model and provide the necessary initial conditions and time-varying boundary conditions to make a model simulation. This can at times be confusing. We encourage the user to look at the [SUMMA test cases](../installation/SUMMA_test_cases.md), which provide working SUMMA setups.

## Master configuration file
<a id="infile_master_configuration"></a>

The master configuration file is provided to SUMMA at run-time as a command-line option. The path to this file needs to be supplied with the `-m` or `--master` command-line flag. The contents of this file orchestrate the remainder of the SUMMA run and are processed by the code in `summa/build/source/hookup/summaFileManager.f90`. The file contents mostly consist of file paths that provide the actual information about the model configuration.

Comments can be added to the file by starting them with `!`. Anything after the `!` will be discarded till the end of the line. You can include as many comments as you want, as they will be stripped (in memory) as SUMMA processes the file.

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

## Output control file
<a id="infile_output_control"></a>

## Local attributes file
<a id="infile_local_attributes"></a>

## Local parameters file
<a id="infile_local_parameters"></a>

## Basin parameters file
<a id="infile_basin_parameters"></a>

## List of forcing files file
<a id="infile_forcing_list"></a>

## Initial conditions file
<a id="infile_initial_conditions"></a>

## Trial parameters file
<a id="infile_trial_parameters"></a>

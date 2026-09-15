# SUMMA base settings
Recall that a SUMMA Grouped Response Unit (GRU) contains at least one but possibly more Hydrologic Response Units (HRUs).

#### basinParamInfo.txt
GRU-level parameters. These control the properties of the aquifer that the HRUs share in certain model setups and the settings for within-GRU routing. The shared aquifer is disabled if all HRUs are modeled as independent soil columns; see setting "settings_summa_connect_HRUs" in the control file. Within-GRU routing is controlled by the model decision `subRouting`. In the default setup in this workflow, this decision is set to `timeDlay` which means that within-GRU routing is active. Note that this means that the mizuRoute setting `doesBasinRoute` should be set to `0` to avoid the water inside a given GRU being routed twice. See:
- https://summa.readthedocs.io/en/latest/input_output/SUMMA_input/#attribute-and-parameter-files

#### localParamInfo.txt
HRU-level parameters. First column is the default value. Second and third columns provide plausible ranges that are currently not used by SUMMA (but must be provided regardless). These columns may be used in the future for built-in parameter sampling or sensitivity analysis. See:
- https://summa.readthedocs.io/en/latest/input_output/SUMMA_input/#attribute-and-parameter-files

#### modelDecisions.txt
Controls the modeling options active within SUMMA. Note that with the current workflow, modeling option `qTopmodl` for decision `groundwatr` should not be used, because the required geospatial parameters are not defined. To use the `qTopmodl` option parameters `tan_slope` and `contourLength` must be set to appropriate values for each HRU in the attributes `.nc` file. As mentioned above, option `timeDlay` for decision `subRouting` should not be used at the same time as mizuRoute's `doesBasinRoute` option. See:
- https://summa.readthedocs.io/en/latest/input_output/SUMMA_input/#model-decisions-file

#### outputControl.txt
Controls which variables SUMMA writes to the output files. See:
- https://summa.readthedocs.io/en/latest/input_output/SUMMA_input/#output-control-file
- https://summa.readthedocs.io/en/latest/input_output/SUMMA_output/

#### `*.TBL`
Contain the lookup tables used for soil parameters (`TBL_SOIL_PARM.TBL`) and vegetation parameters (`TBL_VEGPARM.TBL`). Files `TBL_GENPARM.TBL` and `TBL_MPTABLE.TBL` are legacy files that are no longer used but still need to be provided.
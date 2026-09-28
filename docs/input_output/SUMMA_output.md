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

<a id="outfile_variables"></a>
## Output variables
The tables below list the time-varying variables that can be requested in the [output control file](SUMMA_input#infile_output_control): the meteorological forcing (`forc`), prognostic/state (`prog`), diagnostic (`diag`), flux (`flux`) and basin-average (`bvar`) structures in `build/source/dshare/popMetadat.f90`. Time-constant parameters are documented separately on the [SUMMA parameters](SUMMA_parameters.md) page. The `attr`, `type`, `grid`, `time`, and `indx` structures are also available for output (see the [output control file](SUMMA_input#infile_output_control) section) but are bookkeeping/identification fields rather than model variables, so they are not enumerated here.

<a id="outvar_forc"></a>
### Forcing variables

| Variable | Description | Units |
|---|---|---|
| `time` | time since time reference | seconds since 1990-1-1 0:0:0.0 -0:00 |
| `pptrate` | precipitation rate | kg m-2 s-1 |
| `SWRadAtm` | downward shortwave radiation at the upper boundary | W m-2 |
| `LWRadAtm` | downward longwave radiation at the upper boundary | W m-2 |
| `airtemp` | air temperature at the measurement height | K |
| `windspd` | wind speed at the measurement height | m s-1 |
| `airpres` | air pressure at the the measurement height | Pa |
| `spechum` | specific humidity at the measurement height | g g-1 |

<a id="outvar_prog"></a>
### Prognostic (state) variables

| Variable | Description | Units |
|---|---|---|
| `dt_init` | length of initial time step at start of next data interval | s |
| `scalarCanopyIce` | mass of ice on the vegetation canopy | kg m-2 |
| `scalarCanopyLiq` | mass of liquid water on the vegetation canopy | kg m-2 |
| `scalarCanopyWat` | mass of total water on the vegetation canopy | kg m-2 |
| `scalarCanairTemp` | temperature of the canopy air space | K |
| `scalarCanopyTemp` | temperature of the vegetation canopy | K |
| `spectralSnowAlbedoDiffuse` | diffuse snow albedo for individual spectral bands | - |
| `scalarSnowAlbedo` | snow albedo for the entire spectral band | - |
| `scalarSnowDepth` | total snow depth | m |
| `scalarSWE` | snow water equivalent | kg m-2 |
| `scalarSfcMeltPond` | ponded water caused by melt of the "snow without a layer" | kg m-2 |
| `glacMass4AreaChange` | since updateJulDay glacier layers together mass change | kg m-2 |
| `scalarGlceWE` | glacier ice (not snow) water equivalent change over simulation | kg m-2 |
| `mLayerTemp` | temperature of each layer | K |
| `mLayerVolFracIce` | volumetric fraction of ice in each layer | - |
| `mLayerVolFracLiq` | volumetric fraction of liquid water in each layer | - |
| `mLayerVolFracWat` | volumetric fraction of total water in each layer | - |
| `mLayerMatricHead` | matric head of water in the soil | m |
| `scalarCanairEnthalpy` | enthalpy of the canopy air space | J m-3 |
| `scalarCanopyEnthalpy` | enthalpy of the vegetation canopy | J m-3 |
| `mLayerEnthalpy` | enthalpy of the layers | J m-3 |
| `scalarAquiferStorage` | water required to bring aquifer to the bottom of the soil profile | m |
| `scalarSurfaceTemp` | surface temperature (just a copy of the upper-layer temperature) | K |
| `mLayerDepth` | depth of each layer | m |
| `mLayerHeight` | height of the layer mid-point (top of soil = 0) | m |
| `iLayerHeight` | height of the layer interface (top of soil = 0) | m |
| `DOMarea` | area of the domain | m2 |
| `DOMelev` | elevation of the domain | m |
| `DOMtan_slope` | tan local ground surface slope of the domain | - |
| `DOMaspect` | azimuth in degrees East of North of the domain | degrees |
| `DOMcontourLength` | length of contour at downslope edge of the domain | m |
| `scalarAblFrac` | fraction of the domain that is in a glacier ablation zone | - |

<a id="outvar_diag"></a>
### Diagnostic variables

| Variable | Description | Units |
|---|---|---|
| `scalarCanopyDepth` | canopy depth | m |
| `scalarBulkVolHeatCapVeg` | bulk volumetric heat capacity of vegetation | J m-3 K-1 |
| `scalarCanopyCm` | Cm for canopy vegetation | J kg-1 |
| `scalarCanopyEmissivity` | effective canopy emissivity | - |
| `scalarRootZoneTemp` | average temperature of the root zone | K |
| `scalarLAI` | one-sided leaf area index | m2 m-2 |
| `scalarSAI` | one-sided stem area index | m2 m-2 |
| `scalarExposedLAI` | exposed leaf area index (after burial by snow) | m2 m-2 |
| `scalarExposedSAI` | exposed stem area index (after burial by snow) | m2 m-2 |
| `scalarAdjMeasHeight` | adjusted measurement height for cases snowDepth>mHeight | m |
| `scalarCanopyIceMax` | maximum interception storage capacity for ice | kg m-2 |
| `scalarCanopyLiqMax` | maximum interception storage capacity for liquid water | kg m-2 |
| `scalarGrowingSeasonIndex` | growing season index (0=off, 1=on) | - |
| `mLayerVolHtCapBulk` | volumetric heat capacity in each layer | J m-3 K-1 |
| `mLayerCm` | Cm for each layer | J m-3 |
| `mLayerThermalC` | thermal conductivity at the mid-point of each layer | W m-1 K-1 |
| `iLayerThermalC` | thermal conductivity at the interface of each layer | W m-1 K-1 |
| `scalarCanopyEnthTemp` | temperature component of enthalpy of the vegetation canopy | J m-3 |
| `mLayerEnthTemp` | temperature component of enthalpy of the layers | J m-3 |
| `scalarTotalSnowEnthalpy` | total enthalpy of the snow column | J m-3 |
| `scalarTotalLakeEnthalpy` | total enthalpy of the lake column | J m-3 |
| `scalarTotalSoilEnthalpy` | total enthalpy of the soil column | J m-3 |
| `scalarTotalGlceEnthalpy` | total enthalpy of the glacier ice column | J m-3 |
| `scalarVPair` | vapor pressure of the air above the vegetation canopy | Pa |
| `scalarVP_CanopyAir` | vapor pressure of the canopy air space | Pa |
| `scalarTwetbulb` | wet bulb temperature | K |
| `scalarSnowfallTemp` | temperature of fresh snow | K |
| `scalarNewSnowDensity` | density of fresh snow | kg m-3 |
| `scalarO2air` | atmospheric o2 concentration | Pa |
| `scalarCO2air` | atmospheric co2 concentration | Pa |
| `windspd_x` | wind speed at 10 meter height in x-direction | m s-1 |
| `windspd_y` | wind speed at 10 meter height in y-direction | m s-1 |
| `scalarCosZenith` | cosine of the solar zenith angle | - |
| `scalarFractionDirect` | fraction of direct radiation (0-1) | - |
| `scalarCanopySunlitFraction` | sunlit fraction of canopy | - |
| `scalarCanopySunlitLAI` | sunlit leaf area | - |
| `scalarCanopyShadedLAI` | shaded leaf area | - |
| `spectralAlbGndDirect` | direct  albedo of underlying surface for each spectral band | - |
| `spectralAlbGndDiffuse` | diffuse albedo of underlying surface for each spectral band | - |
| `scalarGroundAlbedo` | albedo of the ground surface | - |
| `scalarLatHeatSubVapCanopy` | latent heat of sublimation/vaporization used for veg canopy | J kg-1 |
| `scalarLatHeatSubVapGround` | latent heat of sublimation/vaporization used for ground surface | J kg-1 |
| `scalarSatVP_CanopyTemp` | saturation vapor pressure at the temperature of vegetation canopy | Pa |
| `scalarSatVP_GroundTemp` | saturation vapor pressure at the temperature of the ground | Pa |
| `scalarZ0Canopy` | roughness length of the canopy | m |
| `scalarWindReductionFactor` | canopy wind reduction factor | - |
| `scalarZeroPlaneDisplacement` | zero plane displacement | m |
| `scalarRiBulkCanopy` | bulk Richardson number for the canopy | - |
| `scalarRiBulkGround` | bulk Richardson number for the ground surface | - |
| `scalarCanopyStabilityCorrection` | stability correction for the canopy | - |
| `scalarGroundStabilityCorrection` | stability correction for the ground surface | - |
| `scalarIntercellularCO2Sunlit` | carbon dioxide partial pressure of leaf interior (sunlit leaves) | Pa |
| `scalarIntercellularCO2Shaded` | carbon dioxide partial pressure of leaf interior (shaded leaves) | Pa |
| `scalarTranspireLim` | aggregate soil moisture and aquifer control on transpiration | - |
| `scalarTranspireLimAqfr` | aquifer storage control on transpiration | - |
| `scalarFoliageNitrogenFactor` | foliage nitrogen concentration (1=saturated) | - |
| `scalarSoilRelHumidity` | relative humidity in the soil pores in the upper-most soil layer | - |
| `mLayerTranspireLim` | soil moist & veg limit on transpiration for each layer | - |
| `mLayerRootDensity` | fraction of roots in each soil layer | - |
| `scalarAquiferRootFrac` | fraction of roots below the soil profile (in the aquifer) | - |
| `scalarFracLiqVeg` | fraction of liquid water on vegetation | - |
| `scalarCanopyWetFraction` | fraction canopy that is wet | - |
| `scalarSnowAge` | non-dimensional snow age | - |
| `scalarGroundSnowFraction` | fraction ground that is covered with snow | - |
| `spectralSnowAlbedoDirect` | direct snow albedo for individual spectral bands | - |
| `mLayerFracLiq` | fraction of liquid water in each snow, lake, or glce layer | - |
| `mLayerThetaResid` | residual volumetric water content in each snow, lake, or glce layer | - |
| `mLayerPoreSpace` | total pore space in each snow, lake, or glce layer | - |
| `mLayerMeltFreeze` | ice content change from melt/freeze in each layer | kg m-3 |
| `spectralFrznWatAlbedo` | albedo of frozen water in each spectral band | - |
| `spectralOpenWatAlbedo` | albedo of open water in each spectral band | - |
| `scalarTotalMassChange` | mass change of all system together (kg m-2 s-1) | kg m-2 s-1 |
| `scalarInfilArea` | fraction of unfrozen area where water can infiltrate | - |
| `scalarSaturatedArea` | fraction of area that is considered saturated | - |
| `scalarFrozenArea` | fraction of area that is considered impermeable due to soil ice | - |
| `scalarSoilControl` | soil control on infiltration for derivative | - |
| `scalarSoilControlBot` | soil control on bottom capillary fluxes for derivative | - |
| `mLayerVolFracAir` | volumetric fraction of air in each layer | - |
| `mLayerTcrit` | critical soil temperature above which all water is unfrozen | K |
| `mLayerCompress` | change in volumetric water content due to compression of soil | s-1 |
| `scalarSoilCompress` | change in total soil storage due to compression of soil matrix | kg m-2 s-1  |
| `mLayerMatricHeadLiq` | matric potential of liquid water | m |
| `scalarTotalSoilLiq` | total mass of liquid water in the soil | kg m-2 |
| `scalarTotalSoilIce` | total mass of ice in the soil | kg m-2 |
| `scalarTotalSoilWat` | total mass of water in the soil | kg m-2 |
| `scalarVGn_m` | van Genuchten "m" parameter | - |
| `numFluxCalls` | number of flux calls | - |
| `wallClockTime` | wall clock time for physics routines | s |
| `meanStepSize` | mean time step size over data window | s |
| `balanceCasNrg` | balance of energy in the canopy air space on data window | W m-3 |
| `balanceVegNrg` | balance of energy in the vegetation on data window | W m-3 |
| `balanceLayerNrg` | balance of energy in each layer on substep | W m-3 |
| `balanceSnowNrg` | balance of energy in the snow on data window | W m-3 |
| `balanceLakeNrg` | balance of energy in the lake on data window | W m-3 |
| `balanceSoilNrg` | balance of energy in the soil on data window | W m-3 |
| `balanceGlceNrg` | balance of energy in the glacier ice on data window | W m-3 |
| `balanceVegMass` | balance of water in the vegetation on data window | kg m-3 s-1 |
| `balanceLayerMass` | balance of water in each layer on substep | kg m-3 s-1 |
| `balanceSnowMass` | balance of water in the snow on data window | kg m-3 s-1 |
| `balanceLakeMass` | balance of water in the lake on data window | kg m-3 s-1 |
| `balanceSoilMass` | balance of water in the soil on data window | kg m-3 s-1 |
| `balanceGlceMass` | balance of water in the glacier ice on data window | kg m-3 s-1 |
| `balanceAqMass` | balance of water in the aquifer on data window | kg m-2 s-1 |
| `numSteps` | number of steps taken by the integrator | - |
| `numResEvals` | number of residual evaluations | - |
| `numLinSolvSetups` | number of linear solver setups | - |
| `numErrTestFails` | number of error test failures | - |
| `kLast` | method order used on the last internal step | - |
| `kCur` | method order to be used on the next internal step | - |
| `hInitUsed` | step size used on the first internal step | s |
| `hLast` | step size used on the last internal step | s |
| `hCur` | step size to be used on the next internal step | s |
| `tCur` | current time reached by the integrator | s |

<a id="outvar_flux"></a>
### Flux variables

| Variable | Description | Units |
|---|---|---|
| `scalarCanairNetNrgFlux` | net energy flux for the canopy air space | W m-2 |
| `scalarCanopyNetNrgFlux` | net energy flux for the vegetation canopy | W m-2 |
| `scalarGroundNetNrgFlux` | net energy flux for the ground surface | W m-2 |
| `scalarCanopyNetLiqFlux` | net liquid water flux for the vegetation canopy | kg m-2 s-1 |
| `scalarRainfall` | computed rainfall rate | kg m-2 s-1 |
| `scalarSnowfall` | computed snowfall rate | kg m-2 s-1 |
| `spectralIncomingDirect` | incoming direct solar radiation in each wave band | W m-2 |
| `spectralIncomingDiffuse` | incoming diffuse solar radiation in each wave band | W m-2 |
| `scalarCanopySunlitPAR` | average absorbed par for sunlit leaves | W m-2 |
| `scalarCanopyShadedPAR` | average absorbed par for shaded leaves | W m-2 |
| `spectralBelowCanopyDirect` | downward direct flux below veg layer for each spectral band | W m-2 |
| `spectralBelowCanopyDiffuse` | downward diffuse flux below veg layer for each spectral band | W m-2 |
| `scalarBelowCanopySolar` | solar radiation transmitted below the canopy | W m-2 |
| `scalarCanopyAbsorbedSolar` | solar radiation absorbed by canopy | W m-2 |
| `scalarGroundAbsorbedSolar` | solar radiation absorbed by ground | W m-2 |
| `scalarLWRadCanopy` | longwave radiation emitted from the canopy | W m-2 |
| `scalarLWRadGround` | longwave radiation emitted at the ground surface | W m-2 |
| `scalarLWRadUbound2Canopy` | downward atmospheric longwave radiation absorbed by the canopy | W m-2 |
| `scalarLWRadUbound2Ground` | downward atmospheric longwave radiation absorbed by the ground | W m-2 |
| `scalarLWRadUbound2Ubound` | atmospheric radiation refl by ground + lost thru upper boundary | W m-2 |
| `scalarLWRadCanopy2Ubound` | longwave radiation emitted from canopy lost thru upper boundary | W m-2 |
| `scalarLWRadCanopy2Ground` | longwave radiation emitted from canopy absorbed by the ground | W m-2 |
| `scalarLWRadCanopy2Canopy` | canopy longwave reflected from ground and absorbed by the canopy | W m-2 |
| `scalarLWRadGround2Ubound` | longwave radiation emitted from ground lost thru upper boundary | W m-2 |
| `scalarLWRadGround2Canopy` | longwave radiation emitted from ground and absorbed by the canopy | W m-2 |
| `scalarLWNetCanopy` | net longwave radiation at the canopy | W m-2 |
| `scalarLWNetGround` | net longwave radiation at the ground surface | W m-2 |
| `scalarLWNetUbound` | net longwave radiation at the upper atmospheric boundary | W m-2 |
| `scalarEddyDiffusCanopyTop` | eddy diffusivity for heat at the top of the canopy | m2 s-1 |
| `scalarFrictionVelocity` | friction velocity (canopy momentum sink) | m s-1 |
| `scalarWindspdCanopyTop` | windspeed at the top of the canopy | m s-1 |
| `scalarWindspdCanopyBottom` | windspeed at the height of the bottom of the canopy | m s-1 |
| `scalarGroundResistance` | below canopy aerodynamic resistance | s m-1 |
| `scalarCanopyResistance` | above canopy aerodynamic resistance | s m-1 |
| `scalarLeafResistance` | mean leaf boundary layer resistance per unit leaf area | s m-1 |
| `scalarSoilResistance` | soil surface resistance | s m-1 |
| `scalarSenHeatTotal` | sensible heat from the canopy air space to the atmosphere | W m-2 |
| `scalarSenHeatCanopy` | sensible heat from the canopy to the canopy air space | W m-2 |
| `scalarSenHeatGround` | sensible heat from the ground (below canopy or non-vegetated) | W m-2 |
| `scalarLatHeatTotal` | latent heat from the canopy air space to the atmosphere | W m-2 |
| `scalarLatHeatCanopyEvap` | evaporation latent heat from the canopy to the canopy air space | W m-2 |
| `scalarLatHeatCanopyTrans` | transpiration latent heat from the canopy to the canopy air space | W m-2 |
| `scalarLatHeatGround` | latent heat from the ground (below canopy or non-vegetated) | W m-2 |
| `scalarCanopyAdvectiveHeatFlux` | heat advected to the canopy with precipitation (snow + rain) | W m-2 |
| `scalarGroundAdvectiveHeatFlux` | heat advected to the ground with throughfall + unloading/drainage | W m-2 |
| `scalarCanopySublimation` | canopy sublimation/frost | kg m-2 s-1 |
| `scalarGroundSublimation` | ground sublimation/frost (below canopy or non-vegetated snow or ice) | kg m-2 s-1 |
| `scalarStomResistSunlit` | stomatal resistance for sunlit leaves | s m-1 |
| `scalarStomResistShaded` | stomatal resistance for shaded leaves | s m-1 |
| `scalarPhotosynthesisSunlit` | sunlit photosynthesis | umolco2 m-2 s-1 |
| `scalarPhotosynthesisShaded` | shaded photosynthesis | umolco2 m-2 s-1 |
| `scalarCanopyTranspiration` | canopy transpiration | kg m-2 s-1 |
| `scalarCanopyEvaporation` | canopy evaporation/condensation | kg m-2 s-1 |
| `scalarGroundEvaporation` | ground evaporation/condensation (below canopy or non-vegetated) | kg m-2 s-1 |
| `mLayerTranspire` | transpiration loss from each soil layer | m s-1 |
| `scalarThroughfallSnow` | snow that reaches the ground without ever touching the canopy | kg m-2 s-1 |
| `scalarThroughfallRain` | rain that reaches the ground without ever touching the canopy | kg m-2 s-1 |
| `scalarCanopySnowUnloading` | unloading of snow from the vegetation canopy | kg m-2 s-1 |
| `scalarCanopyLiqDrainage` | drainage of liquid water from the vegetation canopy | kg m-2 s-1 |
| `iLayerConductiveFlux` | conductive energy flux at layer interfaces | W m-2 |
| `iLayerAdvectiveFlux` | advective energy flux at layer interfaces | W m-2 |
| `iLayerNrgFlux` | energy flux at layer interfaces | W m-2 |
| `mLayerNrgFlux` | net energy flux for each layer within the layer domains | J m-3 s-1 |
| `scalarSnowDrainage` | drainage from the bottom of the snow profile | m s-1 |
| `scalarLakeDrainage` | drainage from the bottom of the lake | m s-1 |
| `scalarGlceMelt` | glacier ice melt | m s-1 |
| `iLayerLiqFluxSnLaGl` | liquid flux at snow lake glce layer interfaces | m s-1 |
| `scalarSurfaceIceMelt` | surface ice melt flux | m s-1 |
| `mLayerLiqFluxSnLaGl` | net liquid water flux for each snow lake glce layer | s-1 |
| `scalarRainPlusMelt` | rain plus melt, used as input to soil before surface runoff | m s-1 |
| `scalarMaxInfilRate` | maximum infiltration rate | m s-1 |
| `scalarInfiltration` | infiltration of water into the soil profile | m s-1 |
| `scalarExfiltration` | exfiltration of water from the top of the soil profile | m s-1 |
| `scalarSurfaceRunoff` | surface runoff | m s-1 |
| `scalarSurfaceRunoff_IE` | infiltration excess surface runoff | m s-1 |
| `scalarSurfaceRunoff_SE` | saturation excess surface runoff | m s-1 |
| `mLayerSatHydCondMP` | saturated hydraulic conductivity of macropores in each layer | m s-1 |
| `mLayerSatHydCond` | saturated hydraulic conductivity in each layer | m s-1 |
| `iLayerSatHydCond` | saturated hydraulic conductivity in each layer interface | m s-1 |
| `mLayerHydCond` | hydraulic conductivity in each layer | m s-1 |
| `iLayerLiqFluxSoil` | liquid flux at soil layer interfaces | m s-1 |
| `mLayerLiqFluxSoil` | net liquid water flux for each soil layer | s-1 |
| `mLayerBaseflow` | baseflow from each soil layer | m s-1 |
| `mLayerColumnInflow` | total inflow to each layer in a given soil column | m3 s-1 |
| `mLayerColumnOutflow` | total outflow from each layer in a given soil column | m3 s-1 |
| `scalarSoilBaseflow` | total baseflow from the soil profile | m s-1 |
| `scalarSoilDrainage` | drainage from the bottom of the soil profile | m s-1 |
| `scalarAquiferRecharge` | recharge to the aquifer | m s-1 |
| `scalarAquiferTranspire` | transpiration loss from the aquifer | m s-1 |
| `scalarAquiferBaseflow` | baseflow from the aquifer | m s-1 |
| `scalarTotalET` | total ET | kg m-2 s-1 |
| `scalarTotalRunoff` | total runoff | m s-1 |
| `scalarGlacierMelt` | glacier system melt (goes into glacier internal reservoir) | m s-1 |
| `scalarNetRadiation` | net radiation | W m-2 |

<a id="outvar_bvar"></a>
### Basin-average variables

| Variable | Description | Units |
|---|---|---|
| `basin__TotalArea` | total basin area | m2 |
| `basin__SurfaceRunoff` | surface runoff | m s-1 |
| `basin__ColumnOutflow` | outflow from all "outlet" HRUs (with no downstream HRU) | m3 s-1 |
| `basin__AquiferStorage` | aquifer storage | m |
| `basin__AquiferRecharge` | recharge to the aquifer | m s-1 |
| `basin__AquiferBaseflow` | baseflow from the aquifer | m s-1 |
| `basin__AquiferTranspire` | transpiration loss from the aquifer | m s-1 |
| `basin__TotalRunoff` | total runoff to channel from all active components | m s-1 |
| `basin__SoilDrainage` | soil drainage | m s-1 |
| `basin__GlacierStorage` | glacier storage | Gt |
| `basin__StorageChange` | change in total basin storage | kg m-2 s-1 |
| `basin__GlacierArea` | glacier area | m2 |
| `updateJulDay` | julian day at which glacier geometry was last updated | day |
| `updateJulDayNext` | next julian day at which glacier geometry will be updated | day |
| `routingRunoffFuture` | runoff in future time steps | m s-1 |
| `routingFractionFuture` | fraction of runoff in future time steps | - |
| `averageInstantRunoff` | instantaneous runoff | m s-1 |
| `averageRoutedRunoff` | routed runoff | m s-1 |
| `glacierAblArea` | per glacier ablation area | m2 |
| `glacierAccArea` | per glacier accumulation area | m2 |
| `glacIceRunoffFuture` | per glacier ice reservoir runoff in future time steps | m s-1 |
| `glacSnowRunoffFuture` | per glacier snow reservoir runoff in future time steps | m s-1 |
| `glacFirnRunoffFuture` | per glacier firn reservoir runoff in future time steps | m s-1 |
| `glacierRoutedRunoff` | lapsed glacier runoff | m s-1 |

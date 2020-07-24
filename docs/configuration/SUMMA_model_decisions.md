# Model Decisions in SUMMA

**TODO**: Add information and background (including references) for all the available model decisions.

Information about the selection of specific model decisions is provided to SUMMA via the [model decisions file](../input_output/SUMMA_input.md#infile_model_decisions).

<a id="soilCatTbl"></a>
##  3. soilCatTbl
Soil-category dataset

| Option | Description |
|---|---|
| STAS | **TODO: Describe STAS <br> [Reference](http://doi.org/)** |
| STAS-RUC | **TODO: Describe STAS-RUC <br> [Reference](http://doi.org/)** |
| ROSETTA | **TODO: Describe ROSETTA <br> [Reference](http://doi.org/)** |


<a id="vegeParTbl"></a>
##  4. vegeParTbl
Vegetation category dataset

| Option | Description |
|---|---|
| USGS | **TODO: Describe USGS <br> [Reference](http://doi.org/)** |
| MODIFIED_IGBP_MODIS_NOAH | **TODO: Describe MODIFIED_IGBP_MODIS_NOAH <br> [Reference](http://doi.org/)** |
| plumberCABLE | **TODO: Describe plumberCABLE <br> [Reference](http://doi.org/)** |
| plumberCHTESSEL | **TODO: Describe plumberCHTESSEL <br> [Reference](http://doi.org/)** |
| plumberSUMMA | **TODO: Describe plumberSUMMA <br> [Reference](http://doi.org/)** |

<a id="soilStress"></a>
##  5. soilStress
Function for the soil moisture control on stomatal resistance

| Option | Description |
|---|---|
| NoahType | **Soil stress is a thresholded linear function of volumetric liquid water content <br> [Reference](http://doi.org/)** |
| CLM_Type | **Soil stress is a thresholded linear function of matric head <br> [CLM, 2010]( https://doi.org/10.1029/2011MS00045)** |
| SiB_Type | **Soil stress is an exponential of the log of matric head <br> [Reference](http://doi.org/)** |

<a id="stomResist"></a>
##  6. stomResist
Function for stomatal resistance

| Option | Description |
|---|---|
| BallBerry | **TODO: Describe BallBerry <br> [Reference](http://doi.org/)** |
| Jarvis | **TODO: Describe Jarvis <br> [Reference](http://doi.org/)** |
| simpleResistance | **TODO: Describe simpleResistance <br> [Reference](http://doi.org/)** |
| BallBerryFlex | **TODO: Describe BallBerryFlex <br> [Reference](http://doi.org/)** |
| BallBerryTest | **TODO: Describe BallBerryTest <br> [Reference](http://doi.org/)** |

<a id="bbTempFunc"></a>
##  7. bbTempFunc
Ball-Berry: leaf temperature controls on photosynthesis + stomatal resistance

| Option | Description |
|---|---|
| q10Func | **TODO: Describe q10Func <br> [Reference](http://doi.org/)** |
| Arrhenius | **TODO: Describe Arrhenius <br> [Reference](http://doi.org/)** |


<a id="bbHumdFunc"></a>
##  8. bbHumdFunc
Ball-Berry: humidity controls on stomatal resistance

| Option | Description |
|---|---|
| humidLeafSurface | **TODO: Describe humidLeafSurface <br> [Reference](http://doi.org/)** |
| scaledHyperbolic | **TODO: Describe scaledHyperbolic <br> [Reference](http://doi.org/)** |


<a id="bbElecFunc"></a>
## 9. bbElecFunc
Ball-Berry: dependence of photosynthesis on PAR

| Option | Description |
|---|---|
| linear | **TODO: Describe linear <br> [Reference](http://doi.org/)** |
| linearJmax | **TODO: Describe linearJmax <br> [Reference](http://doi.org/)** |
| quadraticJmax | **TODO: Describe quadraticJmax <br> [Reference](http://doi.org/)** |

<a id="bbCO2point"></a>
## 10. bbCO2point
Ball-Berry: use of CO2 compensation point to calculate stomatal resistance

| Option | Description |
|---|---|
| origBWB | **TODO: Describe origBWB <br> [Reference](http://doi.org/)** |
| Leuning | **TODO: Describe Leuning <br> [Reference](http://doi.org/)** |


<a id="bbNumerics"></a>
## 11. bbNumerics
Ball-Berry: iterative numerical solution method

| Option | Description |
|---|---|
| NoahMPsolution | **TODO: Describe NoahMPsolution <br> [Reference](http://doi.org/)** |
| newtonRaphson | **TODO: Describe newtonRaphson <br> [Reference](http://doi.org/)** |


<a id="bbAssimFnc"></a>
## 12. bbAssimFnc
Ball-Berry: controls on carbon assimilation

| Option | Description |
|---|---|
| colimitation | **TODO: Describe colimitation <br> [Reference](http://doi.org/)** |
| minFunc | **TODO: Describe minFunc <br> [Reference](http://doi.org/)** |


<a id="bbCanIntg8"></a>
## 13. bbCanIntg8
Ball-Berry: scaling of photosynthesis from the leaf to the canopy

| Option | Description |
|---|---|
| constantScaling | **TODO: Describe constantScaling <br> [Reference](http://doi.org/)** |
| laiScaling | **TODO: Describe laiScaling <br> [Reference](http://doi.org/)** |


<a id="num_method"></a>
## 14. num_method
Numerical method

| Option | Description |
|---|---|
| itertive | **TODO: Describe itertive <br> [Reference](http://doi.org/)** |
| non_iter | **TODO: Describe non_iter <br> [Reference](http://doi.org/)** |
| itersurf | **TODO: Describe itersurf <br> [Reference](http://doi.org/)** |

<a id="fDerivMeth"></a>
## 15. fDerivMeth
Method to calculate flux derivatives

| Option | Description |
|---|---|
| numericl | **TODO: Describe numericl <br> [Reference](http://doi.org/)** |
| analytic | **TODO: Describe analytic <br> [Reference](http://doi.org/)** |


<a id="LAI_method"></a>
## 16. LAI_method
Method to determine LAI and SAI

| Option | Description |
|---|---|
| monTable | **TODO: Describe monTable <br> [Reference](http://doi.org/)** |
| specified | **TODO: Describe specified <br> [Reference](http://doi.org/)** |


<a id="cIntercept"></a>
## 17. cIntercept
Parameterization for canopy interception

| Option | Description |
|---|---|
| sparseCanopy | **TODO: Describe sparseCanopy <br> [Reference](http://doi.org/)** |
| storageFunc | **TODO: Describe storageFunc <br> [Reference](http://doi.org/)** |
| notPopulatedYet | **TODO: Describe notPopulatedYet <br> [Reference](http://doi.org/)** |

<a id="f_Richards"></a>
## 18. f_Richards
Form of Richards' equation

| Option | Description |
|---|---|
| moisture | **TODO: Describe moisture <br> [Reference](http://doi.org/)** |
| mixdform | **TODO: Describe mixdform <br> [Reference](http://doi.org/)** |


<a id="groundwatr"></a>
## 19. groundwatr
Groundwater parameterization

| Option | Description |
|---|---|
| qTopmodl | **TODO: Describe qTopmodl <br> [Reference](http://doi.org/)** |
| bigBuckt | **TODO: Describe bigBuckt <br> [Reference](http://doi.org/)** |
| noXplict | **TODO: Describe noXplict <br> [Reference](http://doi.org/)** |

<a id="hc_profile"></a>
## 20. hc_profile
Hydraulic conductivity profile

| Option | Description |
|---|---|
| constant | **TODO: Describe constant <br> [Reference](http://doi.org/)** |
| pow_prof | **TODO: Describe pow_prof <br> [Reference](http://doi.org/)** |


<a id="bcUpprTdyn"></a>
## 21. bcUpprTdyn
Upper boundary condition for thermodynamics

| Option | Description |
|---|---|
| presTemp | **TODO: Describe presTemp <br> [Reference](http://doi.org/)** |
| nrg_flux | **TODO: Describe nrg_flux <br> [Reference](http://doi.org/)** |
| zeroFlux | **TODO: Describe zeroFlux <br> [Reference](http://doi.org/)** |

<a id="bcLowrTdyn"></a>
## 22. bcLowrTdyn
Lower boundary condition for thermodynamics

| Option | Description |
|---|---|
| presTemp | **TODO: Describe presTemp <br> [Reference](http://doi.org/)** |
| zeroFlux | **TODO: Describe zeroFlux <br> [Reference](http://doi.org/)** |


<a id="bcUpprSoiH"></a>
## 23. bcUpprSoiH
Upper boundary condition for soil hydrology

| Option | Description |
|---|---|
| presHead | **TODO: Describe presHead <br> [Reference](http://doi.org/)** |
| liq_flux | **TODO: Describe liq_flux <br> [Reference](http://doi.org/)** |


<a id="bcLowrSoiH"></a>
## 24. bcLowrSoiH
Lower boundary condition for soil hydrology

| Option | Description |
|---|---|
| presHead | **TODO: Describe presHead <br> [Reference](http://doi.org/)** |
| bottmPsi | **TODO: Describe bottmPsi <br> [Reference](http://doi.org/)** |
| drainage | **TODO: Describe drainage <br> [Reference](http://doi.org/)** |
| zeroFlux | **TODO: Describe zeroFlux <br> [Reference](http://doi.org/)** |

<a id="veg_traits"></a>
## 25. veg_traits
Parameterization for vegetation roughness length and displacement height

| Option | Description |
|---|---|
| Raupach_BLM1994 | **TODO: Describe Raupach_BLM1994 <br> [Raupach, 1994](http://doi.org/10.1007/BF00709229)** |
| CM_QJRMS1988 | **TODO: Describe CM_QJRMS1988 <br> [Choudhury and Monteith, 1988](http://doi.org/10.1002/qj.49711448006)** |
| vegTypeTable | **TODO: Describe vegTypeTable <br> [Reference](http://doi.org/)** |

<a id="rootProfil"></a>
## 26. rootProfil
Parameterization for the rooting profile

| Option | Description |
|---|---|
| powerLaw | **TODO: Describe powerLaw <br> [Reference](http://doi.org/)** |
| doubleExp | **TODO: Describe doubleExp <br> [Reference](http://doi.org/)** |


<a id="canopyEmis"></a>
## 27. canopyEmis
Parameterization for canopy emissivity

| Option | Description |
|---|---|
| simplExp | **TODO: Describe simplExp <br> [Reference](http://doi.org/)** |
| difTrans | **TODO: Describe difTrans <br> [Reference](http://doi.org/)** |


<a id="snowIncept"></a>
## 28. snowIncept
Parameterization for snow interception

| Option | Description |
|---|---|
| stickySnow | **Includes a rapid interception increasae between -3 and 0 C from observations of increased cohesion in warm regions.  <br> [Andreadis et al. 2009](https://doi.org/10.1029/2008WR007042)** |
| lightSnow | **Includes a slight decrease in interception after -3 C from obervations in cold regions. <br> [Hedstom and Pomeroy, 1998](https://doi.org/10.1002/(SICI)1099-1085(199808/09)12:10/11<1611::AID-HYP684>3.0.CO;2-4)** |


<a id="windPrfile"></a>
## 29. windPrfile
Canopy wind profile

| Option | Description |
|---|---|
| exponential | **The wind speed profile through the canopy is an exponential decay function. <br> [Reference](http://doi.org/)** |
| logBelowCanopy | **The wind speed profile through the canopy is a logarithmic decay function. <br> [Reference](http://doi.org/)** |


<a id="astability"></a>
## 30. astability
Stability function

| Option | Description |
|---|---|
| standard | **TODO: Describe standard <br> [Reference](http://doi.org/)** |
| louisinv | **TODO: Describe louisinv <br> [Reference](http://doi.org/)** |
| mahrtexp | **TODO: Describe mahrtexp <br> [Reference](http://doi.org/)** |

<a id="compaction"></a>
## 31. compaction
Compaction routine

| Option | Description |
|---|---|
| consettl | **TODO: Describe consettl <br> [Reference](http://doi.org/)** |
| anderson | **TODO: Describe anderson <br> [Reference](http://doi.org/)** |


<a id="snowLayers"></a>
## 32. snowLayers
Method to combine and sub-divide snow layers

| Option | Description |
|---|---|
| jrdn1991 | **Divides the snowpack into a growing layer system where the number of layers is 100.  <br> [Jordan, 1991](https://apps.dtic.mil/docs/citations/ADA245493)** |
| CLM_2010 | **Divides the snowpack into a 5-layer system. The rules of the layer division/merge can be altered to create <5 layer snowpacks. <br> [Community Land Model, 2010]( https://doi.org/10.1029/2011MS00045)** |


<a id="thCondSnow"></a>
## 33. thCondSnow
Thermal conductivity representation for snow

| Option | Description |
|---|---|
| tyen1965 | **TODO: Describe tyen1965 <br> [Reference](http://doi.org/)** |
| melr1977 | **TODO: Describe melr1977 <br> [Reference](http://doi.org/)** |
| jrdn1991 | **TODO: Describe jrdn1991 <br> [Reference](http://doi.org/)** |
| smnv2000 | **TODO: Describe smnv2000 <br> [Reference](http://doi.org/)** |

<a id="thCondSoil"></a>
## 34. thCondSoil
Thermal conductivity representation for soil

| Option | Description |
|---|---|
| funcSoilWet | **TODO: Describe funcSoilWet <br> [Reference](http://doi.org/)** |
| mixConstit | **TODO: Describe mixConstit <br> [Reference](http://doi.org/)** |
| hanssonVZJ | **TODO: Describe hanssonVZJ <br> [Reference](http://doi.org/)** |

<a id="canopySrad"></a>
## 35. canopySrad
Method for canopy shortwave radiation

| Option | Description |
|---|---|
| noah_mp | **TODO: Describe noah_mp <br> [Reference](http://doi.org/)** |
| CLM_2stream | **TODO: Describe CLM_2stream <br> [Reference](http://doi.org/)** |
| UEB_2stream | **TODO: Describe UEB_2stream <br> [Reference](http://doi.org/)** |
| NL_scatter | **TODO: Describe NL_scatter <br> [Reference](http://doi.org/)** |
| BeersLaw | **TODO: Describe BeersLaw <br> [Reference](http://doi.org/)** |

<a id="alb_method"></a>
## 36. alb_method
Albedo representation

| Option | Description |
|---|---|
| conDecay | **TODO: Describe conDecay <br> [Reference](http://doi.org/)** |
| varDecay | **TODO: Describe varDecay <br> [Reference](http://doi.org/)** |


<a id="spatial_gw"></a>
## 37. spatial_gw
Method for spatial representation of groundwater

| Option | Description |
|---|---|
| localColumn | **TODO: Describe localColumn <br> [Reference](http://doi.org/)** |
| singleBasin | **TODO: Describe singleBasin <br> [Reference](http://doi.org/)** |


<a id="subRouting"></a>
## 38. subRouting
Method for sub-grid routing

| Option | Description |
|---|---|
| timeDlay | **TODO: Describe timeDlay <br> [Reference](http://doi.org/)** |
| qInstant | **TODO: Describe qInstant <br> [Reference](http://doi.org/)** |


<a id="snowDenNew"></a>
## 39. snowDenNew
Method for new snow density

| Option | Description |
|---|---|
| hedAndPom | **An empirical calculation dependant on air temperature. <br> [Hedstom and Pomeroy, 1998](https://doi.org/10.1002/(SICI)1099-1085(199808/09)12:10/11<1611::AID-HYP684>3.0.CO;2-4)** |
| anderson | **TODO: Describe anderson <br> [Reference](http://doi.org/)** |
| pahaut_76 | **An empirical calculation dependant on air temperature and wind speed. <br> [Pahaut, 1976](http://doi.org/)** |
| constDens | **A constant new snow density of 330 kg/m^3 <br> [Reference](http://doi.org/)** |

<a id="snowUnload"></a>
## 40. snowUnload
Method for unloading snow from the canopy

| Option | Description |
|---|---|
| meltDripUnload | **Contains a temperature unloading function where the parameter *snowUnloadingCoeff* controls the exponential unloading rate and *ratioDrip2Unloading* is the ratio of liquid water drip from the canopy to snow unloading. <br> [Hedstom and Pomeroy, 1998](https://doi.org/10.1002/(SICI)1099-1085(199808/09)12:10/11<1611::AID-HYP684>3.0.CO;2-4) <br> [Storck et al. 2002]( https://doi.org/10.1029/2002WR001281)** |
| windUnload | **Contains temperature and wind dependent unloading functions. The rates of temperature and wind unloading are adjustable through parameters *rateTempUnloading* and *rateWindUnloading*. Both functions contain parameter thresholds for the minimum temperature and windspeed required for unloading.  <br> [Roesch et al. 2001](https://doi.org/10.1007/s003820100153)** |


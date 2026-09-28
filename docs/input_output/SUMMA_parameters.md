# SUMMA Parameters

<a id="params_local"></a>
## Local parameters
These are the spatially constant, HRU-level parameters set in the [local parameters file](SUMMA_input#infile_local_parameters) (or overridden per HRU in the [trial parameters file](SUMMA_input#infile_trial_parameters)). They correspond to the `iLookPARAM` members in `build/source/dshare/var_lookup.f90`, with names, descriptions and units as registered in `build/source/dshare/popMetadat.f90`.

| Variable | Description | Units |
|---|---|---|
| `upperBoundHead` | matric head at the upper boundary | m |
| `lowerBoundHead` | matric head at the lower boundary | m |
| `upperBoundTemp` | temperature of the upper boundary | K |
| `lowerBoundTemp` | temperature of the lower boundary | K |
| `tempCritRain` | critical temperature where precipitation is rain | K |
| `tempRangeTimestep` | temperature range over the time step | K |
| `frozenPrecipMultip` | frozen precipitation multiplier | - |
| `snowfrz_scale` | scaling parameter for the freezing curve for snow | K-1 |
| `fixedThermalCond_snow` | temporally constant thermal conductivity for snow | W m-1 K-1 |
| `albedoMax` | maximum snow albedo (single spectral band) | - |
| `albedoMinWinter` | minimum snow albedo during winter (single spectral band) | - |
| `albedoMinSpring` | minimum snow albedo during spring (single spectral band) | - |
| `albedoMaxVisible` | maximum snow albedo in the visible part of the spectrum | - |
| `albedoMinVisible` | minimum snow albedo in the visible part of the spectrum | - |
| `albedoMaxNearIR` | maximum snow albedo in the near infra-red part of the spectrum | - |
| `albedoMinNearIR` | minimum snow albedo in the near infra-red part of the spectrum | - |
| `albedoDecayRate` | albedo decay rate | s |
| `albedoSootLoad` | soot load factor | - |
| `albedoRefresh` | critical mass necessary for albedo refreshment | kg m-2 |
| `albedoFrznWatVisible` | albedo of frozen water in the visible part of the spectrum | - |
| `albedoFrznWatNearIR` | albedo of frozen water in the near infra-red part of the spectrum | - |
| `albedoOpenWatVisible` | albedo of open water in the visible part of the spectrum | - |
| `albedoOpenWatNearIR` | albedo of open water in the near infra-red part of the spectrum | - |
| `radExt_snow` | extinction coefficient for radiation penetration into snowpack | m-1 |
| `directScale` | scaling factor for fractional driect radiaion parameterization | - |
| `Frad_direct` | fraction direct solar radiation | - |
| `Frad_vis` | fraction radiation in visible part of spectrum | - |
| `newSnowDenMin` | minimum new snow density | kg m-3 |
| `newSnowDenMult` | multiplier for new snow density | kg m-3 |
| `newSnowDenScal` | scaling factor for new snow density | K |
| `constSnowDen` | Constant new snow density | kg m-3 |
| `newSnowDenAdd` | Pahaut 1976, additive factor for new snow density | kg m-3 |
| `newSnowDenMultTemp` | Pahaut 1976, multiplier for new snow density for air temperature | kg m-3 K-1 |
| `newSnowDenMultWind` | Pahaut 1976, multiplier for new snow density for wind speed | kg m-7/2 s-1/2 |
| `newSnowDenMultAnd` | Anderson 1976, multiplier for new snow density (Anderson func) | K-1 |
| `newSnowDenBase` | Anderson 1976, base value that is rasied to the (3/2) power | K |
| `densScalGrowth` | density scaling factor for grain growth | kg-1 m3 |
| `tempScalGrowth` | temperature scaling factor for grain growth | K-1 |
| `grainGrowthRate` | rate of grain growth | s-1 |
| `densScalOvrbdn` | density scaling factor for overburden pressure | kg-1 m3 |
| `tempScalOvrbdn` | temperature scaling factor for overburden pressure | K-1 |
| `baseViscosity` | viscosity coefficient at T=T_frz and snow density=0 | kg s m-2 |
| `Fcapil` | capillary retention (fraction of total pore volume) | - |
| `k_snow` | hydraulic conductivity of snow | m s-1 |
| `mw_exp` | exponent for meltwater flow | - |
| `z0Water` | roughness length of open water | m |
| `z0Ice` | roughness length of ice | m |
| `z0Snow` | roughness length of snow | m |
| `z0Soil` | roughness length of bare soil below the canopy | m |
| `z0Canopy` | roughness length of the canopy, only used if decision veg_traits==vegTypeTable | m |
| `zpdFraction` | zero plane displacement / canopy height | - |
| `critRichNumber` | critical value for the bulk Richardson number | - |
| `Louis79_bparam` | parameter in Louis (1979) stability function | - |
| `Louis79_cStar` | parameter in Louis (1979) stability function | - |
| `Mahrt87_eScale` | exponential scaling factor in the Mahrt (1987) stability function | - |
| `leafExchangeCoeff` | turbulent exchange coeff between canopy surface and canopy air | m s-(1/2) |
| `windReductionParam` | canopy wind reduction parameter | - |
| `glacierWindFactor` | wind speed increase to account for glacier katabatic wind profile | - |
| `glacierTempReduction` | air temperature decrease to account for down-glacier katabatic wind | - |
| `Kc25` | Michaelis-Menten constant for CO2 at 25 degrees C | umol mol-1 |
| `Ko25` | Michaelis-Menten constant for O2 at 25 degrees C | mol mol-1 |
| `Kc_qFac` | factor in the q10 function defining temperature controls on Kc | - |
| `Ko_qFac` | factor in the q10 function defining temperature controls on Ko | - |
| `kc_Ha` | activation energy for the Michaelis-Menten constant for CO2 | J mol-1 |
| `ko_Ha` | activation energy for the Michaelis-Menten constant for O2 | J mol-1 |
| `vcmax25_canopyTop` | potential carboxylation rate at 25 degrees C at the canopy top | umol co2 m-2 s-1 |
| `vcmax_qFac` | factor in the q10 function defining temperature controls on vcmax | - |
| `vcmax_Ha` | activation energy in the vcmax function | J mol-1 |
| `vcmax_Hd` | deactivation energy in the vcmax function | J mol-1 |
| `vcmax_Sv` | entropy term in the vcmax function | J mol-1 K-1 |
| `vcmax_Kn` | foliage nitrogen decay coefficient | - |
| `jmax25_scale` | scaling factor to relate jmax25 to vcmax25 | - |
| `jmax_Ha` | activation energy in the jmax function | J mol-1 |
| `jmax_Hd` | deactivation energy in the jmax function | J mol-1 |
| `jmax_Sv` | entropy term in the jmax function | J mol-1 K-1 |
| `fractionJ` | fraction of light lost by other than the chloroplast lamellae | - |
| `quantamYield` | quantam yield | mol e mol-1 q |
| `vpScaleFactor` | vapor pressure scaling factor in stomatal conductance function | Pa |
| `cond2photo_slope` | slope of conductance-photosynthesis relationship | - |
| `minStomatalConductance` | minimum stomatal conductance | umol H2O m-2 s-1 |
| `winterSAI` | stem area index prior to the start of the growing season | m2 m-2 |
| `summerLAI` | maximum leaf area index at the peak of the growing season | m2 m-2 |
| `rootScaleFactor1` | 1st scaling factor (a) in Y = 1 - 0.5*( exp(-aZ) + exp(-bZ) ) | m-1 |
| `rootScaleFactor2` | 2nd scaling factor (b) in Y = 1 - 0.5*( exp(-aZ) + exp(-bZ) ) | m-1 |
| `rootingDepth` | rooting depth | m |
| `rootDistExp` | exponent for the vertical distribution of root density | - |
| `plantWiltPsi` | matric head at wilting point | m |
| `soilStressParam` | parameter in the exponential soil stress function | - |
| `critSoilWilting` | critical vol. liq. water content when plants are wilting | - |
| `critSoilTranspire` | critical vol. liq. water content when transpiration is limited | - |
| `critAquiferTranspire` | critical aquifer storage value when transpiration is limited | m |
| `minStomatalResistance` | minimum stomatal resistance | s m-1 |
| `leafDimension` | characteristic leaf dimension | m |
| `heightCanopyTop` | height of top of the vegetation canopy above ground surface | m |
| `heightCanopyBottom` | height of bottom of the vegetation canopy above ground surface | m |
| `specificHeatVeg` | specific heat of vegetation | J kg-1 K-1 |
| `maxMassVegetation` | maximum mass of vegetation (full foliage) | kg m-2 |
| `throughfallScaleSnow` | scaling factor for throughfall (snow) | - |
| `throughfallScaleRain` | scaling factor for throughfall (rain) | - |
| `refInterceptCapSnow` | reference canopy interception capacity per unit leaf area (snow) | kg m-2 |
| `refInterceptCapRain` | canopy interception capacity per unit leaf area (rain) | kg m-2 |
| `snowUnloadingCoeff` | time constant for unloading of snow from the forest canopy | s-1 |
| `canopyDrainageCoeff` | time constant for drainage of liquid water from the forest canopy | s-1 |
| `ratioDrip2Unloading` | ratio of canopy drip to unloading of snow from the forest canopy | - |
| `canopyWettingFactor` | maximum wetted fraction of the canopy | - |
| `canopyWettingExp` | exponent in canopy wetting function | - |
| `minTempUnloading` | min temp for unloading in windySnow | K |
| `rateTempUnloading` | how quickly to unload due to temperature | K s |
| `minWindUnloading` | min wind speed for unloading in windySnow | m s-1 |
| `rateWindUnloading` | how quickly to unload due to wind | m |
| `soil_dens_intr` | intrinsic soil density | kg m-3 |
| `thCond_soil` | thermal conductivity of soil (includes quartz and other minerals) | W m-1 K-1 |
| `frac_sand` | fraction of sand | - |
| `frac_silt` | fraction of silt | - |
| `frac_clay` | fraction of clay | - |
| `theta_sat` | soil porosity | - |
| `theta_res` | volumetric residual water content | - |
| `vGn_alpha` | van Genuchten "alpha" parameter | m-1 |
| `vGn_n` | van Genuchten "n" parameter | - |
| `k_soil` | saturated hydraulic conductivity | m s-1 |
| `k_macropore` | saturated hydraulic conductivity for macropores | m s-1 |
| `fieldCapacity` | soil field capacity (vol liq water content when baseflow begins) | - |
| `wettingFrontSuction` | Green-Ampt wetting front suction | m |
| `theta_mp` | volumetric liquid water content when macropore flow begins | - |
| `mpExp` | empirical exponent in macropore flow equation | - |
| `kAnisotropic` | anisotropy factor for lateral hydraulic conductivity | - |
| `zScale_TOPMODEL` | TOPMODEL scaling factor used in lower boundary condition for soil | m |
| `compactedDepth` | depth where k_soil reaches the compacted value given by CH78 | m |
| `f_hydCond` | decay rate of hydraulic conductivity with depth, exp profile | m-1 |
| `aquiferBaseflowRate` | baseflow rate when aquifer storage = aquiferScaleFactor | m s-1 |
| `aquiferScaleFactor` | scaling factor for aquifer storage in the big bucket | m |
| `aquiferBaseflowExp` | baseflow exponent | - |
| `qSurfScale` | scaling factor in the surface runoff parameterization | - |
| `specificYield` | specific yield | - |
| `specificStorage` | specific storage coefficient | m-1 |
| `f_impede` | ice impedence factor | - |
| `soilIceScale` | scaling factor for depth of soil ice, used to get frozen fraction | m |
| `soilIceCV` | CV of depth of soil ice, used to get frozen fraction | - |
| `FUSE_Ac_max` | FUSE PRMS max saturated area | - |
| `FUSE_phi_tens` | FUSE PRMS tension storage fraction | - |
| `FUSE_b` | FUSE ARNO/VIC exponent | - |
| `FUSE_lambda` | FUSE TOPMODEL gamma distribution lambda parameter | m |
| `FUSE_chi` | FUSE TOPMODEL gamma distribution chi parameter | - |
| `FUSE_mu` | FUSE TOPMODEL gamma distribution mu parameter | m |
| `FUSE_n` | FUSE TOPMODEL exponent | - |
| `minwind` | minimum wind speed | m s-1 |
| `minstep` | minimum length of the time step homegrown | s |
| `maxstep` | maximum length of the time step (data window) | s |
| `be_steps` | number of equal substeps to dividing the data window for BE | - |
| `wimplicit` | weight assigned to the start-of-step fluxes ,homegrown, not currently used | - |
| `maxiter` | maximum number of iterations homegrown and kinsol | - |
| `relConvTol_liquid` | BE relative convergence tolerance for vol frac liq water homegrown | - |
| `absConvTol_liquid` | BE absolute convergence tolerance for vol frac liq water homegrown | - |
| `relConvTol_matric` | BE relative convergence tolerance for matric head homegrown | - |
| `absConvTol_matric` | BE absolute convergence tolerance for matric head homegrown | m |
| `relConvTol_energy` | BE relative convergence tolerance for energy homegrown | - |
| `absConvTol_energy` | BE absolute convergence tolerance for energy homegrown | J m-3 |
| `relConvTol_aquifr` | BE relative convergence tolerance for aquifer storage homegrown | - |
| `absConvTol_aquifr` | BE absolute convergence tolerance for aquifer storage homegrown | m |
| `relTolTempCas` | IDA relative error tolerance for canopy temperature state variable | - |
| `absTolTempCas` | IDA absolute error tolerance for canopy temperature state variable | - |
| `relTolTempVeg` | IDA relative error tolerance for vegitation temp state var | - |
| `absTolTempVeg` | IDA absolute error tolerance for vegitation temp state var | - |
| `relTolWatVeg` | IDA relative error tolerance for vegitation hydrology | - |
| `absTolWatVeg` | IDA absolute error tolerance for vegitation hydrology | - |
| `relTolTempSoilSnow` | IDA relative error tolerance for layers energy | - |
| `absTolTempSoilSnow` | IDA absolute error tolerance for layers energy | - |
| `relTolWatSnow` | IDA relative error tolerance for snow hydrology | - |
| `absTolWatSnow` | IDA absolute error tolerance for snow hydrology | - |
| `relTolMatric` | IDA relative error tolerance for matric head | - |
| `absTolMatric` | IDA absolute error tolerance for matric head | - |
| `relTolAquifr` | IDA relative error tolerance for aquifer hydrology | - |
| `absTolAquifr` | IDA absolute error tolerance for aquifer hydrology | - |
| `idaMaxInternalSteps` | maximum number of internal steps for IDA before tout | - |
| `idaInitStepSize` | initial step size for IDA | - |
| `idaMinStepSize` | minimum step size for IDA | - |
| `idaMaxStepSize` | maximum step size for IDA | - |
| `idaMaxErrTestFail` | maximum number of error test failures for IDA | - |
| `idaMaxOrder` | maximum order for IDA | - |
| `idaMaxDataWindowSteps` | maximum number of steps for IDA within one data window | - |
| `idaDetectEvents` | flag to turn on event detection in IDA, 0=off, 1=on | - |
| `zmin` | minimum layer depth | m |
| `zmax` | maximum layer depth | m |
| `zminLayer1` | minimum layer depth for the 1st (top) layer | m |
| `zminLayer2` | minimum layer depth for the 2nd layer | m |
| `zminLayer3` | minimum layer depth for the 3rd layer | m |
| `zminLayer4` | minimum layer depth for the 4th layer | m |
| `zminLayer5` | minimum layer depth for the 5th (bottom) layer | m |
| `zmaxLayer1_lower` | maximum layer depth for the 1st (top) layer when only 1 layer | m |
| `zmaxLayer2_lower` | maximum layer depth for the 2nd layer when only 2 layers | m |
| `zmaxLayer3_lower` | maximum layer depth for the 3rd layer when only 3 layers | m |
| `zmaxLayer4_lower` | maximum layer depth for the 4th layer when only 4 layers | m |
| `zmaxLayer1_upper` | maximum layer depth for the 1st (top) layer when > 1 layer | m |
| `zmaxLayer2_upper` | maximum layer depth for the 2nd layer when > 2 layers | m |
| `zmaxLayer3_upper` | maximum layer depth for the 3rd layer when > 3 layers | m |
| `zmaxLayer4_upper` | maximum layer depth for the 4th layer when > 4 layers | m |

<a id="params_basin"></a>
## Basin parameters
These are the spatially constant, GRU-level parameters set in the [basin parameters file](SUMMA_input#infile_basin_parameters). They correspond to the `iLookBPAR` members in `build/source/dshare/var_lookup.f90`.

| Variable | Description | Units |
|---|---|---|
| `basin__aquiferHydCond` | hydraulic conductivity of the aquifer | m s-1 |
| `basin__aquiferScaleFactor` | scaling factor for aquifer storage in the big bucket | m |
| `basin__aquiferBaseflowExp` | baseflow exponent for the big bucket | - |
| `routingGammaShape` | shape parameter in Gamma distribution used for sub-grid routing | - |
| `routingGammaScale` | scale parameter in Gamma distribution used for sub-grid routing | s |
| `glacStor_kIce` | storage coefficient glacier ice reservoir | s |
| `glacStor_kSnow` | storage coefficient glacier snow reservoir | s |
| `glacStor_kFirn` | storage coefficient glacier firn reservoir | s |
| `debrisConc` | englacial debris concentration | kg m-3 |
| `wallErosionRate` | glacier wall erosion rate input for debris advection | mm yr-1 |
| `debrisCritStress` | critical driving stress where debris slides on terminal wedge | Pa |
| `latMoraineWidth` | lateral moraine width or rockfall length | m |


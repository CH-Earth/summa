MODULE var_lookup
 ! defines named variables used to index array elements
 USE nrtype
 implicit none
 private
 ! define missing value
 integer(i4b),parameter     :: imiss = -999      ! used to initialize named variables
 ! ***********************************************************************************************************
 ! (0) define model decisions
 ! ***********************************************************************************************************
 type, public  ::  iLook_decision
  integer(i4b)    :: soilCatTbl       = 1  ! soil-category dateset
  integer(i4b)    :: vegeParTbl       = 2  ! vegetation category dataset
  integer(i4b)    :: num_method       = 3  ! choice of numerical method
  integer(i4b)    :: fDerivMeth       = 4  ! method used to calculate flux derivatives
  integer(i4b)    :: f_Richards       = 5  ! form of richards' equation
  integer(i4b)    :: groundwatr       = 6  ! choice of groundwater parameterization
  integer(i4b)    :: hc_profile       = 7  ! choice of hydraulic conductivity profile
  integer(i4b)    :: bcUpprTdyn       = 8  ! type of upper boundary condition for thermodynamics 
  integer(i4b)    :: bcLowrTdyn       = 9  ! type of lower boundary condition for thermodynamics
  integer(i4b)    :: bcUpprSoiH       = 10 ! type of upper boundary condition for soil hydrology
  integer(i4b)    :: bcLowrSoiH       = 11 ! type of lower boundary condition for soil hydrology
  integer(i4b)    :: astability       = 12 ! choice of stability function
  integer(i4b)    :: alb_method       = 13 ! choice of albedo representation
  integer(i4b)    :: compaction       = 14 ! choice of compaction routine
  integer(i4b)    :: thermlcond       = 15 ! choice of thermal conductivity representation
  integer(i4b)    :: subRouting       = 16 ! choice of method for sub-grid routing
 endtype iLook_decision
 ! ***********************************************************************************************************
 ! (1) define model time
 ! ***********************************************************************************************************
 type, public  ::  iLook_time
  integer(i4b)    :: iyyy             = 1  ! year
  integer(i4b)    :: im               = 2  ! month
  integer(i4b)    :: id               = 3  ! day
  integer(i4b)    :: ih               = 4  ! hour
  integer(i4b)    :: imin             = 5  ! minute
 endtype iLook_time
 ! ***********************************************************************************************************
 ! (2) define model forcing data
 ! ***********************************************************************************************************
 type, public  ::  iLook_force
  integer(i4b)    :: time             = 1  ! time since time reference       (s)
  integer(i4b)    :: pptrate          = 2  ! precipitation rate              (kg m-2 s-1)
  integer(i4b)    :: airtemp          = 3  ! air temperature                 (K)
  integer(i4b)    :: spechum          = 4  ! specific humidity               (g/g)
  integer(i4b)    :: windspd          = 5  ! windspeed                       (m/s)
  integer(i4b)    :: sw_down          = 6  ! downwelling shortwave radiaiton (W m-2)
  integer(i4b)    :: lw_down          = 7  ! downwelling longwave radiation  (W m-2)
  integer(i4b)    :: airpres          = 8  ! pressure                        (Pa)
 endtype iLook_force
 ! ***********************************************************************************************************
 ! (3) define local attributes
 ! ***********************************************************************************************************
 type, public  ::  iLook_attr
  integer(i4b)    :: latitude         = 1  ! latitude    (degrees north)
  integer(i4b)    :: longitude        = 2  ! longitude   (degrees east)
  integer(i4b)    :: elevation        = 3  ! elevation   (m)
 end type iLook_attr
 ! ***********************************************************************************************************
 ! (4) define local classification of veg, soil, etc. 
 ! ***********************************************************************************************************
 type, public  ::  iLook_type
  integer(i4b)    :: vegTypeIndex     = 1  ! index defining vegetation type (-)
  integer(i4b)    :: soilTypeIndex    = 2  ! index defining soil type (-)
  integer(i4b)    :: slopeTypeIndex   = 3  ! index defining slope (-)
 end type iLook_type
 ! ***********************************************************************************************************
 ! (5) define model parameters
 ! ***********************************************************************************************************
 type, public  ::  iLook_param
  ! boundary conditions
  integer(i4b)    :: upperBoundHead       = 1 ! matric head of the upper boundary (m)
  integer(i4b)    :: lowerBoundHead       = 2 ! matric head of the lower boundary (m)
  integer(i4b)    :: upperBoundTheta      = 3 ! volumetric liquid water content of the upper boundary (-)
  integer(i4b)    :: lowerBoundTheta      = 4 ! volumetric liquid water content of the lower boundary (-)
  integer(i4b)    :: upperBoundTemp       = 5 ! temperature of the upper boundary (K)
  integer(i4b)    :: lowerBoundTemp       = 6 ! temperature of the lower boundary (K)
  ! precipitation partitioning
  integer(i4b)    :: tempCritRain         = 7  ! critical temperature where precipitation is rain (K)
  integer(i4b)    :: tempRangeTimestep    = 8  ! temperature range over the time step (K)
  ! freezing curve for snow
  integer(i4b)    :: snowfrz_scale        = 9  ! scaling parameter for the freezing curve for snow (K-1)
  ! snow albedo
  integer(i4b)    :: snw_crit             = 10 ! critical mass necessary for albedo refreshment (kg m-2)
  integer(i4b)    :: alb_fresh            = 11 ! fresh snow albedo (-)
  integer(i4b)    :: alb_dry              = 12 ! minimum snow albedo during winter (-)
  integer(i4b)    :: alb_wet              = 13 ! minimum snow albedo during spring (-)
  integer(i4b)    :: alb_decay            = 14 ! temporal decay factor for snow albedo (s-1)
  integer(i4b)    :: alb_scale            = 15 ! albedo scaling factor (s)
  integer(i4b)    :: soot_load            = 16 ! temporal decay in snow albedo associated with the soot load (days-1)
  ! radiation transfer within snow
  integer(i4b)    :: radExt_snow          = 17 ! extinction coefficient for radiation penetration into the snowpack (m-1)
  integer(i4b)    :: Frad_vis             = 18 ! fraction of radiation in the visible part of the spectrum (-)
  ! new snow density
  integer(i4b)    :: newSnowDenMin        = 19 ! minimum new snow density (kg m-3)  
  integer(i4b)    :: newSnowDenMult       = 20 ! multiplier for new snow density (kg m-3)
  integer(i4b)    :: newSnowDenScal       = 21 ! scaling factor for new snow density (K)
  ! snow compaction
  integer(i4b)    :: densScalGrowth       = 22 ! density scaling factor for grain growth (kg-1 m3)
  integer(i4b)    :: tempScalGrowth       = 23 ! temperature scaling factor for grain growth (K-1)
  integer(i4b)    :: grainGrowthRate      = 24 ! rate of grain growth (s-1)
  integer(i4b)    :: densScalOvrbdn       = 25 ! density scaling factor for overburden pressure (kg-1 m3)
  integer(i4b)    :: tempScalOvrbdn       = 26 ! temperature scaling factor for overburden pressure (K-1)
  integer(i4b)    :: base_visc            = 27 ! viscosity coefficient at T=T_frz and snow density=0  (kg s m-2)
  ! water flow within snow
  integer(i4b)    :: Fcapil               = 28 ! capillary retention as a fraction of the total pore volume (-)
  integer(i4b)    :: k_snow               = 29 ! hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
  integer(i4b)    :: mw_exp               = 30 ! exponent for meltwater flow (-)
  ! turbulent heat fluxes
  integer(i4b)    :: mheight              = 31 ! measurement height (m)
  integer(i4b)    :: zon                  = 32 ! roughness length (m)
  integer(i4b)    :: c_star               = 33 ! parameter in Louis (1979) stability function
  integer(i4b)    :: bparam               = 34 ! parameter in Louis (1979) stability function
  integer(i4b)    :: Mahrt_m              = 35 ! the m parameter from the Mahrt (1987) stability function
  ! vegetation properties
  integer(i4b)    :: rootingDepth         = 36 ! rooting depth (m)
  integer(i4b)    :: rootDistExp          = 37 ! exponent controlling the vertical distribution of root density (-)
  integer(i4b)    :: minStomatalResist    = 38 ! minimum stomatal resistance (s m-1)
  integer(i4b)    :: maxStomatalResist    = 39 ! maximum stomatal resistance (s m-1)
  integer(i4b)    :: plantWiltPsi         = 40 ! critical matric head when stomatal resitance 2 x min (m)
  integer(i4b)    :: plantWiltExp         = 41 ! empirical exponent in plant wilting factor expression (-)
  integer(i4b)    :: critAquiferTranspire = 42 ! critical aquifer storage value when transpiration is limited (m)
  ! soil properties
  integer(i4b)    :: soilAlbedo           = 43 ! soil albedo (-)
  integer(i4b)    :: soil_dens_intr       = 44 ! intrinsic soil density (kg m-3)
  integer(i4b)    :: frac_sand            = 45 ! fraction of sand (-)
  integer(i4b)    :: frac_silt            = 46 ! fraction of silt (-)
  integer(i4b)    :: frac_clay            = 47 ! fraction of clay (-)
  integer(i4b)    :: theta_sat            = 48 ! porosity (-)
  integer(i4b)    :: theta_res            = 49 ! volumetric residual water content (-)
  integer(i4b)    :: vGn_alpha            = 50 ! van Genuchten "alpha" parameter (m-1)
  integer(i4b)    :: vGn_n                = 51 ! van Genuchten "n" parameter (-)
  integer(i4b)    :: k_soil               = 52 ! hydraulic conductivity of soil (m s-1) 
  integer(i4b)    :: kAnisotropic         = 53 ! anisotropy factor for lateral hydraulic conductivity (-)
  integer(i4b)    :: zScale_TOPMODEL      = 54 ! scale factor for TOPMODEL-ish baseflow parameterization (m)
  integer(i4b)    :: compactedDepth       = 55 ! depth where k_soil reaches the compacted value given by CH78 (m)
  integer(i4b)    :: bpar_VIC             = 56 ! b-parameter in the VIC surface runoff parameterization (-)
  integer(i4b)    :: specificYield        = 57 ! specific yield (-)
  integer(i4b)    :: specificStorage      = 58 ! specific storage coefficient (m-1)
  integer(i4b)    :: aquiferScaleFactor   = 59 ! scaling factor for aquifer storage in the big bucket (m)
  integer(i4b)    :: bucketBaseflowExp    = 60 ! baseflow exponent for the big bucket (-)
  integer(i4b)    :: f_impede             = 61 ! ice impedence factor (-)
  ! within-grid routing
  integer(i4b)    :: routingGammaShape    = 62 ! shape parameter in Gamma distribution used for sub-grid routing (-)
  integer(i4b)    :: routingGammaScale    = 63 ! scale parameter in Gamma distribution used for sub-grid routing (s)
  ! algorithmic control parameters
  integer(i4b)    :: minwind              = 64 ! minimum wind speed (m s-1)
  integer(i4b)    :: minstep              = 65 ! minimum length of the time step
  integer(i4b)    :: maxstep              = 66 ! maximum length of the time step
  integer(i4b)    :: wimplicit            = 67 ! weight assigned to the start-of-step fluxes
  integer(i4b)    :: maxiter              = 68 ! maximum number of iteration 
  integer(i4b)    :: relConvTol_liquid    = 69 ! relative convergence tolerance for vol frac liq water (-)
  integer(i4b)    :: absConvTol_liquid    = 70 ! absolute convergence tolerance for vol frac liq water (-)
  integer(i4b)    :: relConvTol_matric    = 71 ! relative convergence tolerance for matric head (-)
  integer(i4b)    :: absConvTol_matric    = 72 ! absolute convergence tolerance for matric head (m)
  integer(i4b)    :: relConvTol_energy    = 73 ! relative convergence tolerance for energy (-)
  integer(i4b)    :: absConvTol_energy    = 74 ! absolute convergence tolerance for energy (J m-3)
  integer(i4b)    :: relConvTol_aquifr    = 75 ! relative convergence tolerance for aquifer storage (-)
  integer(i4b)    :: absConvTol_aquifr    = 76 ! absolute convergence tolerance for aquifer storage (J m-3)
  integer(i4b)    :: zmin                 = 77 ! minimum layer depth (m)
  integer(i4b)    :: zmax                 = 78 ! maximum layer depth (m)
 endtype ilook_param
 ! ***********************************************************************************************************
 ! (6) define model variables
 ! ***********************************************************************************************************
 type, public :: iLook_mvar
  ! define timestep-average fluxes for a few key variables
  integer(i4b)    :: averageMassLiquid            = 1   ! evaporation or dew (kg m-2 s-1)
  integer(i4b)    :: averageMassSolid             = 2   ! sublimation or frost (kg m-2 s-1)
  integer(i4b)    :: averageRainPlusMelt          = 3   ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  integer(i4b)    :: averageSurfaceRunoff         = 4   ! surface runoff (m s-1)
  integer(i4b)    :: averageSoilInflux            = 5   ! influx of water at the top of the soil profile (m s-1)
  integer(i4b)    :: averageSoilBaseflow          = 6   ! total baseflow from throughout the soil profile (m s-1)
  integer(i4b)    :: averageSoilDrainage          = 7   ! drainage from the bottom of the soil profile (m s-1)
  integer(i4b)    :: averageSoilEjection          = 8   ! ejected water from the soil matrix (m s-1)
  integer(i4b)    :: averageAquiferRecharge       = 9   ! recharge to the aquifer (m s-1)
  integer(i4b)    :: averageAquiferBaseflow       = 10  ! baseflow from the aquifer (m s-1)
  ! define scalar variables
  integer(i4b)    :: scalarSwDownVis              = 11  ! downwelling shortwave radiation in visible part of spectrum (W m-2)
  integer(i4b)    :: scalarSwDownNir              = 12  ! downwelling shortwave radiation in near-infrared part of spectrum (W m-2)
  integer(i4b)    :: scalarTwetbulb               = 13  ! wet bulb temperature (K)
  integer(i4b)    :: scalarRainfall               = 14  ! computed rainfall rate (kg m-2 s-1)
  integer(i4b)    :: scalarSnowfall               = 15  ! computed snowfall rate (kg m-2 s-1)
  integer(i4b)    :: scalarSnowfallTemp           = 16  ! temperature of fresh snow (K)
  integer(i4b)    :: scalarSurfaceTemp            = 17  ! temperature of the top layer (K)
  integer(i4b)    :: scalarAlbedo                 = 18  ! albedo of the surface, soil or snow (-)
  integer(i4b)    :: scalarSnowDepth              = 19  ! total snow depth (m)
  integer(i4b)    :: scalarSWE                    = 20  ! snow water equivalent (kg m-2)
  integer(i4b)    :: scalarSfcMeltPond            = 21  ! ponded water caused by melt of the "snow without a layer" (kg m-2)
  integer(i4b)    :: scalarAquiferStorage         = 22  ! relative aquifer storage -- above bottom of the soil profile (m)
  integer(i4b)    :: scalarWaterTableDepth        = 23  ! depth of the water table (m)
  integer(i4b)    :: scalarExCoef                 = 24  ! turbulent exchange coefficient (-)
  integer(i4b)    :: scalarExSen                  = 25  ! exchange factor for sensible heat (W m-2 K-1)
  integer(i4b)    :: scalarExLat                  = 26  ! exchange factor for latent heat (W m-2) 
  integer(i4b)    :: scalarSenHeat                = 27  ! sensible heat flux at the surface (W m-2) 
  integer(i4b)    :: scalarLatHeat                = 28  ! latent heat flux at the surface (W m-2)
  integer(i4b)    :: scalarPotentialET            = 29  ! potential evapotranspiration (kg m-2 s-1)
  integer(i4b)    :: scalarTranspireLim           = 30  ! aggregate soil moist & vegetation limit on transpiration (-)
  integer(i4b)    :: scalarTranspireLimAqfr       = 31  ! soil moist & vegetation limit on transpiration in the aquifer (-)
  integer(i4b)    :: scalarMassSolid              = 32  ! sublimation or frost (kg m-2 s-1)
  integer(i4b)    :: scalarMassLiquid             = 33  ! evaporation or dew (kg m-2 s-1)
  integer(i4b)    :: scalarHeatPrecip             = 34  ! sensible heat of precipitation (W m-2)
  integer(i4b)    :: scalarRainPlusMelt           = 35  ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  integer(i4b)    :: scalarSurfaceRunoff          = 36  ! surface runoff (m s-1)
  integer(i4b)    :: scalarInitAquiferRecharge    = 37  ! recharge to the aquifer at the start of the step (m s-1)
  integer(i4b)    :: scalarAquiferRecharge        = 38  ! recharge to the aquifer (m s-1)
  integer(i4b)    :: scalarInitAquiferTranspire   = 39  ! transpiration from the aquifer at the start of the step (m s-1)
  integer(i4b)    :: scalarAquiferTranspire       = 40  ! transpiration from the aquifer (m s-1)
  integer(i4b)    :: scalarInitAquiferBaseflow    = 41  ! baseflow from the aquifer at the start of the step (m s-1)
  integer(i4b)    :: scalarAquiferBaseflow        = 42  ! baseflow from the aquifer (m s-1)
  integer(i4b)    :: scalarSoilInflux             = 43  ! influx of water at the top of the soil profile (m s-1)
  integer(i4b)    :: scalarSoilBaseflow           = 44  ! total baseflow from throughout the soil profile (m s-1)
  integer(i4b)    :: scalarSoilDrainage           = 45  ! drainage from the bottom of the soil profile (m s-1) 
  integer(i4b)    :: scalarSoilEjection           = 46  ! total water ejected from all soil layers (m s-1) 
  integer(i4b)    :: scalarSoilWatBalError        = 47  ! error in the total soil water balance (kg m-2)
  integer(i4b)    :: scalarAquiferBalError        = 48  ! error in the aquifer water balance (kg m-2)
  integer(i4b)    :: scalarTotalSoilLiq           = 49  ! total mass of liquid water in the soil (kg m-2)
  integer(i4b)    :: scalarTotalSoilIce           = 50  ! total mass of ice in the soil (kg m-2)
  ! define variables at the mid-point of each layer -- domain geometry
  integer(i4b)    :: mLayerDepth                  = 51  ! depth of each layer (m)
  integer(i4b)    :: mLayerHeight                 = 52  ! height at the mid-point of each layer (m)
  integer(i4b)    :: mLayerRootDensity            = 53  ! fraction of roots in each soil layer (-)
  ! define variables at the mid-point of each layer -- coupled energy and mass
  integer(i4b)    :: mLayerTemp                   = 54  ! temperature of each layer (K)
  integer(i4b)    :: mLayerVolFracAir             = 55  ! volumetric fraction of air in each layer (-)
  integer(i4b)    :: mLayerVolFracIce             = 56  ! volumetric fraction of ice water in each layer (-)
  integer(i4b)    :: mLayerVolFracLiq             = 57  ! volumetric fraction of liquid water in each layer (-)
  integer(i4b)    :: mLayerVolHtCapBulk           = 58  ! volumetric heat capacity in each layer (J m-3 K-1)
  integer(i4b)    :: mLayerTcrit                  = 59  ! critical soil temperature above which all water is unfrozen (K)
  integer(i4b)    :: mLayerdTheta_dTk             = 60  ! derivative in volumetric liquid water content wrt temperature (K-1)
  integer(i4b)    :: mLayerThermalC               = 61  ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
  integer(i4b)    :: mLayerRadCondFlux            = 62  ! temporal derivative in energy from radiative and conductive flux (J m-2 s-1)  
  integer(i4b)    :: mLayerMeltFreeze             = 63  ! melt/freeze in each layer (kg m-3)
  integer(i4b)    :: mLayerSatHydCond             = 64  ! saturated hydraulic conductivity in each layer (m s-1)
  integer(i4b)    :: mLayerMatricHead             = 65  ! matric head of water in the soil (m)
  integer(i4b)    :: mLayerdTheta_dPsi            = 66  ! derivative in the soil water characteristic (m-1)
  integer(i4b)    :: mLayerdPsi_dTheta            = 67  ! derivative in the soil water characteristic (m)
  integer(i4b)    :: mLayerThetaResid             = 68  ! residual volumetric water content in each snow layer (-)
  integer(i4b)    :: mLayerPoreSpace              = 69  ! total pore space in each snow layer (-)
  integer(i4b)    :: mLayerInfilFreeze            = 70  ! volumetric ice content increase by freezing infiltrating flux (-)
  integer(i4b)    :: mLayerTranspireLim           = 71  ! soil moist & veg limit on transpiration for each layer (-)
  integer(i4b)    :: mLayerInitTranspire          = 72  ! transpiration loss from each soil layer at the start of the step (kg m-2 s-1)
  integer(i4b)    :: mLayerTranspire              = 73  ! transpiration loss from each soil layer (kg m-2 s-1)
  integer(i4b)    :: mLayerInitEjectWater         = 74  ! water ejected from each soil layer at the start-of-step (m s-1)
  integer(i4b)    :: mLayerEjectWater             = 75  ! water ejected from each soil layer (m s-1)
  integer(i4b)    :: mLayerInitBaseflow           = 76  ! baseflow from each soil layer at the start of the time step (m s-1)
  integer(i4b)    :: mLayerBaseflow               = 77  ! baseflow from each soil layer (m s-1)
  ! define variables at the interface of each layer
  integer(i4b)    :: iLayerHeight                 = 78  ! height of the layer interface; top of soil = 0 (m)
  integer(i4b)    :: iLayerThermalC               = 79  ! thermal conductivity at the interface of each layer (W m-1 K-1)
  integer(i4b)    :: iLayerInitNrgFlux            = 80  ! energy flux at layer interfaces at the start of the time step (W m-2) 
  integer(i4b)    :: iLayerNrgFlux                = 81  ! energy flux at layer interfaces at the end of the time step (W m-2) 
  integer(i4b)    :: iLayerSatHydCond             = 82  ! saturated hydraulic conductivity at each layer interface (m s-1)
  integer(i4b)    :: iLayerInitLiqFluxSnow        = 83  ! liquid flux at snow layer interfaces at the start of the time step (m s-1) 
  integer(i4b)    :: iLayerInitLiqFluxSoil        = 84  ! liquid flux at soil layer interfaces at the start of the time step (m s-1) 
  integer(i4b)    :: iLayerLiqFluxSnow            = 85  ! liquid flux at snow layer interfaces at the end of the time step (m s-1)
  integer(i4b)    :: iLayerLiqFluxSoil            = 86  ! liquid flux at soil layer interfaces at the end of the time step (m s-1) 
  ! define variables for runoff
  integer(i4b)    :: routingRunoffFuture          = 87  ! runoff in future time steps (m s-1)
  integer(i4b)    :: routingFractionFuture        = 88  ! fraction of runoff in future time steps (-)
  integer(i4b)    :: averageInstantRunoff         = 89  ! instantaneous runoff (m s-1)
  integer(i4b)    :: averageRoutedRunoff          = 90  ! routed runoff (m s-1)
  ! define derived variables
  integer(i4b)    :: scalarExNeut                 = 91  ! exchange coefficient in neutral conditions
  integer(i4b)    :: scalarBprime                 = 92  ! stable b parameter in Louis (1979) stability function
  integer(i4b)    :: scalarCparam                 = 93  ! c parameter in Louis (1979) stability function
  integer(i4b)    :: scalarVGn_m                  = 94  ! van Genuchten "m" parameter (-)
  integer(i4b)    :: scalarKappa                  = 95  ! constant in the freezing curve function (m K-1)
  integer(i4b)    :: scalarVolHtCap_air           = 96  ! volumetric heat capacity air         (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_ice           = 97  ! volumetric heat capacity ice         (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_soil          = 98  ! volumetric heat capacity dry soil    (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_water         = 99  ! volumetric heat capacity liquid wat  (J m-3 K-1)
  integer(i4b)    :: scalarLambda_drysoil         = 100 ! thermal conductivity of dry soil     (W m-1)
  integer(i4b)    :: scalarLambda_wetsoil         = 101 ! thermal conductivity of wet soil     (W m-1)
  integer(i4b)    :: scalarVolLatHt_fus           = 102 ! volumetric latent heat of fusion     (J m-3)
  integer(i4b)    :: scalarAquiferRootFrac        = 103 ! fraction of roots below the soil profile (-)
  ! define NOAH-MP vegetation variables
  integer(i4b)    :: scalarCanopyHeight           = 104  ! height of the top of the canopy layer (m)
  integer(i4b)    :: scalarVegetationTemp         = 105  ! vegetation temperature (K)
  integer(i4b)    :: scalarRootZoneTemp           = 106  ! average temperature of the root zone (K)
  integer(i4b)    :: scalarLAI                    = 107  ! one-sided leaf area index (m2 m-2)
  integer(i4b)    :: scalarSAI                    = 108  ! one-sided stem area index (m2 m-2)
  integer(i4b)    :: scalarEffectiveLAI           = 119  ! effective leaf area index after burial by snow (m2 m-2)
  integer(i4b)    :: scalarEffectiveSAI           = 110  ! effective stem area index after burial by snow(m2 m-2)
  integer(i4b)    :: scalarGrowingSeasonIndex     = 111  ! growing season index (0=off, 1=on)
 endtype iLook_mvar

 ! ***********************************************************************************************************
 ! (6) define model indices
 ! ***********************************************************************************************************
 type, public :: iLook_index
  integer(i4b)    :: nSnow                 = 1  ! number of snow layers
  integer(i4b)    :: nSoil                 = 2  ! number of soil layers
  integer(i4b)    :: nLayers               = 3  ! total number of layers
  integer(i4b)    :: midSnowStartIndex     = 4  ! start index of the midSnow vector for a given timestep
  integer(i4b)    :: midSoilStartIndex     = 5  ! start index of the midSoil vector for a given timestep
  integer(i4b)    :: midTotoStartIndex     = 6  ! start index of the midToto vector for a given timestep
  integer(i4b)    :: ifcSnowStartIndex     = 7  ! start index of the ifcSnow vector for a given timestep
  integer(i4b)    :: ifcSoilStartIndex     = 8  ! start index of the ifcSoil vector for a given timestep
  integer(i4b)    :: ifcTotoStartIndex     = 9  ! start index of the ifcToto vector for a given timestep
  integer(i4b)    :: layerType             = 10 ! type of layer (soil or snow)
 endtype iLook_index
 ! ***********************************************************************************************************
 ! (X) define data structures and maximum number of variables of each type
 ! ***********************************************************************************************************
 ! define look-up structures
 type(iLook_decision),public,parameter :: iLookDECISIONS=iLook_decision( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,&
                                                                        11,12,13,14,15,16)
 type(iLook_time),    public,parameter :: iLookTIME     =iLook_time    ( 1, 2, 3, 4, 5)
 type(iLook_force),   public,parameter :: iLookFORCE    =iLook_force   ( 1, 2, 3, 4, 5, 6, 7, 8)
 type(iLook_attr),    public,parameter :: iLookATTR     =iLook_attr    ( 1, 2, 3)
 type(iLook_type),    public,parameter :: iLookTYPE     =iLook_type    ( 1, 2, 3)
 type(iLook_param),   public,parameter :: iLookPARAM    =iLook_param   ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,&
                                                                        11,12,13,14,15,16,17,18,19,20,&
                                                                        21,22,23,24,25,26,27,28,29,30,&
                                                                        31,32,33,34,35,36,37,38,39,40,&
                                                                        41,42,43,44,45,46,47,48,49,50,&
                                                                        51,52,53,54,55,56,57,58,59,60,&
                                                                        61,62,63,64,65,66,67,68,69,70,&
                                                                        71,72,73,74,75,76,77,78)
 type(iLook_mvar),    public,parameter :: iLookMVAR     =ilook_mvar    (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,&
                                                                         21, 22, 23, 24, 25, 26, 27, 28, 29, 30,&
                                                                         31, 32, 33, 34, 35, 36, 37, 38, 39, 40,&
                                                                         41, 42, 43, 44, 45, 46, 47, 48, 49, 50,&
                                                                         51, 52, 53, 54, 55, 56, 57, 58, 59, 60,&
                                                                         61, 62, 63, 64, 65, 66, 67, 68, 69, 70,&
                                                                         71, 72, 73, 74, 75, 76, 77, 78, 79, 80,&
                                                                         81, 82, 83, 84, 85, 86, 87, 88, 89, 90,&
                                                                         91, 92, 93, 94, 95, 96, 97, 98, 99,100,&
                                                                        101,102,103,104,105,106,107,108,109,110,&
                                                                        111)
 type(iLook_index),   public,parameter :: iLookINDEX    =ilook_index   ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10)
 ! define maximum number of variables of each type
 integer(i4b),parameter,public :: maxvarDecisions= 16
 integer(i4b),parameter,public :: maxvarTime     = 5
 integer(i4b),parameter,public :: maxvarForc     = 8
 integer(i4b),parameter,public :: maxvarAttr     = 3
 integer(i4b),parameter,public :: maxvarType     = 3
 integer(i4b),parameter,public :: maxvarMpar     = 78
 integer(i4b),parameter,public :: maxvarMvar     = 111
 integer(i4b),parameter,public :: maxvarIndx     = 10
 ! ***********************************************************************************************************
END MODULE var_lookup

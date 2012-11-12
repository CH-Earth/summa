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
  integer(i4b)    :: num_method       = 1  ! choice of numerical method
  integer(i4b)    :: fDerivMeth       = 2  ! method used to calculate flux derivatives
  integer(i4b)    :: f_Richards       = 3  ! form of richards' equation
  integer(i4b)    :: groundwatr       = 4  ! choice of groundwater parameterization
  integer(i4b)    :: hc_profile       = 5  ! choice of hydraulic conductivity profile
  integer(i4b)    :: bcUpprTdyn       = 6  ! type of upper boundary condition for thermodynamics 
  integer(i4b)    :: bcLowrTdyn       = 7  ! type of lower boundary condition for thermodynamics
  integer(i4b)    :: bcUpprSoiH       = 8  ! type of upper boundary condition for soil hydrology
  integer(i4b)    :: bcLowrSoiH       = 9  ! type of lower boundary condition for soil hydrology
  integer(i4b)    :: astability       = 10  ! choice of stability function
  integer(i4b)    :: compaction       = 11 ! choice of compaction routine
  integer(i4b)    :: thermlcond       = 12 ! choice of thermal conductivity representation
  integer(i4b)    :: alb_method       = 13 ! choice of albedo representation
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
 ! (3) define model parameters
 ! ***********************************************************************************************************
 type, public  ::  iLook_param
  ! boundary conditions
  integer(i4b)    :: upperBoundHead     = 1 ! matric head of the upper boundary (m)
  integer(i4b)    :: lowerBoundHead     = 2 ! matric head of the lower boundary (m)
  integer(i4b)    :: upperBoundTheta    = 3 ! volumetric liquid water content of the upper boundary (-)
  integer(i4b)    :: lowerBoundTheta    = 4 ! volumetric liquid water content of the lower boundary (-)
  integer(i4b)    :: upperBoundTemp     = 5 ! temperature of the upper boundary (K)
  integer(i4b)    :: lowerBoundTemp     = 6 ! temperature of the lower boundary (K)
  ! precipitation partitioning
  integer(i4b)    :: tempCritRain       = 7  ! critical temperature where precipitation is rain (K)
  integer(i4b)    :: tempRangeTimestep  = 8  ! temperature range over the time step (K)
  ! freezing curve for snow
  integer(i4b)    :: snowfrz_scale      = 9  ! scaling parameter for the freezing curve for snow (K-1)
  ! snow albedo
  integer(i4b)    :: snw_crit           = 10 ! critical mass necessary for albedo refreshment (kg m-2)
  integer(i4b)    :: alb_fresh          = 11 ! fresh snow albedo (-)
  integer(i4b)    :: alb_dry            = 12 ! minimum snow albedo during winter (-)
  integer(i4b)    :: alb_wet            = 13 ! minimum snow albedo during spring (-)
  integer(i4b)    :: alb_decay          = 14 ! temporal decay factor for snow albedo (s-1)
  integer(i4b)    :: alb_scale          = 15 ! albedo scaling factor (s)
  integer(i4b)    :: soot_load          = 16 ! temporal decay in snow albedo associated with the soot load (days-1)
  ! radiation transfer within snow
  integer(i4b)    :: rad_ext            = 17 ! extinction coefficient for radiation penetration (m-1)
  integer(i4b)    :: Fabs_vis           = 18 ! fraction of absorbed radiation in the visible part of the spectrum (-)
  ! turbulent heat fluxes
  integer(i4b)    :: mheight            = 19 ! measurement height (m)
  integer(i4b)    :: zon                = 20 ! roughness length (m)
  integer(i4b)    :: c_star             = 21 ! parameter in Louis (1979) stability function
  integer(i4b)    :: bparam             = 22 ! parameter in Louis (1979) stability function
  integer(i4b)    :: Mahrt_m            = 23 ! the m parameter from the Mahrt (1987) stability function
  ! water flow within snow
  integer(i4b)    :: Fcapil             = 24 ! capillary retention as a fraction of the total pore volume (-)
  integer(i4b)    :: k_snow             = 25 ! hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
  integer(i4b)    :: mw_exp             = 26 ! exponent for meltwater flow (-)
  ! soil properties
  integer(i4b)    :: soilAlbedo         = 27 ! soil albedo (-)
  integer(i4b)    :: soil_dens_intr     = 28 ! intrinsic soil density (kg m-3)
  integer(i4b)    :: frac_sand          = 29 ! fraction of sand (-)
  integer(i4b)    :: frac_silt          = 30 ! fraction of silt (-)
  integer(i4b)    :: frac_clay          = 31 ! fraction of clay (-)
  integer(i4b)    :: theta_sat          = 32 ! porosity (-)
  integer(i4b)    :: theta_res          = 33 ! volumetric residual water content (-)
  integer(i4b)    :: vGn_alpha          = 34 ! van Genuchten "alpha" parameter (m-1)
  integer(i4b)    :: vGn_n              = 35 ! van Genuchten "n" parameter (-)
  integer(i4b)    :: k_soil             = 36 ! hydraulic conductivity of soil (m s-1) 
  integer(i4b)    :: kAnisotropic       = 37 ! anisotropy factor for lateral hydraulic conductivity (-)
  integer(i4b)    :: zScale_TOPMODEL    = 38 ! scale factor for TOPMODEL-ish baseflow parameterization (m)
  integer(i4b)    :: compactedDepth     = 39 ! depth where k_soil reaches the compacted value given by CH78 (m)
  integer(i4b)    :: bpar_VIC           = 40 ! b-parameter in the VIC surface runoff parameterization (-)
  integer(i4b)    :: specificYield      = 41 ! specific yield (-)
  integer(i4b)    :: specificStorage    = 42 ! specific storage coefficient (m-1)
  integer(i4b)    :: f_impede           = 43 ! ice impedence factor (-)
  ! vegetation properties
  integer(i4b)    :: rootingDepth       = 44 ! rooting depth (m)
  integer(i4b)    :: rootDistExp        = 45 ! exponent controlling the vertical distribution of root density (-)
  integer(i4b)    :: LAI                = 46 ! leaf area index (m2 m-2) 
  integer(i4b)    :: minStomatalResist  = 47 ! minimum stomatal resistance (s m-1)
  integer(i4b)    :: maxStomatalResist  = 48 ! maximum stomatal resistance (s m-1)
  integer(i4b)    :: plantWiltPsi       = 49 ! critical matric head when stomatal resitance 2 x min (m)
  integer(i4b)    :: plantWiltExp       = 50 ! empirical exponent in plant wilting factor expression (-)
  ! new snow density
  integer(i4b)    :: newSnowDenMin      = 51 ! minimum new snow density (kg m-3)  
  integer(i4b)    :: newSnowDenMult     = 52 ! multiplier for new snow density (kg m-3)
  integer(i4b)    :: newSnowDenScal     = 53 ! scaling factor for new snow density (K)
  ! snow compaction
  integer(i4b)    :: densScalGrowth     = 54 ! density scaling factor for grain growth (kg-1 m3)
  integer(i4b)    :: tempScalGrowth     = 55 ! temperature scaling factor for grain growth (K-1)
  integer(i4b)    :: grainGrowthRate    = 56 ! rate of grain growth (s-1)
  integer(i4b)    :: densScalOvrbdn     = 57 ! density scaling factor for overburden pressure (kg-1 m3)
  integer(i4b)    :: tempScalOvrbdn     = 58 ! temperature scaling factor for overburden pressure (K-1)
  integer(i4b)    :: base_visc          = 59 ! viscosity coefficient at T=T_frz and snow density=0  (kg s m-2)
  ! algorithmic control parameters
  integer(i4b)    :: minwind            = 60 ! minimum wind speed (m s-1)
  integer(i4b)    :: minstep            = 61 ! minimum length of the time step
  integer(i4b)    :: maxstep            = 62 ! maximum length of the time step
  integer(i4b)    :: wimplicit          = 63 ! weight assigned to the start-of-step fluxes
  integer(i4b)    :: maxiter            = 64 ! maximum number of iteration 
  integer(i4b)    :: relConvTol_liquid  = 65 ! relative convergence tolerance for vol frac liq water (-)
  integer(i4b)    :: absConvTol_liquid  = 66 ! absolute convergence tolerance for vol frac liq water (-)
  integer(i4b)    :: relConvTol_matric  = 67 ! relative convergence tolerance for matric head (-)
  integer(i4b)    :: absConvTol_matric  = 68 ! absolute convergence tolerance for matric head (m)
  integer(i4b)    :: relConvTol_energy  = 69 ! relative convergence tolerance for energy (-)
  integer(i4b)    :: absConvTol_energy  = 70 ! absolute convergence tolerance for energy (J m-3)
  integer(i4b)    :: zmin               = 71 ! minimum layer depth (m)
  integer(i4b)    :: zmax               = 72 ! maximum layer depth (m)
 endtype ilook_param
 ! ***********************************************************************************************************
 ! (4) define model variables
 ! ***********************************************************************************************************
 type, public :: iLook_mvar
  ! define timestep-average fluxes for a few key variables
  integer(i4b)    :: averageMassLiquid          = 1  ! evaporation or dew (kg m-2 s-1)
  integer(i4b)    :: averageMassSolid           = 2  ! sublimation or frost (kg m-2 s-1)
  integer(i4b)    :: averageRainPlusMelt        = 3  ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  integer(i4b)    :: averageSurfaceRunoff       = 4  ! surface runoff (m s-1)
  integer(i4b)    :: averageSoilInflux          = 5  ! influx of water at the top of the soil profile (m s-1)
  integer(i4b)    :: averageSoilDrainage        = 6  ! drainage from the bottom of the soil profile (m s-1)
  ! define scalar variables
  integer(i4b)    :: scalarTwetbulb             = 7  ! wet bulb temperature (K)
  integer(i4b)    :: scalarRainfall             = 8  ! computed rainfall rate (kg m-2 s-1)
  integer(i4b)    :: scalarSnowfall             = 9  ! computed snowfall rate (kg m-2 s-1)
  integer(i4b)    :: scalarSnowfallTemp         = 10 ! temperature of fresh snow (K)
  integer(i4b)    :: scalarSurfaceTemp          = 11 ! temperature of the top layer (K)
  integer(i4b)    :: scalarAlbedo               = 12 ! albedo of the surface, soil or snow (-)
  integer(i4b)    :: scalarSnowDepth            = 13 ! total snow depth (m)
  integer(i4b)    :: scalarSWE                  = 14 ! snow water equivalent (kg m-2)
  integer(i4b)    :: scalarSfcMeltPond          = 15 ! ponded water caused by melt of the "snow without a layer" (kg m-2)
  integer(i4b)    :: scalarAquiferStorage       = 16 ! relative aquifer storage -- above bottom of the soil profile (m)
  integer(i4b)    :: scalarWaterTableDepth      = 17 ! depth of the water table (m)
  integer(i4b)    :: scalarExCoef               = 18 ! turbulent exchange coefficient (-)
  integer(i4b)    :: scalarExSen                = 19 ! exchange factor for sensible heat (W m-2 K-1)
  integer(i4b)    :: scalarExLat                = 20 ! exchange factor for latent heat (W m-2) 
  integer(i4b)    :: scalarSenHeat              = 21 ! sensible heat flux at the surface (W m-2) 
  integer(i4b)    :: scalarLatHeat              = 22 ! latent heat flux at the surface (W m-2)
  integer(i4b)    :: scalarPotentialET          = 23 ! potential evapotranspiration (kg m-2 s-1)
  integer(i4b)    :: scalarTranspireLim         = 24 ! aggregate soil moist & vegetation limit on transpiration (-)
  integer(i4b)    :: scalarTranspireLimAqfr     = 25 ! soil moist & vegetation limit on transpiration in the aquifer (-)
  integer(i4b)    :: scalarMassSolid            = 26 ! sublimation or frost (kg m-2 s-1)
  integer(i4b)    :: scalarMassLiquid           = 27 ! evaporation or dew (kg m-2 s-1)
  integer(i4b)    :: scalarHeatPrecip           = 28 ! sensible heat of precipitation (W m-2)
  integer(i4b)    :: scalarRainPlusMelt         = 29 ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  integer(i4b)    :: scalarSurfaceRunoff        = 30 ! surface runoff (m s-1)
  integer(i4b)    :: scalarInitAquiferTranspire = 31 ! transpiration from the aquifer at the start of the step (m s-1)
  integer(i4b)    :: scalarAquiferTranspire     = 32 ! transpiration from the aquifer (m s-1)
  integer(i4b)    :: scalarInitAquiferRcharge   = 33 ! aquifer recharge rate at the start of the step (m s-1) 
  integer(i4b)    :: scalarAquiferRcharge       = 34 ! aquifer recharge rate (m s-1) 
  integer(i4b)    :: scalarInitAquiferBaseflow  = 35 ! baseflow from the aquifer at the start of the step (m s-1)
  integer(i4b)    :: scalarAquiferBaseflow      = 36 ! baseflow from the aquifer (m s-1)
  integer(i4b)    :: scalarSoilInflux           = 37 ! influx of water at the top of the soil profile (m s-1)
  integer(i4b)    :: scalarSoilDrainage         = 38 ! drainage from the bottom of the soil profile (m s-1) 
  integer(i4b)    :: scalarSoilEjection         = 39 ! total water ejected from each soil layer (m s-1) 
  integer(i4b)    :: scalarSoilWatBalError      = 40 ! error in the total soil water balance (kg m-2)
  integer(i4b)    :: scalarAquiferBalError      = 41 ! error in the aquifer water balance (kg m-2)
  integer(i4b)    :: scalarTotalSoilLiq         = 42 ! total mass of liquid water in the soil (kg m-2)
  integer(i4b)    :: scalarTotalSoilIce         = 43 ! total mass of ice in the soil (kg m-2)
  ! define variables at the mid-point of each layer -- domain geometry
  integer(i4b)    :: mLayerDepth                = 44 ! depth of each layer (m)
  integer(i4b)    :: mLayerHeight               = 45 ! height at the mid-point of each layer (m)
  integer(i4b)    :: mLayerRootDensity          = 46 ! fraction of roots in each soil layer (-)
  ! define variables at the mid-point of each layer -- coupled energy and mass
  integer(i4b)    :: mLayerTemp                 = 47 ! temperature of each layer (K)
  integer(i4b)    :: mLayerVolFracAir           = 48 ! volumetric fraction of air in each layer (-)
  integer(i4b)    :: mLayerVolFracIce           = 49 ! volumetric fraction of ice water in each layer (-)
  integer(i4b)    :: mLayerVolFracLiq           = 50 ! volumetric fraction of liquid water in each layer (-)
  integer(i4b)    :: mLayerVolHtCapBulk         = 51 ! volumetric heat capacity in each layer (J m-3 K-1)
  integer(i4b)    :: mLayerTcrit                = 52 ! critical soil temperature above which all water is unfrozen (K)
  integer(i4b)    :: mLayerdTheta_dTk           = 53 ! derivative in volumetric liquid water content wrt temperature (K-1)
  integer(i4b)    :: mLayerThermalC             = 54 ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
  integer(i4b)    :: mLayerRadCondFlux          = 55 ! temporal derivative in energy from radiative and conductive flux (J m-2 s-1)  
  integer(i4b)    :: mLayerMeltFreeze           = 56 ! melt/freeze in each layer (kg m-3)
  integer(i4b)    :: mLayerSatHydCond           = 57 ! saturated hydraulic conductivity in each layer (m s-1)
  integer(i4b)    :: mLayerMatricHead           = 58 ! matric head of water in the soil (m)
  integer(i4b)    :: mLayerdTheta_dPsi          = 59 ! derivative in the soil water characteristic (m-1)
  integer(i4b)    :: mLayerdPsi_dTheta          = 60 ! derivative in the soil water characteristic (m)
  integer(i4b)    :: mLayerThetaResid           = 61 ! residual volumetric water content in each snow layer (-)
  integer(i4b)    :: mLayerPoreSpace            = 62 ! total pore space in each snow layer (-)
  integer(i4b)    :: mLayerInfilFreeze          = 63 ! volumetric ice content increase by freezing infiltrating flux (-)
  integer(i4b)    :: mLayerTranspireLim         = 64 ! soil moist & veg limit on transpiration for each layer (-)
  integer(i4b)    :: mLayerInitTranspire        = 65 ! transpiration loss from each soil layer at the start of the step (kg m-2 s-1)
  integer(i4b)    :: mLayerTranspire            = 66 ! transpiration loss from each soil layer (kg m-2 s-1)
  integer(i4b)    :: mLayerInitEjectWater       = 67 ! water ejected from each soil layer at the start-of-step (m s-1)
  integer(i4b)    :: mLayerEjectWater           = 68 ! water ejected from each soil layer (m s-1)
  ! define variables at the interface of each layer
  integer(i4b)    :: iLayerHeight               = 69 ! height of the layer interface; top of soil = 0 (m)
  integer(i4b)    :: iLayerThermalC             = 70 ! thermal conductivity at the interface of each layer (W m-1 K-1)
  integer(i4b)    :: iLayerInitNrgFlux          = 71 ! energy flux at layer interfaces at the start of the time step (W m-2) 
  integer(i4b)    :: iLayerNrgFlux              = 72 ! energy flux at layer interfaces at the end of the time step (W m-2) 
  integer(i4b)    :: iLayerSatHydCond           = 73 ! saturated hydraulic conductivity at each layer interface (m s-1)
  integer(i4b)    :: iLayerInitLiqFluxSnow      = 74 ! liquid flux at snow layer interfaces at the start of the time step (m s-1) 
  integer(i4b)    :: iLayerInitLiqFluxSoil      = 75 ! liquid flux at soil layer interfaces at the start of the time step (m s-1) 
  integer(i4b)    :: iLayerLiqFluxSnow          = 76 ! liquid flux at snow layer interfaces at the end of the time step (m s-1)
  integer(i4b)    :: iLayerLiqFluxSoil          = 77 ! liquid flux at soil layer interfaces at the end of the time step (m s-1) 
  ! define derived variables
  integer(i4b)    :: scalarExNeut               = 78 ! exchange coefficient in neutral conditions
  integer(i4b)    :: scalarBprime               = 79 ! stable b parameter in Louis (1979) stability function
  integer(i4b)    :: scalarCparam               = 80 ! c parameter in Louis (1979) stability function
  integer(i4b)    :: scalarVGn_m                = 81 ! van Genuchten "m" parameter (-)
  integer(i4b)    :: scalarKappa                = 82 ! constant in the freezing curve function (m K-1)
  integer(i4b)    :: scalarVolHtCap_air         = 83 ! volumetric heat capacity air         (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_ice         = 84 ! volumetric heat capacity ice         (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_soil        = 85 ! volumetric heat capacity dry soil    (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_water       = 86 ! volumetric heat capacity liquid wat  (J m-3 K-1)
  integer(i4b)    :: scalarLambda_drysoil       = 87 ! thermal conductivity of dry soil     (W m-1)
  integer(i4b)    :: scalarLambda_wetsoil       = 88 ! thermal conductivity of wet soil     (W m-1)
  integer(i4b)    :: scalarVolLatHt_fus         = 89 ! volumetric latent heat of fusion     (J m-3)
  integer(i4b)    :: scalarAquiferRootFrac      = 90 ! fraction of roots below the soil profile (-)
 endtype iLook_mvar

 ! ***********************************************************************************************************
 ! (5) define model indices
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
                                                                        11)
 type(iLook_time),    public,parameter :: iLookTIME     =iLook_time    ( 1, 2, 3, 4, 5)
 type(iLook_force),   public,parameter :: iLookFORCE    =iLook_force   ( 1, 2, 3, 4, 5, 6, 7, 8)
 type(iLook_param),   public,parameter :: iLookPARAM    =iLook_param   ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,&
                                                                        11,12,13,14,15,16,17,18,19,20,&
                                                                        21,22,23,24,25,26,27,28,29,30,&
                                                                        31,32,33,34,35,36,37,38,39,40,&
                                                                        41,42,43,44,45,46,47,48,49,50,&
                                                                        51,52,53,54,55,56,57,58,59,60,&
                                                                        61,62,63,64,65,66,67,68,69,70,&
                                                                        71,72)
 type(iLook_mvar),    public,parameter :: iLookMVAR     =ilook_mvar    ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10,&
                                                                        11,12,13,14,15,16,17,18,19,20,&
                                                                        21,22,23,24,25,26,27,28,29,30,&
                                                                        31,32,33,34,35,36,37,38,39,40,&
                                                                        41,42,43,44,45,46,47,48,49,50,&
                                                                        51,52,53,54,55,56,57,58,59,60,&
                                                                        61,62,63,64,65,66,67,68,69,70,&
                                                                        71,72,73,74,75,76,77,78,79,80,&
                                                                        81,82,83,84,85,86,87,88,89,90)
 type(iLook_index),   public,parameter :: iLookINDEX    =ilook_index   ( 1, 2, 3, 4, 5, 6, 7, 8, 9,10)
 ! define maximum number of variables of each type
 integer(i4b),parameter,public :: maxvarDecisions= 11
 integer(i4b),parameter,public :: maxvarTime     = 5
 integer(i4b),parameter,public :: maxvarForc     = 8
 integer(i4b),parameter,public :: maxvarMpar     = 72
 integer(i4b),parameter,public :: maxvarMvar     = 90
 integer(i4b),parameter,public :: maxvarIndx     = 10
 ! ***********************************************************************************************************
END MODULE var_lookup

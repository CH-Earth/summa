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
  ! simulation options
  integer(i4b)    :: simulStart       = 1  ! simulation start time
  integer(i4b)    :: simulFinsh       = 2  ! simulation end time
  ! Noah-MP options
  integer(i4b)    :: soilCatTbl       = 3  ! soil-category dateset
  integer(i4b)    :: vegeParTbl       = 4  ! vegetation category dataset
  integer(i4b)    :: soilStress       = 5  ! choice of function for the soil moisture control on stomatal resistance
  integer(i4b)    :: stomResist       = 6  ! choice of function for stomatal resistance
  ! FUSE options
  integer(i4b)    :: num_method       = 7  ! choice of numerical method
  integer(i4b)    :: fDerivMeth       = 8  ! method used to calculate flux derivatives
  integer(i4b)    :: LAI_method       = 9  ! method used to determine LAI and SAI
  integer(i4b)    :: f_Richards       = 10  ! form of richards' equation
  integer(i4b)    :: groundwatr       = 11 ! choice of groundwater parameterization
  integer(i4b)    :: hc_profile       = 12 ! choice of hydraulic conductivity profile
  integer(i4b)    :: bcUpprTdyn       = 13 ! type of upper boundary condition for thermodynamics 
  integer(i4b)    :: bcLowrTdyn       = 14 ! type of lower boundary condition for thermodynamics
  integer(i4b)    :: bcUpprSoiH       = 15 ! type of upper boundary condition for soil hydrology
  integer(i4b)    :: bcLowrSoiH       = 16 ! type of lower boundary condition for soil hydrology
  integer(i4b)    :: veg_traits       = 17 ! choice of parameterization for vegetation roughness length and displacement height
  integer(i4b)    :: canopyEmis       = 18 ! choice of parameterization for canopy emissivity 
  integer(i4b)    :: astability       = 19 ! choice of stability function
  integer(i4b)    :: alb_method       = 20 ! choice of albedo representation
  integer(i4b)    :: snowLayers       = 21 ! choice of method to combine and sub-divide snow layers
  integer(i4b)    :: compaction       = 22 ! choice of compaction routine
  integer(i4b)    :: thermlcond       = 23 ! choice of thermal conductivity representation
  integer(i4b)    :: spatial_gw       = 24 ! choice of method for spatial representation of groundwater
  integer(i4b)    :: subRouting       = 25 ! choice of method for sub-grid routing
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
  integer(i4b)    :: SWRadAtm         = 6  ! downwelling shortwave radiaiton (W m-2)
  integer(i4b)    :: LWRadAtm         = 7  ! downwelling longwave radiation  (W m-2)
  integer(i4b)    :: airpres          = 8  ! pressure                        (Pa)
 endtype iLook_force
 ! ***********************************************************************************************************
 ! (3) define local attributes
 ! ***********************************************************************************************************
 type, public  ::  iLook_attr
  integer(i4b)    :: latitude         = 1  ! latitude (degrees north)
  integer(i4b)    :: longitude        = 2  ! longitude (degrees east)
  integer(i4b)    :: elevation        = 3  ! elevation (m)
  integer(i4b)    :: tan_slope        = 4  ! tan water table slope, taken as tan local ground surface slope (-)
  integer(i4b)    :: contourLength    = 5  ! length of contour at downslope edge of HRU (m)
  integer(i4b)    :: HRUarea          = 6  ! area of each HRU  (m2)
  integer(i4b)    :: mHeight          = 7  ! measurement height above bare ground (m)
 end type iLook_attr
 ! ***********************************************************************************************************
 ! (4) define local classification of veg, soil, etc. 
 ! ***********************************************************************************************************
 type, public  ::  iLook_type
  integer(i4b)    :: hruIndex         = 1  ! index defining hydrologic response unit (-)
  integer(i4b)    :: vegTypeIndex     = 2  ! index defining vegetation type (-)
  integer(i4b)    :: soilTypeIndex    = 3  ! index defining soil type (-)
  integer(i4b)    :: slopeTypeIndex   = 4  ! index defining slope (-)
  integer(i4b)    :: downHRUindex     = 5  ! index of downslope HRU (0 = basin outlet)
 end type iLook_type
 ! ***********************************************************************************************************
 ! (5) define model parameters
 ! ***********************************************************************************************************
 type, public  ::  iLook_param
  ! boundary conditions
  integer(i4b)    :: upperBoundHead       = 1   ! matric head of the upper boundary (m)
  integer(i4b)    :: lowerBoundHead       = 2   ! matric head of the lower boundary (m)
  integer(i4b)    :: upperBoundTheta      = 3   ! volumetric liquid water content of the upper boundary (-)
  integer(i4b)    :: lowerBoundTheta      = 4   ! volumetric liquid water content of the lower boundary (-)
  integer(i4b)    :: upperBoundTemp       = 5   ! temperature of the upper boundary (K)
  integer(i4b)    :: lowerBoundTemp       = 6   ! temperature of the lower boundary (K)
  ! precipitation partitioning
  integer(i4b)    :: tempCritRain         = 7   ! critical temperature where precipitation is rain (K)
  integer(i4b)    :: tempRangeTimestep    = 8   ! temperature range over the time step (K)
  integer(i4b)    :: frozenPrecipMultip   = 9   ! frozen precipitation multiplier (-)
  ! freezing curve for snow
  integer(i4b)    :: snowfrz_scale        = 10  ! scaling parameter for the freezing curve for snow (K-1)
  ! snow albedo
  integer(i4b)    :: snw_crit             = 11  ! critical mass necessary for albedo refreshment (kg m-2)
  integer(i4b)    :: alb_fresh            = 12  ! fresh snow albedo (-)
  integer(i4b)    :: alb_dry              = 13  ! minimum snow albedo during winter (-)
  integer(i4b)    :: alb_wet              = 14  ! minimum snow albedo during spring (-)
  integer(i4b)    :: alb_decay            = 15  ! temporal decay factor for snow albedo (s-1)
  integer(i4b)    :: alb_scale            = 16  ! albedo scaling factor (s)
  integer(i4b)    :: soot_load            = 17  ! temporal decay in snow albedo associated with the soot load (days-1)
  ! radiation transfer within snow
  integer(i4b)    :: radExt_snow          = 18  ! extinction coefficient for radiation penetration into the snowpack (m-1)
  integer(i4b)    :: Frad_direct          = 19  ! fraction of direct solar radiation (-)
  integer(i4b)    :: Frad_vis             = 20  ! fraction of radiation in the visible part of the spectrum (-)
  ! new snow density
  integer(i4b)    :: newSnowDenMin        = 21  ! minimum new snow density (kg m-3)  
  integer(i4b)    :: newSnowDenMult       = 22  ! multiplier for new snow density (kg m-3)
  integer(i4b)    :: newSnowDenScal       = 23  ! scaling factor for new snow density (K)
  ! snow compaction
  integer(i4b)    :: densScalGrowth       = 24  ! density scaling factor for grain growth (kg-1 m3)
  integer(i4b)    :: tempScalGrowth       = 25  ! temperature scaling factor for grain growth (K-1)
  integer(i4b)    :: grainGrowthRate      = 26  ! rate of grain growth (s-1)
  integer(i4b)    :: densScalOvrbdn       = 27  ! density scaling factor for overburden pressure (kg-1 m3)
  integer(i4b)    :: tempScalOvrbdn       = 28  ! temperature scaling factor for overburden pressure (K-1)
  integer(i4b)    :: base_visc            = 29  ! viscosity coefficient at T=T_frz and snow density=0  (kg s m-2)
  ! water flow within snow
  integer(i4b)    :: Fcapil               = 30  ! capillary retention as a fraction of the total pore volume (-)
  integer(i4b)    :: k_snow               = 31  ! hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
  integer(i4b)    :: mw_exp               = 32  ! exponent for meltwater flow (-)
  ! turbulent heat fluxes
  integer(i4b)    :: z0Snow               = 33  ! roughness length of snow (m)
  integer(i4b)    :: z0Soil               = 34  ! roughness length of bare soil below the canopy (m)
  integer(i4b)    :: z0Canopy             = 35  ! roughness length of the canopy (m)
  integer(i4b)    :: zpdFraction          = 36  ! zero plane displacement / canopy height (-)
  integer(i4b)    :: critRichNumber       = 37  ! critical value for the bulk Richardson number (-)
  integer(i4b)    :: Louis79_bparam       = 38  ! parameter in Louis (1979) stability function (-)
  integer(i4b)    :: Louis79_cStar        = 39  ! parameter in Louis (1979) stability function (-)
  integer(i4b)    :: Mahrt87_eScale       = 40  ! exponential scaling factor in the Mahrt (1987) stability function (-)
  integer(i4b)    :: leafExchangeCoeff    = 41  ! turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
  integer(i4b)    :: windReductionParam   = 42  ! canopy wind reduction parameter (-)
  ! vegetation properties
  integer(i4b)    :: winterSAI            = 43  ! stem area index prior to the start of the growing season (m2 m-2)
  integer(i4b)    :: summerLAI            = 44  ! maximum leaf area index at the peak of the growing season (m2 m-2)
  integer(i4b)    :: rootingDepth         = 45  ! rooting depth (m)
  integer(i4b)    :: rootDistExp          = 46  ! exponent controlling the vertical distribution of root density (-)
  integer(i4b)    :: plantWiltPsi         = 47  ! matric head at wilting point (m)
  integer(i4b)    :: soilStressParam      = 48  ! parameter in the exponential soil stress function
  integer(i4b)    :: critSoilWilting      = 49  ! critical vol. liq. water content when plants are wilting (-) 
  integer(i4b)    :: critSoilTranspire    = 50  ! critical vol. liq. water content when transpiration is limited (-)
  integer(i4b)    :: critAquiferTranspire = 51  ! critical aquifer storage value when transpiration is limited (m)
  integer(i4b)    :: leafDimension        = 52  ! characteristic leaf dimension (m)
  integer(i4b)    :: heightCanopyTop      = 53  ! height of top of the vegetation canopy above ground surface (m)
  integer(i4b)    :: heightCanopyBottom   = 54  ! height of bottom of the vegetation canopy above ground surface (m)
  integer(i4b)    :: maxCanopyIce         = 55  ! maximum mass of ice on the vegetation canopy per unit lai+sai (kg m-2)
  integer(i4b)    :: maxCanopyLiquid      = 56  ! maximum mass of liq water on vegetation canopy per unit lai+sai (kg m-2)
  integer(i4b)    :: specificHeatVeg      = 57  ! specific heat of vegetation (J kg-1 K-1)
  integer(i4b)    :: maxMassVegetation    = 58  ! maximum mass of vegetation (full foliage) (kg m-2)
  ! soil properties
  integer(i4b)    :: soil_dens_intr       = 59  ! intrinsic soil density (kg m-3)
  integer(i4b)    :: frac_sand            = 60  ! fraction of sand (-)
  integer(i4b)    :: frac_silt            = 61  ! fraction of silt (-)
  integer(i4b)    :: frac_clay            = 62  ! fraction of clay (-)
  integer(i4b)    :: theta_sat            = 63  ! porosity (-)
  integer(i4b)    :: theta_res            = 64  ! volumetric residual water content (-)
  integer(i4b)    :: vGn_alpha            = 65  ! van Genuchten "alpha" parameter (m-1)
  integer(i4b)    :: vGn_n                = 66  ! van Genuchten "n" parameter (-)
  integer(i4b)    :: k_soil               = 67  ! hydraulic conductivity of soil (m s-1) 
  integer(i4b)    :: kAnisotropic         = 68  ! anisotropy factor for lateral hydraulic conductivity (-)
  integer(i4b)    :: zScale_TOPMODEL      = 69  ! TOPMODEL scaling factor used in lower boundary condition for soil (m)
  integer(i4b)    :: compactedDepth       = 70  ! depth where k_soil reaches the compacted value given by CH78 (m)
  integer(i4b)    :: aquiferScaleFactor   = 71  ! scaling factor for aquifer storage in the big bucket (m)
  integer(i4b)    :: aquiferBaseflowExp   = 72  ! baseflow exponent (-)
  integer(i4b)    :: bpar_VIC             = 73  ! b-parameter in the VIC surface runoff parameterization (-)
  integer(i4b)    :: specificYield        = 74  ! specific yield (-)
  integer(i4b)    :: specificStorage      = 75  ! specific storage coefficient (m-1)
  integer(i4b)    :: f_impede             = 76  ! ice impedence factor (-)
  ! algorithmic control parameters
  integer(i4b)    :: minwind              = 77  ! minimum wind speed (m s-1)
  integer(i4b)    :: minstep              = 78  ! minimum length of the time step
  integer(i4b)    :: maxstep              = 79  ! maximum length of the time step
  integer(i4b)    :: wimplicit            = 80  ! weight assigned to the start-of-step fluxes
  integer(i4b)    :: maxiter              = 81  ! maximum number of iteration 
  integer(i4b)    :: relConvTol_liquid    = 82  ! relative convergence tolerance for vol frac liq water (-)
  integer(i4b)    :: absConvTol_liquid    = 83  ! absolute convergence tolerance for vol frac liq water (-)
  integer(i4b)    :: relConvTol_matric    = 84  ! relative convergence tolerance for matric head (-)
  integer(i4b)    :: absConvTol_matric    = 85  ! absolute convergence tolerance for matric head (m)
  integer(i4b)    :: relConvTol_energy    = 86  ! relative convergence tolerance for energy (-)
  integer(i4b)    :: absConvTol_energy    = 87  ! absolute convergence tolerance for energy (J m-3)
  integer(i4b)    :: relConvTol_aquifr    = 88  ! relative convergence tolerance for aquifer storage (-)
  integer(i4b)    :: absConvTol_aquifr    = 89  ! absolute convergence tolerance for aquifer storage (J m-3)
  integer(i4b)    :: zmin                 = 90  ! minimum layer depth (m)
  integer(i4b)    :: zmax                 = 91  ! maximum layer depth (m)
  integer(i4b)    :: zminLayer1           = 92  ! minimum layer depth for the 1st (top) layer (m)
  integer(i4b)    :: zminLayer2           = 93  ! minimum layer depth for the 2nd layer (m) 
  integer(i4b)    :: zminLayer3           = 94  ! minimum layer depth for the 3rd layer (m) 
  integer(i4b)    :: zminLayer4           = 95  ! minimum layer depth for the 4th layer (m) 
  integer(i4b)    :: zminLayer5           = 96  ! minimum layer depth for the 5th (bottom) layer (m) 
  integer(i4b)    :: zmaxLayer1_lower     = 97  ! maximum layer depth for the 1st (top) layer when only 1 layer (m) 
  integer(i4b)    :: zmaxLayer2_lower     = 98  ! maximum layer depth for the 2nd layer when only 2 layers (m) 
  integer(i4b)    :: zmaxLayer3_lower     = 99  ! maximum layer depth for the 3rd layer when only 3 layers (m) 
  integer(i4b)    :: zmaxLayer4_lower     = 100 ! maximum layer depth for the 4th layer when only 4 layers (m) 
  integer(i4b)    :: zmaxLayer1_upper     = 101 ! maximum layer depth for the 1st (top) layer when > 1 layer (m) 
  integer(i4b)    :: zmaxLayer2_upper     = 102 ! maximum layer depth for the 2nd layer when > 2 layers (m) 
  integer(i4b)    :: zmaxLayer3_upper     = 103 ! maximum layer depth for the 3rd layer when > 3 layers (m) 
  integer(i4b)    :: zmaxLayer4_upper     = 104 ! maximum layer depth for the 4th layer when > 4 layers (m) 
 endtype ilook_param
 ! ***********************************************************************************************************
 ! (6) define model variables
 ! ***********************************************************************************************************
 type, public :: iLook_mvar
  ! define timestep-average fluxes for a few key variables
  integer(i4b)    :: averageSnowSublimation          = 1   ! snow sublimation/frost - below canopy or non-vegetated (kg m-2 s-1)
  integer(i4b)    :: averageGroundEvaporation        = 2   ! ground evaporation/condensation - below canopy or non-vegetated (kg m-2 s-1)
  integer(i4b)    :: averageRainPlusMelt             = 3   ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  integer(i4b)    :: averageSurfaceRunoff            = 4   ! surface runoff (m s-1)
  integer(i4b)    :: averageSoilInflux               = 5   ! influx of water at the top of the soil profile (m s-1)
  integer(i4b)    :: averageSoilBaseflow             = 6   ! total baseflow from throughout the soil profile (m s-1)
  integer(i4b)    :: averageSoilDrainage             = 7   ! drainage from the bottom of the soil profile (m s-1)
  integer(i4b)    :: averageSoilEjection             = 8   ! ejected water from the soil matrix (m s-1)
  integer(i4b)    :: averageAquiferRecharge          = 9   ! recharge to the aquifer (m s-1)
  integer(i4b)    :: averageAquiferBaseflow          = 10  ! baseflow from the aquifer (m s-1)
  integer(i4b)    :: averageAquiferTranspire         = 11  ! transpiration from the aquifer (m s-1)
  ! define scalar variables -- forcing
  integer(i4b)    :: scalarCosZenith                 = 12  ! cosine of the solar zenith angle (0-1)
  integer(i4b)    :: spectralIncomingDirect          = 13  ! incoming direct solar radiation in each wave band (W m-2)
  integer(i4b)    :: spectralIncomingDiffuse         = 14  ! incoming diffuse solar radiation in each wave band (W m-2)
  integer(i4b)    :: scalarVPair                     = 15  ! vapor pressure of the air above the vegetation canopy (Pa)
  integer(i4b)    :: scalarTwetbulb                  = 16  ! wet bulb temperature (K)
  integer(i4b)    :: scalarRainfall                  = 17  ! computed rainfall rate (kg m-2 s-1)
  integer(i4b)    :: scalarSnowfall                  = 18  ! computed snowfall rate (kg m-2 s-1)
  integer(i4b)    :: scalarSnowfallTemp              = 19  ! temperature of fresh snow (K)
  integer(i4b)    :: scalarO2air                     = 20  ! atmospheric o2 concentration (Pa)
  integer(i4b)    :: scalarCO2air                    = 21  ! atmospheric co2 concentration (Pa)
  ! define scalar variables -- state variables
  integer(i4b)    :: scalarCanopyIce                 = 22  ! mass of ice on the vegetation canopy (kg m-2)
  integer(i4b)    :: scalarCanopyLiq                 = 23  ! mass of liquid water on the vegetation canopy (kg m-2)
  integer(i4b)    :: scalarCanopyTemp                = 24  ! temperature of the vegetation canopy (K)
  integer(i4b)    :: scalarAlbedo                    = 25  ! snow albedo (-)
  integer(i4b)    :: scalarSnowAge                   = 26  ! non-dimensional snow age (-)
  integer(i4b)    :: scalarSnowDepth                 = 27  ! total snow depth (m)
  integer(i4b)    :: scalarSWE                       = 28  ! snow water equivalent (kg m-2)
  integer(i4b)    :: scalarSfcMeltPond               = 29  ! ponded water caused by melt of the "snow without a layer" (kg m-2)
  integer(i4b)    :: scalarAquiferStorage            = 30  ! relative aquifer storage -- above bottom of the soil profile (m)
  integer(i4b)    :: scalarSurfaceTemp               = 31  ! surface temperature (K)
  ! define NOAH-MP vegetation variables -- general
  integer(i4b)    :: scalarGreenVegFraction          = 32  ! green vegetation fraction used to compute LAI (-)
  integer(i4b)    :: scalarBulkVolHeatCapVeg         = 33  ! bulk volumetric heat capacity of vegetation (J m-3 K-1)
  integer(i4b)    :: scalarRootZoneTemp              = 34  ! average temperature of the root zone (K)
  integer(i4b)    :: scalarLAI                       = 35  ! one-sided leaf area index (m2 m-2)
  integer(i4b)    :: scalarSAI                       = 36  ! one-sided stem area index (m2 m-2)
  integer(i4b)    :: scalarExposedLAI                = 37  ! exposed leaf area index after burial by snow (m2 m-2)
  integer(i4b)    :: scalarExposedSAI                = 38  ! exposed stem area index after burial by snow(m2 m-2)
  integer(i4b)    :: scalarGrowingSeasonIndex        = 39  ! growing season index (0=off, 1=on)
  integer(i4b)    :: scalarVP_CanopyAir              = 40  ! vapor pressure of the canopy air space (Pa)
  integer(i4b)    :: scalarTemp_CanopyAir            = 41  ! temperature of the canopy air space (Pa)
  ! define NOAH-MP vegetation variables -- shortwave radiation
  integer(i4b)    :: scalarCanopySunlitFraction      = 42  ! sunlit fraction of canopy (-)
  integer(i4b)    :: scalarCanopySunlitLAI           = 43  ! sunlit leaf area (-)
  integer(i4b)    :: scalarCanopyShadedLAI           = 44  ! shaded leaf area (-)
  integer(i4b)    :: scalarCanopySunlitPAR           = 45  ! average absorbed par for sunlit leaves (w m-2)
  integer(i4b)    :: scalarCanopyShadedPAR           = 46  ! average absorbed par for shaded leaves (w m-2)
  integer(i4b)    :: scalarCanopyAbsorbedSolar       = 47  ! solar radiation absorbed by canopy (W m-2)
  integer(i4b)    :: scalarGroundAbsorbedSolar       = 48  ! solar radiation absorbed by ground (W m-2)
  integer(i4b)    :: scalarTotalReflectedSolar       = 49  ! total reflected solar radiation (W m-2)
  integer(i4b)    :: scalarTotalAbsorbedSolar        = 50  ! total absorbed solar radiation (W m-2)
  integer(i4b)    :: scalarCanopyReflectedSolar      = 51  ! solar radiation reflected from the canopy (W m-2)
  integer(i4b)    :: scalarGroundReflectedSolar      = 52  ! solar radiation reflected from the ground (W m-2) 
  integer(i4b)    :: scalarBetweenCanopyGapFraction  = 53  ! between canopy gap fraction for beam (-)
  integer(i4b)    :: scalarWithinCanopyGapFraction   = 54  ! within canopy gap fraction for beam (-)
  ! define NOAH-MP vegetation variables -- longwave radiation
  integer(i4b)    :: scalarCanopyEmissivity          = 55  ! effective canopy emissivity (-)
  integer(i4b)    :: scalarLWRadCanopy               = 56  ! longwave radiation emitted from the canopy (W m-2)
  integer(i4b)    :: scalarLWRadGround               = 57  ! longwave radiation emitted at the ground surface  (W m-2)
  integer(i4b)    :: scalarLWRadUbound2Canopy        = 58  ! downward atmospheric longwave radiation absorbed by the canopy (W m-2)
  integer(i4b)    :: scalarLWRadUbound2Ground        = 59  ! downward atmospheric longwave radiation absorbed by the ground (W m-2)
  integer(i4b)    :: scalarLWRadUbound2Ubound        = 60  ! atmospheric radiation refl by ground + lost thru upper boundary (W m-2)
  integer(i4b)    :: scalarLWRadCanopy2Ubound        = 61  ! longwave radiation emitted from canopy lost thru upper boundary (W m-2)
  integer(i4b)    :: scalarLWRadCanopy2Ground        = 62  ! longwave radiation emitted from canopy absorbed by the ground (W m-2)
  integer(i4b)    :: scalarLWRadCanopy2Canopy        = 63  ! canopy longwave reflected from ground and absorbed by the canopy (W m-2)
  integer(i4b)    :: scalarLWRadGround2Ubound        = 64  ! longwave radiation emitted from ground lost thru upper boundary (W m-2)
  integer(i4b)    :: scalarLWRadGround2Canopy        = 65  ! longwave radiation emitted from ground and absorbed by the canopy (W m-2)
  integer(i4b)    :: scalarLWNetCanopy               = 66  ! net longwave radiation at the canopy (W m-2)
  integer(i4b)    :: scalarLWNetGround               = 67  ! net longwave radiation at the ground surface (W m-2)
  integer(i4b)    :: scalarLWNetUbound               = 68  ! net longwave radiation at the upper atmospheric boundary (W m-2)
  ! define NOAH-MP vegetation variables -- turbulent heat transfer
  integer(i4b)    :: scalarLatHeatSubVapCanopy       = 69  ! latent heat of sublimation/vaporization used for veg canopy (J kg-1)
  integer(i4b)    :: scalarLatHeatSubVapGround       = 70  ! latent heat of sublimation/vaporization used for ground surface (J kg-1)
  integer(i4b)    :: scalarSatVP_CanopyTemp          = 71  ! saturation vapor pressure at the temperature of vegetation canopy (Pa)
  integer(i4b)    :: scalarSatVP_GroundTemp          = 72  ! saturation vapor pressure at the temperature of the ground (Pa)
  integer(i4b)    :: scalarZ0Canopy                  = 73  ! roughness length of the canopy (m)
  integer(i4b)    :: scalarWindReductionFactor       = 74  ! canopy wind reduction factor (-)
  integer(i4b)    :: scalarZeroPlaneDisplacement     = 75  ! zero plane displacement (m) 
  integer(i4b)    :: scalarRiBulkCanopy              = 76  ! bulk Richardson number for the canopy (-)
  integer(i4b)    :: scalarRiBulkGround              = 77  ! bulk Richardson number for the ground surface (-)
  integer(i4b)    :: scalarCanopyStabilityCorrection = 78  ! stability correction for the canopy (-)
  integer(i4b)    :: scalarGroundStabilityCorrection = 79  ! stability correction for the ground surface (-)
  integer(i4b)    :: scalarEddyDiffusCanopyTop       = 80  ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
  integer(i4b)    :: scalarFrictionVelocity          = 81  ! friction velocity - canopy momentum sink (m s-1)
  integer(i4b)    :: scalarWindspdCanopyTop          = 82  ! windspeed at the top of the canopy (m s-1)
  integer(i4b)    :: scalarWindspdCanopyBottom       = 83  ! windspeed at the height of the bottom of the canopy (m s-1)
  integer(i4b)    :: scalarGroundResistance          = 84  ! below canopy aerodynamic resistance (s m-1) 
  integer(i4b)    :: scalarCanopyResistance          = 85  ! above canopy aerodynamic resistance (s m-1)
  integer(i4b)    :: scalarLeafResistance            = 86  ! mean leaf boundary layer resistance per unit leaf area (s m-1)
  integer(i4b)    :: scalarSoilResistance            = 87  ! soil surface resistance (s m-1)
  integer(i4b)    :: scalarSoilRelHumidity           = 88  ! relative humidity in the soil pores in the upper-most soil layer (-)
  integer(i4b)    :: scalarSenHeatCanopy             = 89  ! sensible heat from the canopy to the canopy air space (W m-2)
  integer(i4b)    :: scalarSenHeatGround             = 90  ! sensible heat from the ground (below canopy or non-vegetated) (W m-2)
  integer(i4b)    :: scalarLatHeatCanopyEvap         = 91  ! evaporation latent heat from the canopy to the canopy air space (W m-2) 
  integer(i4b)    :: scalarLatHeatCanopyTrans        = 92  ! transpiration latent heat from the canopy to the canopy air space (W m-2)
  integer(i4b)    :: scalarLatHeatGround             = 93  ! latent heat from the ground (below canopy or non-vegetated) (W m-2)
  integer(i4b)    :: scalarAdvectiveHeatFlux         = 94  ! heat advected to the surface with precipitation (W m-2)
  integer(i4b)    :: scalarCanopyTranspiration       = 95  ! canopy transpiration (kg m-2 s-1)
  integer(i4b)    :: scalarCanopyEvaporation         = 96  ! canopy evaporation/condensation (kg m-2 s-1)
  integer(i4b)    :: scalarCanopySublimation         = 97  ! canopy sublimation/frost (kg m-2 s-1)
  integer(i4b)    :: scalarGroundEvaporation         = 98  ! ground evaporation/condensation (below canopy or non-vegetated) (kg m-2 s-1)
  integer(i4b)    :: scalarSnowSublimation           = 99  ! snow sublimation/frost (below canopy or non-vegetated) (kg m-2 s-1)
  ! define NOAH-MP vegetation variables -- transpiration
  integer(i4b)    :: scalarTranspireLim              = 100 ! aggregate soil moisture + aquifer storage limit on transpiration (-)
  integer(i4b)    :: scalarTranspireLimAqfr          = 101 ! aquifer storage limit on transpiration (-)
  integer(i4b)    :: scalarFoliageNitrogenFactor     = 102 ! foliage nitrogen concentration, 1=saturated (-)
  integer(i4b)    :: scalarStomResistSunlit          = 103 ! stomatal resistance for sunlit leaves (s m-1)
  integer(i4b)    :: scalarStomResistShaded          = 104 ! stomatal resistance for shaded leaves (s m-1)
  integer(i4b)    :: scalarPhotosynthesisSunlit      = 105 ! sunlit photosynthesis (umolco2 m-2 s-1)
  integer(i4b)    :: scalarPhotosynthesisShaded      = 106 ! shaded photosynthesis (umolco2 m-2 s-1)
  ! define NOAH-MP vegetation variables --hydrology
  integer(i4b)    :: scalarCanopyWetFraction         = 107 ! fraction of canopy that is wet
  integer(i4b)    :: scalarGroundSnowFraction        = 108 ! fraction of ground that is covered with snow (-)
  ! define scalar variables -- soil and aquifer fluxes
  integer(i4b)    :: scalarRainPlusMelt              = 109 ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  integer(i4b)    :: scalarSurfaceRunoff             = 110 ! surface runoff (m s-1)
  integer(i4b)    :: scalarInitAquiferRecharge       = 111 ! recharge to the aquifer at the start of the step (m s-1)
  integer(i4b)    :: scalarAquiferRecharge           = 112 ! recharge to the aquifer (m s-1)
  integer(i4b)    :: scalarInitAquiferTranspire      = 113 ! transpiration from the aquifer at the start of the step (m s-1)
  integer(i4b)    :: scalarAquiferTranspire          = 114 ! transpiration from the aquifer (m s-1)
  integer(i4b)    :: scalarInitAquiferBaseflow       = 115 ! baseflow from the aquifer at the start of the step (m s-1)
  integer(i4b)    :: scalarAquiferBaseflow           = 116 ! baseflow from the aquifer (m s-1)
  ! scalar variables -- sub-step average fluxes for the soil zone
  integer(i4b)    :: scalarSoilInflux                = 117 ! sub-step average: influx of water at the top of the soil profile (m s-1)
  integer(i4b)    :: scalarSoilBaseflow              = 118 ! sub-step average: total baseflow from throughout the soil profile (m s-1)
  integer(i4b)    :: scalarSoilDrainage              = 119 ! sub-step average: drainage from the bottom of the soil profile (m s-1) 
  integer(i4b)    :: scalarSoilEjection              = 120 ! sub-step average: total water ejected from all soil layers (m s-1) 
  integer(i4b)    :: scalarSoilTranspiration         = 121 ! sub-step average: total transpiration from the soil (m s-1)
  ! define scalar variables -- mass balance check
  integer(i4b)    :: scalarSoilWatBalError           = 122 ! error in the total soil water balance (kg m-2)
  integer(i4b)    :: scalarAquiferBalError           = 123 ! error in the aquifer water balance (kg m-2)
  integer(i4b)    :: scalarTotalSoilLiq              = 124 ! total mass of liquid water in the soil (kg m-2)
  integer(i4b)    :: scalarTotalSoilIce              = 125 ! total mass of ice in the soil (kg m-2)
  ! define variables at the mid-point of each layer -- domain geometry
  integer(i4b)    :: mLayerDepth                     = 126 ! depth of each layer (m)
  integer(i4b)    :: mLayerHeight                    = 127 ! height at the mid-point of each layer (m)
  integer(i4b)    :: mLayerRootDensity               = 128 ! fraction of roots in each soil layer (-)
  ! define variables at the mid-point of each layer -- coupled energy and mass
  integer(i4b)    :: mLayerTemp                      = 129 ! temperature of each layer (K)
  integer(i4b)    :: mLayerVolFracAir                = 130 ! volumetric fraction of air in each layer (-)
  integer(i4b)    :: mLayerVolFracIce                = 131 ! volumetric fraction of ice water in each layer (-)
  integer(i4b)    :: mLayerVolFracLiq                = 132 ! volumetric fraction of liquid water in each layer (-)
  integer(i4b)    :: mLayerVolHtCapBulk              = 133 ! volumetric heat capacity in each layer (J m-3 K-1)
  integer(i4b)    :: mLayerTcrit                     = 134 ! critical soil temperature above which all water is unfrozen (K)
  integer(i4b)    :: mLayerdTheta_dTk                = 135 ! derivative in volumetric liquid water content wrt temperature (K-1)
  integer(i4b)    :: mLayerThermalC                  = 136 ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
  integer(i4b)    :: mLayerRadCondFlux               = 137 ! temporal derivative in energy from radiative and conductive flux (J m-2 s-1)  
  integer(i4b)    :: mLayerMeltFreeze                = 138 ! rate of ice content change from melt/freeze in each layer (kg m-3 s-1)
  integer(i4b)    :: mLayerInfilFreeze               = 139 ! rate of ice content change by freezing infiltrating flux (kg m-3 s-1)
  integer(i4b)    :: mLayerSatHydCond                = 140 ! saturated hydraulic conductivity in each layer (m s-1)
  integer(i4b)    :: mLayerMatricHead                = 141 ! matric head of water in the soil (m)
  integer(i4b)    :: mLayerdTheta_dPsi               = 142 ! derivative in the soil water characteristic (m-1)
  integer(i4b)    :: mLayerdPsi_dTheta               = 143 ! derivative in the soil water characteristic (m)
  integer(i4b)    :: mLayerThetaResid                = 144 ! residual volumetric water content in each snow layer (-)
  integer(i4b)    :: mLayerPoreSpace                 = 145 ! total pore space in each snow layer (-)
  integer(i4b)    :: mLayerTranspireLim              = 146 ! soil moist & veg limit on transpiration for each layer (-)
  integer(i4b)    :: mLayerInitTranspire             = 147 ! transpiration loss from each soil layer at the start of the step (kg m-2 s-1)
  integer(i4b)    :: mLayerTranspire                 = 148 ! transpiration loss from each soil layer (kg m-2 s-1)
  integer(i4b)    :: mLayerInitEjectWater            = 149 ! water ejected from each soil layer at the start-of-step (m s-1)
  integer(i4b)    :: mLayerEjectWater                = 150 ! water ejected from each soil layer (m s-1)
  integer(i4b)    :: mLayerInitBaseflow              = 151 ! baseflow from each soil layer at the start of the time step (m s-1)
  integer(i4b)    :: mLayerBaseflow                  = 152 ! baseflow from each soil layer (m s-1)
  integer(i4b)    :: mLayerColumnInflow              = 153 ! total inflow to each layer in a given soil column (m3 s-1)
  integer(i4b)    :: mLayerColumnOutflow             = 154 ! total outflow from each layer in a given soil column (m3 s-1)
  ! define variables at the interface of each layer
  integer(i4b)    :: iLayerHeight                    = 155 ! height of the layer interface; top of soil = 0 (m)
  integer(i4b)    :: iLayerThermalC                  = 156 ! thermal conductivity at the interface of each layer (W m-1 K-1)
  integer(i4b)    :: iLayerConductiveFlux            = 157 ! conductive energy flux at layer interfaces at end of time step (W m-2)
  integer(i4b)    :: iLayerAdvectiveFlux             = 158 ! advective energy flux at layer interfaces at end of time step (W m-2)
  integer(i4b)    :: iLayerInitNrgFlux               = 159 ! energy flux at layer interfaces at the start of the time step (W m-2) 
  integer(i4b)    :: iLayerNrgFlux                   = 160 ! energy flux at layer interfaces at the end of the time step (W m-2) 
  integer(i4b)    :: iLayerSatHydCond                = 161 ! saturated hydraulic conductivity at each layer interface (m s-1)
  integer(i4b)    :: iLayerInitLiqFluxSnow           = 162 ! liquid flux at snow layer interfaces at the start of the time step (m s-1) 
  integer(i4b)    :: iLayerInitLiqFluxSoil           = 163 ! liquid flux at soil layer interfaces at the start of the time step (m s-1) 
  integer(i4b)    :: iLayerLiqFluxSnow               = 164 ! liquid flux at snow layer interfaces at the end of the time step (m s-1)
  integer(i4b)    :: iLayerLiqFluxSoil               = 165 ! liquid flux at soil layer interfaces at the end of the time step (m s-1) 
  ! define derived variables
  integer(i4b)    :: scalarVGn_m                     = 166 ! van Genuchten "m" parameter (-)
  integer(i4b)    :: scalarKappa                     = 167 ! constant in the freezing curve function (m K-1)
  integer(i4b)    :: scalarVolHtCap_air              = 168 ! volumetric heat capacity air         (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_ice              = 169 ! volumetric heat capacity ice         (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_soil             = 170 ! volumetric heat capacity dry soil    (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_water            = 171 ! volumetric heat capacity liquid wat  (J m-3 K-1)
  integer(i4b)    :: scalarLambda_drysoil            = 172 ! thermal conductivity of dry soil     (W m-1)
  integer(i4b)    :: scalarLambda_wetsoil            = 173 ! thermal conductivity of wet soil     (W m-1)
  integer(i4b)    :: scalarVolLatHt_fus              = 174 ! volumetric latent heat of fusion     (J m-3)
  integer(i4b)    :: scalarAquiferRootFrac           = 175 ! fraction of roots below the soil profile (-)
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
 ! (7) define basin-average model parameters
 ! ***********************************************************************************************************
 type, public :: iLook_bpar
  ! baseflow
  integer(i4b)    :: basin__aquiferHydCond       =  1  ! hydraulic conductivity for the aquifer (m s-1)
  integer(i4b)    :: basin__aquiferScaleFactor   =  2  ! scaling factor for aquifer storage in the big bucket (m)
  integer(i4b)    :: basin__aquiferBaseflowExp   =  3  ! baseflow exponent for the big bucket (-)
  ! within-grid routing
  integer(i4b)    :: routingGammaShape           =  4  ! shape parameter in Gamma distribution used for sub-grid routing (-)
  integer(i4b)    :: routingGammaScale           =  5  ! scale parameter in Gamma distribution used for sub-grid routing (s)
 endtype iLook_bpar

 ! ***********************************************************************************************************
 ! (8) define basin-average model variables
 ! ***********************************************************************************************************
 type, public :: iLook_bvar
  ! define derived variables
  integer(i4b)    :: basin__totalArea                =  1 ! total basin area (m2)
  ! define fluxes
  integer(i4b)    :: basin__SurfaceRunoff            =  2 ! surface runoff (m s-1)
  integer(i4b)    :: basin__SoilEjection             =  3 ! ejected water from the soil profile (m s-1)
  integer(i4b)    :: basin__SoilBaseflow             =  4 ! baseflow from the soil profile (m s-1)
  integer(i4b)    :: basin__AquiferStorage           =  5 ! aquifer storage (m s-1)  
  integer(i4b)    :: basin__AquiferRecharge          =  6 ! recharge to the aquifer (m s-1)
  integer(i4b)    :: basin__AquiferBaseflow          =  7 ! baseflow from the aquifer (m s-1)
  integer(i4b)    :: basin__AquiferTranspire         =  8 ! transpiration from the aquifer (m s-1)
  ! define variables for runoff
  integer(i4b)    :: routingRunoffFuture             =  9 ! runoff in future time steps (m s-1)
  integer(i4b)    :: routingFractionFuture           = 10 ! fraction of runoff in future time steps (-)
  integer(i4b)    :: averageInstantRunoff            = 11 ! instantaneous runoff (m s-1)
  integer(i4b)    :: averageRoutedRunoff             = 12 ! routed runoff (m s-1)
 endtype iLook_bvar

 ! ***********************************************************************************************************
 ! (X) define data structures and maximum number of variables of each type
 ! ***********************************************************************************************************
 ! define look-up structures
 type(iLook_decision),public,parameter :: iLookDECISIONS=iLook_decision(  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,&
                                                                         21, 22, 23, 24, 25)
 type(iLook_time),    public,parameter :: iLookTIME     =iLook_time    (  1,  2,  3,  4,  5)
 type(iLook_force),   public,parameter :: iLookFORCE    =iLook_force   (  1,  2,  3,  4,  5,  6,  7,  8)
 type(iLook_attr),    public,parameter :: iLookATTR     =iLook_attr    (  1,  2,  3,  4,  5,  6,  7)
 type(iLook_type),    public,parameter :: iLookTYPE     =iLook_type    (  1,  2,  3,  4,  5)
 type(iLook_param),   public,parameter :: iLookPARAM    =iLook_param   (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,&
                                                                         21, 22, 23, 24, 25, 26, 27, 28, 29, 30,&
                                                                         31, 32, 33, 34, 35, 36, 37, 38, 39, 40,&
                                                                         41, 42, 43, 44, 45, 46, 47, 48, 49, 50,&
                                                                         51, 52, 53, 54, 55, 56, 57, 58, 59, 60,&
                                                                         61, 62, 63, 64, 65, 66, 67, 68, 69, 70,&
                                                                         71, 72, 73, 74, 75, 76, 77, 78, 79, 80,&
                                                                         81, 82, 83, 84, 85, 86, 87, 88, 89, 90,&
                                                                         91, 92, 93, 94, 95, 96, 97, 98, 99,100,&
                                                                        101,102,103, 104)
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
                                                                        111,112,113,114,115,116,117,118,119,120,&
                                                                        121,122,123,124,125,126,127,128,129,130,&
                                                                        131,132,133,134,135,136,137,138,139,140,&
                                                                        141,142,143,144,145,146,147,148,149,150,&
                                                                        151,152,153,154,155,156,157,158,159,160,&
                                                                        161,162,163,164,165,166,167,168,169,170,&
                                                                        171,172,173,174,175)
 type(iLook_index),   public,parameter :: iLookINDEX    =ilook_index   (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10)
 type(iLook_bpar),    public,parameter :: iLookBPAR     =ilook_bpar    (  1,  2,  3,  4,  5)
 type(iLook_bvar),    public,parameter :: iLookBVAR     =ilook_bvar    (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12)
 ! define maximum number of variables of each type
 integer(i4b),parameter,public :: maxvarDecisions= 25
 integer(i4b),parameter,public :: maxvarTime     = 5
 integer(i4b),parameter,public :: maxvarForc     = 8
 integer(i4b),parameter,public :: maxvarAttr     = 7
 integer(i4b),parameter,public :: maxvarType     = 5
 integer(i4b),parameter,public :: maxvarMpar     = 104
 integer(i4b),parameter,public :: maxvarMvar     = 175
 integer(i4b),parameter,public :: maxvarIndx     = 10
 integer(i4b),parameter,public :: maxvarBpar     = 5
 integer(i4b),parameter,public :: maxvarBvar     = 12
 ! ***********************************************************************************************************
END MODULE var_lookup

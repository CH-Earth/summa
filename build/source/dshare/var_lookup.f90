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
  ! Noah-MP options
  integer(i4b)    :: soilCatTbl       = 1  ! soil-category dateset
  integer(i4b)    :: vegeParTbl       = 2  ! vegetation category dataset
  integer(i4b)    :: soilStress       = 3  ! choice of function for the soil moisture control on stomatal resistance
  integer(i4b)    :: stomResist       = 4  ! choice of function for stomatal resistance
  ! FUSE options
  integer(i4b)    :: num_method       = 5  ! choice of numerical method
  integer(i4b)    :: fDerivMeth       = 6  ! method used to calculate flux derivatives
  integer(i4b)    :: LAI_method       = 7  ! method used to determine LAI and SAI
  integer(i4b)    :: f_Richards       = 8  ! form of richards' equation
  integer(i4b)    :: groundwatr       = 9  ! choice of groundwater parameterization
  integer(i4b)    :: hc_profile       = 10 ! choice of hydraulic conductivity profile
  integer(i4b)    :: bcUpprTdyn       = 11 ! type of upper boundary condition for thermodynamics 
  integer(i4b)    :: bcLowrTdyn       = 12 ! type of lower boundary condition for thermodynamics
  integer(i4b)    :: bcUpprSoiH       = 13 ! type of upper boundary condition for soil hydrology
  integer(i4b)    :: bcLowrSoiH       = 14 ! type of lower boundary condition for soil hydrology
  integer(i4b)    :: veg_traits       = 15 ! choice of parameterization for vegetation roughness length and displacement height
  integer(i4b)    :: canopyEmis       = 16 ! choice of parameterization for canopy emissivity 
  integer(i4b)    :: astability       = 17 ! choice of stability function
  integer(i4b)    :: alb_method       = 18 ! choice of albedo representation
  integer(i4b)    :: snowLayers       = 19 ! choice of method to combine and sub-divide snow layers
  integer(i4b)    :: compaction       = 20 ! choice of compaction routine
  integer(i4b)    :: thermlcond       = 21 ! choice of thermal conductivity representation
  integer(i4b)    :: subRouting       = 22 ! choice of method for sub-grid routing
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
  integer(i4b)    :: latitude         = 1  ! latitude    (degrees north)
  integer(i4b)    :: longitude        = 2  ! longitude   (degrees east)
  integer(i4b)    :: elevation        = 3  ! elevation   (m)
  integer(i4b)    :: mHeight          = 4  ! measurement height above bare ground (m)
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
  integer(i4b)    :: upperBoundHead       = 1   ! matric head of the upper boundary (m)
  integer(i4b)    :: lowerBoundHead       = 2   ! matric head of the lower boundary (m)
  integer(i4b)    :: upperBoundTheta      = 3   ! volumetric liquid water content of the upper boundary (-)
  integer(i4b)    :: lowerBoundTheta      = 4   ! volumetric liquid water content of the lower boundary (-)
  integer(i4b)    :: upperBoundTemp       = 5   ! temperature of the upper boundary (K)
  integer(i4b)    :: lowerBoundTemp       = 6   ! temperature of the lower boundary (K)
  ! precipitation partitioning
  integer(i4b)    :: tempCritRain         = 7   ! critical temperature where precipitation is rain (K)
  integer(i4b)    :: tempRangeTimestep    = 8   ! temperature range over the time step (K)
  ! freezing curve for snow
  integer(i4b)    :: snowfrz_scale        = 9   ! scaling parameter for the freezing curve for snow (K-1)
  ! snow albedo
  integer(i4b)    :: snw_crit             = 10  ! critical mass necessary for albedo refreshment (kg m-2)
  integer(i4b)    :: alb_fresh            = 11  ! fresh snow albedo (-)
  integer(i4b)    :: alb_dry              = 12  ! minimum snow albedo during winter (-)
  integer(i4b)    :: alb_wet              = 13  ! minimum snow albedo during spring (-)
  integer(i4b)    :: alb_decay            = 14  ! temporal decay factor for snow albedo (s-1)
  integer(i4b)    :: alb_scale            = 15  ! albedo scaling factor (s)
  integer(i4b)    :: soot_load            = 16  ! temporal decay in snow albedo associated with the soot load (days-1)
  ! radiation transfer within snow
  integer(i4b)    :: radExt_snow          = 17  ! extinction coefficient for radiation penetration into the snowpack (m-1)
  integer(i4b)    :: Frad_direct          = 18  ! fraction of direct solar radiation (-)
  integer(i4b)    :: Frad_vis             = 19  ! fraction of radiation in the visible part of the spectrum (-)
  ! new snow density
  integer(i4b)    :: newSnowDenMin        = 20  ! minimum new snow density (kg m-3)  
  integer(i4b)    :: newSnowDenMult       = 21  ! multiplier for new snow density (kg m-3)
  integer(i4b)    :: newSnowDenScal       = 22  ! scaling factor for new snow density (K)
  ! snow compaction
  integer(i4b)    :: densScalGrowth       = 23  ! density scaling factor for grain growth (kg-1 m3)
  integer(i4b)    :: tempScalGrowth       = 24  ! temperature scaling factor for grain growth (K-1)
  integer(i4b)    :: grainGrowthRate      = 25  ! rate of grain growth (s-1)
  integer(i4b)    :: densScalOvrbdn       = 26  ! density scaling factor for overburden pressure (kg-1 m3)
  integer(i4b)    :: tempScalOvrbdn       = 27  ! temperature scaling factor for overburden pressure (K-1)
  integer(i4b)    :: base_visc            = 28  ! viscosity coefficient at T=T_frz and snow density=0  (kg s m-2)
  ! water flow within snow
  integer(i4b)    :: Fcapil               = 29  ! capillary retention as a fraction of the total pore volume (-)
  integer(i4b)    :: k_snow               = 30  ! hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
  integer(i4b)    :: mw_exp               = 31  ! exponent for meltwater flow (-)
  ! turbulent heat fluxes
  integer(i4b)    :: z0Snow               = 32  ! roughness length of snow (m)
  integer(i4b)    :: z0Soil               = 33  ! roughness length of bare soil below the canopy (m)
  integer(i4b)    :: z0Canopy             = 34  ! roughness length of the canopy (m)
  integer(i4b)    :: zpdFraction          = 35  ! zero plane displacement / canopy height (-)
  integer(i4b)    :: critRichNumber       = 36  ! critical value for the bulk Richardson number (-)
  integer(i4b)    :: Louis79_bparam       = 37  ! parameter in Louis (1979) stability function (-)
  integer(i4b)    :: Louis79_cStar        = 38  ! parameter in Louis (1979) stability function (-)
  integer(i4b)    :: Mahrt87_eScale       = 39  ! exponential scaling factor in the Mahrt (1987) stability function (-)
  integer(i4b)    :: leafExchangeCoeff    = 40  ! turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
  integer(i4b)    :: windReductionParam   = 41  ! canopy wind reduction parameter (-)
  ! vegetation properties
  integer(i4b)    :: winterSAI            = 42  ! stem area index prior to the start of the growing season (m2 m-2)
  integer(i4b)    :: summerLAI            = 43  ! maximum leaf area index at the peak of the growing season (m2 m-2)
  integer(i4b)    :: rootingDepth         = 44  ! rooting depth (m)
  integer(i4b)    :: rootDistExp          = 45  ! exponent controlling the vertical distribution of root density (-)
  integer(i4b)    :: plantWiltPsi         = 46  ! matric head at wilting point (m)
  integer(i4b)    :: soilStressParam      = 47  ! parameter in the exponential soil stress function
  integer(i4b)    :: critSoilWilting      = 48  ! critical vol. liq. water content when plants are wilting (-) 
  integer(i4b)    :: critSoilTranspire    = 49  ! critical vol. liq. water content when transpiration is limited (-)
  integer(i4b)    :: critAquiferTranspire = 50  ! critical aquifer storage value when transpiration is limited (m)
  integer(i4b)    :: leafDimension        = 51  ! characteristic leaf dimension (m)
  integer(i4b)    :: heightCanopyTop      = 52  ! height of top of the vegetation canopy above ground surface (m)
  integer(i4b)    :: heightCanopyBottom   = 53  ! height of bottom of the vegetation canopy above ground surface (m)
  integer(i4b)    :: maxCanopyIce         = 54  ! maximum mass of ice on the vegetation canopy per unit lai+sai (kg m-2)
  integer(i4b)    :: maxCanopyLiquid      = 55  ! maximum mass of liq water on vegetation canopy per unit lai+sai (kg m-2)
  integer(i4b)    :: specificHeatVeg      = 56  ! specific heat of vegetation (J kg-1 K-1)
  integer(i4b)    :: maxMassVegetation    = 57  ! maximum mass of vegetation (full foliage) (kg m-2)
  ! soil properties
  integer(i4b)    :: soil_dens_intr       = 58  ! intrinsic soil density (kg m-3)
  integer(i4b)    :: frac_sand            = 59  ! fraction of sand (-)
  integer(i4b)    :: frac_silt            = 60  ! fraction of silt (-)
  integer(i4b)    :: frac_clay            = 61  ! fraction of clay (-)
  integer(i4b)    :: theta_sat            = 62  ! porosity (-)
  integer(i4b)    :: theta_res            = 63  ! volumetric residual water content (-)
  integer(i4b)    :: vGn_alpha            = 64  ! van Genuchten "alpha" parameter (m-1)
  integer(i4b)    :: vGn_n                = 65  ! van Genuchten "n" parameter (-)
  integer(i4b)    :: k_soil               = 66  ! hydraulic conductivity of soil (m s-1) 
  integer(i4b)    :: kAnisotropic         = 67  ! anisotropy factor for lateral hydraulic conductivity (-)
  integer(i4b)    :: zScale_TOPMODEL      = 68  ! scale factor for TOPMODEL-ish baseflow parameterization (m)
  integer(i4b)    :: compactedDepth       = 69  ! depth where k_soil reaches the compacted value given by CH78 (m)
  integer(i4b)    :: bpar_VIC             = 70  ! b-parameter in the VIC surface runoff parameterization (-)
  integer(i4b)    :: specificYield        = 71  ! specific yield (-)
  integer(i4b)    :: specificStorage      = 72  ! specific storage coefficient (m-1)
  integer(i4b)    :: aquiferScaleFactor   = 73  ! scaling factor for aquifer storage in the big bucket (m)
  integer(i4b)    :: bucketBaseflowExp    = 74  ! baseflow exponent for the big bucket (-)
  integer(i4b)    :: f_impede             = 75  ! ice impedence factor (-)
  ! within-grid routing
  integer(i4b)    :: routingGammaShape    = 76  ! shape parameter in Gamma distribution used for sub-grid routing (-)
  integer(i4b)    :: routingGammaScale    = 77  ! scale parameter in Gamma distribution used for sub-grid routing (s)
  ! algorithmic control parameters
  integer(i4b)    :: minwind              = 78  ! minimum wind speed (m s-1)
  integer(i4b)    :: minstep              = 79  ! minimum length of the time step
  integer(i4b)    :: maxstep              = 80  ! maximum length of the time step
  integer(i4b)    :: wimplicit            = 81  ! weight assigned to the start-of-step fluxes
  integer(i4b)    :: maxiter              = 82  ! maximum number of iteration 
  integer(i4b)    :: relConvTol_liquid    = 83  ! relative convergence tolerance for vol frac liq water (-)
  integer(i4b)    :: absConvTol_liquid    = 84  ! absolute convergence tolerance for vol frac liq water (-)
  integer(i4b)    :: relConvTol_matric    = 85  ! relative convergence tolerance for matric head (-)
  integer(i4b)    :: absConvTol_matric    = 86  ! absolute convergence tolerance for matric head (m)
  integer(i4b)    :: relConvTol_energy    = 87  ! relative convergence tolerance for energy (-)
  integer(i4b)    :: absConvTol_energy    = 88  ! absolute convergence tolerance for energy (J m-3)
  integer(i4b)    :: relConvTol_aquifr    = 89  ! relative convergence tolerance for aquifer storage (-)
  integer(i4b)    :: absConvTol_aquifr    = 90  ! absolute convergence tolerance for aquifer storage (J m-3)
  integer(i4b)    :: zmin                 = 91  ! minimum layer depth (m)
  integer(i4b)    :: zmax                 = 92  ! maximum layer depth (m)
  integer(i4b)    :: zminLayer1           = 93  ! minimum layer depth for the 1st (top) layer (m)
  integer(i4b)    :: zminLayer2           = 94  ! minimum layer depth for the 2nd layer (m) 
  integer(i4b)    :: zminLayer3           = 95  ! minimum layer depth for the 3rd layer (m) 
  integer(i4b)    :: zminLayer4           = 96  ! minimum layer depth for the 4th layer (m) 
  integer(i4b)    :: zminLayer5           = 97  ! minimum layer depth for the 5th (bottom) layer (m) 
  integer(i4b)    :: zmaxLayer1_lower     = 98  ! maximum layer depth for the 1st (top) layer when only 1 layer (m) 
  integer(i4b)    :: zmaxLayer2_lower     = 99  ! maximum layer depth for the 2nd layer when only 2 layers (m) 
  integer(i4b)    :: zmaxLayer3_lower     = 100 ! maximum layer depth for the 3rd layer when only 3 layers (m) 
  integer(i4b)    :: zmaxLayer4_lower     = 101 ! maximum layer depth for the 4th layer when only 4 layers (m) 
  integer(i4b)    :: zmaxLayer1_upper     = 102 ! maximum layer depth for the 1st (top) layer when > 1 layer (m) 
  integer(i4b)    :: zmaxLayer2_upper     = 103 ! maximum layer depth for the 2nd layer when > 2 layers (m) 
  integer(i4b)    :: zmaxLayer3_upper     = 104 ! maximum layer depth for the 3rd layer when > 3 layers (m) 
  integer(i4b)    :: zmaxLayer4_upper     = 105 ! maximum layer depth for the 4th layer when > 4 layers (m) 
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
  ! define scalar variables -- forcing
  integer(i4b)    :: scalarCosZenith                 = 11  ! cosine of the solar zenith angle (0-1)
  integer(i4b)    :: spectralIncomingDirect          = 12  ! incoming direct solar radiation in each wave band (W m-2)
  integer(i4b)    :: spectralIncomingDiffuse         = 13  ! incoming diffuse solar radiation in each wave band (W m-2)
  integer(i4b)    :: scalarVPair                     = 14  ! vapor pressure of the air above the vegetation canopy (Pa)
  integer(i4b)    :: scalarTwetbulb                  = 15  ! wet bulb temperature (K)
  integer(i4b)    :: scalarRainfall                  = 16  ! computed rainfall rate (kg m-2 s-1)
  integer(i4b)    :: scalarSnowfall                  = 17  ! computed snowfall rate (kg m-2 s-1)
  integer(i4b)    :: scalarSnowfallTemp              = 18  ! temperature of fresh snow (K)
  integer(i4b)    :: scalarO2air                     = 19  ! atmospheric o2 concentration (Pa)
  integer(i4b)    :: scalarCO2air                    = 20  ! atmospheric co2 concentration (Pa)
  ! define scalar variables -- state variables
  integer(i4b)    :: scalarCanopyIce                 = 21  ! mass of ice on the vegetation canopy (kg m-2)
  integer(i4b)    :: scalarCanopyLiq                 = 22  ! mass of liquid water on the vegetation canopy (kg m-2)
  integer(i4b)    :: scalarCanopyTemp                = 23  ! temperature of the vegetation canopy (K)
  integer(i4b)    :: scalarAlbedo                    = 24  ! snow albedo (-)
  integer(i4b)    :: scalarSnowAge                   = 25  ! non-dimensional snow age (-)
  integer(i4b)    :: scalarSnowDepth                 = 26  ! total snow depth (m)
  integer(i4b)    :: scalarSWE                       = 27  ! snow water equivalent (kg m-2)
  integer(i4b)    :: scalarSfcMeltPond               = 28  ! ponded water caused by melt of the "snow without a layer" (kg m-2)
  integer(i4b)    :: scalarAquiferStorage            = 29  ! relative aquifer storage -- above bottom of the soil profile (m)
  integer(i4b)    :: scalarWaterTableDepth           = 30  ! depth of the water table (m)
  ! define NOAH-MP vegetation variables -- general
  integer(i4b)    :: scalarGreenVegFraction          = 31  ! green vegetation fraction used to compute LAI (-)
  integer(i4b)    :: scalarBulkVolHeatCapVeg         = 32  ! bulk volumetric heat capacity of vegetation (J m-3 K-1)
  integer(i4b)    :: scalarRootZoneTemp              = 33  ! average temperature of the root zone (K)
  integer(i4b)    :: scalarLAI                       = 34  ! one-sided leaf area index (m2 m-2)
  integer(i4b)    :: scalarSAI                       = 35  ! one-sided stem area index (m2 m-2)
  integer(i4b)    :: scalarExposedLAI                = 36  ! exposed leaf area index after burial by snow (m2 m-2)
  integer(i4b)    :: scalarExposedSAI                = 37  ! exposed stem area index after burial by snow(m2 m-2)
  integer(i4b)    :: scalarGrowingSeasonIndex        = 38  ! growing season index (0=off, 1=on)
  integer(i4b)    :: scalarVP_CanopyAir              = 39  ! vapor pressure of the canopy air space (Pa)
  integer(i4b)    :: scalarTemp_CanopyAir            = 40  ! temperature of the canopy air space (Pa)
  ! define NOAH-MP vegetation variables -- shortwave radiation
  integer(i4b)    :: scalarCanopySunlitFraction      = 41  ! sunlit fraction of canopy (-)
  integer(i4b)    :: scalarCanopySunlitLAI           = 42  ! sunlit leaf area (-)
  integer(i4b)    :: scalarCanopyShadedLAI           = 43  ! shaded leaf area (-)
  integer(i4b)    :: scalarCanopySunlitPAR           = 44  ! average absorbed par for sunlit leaves (w m-2)
  integer(i4b)    :: scalarCanopyShadedPAR           = 45  ! average absorbed par for shaded leaves (w m-2)
  integer(i4b)    :: scalarCanopyAbsorbedSolar       = 46  ! solar radiation absorbed by canopy (W m-2)
  integer(i4b)    :: scalarGroundAbsorbedSolar       = 47  ! solar radiation absorbed by ground (W m-2)
  integer(i4b)    :: scalarTotalReflectedSolar       = 48  ! total reflected solar radiation (W m-2)
  integer(i4b)    :: scalarTotalAbsorbedSolar        = 49  ! total absorbed solar radiation (W m-2)
  integer(i4b)    :: scalarCanopyReflectedSolar      = 50  ! solar radiation reflected from the canopy (W m-2)
  integer(i4b)    :: scalarGroundReflectedSolar      = 51  ! solar radiation reflected from the ground (W m-2) 
  integer(i4b)    :: scalarBetweenCanopyGapFraction  = 52  ! between canopy gap fraction for beam (-)
  integer(i4b)    :: scalarWithinCanopyGapFraction   = 53  ! within canopy gap fraction for beam (-)
  ! define NOAH-MP vegetation variables -- longwave radiation
  integer(i4b)    :: scalarCanopyEmissivity          = 54  ! effective canopy emissivity (-)
  integer(i4b)    :: scalarLWRadCanopy               = 55  ! longwave radiation emitted from the canopy (W m-2)
  integer(i4b)    :: scalarLWRadGround               = 56  ! longwave radiation emitted at the ground surface  (W m-2)
  integer(i4b)    :: scalarLWRadUbound2Canopy        = 57  ! downward atmospheric longwave radiation absorbed by the canopy (W m-2)
  integer(i4b)    :: scalarLWRadUbound2Ground        = 58  ! downward atmospheric longwave radiation absorbed by the ground (W m-2)
  integer(i4b)    :: scalarLWRadUbound2Ubound        = 59  ! atmospheric radiation refl by ground + lost thru upper boundary (W m-2)
  integer(i4b)    :: scalarLWRadCanopy2Ubound        = 60  ! longwave radiation emitted from canopy lost thru upper boundary (W m-2)
  integer(i4b)    :: scalarLWRadCanopy2Ground        = 61  ! longwave radiation emitted from canopy absorbed by the ground (W m-2)
  integer(i4b)    :: scalarLWRadCanopy2Canopy        = 62  ! canopy longwave reflected from ground and absorbed by the canopy (W m-2)
  integer(i4b)    :: scalarLWRadGround2Ubound        = 63  ! longwave radiation emitted from ground lost thru upper boundary (W m-2)
  integer(i4b)    :: scalarLWRadGround2Canopy        = 64  ! longwave radiation emitted from ground and absorbed by the canopy (W m-2)
  integer(i4b)    :: scalarLWNetCanopy               = 65  ! net longwave radiation at the canopy (W m-2)
  integer(i4b)    :: scalarLWNetGround               = 66  ! net longwave radiation at the ground surface (W m-2)
  integer(i4b)    :: scalarLWNetUbound               = 67  ! net longwave radiation at the upper atmospheric boundary (W m-2)
  ! define NOAH-MP vegetation variables -- turbulent heat transfer
  integer(i4b)    :: scalarLatHeatSubVapCanopy       = 68  ! latent heat of sublimation/vaporization used for veg canopy (J kg-1)
  integer(i4b)    :: scalarLatHeatSubVapGround       = 69  ! latent heat of sublimation/vaporization used for ground surface (J kg-1)
  integer(i4b)    :: scalarSatVP_CanopyTemp          = 70  ! saturation vapor pressure at the temperature of vegetation canopy (Pa)
  integer(i4b)    :: scalarSatVP_GroundTemp          = 71  ! saturation vapor pressure at the temperature of the ground (Pa)
  integer(i4b)    :: scalarZ0Canopy                  = 72  ! roughness length of the canopy (m)
  integer(i4b)    :: scalarWindReductionFactor       = 73  ! canopy wind reduction factor (-)
  integer(i4b)    :: scalarZeroPlaneDisplacement     = 74  ! zero plane displacement (m) 
  integer(i4b)    :: scalarRiBulkCanopy              = 75  ! bulk Richardson number for the canopy (-)
  integer(i4b)    :: scalarRiBulkGround              = 76  ! bulk Richardson number for the ground surface (-)
  integer(i4b)    :: scalarCanopyStabilityCorrection = 77  ! stability correction for the canopy (-)
  integer(i4b)    :: scalarGroundStabilityCorrection = 78  ! stability correction for the ground surface (-)
  integer(i4b)    :: scalarEddyDiffusCanopyTop       = 79  ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
  integer(i4b)    :: scalarFrictionVelocity          = 80  ! friction velocity - canopy momentum sink (m s-1)
  integer(i4b)    :: scalarWindspdCanopyTop          = 81  ! windspeed at the top of the canopy (m s-1)
  integer(i4b)    :: scalarWindspdCanopyBottom       = 82  ! windspeed at the height of the bottom of the canopy (m s-1)
  integer(i4b)    :: scalarGroundResistance          = 83  ! below canopy aerodynamic resistance (s m-1) 
  integer(i4b)    :: scalarCanopyResistance          = 84  ! above canopy aerodynamic resistance (s m-1)
  integer(i4b)    :: scalarLeafResistance            = 85  ! mean leaf boundary layer resistance per unit leaf area (s m-1)
  integer(i4b)    :: scalarSoilResistance            = 86  ! soil surface resistance (s m-1)
  integer(i4b)    :: scalarSoilRelHumidity           = 87  ! relative humidity in the soil pores in the upper-most soil layer (-)
  integer(i4b)    :: scalarSenHeatCanopy             = 88  ! sensible heat from the canopy to the canopy air space (W m-2)
  integer(i4b)    :: scalarSenHeatGround             = 89  ! sensible heat from the ground (below canopy or non-vegetated) (W m-2)
  integer(i4b)    :: scalarLatHeatCanopyEvap         = 90  ! evaporation latent heat from the canopy to the canopy air space (W m-2) 
  integer(i4b)    :: scalarLatHeatCanopyTrans        = 91  ! transpiration latent heat from the canopy to the canopy air space (W m-2)
  integer(i4b)    :: scalarLatHeatGround             = 92  ! latent heat from the ground (below canopy or non-vegetated) (W m-2)
  integer(i4b)    :: scalarCanopyTranspiration       = 93  ! canopy transpiration (kg m-2 s-1)
  integer(i4b)    :: scalarCanopyEvaporation         = 94  ! canopy evaporation/condensation (kg m-2 s-1)
  integer(i4b)    :: scalarCanopySublimation         = 95  ! canopy sublimation/frost (kg m-2 s-1)
  integer(i4b)    :: scalarGroundEvaporation         = 96  ! ground evaporation/condensation (below canopy or non-vegetated) (kg m-2 s-1)
  integer(i4b)    :: scalarSnowSublimation           = 97  ! snow sublimation/frost (below canopy or non-vegetated) (kg m-2 s-1)
  ! define NOAH-MP vegetation variables -- transpiration
  integer(i4b)    :: scalarTranspireLim              = 98  ! aggregate soil moisture + aquifer storage limit on transpiration (-)
  integer(i4b)    :: scalarTranspireLimAqfr          = 99  ! aquifer storage limit on transpiration (-)
  integer(i4b)    :: scalarFoliageNitrogenFactor     = 100 ! foliage nitrogen concentration, 1=saturated (-)
  integer(i4b)    :: scalarStomResistSunlit          = 101 ! stomatal resistance for sunlit leaves (s m-1)
  integer(i4b)    :: scalarStomResistShaded          = 102 ! stomatal resistance for shaded leaves (s m-1)
  integer(i4b)    :: scalarPhotosynthesisSunlit      = 103 ! sunlit photosynthesis (umolco2 m-2 s-1)
  integer(i4b)    :: scalarPhotosynthesisShaded      = 104 ! shaded photosynthesis (umolco2 m-2 s-1)
  ! define NOAH-MP vegetation variables --hydrology
  integer(i4b)    :: scalarCanopyWetFraction         = 105 ! fraction of canopy that is wet
  integer(i4b)    :: scalarGroundSnowFraction        = 106 ! fraction of ground that is covered with snow (-)
  ! define scalar variables -- soil and aquifer fluxes
  integer(i4b)    :: scalarRainPlusMelt              = 107 ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  integer(i4b)    :: scalarSurfaceRunoff             = 108 ! surface runoff (m s-1)
  integer(i4b)    :: scalarInitAquiferRecharge       = 109 ! recharge to the aquifer at the start of the step (m s-1)
  integer(i4b)    :: scalarAquiferRecharge           = 110 ! recharge to the aquifer (m s-1)
  integer(i4b)    :: scalarInitAquiferTranspire      = 111 ! transpiration from the aquifer at the start of the step (m s-1)
  integer(i4b)    :: scalarAquiferTranspire          = 112 ! transpiration from the aquifer (m s-1)
  integer(i4b)    :: scalarInitAquiferBaseflow       = 113 ! baseflow from the aquifer at the start of the step (m s-1)
  integer(i4b)    :: scalarAquiferBaseflow           = 114 ! baseflow from the aquifer (m s-1)
  ! scalar variables -- sub-step average fluxes for the soil zone
  integer(i4b)    :: scalarSoilInflux                = 115 ! sub-step average: influx of water at the top of the soil profile (m s-1)
  integer(i4b)    :: scalarSoilBaseflow              = 116 ! sub-step average: total baseflow from throughout the soil profile (m s-1)
  integer(i4b)    :: scalarSoilDrainage              = 117 ! sub-step average: drainage from the bottom of the soil profile (m s-1) 
  integer(i4b)    :: scalarSoilEjection              = 118 ! sub-step average: total water ejected from all soil layers (m s-1) 
  integer(i4b)    :: scalarSoilTranspiration         = 119 ! sub-step average: total transpiration from the soil (m s-1)
  ! define scalar variables -- mass balance check
  integer(i4b)    :: scalarSoilWatBalError           = 120 ! error in the total soil water balance (kg m-2)
  integer(i4b)    :: scalarAquiferBalError           = 121 ! error in the aquifer water balance (kg m-2)
  integer(i4b)    :: scalarTotalSoilLiq              = 122 ! total mass of liquid water in the soil (kg m-2)
  integer(i4b)    :: scalarTotalSoilIce              = 123 ! total mass of ice in the soil (kg m-2)
  ! define variables at the mid-point of each layer -- domain geometry
  integer(i4b)    :: mLayerDepth                     = 124 ! depth of each layer (m)
  integer(i4b)    :: mLayerHeight                    = 125 ! height at the mid-point of each layer (m)
  integer(i4b)    :: mLayerRootDensity               = 126 ! fraction of roots in each soil layer (-)
  ! define variables at the mid-point of each layer -- coupled energy and mass
  integer(i4b)    :: mLayerTemp                      = 127 ! temperature of each layer (K)
  integer(i4b)    :: mLayerVolFracAir                = 128 ! volumetric fraction of air in each layer (-)
  integer(i4b)    :: mLayerVolFracIce                = 129 ! volumetric fraction of ice water in each layer (-)
  integer(i4b)    :: mLayerVolFracLiq                = 130 ! volumetric fraction of liquid water in each layer (-)
  integer(i4b)    :: mLayerVolHtCapBulk              = 131 ! volumetric heat capacity in each layer (J m-3 K-1)
  integer(i4b)    :: mLayerTcrit                     = 132 ! critical soil temperature above which all water is unfrozen (K)
  integer(i4b)    :: mLayerdTheta_dTk                = 133 ! derivative in volumetric liquid water content wrt temperature (K-1)
  integer(i4b)    :: mLayerThermalC                  = 134 ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
  integer(i4b)    :: mLayerRadCondFlux               = 135 ! temporal derivative in energy from radiative and conductive flux (J m-2 s-1)  
  integer(i4b)    :: mLayerMeltFreeze                = 136 ! rate of ice content change from melt/freeze in each layer (kg m-3 s-1)
  integer(i4b)    :: mLayerInfilFreeze               = 137 ! rate of ice content change by freezing infiltrating flux (kg m-3 s-1)
  integer(i4b)    :: mLayerSatHydCond                = 138 ! saturated hydraulic conductivity in each layer (m s-1)
  integer(i4b)    :: mLayerMatricHead                = 139 ! matric head of water in the soil (m)
  integer(i4b)    :: mLayerdTheta_dPsi               = 140 ! derivative in the soil water characteristic (m-1)
  integer(i4b)    :: mLayerdPsi_dTheta               = 141 ! derivative in the soil water characteristic (m)
  integer(i4b)    :: mLayerThetaResid                = 142 ! residual volumetric water content in each snow layer (-)
  integer(i4b)    :: mLayerPoreSpace                 = 143 ! total pore space in each snow layer (-)
  integer(i4b)    :: mLayerTranspireLim              = 144 ! soil moist & veg limit on transpiration for each layer (-)
  integer(i4b)    :: mLayerInitTranspire             = 145 ! transpiration loss from each soil layer at the start of the step (kg m-2 s-1)
  integer(i4b)    :: mLayerTranspire                 = 146 ! transpiration loss from each soil layer (kg m-2 s-1)
  integer(i4b)    :: mLayerInitEjectWater            = 147 ! water ejected from each soil layer at the start-of-step (m s-1)
  integer(i4b)    :: mLayerEjectWater                = 148 ! water ejected from each soil layer (m s-1)
  integer(i4b)    :: mLayerInitBaseflow              = 149 ! baseflow from each soil layer at the start of the time step (m s-1)
  integer(i4b)    :: mLayerBaseflow                  = 150 ! baseflow from each soil layer (m s-1)
  ! define variables at the interface of each layer
  integer(i4b)    :: iLayerHeight                    = 151 ! height of the layer interface; top of soil = 0 (m)
  integer(i4b)    :: iLayerThermalC                  = 152 ! thermal conductivity at the interface of each layer (W m-1 K-1)
  integer(i4b)    :: iLayerInitNrgFlux               = 153 ! energy flux at layer interfaces at the start of the time step (W m-2) 
  integer(i4b)    :: iLayerNrgFlux                   = 154 ! energy flux at layer interfaces at the end of the time step (W m-2) 
  integer(i4b)    :: iLayerSatHydCond                = 155 ! saturated hydraulic conductivity at each layer interface (m s-1)
  integer(i4b)    :: iLayerInitLiqFluxSnow           = 156 ! liquid flux at snow layer interfaces at the start of the time step (m s-1) 
  integer(i4b)    :: iLayerInitLiqFluxSoil           = 157 ! liquid flux at soil layer interfaces at the start of the time step (m s-1) 
  integer(i4b)    :: iLayerLiqFluxSnow               = 158 ! liquid flux at snow layer interfaces at the end of the time step (m s-1)
  integer(i4b)    :: iLayerLiqFluxSoil               = 159 ! liquid flux at soil layer interfaces at the end of the time step (m s-1) 
  ! define variables for runoff
  integer(i4b)    :: routingRunoffFuture             = 160 ! runoff in future time steps (m s-1)
  integer(i4b)    :: routingFractionFuture           = 161 ! fraction of runoff in future time steps (-)
  integer(i4b)    :: averageInstantRunoff            = 162 ! instantaneous runoff (m s-1)
  integer(i4b)    :: averageRoutedRunoff             = 163 ! routed runoff (m s-1)
  ! define derived variables
  integer(i4b)    :: scalarVGn_m                     = 164 ! van Genuchten "m" parameter (-)
  integer(i4b)    :: scalarKappa                     = 165 ! constant in the freezing curve function (m K-1)
  integer(i4b)    :: scalarVolHtCap_air              = 166 ! volumetric heat capacity air         (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_ice              = 167 ! volumetric heat capacity ice         (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_soil             = 168 ! volumetric heat capacity dry soil    (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_water            = 169 ! volumetric heat capacity liquid wat  (J m-3 K-1)
  integer(i4b)    :: scalarLambda_drysoil            = 170 ! thermal conductivity of dry soil     (W m-1)
  integer(i4b)    :: scalarLambda_wetsoil            = 171 ! thermal conductivity of wet soil     (W m-1)
  integer(i4b)    :: scalarVolLatHt_fus              = 172 ! volumetric latent heat of fusion     (J m-3)
  integer(i4b)    :: scalarAquiferRootFrac           = 173 ! fraction of roots below the soil profile (-)
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
 type(iLook_decision),public,parameter :: iLookDECISIONS=iLook_decision(  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,&
                                                                         21)
 type(iLook_time),    public,parameter :: iLookTIME     =iLook_time    (  1,  2,  3,  4,  5)
 type(iLook_force),   public,parameter :: iLookFORCE    =iLook_force   (  1,  2,  3,  4,  5,  6,  7,  8)
 type(iLook_attr),    public,parameter :: iLookATTR     =iLook_attr    (  1,  2,  3,  4)
 type(iLook_type),    public,parameter :: iLookTYPE     =iLook_type    (  1,  2,  3)
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
                                                                        101,102,103,104,105)
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
                                                                        171,172,173)
 type(iLook_index),   public,parameter :: iLookINDEX    =ilook_index   (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10)
 ! define maximum number of variables of each type
 integer(i4b),parameter,public :: maxvarDecisions= 21
 integer(i4b),parameter,public :: maxvarTime     = 5
 integer(i4b),parameter,public :: maxvarForc     = 8
 integer(i4b),parameter,public :: maxvarAttr     = 4
 integer(i4b),parameter,public :: maxvarType     = 3
 integer(i4b),parameter,public :: maxvarMpar     = 105
 integer(i4b),parameter,public :: maxvarMvar     = 173
 integer(i4b),parameter,public :: maxvarIndx     = 10
 ! ***********************************************************************************************************
END MODULE var_lookup

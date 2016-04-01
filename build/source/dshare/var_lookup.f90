! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE var_lookup
 ! defines named variables used to index array elements
 USE nrtype
 implicit none
 private
 ! local variables
 integer(i4b),parameter     :: imiss = -999             ! missing value: used to initialize named variables
 integer(i4b),parameter     :: iLength=storage_size(1)  ! size of an integer


 ! ***********************************************************************************************************
 ! (0) define model decisions
 ! ***********************************************************************************************************
 type, public  ::  iLook_decision
  integer(i4b)    :: simulStart = imiss     ! simulation start time
  integer(i4b)    :: simulFinsh = imiss     ! simulation end time
  integer(i4b)    :: soilCatTbl = imiss     ! soil-category dateset
  integer(i4b)    :: vegeParTbl = imiss     ! vegetation category dataset
  integer(i4b)    :: soilStress = imiss     ! choice of function for the soil moisture control on stomatal resistance
  integer(i4b)    :: stomResist = imiss     ! choice of function for stomatal resistance
  integer(i4b)    :: bbTempFunc = imiss     ! Ball-Berry: leaf temperature controls on photosynthesis + stomatal resistance
  integer(i4b)    :: bbHumdFunc = imiss     ! Ball-Berry: humidity controls on stomatal resistance
  integer(i4b)    :: bbElecFunc = imiss     ! Ball-Berry: dependence of photosynthesis on PAR
  integer(i4b)    :: bbCO2point = imiss     ! Ball-Berry: use of CO2 compensation point to calculate stomatal resistance
  integer(i4b)    :: bbNumerics = imiss     ! Ball-Berry: iterative numerical solution method
  integer(i4b)    :: bbAssimFnc = imiss     ! Ball-Berry: controls on carbon assimilation
  integer(i4b)    :: bbCanIntg8 = imiss     ! Ball-Berry: scaling of photosynthesis from the leaf to the canopy
  integer(i4b)    :: num_method = imiss     ! choice of numerical method
  integer(i4b)    :: fDerivMeth = imiss     ! method used to calculate flux derivatives
  integer(i4b)    :: LAI_method = imiss     ! method used to determine LAI and SAI
  integer(i4b)    :: cIntercept = imiss     ! choice of parameterization for canopy interception
  integer(i4b)    :: f_Richards = imiss     ! form of richards' equation
  integer(i4b)    :: groundwatr = imiss     ! choice of groundwater parameterization
  integer(i4b)    :: hc_profile = imiss     ! choice of hydraulic conductivity profile
  integer(i4b)    :: bcUpprTdyn = imiss     ! type of upper boundary condition for thermodynamics
  integer(i4b)    :: bcLowrTdyn = imiss     ! type of lower boundary condition for thermodynamics
  integer(i4b)    :: bcUpprSoiH = imiss     ! type of upper boundary condition for soil hydrology
  integer(i4b)    :: bcLowrSoiH = imiss     ! type of lower boundary condition for soil hydrology
  integer(i4b)    :: veg_traits = imiss     ! choice of parameterization for vegetation roughness length and displacement height
  integer(i4b)    :: rootProfil = imiss     ! choice of parameterization for the rooting profile
  integer(i4b)    :: canopyEmis = imiss     ! choice of parameterization for canopy emissivity
  integer(i4b)    :: snowIncept = imiss     ! choice of parameterization for snow interception
  integer(i4b)    :: windPrfile = imiss     ! choice of canopy wind profile
  integer(i4b)    :: astability = imiss     ! choice of stability function
  integer(i4b)    :: canopySrad = imiss     ! choice of method for canopy shortwave radiation
  integer(i4b)    :: alb_method = imiss     ! choice of albedo representation
  integer(i4b)    :: snowLayers = imiss     ! choice of method to combine and sub-divide snow layers
  integer(i4b)    :: compaction = imiss     ! choice of compaction routine
  integer(i4b)    :: thCondSnow = imiss     ! choice of thermal conductivity representation for snow
  integer(i4b)    :: thCondSoil = imiss     ! choice of thermal conductivity representation for soil
  integer(i4b)    :: spatial_gw = imiss     ! choice of method for spatial representation of groundwater
  integer(i4b)    :: subRouting = imiss     ! choice of method for sub-grid routing
  integer(i4b)    :: snowDenNew = imiss     ! choice of method for new snow density
 endtype iLook_decision

 ! ***********************************************************************************************************
 ! (1) define model time
 ! ***********************************************************************************************************
 type, public  ::  iLook_time
  integer(i4b)    :: iyyy       = imiss     ! year
  integer(i4b)    :: im         = imiss     ! month
  integer(i4b)    :: id         = imiss     ! day
  integer(i4b)    :: ih         = imiss     ! hour
  integer(i4b)    :: imin       = imiss     ! minute
 endtype iLook_time

 ! ***********************************************************************************************************
 ! (2) define model forcing data
 ! ***********************************************************************************************************
 type, public  ::  iLook_force
  integer(i4b)    :: time       = imiss     ! time since time reference       (s)
  integer(i4b)    :: pptrate    = imiss     ! precipitation rate              (kg m-2 s-1)
  integer(i4b)    :: airtemp    = imiss     ! air temperature                 (K)
  integer(i4b)    :: spechum    = imiss     ! specific humidity               (g/g)
  integer(i4b)    :: windspd    = imiss     ! windspeed                       (m/s)
  integer(i4b)    :: SWRadAtm   = imiss     ! downwelling shortwave radiaiton (W m-2)
  integer(i4b)    :: LWRadAtm   = imiss     ! downwelling longwave radiation  (W m-2)
  integer(i4b)    :: airpres    = imiss     ! pressure                        (Pa)
 endtype iLook_force

 ! ***********************************************************************************************************
 ! (3) define local attributes
 ! ***********************************************************************************************************
 type, public  ::  iLook_attr
  integer(i4b)    :: latitude      = imiss  ! latitude (degrees north)
  integer(i4b)    :: longitude     = imiss  ! longitude (degrees east)
  integer(i4b)    :: elevation     = imiss  ! elevation (m)
  integer(i4b)    :: tan_slope     = imiss  ! tan water table slope, taken as tan local ground surface slope (-)
  integer(i4b)    :: contourLength = imiss  ! length of contour at downslope edge of HRU (m)
  integer(i4b)    :: HRUarea       = imiss  ! area of each HRU  (m2)
  integer(i4b)    :: mHeight       = imiss  ! measurement height above bare ground (m)
 end type iLook_attr

 ! ***********************************************************************************************************
 ! (4) define local classification of veg, soil, etc.
 ! ***********************************************************************************************************
 type, public  ::  iLook_type
  integer(i4b)    :: hruIndex      = imiss  ! index defining hydrologic response unit (-)
  integer(i4b)    :: vegTypeIndex  = imiss  ! index defining vegetation type (-)
  integer(i4b)    :: soilTypeIndex = imiss  ! index defining soil type (-)
  integer(i4b)    :: slopeTypeIndex= imiss  ! index defining slope (-)
  integer(i4b)    :: downHRUindex  = imiss  ! index of downslope HRU (0 = basin outlet)
 end type iLook_type

 ! ***********************************************************************************************************
 ! (5) define model parameters
 ! ***********************************************************************************************************
 type, public  ::  iLook_param
  ! boundary conditions
  integer(i4b)    :: upperBoundHead        = imiss    ! matric head of the upper boundary (m)
  integer(i4b)    :: lowerBoundHead        = imiss    ! matric head of the lower boundary (m)
  integer(i4b)    :: upperBoundTheta       = imiss    ! volumetric liquid water content of the upper boundary (-)
  integer(i4b)    :: lowerBoundTheta       = imiss    ! volumetric liquid water content of the lower boundary (-)
  integer(i4b)    :: upperBoundTemp        = imiss    ! temperature of the upper boundary (K)
  integer(i4b)    :: lowerBoundTemp        = imiss    ! temperature of the lower boundary (K)
  ! precipitation partitioning
  integer(i4b)    :: tempCritRain          = imiss    ! critical temperature where precipitation is rain (K)
  integer(i4b)    :: tempRangeTimestep     = imiss    ! temperature range over the time step (K)
  integer(i4b)    :: frozenPrecipMultip    = imiss    ! frozen precipitation multiplier (-)
  ! freezing curve for snow
  integer(i4b)    :: snowfrz_scale         = imiss    ! scaling parameter for the freezing curve for snow (K-1)
  ! snow albedo
  integer(i4b)    :: albedoMax             = imiss    ! maximum snow albedo for a single spectral band (-)
  integer(i4b)    :: albedoMinWinter       = imiss    ! minimum snow albedo during winter for a single spectral band (-)
  integer(i4b)    :: albedoMinSpring       = imiss    ! minimum snow albedo during spring for a single spectral band (-)
  integer(i4b)    :: albedoMaxVisible      = imiss    ! maximum snow albedo in the visible part of the spectrum (-)
  integer(i4b)    :: albedoMinVisible      = imiss    ! minimum snow albedo in the visible part of the spectrum (-)
  integer(i4b)    :: albedoMaxNearIR       = imiss    ! maximum snow albedo in the near infra-red part of the spectrum (-)
  integer(i4b)    :: albedoMinNearIR       = imiss    ! minimum snow albedo in the near infra-red part of the spectrum (-)
  integer(i4b)    :: albedoDecayRate       = imiss    ! albedo decay rate (s)
  integer(i4b)    :: albedoSootLoad        = imiss    ! soot load factor (-)
  integer(i4b)    :: albedoRefresh         = imiss    ! critical mass necessary for albedo refreshment (kg m-2)
  ! radiation transfer within snow
  integer(i4b)    :: radExt_snow           = imiss    ! extinction coefficient for radiation penetration into the snowpack (m-1)
  integer(i4b)    :: directScale           = imiss    ! scaling factor for fractional driect radiaion parameterization (-)
  integer(i4b)    :: Frad_direct           = imiss    ! maximum fraction of direct solar radiation (-)
  integer(i4b)    :: Frad_vis              = imiss    ! fraction of radiation in the visible part of the spectrum (-)
  ! new snow density
  integer(i4b)    :: newSnowDenMin         = imiss    ! minimum new snow density (kg m-3)
  integer(i4b)    :: newSnowDenMult        = imiss    ! multiplier for new snow density (kg m-3)
  integer(i4b)    :: newSnowDenScal        = imiss    ! scaling factor for new snow density (K)
  integer(i4b)    :: constSnowDen          = imiss    ! constDens, Constant new snow density (kg m-3)
  integer(i4b)    :: newSnowDenAdd         = imiss    ! Pahaut 1976, additive factor for new snow density (kg m-3)
  integer(i4b)    :: newSnowDenMultTemp    = imiss    ! Pahaut 1976, multiplier for new snow density applied to air temperature (kg m-3 K-1)
  integer(i4b)    :: newSnowDenMultWind    = imiss    ! Pahaut 1976, multiplier for new snow density applied to wind speed (kg m-7/2 s-1/2)
  integer(i4b)    :: newSnowDenMultAnd     = imiss    ! Anderson 1976, multiplier for new snow density for Anderson function (K-1)
  integer(i4b)    :: newSnowDenBase        = imiss    ! Anderson 1976, base value that is rasied to the (3/2) power (K)
  ! snow compaction
  integer(i4b)    :: densScalGrowth        = imiss    ! density scaling factor for grain growth (kg-1 m3)
  integer(i4b)    :: tempScalGrowth        = imiss    ! temperature scaling factor for grain growth (K-1)
  integer(i4b)    :: grainGrowthRate       = imiss    ! rate of grain growth (s-1)
  integer(i4b)    :: densScalOvrbdn        = imiss    ! density scaling factor for overburden pressure (kg-1 m3)
  integer(i4b)    :: tempScalOvrbdn        = imiss    ! temperature scaling factor for overburden pressure (K-1)
  integer(i4b)    :: baseViscosity         = imiss    ! viscosity coefficient at T=T_frz and snow density=0  (kg s m-2)
  ! water flow within snow
  integer(i4b)    :: Fcapil                = imiss    ! capillary retention as a fraction of the total pore volume (-)
  integer(i4b)    :: k_snow                = imiss    ! hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
  integer(i4b)    :: mw_exp                = imiss    ! exponent for meltwater flow (-)
  ! turbulent heat fluxes
  integer(i4b)    :: z0Snow                = imiss    ! roughness length of snow (m)
  integer(i4b)    :: z0Soil                = imiss    ! roughness length of bare soil below the canopy (m)
  integer(i4b)    :: z0Canopy              = imiss    ! roughness length of the canopy (m)
  integer(i4b)    :: zpdFraction           = imiss    ! zero plane displacement / canopy height (-)
  integer(i4b)    :: critRichNumber        = imiss    ! critical value for the bulk Richardson number (-)
  integer(i4b)    :: Louis79_bparam        = imiss    ! parameter in Louis (1979) stability function (-)
  integer(i4b)    :: Louis79_cStar         = imiss    ! parameter in Louis (1979) stability function (-)
  integer(i4b)    :: Mahrt87_eScale        = imiss    ! exponential scaling factor in the Mahrt (1987) stability function (-)
  integer(i4b)    :: leafExchangeCoeff     = imiss    ! turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
  integer(i4b)    :: windReductionParam    = imiss    ! canopy wind reduction parameter (-)
  ! stomatal conductance
  integer(i4b)    :: Kc25                  = imiss    ! Michaelis-Menten constant for CO2 at 25 degrees C (umol mol-1)
  integer(i4b)    :: Ko25                  = imiss    ! Michaelis-Menten constant for O2 at 25 degrees C (mol mol-1)
  integer(i4b)    :: Kc_qFac               = imiss    ! factor in the q10 function defining temperature controls on Kc (-)
  integer(i4b)    :: Ko_qFac               = imiss    ! factor in the q10 function defining temperature controls on Ko (-)
  integer(i4b)    :: kc_Ha                 = imiss    ! activation energy for the Michaelis-Menten constant for CO2 (J mol-1)
  integer(i4b)    :: ko_Ha                 = imiss    ! activation energy for the Michaelis-Menten constant for O2 (J mol-1)
  integer(i4b)    :: vcmax25_canopyTop     = imiss    ! potential carboxylation rate at 25 degrees C at the canopy top (umol co2 m-2 s-1)
  integer(i4b)    :: vcmax_qFac            = imiss    ! factor in the q10 function defining temperature controls on vcmax (-)
  integer(i4b)    :: vcmax_Ha              = imiss    ! activation energy in the vcmax function (J mol-1)
  integer(i4b)    :: vcmax_Hd              = imiss    ! deactivation energy in the vcmax function (J mol-1)
  integer(i4b)    :: vcmax_Sv              = imiss    ! entropy term in the vcmax function (J mol-1 K-1)
  integer(i4b)    :: vcmax_Kn              = imiss    ! foliage nitrogen decay coefficient (-)
  integer(i4b)    :: jmax25_scale          = imiss    ! scaling factor to relate jmax25 to vcmax25 (-)
  integer(i4b)    :: jmax_Ha               = imiss    ! activation energy in the jmax function (J mol-1)
  integer(i4b)    :: jmax_Hd               = imiss    ! deactivation energy in the jmax function (J mol-1)
  integer(i4b)    :: jmax_Sv               = imiss    ! entropy term in the jmax function (J mol-1 K-1)
  integer(i4b)    :: fractionJ             = imiss    ! fraction of light lost by other than the chloroplast lamellae (-)
  integer(i4b)    :: quantamYield          = imiss    ! quantam yield (mol e mol-1 q)
  integer(i4b)    :: vpScaleFactor         = imiss    ! vapor pressure scaling factor in stomatal conductance function (Pa)
  integer(i4b)    :: cond2photo_slope      = imiss    ! slope of conductance-photosynthesis relationship (-)
  integer(i4b)    :: minStomatalConductance= imiss    ! minimum stomatal conductance (umol H2O m-2 s-1)
  ! vegetation properties
  integer(i4b)    :: winterSAI             = imiss    ! stem area index prior to the start of the growing season (m2 m-2)
  integer(i4b)    :: summerLAI             = imiss    ! maximum leaf area index at the peak of the growing season (m2 m-2)
  integer(i4b)    :: rootScaleFactor1      = imiss    ! 1st scaling factor (a) in Y = 1 - 0.5*( exp(-aZ) + exp(-bZ) )  (m-1)
  integer(i4b)    :: rootScaleFactor2      = imiss    ! 2nd scaling factor (b) in Y = 1 - 0.5*( exp(-aZ) + exp(-bZ) )  (m-1)
  integer(i4b)    :: rootingDepth          = imiss    ! rooting depth (m)
  integer(i4b)    :: rootDistExp           = imiss    ! exponent controlling the vertical distribution of root density (-)
  integer(i4b)    :: plantWiltPsi          = imiss    ! matric head at wilting point (m)
  integer(i4b)    :: soilStressParam       = imiss    ! parameter in the exponential soil stress function
  integer(i4b)    :: critSoilWilting       = imiss    ! critical vol. liq. water content when plants are wilting (-)
  integer(i4b)    :: critSoilTranspire     = imiss    ! critical vol. liq. water content when transpiration is limited (-)
  integer(i4b)    :: critAquiferTranspire  = imiss    ! critical aquifer storage value when transpiration is limited (m)
  integer(i4b)    :: minStomatalResistance = imiss    ! minimum canopy resistance (s m-1)
  integer(i4b)    :: leafDimension         = imiss    ! characteristic leaf dimension (m)
  integer(i4b)    :: heightCanopyTop       = imiss    ! height of top of the vegetation canopy above ground surface (m)
  integer(i4b)    :: heightCanopyBottom    = imiss    ! height of bottom of the vegetation canopy above ground surface (m)
  integer(i4b)    :: specificHeatVeg       = imiss    ! specific heat of vegetation (J kg-1 K-1)
  integer(i4b)    :: maxMassVegetation     = imiss    ! maximum mass of vegetation (full foliage) (kg m-2)
  integer(i4b)    :: throughfallScaleSnow  = imiss    ! scaling factor for throughfall (snow) (-)
  integer(i4b)    :: throughfallScaleRain  = imiss    ! scaling factor for throughfall (rain) (-)
  integer(i4b)    :: refInterceptCapSnow   = imiss    ! reference canopy interception capacity per unit leaf area (snow) (kg m-2)
  integer(i4b)    :: refInterceptCapRain   = imiss    ! canopy interception capacity per unit leaf area (rain) (kg m-2)
  integer(i4b)    :: snowUnloadingCoeff    = imiss    ! time constant for unloading of snow from the forest canopy (s-1)
  integer(i4b)    :: canopyDrainageCoeff   = imiss    ! time constant for drainage of liquid water from the forest canopy (s-1)
  integer(i4b)    :: ratioDrip2Unloading   = imiss    ! ratio of canopy drip to unloading of snow from the forest canopy (-)
  integer(i4b)    :: canopyWettingFactor   = imiss    ! maximum wetted fraction of the canopy (-)
  integer(i4b)    :: canopyWettingExp      = imiss    ! exponent in canopy wetting function (-)
  ! soil properties
  integer(i4b)    :: soil_dens_intr        = imiss    ! intrinsic soil density (kg m-3)
  integer(i4b)    :: thCond_soil           = imiss    ! thermal conductivity of soil (W m-1 K-1)
  integer(i4b)    :: frac_sand             = imiss    ! fraction of sand (-)
  integer(i4b)    :: frac_silt             = imiss    ! fraction of silt (-)
  integer(i4b)    :: frac_clay             = imiss    ! fraction of clay (-)
  integer(i4b)    :: fieldCapacity         = imiss    ! field capacity (-)
  integer(i4b)    :: wettingFrontSuction   = imiss    ! Green-Ampt wetting front suction (m)
  integer(i4b)    :: theta_mp              = imiss    ! volumetric liquid water content when macropore flow begins (-)
  integer(i4b)    :: theta_sat             = imiss    ! porosity (-)
  integer(i4b)    :: theta_res             = imiss    ! volumetric residual water content (-)
  integer(i4b)    :: vGn_alpha             = imiss    ! van Genuchten "alpha" parameter (m-1)
  integer(i4b)    :: vGn_n                 = imiss    ! van Genuchten "n" parameter (-)
  integer(i4b)    :: mpExp                 = imiss    ! empirical exponent in macropore flow equation (-)
  integer(i4b)    :: k_soil                = imiss    ! hydraulic conductivity of soil (m s-1)
  integer(i4b)    :: k_macropore           = imiss    ! saturated hydraulic conductivity for macropores (m s-1)
  integer(i4b)    :: kAnisotropic          = imiss    ! anisotropy factor for lateral hydraulic conductivity (-)
  integer(i4b)    :: zScale_TOPMODEL       = imiss    ! TOPMODEL scaling factor used in lower boundary condition for soil (m)
  integer(i4b)    :: compactedDepth        = imiss    ! depth where k_soil reaches the compacted value given by CH78 (m)
  integer(i4b)    :: aquiferScaleFactor    = imiss    ! scaling factor for aquifer storage in the big bucket (m)
  integer(i4b)    :: aquiferBaseflowExp    = imiss    ! baseflow exponent (-)
  integer(i4b)    :: qSurfScale            = imiss    ! scaling factor in the surface runoff parameterization (-)
  integer(i4b)    :: specificYield         = imiss    ! specific yield (-)
  integer(i4b)    :: specificStorage       = imiss    ! specific storage coefficient (m-1)
  integer(i4b)    :: f_impede              = imiss    ! ice impedence factor (-)
  integer(i4b)    :: soilIceScale          = imiss    ! scaling factor for depth of soil ice, used to get frozen fraction (m)
  integer(i4b)    :: soilIceCV             = imiss    ! CV of depth of soil ice, used to get frozen fraction (-)
  ! algorithmic control parameters
  integer(i4b)    :: minwind               = imiss    ! minimum wind speed (m s-1)
  integer(i4b)    :: minstep               = imiss    ! minimum length of the time step
  integer(i4b)    :: maxstep               = imiss    ! maximum length of the time step
  integer(i4b)    :: wimplicit             = imiss    ! weight assigned to the start-of-step fluxes
  integer(i4b)    :: maxiter               = imiss    ! maximum number of iteration
  integer(i4b)    :: relConvTol_liquid     = imiss    ! relative convergence tolerance for vol frac liq water (-)
  integer(i4b)    :: absConvTol_liquid     = imiss    ! absolute convergence tolerance for vol frac liq water (-)
  integer(i4b)    :: relConvTol_matric     = imiss    ! relative convergence tolerance for matric head (-)
  integer(i4b)    :: absConvTol_matric     = imiss    ! absolute convergence tolerance for matric head (m)
  integer(i4b)    :: relConvTol_energy     = imiss    ! relative convergence tolerance for energy (-)
  integer(i4b)    :: absConvTol_energy     = imiss    ! absolute convergence tolerance for energy (J m-3)
  integer(i4b)    :: relConvTol_aquifr     = imiss    ! relative convergence tolerance for aquifer storage (-)
  integer(i4b)    :: absConvTol_aquifr     = imiss    ! absolute convergence tolerance for aquifer storage (J m-3)
  integer(i4b)    :: zmin                  = imiss    ! minimum layer depth (m)
  integer(i4b)    :: zmax                  = imiss    ! maximum layer depth (m)
  integer(i4b)    :: zminLayer1            = imiss    ! minimum layer depth for the 1st (top) layer (m)
  integer(i4b)    :: zminLayer2            = imiss    ! minimum layer depth for the 2nd layer (m)
  integer(i4b)    :: zminLayer3            = imiss    ! minimum layer depth for the 3rd layer (m)
  integer(i4b)    :: zminLayer4            = imiss    ! minimum layer depth for the 4th layer (m)
  integer(i4b)    :: zminLayer5            = imiss    ! minimum layer depth for the 5th (bottom) layer (m)
  integer(i4b)    :: zmaxLayer1_lower      = imiss    ! maximum layer depth for the 1st (top) layer when only 1 layer (m)
  integer(i4b)    :: zmaxLayer2_lower      = imiss    ! maximum layer depth for the 2nd layer when only 2 layers (m)
  integer(i4b)    :: zmaxLayer3_lower      = imiss    ! maximum layer depth for the 3rd layer when only 3 layers (m)
  integer(i4b)    :: zmaxLayer4_lower      = imiss    ! maximum layer depth for the 4th layer when only 4 layers (m)
  integer(i4b)    :: zmaxLayer1_upper      = imiss    ! maximum layer depth for the 1st (top) layer when > 1 layer (m)
  integer(i4b)    :: zmaxLayer2_upper      = imiss    ! maximum layer depth for the 2nd layer when > 2 layers (m)
  integer(i4b)    :: zmaxLayer3_upper      = imiss    ! maximum layer depth for the 3rd layer when > 3 layers (m)
  integer(i4b)    :: zmaxLayer4_upper      = imiss    ! maximum layer depth for the 4th layer when > 4 layers (m)
 endtype ilook_param


 ! ***********************************************************************************************************
 ! (6) define model prognostic (state) variables
 ! ***********************************************************************************************************
 type, public :: iLook_prog
  ! variables for time stepping
  integer(i4b)    :: dt_init                     = imiss    ! length of initial time step at start of next data interval (s)
  ! state variables for vegetation
  integer(i4b)    :: scalarCanopyIce             = imiss    ! mass of ice on the vegetation canopy (kg m-2)
  integer(i4b)    :: scalarCanopyLiq             = imiss    ! mass of liquid water on the vegetation canopy (kg m-2)
  integer(i4b)    :: scalarCanairTemp            = imiss    ! temperature of the canopy air space (Pa)
  integer(i4b)    :: scalarCanopyTemp            = imiss    ! temperature of the vegetation canopy (K)
  ! state variables for snow
  integer(i4b)    :: spectralSnowAlbedoDiffuse   = imiss    ! diffuse snow albedo for individual spectral bands (-)
  integer(i4b)    :: scalarSnowAlbedo            = imiss    ! snow albedo for the entire spectral band (-)
  integer(i4b)    :: scalarSnowDepth             = imiss    ! total snow depth (m)
  integer(i4b)    :: scalarSWE                   = imiss    ! snow water equivalent (kg m-2)
  integer(i4b)    :: scalarSfcMeltPond           = imiss    ! ponded water caused by melt of the "snow without a layer" (kg m-2)
  ! state variables for the snow+soil domain
  integer(i4b)    :: mLayerTemp                  = imiss    ! temperature of each layer (K)
  integer(i4b)    :: mLayerVolFracIce            = imiss    ! volumetric fraction of ice water in each layer (-)
  integer(i4b)    :: mLayerVolFracLiq            = imiss    ! volumetric fraction of liquid water in each layer (-)
  integer(i4b)    :: mLayerMatricHead            = imiss    ! matric head of water in the soil (m)
  ! other state variables
  integer(i4b)    :: scalarAquiferStorage        = imiss    ! relative aquifer storage -- above bottom of the soil profile (m)
  integer(i4b)    :: scalarSurfaceTemp           = imiss    ! surface temperature (K)
  ! coordinate variables
  integer(i4b)    :: mLayerDepth                 = imiss    ! depth of each layer (m)
  integer(i4b)    :: mLayerHeight                = imiss    ! height at the mid-point of each layer (m)
  integer(i4b)    :: iLayerHeight                = imiss    ! height of the layer interface; top of soil = 0 (m)
 endtype iLook_prog

 ! ***********************************************************************************************************
 ! (7) define diagnostic variables
 ! ***********************************************************************************************************
 type, public :: iLook_diag
  ! local properties
  integer(i4b)    :: scalarGreenVegFraction          = imiss ! green vegetation fraction used to compute LAI (-)
  integer(i4b)    :: scalarBulkVolHeatCapVeg         = imiss ! bulk volumetric heat capacity of vegetation (J m-3 K-1)
  integer(i4b)    :: scalarCanopyEmissivity          = imiss ! effective canopy emissivity (-)
  integer(i4b)    :: scalarRootZoneTemp              = imiss ! average temperature of the root zone (K)
  integer(i4b)    :: scalarLAI                       = imiss ! one-sided leaf area index (m2 m-2)
  integer(i4b)    :: scalarSAI                       = imiss ! one-sided stem area index (m2 m-2)
  integer(i4b)    :: scalarExposedLAI                = imiss ! exposed leaf area index after burial by snow (m2 m-2)
  integer(i4b)    :: scalarExposedSAI                = imiss ! exposed stem area index after burial by snow(m2 m-2)
  integer(i4b)    :: scalarCanopyIceMax              = imiss ! maximum interception storage capacity for ice (kg m-2)
  integer(i4b)    :: scalarCanopyLiqMax              = imiss ! maximum interception storage capacity for liquid water (kg m-2)
  integer(i4b)    :: scalarGrowingSeasonIndex        = imiss ! growing season index (0=off, 1=on)
  integer(i4b)    :: scalarVolHtCap_air              = imiss ! volumetric heat capacity air (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_ice              = imiss ! volumetric heat capacity ice (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_soil             = imiss ! volumetric heat capacity dry soil (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_water            = imiss ! volumetric heat capacity liquid wat (J m-3 K-1)
  integer(i4b)    :: mLayerVolHtCapBulk              = imiss ! volumetric heat capacity in each layer (J m-3 K-1)
  integer(i4b)    :: scalarLambda_drysoil            = imiss ! thermal conductivity of dry soil     (W m-1 K-1)
  integer(i4b)    :: scalarLambda_wetsoil            = imiss ! thermal conductivity of wet soil     (W m-1 K-1)
  integer(i4b)    :: mLayerThermalC                  = imiss ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
  integer(i4b)    :: iLayerThermalC                  = imiss ! thermal conductivity at the interface of each layer (W m-1 K-1)
  ! forcing
  integer(i4b)    :: scalarVPair                     = imiss ! vapor pressure of the air above the vegetation canopy (Pa)
  integer(i4b)    :: scalarVP_CanopyAir              = imiss ! vapor pressure of the canopy air space (Pa)
  integer(i4b)    :: scalarTwetbulb                  = imiss ! wet bulb temperature (K)
  integer(i4b)    :: scalarSnowfallTemp              = imiss ! temperature of fresh snow (K)
  integer(i4b)    :: scalarNewSnowDensity            = imiss ! density of fresh snow (kg m-3)
  integer(i4b)    :: scalarO2air                     = imiss ! atmospheric o2 concentration (Pa)
  integer(i4b)    :: scalarCO2air                    = imiss ! atmospheric co2 concentration (Pa)
  ! shortwave radiation
  integer(i4b)    :: scalarCosZenith                 = imiss ! cosine of the solar zenith angle (0-1)
  integer(i4b)    :: scalarFractionDirect            = imiss ! fraction of direct radiation (0-1)
  integer(i4b)    :: scalarCanopySunlitFraction      = imiss ! sunlit fraction of canopy (-)
  integer(i4b)    :: scalarCanopySunlitLAI           = imiss ! sunlit leaf area (-)
  integer(i4b)    :: scalarCanopyShadedLAI           = imiss ! shaded leaf area (-)
  integer(i4b)    :: spectralAlbGndDirect            = imiss ! direct  albedo of underlying surface for each spectral band (-)
  integer(i4b)    :: spectralAlbGndDiffuse           = imiss ! diffuse albedo of underlying surface for each spectral band (-)
  integer(i4b)    :: scalarGroundAlbedo              = imiss ! albedo of the ground surface (-)
  ! turbulent heat transfer
  integer(i4b)    :: scalarLatHeatSubVapCanopy       = imiss ! latent heat of sublimation/vaporization used for veg canopy (J kg-1)
  integer(i4b)    :: scalarLatHeatSubVapGround       = imiss ! latent heat of sublimation/vaporization used for ground surface (J kg-1)
  integer(i4b)    :: scalarSatVP_CanopyTemp          = imiss ! saturation vapor pressure at the temperature of vegetation canopy (Pa)
  integer(i4b)    :: scalarSatVP_GroundTemp          = imiss ! saturation vapor pressure at the temperature of the ground (Pa)
  integer(i4b)    :: scalarZ0Canopy                  = imiss ! roughness length of the canopy (m)
  integer(i4b)    :: scalarWindReductionFactor       = imiss ! canopy wind reduction factor (-)
  integer(i4b)    :: scalarZeroPlaneDisplacement     = imiss ! zero plane displacement (m)
  integer(i4b)    :: scalarRiBulkCanopy              = imiss ! bulk Richardson number for the canopy (-)
  integer(i4b)    :: scalarRiBulkGround              = imiss ! bulk Richardson number for the ground surface (-)
  integer(i4b)    :: scalarCanopyStabilityCorrection = imiss ! stability correction for the canopy (-)
  integer(i4b)    :: scalarGroundStabilityCorrection = imiss ! stability correction for the ground surface (-)
  ! evapotranspiration
  integer(i4b)    :: scalarIntercellularCO2Sunlit    = imiss ! carbon dioxide partial pressure of leaf interior (sunlit leaves) (Pa)
  integer(i4b)    :: scalarIntercellularCO2Shaded    = imiss ! carbon dioxide partial pressure of leaf interior (shaded leaves) (Pa)
  integer(i4b)    :: scalarTranspireLim              = imiss ! aggregate soil moisture + aquifer storage limit on transpiration (-)
  integer(i4b)    :: scalarTranspireLimAqfr          = imiss ! aquifer storage limit on transpiration (-)
  integer(i4b)    :: scalarFoliageNitrogenFactor     = imiss ! foliage nitrogen concentration, 1=saturated (-)
  integer(i4b)    :: scalarSoilRelHumidity           = imiss ! relative humidity in the soil pores in the upper-most soil layer (-)
  integer(i4b)    :: mLayerTranspireLim              = imiss ! soil moist & veg limit on transpiration for each layer (-)
  integer(i4b)    :: mLayerRootDensity               = imiss ! fraction of roots in each soil layer (-)
  integer(i4b)    :: scalarAquiferRootFrac           = imiss ! fraction of roots below the soil profile (-)
  ! canopy hydrology
  integer(i4b)    :: scalarFracLiqVeg                = imiss ! fraction of liquid water on vegetation (-)
  integer(i4b)    :: scalarCanopyWetFraction         = imiss ! fraction of canopy that is wet
  ! snow hydrology
  integer(i4b)    :: scalarSnowAge                   = imiss ! non-dimensional snow age (-)
  integer(i4b)    :: scalarGroundSnowFraction        = imiss ! fraction of ground that is covered with snow (-)
  integer(i4b)    :: spectralSnowAlbedoDirect        = imiss ! direct snow albedo for individual spectral bands (-)
  integer(i4b)    :: scalarFracLiqSnow               = imiss ! fraction of liquid water in each snow layer (-)
  integer(i4b)    :: mLayerThetaResid                = imiss ! residual volumetric water content in each snow layer (-)
  integer(i4b)    :: mLayerPoreSpace                 = imiss ! total pore space in each snow layer (-)
  integer(i4b)    :: mLayerMeltFreeze                = imiss ! change in ice content due to melt/freeze in each layer (kg m-3)
  ! soil hydrology
  integer(i4b)    :: scalarInfilArea                 = imiss ! fraction of unfrozen area where water can infiltrate (-)
  integer(i4b)    :: scalarFrozenArea                = imiss ! fraction of area that is considered impermeable due to soil ice (-)
  integer(i4b)    :: scalarSoilControl               = imiss ! soil control on infiltration: 1=controlling; 0=not (-)
  integer(i4b)    :: mLayerVolFracAir                = imiss ! volumetric fraction of air in each layer (-)
  integer(i4b)    :: mLayerTcrit                     = imiss ! critical soil temperature above which all water is unfrozen (K)
  integer(i4b)    :: mLayerCompress                  = imiss ! change in volumetric water content due to compression of soil (-)
  integer(i4b)    :: scalarSoilCompress              = imiss ! change in total soil storage due to compression of the soil matrix (kg m-2)
  ! mass balance check
  integer(i4b)    :: scalarSoilWatBalError           = imiss ! error in the total soil water balance (kg m-2)
  integer(i4b)    :: scalarAquiferBalError           = imiss ! error in the aquifer water balance (kg m-2)
  integer(i4b)    :: scalarTotalSoilLiq              = imiss ! total mass of liquid water in the soil (kg m-2)
  integer(i4b)    :: scalarTotalSoilIce              = imiss ! total mass of ice in the soil (kg m-2)
  ! variable shortcuts
  integer(i4b)    :: scalarVGn_m                     = imiss ! van Genuchten "m" parameter (-)
  integer(i4b)    :: scalarKappa                     = imiss ! constant in the freezing curve function (m K-1)
  integer(i4b)    :: scalarVolLatHt_fus              = imiss ! volumetric latent heat of fusion     (J m-3)
 endtype iLook_diag

 ! ***********************************************************************************************************
 ! (8) define model fluxes
 ! ***********************************************************************************************************
 type, public :: iLook_flux
  ! net energy and mass fluxes for the vegetation domain
  integer(i4b)    :: scalarCanairNetNrgFlux          = imiss ! net energy flux for the canopy air space (W m-2)
  integer(i4b)    :: scalarCanopyNetNrgFlux          = imiss ! net energy flux for the vegetation canopy (W m-2)
  integer(i4b)    :: scalarGroundNetNrgFlux          = imiss ! net energy flux for the ground surface (W m-2)
  integer(i4b)    :: scalarCanopyNetLiqFlux          = imiss ! net liquid water flux for the vegetation canopy (kg m-2 s-1)
  ! forcing
  integer(i4b)    :: scalarRainfall                  = imiss ! computed rainfall rate (kg m-2 s-1)
  integer(i4b)    :: scalarSnowfall                  = imiss ! computed snowfall rate (kg m-2 s-1)
  ! shortwave radiation
  integer(i4b)    :: spectralIncomingDirect          = imiss ! incoming direct solar radiation in each wave band (W m-2)
  integer(i4b)    :: spectralIncomingDiffuse         = imiss ! incoming diffuse solar radiation in each wave band (W m-2)
  integer(i4b)    :: scalarCanopySunlitPAR           = imiss ! average absorbed par for sunlit leaves (W m-2)
  integer(i4b)    :: scalarCanopyShadedPAR           = imiss ! average absorbed par for shaded leaves (W m-2)
  integer(i4b)    :: spectralBelowCanopyDirect       = imiss ! downward direct flux below veg layer for each spectral band  (W m-2)
  integer(i4b)    :: spectralBelowCanopyDiffuse      = imiss ! downward diffuse flux below veg layer for each spectral band (W m-2)
  integer(i4b)    :: scalarBelowCanopySolar          = imiss ! solar radiation transmitted below the canopy (W m-2)
  integer(i4b)    :: scalarCanopyAbsorbedSolar       = imiss ! solar radiation absorbed by canopy (W m-2)
  integer(i4b)    :: scalarGroundAbsorbedSolar       = imiss ! solar radiation absorbed by ground (W m-2)
  ! longwave radiation
  integer(i4b)    :: scalarLWRadCanopy               = imiss ! longwave radiation emitted from the canopy (W m-2)
  integer(i4b)    :: scalarLWRadGround               = imiss ! longwave radiation emitted at the ground surface  (W m-2)
  integer(i4b)    :: scalarLWRadUbound2Canopy        = imiss ! downward atmospheric longwave radiation absorbed by the canopy (W m-2)
  integer(i4b)    :: scalarLWRadUbound2Ground        = imiss ! downward atmospheric longwave radiation absorbed by the ground (W m-2)
  integer(i4b)    :: scalarLWRadUbound2Ubound        = imiss ! atmospheric radiation refl by ground + lost thru upper boundary (W m-2)
  integer(i4b)    :: scalarLWRadCanopy2Ubound        = imiss ! longwave radiation emitted from canopy lost thru upper boundary (W m-2)
  integer(i4b)    :: scalarLWRadCanopy2Ground        = imiss ! longwave radiation emitted from canopy absorbed by the ground (W m-2)
  integer(i4b)    :: scalarLWRadCanopy2Canopy        = imiss ! canopy longwave reflected from ground and absorbed by the canopy (W m-2)
  integer(i4b)    :: scalarLWRadGround2Ubound        = imiss ! longwave radiation emitted from ground lost thru upper boundary (W m-2)
  integer(i4b)    :: scalarLWRadGround2Canopy        = imiss ! longwave radiation emitted from ground and absorbed by the canopy (W m-2)
  integer(i4b)    :: scalarLWNetCanopy               = imiss ! net longwave radiation at the canopy (W m-2)
  integer(i4b)    :: scalarLWNetGround               = imiss ! net longwave radiation at the ground surface (W m-2)
  integer(i4b)    :: scalarLWNetUbound               = imiss ! net longwave radiation at the upper atmospheric boundary (W m-2)
  ! turbulent heat transfer
  integer(i4b)    :: scalarEddyDiffusCanopyTop       = imiss ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
  integer(i4b)    :: scalarFrictionVelocity          = imiss ! friction velocity - canopy momentum sink (m s-1)
  integer(i4b)    :: scalarWindspdCanopyTop          = imiss ! windspeed at the top of the canopy (m s-1)
  integer(i4b)    :: scalarWindspdCanopyBottom       = imiss ! windspeed at the height of the bottom of the canopy (m s-1)
  integer(i4b)    :: scalarGroundResistance          = imiss ! below canopy aerodynamic resistance (s m-1)
  integer(i4b)    :: scalarCanopyResistance          = imiss ! above canopy aerodynamic resistance (s m-1)
  integer(i4b)    :: scalarLeafResistance            = imiss ! mean leaf boundary layer resistance per unit leaf area (s m-1)
  integer(i4b)    :: scalarSoilResistance            = imiss ! soil surface resistance (s m-1)
  integer(i4b)    :: scalarSenHeatTotal              = imiss ! sensible heat from the canopy air space to the atmosphere (W m-2)
  integer(i4b)    :: scalarSenHeatCanopy             = imiss ! sensible heat from the canopy to the canopy air space (W m-2)
  integer(i4b)    :: scalarSenHeatGround             = imiss ! sensible heat from the ground (below canopy or non-vegetated) (W m-2)
  integer(i4b)    :: scalarLatHeatTotal              = imiss ! latent heat from the canopy air space to the atmosphere (W m-2)
  integer(i4b)    :: scalarLatHeatCanopyEvap         = imiss ! evaporation latent heat from the canopy to the canopy air space (W m-2)
  integer(i4b)    :: scalarLatHeatCanopyTrans        = imiss ! transpiration latent heat from the canopy to the canopy air space (W m-2)
  integer(i4b)    :: scalarLatHeatGround             = imiss ! latent heat from the ground (below canopy or non-vegetated) (W m-2)
  integer(i4b)    :: scalarCanopyAdvectiveHeatFlux   = imiss ! heat advected to the canopy surface with rain + snow (W m-2)
  integer(i4b)    :: scalarGroundAdvectiveHeatFlux   = imiss ! heat advected to the ground surface with throughfall and unloading/drainage (W m-2)
  integer(i4b)    :: scalarCanopySublimation         = imiss ! canopy sublimation/frost (kg m-2 s-1)
  integer(i4b)    :: scalarSnowSublimation           = imiss ! snow sublimation/frost (below canopy or non-vegetated) (kg m-2 s-1)
  ! liquid water fluxes associated with evapotranspiration
  integer(i4b)    :: scalarStomResistSunlit          = imiss ! stomatal resistance for sunlit leaves (s m-1)
  integer(i4b)    :: scalarStomResistShaded          = imiss ! stomatal resistance for shaded leaves (s m-1)
  integer(i4b)    :: scalarPhotosynthesisSunlit      = imiss ! sunlit photosynthesis (umolco2 m-2 s-1)
  integer(i4b)    :: scalarPhotosynthesisShaded      = imiss ! shaded photosynthesis (umolco2 m-2 s-1)
  integer(i4b)    :: scalarCanopyTranspiration       = imiss ! canopy transpiration (kg m-2 s-1)
  integer(i4b)    :: scalarCanopyEvaporation         = imiss ! canopy evaporation/condensation (kg m-2 s-1)
  integer(i4b)    :: scalarGroundEvaporation         = imiss ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
  integer(i4b)    :: mLayerTranspire                 = imiss ! transpiration loss from each soil layer (kg m-2 s-1)
  ! liquid and solid water fluxes through the canopy
  integer(i4b)    :: scalarThroughfallSnow           = imiss ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
  integer(i4b)    :: scalarThroughfallRain           = imiss ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
  integer(i4b)    :: scalarCanopySnowUnloading       = imiss ! unloading of snow from the vegetion canopy (kg m-2 s-1)
  integer(i4b)    :: scalarCanopyLiqDrainage         = imiss ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
  integer(i4b)    :: scalarCanopyMeltFreeze          = imiss ! melt/freeze of water stored in the canopy (kg m-2 s-1)
  ! energy fluxes and for the snow and soil domains
  integer(i4b)    :: iLayerConductiveFlux            = imiss ! conductive energy flux at layer interfaces (W m-2)
  integer(i4b)    :: iLayerAdvectiveFlux             = imiss ! advective energy flux at layer interfaces (W m-2)
  integer(i4b)    :: iLayerNrgFlux                   = imiss ! energy flux at layer interfaces (W m-2)
  integer(i4b)    :: mLayerNrgFlux                   = imiss ! net energy flux for each layer in the snow+soil domain (J m-3 s-1)
  ! liquid water fluxes for the snow domain
  integer(i4b)    :: iLayerLiqFluxSnow               = imiss ! liquid flux at snow layer interfaces (m s-1)
  integer(i4b)    :: mLayerLiqFluxSnow               = imiss ! net liquid water flux for each snow layer (s-1)
  ! liquid water fluxes for the soil domain
  integer(i4b)    :: scalarRainPlusMelt              = imiss ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  integer(i4b)    :: scalarMaxInfilRate              = imiss ! maximum infiltration rate (m s-1)
  integer(i4b)    :: scalarInfiltration              = imiss ! infiltration of water into the soil profile (m s-1)
  integer(i4b)    :: scalarExfiltration              = imiss ! exfiltration of water from the top of the soil profile (m s-1)
  integer(i4b)    :: scalarSurfaceRunoff             = imiss ! surface runoff (m s-1)
  integer(i4b)    :: mLayerSatHydCondMP              = imiss ! saturated hydraulic conductivity of macropores in each layer (m s-1)
  integer(i4b)    :: mLayerSatHydCond                = imiss ! saturated hydraulic conductivity in each layer (m s-1)
  integer(i4b)    :: iLayerSatHydCond                = imiss ! saturated hydraulic conductivity at each layer interface (m s-1)
  integer(i4b)    :: mLayerHydCond                   = imiss ! hydraulic conductivity in each soil layer (m s-1)
  integer(i4b)    :: iLayerLiqFluxSoil               = imiss ! liquid flux at soil layer interfaces (m s-1)
  integer(i4b)    :: mLayerLiqFluxSoil               = imiss ! net liquid water flux for each soil layer (s-1)
  integer(i4b)    :: mLayerBaseflow                  = imiss ! baseflow from each soil layer (m s-1)
  integer(i4b)    :: mLayerColumnInflow              = imiss ! total inflow to each layer in a given soil column (m3 s-1)
  integer(i4b)    :: mLayerColumnOutflow             = imiss ! total outflow from each layer in a given soil column (m3 s-1)
  integer(i4b)    :: scalarSoilBaseflow              = imiss ! total baseflow from throughout the soil profile (m s-1)
  integer(i4b)    :: scalarSoilDrainage              = imiss ! drainage from the bottom of the soil profile (m s-1)
  integer(i4b)    :: scalarAquiferRecharge           = imiss ! recharge to the aquifer (m s-1)
  integer(i4b)    :: scalarAquiferTranspire          = imiss ! transpiration from the aquifer (m s-1)
  integer(i4b)    :: scalarAquiferBaseflow           = imiss ! baseflow from the aquifer (m s-1)
 endtype iLook_flux

 ! ***********************************************************************************************************
 ! (9) define derivatives
 ! ***********************************************************************************************************
 type, public :: iLook_deriv
  ! derivatives in net vegetation energy fluxes w.r.t. relevant state variables
  integer(i4b)    :: dCanairNetFlux_dCanairTemp      = imiss ! derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
  integer(i4b)    :: dCanairNetFlux_dCanopyTemp      = imiss ! derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
  integer(i4b)    :: dCanairNetFlux_dGroundTemp      = imiss ! derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
  integer(i4b)    :: dCanopyNetFlux_dCanairTemp      = imiss ! derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
  integer(i4b)    :: dCanopyNetFlux_dCanopyTemp      = imiss ! derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
  integer(i4b)    :: dCanopyNetFlux_dGroundTemp      = imiss ! derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
  integer(i4b)    :: dCanopyNetFlux_dCanLiq          = imiss ! derivative in net canopy fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
  integer(i4b)    :: dGroundNetFlux_dCanairTemp      = imiss ! derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
  integer(i4b)    :: dGroundNetFlux_dCanopyTemp      = imiss ! derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
  integer(i4b)    :: dGroundNetFlux_dGroundTemp      = imiss ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
  integer(i4b)    :: dGroundNetFlux_dCanLiq          = imiss ! derivative in net ground fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
  ! derivatives in evaporative fluxes w.r.t. relevant state variables
  integer(i4b)    :: dCanopyEvaporation_dTCanair     = imiss ! derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
  integer(i4b)    :: dCanopyEvaporation_dTCanopy     = imiss ! derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
  integer(i4b)    :: dCanopyEvaporation_dTGround     = imiss ! derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
  integer(i4b)    :: dCanopyEvaporation_dCanLiq      = imiss ! derivative in canopy evaporation w.r.t. canopy liquid water content (s-1)
  integer(i4b)    :: dGroundEvaporation_dTCanair     = imiss ! derivative in ground evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
  integer(i4b)    :: dGroundEvaporation_dTCanopy     = imiss ! derivative in ground evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
  integer(i4b)    :: dGroundEvaporation_dTGround     = imiss ! derivative in ground evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
  integer(i4b)    :: dGroundEvaporation_dCanLiq      = imiss ! derivative in ground evaporation w.r.t. canopy liquid water content (s-1)
  ! derivatives in canopy water w.r.t canopy temperature
  integer(i4b)    :: dTheta_dTkCanopy                = imiss ! derivative of volumetric liquid water content w.r.t. temperature (K-1)
  integer(i4b)    :: dCanLiq_dTcanopy                = imiss ! derivative of canopy liquid storage w.r.t. temperature (kg m-2 K-1)
  ! derivatives in canopy liquid fluxes w.r.t. canopy water
  integer(i4b)    :: scalarCanopyLiqDeriv            = imiss ! derivative in (throughfall + canopy drainage) w.r.t. canopy liquid water (s-1)
  integer(i4b)    :: scalarThroughfallRainDeriv      = imiss ! derivative in throughfall w.r.t. canopy liquid water (s-1)
  integer(i4b)    :: scalarCanopyLiqDrainageDeriv    = imiss ! derivative in canopy drainage w.r.t. canopy liquid water (s-1)
  ! derivatives in energy fluxes at the interface of snow+soil layers w.r.t. temperature in layers above and below
  integer(i4b)    :: dNrgFlux_dTempAbove             = imiss ! derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
  integer(i4b)    :: dNrgFlux_dTempBelow             = imiss ! derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
  ! derivative in liquid water fluxes at the interface of snow layers w.r.t. volumetric liquid water content in the layer above
  integer(i4b)    :: iLayerLiqFluxSnowDeriv          = imiss ! derivative in vertical liquid water flux at layer interfaces (m s-1)
  ! derivative in liquid water fluxes for the soil domain w.r.t hydrology state variables
  integer(i4b)    :: dVolTot_dPsi0                   = imiss ! derivative in total water content w.r.t. total water matric potential (m-1)
  integer(i4b)    :: dq_dHydStateAbove               = imiss ! change in the flux in layer interfaces w.r.t. state variables in the layer above
  integer(i4b)    :: dq_dHydStateBelow               = imiss ! change in the flux in layer interfaces w.r.t. state variables in the layer below
  integer(i4b)    :: mLayerdTheta_dPsi               = imiss ! derivative in the soil water characteristic w.r.t. psi (m-1)
  integer(i4b)    :: mLayerdPsi_dTheta               = imiss ! derivative in the soil water characteristic w.r.t. theta (m)
  integer(i4b)    :: dCompress_dPsi                  = imiss ! derivative in compressibility w.r.t matric head (m-1)
  ! derivative in liquid water fluxes for the soil domain w.r.t energy state variables
  integer(i4b)    :: dq_dNrgStateAbove               = imiss ! change in the flux in layer interfaces w.r.t. state variables in the layer above
  integer(i4b)    :: dq_dNrgStateBelow               = imiss ! change in the flux in layer interfaces w.r.t. state variables in the layer below
  integer(i4b)    :: mLayerdTheta_dTk                = imiss ! derivative of volumetric liquid water content w.r.t. temperature (K-1)
  integer(i4b)    :: dPsiLiq_dTemp                   = imiss ! derivative in the liquid water matric potential w.r.t. temperature (m K-1)
 endtype iLook_deriv

 ! ***********************************************************************************************************
 ! (10) define model indices
 ! ***********************************************************************************************************
 type, public :: iLook_index
  ! number of state variables of different type
  integer(i4b)    :: nVegNrg                = imiss ! number of energy state variables for vegetation
  integer(i4b)    :: nVegMass               = imiss ! number of hydrology state variables for vegetation (mass of water)
  integer(i4b)    :: nVegState              = imiss ! number of vegetation state variables
  integer(i4b)    :: nNrgState              = imiss ! number of energy state variables
  integer(i4b)    :: nWatState              = imiss ! number of "total water" state variables (volumetric total water content)
  integer(i4b)    :: nMatState              = imiss ! number of matric head state variables
  integer(i4b)    :: nMassState             = imiss ! number of hydrology state variables (mass of water)
  integer(i4b)    :: nState                 = imiss ! total number of model state variables
  ! number of model layers and layer indices= imiss
  integer(i4b)    :: nSnow                  = imiss ! number of snow layers
  integer(i4b)    :: nSoil                  = imiss ! number of soil layers
  integer(i4b)    :: nLayers                = imiss ! total number of layers
  integer(i4b)    :: layerType              = imiss ! type of layer (soil or snow)
  ! indices of model state variables
  integer(i4b)    :: ixCasNrg               = imiss ! index of the canopy air space state variable
  integer(i4b)    :: ixVegNrg               = imiss ! index of the canopy energy state variable
  integer(i4b)    :: ixVegWat               = imiss ! index of the canopy total water state variable
  integer(i4b)    :: ixTopNrg               = imiss ! index of the upper-most energy state variable in the snow-soil subdomain
  integer(i4b)    :: ixTopWat               = imiss ! index of the upper-most total water state variable in the snow-soil subdomain
  integer(i4b)    :: ixTopMat               = imiss ! index of the upper-most matric head state variable in the soil subdomain
  integer(i4b)    :: ixSnowSoilNrg          = imiss ! indices for energy state variables in the snow-soil subdomain
  integer(i4b)    :: ixSnowSoilWat          = imiss ! indices for total water state variables in the snow-soil subdomain
  integer(i4b)    :: ixSnowOnlyNrg          = imiss ! indices for energy state variables in the snow subdomain
  integer(i4b)    :: ixSnowOnlyWat          = imiss ! indices for total water state variables in the snow subdomain
  integer(i4b)    :: ixSoilOnlyNrg          = imiss ! indices for energy state variables in the soil subdomain
  integer(i4b)    :: ixSoilOnlyHyd          = imiss ! indices for hydrology state variables in the soil subdomain
  ! type of model state variables
  integer(i4b)    :: ixStateType            = imiss ! indices defining the type of the state (ixNrgState, ixWatState, ixMatState...)
  integer(i4b)    :: ixAllState             = imiss ! list of indices for all model state variables
  integer(i4b)    :: ixSoilState            = imiss ! list of indices for all soil layers
  integer(i4b)    :: ixLayerState           = imiss ! list of indices for all model layers
  integer(i4b)    :: ixNrgOnly              = imiss ! list of indices for all energy states
  integer(i4b)    :: ixWatOnly              = imiss ! list of indices for all "total water" state variables (volumetric total water content)
  integer(i4b)    :: ixMatOnly              = imiss ! list of indices for matric head state variables
  integer(i4b)    :: ixMassOnly             = imiss ! list of indices for hydrology state variables (mass of water)
  ! indices for the model output files
  integer(i4b)    :: midSnowStartIndex      = imiss ! start index of the midSnow vector for a given timestep
  integer(i4b)    :: midSoilStartIndex      = imiss ! start index of the midSoil vector for a given timestep
  integer(i4b)    :: midTotoStartIndex      = imiss ! start index of the midToto vector for a given timestep
  integer(i4b)    :: ifcSnowStartIndex      = imiss ! start index of the ifcSnow vector for a given timestep
  integer(i4b)    :: ifcSoilStartIndex      = imiss ! start index of the ifcSoil vector for a given timestep
  integer(i4b)    :: ifcTotoStartIndex      = imiss ! start index of the ifcToto vector for a given timestep
 endtype iLook_index

 ! ***********************************************************************************************************
 ! (11) define basin-average model parameters
 ! ***********************************************************************************************************
 type, public :: iLook_bpar
  ! baseflow
  integer(i4b)    :: basin__aquiferHydCond      = imiss ! hydraulic conductivity for the aquifer (m s-1)
  integer(i4b)    :: basin__aquiferScaleFactor  = imiss ! scaling factor for aquifer storage in the big bucket (m)
  integer(i4b)    :: basin__aquiferBaseflowExp  = imiss ! baseflow exponent for the big bucket (-)
  ! within-grid routing
  integer(i4b)    :: routingGammaShape          = imiss ! shape parameter in Gamma distribution used for sub-grid routing (-)
  integer(i4b)    :: routingGammaScale          = imiss ! scale parameter in Gamma distribution used for sub-grid routing (s)
 endtype iLook_bpar

 ! ***********************************************************************************************************
 ! (12) define basin-average model variables
 ! ***********************************************************************************************************
 type, public :: iLook_bvar
  ! define derived variables
  integer(i4b)    :: basin__totalArea           = imiss ! total basin area (m2)
  ! define fluxes
  integer(i4b)    :: basin__SurfaceRunoff       = imiss ! surface runoff (m s-1)
  integer(i4b)    :: basin__ColumnOutflow       = imiss ! outflow from all "outlet" HRUs (those with no downstream HRU)
  integer(i4b)    :: basin__AquiferStorage      = imiss ! aquifer storage (m s-1)
  integer(i4b)    :: basin__AquiferRecharge     = imiss ! recharge to the aquifer (m s-1)
  integer(i4b)    :: basin__AquiferBaseflow     = imiss ! baseflow from the aquifer (m s-1)
  integer(i4b)    :: basin__AquiferTranspire    = imiss ! transpiration from the aquifer (m s-1)
  ! define variables for runoff
  integer(i4b)    :: routingRunoffFuture        = imiss ! runoff in future time steps (m s-1)
  integer(i4b)    :: routingFractionFuture      = imiss ! fraction of runoff in future time steps (-)
  integer(i4b)    :: averageInstantRunoff       = imiss ! instantaneous runoff (m s-1)
  integer(i4b)    :: averageRoutedRunoff        = imiss ! routed runoff (m s-1)
 endtype iLook_bvar

 ! ***********************************************************************************************************
 ! (X) define data structures and maximum number of variables of each type
 ! ***********************************************************************************************************

 ! named variables: model decisions
 type(iLook_decision),public,parameter :: iLookDECISIONS=iLook_decision(  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,&
                                                                         21, 22, 23, 24, 25, 26, 27, 28, 29, 30,&
                                                                         31, 32, 33, 34, 35, 36, 37, 38, 39)

 ! named variables: model time
 type(iLook_time),    public,parameter :: iLookTIME     =iLook_time    (  1,  2,  3,  4,  5)

 ! named variables: model forcing data
 type(iLook_force),   public,parameter :: iLookFORCE    =iLook_force   (  1,  2,  3,  4,  5,  6,  7,  8)

 ! named variables: model attributes
 type(iLook_attr),    public,parameter :: iLookATTR     =iLook_attr    (  1,  2,  3,  4,  5,  6,  7)

 ! named variables: soil and vegetation types
 type(iLook_type),    public,parameter :: iLookTYPE     =iLook_type    (  1,  2,  3,  4,  5)

 ! named variables: model parameters
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
                                                                        101,102,103,104,105,106,107,108,109,110,&
                                                                        111,112,113,114,115,116,117,118,119,120,&
                                                                        121,122,123,124,125,126,127,128,129,130,&
                                                                        131,132,133,134,135,136,137,138,139,140,&
                                                                        141,142,143,144,145,146,147,148,149,150,&
                                                                        151,152,153)

 ! named variables: model prognostic (state) variables
 type(iLook_prog),   public,parameter  :: iLookPROG     =iLook_prog    (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12, 13, 14, 15, 16, 17, 18, 19)
 ! named variables: model diagnostic variables
 type(iLook_diag),    public,parameter :: iLookDIAG     =iLook_diag    (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,&
                                                                         21, 22, 23, 24, 25, 26, 27, 28, 29, 30,&
                                                                         31, 32, 33, 34, 35, 36, 37, 38, 39, 40,&
                                                                         41, 42, 43, 44, 45, 46, 47, 48, 49, 50,&
                                                                         51, 52, 53, 54, 55, 56, 57, 58, 59, 60,&
                                                                         61, 62, 63, 64, 65, 66, 67, 68, 69, 70,&
                                                                         71, 72, 73, 74, 75, 76, 77, 78)
 ! named variables: model fluxes
 type(iLook_flux),    public,parameter :: iLookFLUX     =iLook_flux    (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,&
                                                                         21, 22, 23, 24, 25, 26, 27, 28, 29, 30,&
                                                                         31, 32, 33, 34, 35, 36, 37, 38, 39, 40,&
                                                                         41, 42, 43, 44, 45, 46, 47, 48, 49, 50,&
                                                                         51, 52, 53, 54, 55, 56, 57, 58, 59, 60,&
                                                                         61, 62, 63, 64, 65, 66, 67, 68, 69, 70,&
                                                                         71, 72, 73, 74, 75, 76, 77, 78, 79, 80,&
                                                                         81, 82, 83, 84, 85)

 ! named variables: derivatives in model fluxes w.r.t. relevant state variables
 type(iLook_deriv),   public,parameter :: iLookDERIV    =iLook_deriv   (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,&
                                                                         21, 22, 23, 24, 25, 26, 27, 28, 29, 30,&
                                                                         31, 32, 33, 34, 35, 36, 37)

 type(iLook_index),   public,parameter :: iLookINDEX    =ilook_index   (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,&
                                                                         21, 22, 23, 24, 25, 26, 27, 28, 29, 30,&
                                                                         31, 32, 33, 34, 35, 36, 37, 38)

 ! named variables: basin-average parameters
 type(iLook_bpar),    public,parameter :: iLookBPAR     =ilook_bpar    (  1,  2,  3,  4,  5)

 ! named variables: basin-average variables
 type(iLook_bvar),    public,parameter :: iLookBVAR     =ilook_bvar    (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11)
 ! define maximum number of variables of each type
 integer(i4b),parameter,public :: maxvarDecisions = storage_size(iLookDECISIONS)/iLength
 integer(i4b),parameter,public :: maxvarTime      = storage_size(iLookTIME)/iLength
 integer(i4b),parameter,public :: maxvarForc      = storage_size(iLookFORCE)/iLength
 integer(i4b),parameter,public :: maxvarAttr      = storage_size(iLookATTR)/iLength
 integer(i4b),parameter,public :: maxvarType      = storage_size(iLookTYPE)/iLength
 integer(i4b),parameter,public :: maxvarMpar      = storage_size(iLookPARAM)/iLength
 integer(i4b),parameter,public :: maxvarProg      = storage_size(iLookPROG)/iLength
 integer(i4b),parameter,public :: maxvarDiag      = storage_size(iLookDIAG)/iLength
 integer(i4b),parameter,public :: maxvarFlux      = storage_size(iLookFLUX)/iLength
 integer(i4b),parameter,public :: maxvarDeriv     = storage_size(iLookDERIV)/iLength
 integer(i4b),parameter,public :: maxvarIndx      = storage_size(iLookINDEX)/iLength
 integer(i4b),parameter,public :: maxvarBpar      = storage_size(iLookBPAR)/iLength
 integer(i4b),parameter,public :: maxvarBvar      = storage_size(iLookBVAR)/iLength

 ! ***********************************************************************************************************
 ! (Y) define ancillary look-up structures
 ! ***********************************************************************************************************

 integer(i4b), public          :: childFLUX_MEAN(maxvarFlux)  ! index of the child data structure: mean flux


END MODULE var_lookup

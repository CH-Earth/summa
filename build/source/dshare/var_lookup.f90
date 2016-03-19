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
  integer(i4b)    :: simulStart       ! simulation start time
  integer(i4b)    :: simulFinsh       ! simulation end time
  integer(i4b)    :: soilCatTbl       ! soil-category dateset
  integer(i4b)    :: vegeParTbl       ! vegetation category dataset
  integer(i4b)    :: soilStress       ! choice of function for the soil moisture control on stomatal resistance
  integer(i4b)    :: stomResist       ! choice of function for stomatal resistance
  integer(i4b)    :: bbTempFunc       ! Ball-Berry: leaf temperature controls on photosynthesis + stomatal resistance
  integer(i4b)    :: bbHumdFunc       ! Ball-Berry: humidity controls on stomatal resistance
  integer(i4b)    :: bbElecFunc       ! Ball-Berry: dependence of photosynthesis on PAR
  integer(i4b)    :: bbCO2point       ! Ball-Berry: use of CO2 compensation point to calculate stomatal resistance
  integer(i4b)    :: bbNumerics       ! Ball-Berry: iterative numerical solution method
  integer(i4b)    :: bbAssimFnc       ! Ball-Berry: controls on carbon assimilation
  integer(i4b)    :: bbCanIntg8       ! Ball-Berry: scaling of photosynthesis from the leaf to the canopy
  integer(i4b)    :: num_method       ! choice of numerical method
  integer(i4b)    :: fDerivMeth       ! method used to calculate flux derivatives
  integer(i4b)    :: LAI_method       ! method used to determine LAI and SAI
  integer(i4b)    :: cIntercept       ! choice of parameterization for canopy interception
  integer(i4b)    :: f_Richards       ! form of richards' equation
  integer(i4b)    :: groundwatr       ! choice of groundwater parameterization
  integer(i4b)    :: hc_profile       ! choice of hydraulic conductivity profile
  integer(i4b)    :: bcUpprTdyn       ! type of upper boundary condition for thermodynamics
  integer(i4b)    :: bcLowrTdyn       ! type of lower boundary condition for thermodynamics
  integer(i4b)    :: bcUpprSoiH       ! type of upper boundary condition for soil hydrology
  integer(i4b)    :: bcLowrSoiH       ! type of lower boundary condition for soil hydrology
  integer(i4b)    :: veg_traits       ! choice of parameterization for vegetation roughness length and displacement height
  integer(i4b)    :: rootProfil       ! choice of parameterization for the rooting profile
  integer(i4b)    :: canopyEmis       ! choice of parameterization for canopy emissivity
  integer(i4b)    :: snowIncept       ! choice of parameterization for snow interception
  integer(i4b)    :: windPrfile       ! choice of canopy wind profile
  integer(i4b)    :: astability       ! choice of stability function
  integer(i4b)    :: canopySrad       ! choice of method for canopy shortwave radiation
  integer(i4b)    :: alb_method       ! choice of albedo representation
  integer(i4b)    :: snowLayers       ! choice of method to combine and sub-divide snow layers
  integer(i4b)    :: compaction       ! choice of compaction routine
  integer(i4b)    :: thCondSnow       ! choice of thermal conductivity representation for snow
  integer(i4b)    :: thCondSoil       ! choice of thermal conductivity representation for soil
  integer(i4b)    :: spatial_gw       ! choice of method for spatial representation of groundwater
  integer(i4b)    :: subRouting       ! choice of method for sub-grid routing
 endtype iLook_decision

 ! ***********************************************************************************************************
 ! (1) define model time
 ! ***********************************************************************************************************
 type, public  ::  iLook_time
  integer(i4b)    :: iyyy             ! year
  integer(i4b)    :: im               ! month
  integer(i4b)    :: id               ! day
  integer(i4b)    :: ih               ! hour
  integer(i4b)    :: imin             ! minute
 endtype iLook_time

 ! ***********************************************************************************************************
 ! (2) define model forcing data
 ! ***********************************************************************************************************
 type, public  ::  iLook_force
  integer(i4b)    :: time             ! time since time reference       (s)
  integer(i4b)    :: pptrate          ! precipitation rate              (kg m-2 s-1)
  integer(i4b)    :: airtemp          ! air temperature                 (K)
  integer(i4b)    :: spechum          ! specific humidity               (g/g)
  integer(i4b)    :: windspd          ! windspeed                       (m/s)
  integer(i4b)    :: SWRadAtm         ! downwelling shortwave radiaiton (W m-2)
  integer(i4b)    :: LWRadAtm         ! downwelling longwave radiation  (W m-2)
  integer(i4b)    :: airpres          ! pressure                        (Pa)
 endtype iLook_force

 ! ***********************************************************************************************************
 ! (3) define local attributes
 ! ***********************************************************************************************************
 type, public  ::  iLook_attr
  integer(i4b)    :: latitude         ! latitude (degrees north)
  integer(i4b)    :: longitude        ! longitude (degrees east)
  integer(i4b)    :: elevation        ! elevation (m)
  integer(i4b)    :: tan_slope        ! tan water table slope, taken as tan local ground surface slope (-)
  integer(i4b)    :: contourLength    ! length of contour at downslope edge of HRU (m)
  integer(i4b)    :: HRUarea          ! area of each HRU  (m2)
  integer(i4b)    :: mHeight          ! measurement height above bare ground (m)
 end type iLook_attr

 ! ***********************************************************************************************************
 ! (4) define local classification of veg, soil, etc.
 ! ***********************************************************************************************************
 type, public  ::  iLook_type
  integer(i4b)    :: hruIndex         ! index defining hydrologic response unit (-)
  integer(i4b)    :: vegTypeIndex     ! index defining vegetation type (-)
  integer(i4b)    :: soilTypeIndex    ! index defining soil type (-)
  integer(i4b)    :: slopeTypeIndex   ! index defining slope (-)
  integer(i4b)    :: downHRUindex     ! index of downslope HRU (0 = basin outlet)
 end type iLook_type

 ! ***********************************************************************************************************
 ! (5) define model parameters
 ! ***********************************************************************************************************
 type, public  ::  iLook_param
  ! boundary conditions
  integer(i4b)    :: upperBoundHead             ! matric head of the upper boundary (m)
  integer(i4b)    :: lowerBoundHead             ! matric head of the lower boundary (m)
  integer(i4b)    :: upperBoundTheta            ! volumetric liquid water content of the upper boundary (-)
  integer(i4b)    :: lowerBoundTheta            ! volumetric liquid water content of the lower boundary (-)
  integer(i4b)    :: upperBoundTemp             ! temperature of the upper boundary (K)
  integer(i4b)    :: lowerBoundTemp             ! temperature of the lower boundary (K)
  ! precipitation partitioning
  integer(i4b)    :: tempCritRain               ! critical temperature where precipitation is rain (K)
  integer(i4b)    :: tempRangeTimestep          ! temperature range over the time step (K)
  integer(i4b)    :: frozenPrecipMultip         ! frozen precipitation multiplier (-)
  ! freezing curve for snow
  integer(i4b)    :: snowfrz_scale              ! scaling parameter for the freezing curve for snow (K-1)
  ! snow albedo
  integer(i4b)    :: albedoMax                  ! maximum snow albedo for a single spectral band (-)
  integer(i4b)    :: albedoMinWinter            ! minimum snow albedo during winter for a single spectral band (-)
  integer(i4b)    :: albedoMinSpring            ! minimum snow albedo during spring for a single spectral band (-)
  integer(i4b)    :: albedoMaxVisible           ! maximum snow albedo in the visible part of the spectrum (-)
  integer(i4b)    :: albedoMinVisible           ! minimum snow albedo in the visible part of the spectrum (-)
  integer(i4b)    :: albedoMaxNearIR            ! maximum snow albedo in the near infra-red part of the spectrum (-)
  integer(i4b)    :: albedoMinNearIR            ! minimum snow albedo in the near infra-red part of the spectrum (-)
  integer(i4b)    :: albedoDecayRate            ! albedo decay rate (s)
  integer(i4b)    :: albedoSootLoad             ! soot load factor (-)
  integer(i4b)    :: albedoRefresh              ! critical mass necessary for albedo refreshment (kg m-2)
  ! radiation transfer within snow
  integer(i4b)    :: radExt_snow                ! extinction coefficient for radiation penetration into the snowpack (m-1)
  integer(i4b)    :: directScale                ! scaling factor for fractional driect radiaion parameterization (-)
  integer(i4b)    :: Frad_direct                ! maximum fraction of direct solar radiation (-)
  integer(i4b)    :: Frad_vis                   ! fraction of radiation in the visible part of the spectrum (-)
  ! new snow density
  integer(i4b)    :: newSnowDenMin              ! minimum new snow density (kg m-3)
  integer(i4b)    :: newSnowDenMult             ! multiplier for new snow density (kg m-3)
  integer(i4b)    :: newSnowDenScal             ! scaling factor for new snow density (K)
  ! snow compaction
  integer(i4b)    :: densScalGrowth             ! density scaling factor for grain growth (kg-1 m3)
  integer(i4b)    :: tempScalGrowth             ! temperature scaling factor for grain growth (K-1)
  integer(i4b)    :: grainGrowthRate            ! rate of grain growth (s-1)
  integer(i4b)    :: densScalOvrbdn             ! density scaling factor for overburden pressure (kg-1 m3)
  integer(i4b)    :: tempScalOvrbdn             ! temperature scaling factor for overburden pressure (K-1)
  integer(i4b)    :: base_visc                  ! viscosity coefficient at T=T_frz and snow density=0  (kg s m-2)
  ! water flow within snow
  integer(i4b)    :: Fcapil                     ! capillary retention as a fraction of the total pore volume (-)
  integer(i4b)    :: k_snow                     ! hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
  integer(i4b)    :: mw_exp                     ! exponent for meltwater flow (-)
  ! turbulent heat fluxes
  integer(i4b)    :: z0Snow                     ! roughness length of snow (m)
  integer(i4b)    :: z0Soil                     ! roughness length of bare soil below the canopy (m)
  integer(i4b)    :: z0Canopy                   ! roughness length of the canopy (m)
  integer(i4b)    :: zpdFraction                ! zero plane displacement / canopy height (-)
  integer(i4b)    :: critRichNumber             ! critical value for the bulk Richardson number (-)
  integer(i4b)    :: Louis79_bparam             ! parameter in Louis (1979) stability function (-)
  integer(i4b)    :: Louis79_cStar              ! parameter in Louis (1979) stability function (-)
  integer(i4b)    :: Mahrt87_eScale             ! exponential scaling factor in the Mahrt (1987) stability function (-)
  integer(i4b)    :: leafExchangeCoeff          ! turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
  integer(i4b)    :: windReductionParam         ! canopy wind reduction parameter (-)
  ! stomatal conductance
  integer(i4b)    :: Kc25                       ! Michaelis-Menten constant for CO2 at 25 degrees C (umol mol-1)
  integer(i4b)    :: Ko25                       ! Michaelis-Menten constant for O2 at 25 degrees C (mol mol-1)
  integer(i4b)    :: Kc_qFac                    ! factor in the q10 function defining temperature controls on Kc (-)
  integer(i4b)    :: Ko_qFac                    ! factor in the q10 function defining temperature controls on Ko (-)
  integer(i4b)    :: kc_Ha                      ! activation energy for the Michaelis-Menten constant for CO2 (J mol-1)
  integer(i4b)    :: ko_Ha                      ! activation energy for the Michaelis-Menten constant for O2 (J mol-1)
  integer(i4b)    :: vcmax25_canopyTop          ! potential carboxylation rate at 25 degrees C at the canopy top (umol co2 m-2 s-1)
  integer(i4b)    :: vcmax_qFac                 ! factor in the q10 function defining temperature controls on vcmax (-)
  integer(i4b)    :: vcmax_Ha                   ! activation energy in the vcmax function (J mol-1)
  integer(i4b)    :: vcmax_Hd                   ! deactivation energy in the vcmax function (J mol-1)
  integer(i4b)    :: vcmax_Sv                   ! entropy term in the vcmax function (J mol-1 K-1)
  integer(i4b)    :: vcmax_Kn                   ! foliage nitrogen decay coefficient (-)
  integer(i4b)    :: jmax25_scale               ! scaling factor to relate jmax25 to vcmax25 (-)
  integer(i4b)    :: jmax_Ha                    ! activation energy in the jmax function (J mol-1)
  integer(i4b)    :: jmax_Hd                    ! deactivation energy in the jmax function (J mol-1)
  integer(i4b)    :: jmax_Sv                    ! entropy term in the jmax function (J mol-1 K-1)
  integer(i4b)    :: fractionJ                  ! fraction of light lost by other than the chloroplast lamellae (-)
  integer(i4b)    :: quantamYield               ! quantam yield (mol e mol-1 q)
  integer(i4b)    :: vpScaleFactor              ! vapor pressure scaling factor in stomatal conductance function (Pa)
  integer(i4b)    :: cond2photo_slope           ! slope of conductance-photosynthesis relationship (-)
  integer(i4b)    :: minStomatalConductance     ! minimum stomatal conductance (umol H2O m-2 s-1)
  ! vegetation properties
  integer(i4b)    :: winterSAI                  ! stem area index prior to the start of the growing season (m2 m-2)
  integer(i4b)    :: summerLAI                  ! maximum leaf area index at the peak of the growing season (m2 m-2)
  integer(i4b)    :: rootScaleFactor1           ! 1st scaling factor (a) in Y = 1 - 0.5*( exp(-aZ) + exp(-bZ) )  (m-1)
  integer(i4b)    :: rootScaleFactor2           ! 2nd scaling factor (b) in Y = 1 - 0.5*( exp(-aZ) + exp(-bZ) )  (m-1)
  integer(i4b)    :: rootingDepth               ! rooting depth (m)
  integer(i4b)    :: rootDistExp                ! exponent controlling the vertical distribution of root density (-)
  integer(i4b)    :: plantWiltPsi               ! matric head at wilting point (m)
  integer(i4b)    :: soilStressParam            ! parameter in the exponential soil stress function
  integer(i4b)    :: critSoilWilting            ! critical vol. liq. water content when plants are wilting (-)
  integer(i4b)    :: critSoilTranspire          ! critical vol. liq. water content when transpiration is limited (-)
  integer(i4b)    :: critAquiferTranspire       ! critical aquifer storage value when transpiration is limited (m)
  integer(i4b)    :: minStomatalResistance      ! minimum canopy resistance (s m-1)
  integer(i4b)    :: leafDimension              ! characteristic leaf dimension (m)
  integer(i4b)    :: heightCanopyTop            ! height of top of the vegetation canopy above ground surface (m)
  integer(i4b)    :: heightCanopyBottom         ! height of bottom of the vegetation canopy above ground surface (m)
  integer(i4b)    :: specificHeatVeg            ! specific heat of vegetation (J kg-1 K-1)
  integer(i4b)    :: maxMassVegetation          ! maximum mass of vegetation (full foliage) (kg m-2)
  integer(i4b)    :: throughfallScaleSnow       ! scaling factor for throughfall (snow) (-)
  integer(i4b)    :: throughfallScaleRain       ! scaling factor for throughfall (rain) (-)
  integer(i4b)    :: refInterceptCapSnow        ! reference canopy interception capacity per unit leaf area (snow) (kg m-2)
  integer(i4b)    :: refInterceptCapRain        ! canopy interception capacity per unit leaf area (rain) (kg m-2)
  integer(i4b)    :: snowUnloadingCoeff         ! time constant for unloading of snow from the forest canopy (s-1)
  integer(i4b)    :: canopyDrainageCoeff        ! time constant for drainage of liquid water from the forest canopy (s-1)
  integer(i4b)    :: ratioDrip2Unloading        ! ratio of canopy drip to unloading of snow from the forest canopy (-)
  integer(i4b)    :: canopyWettingFactor        ! maximum wetted fraction of the canopy (-)
  integer(i4b)    :: canopyWettingExp           ! exponent in canopy wetting function (-)
  ! soil properties
  integer(i4b)    :: soil_dens_intr             ! intrinsic soil density (kg m-3)
  integer(i4b)    :: thCond_soil                ! thermal conductivity of soil (W m-1 K-1)
  integer(i4b)    :: frac_sand                  ! fraction of sand (-)
  integer(i4b)    :: frac_silt                  ! fraction of silt (-)
  integer(i4b)    :: frac_clay                  ! fraction of clay (-)
  integer(i4b)    :: fieldCapacity              ! field capacity (-)
  integer(i4b)    :: wettingFrontSuction        ! Green-Ampt wetting front suction (m)
  integer(i4b)    :: theta_mp                   ! volumetric liquid water content when macropore flow begins (-)
  integer(i4b)    :: theta_sat                  ! porosity (-)
  integer(i4b)    :: theta_res                  ! volumetric residual water content (-)
  integer(i4b)    :: vGn_alpha                  ! van Genuchten "alpha" parameter (m-1)
  integer(i4b)    :: vGn_n                      ! van Genuchten "n" parameter (-)
  integer(i4b)    :: mpExp                      ! empirical exponent in macropore flow equation (-)
  integer(i4b)    :: k_soil                     ! hydraulic conductivity of soil (m s-1)
  integer(i4b)    :: k_macropore                ! saturated hydraulic conductivity for macropores (m s-1)
  integer(i4b)    :: kAnisotropic               ! anisotropy factor for lateral hydraulic conductivity (-)
  integer(i4b)    :: zScale_TOPMODEL            ! TOPMODEL scaling factor used in lower boundary condition for soil (m)
  integer(i4b)    :: compactedDepth             ! depth where k_soil reaches the compacted value given by CH78 (m)
  integer(i4b)    :: aquiferScaleFactor         ! scaling factor for aquifer storage in the big bucket (m)
  integer(i4b)    :: aquiferBaseflowExp         ! baseflow exponent (-)
  integer(i4b)    :: qSurfScale                 ! scaling factor in the surface runoff parameterization (-)
  integer(i4b)    :: specificYield              ! specific yield (-)
  integer(i4b)    :: specificStorage            ! specific storage coefficient (m-1)
  integer(i4b)    :: f_impede                   ! ice impedence factor (-)
  integer(i4b)    :: soilIceScale               ! scaling factor for depth of soil ice, used to get frozen fraction (m)
  integer(i4b)    :: soilIceCV                  ! CV of depth of soil ice, used to get frozen fraction (-)
  ! algorithmic control parameters
  integer(i4b)    :: minwind                    ! minimum wind speed (m s-1)
  integer(i4b)    :: minstep                    ! minimum length of the time step
  integer(i4b)    :: maxstep                    ! maximum length of the time step
  integer(i4b)    :: wimplicit                  ! weight assigned to the start-of-step fluxes
  integer(i4b)    :: maxiter                    ! maximum number of iteration
  integer(i4b)    :: relConvTol_liquid          ! relative convergence tolerance for vol frac liq water (-)
  integer(i4b)    :: absConvTol_liquid          ! absolute convergence tolerance for vol frac liq water (-)
  integer(i4b)    :: relConvTol_matric          ! relative convergence tolerance for matric head (-)
  integer(i4b)    :: absConvTol_matric          ! absolute convergence tolerance for matric head (m)
  integer(i4b)    :: relConvTol_energy          ! relative convergence tolerance for energy (-)
  integer(i4b)    :: absConvTol_energy          ! absolute convergence tolerance for energy (J m-3)
  integer(i4b)    :: relConvTol_aquifr          ! relative convergence tolerance for aquifer storage (-)
  integer(i4b)    :: absConvTol_aquifr          ! absolute convergence tolerance for aquifer storage (J m-3)
  integer(i4b)    :: zmin                       ! minimum layer depth (m)
  integer(i4b)    :: zmax                       ! maximum layer depth (m)
  integer(i4b)    :: zminLayer1                 ! minimum layer depth for the 1st (top) layer (m)
  integer(i4b)    :: zminLayer2                 ! minimum layer depth for the 2nd layer (m)
  integer(i4b)    :: zminLayer3                 ! minimum layer depth for the 3rd layer (m)
  integer(i4b)    :: zminLayer4                 ! minimum layer depth for the 4th layer (m)
  integer(i4b)    :: zminLayer5                 ! minimum layer depth for the 5th (bottom) layer (m)
  integer(i4b)    :: zmaxLayer1_lower           ! maximum layer depth for the 1st (top) layer when only 1 layer (m)
  integer(i4b)    :: zmaxLayer2_lower           ! maximum layer depth for the 2nd layer when only 2 layers (m)
  integer(i4b)    :: zmaxLayer3_lower           ! maximum layer depth for the 3rd layer when only 3 layers (m)
  integer(i4b)    :: zmaxLayer4_lower           ! maximum layer depth for the 4th layer when only 4 layers (m)
  integer(i4b)    :: zmaxLayer1_upper           ! maximum layer depth for the 1st (top) layer when > 1 layer (m)
  integer(i4b)    :: zmaxLayer2_upper           ! maximum layer depth for the 2nd layer when > 2 layers (m)
  integer(i4b)    :: zmaxLayer3_upper           ! maximum layer depth for the 3rd layer when > 3 layers (m)
  integer(i4b)    :: zmaxLayer4_upper           ! maximum layer depth for the 4th layer when > 4 layers (m)
 endtype ilook_param


 ! ***********************************************************************************************************
 ! (6) define model prognostic (state) variables
 ! ***********************************************************************************************************
 type, public :: iLook_prog
  ! define variables for time stepping
  integer(i4b)    :: dt_init                          ! length of initial time step at start of next data interval (s)
  ! define state variables for vegetation
  integer(i4b)    :: scalarCanopyIce                  ! mass of ice on the vegetation canopy (kg m-2)
  integer(i4b)    :: scalarCanopyLiq                  ! mass of liquid water on the vegetation canopy (kg m-2)
  integer(i4b)    :: scalarCanairTemp                 ! temperature of the canopy air space (Pa)
  integer(i4b)    :: scalarCanopyTemp                 ! temperature of the vegetation canopy (K)
  ! define state variables for snow
  integer(i4b)    :: spectralSnowAlbedoDiffuse        ! diffuse snow albedo for individual spectral bands (-)
  integer(i4b)    :: scalarSnowAlbedo                 ! snow albedo for the entire spectral band (-)
  integer(i4b)    :: scalarSnowDepth                  ! total snow depth (m)
  integer(i4b)    :: scalarSWE                        ! snow water equivalent (kg m-2)
  integer(i4b)    :: scalarSfcMeltPond                ! ponded water caused by melt of the "snow without a layer" (kg m-2)
  ! define state variables for the snow+soil domain
  integer(i4b)    :: mLayerTemp                       ! temperature of each layer (K)
  integer(i4b)    :: mLayerVolFracIce                 ! volumetric fraction of ice water in each layer (-)
  integer(i4b)    :: mLayerVolFracLiq                 ! volumetric fraction of liquid water in each layer (-)
  integer(i4b)    :: mLayerMatricHead                 ! matric head of water in the soil (m)
  ! define other state variables
  integer(i4b)    :: scalarAquiferStorage             ! relative aquifer storage -- above bottom of the soil profile (m)
  integer(i4b)    :: scalarSurfaceTemp                ! surface temperature (K)
  ! define coordinate variables
  integer(i4b)    :: mLayerDepth                      ! depth of each layer (m)
  integer(i4b)    :: mLayerHeight                     ! height at the mid-point of each layer (m)
  integer(i4b)    :: iLayerHeight                     ! height of the layer interface; top of soil = 0 (m)
 endtype iLook_prog

 ! ***********************************************************************************************************
 ! (7) define diagnostic variables
 ! ***********************************************************************************************************
 type, public :: iLook_diag
  ! local properties
  integer(i4b)    :: scalarGreenVegFraction           ! green vegetation fraction used to compute LAI (-)
  integer(i4b)    :: scalarBulkVolHeatCapVeg          ! bulk volumetric heat capacity of vegetation (J m-3 K-1)
  integer(i4b)    :: scalarCanopyEmissivity           ! effective canopy emissivity (-)
  integer(i4b)    :: scalarRootZoneTemp               ! average temperature of the root zone (K)
  integer(i4b)    :: scalarLAI                        ! one-sided leaf area index (m2 m-2)
  integer(i4b)    :: scalarSAI                        ! one-sided stem area index (m2 m-2)
  integer(i4b)    :: scalarExposedLAI                 ! exposed leaf area index after burial by snow (m2 m-2)
  integer(i4b)    :: scalarExposedSAI                 ! exposed stem area index after burial by snow(m2 m-2)
  integer(i4b)    :: scalarCanopyIceMax               ! maximum interception storage capacity for ice (kg m-2)
  integer(i4b)    :: scalarCanopyLiqMax               ! maximum interception storage capacity for liquid water (kg m-2)
  integer(i4b)    :: scalarGrowingSeasonIndex         ! growing season index (0=off, 1=on)
  integer(i4b)    :: scalarVolHtCap_air               ! volumetric heat capacity air (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_ice               ! volumetric heat capacity ice (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_soil              ! volumetric heat capacity dry soil (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_water             ! volumetric heat capacity liquid wat (J m-3 K-1)
  integer(i4b)    :: mLayerVolHtCapBulk               ! volumetric heat capacity in each layer (J m-3 K-1)
  integer(i4b)    :: scalarLambda_drysoil             ! thermal conductivity of dry soil     (W m-1 K-1)
  integer(i4b)    :: scalarLambda_wetsoil             ! thermal conductivity of wet soil     (W m-1 K-1)
  integer(i4b)    :: mLayerThermalC                   ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
  integer(i4b)    :: iLayerThermalC                   ! thermal conductivity at the interface of each layer (W m-1 K-1)
  ! forcing
  integer(i4b)    :: scalarVPair                      ! vapor pressure of the air above the vegetation canopy (Pa)
  integer(i4b)    :: scalarVP_CanopyAir               ! vapor pressure of the canopy air space (Pa)
  integer(i4b)    :: scalarTwetbulb                   ! wet bulb temperature (K)
  integer(i4b)    :: scalarSnowfallTemp               ! temperature of fresh snow (K)
  integer(i4b)    :: scalarNewSnowDensity             ! density of fresh snow (kg m-3)
  integer(i4b)    :: scalarO2air                      ! atmospheric o2 concentration (Pa)
  integer(i4b)    :: scalarCO2air                     ! atmospheric co2 concentration (Pa)
  ! shortwave radiation
  integer(i4b)    :: scalarCosZenith                  ! cosine of the solar zenith angle (0-1)
  integer(i4b)    :: scalarFractionDirect             ! fraction of direct radiation (0-1)
  integer(i4b)    :: scalarCanopySunlitFraction       ! sunlit fraction of canopy (-)
  integer(i4b)    :: scalarCanopySunlitLAI            ! sunlit leaf area (-)
  integer(i4b)    :: scalarCanopyShadedLAI            ! shaded leaf area (-)
  integer(i4b)    :: spectralAlbGndDirect             ! direct  albedo of underlying surface for each spectral band (-)
  integer(i4b)    :: spectralAlbGndDiffuse            ! diffuse albedo of underlying surface for each spectral band (-)
  integer(i4b)    :: scalarGroundAlbedo               ! albedo of the ground surface (-)
  ! turbulent heat transfer
  integer(i4b)    :: scalarLatHeatSubVapCanopy        ! latent heat of sublimation/vaporization used for veg canopy (J kg-1)
  integer(i4b)    :: scalarLatHeatSubVapGround        ! latent heat of sublimation/vaporization used for ground surface (J kg-1)
  integer(i4b)    :: scalarSatVP_CanopyTemp           ! saturation vapor pressure at the temperature of vegetation canopy (Pa)
  integer(i4b)    :: scalarSatVP_GroundTemp           ! saturation vapor pressure at the temperature of the ground (Pa)
  integer(i4b)    :: scalarZ0Canopy                   ! roughness length of the canopy (m)
  integer(i4b)    :: scalarWindReductionFactor        ! canopy wind reduction factor (-)
  integer(i4b)    :: scalarZeroPlaneDisplacement      ! zero plane displacement (m)
  integer(i4b)    :: scalarRiBulkCanopy               ! bulk Richardson number for the canopy (-)
  integer(i4b)    :: scalarRiBulkGround               ! bulk Richardson number for the ground surface (-)
  integer(i4b)    :: scalarCanopyStabilityCorrection  ! stability correction for the canopy (-)
  integer(i4b)    :: scalarGroundStabilityCorrection  ! stability correction for the ground surface (-)
  ! evapotranspiration
  integer(i4b)    :: scalarIntercellularCO2Sunlit     ! carbon dioxide partial pressure of leaf interior (sunlit leaves) (Pa)
  integer(i4b)    :: scalarIntercellularCO2Shaded     ! carbon dioxide partial pressure of leaf interior (shaded leaves) (Pa)
  integer(i4b)    :: scalarTranspireLim               ! aggregate soil moisture + aquifer storage limit on transpiration (-)
  integer(i4b)    :: scalarTranspireLimAqfr           ! aquifer storage limit on transpiration (-)
  integer(i4b)    :: scalarFoliageNitrogenFactor      ! foliage nitrogen concentration, 1=saturated (-)
  integer(i4b)    :: scalarSoilRelHumidity            ! relative humidity in the soil pores in the upper-most soil layer (-)
  integer(i4b)    :: mLayerTranspireLim               ! soil moist & veg limit on transpiration for each layer (-)
  integer(i4b)    :: mLayerRootDensity                ! fraction of roots in each soil layer (-)
  integer(i4b)    :: scalarAquiferRootFrac            ! fraction of roots below the soil profile (-)
  ! canopy hydrology
  integer(i4b)    :: scalarFracLiqVeg                 ! fraction of liquid water on vegetation (-)
  integer(i4b)    :: scalarCanopyWetFraction          ! fraction of canopy that is wet
  ! snow hydrology
  integer(i4b)    :: scalarSnowAge                    ! non-dimensional snow age (-)
  integer(i4b)    :: scalarGroundSnowFraction         ! fraction of ground that is covered with snow (-)
  integer(i4b)    :: spectralSnowAlbedoDirect         ! direct snow albedo for individual spectral bands (-)
  integer(i4b)    :: scalarFracLiqSnow                ! fraction of liquid water in each snow layer (-)
  integer(i4b)    :: mLayerThetaResid                 ! residual volumetric water content in each snow layer (-)
  integer(i4b)    :: mLayerPoreSpace                  ! total pore space in each snow layer (-)
  integer(i4b)    :: mLayerMeltFreeze                 ! change in ice content due to melt/freeze in each layer (kg m-3)
  ! soil hydrology
  integer(i4b)    :: scalarInfilArea                  ! fraction of unfrozen area where water can infiltrate (-)
  integer(i4b)    :: scalarFrozenArea                 ! fraction of area that is considered impermeable due to soil ice (-)
  integer(i4b)    :: scalarSoilControl                ! soil control on infiltration: 1=controlling; 0=not (-)
  integer(i4b)    :: mLayerVolFracAir                 ! volumetric fraction of air in each layer (-)
  integer(i4b)    :: mLayerTcrit                      ! critical soil temperature above which all water is unfrozen (K)
  integer(i4b)    :: mLayerCompress                   ! change in volumetric water content due to compression of soil (-)
  integer(i4b)    :: scalarSoilCompress               ! change in total soil storage due to compression of the soil matrix (kg m-2)
  ! mass balance check
  integer(i4b)    :: scalarSoilWatBalError            ! error in the total soil water balance (kg m-2)
  integer(i4b)    :: scalarAquiferBalError            ! error in the aquifer water balance (kg m-2)
  integer(i4b)    :: scalarTotalSoilLiq               ! total mass of liquid water in the soil (kg m-2)
  integer(i4b)    :: scalarTotalSoilIce               ! total mass of ice in the soil (kg m-2)
  ! variable shortcuts
  integer(i4b)    :: scalarVGn_m                      ! van Genuchten "m" parameter (-)
  integer(i4b)    :: scalarKappa                      ! constant in the freezing curve function (m K-1)
  integer(i4b)    :: scalarVolLatHt_fus               ! volumetric latent heat of fusion     (J m-3)
 endtype iLook_diag

 ! ***********************************************************************************************************
 ! (8) define model fluxes
 ! ***********************************************************************************************************
 type, public :: iLook_flux
  ! net energy and mass fluxes for the vegetation domain
  integer(i4b)    :: scalarCanairNetNrgFlux           ! net energy flux for the canopy air space (W m-2)
  integer(i4b)    :: scalarCanopyNetNrgFlux           ! net energy flux for the vegetation canopy (W m-2)
  integer(i4b)    :: scalarGroundNetNrgFlux           ! net energy flux for the ground surface (W m-2)
  integer(i4b)    :: scalarCanopyNetLiqFlux           ! net liquid water flux for the vegetation canopy (kg m-2 s-1)
  ! forcing
  integer(i4b)    :: scalarRainfall                   ! computed rainfall rate (kg m-2 s-1)
  integer(i4b)    :: scalarSnowfall                   ! computed snowfall rate (kg m-2 s-1)
  ! shortwave radiation
  integer(i4b)    :: spectralIncomingDirect           ! incoming direct solar radiation in each wave band (W m-2)
  integer(i4b)    :: spectralIncomingDiffuse          ! incoming diffuse solar radiation in each wave band (W m-2)
  integer(i4b)    :: scalarCanopySunlitPAR            ! average absorbed par for sunlit leaves (W m-2)
  integer(i4b)    :: scalarCanopyShadedPAR            ! average absorbed par for shaded leaves (W m-2)
  integer(i4b)    :: spectralBelowCanopyDirect        ! downward direct flux below veg layer for each spectral band  (W m-2)
  integer(i4b)    :: spectralBelowCanopyDiffuse       ! downward diffuse flux below veg layer for each spectral band (W m-2)
  integer(i4b)    :: scalarBelowCanopySolar           ! solar radiation transmitted below the canopy (W m-2)
  integer(i4b)    :: scalarCanopyAbsorbedSolar        ! solar radiation absorbed by canopy (W m-2)
  integer(i4b)    :: scalarGroundAbsorbedSolar        ! solar radiation absorbed by ground (W m-2)
  ! longwave radiation
  integer(i4b)    :: scalarLWRadCanopy                ! longwave radiation emitted from the canopy (W m-2)
  integer(i4b)    :: scalarLWRadGround                ! longwave radiation emitted at the ground surface  (W m-2)
  integer(i4b)    :: scalarLWRadUbound2Canopy         ! downward atmospheric longwave radiation absorbed by the canopy (W m-2)
  integer(i4b)    :: scalarLWRadUbound2Ground         ! downward atmospheric longwave radiation absorbed by the ground (W m-2)
  integer(i4b)    :: scalarLWRadUbound2Ubound         ! atmospheric radiation refl by ground + lost thru upper boundary (W m-2)
  integer(i4b)    :: scalarLWRadCanopy2Ubound         ! longwave radiation emitted from canopy lost thru upper boundary (W m-2)
  integer(i4b)    :: scalarLWRadCanopy2Ground         ! longwave radiation emitted from canopy absorbed by the ground (W m-2)
  integer(i4b)    :: scalarLWRadCanopy2Canopy         ! canopy longwave reflected from ground and absorbed by the canopy (W m-2)
  integer(i4b)    :: scalarLWRadGround2Ubound         ! longwave radiation emitted from ground lost thru upper boundary (W m-2)
  integer(i4b)    :: scalarLWRadGround2Canopy         ! longwave radiation emitted from ground and absorbed by the canopy (W m-2)
  integer(i4b)    :: scalarLWNetCanopy                ! net longwave radiation at the canopy (W m-2)
  integer(i4b)    :: scalarLWNetGround                ! net longwave radiation at the ground surface (W m-2)
  integer(i4b)    :: scalarLWNetUbound                ! net longwave radiation at the upper atmospheric boundary (W m-2)
  ! turbulent heat transfer
  integer(i4b)    :: scalarEddyDiffusCanopyTop        ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
  integer(i4b)    :: scalarFrictionVelocity           ! friction velocity - canopy momentum sink (m s-1)
  integer(i4b)    :: scalarWindspdCanopyTop           ! windspeed at the top of the canopy (m s-1)
  integer(i4b)    :: scalarWindspdCanopyBottom        ! windspeed at the height of the bottom of the canopy (m s-1)
  integer(i4b)    :: scalarGroundResistance           ! below canopy aerodynamic resistance (s m-1)
  integer(i4b)    :: scalarCanopyResistance           ! above canopy aerodynamic resistance (s m-1)
  integer(i4b)    :: scalarLeafResistance             ! mean leaf boundary layer resistance per unit leaf area (s m-1)
  integer(i4b)    :: scalarSoilResistance             ! soil surface resistance (s m-1)
  integer(i4b)    :: scalarSenHeatTotal               ! sensible heat from the canopy air space to the atmosphere (W m-2)
  integer(i4b)    :: scalarSenHeatCanopy              ! sensible heat from the canopy to the canopy air space (W m-2)
  integer(i4b)    :: scalarSenHeatGround              ! sensible heat from the ground (below canopy or non-vegetated) (W m-2)
  integer(i4b)    :: scalarLatHeatTotal               ! latent heat from the canopy air space to the atmosphere (W m-2)
  integer(i4b)    :: scalarLatHeatCanopyEvap          ! evaporation latent heat from the canopy to the canopy air space (W m-2)
  integer(i4b)    :: scalarLatHeatCanopyTrans         ! transpiration latent heat from the canopy to the canopy air space (W m-2)
  integer(i4b)    :: scalarLatHeatGround              ! latent heat from the ground (below canopy or non-vegetated) (W m-2)
  integer(i4b)    :: scalarCanopyAdvectiveHeatFlux    ! heat advected to the canopy surface with rain + snow (W m-2)
  integer(i4b)    :: scalarGroundAdvectiveHeatFlux    ! heat advected to the ground surface with throughfall and unloading/drainage (W m-2)
  integer(i4b)    :: scalarCanopySublimation          ! canopy sublimation/frost (kg m-2 s-1)
  integer(i4b)    :: scalarSnowSublimation            ! snow sublimation/frost (below canopy or non-vegetated) (kg m-2 s-1)
  ! liquid water fluxes associated with evapotranspiration
  integer(i4b)    :: scalarStomResistSunlit           ! stomatal resistance for sunlit leaves (s m-1)
  integer(i4b)    :: scalarStomResistShaded           ! stomatal resistance for shaded leaves (s m-1)
  integer(i4b)    :: scalarPhotosynthesisSunlit       ! sunlit photosynthesis (umolco2 m-2 s-1)
  integer(i4b)    :: scalarPhotosynthesisShaded       ! shaded photosynthesis (umolco2 m-2 s-1)
  integer(i4b)    :: scalarCanopyTranspiration        ! canopy transpiration (kg m-2 s-1)
  integer(i4b)    :: scalarCanopyEvaporation          ! canopy evaporation/condensation (kg m-2 s-1)
  integer(i4b)    :: scalarGroundEvaporation          ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
  integer(i4b)    :: mLayerTranspire                  ! transpiration loss from each soil layer (kg m-2 s-1)
  ! liquid and solid water fluxes through the canopy
  integer(i4b)    :: scalarThroughfallSnow            ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
  integer(i4b)    :: scalarThroughfallRain            ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
  integer(i4b)    :: scalarCanopySnowUnloading        ! unloading of snow from the vegetion canopy (kg m-2 s-1)
  integer(i4b)    :: scalarCanopyLiqDrainage          ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
  integer(i4b)    :: scalarCanopyMeltFreeze           ! melt/freeze of water stored in the canopy (kg m-2 s-1)
  ! energy fluxes and for the snow and soil domains
  integer(i4b)    :: iLayerConductiveFlux             ! conductive energy flux at layer interfaces (W m-2)
  integer(i4b)    :: iLayerAdvectiveFlux              ! advective energy flux at layer interfaces (W m-2)
  integer(i4b)    :: iLayerNrgFlux                    ! energy flux at layer interfaces (W m-2)
  integer(i4b)    :: mLayerNrgFlux                    ! net energy flux for each layer in the snow+soil domain (J m-3 s-1)
  ! liquid water fluxes for the snow domain
  integer(i4b)    :: iLayerLiqFluxSnow                ! liquid flux at snow layer interfaces (m s-1)
  integer(i4b)    :: mLayerLiqFluxSnow                ! net liquid water flux for each snow layer (s-1)
  ! liquid water fluxes for the soil domain
  integer(i4b)    :: scalarRainPlusMelt               ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  integer(i4b)    :: scalarInfiltration               ! infiltration of water into the soil profile (m s-1)
  integer(i4b)    :: scalarExfiltration               ! exfiltration of water from the top of the soil profile (m s-1)
  integer(i4b)    :: scalarSurfaceRunoff              ! surface runoff (m s-1)
  integer(i4b)    :: mLayerSatHydCondMP               ! saturated hydraulic conductivity of macropores in each layer (m s-1)
  integer(i4b)    :: mLayerSatHydCond                 ! saturated hydraulic conductivity in each layer (m s-1)
  integer(i4b)    :: iLayerSatHydCond                 ! saturated hydraulic conductivity at each layer interface (m s-1)
  integer(i4b)    :: mLayerHydCond                    ! hydraulic conductivity in each soil layer (m s-1)
  integer(i4b)    :: iLayerLiqFluxSoil                ! liquid flux at soil layer interfaces (m s-1)
  integer(i4b)    :: mLayerLiqFluxSoil                ! net liquid water flux for each soil layer (s-1)
  integer(i4b)    :: mLayerBaseflow                   ! baseflow from each soil layer (m s-1)
  integer(i4b)    :: mLayerColumnInflow               ! total inflow to each layer in a given soil column (m3 s-1)
  integer(i4b)    :: mLayerColumnOutflow              ! total outflow from each layer in a given soil column (m3 s-1)
  integer(i4b)    :: scalarSoilBaseflow               ! total baseflow from throughout the soil profile (m s-1)
  integer(i4b)    :: scalarSoilDrainage               ! drainage from the bottom of the soil profile (m s-1)
  integer(i4b)    :: scalarAquiferRecharge            ! recharge to the aquifer (m s-1)
  integer(i4b)    :: scalarAquiferTranspire           ! transpiration from the aquifer (m s-1)
  integer(i4b)    :: scalarAquiferBaseflow            ! baseflow from the aquifer (m s-1)
 endtype iLook_flux

 ! ***********************************************************************************************************
 ! (9) define derivatives
 ! ***********************************************************************************************************
 type, public :: iLook_deriv
  ! derivatives in net vegetation energy fluxes w.r.t. relevant state variables
  integer(i4b)    :: dCanairNetFlux_dCanairTemp   ! derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
  integer(i4b)    :: dCanairNetFlux_dCanopyTemp   ! derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
  integer(i4b)    :: dCanairNetFlux_dGroundTemp   ! derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
  integer(i4b)    :: dCanopyNetFlux_dCanairTemp   ! derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
  integer(i4b)    :: dCanopyNetFlux_dCanopyTemp   ! derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
  integer(i4b)    :: dCanopyNetFlux_dGroundTemp   ! derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
  integer(i4b)    :: dCanopyNetFlux_dCanLiq       ! derivative in net canopy fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
  integer(i4b)    :: dGroundNetFlux_dCanairTemp   ! derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
  integer(i4b)    :: dGroundNetFlux_dCanopyTemp   ! derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
  integer(i4b)    :: dGroundNetFlux_dGroundTemp   ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
  integer(i4b)    :: dGroundNetFlux_dCanLiq       ! derivative in net ground fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
  ! derivatives in evaporative fluxes w.r.t. relevant state variables
  integer(i4b)    :: dCanopyEvaporation_dTCanair  ! derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
  integer(i4b)    :: dCanopyEvaporation_dTCanopy  ! derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
  integer(i4b)    :: dCanopyEvaporation_dTGround  ! derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
  integer(i4b)    :: dCanopyEvaporation_dCanLiq   ! derivative in canopy evaporation w.r.t. canopy liquid water content (s-1)
  integer(i4b)    :: dGroundEvaporation_dTCanair  ! derivative in ground evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
  integer(i4b)    :: dGroundEvaporation_dTCanopy  ! derivative in ground evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
  integer(i4b)    :: dGroundEvaporation_dTGround  ! derivative in ground evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
  integer(i4b)    :: dGroundEvaporation_dCanLiq   ! derivative in ground evaporation w.r.t. canopy liquid water content (s-1)
  ! derivatives in canopy water w.r.t canopy temperature
  integer(i4b)    :: dTheta_dTkCanopy             ! derivative of volumetric liquid water content w.r.t. temperature (K-1)
  integer(i4b)    :: dCanLiq_dTcanopy             ! derivative of canopy liquid storage w.r.t. temperature (kg m-2 K-1)
  ! derivatives in canopy liquid fluxes w.r.t. canopy water
  integer(i4b)    :: scalarCanopyLiqDeriv         ! derivative in (throughfall + canopy drainage) w.r.t. canopy liquid water (s-1)
  integer(i4b)    :: scalarThroughfallRainDeriv   ! derivative in throughfall w.r.t. canopy liquid water (s-1)
  integer(i4b)    :: scalarCanopyLiqDrainageDeriv ! derivative in canopy drainage w.r.t. canopy liquid water (s-1)
  ! derivatives in energy fluxes at the interface of snow+soil layers w.r.t. temperature in layers above and below
  integer(i4b)    :: dNrgFlux_dTempAbove          ! derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
  integer(i4b)    :: dNrgFlux_dTempBelow          ! derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
  ! derivative in liquid water fluxes at the interface of snow layers w.r.t. volumetric liquid water content in the layer above
  integer(i4b)    :: iLayerLiqFluxSnowDeriv       ! derivative in vertical liquid water flux at layer interfaces (m s-1)
  ! derivative in liquid water fluxes for the soil domain w.r.t hydrology state variables
  integer(i4b)    :: dVolTot_dPsi0                ! derivative in total water content w.r.t. total water matric potential (m-1)
  integer(i4b)    :: dq_dHydStateAbove            ! change in the flux in layer interfaces w.r.t. state variables in the layer above
  integer(i4b)    :: dq_dHydStateBelow            ! change in the flux in layer interfaces w.r.t. state variables in the layer below
  integer(i4b)    :: mLayerdTheta_dPsi            ! derivative in the soil water characteristic w.r.t. psi (m-1)
  integer(i4b)    :: mLayerdPsi_dTheta            ! derivative in the soil water characteristic w.r.t. theta (m)
  integer(i4b)    :: dCompress_dPsi               ! derivative in compressibility w.r.t matric head (m-1)
  ! derivative in liquid water fluxes for the soil domain w.r.t energy state variables
  integer(i4b)    :: dq_dNrgStateAbove            ! change in the flux in layer interfaces w.r.t. state variables in the layer above
  integer(i4b)    :: dq_dNrgStateBelow            ! change in the flux in layer interfaces w.r.t. state variables in the layer below
  integer(i4b)    :: mLayerdTheta_dTk             ! derivative of volumetric liquid water content w.r.t. temperature (K-1)
  integer(i4b)    :: dPsiLiq_dTemp                ! derivative in the liquid water matric potential w.r.t. temperature (m K-1)
 endtype iLook_deriv

 ! ***********************************************************************************************************
 ! (10) define model indices
 ! ***********************************************************************************************************
 type, public :: iLook_index
  ! number of state variables of different type
  integer(i4b)    :: nVegNrg                 ! number of energy state variables for vegetation
  integer(i4b)    :: nVegMass                ! number of hydrology state variables for vegetation (mass of water)
  integer(i4b)    :: nVegState               ! number of vegetation state variables
  integer(i4b)    :: nNrgState               ! number of energy state variables
  integer(i4b)    :: nWatState               ! number of "total water" state variables (volumetric total water content) 
  integer(i4b)    :: nMatState               ! number of matric head state variables
  integer(i4b)    :: nMassState              ! number of hydrology state variables (mass of water)
  integer(i4b)    :: nState                  ! total number of model state variables
  ! number of model layers, and layer indices
  integer(i4b)    :: nSnow                   ! number of snow layers
  integer(i4b)    :: nSoil                   ! number of soil layers
  integer(i4b)    :: nLayers                 ! total number of layers
  integer(i4b)    :: layerType               ! type of layer (soil or snow)
  ! indices of model state variables
  integer(i4b)    :: ixCasNrg                ! index of the canopy air space state variable
  integer(i4b)    :: ixVegNrg                ! index of the canopy energy state variable
  integer(i4b)    :: ixVegWat                ! index of the canopy total water state variable
  integer(i4b)    :: ixTopNrg                ! index of the upper-most energy state variable in the snow-soil subdomain
  integer(i4b)    :: ixTopWat                ! index of the upper-most total water state variable in the snow-soil subdomain
  integer(i4b)    :: ixTopMat                ! index of the upper-most matric head state variable in the soil subdomain
  integer(i4b)    :: ixSnowSoilNrg           ! indices for energy state variables in the snow-soil subdomain
  integer(i4b)    :: ixSnowSoilWat           ! indices for total water state variables in the snow-soil subdomain
  integer(i4b)    :: ixSnowOnlyNrg           ! indices for energy state variables in the snow subdomain
  integer(i4b)    :: ixSnowOnlyWat           ! indices for total water state variables in the snow subdomain
  integer(i4b)    :: ixSoilOnlyNrg           ! indices for energy state variables in the soil subdomain
  integer(i4b)    :: ixSoilOnlyHyd           ! indices for hydrology state variables in the soil subdomain
  ! type of model state variables
  integer(i4b)    :: ixStateType             ! indices defining the type of the state (ixNrgState, ixWatState, ixMatState...)
  integer(i4b)    :: ixAllState              ! list of indices for all model state variables
  integer(i4b)    :: ixSoilState             ! list of indices for all soil layers
  integer(i4b)    :: ixLayerState            ! list of indices for all model layers
  integer(i4b)    :: ixNrgOnly               ! list of indices for all energy states
  integer(i4b)    :: ixWatOnly               ! list of indices for all "total water" state variables (volumetric total water content)
  integer(i4b)    :: ixMatOnly               ! list of indices for matric head state variables
  integer(i4b)    :: ixMassOnly              ! list of indices for hydrology state variables (mass of water)
  ! indices for the model output files
  integer(i4b)    :: midSnowStartIndex       ! start index of the midSnow vector for a given timestep
  integer(i4b)    :: midSoilStartIndex       ! start index of the midSoil vector for a given timestep
  integer(i4b)    :: midTotoStartIndex       ! start index of the midToto vector for a given timestep
  integer(i4b)    :: ifcSnowStartIndex       ! start index of the ifcSnow vector for a given timestep
  integer(i4b)    :: ifcSoilStartIndex       ! start index of the ifcSoil vector for a given timestep
  integer(i4b)    :: ifcTotoStartIndex       ! start index of the ifcToto vector for a given timestep
 endtype iLook_index

 ! ***********************************************************************************************************
 ! (11) define basin-average model parameters
 ! ***********************************************************************************************************
 type, public :: iLook_bpar
  ! baseflow
  integer(i4b)    :: basin__aquiferHydCond        ! hydraulic conductivity for the aquifer (m s-1)
  integer(i4b)    :: basin__aquiferScaleFactor    ! scaling factor for aquifer storage in the big bucket (m)
  integer(i4b)    :: basin__aquiferBaseflowExp    ! baseflow exponent for the big bucket (-)
  ! within-grid routing
  integer(i4b)    :: routingGammaShape            ! shape parameter in Gamma distribution used for sub-grid routing (-)
  integer(i4b)    :: routingGammaScale            ! scale parameter in Gamma distribution used for sub-grid routing (s)
 endtype iLook_bpar

 ! ***********************************************************************************************************
 ! (12) define basin-average model variables
 ! ***********************************************************************************************************
 type, public :: iLook_bvar
  ! define derived variables
  integer(i4b)    :: basin__totalArea                ! total basin area (m2)
  ! define fluxes
  integer(i4b)    :: basin__SurfaceRunoff            ! surface runoff (m s-1)
  integer(i4b)    :: basin__ColumnOutflow            ! outflow from all "outlet" HRUs (those with no downstream HRU)
  integer(i4b)    :: basin__AquiferStorage           ! aquifer storage (m s-1)
  integer(i4b)    :: basin__AquiferRecharge          ! recharge to the aquifer (m s-1)
  integer(i4b)    :: basin__AquiferBaseflow          ! baseflow from the aquifer (m s-1)
  integer(i4b)    :: basin__AquiferTranspire         ! transpiration from the aquifer (m s-1)
  ! define variables for runoff
  integer(i4b)    :: routingRunoffFuture             ! runoff in future time steps (m s-1)
  integer(i4b)    :: routingFractionFuture           ! fraction of runoff in future time steps (-)
  integer(i4b)    :: averageInstantRunoff            ! instantaneous runoff (m s-1)
  integer(i4b)    :: averageRoutedRunoff             ! routed runoff (m s-1)
 endtype iLook_bvar

 ! ***********************************************************************************************************
 ! (X) define data structures and maximum number of variables of each type
 ! ***********************************************************************************************************

 ! named variables: model decisions
 type(iLook_decision),public,parameter :: iLookDECISIONS=iLook_decision(  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,&
                                                                         21, 22, 23, 24, 25, 26, 27, 28, 29, 30,&
                                                                         31, 32, 33, 34, 35, 36, 37, 38)

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
                                                                        141,142,143,144,145,146,147)

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
                                                                         81, 82, 83, 84)

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

END MODULE var_lookup

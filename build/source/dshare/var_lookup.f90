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
 ! define missing value
 integer(i4b),parameter     :: imiss = -999      ! used to initialize named variables
 ! ***********************************************************************************************************
 ! (0) define model decisions
 ! ***********************************************************************************************************
 type, public  ::  iLook_decision
  integer(i4b)    :: simulStart       = 1  ! simulation start time
  integer(i4b)    :: simulFinsh       = 2  ! simulation end time
  integer(i4b)    :: soilCatTbl       = 3  ! soil-category dateset
  integer(i4b)    :: vegeParTbl       = 4  ! vegetation category dataset
  integer(i4b)    :: soilStress       = 5  ! choice of function for the soil moisture control on stomatal resistance
  integer(i4b)    :: stomResist       = 6  ! choice of function for stomatal resistance
  integer(i4b)    :: bbTempFunc       = 7  ! Ball-Berry: leaf temperature controls on photosynthesis + stomatal resistance
  integer(i4b)    :: bbHumdFunc       = 8  ! Ball-Berry: humidity controls on stomatal resistance
  integer(i4b)    :: bbElecFunc       = 9  ! Ball-Berry: dependence of photosynthesis on PAR
  integer(i4b)    :: bbCO2point       = 10 ! Ball-Berry: use of CO2 compensation point to calculate stomatal resistance
  integer(i4b)    :: num_method       = 11 ! choice of numerical method
  integer(i4b)    :: fDerivMeth       = 12 ! method used to calculate flux derivatives
  integer(i4b)    :: LAI_method       = 13 ! method used to determine LAI and SAI
  integer(i4b)    :: cIntercept       = 14 ! choice of parameterization for canopy interception
  integer(i4b)    :: f_Richards       = 15 ! form of richards' equation
  integer(i4b)    :: groundwatr       = 16 ! choice of groundwater parameterization
  integer(i4b)    :: hc_profile       = 17 ! choice of hydraulic conductivity profile
  integer(i4b)    :: bcUpprTdyn       = 18 ! type of upper boundary condition for thermodynamics
  integer(i4b)    :: bcLowrTdyn       = 19 ! type of lower boundary condition for thermodynamics
  integer(i4b)    :: bcUpprSoiH       = 20 ! type of upper boundary condition for soil hydrology
  integer(i4b)    :: bcLowrSoiH       = 21 ! type of lower boundary condition for soil hydrology
  integer(i4b)    :: veg_traits       = 22 ! choice of parameterization for vegetation roughness length and displacement height
  integer(i4b)    :: rootProfil       = 23 ! choice of parameterization for the rooting profile
  integer(i4b)    :: canopyEmis       = 24 ! choice of parameterization for canopy emissivity
  integer(i4b)    :: snowIncept       = 25 ! choice of parameterization for snow interception
  integer(i4b)    :: windPrfile       = 26 ! choice of canopy wind profile
  integer(i4b)    :: astability       = 27 ! choice of stability function
  integer(i4b)    :: canopySrad       = 28 ! choice of method for canopy shortwave radiation
  integer(i4b)    :: alb_method       = 29 ! choice of albedo representation
  integer(i4b)    :: snowLayers       = 30 ! choice of method to combine and sub-divide snow layers
  integer(i4b)    :: compaction       = 31 ! choice of compaction routine
  integer(i4b)    :: thCondSnow       = 32 ! choice of thermal conductivity representation for snow
  integer(i4b)    :: thCondSoil       = 33 ! choice of thermal conductivity representation for soil
  integer(i4b)    :: spatial_gw       = 34 ! choice of method for spatial representation of groundwater
  integer(i4b)    :: subRouting       = 35 ! choice of method for sub-grid routing
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
  integer(i4b)    :: upperBoundHead            = 1   ! matric head of the upper boundary (m)
  integer(i4b)    :: lowerBoundHead            = 2   ! matric head of the lower boundary (m)
  integer(i4b)    :: upperBoundTheta           = 3   ! volumetric liquid water content of the upper boundary (-)
  integer(i4b)    :: lowerBoundTheta           = 4   ! volumetric liquid water content of the lower boundary (-)
  integer(i4b)    :: upperBoundTemp            = 5   ! temperature of the upper boundary (K)
  integer(i4b)    :: lowerBoundTemp            = 6   ! temperature of the lower boundary (K)
  ! precipitation partitioning
  integer(i4b)    :: tempCritRain              = 7   ! critical temperature where precipitation is rain (K)
  integer(i4b)    :: tempRangeTimestep         = 8   ! temperature range over the time step (K)
  integer(i4b)    :: frozenPrecipMultip        = 9   ! frozen precipitation multiplier (-)
  ! freezing curve for snow
  integer(i4b)    :: snowfrz_scale             = 10  ! scaling parameter for the freezing curve for snow (K-1)
  ! snow albedo
  integer(i4b)    :: albedoMax                 = 11  ! maximum snow albedo for a single spectral band (-)
  integer(i4b)    :: albedoMinWinter           = 12  ! minimum snow albedo during winter for a single spectral band (-)
  integer(i4b)    :: albedoMinSpring           = 13  ! minimum snow albedo during spring for a single spectral band (-)
  integer(i4b)    :: albedoMaxVisible          = 14  ! maximum snow albedo in the visible part of the spectrum (-)
  integer(i4b)    :: albedoMinVisible          = 15  ! minimum snow albedo in the visible part of the spectrum (-)
  integer(i4b)    :: albedoMaxNearIR           = 16  ! maximum snow albedo in the near infra-red part of the spectrum (-)
  integer(i4b)    :: albedoMinNearIR           = 17  ! minimum snow albedo in the near infra-red part of the spectrum (-)
  integer(i4b)    :: albedoDecayRate           = 18  ! albedo decay rate (s)
  integer(i4b)    :: albedoSootLoad            = 19  ! soot load factor (-)
  integer(i4b)    :: albedoRefresh             = 20  ! critical mass necessary for albedo refreshment (kg m-2)
  ! radiation transfer within snow
  integer(i4b)    :: radExt_snow               = 21  ! extinction coefficient for radiation penetration into the snowpack (m-1)
  integer(i4b)    :: directScale               = 22  ! scaling factor for fractional driect radiaion parameterization (-)
  integer(i4b)    :: Frad_direct               = 23  ! maximum fraction of direct solar radiation (-)
  integer(i4b)    :: Frad_vis                  = 24  ! fraction of radiation in the visible part of the spectrum (-)
  ! new snow density
  integer(i4b)    :: newSnowDenMin             = 25  ! minimum new snow density (kg m-3)
  integer(i4b)    :: newSnowDenMult            = 26  ! multiplier for new snow density (kg m-3)
  integer(i4b)    :: newSnowDenScal            = 27  ! scaling factor for new snow density (K)
  ! snow compaction
  integer(i4b)    :: densScalGrowth            = 28  ! density scaling factor for grain growth (kg-1 m3)
  integer(i4b)    :: tempScalGrowth            = 29  ! temperature scaling factor for grain growth (K-1)
  integer(i4b)    :: grainGrowthRate           = 30  ! rate of grain growth (s-1)
  integer(i4b)    :: densScalOvrbdn            = 31  ! density scaling factor for overburden pressure (kg-1 m3)
  integer(i4b)    :: tempScalOvrbdn            = 32  ! temperature scaling factor for overburden pressure (K-1)
  integer(i4b)    :: base_visc                 = 33  ! viscosity coefficient at T=T_frz and snow density=0  (kg s m-2)
  ! water flow within snow
  integer(i4b)    :: Fcapil                    = 34  ! capillary retention as a fraction of the total pore volume (-)
  integer(i4b)    :: k_snow                    = 35  ! hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
  integer(i4b)    :: mw_exp                    = 36  ! exponent for meltwater flow (-)
  ! turbulent heat fluxes
  integer(i4b)    :: z0Snow                    = 37  ! roughness length of snow (m)
  integer(i4b)    :: z0Soil                    = 38  ! roughness length of bare soil below the canopy (m)
  integer(i4b)    :: z0Canopy                  = 39  ! roughness length of the canopy (m)
  integer(i4b)    :: zpdFraction               = 40  ! zero plane displacement / canopy height (-)
  integer(i4b)    :: critRichNumber            = 41  ! critical value for the bulk Richardson number (-)
  integer(i4b)    :: Louis79_bparam            = 42  ! parameter in Louis (1979) stability function (-)
  integer(i4b)    :: Louis79_cStar             = 43  ! parameter in Louis (1979) stability function (-)
  integer(i4b)    :: Mahrt87_eScale            = 44  ! exponential scaling factor in the Mahrt (1987) stability function (-)
  integer(i4b)    :: leafExchangeCoeff         = 45  ! turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
  integer(i4b)    :: windReductionParam        = 46  ! canopy wind reduction parameter (-)
  ! stomatal conductance
  integer(i4b)    :: Kc25                      = 47  ! Michaelis-Menten constant for CO2 at 25 degrees C (umol mol-1)
  integer(i4b)    :: Ko25                      = 48  ! Michaelis-Menten constant for O2 at 25 degrees C (mol mol-1)
  integer(i4b)    :: Kc_qFac                   = 49  ! factor in the q10 function defining temperature controls on Kc (-)
  integer(i4b)    :: Ko_qFac                   = 50  ! factor in the q10 function defining temperature controls on Ko (-)
  integer(i4b)    :: kc_Ha                     = 51  ! activation energy for the Michaelis-Menten constant for CO2 (J mol-1)
  integer(i4b)    :: ko_Ha                     = 52  ! activation energy for the Michaelis-Menten constant for O2 (J mol-1)
  integer(i4b)    :: vcmax25                   = 53  ! potential carboxylation rate at 25 degrees C (umol co2 m-2 s-1)
  integer(i4b)    :: vcmax_qFac                = 54  ! factor in the q10 function defining temperature controls on vcmax (-)
  integer(i4b)    :: vcmax_Ha                  = 55  ! activation energy in the vcmax function (J mol-1)
  integer(i4b)    :: vcmax_Hd                  = 56  ! deactivation energy in the vcmax function (J mol-1)
  integer(i4b)    :: vcmax_Sv                  = 57  ! entropy term in the vcmax function (J mol-1 K-1)
  integer(i4b)    :: jmax25_scale              = 58  ! scaling factor to relate jmax25 to vcmax25 (-)
  integer(i4b)    :: jmax_Ha                   = 59  ! activation energy in the jmax function (J mol-1)
  integer(i4b)    :: jmax_Hd                   = 60  ! deactivation energy in the jmax function (J mol-1)
  integer(i4b)    :: jmax_Sv                   = 61  ! entropy term in the jmax function (J mol-1 K-1)
  integer(i4b)    :: fractionJ                 = 62  ! fraction of light lost by other than the chloroplast lamellae (-)
  integer(i4b)    :: quantamYield              = 63  ! quantam yield (mol e mol-1 q)
  integer(i4b)    :: vpScaleFactor             = 64  ! vapor pressure scaling factor in stomatal conductance function (Pa)
  integer(i4b)    :: cond2photo_slope          = 65  ! slope of conductance-photosynthesis relationship (-)
  integer(i4b)    :: minStomatalConductance    = 66  ! minimum stomatal conductance (umol H2O m-2 s-1)
  ! vegetation properties
  integer(i4b)    :: winterSAI                 = 67  ! stem area index prior to the start of the growing season (m2 m-2)
  integer(i4b)    :: summerLAI                 = 68  ! maximum leaf area index at the peak of the growing season (m2 m-2)
  integer(i4b)    :: rootScaleFactor1          = 69  ! 1st scaling factor (a) in Y = 1 - 0.5*( exp(-aZ) + exp(-bZ) )  (m-1)
  integer(i4b)    :: rootScaleFactor2          = 70  ! 2nd scaling factor (b) in Y = 1 - 0.5*( exp(-aZ) + exp(-bZ) )  (m-1)
  integer(i4b)    :: rootingDepth              = 71  ! rooting depth (m)
  integer(i4b)    :: rootDistExp               = 72  ! exponent controlling the vertical distribution of root density (-)
  integer(i4b)    :: plantWiltPsi              = 73  ! matric head at wilting point (m)
  integer(i4b)    :: soilStressParam           = 74  ! parameter in the exponential soil stress function
  integer(i4b)    :: critSoilWilting           = 75  ! critical vol. liq. water content when plants are wilting (-)
  integer(i4b)    :: critSoilTranspire         = 76  ! critical vol. liq. water content when transpiration is limited (-)
  integer(i4b)    :: critAquiferTranspire      = 77  ! critical aquifer storage value when transpiration is limited (m)
  integer(i4b)    :: minStomatalResistance     = 78  ! minimum canopy resistance (s m-1)
  integer(i4b)    :: leafDimension             = 79  ! characteristic leaf dimension (m)
  integer(i4b)    :: heightCanopyTop           = 80  ! height of top of the vegetation canopy above ground surface (m)
  integer(i4b)    :: heightCanopyBottom        = 81  ! height of bottom of the vegetation canopy above ground surface (m)
  integer(i4b)    :: specificHeatVeg           = 82  ! specific heat of vegetation (J kg-1 K-1)
  integer(i4b)    :: maxMassVegetation         = 83  ! maximum mass of vegetation (full foliage) (kg m-2)
  integer(i4b)    :: throughfallScaleSnow      = 84  ! scaling factor for throughfall (snow) (-)
  integer(i4b)    :: throughfallScaleRain      = 85  ! scaling factor for throughfall (rain) (-)
  integer(i4b)    :: refInterceptCapSnow       = 86  ! reference canopy interception capacity per unit leaf area (snow) (kg m-2)
  integer(i4b)    :: refInterceptCapRain       = 87  ! canopy interception capacity per unit leaf area (rain) (kg m-2)
  integer(i4b)    :: snowUnloadingCoeff        = 88  ! time constant for unloading of snow from the forest canopy (s-1)
  integer(i4b)    :: canopyDrainageCoeff       = 89  ! time constant for drainage of liquid water from the forest canopy (s-1)
  integer(i4b)    :: ratioDrip2Unloading       = 90  ! ratio of canopy drip to unloading of snow from the forest canopy (-)
  integer(i4b)    :: canopyWettingFactor       = 91  ! maximum wetted fraction of the canopy (-)
  integer(i4b)    :: canopyWettingExp          = 92  ! exponent in canopy wetting function (-)
  ! soil properties
  integer(i4b)    :: soil_dens_intr            = 93  ! intrinsic soil density (kg m-3)
  integer(i4b)    :: thCond_soil               = 94  ! thermal conductivity of soil (W m-1 K-1)
  integer(i4b)    :: frac_sand                 = 95  ! fraction of sand (-)
  integer(i4b)    :: frac_silt                 = 96  ! fraction of silt (-)
  integer(i4b)    :: frac_clay                 = 97  ! fraction of clay (-)
  integer(i4b)    :: fieldCapacity             = 98  ! field capacity (-)
  integer(i4b)    :: wettingFrontSuction       = 99  ! Green-Ampt wetting front suction (m)
  integer(i4b)    :: theta_mp                  = 100 ! volumetric liquid water content when macropore flow begins (-)
  integer(i4b)    :: theta_sat                 = 101 ! porosity (-)
  integer(i4b)    :: theta_res                 = 102 ! volumetric residual water content (-)
  integer(i4b)    :: vGn_alpha                 = 103 ! van Genuchten "alpha" parameter (m-1)
  integer(i4b)    :: vGn_n                     = 104 ! van Genuchten "n" parameter (-)
  integer(i4b)    :: mpExp                     = 105 ! empirical exponent in macropore flow equation (-)
  integer(i4b)    :: k_soil                    = 106 ! hydraulic conductivity of soil (m s-1)
  integer(i4b)    :: k_macropore               = 107 ! saturated hydraulic conductivity for macropores (m s-1)
  integer(i4b)    :: kAnisotropic              = 108 ! anisotropy factor for lateral hydraulic conductivity (-)
  integer(i4b)    :: zScale_TOPMODEL           = 109 ! TOPMODEL scaling factor used in lower boundary condition for soil (m)
  integer(i4b)    :: compactedDepth            = 110 ! depth where k_soil reaches the compacted value given by CH78 (m)
  integer(i4b)    :: aquiferScaleFactor        = 111 ! scaling factor for aquifer storage in the big bucket (m)
  integer(i4b)    :: aquiferBaseflowExp        = 112 ! baseflow exponent (-)
  integer(i4b)    :: qSurfScale                = 113 ! scaling factor in the surface runoff parameterization (-)
  integer(i4b)    :: specificYield             = 114 ! specific yield (-)
  integer(i4b)    :: specificStorage           = 115 ! specific storage coefficient (m-1)
  integer(i4b)    :: f_impede                  = 116 ! ice impedence factor (-)
  integer(i4b)    :: soilIceScale              = 117 ! scaling factor for depth of soil ice, used to get frozen fraction (m)
  integer(i4b)    :: soilIceCV                 = 118 ! CV of depth of soil ice, used to get frozen fraction (-)
  ! algorithmic control parameters
  integer(i4b)    :: minwind                   = 119 ! minimum wind speed (m s-1)
  integer(i4b)    :: minstep                   = 120 ! minimum length of the time step
  integer(i4b)    :: maxstep                   = 121 ! maximum length of the time step
  integer(i4b)    :: wimplicit                 = 122 ! weight assigned to the start-of-step fluxes
  integer(i4b)    :: maxiter                   = 123 ! maximum number of iteration
  integer(i4b)    :: relConvTol_liquid         = 124 ! relative convergence tolerance for vol frac liq water (-)
  integer(i4b)    :: absConvTol_liquid         = 125 ! absolute convergence tolerance for vol frac liq water (-)
  integer(i4b)    :: relConvTol_matric         = 126 ! relative convergence tolerance for matric head (-)
  integer(i4b)    :: absConvTol_matric         = 127 ! absolute convergence tolerance for matric head (m)
  integer(i4b)    :: relConvTol_energy         = 128 ! relative convergence tolerance for energy (-)
  integer(i4b)    :: absConvTol_energy         = 129 ! absolute convergence tolerance for energy (J m-3)
  integer(i4b)    :: relConvTol_aquifr         = 130 ! relative convergence tolerance for aquifer storage (-)
  integer(i4b)    :: absConvTol_aquifr         = 131 ! absolute convergence tolerance for aquifer storage (J m-3)
  integer(i4b)    :: zmin                      = 132 ! minimum layer depth (m)
  integer(i4b)    :: zmax                      = 133 ! maximum layer depth (m)
  integer(i4b)    :: zminLayer1                = 134 ! minimum layer depth for the 1st (top) layer (m)
  integer(i4b)    :: zminLayer2                = 135 ! minimum layer depth for the 2nd layer (m)
  integer(i4b)    :: zminLayer3                = 136 ! minimum layer depth for the 3rd layer (m)
  integer(i4b)    :: zminLayer4                = 137 ! minimum layer depth for the 4th layer (m)
  integer(i4b)    :: zminLayer5                = 138 ! minimum layer depth for the 5th (bottom) layer (m)
  integer(i4b)    :: zmaxLayer1_lower          = 139 ! maximum layer depth for the 1st (top) layer when only 1 layer (m)
  integer(i4b)    :: zmaxLayer2_lower          = 140 ! maximum layer depth for the 2nd layer when only 2 layers (m)
  integer(i4b)    :: zmaxLayer3_lower          = 141 ! maximum layer depth for the 3rd layer when only 3 layers (m)
  integer(i4b)    :: zmaxLayer4_lower          = 142 ! maximum layer depth for the 4th layer when only 4 layers (m)
  integer(i4b)    :: zmaxLayer1_upper          = 143 ! maximum layer depth for the 1st (top) layer when > 1 layer (m)
  integer(i4b)    :: zmaxLayer2_upper          = 144 ! maximum layer depth for the 2nd layer when > 2 layers (m)
  integer(i4b)    :: zmaxLayer3_upper          = 145 ! maximum layer depth for the 3rd layer when > 3 layers (m)
  integer(i4b)    :: zmaxLayer4_upper          = 146 ! maximum layer depth for the 4th layer when > 4 layers (m)
 endtype ilook_param
 ! ***********************************************************************************************************
 ! (6) define model variables
 ! ***********************************************************************************************************
 type, public :: iLook_mvar
  ! define timestep-average fluxes for a few key variables
  integer(i4b)    :: totalSoilCompress               = 1   ! change in total soil storage due to compression of the soil matrix (kg m-2)
  integer(i4b)    :: averageThroughfallSnow          = 2   ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
  integer(i4b)    :: averageThroughfallRain          = 3   ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
  integer(i4b)    :: averageCanopySnowUnloading      = 4   ! unloading of snow from the vegetion canopy (kg m-2 s-1)
  integer(i4b)    :: averageCanopyLiqDrainage        = 5   ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
  integer(i4b)    :: averageCanopyMeltFreeze         = 6   ! melt/freeze of water stored in the canopy (kg m-2 s-1)
  integer(i4b)    :: averageCanopyTranspiration      = 7   ! canopy transpiration (kg m-2 s-1)
  integer(i4b)    :: averageCanopyEvaporation        = 8   ! canopy evaporation/condensation (kg m-2 s-1)
  integer(i4b)    :: averageCanopySublimation        = 9   ! canopy sublimation/frost (kg m-2 s-1)
  integer(i4b)    :: averageSnowSublimation          = 10  ! snow sublimation/frost - below canopy or non-vegetated (kg m-2 s-1)
  integer(i4b)    :: averageGroundEvaporation        = 11  ! ground evaporation/condensation - below canopy or non-vegetated (kg m-2 s-1)
  integer(i4b)    :: averageRainPlusMelt             = 12  ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  integer(i4b)    :: averageSurfaceRunoff            = 13  ! surface runoff (m s-1)
  integer(i4b)    :: averageSoilInflux               = 14  ! influx of water at the top of the soil profile (m s-1)
  integer(i4b)    :: averageSoilBaseflow             = 15  ! total baseflow from throughout the soil profile (m s-1)
  integer(i4b)    :: averageSoilDrainage             = 16  ! drainage from the bottom of the soil profile (m s-1)
  integer(i4b)    :: averageAquiferRecharge          = 17  ! recharge to the aquifer (m s-1)
  integer(i4b)    :: averageAquiferBaseflow          = 18  ! baseflow from the aquifer (m s-1)
  integer(i4b)    :: averageAquiferTranspire         = 19  ! transpiration from the aquifer (m s-1)
  integer(i4b)    :: averageColumnOutflow            = 20  ! outflow from each layer in the soil profile (m3 s-1)
  ! define scalar variables -- forcing
  integer(i4b)    :: scalarCosZenith                 = 21  ! cosine of the solar zenith angle (0-1)
  integer(i4b)    :: scalarFractionDirect            = 22  ! fraction of direct radiation (0-1)
  integer(i4b)    :: spectralIncomingDirect          = 23  ! incoming direct solar radiation in each wave band (W m-2)
  integer(i4b)    :: spectralIncomingDiffuse         = 24  ! incoming diffuse solar radiation in each wave band (W m-2)
  integer(i4b)    :: scalarVPair                     = 25  ! vapor pressure of the air above the vegetation canopy (Pa)
  integer(i4b)    :: scalarTwetbulb                  = 26  ! wet bulb temperature (K)
  integer(i4b)    :: scalarRainfall                  = 27  ! computed rainfall rate (kg m-2 s-1)
  integer(i4b)    :: scalarSnowfall                  = 28  ! computed snowfall rate (kg m-2 s-1)
  integer(i4b)    :: scalarSnowfallTemp              = 29  ! temperature of fresh snow (K)
  integer(i4b)    :: scalarNewSnowDensity            = 30  ! density of fresh snow, should snow be falling in this time step (kg m-3)
  integer(i4b)    :: scalarO2air                     = 31  ! atmospheric o2 concentration (Pa)
  integer(i4b)    :: scalarCO2air                    = 32  ! atmospheric co2 concentration (Pa)
  ! define scalar variables -- state variables
  integer(i4b)    :: scalarCanopyIce                 = 33  ! mass of ice on the vegetation canopy (kg m-2)
  integer(i4b)    :: scalarCanopyLiq                 = 34  ! mass of liquid water on the vegetation canopy (kg m-2)
  integer(i4b)    :: scalarCanairTemp                = 35  ! temperature of the canopy air space (Pa)
  integer(i4b)    :: scalarCanopyTemp                = 36  ! temperature of the vegetation canopy (K)
  integer(i4b)    :: scalarSnowAge                   = 37  ! non-dimensional snow age (-)
  integer(i4b)    :: scalarSnowAlbedo                = 38  ! snow albedo for the entire spectral band (-)
  integer(i4b)    :: spectralSnowAlbedoDirect        = 39  ! direct snow albedo for individual spectral bands (-)
  integer(i4b)    :: spectralSnowAlbedoDiffuse       = 40  ! diffuse snow albedo for individual spectral bands (-)
  integer(i4b)    :: scalarSnowDepth                 = 41  ! total snow depth (m)
  integer(i4b)    :: scalarSWE                       = 42  ! snow water equivalent (kg m-2)
  integer(i4b)    :: scalarSfcMeltPond               = 43  ! ponded water caused by melt of the "snow without a layer" (kg m-2)
  integer(i4b)    :: scalarAquiferStorage            = 44  ! relative aquifer storage -- above bottom of the soil profile (m)
  integer(i4b)    :: scalarSurfaceTemp               = 45  ! surface temperature (K)
  ! define NOAH-MP vegetation variables -- general
  integer(i4b)    :: scalarGreenVegFraction          = 46  ! green vegetation fraction used to compute LAI (-)
  integer(i4b)    :: scalarBulkVolHeatCapVeg         = 47  ! bulk volumetric heat capacity of vegetation (J m-3 K-1)
  integer(i4b)    :: scalarRootZoneTemp              = 48  ! average temperature of the root zone (K)
  integer(i4b)    :: scalarLAI                       = 49  ! one-sided leaf area index (m2 m-2)
  integer(i4b)    :: scalarSAI                       = 50  ! one-sided stem area index (m2 m-2)
  integer(i4b)    :: scalarExposedLAI                = 51  ! exposed leaf area index after burial by snow (m2 m-2)
  integer(i4b)    :: scalarExposedSAI                = 52  ! exposed stem area index after burial by snow(m2 m-2)
  integer(i4b)    :: scalarCanopyIceMax              = 53  ! maximum interception storage capacity for ice (kg m-2)
  integer(i4b)    :: scalarCanopyLiqMax              = 54  ! maximum interception storage capacity for liquid water (kg m-2)
  integer(i4b)    :: scalarGrowingSeasonIndex        = 55  ! growing season index (0=off, 1=on)
  integer(i4b)    :: scalarVP_CanopyAir              = 56  ! vapor pressure of the canopy air space (Pa)
  ! define NOAH-MP vegetation variables -- shortwave radiation
  integer(i4b)    :: scalarCanopySunlitFraction      = 57  ! sunlit fraction of canopy (-)
  integer(i4b)    :: scalarCanopySunlitLAI           = 58  ! sunlit leaf area (-)
  integer(i4b)    :: scalarCanopyShadedLAI           = 59  ! shaded leaf area (-)
  integer(i4b)    :: scalarCanopySunlitPAR           = 60  ! average absorbed par for sunlit leaves (w m-2)
  integer(i4b)    :: scalarCanopyShadedPAR           = 61  ! average absorbed par for shaded leaves (w m-2)
  integer(i4b)    :: spectralBelowCanopyDirect       = 62  ! downward direct flux below veg layer for each spectral band  W m-2)
  integer(i4b)    :: spectralBelowCanopyDiffuse      = 63  ! downward diffuse flux below veg layer for each spectral band (W m-2)
  integer(i4b)    :: scalarBelowCanopySolar          = 64  ! solar radiation transmitted below the canopy (W m-2)
  integer(i4b)    :: spectralAlbGndDirect            = 65  ! direct  albedo of underlying surface for each spectral band (-)
  integer(i4b)    :: spectralAlbGndDiffuse           = 66  ! diffuse albedo of underlying surface for each spectral band (-)
  integer(i4b)    :: scalarGroundAlbedo              = 67  ! albedo of the ground surface (-)
  integer(i4b)    :: scalarCanopyAbsorbedSolar       = 68  ! solar radiation absorbed by canopy (W m-2)
  integer(i4b)    :: scalarGroundAbsorbedSolar       = 69  ! solar radiation absorbed by ground (W m-2)
  ! define NOAH-MP vegetation variables -- longwave radiation
  integer(i4b)    :: scalarCanopyEmissivity          = 70  ! effective canopy emissivity (-)
  integer(i4b)    :: scalarLWRadCanopy               = 71  ! longwave radiation emitted from the canopy (W m-2)
  integer(i4b)    :: scalarLWRadGround               = 72  ! longwave radiation emitted at the ground surface  (W m-2)
  integer(i4b)    :: scalarLWRadUbound2Canopy        = 73  ! downward atmospheric longwave radiation absorbed by the canopy (W m-2)
  integer(i4b)    :: scalarLWRadUbound2Ground        = 74  ! downward atmospheric longwave radiation absorbed by the ground (W m-2)
  integer(i4b)    :: scalarLWRadUbound2Ubound        = 75  ! atmospheric radiation refl by ground + lost thru upper boundary (W m-2)
  integer(i4b)    :: scalarLWRadCanopy2Ubound        = 76  ! longwave radiation emitted from canopy lost thru upper boundary (W m-2)
  integer(i4b)    :: scalarLWRadCanopy2Ground        = 77  ! longwave radiation emitted from canopy absorbed by the ground (W m-2)
  integer(i4b)    :: scalarLWRadCanopy2Canopy        = 78  ! canopy longwave reflected from ground and absorbed by the canopy (W m-2)
  integer(i4b)    :: scalarLWRadGround2Ubound        = 79  ! longwave radiation emitted from ground lost thru upper boundary (W m-2)
  integer(i4b)    :: scalarLWRadGround2Canopy        = 80  ! longwave radiation emitted from ground and absorbed by the canopy (W m-2)
  integer(i4b)    :: scalarLWNetCanopy               = 81  ! net longwave radiation at the canopy (W m-2)
  integer(i4b)    :: scalarLWNetGround               = 82  ! net longwave radiation at the ground surface (W m-2)
  integer(i4b)    :: scalarLWNetUbound               = 83  ! net longwave radiation at the upper atmospheric boundary (W m-2)
  ! define NOAH-MP vegetation variables -- turbulent heat transfer
  integer(i4b)    :: scalarLatHeatSubVapCanopy       = 84  ! latent heat of sublimation/vaporization used for veg canopy (J kg-1)
  integer(i4b)    :: scalarLatHeatSubVapGround       = 85  ! latent heat of sublimation/vaporization used for ground surface (J kg-1)
  integer(i4b)    :: scalarSatVP_CanopyTemp          = 86  ! saturation vapor pressure at the temperature of vegetation canopy (Pa)
  integer(i4b)    :: scalarSatVP_GroundTemp          = 87  ! saturation vapor pressure at the temperature of the ground (Pa)
  integer(i4b)    :: scalarZ0Canopy                  = 88  ! roughness length of the canopy (m)
  integer(i4b)    :: scalarWindReductionFactor       = 89  ! canopy wind reduction factor (-)
  integer(i4b)    :: scalarZeroPlaneDisplacement     = 90  ! zero plane displacement (m)
  integer(i4b)    :: scalarRiBulkCanopy              = 91  ! bulk Richardson number for the canopy (-)
  integer(i4b)    :: scalarRiBulkGround              = 92  ! bulk Richardson number for the ground surface (-)
  integer(i4b)    :: scalarCanopyStabilityCorrection = 93  ! stability correction for the canopy (-)
  integer(i4b)    :: scalarGroundStabilityCorrection = 94  ! stability correction for the ground surface (-)
  integer(i4b)    :: scalarEddyDiffusCanopyTop       = 95  ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
  integer(i4b)    :: scalarFrictionVelocity          = 96  ! friction velocity - canopy momentum sink (m s-1)
  integer(i4b)    :: scalarWindspdCanopyTop          = 97  ! windspeed at the top of the canopy (m s-1)
  integer(i4b)    :: scalarWindspdCanopyBottom       = 98  ! windspeed at the height of the bottom of the canopy (m s-1)
  integer(i4b)    :: scalarGroundResistance          = 99  ! below canopy aerodynamic resistance (s m-1)
  integer(i4b)    :: scalarCanopyResistance          = 100 ! above canopy aerodynamic resistance (s m-1)
  integer(i4b)    :: scalarLeafResistance            = 101 ! mean leaf boundary layer resistance per unit leaf area (s m-1)
  integer(i4b)    :: scalarSoilResistance            = 102 ! soil surface resistance (s m-1)
  integer(i4b)    :: scalarSoilRelHumidity           = 103 ! relative humidity in the soil pores in the upper-most soil layer (-)
  integer(i4b)    :: scalarSenHeatTotal              = 104 ! sensible heat from the canopy air space to the atmosphere (W m-2)
  integer(i4b)    :: scalarSenHeatCanopy             = 105 ! sensible heat from the canopy to the canopy air space (W m-2)
  integer(i4b)    :: scalarSenHeatGround             = 106 ! sensible heat from the ground (below canopy or non-vegetated) (W m-2)
  integer(i4b)    :: scalarLatHeatTotal              = 107 ! latent heat from the canopy air space to the atmosphere (W m-2)
  integer(i4b)    :: scalarLatHeatCanopyEvap         = 108 ! evaporation latent heat from the canopy to the canopy air space (W m-2)
  integer(i4b)    :: scalarLatHeatCanopyTrans        = 109 ! transpiration latent heat from the canopy to the canopy air space (W m-2)
  integer(i4b)    :: scalarLatHeatGround             = 110 ! latent heat from the ground (below canopy or non-vegetated) (W m-2)
  integer(i4b)    :: scalarCanopyAdvectiveHeatFlux   = 111 ! heat advected to the canopy surface with rain + snow (W m-2)
  integer(i4b)    :: scalarGroundAdvectiveHeatFlux   = 112 ! heat advected to the ground surface with throughfall and unloading/drainage (W m-2)
  integer(i4b)    :: scalarCanopyTranspiration       = 113 ! canopy transpiration (kg m-2 s-1)
  integer(i4b)    :: scalarCanopyEvaporation         = 114 ! canopy evaporation/condensation (kg m-2 s-1)
  integer(i4b)    :: scalarCanopySublimation         = 115 ! canopy sublimation/frost (kg m-2 s-1)
  integer(i4b)    :: scalarGroundEvaporation         = 116 ! ground evaporation/condensation (below canopy or non-vegetated) (kg m-2 s-1)
  integer(i4b)    :: scalarSnowSublimation           = 117 ! snow sublimation/frost (below canopy or non-vegetated) (kg m-2 s-1)
  ! define NOAH-MP vegetation variables -- transpiration
  integer(i4b)    :: scalarTranspireLim              = 118 ! aggregate soil moisture + aquifer storage limit on transpiration (-)
  integer(i4b)    :: scalarTranspireLimAqfr          = 119 ! aquifer storage limit on transpiration (-)
  integer(i4b)    :: scalarFoliageNitrogenFactor     = 120 ! foliage nitrogen concentration, 1=saturated (-)
  integer(i4b)    :: scalarStomResistSunlit          = 121 ! stomatal resistance for sunlit leaves (s m-1)
  integer(i4b)    :: scalarStomResistShaded          = 122 ! stomatal resistance for shaded leaves (s m-1)
  integer(i4b)    :: scalarPhotosynthesisSunlit      = 123 ! sunlit photosynthesis (umolco2 m-2 s-1)
  integer(i4b)    :: scalarPhotosynthesisShaded      = 124 ! shaded photosynthesis (umolco2 m-2 s-1)
  ! define vegetation variables -- canopy water
  integer(i4b)    :: scalarCanopyWetFraction         = 125 ! fraction of canopy that is wet
  integer(i4b)    :: scalarGroundSnowFraction        = 126 ! fraction of ground that is covered with snow (-)
  integer(i4b)    :: scalarThroughfallSnow           = 127 ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
  integer(i4b)    :: scalarThroughfallRain           = 128 ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
  integer(i4b)    :: scalarCanopySnowUnloading       = 129 ! unloading of snow from the vegetion canopy (kg m-2 s-1)
  integer(i4b)    :: scalarCanopyLiqDrainage         = 130 ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
  integer(i4b)    :: scalarCanopyMeltFreeze          = 131 ! melt/freeze of water stored in the canopy (kg m-2 s-1)
  ! define scalar variables -- soil and aquifer fluxes
  integer(i4b)    :: scalarRainPlusMelt              = 132 ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  integer(i4b)    :: scalarInfilArea                 = 133 ! fraction of unfrozen area where water can infiltrate (-)
  integer(i4b)    :: scalarFrozenArea                = 134 ! fraction of area that is considered impermeable due to soil ice (-)
  integer(i4b)    :: scalarInfiltration              = 135 ! infiltration of water into the soil profile (m s-1)
  integer(i4b)    :: scalarExfiltration              = 136 ! exfiltration of water from the top of the soil profile (m s-1)
  integer(i4b)    :: scalarSurfaceRunoff             = 137 ! surface runoff (m s-1)
  integer(i4b)    :: scalarInitAquiferRecharge       = 138 ! recharge to the aquifer at the start of the step (m s-1)
  integer(i4b)    :: scalarAquiferRecharge           = 139 ! recharge to the aquifer (m s-1)
  integer(i4b)    :: scalarInitAquiferTranspire      = 140 ! transpiration from the aquifer at the start of the step (m s-1)
  integer(i4b)    :: scalarAquiferTranspire          = 141 ! transpiration from the aquifer (m s-1)
  integer(i4b)    :: scalarInitAquiferBaseflow       = 142 ! baseflow from the aquifer at the start of the step (m s-1)
  integer(i4b)    :: scalarAquiferBaseflow           = 143 ! baseflow from the aquifer (m s-1)
  ! scalar variables -- sub-step average fluxes for the soil zone
  integer(i4b)    :: scalarSoilInflux                = 144 ! influx of water at the top of the soil profile (m s-1)
  integer(i4b)    :: scalarSoilCompress              = 145 ! change in total soil storage due to compression of the soil matrix (kg m-2)
  integer(i4b)    :: scalarSoilBaseflow              = 146 ! sub-step average: total baseflow from throughout the soil profile (m s-1)
  integer(i4b)    :: scalarSoilDrainage              = 147 ! sub-step average: drainage from the bottom of the soil profile (m s-1)
  integer(i4b)    :: scalarSoilTranspiration         = 148 ! sub-step average: total transpiration from the soil (m s-1)
  ! define scalar variables -- mass balance check
  integer(i4b)    :: scalarSoilWatBalError           = 149 ! error in the total soil water balance (kg m-2)
  integer(i4b)    :: scalarAquiferBalError           = 150 ! error in the aquifer water balance (kg m-2)
  integer(i4b)    :: scalarTotalSoilLiq              = 151 ! total mass of liquid water in the soil (kg m-2)
  integer(i4b)    :: scalarTotalSoilIce              = 152 ! total mass of ice in the soil (kg m-2)
  ! define variables at the mid-point of each layer -- domain geometry
  integer(i4b)    :: mLayerDepth                     = 153 ! depth of each layer (m)
  integer(i4b)    :: mLayerHeight                    = 154 ! height at the mid-point of each layer (m)
  integer(i4b)    :: mLayerRootDensity               = 155 ! fraction of roots in each soil layer (-)
  ! define variables at the mid-point of each layer -- coupled energy and mass
  integer(i4b)    :: mLayerTemp                      = 156 ! temperature of each layer (K)
  integer(i4b)    :: mLayerVolFracAir                = 157 ! volumetric fraction of air in each layer (-)
  integer(i4b)    :: mLayerVolFracIce                = 158 ! volumetric fraction of ice water in each layer (-)
  integer(i4b)    :: mLayerVolFracLiq                = 159 ! volumetric fraction of liquid water in each layer (-)
  integer(i4b)    :: mLayerVolHtCapBulk              = 160 ! volumetric heat capacity in each layer (J m-3 K-1)
  integer(i4b)    :: mLayerTcrit                     = 161 ! critical soil temperature above which all water is unfrozen (K)
  integer(i4b)    :: mLayerdTheta_dTk                = 162 ! derivative in volumetric liquid water content wrt temperature (K-1)
  integer(i4b)    :: mLayerThermalC                  = 163 ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
  integer(i4b)    :: mLayerRadCondFlux               = 164 ! temporal derivative in energy from radiative and conductive flux (J m-2 s-1)
  integer(i4b)    :: mLayerMeltFreeze                = 165 ! rate of ice content change from melt/freeze in each layer (kg m-3 s-1)
  integer(i4b)    :: mLayerInfilFreeze               = 166 ! rate of ice content change by freezing infiltrating flux (kg m-3 s-1)
  integer(i4b)    :: mLayerSatHydCond                = 167 ! saturated hydraulic conductivity in each layer (m s-1)
  integer(i4b)    :: mLayerSatHydCondMP              = 168 ! saturated hydraulic conductivity of macropores in each layer (m s-1)
  integer(i4b)    :: mLayerMatricHead                = 169 ! matric head of water in the soil (m)
  integer(i4b)    :: mLayerdTheta_dPsi               = 170 ! derivative in the soil water characteristic (m-1)
  integer(i4b)    :: mLayerdPsi_dTheta               = 171 ! derivative in the soil water characteristic (m)
  integer(i4b)    :: mLayerThetaResid                = 172 ! residual volumetric water content in each snow layer (-)
  integer(i4b)    :: mLayerPoreSpace                 = 173 ! total pore space in each snow layer (-)
  integer(i4b)    :: mLayerCompress                  = 174 ! change in volumetric water content due to compression of soil (-)
  integer(i4b)    :: mLayerTranspireLim              = 175 ! soil moist & veg limit on transpiration for each layer (-)
  integer(i4b)    :: mLayerInitTranspire             = 176 ! transpiration loss from each soil layer at the start of the step (kg m-2 s-1)
  integer(i4b)    :: mLayerTranspire                 = 177 ! transpiration loss from each soil layer (kg m-2 s-1)
  integer(i4b)    :: mLayerInitQMacropore            = 178 ! liquid flux from micropores to macropores at the start-of-step (m s-1)
  integer(i4b)    :: mLayerQMacropore                = 179 ! liquid flux from micropores to macropores (m s-1)
  integer(i4b)    :: mLayerInitBaseflow              = 180 ! baseflow from each soil layer at the start of the time step (m s-1)
  integer(i4b)    :: mLayerBaseflow                  = 181 ! baseflow from each soil layer (m s-1)
  integer(i4b)    :: mLayerColumnInflow              = 182 ! total inflow to each layer in a given soil column (m3 s-1)
  integer(i4b)    :: mLayerColumnOutflow             = 183 ! total outflow from each layer in a given soil column (m3 s-1)
  ! define variables at the interface of each layer
  integer(i4b)    :: iLayerHeight                    = 184 ! height of the layer interface; top of soil = 0 (m)
  integer(i4b)    :: iLayerThermalC                  = 185 ! thermal conductivity at the interface of each layer (W m-1 K-1)
  integer(i4b)    :: iLayerConductiveFlux            = 186 ! conductive energy flux at layer interfaces at end of time step (W m-2)
  integer(i4b)    :: iLayerAdvectiveFlux             = 187 ! advective energy flux at layer interfaces at end of time step (W m-2)
  integer(i4b)    :: iLayerInitNrgFlux               = 188 ! energy flux at layer interfaces at the start of the time step (W m-2)
  integer(i4b)    :: iLayerNrgFlux                   = 189 ! energy flux at layer interfaces at the end of the time step (W m-2)
  integer(i4b)    :: iLayerSatHydCond                = 190 ! saturated hydraulic conductivity at each layer interface (m s-1)
  integer(i4b)    :: iLayerInitLiqFluxSnow           = 191 ! liquid flux at snow layer interfaces at the start of the time step (m s-1)
  integer(i4b)    :: iLayerInitLiqFluxSoil           = 192 ! liquid flux at soil layer interfaces at the start of the time step (m s-1)
  integer(i4b)    :: iLayerInitFluxReversal          = 193 ! start of step liquid flux at soil layer interfaces from impedance (m s-1)
  integer(i4b)    :: iLayerLiqFluxSnow               = 194 ! liquid flux at snow layer interfaces at the end of the time step (m s-1)
  integer(i4b)    :: iLayerLiqFluxSoil               = 195 ! liquid flux at soil layer interfaces at the end of the time step (m s-1)
  integer(i4b)    :: iLayerFluxReversal              = 196 ! end of step liquid flux at soil layer interfaces from impedance (m s-1)
  ! define variables for time stepping
  integer(i4b)    :: dt_init                         = 197 ! length of initial time step at start of next data interval (s)
  ! define derived variable
  integer(i4b)    :: scalarVGn_m                     = 198 ! van Genuchten "m" parameter (-)
  integer(i4b)    :: scalarKappa                     = 199 ! constant in the freezing curve function (m K-1)
  integer(i4b)    :: scalarVolHtCap_air              = 200 ! volumetric heat capacity air         (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_ice              = 201 ! volumetric heat capacity ice         (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_soil             = 202 ! volumetric heat capacity dry soil    (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_water            = 203 ! volumetric heat capacity liquid wat  (J m-3 K-1)
  integer(i4b)    :: scalarLambda_drysoil            = 204 ! thermal conductivity of dry soil     (W m-1)
  integer(i4b)    :: scalarLambda_wetsoil            = 205 ! thermal conductivity of wet soil     (W m-1)
  integer(i4b)    :: scalarVolLatHt_fus              = 206 ! volumetric latent heat of fusion     (J m-3)
  integer(i4b)    :: scalarAquiferRootFrac           = 207 ! fraction of roots below the soil profile (-)
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
  integer(i4b)    :: basin__ColumnOutflow            =  3 ! outflow from all "outlet" HRUs (those with no downstream HRU)
  integer(i4b)    :: basin__AquiferStorage           =  4 ! aquifer storage (m s-1)
  integer(i4b)    :: basin__AquiferRecharge          =  5 ! recharge to the aquifer (m s-1)
  integer(i4b)    :: basin__AquiferBaseflow          =  6 ! baseflow from the aquifer (m s-1)
  integer(i4b)    :: basin__AquiferTranspire         =  7 ! transpiration from the aquifer (m s-1)
  ! define variables for runoff
  integer(i4b)    :: routingRunoffFuture             =  8 ! runoff in future time steps (m s-1)
  integer(i4b)    :: routingFractionFuture           =  9 ! fraction of runoff in future time steps (-)
  integer(i4b)    :: averageInstantRunoff            = 10 ! instantaneous runoff (m s-1)
  integer(i4b)    :: averageRoutedRunoff             = 11 ! routed runoff (m s-1)
 endtype iLook_bvar

 ! ***********************************************************************************************************
 ! (X) define data structures and maximum number of variables of each type
 ! ***********************************************************************************************************
 ! define look-up structures
 type(iLook_decision),public,parameter :: iLookDECISIONS=iLook_decision(  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,&
                                                                         21, 22, 23, 24, 25, 26, 27, 28, 29, 30,&
                                                                         31, 32, 33, 34, 35)
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
                                                                        101,102,103,104,105,106,107,108,109,110,&
                                                                        111,112,113,114,115,116,117,118,119,120,&
                                                                        121,122,123,124,125,126,127,128,129,130,&
                                                                        131,132,133,134,135,136,137,138,139,140,&
                                                                        141,142,143,144,145,146)
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
                                                                        171,172,173,174,175,176,177,178,179,180,&
                                                                        181,182,183,184,185,186,187,188,189,190,&
                                                                        191,192,193,194,195,196,197,198,199,200,&
                                                                        201,202,203,204,205,206,207)
 type(iLook_index),   public,parameter :: iLookINDEX    =ilook_index   (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10)
 type(iLook_bpar),    public,parameter :: iLookBPAR     =ilook_bpar    (  1,  2,  3,  4,  5)
 type(iLook_bvar),    public,parameter :: iLookBVAR     =ilook_bvar    (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11)
 ! define maximum number of variables of each type
 integer(i4b),parameter,public :: maxvarDecisions= 35
 integer(i4b),parameter,public :: maxvarTime     = 5
 integer(i4b),parameter,public :: maxvarForc     = 8
 integer(i4b),parameter,public :: maxvarAttr     = 7
 integer(i4b),parameter,public :: maxvarType     = 5
 integer(i4b),parameter,public :: maxvarMpar     = 146
 integer(i4b),parameter,public :: maxvarMvar     = 207
 integer(i4b),parameter,public :: maxvarIndx     = 10
 integer(i4b),parameter,public :: maxvarBpar     = 5
 integer(i4b),parameter,public :: maxvarBvar     = 11


END MODULE var_lookup

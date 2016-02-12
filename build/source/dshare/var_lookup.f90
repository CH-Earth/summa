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
  integer(i4b)    :: bbNumerics       = 11 ! Ball-Berry: iterative numerical solution method
  integer(i4b)    :: bbAssimFnc       = 12 ! Ball-Berry: controls on carbon assimilation
  integer(i4b)    :: bbCanIntg8       = 13 ! Ball-Berry: scaling of photosynthesis from the leaf to the canopy
  integer(i4b)    :: num_method       = 14 ! choice of numerical method
  integer(i4b)    :: fDerivMeth       = 15 ! method used to calculate flux derivatives
  integer(i4b)    :: LAI_method       = 16 ! method used to determine LAI and SAI
  integer(i4b)    :: cIntercept       = 17 ! choice of parameterization for canopy interception
  integer(i4b)    :: f_Richards       = 18 ! form of richards' equation
  integer(i4b)    :: groundwatr       = 19 ! choice of groundwater parameterization
  integer(i4b)    :: hc_profile       = 20 ! choice of hydraulic conductivity profile
  integer(i4b)    :: bcUpprTdyn       = 21 ! type of upper boundary condition for thermodynamics
  integer(i4b)    :: bcLowrTdyn       = 22 ! type of lower boundary condition for thermodynamics
  integer(i4b)    :: bcUpprSoiH       = 23 ! type of upper boundary condition for soil hydrology
  integer(i4b)    :: bcLowrSoiH       = 24 ! type of lower boundary condition for soil hydrology
  integer(i4b)    :: veg_traits       = 25 ! choice of parameterization for vegetation roughness length and displacement height
  integer(i4b)    :: rootProfil       = 26 ! choice of parameterization for the rooting profile
  integer(i4b)    :: canopyEmis       = 27 ! choice of parameterization for canopy emissivity
  integer(i4b)    :: snowIncept       = 28 ! choice of parameterization for snow interception
  integer(i4b)    :: windPrfile       = 29 ! choice of canopy wind profile
  integer(i4b)    :: astability       = 30 ! choice of stability function
  integer(i4b)    :: canopySrad       = 31 ! choice of method for canopy shortwave radiation
  integer(i4b)    :: alb_method       = 32 ! choice of albedo representation
  integer(i4b)    :: snowLayers       = 33 ! choice of method to combine and sub-divide snow layers
  integer(i4b)    :: compaction       = 34 ! choice of compaction routine
  integer(i4b)    :: thCondSnow       = 35 ! choice of thermal conductivity representation for snow
  integer(i4b)    :: thCondSoil       = 36 ! choice of thermal conductivity representation for soil
  integer(i4b)    :: spatial_gw       = 37 ! choice of method for spatial representation of groundwater
  integer(i4b)    :: subRouting       = 38 ! choice of method for sub-grid routing
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
  integer(i4b)    :: vcmax25_canopyTop         = 53  ! potential carboxylation rate at 25 degrees C at the canopy top (umol co2 m-2 s-1)
  integer(i4b)    :: vcmax_qFac                = 54  ! factor in the q10 function defining temperature controls on vcmax (-)
  integer(i4b)    :: vcmax_Ha                  = 55  ! activation energy in the vcmax function (J mol-1)
  integer(i4b)    :: vcmax_Hd                  = 56  ! deactivation energy in the vcmax function (J mol-1)
  integer(i4b)    :: vcmax_Sv                  = 57  ! entropy term in the vcmax function (J mol-1 K-1)
  integer(i4b)    :: vcmax_Kn                  = 58  ! foliage nitrogen decay coefficient (-)
  integer(i4b)    :: jmax25_scale              = 59  ! scaling factor to relate jmax25 to vcmax25 (-)
  integer(i4b)    :: jmax_Ha                   = 60  ! activation energy in the jmax function (J mol-1)
  integer(i4b)    :: jmax_Hd                   = 61  ! deactivation energy in the jmax function (J mol-1)
  integer(i4b)    :: jmax_Sv                   = 62  ! entropy term in the jmax function (J mol-1 K-1)
  integer(i4b)    :: fractionJ                 = 63  ! fraction of light lost by other than the chloroplast lamellae (-)
  integer(i4b)    :: quantamYield              = 64  ! quantam yield (mol e mol-1 q)
  integer(i4b)    :: vpScaleFactor             = 65  ! vapor pressure scaling factor in stomatal conductance function (Pa)
  integer(i4b)    :: cond2photo_slope          = 66  ! slope of conductance-photosynthesis relationship (-)
  integer(i4b)    :: minStomatalConductance    = 67  ! minimum stomatal conductance (umol H2O m-2 s-1)
  ! vegetation properties
  integer(i4b)    :: winterSAI                 = 68  ! stem area index prior to the start of the growing season (m2 m-2)
  integer(i4b)    :: summerLAI                 = 69  ! maximum leaf area index at the peak of the growing season (m2 m-2)
  integer(i4b)    :: rootScaleFactor1          = 70  ! 1st scaling factor (a) in Y = 1 - 0.5*( exp(-aZ) + exp(-bZ) )  (m-1)
  integer(i4b)    :: rootScaleFactor2          = 71  ! 2nd scaling factor (b) in Y = 1 - 0.5*( exp(-aZ) + exp(-bZ) )  (m-1)
  integer(i4b)    :: rootingDepth              = 72  ! rooting depth (m)
  integer(i4b)    :: rootDistExp               = 73  ! exponent controlling the vertical distribution of root density (-)
  integer(i4b)    :: plantWiltPsi              = 74  ! matric head at wilting point (m)
  integer(i4b)    :: soilStressParam           = 75  ! parameter in the exponential soil stress function
  integer(i4b)    :: critSoilWilting           = 76  ! critical vol. liq. water content when plants are wilting (-)
  integer(i4b)    :: critSoilTranspire         = 77  ! critical vol. liq. water content when transpiration is limited (-)
  integer(i4b)    :: critAquiferTranspire      = 78  ! critical aquifer storage value when transpiration is limited (m)
  integer(i4b)    :: minStomatalResistance     = 79  ! minimum canopy resistance (s m-1)
  integer(i4b)    :: leafDimension             = 80  ! characteristic leaf dimension (m)
  integer(i4b)    :: heightCanopyTop           = 81  ! height of top of the vegetation canopy above ground surface (m)
  integer(i4b)    :: heightCanopyBottom        = 82  ! height of bottom of the vegetation canopy above ground surface (m)
  integer(i4b)    :: specificHeatVeg           = 83  ! specific heat of vegetation (J kg-1 K-1)
  integer(i4b)    :: maxMassVegetation         = 84  ! maximum mass of vegetation (full foliage) (kg m-2)
  integer(i4b)    :: throughfallScaleSnow      = 85  ! scaling factor for throughfall (snow) (-)
  integer(i4b)    :: throughfallScaleRain      = 86  ! scaling factor for throughfall (rain) (-)
  integer(i4b)    :: refInterceptCapSnow       = 87  ! reference canopy interception capacity per unit leaf area (snow) (kg m-2)
  integer(i4b)    :: refInterceptCapRain       = 88  ! canopy interception capacity per unit leaf area (rain) (kg m-2)
  integer(i4b)    :: snowUnloadingCoeff        = 89  ! time constant for unloading of snow from the forest canopy (s-1)
  integer(i4b)    :: canopyDrainageCoeff       = 90  ! time constant for drainage of liquid water from the forest canopy (s-1)
  integer(i4b)    :: ratioDrip2Unloading       = 91  ! ratio of canopy drip to unloading of snow from the forest canopy (-)
  integer(i4b)    :: canopyWettingFactor       = 92  ! maximum wetted fraction of the canopy (-)
  integer(i4b)    :: canopyWettingExp          = 93  ! exponent in canopy wetting function (-)
  ! soil properties
  integer(i4b)    :: soil_dens_intr            = 94  ! intrinsic soil density (kg m-3)
  integer(i4b)    :: thCond_soil               = 95  ! thermal conductivity of soil (W m-1 K-1)
  integer(i4b)    :: frac_sand                 = 96  ! fraction of sand (-)
  integer(i4b)    :: frac_silt                 = 97  ! fraction of silt (-)
  integer(i4b)    :: frac_clay                 = 98  ! fraction of clay (-)
  integer(i4b)    :: fieldCapacity             = 99  ! field capacity (-)
  integer(i4b)    :: wettingFrontSuction       = 100 ! Green-Ampt wetting front suction (m)
  integer(i4b)    :: theta_mp                  = 101 ! volumetric liquid water content when macropore flow begins (-)
  integer(i4b)    :: theta_sat                 = 102 ! porosity (-)
  integer(i4b)    :: theta_res                 = 103 ! volumetric residual water content (-)
  integer(i4b)    :: vGn_alpha                 = 104 ! van Genuchten "alpha" parameter (m-1)
  integer(i4b)    :: vGn_n                     = 105 ! van Genuchten "n" parameter (-)
  integer(i4b)    :: mpExp                     = 106 ! empirical exponent in macropore flow equation (-)
  integer(i4b)    :: k_soil                    = 107 ! hydraulic conductivity of soil (m s-1)
  integer(i4b)    :: k_macropore               = 108 ! saturated hydraulic conductivity for macropores (m s-1)
  integer(i4b)    :: kAnisotropic              = 109 ! anisotropy factor for lateral hydraulic conductivity (-)
  integer(i4b)    :: zScale_TOPMODEL           = 110 ! TOPMODEL scaling factor used in lower boundary condition for soil (m)
  integer(i4b)    :: compactedDepth            = 111 ! depth where k_soil reaches the compacted value given by CH78 (m)
  integer(i4b)    :: aquiferScaleFactor        = 112 ! scaling factor for aquifer storage in the big bucket (m)
  integer(i4b)    :: aquiferBaseflowExp        = 113 ! baseflow exponent (-)
  integer(i4b)    :: qSurfScale                = 114 ! scaling factor in the surface runoff parameterization (-)
  integer(i4b)    :: specificYield             = 115 ! specific yield (-)
  integer(i4b)    :: specificStorage           = 116 ! specific storage coefficient (m-1)
  integer(i4b)    :: f_impede                  = 117 ! ice impedence factor (-)
  integer(i4b)    :: soilIceScale              = 118 ! scaling factor for depth of soil ice, used to get frozen fraction (m)
  integer(i4b)    :: soilIceCV                 = 119 ! CV of depth of soil ice, used to get frozen fraction (-)
  ! algorithmic control parameters
  integer(i4b)    :: minwind                   = 120 ! minimum wind speed (m s-1)
  integer(i4b)    :: minstep                   = 121 ! minimum length of the time step
  integer(i4b)    :: maxstep                   = 122 ! maximum length of the time step
  integer(i4b)    :: wimplicit                 = 123 ! weight assigned to the start-of-step fluxes
  integer(i4b)    :: maxiter                   = 124 ! maximum number of iteration
  integer(i4b)    :: relConvTol_liquid         = 125 ! relative convergence tolerance for vol frac liq water (-)
  integer(i4b)    :: absConvTol_liquid         = 126 ! absolute convergence tolerance for vol frac liq water (-)
  integer(i4b)    :: relConvTol_matric         = 127 ! relative convergence tolerance for matric head (-)
  integer(i4b)    :: absConvTol_matric         = 128 ! absolute convergence tolerance for matric head (m)
  integer(i4b)    :: relConvTol_energy         = 129 ! relative convergence tolerance for energy (-)
  integer(i4b)    :: absConvTol_energy         = 130 ! absolute convergence tolerance for energy (J m-3)
  integer(i4b)    :: relConvTol_aquifr         = 131 ! relative convergence tolerance for aquifer storage (-)
  integer(i4b)    :: absConvTol_aquifr         = 132 ! absolute convergence tolerance for aquifer storage (J m-3)
  integer(i4b)    :: zmin                      = 133 ! minimum layer depth (m)
  integer(i4b)    :: zmax                      = 134 ! maximum layer depth (m)
  integer(i4b)    :: zminLayer1                = 135 ! minimum layer depth for the 1st (top) layer (m)
  integer(i4b)    :: zminLayer2                = 136 ! minimum layer depth for the 2nd layer (m)
  integer(i4b)    :: zminLayer3                = 137 ! minimum layer depth for the 3rd layer (m)
  integer(i4b)    :: zminLayer4                = 138 ! minimum layer depth for the 4th layer (m)
  integer(i4b)    :: zminLayer5                = 139 ! minimum layer depth for the 5th (bottom) layer (m)
  integer(i4b)    :: zmaxLayer1_lower          = 140 ! maximum layer depth for the 1st (top) layer when only 1 layer (m)
  integer(i4b)    :: zmaxLayer2_lower          = 141 ! maximum layer depth for the 2nd layer when only 2 layers (m)
  integer(i4b)    :: zmaxLayer3_lower          = 142 ! maximum layer depth for the 3rd layer when only 3 layers (m)
  integer(i4b)    :: zmaxLayer4_lower          = 143 ! maximum layer depth for the 4th layer when only 4 layers (m)
  integer(i4b)    :: zmaxLayer1_upper          = 144 ! maximum layer depth for the 1st (top) layer when > 1 layer (m)
  integer(i4b)    :: zmaxLayer2_upper          = 145 ! maximum layer depth for the 2nd layer when > 2 layers (m)
  integer(i4b)    :: zmaxLayer3_upper          = 146 ! maximum layer depth for the 3rd layer when > 3 layers (m)
  integer(i4b)    :: zmaxLayer4_upper          = 147 ! maximum layer depth for the 4th layer when > 4 layers (m)
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
  integer(i4b)    :: scalarIntercellularCO2Sunlit    = 125 ! carbon dioxide partial pressure of leaf interior (sunlit leaves) (Pa)
  integer(i4b)    :: scalarIntercellularCO2Shaded    = 126 ! carbon dioxide partial pressure of leaf interior (shaded leaves) (Pa)
  ! define vegetation variables -- canopy water
  integer(i4b)    :: scalarCanopyWetFraction         = 127 ! fraction of canopy that is wet
  integer(i4b)    :: scalarGroundSnowFraction        = 128 ! fraction of ground that is covered with snow (-)
  integer(i4b)    :: scalarThroughfallSnow           = 129 ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
  integer(i4b)    :: scalarThroughfallRain           = 130 ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
  integer(i4b)    :: scalarCanopySnowUnloading       = 131 ! unloading of snow from the vegetion canopy (kg m-2 s-1)
  integer(i4b)    :: scalarCanopyLiqDrainage         = 132 ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
  integer(i4b)    :: scalarCanopyMeltFreeze          = 133 ! melt/freeze of water stored in the canopy (kg m-2 s-1)
  ! define scalar variables -- soil and aquifer fluxes
  integer(i4b)    :: scalarRainPlusMelt              = 134 ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  integer(i4b)    :: scalarInfilArea                 = 135 ! fraction of unfrozen area where water can infiltrate (-)
  integer(i4b)    :: scalarFrozenArea                = 136 ! fraction of area that is considered impermeable due to soil ice (-)
  integer(i4b)    :: scalarInfiltration              = 137 ! infiltration of water into the soil profile (m s-1)
  integer(i4b)    :: scalarExfiltration              = 138 ! exfiltration of water from the top of the soil profile (m s-1)
  integer(i4b)    :: scalarSurfaceRunoff             = 139 ! surface runoff (m s-1)
  integer(i4b)    :: scalarInitAquiferRecharge       = 140 ! recharge to the aquifer at the start of the step (m s-1)
  integer(i4b)    :: scalarAquiferRecharge           = 141 ! recharge to the aquifer (m s-1)
  integer(i4b)    :: scalarInitAquiferTranspire      = 142 ! transpiration from the aquifer at the start of the step (m s-1)
  integer(i4b)    :: scalarAquiferTranspire          = 143 ! transpiration from the aquifer (m s-1)
  integer(i4b)    :: scalarInitAquiferBaseflow       = 144 ! baseflow from the aquifer at the start of the step (m s-1)
  integer(i4b)    :: scalarAquiferBaseflow           = 145 ! baseflow from the aquifer (m s-1)
  ! scalar variables -- sub-step average fluxes for the soil zone
  integer(i4b)    :: scalarSoilInflux                = 146 ! influx of water at the top of the soil profile (m s-1)
  integer(i4b)    :: scalarSoilCompress              = 147 ! change in total soil storage due to compression of the soil matrix (kg m-2)
  integer(i4b)    :: scalarSoilBaseflow              = 148 ! sub-step average: total baseflow from throughout the soil profile (m s-1)
  integer(i4b)    :: scalarSoilDrainage              = 149 ! sub-step average: drainage from the bottom of the soil profile (m s-1)
  integer(i4b)    :: scalarSoilTranspiration         = 150 ! sub-step average: total transpiration from the soil (m s-1)
  ! define scalar variables -- mass balance check
  integer(i4b)    :: scalarSoilWatBalError           = 151 ! error in the total soil water balance (kg m-2)
  integer(i4b)    :: scalarAquiferBalError           = 152 ! error in the aquifer water balance (kg m-2)
  integer(i4b)    :: scalarTotalSoilLiq              = 153 ! total mass of liquid water in the soil (kg m-2)
  integer(i4b)    :: scalarTotalSoilIce              = 154 ! total mass of ice in the soil (kg m-2)
  ! define variables at the mid-point of each layer -- domain geometry
  integer(i4b)    :: mLayerDepth                     = 155 ! depth of each layer (m)
  integer(i4b)    :: mLayerHeight                    = 156 ! height at the mid-point of each layer (m)
  integer(i4b)    :: mLayerRootDensity               = 157 ! fraction of roots in each soil layer (-)
  ! define variables at the mid-point of each layer -- coupled energy and mass
  integer(i4b)    :: mLayerTemp                      = 158 ! temperature of each layer (K)
  integer(i4b)    :: mLayerVolFracAir                = 159 ! volumetric fraction of air in each layer (-)
  integer(i4b)    :: mLayerVolFracIce                = 160 ! volumetric fraction of ice water in each layer (-)
  integer(i4b)    :: mLayerVolFracLiq                = 161 ! volumetric fraction of liquid water in each layer (-)
  integer(i4b)    :: mLayerVolHtCapBulk              = 162 ! volumetric heat capacity in each layer (J m-3 K-1)
  integer(i4b)    :: mLayerTcrit                     = 163 ! critical soil temperature above which all water is unfrozen (K)
  integer(i4b)    :: mLayerdTheta_dTk                = 164 ! derivative in volumetric liquid water content wrt temperature (K-1)
  integer(i4b)    :: mLayerThermalC                  = 165 ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
  integer(i4b)    :: mLayerRadCondFlux               = 166 ! temporal derivative in energy from radiative and conductive flux (J m-2 s-1)
  integer(i4b)    :: mLayerMeltFreeze                = 167 ! rate of ice content change from melt/freeze in each layer (kg m-3 s-1)
  integer(i4b)    :: mLayerInfilFreeze               = 168 ! rate of ice content change by freezing infiltrating flux (kg m-3 s-1)
  integer(i4b)    :: mLayerSatHydCond                = 169 ! saturated hydraulic conductivity in each layer (m s-1)
  integer(i4b)    :: mLayerSatHydCondMP              = 170 ! saturated hydraulic conductivity of macropores in each layer (m s-1)
  integer(i4b)    :: mLayerMatricHead                = 171 ! matric head of water in the soil (m)
  integer(i4b)    :: mLayerdTheta_dPsi               = 172 ! derivative in the soil water characteristic (m-1)
  integer(i4b)    :: mLayerdPsi_dTheta               = 173 ! derivative in the soil water characteristic (m)
  integer(i4b)    :: mLayerThetaResid                = 174 ! residual volumetric water content in each snow layer (-)
  integer(i4b)    :: mLayerPoreSpace                 = 175 ! total pore space in each snow layer (-)
  integer(i4b)    :: mLayerCompress                  = 176 ! change in volumetric water content due to compression of soil (-)
  integer(i4b)    :: mLayerTranspireLim              = 177 ! soil moist & veg limit on transpiration for each layer (-)
  integer(i4b)    :: mLayerInitTranspire             = 178 ! transpiration loss from each soil layer at the start of the step (kg m-2 s-1)
  integer(i4b)    :: mLayerTranspire                 = 179 ! transpiration loss from each soil layer (kg m-2 s-1)
  integer(i4b)    :: mLayerInitQMacropore            = 180 ! liquid flux from micropores to macropores at the start-of-step (m s-1)
  integer(i4b)    :: mLayerQMacropore                = 181 ! liquid flux from micropores to macropores (m s-1)
  integer(i4b)    :: mLayerInitBaseflow              = 182 ! baseflow from each soil layer at the start of the time step (m s-1)
  integer(i4b)    :: mLayerBaseflow                  = 183 ! baseflow from each soil layer (m s-1)
  integer(i4b)    :: mLayerColumnInflow              = 184 ! total inflow to each layer in a given soil column (m3 s-1)
  integer(i4b)    :: mLayerColumnOutflow             = 185 ! total outflow from each layer in a given soil column (m3 s-1)
  ! define variables at the interface of each layer
  integer(i4b)    :: iLayerHeight                    = 186 ! height of the layer interface; top of soil = 0 (m)
  integer(i4b)    :: iLayerThermalC                  = 187 ! thermal conductivity at the interface of each layer (W m-1 K-1)
  integer(i4b)    :: iLayerConductiveFlux            = 188 ! conductive energy flux at layer interfaces at end of time step (W m-2)
  integer(i4b)    :: iLayerAdvectiveFlux             = 189 ! advective energy flux at layer interfaces at end of time step (W m-2)
  integer(i4b)    :: iLayerInitNrgFlux               = 190 ! energy flux at layer interfaces at the start of the time step (W m-2)
  integer(i4b)    :: iLayerNrgFlux                   = 191 ! energy flux at layer interfaces at the end of the time step (W m-2)
  integer(i4b)    :: iLayerSatHydCond                = 192 ! saturated hydraulic conductivity at each layer interface (m s-1)
  integer(i4b)    :: iLayerInitLiqFluxSnow           = 193 ! liquid flux at snow layer interfaces at the start of the time step (m s-1)
  integer(i4b)    :: iLayerInitLiqFluxSoil           = 194 ! liquid flux at soil layer interfaces at the start of the time step (m s-1)
  integer(i4b)    :: iLayerInitFluxReversal          = 195 ! start of step liquid flux at soil layer interfaces from impedance (m s-1)
  integer(i4b)    :: iLayerLiqFluxSnow               = 196 ! liquid flux at snow layer interfaces at the end of the time step (m s-1)
  integer(i4b)    :: iLayerLiqFluxSoil               = 197 ! liquid flux at soil layer interfaces at the end of the time step (m s-1)
  integer(i4b)    :: iLayerFluxReversal              = 198 ! end of step liquid flux at soil layer interfaces from impedance (m s-1)
  ! define variables for time stepping
  integer(i4b)    :: dt_init                         = 199 ! length of initial time step at start of next data interval (s)
  ! define derived variable
  integer(i4b)    :: scalarVGn_m                     = 200 ! van Genuchten "m" parameter (-)
  integer(i4b)    :: scalarKappa                     = 201 ! constant in the freezing curve function (m K-1)
  integer(i4b)    :: scalarVolHtCap_air              = 202 ! volumetric heat capacity air         (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_ice              = 203 ! volumetric heat capacity ice         (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_soil             = 204 ! volumetric heat capacity dry soil    (J m-3 K-1)
  integer(i4b)    :: scalarVolHtCap_water            = 205 ! volumetric heat capacity liquid wat  (J m-3 K-1)
  integer(i4b)    :: scalarLambda_drysoil            = 206 ! thermal conductivity of dry soil     (W m-1)
  integer(i4b)    :: scalarLambda_wetsoil            = 207 ! thermal conductivity of wet soil     (W m-1)
  integer(i4b)    :: scalarVolLatHt_fus              = 208 ! volumetric latent heat of fusion     (J m-3)
  integer(i4b)    :: scalarAquiferRootFrac           = 209 ! fraction of roots below the soil profile (-)
 endtype iLook_mvar

 ! ***********************************************************************************************************
 ! (7) define model fluxes
 ! ***********************************************************************************************************
 type, public :: iLook_flux
  ! net energy and mass fluxes for the vegetation domain
  integer(i4b)    :: canairNetNrgFlux             ! net energy flux for the canopy air space (W m-2)
  integer(i4b)    :: canopyNetNrgFlux             ! net energy flux for the vegetation canopy (W m-2)
  integer(i4b)    :: groundNetNrgFlux             ! net energy flux for the ground surface (W m-2)
  integer(i4b)    :: canopyNetLiqFlux             ! net liquid water flux for the vegetation canopy (kg m-2 s-1)
  ! liquid water fluxes associated with evapotranspiration
  integer(i4b)    :: scalarCanopyTranspiration    ! canopy transpiration (kg m-2 s-1)
  integer(i4b)    :: scalarCanopyEvaporation      ! canopy evaporation/condensation (kg m-2 s-1)
  integer(i4b)    :: scalarGroundEvaporation      ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
  ! energy fluxes and derivatives for the snow and soil domains
  integer(i4b)    :: ssdNetNrgFlux                ! net energy flux for each layer in the snow+soil domain (J m-3 s-1)
  ! liquid water fluxes and derivatives for the soil domain
  integer(i4b)    :: xMaxInfilRate                ! maximum infiltration rate (m s-1)
  integer(i4b)    :: scalarSurfaceInfiltration    ! surface infiltration rate (m s-1) -- only computed for iter==1
  integer(i4b)    :: mLayerTranspire              ! transpiration loss from each soil layer (m s-1)
  integer(i4b)    :: mLayerHydCond                ! hydraulic conductivity in each soil layer (m s-1)
  integer(i4b)    :: snowNetLiqFlux               ! net liquid water flux for each snow layer (s-1)
  integer(i4b)    :: soilNetLiqFlux               ! net liquid water flux for each soil layer (s-1)
 endtype iLook_flux

 ! ***********************************************************************************************************
 ! (8) define derivatives
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
 ! (8) define diagnostic variables
 ! ***********************************************************************************************************
 type, public :: iLook_diag
  integer(i4b)    :: fracLiqVeg                   ! fraction of liquid water on vegetation (-)
  integer(i4b)    :: fracLiqSnow                  ! fraction of liquid water in each snow layer (-)
  integer(i4b)    :: scalarInfilArea              ! fraction of unfrozen area where water can infiltrate (-)
  integer(i4b)    :: scalarFrozenArea             ! fraction of area that is considered impermeable due to soil ice (-)
  integer(i4b)    :: soilControl                  ! soil control on infiltration (-)
 endtype iLook_diag

 ! ***********************************************************************************************************
 ! (9) define model indices
 ! ***********************************************************************************************************
 type, public :: iLook_index
  ! number of state variables of different type
  integer(i4b)    :: nVegNrg               = 1    ! number of energy state variables for vegetation
  integer(i4b)    :: nVegMass              = 2    ! number of hydrology state variables for vegetation (mass of water)
  integer(i4b)    :: nVegState             = 3    ! number of vegetation state variables
  integer(i4b)    :: nNrgState             = 4    ! number of energy state variables
  integer(i4b)    :: nWatState             = 5    ! number of "total water" state variables (volumetric total water content) 
  integer(i4b)    :: nMatState             = 6    ! number of matric head state variables
  integer(i4b)    :: nMassState            = 7    ! number of hydrology state variables (mass of water)
  integer(i4b)    :: nState                = 8    ! total number of model state variables
  ! number of model layers, and layer indices
  integer(i4b)    :: nSnow                 = 9    ! number of snow layers
  integer(i4b)    :: nSoil                 = 10   ! number of soil layers
  integer(i4b)    :: nLayers               = 11   ! total number of layers
  integer(i4b)    :: layerType             = 12   ! type of layer (soil or snow)
  ! indices of model state variables
  integer(i4b)    :: ixCasNrg              = 13   ! index of the canopy air space state variable
  integer(i4b)    :: ixVegNrg              = 14   ! index of the canopy energy state variable
  integer(i4b)    :: ixVegWat              = 15   ! index of the canopy total water state variable
  integer(i4b)    :: ixTopNrg              = 16   ! index of the upper-most energy state variable in the snow-soil subdomain
  integer(i4b)    :: ixTopWat              = 17   ! index of the upper-most total water state variable in the snow-soil subdomain
  integer(i4b)    :: ixTopMat              = 18   ! index of the upper-most matric head state variable in the soil subdomain
  integer(i4b)    :: ixSnowSoilNrg         = 19   ! indices for energy state variables in the snow-soil subdomain
  integer(i4b)    :: ixSnowSoilWat         = 20   ! indices for total water state variables in the snow-soil subdomain
  integer(i4b)    :: ixSnowOnlyNrg         = 21   ! indices for energy state variables in the snow subdomain
  integer(i4b)    :: ixSnowOnlyWat         = 22   ! indices for total water state variables in the snow subdomain
  integer(i4b)    :: ixSoilOnlyNrg         = 23   ! indices for energy state variables in the soil subdomain
  integer(i4b)    :: ixSoilOnlyHyd         = 24   ! indices for hydrology state variables in the soil subdomain
  ! type of model state variables
  integer(i4b)    :: ixStateType           = 25   ! indices defining the type of the state (ixNrgState, ixWatState, ixMatState...)
  integer(i4b)    :: ixAllState            = 26   ! list of indices for all model state variables
  integer(i4b)    :: ixSoilState           = 27   ! list of indices for all soil layers
  integer(i4b)    :: ixLayerState          = 28   ! list of indices for all model layers
  integer(i4b)    :: ixNrgOnly             = 29   ! list of indices for all energy states
  integer(i4b)    :: ixWatOnly             = 30   ! list of indices for all "total water" state variables (volumetric total water content)
  integer(i4b)    :: ixMatOnly             = 31   ! list of indices for matric head state variables
  integer(i4b)    :: ixMassOnly            = 32   ! list of indices for hydrology state variables (mass of water)
  ! indices for the model output files
  integer(i4b)    :: midSnowStartIndex     = 33   ! start index of the midSnow vector for a given timestep
  integer(i4b)    :: midSoilStartIndex     = 34   ! start index of the midSoil vector for a given timestep
  integer(i4b)    :: midTotoStartIndex     = 35   ! start index of the midToto vector for a given timestep
  integer(i4b)    :: ifcSnowStartIndex     = 36   ! start index of the ifcSnow vector for a given timestep
  integer(i4b)    :: ifcSoilStartIndex     = 37   ! start index of the ifcSoil vector for a given timestep
  integer(i4b)    :: ifcTotoStartIndex     = 38   ! start index of the ifcToto vector for a given timestep
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
                                                                         31, 32, 33, 34, 35, 36, 37, 38)
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
                                                                        141,142,143,144,145,146,147)
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
                                                                        201,202,203,204,205,206,207,208,209)
 type(iLook_index),   public,parameter :: iLookINDEX    =ilook_index   (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11, 12, 13, 14, 15, 16, 17, 18, 19, 20,&
                                                                         21, 22, 23, 24, 25, 26, 27, 28, 29, 30,&
                                                                         31, 32, 33, 34, 35, 36, 37, 38)
 type(iLook_bpar),    public,parameter :: iLookBPAR     =ilook_bpar    (  1,  2,  3,  4,  5)
 type(iLook_bvar),    public,parameter :: iLookBVAR     =ilook_bvar    (  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,&
                                                                         11)
 ! define maximum number of variables of each type
 integer(i4b),parameter,public :: maxvarDecisions= 38
 integer(i4b),parameter,public :: maxvarTime     = 5
 integer(i4b),parameter,public :: maxvarForc     = 8
 integer(i4b),parameter,public :: maxvarAttr     = 7
 integer(i4b),parameter,public :: maxvarType     = 5
 integer(i4b),parameter,public :: maxvarMpar     = 147
 integer(i4b),parameter,public :: maxvarMvar     = 209
 integer(i4b),parameter,public :: maxvarIndx     = 38
 integer(i4b),parameter,public :: maxvarBpar     = 5
 integer(i4b),parameter,public :: maxvarBvar     = 11


END MODULE var_lookup

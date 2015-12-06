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

module get_ixname_module
! used to get the index of a named variable
USE nrtype                                          ! variable types, etc.
implicit none
private
public::get_ixdecisions
public::get_ixTime
public::get_ixAttr
public::get_ixType
public::get_ixForce
public::get_ixParam
public::get_ixMvar
public::get_ixIndex
public::get_ixBpar
public::get_ixBvar
contains

 ! *******************************************************************************************************************
 ! public function get_ixdecisions: get the index of the named variables for the model decisions
 ! *******************************************************************************************************************
 function get_ixdecisions(varName)
 USE var_lookup,only:iLookDECISIONS                  ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixdecisions         ! index of the named variable
 ! define local variables
 integer(i4b), parameter  :: imiss = -999            ! missing value
 ! get the index of the named variables
 select case(trim(varName))
  ! simulation options
  case('simulStart'      ); get_ixdecisions=iLookDECISIONS%simulStart  ! ( 1) simulation start time
  case('simulFinsh'      ); get_ixdecisions=iLookDECISIONS%simulFinsh  ! ( 2) simulation end time
  ! Noah-MP decisions
  case('soilCatTbl'      ); get_ixdecisions=iLookDECISIONS%soilCatTbl  ! ( 3) soil-category dateset
  case('vegeParTbl'      ); get_ixdecisions=iLookDECISIONS%vegeParTbl  ! ( 4) vegetation category dataset
  case('soilStress'      ); get_ixdecisions=iLookDECISIONS%soilStress  ! ( 5) choice of function for the soil moisture control on stomatal resistance
  case('stomResist'      ); get_ixdecisions=iLookDECISIONS%stomResist  ! ( 6) choice of function for stomatal resistance
  ! SUMMA decisions
  case('num_method'      ); get_ixdecisions=iLookDECISIONS%num_method  ! ( 7) choice of numerical method
  case('fDerivMeth'      ); get_ixdecisions=iLookDECISIONS%fDerivMeth  ! ( 8) choice of method to calculate flux derivatives
  case('LAI_method'      ); get_ixdecisions=iLookDECISIONS%LAI_method  ! ( 9) choice of method to determine LAI and SAI
  case('f_Richards'      ); get_ixdecisions=iLookDECISIONS%f_Richards  ! (10) form of Richards' equation
  case('groundwatr'      ); get_ixdecisions=iLookDECISIONS%groundwatr  ! (11) choice of groundwater parameterization
  case('hc_profile'      ); get_ixdecisions=iLookDECISIONS%hc_profile  ! (12) choice of hydraulic conductivity profile
  case('bcUpprTdyn'      ); get_ixdecisions=iLookDECISIONS%bcUpprTdyn  ! (13) type of upper boundary condition for thermodynamics
  case('bcLowrTdyn'      ); get_ixdecisions=iLookDECISIONS%bcLowrTdyn  ! (14) type of lower boundary condition for thermodynamics
  case('bcUpprSoiH'      ); get_ixdecisions=iLookDECISIONS%bcUpprSoiH  ! (15) type of upper boundary condition for soil hydrology
  case('bcLowrSoiH'      ); get_ixdecisions=iLookDECISIONS%bcLowrSoiH  ! (16) type of lower boundary condition for soil hydrology
  case('veg_traits'      ); get_ixdecisions=iLookDECISIONS%veg_traits  ! (17) choice of parameterization for vegetation roughness length and displacement height
  case('canopyEmis'      ); get_ixdecisions=iLookDECISIONS%canopyEmis  ! (18) choice of parameterization for canopy emissivity
  case('snowIncept'      ); get_ixdecisions=iLookDECISIONS%snowIncept  ! (19) choice of parameterization for snow interception
  case('windPrfile'      ); get_ixdecisions=iLookDECISIONS%windPrfile  ! (20) choice of canopy wind profile
  case('astability'      ); get_ixdecisions=iLookDECISIONS%astability  ! (21) choice of stability function
  case('compaction'      ); get_ixdecisions=iLookDECISIONS%compaction  ! (22) choice of compaction routine
  case('snowLayers'      ); get_ixdecisions=iLookDECISIONS%snowLayers  ! (23) choice of method to combine and sub-divide snow layers
  case('thCondSnow'      ); get_ixdecisions=iLookDECISIONS%thCondSnow  ! (24) choice of thermal conductivity representation for snow
  case('thCondSoil'      ); get_ixdecisions=iLookDECISIONS%thCondSoil  ! (25) choice of thermal conductivity representation for soil
  case('canopySrad'      ); get_ixdecisions=iLookDECISIONS%canopySrad  ! (26) choice of method for canopy shortwave radiation
  case('alb_method'      ); get_ixdecisions=iLookDECISIONS%alb_method  ! (27) choice of albedo representation
  case('spatial_gw'      ); get_ixdecisions=iLookDECISIONS%spatial_gw  ! (28) choice of method for spatial representation of groundwater
  case('subRouting'      ); get_ixdecisions=iLookDECISIONS%subRouting  ! (29) choice of method for sub-grid routing
  ! get to here if cannot find the variable
  case default
   get_ixdecisions = imiss
 endselect
 end function get_ixdecisions


 ! *******************************************************************************************************************
 ! public function get_ixtime: get the index of the named variables for the model time
 ! *******************************************************************************************************************
 function get_ixtime(varName)
 USE var_lookup,only:iLookTIME                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixtime              ! index of the named variable
 ! define local variables
 integer(i4b), parameter  :: imiss = -999            ! missing value
 ! get the index of the named variables
 select case(trim(varName))
  case('iyyy'            ); get_ixtime = iLookTIME%iyyy             ! year
  case('im'              ); get_ixtime = iLookTIME%im               ! month
  case('id'              ); get_ixtime = iLookTIME%id               ! day
  case('ih'              ); get_ixtime = iLookTIME%ih               ! hour
  case('imin'            ); get_ixtime = iLookTIME%imin             ! minute
  case('isec'            ); get_ixtime = iLookTIME%isec             ! second
  ! get to here if cannot find the variable
  case default
   get_ixtime = imiss
 endselect
 end function get_ixtime


 ! *******************************************************************************************************************
 ! public function get_ixforce: get the index of the named variables for the model forcing data
 ! *******************************************************************************************************************
 function get_ixforce(varName)
 USE var_lookup,only:iLookFORCE                      ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixforce             ! index of the named variable
 ! define local variables
 integer(i4b), parameter  :: imiss = -999            ! missing value
 ! get the index of the named variables
 select case(trim(varName))
  case('time'            ); get_ixforce = iLookFORCE%time             ! time since time reference       (s)
  case('pptrate'         ); get_ixforce = iLookFORCE%pptrate          ! precipitation rate              (kg m-2 s-1)
  case('airtemp'         ); get_ixforce = iLookFORCE%airtemp          ! air temperature                 (K)
  case('spechum'         ); get_ixforce = iLookFORCE%spechum          ! specific humidity               (g/g)
  case('windspd'         ); get_ixforce = iLookFORCE%windspd          ! windspeed                       (m/s)
  case('SWRadAtm'        ); get_ixforce = iLookFORCE%SWRadAtm         ! downwelling shortwave radiaiton (W m-2)
  case('LWRadAtm'        ); get_ixforce = iLookFORCE%LWRadAtm         ! downwelling longwave radiation  (W m-2)
  case('airpres'         ); get_ixforce = iLookFORCE%airpres          ! pressure                        (Pa)
  ! get to here if cannot find the variable
  case default
   get_ixforce = imiss
 endselect
 end function get_ixforce


 ! *******************************************************************************************************************
 ! public function get_ixAttr: get the index of the named variables for the site characteristics
 ! *******************************************************************************************************************
 function get_ixAttr(varName)
 USE var_lookup,only:iLookATTR                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixAttr              ! index of the named variable
 ! define local variables
 integer(i4b), parameter  :: imiss = -999            ! missing value
 ! get the index of the named variables
 select case(trim(varName))
  case('latitude'      ); get_ixAttr = iLookATTR%latitude       ! latitude (degrees north)
  case('longitude'     ); get_ixAttr = iLookATTR%longitude      ! longitude (degrees east)
  case('elevation'     ); get_ixAttr = iLookATTR%elevation      ! elevation (m)
  case('tan_slope'     ); get_ixAttr = iLookATTR%tan_slope      ! tan water table slope, taken as tan local ground surface slope (-)
  case('contourLength' ); get_ixAttr = iLookATTR%contourLength  ! length of contour at downslope edge of HRU (m)
  case('HRUarea'       ); get_ixAttr = iLookATTR%HRUarea        ! area of each HRU (m2)
  case('mHeight'       ); get_ixAttr = iLookATTR%mHeight        ! measurement height above bare ground (m)
  ! get to here if cannot find the variable
  case default
   get_ixAttr = imiss
 endselect
 end function get_ixAttr


 ! *******************************************************************************************************************
 ! public function get_ixType: get the index of the named variables for the local classification of veg, soil, etc.
 ! *******************************************************************************************************************
 function get_ixType(varName)
 USE var_lookup,only:iLookTYPE                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixType              ! index of the named variable
 ! define local variables
 integer(i4b), parameter  :: imiss = -999            ! missing value
 ! get the index of the named variables
 select case(trim(varName))
  case('hruIndex'       ); get_ixType = iLookTYPE%hruIndex           ! index defining HRU index
  case('vegTypeIndex'   ); get_ixType = iLookTYPE%vegTypeIndex       ! index defining vegetation type
  case('soilTypeIndex'  ); get_ixType = iLookTYPE%soilTypeIndex      ! index defining soil type
  case('slopeTypeIndex' ); get_ixType = iLookTYPE%slopeTypeIndex     ! index defining slope
  case('downHRUindex'   ); get_ixType = iLookTYPE%downHRUindex       ! index of downslope HRU (0 = basin outlet)
  ! get to here if cannot find the variable
  case default
   get_ixType = imiss
 endselect
 end function get_ixType


 ! *******************************************************************************************************************
 ! public function get_ixparam: get the index of the named variables for the model parameters
 ! *******************************************************************************************************************
 function get_ixparam(varName)
 USE var_lookup,only:iLookPARAM                      ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixparam             ! index of the named variable
 ! define local variables
 integer(i4b), parameter  :: imiss = -999            ! missing value
 ! get the index of the named variables
 select case(trim(varName))
  ! boundary conditions
  case('upperBoundHead'           ); get_ixparam = iLookPARAM%upperBoundHead         ! matric head of the upper boundary (m)
  case('lowerBoundHead'           ); get_ixparam = iLookPARAM%lowerBoundHead         ! matric head of the lower boundary (m)
  case('upperBoundTheta'          ); get_ixparam = iLookPARAM%upperBoundTheta        ! volumetric liquid water content at the upper boundary (-)
  case('lowerBoundTheta'          ); get_ixparam = iLookPARAM%lowerBoundTheta        ! volumetric liquid water content at the lower boundary (-)
  case('upperBoundTemp'           ); get_ixparam = iLookPARAM%upperBoundTemp         ! temperature of the upper boundary (K)
  case('lowerBoundTemp'           ); get_ixparam = iLookPARAM%lowerBoundTemp         ! temperature of the lower boundary (K)
  ! precipitation partitioning
  case('tempCritRain'             ); get_ixparam = iLookPARAM%tempCritRain           ! critical temperature where precipitation is rain (K)
  case('tempRangeTimestep'        ); get_ixparam = iLookPARAM%tempRangeTimestep      ! temperature range over the time step (K)
  case('frozenPrecipMultip'       ); get_ixparam = iLookPARAM%frozenPrecipMultip     ! frozen precipitation multiplier (-)
  ! freezing curve for snow
  case('snowfrz_scale'            ); get_ixparam = iLookPARAM%snowfrz_scale          ! scaling parameter for the freezing curve for snow (K-1)
  ! snow albedo
  case('albedoMax'                ); get_ixparam = iLookPARAM%albedoMax              ! maximum snow albedo for a single spectral band (-)
  case('albedoMinWinter'          ); get_ixparam = iLookPARAM%albedoMinWinter        ! minimum snow albedo during winter for a single spectral band (-)
  case('albedoMinSpring'          ); get_ixparam = iLookPARAM%albedoMinSpring        ! minimum snow albedo during spring for a single spectral band (-)
  case('albedoMaxVisible'         ); get_ixparam = iLookPARAM%albedoMaxVisible       ! maximum snow albedo in the visible part of the spectrum (-)
  case('albedoMinVisible'         ); get_ixparam = iLookPARAM%albedoMinVisible       ! minimum snow albedo in the visible part of the spectrum (-)
  case('albedoMaxNearIR'          ); get_ixparam = iLookPARAM%albedoMaxNearIR        ! maximum snow albedo in the near infra-red part of the spectrum (-)
  case('albedoMinNearIR'          ); get_ixparam = iLookPARAM%albedoMinNearIR        ! minimum snow albedo in the near infra-red part of the spectrum (-)
  case('albedoDecayRate'          ); get_ixparam = iLookPARAM%albedoDecayRate        ! albedo decay rate (s)
  case('albedoSootLoad'           ); get_ixparam = iLookPARAM%albedoSootLoad         ! soot load factor (-)
  case('albedoRefresh'            ); get_ixparam = iLookPARAM%albedoRefresh          ! critical mass necessary for albedo refreshment (kg m-2)
  ! radiation transfer
  case('radExt_snow'              ); get_ixparam = iLookPARAM%radExt_snow            ! extinction coefficient for radiation penetration within the snowpack (m-1)
  case('directScale'              ); get_ixparam = iLookPARAM%directScale            ! scaling factor for fractional driect radiaion parameterization (-)
  case('Frad_direct'              ); get_ixparam = iLookPARAM%Frad_direct            ! maximum fraction of direct radiation (-)
  case('Frad_vis'                 ); get_ixparam = iLookPARAM%Frad_vis               ! fraction of radiation in the visible part of the spectrum (-)
  ! new snow density
  case('newSnowDenMin'            ); get_ixparam = iLookPARAM%newSnowDenMin          ! minimum new snow density (kg m-3)
  case('newSnowDenMult'           ); get_ixparam = iLookPARAM%newSnowDenMult         ! multiplier for new snow density (kg m-3)
  case('newSnowDenScal'           ); get_ixparam = iLookPARAM%newSnowDenScal         ! scaling factor for new snow density (K)
  ! snow compaction
  case('densScalGrowth'           ); get_ixparam = iLookPARAM%densScalGrowth         ! density scaling factor for grain growth (kg-1 m3)
  case('tempScalGrowth'           ); get_ixparam = iLookPARAM%tempScalGrowth         ! temperature scaling factor for grain growth (K-1)
  case('grainGrowthRate'          ); get_ixparam = iLookPARAM%grainGrowthRate        ! rate of grain growth (s-1)
  case('densScalOvrbdn'           ); get_ixparam = iLookPARAM%densScalOvrbdn         ! density scaling factor for overburden pressure (kg-1 m3)
  case('tempScalOvrbdn'           ); get_ixparam = iLookPARAM%tempScalOvrbdn         ! temperature scaling factor for overburden pressure (K-1)
  case('base_visc'                ); get_ixparam = iLookPARAM%base_visc              ! viscosity coefficient at T=T_frz and snow density=0  (kg s m-2)
  ! water flow through snow
  case('Fcapil'                   ); get_ixparam = iLookPARAM%Fcapil                 ! capillary retention as a fraction of the total pore volume (-)
  case('k_snow'                   ); get_ixparam = iLookPARAM%k_snow                 ! hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
  case('mw_exp'                   ); get_ixparam = iLookPARAM%mw_exp                 ! exponent for meltwater flow (-)
  ! turbulent heat fluxes
  case('z0Snow'                   ); get_ixparam = iLookPARAM%z0Snow                 ! roughness length of snow (m)
  case('z0Soil'                   ); get_ixparam = iLookPARAM%z0Soil                 ! roughness length of bare soil below the canopy (m)
  case('z0Canopy'                 ); get_ixparam = iLookPARAM%z0Canopy               ! roughness length of the canopy (m)
  case('zpdFraction'              ); get_ixparam = iLookPARAM%zpdFraction            ! zero plane displacement / canopy height (-)
  case('critRichNumber'           ); get_ixparam = iLookPARAM%critRichNumber         ! critical value for the bulk Richardson number (-)
  case('Louis79_bparam'           ); get_ixparam = iLookPARAM%Louis79_bparam         ! parameter in Louis (1979) stability function (-)
  case('Louis79_cStar'            ); get_ixparam = iLookPARAM%Louis79_cStar          ! parameter in Louis (1979) stability function (-)
  case('Mahrt87_eScale'           ); get_ixparam = iLookPARAM%Mahrt87_eScale         ! exponential scaling factor in the Mahrt (1987) stability function (-)
  case('leafExchangeCoeff'        ); get_ixparam = iLookPARAM%leafExchangeCoeff      ! turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
  case('windReductionParam'       ); get_ixparam = iLookPARAM%windReductionParam     ! canopy wind reduction parameter (-)
  ! vegetation properties
  case('winterSAI'                ); get_ixparam = iLookPARAM%winterSAI              ! stem area index prior to the start of the growing season (m2 m-2)
  case('summerLAI'                ); get_ixparam = iLookPARAM%summerLAI              ! maximum leaf area index at the peak of the growing season (m2 m-2)
  case('rootingDepth'             ); get_ixparam = iLookPARAM%rootingDepth           ! rooting depth (m)
  case('rootDistExp'              ); get_ixparam = iLookPARAM%rootDistExp            ! exponent for the vertical distriution of root density (-)
  case('plantWiltPsi'             ); get_ixparam = iLookPARAM%plantWiltPsi           ! matric head at wilting point (m)
  case('soilStressParam'          ); get_ixparam = iLookPARAM%soilStressParam        ! parameter in the exponential soil stress function
  case('critSoilWilting'          ); get_ixparam = iLookPARAM%critSoilWilting        ! critical vol. liq. water content when plants are wilting (-)
  case('critSoilTranspire'        ); get_ixparam = iLookPARAM%critSoilTranspire      ! critical vol. liq. water content when transpiration is limited (-)
  case('critAquiferTranspire'     ); get_ixparam = iLookPARAM%critAquiferTranspire   ! critical aquifer storage value when transpiration is limited (m)
  case('minStomatalResistance'    ); get_ixparam = iLookPARAM%minStomatalResistance  ! minimum canopy resistance (s m-1)
  case('leafDimension'            ); get_ixparam = iLookPARAM%leafDimension          ! characteristic leaf dimension (m)
  case('heightCanopyTop'          ); get_ixparam = iLookPARAM%heightCanopyTop        ! height of top of the vegetation canopy above ground surface (m)
  case('heightCanopyBottom'       ); get_ixparam = iLookPARAM%heightCanopyBottom     ! height of bottom of the vegetation canopy above ground surface (m)
  case('specificHeatVeg'          ); get_ixparam = iLookPARAM%specificHeatVeg        ! specific heat of vegetation (J kg-1 K-1)
  case('maxMassVegetation'        ); get_ixparam = iLookPARAM%maxMassVegetation      ! maximum mass of vegetation (full foliage) (kg m-2)
  case('throughfallScaleSnow'     ); get_ixparam = iLookPARAM%throughfallScaleSnow   ! scaling factor for throughfall (snow) (-)
  case('throughfallScaleRain'     ); get_ixparam = iLookPARAM%throughfallScaleRain   ! scaling factor for throughfall (rain) (-)
  case('refInterceptCapSnow'      ); get_ixparam = iLookPARAM%refInterceptCapSnow    ! reference canopy interception capacity per unit leaf area (snow) (kg m-2)
  case('refInterceptCapRain'      ); get_ixparam = iLookPARAM%refInterceptCapRain    ! canopy interception capacity per unit leaf area (rain) (kg m-2)
  case('snowUnloadingCoeff'       ); get_ixparam = iLookPARAM%snowUnloadingCoeff     ! time constant for unloading of snow from the forest canopy (s-1)
  case('canopyDrainageCoeff'      ); get_ixparam = iLookPARAM%canopyDrainageCoeff    ! time constant for drainage of liquid water from the forest canopy (s-1)
  case('ratioDrip2Unloading'      ); get_ixparam = iLookPARAM%ratioDrip2Unloading    ! ratio of canopy drip to unloading of snow from the forest canopy (-)
  ! soil properties
  case('soil_dens_intr'           ); get_ixparam = iLookPARAM%soil_dens_intr         ! intrinsic soil density (kg m-3)
  case('thCond_soil'              ); get_ixparam = iLookPARAM%thCond_soil            ! thermal conductivity of soil (W m-1 K-1)
  case('frac_sand'                ); get_ixparam = iLookPARAM%frac_sand              ! fraction of sand (-)
  case('frac_silt'                ); get_ixparam = iLookPARAM%frac_silt              ! fraction of silt (-)
  case('frac_clay'                ); get_ixparam = iLookPARAM%frac_clay              ! fraction of clay (-)
  case('fieldCapacity'            ); get_ixparam = iLookPARAM%fieldCapacity          ! field capacity (-)
  case('wettingFrontSuction'      ); get_ixparam = iLookPARAM%wettingFrontSuction    ! Green-Ampt wetting front suction (m)
  case('theta_mp'                 ); get_ixparam = iLookPARAM%theta_mp               ! volumetric liquid water content when macropore flow begins (-)
  case('theta_sat'                ); get_ixparam = iLookPARAM%theta_sat              ! soil porosity (-)
  case('theta_res'                ); get_ixparam = iLookPARAM%theta_res              ! volumetric residual water content (-)
  case('vGn_alpha'                ); get_ixparam = iLookPARAM%vGn_alpha              ! van Genuchten "alpha" parameter (m-1)
  case('vGn_n'                    ); get_ixparam = iLookPARAM%vGn_n                  ! van Genuchten "n" parameter (-)
  case('mpExp'                    ); get_ixparam = iLookPARAM%mpExp                  ! empirical exponent in macropore flow equation (-)
  case('k_soil'                   ); get_ixparam = iLookPARAM%k_soil                 ! saturated hydraulic conductivity (m s-1)
  case('k_macropore'              ); get_ixparam = iLookPARAM%k_macropore            ! saturated hydraulic conductivity for the macropores (m s-1)
  case('kAnisotropic'             ); get_ixparam = iLookPARAM%kAnisotropic           ! anisotropy factor for lateral hydraulic conductivity (-)
  case('zScale_TOPMODEL'          ); get_ixparam = iLookPARAM%zScale_TOPMODEL        ! TOPMODEL scaling factor used in lower boundary condition for soil (m)
  case('compactedDepth'           ); get_ixparam = iLookPARAM%compactedDepth         ! depth where k_soil reaches the compacted value given by CH78 (m)
  case('aquiferScaleFactor'       ); get_ixparam = iLookPARAM%aquiferScaleFactor     ! scaling factor for aquifer storage in the big bucket (m)
  case('aquiferBaseflowExp'       ); get_ixparam = iLookPARAM%aquiferBaseflowExp     ! baseflow exponent (-)
  case('qSurfScale'               ); get_ixparam = iLookPARAM%qSurfScale             ! scaling factor in the surface runoff parameterization (-)
  case('specificYield'            ); get_ixparam = iLookPARAM%specificYield          ! specific yield (-)
  case('specificStorage'          ); get_ixparam = iLookPARAM%specificStorage        ! specific storage coefficient (m-1)
  case('f_impede'                 ); get_ixparam = iLookPARAM%f_impede               ! ice impedence factor (-)
  case('soilIceScale'             ); get_ixparam = iLookPARAM%soilIceScale           ! scaling factor for depth of soil ice, used to get frozen fraction (m)
  case('soilIceCV'                ); get_ixparam = iLookPARAM%soilIceCV              ! CV of depth of soil ice, used to get frozen fraction (-)
  ! algorithmic control parameters
  case('minwind'                  ); get_ixparam = iLookPARAM%minwind                ! minimum wind speed (m s-1)
  case('minstep'                  ); get_ixparam = iLookPARAM%minstep                ! minimum length of the time step
  case('maxstep'                  ); get_ixparam = iLookPARAM%maxstep                ! maximum length of the time step
  case('wimplicit'                ); get_ixparam = iLookPARAM%wimplicit              ! weight assigned to start-of-step fluxes
  case('maxiter'                  ); get_ixparam = iLookPARAM%maxiter                ! maximum number of iterations
  case('relConvTol_liquid'        ); get_ixparam = iLookPARAM%relConvTol_liquid      ! relative convergence tolerance for vol frac liq water (-)
  case('absConvTol_liquid'        ); get_ixparam = iLookPARAM%absConvTol_liquid      ! absolute convergence tolerance for vol frac liq water (-)
  case('relConvTol_matric'        ); get_ixparam = iLookPARAM%relConvTol_matric      ! relative convergence tolerance for matric head (-)
  case('absConvTol_matric'        ); get_ixparam = iLookPARAM%absConvTol_matric      ! absolute convergence tolerance for matric head (m)
  case('relConvTol_energy'        ); get_ixparam = iLookPARAM%relConvTol_energy      ! relative convergence tolerance for energy (-)
  case('absConvTol_energy'        ); get_ixparam = iLookPARAM%absConvTol_energy      ! absolute convergence tolerance for energy (J m-3)
  case('relConvTol_aquifr'        ); get_ixparam = iLookPARAM%relConvTol_aquifr      ! relative convergence tolerance for aquifer storage (-)
  case('absConvTol_aquifr'        ); get_ixparam = iLookPARAM%absConvTol_aquifr      ! absolute convergence tolerance for aquifer storage (m)
  case('zmin'                     ); get_ixparam = iLookPARAM%zmin                   ! minimum layer depth (m)
  case('zmax'                     ); get_ixparam = iLookPARAM%zmax                   ! maximum layer depth (m)
  case('zminLayer1'               ); get_ixparam = iLookPARAM%zminLayer1             ! minimum layer depth for the 1st (top) layer (m)
  case('zminLayer2'               ); get_ixparam = iLookPARAM%zminLayer2             ! minimum layer depth for the 2nd layer (m)
  case('zminLayer3'               ); get_ixparam = iLookPARAM%zminLayer3             ! minimum layer depth for the 3rd layer (m)
  case('zminLayer4'               ); get_ixparam = iLookPARAM%zminLayer4             ! minimum layer depth for the 4th layer (m)
  case('zminLayer5'               ); get_ixparam = iLookPARAM%zminLayer5             ! minimum layer depth for the 5th (bottom) layer (m)
  case('zmaxLayer1_lower'         ); get_ixparam = iLookPARAM%zmaxLayer1_lower       ! maximum layer depth for the 1st (top) layer when only 1 layer (m)
  case('zmaxLayer2_lower'         ); get_ixparam = iLookPARAM%zmaxLayer2_lower       ! maximum layer depth for the 2nd layer when only 2 layers (m)
  case('zmaxLayer3_lower'         ); get_ixparam = iLookPARAM%zmaxLayer3_lower       ! maximum layer depth for the 3rd layer when only 3 layers (m)
  case('zmaxLayer4_lower'         ); get_ixparam = iLookPARAM%zmaxLayer4_lower       ! maximum layer depth for the 4th layer when only 4 layers (m)
  case('zmaxLayer1_upper'         ); get_ixparam = iLookPARAM%zmaxLayer1_upper       ! maximum layer depth for the 1st (top) layer when > 1 layer (m)
  case('zmaxLayer2_upper'         ); get_ixparam = iLookPARAM%zmaxLayer2_upper       ! maximum layer depth for the 2nd layer when > 2 layers (m)
  case('zmaxLayer3_upper'         ); get_ixparam = iLookPARAM%zmaxLayer3_upper       ! maximum layer depth for the 3rd layer when > 3 layers (m)
  case('zmaxLayer4_upper'         ); get_ixparam = iLookPARAM%zmaxLayer4_upper       ! maximum layer depth for the 4th layer when > 4 layers (m)
  ! get to here if cannot find the variable
  case default
   get_ixparam = imiss
 endselect
 end function get_ixparam


 ! *******************************************************************************************************************
 ! public function get_ixmvar: get the index of the named variables for the model variables
 ! *******************************************************************************************************************
 function get_ixmvar(varName)
 USE var_lookup,only:iLookMVAR                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! parameter name
 integer(i4b)             :: get_ixmvar              ! index of the named variable
 ! define local variables
 integer(i4b), parameter  :: imiss = -999            ! missing value
 ! get the index of the named variables
 select case(trim(varName))
  ! define timestep-average fluxes for a few key variables
  case('totalSoilCompress'              ); get_ixmvar = iLookMVAR%totalSoilCompress                ! change in total soil storage due to compression of the soil matrix (kg m-2)
  case('averageThroughfallSnow'         ); get_ixmvar = iLookMVAR%averageThroughfallSnow           ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
  case('averageThroughfallRain'         ); get_ixmvar = iLookMVAR%averageThroughfallRain           ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
  case('averageCanopySnowUnloading'     ); get_ixmvar = iLookMVAR%averageCanopySnowUnloading       ! unloading of snow from the vegetion canopy (kg m-2 s-1)
  case('averageCanopyLiqDrainage'       ); get_ixmvar = iLookMVAR%averageCanopyLiqDrainage         ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
  case('averageCanopyMeltFreeze'        ); get_ixmvar = iLookMVAR%averageCanopyMeltFreeze          ! melt/freeze of water stored in the canopy (kg m-2 s-1)
  case('averageCanopyTranspiration'     ); get_ixmvar = iLookMVAR%averageCanopyTranspiration       ! canopy transpiration (kg m-2 s-1)
  case('averageCanopyEvaporation'       ); get_ixmvar = iLookMVAR%averageCanopyEvaporation         ! canopy evaporation/condensation (kg m-2 s-1)
  case('averageCanopySublimation'       ); get_ixmvar = iLookMVAR%averageCanopySublimation         ! canopy sublimation/frost (kg m-2 s-1)
  case('averageSnowSublimation'         ); get_ixmvar = iLookMVAR%averageSnowSublimation           ! snow sublimation/frost - below canopy or non-vegetated (kg m-2 s-1)
  case('averageGroundEvaporation'       ); get_ixmvar = iLookMVAR%averageGroundEvaporation         ! ground evaporation/condensation - below canopy or non-vegetated (kg m-2 s-1)
  case('averageRainPlusMelt'            ); get_ixmvar = iLookMVAR%averageRainPlusMelt              ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  case('averageSurfaceRunoff'           ); get_ixmvar = iLookMVAR%averageSurfaceRunoff             ! surface runoff (m s-1)
  case('averageSoilInflux'              ); get_ixmvar = iLookMVAR%averageSoilInflux                ! influx of water at the top of the soil profile (m s-1)
  case('averageSoilBaseflow'            ); get_ixmvar = iLookMVAR%averageSoilBaseflow              ! total baseflow from throughout the soil profile (m s-1)
  case('averageSoilDrainage'            ); get_ixmvar = iLookMVAR%averageSoilDrainage              ! drainage from the bottom of the soil profile (m s-1)
  case('averageAquiferRecharge'         ); get_ixmvar = iLookMVAR%averageAquiferRecharge           ! recharge to the aquifer (m s-1)
  case('averageAquiferBaseflow'         ); get_ixmvar = iLookMVAR%averageAquiferBaseflow           ! baseflow from the aquifer (m s-1)
  case('averageAquiferTranspire'        ); get_ixmvar = iLookMVAR%averageAquiferTranspire          ! transpiration from the aquifer (m s-1)
  case('averageColumnOutflow'           ); get_ixmvar = iLookMVAR%averageColumnOutflow             ! outflow from each layer in the soil profile (m3 s-1)
  ! scalar variables -- forcing
  case('scalarCosZenith'                ); get_ixmvar = iLookMVAR%scalarCosZenith                  ! cosine of the solar zenith angle (0-1)
  case('scalarFractionDirect'           ); get_ixmvar = iLookMVAR%scalarFractionDirect             ! fraction of direct radiation (0-1)
  case('spectralIncomingDirect'         ); get_ixmvar = iLookMVAR%spectralIncomingDirect           ! incoming direct solar radiation in each wave band (W m-2)
  case('spectralIncomingDiffuse'        ); get_ixmvar = iLookMVAR%spectralIncomingDiffuse          ! incoming diffuse solar radiation in each wave band (W m-2)
  case('scalarVPair'                    ); get_ixmvar = iLookMVAR%scalarVPair                      ! vapor pressure of the air above the vegetation canopy (Pa)
  case('scalarTwetbulb'                 ); get_ixmvar = iLookMVAR%scalarTwetbulb                   ! wetbulb temperature (K)
  case('scalarRainfall'                 ); get_ixmvar = iLookMVAR%scalarRainfall                   ! computed rainfall rate (kg m-2 s-1)
  case('scalarSnowfall'                 ); get_ixmvar = iLookMVAR%scalarSnowfall                   ! computed snowfall rate (kg m-2 s-1)
  case('scalarSnowfallTemp'             ); get_ixmvar = iLookMVAR%scalarSnowfallTemp               ! temperature of fresh snow (K)
  case('scalarNewSnowDensity'           ); get_ixmvar = iLookMVAR%scalarNewSnowDensity             ! density of fresh snow, should snow be falling in this time step (kg m-3)
  case('scalarO2air'                    ); get_ixmvar = iLookMVAR%scalarO2air                      ! atmospheric o2 concentration (Pa)
  case('scalarCO2air'                   ); get_ixmvar = iLookMVAR%scalarCO2air                     ! atmospheric co2 concentration (Pa)
  ! scalar variables -- state variables
  case('scalarCanopyIce'                ); get_ixmvar = iLookMVAR%scalarCanopyIce                  ! mass of ice on the vegetation canopy (kg m-2)
  case('scalarCanopyLiq'                ); get_ixmvar = iLookMVAR%scalarCanopyLiq                  ! mass of liquid water on the vegetation canopy (kg m-2)
  case('scalarCanairTemp'               ); get_ixmvar = iLookMVAR%scalarCanairTemp                 ! temperature of the canopy air space (K)
  case('scalarCanopyTemp'               ); get_ixmvar = iLookMVAR%scalarCanopyTemp                 ! temperature of the vegetation canopy (K)
  case('scalarSnowAge'                  ); get_ixmvar = iLookMVAR%scalarSnowAge                    ! non-dimensional snow age (-)
  case('scalarSnowAlbedo'               ); get_ixmvar = iLookMVAR%scalarSnowAlbedo                 ! snow albedo for the entire spectral band (-)
  case('spectralSnowAlbedoDirect'       ); get_ixmvar = iLookMVAR%spectralSnowAlbedoDirect         ! direct snow albedo for individual spectral bands (-)
  case('spectralSnowAlbedoDiffuse'      ); get_ixmvar = iLookMVAR%spectralSnowAlbedoDiffuse        ! diffuse snow albedo for individual spectral bands (-)
  case('scalarSnowDepth'                ); get_ixmvar = iLookMVAR%scalarSnowDepth                  ! total snow depth (m)
  case('scalarSWE'                      ); get_ixmvar = iLookMVAR%scalarSWE                        ! snow water equivalent (kg m-2)
  case('scalarSfcMeltPond'              ); get_ixmvar = iLookMVAR%scalarSfcMeltPond                ! ponded water caused by melt of the "snow without a layer" (kg m-2)
  case('scalarAquiferStorage'           ); get_ixmvar = iLookMVAR%scalarAquiferStorage             ! relative aquifer storage -- above bottom of the soil profile (m)
  case('scalarSurfaceTemp'              ); get_ixmvar = iLookMVAR%scalarSurfaceTemp                ! surface temperature (K)
  ! NOAH-MP vegetation variables (general)
  case('scalarGreenVegFraction'         ); get_ixmvar = iLookMVAR%scalarGreenVegFraction           ! green vegetation fraction used to compute LAI (-)
  case('scalarBulkVolHeatCapVeg'        ); get_ixmvar = iLookMVAR%scalarBulkVolHeatCapVeg          ! bulk volumetric heat capacity of vegetation (J m-3 K-1)
  case('scalarRootZoneTemp'             ); get_ixmvar = iLookMVAR%scalarRootZoneTemp               ! average temperature of the root zone (K)
  case('scalarLAI'                      ); get_ixmvar = iLookMVAR%scalarLAI                        ! one-sided leaf area index (m2 m-2)
  case('scalarSAI'                      ); get_ixmvar = iLookMVAR%scalarSAI                        ! one-sided stem area index (m2 m-2)
  case('scalarExposedLAI'               ); get_ixmvar = iLookMVAR%scalarExposedLAI                 ! exposed leaf area index after burial by snow (m2 m-2)
  case('scalarExposedSAI'               ); get_ixmvar = iLookMVAR%scalarExposedSAI                 ! exposed stem area index after burial by snow (m2 m-2)
  case('scalarCanopyIceMax'             ); get_ixmvar = iLookMVAR%scalarCanopyIceMax               ! maximum interception storage capacity for ice (kg m-2)
  case('scalarCanopyLiqMax'             ); get_ixmvar = iLookMVAR%scalarCanopyLiqMax               ! maximum interception storage capacity for liquid water (kg m-2)
  case('scalarGrowingSeasonIndex'       ); get_ixmvar = iLookMVAR%scalarGrowingSeasonIndex         ! growing season index (0=off, 1=on)
  case('scalarVP_CanopyAir'             ); get_ixmvar = iLookMVAR%scalarVP_CanopyAir               ! vapor pressure of the canopy air space (Pa)
  ! NOAH-MP vegetation variables (shortwave radiation)
  case('scalarCanopySunlitFraction'     ); get_ixmvar = iLookMVAR%scalarCanopySunlitFraction       ! sunlit fraction of canopy (-)
  case('scalarCanopySunlitLAI'          ); get_ixmvar = iLookMVAR%scalarCanopySunlitLAI            ! sunlit leaf area (-)
  case('scalarCanopyShadedLAI'          ); get_ixmvar = iLookMVAR%scalarCanopyShadedLAI            ! shaded leaf area (-)
  case('scalarCanopySunlitPAR'          ); get_ixmvar = iLookMVAR%scalarCanopySunlitPAR            ! average absorbed par for sunlit leaves (w m-2)
  case('scalarCanopyShadedPAR'          ); get_ixmvar = iLookMVAR%scalarCanopyShadedPAR            ! average absorbed par for shaded leaves (w m-2)
  case('spectralBelowCanopyDirect'      ); get_ixmvar = iLookMVAR%spectralBelowCanopyDirect        ! downward direct flux below veg layer for each spectral band  W m-2)
  case('spectralBelowCanopyDiffuse'     ); get_ixmvar = iLookMVAR%spectralBelowCanopyDiffuse       ! downward diffuse flux below veg layer for each spectral band (W m-2)
  case('scalarBelowCanopySolar'         ); get_ixmvar = iLookMVAR%scalarBelowCanopySolar           ! solar radiation transmitted below the canopy (W m-2)
  case('spectralAlbGndDirect'           ); get_ixmvar = iLookMVAR%spectralAlbGndDirect             ! direct  albedo of underlying surface for each spectral band (-)
  case('spectralAlbGndDiffuse'          ); get_ixmvar = iLookMVAR%spectralAlbGndDiffuse            ! diffuse albedo of underlying surface for each spectral band (-)
  case('scalarGroundAlbedo'             ); get_ixmvar = iLookMVAR%scalarGroundAlbedo               ! albedo of the ground surface (-)
  case('scalarCanopyAbsorbedSolar'      ); get_ixmvar = iLookMVAR%scalarCanopyAbsorbedSolar        ! solar radiation absorbed by canopy (W m-2)
  case('scalarGroundAbsorbedSolar'      ); get_ixmvar = iLookMVAR%scalarGroundAbsorbedSolar        ! solar radiation absorbed by ground (W m-2)
  ! NOAH-MP vegetation variables (longwave radiation)
  case('scalarCanopyEmissivity'         ); get_ixmvar = iLookMVAR%scalarCanopyEmissivity           ! effective canopy emissivity (-)
  case('scalarLWRadCanopy'              ); get_ixmvar = iLookMVAR%scalarLWRadCanopy                ! longwave radiation emitted from the canopy (W m-2)
  case('scalarLWRadGround'              ); get_ixmvar = iLookMVAR%scalarLWRadGround                ! longwave radiation emitted at the ground surface  (W m-2)
  case('scalarLWRadUbound2Canopy'       ); get_ixmvar = iLookMVAR%scalarLWRadUbound2Canopy         ! downward atmospheric longwave radiation absorbed by the canopy (W m-2)
  case('scalarLWRadUbound2Ground'       ); get_ixmvar = iLookMVAR%scalarLWRadUbound2Ground         ! downward atmospheric longwave radiation absorbed by the ground (W m-2)
  case('scalarLWRadUbound2Ubound'       ); get_ixmvar = iLookMVAR%scalarLWRadUbound2Ubound         ! atmospheric radiation refl by ground + lost thru upper boundary (W m-2)
  case('scalarLWRadCanopy2Ubound'       ); get_ixmvar = iLookMVAR%scalarLWRadCanopy2Ubound         ! longwave radiation emitted from canopy lost thru upper boundary (W m-2)
  case('scalarLWRadCanopy2Ground'       ); get_ixmvar = iLookMVAR%scalarLWRadCanopy2Ground         ! longwave radiation emitted from canopy absorbed by the ground (W m-2)
  case('scalarLWRadCanopy2Canopy'       ); get_ixmvar = iLookMVAR%scalarLWRadCanopy2Canopy         ! canopy longwave reflected from ground and absorbed by the canopy (W m-2)
  case('scalarLWRadGround2Ubound'       ); get_ixmvar = iLookMVAR%scalarLWRadGround2Ubound         ! longwave radiation emitted from ground lost thru upper boundary (W m-2)
  case('scalarLWRadGround2Canopy'       ); get_ixmvar = iLookMVAR%scalarLWRadGround2Canopy         ! longwave radiation emitted from ground and absorbed by the canopy (W m-2)
  case('scalarLWNetCanopy'              ); get_ixmvar = iLookMVAR%scalarLWNetCanopy                ! net longwave radiation at the canopy (W m-2)
  case('scalarLWNetGround'              ); get_ixmvar = iLookMVAR%scalarLWNetGround                ! net longwave radiation at the ground surface (W m-2)
  case('scalarLWNetUbound'              ); get_ixmvar = iLookMVAR%scalarLWNetUbound                ! net longwave radiation at the upper atmospheric boundary (W m-2)
  ! NOAH-MP vegetation variables (turbulent heat transfer)
  case('scalarLatHeatSubVapCanopy'      ); get_ixmvar = iLookMVAR%scalarLatHeatSubVapCanopy        ! latent heat of sublimation/vaporization used for veg canopy (J kg-1)
  case('scalarLatHeatSubVapGround'      ); get_ixmvar = iLookMVAR%scalarLatHeatSubVapGround        ! latent heat of sublimation/vaporization used for ground surface (J kg-1)
  case('scalarSatVP_CanopyTemp'         ); get_ixmvar = iLookMVAR%scalarSatVP_CanopyTemp           ! saturation vapor pressure at the temperature of vegetation canopy (Pa)
  case('scalarSatVP_GroundTemp'         ); get_ixmvar = iLookMVAR%scalarSatVP_GroundTemp           ! saturation vapor pressure at the temperature of the ground (Pa)
  case('scalarZ0Canopy'                 ); get_ixmvar = iLookMVAR%scalarZ0Canopy                   ! roughness length of the canopy (m)
  case('scalarWindReductionFactor'      ); get_ixmvar = iLookMVAR%scalarWindReductionFactor        ! canopy wind reduction factor (-)
  case('scalarZeroPlaneDisplacement'    ); get_ixmvar = iLookMVAR%scalarZeroPlaneDisplacement      ! zero plane displacement (m)
  case('scalarRiBulkCanopy'             ); get_ixmvar = iLookMVAR%scalarRiBulkCanopy               ! bulk Richardson number for the canopy (-)
  case('scalarRiBulkGround'             ); get_ixmvar = iLookMVAR%scalarRiBulkGround               ! bulk Richardson number for the ground surface (-)
  case('scalarCanopyStabilityCorrection'); get_ixmvar = iLookMVAR%scalarCanopyStabilityCorrection  ! stability correction for the canopy (-)
  case('scalarGroundStabilityCorrection'); get_ixmvar = iLookMVAR%scalarGroundStabilityCorrection  ! stability correction for the ground surface (-)
  case('scalarEddyDiffusCanopyTop'      ); get_ixmvar = iLookMVAR%scalarEddyDiffusCanopyTop        ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
  case('scalarFrictionVelocity'         ); get_ixmvar = iLookMVAR%scalarFrictionVelocity           ! friction velocity - canopy momentum sink (m s-1)
  case('scalarWindspdCanopyTop'         ); get_ixmvar = iLookMVAR%scalarWindspdCanopyTop           ! windspeed at the top of the canopy (m s-1)
  case('scalarWindspdCanopyBottom'      ); get_ixmvar = iLookMVAR%scalarWindspdCanopyBottom        ! windspeed at the height of the bottom of the canopy (m s-1)
  case('scalarGroundResistance'         ); get_ixmvar = iLookMVAR%scalarGroundResistance           ! below canopy aerodynamic resistance (s m-1)
  case('scalarCanopyResistance'         ); get_ixmvar = iLookMVAR%scalarCanopyResistance           ! above canopy aerodynamic resistance (s m-1)
  case('scalarLeafResistance'           ); get_ixmvar = iLookMVAR%scalarLeafResistance             ! mean leaf boundary layer resistance per unit leaf area (s m-1)
  case('scalarSoilResistance'           ); get_ixmvar = iLookMVAR%scalarSoilResistance             ! soil surface resistance (s m-1)
  case('scalarSoilRelHumidity'          ); get_ixmvar = iLookMVAR%scalarSoilRelHumidity            ! relative humidity in the soil pores in the upper-most soil layer (-)
  case('scalarSenHeatTotal'             ); get_ixmvar = iLookMVAR%scalarSenHeatTotal               ! sensible heat from the canopy air space to the atmosphere (W m-2)
  case('scalarSenHeatCanopy'            ); get_ixmvar = iLookMVAR%scalarSenHeatCanopy              ! sensible heat from the canopy to the canopy air space (W m-2)
  case('scalarSenHeatGround'            ); get_ixmvar = iLookMVAR%scalarSenHeatGround              ! sensible heat from the ground (below canopy or non-vegetated) (W m-2)
  case('scalarLatHeatTotal'             ); get_ixmvar = iLookMVAR%scalarLatHeatTotal               ! latent heat from the canopy air space to the atmosphere (W m-2)
  case('scalarLatHeatCanopyEvap'        ); get_ixmvar = iLookMVAR%scalarLatHeatCanopyEvap          ! evaporation latent heat from the canopy to the canopy air space (W m-2)
  case('scalarLatHeatCanopyTrans'       ); get_ixmvar = iLookMVAR%scalarLatHeatCanopyTrans         ! transpiration latent heat from the canopy to the canopy air space (W m-2)
  case('scalarLatHeatGround'            ); get_ixmvar = iLookMVAR%scalarLatHeatGround              ! latent heat from the ground (below canopy or non-vegetated) (W m-2)
  case('scalarCanopyAdvectiveHeatFlux'  ); get_ixmvar = iLookMVAR%scalarCanopyAdvectiveHeatFlux    ! heat advected to the canopy surface with rain + snow (W m-2)
  case('scalarGroundAdvectiveHeatFlux'  ); get_ixmvar = iLookMVAR%scalarGroundAdvectiveHeatFlux    ! heat advected to the ground surface with throughfall and unloading/drainage (W m-2)
  case('scalarCanopyTranspiration'      ); get_ixmvar = iLookMVAR%scalarCanopyTranspiration        ! canopy transpiration (kg m-2 s-1)
  case('scalarCanopyEvaporation'        ); get_ixmvar = iLookMVAR%scalarCanopyEvaporation          ! canopy evaporation/condensation (kg m-2 s-1)
  case('scalarCanopySublimation'        ); get_ixmvar = iLookMVAR%scalarCanopySublimation          ! canopy sublimation/frost (kg m-2 s-1)
  case('scalarGroundEvaporation'        ); get_ixmvar = iLookMVAR%scalarGroundEvaporation          ! ground evaporation/condensation (below canopy or non-vegetated) (kg m-2 s-1)
  case('scalarSnowSublimation'          ); get_ixmvar = iLookMVAR%scalarSnowSublimation            ! snow sublimation/frost (below canopy or non-vegetated) (kg m-2 s-1)
  ! NOAH-MP vegetation variables (transpiration)
  case('scalarTranspireLim'             ); get_ixmvar = iLookMVAR%scalarTranspireLim               ! aggregate soil moisture and aquifer storage limit on transpiration (-)
  case('scalarTranspireLimAqfr'         ); get_ixmvar = iLookMVAR%scalarTranspireLimAqfr           ! aquifer storage limit on transpiration (-)
  case('scalarFoliageNitrogenFactor'    ); get_ixmvar = iLookMVAR%scalarFoliageNitrogenFactor      ! foliage nitrogen concentration, 1=saturated (-)
  case('scalarStomResistSunlit'         ); get_ixmvar = iLookMVAR%scalarStomResistSunlit           ! stomatal resistance for sunlit leaves (s m-1)
  case('scalarStomResistShaded'         ); get_ixmvar = iLookMVAR%scalarStomResistShaded           ! stomatal resistance for shaded leaves (s m-1)
  case('scalarPhotosynthesisSunlit'     ); get_ixmvar = iLookMVAR%scalarPhotosynthesisSunlit       ! sunlit photosynthesis (umolco2 m-2 s-1)
  case('scalarPhotosynthesisShaded'     ); get_ixmvar = iLookMVAR%scalarPhotosynthesisShaded       ! shaded photosynthesis (umolco2 m-2 s-1)
  ! NOAH-MP vegetation variables (canopy water)
  case('scalarCanopyWetFraction'        ); get_ixmvar = iLookMVAR%scalarCanopyWetFraction          ! fraction of canopy that is wet
  case('scalarGroundSnowFraction'       ); get_ixmvar = iLookMVAR%scalarGroundSnowFraction         ! fraction of ground that is covered with snow (-)
  case('scalarThroughfallSnow'          ); get_ixmvar = iLookMVAR%scalarThroughfallSnow            ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
  case('scalarThroughfallRain'          ); get_ixmvar = iLookMVAR%scalarThroughfallRain            ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
  case('scalarCanopySnowUnloading'      ); get_ixmvar = iLookMVAR%scalarCanopySnowUnloading        ! unloading of snow from the vegetion canopy (kg m-2 s-1)
  case('scalarCanopyLiqDrainage'        ); get_ixmvar = iLookMVAR%scalarCanopyLiqDrainage          ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
  case('scalarCanopyMeltFreeze'         ); get_ixmvar = iLookMVAR%scalarCanopyMeltFreeze           ! melt/freeze of water stored in the canopy (kg m-2 s-1)
  ! scalar variables -- soil and aquifer fluxes
  case('scalarRainPlusMelt'             ); get_ixmvar = iLookMVAR%scalarRainPlusMelt               ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  case('scalarInfilArea'                ); get_ixmvar = iLookMVAR%scalarInfilArea                  ! fraction of unfrozen area where water can infiltrate (-)
  case('scalarFrozenArea'               ); get_ixmvar = iLookMVAR%scalarFrozenArea                 ! fraction of area that is considered impermeable due to soil ice (-)
  case('scalarInfiltration'             ); get_ixmvar = iLookMVAR%scalarInfiltration               ! infiltration of water into the soil profile (m s-1)
  case('scalarExfiltration'             ); get_ixmvar = iLookMVAR%scalarExfiltration               ! exfiltration of water from the top of the soil profile (m s-1)
  case('scalarSurfaceRunoff'            ); get_ixmvar = iLookMVAR%scalarSurfaceRunoff              ! surface runoff (m s-1)
  case('scalarInitAquiferRecharge'      ); get_ixmvar = iLookMVAR%scalarInitAquiferRecharge        ! recharge to the aquifer -- at the start of the step (m s-1)
  case('scalarAquiferRecharge'          ); get_ixmvar = iLookMVAR%scalarAquiferRecharge            ! recharge to the aquifer (m s-1)
  case('scalarInitAquiferTranspire'     ); get_ixmvar = iLookMVAR%scalarInitAquiferTranspire       ! transpiration from the aquifer (m s-1)
  case('scalarAquiferTranspire'         ); get_ixmvar = iLookMVAR%scalarAquiferTranspire           ! transpiration from the aquifer (m s-1)
  case('scalarInitAquiferBaseflow'      ); get_ixmvar = iLookMVAR%scalarInitAquiferBaseflow        ! baseflow from the aquifer (m s-1)
  case('scalarAquiferBaseflow'          ); get_ixmvar = iLookMVAR%scalarAquiferBaseflow            ! baseflow from the aquifer (m s-1)
  ! scalar variables -- sub-step average fluxes for the soil zone
  case('scalarSoilInflux'               ); get_ixmvar = iLookMVAR%scalarSoilInflux                 ! sub-step average: influx of water at the top of the soil profile (m s-1)
  case('scalarSoilCompress'             ); get_ixmvar = iLookMVAR%scalarSoilCompress               ! change in total soil storage due to compression of the soil matrix (kg m-2)
  case('scalarSoilBaseflow'             ); get_ixmvar = iLookMVAR%scalarSoilBaseflow               ! sub-step average: total baseflow from throughout the soil profile (m s-1)
  case('scalarSoilDrainage'             ); get_ixmvar = iLookMVAR%scalarSoilDrainage               ! sub-step average: drainage from the bottom of the soil profile (m s-1)
  case('scalarSoilTranspiration'        ); get_ixmvar = iLookMVAR%scalarSoilTranspiration          ! sub-step average: total transpiration from the soil (m s-1)
  ! scalar variables -- mass balance check
  case('scalarSoilWatBalError'          ); get_ixmvar = iLookMVAR%scalarSoilWatBalError            ! error in the total soil water balance (kg m-2)
  case('scalarAquiferBalError'          ); get_ixmvar = iLookMVAR%scalarAquiferBalError            ! error in the aquifer water balance (kg m-2)
  case('scalarTotalSoilLiq'             ); get_ixmvar = iLookMVAR%scalarTotalSoilLiq               ! total mass of liquid water in the soil (kg m-2)
  case('scalarTotalSoilIce'             ); get_ixmvar = iLookMVAR%scalarTotalSoilIce               ! total mass of ice in the soil (kg m-2)
  ! variables at the mid-point of each layer -- domain geometry
  case('mLayerDepth'                    ); get_ixmvar = iLookMVAR%mLayerDepth                      ! depth of each layer (m)
  case('mLayerHeight'                   ); get_ixmvar = iLookMVAR%mLayerHeight                     ! height at the midpoint of each layer (m)
  case('mLayerRootDensity'              ); get_ixmvar = iLookMVAR%mLayerRootDensity                ! fraction of roots in each soil layer (-)
  ! variables at the mid-point of each layer -- coupled energy and mass
  case('mLayerTemp'                     ); get_ixmvar = iLookMVAR%mLayerTemp                       ! temperature of each layer (K)
  case('mLayerVolFracAir'               ); get_ixmvar = iLookMVAR%mLayerVolFracAir                 ! volumetric fraction of air in each layer (-)
  case('mLayerVolFracIce'               ); get_ixmvar = iLookMVAR%mLayerVolFracIce                 ! volumetric fraction of icein each layer (-)
  case('mLayerVolFracLiq'               ); get_ixmvar = iLookMVAR%mLayerVolFracLiq                 ! volumetric fraction of liquid water in each layer (-)
  case('mLayerVolHtCapBulk'             ); get_ixmvar = iLookMVAR%mLayerVolHtCapBulk               ! volumetric heat capacity in each layer (J m-3 K-1)
  case('mLayerTcrit'                    ); get_ixmvar = iLookMVAR%mLayerTcrit                      ! critical soil temperature above which all water is unfrozen (K)
  case('mLayerdTheta_dTk'               ); get_ixmvar = iLookMVAR%mLayerdTheta_dTk                 ! analytical derivative in the freezing curve (K-1)
  case('mLayerThermalC'                 ); get_ixmvar = iLookMVAR%mLayerThermalC                   ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
  case('mLayerRadCondFlux'              ); get_ixmvar = iLookMVAR%mLayerRadCondFlux                ! temporal derivative in energy from radiative and conductive flux (J m-2 s-1)
  case('mLayerMeltFreeze'               ); get_ixmvar = iLookMVAR%mLayerMeltFreeze                 ! rate of ice content change from melt/freeze in each layer (kg m-3 s-1)
  case('mLayerInfilFreeze'              ); get_ixmvar = iLookMVAR%mLayerInfilFreeze                ! rate of ice content change by freezing infiltrating flux (kg m-3 s-1)
  case('mLayerSatHydCond'               ); get_ixmvar = iLookMVAR%mLayerSatHydCond                 ! saturated hydraulic conductivity in each layer (m s-1)
  case('mLayerSatHydCondMP'             ); get_ixmvar = iLookMVAR%mLayerSatHydCondMP               ! saturated hydraulic conductivity of macropores in each layer (m s-1)
  case('mLayerMatricHead'               ); get_ixmvar = iLookMVAR%mLayerMatricHead                 ! matric head of water in the soil (m)
  case('mLayerdTheta_dPsi'              ); get_ixmvar = iLookMVAR%mLayerdTheta_dPsi                ! analytical derivative in the soil water characteristic w.r.t. psi (m-1)
  case('mLayerdPsi_dTheta'              ); get_ixmvar = iLookMVAR%mLayerdPsi_dTheta                ! analytical derivative in the soil water characteristic w.r.t. theta (m)
  case('mLayerThetaResid'               ); get_ixmvar = iLookMVAR%mLayerThetaResid                 ! residual volumetric water content in each snow layer (-)
  case('mLayerPoreSpace'                ); get_ixmvar = iLookMVAR%mLayerPoreSpace                  ! total pore space in each snow layer (-)
  case('mLayerCompress'                 ); get_ixmvar = iLookMVAR%mLayerCompress                   ! change in volumetric water content due to compression of soil (-)
  case('mLayerTranspireLim'             ); get_ixmvar = iLookMVAR%mLayerTranspireLim               ! moisture avail factor limiting transpiration in each layer (-)
  case('mLayerInitTranspire'            ); get_ixmvar = iLookMVAR%mLayerInitTranspire              ! transpiration loss from each soil layer at the start of the step (kg m-2 s-1)
  case('mLayerTranspire'                ); get_ixmvar = iLookMVAR%mLayerTranspire                  ! transpiration loss from each soil layer (kg m-2 s-1)
  case('mLayerInitQMacropore'           ); get_ixmvar = iLookMVAR%mLayerInitQMacropore             ! liquid flux from micropores to macropores at the start-of-step (m s-1)
  case('mLayerQMacropore'               ); get_ixmvar = iLookMVAR%mLayerQMacropore                 ! liquid flux from micropores to macropores (m s-1)
  case('mLayerInitBaseflow'             ); get_ixmvar = iLookMVAR%mLayerInitBaseflow               ! baseflow from each soil layer at the start of the time step (m s-1)
  case('mLayerBaseflow'                 ); get_ixmvar = iLookMVAR%mLayerBaseflow                   ! baseflow from each soil layer (m s-1)
  case('mLayerColumnInflow'             ); get_ixmvar = iLookMVAR%mLayerColumnInflow               ! total inflow to each layer in a given soil column (m3 s-1)
  case('mLayerColumnOutflow'            ); get_ixmvar = iLookMVAR%mLayerColumnOutflow              ! total outflow from each layer in a given soil column (m3 s-1)
  ! variables at the interface of each layer
  case('iLayerHeight'                   ); get_ixmvar = iLookMVAR%iLayerHeight                     ! height at the interface of each layer (m)
  case('iLayerThermalC'                 ); get_ixmvar = iLookMVAR%iLayerThermalC                   ! thermal conductivity at the interface of each layer (W m-1 K-1)
  case('iLayerConductiveFlux'           ); get_ixmvar = iLookMVAR%iLayerConductiveFlux             ! conductive energy flux at layer interfaces at end of time step (W m-2)
  case('iLayerAdvectiveFlux'            ); get_ixmvar = iLookMVAR%iLayerAdvectiveFlux              ! advective energy flux at layer interfaces at end of time step (W m-2)
  case('iLayerInitNrgFlux'              ); get_ixmvar = iLookMVAR%iLayerInitNrgFlux                ! energy flux at layer interfaces at the start of the time step (W m-2)
  case('iLayerNrgFlux'                  ); get_ixmvar = iLookMVAR%iLayerNrgFlux                    ! energy flux at layer interfaces at the end of the time step (W m-2)
  case('iLayerSatHydCond'               ); get_ixmvar = iLookMVAR%iLayerSatHydCond                 ! saturated hydraulic conductivity in each layer (m s-1)
  case('iLayerInitLiqFluxSnow'          ); get_ixmvar = iLookMVAR%iLayerInitLiqFluxSnow            ! liquid flux at snow layer interfaces at the start of the time step (m s-1)
  case('iLayerInitLiqFluxSoil'          ); get_ixmvar = iLookMVAR%iLayerInitLiqFluxSoil            ! liquid flux at soil layer interfaces at the start of the time step (m s-1)
  case('iLayerInitFluxReversal'         ); get_ixmvar = iLookMVAR%iLayerInitFluxReversal           ! start of step liquid flux at soil layer interfaces from impedance (m s-1)
  case('iLayerLiqFluxSnow'              ); get_ixmvar = iLookMVAR%iLayerLiqFluxSnow                ! liquid flux at snow layer interfaces at the end of the time step (m s-1)
  case('iLayerLiqFluxSoil'              ); get_ixmvar = iLookMVAR%iLayerLiqFluxSoil                ! liquid flux at soil layer interfaces at the end of the time step (m s-1)
  case('iLayerFluxReversal'             ); get_ixmvar = iLookMVAR%iLayerFluxReversal               ! end of step liquid flux at soil layer interfaces from impedance (m s-1)
  ! time stepping variables
  case('dt_init'                        ); get_ixmvar = iLookMVAR%dt_init                          ! length of initial time step at start of next data interval (s)
  ! "short-cut" variables
  case('scalarVGn_m'                    ); get_ixmvar = iLookMVAR%scalarVGn_m                      ! van Genuchten "m" parameter (-)
  case('scalarKappa'                    ); get_ixmvar = iLookMVAR%scalarKappa                      ! constant in the freezing curve function (m K-1)
  case('scalarVolHtCap_air'             ); get_ixmvar = iLookMVAR%scalarVolHtCap_air               ! volumetric heat capacity air         (J m-3 K-1)
  case('scalarVolHtCap_ice'             ); get_ixmvar = iLookMVAR%scalarVolHtCap_ice               ! volumetric heat capacity ice         (J m-3 K-1)
  case('scalarVolHtCap_soil'            ); get_ixmvar = iLookMVAR%scalarVolHtCap_soil              ! volumetric heat capacity dry soil    (J m-3 K-1)
  case('scalarVolHtCap_water'           ); get_ixmvar = iLookMVAR%scalarVolHtCap_water             ! volumetric heat capacity liquid wat  (J m-3 K-1)
  case('scalarLambda_drysoil'           ); get_ixmvar = iLookMVAR%scalarLambda_drysoil             ! thermal conductivity of dry soil     (W m-1)
  case('scalarLambda_wetsoil'           ); get_ixmvar = iLookMVAR%scalarLambda_wetsoil             ! thermal conductivity of wet soil     (W m-1)
  case('scalarVolLatHt_fus'             ); get_ixmvar = iLookMVAR%scalarVolLatHt_fus               ! volumetric latent heat of fusion     (J m-3)
  case('scalarAquiferRootFrac'          ); get_ixmvar = iLookMVAR%scalarAquiferRootFrac            ! fraction of roots below the soil profile (-)
  ! get to here if cannot find the variable
  case default
   get_ixmvar = imiss
 endselect
 end function get_ixmvar


 ! *******************************************************************************************************************
 ! public function get_ixindex: get the index of the named variables for the model indices
 ! *******************************************************************************************************************
 function get_ixindex(varName)
 USE var_lookup,only:iLookINDEX                      ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! parameter name
 integer(i4b)             :: get_ixINDEX             ! index of the named variable
 ! define local variables
 integer(i4b), parameter  :: imiss = -999            ! missing value
 ! get the index of the named variables
 select case(trim(varName))
  case('nSnow'            ); get_ixindex = iLookINDEX%nSnow             ! number of snow layers
  case('nSoil'            ); get_ixindex = iLookINDEX%nSoil             ! number of soil layers
  case('nLayers'          ); get_ixindex = iLookINDEX%nLayers           ! total number of layers
  case('midSnowStartIndex'); get_ixindex = iLookINDEX%midSnowStartIndex ! start index of the midSnow vector for a given timestep
  case('midSoilStartIndex'); get_ixindex = iLookINDEX%midSoilStartIndex ! start index of the midSoil vector for a given timestep
  case('midTotoStartIndex'); get_ixindex = iLookINDEX%midTotoStartIndex ! start index of the midToto vector for a given timestep
  case('ifcSnowStartIndex'); get_ixindex = iLookINDEX%ifcSnowStartIndex ! start index of the ifcSnow vector for a given timestep
  case('ifcSoilStartIndex'); get_ixindex = iLookINDEX%ifcSoilStartIndex ! start index of the ifcSoil vector for a given timestep
  case('ifcTotoStartIndex'); get_ixindex = iLookINDEX%ifcTotoStartIndex ! start index of the ifcToto vector for a given timestep
  case('layerType'        ); get_ixindex = iLookINDEX%layerType         ! type of layer (soil or snow)
  ! get to here if cannot find the variable
  case default
   get_ixindex = imiss
 endselect
 end function get_ixindex


 ! *******************************************************************************************************************
 ! public function get_ixbpar: get the index of the named variables for the basin-average variables
 ! *******************************************************************************************************************
 function get_ixbpar(varName)
 USE var_lookup,only:iLookBPAR                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! parameter name
 integer(i4b)             :: get_ixbpar              ! index of the named variable
 ! define local variables
 integer(i4b), parameter  :: imiss = -999            ! missing value
 ! get the index of the named variables
 select case(trim(varName))
  ! baseflow
  case('basin__aquiferHydCond'    ); get_ixbpar = iLookBPAR%basin__aquiferHydCond     ! hydraulic conductivity of the basin aquifer (m s-1)
  case('basin__aquiferScaleFactor'); get_ixbpar = iLookBPAR%basin__aquiferScaleFactor ! scaling factor for aquifer storage in the big bucket (m)
  case('basin__aquiferBaseflowExp'); get_ixbpar = iLookBPAR%basin__aquiferBaseflowExp ! baseflow exponent for the big bucket (-)
  ! sub-grid routing
  case('routingGammaShape'        ); get_ixbpar = iLookBPAR%routingGammaShape         ! shape parameter in Gamma distribution used for sub-grid routing (-)
  case('routingGammaScale'        ); get_ixbpar = iLookBPAR%routingGammaScale         ! scale parameter in Gamma distribution used for sub-grid routing (s)
  ! get to here if cannot find the variable
  case default
   get_ixbpar = imiss
 endselect
 end function get_ixbpar


 ! *******************************************************************************************************************
 ! public function get_ixbvar: get the index of the named variables for the basin-average variables
 ! *******************************************************************************************************************
 function get_ixbvar(varName)
 USE var_lookup,only:iLookBVAR                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! parameter name
 integer(i4b)             :: get_ixbvar              ! index of the named variable
 ! define local variables
 integer(i4b), parameter  :: imiss = -999            ! missing value
 ! get the index of the named variables
 select case(trim(varName))
  ! derived variables
  case('basin__totalArea'              ); get_ixbvar = iLookBVAR%basin__totalArea                ! total basin area (m2)
  ! scalar variables -- basin-average runoff and aquifer fluxes
  case('basin__SurfaceRunoff'          ); get_ixbvar = iLookBVAR%basin__SurfaceRunoff            ! surface runoff (m s-1)
  case('basin__ColumnOutflow'          ); get_ixbvar = iLookBVAR%basin__ColumnOutflow            ! outflow from all "outlet" HRUs (those with no downstream HRU)
  case('basin__AquiferStorage'         ); get_ixbvar = iLookBVAR%basin__AquiferStorage           ! aquifer storage (m s-1)
  case('basin__AquiferRecharge'        ); get_ixbvar = iLookBVAR%basin__AquiferRecharge          ! recharge to the aquifer (m s-1)
  case('basin__AquiferBaseflow'        ); get_ixbvar = iLookBVAR%basin__AquiferBaseflow          ! baseflow from the aquifer (m s-1)
  case('basin__AquiferTranspire'       ); get_ixbvar = iLookBVAR%basin__AquiferTranspire         ! transpiration from the aquifer (m s-1)
  ! variables to compute runoff
  case('routingRunoffFuture'           ); get_ixbvar = iLookBVAR%routingRunoffFuture             ! runoff in future time steps (m s-1)
  case('routingFractionFuture'         ); get_ixbvar = iLookBVAR%routingFractionFuture           ! fraction of runoff in future time steps (-)
  case('averageInstantRunoff'          ); get_ixbvar = iLookBVAR%averageInstantRunoff            ! instantaneous runoff (m s-1)
  case('averageRoutedRunoff'           ); get_ixbvar = iLookBVAR%averageRoutedRunoff             ! routed runoff (m s-1)
  ! get to here if cannot find the variable
  case default
   get_ixbvar = imiss
 endselect
 end function get_ixbvar


end module get_ixname_module

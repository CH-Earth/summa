! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
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
USE nrtype, integerMissing=>nr_integerMissing
implicit none
private
public::get_ixdecisions
public::get_ixTime
public::get_ixAttr
public::get_ixType
public::get_ixId
public::get_ixForce
public::get_ixParam
public::get_ixProg
public::get_ixDiag
public::get_ixFlux
public::get_ixDeriv
public::get_ixIndex
public::get_ixBpar
public::get_ixBvar
public::get_ixVarType
public::get_varTypeName
public::get_ixUnknown
public::get_ixFreq
public::get_ixStat
public::get_freqName
public::get_statName
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
 ! get the index of the named variables
 select case(trim(varName))
  case('soilCatTbl'      ); get_ixdecisions=iLookDECISIONS%soilCatTbl  ! soil-category dateset
  case('vegeParTbl'      ); get_ixdecisions=iLookDECISIONS%vegeParTbl  ! vegetation category dataset
  case('soilStress'      ); get_ixdecisions=iLookDECISIONS%soilStress  ! choice of function for the soil moisture control on stomatal resistance
  case('stomResist'      ); get_ixdecisions=iLookDECISIONS%stomResist  ! choice of function for stomatal resistance
  case('bbTempFunc'      ); get_ixdecisions=iLookDECISIONS%bbTempFunc  ! Ball-Berry: leaf temperature controls on photosynthesis + stomatal resistance
  case('bbHumdFunc'      ); get_ixdecisions=iLookDECISIONS%bbHumdFunc  ! Ball-Berry: humidity controls on stomatal resistance
  case('bbElecFunc'      ); get_ixdecisions=iLookDECISIONS%bbElecFunc  ! Ball-Berry: dependence of photosynthesis on PAR
  case('bbCO2point'      ); get_ixdecisions=iLookDECISIONS%bbCO2point  ! Ball-Berry: use of CO2 compensation point to calculate stomatal resistance
  case('bbNumerics'      ); get_ixdecisions=iLookDECISIONS%bbNumerics  ! Ball-Berry: iterative numerical solution method
  case('bbAssimFnc'      ); get_ixdecisions=iLookDECISIONS%bbAssimFnc  ! Ball-Berry: controls on carbon assimilation
  case('bbCanIntg8'      ); get_ixdecisions=iLookDECISIONS%bbCanIntg8  ! Ball-Berry: scaling of photosynthesis from the leaf to the canopy
  case('num_method'      ); get_ixdecisions=iLookDECISIONS%num_method  ! choice of numerical method
  case('fDerivMeth'      ); get_ixdecisions=iLookDECISIONS%fDerivMeth  ! choice of method to calculate flux derivatives
  case('LAI_method'      ); get_ixdecisions=iLookDECISIONS%LAI_method  ! choice of method to determine LAI and SAI
  case('cIntercept'      ); get_ixdecisions=iLookDECISIONS%cIntercept  ! choice of parameterization for canopy interception
  case('f_Richards'      ); get_ixdecisions=iLookDECISIONS%f_Richards  ! form of Richards' equation
  case('groundwatr'      ); get_ixdecisions=iLookDECISIONS%groundwatr  ! choice of groundwater parameterization
  case('hc_profile'      ); get_ixdecisions=iLookDECISIONS%hc_profile  ! choice of hydraulic conductivity profile
  case('bcUpprTdyn'      ); get_ixdecisions=iLookDECISIONS%bcUpprTdyn  ! type of upper boundary condition for thermodynamics
  case('bcLowrTdyn'      ); get_ixdecisions=iLookDECISIONS%bcLowrTdyn  ! type of lower boundary condition for thermodynamics
  case('bcUpprSoiH'      ); get_ixdecisions=iLookDECISIONS%bcUpprSoiH  ! type of upper boundary condition for soil hydrology
  case('bcLowrSoiH'      ); get_ixdecisions=iLookDECISIONS%bcLowrSoiH  ! type of lower boundary condition for soil hydrology
  case('veg_traits'      ); get_ixdecisions=iLookDECISIONS%veg_traits  ! choice of parameterization for vegetation roughness length and displacement height
  case('rootProfil'      ); get_ixdecisions=iLookDECISIONS%rootProfil  ! choice of parameterization for the rooting profile
  case('canopyEmis'      ); get_ixdecisions=iLookDECISIONS%canopyEmis  ! choice of parameterization for canopy emissivity
  case('snowIncept'      ); get_ixdecisions=iLookDECISIONS%snowIncept  ! choice of parameterization for snow interception
  case('windPrfile'      ); get_ixdecisions=iLookDECISIONS%windPrfile  ! choice of canopy wind profile
  case('astability'      ); get_ixdecisions=iLookDECISIONS%astability  ! choice of stability function
  case('compaction'      ); get_ixdecisions=iLookDECISIONS%compaction  ! choice of compaction routine
  case('snowLayers'      ); get_ixdecisions=iLookDECISIONS%snowLayers  ! choice of method to combine and sub-divide snow layers
  case('thCondSnow'      ); get_ixdecisions=iLookDECISIONS%thCondSnow  ! choice of thermal conductivity representation for snow
  case('thCondSoil'      ); get_ixdecisions=iLookDECISIONS%thCondSoil  ! choice of thermal conductivity representation for soil
  case('canopySrad'      ); get_ixdecisions=iLookDECISIONS%canopySrad  ! choice of method for canopy shortwave radiation
  case('alb_method'      ); get_ixdecisions=iLookDECISIONS%alb_method  ! choice of albedo representation
  case('spatial_gw'      ); get_ixdecisions=iLookDECISIONS%spatial_gw  ! choice of method for spatial representation of groundwater
  case('subRouting'      ); get_ixdecisions=iLookDECISIONS%subRouting  ! choice of method for sub-grid routing
  case('snowDenNew'      ); get_ixdecisions=iLookDECISIONS%snowDenNew  ! choice of method for new snow density
  case('snowUnload'      ); get_ixdecisions=iLookDECISIONS%snowUnload  ! choice of parameterization for snow unloading from canopy
  case('nrgConserv'      ); get_ixdecisions=iLookDECISIONS%nrgConserv  ! choice of variable in either energy backward Euler residual or IDA state variable
  case('aquiferIni'      ); get_ixdecisions=iLookDECISIONS%aquiferIni  ! choice of full or empty aquifer at start
  ! get to here if cannot find the variable
  case default
   get_ixdecisions = integerMissing
 end select
 end function get_ixdecisions


 ! *******************************************************************************************************************
 ! public function get_ixTime: get the index of the named variables for the model time
 ! *******************************************************************************************************************
 function get_ixTime(varName)
 USE var_lookup,only:iLookTIME                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixTime              ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  case('iyyy'            ); get_ixTime = iLookTIME%iyyy             ! year
  case('im'              ); get_ixTime = iLookTIME%im               ! month
  case('id'              ); get_ixTime = iLookTIME%id               ! day
  case('ih'              ); get_ixTime = iLookTIME%ih               ! hour
  case('imin'            ); get_ixTime = iLookTIME%imin             ! minute
  case('ih_tz'           ); get_ixTime = iLookTIME%ih_tz            ! hour for time zone offset
  case('imin_tz'         ); get_ixTime = iLookTIME%imin_tz          ! minute for time zone offset
  ! get to here if cannot find the variable
  case default
   get_ixTime = integerMissing
 end select
 end function get_ixTime


 ! *******************************************************************************************************************
 ! public function get_ixForce: get the index of the named variables for the model forcing data
 ! *******************************************************************************************************************
 function get_ixForce(varName)
 USE var_lookup,only:iLookFORCE                      ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixForce             ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  case('time'            ); get_ixForce = iLookFORCE%time             ! time since time reference       (s)
  case('pptrate'         ); get_ixForce = iLookFORCE%pptrate          ! precipitation rate              (kg m-2 s-1)
  case('airtemp'         ); get_ixForce = iLookFORCE%airtemp          ! air temperature                 (K)
  case('spechum'         ); get_ixForce = iLookFORCE%spechum          ! specific humidity               (g/g)
  case('windspd'         ); get_ixForce = iLookFORCE%windspd          ! windspeed                       (m/s)
  case('SWRadAtm'        ); get_ixForce = iLookFORCE%SWRadAtm         ! downwelling shortwave radiaiton (W m-2)
  case('LWRadAtm'        ); get_ixForce = iLookFORCE%LWRadAtm         ! downwelling longwave radiation  (W m-2)
  case('airpres'         ); get_ixForce = iLookFORCE%airpres          ! pressure                        (Pa)
  ! get to here if cannot find the variable
  case default
   get_ixForce = integerMissing
 end select
 end function get_ixForce


 ! *******************************************************************************************************************
 ! public function get_ixAttr: get the index of the named variables for the site characteristics
 ! *******************************************************************************************************************
 function get_ixAttr(varName)
 USE var_lookup,only:iLookATTR                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixAttr              ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  case('latitude'      ); get_ixAttr = iLookATTR%latitude       ! latitude (degrees north)
  case('longitude'     ); get_ixAttr = iLookATTR%longitude      ! longitude (degrees east)
  case('elevation'     ); get_ixAttr = iLookATTR%elevation      ! elevation (m)
  case('tan_slope'     ); get_ixAttr = iLookATTR%tan_slope      ! tan water table slope, taken as tan local ground surface slope (-)
  case('contourLength' ); get_ixAttr = iLookATTR%contourLength  ! length of contour at downslope edge of HRU (m)
  case('HRUarea'       ); get_ixAttr = iLookATTR%HRUarea        ! area of each HRU (m2)
  case('mHeight'       ); get_ixAttr = iLookATTR%mHeight        ! measurement height above bare ground (m)
  case('aspect'        ); get_ixAttr = iLookATTR%aspect         ! azimuth in degrees East of North (degrees)
  ! get to here if cannot find the variable
  case default
   get_ixAttr = integerMissing
 end select
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
 ! get the index of the named variables
 select case(trim(varName))
  case('vegTypeIndex'   ); get_ixType = iLookTYPE%vegTypeIndex       ! index defining vegetation type
  case('soilTypeIndex'  ); get_ixType = iLookTYPE%soilTypeIndex      ! index defining soil type
  case('slopeTypeIndex' ); get_ixType = iLookTYPE%slopeTypeIndex     ! index defining slope
  case('downHRUindex'   ); get_ixType = iLookTYPE%downHRUindex       ! index of downslope HRU (0 = basin outlet)
  ! get to here if cannot find the variable
  case default
   get_ixType = integerMissing
 end select
 end function get_ixType


 ! *******************************************************************************************************************
 ! public function get_ixId: get the index of the named variables for hru and gru IDs and related information
 ! *******************************************************************************************************************
 function get_ixId(varName)
 USE var_lookup,only:iLookID                         ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixId                ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  case('hruId'          ); get_ixId = iLookID%hruId              ! id defining HRU index
  ! get to here if cannot find the variable
  case default
   get_ixId = integerMissing
 end select
 end function get_ixId


 ! *******************************************************************************************************************
 ! public function get_ixParam: get the index of the named variables for the model parameters
 ! *******************************************************************************************************************
 function get_ixParam(varName)
 USE var_lookup,only:iLookPARAM                      ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixParam             ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  ! boundary conditions
  case('upperBoundHead'           ); get_ixParam = iLookPARAM%upperBoundHead         ! matric head of the upper boundary (m)
  case('lowerBoundHead'           ); get_ixParam = iLookPARAM%lowerBoundHead         ! matric head of the lower boundary (m)
  case('upperBoundTheta'          ); get_ixParam = iLookPARAM%upperBoundTheta        ! volumetric liquid water content at the upper boundary (-)
  case('lowerBoundTheta'          ); get_ixParam = iLookPARAM%lowerBoundTheta        ! volumetric liquid water content at the lower boundary (-)
  case('upperBoundTemp'           ); get_ixParam = iLookPARAM%upperBoundTemp         ! temperature of the upper boundary (K)
  case('lowerBoundTemp'           ); get_ixParam = iLookPARAM%lowerBoundTemp         ! temperature of the lower boundary (K)
  ! precipitation partitioning
  case('tempCritRain'             ); get_ixParam = iLookPARAM%tempCritRain           ! critical temperature where precipitation is rain (K)
  case('tempRangeTimestep'        ); get_ixParam = iLookPARAM%tempRangeTimestep      ! temperature range over the time step (K)
  case('frozenPrecipMultip'       ); get_ixParam = iLookPARAM%frozenPrecipMultip     ! frozen precipitation multiplier (-)
  ! freezing curve for snow
  case('snowfrz_scale'            ); get_ixParam = iLookPARAM%snowfrz_scale          ! scaling parameter for the freezing curve for snow (K-1)
  case('fixedThermalCond_snow'    ); get_ixParam = iLookPARAM%fixedThermalCond_snow  ! temporally constant thermal conductivity for snow (W m-1 K-1)
  ! snow albedo
  case('albedoMax'                ); get_ixParam = iLookPARAM%albedoMax              ! maximum snow albedo for a single spectral band (-)
  case('albedoMinWinter'          ); get_ixParam = iLookPARAM%albedoMinWinter        ! minimum snow albedo during winter for a single spectral band (-)
  case('albedoMinSpring'          ); get_ixParam = iLookPARAM%albedoMinSpring        ! minimum snow albedo during spring for a single spectral band (-)
  case('albedoMaxVisible'         ); get_ixParam = iLookPARAM%albedoMaxVisible       ! maximum snow albedo in the visible part of the spectrum (-)
  case('albedoMinVisible'         ); get_ixParam = iLookPARAM%albedoMinVisible       ! minimum snow albedo in the visible part of the spectrum (-)
  case('albedoMaxNearIR'          ); get_ixParam = iLookPARAM%albedoMaxNearIR        ! maximum snow albedo in the near infra-red part of the spectrum (-)
  case('albedoMinNearIR'          ); get_ixParam = iLookPARAM%albedoMinNearIR        ! minimum snow albedo in the near infra-red part of the spectrum (-)
  case('albedoDecayRate'          ); get_ixParam = iLookPARAM%albedoDecayRate        ! albedo decay rate (s)
  case('albedoSootLoad'           ); get_ixParam = iLookPARAM%albedoSootLoad         ! soot load factor (-)
  case('albedoRefresh'            ); get_ixParam = iLookPARAM%albedoRefresh          ! critical mass necessary for albedo refreshment (kg m-2)
  ! radiation transfer
  case('radExt_snow'              ); get_ixParam = iLookPARAM%radExt_snow            ! extinction coefficient for radiation penetration within the snowpack (m-1)
  case('directScale'              ); get_ixParam = iLookPARAM%directScale            ! scaling factor for fractional driect radiaion parameterization (-)
  case('Frad_direct'              ); get_ixParam = iLookPARAM%Frad_direct            ! maximum fraction of direct radiation (-)
  case('Frad_vis'                 ); get_ixParam = iLookPARAM%Frad_vis               ! fraction of radiation in the visible part of the spectrum (-)
  ! new snow density
  case('newSnowDenMin'            ); get_ixParam = iLookPARAM%newSnowDenMin          ! minimum new snow density (kg m-3)
  case('newSnowDenMult'           ); get_ixParam = iLookPARAM%newSnowDenMult         ! multiplier for new snow density (kg m-3)
  case('newSnowDenScal'           ); get_ixParam = iLookPARAM%newSnowDenScal         ! scaling factor for new snow density (K)
  case('constSnowDen'             ); get_ixParam = iLookPARAM%constSnowDen           ! Constant new snow density (kg m-3)
  case('newSnowDenAdd'            ); get_ixParam = iLookPARAM%newSnowDenAdd          ! Pahaut 1976, additive factor for new snow density (kg m-3)
  case('newSnowDenMultTemp'       ); get_ixParam = iLookPARAM%newSnowDenMultTemp     ! Pahaut 1976, multiplier for new snow density applied to air temperature (kg m-3 K-1)
  case('newSnowDenMultWind'       ); get_ixParam = iLookPARAM%newSnowDenMultWind     ! Pahaut 1976, multiplier for new snow density applied to wind speed (kg m-7/2 s-1/2)
  case('newSnowDenMultAnd'        ); get_ixParam = iLookPARAM%newSnowDenMultAnd      ! Anderson 1976, multiplier for new snow density for Anderson function (K-1)
  case('newSnowDenBase'           ); get_ixParam = iLookPARAM%newSnowDenBase         ! Anderson 1976, base value that is rasied to the (3/2) power (K)
  ! snow compaction
  case('densScalGrowth'           ); get_ixParam = iLookPARAM%densScalGrowth         ! density scaling factor for grain growth (kg-1 m3)
  case('tempScalGrowth'           ); get_ixParam = iLookPARAM%tempScalGrowth         ! temperature scaling factor for grain growth (K-1)
  case('grainGrowthRate'          ); get_ixParam = iLookPARAM%grainGrowthRate        ! rate of grain growth (s-1)
  case('densScalOvrbdn'           ); get_ixParam = iLookPARAM%densScalOvrbdn         ! density scaling factor for overburden pressure (kg-1 m3)
  case('tempScalOvrbdn'           ); get_ixParam = iLookPARAM%tempScalOvrbdn         ! temperature scaling factor for overburden pressure (K-1)
  case('baseViscosity'            ); get_ixParam = iLookPARAM%baseViscosity          ! viscosity coefficient at T=T_frz and snow density=0  (kg s m-2)
  ! water flow through snow
  case('Fcapil'                   ); get_ixParam = iLookPARAM%Fcapil                 ! capillary retention as a fraction of the total pore volume (-)
  case('k_snow'                   ); get_ixParam = iLookPARAM%k_snow                 ! hydraulic conductivity of snow (m s-1), 0.0055 = approx. 20 m/hr, from UEB
  case('mw_exp'                   ); get_ixParam = iLookPARAM%mw_exp                 ! exponent for meltwater flow (-)
  ! turbulent heat fluxes
  case('z0Snow'                   ); get_ixParam = iLookPARAM%z0Snow                 ! roughness length of snow (m)
  case('z0Soil'                   ); get_ixParam = iLookPARAM%z0Soil                 ! roughness length of bare soil below the canopy (m)
  case('z0Canopy'                 ); get_ixParam = iLookPARAM%z0Canopy               ! roughness length of the canopy (m)
  case('zpdFraction'              ); get_ixParam = iLookPARAM%zpdFraction            ! zero plane displacement / canopy height (-)
  case('critRichNumber'           ); get_ixParam = iLookPARAM%critRichNumber         ! critical value for the bulk Richardson number (-)
  case('Louis79_bparam'           ); get_ixParam = iLookPARAM%Louis79_bparam         ! parameter in Louis (1979) stability function (-)
  case('Louis79_cStar'            ); get_ixParam = iLookPARAM%Louis79_cStar          ! parameter in Louis (1979) stability function (-)
  case('Mahrt87_eScale'           ); get_ixParam = iLookPARAM%Mahrt87_eScale         ! exponential scaling factor in the Mahrt (1987) stability function (-)
  case('leafExchangeCoeff'        ); get_ixParam = iLookPARAM%leafExchangeCoeff      ! turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
  case('windReductionParam'       ); get_ixParam = iLookPARAM%windReductionParam     ! canopy wind reduction parameter (-)
  ! stomatal conductance
  case('Kc25'                     ); get_ixParam = iLookPARAM%Kc25                   ! Michaelis-Menten constant for CO2 at 25 degrees C (umol mol-1)
  case('Ko25'                     ); get_ixParam = iLookPARAM%Ko25                   ! Michaelis-Menten constant for O2 at 25 degrees C (mol mol-1)
  case('Kc_qFac'                  ); get_ixParam = iLookPARAM%Kc_qFac                ! factor in the q10 function defining temperature controls on Kc (-)
  case('Ko_qFac'                  ); get_ixParam = iLookPARAM%Ko_qFac                ! factor in the q10 function defining temperature controls on Ko (-)
  case('kc_Ha'                    ); get_ixParam = iLookPARAM%kc_Ha                  ! activation energy for the Michaelis-Menten constant for CO2 (J mol-1)
  case('ko_Ha'                    ); get_ixParam = iLookPARAM%ko_Ha                  ! activation energy for the Michaelis-Menten constant for CO2 (J mol-1)
  case('vcmax25_canopyTop'        ); get_ixParam = iLookPARAM%vcmax25_canopyTop      ! potential carboxylation rate at 25 degrees C at the canopy top (umol co2 m-2 s-1)
  case('vcmax_qFac'               ); get_ixParam = iLookPARAM%vcmax_qFac             ! factor in the q10 function defining temperature controls on vcmax (-)
  case('vcmax_Ha'                 ); get_ixParam = iLookPARAM%vcmax_Ha               ! activation energy in the vcmax function (J mol-1)
  case('vcmax_Hd'                 ); get_ixParam = iLookPARAM%vcmax_Hd               ! deactivation energy in the vcmax function (J mol-1)
  case('vcmax_Sv'                 ); get_ixParam = iLookPARAM%vcmax_Sv               ! entropy term in the vcmax function (J mol-1 K-1)
  case('vcmax_Kn'                 ); get_ixParam = iLookPARAM%vcmax_Kn               ! foliage nitrogen decay coefficient (-)
  case('jmax25_scale'             ); get_ixParam = iLookPARAM%jmax25_scale           ! scaling factor to relate jmax25 to vcmax25 (-)
  case('jmax_Ha'                  ); get_ixParam = iLookPARAM%jmax_Ha                ! activation energy in the jmax function (J mol-1)
  case('jmax_Hd'                  ); get_ixParam = iLookPARAM%jmax_Hd                ! deactivation energy in the jmax function (J mol-1)
  case('jmax_Sv'                  ); get_ixParam = iLookPARAM%jmax_Sv                ! entropy term in the jmax function (J mol-1 K-1)
  case('fractionJ'                ); get_ixParam = iLookPARAM%fractionJ              ! fraction of light lost by other than the chloroplast lamellae (-)
  case('quantamYield'             ); get_ixParam = iLookPARAM%quantamYield           ! quantam yield (mol e mol-1 quanta)
  case('vpScaleFactor'            ); get_ixParam = iLookPARAM%vpScaleFactor          ! vapor pressure scaling factor in stomatal conductance function (Pa)
  case('cond2photo_slope'         ); get_ixParam = iLookPARAM%cond2photo_slope       ! slope of conductance-photosynthesis relationship (-)
  case('minStomatalConductance'   ); get_ixParam = iLookPARAM%minStomatalConductance ! minimum stomatal conductance (umol H2O m-2 s-1)
  ! vegetation properties
  case('winterSAI'                ); get_ixParam = iLookPARAM%winterSAI              ! stem area index prior to the start of the growing season (m2 m-2)
  case('summerLAI'                ); get_ixParam = iLookPARAM%summerLAI              ! maximum leaf area index at the peak of the growing season (m2 m-2)
  case('rootScaleFactor1'         ); get_ixParam = iLookPARAM%rootScaleFactor1       ! 1st scaling factor (a) in Y = 1 - 0.5*( exp(-aZ) + exp(-bZ) )   (m-1)
  case('rootScaleFactor2'         ); get_ixParam = iLookPARAM%rootScaleFactor2       ! 2nd scaling factor (b) in Y = 1 - 0.5*( exp(-aZ) + exp(-bZ) )   (m-1)
  case('rootingDepth'             ); get_ixParam = iLookPARAM%rootingDepth           ! rooting depth (m)
  case('rootDistExp'              ); get_ixParam = iLookPARAM%rootDistExp            ! exponent for the vertical distriution of root density (-)
  case('plantWiltPsi'             ); get_ixParam = iLookPARAM%plantWiltPsi           ! matric head at wilting point (m)
  case('soilStressParam'          ); get_ixParam = iLookPARAM%soilStressParam        ! parameter in the exponential soil stress function
  case('critSoilWilting'          ); get_ixParam = iLookPARAM%critSoilWilting        ! critical vol. liq. water content when plants are wilting (-)
  case('critSoilTranspire'        ); get_ixParam = iLookPARAM%critSoilTranspire      ! critical vol. liq. water content when transpiration is limited (-)
  case('critAquiferTranspire'     ); get_ixParam = iLookPARAM%critAquiferTranspire   ! critical aquifer storage value when transpiration is limited (m)
  case('minStomatalResistance'    ); get_ixParam = iLookPARAM%minStomatalResistance  ! minimum canopy resistance (s m-1)
  case('leafDimension'            ); get_ixParam = iLookPARAM%leafDimension          ! characteristic leaf dimension (m)
  case('heightCanopyTop'          ); get_ixParam = iLookPARAM%heightCanopyTop        ! height of top of the vegetation canopy above ground surface (m)
  case('heightCanopyBottom'       ); get_ixParam = iLookPARAM%heightCanopyBottom     ! height of bottom of the vegetation canopy above ground surface (m)
  case('specificHeatVeg'          ); get_ixParam = iLookPARAM%specificHeatVeg        ! specific heat of vegetation (J kg-1 K-1)
  case('maxMassVegetation'        ); get_ixParam = iLookPARAM%maxMassVegetation      ! maximum mass of vegetation (full foliage) (kg m-2)
  case('throughfallScaleSnow'     ); get_ixParam = iLookPARAM%throughfallScaleSnow   ! scaling factor for throughfall (snow) (-)
  case('throughfallScaleRain'     ); get_ixParam = iLookPARAM%throughfallScaleRain   ! scaling factor for throughfall (rain) (-)
  case('refInterceptCapSnow'      ); get_ixParam = iLookPARAM%refInterceptCapSnow    ! reference canopy interception capacity per unit leaf area (snow) (kg m-2)
  case('refInterceptCapRain'      ); get_ixParam = iLookPARAM%refInterceptCapRain    ! canopy interception capacity per unit leaf area (rain) (kg m-2)
  case('snowUnloadingCoeff'       ); get_ixParam = iLookPARAM%snowUnloadingCoeff     ! time constant for unloading of snow from the forest canopy (s-1)
  case('canopyDrainageCoeff'      ); get_ixParam = iLookPARAM%canopyDrainageCoeff    ! time constant for drainage of liquid water from the forest canopy (s-1)
  case('ratioDrip2Unloading'      ); get_ixParam = iLookPARAM%ratioDrip2Unloading    ! ratio of canopy drip to unloading of snow from the forest canopy (-)
  case('canopyWettingFactor'      ); get_ixParam = iLookPARAM%canopyWettingFactor    ! maximum wetted fraction of the canopy (-)
  case('canopyWettingExp'         ); get_ixParam = iLookPARAM%canopyWettingExp       ! exponent in canopy wetting function (-)
  case('minTempUnloading'         ); get_ixParam = iLookPARAM%minTempUnloading       ! min temp for unloading in windySnow (K)
  case('rateTempUnloading'        ); get_ixParam = iLookPARAM%rateTempUnloading      ! how quickly to unload due to temperature (K s)
  case('minWindUnloading'         ); get_ixParam = iLookPARAM%minWindUnloading       ! min wind speed for unloading in windySnow (m s-1)
  case('rateWindUnloading'        ); get_ixParam = iLookPARAM%rateWindUnloading      ! how quickly to unload due to wind (m)
  ! soil properties
  case('soil_dens_intr'           ); get_ixParam = iLookPARAM%soil_dens_intr         ! intrinsic soil density (kg m-3)
  case('thCond_soil'              ); get_ixParam = iLookPARAM%thCond_soil            ! thermal conductivity of soil (W m-1 K-1)
  case('frac_sand'                ); get_ixParam = iLookPARAM%frac_sand              ! fraction of sand (-)
  case('frac_silt'                ); get_ixParam = iLookPARAM%frac_silt              ! fraction of silt (-)
  case('frac_clay'                ); get_ixParam = iLookPARAM%frac_clay              ! fraction of clay (-)
  case('fieldCapacity'            ); get_ixParam = iLookPARAM%fieldCapacity          ! field capacity (-)
  case('wettingFrontSuction'      ); get_ixParam = iLookPARAM%wettingFrontSuction    ! Green-Ampt wetting front suction (m)
  case('theta_mp'                 ); get_ixParam = iLookPARAM%theta_mp               ! volumetric liquid water content when macropore flow begins (-)
  case('theta_sat'                ); get_ixParam = iLookPARAM%theta_sat              ! soil porosity (-)
  case('theta_res'                ); get_ixParam = iLookPARAM%theta_res              ! volumetric residual water content (-)
  case('vGn_alpha'                ); get_ixParam = iLookPARAM%vGn_alpha              ! van Genuchten "alpha" parameter (m-1)
  case('vGn_n'                    ); get_ixParam = iLookPARAM%vGn_n                  ! van Genuchten "n" parameter (-)
  case('mpExp'                    ); get_ixParam = iLookPARAM%mpExp                  ! empirical exponent in macropore flow equation (-)
  case('k_soil'                   ); get_ixParam = iLookPARAM%k_soil                 ! saturated hydraulic conductivity (m s-1)
  case('k_macropore'              ); get_ixParam = iLookPARAM%k_macropore            ! saturated hydraulic conductivity for the macropores (m s-1)
  case('kAnisotropic'             ); get_ixParam = iLookPARAM%kAnisotropic           ! anisotropy factor for lateral hydraulic conductivity (-)
  case('zScale_TOPMODEL'          ); get_ixParam = iLookPARAM%zScale_TOPMODEL        ! TOPMODEL scaling factor used in lower boundary condition for soil (m)
  case('compactedDepth'           ); get_ixParam = iLookPARAM%compactedDepth         ! depth where k_soil reaches the compacted value given by CH78 (m)
  case('aquiferBaseflowRate'      ); get_ixParam = iLookPARAM%aquiferBaseflowRate    ! baseflow rate when aquifer storage = aquiferScaleFactor (m s-1)
  case('aquiferScaleFactor'       ); get_ixParam = iLookPARAM%aquiferScaleFactor     ! scaling factor for aquifer storage in the big bucket (m)
  case('aquiferBaseflowExp'       ); get_ixParam = iLookPARAM%aquiferBaseflowExp     ! baseflow exponent (-)
  case('qSurfScale'               ); get_ixParam = iLookPARAM%qSurfScale             ! scaling factor in the surface runoff parameterization (-)
  case('specificYield'            ); get_ixParam = iLookPARAM%specificYield          ! specific yield (-)
  case('specificStorage'          ); get_ixParam = iLookPARAM%specificStorage        ! specific storage coefficient (m-1)
  case('f_impede'                 ); get_ixParam = iLookPARAM%f_impede               ! ice impedence factor (-)
  case('soilIceScale'             ); get_ixParam = iLookPARAM%soilIceScale           ! scaling factor for depth of soil ice, used to get frozen fraction (m)
  case('soilIceCV'                ); get_ixParam = iLookPARAM%soilIceCV              ! CV of depth of soil ice, used to get frozen fraction (-)
  ! algorithmic control parameters
  case('minwind'                  ); get_ixParam = iLookPARAM%minwind                ! minimum wind speed (m s-1)
  case('minstep'                  ); get_ixParam = iLookPARAM%minstep                ! minimum length of the time step homegrown, not currently used
  case('maxstep'                  ); get_ixParam = iLookPARAM%maxstep                ! maximum length of the time step homegrown
  case('be_steps'                 ); get_ixParam = iLookPARAM%be_steps               ! minimum number of substeps to take in a maxstep homegrown
  case('wimplicit'                ); get_ixParam = iLookPARAM%wimplicit              ! weight assigned to start-of-step fluxes homegrown, not currently used
  case('maxiter'                  ); get_ixParam = iLookPARAM%maxiter                ! maximum number of iterations homegrown and kinsol
  case('relConvTol_liquid'        ); get_ixParam = iLookPARAM%relConvTol_liquid      ! relative convergence tolerance for vol frac liq water (-) homegrown
  case('absConvTol_liquid'        ); get_ixParam = iLookPARAM%absConvTol_liquid      ! absolute convergence tolerance for vol frac liq water (-) homegrown
  case('relConvTol_matric'        ); get_ixParam = iLookPARAM%relConvTol_matric      ! relative convergence tolerance for matric head (-) homegrown
  case('absConvTol_matric'        ); get_ixParam = iLookPARAM%absConvTol_matric      ! absolute convergence tolerance for matric head (m) homegrown
  case('relConvTol_energy'        ); get_ixParam = iLookPARAM%relConvTol_energy      ! relative convergence tolerance for energy (-) homegrown
  case('absConvTol_energy'        ); get_ixParam = iLookPARAM%absConvTol_energy      ! absolute convergence tolerance for energy (J m-3) homegrown
  case('relConvTol_aquifr'        ); get_ixParam = iLookPARAM%relConvTol_aquifr      ! relative convergence tolerance for aquifer storage (-) homegrown
  case('absConvTol_aquifr'        ); get_ixParam = iLookPARAM%absConvTol_aquifr      ! absolute convergence tolerance for aquifer storage (m) homegrown
  case('relTolTempCas'            ); get_ixParam = iLookPARAM%relTolTempCas          ! relative error tolerance for canopy temperature state variable
  case('absTolTempCas'            ); get_ixParam = iLookPARAM%absTolTempCas          ! absolute error tolerance for canopy temperature state variable
  case('relTolTempVeg'            ); get_ixParam = iLookPARAM%relTolTempVeg          ! relative error tolerance for vegitation temp state var
  case('absTolTempVeg'            ); get_ixParam = iLookPARAM%absTolTempVeg          ! absolute error tolerance for vegitation temp state var
  case('relTolWatVeg'             ); get_ixParam = iLookPARAM%relTolWatVeg           ! relative error tolerance for vegitation hydrology
  case('absTolWatVeg'             ); get_ixParam = iLookPARAM%absTolWatVeg           ! absolute error tolerance for vegitation hydrology
  case('relTolTempSoilSnow'       ); get_ixParam = iLookPARAM%relTolTempSoilSnow     ! relative error tolerance for snow+soil energy
  case('absTolTempSoilSnow'       ); get_ixParam = iLookPARAM%absTolTempSoilSnow     ! absolute error tolerance for snow+soil energy
  case('relTolWatSnow'            ); get_ixParam = iLookPARAM%relTolWatSnow          ! relative error tolerance for snow hydrology
  case('absTolWatSnow'            ); get_ixParam = iLookPARAM%absTolWatSnow          ! absolute error tolerance for snow hydrology
  case('relTolMatric'             ); get_ixParam = iLookPARAM%relTolMatric           ! relative error tolerance for matric head
  case('absTolMatric'             ); get_ixParam = iLookPARAM%absTolMatric           ! absolute error tolerance for matric head
  case('relTolAquifr'             ); get_ixParam = iLookPARAM%relTolAquifr           ! relative error tolerance for aquifer hydrology
  case('absTolAquifr'             ); get_ixParam = iLookPARAM%absTolAquifr           ! absolute error tolerance for aquifer hydrology
  case('idaMaxOrder'              ); get_ixParam = iLookPARAM%idaMaxOrder            ! maximum order for IDA  
  case('idaMaxInternalSteps'      ); get_ixParam = iLookPARAM%idaMaxInternalSteps    ! maximum number of internal steps for IDA before tout 
  case('idaInitStepSize'          ); get_ixParam = iLookPARAM%idaInitStepSize        ! initial step size for IDA 
  case('idaMinStepSize'           ); get_ixParam = iLookPARAM%idaMinStepSize         ! minimum step size for IDA
  case('idaMaxStepSize'           ); get_ixParam = iLookPARAM%idaMaxStepSize         ! maximum step size for IDA
  case('idaMaxErrTestFail'        ); get_ixParam = iLookPARAM%idaMaxErrTestFail      ! maximum number of error test failures for IDA
  case('zmin'                     ); get_ixParam = iLookPARAM%zmin                   ! minimum layer depth (m)
  case('zmax'                     ); get_ixParam = iLookPARAM%zmax                   ! maximum layer depth (m)
  case('zminLayer1'               ); get_ixParam = iLookPARAM%zminLayer1             ! minimum layer depth for the 1st (top) layer (m)
  case('zminLayer2'               ); get_ixParam = iLookPARAM%zminLayer2             ! minimum layer depth for the 2nd layer (m)
  case('zminLayer3'               ); get_ixParam = iLookPARAM%zminLayer3             ! minimum layer depth for the 3rd layer (m)
  case('zminLayer4'               ); get_ixParam = iLookPARAM%zminLayer4             ! minimum layer depth for the 4th layer (m)
  case('zminLayer5'               ); get_ixParam = iLookPARAM%zminLayer5             ! minimum layer depth for the 5th (bottom) layer (m)
  case('zmaxLayer1_lower'         ); get_ixParam = iLookPARAM%zmaxLayer1_lower       ! maximum layer depth for the 1st (top) layer when only 1 layer (m)
  case('zmaxLayer2_lower'         ); get_ixParam = iLookPARAM%zmaxLayer2_lower       ! maximum layer depth for the 2nd layer when only 2 layers (m)
  case('zmaxLayer3_lower'         ); get_ixParam = iLookPARAM%zmaxLayer3_lower       ! maximum layer depth for the 3rd layer when only 3 layers (m)
  case('zmaxLayer4_lower'         ); get_ixParam = iLookPARAM%zmaxLayer4_lower       ! maximum layer depth for the 4th layer when only 4 layers (m)
  case('zmaxLayer1_upper'         ); get_ixParam = iLookPARAM%zmaxLayer1_upper       ! maximum layer depth for the 1st (top) layer when > 1 layer (m)
  case('zmaxLayer2_upper'         ); get_ixParam = iLookPARAM%zmaxLayer2_upper       ! maximum layer depth for the 2nd layer when > 2 layers (m)
  case('zmaxLayer3_upper'         ); get_ixParam = iLookPARAM%zmaxLayer3_upper       ! maximum layer depth for the 3rd layer when > 3 layers (m)
  case('zmaxLayer4_upper'         ); get_ixParam = iLookPARAM%zmaxLayer4_upper       ! maximum layer depth for the 4th layer when > 4 layers (m)
  ! get to here if cannot find the variable
  case default
   get_ixParam = integerMissing
 end select
 end function get_ixParam


 ! *******************************************************************************************************************
 ! public function get_ixProg: get the index of the named variables for the prognostic (state) variables
 ! *******************************************************************************************************************
 function get_ixProg(varName)
 USE var_lookup,only:iLookPROG                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixProg              ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  ! variables for time stepping
  case('dt_init'                        ); get_ixProg = iLookPROG%dt_init                          ! length of initial time step at start of next data interval (s)
  ! state variables for vegetation
  case('scalarCanopyIce'                ); get_ixProg = iLookPROG%scalarCanopyIce                  ! mass of ice on the vegetation canopy (kg m-2)
  case('scalarCanopyLiq'                ); get_ixProg = iLookPROG%scalarCanopyLiq                  ! mass of liquid water on the vegetation canopy (kg m-2)
  case('scalarCanopyWat'                ); get_ixProg = iLookPROG%scalarCanopyWat                  ! mass of total water on the vegetation canopy (kg m-2)
  case('scalarCanairTemp'               ); get_ixProg = iLookPROG%scalarCanairTemp                 ! temperature of the canopy air space (K)
  case('scalarCanopyTemp'               ); get_ixProg = iLookPROG%scalarCanopyTemp                 ! temperature of the vegetation canopy (K)
  ! state variables for snow
  case('spectralSnowAlbedoDiffuse'      ); get_ixProg = iLookPROG%spectralSnowAlbedoDiffuse        ! diffuse snow albedo for individual spectral bands (-)
  case('scalarSnowAlbedo'               ); get_ixProg = iLookPROG%scalarSnowAlbedo                 ! snow albedo for the entire spectral band (-)
  case('scalarSnowDepth'                ); get_ixProg = iLookPROG%scalarSnowDepth                  ! total snow depth (m)
  case('scalarSWE'                      ); get_ixProg = iLookPROG%scalarSWE                        ! snow water equivalent (kg m-2)
  case('scalarSfcMeltPond'              ); get_ixProg = iLookPROG%scalarSfcMeltPond                ! ponded water caused by melt of the "snow without a layer" (kg m-2)
  ! state variables for the snow+soil domain
  case('mLayerTemp'                     ); get_ixProg = iLookPROG%mLayerTemp                       ! temperature of each layer (K)
  case('mLayerVolFracIce'               ); get_ixProg = iLookPROG%mLayerVolFracIce                 ! volumetric fraction of icein each layer (-)
  case('mLayerVolFracLiq'               ); get_ixProg = iLookPROG%mLayerVolFracLiq                 ! volumetric fraction of liquid water in each layer (-)
  case('mLayerVolFracWat'               ); get_ixProg = iLookPROG%mLayerVolFracWat                 ! volumetric fraction of total water in each layer (-)
  case('mLayerMatricHead'               ); get_ixProg = iLookPROG%mLayerMatricHead                 ! matric head of water in the soil (m)
  ! other state variables
  case('scalarAquiferStorage'           ); get_ixProg = iLookPROG%scalarAquiferStorage             ! relative aquifer storage -- above bottom of the soil profile (m)
  case('scalarSurfaceTemp'              ); get_ixProg = iLookPROG%scalarSurfaceTemp                ! surface temperature (K)
  ! coordinate variables
  case('mLayerDepth'                    ); get_ixProg = iLookPROG%mLayerDepth                      ! depth of each layer (m)
  case('mLayerHeight'                   ); get_ixProg = iLookPROG%mLayerHeight                     ! height at the midpoint of each layer (m)
  case('iLayerHeight'                   ); get_ixProg = iLookPROG%iLayerHeight                     ! height at the interface of each layer (m)
  ! get to here if cannot find the variable
  case default
   get_ixProg = integerMissing
 end select
 end function get_ixProg


 ! *******************************************************************************************************************
 ! public function get_ixDiag: get the index of the named variables for the diagnostic variables
 ! *******************************************************************************************************************
 function get_ixDiag(varName)
 USE var_lookup,only:iLookDIAG                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixDiag              ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  ! local properties
  case('scalarCanopyDepth'              ); get_ixDiag = iLookDIAG%scalarCanopyDepth                ! canopy depth (m)
  case('scalarGreenVegFraction'         ); get_ixDiag = iLookDIAG%scalarGreenVegFraction           ! green vegetation fraction used to compute LAI (-)
  case('scalarBulkVolHeatCapVeg'        ); get_ixDiag = iLookDIAG%scalarBulkVolHeatCapVeg          ! bulk volumetric heat capacity of vegetation (J m-3 K-1)
  case('scalarCanopyCm'                 ); get_ixDiag = iLookDIAG%scalarCanopyCm                   ! Cm of canopy (J kg-1 K-1)
  case('scalarCanopyEmissivity'         ); get_ixDiag = iLookDIAG%scalarCanopyEmissivity           ! effective canopy emissivity (-)
  case('scalarRootZoneTemp'             ); get_ixDiag = iLookDIAG%scalarRootZoneTemp               ! average temperature of the root zone (K)
  case('scalarLAI'                      ); get_ixDiag = iLookDIAG%scalarLAI                        ! one-sided leaf area index (m2 m-2)
  case('scalarSAI'                      ); get_ixDiag = iLookDIAG%scalarSAI                        ! one-sided stem area index (m2 m-2)
  case('scalarExposedLAI'               ); get_ixDiag = iLookDIAG%scalarExposedLAI                 ! exposed leaf area index after burial by snow (m2 m-2)
  case('scalarExposedSAI'               ); get_ixDiag = iLookDIAG%scalarExposedSAI                 ! exposed stem area index after burial by snow (m2 m-2)
  case('scalarAdjMeasHeight'            ); get_ixDiag = iLookDIAG%scalarAdjMeasHeight              ! adjusted measurement height for cases snowDepth>mHeight (m)
  case('scalarCanopyIceMax'             ); get_ixDiag = iLookDIAG%scalarCanopyIceMax               ! maximum interception storage capacity for ice (kg m-2)
  case('scalarCanopyLiqMax'             ); get_ixDiag = iLookDIAG%scalarCanopyLiqMax               ! maximum interception storage capacity for liquid water (kg m-2)
  case('scalarGrowingSeasonIndex'       ); get_ixDiag = iLookDIAG%scalarGrowingSeasonIndex         ! growing season index (0=off, 1=on)
  case('scalarVolHtCap_air'             ); get_ixDiag = iLookDIAG%scalarVolHtCap_air               ! volumetric heat capacity air (J m-3 K-1)
  case('scalarVolHtCap_ice'             ); get_ixDiag = iLookDIAG%scalarVolHtCap_ice               ! volumetric heat capacity ice (J m-3 K-1)
  case('scalarVolHtCap_soil'            ); get_ixDiag = iLookDIAG%scalarVolHtCap_soil              ! volumetric heat capacity dry soil (J m-3 K-1)
  case('scalarVolHtCap_water'           ); get_ixDiag = iLookDIAG%scalarVolHtCap_water             ! volumetric heat capacity liquid wat (J m-3 K-1)
  case('mLayerVolHtCapBulk'             ); get_ixDiag = iLookDIAG%mLayerVolHtCapBulk               ! volumetric heat capacity in each layer (J m-3 K-1)
  case('mLayerCm'                       ); get_ixDiag = iLookDIAG%mLayerCm                         ! Cm of each layer (J kg-1 K-1)
  case('scalarLambda_drysoil'           ); get_ixDiag = iLookDIAG%scalarLambda_drysoil             ! thermal conductivity of dry soil     (W m-1)
  case('scalarLambda_wetsoil'           ); get_ixDiag = iLookDIAG%scalarLambda_wetsoil             ! thermal conductivity of wet soil     (W m-1)
  case('mLayerThermalC'                 ); get_ixDiag = iLookDIAG%mLayerThermalC                   ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
  case('iLayerThermalC'                 ); get_ixDiag = iLookDIAG%iLayerThermalC                   ! thermal conductivity at the interface of each layer (W m-1 K-1)
  ! enthalpy
  case('scalarCanairEnthalpy'           ); get_ixDiag = iLookDIAG%scalarCanairEnthalpy             ! enthalpy of the canopy air space (J m-3)
  case('scalarCanopyEnthTemp'           ); get_ixDiag = iLookDIAG%scalarCanopyEnthTemp             ! temperature component of enthalpy of the vegetation canopy (J m-3)
  case('scalarCanopyEnthalpy'           ); get_ixDiag = iLookDIAG%scalarCanopyEnthalpy             ! enthalpy of the vegetation canopy (J m-3)
  case('mLayerEnthTemp'                 ); get_ixDiag = iLookDIAG%mLayerEnthTemp                   ! temperature component of enthalpy of the snow+soil layers (J m-3)
  case('mLayerEnthalpy'                 ); get_ixDiag = iLookDIAG%mLayerEnthalpy                   ! enthalpy of the snow+soil layers (J m-3)
  case('scalarTotalSoilEnthalpy'        ); get_ixDiag = iLookDIAG%scalarTotalSoilEnthalpy          ! total enthalpy of the soil column (J m-3)
  case('scalarTotalSnowEnthalpy'        ); get_ixDiag = iLookDIAG%scalarTotalSnowEnthalpy          ! total enthalpy of the snow column (J m-3)   
  ! forcing
  case('scalarVPair'                    ); get_ixDiag = iLookDIAG%scalarVPair                      ! vapor pressure of the air above the vegetation canopy (Pa)
  case('scalarVP_CanopyAir'             ); get_ixDiag = iLookDIAG%scalarVP_CanopyAir               ! vapor pressure of the canopy air space (Pa)
  case('scalarTwetbulb'                 ); get_ixDiag = iLookDIAG%scalarTwetbulb                   ! wetbulb temperature (K)
  case('scalarSnowfallTemp'             ); get_ixDiag = iLookDIAG%scalarSnowfallTemp               ! temperature of fresh snow (K)
  case('scalarNewSnowDensity'           ); get_ixDiag = iLookDIAG%scalarNewSnowDensity             ! density of fresh snow, should snow be falling in this time step (kg m-3)
  case('scalarO2air'                    ); get_ixDiag = iLookDIAG%scalarO2air                      ! atmospheric o2 concentration (Pa)
  case('scalarCO2air'                   ); get_ixDiag = iLookDIAG%scalarCO2air                     ! atmospheric co2 concentration (Pa)
  case('windspd_x'                      ); get_ixDiag = iLookDIAG%windspd_x                        ! wind speed at 10 meter height in x-direction (m s-1)
  case('windspd_y'                      ); get_ixDiag = iLookDIAG%windspd_y                        ! wind speed at 10 meter height in y-direction (m s-1)
  ! shortwave radiation
  case('scalarCosZenith'                ); get_ixDiag = iLookDIAG%scalarCosZenith                  ! cosine of the solar zenith angle (0-1)
  case('scalarFractionDirect'           ); get_ixDiag = iLookDIAG%scalarFractionDirect             ! fraction of direct radiation (0-1)
  case('scalarCanopySunlitFraction'     ); get_ixDiag = iLookDIAG%scalarCanopySunlitFraction       ! sunlit fraction of canopy (-)
  case('scalarCanopySunlitLAI'          ); get_ixDiag = iLookDIAG%scalarCanopySunlitLAI            ! sunlit leaf area (-)
  case('scalarCanopyShadedLAI'          ); get_ixDiag = iLookDIAG%scalarCanopyShadedLAI            ! shaded leaf area (-)
  case('spectralAlbGndDirect'           ); get_ixDiag = iLookDIAG%spectralAlbGndDirect             ! direct  albedo of underlying surface for each spectral band (-)
  case('spectralAlbGndDiffuse'          ); get_ixDiag = iLookDIAG%spectralAlbGndDiffuse            ! diffuse albedo of underlying surface for each spectral band (-)
  case('scalarGroundAlbedo'             ); get_ixDiag = iLookDIAG%scalarGroundAlbedo               ! albedo of the ground surface (-)
  ! turbulent heat transfer
  case('scalarLatHeatSubVapCanopy'      ); get_ixDiag = iLookDIAG%scalarLatHeatSubVapCanopy        ! latent heat of sublimation/vaporization used for veg canopy (J kg-1)
  case('scalarLatHeatSubVapGround'      ); get_ixDiag = iLookDIAG%scalarLatHeatSubVapGround        ! latent heat of sublimation/vaporization used for ground surface (J kg-1)
  case('scalarSatVP_CanopyTemp'         ); get_ixDiag = iLookDIAG%scalarSatVP_CanopyTemp           ! saturation vapor pressure at the temperature of vegetation canopy (Pa)
  case('scalarSatVP_GroundTemp'         ); get_ixDiag = iLookDIAG%scalarSatVP_GroundTemp           ! saturation vapor pressure at the temperature of the ground (Pa)
  case('scalarZ0Canopy'                 ); get_ixDiag = iLookDIAG%scalarZ0Canopy                   ! roughness length of the canopy (m)
  case('scalarWindReductionFactor'      ); get_ixDiag = iLookDIAG%scalarWindReductionFactor        ! canopy wind reduction factor (-)
  case('scalarZeroPlaneDisplacement'    ); get_ixDiag = iLookDIAG%scalarZeroPlaneDisplacement      ! zero plane displacement (m)
  case('scalarRiBulkCanopy'             ); get_ixDiag = iLookDIAG%scalarRiBulkCanopy               ! bulk Richardson number for the canopy (-)
  case('scalarRiBulkGround'             ); get_ixDiag = iLookDIAG%scalarRiBulkGround               ! bulk Richardson number for the ground surface (-)
  case('scalarCanopyStabilityCorrection'); get_ixDiag = iLookDIAG%scalarCanopyStabilityCorrection  ! stability correction for the canopy (-)
  case('scalarGroundStabilityCorrection'); get_ixDiag = iLookDIAG%scalarGroundStabilityCorrection  ! stability correction for the ground surface (-)
  ! evapotranspiration
  case('scalarIntercellularCO2Sunlit'   ); get_ixDiag = iLookDIAG%scalarIntercellularCO2Sunlit     ! carbon dioxide partial pressure of leaf interior (sunlit leaves) (Pa)
  case('scalarIntercellularCO2Shaded'   ); get_ixDiag = iLookDIAG%scalarIntercellularCO2Shaded     ! carbon dioxide partial pressure of leaf interior (shaded leaves) (Pa)
  case('scalarTranspireLim'             ); get_ixDiag = iLookDIAG%scalarTranspireLim               ! aggregate soil moisture and aquifer storage limit on transpiration (-)
  case('scalarTranspireLimAqfr'         ); get_ixDiag = iLookDIAG%scalarTranspireLimAqfr           ! aquifer storage limit on transpiration (-)
  case('scalarFoliageNitrogenFactor'    ); get_ixDiag = iLookDIAG%scalarFoliageNitrogenFactor      ! foliage nitrogen concentration, 1=saturated (-)
  case('scalarSoilRelHumidity'          ); get_ixDiag = iLookDIAG%scalarSoilRelHumidity            ! relative humidity in the soil pores in the upper-most soil layer (-)
  case('mLayerTranspireLim'             ); get_ixDiag = iLookDIAG%mLayerTranspireLim               ! moisture avail factor limiting transpiration in each layer (-)
  case('mLayerRootDensity'              ); get_ixDiag = iLookDIAG%mLayerRootDensity                ! fraction of roots in each soil layer (-)
  case('scalarAquiferRootFrac'          ); get_ixDiag = iLookDIAG%scalarAquiferRootFrac            ! fraction of roots below the soil profile (-)
  ! canopy hydrology
  case('scalarFracLiqVeg'               ); get_ixDiag = iLookDIAG%scalarFracLiqVeg                 ! fraction of liquid water on vegetation (-)
  case('scalarCanopyWetFraction'        ); get_ixDiag = iLookDIAG%scalarCanopyWetFraction          ! fraction of canopy that is wet
  ! snow hydrology
  case('scalarSnowAge'                  ); get_ixDiag = iLookDIAG%scalarSnowAge                    ! non-dimensional snow age (-)
  case('scalarGroundSnowFraction'       ); get_ixDiag = iLookDIAG%scalarGroundSnowFraction         ! fraction of ground that is covered with snow (-)
  case('spectralSnowAlbedoDirect'       ); get_ixDiag = iLookDIAG%spectralSnowAlbedoDirect         ! direct snow albedo for individual spectral bands (-)
  case('mLayerFracLiqSnow'              ); get_ixDiag = iLookDIAG%mLayerFracLiqSnow                ! fraction of liquid water in each snow layer (-)
  case('mLayerThetaResid'               ); get_ixDiag = iLookDIAG%mLayerThetaResid                 ! residual volumetric water content in each snow layer (-)
  case('mLayerPoreSpace'                ); get_ixDiag = iLookDIAG%mLayerPoreSpace                  ! total pore space in each snow layer (-)
  case('mLayerMeltFreeze'               ); get_ixDiag = iLookDIAG%mLayerMeltFreeze                 ! ice content change from melt/freeze in each layer (kg m-3)
  ! soil hydrology
  case('scalarInfilArea'                ); get_ixDiag = iLookDIAG%scalarInfilArea                  ! fraction of unfrozen area where water can infiltrate (-)
  case('scalarFrozenArea'               ); get_ixDiag = iLookDIAG%scalarFrozenArea                 ! fraction of area that is considered impermeable due to soil ice (-)
  case('scalarSoilControl'              ); get_ixDiag = iLookDIAG%scalarSoilControl                ! soil control on infiltration: 1=controlling; 0=not (-)
  case('mLayerVolFracAir'               ); get_ixDiag = iLookDIAG%mLayerVolFracAir                 ! volumetric fraction of air in each layer (-)
  case('mLayerTcrit'                    ); get_ixDiag = iLookDIAG%mLayerTcrit                      ! critical soil temperature above which all water is unfrozen (K)
  case('mLayerCompress'                 ); get_ixDiag = iLookDIAG%mLayerCompress                   ! change in volumetric water content due to compression of soil (s-1)
  case('scalarSoilCompress'             ); get_ixDiag = iLookDIAG%scalarSoilCompress               ! change in total soil storage due to compression of the soil matrix (kg m-2 s-1)
  case('mLayerMatricHeadLiq'            ); get_ixDiag = iLookDIAG%mLayerMatricHeadLiq              ! matric potential of liquid water (m)
  ! mass balance check
  case('scalarTotalSoilLiq'             ); get_ixDiag = iLookDIAG%scalarTotalSoilLiq               ! total mass of liquid water in the soil (kg m-2)
  case('scalarTotalSoilIce'             ); get_ixDiag = iLookDIAG%scalarTotalSoilIce               ! total mass of ice in the soil (kg m-2)
  case('scalarTotalSoilWat'             ); get_ixDiag = iLookDIAG%scalarTotalSoilWat               ! total mass of water in the soil (kg m-2)
  ! variable shortcuts
  case('scalarVGn_m'                    ); get_ixDiag = iLookDIAG%scalarVGn_m                      ! van Genuchten "m" parameter (-)
  case('scalarKappa'                    ); get_ixDiag = iLookDIAG%scalarKappa                      ! constant in the freezing curve function (m K-1)
  case('scalarVolLatHt_fus'             ); get_ixDiag = iLookDIAG%scalarVolLatHt_fus               ! volumetric latent heat of fusion     (J m-3)
  ! timing information
  case('numFluxCalls'                   ); get_ixDiag = iLookDIAG%numFluxCalls                     ! number of flux calls (-)
  case('wallClockTime'                  ); get_ixDiag = iLookDIAG%wallClockTime                    ! wall clock time (s)
  case('meanStepSize'                   ); get_ixDiag = iLookDIAG%meanStepSize                     ! mean time step size (s) over data window
  ! balances
  case('balanceCasNrg'                  ); get_ixDiag = iLookDIAG%balanceCasNrg                    ! balance of energy in the canopy air space (W m-3)
  case('balanceVegNrg'                  ); get_ixDiag = iLookDIAG%balanceVegNrg                    ! balance of energy in the vegetation canopy (W m-3)
  case('balanceLayerNrg'                ); get_ixDiag = iLookDIAG%balanceLayerNrg                  ! balance of energy in each snow+soil layer (W m-3)
  case('balanceSnowNrg'                 ); get_ixDiag = iLookDIAG%balanceSnowNrg                   ! balance of energy in the snow (W m-3)
  case('balanceSoilNrg'                 ); get_ixDiag = iLookDIAG%balanceSoilNrg                   ! balance of energy in the soil (W m-3)
  case('balanceVegMass'                 ); get_ixDiag = iLookDIAG%balanceVegMass                   ! balance of water in the vegetation canopy (kg m-2 s-1)
  case('balanceLayerMass'               ); get_ixDiag = iLookDIAG%balanceLayerMass                 ! balance of water in each snow+soil layer (kg m-2 s-1)
  case('balanceSnowMass'                ); get_ixDiag = iLookDIAG%balanceSnowMass                  ! balance of water in the snow (kg m-2 s-1)
  case('balanceSoilMass'                ); get_ixDiag = iLookDIAG%balanceSoilMass                  ! balance of water in the soil (kg m-2 s-1)
  case('balanceAqMass'                  ); get_ixDiag = iLookDIAG%balanceAqMass                    ! balance of water in the aquifer (kg m-2 s-1)
  ! sundials integrator stats
  case('numSteps'                       ); get_ixDiag = iLookDIAG%numSteps
  case('numResEvals'                    ); get_ixDiag = iLookDIAG%numResEvals
  case('numLinSolvSetups'               ); get_ixDiag = iLookDIAG%numLinSolvSetups
  case('numErrTestFails'                ); get_ixDiag = iLookDIAG%numErrTestFails
  case('kLast'                          ); get_ixDiag = iLookDIAG%kLast
  case('kCur'                           ); get_ixDiag = iLookDIAG%kCur
  case('hInitUsed'                      ); get_ixDiag = iLookDIAG%hInitUsed
  case('hLast'                          ); get_ixDiag = iLookDIAG%hLast
  case('hCur'                           ); get_ixDiag = iLookDIAG%hCur
  case('tCur'                           ); get_ixDiag = iLookDIAG%tCur
  ! get to here if cannot find the variable
  case default
   get_ixDiag = integerMissing
 end select
 end function get_ixDiag


 ! *******************************************************************************************************************
 ! public function get_ixDiag: get the index of the named variables for the fluxes
 ! *******************************************************************************************************************
 function get_ixFlux(varName)
 USE var_lookup,only:iLookFLUX                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixFlux              ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  ! net energy and mass fluxes for the vegetation domain
  case('scalarCanairNetNrgFlux'         ); get_ixFlux = iLookFLUX%scalarCanairNetNrgFlux           ! net energy flux for the canopy air space (W m-2)
  case('scalarCanopyNetNrgFlux'         ); get_ixFlux = iLookFLUX%scalarCanopyNetNrgFlux           ! net energy flux for the vegetation canopy (W m-2)
  case('scalarGroundNetNrgFlux'         ); get_ixFlux = iLookFLUX%scalarGroundNetNrgFlux           ! net energy flux for the ground surface (W m-2)
  case('scalarCanopyNetLiqFlux'         ); get_ixFlux = iLookFLUX%scalarCanopyNetLiqFlux           ! net liquid water flux for the vegetation canopy (kg m-2 s-1)
  ! forcing
  case('scalarRainfall'                 ); get_ixFlux = iLookFLUX%scalarRainfall                   ! computed rainfall rate (kg m-2 s-1)
  case('scalarSnowfall'                 ); get_ixFlux = iLookFLUX%scalarSnowfall                   ! computed snowfall rate (kg m-2 s-1)
  ! shortwave radiation
  case('spectralIncomingDirect'         ); get_ixFlux = iLookFLUX%spectralIncomingDirect           ! incoming direct solar radiation in each wave band (W m-2)
  case('spectralIncomingDiffuse'        ); get_ixFlux = iLookFLUX%spectralIncomingDiffuse          ! incoming diffuse solar radiation in each wave band (W m-2)
  case('scalarCanopySunlitPAR'          ); get_ixFlux = iLookFLUX%scalarCanopySunlitPAR            ! average absorbed par for sunlit leaves (w m-2)
  case('scalarCanopyShadedPAR'          ); get_ixFlux = iLookFLUX%scalarCanopyShadedPAR            ! average absorbed par for shaded leaves (w m-2)
  case('spectralBelowCanopyDirect'      ); get_ixFlux = iLookFLUX%spectralBelowCanopyDirect        ! downward direct flux below veg layer for each spectral band  W m-2)
  case('spectralBelowCanopyDiffuse'     ); get_ixFlux = iLookFLUX%spectralBelowCanopyDiffuse       ! downward diffuse flux below veg layer for each spectral band (W m-2)
  case('scalarBelowCanopySolar'         ); get_ixFlux = iLookFLUX%scalarBelowCanopySolar           ! solar radiation transmitted below the canopy (W m-2)
  case('scalarCanopyAbsorbedSolar'      ); get_ixFlux = iLookFLUX%scalarCanopyAbsorbedSolar        ! solar radiation absorbed by canopy (W m-2)
  case('scalarGroundAbsorbedSolar'      ); get_ixFlux = iLookFLUX%scalarGroundAbsorbedSolar        ! solar radiation absorbed by ground (W m-2)
  ! longwave radiation
  case('scalarLWRadCanopy'              ); get_ixFlux = iLookFLUX%scalarLWRadCanopy                ! longwave radiation emitted from the canopy (W m-2)
  case('scalarLWRadGround'              ); get_ixFlux = iLookFLUX%scalarLWRadGround                ! longwave radiation emitted at the ground surface  (W m-2)
  case('scalarLWRadUbound2Canopy'       ); get_ixFlux = iLookFLUX%scalarLWRadUbound2Canopy         ! downward atmospheric longwave radiation absorbed by the canopy (W m-2)
  case('scalarLWRadUbound2Ground'       ); get_ixFlux = iLookFLUX%scalarLWRadUbound2Ground         ! downward atmospheric longwave radiation absorbed by the ground (W m-2)
  case('scalarLWRadUbound2Ubound'       ); get_ixFlux = iLookFLUX%scalarLWRadUbound2Ubound         ! atmospheric radiation refl by ground + lost thru upper boundary (W m-2)
  case('scalarLWRadCanopy2Ubound'       ); get_ixFlux = iLookFLUX%scalarLWRadCanopy2Ubound         ! longwave radiation emitted from canopy lost thru upper boundary (W m-2)
  case('scalarLWRadCanopy2Ground'       ); get_ixFlux = iLookFLUX%scalarLWRadCanopy2Ground         ! longwave radiation emitted from canopy absorbed by the ground (W m-2)
  case('scalarLWRadCanopy2Canopy'       ); get_ixFlux = iLookFLUX%scalarLWRadCanopy2Canopy         ! canopy longwave reflected from ground and absorbed by the canopy (W m-2)
  case('scalarLWRadGround2Ubound'       ); get_ixFlux = iLookFLUX%scalarLWRadGround2Ubound         ! longwave radiation emitted from ground lost thru upper boundary (W m-2)
  case('scalarLWRadGround2Canopy'       ); get_ixFlux = iLookFLUX%scalarLWRadGround2Canopy         ! longwave radiation emitted from ground and absorbed by the canopy (W m-2)
  case('scalarLWNetCanopy'              ); get_ixFlux = iLookFLUX%scalarLWNetCanopy                ! net longwave radiation at the canopy (W m-2)
  case('scalarLWNetGround'              ); get_ixFlux = iLookFLUX%scalarLWNetGround                ! net longwave radiation at the ground surface (W m-2)
  case('scalarLWNetUbound'              ); get_ixFlux = iLookFLUX%scalarLWNetUbound                ! net longwave radiation at the upper atmospheric boundary (W m-2)
  ! turbulent heat transfer
  case('scalarEddyDiffusCanopyTop'      ); get_ixFlux = iLookFLUX%scalarEddyDiffusCanopyTop        ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
  case('scalarFrictionVelocity'         ); get_ixFlux = iLookFLUX%scalarFrictionVelocity           ! friction velocity - canopy momentum sink (m s-1)
  case('scalarWindspdCanopyTop'         ); get_ixFlux = iLookFLUX%scalarWindspdCanopyTop           ! windspeed at the top of the canopy (m s-1)
  case('scalarWindspdCanopyBottom'      ); get_ixFlux = iLookFLUX%scalarWindspdCanopyBottom        ! windspeed at the height of the bottom of the canopy (m s-1)
  case('scalarGroundResistance'         ); get_ixFlux = iLookFLUX%scalarGroundResistance           ! below canopy aerodynamic resistance (s m-1)
  case('scalarCanopyResistance'         ); get_ixFlux = iLookFLUX%scalarCanopyResistance           ! above canopy aerodynamic resistance (s m-1)
  case('scalarLeafResistance'           ); get_ixFlux = iLookFLUX%scalarLeafResistance             ! mean leaf boundary layer resistance per unit leaf area (s m-1)
  case('scalarSoilResistance'           ); get_ixFlux = iLookFLUX%scalarSoilResistance             ! soil surface resistance (s m-1)
  case('scalarSenHeatTotal'             ); get_ixFlux = iLookFLUX%scalarSenHeatTotal               ! sensible heat from the canopy air space to the atmosphere (W m-2)
  case('scalarSenHeatCanopy'            ); get_ixFlux = iLookFLUX%scalarSenHeatCanopy              ! sensible heat from the canopy to the canopy air space (W m-2)
  case('scalarSenHeatGround'            ); get_ixFlux = iLookFLUX%scalarSenHeatGround              ! sensible heat from the ground (below canopy or non-vegetated) (W m-2)
  case('scalarLatHeatTotal'             ); get_ixFlux = iLookFLUX%scalarLatHeatTotal               ! latent heat from the canopy air space to the atmosphere (W m-2)
  case('scalarLatHeatCanopyEvap'        ); get_ixFlux = iLookFLUX%scalarLatHeatCanopyEvap          ! evaporation latent heat from the canopy to the canopy air space (W m-2)
  case('scalarLatHeatCanopyTrans'       ); get_ixFlux = iLookFLUX%scalarLatHeatCanopyTrans         ! transpiration latent heat from the canopy to the canopy air space (W m-2)
  case('scalarLatHeatGround'            ); get_ixFlux = iLookFLUX%scalarLatHeatGround              ! latent heat from the ground (below canopy or non-vegetated) (W m-2)
  case('scalarCanopyAdvectiveHeatFlux'  ); get_ixFlux = iLookFLUX%scalarCanopyAdvectiveHeatFlux    ! heat advected to the canopy surface with rain + snow (W m-2)
  case('scalarGroundAdvectiveHeatFlux'  ); get_ixFlux = iLookFLUX%scalarGroundAdvectiveHeatFlux    ! heat advected to the ground surface with throughfall and unloading/drainage (W m-2)
  case('scalarCanopySublimation'        ); get_ixFlux = iLookFLUX%scalarCanopySublimation          ! canopy sublimation/frost (kg m-2 s-1)
  case('scalarSnowSublimation'          ); get_ixFlux = iLookFLUX%scalarSnowSublimation            ! snow sublimation/frost (below canopy or non-vegetated) (kg m-2 s-1)
  ! liquid water fluxes associated with evapotranspiration
  case('scalarStomResistSunlit'         ); get_ixFlux = iLookFLUX%scalarStomResistSunlit           ! stomatal resistance for sunlit leaves (s m-1)
  case('scalarStomResistShaded'         ); get_ixFlux = iLookFLUX%scalarStomResistShaded           ! stomatal resistance for shaded leaves (s m-1)
  case('scalarPhotosynthesisSunlit'     ); get_ixFlux = iLookFLUX%scalarPhotosynthesisSunlit       ! sunlit photosynthesis (umolco2 m-2 s-1)
  case('scalarPhotosynthesisShaded'     ); get_ixFlux = iLookFLUX%scalarPhotosynthesisShaded       ! shaded photosynthesis (umolco2 m-2 s-1)
  case('scalarCanopyTranspiration'      ); get_ixFlux = iLookFLUX%scalarCanopyTranspiration        ! canopy transpiration (kg m-2 s-1)
  case('scalarCanopyEvaporation'        ); get_ixFlux = iLookFLUX%scalarCanopyEvaporation          ! canopy evaporation/condensation (kg m-2 s-1)
  case('scalarGroundEvaporation'        ); get_ixFlux = iLookFLUX%scalarGroundEvaporation          ! ground evaporation/condensation (below canopy or non-vegetated) (kg m-2 s-1)
  case('mLayerTranspire'                ); get_ixFlux = iLookFLUX%mLayerTranspire                  ! transpiration loss from each soil layer (kg m-2 s-1)
  ! liquid and solid water fluxes through the canopy
  case('scalarThroughfallSnow'          ); get_ixFlux = iLookFLUX%scalarThroughfallSnow            ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
  case('scalarThroughfallRain'          ); get_ixFlux = iLookFLUX%scalarThroughfallRain            ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
  case('scalarCanopySnowUnloading'      ); get_ixFlux = iLookFLUX%scalarCanopySnowUnloading        ! unloading of snow from the vegetion canopy (kg m-2 s-1)
  case('scalarCanopyLiqDrainage'        ); get_ixFlux = iLookFLUX%scalarCanopyLiqDrainage          ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
  case('scalarCanopyMeltFreeze'         ); get_ixFlux = iLookFLUX%scalarCanopyMeltFreeze           ! melt/freeze of water stored in the canopy (kg m-2 s-1)
  ! energy fluxes and for the snow and soil domains
  case('iLayerConductiveFlux'           ); get_ixFlux = iLookFLUX%iLayerConductiveFlux             ! conductive energy flux at layer interfaces at end of time step (W m-2)
  case('iLayerAdvectiveFlux'            ); get_ixFlux = iLookFLUX%iLayerAdvectiveFlux              ! advective energy flux at layer interfaces at end of time step (W m-2)
  case('iLayerNrgFlux'                  ); get_ixFlux = iLookFLUX%iLayerNrgFlux                    ! energy flux at layer interfaces at the end of the time step (W m-2)
  case('mLayerNrgFlux'                  ); get_ixFlux = iLookFLUX%mLayerNrgFlux                    ! net energy flux for each layer in the snow+soil domain (J m-3 s-1)
  ! liquid water fluxes for the snow domain
  case('scalarSnowDrainage'             ); get_ixFlux = iLookFLUX%scalarSnowDrainage               ! drainage from the bottom of the snow profile (m s-1)
  case('iLayerLiqFluxSnow'              ); get_ixFlux = iLookFLUX%iLayerLiqFluxSnow                ! liquid flux at snow layer interfaces at the end of the time step (m s-1)
  case('mLayerLiqFluxSnow'              ); get_ixFlux = iLookFLUX%mLayerLiqFluxSnow                ! net liquid water flux for each snow layer (s-1)
  ! liquid water fluxes for the soil domain
  case('scalarRainPlusMelt'             ); get_ixFlux = iLookFLUX%scalarRainPlusMelt               ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  case('scalarMaxInfilRate'             ); get_ixFlux = iLookFLUX%scalarMaxInfilRate               ! maximum infiltration rate (m s-1)
  case('scalarInfiltration'             ); get_ixFlux = iLookFLUX%scalarInfiltration               ! infiltration of water into the soil profile (m s-1)
  case('scalarExfiltration'             ); get_ixFlux = iLookFLUX%scalarExfiltration               ! exfiltration of water from the top of the soil profile (m s-1)
  case('scalarSurfaceRunoff'            ); get_ixFlux = iLookFLUX%scalarSurfaceRunoff              ! surface runoff (m s-1)
  case('mLayerSatHydCondMP'             ); get_ixFlux = iLookFLUX%mLayerSatHydCondMP               ! saturated hydraulic conductivity of macropores in each layer (m s-1)
  case('mLayerSatHydCond'               ); get_ixFlux = iLookFLUX%mLayerSatHydCond                 ! saturated hydraulic conductivity in each layer (m s-1)
  case('iLayerSatHydCond'               ); get_ixFlux = iLookFLUX%iLayerSatHydCond                 ! saturated hydraulic conductivity in each layer interface (m s-1)
  case('mLayerHydCond'                  ); get_ixFlux = iLookFLUX%mLayerHydCond                    ! hydraulic conductivity in each layer (m s-1)
  case('iLayerLiqFluxSoil'              ); get_ixFlux = iLookFLUX%iLayerLiqFluxSoil                ! liquid flux at soil layer interfaces at the end of the time step (m s-1)
  case('mLayerLiqFluxSoil'              ); get_ixFlux = iLookFLUX%mLayerLiqFluxSoil                ! net liquid water flux for each soil layer (s-1)
  case('mLayerBaseflow'                 ); get_ixFlux = iLookFLUX%mLayerBaseflow                   ! baseflow from each soil layer (m s-1)
  case('mLayerColumnInflow'             ); get_ixFlux = iLookFLUX%mLayerColumnInflow               ! total inflow to each layer in a given soil column (m3 s-1)
  case('mLayerColumnOutflow'            ); get_ixFlux = iLookFLUX%mLayerColumnOutflow              ! total outflow from each layer in a given soil column (m3 s-1)
  case('scalarSoilBaseflow'             ); get_ixFlux = iLookFLUX%scalarSoilBaseflow               ! total baseflow from throughout the soil profile (m s-1)
  case('scalarSoilDrainage'             ); get_ixFlux = iLookFLUX%scalarSoilDrainage               ! drainage from the bottom of the soil profile (m s-1)
  case('scalarAquiferRecharge'          ); get_ixFlux = iLookFLUX%scalarAquiferRecharge            ! recharge to the aquifer (m s-1)
  case('scalarAquiferTranspire'         ); get_ixFlux = iLookFLUX%scalarAquiferTranspire           ! transpiration from the aquifer (m s-1)
  case('scalarAquiferBaseflow'          ); get_ixFlux = iLookFLUX%scalarAquiferBaseflow            ! baseflow from the aquifer (m s-1)
  ! derived variables
  case('scalarTotalET'                  ); get_ixFlux = iLookFLUX%scalarTotalET                    ! total ET (kg m-2 s-1)
  case('scalarTotalRunoff'              ); get_ixFlux = iLookFLUX%scalarTotalRunoff                ! total runoff (m s-1)
  case('scalarNetRadiation'             ); get_ixFlux = iLookFLUX%scalarNetRadiation               ! net radiation (W m-2)
  ! return missing if variable not found
  case default
   get_ixFlux = integerMissing
 end select
 end function get_ixFlux


 ! *******************************************************************************************************************
 ! public function get_ixDeriv: get the index of the named variables for the model derivatives
 ! *******************************************************************************************************************
 function get_ixDeriv(varName)
 USE var_lookup,only:iLookDERIV                      ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! parameter name
 integer(i4b)             :: get_ixDeriv             ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  ! derivatives in net vegetation energy fluxes w.r.t. relevant state variables
  case('dCanairNetFlux_dCanairTemp'     ); get_ixDeriv = iLookDERIV%dCanairNetFlux_dCanairTemp     ! derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
  case('dCanairNetFlux_dCanopyTemp'     ); get_ixDeriv = iLookDERIV%dCanairNetFlux_dCanopyTemp     ! derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
  case('dCanairNetFlux_dGroundTemp'     ); get_ixDeriv = iLookDERIV%dCanairNetFlux_dGroundTemp     ! derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
  case('dCanopyNetFlux_dCanairTemp'     ); get_ixDeriv = iLookDERIV%dCanopyNetFlux_dCanairTemp     ! derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
  case('dCanopyNetFlux_dCanopyTemp'     ); get_ixDeriv = iLookDERIV%dCanopyNetFlux_dCanopyTemp     ! derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
  case('dCanopyNetFlux_dGroundTemp'     ); get_ixDeriv = iLookDERIV%dCanopyNetFlux_dGroundTemp     ! derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
  case('dCanopyNetFlux_dCanWat'         ); get_ixDeriv = iLookDERIV%dCanopyNetFlux_dCanWat         ! derivative in net canopy fluxes w.r.t. canopy total water content (J kg-1 s-1)
  case('dGroundNetFlux_dCanairTemp'     ); get_ixDeriv = iLookDERIV%dGroundNetFlux_dCanairTemp     ! derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
  case('dGroundNetFlux_dCanopyTemp'     ); get_ixDeriv = iLookDERIV%dGroundNetFlux_dCanopyTemp     ! derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
  case('dGroundNetFlux_dGroundTemp'     ); get_ixDeriv = iLookDERIV%dGroundNetFlux_dGroundTemp     ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
  case('dGroundNetFlux_dCanWat'         ); get_ixDeriv = iLookDERIV%dGroundNetFlux_dCanWat         ! derivative in net ground fluxes w.r.t. canopy total water content (J kg-1 s-1)
  ! derivatives in evaporative fluxes w.r.t. relevant state variables
  case('dCanopyEvaporation_dTCanair'    ); get_ixDeriv = iLookDERIV%dCanopyEvaporation_dTCanair    ! derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
  case('dCanopyEvaporation_dTCanopy'    ); get_ixDeriv = iLookDERIV%dCanopyEvaporation_dTCanopy    ! derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
  case('dCanopyEvaporation_dTGround'    ); get_ixDeriv = iLookDERIV%dCanopyEvaporation_dTGround    ! derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
  case('dCanopyEvaporation_dCanWat'     ); get_ixDeriv = iLookDERIV%dCanopyEvaporation_dCanWat     ! derivative in canopy evaporation w.r.t. canopy total water content (s-1)
  case('dGroundEvaporation_dTCanair'    ); get_ixDeriv = iLookDERIV%dGroundEvaporation_dTCanair    ! derivative in ground evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
  case('dGroundEvaporation_dTCanopy'    ); get_ixDeriv = iLookDERIV%dGroundEvaporation_dTCanopy    ! derivative in ground evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
  case('dGroundEvaporation_dTGround'    ); get_ixDeriv = iLookDERIV%dGroundEvaporation_dTGround    ! derivative in ground evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
  case('dGroundEvaporation_dCanWat'     ); get_ixDeriv = iLookDERIV%dGroundEvaporation_dCanWat     ! derivative in ground evaporation w.r.t. canopy total water content (s-1)
  ! derivatives in transpiration
  case('dCanopyTrans_dTCanair'          ); get_ixDeriv = iLookDERIV%dCanopyTrans_dTCanair          ! derivative in canopy transpiration w.r.t. canopy air temperature (kg m-2 s-1 K-1)
  case('dCanopyTrans_dTCanopy'          ); get_ixDeriv = iLookDERIV%dCanopyTrans_dTCanopy          ! derivative in canopy transpiration w.r.t. canopy temperature (kg m-2 s-1 K-1)
  case('dCanopyTrans_dTGround'          ); get_ixDeriv = iLookDERIV%dCanopyTrans_dTGround          ! derivative in canopy transpiration w.r.t. ground temperature (kg m-2 s-1 K-1)
  case('dCanopyTrans_dCanWat'           ); get_ixDeriv = iLookDERIV%dCanopyTrans_dCanWat           ! derivative in canopy transpiration w.r.t. canopy total water content (s-1)
  ! derivatives in canopy water w.r.t canopy temperature
  case('dTheta_dTkCanopy'               ); get_ixDeriv = iLookDERIV%dTheta_dTkCanopy               ! derivative of volumetric liquid water content w.r.t. temperature (K-1)
  case('d2Theta_dTkCanopy2'             ); get_ixDeriv = iLookDERIV%d2Theta_dTkCanopy2             ! second derivative of volumetric liquid water content w.r.t. temperature
  case('dCanLiq_dTcanopy'               ); get_ixDeriv = iLookDERIV%dCanLiq_dTcanopy               ! derivative of canopy liquid storage w.r.t. temperature (kg m-2 K-1)
  case('dFracLiqVeg_dTkCanopy'          ); get_ixDeriv = iLookDERIV%dFracLiqVeg_dTkCanopy          ! derivative in fraction of (throughfall + drainage)  w.r.t. temperature
  ! derivatives in canopy liquid fluxes w.r.t. canopy water
  case('scalarCanopyLiqDeriv'           ); get_ixDeriv = iLookDERIV%scalarCanopyLiqDeriv           ! derivative in (throughfall + canopy drainage) w.r.t. canopy liquid water (s-1)
  case('scalarThroughfallRainDeriv'     ); get_ixDeriv = iLookDERIV%scalarThroughfallRainDeriv     ! derivative in throughfall w.r.t. canopy liquid water (s-1)
  case('scalarCanopyLiqDrainageDeriv'   ); get_ixDeriv = iLookDERIV%scalarCanopyLiqDrainageDeriv   ! derivative in canopy drainage w.r.t. canopy liquid water (s-1)
  ! energy derivatives that might be treated as constant if heat capacity and thermal conductivity not updated
  case('dVolHtCapBulk_dPsi0'            ); get_ixDeriv = iLookDERIV%dVolHtCapBulk_dPsi0            ! derivative in bulk heat capacity w.r.t. matric potential
  case('dVolHtCapBulk_dTheta'           ); get_ixDeriv = iLookDERIV%dVolHtCapBulk_dTheta           ! derivative in bulk heat capacity w.r.t. volumetric water content
  case('dVolHtCapBulk_dCanWat'          ); get_ixDeriv = iLookDERIV%dVolHtCapBulk_dCanWat          ! derivative in bulk heat capacity w.r.t. canopy volumetric water content
  case('dVolHtCapBulk_dTk'              ); get_ixDeriv = iLookDERIV%dVolHtCapBulk_dTk              ! derivative in bulk heat capacity w.r.t. temperature
  case('dVolHtCapBulk_dTkCanopy'        ); get_ixDeriv = iLookDERIV%dVolHtCapBulk_dTkCanopy        ! derivative in bulk heat capacity w.r.t. canopy temperature
  case('dThermalC_dTempAbove'           ); get_ixDeriv = iLookDERIV%dThermalC_dTempAbove           ! derivative in the thermal conductivity w.r.t. energy state in the layer above
  case('dThermalC_dTempBelow'           ); get_ixDeriv = iLookDERIV%dThermalC_dTempBelow           ! derivative in the thermal conductivity w.r.t. energy state in the layer above
  case('dThermalC_dWatAbove'            ); get_ixDeriv = iLookDERIV%dThermalC_dWatAbove            ! derivative in the thermal conductivity w.r.t. water state in the layer above
  case('dThermalC_dWatBelow'            ); get_ixDeriv = iLookDERIV%dThermalC_dWatBelow            ! derivative in the thermal conductivity w.r.t. water state in the layer above
  case('dNrgFlux_dTempAbove'            ); get_ixDeriv = iLookDERIV%dNrgFlux_dTempAbove            ! derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
  case('dNrgFlux_dTempBelow'            ); get_ixDeriv = iLookDERIV%dNrgFlux_dTempBelow            ! derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
  ! energy derivatives that might be treated as constant if Cm not updated
  case('dCm_dTk'                        ); get_ixDeriv = iLookDERIV%dCm_dTk                        ! derivative in Cm w.r.t. temperature (J kg K-2)
  case('dCm_dTkCanopy'                  ); get_ixDeriv = iLookDERIV%dCm_dTkCanopy                  ! derivative in Cm w.r.t. canopy temperature (J kg K-2)
  ! derivatives in energy fluxes at the interface of snow+soil layers w.r.t. water state in layers above and below
  case('dNrgFlux_dWatAbove'             ); get_ixDeriv = iLookDERIV%dNrgFlux_dWatAbove             ! derivatives in the flux w.r.t. water state temperature in the layer above
  case('dNrgFlux_dWatBelow'             ); get_ixDeriv = iLookDERIV%dNrgFlux_dWatBelow             ! derivatives in the flux w.r.t. water state in the layer below
  ! derivative in liquid water fluxes at the interface of snow layers w.r.t. volumetric liquid water content in the layer above
  case('iLayerLiqFluxSnowDeriv'         ); get_ixDeriv = iLookDERIV%iLayerLiqFluxSnowDeriv         ! derivative in vertical liquid water flux at layer interfaces (m s-1)
  ! derivative in liquid water fluxes for the soil domain w.r.t hydrology state variables
  case('dVolTot_dPsi0'                  ); get_ixDeriv = iLookDERIV%dVolTot_dPsi0                  ! derivative in total water content w.r.t. total water matric potential (m-1)
  case('d2VolTot_dPsi02'                ); get_ixDeriv = iLookDERIV%d2VolTot_dPsi02                ! second derivative in total water content w.r.t. total water matric potential
  case('dq_dHydStateAbove'              ); get_ixDeriv = iLookDERIV%dq_dHydStateAbove              ! change in the flux in layer interfaces w.r.t. state variables in the layer above
  case('dq_dHydStateBelow'              ); get_ixDeriv = iLookDERIV%dq_dHydStateBelow              ! change in the flux in layer interfaces w.r.t. state variables in the layer below
  case('dq_dHydStateLayerSurfVec'       ); get_ixDeriv = iLookDERIV%dq_dHydStateLayerSurfVec       ! change in the flux in soil surface interface w.r.t. state variables in layer above and below
  case('mLayerdTheta_dPsi'              ); get_ixDeriv = iLookDERIV%mLayerdTheta_dPsi              ! derivative in the soil water characteristic w.r.t. psi (m-1)
  case('mLayerdPsi_dTheta'              ); get_ixDeriv = iLookDERIV%mLayerdPsi_dTheta              ! derivative in the soil water characteristic w.r.t. theta (m)
  case('dCompress_dPsi'                 ); get_ixDeriv = iLookDERIV%dCompress_dPsi                 ! derivative in compressibility w.r.t matric head (m-1)
  ! derivative in baseflow flux w.r.t. aquifer storage
  case('dBaseflow_dAquifer'             ); get_ixDeriv = iLookDERIV%dBaseflow_dAquifer             ! derivative in baseflow flux w.r.t. aquifer storage (s-1)
  ! derivative in liquid water fluxes for the soil domain w.r.t energy state variables
  case('dq_dNrgStateAbove'              ); get_ixDeriv = iLookDERIV%dq_dNrgStateAbove              ! change in the flux in layer interfaces w.r.t. state variables in the layer above
  case('dq_dNrgStateBelow'              ); get_ixDeriv = iLookDERIV%dq_dNrgStateBelow              ! change in the flux in layer interfaces w.r.t. state variables in the layer below
  case('dq_dNrgStateLayerSurfVec'       ); get_ixDeriv = iLookDERIV%dq_dNrgStateLayerSurfVec       ! change in the flux in soil surface interface w.r.t. state variables in layer above and below
  case('dPsiLiq_dTemp'                  ); get_ixDeriv = iLookDERIV%dPsiLiq_dTemp                  ! derivative in the liquid water matric potential w.r.t. temperature (m K-1)
  case('dPsiLiq_dPsi0'                  ); get_ixDeriv = iLookDERIV%dPsiLiq_dPsi0                  ! derivative in liquid matric potential w.r.t. total  matric potential (-)
 ! derivatives in soil transpiration w.r.t. canopy state variables
  case('mLayerdTrans_dTCanair'          ); get_ixDeriv = iLookDERIV%mLayerdTrans_dTCanair          ! derivatives in the soil layer transpiration flux w.r.t. canopy air temperature
  case('mLayerdTrans_dTCanopy'          ); get_ixDeriv = iLookDERIV%mLayerdTrans_dTCanopy          ! derivatives in the soil layer transpiration flux w.r.t. canopy temperature
  case('mLayerdTrans_dTGround'          ); get_ixDeriv = iLookDERIV%mLayerdTrans_dTGround          ! derivatives in the soil layer transpiration flux w.r.t. ground temperature
  case('mLayerdTrans_dCanWat'           ); get_ixDeriv = iLookDERIV%mLayerdTrans_dCanWat           ! derivatives in the soil layer transpiration flux w.r.t. canopy total water
 ! derivatives in aquifer transpiration w.r.t. canopy state variables
  case('dAquiferTrans_dTCanair'         ); get_ixDeriv = iLookDERIV%dAquiferTrans_dTCanair         ! derivative in the aquifer transpiration flux w.r.t. canopy air temperature
  case('dAquiferTrans_dTCanopy'         ); get_ixDeriv = iLookDERIV%dAquiferTrans_dTCanopy         ! derivative in the aquifer transpiration flux w.r.t. canopy temperature
  case('dAquiferTrans_dTGround'         ); get_ixDeriv = iLookDERIV%dAquiferTrans_dTGround         ! derivative in the aquifer transpiration flux w.r.t. ground temperature
  case('dAquiferTrans_dCanWat'          ); get_ixDeriv = iLookDERIV%dAquiferTrans_dCanWat          ! derivative in the aquifer transpiration flux w.r.t. canopy total water
 ! derivative in liquid water fluxes for the soil and snow domain w.r.t temperature
  case('dFracLiqSnow_dTk'               ); get_ixDeriv = iLookDERIV%dFracLiqSnow_dTk               ! derivative in fraction of liquid snow w.r.t. temperature
  case('mLayerdTheta_dTk'               ); get_ixDeriv = iLookDERIV%mLayerdTheta_dTk               ! derivative of volumetric liquid water content w.r.t. temperature (K-1)
  case('mLayerd2Theta_dTk2'             ); get_ixDeriv = iLookDERIV%mLayerd2Theta_dTk2             ! second derivative of volumetric liquid water content w.r.t. temperature
  ! derivatives in time
  case('mLayerdTemp_dt'                 ); get_ixDeriv = iLookDERIV%mLayerdTemp_dt                 ! timestep change in layer temperature
  case('scalarCanopydTemp_dt'           ); get_ixDeriv = iLookDERIV%scalarCanopydTemp_dt           ! timestep change in canopy temperature
  case('mLayerdWat_dt'                  ); get_ixDeriv = iLookDERIV%mLayerdWat_dt                  ! timestep change in layer volumetric fraction of total water
  case('scalarCanopydWat_dt'            ); get_ixDeriv = iLookDERIV%scalarCanopydWat_dt            ! timestep change in canopy water content
  ! derivatives of temperature if enthalpy is the state variable
  case('dCanairTemp_dEnthalpy'          ); get_ixDeriv = iLookDERIV%dCanairTemp_dEnthalpy          ! derivative of canopy air temperature w.r.t. enthalpy
  case('dCanopyTemp_dEnthalpy'          ); get_ixDeriv = iLookDERIV%dCanopyTemp_dEnthalpy          ! derivative of canopy temperature w.r.t. enthalpy 
  case('dTemp_dEnthalpy'                ); get_ixDeriv = iLookDERIV%dTemp_dEnthalpy                ! derivative of temperature w.r.t. enthalpy      
  case('dCanopyTemp_dCanWat'            ); get_ixDeriv = iLookDERIV%dCanopyTemp_dCanWat            ! derivative of canopy temperature w.r.t. volumetric water content  
  case('dTemp_dTheta'                   ); get_ixDeriv = iLookDERIV%dTemp_dTheta                   ! derivative of temperature w.r.t. volumetric water content         
  case('dTemp_dPsi0'                    ); get_ixDeriv = iLookDERIV%dTemp_dPsi0                    ! derivative of temperature w.r.t. total water matric potential         

  case default
   get_ixDeriv = integerMissing
 end select
 end function get_ixDeriv


 ! *******************************************************************************************************************
 ! public function get_ixIndex: get the index of the named variables for the model indices
 ! *******************************************************************************************************************
 function get_ixIndex(varName)
 USE var_lookup,only:iLookINDEX                      ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! parameter name
 integer(i4b)             :: get_ixINDEX             ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  ! number of model layers, and layer indices
  case('nSnow'                ); get_ixINDEX = iLookINDEX%nSnow                 ! number of snow layers                                                    (-)
  case('nSoil'                ); get_ixINDEX = iLookINDEX%nSoil                 ! number of soil layers                                                    (-)
  case('nLayers'              ); get_ixINDEX = iLookINDEX%nLayers               ! total number of layers                                                   (-)
  case('layerType'            ); get_ixINDEX = iLookINDEX%layerType             ! index defining type of layer (snow or soil)                              (-)
  ! number of state variables of different type
  case('nCasNrg'              ); get_ixINDEX = iLookINDEX%nCasNrg               ! number of energy state variables for the canopy air space domain         (-)
  case('nVegNrg'              ); get_ixINDEX = iLookINDEX%nVegNrg               ! number of energy state variables for the vegetation canopy               (-)
  case('nVegMass'             ); get_ixINDEX = iLookINDEX%nVegMass              ! number of hydrology states for vegetation (mass of water)                (-)
  case('nVegState'            ); get_ixINDEX = iLookINDEX%nVegState             ! number of vegetation state variables                                     (-)
  case('nNrgState'            ); get_ixINDEX = iLookINDEX%nNrgState             ! number of energy state variables                                         (-)
  case('nWatState'            ); get_ixINDEX = iLookINDEX%nWatState             ! number of "total water" states (vol. total water content)                (-)
  case('nMatState'            ); get_ixINDEX = iLookINDEX%nMatState             ! number of matric head state variables                                    (-)
  case('nMassState'           ); get_ixINDEX = iLookINDEX%nMassState            ! number of hydrology state variables (mass of water)                      (-)
  case('nState'               ); get_ixINDEX = iLookINDEX%nState                ! total number of model state variables                                    (-)
  ! number of state variables within different domains in the snow+soil system  !
  case('nSnowSoilNrg'         ); get_ixINDEX = iLookINDEX%nSnowSoilNrg          ! number of energy states in the snow+soil domain                          (-)
  case('nSnowOnlyNrg'         ); get_ixINDEX = iLookINDEX%nSnowOnlyNrg          ! number of energy states in the snow domain                               (-)
  case('nSoilOnlyNrg'         ); get_ixINDEX = iLookINDEX%nSoilOnlyNrg          ! number of energy states in the soil domain                               (-)
  case('nSnowSoilHyd'         ); get_ixINDEX = iLookINDEX%nSnowSoilHyd          ! number of hydrology states in the snow+soil domain                       (-)
  case('nSnowOnlyHyd'         ); get_ixINDEX = iLookINDEX%nSnowOnlyHyd          ! number of hydrology states in the snow domain                            (-)
  case('nSoilOnlyHyd'         ); get_ixINDEX = iLookINDEX%nSoilOnlyHyd          ! number of hydrology states in the soil domain                            (-)
  ! type of model state variables
  case('ixControlVolume'      ); get_ixINDEX = iLookINDEX%ixControlVolume       ! index of the control volume for different domains (veg, snow, soil)      (-)
  case('ixDomainType'         ); get_ixINDEX = iLookINDEX%ixDomainType          ! index of the type of domain (iname_veg, iname_snow, iname_soil)          (-)
  case('ixStateType'          ); get_ixINDEX = iLookINDEX%ixStateType           ! index of the type of every state variable (iname_nrgCanair, ...)         (-)
  case('ixHydType'            ); get_ixINDEX = iLookINDEX%ixHydType             ! index of the type of hydrology states in snow+soil domain                (-)
  ! type of model state variables (state subset)
  case('ixDomainType_subset'  ); get_ixINDEX = iLookINDEX%ixDomainType_subset   ! [state subset] id of domain for desired model state variables            (-)
  case('ixStateType_subset'   ); get_ixINDEX = iLookINDEX%ixStateType_subset    ! [state subset] type of desired model state variables                     (-)
  ! mapping between state subset and the full state vector
  case('ixMapFull2Subset'     ); get_ixINDEX = iLookINDEX%ixMapFull2Subset      ! list of indices of the state subset in the full state vector             (-)
  case('ixMapSubset2Full'     ); get_ixINDEX = iLookINDEX%ixMapSubset2Full      ! list of indices of the full state vector in the state subset             (-)
  ! indices of model specific state variables
  case('ixCasNrg'             ); get_ixINDEX = iLookINDEX%ixCasNrg              ! index of canopy air space energy state variable                          (-)
  case('ixVegNrg'             ); get_ixINDEX = iLookINDEX%ixVegNrg              ! index of canopy energy state variable                                    (-)
  case('ixVegHyd'             ); get_ixINDEX = iLookINDEX%ixVegHyd              ! index of canopy hydrology state variable (mass)                          (-)
  case('ixTopNrg'             ); get_ixINDEX = iLookINDEX%ixTopNrg              ! index of upper-most energy state in the snow+soil subdomain              (-)
  case('ixTopHyd'             ); get_ixINDEX = iLookINDEX%ixTopHyd              ! index of upper-most hydrology state in the snow+soil subdomain           (-)
  case('ixAqWat'              ); get_ixINDEX = iLookINDEX%ixAqWat               ! index of storage of water in the aquifer                                 (-)
  ! vectors of indices for specific state types
  case('ixNrgOnly'            ); get_ixINDEX = iLookINDEX%ixNrgOnly             ! indices IN THE STATE SUBSET for all energy states                        (-)
  case('ixHydOnly'            ); get_ixINDEX = iLookINDEX%ixHydOnly             ! indices IN THE STATE SUBSET for hydrology states in the snow+soil domain (-)
  case('ixMatOnly'            ); get_ixINDEX = iLookINDEX%ixMatOnly             ! indices IN THE STATE SUBSET for matric head state variables              (-)
  case('ixMassOnly'           ); get_ixINDEX = iLookINDEX%ixMassOnly            ! indices IN THE STATE SUBSET for hydrology states (mass of water)         (-)
  ! vectors of indicesfor specific state types within specific sub-domains
  case('ixSnowSoilNrg'        ); get_ixINDEX = iLookINDEX%ixSnowSoilNrg         ! indices IN THE STATE SUBSET for energy states in the snow+soil domain    (-)
  case('ixSnowOnlyNrg'        ); get_ixINDEX = iLookINDEX%ixSnowOnlyNrg         ! indices IN THE STATE SUBSET for energy states in the snow domain         (-)
  case('ixSoilOnlyNrg'        ); get_ixINDEX = iLookINDEX%ixSoilOnlyNrg         ! indices IN THE STATE SUBSET for energy states in the soil domain         (-)
  case('ixSnowSoilHyd'        ); get_ixINDEX = iLookINDEX%ixSnowSoilHyd         ! indices IN THE STATE SUBSET for hydrology states in the snow+soil domain (-)
  case('ixSnowOnlyHyd'        ); get_ixINDEX = iLookINDEX%ixSnowOnlyHyd         ! indices IN THE STATE SUBSET for hydrology states in the snow domain      (-)
  case('ixSoilOnlyHyd'        ); get_ixINDEX = iLookINDEX%ixSoilOnlyHyd         ! indices IN THE STATE SUBSET for hydrology states in the soil domain      (-)
  ! vectors of indices for specfic state types within specific sub-domains
  case('ixNrgCanair'          ); get_ixINDEX = iLookINDEX%ixNrgCanair           ! indices IN THE FULL VECTOR for energy states in canopy air space domain (-)
  case('ixNrgCanopy'          ); get_ixINDEX = iLookINDEX%ixNrgCanopy           ! indices IN THE FULL VECTOR for energy states in the canopy domain       (-)
  case('ixHydCanopy'          ); get_ixINDEX = iLookINDEX%ixHydCanopy           ! indices IN THE FULL VECTOR for hydrology states in the canopy domain    (-)
  case('ixNrgLayer'           ); get_ixINDEX = iLookINDEX%ixNrgLayer            ! indices IN THE FULL VECTOR for energy states in the snow+soil domain    (-)
  case('ixHydLayer'           ); get_ixINDEX = iLookINDEX%ixHydLayer            ! indices IN THE FULL VECTOR for hydrology states in the snow+soil domain (-)
  case('ixWatAquifer'         ); get_ixINDEX = iLookINDEX%ixWatAquifer          ! indices IN THE FULL VECTOR for storage of water in the aquifer          (-)
  ! vectors of indices for specific state types IN SPECIFIC SUB-DOMAINS
  case('ixVolFracWat'         ); get_ixINDEX = iLookINDEX%ixVolFracWat          ! indices IN THE SNOW+SOIL VECTOR for hyd states                          (-)
  case('ixMatricHead'         ); get_ixINDEX = iLookINDEX%ixMatricHead          ! indices IN THE SOIL VECTOR for hyd states                               (-)
  ! indices within state vectors
  case('ixAllState'           ); get_ixINDEX = iLookINDEX%ixAllState            ! list of indices for all model state variables                           (-)
  case('ixSoilState'          ); get_ixINDEX = iLookINDEX%ixSoilState           ! list of indices for all soil layers                                     (-)
  case('ixLayerState'         ); get_ixINDEX = iLookINDEX%ixLayerState          ! list of indices for all model layers                                    (-)
  case('ixLayerActive'        ); get_ixINDEX = iLookINDEX%ixLayerActive         ! list of indices for all active model layers                             (-)
  ! number of trials
  case('numberFluxCalc'       ); get_ixINDEX = iLookINDEX%numberFluxCalc        ! number of flux calculations                                             (-)
  case('numberStateSplit'     ); get_ixINDEX = iLookINDEX%numberStateSplit      ! number of state splitting solutions                                     (-)
  case('numberDomainSplitNrg' ); get_ixINDEX = iLookINDEX%numberDomainSplitNrg  ! number of domain splitting solutions for energy                         (-)
  case('numberDomainSplitMass'); get_ixINDEX = iLookINDEX%numberDomainSplitMass ! number of domain splitting solutions for mass                           (-)
  case('numberScalarSolutions'); get_ixINDEX = iLookINDEX%numberScalarSolutions ! number of scalar solutions                                              (-)
  ! default
  case default
   get_ixIndex = integerMissing
 end select
 end function get_ixIndex


 ! *******************************************************************************************************************
 ! public function get_ixBpar: get the index of the named variables for the basin-average variables
 ! *******************************************************************************************************************
 function get_ixBpar(varName)
 USE var_lookup,only:iLookBPAR                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! parameter name
 integer(i4b)             :: get_ixBpar              ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  ! baseflow
  case('basin__aquiferHydCond'    ); get_ixBpar = iLookBPAR%basin__aquiferHydCond     ! hydraulic conductivity of the basin aquifer (m s-1)
  case('basin__aquiferScaleFactor'); get_ixBpar = iLookBPAR%basin__aquiferScaleFactor ! scaling factor for aquifer storage in the big bucket (m)
  case('basin__aquiferBaseflowExp'); get_ixBpar = iLookBPAR%basin__aquiferBaseflowExp ! baseflow exponent for the big bucket (-)
  ! sub-grid routing
  case('routingGammaShape'        ); get_ixBpar = iLookBPAR%routingGammaShape         ! shape parameter in Gamma distribution used for sub-grid routing (-)
  case('routingGammaScale'        ); get_ixBpar = iLookBPAR%routingGammaScale         ! scale parameter in Gamma distribution used for sub-grid routing (s)
  ! get to here if cannot find the variable
  case default
   get_ixBpar = integerMissing
 end select
 end function get_ixBpar


 ! *******************************************************************************************************************
 ! public function get_ixBvar: get the index of the named variables for the basin-average variables
 ! *******************************************************************************************************************
 function get_ixBvar(varName)
 USE var_lookup,only:iLookBVAR                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! parameter name
 integer(i4b)             :: get_ixBvar              ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  ! derived variables
  case('basin__TotalArea'              ); get_ixBvar = iLookBVAR%basin__totalArea                ! total basin area (m2)
  ! scalar variables -- basin-average runoff and aquifer fluxes
  case('basin__SurfaceRunoff'          ); get_ixBvar = iLookBVAR%basin__SurfaceRunoff            ! surface runoff (m s-1)
  case('basin__ColumnOutflow'          ); get_ixBvar = iLookBVAR%basin__ColumnOutflow            ! outflow from all "outlet" HRUs (those with no downstream HRU)
  case('basin__AquiferStorage'         ); get_ixBvar = iLookBVAR%basin__AquiferStorage           ! aquifer storage (m s-1)
  case('basin__AquiferRecharge'        ); get_ixBvar = iLookBVAR%basin__AquiferRecharge          ! recharge to the aquifer (m s-1)
  case('basin__AquiferBaseflow'        ); get_ixBvar = iLookBVAR%basin__AquiferBaseflow          ! baseflow from the aquifer (m s-1)
  case('basin__AquiferTranspire'       ); get_ixBvar = iLookBVAR%basin__AquiferTranspire         ! transpiration from the aquifer (m s-1)
  case('basin__TotalRunoff'            ); get_ixBvar = iLookBVAR%basin__TotalRunoff              ! total runoff to channel from all active components (m s-1)
  case('basin__SoilDrainage'           ); get_ixBvar = iLookBVAR%basin__SoilDrainage             ! soil drainage (m s-1)
  ! variables to compute runoff
  case('routingRunoffFuture'           ); get_ixBvar = iLookBVAR%routingRunoffFuture             ! runoff in future time steps (m s-1)
  case('routingFractionFuture'         ); get_ixBvar = iLookBVAR%routingFractionFuture           ! fraction of runoff in future time steps (-)
  case('averageInstantRunoff'          ); get_ixBvar = iLookBVAR%averageInstantRunoff            ! instantaneous runoff (m s-1)
  case('averageRoutedRunoff'           ); get_ixBvar = iLookBVAR%averageRoutedRunoff             ! routed runoff (m s-1)
  ! get to here if cannot find the variable
  case default
   get_ixBvar = integerMissing
 end select
 end function get_ixBvar

 ! *********************************************************************************************************
 ! public function get_ixVarType: get the index of the named variable type
 ! *********************************************************************************************************
 function get_ixVarType(varType)
 USE var_lookup,only:iLookVarType                    ! indices of the named variable types
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varType                 ! variable type name
 integer(i4b)             :: get_ixVarType          ! index of the named variable type list
 ! get the index of the named variables
 select case(trim(varType))
  case('scalarv'); get_ixVarType = iLookVarType%scalarv
  case('wLength'); get_ixVarType = iLookVarType%wLength
  case('midSnow'); get_ixVarType = iLookVarType%midSnow
  case('midSoil'); get_ixVarType = iLookVarType%midSoil
  case('midToto'); get_ixVarType = iLookVarType%midToto
  case('ifcSnow'); get_ixVarType = iLookVarType%ifcSnow
  case('ifcSoil'); get_ixVarType = iLookVarType%ifcSoil
  case('ifcToto'); get_ixVarType = iLookVarType%ifcToto
  case('parSoil'); get_ixVarType = iLookVarType%parSoil
  case('routing'); get_ixVarType = iLookVarType%routing
  case('unknown'); get_ixVarType = iLookVarType%unknown
  ! get to here if cannot find the variable
  case default
   get_ixVarType = integerMissing
 end select
 end function get_ixVarType

 ! ****************************************************************************************************************
 ! public function get_varTypeName: get the index of the named variable type
 ! ****************************************************************************************************************
 function get_varTypeName(varType)
 USE var_lookup,only:iLookVarType                    ! indices of the named variable types
 implicit none
 ! define dummy variables
 integer(i4b), intent(in) :: varType                 ! variable type name
 character(LEN=7)         :: get_varTypeName         ! index of the named variable type list
 ! get the index of the named variables
 select case(varType)
  case(iLookVarType%scalarv);get_varTypeName='scalarv'
  case(iLookVarType%wLength);get_varTypeName='wLength'
  case(iLookVarType%midSnow);get_varTypeName='midSnow'
  case(iLookVarType%midSoil);get_varTypeName='midSoil'
  case(iLookVarType%midToto);get_varTypeName='midToto'
  case(iLookVarType%ifcSnow);get_varTypeName='ifcSnow'
  case(iLookVarType%ifcSoil);get_varTypeName='ifcSoil'
  case(iLookVarType%ifcToto);get_varTypeName='ifcToto'
  case(iLookVarType%parSoil);get_varTypeName='parSoil'
  case(iLookVarType%routing);get_varTypeName='routing'
  case(iLookVarType%unknown);get_varTypeName='unknown'
  ! get to here if cannot find the variable
  case default
   get_VarTypeName = 'missing'
 end select
 end function get_VarTypeName

 ! *******************************************************************************************************************
 ! public subroutine get_ixUnknown: get the index of the named variable type from ANY structure, as well as the
 ! structure that it was found in
 ! *******************************************************************************************************************
 subroutine get_ixUnknown(varName,typeName,vDex,err,message)
 USE nrtype
 USE globalData,only:structInfo        ! information on the data structures
 implicit none

 ! dummies
 character(*),intent(in)  :: varName   ! variable name
 character(*),intent(out) :: typeName  ! variable type (structure) name
 integer(i4b),intent(out) :: vDex      ! variable index in structure
 integer(i4b),intent(out) :: err       ! error code
 character(*),intent(out) :: message   ! error message

 ! internals
 integer(i4b)             :: iStruc    ! index for looping through structure types

 ! error init
 err=0
 message='get_ixUnknown/'

 ! loop through all structure types to find the one with the given variable name
 ! poll variable index plus return which structure it was found in
 do iStruc = 1,size(structInfo)
  select case(trim(structInfo(iStruc)%structName))
   case ('time' );  vDex = get_ixTime(trim(varName))
   case ('forc' );  vDex = get_ixForce(trim(varName))
   case ('attr' );  vDex = get_ixAttr(trim(varName))
   case ('type' );  vDex = get_ixType(trim(varName))
   case ('id'   );  vDex = get_ixId(trim(varName))
   case ('mpar' );  vDex = get_ixParam(trim(varName))
   case ('indx' );  vDex = get_ixIndex(trim(varName))
   case ('prog' );  vDex = get_ixProg(trim(varName))
   case ('diag' );  vDex = get_ixDiag(trim(varName))
   case ('flux' );  vDex = get_ixFlux(trim(varName))
   case ('bpar' );  vDex = get_ixBpar(trim(varName))
   case ('bvar' );  vDex = get_ixBvar(trim(varName))
   case ('deriv');  vDex = get_ixDeriv(trim(varName))
   case ('lookup'); vDex = get_ixLookup(trim(varName))
  end select
  if (vDex>0) then; typeName=trim(structInfo(iStruc)%structName); return; end if
 end do

 ! 404
 err=20;message=trim(message)//'variable '//trim(varName)//' is not found in any structure'; return

 end subroutine get_ixUnknown

 ! *******************************************************************************************************************
 ! public function get_ixFreq: get the index of the named variables for the output frequencies
 ! *******************************************************************************************************************
 function get_ixLookup(varName)
 USE var_lookup,only:iLookLOOKUP                     ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixLookup            ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  case('temperature'); get_ixLookup = iLookLOOKUP%temperature     ! temperature (K)
  case('psiLiq_int' ); get_ixLookup = iLookLOOKUP%psiLiq_int      ! integral of mLayerPsiLiq from Tfreeze to Tk (K)
  case('deriv2'     ); get_ixLookup = iLookLOOKUP%deriv2          ! secind derivative of the interpolating function
  ! get to here if cannot find the variable
  case default
   get_ixLookup = integerMissing
 end select
 end function get_ixLookup

 ! *******************************************************************************************************************
 ! public function get_ixFreq: get the index of the named variables for the output frequencies
 ! *******************************************************************************************************************
 function get_ixFreq(varName)
 USE var_lookup,only:iLookFREQ                       ! indices of the named variables
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName                 ! variable name
 integer(i4b)             :: get_ixFreq              ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  case('day'     ); get_ixFreq = iLookFREQ%day      ! daily aggregation
  case('month'   ); get_ixFreq = iLookFREQ%month    ! monthly aggregation
  case('annual'  ); get_ixFreq = iLookFREQ%annual   ! yearly (annual) aggregation
  case('timestep'); get_ixFreq = iLookFREQ%timestep ! timestep-level output (no temporal aggregation)
  ! get to here if cannot find the variable
  case default
   get_ixFreq = integerMissing
 end select
 end function get_ixFreq

 ! ***************************************************************************************************************
 ! public function get_ixStat: get the named variables for the statistics
 ! ***************************************************************************************************************
 function get_ixStat(varName)
 USE var_lookup,only:iLookSTAT                   ! indices of the possible output statistics
 implicit none
 ! define dummy variables
 character(*), intent(in) :: varName             ! variable name
 integer(i4b)             :: get_ixStat          ! index of the named variable
 ! get the index of the named variables
 select case(trim(varName))
  case('total'   ); get_ixStat = iLookSTAT%totl
  case('instant' ); get_ixStat = iLookSTAT%inst
  case('mean'    ); get_ixStat = iLookSTAT%mean
  case('variance'); get_ixStat = iLookSTAT%vari
  case('minimum' ); get_ixStat = iLookSTAT%mini
  case('maximum' ); get_ixStat = iLookSTAT%maxi
  case('mode'    ); get_ixStat = iLookSTAT%mode
  ! get to here if cannot find the variable
  case default
   get_ixStat = integerMissing
 end select
 end function get_ixStat

 ! ***************************************************************************************************************
 ! public function get_freqName: get the name of the output frequency type
 ! ***************************************************************************************************************
 function get_freqName(ifreq)
 USE var_lookup,only:iLookFREQ                   ! indices of the possible output frequencies
 implicit none
 ! define dummy variables
 integer(i4b), intent(in) :: ifreq               ! output frequency index
 character(LEN=10)        :: get_freqName        ! name of the output frequency
 ! get the index of the named variables
 select case(ifreq)
  case(iLookFREQ%day);      get_freqName='day'
  case(iLookFREQ%month);    get_freqName='month'
  case(iLookFREQ%annual);   get_freqName='annual'
  case(iLookFREQ%timestep); get_freqName='timestep'
  ! get to here if cannot find the variable
  case default
   get_freqName = 'unknown'
 end select
 end function get_freqName

 ! ***************************************************************************************************************
 ! public function get_statName: get the name of the output statistics type
 ! ***************************************************************************************************************
 function get_statName(istat)
 USE var_lookup,only:iLookSTAT                   ! indices of the possible output statistics
 implicit none
 ! define dummy variables
 integer(i4b), intent(in) :: istat               ! stat type name
 character(LEN=10)         :: get_statName       ! name of the statistic
 ! get the index of the named variables
 select case(istat)
  case(iLookSTAT%totl);get_statName='total'
  case(iLookSTAT%inst);get_statName='instant'
  case(iLookSTAT%mean);get_statName='mean'
  case(iLookSTAT%vari);get_statName='variance'
  case(iLookSTAT%mini);get_statName='minimum'
  case(iLookSTAT%maxi);get_statName='maximum'
  case(iLookSTAT%mode);get_statName='mode'
  ! get to here if cannot find the variable
  case default
   get_statName = 'unknown'
 end select
 end function get_statName

end module get_ixname_module

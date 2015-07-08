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

module mDecisions_module
USE nrtype
implicit none
private
public::mDecisions
! -----------------------------------------------------------------------------------------------------------
! ***** define look-up values for different Noah-MP decisions *****
! -----------------------------------------------------------------------------------------------------------
! look-up values for the choice of function for the soil moisture control on stomatal resistance
integer(i4b),parameter,public :: NoahType          = 1    ! thresholded linear function of volumetric liquid water content
integer(i4b),parameter,public :: CLM_Type          = 2    ! thresholded linear function of matric head
integer(i4b),parameter,public :: SiB_Type          = 3    ! exponential of the log of matric head
! look-up values for the choice of stomatal resistance formulation
integer(i4b),parameter,public :: BallBerry         = 1    ! Ball-Berry
integer(i4b),parameter,public :: Jarvis            = 2    ! Jarvis
integer(i4b),parameter,public :: simpleResistance  = 3    ! simple resistance formulation
! -----------------------------------------------------------------------------------------------------------
! ***** define look-up values for different SUMMA model decisions *****
! -----------------------------------------------------------------------------------------------------------
! look-up values for the choice of numerical method
integer(i4b),parameter,public :: iterative            =  11    ! iterative
integer(i4b),parameter,public :: nonIterative         =  12    ! non-iterative
integer(i4b),parameter,public :: iterSurfEnergyBal    =  13    ! iterate only on the surface energy balance
! look-up values for method used to compute derivative
integer(i4b),parameter,public :: numerical            =  21    ! numerical solution
integer(i4b),parameter,public :: analytical           =  22    ! analytical solution
! look-up values for method used to determine LAI and SAI
integer(i4b),parameter,public :: monthlyTable         =  31    ! LAI/SAI taken directly from a monthly table for different vegetation classes
integer(i4b),parameter,public :: specified            =  32    ! LAI/SAI computed from green vegetation fraction and winterSAI and summerLAI parameters
! look-up values for the form of Richards' equation
integer(i4b),parameter,public :: moisture             =  41    ! moisture-based form of Richards' equation
integer(i4b),parameter,public :: mixdform             =  42    ! mixed form of Richards' equation
! look-up values for the choice of groundwater parameterization
integer(i4b),parameter,public :: qbaseTopmodel        =  51    ! TOPMODEL-ish baseflow parameterization
integer(i4b),parameter,public :: bigBucket            =  52    ! a big bucket (lumped aquifer model)
integer(i4b),parameter,public :: noExplicit           =  53    ! no explicit groundwater parameterization
! look-up values for the choice of hydraulic conductivity profile
integer(i4b),parameter,public :: constant             =  61    ! constant hydraulic conductivity with depth
integer(i4b),parameter,public :: powerLaw_profile     =  62    ! power-law profile
! look-up values for the choice of boundary conditions for thermodynamics
integer(i4b),parameter,public :: prescribedTemp       =  71    ! prescribed temperature
integer(i4b),parameter,public :: energyFlux           =  72    ! energy flux
integer(i4b),parameter,public :: zeroFlux             =  73    ! zero flux
! look-up values for the choice of boundary conditions for hydrology
integer(i4b),parameter,public :: liquidFlux           =  81    ! liquid water flux
integer(i4b),parameter,public :: prescribedHead       =  82    ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
integer(i4b),parameter,public :: funcBottomHead       =  83    ! function of matric head in the lower-most layer
integer(i4b),parameter,public :: freeDrainage         =  84    ! free drainage
! look-up values for the choice of parameterization for vegetation roughness length and displacement height
integer(i4b),parameter,public :: Raupach_BLM1994      =  91    ! Raupach (BLM 1994) "Simplified expressions..."
integer(i4b),parameter,public :: CM_QJRMS1998         =  92    ! Choudhury and Monteith (QJRMS 1998) "A four layer model for the heat budget..."
integer(i4b),parameter,public :: vegTypeTable         =  93    ! constant parameters dependent on the vegetation type
! look-up values for the choice of parameterization for canopy emissivity
integer(i4b),parameter,public :: simplExp             = 101    ! simple exponential function
integer(i4b),parameter,public :: difTrans             = 102    ! parameterized as a function of diffuse transmissivity
! look-up values for the choice of parameterization for snow interception
integer(i4b),parameter,public :: stickySnow           = 111    ! maximum interception capacity an increasing function of temerature
integer(i4b),parameter,public :: lightSnow            = 112    ! maximum interception capacity an inverse function of new snow densit
! look-up values for the choice of wind profile
integer(i4b),parameter,public :: exponential          = 121    ! exponential wind profile extends to the surface
integer(i4b),parameter,public :: logBelowCanopy       = 122    ! logarithmic profile below the vegetation canopy
! look-up values for the choice of stability function
integer(i4b),parameter,public :: standard             = 131    ! standard MO similarity, a la Anderson (1976)
integer(i4b),parameter,public :: louisInversePower    = 132    ! Louis (1979) inverse power function
integer(i4b),parameter,public :: mahrtExponential     = 133    ! Mahrt (1987) exponential
! look-up values for the choice of canopy shortwave radiation method
integer(i4b),parameter,public :: noah_mp              = 141    ! full Noah-MP implementation (including albedo)
integer(i4b),parameter,public :: CLM_2stream          = 142    ! CLM 2-stream model (see CLM documentation)
integer(i4b),parameter,public :: UEB_2stream          = 143    ! UEB 2-stream model (Mahat and Tarboton, WRR 2011)
integer(i4b),parameter,public :: NL_scatter           = 144    ! Simplified method Nijssen and Lettenmaier (JGR 1999)
integer(i4b),parameter,public :: BeersLaw             = 145    ! Beer's Law (as implemented in VIC)
! look-up values for the choice of albedo representation
integer(i4b),parameter,public :: constantDecay        = 151    ! constant decay (e.g., VIC, CLASS)
integer(i4b),parameter,public :: variableDecay        = 152    ! variable decay (e.g., BATS approach, with destructive metamorphism + soot content)
! look-up values for the choice of compaction routine
integer(i4b),parameter,public :: constantSettlement   = 161    ! constant settlement rate
integer(i4b),parameter,public :: andersonEmpirical    = 162    ! semi-empirical method of Anderson (1976)
! look-up values for the choice of method to combine and sub-divide snow layers
integer(i4b),parameter,public :: sameRulesAllLayers   = 171    ! same combination/sub-division rules applied to all layers
integer(i4b),parameter,public :: rulesDependLayerIndex= 172    ! combination/sub-dividion rules depend on layer index
! look-up values for the choice of thermal conductivity representation for snow
integer(i4b),parameter,public :: Yen1965              = 181    ! Yen (1965)
integer(i4b),parameter,public :: Mellor1977           = 182    ! Mellor (1977)
integer(i4b),parameter,public :: Jordan1991           = 183    ! Jordan (1991)
integer(i4b),parameter,public :: Smirnova2000         = 184    ! Smirnova et al. (2000)
! look-up values for the choice of thermal conductivityi representation for soil
integer(i4b),parameter,public :: funcSoilWet          = 191    ! function of soil wetness
integer(i4b),parameter,public :: mixConstit           = 192    ! mixture of constituents
integer(i4b),parameter,public :: hanssonVZJ           = 193    ! test case for the mizoguchi lab experiment, Hansson et al. VZJ 2004
! look-up values for the choice of method for the spatial representation of groundwater
integer(i4b),parameter,public :: localColumn          = 201    ! separate groundwater representation in each local soil column
integer(i4b),parameter,public :: singleBasin          = 202    ! single groundwater store over the entire basin
! look-up values for the choice of sub-grid routing method
integer(i4b),parameter,public :: timeDelay            = 211    ! time-delay histogram
integer(i4b),parameter,public :: qInstant             = 212    ! instantaneous routing
! -----------------------------------------------------------------------------------------------------------
contains


 ! ************************************************************************************************
 ! public subroutine mDecisions: save model decisions as named integers
 ! ************************************************************************************************
 subroutine mDecisions(err,message)
 ! model time structures
 USE multiconst,only:secprday               ! number of seconds in a day
 USE data_struc,only:time_meta              ! time metadata
 USE var_lookup,only:iLookTIME              ! named variables that identify indices in the time structures
 USE data_struc,only:startTime,finshTime    ! start/end time of simulation
 USE data_struc,only:dJulianStart           ! julian day of start time of simulation
 USE data_struc,only:dJulianFinsh           ! julian day of end time of simulation
 USE data_struc,only:data_step              ! length of data step (s)
 USE data_struc,only:numtim                 ! number of time steps in the simulation
 ! model decision structures
 USE data_struc,only:model_decisions        ! model decision structure
 USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure
 ! Noah-MP decision structures
 USE noahmp_globals,only:DVEG               ! decision for dynamic vegetation
 USE noahmp_globals,only:OPT_RAD            ! decision for canopy radiation
 USE noahmp_globals,only:OPT_ALB            ! decision for snow albedo
 ! time utility programs
 USE time_utils_module,only:extractTime     ! extract time info from units string
 USE time_utils_module,only:compjulday      ! compute the julian day
 implicit none
 ! define output
 integer(i4b),intent(out)             :: err            ! error code
 character(*),intent(out)             :: message        ! error message
 ! define local variables
 character(len=256)                   :: cmessage       ! error message for downwind routine
 integer(i4b)                         :: nAtt           ! number of attributes in the time structures
 real(dp)                             :: dsec           ! second
 ! initialize error control
 err=0; message='mDecisions/'

 ! -------------------------------------------------------------------------------------------------
 ! -------------------------------------------------------------------------------------------------

 ! read information from model decisions file, and populate model decisions structure
 call readoption(err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! -------------------------------------------------------------------------------------------------

 ! allocate space for start/end time structures
 nAtt = size(time_meta)  ! number of attributes in the time structure
 allocate(startTime,finshTime, stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the time structures'; return; endif
 allocate(startTime%var(nAtt),finshTime%var(nAtt), stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the time structure components'; return; endif

 ! put simulation start time information into the time structures
 call extractTime(model_decisions(iLookDECISIONS%simulStart)%cDecision,  & ! date-time string
                  startTime%var(iLookTIME%iyyy),                         & ! year
                  startTime%var(iLookTIME%im),                           & ! month
                  startTime%var(iLookTIME%id),                           & ! day
                  startTime%var(iLookTIME%ih),                           & ! hour
                  startTime%var(iLookTIME%imin),                         & ! minute
                  dsec,                                                  & ! second
                  err,cmessage)                                            ! error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! put simulation end time information into the time structures
 call extractTime(model_decisions(iLookDECISIONS%simulFinsh)%cDecision,  & ! date-time string
                  finshTime%var(iLookTIME%iyyy),                         & ! year
                  finshTime%var(iLookTIME%im),                           & ! month
                  finshTime%var(iLookTIME%id),                           & ! day
                  finshTime%var(iLookTIME%ih),                           & ! hour
                  finshTime%var(iLookTIME%imin),                         & ! minute
                  dsec,                                                  & ! second
                  err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! compute the julian date (fraction of day) for the start of the simulation
 call compjulday(&
                 startTime%var(iLookTIME%iyyy),                         & ! year
                 startTime%var(iLookTIME%im),                           & ! month
                 startTime%var(iLookTIME%id),                           & ! day
                 startTime%var(iLookTIME%ih),                           & ! hour
                 startTime%var(iLookTIME%imin),                         & ! minute
                 0._dp,                                                 & ! second
                 dJulianStart,                                          & ! julian date for the start of the simulation
                 err, cmessage)                                           ! error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! compute the julian date (fraction of day) for the end of the simulation
 call compjulday(&
                 finshTime%var(iLookTIME%iyyy),                         & ! year
                 finshTime%var(iLookTIME%im),                           & ! month
                 finshTime%var(iLookTIME%id),                           & ! day
                 finshTime%var(iLookTIME%ih),                           & ! hour
                 finshTime%var(iLookTIME%imin),                         & ! minute
                 0._dp,                                                 & ! second
                 dJulianFinsh,                                          & ! julian date for the end of the simulation
                 err, cmessage)                                           ! error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! check that simulation end time is > start time
 if(dJulianFinsh < dJulianStart)then; err=20; message=trim(message)//'end time of simulation occurs before start time'; return; endif

 ! compute the number of time steps
 numtim = nint( (dJulianFinsh - dJulianStart)*secprday/data_step ) + 1

 ! -------------------------------------------------------------------------------------------------

 ! set Noah-MP options
 DVEG=3      ! option for dynamic vegetation
 OPT_RAD=3   ! option for canopy radiation
 OPT_ALB=2   ! option for snow albedo

 ! identify the choice of function for the soil moisture control on stomatal resistance
 select case(trim(model_decisions(iLookDECISIONS%soilStress)%cDecision))
  case('NoahType'); model_decisions(iLookDECISIONS%soilStress)%iDecision = NoahType             ! thresholded linear function of volumetric liquid water content
  case('CLM_Type'); model_decisions(iLookDECISIONS%soilStress)%iDecision = CLM_Type             ! thresholded linear function of matric head
  case('SiB_Type'); model_decisions(iLookDECISIONS%soilStress)%iDecision = SiB_Type             ! exponential of the log of matric head
  case default
   err=10; message=trim(message)//"unknown numerical [option="//trim(model_decisions(iLookDECISIONS%soilStress)%cDecision)//"]"; return
 end select

 ! identify the choice of function for stomatal resistance
 select case(trim(model_decisions(iLookDECISIONS%stomResist)%cDecision))
  case('BallBerry'          ); model_decisions(iLookDECISIONS%stomResist)%iDecision = BallBerry           ! Ball-Berry
  case('Jarvis'             ); model_decisions(iLookDECISIONS%stomResist)%iDecision = Jarvis              ! Jarvis
  case('simpleResistance'   ); model_decisions(iLookDECISIONS%stomResist)%iDecision = simpleResistance    ! simple resistance formulation
  case default
   err=10; message=trim(message)//"unknown numerical [option="//trim(model_decisions(iLookDECISIONS%stomResist)%cDecision)//"]"; return
 end select

 ! -------------------------------------------------------------------------------------------------

 ! identify the numerical method
 select case(trim(model_decisions(iLookDECISIONS%num_method)%cDecision))
  case('itertive'); model_decisions(iLookDECISIONS%num_method)%iDecision = iterative           ! iterative
  case('non_iter'); model_decisions(iLookDECISIONS%num_method)%iDecision = nonIterative        ! non-iterative
  case('itersurf'); model_decisions(iLookDECISIONS%num_method)%iDecision = iterSurfEnergyBal   ! iterate only on the surface energy balance
  case default
   err=10; message=trim(message)//"unknown numerical [option="//trim(model_decisions(iLookDECISIONS%num_method)%cDecision)//"]"; return
 end select

 ! identify the method used to calculate flux derivatives
 select case(trim(model_decisions(iLookDECISIONS%fDerivMeth)%cDecision))
  case('numericl'); model_decisions(iLookDECISIONS%fDerivMeth)%iDecision = numerical           ! numerical
  case('analytic'); model_decisions(iLookDECISIONS%fDerivMeth)%iDecision = analytical          ! analytical
  case default
   err=10; message=trim(message)//"unknown method used to calculate flux derivatives [option="//trim(model_decisions(iLookDECISIONS%fDerivMeth)%cDecision)//"]"; return
 end select

 ! identify the method used to determine LAI and SAI
 select case(trim(model_decisions(iLookDECISIONS%LAI_method)%cDecision))
  case('monTable');  model_decisions(iLookDECISIONS%LAI_method)%iDecision = monthlyTable       ! LAI/SAI taken directly from a monthly table for different vegetation classes
  case('specified'); model_decisions(iLookDECISIONS%LAI_method)%iDecision = specified          ! LAI/SAI computed from green vegetation fraction and winterSAI and summerLAI parameters
  case default
   err=10; message=trim(message)//"unknown method to determine LAI and SAI [option="//trim(model_decisions(iLookDECISIONS%LAI_method)%cDecision)//"]"; return
 end select

 ! identify the form of Richards' equation
 select case(trim(model_decisions(iLookDECISIONS%f_Richards)%cDecision))
  case('moisture'); model_decisions(iLookDECISIONS%f_Richards)%iDecision = moisture            ! moisture-based form
  case('mixdform'); model_decisions(iLookDECISIONS%f_Richards)%iDecision = mixdform            ! mixed form
  case default
   err=10; message=trim(message)//"unknown form of Richards' equation [option="//trim(model_decisions(iLookDECISIONS%f_Richards)%cDecision)//"]"; return
 end select

 ! identify the groundwater parameterization
 select case(trim(model_decisions(iLookDECISIONS%groundwatr)%cDecision))
  case('qTopmodl'); model_decisions(iLookDECISIONS%groundwatr)%iDecision = qbaseTopmodel       ! TOPMODEL-ish baseflow parameterization
  case('bigBuckt'); model_decisions(iLookDECISIONS%groundwatr)%iDecision = bigBucket           ! a big bucket (lumped aquifer model)
  case('noXplict'); model_decisions(iLookDECISIONS%groundwatr)%iDecision = noExplicit          ! no explicit groundwater parameterization
  case default
   err=10; message=trim(message)//"unknown groundwater parameterization [option="//trim(model_decisions(iLookDECISIONS%groundwatr)%cDecision)//"]"; return
 end select

 ! identify the hydraulic conductivity profile
 select case(trim(model_decisions(iLookDECISIONS%hc_profile)%cDecision))
  case('constant'); model_decisions(iLookDECISIONS%hc_profile)%iDecision = constant            ! constant hydraulic conductivity with depth
  case('pow_prof'); model_decisions(iLookDECISIONS%hc_profile)%iDecision = powerLaw_profile    ! power-law profile
  case default
   err=10; message=trim(message)//"unknown hydraulic conductivity profile [option="//trim(model_decisions(iLookDECISIONS%hc_profile)%cDecision)//"]"; return
 end select

 ! identify the upper boundary conditions for thermodynamics
 select case(trim(model_decisions(iLookDECISIONS%bcUpprTdyn)%cDecision))
  case('presTemp'); model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision = prescribedTemp      ! prescribed temperature
  case('nrg_flux'); model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision = energyFlux          ! energy flux
  case('zeroFlux'); model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision = zeroFlux            ! zero flux
  case default
   err=10; message=trim(message)//"unknown upper boundary conditions for thermodynamics [option="//trim(model_decisions(iLookDECISIONS%bcUpprTdyn)%cDecision)//"]"; return
 end select

 ! identify the lower boundary conditions for thermodynamics
 select case(trim(model_decisions(iLookDECISIONS%bcLowrTdyn)%cDecision))
  case('presTemp'); model_decisions(iLookDECISIONS%bcLowrTdyn)%iDecision = prescribedTemp      ! prescribed temperature
  case('zeroFlux'); model_decisions(iLookDECISIONS%bcLowrTdyn)%iDecision = zeroFlux            ! zero flux
  case default
   err=10; message=trim(message)//"unknown lower boundary conditions for thermodynamics [option="//trim(model_decisions(iLookDECISIONS%bcLowrTdyn)%cDecision)//"]"; return
 end select

 ! identify the upper boundary conditions for soil hydrology
 select case(trim(model_decisions(iLookDECISIONS%bcUpprSoiH)%cDecision))
  case('presHead'); model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision = prescribedHead      ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
  case('liq_flux'); model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision = liquidFlux          ! liquid water flux
  case default
   err=10; message=trim(message)//"unknown upper boundary conditions for soil hydrology [option="//trim(model_decisions(iLookDECISIONS%bcUpprSoiH)%cDecision)//"]"; return
 end select

 ! identify the lower boundary conditions for soil hydrology
 select case(trim(model_decisions(iLookDECISIONS%bcLowrSoiH)%cDecision))
  case('presHead'); model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision = prescribedHead      ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
  case('bottmPsi'); model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision = funcBottomHead      ! function of matric head in the lower-most layer
  case('drainage'); model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision = freeDrainage        ! free drainage
  case('zeroFlux'); model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision = zeroFlux            ! zero flux
  case default
   err=10; message=trim(message)//"unknown lower boundary conditions for soil hydrology [option="//trim(model_decisions(iLookDECISIONS%bcLowrSoiH)%cDecision)//"]"; return
 end select

 ! identify the choice of parameterization for vegetation roughness length and displacement height
 select case(trim(model_decisions(iLookDECISIONS%veg_traits)%cDecision))
  case('Raupach_BLM1994'); model_decisions(iLookDECISIONS%veg_traits)%iDecision = Raupach_BLM1994  ! Raupach (BLM 1994) "Simplified expressions..."
  case('CM_QJRMS1998'   ); model_decisions(iLookDECISIONS%veg_traits)%iDecision = CM_QJRMS1998     ! Choudhury and Monteith (QJRMS 1998) "A four layer model for the heat budget..."
  case('vegTypeTable'   ); model_decisions(iLookDECISIONS%veg_traits)%iDecision = vegTypeTable     ! constant parameters dependent on the vegetation type
  case default
   err=10; message=trim(message)//"unknown parameterization for vegetation roughness length and displacement height [option="//trim(model_decisions(iLookDECISIONS%veg_traits)%cDecision)//"]"; return
 end select

 ! identify the choice of parameterization for canopy emissivity
 select case(trim(model_decisions(iLookDECISIONS%canopyEmis)%cDecision))
  case('simplExp'); model_decisions(iLookDECISIONS%canopyEmis)%iDecision = simplExp            ! simple exponential function
  case('difTrans'); model_decisions(iLookDECISIONS%canopyEmis)%iDecision = difTrans            ! parameterized as a function of diffuse transmissivity
  case default
   err=10; message=trim(message)//"unknown parameterization for canopy emissivity [option="//trim(model_decisions(iLookDECISIONS%canopyEmis)%cDecision)//"]"; return
 end select

 ! choice of parameterization for snow interception
 select case(trim(model_decisions(iLookDECISIONS%snowIncept)%cDecision))
  case('stickySnow'); model_decisions(iLookDECISIONS%snowIncept)%iDecision = stickySnow        ! maximum interception capacity an increasing function of temerature
  case('lightSnow' ); model_decisions(iLookDECISIONS%snowIncept)%iDecision = lightSnow         ! maximum interception capacity an inverse function of new snow density
  case default
   err=10; message=trim(message)//"unknown option for snow interception capacity[option="//trim(model_decisions(iLookDECISIONS%snowIncept)%cDecision)//"]"; return
 end select

 ! identify the choice of wind profile
 select case(trim(model_decisions(iLookDECISIONS%windPrfile)%cDecision))
  case('exponential'   ); model_decisions(iLookDECISIONS%windPrfile)%iDecision = exponential      ! exponential wind profile extends to the surface
  case('logBelowCanopy'); model_decisions(iLookDECISIONS%windPrfile)%iDecision = logBelowCanopy   ! logarithmic profile below the vegetation canopy
  case default
   err=10; message=trim(message)//"unknown option for choice of wind profile[option="//trim(model_decisions(iLookDECISIONS%windPrfile)%cDecision)//"]"; return
 end select

 ! identify the choice of atmospheric stability function
 select case(trim(model_decisions(iLookDECISIONS%astability)%cDecision))
  case('standard'); model_decisions(iLookDECISIONS%astability)%iDecision = standard            ! standard MO similarity, a la Anderson (1976)
  case('louisinv'); model_decisions(iLookDECISIONS%astability)%iDecision = louisInversePower   ! Louis (1979) inverse power function
  case('mahrtexp'); model_decisions(iLookDECISIONS%astability)%iDecision = mahrtExponential    ! Mahrt (1987) exponential
  case default
   err=10; message=trim(message)//"unknown stability function [option="//trim(model_decisions(iLookDECISIONS%astability)%cDecision)//"]"; return
 end select

 ! choice of canopy shortwave radiation method
 select case(trim(model_decisions(iLookDECISIONS%canopySrad)%cDecision))
  case('noah_mp'    ); model_decisions(iLookDECISIONS%canopySrad)%iDecision = noah_mp          ! full Noah-MP implementation (including albedo)
  case('CLM_2stream'); model_decisions(iLookDECISIONS%canopySrad)%iDecision = CLM_2stream      ! CLM 2-stream model (see CLM documentation)
  case('UEB_2stream'); model_decisions(iLookDECISIONS%canopySrad)%iDecision = UEB_2stream      ! UEB 2-stream model (Mahat and Tarboton, WRR 2011)
  case('NL_scatter' ); model_decisions(iLookDECISIONS%canopySrad)%iDecision = NL_scatter       ! Simplified method Nijssen and Lettenmaier (JGR 1999)
  case('BeersLaw'   ); model_decisions(iLookDECISIONS%canopySrad)%iDecision = BeersLaw         ! Beer's Law (as implemented in VIC)
  case default
   err=10; message=trim(message)//"unknown canopy radiation method [option="//trim(model_decisions(iLookDECISIONS%canopySrad)%cDecision)//"]"; return
 end select

 ! choice of albedo representation
 select case(trim(model_decisions(iLookDECISIONS%alb_method)%cDecision))
  case('conDecay'); model_decisions(iLookDECISIONS%alb_method)%iDecision = constantDecay       ! constant decay (e.g., VIC, CLASS)
  case('varDecay'); model_decisions(iLookDECISIONS%alb_method)%iDecision = variableDecay       ! variable decay (e.g., BATS approach, with destructive metamorphism + soot content)
  case default
   err=10; message=trim(message)//"unknown option for snow albedo [option="//trim(model_decisions(iLookDECISIONS%alb_method)%cDecision)//"]"; return
 end select

 ! choice of snow compaction routine
 select case(trim(model_decisions(iLookDECISIONS%compaction)%cDecision))
  case('consettl'); model_decisions(iLookDECISIONS%compaction)%iDecision = constantSettlement  ! constant settlement rate
  case('anderson'); model_decisions(iLookDECISIONS%compaction)%iDecision = andersonEmpirical   ! semi-empirical method of Anderson (1976)
  case default
   err=10; message=trim(message)//"unknown option for snow compaction [option="//trim(model_decisions(iLookDECISIONS%compaction)%cDecision)//"]"; return
 end select

 ! choice of method to combine and sub-divide snow layers
 select case(trim(model_decisions(iLookDECISIONS%snowLayers)%cDecision))
  case('jrdn1991'); model_decisions(iLookDECISIONS%snowLayers)%iDecision = sameRulesAllLayers    ! SNTHERM option: same combination/sub-dividion rules applied to all layers
  case('CLM_2010'); model_decisions(iLookDECISIONS%snowLayers)%iDecision = rulesDependLayerIndex ! CLM option: combination/sub-dividion rules depend on layer index
  case default
   err=10; message=trim(message)//"unknown option for combination/sub-division of snow layers [option="//trim(model_decisions(iLookDECISIONS%snowLayers)%cDecision)//"]"; return
 end select

 ! choice of thermal conductivity representation for snow
 select case(trim(model_decisions(iLookDECISIONS%thCondSnow)%cDecision))
  case('tyen1965'); model_decisions(iLookDECISIONS%thCondSnow)%iDecision = Yen1965             ! Yen (1965)
  case('melr1977'); model_decisions(iLookDECISIONS%thCondSnow)%iDecision = Mellor1977          ! Mellor (1977)
  case('jrdn1991'); model_decisions(iLookDECISIONS%thCondSnow)%iDecision = Jordan1991          ! Jordan (1991)
  case('smnv2000'); model_decisions(iLookDECISIONS%thCondSnow)%iDecision = Smirnova2000        ! Smirnova et al. (2000)
  case default
   err=10; message=trim(message)//"unknown option for thermal conductivity of snow [option="//trim(model_decisions(iLookDECISIONS%thCondSnow)%cDecision)//"]"; return
 end select

 ! choice of thermal conductivity representation for soil
 select case(trim(model_decisions(iLookDECISIONS%thCondSoil)%cDecision))
  case('funcSoilWet'); model_decisions(iLookDECISIONS%thCondSoil)%iDecision = funcSoilWet      ! function of soil wetness 
  case('mixConstit' ); model_decisions(iLookDECISIONS%thCondSoil)%iDecision = mixConstit       ! mixture of constituents
  case('hanssonVZJ' ); model_decisions(iLookDECISIONS%thCondSoil)%iDecision = hanssonVZJ       ! test case for the mizoguchi lab experiment, Hansson et al. VZJ 2004
  case default
   err=10; message=trim(message)//"unknown option for thermal conductivity of soil [option="//trim(model_decisions(iLookDECISIONS%thCondSoil)%cDecision)//"]"; return
 end select

 ! choice of method for the spatial representation of groundwater
 select case(trim(model_decisions(iLookDECISIONS%spatial_gw)%cDecision))
  case('localColumn'); model_decisions(iLookDECISIONS%spatial_gw)%iDecision = localColumn       ! separate groundwater in each local soil column
  case('singleBasin'); model_decisions(iLookDECISIONS%spatial_gw)%iDecision = singleBasin       ! single groundwater store over the entire basin
  case default
   err=10; message=trim(message)//"unknown option for spatial representation of groundwater [option="//trim(model_decisions(iLookDECISIONS%spatial_gw)%cDecision)//"]"; return
 end select

 ! choice of routing method
 select case(trim(model_decisions(iLookDECISIONS%subRouting)%cDecision))
  case('timeDlay'); model_decisions(iLookDECISIONS%subRouting)%iDecision = timeDelay           ! time-delay histogram
  case('qInstant'); model_decisions(iLookDECISIONS%subRouting)%iDecision = qInstant            ! instantaneous routing
  case default
   err=10; message=trim(message)//"unknown option for sub-grid routing [option="//trim(model_decisions(iLookDECISIONS%subRouting)%cDecision)//"]"; return
 end select

 ! -----------------------------------------------------------------------------------------------------------------------------------------------
 ! check for consistency among options
 ! -----------------------------------------------------------------------------------------------------------------------------------------------

 ! check there is prescribedHead for soil hydrology when zeroFlux or prescribedTemp for thermodynamics
 !select case(model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision)
 ! case(prescribedTemp,zeroFlux)
 !  if(model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision /= prescribedHead)then
 !   message=trim(message)//'upper boundary condition for soil hydology must be presHead with presTemp and zeroFlux options for thermodynamics'
 !   err=20; return
 !  endif
 !end select

 ! check there is prescribedTemp or zeroFlux for thermodynamics when using prescribedHead for soil hydrology
 !select case(model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision)
 ! case(prescribedHead)
 !  ! check that upper boundary condition for thermodynamics is presTemp or zeroFlux
 !  select case(model_decisions(iLookDECISIONS%bcUpprTdyn)%iDecision)
 !   case(prescribedTemp,zeroFlux) ! do nothing: this is OK
 !   case default
 !    message=trim(message)//'upper boundary condition for thermodynamics must be presTemp or zeroFlux with presHead option for soil hydology'
 !    err=20; return
 !  end select
 !end select

 ! check zero flux lower boundary for topmodel baseflow option
 select case(model_decisions(iLookDECISIONS%groundwatr)%iDecision)
  case(qbaseTopmodel)
   if(model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision /= zeroFlux)then
    message=trim(message)//'lower boundary condition for soil hydology must be zeroFlux with qbaseTopmodel option for groundwater'
    err=20; return
   endif
 end select

 ! check power-law profile is selected when using topmodel baseflow option
 select case(model_decisions(iLookDECISIONS%groundwatr)%iDecision)
  case(qbaseTopmodel)
   if(model_decisions(iLookDECISIONS%hc_profile)%iDecision /= powerLaw_profile)then
    message=trim(message)//'power-law transmissivity profile must be selected when using topmodel baseflow option'
    err=20; return
   endif
 end select

 ! check bigBucket groundwater option is used when for spatial groundwater is singleBasin
 if(model_decisions(iLookDECISIONS%spatial_gw)%iDecision == singleBasin)then
  if(model_decisions(iLookDECISIONS%groundwatr)%iDecision /= bigBucket)then
   message=trim(message)//'groundwater parameterization must be bigBucket when using singleBasin for spatial_gw'
   err=20; return
  endif
 endif

 ! ensure that the LAI seaonality option is switched off (this was a silly idea, in retrospect)
 !if(model_decisions(iLookDECISIONS%LAI_method)%iDecision == specified)then
 ! message=trim(message)//'parameterization of LAI in terms of seasonal cycle of green veg fraction was a silly idea '&
 !                      //' -- the LAI_method option ["specified"] is no longer supported'
 ! err=20; return
 !endif

 end subroutine mDecisions


 ! ************************************************************************************************
 ! private subroutine readoption: read information from model decisions file
 ! ************************************************************************************************
 subroutine readoption(err,message)
 ! used to read information from model decisions file
 USE ascii_util_module,only:file_open       ! open file
 USE ascii_util_module,only:get_vlines      ! get a vector of non-comment lines
 USE summaFileManager,only:SETNGS_PATH      ! path for metadata files
 USE summaFileManager,only:M_DECISIONS      ! definition of modeling options
 USE get_ixname_module,only:get_ixdecisions ! identify index of named variable
 USE data_struc,only:model_decisions        ! model decision structure
 implicit none
 ! define output
 integer(i4b),intent(out)             :: err            ! error code
 character(*),intent(out)             :: message        ! error message
 ! define local variables
 character(len=256)                   :: cmessage       ! error message for downwind routine
 character(LEN=256)                   :: infile         ! input filename
 integer(i4b),parameter               :: unt=99         ! DK: need to either define units globally, or use getSpareUnit
 character(LEN=256),allocatable       :: charline(:)    ! vector of character strings
 integer(i4b)                         :: nDecisions     ! number of model decisions
 integer(i4b)                         :: iDecision      ! index of model decisions
 character(len=32)                    :: decision       ! name of model decision
 character(len=32)                    :: option         ! option for model decision
 integer(i4b)                         :: iVar           ! index of the decision in the data structure
 ! Start procedure here
 err=0; message='readoption/'
 ! build filename
 infile = trim(SETNGS_PATH)//trim(M_DECISIONS)
 write(*,'(2(a,1x))') 'decisions file = ', trim(infile)
 ! open file
 call file_open(trim(infile),unt,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! get a list of character strings from non-comment lines
 call get_vlines(unt,charline,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! close the file unit
 close(unt)
 ! get the number of model decisions
 nDecisions = size(charline)
 ! allocate space for the model decisions
 if(associated(model_decisions)) deallocate(model_decisions)
 allocate(model_decisions(nDecisions),stat=err)
 if(err/=0)then;err=30;message=trim(message)//"problemAllocateModelDecisions"; return; endif
 ! populate the model decisions structure
 do iDecision=1,nDecisions
  ! extract name of decision and the decision selected
  read(charline(iDecision),*,iostat=err) option, decision
  if (err/=0) then; err=30; message=trim(message)//"errorReadLine"; return; endif
  ! get the index of the decision in the data structure
  iVar = get_ixdecisions(trim(option))
  if(iVar<=0)then; err=40; message=trim(message)//"cannotFindDecisionIndex[name='"//trim(option)//"']"; return; endif
  ! populate the model decisions structure
  model_decisions(iVar)%cOption   = trim(option)
  model_decisions(iVar)%cDecision = trim(decision)
 end do
 end subroutine readoption


end module mDecisions_module

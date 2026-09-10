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

module run_oneHRU_module

! numerical recipes data types
USE nr_type

! data types
USE data_types,only:&
               dom_info,                 & ! domain info                    (i4b)
               dom_d,                    & ! x%dom(:)                       (rkind)
               dom_double,               & ! x%dom(:)%var(:)                (rkind)
               dom_intVec,               & ! x%dom(:)%var(:)%dat            (i4b)
               dom_doubleVec,            & ! x%dom(:)%var(:)%dat            (rkind)
               dom_z_vLookup,            & ! x%dom(:)%z(:)%var(:)%lookup(:) (rkind)
               var_i,                    & ! x%var(:)                (i4b)
               var_d,                    & ! x%var(:)                (rkind)
               var_dlength                 ! x%var(:)%dat            (rkind)

! access vegetation data
USE globalData,only:greenVegFrac_monthly   ! fraction of green vegetation in each month (0-1)
USE globalData,only:overwriteRSMIN         ! flag to overwrite RSMIN

! access domain types
USE globalData,only:upland                 ! horizontal domain type for upland areas
USE globalData,only:glacCln1               ! first horizontal domain type for glacier clean areas
USE globalData,only:glacCln2               ! second horizontal domain type for glacier clean areas
USE globalData,only:glacDbr                ! horizontal domain type for glacier debris areas
USE globalData,only:wetland                ! horizontal domain type for wetland areas

! provide access to Noah-MP constants
USE module_sf_noahmplsm,only:isWater       ! parameter for water land cover type

! provide access to the named variables that describe elements of parameter structures
USE var_lookup,only:iLookTYPE              ! look-up values for classification of veg, soils etc.
USE var_lookup,only:iLookATTR              ! look-up values for local attributes
USE var_lookup,only:iLookPARAM             ! look-up values for local column model parameters

! global metadata
USE globalData,only:flux_meta              ! metadata on the model fluxes

! provide access to the named variables that describe elements of variable structures
USE var_lookup,only:iLookPROG              ! look-up values for local column model prognostic (state) variables
USE var_lookup,only:iLookDIAG              ! look-up values for local column model diagnostic variables
USE var_lookup,only:iLookINDEX             ! look-up values for local column index variables

! provide access to model decisions
USE globalData,only:model_decisions        ! model decision structure
USE var_lookup,only:iLookDECISIONS         ! look-up values for model decisions

! these are needed because we cannot access them in modules locally if we might use those modules with Actors
USE globalData,only:fracJulDay             ! fractional julian days since the start of year
USE globalData,only:yearLength             ! number of days in the current year
USE globalData,only:tmZoneOffsetFracDay    ! time zone offset in fractional days

! provide access to the named variables that describe model decisions
USE mDecisions_module,only:        &       ! look-up values for LAI decisions
                      monthlyTable,&       ! LAI/SAI taken directly from a monthly table for different vegetation classes
                      specified            ! LAI/SAI computed from green vegetation fraction and winterSAI and summerLAI parameters

! ----- global variables that are modified ------------------------------------------------------------------------------------------

! Noah-MP parameters
USE NOAHMP_VEG_PARAMETERS,only:SAIM,LAIM   ! 2-d tables for stem area index and leaf area index (vegType,month)
USE NOAHMP_VEG_PARAMETERS,only:HVT,HVB     ! height at the top and bottom of vegetation (vegType)
USE noahmp_globals,only:RSMIN              ! minimum stomatal resistance (vegType)

! urban vegetation category (could be local)
USE globalData,only:urbanVegCategory       ! vegetation category for urban areas

implicit none
private
public::run_oneHRU
contains


! ************************************************************************************************
! public subroutine run_oneHRU: simulation for a single HRU
! ************************************************************************************************
subroutine run_oneHRU(&
                      ! model control
                      hru_nc,              & ! intent(in):    hru index in netcdf
                      hruId,               & ! intent(in):    hruId
                      dt_init,             & ! intent(inout): used to initialize the length of the sub-step for each HRU
                      computeVegFlux,      & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
                      ndom,                & ! intent(in):    number of domains
                      domInfo,             & ! intent(inout): domain type and number of layers
                      ! data structures (input)
                      typeData,            & ! intent(in):    local classification of soil veg etc. for each HRU
                      attrData,            & ! intent(in):    local attributes for each HRU
                      lookupData,          & ! intent(in):    local lookup tables for each HRU
                      bvarData,            & ! intent(in):    basin-average variables
                      ! data structures (input-output)
                      mparData,            & ! intent(inout): local model parameters
                      indxData,            & ! intent(inout): model indices
                      forcData,            & ! intent(inout): model forcing data
                      progData,            & ! intent(inout): prognostic variables for a local HRU
                      diagData,            & ! intent(inout): diagnostic variables for a local HRU
                      fluxData,            & ! intent(inout): model fluxes for a local HRU
                      ! error control
                      err,message)           ! intent(out):   error control
  ! ----- define downstream subroutines -----------------------------------------------------------------------------------
  USE module_sf_noahmplsm,only:redprm          ! module to assign more Noah-MP parameters
  USE derivforce_module,only:derivforce        ! module to compute derived forcing data
  USE coupled_em_module,only:coupled_em        ! module to run the coupled energy and mass model
  implicit none
  ! ----- define dummy variables ------------------------------------------------------------------------------------------
  ! model control
  integer(i4b)       , intent(in)    :: hru_nc              ! hru index in netcdf
  integer(i8b)       , intent(in)    :: hruId               ! hruId
  type(dom_d)        , intent(inout) :: dt_init             ! used to initialize the length of the sub-step for each HRU
  logical(lgt)       , intent(inout) :: computeVegFlux      ! flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
  integer(i4b)       , intent(in)    :: ndom                ! number of domains
  type(dom_info)     , intent(inout) :: domInfo(:)          ! domain type and number of layers
  ! data structures (input)
  type(var_i)        , intent(in)    :: typeData            ! x%var(:)                 -- local classification of soil veg etc. for each HRU
  type(var_d)        , intent(in)    :: attrData            ! x%var(:)                 -- local attributes for each HRU
  type(dom_z_vLookup), intent(in)    :: lookupData          ! x%dom(:)%z(:)%var(:)%lookup(:) -- local lookup tables for each HRU
  type(var_dlength)  , intent(in)    :: bvarData            ! x%var(:)%dat             -- basin-average variables
  ! data structures (input-output)
  type(dom_doubleVec), intent(in)    :: mparData            ! x%dom(:)%var(:)%dat -- local (HRU domain) model parameters
  type(dom_intVec)   , intent(inout) :: indxData            ! x%dom(:)%var(:)%dat -- model indices
  type(var_d)        , intent(inout) :: forcData            ! x%var(:)            -- model forcing data
  type(dom_doubleVec), intent(inout) :: progData            ! x%dom(:)%var(:)%dat -- model prognostic (state) variables
  type(dom_doubleVec), intent(inout) :: diagData            ! x%dom(:)%var(:)%dat -- model diagnostic variables
  type(dom_doubleVec), intent(inout) :: fluxData            ! x%dom(:)%var(:)%dat -- model fluxes
  ! error control
  integer(i4b)       , intent(out)   :: err                 ! error code
  character(*)       , intent(out)   :: message             ! error message
  ! ----- define local variables ------------------------------------------------------------------------------------------
  integer(i4b)                      :: i, iVar             ! loop index
  logical(lgt)                      :: use_computeVegFlux  ! computeVegFlux flag for the current domain
  logical(lgt)                      :: is_glac             ! flag to indicate if is a glacier domain
  type(var_d)                       :: forcData0           ! original forcings
  character(len=256)                :: cmessage            ! error message
  ! ----------------------------------------------------------------------------------------------------------------------------------------------
  ! save original forcings
  forcData0 = forcData 

  ! loop over each domain
  do i = 1, ndom
    ! initialize error control
    err=0; write(message, '(A21,I0,A10,I0,A13,I0,A2)' ) 'run_oneHRU (hru_nc = ',hru_nc,', hruId = ',hruId,', dom_type = ',domInfo(i)%dom_type,')/'

    is_glac = .false. ! reset the flag for the next domain
    forcData = forcData0 ! reset the forcings for the next domain

    ! initialize the number of flux calls
    diagData%dom(i)%var(iLookDIAG%numFluxCalls)%dat(1) = 0._rkind

    ! if water pixel or if the fraction of the domain is zero, do not run the model
    if ( typeData%var(iLookTYPE%vegTypeIndex)==isWater .or. progData%dom(i)%var(iLookPROG%DOMarea)%dat(1) <= 0._rkind )then
      ! Set wall_clock time to zero so it does not get a random value
      diagData%dom(i)%var(iLookDIAG%wallClockTime)%dat(1) = 0._rkind

      ! set fluxes to zero
      do iVar=1,size(flux_meta)
        fluxData%dom(i)%var(iVar)%dat(:) = 0
      end do

    else ! not a water pixel and area of the domain is greater than zero
      if (domInfo(i)%dom_type == upland) then

        ! populate parameters in Noah-MP modules
        ! Passing an unrealistically large number of soil layers in order to pass the check for NROOT, that is done to avoid making any changes to Noah-MP code.
        !  --> NROOT from Noah-MP veg tables (as read here) is not used in SUMMA
        call REDPRM(typeData%var(iLookTYPE%vegTypeIndex),      & ! vegetation type index
                    typeData%var(iLookTYPE%soilTypeIndex),     & ! soil type
                    typeData%var(iLookTYPE%slopeTypeIndex),    & ! slope type index
                    10000_i4b,                                 & ! number of soil layers
                    urbanVegCategory)                            ! vegetation category for urban areas

        ! overwrite the minimum resistance
        if(overwriteRSMIN) RSMIN = mparData%dom(i)%var(iLookPARAM%minStomatalResistance)%dat(1)
        ! overwrite the vegetation height
        HVT(typeData%var(iLookTYPE%vegTypeIndex)) = mparData%dom(i)%var(iLookPARAM%heightCanopyTop)%dat(1)
        HVB(typeData%var(iLookTYPE%vegTypeIndex)) = mparData%dom(i)%var(iLookPARAM%heightCanopyBottom)%dat(1)

        ! overwrite the tables for LAI and SAI
        if(model_decisions(iLookDECISIONS%LAI_method)%iDecision == specified)then
          SAIM(typeData%var(iLookTYPE%vegTypeIndex),:) = mparData%dom(i)%var(iLookPARAM%winterSAI)%dat(1)
          LAIM(typeData%var(iLookTYPE%vegTypeIndex),:) = mparData%dom(i)%var(iLookPARAM%summerLAI)%dat(1)*greenVegFrac_monthly
        end if

        use_computeVegFlux = computeVegFlux

      elseif ( domInfo(i)%dom_type == glacCln1 .or. domInfo(i)%dom_type == glacCln2 .or. domInfo(i)%dom_type == glacDbr .or. &
               domInfo(i)%dom_type == wetland )then ! don't need vegetation parameters for glaciers
        use_computeVegFlux = .false.
        if (domInfo(i)%dom_type == glacCln1 .or. domInfo(i)%dom_type == glacCln2 .or. domInfo(i)%dom_type == glacDbr)&
            is_glac = .true.
      else
        err=20; message=trim(message)//'domain type not recognized';return
      endif

      ! ----- hru forcing ----------------------------------------------------------------------------------------------------

      ! compute derived forcing variables, different for each domain type
      call derivforce(&
                      is_glac,              & ! intent(in):    flag to indicate if is a glacier domain
                      forcData,             & ! intent(inout): vector of model forcing data
                      attrData,             & ! intent(in):    vector of model attributes
                      mparData%dom(i),      & ! intent(in):    data structure of model parameters
                      progData%dom(i),      & ! intent(in):    data structure of model prognostic variables
                      diagData%dom(i),      & ! intent(inout): data structure of model diagnostic variables
                      fluxData%dom(i),      & ! intent(inout): data structure of model fluxes
                      tmZoneOffsetFracDay,  & ! intent(inout): time zone offset in fractional days
                      err,cmessage)           ! intent(out):   error control
      if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
      ! ----- run the model --------------------------------------------------------------------------------------------------

      ! run the model for a single HRU
      call coupled_em(&
                     ! model control
                     hruId,               & ! intent(in):    hruId
                     dt_init%dom(i),      & ! intent(inout): initial time step
                     1,                   & ! intent(in):    used to adjust the length of the timestep with failure in Actors (non-Actors here, always 1)
                     use_computeVegFlux,  & ! intent(inout): flag to indicate if we are computing fluxes over vegetation
                     fracJulDay,          & ! intent(in):    fractional julian days since the start of year
                     yearLength,          & ! intent(in):    number of days in the current year
                     ! data structures (input)
                     typeData,            & ! intent(in):    local classification of soil veg etc. for each HRU
                     attrData,            & ! intent(in):    local attributes for each HRU
                     forcData,            & ! intent(in):    model forcing data
                     mparData%dom(i),     & ! intent(in):    model parameters
                     bvarData,            & ! intent(in):    basin-average model variables
                     lookupData%dom(i),   & ! intent(in):    lookup tables
                     ! data structures (input-output)
                     indxData%dom(i),     & ! intent(inout): model indices
                     progData%dom(i),     & ! intent(inout): model prognostic variables for a local HRU
                     diagData%dom(i),     & ! intent(inout): model diagnostic variables for a local HRU
                     fluxData%dom(i),     & ! intent(inout): model fluxes for a local HRU
                     ! error control
                     err,cmessage)          ! intent(out):   error control
      if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

      if(domInfo(i)%dom_type == upland)then
        computeVegFlux = use_computeVegFlux ! update the flag for the next domain on upland areas
      end if
    endif ! not a water pixel and area of the domain is greater than zero

    ! update the number of layers regardless of whether the model was run
    domInfo(i)%nSnow   = indxData%dom(i)%var(iLookINDEX%nSnow)%dat(1)    ! number of snow layers
    domInfo(i)%nLake   = indxData%dom(i)%var(iLookINDEX%nLake)%dat(1)    ! number of lake layers
    domInfo(i)%nSoil   = indxData%dom(i)%var(iLookINDEX%nSoil)%dat(1)    ! number of soil layers
    domInfo(i)%nGlce   = indxData%dom(i)%var(iLookINDEX%nGlce)%dat(1)    ! number of glacier ice layers
  end do

end subroutine run_oneHRU

end module run_oneHRU_module

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

module run_oneGRU_module

! numerical recipes data types
USE nr_type

! constants
USE globalData,only: yes,no             ! .true. and .false.
USE globalData,only: data_step          ! length of data step (s)

! define data types
USE data_types,only:&
                    ! GRU-to-HRU mapping
                    gru2hru_map,       & ! HRU info
                    ! no spatial dimension
                    var_i,             & ! x%var(:)            (i4b)
                    var_d,             & ! x%var(:)            (rkind)
                    var_ilength,       & ! x%var(:)%dat        (i4b)
                    var_dlength,       & ! x%var(:)%dat        (rkind)
                    ! no variable dimension
                    hru_i,             & ! x%hru(:)            (i4b)
                    hru_dom_d,         & ! x%hru(:)%dom(:)     (rkind)
                    ! hru dimension
                    hru_int,           & ! x%hru(:)%var(:)     (i4b)
                    hru_int8,          & ! x%hru(:)%var(:)     (i8b)
                    hru_double,        & ! x%hru(:)%var(:)     (rkind)
                    hru_intVec,        & ! x%hru(:)%var(:)%dat (i4b)
                    !hru+dom dimension
                    hru_dom_intVec,    & ! x%hru(:)%dom(:)%var(:)%dat (i4b)
                    hru_dom_double,    & ! x%hru(:)%dom(:)%var(:)     (rkind)
                    hru_dom_doubleVec, & ! x%hru(:)%dom(:)%var(:)%dat (rkind)
                    ! hru+dom+z dimension
                    hru_dom_z_vLookup, & ! x%hru(:)%z(:)%var(:)%lookup(:)
                    ! grid dimension
                    grid_double          ! x%grid(:)%var(:)%dat2(:,:) (dp)


! provide access to the named variables that describe elements of parameter structures
USE var_lookup,only:iLookTYPE          ! look-up values for classification of veg, soils etc.
USE var_lookup,only:iLookID            ! look-up values for hru and gru IDs
USE var_lookup,only:iLookATTR          ! look-up values for local attributes
USE var_lookup,only:iLookFLUX          ! look-up values for local column model fluxes
USE var_lookup,only:iLookDIAG          ! look-up values model diagnostic variables
USE var_lookup,only:iLookBPAR          ! look-up values for basin-average model parameters
USE var_lookup,only:iLookBVAR          ! look-up values for basin-average model variables
USE var_lookup,only:iLookTIME          ! look-up values for model time data
USE var_lookup,only:iLookPROG          ! look-up values for model prognostic (state) variables
USE var_lookup,only:iLookPARAM         ! look-up values for model parameters

! provide access to model decisions
USE globalData,only:model_decisions    ! model decision structure
USE var_lookup,only:iLookDECISIONS     ! look-up values for model decisions

! access missing values
USE globalData,only:realMissing        ! missing real number

! access domain types
USE globalData,only:upland             ! horizontal domain type for upland areas
USE globalData,only:glacCln1           ! first horizontal domain type for glacier clean areas
USE globalData,only:glacCln2           ! second horizontal domain type for glacier clean areas
USE globalData,only:glacDbr            ! horizontal domain type for glacier debris areas
USE globalData,only:wetland            ! horizontal domain type for wetland areas

! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:       &
 qbaseTopmodel,                   & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                       & ! a big bucket (lumped aquifer model)
 noExplicit                         ! no explicit groundwater parameterization

! look-up values for the choice of method for the spatial representation of groundwater
USE mDecisions_module,only:       &
 localColumn,                     & ! separate groundwater representation in each local soil column
 singleBasin                        ! single groundwater store over the entire basin

 implicit none
private
public::run_oneGRU
contains

! ************************************************************************************************
! public subroutine run_oneGRU: simulation for a single GRU
! ************************************************************************************************
subroutine run_oneGRU(&
                      ! model control
                      gruInfo,            & ! intent(inout): HRU information for given GRU (# HRUs, #layers)
                      dt_init,            & ! intent(inout): used to initialize the length of the sub-step for each HRU
                      ixComputeVegFlux,   & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
                      ! data structures (input)
                      timeVec,            & ! intent(in):    model time data
                      typeHRU,            & ! intent(in):    local classification of soil veg etc. for each HRU
                      idHRU,              & ! intent(in):    local values of hru and gru IDs
                      attrHRU,            & ! intent(in):    local attributes for each HRU
                      lookupHRU,          & ! intent(in):    local lookup tables for each HRU
                      ! data structures (input-output)
                      mparHRU,            & ! intent(in):    local model parameters
                      bparData,           & ! intent(in):    basin model parameters
                      indxHRU,            & ! intent(inout): model indices
                      forcHRU,            & ! intent(inout): model forcing data
                      progHRU,            & ! intent(inout): prognostic variables for a local HRU
                      diagHRU,            & ! intent(inout): diagnostic variables for a local HRU
                      fluxHRU,            & ! intent(inout): model fluxes for a local HRU
                      bvarData,           & ! intent(inout): basin-average variables
                      gridData,           & ! intent(inout): basin glacier grids, may be null
                      ! error control
                      elapsedUpdateArea,  & ! intent(inout): elapsed time for updating glacier and wetland area for all GRUs (s)
                      err,message)          ! intent(out):   error control
  ! ----- define downstream subroutines -----------------------------------------------------------------------------------
  USE run_oneHRU_module,only:run_oneHRU                       ! module to run for one HRU
  USE qTimeDelay_module,only:qGlacier                         ! route water through glacier (time lapse)
  USE qTimeDelay_module,only:qOverland                        ! route water through an "unresolved" river network
  USE glacAreaChange_module,only:time_updateGlacArea          ! check if glacier area needs to be updated
  USE glacAreaChange_module,only:glacAreaChange               ! change glacier area with ice flow model
  USE glacAreaChange_module,only:updateGlacDomain             ! change glacier domain area, elevation, layering
  USE time_utils_module,only:elapsedSec                       ! calculate the elapsed time
  ! ----- define dummy variables ------------------------------------------------------------------------------------------
  implicit none
  ! model control
  type(gru2hru_map)       , intent(inout) :: gruInfo              ! HRU information for given GRU (# HRUs, #layers)
  type(hru_dom_d)         , intent(inout) :: dt_init              ! used to initialize the length of the sub-step for each domain
  type(hru_i)             , intent(inout) :: ixComputeVegFlux     ! flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
  ! data structures (input)
  type(var_i)             , intent(in)    :: timeVec              ! x%var(:)                               -- model time data
  type(hru_int)           , intent(in)    :: typeHRU              ! x%hru(:)%var(:)                        -- local classification of soil veg etc. for each HRU
  type(hru_int8)          , intent(in)    :: idHRU                ! x%hru(:)%var(:)                        -- local values of hru and gru IDs
  type(hru_double)        , intent(in)    :: attrHRU              ! x%hru(:)%var(:)                        -- local attributes for each HRU
  type(hru_dom_z_vLookup) , intent(in)    :: lookupHRU            ! x%hru(:)%dom(:)%z(:)%var(:)%lookup(:) -- lookup values for each HRU
  ! data structures (input-output)
  type(hru_dom_doubleVec) , intent(in)    :: mparHRU              ! x%hru(:)%dom(:)%var(:)%dat   -- local (HRU) model parameters
  type(var_d)             , intent(in)    :: bparData             ! x%var                        -- basin-average parameters
  type(hru_dom_intVec)    , intent(inout) :: indxHRU              ! x%hru(:)%dom(:)%var(:)%dat   -- model indices
  type(hru_double)        , intent(inout) :: forcHRU              ! x%hru(:)%dom(:)%var(:)       -- model forcing data
  type(hru_dom_doubleVec) , intent(inout) :: progHRU              ! x%hru(:)%dom(:)%var(:)%dat   -- model prognostic (state) variables
  type(hru_dom_doubleVec) , intent(inout) :: diagHRU              ! x%hru(:)%dom(:)%var(:)%dat   -- model diagnostic variables
  type(hru_dom_doubleVec) , intent(inout) :: fluxHRU              ! x%hru(:)%dom(:)%var(:)%dat   -- model fluxes
  type(var_dlength)       , intent(inout) :: bvarData             ! x%var(:)%dat                 -- basin-average variables
  type(grid_double)       , intent(inout) :: gridData             ! x%grid(:)%var(:)%dat2(:,:)   -- basin grids, currently used for glaciers only
  ! error control
  real(rkind)             , intent(inout) :: elapsedUpdateArea    ! time for updating glacier and wetland area for all GRUs (s)
  integer(i4b)            , intent(out)   :: err                  ! error code
  character(*)            , intent(out)   :: message              ! error message
  ! ----- define local variables ------------------------------------------------------------------------------------------
  character(len=512)                  :: cmessage                       ! error message
  integer(i4b)                        :: iHRU,jHRU,kHRU                 ! HRU indices
  integer(i4b)                        :: iDOM                           ! domain index
  integer(i4b)                        :: iUpland                        ! index of the upland domain in an HRU
  integer(i4b)                        :: iSeq                           ! position in the cascade processing order
  integer(i4b)                        :: nOrder                         ! number of HRUs placed in the processing order
  integer(i4b), allocatable           :: downIdx(:)                     ! index of the downslope HRU (0 = GRU outlet)
  integer(i4b), allocatable           :: inDegree(:)                    ! number of HRUs draining into a given HRU
  integer(i4b), allocatable           :: hruOrder(:)                    ! HRU indices in cascade order, upslope before downslope
  logical(lgt)                        :: runHRU                         ! flag to run the HRU (it has area)
  logical(lgt)                        :: computeVegFluxFlag             ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
  real(rkind)                         :: fracDOM                        ! fractional area of a given HRU domain in GRU (-)
  real(rkind)                         :: ablMelt                        ! melt from the ablation part of a glacier domain (m s-1)
  real(rkind)                         :: glacIceMelt                    ! glacier ice reservoir melt (m s-1)
  real(rkind)                         :: glacSnowMelt                   ! glacier snow reservoir melt (m s-1)
  real(rkind)                         :: glacFirnMelt                   ! glacier firn reservoir melt (m s-1)
  ! glacier area update
  logical(lgt)                        :: updateGlacArea                 ! flag to update glacier area this time step
  logical(lgt)                        :: updateLakeArea                 ! flag to update wetland area this time step
  logical(lgt)                        :: checkedUpdateTime              ! flag that the glacier update time has been checked for this GRU
  logical(lgt)                        :: hasGlacier                     ! flag that the current HRU has a glacier domain
  real(rkind)                         :: sec_since_last_update          ! seconds since last glacier area update
  integer(i4b)                        :: nglacDOM,iglacDOM              ! number of glacier domains in the GRU, and index
  integer(i4b)                        :: nglacHRU,iglacHRU              ! number of HRUs with a glacier domain in the GRU, and index
  integer(i4b), allocatable           :: nclean(:)                      ! number of clean glacier domains in each glacier HRU
  integer(i4b), allocatable           :: ndebris(:)                     ! number of debris glacier domains in each glacier HRU
  integer(i8b), allocatable           :: glac_hru(:)                    ! HRU index of each glacier domain
  real(rkind), allocatable            :: glac_area(:)                   ! area of each glacier domain (m2)
  real(rkind), allocatable            :: glac_elev(:)                   ! elevation of each glacier domain (m)
  real(rkind), allocatable            :: glac_tan_slope(:)              ! tan local ground surface slope of each glacier domain (m/m)
  real(rkind), allocatable            :: glac_aspect(:)                 ! azimuth in degrees East of North of each glacier domain (degrees)
  real(rkind), allocatable            :: glac_contourLength(:)          ! length of contour at downslope edge of each glacier domain (m)
  real(rkind), allocatable            :: glac_debris_thick(:)           ! debris thickness of each glacier domain (m)
  real(rkind), allocatable            :: glac_ablFrac(:)                ! ablation fraction of each glacier domain (-)
  real(rkind), allocatable            :: massChange(:)                  ! glacier water equivalent change of each glacier domain since the last update (kg m-2)
  real(rkind), allocatable            :: iden_soil_mean(:)              ! depth-weighted mean debris density of each glacier domain (kg m-3)
  real(rkind), allocatable            :: theta_sat_mean(:)              ! depth-weighted mean debris porosity of each glacier domain (-)
  real(rkind)                         :: soil_thick                     ! debris (soil) thickness of a debris domain (m)
  real(rkind)                         :: remaining_area                 ! HRU area not taken by glacier or wetland domains (m2)
  real(rkind)                         :: remaining_elev                 ! area-weighted elevation of the remaining area (m m2)
  real(rkind)                         :: remaining_tan_slope            ! area-weighted tan slope of the remaining area (m2)
  real(rkind)                         :: remaining_aspect_sin           ! area-weighted sine of the aspect of the remaining area (m2)
  real(rkind)                         :: remaining_aspect_cos           ! area-weighted cosine of the aspect of the remaining area (m2)
  integer(i4b),dimension(8)           :: startUpdateArea,endUpdateArea  ! time at start and end of updating glacier and wetland area
  real(rkind),parameter               :: deg2rad=PI_D/180._rkind        ! convert degrees to radians
  real(rkind),parameter               :: rad2deg=180._rkind/PI_D        ! convert radians to degrees
  real(rkind),parameter               :: aspect_tol=1.e-12_rkind        ! tolerance for undefined circular mean
  ! ----------------------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; write(message, '(A21,I0,A10,I0,A2)' ) 'run_oneGRU (gru_nc = ',gruInfo%gru_nc,', gruId = ',gruInfo%gru_id,')/'

  if(model_decisions(iLookDECISIONS%spatial_gw)%iDecision == singleBasin)then
    message=trim(message)//'multi_driver/bigBucket groundwater code not transferred from old code base yet'
    err=20; return
  endif

  ! ----- basin initialization --------------------------------------------------------------------------------------------
  associate(bvar => bvarData%var)
    ! runoff variables
    bvar(iLookBVAR%basin__SurfaceRunoff)%dat(1)    = 0._rkind ! surface runoff (m s-1)
    bvar(iLookBVAR%basin__SoilDrainage)%dat(1)     = 0._rkind ! soil drainage (m s-1)
    bvar(iLookBVAR%basin__ColumnOutflow)%dat(1)    = 0._rkind ! outflow from all "outlet" HRUs (those with no downstream HRU)
    bvar(iLookBVAR%basin__TotalRunoff)%dat(1)      = 0._rkind ! total runoff to the channel from all active components (m s-1)
    ! baseflow variables
    bvar(iLookBVAR%basin__AquiferRecharge)%dat(1)  = 0._rkind ! recharge to the aquifer (m s-1)
    bvar(iLookBVAR%basin__AquiferBaseflow)%dat(1)  = 0._rkind ! baseflow from the aquifer (m s-1)
    bvar(iLookBVAR%basin__AquiferTranspire)%dat(1) = 0._rkind ! transpiration loss from the aquifer (m s-1)
    ! storage change and glacier variables
    bvar(iLookBVAR%basin__StorageChange)%dat(1)    = 0._rkind ! change in total basin storage (kg m-2 s-1)
    bvar(iLookBVAR%basin__GlacierArea)%dat(1)      = 0._rkind ! glacier area (m2)
  end associate
  glacIceMelt    = 0._rkind
  glacSnowMelt   = 0._rkind
  glacFirnMelt   = 0._rkind
  updateGlacArea = .false.
  updateLakeArea = .false.

  ! ----- initialize the column inflows, count the glacier domains, and check if the glacier area is updated this step ------
  nglacDOM = 0
  nglacHRU = 0
  checkedUpdateTime = .false.
  do iHRU=1,gruInfo%hruCount
    hasGlacier = .false.
    do iDOM = 1, gruInfo%hruInfo(iHRU)%domCount
      associate(typeDOM => gruInfo%hruInfo(iHRU)%domInfo(iDOM)%dom_type, &
                DOMarea => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1))
        if(typeDOM==wetland)then; err=20; message=trim(message)//'ERROR:  wetland fluxes not yet implemented'; return; endif
        fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%mLayerColumnInflow)%dat(:) = 0._rkind
        if(typeDOM/=glacCln1 .and. typeDOM/=glacCln2 .and. typeDOM/=glacDbr) cycle
        nglacDOM = nglacDOM + 1
        if(.not.hasGlacier) nglacHRU = nglacHRU + 1
        hasGlacier = .true.
        ! the glacier area is updated on the first of the lowest mass month of the year (October in the Northern Hemisphere,
        !  April in the Southern, January at low latitudes); check once per GRU, from the first glacier domain with area
        if(DOMarea>0._rkind .and. .not.checkedUpdateTime)then
          call time_updateGlacArea(&
                      ! input
                      timeVec%var(iLookTIME%iyyy),timeVec%var(iLookTIME%im),timeVec%var(iLookTIME%id), timeVec%var(iLookTIME%ih),timeVec%var(iLookTIME%imin), & ! intent(in): current model time
                      attrHRU%hru(iHRU)%var(iLookATTR%latitude),        & ! intent(in): latitude of HRU (degrees)
                      ! output
                      bvarData%var(iLookBVAR%updateJulDay)%dat(1),      & ! intent(inout): julian day of last glacier area update (fraction of day)
                      bvarData%var(iLookBVAR%updateJulDayNext)%dat(1),  & ! intent(inout): julian day of next glacier area update (fraction of day)
                      updateGlacArea,                                   & ! intent(inout): flag to update glacier area this time step
                      sec_since_last_update,                            & ! intent(out):   seconds since last glacier area update
                      ! error control
                      err, cmessage)                                       ! intent(out):   error control
          if(err/=0)then; err=30; message=trim(message)//trim(cmessage); return; endif
          checkedUpdateTime = .true.
        endif
      end associate
    enddo ! (looping through domains)
  enddo ! (looping through HRUs)
  if(updateGlacArea)then
    allocate(glac_hru(nglacDOM),glac_area(nglacDOM),glac_elev(nglacDOM),glac_tan_slope(nglacDOM),glac_aspect(nglacDOM), &
             glac_contourLength(nglacDOM),glac_debris_thick(nglacDOM),glac_ablFrac(nglacDOM),massChange(nglacDOM), &
             iden_soil_mean(nglacDOM),theta_sat_mean(nglacDOM),nclean(nglacHRU),ndebris(nglacHRU))
  endif

  ! ----- order the HRUs so that an HRU is run after everything that drains into it -----------------------------------------
  allocate(downIdx(gruInfo%hruCount), inDegree(gruInfo%hruCount), hruOrder(gruInfo%hruCount), stat=err)
  if(err/=0)then; message=trim(message)//'problem allocating cascade ordering arrays'; return; endif
  downIdx(:) = 0; inDegree(:) = 0
  do iHRU=1,gruInfo%hruCount
    dsHRU: do jHRU=1,gruInfo%hruCount
      if(typeHRU%hru(iHRU)%var(iLookTYPE%downHRUindex) == idHRU%hru(jHRU)%var(iLookID%hruId))then
        downIdx(iHRU) = jHRU                  ! first match wins, as before
        inDegree(jHRU) = inDegree(jHRU) + 1
        exit dsHRU
      endif
    enddo dsHRU
  enddo
  ! repeatedly take an HRU nothing drains into, then remove its own contribution
  nOrder = 0
  do iHRU=1,gruInfo%hruCount
    if(inDegree(iHRU)==0)then; nOrder = nOrder + 1; hruOrder(nOrder) = iHRU; endif
  enddo
  iSeq = 0
  do while(iSeq < nOrder)
    iSeq = iSeq + 1
    kHRU = downIdx(hruOrder(iSeq))
    if(kHRU > 0)then
      inDegree(kHRU) = inDegree(kHRU) - 1
      if(inDegree(kHRU)==0)then; nOrder = nOrder + 1; hruOrder(nOrder) = kHRU; endif
    endif
  end do
  if(nOrder /= gruInfo%hruCount)then
    err=20; message=trim(message)//'the downHRUindex cascade network contains a loop, so the HRUs cannot be ordered upslope to &
      &downslope (check downHRUindex in the attributes file)'; return
  endif

  ! ********** RUN FOR ONE HRU ********************************************************************************************
  do iSeq=1,gruInfo%hruCount
    iHRU = hruOrder(iSeq)

    ! skip HRUs with no area
    runHRU = .false.
    do iDOM = 1, gruInfo%hruInfo(iHRU)%domCount
      if(progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1)>0._rkind) runHRU = .true.
    enddo
    if(.not.runHRU) cycle

    ! ----- run the model --------------------------------------------------------------------------------------------------
    computeVegFluxFlag = (ixComputeVegFlux%hru(iHRU) == yes)
    call run_oneHRU(&
                   ! model control
                   gruInfo%hruInfo(iHRU)%hru_nc,   & ! intent(in):    hru count Id
                   gruInfo%hruInfo(iHRU)%hru_id,   & ! intent(in):    hruId
                   dt_init%hru(iHRU),              & ! intent(inout): initial time step
                   computeVegFluxFlag,             & ! intent(inout): flag to indicate if we are computing fluxes over vegetation (false=no, true=yes)
                   gruInfo%hruInfo(iHRU)%domCount, & ! intent(in):    total number of domains
                   gruInfo%hruInfo(iHRU)%domInfo,  & ! intent(inout): domain type and layer information
                   ! data structures (input)
                   typeHRU%hru(iHRU),              & ! intent(in):    local classification of soil veg etc. for each HRU
                   attrHRU%hru(iHRU),              & ! intent(in):    local attributes for each HRU
                   lookupHRU%hru(iHRU),            & ! intent(in):    local lookup tables for each HRU
                   bvarData,                       & ! intent(in):    basin-average model variables
                   ! data structures (input-output)
                   mparHRU%hru(iHRU),              & ! intent(in):    model parameters
                   indxHRU%hru(iHRU),              & ! intent(inout): model indices
                   forcHRU%hru(iHRU),              & ! intent(inout): model forcing data
                   progHRU%hru(iHRU),              & ! intent(inout): model prognostic variables for a local HRU
                   diagHRU%hru(iHRU),              & ! intent(inout): model diagnostic variables for a local HRU
                   fluxHRU%hru(iHRU),              & ! intent(inout): model fluxes for a local HRU
                   ! error control
                   err,cmessage)                      ! intent(out):   error control
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
    ixComputeVegFlux%hru(iHRU) = merge(yes, no, computeVegFluxFlag)

    ! ----- lateral flow to the downslope HRU, and area-weighted basin (GRU) fluxes ------------------------------------------
    kHRU = downIdx(iHRU) ! the downslope HRU, found once above with the cascade ordering
    do iDOM = 1, gruInfo%hruInfo(iHRU)%domCount
      associate(typeDOM   => gruInfo%hruInfo(iHRU)%domInfo(iDOM)%dom_type, &
                flux      => fluxHRU%hru(iHRU)%dom(iDOM)%var, &
                prog      => progHRU%hru(iHRU)%dom(iDOM)%var, &
                diag      => diagHRU%hru(iHRU)%dom(iDOM)%var, &
                bvar      => bvarData%var, &
                DOMarea   => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1), &
                totalArea => bvarData%var(iLookBVAR%basin__totalArea)%dat(1))
        if(DOMarea==0._rkind) cycle ! skip domains with no area
        fracDOM = DOMarea/totalArea

        ! upland outflow goes to the downslope HRU (m3 s-1), or to the basin (GRU) column outflow if there is none
        if(typeDOM==upland)then
          if(kHRU > 0)then
            fluxHRU%hru(kHRU)%dom(1)%var(iLookFLUX%mLayerColumnInflow)%dat(:) = fluxHRU%hru(kHRU)%dom(1)%var(iLookFLUX%mLayerColumnInflow)%dat(:) + flux(iLookFLUX%mLayerColumnOutflow)%dat(:)
          else
            bvar(iLookBVAR%basin__ColumnOutflow)%dat(1) = bvar(iLookBVAR%basin__ColumnOutflow)%dat(1) + sum(flux(iLookFLUX%mLayerColumnOutflow)%dat(:))
          endif
        endif

        bvar(iLookBVAR%basin__StorageChange)%dat(1) = bvar(iLookBVAR%basin__StorageChange)%dat(1) + diag(iLookDIAG%scalarTotalMassChange)%dat(1)*fracDOM
        if(typeDOM==upland)then
          bvar(iLookBVAR%basin__SurfaceRunoff)%dat(1) = bvar(iLookBVAR%basin__SurfaceRunoff)%dat(1) + flux(iLookFLUX%scalarSurfaceRunoff)%dat(1)*fracDOM
          bvar(iLookBVAR%basin__SoilDrainage)%dat(1)  = bvar(iLookBVAR%basin__SoilDrainage)%dat(1)  + flux(iLookFLUX%scalarSoilDrainage)%dat(1) *fracDOM
          ! aquifer fluxes, only if the aquifer is computed for each HRU (singleBasin is computed later; glaciers have no groundwater)
          if(model_decisions(iLookDECISIONS%spatial_gw)%iDecision == localColumn .and. model_decisions(iLookDECISIONS%groundwatr)%iDecision == bigBucket)then
            bvar(iLookBVAR%basin__AquiferRecharge)%dat(1)  = bvar(iLookBVAR%basin__AquiferRecharge)%dat(1)  + flux(iLookFLUX%scalarAquiferRecharge)%dat(1) *fracDOM
            bvar(iLookBVAR%basin__AquiferTranspire)%dat(1) = bvar(iLookBVAR%basin__AquiferTranspire)%dat(1) + flux(iLookFLUX%scalarAquiferTranspire)%dat(1)*fracDOM
            bvar(iLookBVAR%basin__AquiferBaseflow)%dat(1)  = bvar(iLookBVAR%basin__AquiferBaseflow)%dat(1)  + flux(iLookFLUX%scalarAquiferBaseflow)%dat(1) *fracDOM
          endif
        else if(typeDOM==glacCln1 .or. typeDOM==glacCln2 .or. typeDOM==glacDbr)then
          ! glacier melt (m s-1) into the firn reservoir from the accumulation zone, and the snow or ice reservoir from the ablation zone
          ! NOTE: assumes either one HRU per GRU with many glaciers, or one glacier per GRU with many HRUs;
          !       glaciers that are absent from a particular glacier HRU are not captured
          associate(glacierMelt => flux(iLookFLUX%scalarGlacierMelt)%dat(1), ablFrac => prog(iLookPROG%scalarAblFrac)%dat(1))
            glacFirnMelt = glacFirnMelt + glacierMelt*fracDOM*(1._rkind - ablFrac) ! no debris in the accumulation zone for lateral flow
            ablMelt = (glacierMelt + sum(flux(iLookFLUX%mLayerColumnOutflow)%dat(:))/totalArea)*fracDOM*ablFrac
          end associate
          if(prog(iLookPROG%scalarSnowDepth)%dat(1)>0._rkind)then
            glacSnowMelt = glacSnowMelt + ablMelt
          else
            glacIceMelt  = glacIceMelt  + ablMelt
          endif
          bvar(iLookBVAR%basin__GlacierArea)%dat(1)    = bvar(iLookBVAR%basin__GlacierArea)%dat(1) + DOMarea ! m2
          bvar(iLookBVAR%basin__GlacierStorage)%dat(1) = bvar(iLookBVAR%basin__GlacierStorage)%dat(1) &
                                                         + diag(iLookDIAG%scalarTotalMassChange)%dat(1)*data_step*DOMarea*1.e-12_rkind ! Gt (km3 of water equivalent)
        endif ! (if domain type)
      end associate
    enddo ! (looping through domains)
  enddo  ! (looping through HRUs)
  ! ********** END LOOP THROUGH HRUS **************************************************************************************

  ! ----- collect the state of each glacier domain for the area update ----------------------------------------------------
  if(updateGlacArea)then
    glac_area          = 0._rkind
    glac_elev          = realMissing
    glac_tan_slope     = realMissing
    glac_aspect        = realMissing
    glac_contourLength = 0._rkind
    glac_debris_thick  = 0._rkind
    massChange         = 0._rkind
    iden_soil_mean     = 0._rkind
    theta_sat_mean     = 0._rkind
    nclean             = 0
    ndebris            = 0
    iglacHRU = 0
    iglacDOM = 0
    do iHRU=1,gruInfo%hruCount
      hasGlacier = .false.
      do iDOM = 1, gruInfo%hruInfo(iHRU)%domCount
        associate(typeDOM => gruInfo%hruInfo(iHRU)%domInfo(iDOM)%dom_type, &
                  nSnow   => gruInfo%hruInfo(iHRU)%domInfo(iDOM)%nSnow, &
                  nLake   => gruInfo%hruInfo(iHRU)%domInfo(iDOM)%nLake, &
                  nSoil   => gruInfo%hruInfo(iHRU)%domInfo(iDOM)%nSoil, &
                  prog    => progHRU%hru(iHRU)%dom(iDOM)%var, &
                  mpar    => mparHRU%hru(iHRU)%dom(iDOM)%var)
          if(typeDOM/=glacCln1 .and. typeDOM/=glacCln2 .and. typeDOM/=glacDbr) cycle
          iglacDOM = iglacDOM + 1
          if(.not.hasGlacier) iglacHRU = iglacHRU + 1
          hasGlacier = .true.
          glac_hru(iglacDOM) = iHRU
          if(typeDOM==glacDbr)then
            ndebris(iglacHRU) = ndebris(iglacHRU) + 1
          else
            nclean(iglacHRU) = nclean(iglacHRU) + 1
          endif
          if(prog(iLookPROG%DOMarea)%dat(1)<=0._rkind) cycle ! a domain with no area keeps the missing values
          glac_area(iglacDOM)          = prog(iLookPROG%DOMarea)%dat(1)
          glac_elev(iglacDOM)          = prog(iLookPROG%DOMelev)%dat(1)
          glac_tan_slope(iglacDOM)     = prog(iLookPROG%DOMtan_slope)%dat(1)
          glac_aspect(iglacDOM)        = prog(iLookPROG%DOMaspect)%dat(1)
          glac_contourLength(iglacDOM) = prog(iLookPROG%DOMcontourLength)%dat(1)
          massChange(iglacDOM)         = prog(iLookPROG%glacMass4AreaChange)%dat(1)
          ! the debris of a debris domain is its soil column: thickness, and depth-weighted density and porosity
          if(typeDOM==glacDbr)then
            associate(soilDepth => prog(iLookPROG%mLayerDepth)%dat(nSnow+nLake+1:nSnow+nLake+nSoil))
              soil_thick = sum(soilDepth)
              glac_debris_thick(iglacDOM) = soil_thick
              iden_soil_mean(iglacDOM)    = sum(mpar(iLookPARAM%soil_dens_intr)%dat(1:nSoil)*soilDepth)/soil_thick
              theta_sat_mean(iglacDOM)    = sum(mpar(iLookPARAM%theta_sat)%dat(1:nSoil)*soilDepth)/soil_thick
            end associate
          endif
        end associate
      enddo ! (looping through domains)
    enddo ! (looping through HRUs)
  endif ! (if need to update glacier area)

  ! ----- basin runoff and routing ----------------------------------------------------------------------------------------
  ! lapse glacier fluxes to the basin by routing through each glacier
  call qGlacier(&
                ! input
                bparData%var(iLookBPAR%glacStor_kIce),              & ! intent(in):    storage coefficient ice reservoir (hours)
                bparData%var(iLookBPAR%glacStor_kSnow),             & ! intent(in):    storage coefficient snow reservoir (hours)
                bparData%var(iLookBPAR%glacStor_kFirn),             & ! intent(in):    storage coefficient firn reservoir (hours)
                glacIceMelt,                                        & ! intent(in):    total melt into ice reservoirs (m s-1)
                glacSnowMelt,                                       & ! intent(in):    total melt into snow reservoirs (m s-1)
                glacFirnMelt,                                       & ! intent(in):    total melt into firn reservoirs (m s-1)
                bvarData%var(iLookBVAR%glacierAblArea)%dat,         & ! intent(in):    per glacier ablation area (m2)
                bvarData%var(iLookBVAR%glacierAccArea)%dat,         & ! intent(in):    per glacier accumulation area (m2)
                gruInfo%nGlac,                                      & ! intent(in):    number of glaciers in GRU
                ! output
                bvarData%var(iLookBVAR%glacIceRunoffFuture)%dat,    & ! intent(inout): per glacier ice reservoir runoff in future time steps (m s-1)
                bvarData%var(iLookBVAR%glacSnowRunoffFuture)%dat,   & ! intent(inout): per glacier snow reservoir runoff in future time steps (m s-1)
                bvarData%var(iLookBVAR%glacFirnRunoffFuture)%dat,   & ! intent(inout): per glacier firn reservoir runoff in future time steps (m s-1)
                bvarData%var(iLookBVAR%glacierRoutedRunoff)%dat(1), & ! intent(out):   routed glacier runoff (m s-1)
                err,cmessage)              ! error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

  associate(bvar => bvarData%var, totalArea => bvarData%var(iLookBVAR%basin__totalArea)%dat(1))
    ! total runoff: with a deep aquifer the column outflow is zero; without one, either the column outflow
    !  (shallow groundwater) or the soil drainage is zero
    if(model_decisions(iLookDECISIONS%groundwatr)%iDecision == bigBucket)then
      bvar(iLookBVAR%basin__TotalRunoff)%dat(1) = bvar(iLookBVAR%basin__SurfaceRunoff)%dat(1) + bvar(iLookBVAR%basin__ColumnOutflow)%dat(1)/totalArea + bvar(iLookBVAR%basin__AquiferBaseflow)%dat(1)
    else
      bvar(iLookBVAR%basin__TotalRunoff)%dat(1) = bvar(iLookBVAR%basin__SurfaceRunoff)%dat(1) + bvar(iLookBVAR%basin__ColumnOutflow)%dat(1)/totalArea + bvar(iLookBVAR%basin__SoilDrainage)%dat(1)
    endif

    call qOverland(&
                   ! input
                   model_decisions(iLookDECISIONS%subRouting)%iDecision, & ! intent(in):    index for routing method
                   bvar(iLookBVAR%basin__TotalRunoff)%dat(1),            & ! intent(in):    total runoff to the channel from all active components (m s-1)
                   bvar(iLookBVAR%routingFractionFuture)%dat,            & ! intent(in):    fraction of runoff in future time steps (m s-1)
                   bvar(iLookBVAR%routingRunoffFuture)%dat,              & ! intent(inout): runoff in future time steps (m s-1)
                   ! output
                   bvar(iLookBVAR%averageInstantRunoff)%dat(1),          & ! intent(out):   instantaneous runoff (m s-1)
                   bvar(iLookBVAR%averageRoutedRunoff)%dat(1),           & ! intent(out):   routed runoff (m s-1)
                   err,cmessage)                                           ! intent(out):   error control
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

    ! add glacier runoff to overland runoff
    bvar(iLookBVAR%averageInstantRunoff)%dat(1) = bvar(iLookBVAR%averageInstantRunoff)%dat(1) + glacIceMelt + glacSnowMelt + glacFirnMelt
    bvar(iLookBVAR%averageRoutedRunoff)%dat(1)  = bvar(iLookBVAR%averageRoutedRunoff)%dat(1)  + bvar(iLookBVAR%glacierRoutedRunoff)%dat(1)
  end associate

  ! ----- update the glacier area, and the glacier and upland domains ------------------------------------------------------
  call date_and_time(values=startUpdateArea)
  if(updateGlacArea)then
    call glacAreaChange(&
                  ! model control
                  sec_since_last_update,                      & ! intent(in):    seconds since last glacier area update
                  nglacHRU,                                   & ! intent(in):    number of HRUs that have a glacier domain
                  nglacDOM,                                   & ! intent(in):    number of domains that have glaciers
                  ndebris,                                    & ! intent(in):    number of debris domains in each HRU
                  nclean,                                     & ! intent(in):    number of clean domains in each HRU
                  glac_hru,                                   & ! intent(in):    HRU index of glacier domain
                  ! glacier topography
                  gruInfo%nGlac,                              & ! intent(in):    number of glaciers in GRU
                  gruInfo%glacInfo,                           & ! intent(in):    information for each glacier
                  gruInfo%gridInfo,                           & ! intent(in):    grid information for each grid
                  gridData,                                   & ! intent(inout): grid data for each grid
                  ! mass balance per glacier domain
                  massChange,                                 & ! intent(in):    glacier water equivalent change since updateJulDay (kg m-2)
                  glac_elev,                                  & ! intent(inout): elevation of each glacier domain (m)
                  glac_tan_slope,                             & ! intent(inout): tan local ground surface slope of each glacier domain (m/m)
                  glac_aspect,                                & ! intent(inout): azimuth in degrees East of North of each glacier domain (degrees)
                  glac_contourLength,                         & ! intent(inout): length of contour at downslope edge of each glacier domain (m)
                  ! debris
                  glac_debris_thick,                          & ! intent(inout): debris thickness of each glacier domain (m)
                  iden_soil_mean,                             & ! intent(in):    mean soil density (kg m-3)
                  theta_sat_mean,                             & ! intent(in):    mean soil porosity (-)
                  bparData%var(iLookBPAR%debrisConc),         & ! intent(in):    englacial debris concentration (kg m-3)
                  bparData%var(iLookBPAR%wallErosionRate),    & ! intent(in):    glacier wall erosion rate input for debris advection (mm yr-1)
                  bparData%var(iLookBPAR%debrisCritStress),   & ! intent(in):    critical driving stress where debris slides on terminal wedge (Pa)
                  bparData%var(iLookBPAR%latMoraineWidth),    & ! intent(in):    lateral moraine width or rockfall length (m)
                  ! area
                  bvarData%var(iLookBVAR%glacierAblArea)%dat, & ! intent(inout): per glacier ablation area (m2)
                  bvarData%var(iLookBVAR%glacierAccArea)%dat, & ! intent(inout): per glacier accumulation area (m2)
                  glac_area,                                  & ! intent(inout): area of each glacier domain (m2)
                  glac_ablFrac,                               & ! intent(out):   fraction of glacier domain that is ablation area
                  ! error handling
                  err, cmessage)                                ! intent(out):   error control
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

    ! update the glacier domains and their layers in each HRU
    iglacDOM = 0
    do iHRU=1,gruInfo%hruCount
      do iDOM = 1, gruInfo%hruInfo(iHRU)%domCount
        associate(domInfo => gruInfo%hruInfo(iHRU)%domInfo(iDOM))
          if(domInfo%dom_type/=glacCln1 .and. domInfo%dom_type/=glacCln2 .and. domInfo%dom_type/=glacDbr) cycle
          iglacDOM = iglacDOM + 1
          call updateGlacDomain(&
                      ! input
                      iglacDOM,                                  & ! intent(inout): glacier domain index
                      glac_elev,                                 & ! intent(in):    elevation of each glacier domain (m) per HRU
                      glac_area,                                 & ! intent(in):    area of each glacier domain (m2)
                      glac_tan_slope,                            & ! intent(in):    tan local ground surface slope of the domain (m/m)
                      glac_aspect,                               & ! intent(in):    azimuth in degrees East of North of the domain (degrees)
                      glac_contourLength,                        & ! intent(in):    length of contour at downslope edge of the domain (m)
                      glac_ablFrac,                              & ! intent(in):    fraction of glacier area that is ablation area
                      glac_debris_thick,                         & ! intent(in):    debris thickness of each glacier domain (m) per HRU
                      domInfo%dom_type,                          & ! intent(in):    domain type
                      domInfo%nSnow,                             & ! intent(in):    number of snow layers
                      domInfo%nLake,                             & ! intent(in):    number of lake layers
                      domInfo%nSoil,                             & ! intent(in):    number of soil layers
                      domInfo%nGlce,                             & ! intent(in):    number of glacier ice layers
                      ! data structures
                      mparHRU%hru(iHRU)%dom(iDOM),               & ! intent(in):    model parameters
                      indxHRU%hru(iHRU)%dom(iDOM),               & ! intent(in):    model indices
                      progHRU%hru(iHRU)%dom(iDOM),               & ! intent(inout): model prognostic variables
                      diagHRU%hru(iHRU)%dom(iDOM),               & ! intent(inout): model diagnostic variables
                      fluxHRU%hru(iHRU)%dom(iDOM),               & ! intent(inout): model fluxes
                      ! error handling
                      err, cmessage)                               ! intent(out):   error control
          if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
        end associate
      enddo ! (looping through domains)
    enddo ! (looping through HRUs)
    deallocate(glac_hru,glac_area,glac_elev,glac_tan_slope,glac_aspect,glac_contourLength,glac_debris_thick, &
               glac_ablFrac,massChange,iden_soil_mean,theta_sat_mean,nclean,ndebris)
  endif ! (if updateGlacArea)

  ! give the upland domain of each HRU the area, elevation, slope and aspect not taken by the other domains
  ! NOTE: contour length is not updated as we do not know how much of the HRU contour length belongs to the glacier/lake
  if(updateGlacArea .or. updateLakeArea)then
    do iHRU=1,gruInfo%hruCount
      associate(attr => attrHRU%hru(iHRU)%var)
        ! start from the HRU attributes and remove the area-weighted contribution of each non-upland domain
        remaining_area       = attr(iLookATTR%HRUarea)
        remaining_elev       = attr(iLookATTR%HRUarea)*attr(iLookATTR%elevation)
        remaining_tan_slope  = attr(iLookATTR%HRUarea)*attr(iLookATTR%tan_slope)
        remaining_aspect_sin = attr(iLookATTR%HRUarea)*sin(attr(iLookATTR%aspect)*deg2rad)
        remaining_aspect_cos = attr(iLookATTR%HRUarea)*cos(attr(iLookATTR%aspect)*deg2rad)
        iUpland = 0
        do iDOM = 1, gruInfo%hruInfo(iHRU)%domCount
          associate(typeDOM => gruInfo%hruInfo(iHRU)%domInfo(iDOM)%dom_type, &
                    DOMarea => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1), &
                    DOMelev => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMelev)%dat(1), &
                    DOMtan_slope => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMtan_slope)%dat(1), &
                    DOMaspect => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMaspect)%dat(1))
            if(typeDOM==upland)then
              iUpland = iDOM
            elseif(DOMarea>0._rkind)then
              remaining_area       = remaining_area       - DOMarea
              remaining_elev       = remaining_elev       - DOMarea*DOMelev
              remaining_tan_slope  = remaining_tan_slope  - DOMarea*DOMtan_slope
              remaining_aspect_sin = remaining_aspect_sin - DOMarea*sin(DOMaspect*deg2rad)
              remaining_aspect_cos = remaining_aspect_cos - DOMarea*cos(DOMaspect*deg2rad)
            endif
          end associate
        enddo
        if(iUpland==0) cycle
        associate(DOMarea => progHRU%hru(iHRU)%dom(iUpland)%var(iLookPROG%DOMarea)%dat(1), &
                  DOMelev => progHRU%hru(iHRU)%dom(iUpland)%var(iLookPROG%DOMelev)%dat(1), &
                  DOMtan_slope => progHRU%hru(iHRU)%dom(iUpland)%var(iLookPROG%DOMtan_slope)%dat(1), &
                  DOMaspect => progHRU%hru(iHRU)%dom(iUpland)%var(iLookPROG%DOMaspect)%dat(1), &
                  DOMcontourLength => progHRU%hru(iHRU)%dom(iUpland)%var(iLookPROG%DOMcontourLength)%dat(1))
          if(remaining_area>0._rkind)then
            ! the upland domain inherits the HRU attributes, re-derived by area weighting if other domains took part of the HRU
            DOMarea          = remaining_area
            DOMelev          = attr(iLookATTR%elevation)
            DOMtan_slope     = attr(iLookATTR%tan_slope)
            DOMaspect        = attr(iLookATTR%aspect)
            DOMcontourLength = attr(iLookATTR%contourLength) ! could be improved in the future
            if(remaining_area /= attr(iLookATTR%HRUarea))then
              DOMelev      = remaining_elev/remaining_area
              DOMtan_slope = remaining_tan_slope/remaining_area
              if(DOMaspect /= realMissing)then ! aspect is optional, realMissing when absent
                DOMaspect = 0._rkind
                if(remaining_aspect_sin**2 + remaining_aspect_cos**2 > aspect_tol) &
                  DOMaspect = modulo(atan2(remaining_aspect_sin,remaining_aspect_cos)*rad2deg,360._rkind)
              endif
            endif
          else
            DOMarea          = 0._rkind
            DOMelev          = realMissing
            DOMtan_slope     = realMissing
            DOMaspect        = realMissing
            DOMcontourLength = 0._rkind
          endif
        end associate
      end associate
    enddo ! (looping through HRUs)
  endif ! (if updated glacier or wetland area)
  call date_and_time(values=endUpdateArea)
  elapsedUpdateArea = elapsedUpdateArea + elapsedSec(startUpdateArea,endUpdateArea)

  deallocate(downIdx,inDegree,hruOrder)

end subroutine run_oneGRU

end module run_oneGRU_module

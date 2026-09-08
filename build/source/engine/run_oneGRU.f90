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

! physical constants
USE multiconst,only: iden_ice           ! intrinsic density of ice (kg m-3)

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
USE var_lookup,only:iLookINDEX         ! look-up values for local column index variables
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

! provide access to the named variables that describe model decisions
USE mDecisions_module,only:&           ! look-up values for the choice of method for the spatial representation of groundwater
 localColumn, &                        ! separate groundwater representation in each local soil column
 singleBasin, &                        ! single groundwater store over the entire basin
 bigBucket                             ! a big bucket (lumped aquifer model)
! -----------------------------------------------------------------------------------------------------------------------------------
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
  integer(i4b)                        :: iHRU                           ! HRU index
  integer(i4b)                        :: jHRU,kHRU                      ! index of the hydrologic response unit
  integer(i4b)                        :: iDOM                           ! domain index
  real(rkind)                         :: fracDOM                        ! fractional area of a given HRU domain in GRU (-)
  integer(i4b)                        :: nglacDOM                       ! number of glacier domains in the GRU
  integer(i4b)                        :: nglacHRU                       ! number of glacier HRUs in the GRU
  integer(i4b)                        :: iglacDOM                       ! glacier domain index
  integer(i4b)                        :: iglacHRU                       ! glacier HRU index
  real(rkind), allocatable            :: glac_elev(:)                   ! elevation of each glacier domain (m)
  real(rkind), allocatable            :: glac_tan_slope(:)              ! tan local ground surface slope of the domain (m/m)
  real(rkind), allocatable            :: glac_aspect(:)                 ! azimuth in degrees East of North of the domain (degrees)
  real(rkind), allocatable            :: glac_contourLength(:)          ! length of contour at downslope edge of the domain (m)
  real(rkind), allocatable            :: glac_debris_thick(:)           ! debris thickness of each glacier domain (m)
  real(rkind), allocatable            :: massChange(:)                  ! since last update mean rate glacier water equivalent change (kg m-2 s-1)
  integer(i8b), allocatable           :: glac_hru(:)                    ! HRU index of the each glacier cell
  real(rkind), allocatable            :: glac_ablFrac(:)                ! ablation fraction of each glacier domain
  real(rkind), allocatable            :: glac_area(:)                   ! area of each glacier domain (m2)
  real(rkind), allocatable            :: iden_soil_mean(:)              ! mean soil identity of each glacier domain
  real(rkind), allocatable            :: theta_sat_mean(:)              ! mean saturated water content of each glacier domain
  integer(i4b), allocatable           :: nclean(:)                      ! number of clean glacier HRUs in each glacier domain
  integer(i4b), allocatable           :: ndebris(:)                     ! number of debris glacier HRUs in each glacier domain
  logical(lgt)                        :: computeVegFluxFlag             ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
  logical(lgt)                        :: updateGlacArea                 ! flag to update glacier area
  logical(lgt)                        :: updateLakeArea                 ! flag to update wetland area
  logical(lgt)                        :: has_glacier                    ! flag to indicate if glaciers are present in HRU
  real(rkind)                         :: remaining_area                 ! remaining area to be distributed
  real(rkind)                         :: remaining_elev                 ! remaining elevation to be distributed
  real(rkind)                         :: remaining_tan_slope            ! remaining tan slope to be distributed (area-weighted)
  real(rkind)                         :: remaining_aspect_sin           ! remaining sine component for circular aspect mean
  real(rkind)                         :: remaining_aspect_cos           ! remaining cosine component for circular aspect mean
  logical(lgt)                        :: runHRU                         ! flag to run the HRU
  logical(lgt)                        :: check_updateGlacArea           ! flag to check if glacier area needs to be updated
  real(rkind)                         :: glacIceMelt                    ! glacier ice reservoir melt (m3 s-1)
  real(rkind)                         :: glacSnowMelt                   ! glacier snow reservoir melt (m3 s-1)
  real(rkind)                         :: glacFirnMelt                   ! glacier firn reservoir melt (m3 s-1)
  real(rkind)                         :: sec_since_last_update          ! seconds since last update
  real(rkind)                         :: soil_thick                     ! depth of soil== debris in debris domain of glacier HRU
  integer(i4b)                        :: nSnow                          ! number of snow layers in debris domain
  integer(i4b)                        :: nLake                          ! number of lake layers in debris domain (should be 0)
  integer(i4b)                        :: nSoil                          ! number of soil layers in debris domain
  real(rkind),parameter               :: deg2rad=PI_D/180._rkind        ! convert degrees to radians
  real(rkind),parameter               :: rad2deg=180._rkind/PI_D        ! convert radians to degrees
  real(rkind),parameter               :: aspect_tol=1.e-12_rkind        ! tolerance for undefined circular mean
  integer(i4b),dimension(8)           :: startUpdateArea, endUpdateArea ! time at end of updating glacier and wetland area
  ! ----------------------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; write(message, '(A21,I0,A10,I0,A2)' ) 'run_oneGRU (gru_nc = ',gruInfo%gru_nc,', gruId = ',gruInfo%gru_id,')/'

  ! ----- basin initialization --------------------------------------------------------------------------------------------
  ! initialize runoff variables
  bvarData%var(iLookBVAR%basin__SurfaceRunoff)%dat(1)    = 0._rkind  ! surface runoff (m s-1)
  bvarData%var(iLookBVAR%basin__SoilDrainage)%dat(1)     = 0._rkind  ! soil drainage (m s-1)
  bvarData%var(iLookBVAR%basin__ColumnOutflow)%dat(1)    = 0._rkind  ! outflow from all "outlet" HRUs (those with no downstream HRU)
  bvarData%var(iLookBVAR%basin__TotalRunoff)%dat(1)      = 0._rkind  ! total runoff to the channel from all active components (m s-1)

  ! initialize baseflow variables
  bvarData%var(iLookBVAR%basin__AquiferRecharge)%dat(1)  = 0._rkind ! recharge to the aquifer (m s-1)
  bvarData%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)  = 0._rkind ! baseflow from the aquifer (m s-1)
  bvarData%var(iLookBVAR%basin__AquiferTranspire)%dat(1) = 0._rkind ! transpiration loss from the aquifer (m s-1)

  ! initialize storage change variable
  bvarData%var(iLookBVAR%basin__StorageChange)%dat(1)    = 0._rkind ! change in total basin storage (kg m-2 s-1)

  ! initialize glacier variables
  glacIceMelt  = 0._rkind ! glacier ice reservoir melt (m3 s-1)
  glacSnowMelt = 0._rkind ! glacier snow reservoir melt (m3 s-1)
  glacFirnMelt = 0._rkind ! glacier firn reservoir melt (m3 s-1)
  updateGlacArea = .false. ! initialize flag to update glacier area
  updateLakeArea = .false. ! initialize flag to update wetland area
  bvarData%var(iLookBVAR%basin__GlacierArea)%dat(1) = 0._rkind ! glacier area (m2)
  nglacDOM = 0 ! initialize number of glacier domains in the GRU
  nglacHRU = 0 ! initialize number of glacier HRUs in the GRU

  ! initialize total inflow for each layer in a soil column and glacier size allocation
  check_updateGlacArea = .true.
  do iHRU=1,gruInfo%hruCount
    do iDOM = 1, gruInfo%hruInfo(iHRU)%domCount
      associate(typeDOM => gruInfo%hruInfo(iHRU)%domInfo(iDOM)%dom_type, &
                DOMarea => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1) )
        if(typeDOM==wetland)then; err=20; message=trim(message)//'ERROR:  wetland fluxes not yet implemented'; return; endif

        fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%mLayerColumnInflow)%dat(:) = 0._rkind
        if(DOMarea==0._rkind) cycle ! skip domains with no area
        if(typeDOM==glacCln1 .or. typeDOM==glacCln2 .or. typeDOM==glacDbr)then        
          if(check_updateGlacArea)then ! find updateJulDay, update glacier area every first of lowest mass month of mod year (choose October North Hemisphere, April South, January low latitudes)
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
            check_updateGlacArea = .false. ! only check this once for the GRU
          endif ! (checking when to update glacier area)
        endif ! (if glacier domain)
      end associate
    enddo ! (looping through domains)
  enddo ! (looping through HRUs)

  ! allocate space for glacier area change module variables
  if(updateGlacArea)then
    do iHRU=1,gruInfo%hruCount
      has_glacier = .false. ! initialize flag to indicate if glaciers are present in HRU
      do iDOM = 1, gruInfo%hruInfo(iHRU)%domCount
        associate(typeDOM => gruInfo%hruInfo(iHRU)%domInfo(iDOM)%dom_type)
          if(typeDOM==glacCln1 .or. typeDOM==glacCln2 .or. typeDOM==glacDbr)then
            nglacDOM = nglacDOM + 1
            if(.not.has_glacier)then
              has_glacier = .true. ! set flag to indicate if glaciers are present in HRU
              nglacHRU = nglacHRU + 1
            endif
          endif
        end associate
      enddo
    enddo
    allocate(glac_elev(nglacDOM),glac_debris_thick(nglacDOM),glac_area(nglacDOM),glac_ablFrac(nglacDOM),massChange(nglacDOM), &
             glac_hru(nglacDOM),iden_soil_mean(nglacDOM),theta_sat_mean(nglacDOM),nclean(nglacHRU),ndebris(nglacHRU), &
             glac_tan_slope(nglacDOM),glac_aspect(nglacDOM),glac_contourLength(nglacDOM))
  endif

  ! ********** RUN FOR ONE HRU ********************************************************************************************
  do iHRU=1,gruInfo%hruCount
    
    ! skip HRUs with no area
    runHRU = .false.
    do iDOM = 1, gruInfo%hruInfo(iHRU)%domCount
      if(progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1)>0._rkind) runHRU = .true.
    enddo
    if(.not. runHRU) cycle

    computeVegFluxFlag = (ixComputeVegFlux%hru(iHRU) == yes)  ! initialize the flag to compute the vegetation flux
    ! ----- run the model --------------------------------------------------------------------------------------------------

    ! simulation for a single HRU
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

    ! save the flag for computing the vegetation fluxes
    if(computeVegFluxFlag)       ixComputeVegFlux%hru(iHRU) = yes
    if(.not. computeVegFluxFlag) ixComputeVegFlux%hru(iHRU) = no

    ! ----- compute fluxes across HRUs --------------------------------------------------------------------------------------------------
    ! identify lateral connectivity
    ! (Note:  for efficiency, this could this be done as a setup task, not every timestep)
    kHRU = 0
    ! identify the downslope HRU
    dsHRU: do jHRU=1,gruInfo%hruCount
      if(typeHRU%hru(iHRU)%var(iLookTYPE%downHRUindex) == idHRU%hru(jHRU)%var(iLookID%hruId))then
        if(kHRU==0)then  ! check there is a unique match
          kHRU=jHRU
          exit dsHRU
        endif  ! (check there is a unique match)
      endif  ! (if identified a downslope HRU)
    enddo dsHRU
    
    do iDOM = 1, gruInfo%hruInfo(iHRU)%domCount
      if(progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1)==0._rkind) cycle ! skip domains with no area
      associate(typeDOM => gruInfo%hruInfo(iHRU)%domInfo(iDOM)%dom_type, &
                DOMarea => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1), &
                DOMelev => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMelev)%dat(1), &
                DOMtan_slope => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMtan_slope)%dat(1), &
                DOMaspect => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMaspect)%dat(1), &
                DOMcontourLength => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMcontourLength)%dat(1), &
                totalArea => bvarData%var(iLookBVAR%basin__totalArea)%dat(1) )

        ! identify the area covered by the current domain
        fracDOM = DOMarea / totalArea

        ! if lateral flows are active, add inflow to the downslope HRU
        if (typeDOM==upland)then
          if(kHRU > 0)then  ! if there is a downslope HRU, add to upland domain outflow to inflow (m3 s-1)
            fluxHRU%hru(kHRU)%dom(1)%var(iLookFLUX%mLayerColumnInflow)%dat(:) = fluxHRU%hru(kHRU)%dom(1)%var(iLookFLUX%mLayerColumnInflow)%dat(:)  + fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%mLayerColumnOutflow)%dat(:)
          else ! otherwise just increment basin (GRU) column outflow (m3 s-1) with the hru fraction
            bvarData%var(iLookBVAR%basin__ColumnOutflow)%dat(1) = bvarData%var(iLookBVAR%basin__ColumnOutflow)%dat(1) + sum(fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%mLayerColumnOutflow)%dat(:))
          endif
        endif ! (if upland domain)

        ! ----- calculate weighted basin (GRU) fluxes --------------------------------------------------------------------------------------
        bvarData%var(iLookBVAR%basin__StorageChange)%dat(1)  = bvarData%var(iLookBVAR%basin__StorageChange)%dat(1) + diagHRU%hru(iHRU)%dom(iDOM)%var(iLookDIAG%scalarTotalMassChange)%dat(1)*fracDOM
        if(typeDOM==upland)then
           ! increment basin surface runoff (m s-1)
          bvarData%var(iLookBVAR%basin__SurfaceRunoff)%dat(1) = bvarData%var(iLookBVAR%basin__SurfaceRunoff)%dat(1) + fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%scalarSurfaceRunoff)%dat(1)*fracDOM

          ! increment basin soil drainage (m s-1)
          bvarData%var(iLookBVAR%basin__SoilDrainage)%dat(1)  = bvarData%var(iLookBVAR%basin__SoilDrainage)%dat(1)  + fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%scalarSoilDrainage)%dat(1) *fracDOM

          ! increment aquifer variables -- ONLY if aquifer baseflow is computed individually for each HRU and aquifer is run
          ! NOTE: groundwater computed later for singleBasin
          ! NOTE: no groundwater for glacier
          if(model_decisions(iLookDECISIONS%spatial_gw)%iDecision == localColumn .and. model_decisions(iLookDECISIONS%groundwatr)%iDecision == bigBucket)then
            bvarData%var(iLookBVAR%basin__AquiferRecharge)%dat(1)  = bvarData%var(iLookBVAR%basin__AquiferRecharge)%dat(1)  + fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%scalarAquiferRecharge)%dat(1) *fracDOM
            bvarData%var(iLookBVAR%basin__AquiferTranspire)%dat(1) = bvarData%var(iLookBVAR%basin__AquiferTranspire)%dat(1) + fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%scalarAquiferTranspire)%dat(1)*fracDOM
            bvarData%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)  = bvarData%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)  + fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%scalarAquiferBaseflow)%dat(1) *fracDOM
          endif
        else if(typeDOM==glacCln1 .or. typeDOM==glacCln2 .or. typeDOM==glacDbr)then ! collect glacier ablation and accumulation melt m s-1
          ! This logic makes sense if assuming multiple glaciers in each HRU and one HRU per GRU, or one glacier in each GRU with multiple HRUs
          ! If some glaciers are not in a particular glacier HRU, this logic will not capture that
          glacFirnMelt = glacFirnMelt + fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%scalarGlacierMelt)%dat(1) *fracDOM * (1.0_rkind - progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarAblFrac)%dat(1)) ! no debris in accumulation zone for lateral flow
          if(progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSnowDepth)%dat(1)>0._rkind)then
            glacSnowMelt = glacSnowMelt + (fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%scalarGlacierMelt)%dat(1) + sum(fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%mLayerColumnOutflow)%dat(:))/totalArea) &
                          *fracDOM * progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarAblFrac)%dat(1)
          else
            glacIceMelt  = glacIceMelt  + (fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%scalarGlacierMelt)%dat(1) + sum(fluxHRU%hru(iHRU)%dom(iDOM)%var(iLookFLUX%mLayerColumnOutflow)%dat(:))/totalArea) &
                          *fracDOM * progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarAblFrac)%dat(1)
          endif
          ! increment basin glacier area (m2)
          bvarData%var(iLookBVAR%basin__GlacierArea)%dat(1) = bvarData%var(iLookBVAR%basin__GlacierArea)%dat(1) + DOMarea
          ! increment Gt of glacier storage (Gt = km3 of water equivalent)
          bvarData%var(iLookBVAR%basin__GlacierStorage)%dat(1) = bvarData%var(iLookBVAR%basin__GlacierStorage)%dat(1) &
                                                                 + diagHRU%hru(iHRU)%dom(iDOM)%var(iLookDIAG%scalarTotalMassChange)%dat(1)*data_step * DOMarea*1.e-12_rkind
        endif ! (if domain type)
      end associate
    enddo ! (looping through domains)

    ! averaging more fluxes (and/or states) can be added to this section as desired
  enddo  ! (looping through HRUs)

  ! if a year passed from last glacier area update, collect fluxes so that the glacier area can be updated
  if(updateGlacArea)then
    iglacHRU = 0 ! initialize number of glacier HRUs in the GRU
    iglacDOM = 0 ! initialize number of glacier domains in the GRU
    iden_soil_mean = 0._rkind ! initialize mean soil(debris) density of each glacier domain
    theta_sat_mean = 0._rkind ! initialize mean soil(debris) porosity of each glacier domain
    nclean = 0 ! initialize number of clean glacier domains in each HRU
    ndebris = 0 ! initialize number of debris glacier domains in each HRU
    do iHRU=1,gruInfo%hruCount
      has_glacier = .false. ! initialize flag to indicate if glaciers are present in HRU
      do iDOM = 1, gruInfo%hruInfo(iHRU)%domCount
        associate(typeDOM => gruInfo%hruInfo(iHRU)%domInfo(iDOM)%dom_type, &
                  DOMarea => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1), &
                  DOMelev => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMelev)%dat(1), &
                  DOMtan_slope => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMtan_slope)%dat(1), &
                  DOMaspect => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMaspect)%dat(1), &
                  DOMcontourLength => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMcontourLength)%dat(1), &
                  nSnow => gruInfo%hruInfo(iHRU)%domInfo(iDOM)%nSnow, &
                  nLake => gruInfo%hruInfo(iHRU)%domInfo(iDOM)%nLake, &
                  nSoil => gruInfo%hruInfo(iHRU)%domInfo(iDOM)%nSoil, &
                  mLayerDepth => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerDepth)%dat(:))
          if(typeDOM==glacCln1 .or. typeDOM==glacCln2 .or. typeDOM==glacDbr)then
            ! average layers mass change for each domain over time since last update
            iglacDOM = iglacDOM + 1
            glac_hru(iglacDOM) = iHRU
            if(.not.has_glacier)then 
              has_glacier = .true. ! set flag to indicate if glaciers are present in HRU
              iglacHRU = iglacHRU + 1 
            endif
            if(typeDOM==glacCln1 .or. typeDOM==glacCln2) nclean(iglacHRU) = nclean(iglacHRU) + 1
            if(typeDOM==glacDbr) ndebris(iglacHRU) = ndebris(iglacHRU) + 1
            if(DOMarea>0._rkind)then 
              glac_elev(iglacDOM) = DOMelev
              glac_area(iglacDOM) = DOMarea
              glac_tan_slope(iglacDOM) = DOMtan_slope
              glac_aspect(iglacDOM) = DOMaspect
              glac_contourLength(iglacDOM) = DOMcontourLength
              massChange(iglacDOM) = progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%glacMass4AreaChange)%dat(1)
              ! debris thickness is soil thickness in debris domain
              if(typeDOM==glacDbr)then
                soil_thick = sum(mLayerDepth(nSnow+nLake+1:nSnow+nLake+nSoil))
                iden_soil_mean(iglacDOM) = iden_soil_mean(iglacDOM) + sum(mparHRU%hru(iHRU)%dom(iDOM)%var(iLookPARAM%soil_dens_intr)%dat(1:nSoil) &
                                          *mLayerDepth(nSnow+nLake+1:nSnow+nLake+nSoil)) /soil_thick
                theta_sat_mean(iglacDOM) = theta_sat_mean(iglacDOM) + sum(mparHRU%hru(iHRU)%dom(iDOM)%var(iLookPARAM%theta_sat)%dat(1:nSoil) &
                                          *mLayerDepth(nSnow+nLake+1:nSnow+nLake+nSoil)) /soil_thick
                glac_debris_thick(iglacDOM) = soil_thick
              else
                glac_debris_thick(iglacDOM) = 0._rkind
              endif
            else ! fill in missing values
              glac_elev(iglacDOM) = realMissing
              glac_tan_slope(iglacDOM) = realMissing
              glac_aspect(iglacDOM) = realMissing
              glac_contourLength(iglacDOM) = 0._rkind
              glac_area(iglacDOM) = 0._rkind
              massChange(iglacDOM) = 0._rkind
              glac_debris_thick(iglacDOM) = 0._rkind
              iden_soil_mean(iglacDOM) = 0._rkind
              theta_sat_mean(iglacDOM) = 0._rkind
            endif ! (if domain has area)
          endif ! (if glacier domain)
        end associate
      enddo ! (looping through domains)
    enddo ! (looping through HRUs)
  endif ! (if need to update glacier area)

  ! ********** END LOOP THROUGH HRUS **************************************************************************************
  ! lapse glacier fluxes to the basin by routing through each glacier
  call qGlacier(&
                ! input
                bparData%var(iLookBPAR%glacStor_kIce),              & ! intent(in):    storage coefficient ice reservoir (hours)
                bparData%var(iLookBPAR%glacStor_kFirn),             & ! intent(in):    storage coefficient snow reservoir (hours)
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
 
  ! perform the routing
  associate(totalArea => bvarData%var(iLookBVAR%basin__totalArea)%dat(1) )
  
    ! compute water balance for the basin aquifer
    if(model_decisions(iLookDECISIONS%spatial_gw)%iDecision == singleBasin)then
      message=trim(message)//'multi_driver/bigBucket groundwater code not transferred from old code base yet'
      err=20; return
    endif

    ! calculate total runoff depending on whether aquifer is connected
    if(model_decisions(iLookDECISIONS%groundwatr)%iDecision == bigBucket)then
      ! deep aquifer (column outflow will be zero)
      bvarData%var(iLookBVAR%basin__TotalRunoff)%dat(1) = bvarData%var(iLookBVAR%basin__SurfaceRunoff)%dat(1) + bvarData%var(iLookBVAR%basin__ColumnOutflow)%dat(1)/totalArea + bvarData%var(iLookBVAR%basin__AquiferBaseflow)%dat(1)
    else
      ! no deep aquifer (may have column outflow from shallow groundwater then soil drainage will be zero, else the converse is true)
      bvarData%var(iLookBVAR%basin__TotalRunoff)%dat(1) = bvarData%var(iLookBVAR%basin__SurfaceRunoff)%dat(1) + bvarData%var(iLookBVAR%basin__ColumnOutflow)%dat(1)/totalArea + bvarData%var(iLookBVAR%basin__SoilDrainage)%dat(1)
    endif
    
    call qOverland(&
                   ! input
                   model_decisions(iLookDECISIONS%subRouting)%iDecision, & ! intent(in):    index for routing method
                   bvarData%var(iLookBVAR%basin__TotalRunoff)%dat(1),    & ! intent(in):    total runoff to the channel from all active components (m s-1)
                   bvarData%var(iLookBVAR%routingFractionFuture)%dat,    & ! intent(in):    fraction of runoff in future time steps (m s-1)
                   bvarData%var(iLookBVAR%routingRunoffFuture)%dat,      & ! intent(inout): runoff in future time steps (m s-1)
                   ! output
                   bvarData%var(iLookBVAR%averageInstantRunoff)%dat(1),  & ! intent(out):   instantaneous runoff (m s-1)
                   bvarData%var(iLookBVAR%averageRoutedRunoff)%dat(1),   & ! intent(out):   routed runoff (m s-1)
                   err,cmessage)                                           ! intent(out):   error control
    if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

    ! add glacier runoff to overland runoff
    bvarData%var(iLookBVAR%averageInstantRunoff)%dat(1) = bvarData%var(iLookBVAR%averageInstantRunoff)%dat(1) + glacIceMelt + glacSnowMelt + glacFirnMelt
    bvarData%var(iLookBVAR%averageRoutedRunoff)%dat(1) = bvarData%var(iLookBVAR%averageRoutedRunoff)%dat(1) + bvarData%var(iLookBVAR%glacierRoutedRunoff)%dat(1)

  end associate

  call date_and_time(values=startUpdateArea)
  ! Need to update the glacier area
  if(updateGlacArea)then
    ! need to save length, bottom topo, and elevation of glaciers from the end of previous update for this GRU in file associated with gruInfo%gru_id
    ! need to associate each glacier with an HRU and domain
    call glacAreaChange(&
                  ! model control
                  sec_since_last_update,                      & ! intent(in):    seconds since last glacier area update
                  nglacHRU,                                   & ! intent(in):    number of HRUs that have a glacier domain
                  nglacDOM,                                   & ! intent(in):    number of domains that have glaciers
                  ndebris,                                    & ! intent(in):    number of debris domains in each HRU
                  nclean,                                     & ! intent(in):    number of clean domains in each HRU
                  glac_hru,                                   & ! intent(in):    HRU index of glacier domain
                  ! glacier topography
                  gruInfo%nGlac,                              & ! intent(inout): number of glaciers in GRU
                  gruInfo%glacInfo,                           & ! intent(inout): information for each glacier
                  gruInfo%gridInfo,                           & ! intent(in):    grid information for each grid
                  gridData,                                   & ! intent(inout): grid data for each grid
                  ! mass balance per glacier domain
                  massChange,                                 & ! intent(in):    since updateJulDay rate glacier water equivalent change (kg m-2 s-1)
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

    ! update the glacier domains and domain layers in each HRU
    iglacDOM = 0 ! initialize glacier domain index
    do iHRU=1,gruInfo%hruCount
      do iDOM = 1, gruInfo%hruInfo(iHRU)%domCount
        associate(typeDOM => gruInfo%hruInfo(iHRU)%domInfo(iDOM)%dom_type)
          if(typeDOM==glacCln1 .or. typeDOM==glacCln2 .or. typeDOM==glacDbr)then
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
                        typeDOM,                                   & ! intent(in):    domain type
                        gruInfo%hruInfo(iHRU)%domInfo(iDOM)%nSnow, & ! intent(in):    number of snow layers
                        gruInfo%hruInfo(iHRU)%domInfo(iDOM)%nLake, & ! intent(in):    number of lake layers
                        gruInfo%hruInfo(iHRU)%domInfo(iDOM)%nSoil, & ! intent(in):    number of soil layers
                        gruInfo%hruInfo(iHRU)%domInfo(iDOM)%nGlce, & ! intent(in):    number of glacier ice layers
                        ! data structures
                        mparHRU%hru(iHRU)%dom(iDOM),               & ! intent(in):    model parameters
                        indxHRU%hru(iHRU)%dom(iDOM),               & ! intent(in):    model indices
                        progHRU%hru(iHRU)%dom(iDOM),               & ! intent(inout): model prognostic variables
                        diagHRU%hru(iHRU)%dom(iDOM),               & ! intent(inout): model diagnostic variables
                        fluxHRU%hru(iHRU)%dom(iDOM),               & ! intent(inout): model fluxes
                        ! error handling
                        err, cmessage)                               ! intent(out):   error control
            if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
          endif ! (if glacier domain)
        end associate
      enddo ! (looping through domains)
    enddo ! (looping through HRUs)
    deallocate(glac_elev,glac_debris_thick,glac_area,glac_ablFrac,massChange,glac_hru,iden_soil_mean, &
               theta_sat_mean,nclean,ndebris,glac_tan_slope,glac_aspect,glac_contourLength)
  endif ! (if updateGlacArea)

  if(updateGlacArea .or. updateLakeArea)then
    do iHRU=1,gruInfo%hruCount
      ! update the upland coordinate variables for the HRU based on the new glacier and wetland areas
      ! NOTE: contour length is not updated as we do not know how much of the original contour length is associated with the glacier/lake 
      remaining_area = attrHRU%hru(iHRU)%var(iLookATTR%HRUarea)
      remaining_elev = attrHRU%hru(iHRU)%var(iLookATTR%HRUarea)*attrHRU%hru(iHRU)%var(iLookATTR%elevation)
      remaining_tan_slope = attrHRU%hru(iHRU)%var(iLookATTR%HRUarea)*attrHRU%hru(iHRU)%var(iLookATTR%tan_slope)
      remaining_aspect_sin = attrHRU%hru(iHRU)%var(iLookATTR%HRUarea)*sin(attrHRU%hru(iHRU)%var(iLookATTR%aspect)*deg2rad)
      remaining_aspect_cos = attrHRU%hru(iHRU)%var(iLookATTR%HRUarea)*cos(attrHRU%hru(iHRU)%var(iLookATTR%aspect)*deg2rad)
      do iDOM = 1, gruInfo%hruInfo(iHRU)%domCount
        associate(typeDOM => gruInfo%hruInfo(iHRU)%domInfo(iDOM)%dom_type, &
                  DOMarea => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1), &
                  DOMelev => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMelev)%dat(1), &
                  DOMtan_slope => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMtan_slope)%dat(1), &
                  DOMaspect => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMaspect)%dat(1) )
          if(typeDOM.ne.upland .and. DOMarea>0._rkind)then
            remaining_area = remaining_area - DOMarea
            remaining_elev = remaining_elev - DOMarea * DOMelev
            remaining_tan_slope = remaining_tan_slope - DOMarea * DOMtan_slope
            remaining_aspect_sin = remaining_aspect_sin - DOMarea*sin(DOMaspect*deg2rad)
            remaining_aspect_cos = remaining_aspect_cos - DOMarea*cos(DOMaspect*deg2rad)
          endif
        end associate
      enddo
      do iDOM = 1, gruInfo%hruInfo(iHRU)%domCount
        associate(typeDOM => gruInfo%hruInfo(iHRU)%domInfo(iDOM)%dom_type, &
                  DOMarea => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMarea)%dat(1), &
                  DOMelev => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMelev)%dat(1), &
                  DOMtan_slope => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMtan_slope)%dat(1), &
                  DOMaspect => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMaspect)%dat(1), &
                  DOMcontourLength => progHRU%hru(iHRU)%dom(iDOM)%var(iLookPROG%DOMcontourLength)%dat(1) )
          if(typeDOM==upland)then
            DOMarea = remaining_area
            if(remaining_area>0._rkind)then 
              DOMelev = remaining_elev/remaining_area
              DOMtan_slope = remaining_tan_slope/remaining_area
              if(remaining_aspect_sin**2 + remaining_aspect_cos**2 > aspect_tol)then
                DOMaspect = modulo(atan2(remaining_aspect_sin,remaining_aspect_cos)*rad2deg,360._rkind)
              else
                DOMaspect = 0._rkind
              endif
              DOMcontourLength = attrHRU%hru(iHRU)%var(iLookATTR%contourLength) ! for now, just keep upchangint at the HRU contour length, but could be improved in the future
            else
              DOMelev = realMissing
              DOMarea = 0._rkind
              DOMtan_slope = realMissing
              DOMaspect = realMissing
              DOMcontourLength = 0._rkind
            endif
          endif ! (if upland domain)
        end associate
      enddo ! (looping through domains)
    enddo ! (looping through HRUs)
  endif ! (if updated glacier or wetland area)
  call date_and_time(values=endUpdateArea)

  ! aggregate the elapsed time for the update area routines
   elapsedUpdateArea = elapsedUpdateArea + elapsedSec(startUpdateArea,endUpdateArea)

end subroutine run_oneGRU

end module run_oneGRU_module

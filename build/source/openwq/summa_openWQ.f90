! Copyright 2020, Diogo Costa (diogo.pinhodacosta@canada.ca)
! This file is part of OpenWQ model.

! This program, openWQ, is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! ==============================================================================
! Coupling of OpenWQ to SUMMA and to the internally coupled mizuRoute
! ==============================================================================
! One OpenWQ instance carries the solutes of the land stores of SUMMA and, when
! mizuRoute is active, of the river reaches. Every water flux of the two models
! that moves liquid water between stores has a solute counterpart here.
!
! Land: one cell column per (HRU, domain) pair. A domain is a stack of
!   canopy -> snow layers -> lake layers -> soil layers -> (glacier ice) -> aquifer
! of which only the stores that exist in the domain are used.
!
!   precipitation     -> canopy, and -> top store (snow, else lake, else RUNOFF)
!   canopy            -> top store (drainage and unloading)
!   snow layers       -> next snow layer, base -> lake top, else RUNOFF
!   lake layers       -> next lake layer, base -> RUNOFF
!   RUNOFF            -> first soil layer (infiltration), -> stream (surface runoff)
!   soil layers       -> next soil layer (either direction)
!   soil layers       -> same layer of the downslope HRU, else stream (lateral outflow)
!   last soil layer   -> aquifer, else stream (drainage)
!   aquifer           -> stream (baseflow)
!   glacier ice melt  -> stream (external water, EWF "GLACIER_ICE_MELT")
!
! RUNOFF is the liquid water that reaches the surface in the step (rain plus melt).
! Evaporation, transpiration and sublimation carry no solute.
!
! "Stream" is the RUNOFF_TO_STREAM pool of the column, which holds the solute
! while SUMMA delays the runoff of the GRU (routing through the unresolved
! network and the glacier reservoirs). It releases the same fraction of its
! content as SUMMA releases of the water it holds:
!   with mizuRoute    -> the reaches that receive the runoff of the GRU
!   without mizuRoute -> out of the domain
! Wetland domains are not connected to the basin runoff in SUMMA, so what
! leaves them leaves the domain.
!
! River (mizuRoute): one cell per reach.
!   reach -> downstream reach, or out of the domain at an outlet
! with the reach mixing volume = storage + upstream inflow + lateral inflow.
!
! Compartment, export and dependency indices are those of OpenWQ_hydrolink.h.
! ==============================================================================

module summa_openwq

  USE nr_type
  USE, intrinsic :: iso_c_binding, only: c_long_long, c_int
  USE openWQ, only: CLASSWQ_openwq
  USE build_options, only: mizuroute_active

  implicit none
  private

  public :: openwq_init
  public :: openwq_run_time_start
  public :: openwq_run_space_step
  public :: openwq_run_time_end
  public :: openwq_finalize

  type(CLASSWQ_openwq), save, public :: openwq_obj

  ! compartment indices (0-based)
  integer(i4b), parameter :: canopy_cmp  = 0
  integer(i4b), parameter :: snow_cmp    = 1
  integer(i4b), parameter :: runoff_cmp  = 2
  integer(i4b), parameter :: soil_cmp    = 3
  integer(i4b), parameter :: aquifer_cmp = 4
  integer(i4b), parameter :: stream_cmp  = 5
  integer(i4b), save      :: lake_cmp    = -1
  integer(i4b), save      :: river_cmp   = -1
  integer(i4b), parameter :: out_cmp     = -1   ! out of the domain

  ! flux-concentration exports (0-based)
  integer(i4b), parameter :: runoffVol_exp     = 0
  integer(i4b), parameter :: routedRunoff_exp  = 1
  integer(i4b), parameter :: totalRunoff_exp   = 2
  integer(i4b), parameter :: reachOutflow_exp  = 3

  ! a pool with less than this depth of water over the column is reported as holding no water (m)
  real(rkind), parameter  :: minPoolDepth = 1.e-9_rkind

  ! land columns
  integer(i4b), save              :: nCol = 0            ! number of (HRU, domain) columns
  integer(i4b), save              :: nSnowMax = 0        ! snow cells per column
  integer(i4b), save              :: nSoilMax = 0        ! soil cells per column
  integer(i4b), save              :: nLakeMax = 0        ! lake cells per column
  integer(i4b), save, allocatable :: gruFirstCol(:)      ! first column of each GRU
  integer(i4b), save, allocatable :: gruLastCol(:)       ! last column of each GRU
  integer(i4b), save, allocatable :: colDown(:)          ! column of the downslope HRU (0 = none)

  ! state at the start of the time step
  logical(lgt), save, allocatable :: colActive(:)        ! column is simulated by SUMMA
  real(rkind),  save, allocatable :: colArea(:)          ! area (m2)
  integer(i4b), save, allocatable :: nSnowStart(:)       ! number of snow layers
  real(rkind),  save, allocatable :: volCanopy(:)        ! water volumes (m3)
  real(rkind),  save, allocatable :: volAquifer(:)
  real(rkind),  save, allocatable :: volSnow(:,:)
  real(rkind),  save, allocatable :: volLake(:,:)
  real(rkind),  save, allocatable :: volSoil(:,:)

  ! water held by SUMMA between its release by the columns and its delivery to the stream (m3)
  real(rkind),  save, allocatable :: gruRoutingStore(:)

  ! river network
  logical(lgt), save              :: riverCoupled = .false.
  integer(i4b), save              :: nReach = 0

  ! running water check of the coupling (m3), reported by openwq_finalize
  real(rkind),  save              :: cumRoutedLand = 0._rkind   ! runoff routed by SUMMA out of its GRUs
  real(rkind),  save              :: cumLateral    = 0._rkind   ! lateral inflow received by the reaches
  real(rkind),  save              :: cumOutflow    = 0._rkind   ! outflow of the outlet reaches
  real(rkind),  save              :: cumResidual   = 0._rkind   ! |start + inflows - outflows - end| summed over reaches and steps
  real(rkind),  save              :: cumStorageEnd = 0._rkind   ! reach storage at the last step

contains

  ! ============================================================================
  ! Declare the OpenWQ domain: land columns and, with mizuRoute, river reaches
  ! ============================================================================
  subroutine openwq_init(summa1_struc, err)
    USE summa_type, only: summa1_type_dec
    USE globalData, only: gru_struc
    USE globalData, only: maxSnowLayers, maxSoilLayers, maxLakeLayers
    USE var_lookup, only: iLookTYPE, iLookID
#ifdef MIZUROUTE_ACTIVE
    USE mizuroute_coupling, only: init_mizuroute_wq_from_summa
#endif
    implicit none
    type(summa1_type_dec), intent(inout) :: summa1_struc
    integer(i4b),          intent(out)   :: err
    ! local variables
    integer(i4b)                      :: iGRU, iHRU, jHRU, iDOM, iCol
    integer(i4b)                      :: hasGlacier
    integer(i4b), allocatable         :: hruFirstCol(:)
    integer(c_long_long), allocatable :: colId(:), reachId(:)
    integer(c_int), allocatable       :: colDom(:)
    character(len=256)                :: message

    err = 0
    ! one OpenWQ instance per executable
    if(allocated(gruFirstCol))then
      write(*,'(a)') 'openwq_init/OpenWQ is already initialized: only one SUMMA instance can be coupled to OpenWQ'
      err = 20; return
    endif
    openwq_obj = CLASSWQ_openwq()

    ! ----- land columns, in the order of the SUMMA loops (GRU, HRU, domain) -----
    nCol = 0
    do iGRU = 1, summa1_struc%nGRU_local
      do iHRU = 1, gru_struc(iGRU)%hruCount
        nCol = nCol + gru_struc(iGRU)%hruInfo(iHRU)%domCount
      end do
    end do
    nSnowMax = max(maxSnowLayers, 1)
    nSoilMax = max(maxSoilLayers, 1)
    nLakeMax = max(maxLakeLayers, 0)

    allocate(gruFirstCol(summa1_struc%nGRU_local), gruLastCol(summa1_struc%nGRU_local), colDown(nCol), &
             colId(nCol), colDom(nCol), colActive(nCol), colArea(nCol), nSnowStart(nCol),              &
             volCanopy(nCol), volAquifer(nCol), volSnow(nSnowMax,nCol), volLake(max(nLakeMax,1),nCol), &
             volSoil(nSoilMax,nCol), gruRoutingStore(summa1_struc%nGRU_local), stat=err)
    if(err/=0) return
    colDown = 0; colActive = .false.; colArea = 0._rkind; nSnowStart = 0
    volCanopy = 0._rkind; volAquifer = 0._rkind; volSnow = 0._rkind; volLake = 0._rkind; volSoil = 0._rkind
    gruRoutingStore = 0._rkind

    hasGlacier = 0
    iCol = 0
    do iGRU = 1, summa1_struc%nGRU_local
      allocate(hruFirstCol(gru_struc(iGRU)%hruCount))
      gruFirstCol(iGRU) = iCol + 1
      do iHRU = 1, gru_struc(iGRU)%hruCount
        hruFirstCol(iHRU) = iCol + 1
        do iDOM = 1, gru_struc(iGRU)%hruInfo(iHRU)%domCount
          iCol = iCol + 1
          colId(iCol)  = gru_struc(iGRU)%hruInfo(iHRU)%hru_id
          colDom(iCol) = merge(iDOM, 0, gru_struc(iGRU)%hruInfo(iHRU)%domCount > 1)
          if(gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nGlce > 0) hasGlacier = 1
        end do
      end do
      gruLastCol(iGRU) = iCol
      ! lateral soil outflow goes to the first (upland) domain of the downslope HRU
      do iHRU = 1, gru_struc(iGRU)%hruCount
        do jHRU = 1, gru_struc(iGRU)%hruCount
          if(summa1_struc%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%downHRUindex) == &
             summa1_struc%idStruct%gru(iGRU)%hru(jHRU)%var(iLookID%hruId))then
            colDown(hruFirstCol(iHRU)) = hruFirstCol(jHRU)
            exit
          endif
        end do
      end do
      deallocate(hruFirstCol)
    end do

    ! ----- river reaches -----
    riverCoupled = .false.
    nReach = 0
#ifdef MIZUROUTE_ACTIVE
    if(mizuroute_active)then
      if(summa1_struc%config%use_mizuroute)then
        call init_mizuroute_wq_from_summa(summa1_struc, err, message)
        if(err/=0)then; write(*,'(a)') 'openwq_init/'//trim(message); return; endif
        riverCoupled = .true.
        nReach = size(summa1_struc%mizu_domain%river_network%driver%seg_id)
        allocate(reachId(nReach))
        reachId(:) = int(summa1_struc%mizu_domain%river_network%driver%seg_id(:), c_long_long)
      endif
    endif
#endif
    if(.not.allocated(reachId)) allocate(reachId(1), source=0_c_long_long)

    ! compartments that only exist in some configurations follow the fixed ones
    lake_cmp  = -1
    river_cmp = -1
    if(nLakeMax > 0) lake_cmp = stream_cmp + 1
    if(riverCoupled) river_cmp = max(stream_cmp, lake_cmp) + 1

    err = openwq_obj%decl( &
      nCol,                & ! land columns
      1,                   & ! canopy layers
      nSnowMax,            & ! snow layers
      nSoilMax,            & ! soil layers
      1,                   & ! runoff layers
      1,                   & ! aquifer layers
      nLakeMax,            & ! lake layers (0 = no lake compartment)
      1,                   & ! cells in the y direction
      colId,               & ! HRU id of each column
      colDom,              & ! domain of each column (0 if the HRU has one domain)
      hasGlacier,          & ! 1 if any domain has glacier ice
      nReach,              & ! river reaches (0 = no river compartment)
      reachId)               ! reach ids

    if(riverCoupled) call check_river_mapping(summa1_struc)

  end subroutine openwq_init


  ! ============================================================================
  ! Compare the delivery map with the land it comes from, once at initialization
  ! ============================================================================
  ! The area a GRU drains through the map (its reach inflow per unit runoff) must
  ! equal its SUMMA area, or the reaches receive more or less water than the
  ! land delivers and the stream concentrations are off by that ratio.
  subroutine check_river_mapping(summa1_struc)
    USE summa_type, only: summa1_type_dec
    USE var_lookup, only: iLookBVAR
    implicit none
    type(summa1_type_dec), intent(in) :: summa1_struc
#ifdef MIZUROUTE_ACTIVE
    integer(i4b) :: iGRU, nBad, nMethods
    real(rkind)  :: areaMap, areaGRU, worst

    associate(wq => summa1_struc%mizu_domain%river_network%driver%wq)
    nBad = 0; worst = 0._rkind
    do iGRU = 1, summa1_struc%nGRU_local
      areaGRU = summa1_struc%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1)
      areaMap = sum(wq%map_area(wq%map_start(iGRU):wq%map_start(iGRU+1)-1))
      if(areaGRU <= 0._rkind .or. areaMap <= 0._rkind) cycle
      worst = max(worst, abs(areaMap/areaGRU - 1._rkind))
      if(abs(areaMap/areaGRU - 1._rkind) > 0.01_rkind) nBad = nBad + 1
    end do
    nMethods = size(summa1_struc%mizu_domain%river_network%driver%method)
    write(*,'(a,i0,a,i0,a)') ' OpenWQ river coupling: ', summa1_struc%nGRU_local, ' GRU(s) deliver to ', nReach, ' reach(es)'
    if(nBad > 0)then
      write(*,'(a,i0,a,f6.1,a)') ' WARNING: OpenWQ river coupling: the river-network area of ', nBad, &
        ' GRU(s) differs from the SUMMA area by more than 1% (worst ', 100._rkind*worst, &
        '%); the reaches then receive a different volume of water than the land delivers'
    endif
    if(nMethods > 1) write(*,'(a)') ' OpenWQ river coupling: several routing methods are active; water quality follows the first one'
    end associate
#endif
  end subroutine check_river_mapping


  ! ============================================================================
  ! Pass the water volumes and dependencies at the start of the time step
  ! ============================================================================
  subroutine openwq_run_time_start(summa1_struc)
    USE summa_type, only: summa1_type_dec
    USE globalData, only: gru_struc
    USE globalData, only: realMissing
    USE globalData, only: model_decisions
    USE mDecisions_module, only: bigBucket, singleBasin
    USE var_lookup, only: iLookPROG, iLookINDEX, iLookTYPE, iLookFORCE, iLookBVAR, iLookDECISIONS
    USE multiconst, only: iden_water
    USE module_sf_noahmplsm, only: isWater
    implicit none
    type(summa1_type_dec), intent(in) :: summa1_struc
    ! local variables
    integer(i4b)             :: iGRU, iHRU, iDOM, iCol, iLayer
    integer(i4b)             :: nSnow, nLake, nSoil
    integer(i4b)             :: simtime(5)
    integer(i4b)             :: err
    real(rkind)              :: airTemp_K, SWrad_Wm2
    real(rkind)              :: soilTemp_K(nSoilMax), soilMoist(nSoilMax)
    logical(lgt)             :: basinAquifer

    call get_simtime(summa1_struc, simtime)

    if(riverCoupled) call set_reach_state(summa1_struc)

    basinAquifer = (model_decisions(iLookDECISIONS%groundwatr)%iDecision == bigBucket .and. &
                    model_decisions(iLookDECISIONS%spatial_gw)%iDecision == singleBasin)

    iCol = 0
    do iGRU = 1, summa1_struc%nGRU_local
      do iHRU = 1, gru_struc(iGRU)%hruCount
        do iDOM = 1, gru_struc(iGRU)%hruInfo(iHRU)%domCount
          iCol = iCol + 1

          associate( &
            prog    => summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%dom(iDOM)%var, &
            indx    => summa1_struc%indxStruct%gru(iGRU)%hru(iHRU)%dom(iDOM)%var, &
            forc    => summa1_struc%forcStruct%gru(iGRU)%hru(iHRU)%var,           &
            domInfo => gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM))

          colArea(iCol)   = prog(iLookPROG%DOMarea)%dat(1)
          colActive(iCol) = (summa1_struc%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex) /= isWater) &
                            .and. (colArea(iCol) > 0._rkind) .and. (colArea(iCol) /= realMissing)
          if(.not.colActive(iCol)) colArea(iCol) = 0._rkind

          nSnow = 0; nLake = 0; nSoil = 0
          volCanopy(iCol)  = 0._rkind
          volAquifer(iCol) = 0._rkind
          volSnow(:,iCol)  = 0._rkind
          volLake(:,iCol)  = 0._rkind
          volSoil(:,iCol)  = 0._rkind
          soilTemp_K(:)    = 0._rkind
          soilMoist(:)     = 0._rkind
          airTemp_K        = forc(iLookFORCE%airtemp)
          SWrad_Wm2        = forc(iLookFORCE%SWRadAtm)
          if(SWrad_Wm2 == realMissing) SWrad_Wm2 = 0._rkind

          if(colActive(iCol))then

            nSnow = min(max(0, indx(iLookINDEX%nSnow)%dat(1)), nSnowMax)
            nLake = min(max(0, domInfo%nLake), nLakeMax)
            nSoil = min(max(0, domInfo%nSoil), nSoilMax)

            ! canopy (kg m-2) and aquifer (m)
            if(prog(iLookPROG%scalarCanopyWat)%dat(1) /= realMissing) &
              volCanopy(iCol) = prog(iLookPROG%scalarCanopyWat)%dat(1) * colArea(iCol) / iden_water
            if(prog(iLookPROG%scalarAquiferStorage)%dat(1) /= realMissing .and. .not.basinAquifer) &
              volAquifer(iCol) = prog(iLookPROG%scalarAquiferStorage)%dat(1) * colArea(iCol)

            ! layers are ordered snow, lake, soil, glacier ice
            do iLayer = 1, nSnow
              volSnow(iLayer,iCol) = layer_volume(iLayer)
            end do
            do iLayer = 1, nLake
              volLake(iLayer,iCol) = layer_volume(nSnow + iLayer)
            end do
            do iLayer = 1, nSoil
              volSoil(iLayer,iCol)  = layer_volume(nSnow + nLake + iLayer)
              soilTemp_K(iLayer)    = prog(iLookPROG%mLayerTemp)%dat(nSnow + nLake + iLayer)
              soilMoist(iLayer)     = prog(iLookPROG%mLayerVolFracLiq)%dat(nSnow + nLake + iLayer)
            end do

          endif

          ! a basin-wide aquifer is carried by the first column of the GRU
          if(basinAquifer .and. iCol == gruFirstCol(iGRU)) &
            volAquifer(iCol) = summa1_struc%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferStorage)%dat(1) &
                               * summa1_struc%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1)

          nSnowStart(iCol) = nSnow

          err = openwq_obj%openwq_run_time_start( &
            iCol - 1,                 & ! 0-based column index
            nSnow,                    &
            nLake,                    &
            nSoil,                    &
            simtime,                  &
            soilMoist,                & ! volumetric liquid water content of the soil layers (-)
            soilTemp_K,               & ! soil temperature (K)
            airTemp_K,                & ! air temperature (K)
            SWrad_Wm2,                & ! incoming shortwave radiation (W m-2)
            volSnow(:,iCol),          &
            volLake(:,iCol),          &
            volCanopy(iCol),          &
            volSoil(:,iCol),          &
            volAquifer(iCol),         &
            colArea(iCol))

          end associate

        end do
      end do
    end do

  contains

    ! total water (liquid + ice) of a layer (m3)
    real(rkind) function layer_volume(ix)
      integer(i4b), intent(in) :: ix
      layer_volume = 0._rkind
      associate(prog => summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%dom(iDOM)%var)
        if(prog(iLookPROG%mLayerVolFracWat)%dat(ix) /= realMissing) &
          layer_volume = prog(iLookPROG%mLayerVolFracWat)%dat(ix) * prog(iLookPROG%mLayerDepth)%dat(ix) * colArea(iCol)
      end associate
    end function layer_volume

  end subroutine openwq_run_time_start


  ! ============================================================================
  ! Reach storage and dependencies at the start of the time step
  ! ============================================================================
  ! Air temperature and radiation of a reach are the area-weighted values of the GRUs draining to it.
  subroutine set_reach_state(summa1_struc)
    USE summa_type, only: summa1_type_dec
    USE globalData, only: gru_struc
    USE globalData, only: realMissing
    USE var_lookup, only: iLookFORCE, iLookATTR
    implicit none
    type(summa1_type_dec), intent(in) :: summa1_struc
#ifdef MIZUROUTE_ACTIVE
    integer(i4b)             :: iGRU, iHRU, iMap, iReach, err
    real(rkind)              :: areaGRU, tairGRU, swGRU, areaHRU, sw
    real(rkind)              :: reachVol(nReach), reachTair(nReach), reachSW(nReach), reachArea(nReach)

    associate(wq => summa1_struc%mizu_domain%river_network%driver%wq)

    reachVol(:)  = max(wq%vol_end(:), 0._rkind)
    reachTair(:) = 0._rkind
    reachSW(:)   = 0._rkind
    reachArea(:) = 0._rkind

    do iGRU = 1, summa1_struc%nGRU_local
      areaGRU = 0._rkind; tairGRU = 0._rkind; swGRU = 0._rkind
      do iHRU = 1, gru_struc(iGRU)%hruCount
        areaHRU = summa1_struc%attrStruct%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea)
        sw      = summa1_struc%forcStruct%gru(iGRU)%hru(iHRU)%var(iLookFORCE%SWRadAtm)
        if(sw == realMissing) sw = 0._rkind
        areaGRU = areaGRU + areaHRU
        tairGRU = tairGRU + areaHRU * summa1_struc%forcStruct%gru(iGRU)%hru(iHRU)%var(iLookFORCE%airtemp)
        swGRU   = swGRU   + areaHRU * sw
      end do
      if(areaGRU <= 0._rkind) cycle
      tairGRU = tairGRU / areaGRU
      swGRU   = swGRU   / areaGRU
      do iMap = wq%map_start(iGRU), wq%map_start(iGRU+1) - 1
        iReach = wq%map_reach(iMap)
        reachArea(iReach) = reachArea(iReach) + wq%map_area(iMap)
        reachTair(iReach) = reachTair(iReach) + wq%map_area(iMap) * tairGRU
        reachSW(iReach)   = reachSW(iReach)   + wq%map_area(iMap) * swGRU
      end do
    end do

    ! a reach with no contributing area gets the freezing point and no radiation
    where(reachArea > 0._rkind)
      reachTair = reachTair / reachArea
      reachSW   = reachSW   / reachArea
    elsewhere
      reachTair = 273.15_rkind
      reachSW   = 0._rkind
    end where

    err = openwq_obj%openwq_set_reach_state(nReach, reachVol, reachTair, reachSW, reachArea)

    end associate
#endif
  end subroutine set_reach_state


  ! ============================================================================
  ! Pass the water fluxes of the time step
  ! ============================================================================
  subroutine openwq_run_space_step(summa1_struc)
    USE summa_type, only: summa1_type_dec
    USE globalData, only: gru_struc
    USE globalData, only: realMissing
    USE globalData, only: model_decisions
    USE globalData, only: upland, wetland
    USE mDecisions_module, only: bigBucket, singleBasin
    USE var_lookup, only: iLookPROG, iLookFLUX, iLookINDEX, iLookBVAR, iLookDECISIONS
    USE multiconst, only: iden_water
    implicit none
    type(summa1_type_dec), intent(in) :: summa1_struc
    ! local variables
    integer(i4b)             :: iGRU, iHRU, iDOM, iCol, jCol, iLayer
    integer(i4b)             :: nSnow, nLake, nSoil, nGlce
    integer(i4b)             :: simtime(5)
    integer(i4b)             :: err
    integer(i4b)             :: topCmp, belowSnowCmp      ! compartments below the canopy and below the snowpack
    integer(i4b)             :: toStream                  ! 1 if the water leaving the column goes to the stream
    integer(i4b)             :: aqCol                     ! column that carries the aquifer
    logical(lgt)             :: hasAquifer, basinAquifer
    real(rkind)              :: dt                        ! time step (s)
    real(rkind)              :: toVol                     ! (m s-1) to (m3 per step)
    real(rkind)              :: toVolMass                 ! (kg m-2 s-1) to (m3 per step)
    real(rkind)              :: q, qCanopy, poolVol, totalArea

    call get_simtime(summa1_struc, simtime)
    dt = summa1_struc%data_step

    hasAquifer   = (model_decisions(iLookDECISIONS%groundwatr)%iDecision == bigBucket)
    basinAquifer = (hasAquifer .and. model_decisions(iLookDECISIONS%spatial_gw)%iDecision == singleBasin)

    ! ----- land columns -----
    iCol = 0
    do iGRU = 1, summa1_struc%nGRU_local
      do iHRU = 1, gru_struc(iGRU)%hruCount
        do iDOM = 1, gru_struc(iGRU)%hruInfo(iHRU)%domCount
          iCol = iCol + 1
          if(.not.colActive(iCol)) cycle

          associate( &
            flux    => summa1_struc%fluxStruct%gru(iGRU)%hru(iHRU)%dom(iDOM)%var, &
            prog    => summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%dom(iDOM)%var, &
            indx    => summa1_struc%indxStruct%gru(iGRU)%hru(iHRU)%dom(iDOM)%var, &
            domInfo => gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM))

          nSnow = min(max(0, indx(iLookINDEX%nSnow)%dat(1)), nSnowMax)
          nLake = min(max(0, domInfo%nLake), nLakeMax)
          nSoil = min(max(0, domInfo%nSoil), nSoilMax)
          nGlce = max(0, domInfo%nGlce)

          toVol     = colArea(iCol) * dt
          toVolMass = colArea(iCol) * dt / iden_water
          toStream  = merge(0, 1, domInfo%dom_type == wetland)

          ! store below the snowpack, and store that receives the water below the canopy
          belowSnowCmp = merge(lake_cmp, runoff_cmp, nLake > 0)
          topCmp       = merge(snow_cmp, belowSnowCmp, nSnow > 0)

          ! ----- 1. precipitation and canopy -----
          ! intercepted precipitation: negative when the throughfall includes water released by the canopy
          qCanopy = 0._rkind
          if(prog(iLookPROG%scalarCanopyWat)%dat(1) /= realMissing) &
            qCanopy = ( (flux(iLookFLUX%scalarRainfall)%dat(1) - flux(iLookFLUX%scalarThroughfallRain)%dat(1)) &
                      + (flux(iLookFLUX%scalarSnowfall)%dat(1) - flux(iLookFLUX%scalarThroughfallSnow)%dat(1)) ) * toVolMass
          q = (flux(iLookFLUX%scalarThroughfallRain)%dat(1) + flux(iLookFLUX%scalarThroughfallSnow)%dat(1)) * toVolMass
          call water_in('PRECIP', canopy_cmp, iCol, 1, qCanopy)
          call water_in('PRECIP', topCmp, iCol, 1, q + min(qCanopy, 0._rkind))
          call move_water(canopy_cmp, iCol, 1, topCmp, iCol, 1, -qCanopy, volCanopy(iCol), 0)
          ! drainage and unloading
          if(prog(iLookPROG%scalarCanopyWat)%dat(1) /= realMissing)then
            q = (flux(iLookFLUX%scalarCanopySnowUnloading)%dat(1) + flux(iLookFLUX%scalarCanopyLiqDrainage)%dat(1)) * toVolMass
            call move_water(canopy_cmp, iCol, 1, topCmp, iCol, 1, q, volCanopy(iCol), 0)
          endif

          ! ----- 2. snow -----
          ! interface i of the snow-lake-glacier flux is the base of layer i of the whole stack
          associate(iFluxSnLaGl => flux(iLookFLUX%iLayerLiqFluxSnLaGl)%dat)
          if(nSnow > 0)then
            do iLayer = 1, nSnow - 1
              q = iFluxSnLaGl(lbound(iFluxSnLaGl,1) + iLayer) * toVol
              call move_between(snow_cmp, iLayer, volSnow(iLayer,iCol), snow_cmp, iLayer + 1, volSnow(iLayer+1,iCol), q)
            end do
            ! drainage from the base of the snowpack
            q = flux(iLookFLUX%scalarSnowDrainage)%dat(1) * toVol
            call move_water(snow_cmp, iCol, nSnow, belowSnowCmp, iCol, 1, q, volSnow(nSnow,iCol), 0)
          endif
          ! solute of the snow layers that no longer exist follows the water: to the deepest layer left, or below the snowpack
          do iLayer = nSnow + 1, nSnowStart(iCol)
            if(nSnow > 0)then
              call move_water(snow_cmp, iCol, iLayer, snow_cmp, iCol, nSnow, 1._rkind, 1._rkind, 0)
            else
              call move_water(snow_cmp, iCol, iLayer, belowSnowCmp, iCol, 1, 1._rkind, 1._rkind, 0)
            endif
          end do

          ! ----- 3. lake -----
          if(nLake > 0)then
            do iLayer = 1, nLake - 1
              q = iFluxSnLaGl(lbound(iFluxSnLaGl,1) + nSnow + iLayer) * toVol
              call move_between(lake_cmp, iLayer, volLake(iLayer,iCol), lake_cmp, iLayer + 1, volLake(iLayer+1,iCol), q)
            end do
            ! drainage from the base of the lake
            q = flux(iLookFLUX%scalarLakeDrainage)%dat(1) * toVol
            call move_water(lake_cmp, iCol, nLake, runoff_cmp, iCol, 1, q, volLake(nLake,iCol), 0)
          endif
          end associate

          ! ----- 4. water at the surface (rain plus melt) -----
          poolVol = max(flux(iLookFLUX%scalarRainPlusMelt)%dat(1), 0._rkind) * toVol
          err = openwq_obj%openwq_set_watervol(runoff_cmp, iCol, 1, 1, merge(poolVol, 0._rkind, poolVol >= minPoolDepth*colArea(iCol)))
          ! infiltration
          if(nSoil > 0)then
            q = min(max(flux(iLookFLUX%scalarInfiltration)%dat(1), 0._rkind) * toVol, poolVol)
            call move_water(runoff_cmp, iCol, 1, soil_cmp, iCol, 1, q, poolVol, 0)
          endif
          ! surface runoff
          q = min(max(flux(iLookFLUX%scalarSurfaceRunoff)%dat(1), 0._rkind) * toVol, poolVol)
          call move_water(runoff_cmp, iCol, 1, out_cmp, -1, -1, q, poolVol, toStream)
          err = openwq_obj%openwq_set_fluxvol(runoffVol_exp, iCol, 1, 1, q)

          ! ----- 5. soil -----
          if(nSoil > 0)then
            ! vertical flux between layers (interface i is the base of soil layer i)
            associate(iFluxSoil => flux(iLookFLUX%iLayerLiqFluxSoil)%dat)
            do iLayer = 1, nSoil - 1
              q = iFluxSoil(lbound(iFluxSoil,1) + iLayer) * toVol
              call move_between(soil_cmp, iLayer, volSoil(iLayer,iCol), soil_cmp, iLayer + 1, volSoil(iLayer+1,iCol), q)
            end do
            end associate
            ! lateral outflow (m3 s-1, includes the exfiltration of the first layer)
            jCol = 0
            if(domInfo%dom_type == upland) jCol = colDown(iCol)
            do iLayer = 1, nSoil
              q = flux(iLookFLUX%mLayerColumnOutflow)%dat(iLayer) * dt
              if(jCol > 0)then
                call move_water(soil_cmp, iCol, iLayer, soil_cmp, jCol, iLayer, q, volSoil(iLayer,iCol), 0)
              else
                call move_water(soil_cmp, iCol, iLayer, out_cmp, -1, -1, q, volSoil(iLayer,iCol), toStream)
              endif
            end do
            ! drainage from the base of the soil
            q = flux(iLookFLUX%scalarSoilDrainage)%dat(1) * toVol
            if(hasAquifer .and. domInfo%dom_type == upland)then
              aqCol = merge(gruFirstCol(iGRU), iCol, basinAquifer)
              call move_water(soil_cmp, iCol, nSoil, aquifer_cmp, aqCol, 1, q, volSoil(nSoil,iCol), 0)
            else
              call move_water(soil_cmp, iCol, nSoil, out_cmp, -1, -1, q, volSoil(nSoil,iCol), toStream)
            endif
          endif

          ! ----- 6. aquifer of the column -----
          if(hasAquifer .and. .not.basinAquifer .and. domInfo%dom_type == upland)then
            q = flux(iLookFLUX%scalarAquiferBaseflow)%dat(1) * toVol
            call move_water(aquifer_cmp, iCol, 1, out_cmp, -1, -1, q, volAquifer(iCol), toStream)
          endif

          ! ----- 7. glacier ice melt (negative scalarGlceMelt is melt water leaving the ice) -----
          if(nGlce > 0)then
            q = -flux(iLookFLUX%scalarGlceMelt)%dat(1) * toVol
            call water_in('GLACIER_ICE_MELT', stream_cmp, iCol, 1, q)
          endif

          ! volumes of the runoff exports
          err = openwq_obj%openwq_set_fluxvol(routedRunoff_exp, iCol, 1, 1, &
                max(summa1_struc%bvarStruct%gru(iGRU)%var(iLookBVAR%averageRoutedRunoff)%dat(1), 0._rkind) * toVol)
          err = openwq_obj%openwq_set_fluxvol(totalRunoff_exp, iCol, 1, 1, &
                max(flux(iLookFLUX%scalarTotalRunoff)%dat(1), 0._rkind) * toVol)

          end associate

        end do ! domains
      end do ! HRUs

      ! ----- aquifer of the basin -----
      if(basinAquifer)then
        jCol      = gruFirstCol(iGRU)
        totalArea = summa1_struc%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1)
        q = summa1_struc%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__AquiferBaseflow)%dat(1) * totalArea * dt
        call move_water(aquifer_cmp, jCol, 1, out_cmp, -1, -1, q, volAquifer(jCol), 1)
      endif

      ! ----- water on its way to the stream -----
      call release_to_stream(iGRU)

    end do ! GRUs

    ! ----- river reaches -----
    if(riverCoupled) call route_reaches()

  contains

    ! water flux from a source cell (flux and source volume in m3); recipient out_cmp leaves the domain
    subroutine move_water(sCmp, sCol, sLayer, rCmp, rCol, rLayer, wflux, sVol, outToStream)
      integer(i4b), intent(in) :: sCmp, sCol, sLayer, rCmp, rCol, rLayer
      real(rkind),  intent(in) :: wflux, sVol
      integer(i4b), intent(in) :: outToStream
      integer(i4b)             :: ierr
      real(rkind)              :: vol
      if(.not.(wflux > 0._rkind)) return
      ! a source that holds no water at the start of the step passes the flux through
      vol = sVol
      if(.not.(vol > 0._rkind)) vol = wflux
      if(rCmp == out_cmp)then
        ierr = openwq_obj%openwq_run_space(simtime, sCmp, sCol, 1, sLayer, out_cmp, -1, -1, -1, wflux, vol, outToStream)
      else
        ierr = openwq_obj%openwq_run_space(simtime, sCmp, sCol, 1, sLayer, rCmp, rCol, 1, rLayer, wflux, vol, 0)
      endif
    end subroutine move_water

    ! vertical flux between two layers of the current column; positive is from A to B
    subroutine move_between(aCmp, aLayer, aVol, bCmp, bLayer, bVol, wflux)
      integer(i4b), intent(in) :: aCmp, aLayer, bCmp, bLayer
      real(rkind),  intent(in) :: aVol, bVol, wflux
      if(wflux > 0._rkind)then
        call move_water(aCmp, iCol, aLayer, bCmp, iCol, bLayer, wflux, aVol, 0)
      else
        call move_water(bCmp, iCol, bLayer, aCmp, iCol, aLayer, -wflux, bVol, 0)
      endif
    end subroutine move_between

    ! water entering the domain from an external source (m3)
    subroutine water_in(ewfName, rCmp, rCol, rLayer, wflux)
      character(*), intent(in) :: ewfName
      integer(i4b), intent(in) :: rCmp, rCol, rLayer
      real(rkind),  intent(in) :: wflux
      integer(i4b)             :: ierr
      if(.not.(wflux > 0._rkind)) return
      ierr = openwq_obj%openwq_run_space_in(simtime, ewfName, rCmp, rCol, 1, rLayer, wflux)
    end subroutine water_in

    ! release the RUNOFF_TO_STREAM pools of a GRU in step with the runoff that SUMMA routes
    subroutine release_to_stream(jGRU)
      integer(i4b), intent(in) :: jGRU
      integer(i4b)             :: kCol, ierr
      real(rkind)              :: areaGRU, instantVol, routedVol, mixVol, fracOut, colVol

      associate(bvar => summa1_struc%bvarStruct%gru(jGRU)%var)
      areaGRU    = bvar(iLookBVAR%basin__totalArea)%dat(1)
      instantVol = max(bvar(iLookBVAR%averageInstantRunoff)%dat(1), 0._rkind) * areaGRU * dt
      routedVol  = max(bvar(iLookBVAR%averageRoutedRunoff)%dat(1),  0._rkind) * areaGRU * dt
      end associate

      cumRoutedLand = cumRoutedLand + routedVol

      ! water held at the start of the step plus the runoff of the step, of which routedVol leaves
      mixVol  = gruRoutingStore(jGRU) + instantVol
      fracOut = 0._rkind
      if(mixVol > 0._rkind) fracOut = min(routedVol / mixVol, 1._rkind)
      gruRoutingStore(jGRU) = max(mixVol - routedVol, 0._rkind)

      do kCol = gruFirstCol(jGRU), gruLastCol(jGRU)
        if(.not.colActive(kCol)) cycle
        colVol = 0._rkind
        if(areaGRU > 0._rkind) colVol = mixVol * colArea(kCol) / areaGRU
        ierr = openwq_obj%openwq_set_watervol(stream_cmp, kCol, 1, 1, merge(colVol, 0._rkind, colVol >= minPoolDepth*colArea(kCol)))
        if(fracOut > 0._rkind) call deliver(jGRU, kCol, fracOut)
      end do

    end subroutine release_to_stream

    ! move the fraction fracOut of a RUNOFF_TO_STREAM pool to the reaches fed by the GRU, or out of the domain
    subroutine deliver(jGRU, kCol, fracOut)
      integer(i4b), intent(in) :: jGRU, kCol
      real(rkind),  intent(in) :: fracOut
      logical(lgt)             :: delivered
#ifdef MIZUROUTE_ACTIVE
      integer(i4b)             :: iMap
      real(rkind)              :: areaMap
#endif
      delivered = .false.
#ifdef MIZUROUTE_ACTIVE
      if(riverCoupled)then
        associate(wq => summa1_struc%mizu_domain%river_network%driver%wq)
        areaMap = sum(wq%map_area(wq%map_start(jGRU):wq%map_start(jGRU+1)-1))
        if(areaMap > 0._rkind)then
          do iMap = wq%map_start(jGRU), wq%map_start(jGRU+1) - 1
            call move_water(stream_cmp, kCol, 1, river_cmp, wq%map_reach(iMap), 1, &
                            fracOut * wq%map_area(iMap) / areaMap, 1._rkind, 0)
          end do
          delivered = .true.
        endif
        end associate
      endif
#endif
      if(.not.delivered) call move_water(stream_cmp, kCol, 1, out_cmp, -1, -1, fracOut, 1._rkind, 0)
    end subroutine deliver

    ! advect each reach to its downstream reach with the water budget of mizuRoute for this step
    subroutine route_reaches()
#ifdef MIZUROUTE_ACTIVE
      integer(i4b)             :: iOrder, iReach, jReach, ierr
      real(rkind)              :: mixVol, outVol, wmVol, storageEnd
      associate(wq => summa1_struc%mizu_domain%river_network%driver%wq)
      storageEnd = 0._rkind
      do iOrder = 1, nReach
        iReach = wq%order(iOrder)
        ! water management: an injection (negative) is solute-free water, an abstraction leaves the domain
        wmVol  = wq%vol_wm(iReach)
        mixVol = max(wq%vol_start(iReach), 0._rkind) + max(wq%vol_upstream(iReach), 0._rkind) &
               + max(wq%vol_lateral(iReach), 0._rkind) + max(-wmVol, 0._rkind)
        outVol = max(wq%vol_outflow(iReach), 0._rkind)
        ierr = openwq_obj%openwq_set_watervol(river_cmp, iReach, 1, 1, mixVol)
        jReach = wq%down_index(iReach)
        if(jReach >= 1 .and. jReach <= nReach)then
          call move_water(river_cmp, iReach, 1, river_cmp, jReach, 1, outVol, mixVol, 0)
        else
          call move_water(river_cmp, iReach, 1, out_cmp, -1, -1, outVol, mixVol, 0)
          cumOutflow = cumOutflow + outVol
        endif
        if(wmVol > 0._rkind) call move_water(river_cmp, iReach, 1, out_cmp, -1, -1, wmVol, mixVol, 0)
        ierr = openwq_obj%openwq_set_fluxvol(reachOutflow_exp, iReach, 1, 1, outVol)
        ! running water check: what mizuRoute stored must be what entered minus what left
        cumLateral  = cumLateral + wq%vol_lateral(iReach)
        cumResidual = cumResidual + abs(wq%vol_start(iReach) + wq%vol_upstream(iReach) + wq%vol_lateral(iReach) &
                                        - max(wmVol, 0._rkind) + max(-wmVol, 0._rkind) - wq%vol_outflow(iReach) - wq%vol_end(iReach))
        storageEnd  = storageEnd + wq%vol_end(iReach)
      end do
      cumStorageEnd = storageEnd
      end associate
#endif
    end subroutine route_reaches

  end subroutine openwq_run_space_step


  ! ============================================================================
  ! Solve the time step and write the output
  ! ============================================================================
  subroutine openwq_run_time_end(summa1_struc)
    USE summa_type, only: summa1_type_dec
    implicit none
    type(summa1_type_dec), intent(in) :: summa1_struc
    integer(i4b)                      :: simtime(5)
    integer(i4b)                      :: err

    call get_simtime(summa1_struc, simtime)
    err = openwq_obj%openwq_run_time_end(simtime)

  end subroutine openwq_run_time_end


  ! ============================================================================
  ! Report the water check of the coupling at the end of the run
  ! ============================================================================
  subroutine openwq_finalize()
    implicit none
    if(.not.riverCoupled) return
    write(*,'(a)')          ' OpenWQ river coupling, water check over the run (m3):'
    write(*,'(a,es12.4)')   '   runoff routed out of the GRUs by SUMMA   ', cumRoutedLand
    write(*,'(a,es12.4,a,f7.2,a)') '   lateral inflow received by the reaches    ', cumLateral, &
      '  (', 100._rkind*(cumLateral/max(cumRoutedLand,tiny(1._rkind)) - 1._rkind), ' % of the runoff)'
    write(*,'(a,es12.4)')   '   outflow of the outlet reach(es)           ', cumOutflow
    write(*,'(a,es12.4)')   '   reach storage at the end                  ', cumStorageEnd
    write(*,'(a,es12.4,a,es9.2,a)') '   reach budget residual, summed |error|     ', cumResidual, &
      '  (', cumResidual/max(cumLateral,tiny(1._rkind)), ' of the lateral inflow)'
  end subroutine openwq_finalize


  ! year, month, day, hour, minute of the current time step
  subroutine get_simtime(summa1_struc, simtime)
    USE summa_type, only: summa1_type_dec
    USE var_lookup, only: iLookTIME
    implicit none
    type(summa1_type_dec), intent(in)  :: summa1_struc
    integer(i4b),          intent(out) :: simtime(5)
    simtime(1) = summa1_struc%timeStruct%var(iLookTIME%iyyy)
    simtime(2) = summa1_struc%timeStruct%var(iLookTIME%im)
    simtime(3) = summa1_struc%timeStruct%var(iLookTIME%id)
    simtime(4) = summa1_struc%timeStruct%var(iLookTIME%ih)
    simtime(5) = summa1_struc%timeStruct%var(iLookTIME%imin)
  end subroutine get_simtime

end module summa_openwq

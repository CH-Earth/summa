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

module glacAreaChange_module

! data types
USE nr_type
USE data_types,only:&
                    var_ilength,     & ! x%var(:)%dat (i4b)
                    var_dlength,     & ! x%var(:)%dat (rkind)
                    glac_info,       & ! glacier information data structure
                    grid_info,       & ! glacier grid info data structure
                    grid_double        ! x%gru(:)%grid(:)%var(:)%dat2(:,:) (dp)

! access missing values
USE globalData,only:integerMissing     ! missing integer number
USE globalData,only:realMissing        ! missing real number

! constants
USE globalData,only:verySmall          ! a small number
USE globalData,only:thick4area         ! an arbitrary small threshold for glacier thickness to be considered as glacier area
USE globalData,only:dJulianStart       ! julian day of start time of simulation
USE globalData,only:data_step          ! length of time steps for the outermost timeloop
USE multiconst,only:&
                    secprday,        & ! seconds per day
                    gravity,         & ! gravitational acceleration    (m s-2)
                    iden_water,      & ! intrinsic density of water    (kg m-3)
                    iden_ice           ! intrinsic density of ice      (kg m-3)

! define data types
USE var_lookup,only:iLookGRID          ! named variables for the glacier grid information
USE var_lookup,only:iLookPROG          ! named variables for the prognostic variables

! access domain types
USE globalData,only:upland             ! horizontal domain type for upland areas
USE globalData,only:glacCln1           ! first horizontal domain type for glacier clean areas
USE globalData,only:glacCln2           ! second horizontal domain type for glacier clean areas
USE globalData,only:glacDbr            ! horizontal domain type for glacier debris areas
USE globalData,only:wetland            ! horizontal domain type for wetland areas                

implicit none

! shared constants
real(rkind),parameter :: flat_threshold=1.e-6_rkind  ! tan slope below which terrain is treated as flat for aspect
real(rkind),parameter :: deg2rad=PI_D/180._rkind     ! factor to convert degrees to radians
real(rkind),parameter :: rad2deg=180._rkind/PI_D     ! factor to convert radians to degrees

! privacy
private::run_flowModel,run_debrisModel,diffusion_MUSCL,advection_MUSCL,static_emergElev_latRockfall
private::fit_massBalanceCurve,cell_slopeAspect,argsort_desc
private::superbee,flux,SIA,midpt,pluss,minus
public::glacAreaChange
public::time_updateGlacArea
public::updateGlacDomain
contains


! ************************************************************************************************
! public subroutine glacAreaChange: get new glacier area and elevation
! ************************************************************************************************
subroutine glacAreaChange(&
                    ! model control
                    t_total,                 & ! intent(in):    total time to run (s), time since last update
                    nHRU,                    & ! intent(in):    number of HRUs that have a glacier domain
                    nDOM,                    & ! intent(in):    number of glacier domains
                    ndebris,                 & ! intent(in):    number of debris domains in each HRU
                    nclean,                  & ! intent(in):    number of clean domains in each HRU
                    hruInd,                  & ! intent(in):    hruInd of each glacier domain
                    ! glacier topography
                    nGlac,                   & ! intent(in):    number of glaciers in GRU
                    glacInfo,                & ! intent(in):    information for each glacier
                    gridInfo,                & ! intent(in):    information for each grid
                    gridData,                & ! intent(inout): grid data for each glacier
                    ! mass balance
                    dom_massChange,          & ! intent(in):    change in glacier water equivalent (kg m-2) in each glacier domain since last update
                    dom_elev,                & ! intent(inout): elevation in each glacier domain (m)
                    dom_tan_slope,           & ! intent(inout): tan local ground surface slope in each glacier domain (m/m)
                    dom_aspect,              & ! intent(inout): azimuth in degrees East of North of each glacier domain (degrees)
                    dom_contourLength,       & ! intent(inout): length of contour at downslope edge of each glacier domain (m)
                    ! debris
                    dom_debris_thick,        & ! intent(inout): debris thickness in each glacier domain (m)
                    iden_soil_mean,          & ! intent(in):    mean intrinsic density of soil in each glacier domain (kg m-3)
                    theta_sat_mean,          & ! intent(in):    mean soil porosity in each glacier domain (-)
                    debrisConc,              & ! intent(in):    englacial debris concentration (kg m-3)
                    wallErosionRate,         & ! intent(in):    glacier wall erosion rate input for debris advection (mm yr-1)
                    debrisCritStress,        & ! intent(in):    critical driving stress where debris slides on terminal wedge (Pa)
                    latMoraineWidth,         & ! intent(in):    lateral moraine width (rockfall length) (m)
                    ! area
                    glacierAblArea,          & ! intent(inout): per glacier ablation area (m2)
                    glacierAccArea,          & ! intent(inout): per glacier accumulation area (m2)
                    dom_area,                & ! intent(inout): area of each domain (m2)
                    dom_ablFrac,             & ! intent(out):   per domain ablation fraction (-)
                    ! error handling
                    err, message)              ! intent(out):   error control
  ! ---------------------------------------------------------------------------------------------
  implicit none
  ! model control
  real(rkind), intent(in)            :: t_total                         ! total time to run (s), time since last update
  integer(i4b), intent(in)           :: nHRU                            ! number of HRUs that have a glacier domain
  integer(i4b), intent(in)           :: nDOM                            ! number of glacier domains
  integer(i4b), intent(in)           :: ndebris(:)                      ! number of debris domains in each HRU
  integer(i4b), intent(in)           :: nclean(:)                       ! number of clean domains in each HRU
  integer(i8b), intent(in)           :: hruInd(:)                       ! hruInd of each glacier domain
  ! glacier topograpy
  integer(i4b), intent(in)           :: nGlac                           ! number of glaciers in GRU
  type(glac_info), intent(in)        :: glacInfo(:)                     ! information for each glacier
  type(grid_info), intent(in)        :: gridInfo(:)                     ! information for each grid
  type(grid_double), intent(inout)   :: gridData                        ! data for each grid
  ! mass balance, realMissing value if domain is missing (i.e. glacier does not have one of ablation or accumulation)
  real(rkind), intent(in)            :: dom_massChange(:)               ! change in glacier water equivalent (kg m-2) in each glacier domain
  real(rkind), intent(inout)         :: dom_elev(:)                     ! elevation of each glacier domain (m)
  real(rkind), intent(inout)         :: dom_tan_slope(:)                ! tan local ground surface slope in each glacier domain (m/m)
  real(rkind), intent(inout)         :: dom_aspect(:)                   ! azimuth in degrees East of North of each glacier domain (degrees)
  real(rkind), intent(inout)         :: dom_contourLength(:)            ! length of contour at downslope edge of each glacier domain (m)
  ! debris
  real(rkind), intent(inout)         :: dom_debris_thick(:)             ! debris thickness in glacier domain (m)
  real(rkind), intent(in)            :: iden_soil_mean(:)               ! mean intrinsic density of soil (kg m-3)
  real(rkind), intent(in)            :: theta_sat_mean(:)               ! mean soil porosity (-)
  real(rkind), intent(in)            :: debrisConc                      ! englacial debris concentration (kg m-3)
  real(rkind), intent(in)            :: wallErosionRate                 ! glacier wall erosion rate input for debris advection (mm yr-1)
  real(rkind), intent(in)            :: debrisCritStress                ! critical driving stress where debris slides on terminal wedge (Pa)
  real(rkind), intent(in)            :: latMoraineWidth                 ! lateral moraine width (rockfall length) (m)
  ! area
  real(rkind), intent(inout)         :: glacierAblArea(nGlac)           ! per glacier ablation area (m2)
  real(rkind), intent(inout)         :: glacierAccArea(nGlac)           ! per glacier accumulation area (m2)
  real(rkind), intent(inout)         :: dom_area(:)                     ! area of each glacier domain (m2)
  real(rkind), intent(out)           :: dom_ablFrac(nDOM)               ! per domain ablation fraction (-)
  integer(i4b),intent(out)           :: err                             ! error code
  character(*),intent(out)           :: message                         ! error message
  ! locals: mass balance vs elevation, one column per debris condition (1 = clean, 2 = debris covered)
  logical(lgt)                       :: isValid(nDOM,2)                 ! domain has area and the debris condition of the column
  integer(i4b)                       :: validCount(2)                   ! number of valid points for each debris condition
  integer(i4b)                       :: maxCount                        ! maximum number of valid points
  real(rkind), allocatable           :: validElev(:,:)                  ! elevation of the valid points (m)
  real(rkind), allocatable           :: validMassChange(:,:)            ! mass change of the valid points (kg m-2)
  real(rkind), allocatable           :: slope(:,:), intercept(:,:)      ! piecewise linear mass change as a function of elevation
  real(rkind)                        :: ELA_elev(2)                     ! ELA elevation for clean and debris ablation (m)
  real(rkind)                        :: ELA_use                         ! ELA used for domain calculations (m)
  real(rkind)                        :: ELA_use_glac                    ! ELA_use clipped to the surface range of the current glacier (m)
  real(rkind)                        :: iden_soil, theta_sat            ! mean soil properties of the debris-covered domains
  ! locals: glacier grid
  integer(i4b), allocatable          :: glacid_to_index(:)              ! index in gridInfo of each glacier
  integer(i4b)                       :: nx, ny                          ! number of grid cells in x and y directions
  real(rkind)                        :: dx, dy                          ! grid cell size in x and y directions (m)
  real(rkind), allocatable           :: surface(:,:)                    ! surface elevation of each cell (m)
  real(rkind), allocatable           :: bed(:,:)                        ! bed elevation of each cell (m)
  real(rkind), allocatable           :: hgt(:,:)                        ! ice thickness of each cell (m)
  real(rkind), allocatable           :: debris(:,:)                     ! debris thickness of each cell (m)
  real(rkind), allocatable           :: cell_tan_slope(:,:)             ! tan local ground surface slope of each cell (m/m)
  real(rkind), allocatable           :: cell_aspect(:,:)                ! aspect of each cell as compass bearing (degrees)
  integer(i4b), allocatable          :: cell2hru(:,:)                   ! HRU index of each cell
  integer(i4b), allocatable          :: glacierMask(:,:)                ! 1-0 mask of cells the glacier can occupy
  logical(lgt), allocatable          :: glacClnMask(:,:)                ! cells in the clean glacier area of the current HRU
  logical(lgt), allocatable          :: glacDbrMask(:,:)                ! cells in the debris-covered glacier area of the current HRU
  logical(lgt), allocatable          :: glacHiMask(:,:)                 ! cells given to the higher of two clean domains
  integer(i4b), allocatable          :: sortedIndices(:)                ! cell indices in order of descending surface elevation
  real(rkind)                        :: volume                          ! volume of each glacier (km3)
  real(rkind)                        :: totVolume                       ! total volume of all glaciers (km3), not currently used
  ! locals: domain sums
  real(rkind)                        :: elev0(nDOM)                     ! elevation of each glacier domain at entry (m)
  real(rkind)                        :: area0(nDOM)                     ! area of each glacier domain at entry (m2)
  real(rkind)                        :: aspect_sin_sum(nDOM)            ! running sine sum for area-weighted domain aspect averaging
  real(rkind)                        :: aspect_cos_sum(nDOM)            ! running cosine sum for area-weighted domain aspect averaging
  real(rkind)                        :: areaAdd                         ! area to give to the higher of two clean domains (m2)
  real(rkind)                        :: cumulativeArea                  ! area given so far to the higher clean domain (m2)
  integer(i4b)                       :: hiInd, loInd                    ! domain indices of the higher and lower clean domains
  integer(i4b)                       :: i, j, k, n, dbr                 ! loop indices
  integer(i4b)                       :: iGlac, iGrid, iDOM, iHRU        ! loop indices
  character(len=256)                 :: cmessage                        ! error message of downwind routine
  real(rkind),parameter              :: min_thickness=0.02_rkind        ! minimum thickness of debris cover to be considered as debris cover (m)
  logical(lgt),parameter             :: printFlag=.false.               ! flag to print details for debugging

  ! ----------------------------------------------------------------------------------------------
  ! initialize
  err=0; message='glacAreaChange/'
  totVolume = 0._rkind
  elev0 = dom_elev
  area0 = dom_area

  ! ----------------------------------------------------------------------------------------------
  ! Fit mass change as a piecewise linear function of elevation, once through the clean domains
  !  (column 1) and once through the debris-covered domains (column 2). Domains with no area have
  !  realMissing elevation and are excluded; a column is empty if the GRU has no domains of that kind.
  ! ----------------------------------------------------------------------------------------------
  isValid(:,1) = dom_elev/=realMissing .and. dom_debris_thick==0._rkind
  isValid(:,2) = dom_elev/=realMissing .and. dom_debris_thick> 0._rkind
  if(sum(nclean)==0)  isValid(:,1) = .false.
  if(sum(ndebris)==0) isValid(:,2) = .false.
  validCount = count(isValid, dim=1)
  maxCount = maxval(validCount)

  ! mean soil properties of the debris-covered domains (usually HRU independent, so a mean of means is adequate)
  iden_soil = 0._rkind
  theta_sat = 0._rkind
  if(validCount(2)>0)then
    iden_soil = sum(iden_soil_mean, mask=isValid(:,2))/validCount(2)
    theta_sat = sum(theta_sat_mean, mask=isValid(:,2))/validCount(2)
  endif

  allocate(validElev(maxCount,2), validMassChange(maxCount,2), slope(maxCount+1,2), intercept(maxCount+1,2))
  validElev = realMissing
  validMassChange = realMissing
  slope = 0._rkind
  intercept = 0._rkind
  ELA_elev = realMissing
  do dbr = 1,2
    if(validCount(dbr)==0) cycle
    validElev(1:validCount(dbr),dbr)       = pack(dom_elev,       isValid(:,dbr))
    validMassChange(1:validCount(dbr),dbr) = pack(dom_massChange, isValid(:,dbr))
    call fit_massBalanceCurve(validCount(dbr), validElev(:,dbr), validMassChange(:,dbr), slope(:,dbr), intercept(:,dbr), ELA_elev(dbr))
    if(printFlag)then
      do i = 1, validCount(dbr)
        write(*,'(a,i1,a,2(1x,f7.1),a,2(1x,e9.2))') "Debris presence ", dbr-1, "   Valid sorted elevation, mass change =", &
              validElev(i,dbr), validMassChange(i,dbr), "   slope, intercept = ", slope(i,dbr), intercept(i,dbr)
      enddo
      write(*,'(a,i1,a,2(1x,e9.2),a,i8)') "Debris presence ", dbr-1, "   Extrapolation below: slope, intercept =", &
              slope(validCount(dbr)+1,dbr), intercept(validCount(dbr)+1,dbr), "   ELA elevation = ", int(ELA_elev(dbr))
    endif
  enddo

  ! use the clean ELA, falling back to the debris ELA if there is no clean ablation
  ELA_use = ELA_elev(1)
  if(ELA_use<0._rkind) ELA_use = ELA_elev(2)

  if(printFlag)then
    do i = 1, nDOM
      write(*,'(a,1x,f7.1,1x,f4.1,1x,f5.2,1x,f7.1)') "Original domain elevation (m), dom_area (km2), debris depth (m), mass change (kg m-2) =",&
            dom_elev(i), dom_area(i)*1.e-6_rkind, dom_debris_thick(i), dom_massChange(i)
    enddo
  endif

  ! ----------------------------------------------------------------------------------------------
  ! Run the flow model for each glacier and rebuild the domains from the new grids
  ! ----------------------------------------------------------------------------------------------
  ! reset the domain variables; they accumulate area-weighted sums over the grids below
  dom_area          = 0._rkind
  dom_ablFrac       = 0._rkind
  dom_elev          = 0._rkind
  dom_tan_slope     = 0._rkind
  dom_aspect        = 0._rkind
  dom_contourLength = 0._rkind
  dom_debris_thick  = 0._rkind
  aspect_sin_sum    = 0._rkind
  aspect_cos_sum    = 0._rkind

  ! map each glacier to its grid
  allocate(glacid_to_index(nGlac))
  do iGlac = 1, nGlac
    glacid_to_index(iGlac) = findloc(gridInfo(:)%grid_id, glacInfo(iGlac)%glac_id, dim=1)
    if(glacid_to_index(iGlac)==0)then
      write(cmessage,'(i0)') glacInfo(iGlac)%glac_id
      message=trim(message)//'no grid found for glacier id '//trim(cmessage)
      err=20; return
    endif
  enddo

  do iGlac = 1, nGlac
    if(glacierAblArea(iGlac) + glacierAccArea(iGlac) <= 0._rkind)then
      if(printFlag) write(*,'(a,i2,a)') ">GLACIER ",iGlac, " SKIP: no glacier area"
      cycle
    endif

    ! set up the glacier grid
    iGrid = glacid_to_index(iGlac)
    nx = gridInfo(iGrid)%nx
    ny = gridInfo(iGrid)%ny
    dx = gridInfo(iGrid)%dx
    dy = gridInfo(iGrid)%dy
    allocate(surface(nx,ny), bed(nx,ny), hgt(nx,ny), debris(nx,ny), cell_tan_slope(nx,ny), cell_aspect(nx,ny), &
             cell2hru(nx,ny), glacierMask(nx,ny), glacClnMask(nx,ny), glacDbrMask(nx,ny), glacHiMask(nx,ny))
    surface     = gridData%grid(iGrid)%var(iLookGRID%surface_elev)%dat2(1:nx,1:ny)
    bed         = gridData%grid(iGrid)%var(iLookGRID%bed_elev)%dat2(1:nx,1:ny)
    debris      = gridData%grid(iGrid)%var(iLookGRID%debris_thick)%dat2(1:nx,1:ny)
    cell2hru    = int(gridData%grid(iGrid)%var(iLookGRID%cell2hru)%dat2(1:nx,1:ny))
    glacierMask = int(gridData%grid(iGrid)%var(iLookGRID%glacierMask)%dat2(1:nx,1:ny))
    if(printFlag)then
      volume = sum(surface - bed) * dx * dy * 1.e-9_rkind ! km3
      write(*,'(a,i2,a,4(1x,f8.3))') "<GLACIER ",iGlac, " START: ablation, accumulation, total area (km2), volume (km3) =", &
        glacierAblArea(iGlac)*1.e-6_rkind, glacierAccArea(iGlac)*1.e-6_rkind, (glacierAblArea(iGlac) + glacierAccArea(iGlac))*1.e-6_rkind, volume
    endif

    ! keep the ELA within the surface range of this glacier so it can have a debris emergence area
    ELA_use_glac = min(max(ELA_use, minval(surface)-verySmall), maxval(surface)+verySmall)

    ! run the flow model, which updates surface and debris
    call run_flowModel(t_total, debris, surface, bed, glacierMask, slope, intercept, validElev, validCount, maxCount, debrisConc, &
                       wallErosionRate, debrisCritStress, latMoraineWidth, iden_soil, theta_sat, ELA_use_glac, nx, ny, dx, dy, volume, printFlag, &
                       err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    totVolume = totVolume + volume ! includes debris, not currently used

    ! new glacier geometry
    hgt = surface - bed
    call cell_slopeAspect(surface, glacierMask==1_i4b .and. hgt>thick4area, nx, ny, dx, dy, cell_tan_slope, cell_aspect)
    glacierAccArea(iGlac) = count(hgt>thick4area .and. surface>=ELA_use_glac)*dx*dy
    glacierAblArea(iGlac) = count(hgt>thick4area .and. surface< ELA_use_glac)*dx*dy
    if(printFlag) write(*,'(a,i2,a,4(1x,f8.3))') ">GLACIER ",iGlac, " END  : ablation, accumulation, total area (km2), volume (km3) =", &
        glacierAblArea(iGlac)*1.e-6_rkind, glacierAccArea(iGlac)*1.e-6_rkind, (glacierAblArea(iGlac) + glacierAccArea(iGlac))*1.e-6_rkind, volume

    ! cells in descending elevation order, needed to split an HRU's clean area between two clean domains
    if(any(nclean>1))then
      allocate(sortedIndices(nx*ny))
      call argsort_desc(reshape(surface,[nx*ny]), sortedIndices)
    endif

    ! Rebuild the domains of each HRU from the grid. Domains are ordered HRU 1: clean1, clean2, debris;
    !  HRU 2: clean1, clean2, debris; etc., with a clean or debris domain absent if the HRU does not have one.
    n = 0 ! index of the last domain filled
    do iHRU = 1, nHRU
      glacClnMask = cell2hru==hruInd(n+1) .and. hgt>thick4area .and. debris< min_thickness
      glacDbrMask = cell2hru==hruInd(n+1) .and. hgt>thick4area .and. debris>=min_thickness
      ! if the HRU lacks a domain, its area goes to the other kind (lowering the mean debris thickness, or making it zero)
      if(nclean(iHRU)==0)  glacDbrMask = glacDbrMask .or. glacClnMask
      if(ndebris(iHRU)==0) glacClnMask = glacClnMask .or. glacDbrMask

      if(nclean(iHRU)>1)then
        ! two clean domains: keep their area ratio (so both persist), giving the highest cells to the higher domain
        n = n+2
        hiInd = n
        loInd = n-1
        if(elev0(n-1)>elev0(n))then
          hiInd = n-1
          loInd = n
        endif
        areaAdd = 0._rkind
        if(area0(n-1)+area0(n)>0._rkind) areaAdd = area0(hiInd)/(area0(n-1)+area0(n)) * count(glacClnMask) *dx*dy
        glacHiMask = .false.
        cumulativeArea = 0._rkind
        do i = 1, nx*ny
          if(cumulativeArea >= areaAdd) exit
          j = mod(sortedIndices(i)-1, nx) + 1
          k = (sortedIndices(i)-1)/nx + 1
          if(glacClnMask(j,k))then
            glacHiMask(j,k) = .true.
            cumulativeArea = cumulativeArea + dx*dy
          endif
        enddo
        call addCellsToDomain(hiInd, glacHiMask)
        call addCellsToDomain(loInd, glacClnMask .and. .not.glacHiMask)
      elseif(nclean(iHRU)==1)then
        n = n+1
        call addCellsToDomain(n, glacClnMask)
      endif

      if(ndebris(iHRU)>0)then
        n = n+1 ! only one debris domain is currently possible
        call addCellsToDomain(n, glacDbrMask)
        dom_debris_thick(n) = dom_debris_thick(n) + sum(debris, mask=glacDbrMask)*dx*dy
      endif
    enddo ! HRU loop

    ! write the new surface and debris back to the grid
    gridData%grid(iGrid)%var(iLookGRID%surface_elev)%dat2(1:nx,1:ny) = surface
    gridData%grid(iGrid)%var(iLookGRID%debris_thick)%dat2(1:nx,1:ny) = debris
    deallocate(surface, bed, hgt, debris, cell_tan_slope, cell_aspect, cell2hru, glacierMask, glacClnMask, glacDbrMask, glacHiMask)
    if(allocated(sortedIndices)) deallocate(sortedIndices)
  enddo ! glacier loop

  ! convert the area-weighted sums to domain means; a domain with no area gets missing values
  do iDOM = 1,nDOM
    if(dom_area(iDOM)>0._rkind)then
      dom_elev(iDOM)          = dom_elev(iDOM)/dom_area(iDOM)
      dom_tan_slope(iDOM)     = dom_tan_slope(iDOM)/dom_area(iDOM)
      dom_ablFrac(iDOM)       = dom_ablFrac(iDOM)/dom_area(iDOM)
      dom_debris_thick(iDOM)  = dom_debris_thick(iDOM)/dom_area(iDOM)
      dom_aspect(iDOM)        = 0._rkind
      if(aspect_sin_sum(iDOM)**2 + aspect_cos_sum(iDOM)**2 > 0._rkind) &
        dom_aspect(iDOM) = modulo(atan2(aspect_sin_sum(iDOM), aspect_cos_sum(iDOM))*rad2deg, 360._rkind)
    else
      dom_elev(iDOM)          = realMissing
      dom_tan_slope(iDOM)     = realMissing
      dom_aspect(iDOM)        = realMissing
      dom_area(iDOM)          = 0._rkind
      dom_ablFrac(iDOM)       = 0._rkind
      dom_contourLength(iDOM) = 0._rkind
      dom_debris_thick(iDOM)  = 0._rkind
    endif
  enddo
  if(printFlag)then
    do i = 1, nDOM
      write(*,'(a,2(1x,f8.2),1x,f5.2)') "After area change domain elevation (m), area (km2), debris depth (m) =", &
            dom_elev(i), dom_area(i)*1.e-6_rkind, dom_debris_thick(i)
    enddo
  endif

  deallocate(glacid_to_index, validElev, validMassChange, slope, intercept)

contains

  ! add the cells selected by mask to the running area-weighted sums of domain iDOM
  subroutine addCellsToDomain(iDOM, mask)
    implicit none
    integer(i4b), intent(in) :: iDOM                 ! domain index
    logical(lgt), intent(in) :: mask(nx,ny)          ! cells to add
    real(rkind)              :: area                 ! area of the cells (m2)
    area = count(mask)*dx*dy
    dom_area(iDOM)          = dom_area(iDOM)          + area
    dom_elev(iDOM)          = dom_elev(iDOM)          + sum(surface, mask=mask)*dx*dy
    dom_tan_slope(iDOM)     = dom_tan_slope(iDOM)     + sum(cell_tan_slope, mask=mask)*dx*dy
    aspect_sin_sum(iDOM)    = aspect_sin_sum(iDOM)    + sum(sin(cell_aspect*deg2rad), mask=mask .and. cell_tan_slope>=flat_threshold)*dx*dy
    aspect_cos_sum(iDOM)    = aspect_cos_sum(iDOM)    + sum(cos(cell_aspect*deg2rad), mask=mask .and. cell_tan_slope>=flat_threshold)*dx*dy
    dom_contourLength(iDOM) = dom_contourLength(iDOM) + sqrt(area)
    dom_ablFrac(iDOM)       = dom_ablFrac(iDOM)       + count(mask .and. surface<ELA_use_glac)*dx*dy
  end subroutine addCellsToDomain

end subroutine glacAreaChange


! ************************************************************************************************
! private subroutine fit_massBalanceCurve: piecewise linear mass change as a function of elevation
!   through the domain points, and the ELA where it crosses zero
! ************************************************************************************************
subroutine fit_massBalanceCurve(nPts, elev, massChange, slope, intercept, ELA)
  implicit none
  integer(i4b), intent(inout) :: nPts           ! number of points; reduced on return if points shared an elevation
  real(rkind), intent(inout)  :: elev(:)        ! elevation of each point (m); sorted in descending order on return
  real(rkind), intent(inout)  :: massChange(:)  ! mass change of each point (kg m-2); sorted with elev, averaged over shared elevations
  real(rkind), intent(out)    :: slope(:)       ! slope of segment i, which joins points i-1 and i; segments 1 and nPts+1 extrapolate above and below (kg m-2 m-1)
  real(rkind), intent(out)    :: intercept(:)   ! intercept of segment i (kg m-2)
  real(rkind), intent(out)    :: ELA            ! elevation where mass change crosses zero (m)
  ! locals
  integer(i4b)                :: idx(nPts)      ! order of the points
  integer(i4b)                :: i, j, m        ! indices
  integer(i4b)                :: ind            ! segment where mass change changes sign
  real(rkind), parameter      :: allAbl= 1.e6_rkind ! ELA when every point is ablation (m)
  real(rkind), parameter      :: allAcc=-1.e6_rkind ! ELA when every point is accumulation (m)

  ! sort the points by descending elevation
  call argsort_desc(elev(1:nPts), idx)
  elev(1:nPts) = elev(idx)
  massChange(1:nPts) = massChange(idx)

  ! average the points that share an elevation
  m = 0
  i = 1
  do while(i <= nPts)
    j = i
    do while(j < nPts)
      if(elev(j+1) /= elev(i)) exit
      j = j + 1
    enddo
    m = m + 1
    elev(m) = elev(i)
    massChange(m) = sum(massChange(i:j))/(j-i+1)
    i = j + 1
  enddo
  nPts = m

  slope = 0._rkind
  intercept = 0._rkind
  if(nPts == 1)then
    intercept(1) = massChange(1)
    if(massChange(1) > 0._rkind)then
      ELA = allAcc
    else
      ELA = allAbl
    endif
    return
  endif

  ! segment i joins points i-1 and i; ind is the segment where mass change changes sign
  ind = 0
  do i = 1, nPts-1
    slope(i+1) = (massChange(i+1) - massChange(i))/(elev(i+1) - elev(i))
    intercept(i+1) = massChange(i) - slope(i+1)*elev(i)
    if(massChange(i) >= 0._rkind) ind = i+1
  enddo
  if(ind == nPts .and. massChange(nPts) >= 0._rkind) ind = nPts+1 ! every point is accumulation
  if(ind == 0) ind = 1                                            ! every point is ablation

  ! extrapolate above and below with the end segments, flattened so a mass balance inversion is not propagated
  slope(1) = slope(2)
  intercept(1) = intercept(2)
  if(slope(1) < 0._rkind)then
    slope(1) = 0._rkind
    intercept(1) = massChange(1)
  endif
  slope(nPts+1) = slope(nPts)
  intercept(nPts+1) = intercept(nPts)
  if(slope(nPts+1) < 0._rkind)then
    slope(nPts+1) = 0._rkind
    intercept(nPts+1) = massChange(nPts)
  endif

  ! ELA is where the sign-change segment crosses zero
  if(slope(ind) /= 0._rkind)then
    ELA = -intercept(ind)/slope(ind)
  elseif(ind == 1)then
    ELA = allAbl
  elseif(ind == nPts+1)then
    ELA = allAcc
  else
    ELA = elev(ind)
  endif

end subroutine fit_massBalanceCurve


! ************************************************************************************************
! private subroutine cell_slopeAspect: tan slope and aspect of each grid cell from the surface
!   elevation, using one-sided differences at the edges and centered differences in the interior
! ************************************************************************************************
subroutine cell_slopeAspect(surface, onGlacier, nx, ny, dx, dy, tan_slope, aspect)
  implicit none
  integer(i4b), intent(in)    :: nx, ny             ! number of grid cells in x and y directions
  real(rkind), intent(in)     :: surface(nx,ny)     ! surface elevation (m)
  logical(lgt), intent(in)    :: onGlacier(nx,ny)   ! cells with glacier ice
  real(rkind), intent(in)     :: dx, dy             ! grid cell size (m)
  real(rkind), intent(out)    :: tan_slope(nx,ny)   ! tan local surface slope (m/m), zero off the glacier
  real(rkind), intent(out)    :: aspect(nx,ny)      ! compass bearing of the downslope direction (degrees), zero where flat or off the glacier
  ! locals
  real(rkind)                 :: dzdx(nx,ny)        ! x-gradient of surface elevation (m/m)
  real(rkind)                 :: dzdy(nx,ny)        ! y-gradient of surface elevation (m/m)

  dzdx = 0._rkind
  dzdy = 0._rkind
  if(nx > 1)then
    dzdx(1,:)  = (surface(2,:) - surface(1,:)) / dx
    dzdx(nx,:) = (surface(nx,:) - surface(nx-1,:)) / dx
    if(nx > 2) dzdx(2:nx-1,:) = (surface(3:nx,:) - surface(1:nx-2,:)) / (2._rkind*dx)
  endif
  if(ny > 1)then
    dzdy(:,1)  = (surface(:,2) - surface(:,1)) / dy
    dzdy(:,ny) = (surface(:,ny) - surface(:,ny-1)) / dy
    if(ny > 2) dzdy(:,2:ny-1) = (surface(:,3:ny) - surface(:,1:ny-2)) / (2._rkind*dy)
  endif

  tan_slope = sqrt(dzdx**2_i4b + dzdy**2_i4b)
  aspect = 0._rkind
  where(.not.onGlacier)
    tan_slope = 0._rkind
  elsewhere(tan_slope >= flat_threshold)
    aspect = modulo(90._rkind - atan2(-dzdx, dzdy)*rad2deg, 360._rkind)
  end where

end subroutine cell_slopeAspect


! ************************************************************************************************
! private subroutine argsort_desc: indices that order x in descending order (stable insertion sort)
! ************************************************************************************************
subroutine argsort_desc(x, idx)
  implicit none
  real(rkind), intent(in)     :: x(:)               ! values to order
  integer(i4b), intent(out)   :: idx(size(x))       ! indices of x in descending order of value
  ! locals
  integer(i4b)                :: i, j, tmp

  idx = [(i, i=1, size(x))]
  do i = 2, size(x)
    tmp = idx(i)
    j = i - 1
    do while(j >= 1)
      if(x(idx(j)) >= x(tmp)) exit
      idx(j+1) = idx(j)
      j = j - 1
    enddo
    idx(j+1) = tmp
  enddo

end subroutine argsort_desc

! ************************************************************************************************
! private subroutine run_flowModel: set up Shallow Ice Approximation (SIA) diffusion flow model 
!   and run for time period
!   Follows the implementation of Jarosch et al., 2013
! ************************************************************************************************
subroutine run_flowModel(t_total, debris, S, B, glacierMask, slope, intercept, validElev, validCount, maxCount, debrisConc, &
                         wallErosionRate, debrisCritStress, latMoraineWidth, iden_soil, theta_sat, ELA, nx, ny, dx, dy, volume, printFlag, &
                         err, message)
  implicit none
  ! Arguments
  real(rkind), intent(in) :: t_total, dx, dy, B(nx,ny)
  real(rkind), intent(inout) :: debris(nx,ny), S(nx,ny)
  real(rkind), intent(in) :: slope(maxCount+1,2), intercept(maxCount+1,2), validElev(maxCount,2), debrisConc
  real(rkind), intent(in) :: wallErosionRate, debrisCritStress, latMoraineWidth, iden_soil, theta_sat, ELA
  integer(i4b), intent(in) :: validCount(2), maxCount
  integer(i4b), intent(in) :: ny, nx, glacierMask(nx,ny)
  logical(lgt), intent(in) :: printFlag
  real(rkind), intent(out) :: volume
  integer(i4b), intent(out) :: err                ! error code
  character(*), intent(out) :: message             ! error message
  ! Local variables
  real(rkind) :: dt, t, max_dt, min_dt, deltat, debris_half_dt, div_q(nx,ny), dt_cfl, meanS
  real(rkind) :: gamma, m_dot(nx,ny), H(nx,ny), lat_rockfall(nx,ny)
  integer(i4b) :: emergenceMask(nx,ny)
  integer(i4b),parameter :: n=3 ! Glen's flow law exponent
  real(rkind),parameter :: A=2.4e-24 ! Modern Glen parameter
  real(rkind),parameter :: cfl= 0.124 ! Courant-Friedrichs-Lewy condition
  real(rkind) :: drivingStress(nx,ny)
  real(rkind) :: S_stepStart(nx,ny)
  real(rkind) :: Sklp(nx,ny), Sklm(nx,ny), Skl(nx,ny), Skpl(nx,ny), Skml(nx,ny)
  real(rkind) :: S_l_up(nx,ny), S_l_dn(nx,ny), S_k_up(nx,ny), S_k_dn(nx,ny)
  real(rkind) :: slope_l(nx,ny), slope_k(nx,ny), slope_g(nx,ny)
  integer(i4b) :: i, j, isteps
  integer(i4b) :: l(ny), lp(ny), lm(ny), lpp(ny), lmm(ny), k(nx), kp(nx), km(nx), kpp(nx), kmm(nx)

  err=0; message='run_flowModel/'
  gamma = 2._rkind * A * (iden_ice * gravity)**n / (n + 2_i4b)
  max_dt = 31._rkind * secprday ! max timestep in seconds, a month
  min_dt = 3600._rkind ! min timestep in seconds, 1 hour
  t = 0._rkind
  isteps = 0 ! counter for debugging print, only print once a max_dt
  volume = 0._rkind

  ! y direction indices
  l = [(i, i=1, ny)]
  lp = [(i, i=2, ny), ny]
  lpp = [(i, i=3, ny), ny, ny]
  lm = [1, (i, i=1, ny-1)]
  lmm = [1, 1, (i, i=1, ny-2)]
  if(ny == 1)then
    lpp = [1]
    lmm = [1]
  endif
  
  ! x direction indices
  k = [(i, i=1, nx)]
  kp = [(i, i=2, nx), nx]
  kpp = [(i, i=3, nx), nx, nx]
  km = [1, (i, i=1, nx-1)]
  kmm = [1, 1, (i, i=1, nx-2)]
  if(nx == 1)then
    kpp = [1]
    kmm = [1]
  endif

  ! one-time debris geometry setup from initial S and ELA for this yearly run
  if(sum(debris)>0._rkind)then
    call static_emergElev_latRockfall(S, B, ELA, latMoraineWidth, debrisConc, wallErosionRate, nx, ny, dx, dy, &
                                      emergenceMask, lat_rockfall, l, lp, lm, k, kp, km)
  else
    emergenceMask = 0_i4b
    lat_rockfall = 0._rkind
  endif

  ! time loop
  do while (t < t_total)
    dt = t_total - t
    H = merge(S-B, 0._rkind, glacierMask==1_i4b)
    ! debugging print
    if(printFlag)then 
      ! only print once a max_dt
      if(t > isteps*max_dt)then
        isteps = isteps + 1
        if(count(H>thick4area)>0_i4b)&
            write(*,'(a,f4.2,a,2(1x,f6.1),a,2(1x,f6.2))') "  time (yr) = ", t/secprday/365.25_rkind, " mean max H (m) =", &
                  sum(H)/count(H>thick4area), maxval(H), ", mean max debris (m) =", sum(debris)/count(H>thick4area), maxval(debris)
      endif
    endif
    if(count(H>thick4area)==0_i4b)then
      if(printFlag) write(*,'(a,f4.2,a)') "  time (yr) = ", t/secprday/365.25_rkind, " no glacier present"
      debris = 0._rkind
      S = B
      return
    endif

    ! get mass balance rate over surface (m s-1) and remove debris from glacier surface for flow calculation
    m_dot = massBalance(S, debris, glacierMask, slope, intercept, t_total, validElev, validCount, maxCount, nx, ny)
    S = S - debris
    S_stepStart = S

    ! call a step of the diffusion scheme
    call diffusion_MUSCL(S, B, glacierMask, gamma, n, cfl, max_dt, nx, ny, dx, dy, div_q, dt_cfl,&
                         l, lp, lm, lpp, lmm, k, kp, km, kpp, kmm)

    ! update time step
    deltat = min(dt_cfl, dt, max_dt)
    ! Apply a floor only when it does not violate CFL stability.
    if(dt_cfl >= min_dt .and. deltat < min_dt) deltat = min_dt
    debris_half_dt = 0.5_rkind*deltat
    t = t + deltat

    ! Strang splitting: half debris step on starting geometry if debris and ablation zone are present
    if(sum(debris)>0._rkind .and. any(m_dot<0._rkind))then
      call run_debrisModel(S_stepStart, B, debris, gamma, n, debris_half_dt, m_dot, emergenceMask, lat_rockfall,&
               debrisConc, theta_sat, iden_soil, nx, ny, dx, dy, l, lp, lm, k, kp, km)
    else
      debris = 0._rkind
    endif

    ! update S (with time discretization S_t+1 - S_t /dt = div(q_t)+ m_dot_t, so S_t+1 = S_t + (m_dot + div(q))*dt)
    S = S + (m_dot + div_q) * deltat

    ! check that the glacier is in boundaries, fix small violations, how small is arbitrary
    if(any((S - B) > verySmall .and. glacierMask==0_i4b))then
      if(any((S - B) > 10._rkind .and. glacierMask==0_i4b))then
        message=trim(message)//'glacier exceeds boundaries in flow model'
        err=20; return
      endif
      S = merge(B, S, (S - B) > verySmall .and. glacierMask==0_i4b)
    endif
    ! check that glacier surface is not infinite (unstable), bring down to mean glacier height
    if(any(((S - B) > 1.e6_rkind .or. isnan(S-B)) .and. glacierMask==1_i4b))then
      meanS = sum(merge(S, 0._rkind, glacierMask==1_i4b .and. (S-B)<1.e6_rkind)) / count(glacierMask==1_i4b .and. (S-B)<1.e6_rkind)
      S = merge(S, meanS, ((S - B)) <= 1.e6_rkind)
    endif

    ! refresh m_dot on the updated surface
    m_dot = massBalance(S + debris, debris, glacierMask, slope, intercept, t_total, validElev, validCount, maxCount, nx, ny)

    ! Strang splitting: half debris step on updated geometry if debris, ablation zone, and glacier are present
    if(sum(debris)>0._rkind .and. any(m_dot<0._rkind) .and. count(S-B>thick4area)>0_i4b)then
      call run_debrisModel(S, B, debris, gamma, n, debris_half_dt, m_dot, emergenceMask, lat_rockfall,&
               debrisConc, theta_sat, iden_soil, nx, ny, dx, dy, l, lp, lm, k, kp, km)

      ! calculate glacier ice surface slope with a finite difference, add debris back to surface for slope calculation
      Sklp  = S(k ,lp) + debris(k ,lp)
      Sklm  = S(k ,lm) + debris(k ,lm)
      Skl   = S(k ,l ) + debris(k ,l )
      S_l_up = midpt(Sklp, Skl)
      S_l_dn = midpt(Skl, Sklm)
      slope_l= (S_l_up - S_l_dn) / dy ! (m m-1)

      Skpl  = S(kp,l ) + debris(kp,l)
      Skml  = S(km,l ) + debris(km,l)
      Skl   = S(k ,l ) + debris(k ,l)
      S_k_up = midpt(Skpl, Skl)
      S_k_dn = midpt(Skl, Skml)
      slope_k= (S_k_up - S_k_dn) / dx ! (m m-1)

      ! remove debris downslope of where slope is too steep to hold debris, > yield stress following Mayer and Licciulli (2021)
      slope_g = merge(sqrt(slope_k**2_i4b + slope_l**2_i4b), 0._rkind, S>B) ! Slope in m/m, set to zero where no glacier
      drivingStress = (iden_soil*(1._rkind - theta_sat)) * gravity * debris * (abs(slope_g)/sqrt(slope_g**2_i4b + 1_i4b)) ! driving stress in Pa
      do j = 1, ny
        do i = 1, nx
          if(drivingStress(i,j) > debrisCritStress)then ! yield stress
            if(slope_k(i,j) > 0._rkind .and. slope_l(i,j) > 0._rkind)then
              debris(i:nx,j:ny) = 0._rkind
              exit
            elseif(slope_k(i,j) > 0._rkind .and. slope_l(i,j) < 0._rkind)then
              debris(i:nx, 1:j) = 0._rkind
              exit
            elseif(slope_l(i,j) < 0._rkind .and. slope_k(i,j) > 0._rkind)then
              debris(1:i ,j:ny) = 0._rkind
              exit
            elseif(slope_k(i,j) < 0._rkind .and. slope_k(i,j) < 0._rkind)then
              debris(1:i , 1:j) = 0._rkind
              exit
            endif
          endif
        enddo
      enddo
    else
      debris = 0._rkind
    endif

    ! add debris back to surface
    S = S + debris
  enddo ! end of time loop

  H = merge(S-B, 0._rkind, glacierMask==1_i4b)
  ! debugging print
  if(printFlag)then 
    if(count(H>thick4area)>0_i4b)&
        write(*,'(a,f4.2,a,2(1x,f6.1),a,2(1x,f6.2))') "  time (yr) = ", t/secprday/365.25_rkind, " mean max H (m) =", &
              sum(H)/count(H>thick4area), maxval(H), ", mean max debris (m) =", sum(debris)/count(H>thick4area), maxval(debris)
  endif
  if(count(H>thick4area)==0_i4b)then
    if(printFlag) write(*,'(a,f4.2,a)') "  time (yr) = ", t/secprday/365.25_rkind, " no glacier present"
    S = B
    return
  endif

  ! calculate volume of glacier (includes debris)
  volume = sum(S-B) * dx * dy * 1.e-9_rkind ! km3

end subroutine run_flowModel


! ************************************************************************************************
! subroutine diffusion_MUSCL: mass conserving scheme for SIA flow
! ************************************************************************************************
subroutine diffusion_MUSCL(S, B, mask, gamma, n, cfl, max_dt, nx, ny, dx, dy, div_q, dt_cfl, &
                          l, lp, lm, lpp, lmm, k, kp, km, kpp, kmm)
  implicit none
  ! Arguments
  real(rkind), intent(in) :: S(nx,ny), B(nx,ny), dx, dy, gamma, cfl, max_dt
  integer(i4b), intent(in) :: mask(nx,ny), nx, ny, n
  integer(i4b), intent(in) :: l(ny), lp(ny), lm(ny), lpp(ny), lmm(ny)
  integer(i4b), intent(in) :: k(nx), kp(nx), km(nx), kpp(nx), kmm(nx)
  real(rkind), intent(out) :: div_q(nx,ny), dt_cfl
  ! Local variables
  real(rkind) :: H(nx,ny)
  real(rkind) :: Sklp(nx,ny), Sklm(nx,ny), Skplp(nx,ny), Skplm(nx,ny), Skpl(nx,ny), Skl(nx,ny)
  real(rkind) :: Skmlp(nx,ny), Skmlm(nx,ny), Skml(nx,ny)
  real(rkind) :: Hkpl(nx,ny), Hkppl(nx,ny), Hkml(nx,ny), Hkmml(nx,ny), Hkl(nx,ny), Hklp(nx,ny)
  real(rkind) :: Hklpp(nx,ny), Hklm(nx,ny), Hklmm(nx,ny)
  real(rkind) :: H_l_up_m(nx,ny), H_l_up_p(nx,ny)
  real(rkind) :: H_l_dn_m(nx,ny), H_l_dn_p(nx,ny)
  real(rkind) :: f_l_p(nx,ny), f_l_m(nx,ny)
  real(rkind) :: D_l_up_m(nx,ny), D_l_up_p(nx,ny), D_l_up_min(nx,ny), D_l_up_max(nx,ny)
  real(rkind) :: D_l_dn_m(nx,ny), D_l_dn_p(nx,ny), D_l_dn_min(nx,ny), D_l_dn_max(nx,ny)
  real(rkind) :: D_l_up(nx,ny), D_l_dn(nx,ny)
  real(rkind) :: H_k_up_m(nx,ny), H_k_up_p(nx,ny), H_k_dn_m(nx,ny), H_k_dn_p(nx,ny)
  real(rkind) :: f_k_p(nx,ny), f_k_m(nx,ny)
  real(rkind) :: D_k_up_m(nx,ny), D_k_up_p(nx,ny), D_k_up_min(nx,ny), D_k_up_max(nx,ny)
  real(rkind) :: D_k_dn_m(nx,ny), D_k_dn_p(nx,ny), D_k_dn_min(nx,ny), D_k_dn_max(nx,ny)
  real(rkind) :: D_k_up(nx,ny), D_k_dn(nx,ny)
  real(rkind) :: divisor, div_k(nx,ny), div_l(nx,ny)
  integer(i4b) :: i, j

  ! Keep diffusion thickness nonnegative and zero outside glacier mask.
  H = merge(max(S - B, 0._rkind), 0._rkind, mask==1_i4b) ! H = ice thickness, S = surface height, B = bed topography

  ! indices
  Sklp  = S(k  ,lp )
  Sklm  = S(k  ,lm )
  Skplp = S(kp ,lp )
  Skplm = S(kp ,lm )
  Skpl  = S(kp ,l  )
  Skl   = S(k  ,l  )
  Skmlp = S(km ,lp )
  Skmlm = S(km ,lm )
  Skml  = S(km ,l  )
  !
  Hkpl  = H(kp ,l  )
  Hkppl = H(kpp,l  )
  Hkml  = H(km ,l  )
  Hkmml = H(kmm,l  )
  Hkl   = H(k  ,l  )
  Hklp  = H(k  ,lp )
  Hklpp = H(k  ,lpp)
  Hklm  = H(k  ,lm )
  Hklmm = H(k  ,lmm)

  ! calculate l+1/2 index
  H_l_up_m = minus(Hklm, Hkl, Hklp)
  H_l_up_p = pluss(Hkl, Hklp, Hklpp)
  H_l_up_m = max(H_l_up_m, 0._rkind)
  H_l_up_p = max(H_l_up_p, 0._rkind)

  ! calculate l-1/2 index
  H_l_dn_m = minus(Hklmm, Hklm, Hkl)
  H_l_dn_p = pluss(Hklm, Hkl, Hklp)
  H_l_dn_m = max(H_l_dn_m, 0._rkind)
  H_l_dn_p = max(H_l_dn_p, 0._rkind)

  ! calculate l flux
  f_l_p = flux(Skpl, Skml, Skplp, Skmlp, Sklp, Skl, dy, dx, n)
  f_l_m  = flux(Skpl, Skml, Skplm, Skmlm, Skl, Sklm, dy, dx, n)

  ! calculate l Diffusivity
  D_l_up_m = gamma * H_l_up_m**(n+2_i4b) * f_l_p  ! equation 30 Jarosh 2013
  D_l_up_p = gamma * H_l_up_p**(n+2_i4b) * f_l_p  ! equation 30 Jarosh 2013
  D_l_up_min = min(D_l_up_m, D_l_up_p)            ! equation 31 Jarosh 2013
  D_l_up_max = max(D_l_up_m, D_l_up_p)            ! equation 32 Jarosh 2013
  !
  D_l_dn_m = gamma * H_l_dn_m**(n+2_i4b) * f_l_m  ! equation 30 Jarosh 2013
  D_l_dn_p = gamma * H_l_dn_p**(n+2_i4b) * f_l_m  ! equation 30 Jarosh 2013
  D_l_dn_min = min(D_l_dn_m, D_l_dn_p)            ! equation 31 Jarosh 2013
  D_l_dn_max = max(D_l_dn_m, D_l_dn_p)            ! equation 32 Jarosh 2013

  ! equation 33 Jarosh 2013
  D_l_up = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if(Sklp(j,i) <= Skl(j,i) .and. H_l_up_m(j,i) <= H_l_up_p(j,i))then
        D_l_up(j,i) = D_l_up_min(j,i)
      else if(Sklp(j,i) <= Skl(j,i) .and. H_l_up_m(j,i) > H_l_up_p(j,i))then
        D_l_up(j,i) = D_l_up_max(j,i)
      else if(Sklp(j,i) > Skl(j,i) .and. H_l_up_m(j,i) <= H_l_up_p(j,i))then
        D_l_up(j,i) = D_l_up_max(j,i)
      else if(Sklp(j,i) > Skl(j,i) .and. H_l_up_m(j,i) > H_l_up_p(j,i))then
        D_l_up(j,i) = D_l_up_min(j,i)
      endif
    enddo
  enddo
  D_l_dn = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if(Skl(j,i) <= Sklm(j,i) .and. H_l_dn_m(j,i) <= H_l_dn_p(j,i))then
        D_l_dn(j,i) = D_l_dn_min(j,i)
      else if(Skl(j,i) <= Sklm(j,i) .and. H_l_dn_m(j,i) > H_l_dn_p(j,i))then
        D_l_dn(j,i) = D_l_dn_max(j,i)
      else if(Skl(j,i) > Sklm(j,i) .and. H_l_dn_m(j,i) <= H_l_dn_p(j,i))then
        D_l_dn(j,i) = D_l_dn_max(j,i)
      else if(Skl(j,i) > Sklm(j,i) .and. H_l_dn_m(j,i) > H_l_dn_p(j,i))then
        D_l_dn(j,i) = D_l_dn_min(j,i)
      endif
    enddo
  enddo
  ! enforce zero diffusion outside the mask
  D_l_up = merge(0._rkind, D_l_up, mask==0)
  D_l_dn = merge(0._rkind, D_l_dn, mask==0)

  ! calculate k+1/2 index
  H_k_up_m = minus(Hkml, Hkl, Hkpl)
  H_k_up_p = pluss(Hkl, Hkpl, Hkppl)
  H_k_up_m = max(H_k_up_m, 0._rkind)
  H_k_up_p = max(H_k_up_p, 0._rkind)

  ! calculate k-1/2 index
  H_k_dn_m = minus(Hkmml, Hkml, Hkl)
  H_k_dn_p = pluss(Hkml, Hkl, Hkpl)
  H_k_dn_m = max(H_k_dn_m, 0._rkind)
  H_k_dn_p = max(H_k_dn_p, 0._rkind)

  ! calculate k flux
  f_k_p = flux(Sklp, Sklm, Skplp, Skplm, Skpl, Skl, dx, dy, n)
  f_k_m = flux(Sklp, Sklm, Skmlp, Skmlm, Skl, Skml, dx, dy, n)

  ! calculate k Diffusivity
  D_k_up_m = gamma * H_k_up_m**(n+2_i4b) * f_k_p  ! equation 30 Jarosh 2013
  D_k_up_p = gamma * H_k_up_p**(n+2_i4b) * f_k_p  ! equation 30 Jarosh 2013
  D_k_up_min = min(D_k_up_m, D_k_up_p)            ! equation 31 Jarosh 2013
  D_k_up_max = max(D_k_up_m, D_k_up_p)            ! equation 32 Jarosh 2013
  !
  D_k_dn_m = gamma * H_k_dn_m**(n+2_i4b) * f_k_m  ! equation 30 Jarosh 2013
  D_k_dn_p = gamma * H_k_dn_p**(n+2_i4b) * f_k_m  ! equation 30 Jarosh 2013
  D_k_dn_min = min(D_k_dn_m, D_k_dn_p)            ! equation 31 Jarosh 2013
  D_k_dn_max = max(D_k_dn_m, D_k_dn_p)            ! equation 32 Jarosh 2013

  ! equation 33 Jarosh 2013
  D_k_up = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if(Skpl(j,i) <= Skl(j,i) .and. H_k_up_m(j,i) <= H_k_up_p(j,i))then
        D_k_up(j,i) = D_k_up_min(j,i)
      else if(Skpl(j,i) <= Skl(j,i) .and. H_k_up_m(j,i) > H_k_up_p(j,i))then
        D_k_up(j,i) = D_k_up_max(j,i)
      else if(Skpl(j,i) > Skl(j,i) .and. H_k_up_m(j,i) <= H_k_up_p(j,i))then
        D_k_up(j,i) = D_k_up_max(j,i)
      else if(Skpl(j,i) > Skl(j,i) .and. H_k_up_m(j,i) > H_k_up_p(j,i))then
        D_k_up(j,i) = D_k_up_min(j,i)
      endif
    enddo
  enddo
  D_k_dn = 0._rkind
  do j = 1, nx
    do i = 1, ny
      if(Skl(j,i) <= Skml(j,i) .and. H_k_dn_m(j,i) <= H_k_dn_p(j,i))then
        D_k_dn(j,i) = D_k_dn_min(j,i)
      else if(Skl(j,i) <= Skml(j,i) .and. H_k_dn_m(j,i) > H_k_dn_p(j,i))then
        D_k_dn(j,i) = D_k_dn_max(j,i)
      else if(Skl(j,i) > Skml(j,i) .and. H_k_dn_m(j,i) <= H_k_dn_p(j,i))then
        D_k_dn(j,i) = D_k_dn_max(j,i)
      else if(Skl(j,i) > Skml(j,i) .and. H_k_dn_m(j,i) > H_k_dn_p(j,i))then
        D_k_dn(j,i) = D_k_dn_min(j,i)
      endif
    enddo
  enddo
  ! enforce zero diffusion outside the mask
  D_k_up = merge(0._rkind, D_k_up, mask==0)
  D_k_dn = merge(0._rkind, D_k_dn, mask==0)

  ! calculate delta t and t
  divisor = max(maxval(abs(D_k_up)), maxval(abs(D_k_dn)), maxval(abs(D_l_up)), maxval(abs(D_l_dn)))
  if(divisor == 0._rkind)then
    dt_cfl = max_dt
  else
    dt_cfl = cfl * min(dy**2, dx**2) / divisor
  endif

  ! Calculate the time step values
  div_l = SIA(D_l_up, Sklp, Skl, D_l_dn, Sklm, dy) ! equation 36 Jarosh 2013
  div_k = SIA(D_k_up, Skpl, Skl, D_k_dn, Skml, dx) ! equation 36 Jarosh 2013
  div_q = div_k + div_l

end subroutine diffusion_MUSCL


! ************************************************************************************************
! private function mass balance: rate distribution after suface height changes
! ************************************************************************************************
function massBalance(S, debris, glacierMask, slope, intercept, t_total, validElev, validCount, maxCount, nx, ny)
  implicit none
  ! Arguments
  real(rkind), intent(in) :: S(nx,ny), debris(nx,ny), t_total
  integer(i4b), intent(in) :: glacierMask(nx,ny), maxCount, nx, ny, validCount(2)
  real(rkind), intent(in) :: slope(maxCount+1,2), intercept(maxCount+1,2), validElev(maxCount,2)
  ! Returns
  real(rkind) :: massBalance(nx,ny)
  ! Local variables
  integer(i4b) :: ind, dbr, i, j, k
  
  ! distribute mass balance over surface, using all points in GRU
  do k = 1, nx
    do j = 1, ny
      ! Set mass balance to zero if not on glacier
      if(glacierMask(k,j)==0_i4b)then
        massBalance(k,j) = 0._rkind
      else
        ! set debris
        dbr = 0
        if(debris(k,j) > 0._rkind) dbr = 1
        ! find the index of the elevation in the validElev array
        ind = validCount(dbr+1)+1 ! elevation above the valid elevation, extrapolate up
        do i = 1, validCount(dbr+1)
          if(S(k,j) <= validElev(i,dbr+1)) ind = i
        enddo
        massBalance(k,j) = (slope(ind,dbr+1) * S(k,j) + intercept(ind,dbr+1))/ t_total
      endif
    enddo
  enddo

  ! convert mass balance to m s-1 from kg m-2 s-1
  massBalance = massBalance / iden_water

end function massBalance


! ************************************************************************************************
! private subroutine run_debrisModel: set up englacial debris advection transport model and run
!   for one sub-time step. Debris input is from englacial debris emergence and lateral moraine rockfall,
!   and output from terminus slumping and ice velocity.
!   Follows the implementation of Mayer and Licciulli (2021), from Anderson and Anderson (2016).
! ************************************************************************************************
subroutine run_debrisModel(S, B, debris, gamma, n, t_total, m_dot, emergenceMask, lat_rockfall, &
                           debrisConc, theta_sat, iden_soil, nx, ny, dx, dy, l, lp, lm, k, kp, km)
  implicit none
  ! Arguments
  real(rkind), intent(in) :: S(nx,ny), B(nx,ny), gamma, lat_rockfall(nx,ny)
  real(rkind), intent(in) :: t_total, m_dot(nx,ny), debrisConc, theta_sat, iden_soil
  real(rkind), intent(in) :: dx, dy
  integer(i4b), intent(in) :: n, nx, ny, emergenceMask(nx,ny)
  integer(i4b), intent(in) :: l(ny), lp(ny), lm(ny), k(nx), kp(nx), km(nx)
  real(rkind), intent(inout) :: debris(nx,ny)
  ! Local variables
  real(rkind),parameter :: cfl= 0.5 ! Courant-Friedrichs-Lewy condition
  real(rkind) :: Sklp(nx,ny), Sklm(nx,ny), Skl(nx,ny), Skpl(nx,ny), Skml(nx,ny)
  real(rkind) :: S_l_up(nx,ny), S_l_dn(nx,ny), S_k_up(nx,ny), S_k_dn(nx,ny)
  real(rkind) :: slope_l(nx,ny), slope_k(nx,ny), u(nx,ny), v(nx,ny)
  real(rkind) :: influx(nx,ny), englacial_emerg(nx, ny)
  real(rkind) :: min_dt, div_uD(nx,ny), t, max_dt, dt, dt_cfl, deltat
  integer(i4b) :: mask(nx,ny)

  max_dt = 7._rkind * secprday ! maximum timestep in seconds, 1 week
  min_dt = 3600._rkind ! min timestep in seconds, 1 hour
  t = 0._rkind

  ! calculate glacier ice surface slope with a finite difference
  Sklp  = S(k ,lp)
  Sklm  = S(k ,lm) 
  Skl   = S(k ,l )
  S_l_up = midpt(Sklp, Skl)
  S_l_dn = midpt(Skl, Sklm)
  slope_l= (S_l_up - S_l_dn) / dy ! (m m-1)

  Skpl  = S(kp,l )
  Skml  = S(km,l )
  Skl   = S(k ,l )
  S_k_up = midpt(Skpl, Skl)
  S_k_dn = midpt(Skl, Skml)
  slope_k= (S_k_up - S_k_dn) / dx ! (m m-1)

  ! Use one-time masks/rockfall map computed from initial S for this yearly run.
  englacial_emerg = -merge(m_dot, 0._rkind, emergenceMask==1_i4b) * debrisConc ! (kg m-2 s-1)

  ! Sum for total debris influx in m/s
  influx = lat_rockfall + englacial_emerg/(iden_soil*(1._rkind - theta_sat)) 

  ! Calculate surface velocity using SIA, will be zero where no glacier
  u = gamma * (n + 2_i4b)/(n + 1_i4b) * abs(slope_l)**(n-1_i4b) * slope_l *(S-B)**(n + 1_i4b) ! (m s-1)
  v = gamma * (n + 2_i4b)/(n + 1_i4b) * abs(slope_k)**(n-1_i4b) * slope_k *(S-B)**(n + 1_i4b) ! (m s-1)

  mask = merge(1_i4b, 0_i4b, S>B) ! mask for glacier on subgrid
  
  ! ** movement with grid **
  do while (t < t_total)
    dt = t_total - t

    ! Call a step of the advection scheme
    call advection_MUSCL(u, v, debris, mask, cfl, max_dt, nx, ny, dx, dy, div_uD, dt_cfl, &
                         l, lp, lm, k, kp, km)
 
    ! update time step
    deltat = min(dt_cfl, dt)
    if(deltat > max_dt) deltat = max_dt
    if(deltat < min_dt) deltat = min_dt
    t = t + deltat 

    ! update debris thickness (with time discretization D_t+1 - D_t /dt  = div(uD_t) + influx_t, so D_t+1 = D_t + (div(uD_t) + influx)*dt)
    debris = debris + (div_uD + influx)* deltat
    debris = merge(debris, 0._rkind, debris > 0._rkind) ! prevent negative debris
  enddo

end subroutine run_debrisModel


! ************************************************************************************************
! private subroutine static_emergElev_latRockfall: build one-time debris geometry fields from initial 
!   yearly glacier surface S and ELA (for simplification, we assume these fields do not change over the year)
! ************************************************************************************************
subroutine static_emergElev_latRockfall(S, B, ELA, latMoraineWidth, debrisConc, wallErosionRate, nx, ny, dx, dy, &
                                        emergenceMask, lat_rockfall, l, lp, lm, k, kp, km)
  implicit none
  ! Arguments
  real(rkind), intent(in) :: S(nx,ny), B(nx,ny), ELA, latMoraineWidth, debrisConc, wallErosionRate, dx, dy
  integer(i4b), intent(in) :: nx, ny, l(ny), lp(ny), lm(ny), k(nx), kp(nx), km(nx)
  real(rkind), intent(out) :: lat_rockfall(nx,ny)
  integer(i4b), intent(out) :: emergenceMask(nx,ny)
  ! Local variables
  real(rkind) :: zeroMask(nx,ny), slope0(nx,ny), slope(nx,ny)
  real(rkind) :: distance(nx,ny), edgeMatrix(nx,ny)
  real(rkind) :: Sklp(nx,ny), Sklm(nx,ny), Skl(nx,ny), Skpl(nx,ny), Skml(nx,ny)
  real(rkind) :: S_l_up(nx,ny), S_l_dn(nx,ny), S_k_up(nx,ny), S_k_dn(nx,ny)
  real(rkind) :: slope_l(nx,ny), slope_k(nx,ny)
  real(rkind) :: rockface_len(nx,ny), rockface_len_local(nx,ny), sin_slope(nx,ny)
  real(rkind) :: local_drop, src_sum, mean_slope, topElev, emergenceElev
  integer(i4b) :: i, j, iter, nsrc, countGlac
  logical(lgt) :: updated
  logical(lgt) :: isEdge                ! .true. if a glacier cell borders a non-glacier cell

  emergenceMask = 0_i4b
  lat_rockfall = 0._rkind
  if (latMoraineWidth <= verySmall .and. debrisConc <= 0._rkind) return ! if no debris, skip rest of subroutine

  ! Glacier slopes for emergence elevation and rockface conversion
  Sklp  = S(k ,lp)
  Sklm  = S(k ,lm)
  Skl   = S(k ,l )
  S_l_up = midpt(Sklp, Skl)
  S_l_dn = midpt(Skl, Sklm)
  slope_l = (S_l_up - S_l_dn) / dy

  Skpl  = S(kp,l )
  Skml  = S(km,l )
  Skl   = S(k ,l )
  S_k_up = midpt(Skpl, Skl)
  S_k_dn = midpt(Skl, Skml)
  slope_k = (S_k_up - S_k_dn) / dx

  slope0 = sqrt(slope_k**2_i4b + slope_l**2_i4b)
  slope = merge(slope0, 0._rkind, S>B)

  ! debris-laden ice trajectories are assumed to follow the glacier surface slope, approximated as symmetric about the ELA
  if(debrisConc > 0._rkind)then
    countGlac = count(S>B)
    if(countGlac > 0_i4b)then
      mean_slope = sum(slope)/countGlac
      topElev = maxval(merge(S, 0._rkind, S>B))
      emergenceElev = 2._rkind*ELA - topElev + latMoraineWidth * mean_slope
      if(emergenceElev > ELA) emergenceElev = ELA
      emergenceMask = merge(1_i4b, 0_i4b, S<emergenceElev .and. S>B)
    endif
  else
    emergenceMask = 0_i4b
  endif

  ! Add rockfall along the sides of the glacier below the ELA
  if(latMoraineWidth>verySmall)then

    ! calculation of spatial mask distance to side for lateral moraine rockfall
    !  use a two-pass chamfer distance transform algorithm, approximate Euclidean distance
    zeroMask = merge(1_i4b, 0_i4b, S<=B+verySmall)
    distance = merge(0._rkind, 1.e6_rkind, zeroMask==1_i4b)

    ! First pass: top-left to bottom-right
    do i = 1, nx
      do j = 1, ny
        if(i > 1) distance(i, j) = min(distance(i, j), distance(i-1, j) + dx)
        if(j > 1) distance(i, j) = min(distance(i, j), distance(i, j-1) + dy)
      enddo
    enddo
    ! Second pass: bottom-right to top-left
    do i = nx, 1, -1
      do j = ny, 1, -1
        if(i < nx) distance(i, j) = min(distance(i, j), distance(i+1, j) + dx)
        if(j < ny) distance(i, j) = min(distance(i, j), distance(i, j+1) + dy)
      enddo
    enddo
    distance = merge(0._rkind, distance, zeroMask==1_i4b)

    ! Build edge matrix for rockfall, with -1 for non-glacier, 0 for glacier edge, and 1 for interior glacier
    edgeMatrix = merge(-1_i4b, 1_i4b, zeroMask==1_i4b)
    do i = 1, nx
      do j = 1, ny
        if(zeroMask(i,j)==0_i4b)then
          ! NOTE: .and. is not guaranteed to short-circuit in Fortran, so the bounds test
          !       and the array access have to be in separate if statements
          isEdge = .false.
          if(i>1 )then; if(zeroMask(i-1,j)==1_i4b) isEdge = .true.; endif
          if(i<nx)then; if(zeroMask(i+1,j)==1_i4b) isEdge = .true.; endif
          if(j>1 )then; if(zeroMask(i,j-1)==1_i4b) isEdge = .true.; endif
          if(j<ny)then; if(zeroMask(i,j+1)==1_i4b) isEdge = .true.; endif
          if(isEdge) edgeMatrix(i,j) = 0_i4b
        endif
      enddo
    enddo 

    ! get rockface length
    sin_slope = 0._rkind
    rockface_len = 0._rkind
    do j = 1, ny
      do i = 1, nx
        if(edgeMatrix(i,j)==0_i4b .and. S(i,j)<ELA .and. S(i,j)>B(i,j))then
          local_drop = 0._rkind
          if(i>1 )then; if(edgeMatrix(i-1,j)==-1_i4b) local_drop = max(local_drop, 2._rkind*(S(i,j) - B(i-1,j))); endif
          if(i<nx)then; if(edgeMatrix(i+1,j)==-1_i4b) local_drop = max(local_drop, 2._rkind*(S(i,j) - B(i+1,j))); endif
          if(j>1 )then; if(edgeMatrix(i,j-1)==-1_i4b) local_drop = max(local_drop, 2._rkind*(S(i,j) - B(i,j-1))); endif
          if(j<ny)then; if(edgeMatrix(i,j+1)==-1_i4b) local_drop = max(local_drop, 2._rkind*(S(i,j) - B(i,j+1))); endif
          if(local_drop>0._rkind .and. slope0(i,j)>verySmall)then
            sin_slope(i,j) = max(2._rkind*slope0(i,j)/sqrt((2._rkind*slope0(i,j))**2_i4b + 1._rkind), verySmall)
            rockface_len(i,j) = local_drop/sin_slope(i,j)
          endif
        endif
      enddo
    enddo

    ! Iteratively fill in rockface length for glacier points where rockfall can reach.
    rockface_len_local = rockface_len
    do iter = 1, nx + ny
      updated = .false.
      do j = 1, ny
        do i = 1, nx
          if(rockface_len_local(i,j)<=0._rkind .and. edgeMatrix(i,j)==1_i4b .and. distance(i,j)<=latMoraineWidth .and. &
             S(i,j)<ELA .and. S(i,j)>B(i,j))then
            src_sum = 0._rkind
            nsrc = 0
            if(i>1 )then
              if(distance(i-1,j)<distance(i,j) .and. rockface_len_local(i-1,j)>0._rkind)then
                src_sum = src_sum + rockface_len_local(i-1,j)
                nsrc = nsrc + 1
              endif
            endif
            if(i<nx)then
              if(distance(i+1,j)<distance(i,j) .and. rockface_len_local(i+1,j)>0._rkind)then
                src_sum = src_sum + rockface_len_local(i+1,j)
                nsrc = nsrc + 1
              endif
            endif
            if(j>1 )then
              if(distance(i,j-1)<distance(i,j) .and. rockface_len_local(i,j-1)>0._rkind)then
                src_sum = src_sum + rockface_len_local(i,j-1)
                nsrc = nsrc + 1
              endif
            endif
            if(j<ny)then
              if(distance(i,j+1)<distance(i,j) .and. rockface_len_local(i,j+1)>0._rkind)then
                src_sum = src_sum + rockface_len_local(i,j+1)
                nsrc = nsrc + 1
              endif
            endif
            if(nsrc>0)then
              rockface_len_local(i,j) = src_sum/real(nsrc, rkind)
              updated = .true.
            endif
          endif
        enddo
      enddo
      if(.not.updated) exit
    enddo
    lat_rockfall = merge(rockface_len_local*wallErosionRate/1000._rkind/secprday/365.25_rkind/latMoraineWidth, &
                      0._rkind, (edgeMatrix==0_i4b .or. (edgeMatrix==1_i4b .and. distance<=latMoraineWidth)) .and. S<ELA .and. S>B)
  else
    lat_rockfall = 0._rkind
  endif

end subroutine static_emergElev_latRockfall


! ************************************************************************************************
! private function advection_MUSCL: mass conserving advection scheme for debris movement with ice velocity
! ************************************************************************************************
subroutine advection_MUSCL(u, v, D, mask, cfl, max_dt, nx, ny, dx, dy, div_uD, dt_cfl,&
                          l, lp, lm, k, kp, km)
  implicit none
  ! Arguments
  real(rkind), intent(in) :: u(nx,ny), v(nx,ny), D(nx,ny), cfl, max_dt, dx, dy
  integer(i4b), intent(in) :: mask(nx,ny), nx, ny, l(ny), lp(ny), lm(ny), k(nx), kp(nx), km(nx)
  real(rkind), intent(out) :: div_uD(nx,ny), dt_cfl
  ! Local variables
  real(rkind) :: D_l_up(nx,ny), D_l_dn(nx,ny), D_k_up(nx,ny), D_k_dn(nx,ny)
  real(rkind) :: u_l_up(nx,ny), u_l_dn(nx,ny)
  real(rkind) :: v_k_up(nx,ny), v_k_dn(nx,ny)
  real(rkind) :: div_k(nx,ny), div_l(nx,ny)
  real(rkind) :: divisor, dt_u, dt_v

  ! MUSCL slope reconstruction for D, u, and v
  D_l_up = pluss(D(k,lm), D(k,l), D(k,lp)) ! D at l+1/2
  D_l_dn = minus(D(k,lm), D(k,l), D(k,lp)) ! D at l-1/2
  D_k_up = pluss(D(km,l), D(k,l), D(kp,l)) ! D at k+1/2
  D_k_dn = minus(D(km,l), D(k,l), D(kp,l)) ! D at k-1/2
  !
  u_l_up = pluss(u(k,lm), u(k,l), u(k,lp)) ! u at l+1/2
  u_l_dn = minus(u(k,lm), u(k,l), u(k,lp)) ! u at l-1/2
  !
  v_k_up = pluss(v(km,l), v(k,l), v(kp,l)) ! v at k+1/2
  v_k_dn = minus(v(km,l), v(k,l), v(kp,l)) ! v at k-1/2

  ! Calculate divergence of fluxes
  div_l = (u_l_dn * D_l_dn - u_l_up * D_l_up) / dy
  div_k = (v_k_dn * D_k_dn - v_k_up * D_k_up) / dx

  ! enforce zero advection outside the mask
  div_l = merge(0._rkind, div_l, mask==0_i4b)
  div_k = merge(0._rkind, div_k, mask==0_i4b)
  div_uD = div_k + div_l ! change in debris thickness with time

  ! calculate delta t from face speeds used in the advection fluxes
  dt_u = max_dt
  dt_v = max_dt
  divisor = max(maxval(abs(u_l_up)), maxval(abs(u_l_dn)))
  if(divisor > 0._rkind) dt_u = cfl * dy / divisor
  divisor = max(maxval(abs(v_k_up)), maxval(abs(v_k_dn)))
  if(divisor > 0._rkind) dt_v = cfl * dx / divisor
  dt_cfl = min(dt_u, dt_v)
  if(dt_cfl <= 0._rkind)then
    dt_cfl = max_dt
  endif

end subroutine advection_MUSCL

! ************************************************************************************************
! private functions to assist in calculations
! ************************************************************************************************
function superbee(r)
  implicit none
  real(rkind), intent(in) :: r(:,:)
  real(rkind), allocatable :: superbee(:,:)

  allocate(superbee(size(r,1), size(r,2)))
  superbee = max(0._rkind, min(2._rkind * r, 1._rkind), min(r, 2._rkind))
end function superbee

function flux(Skpl, Skml, Skplp, Skmlp, Sklp, Skl, dx, dy, n)
  implicit none
  real(rkind), intent(in) :: Skpl(:,:), Skml(:,:), Skplp(:,:), Skmlp(:,:), Sklp(:,:), Skl(:,:), dx, dy
  integer(i4b), intent(in) :: n
  real(rkind), allocatable :: flux(:,:)

  allocate(flux(size(Skpl,1), size(Skpl,2)))
  flux = ( (Skpl - Skml + Skplp - Skmlp)**2_i4b / (4._rkind * dx)**2_i4b & 
                 + (Sklp - Skl)**2_i4b / dy**2_i4b )**((n - 1.0) / 2._rkind)
end function flux

function SIA(Dup, Sup, S, Ddn, Sdn, d)
  implicit none
  real(rkind), intent(in) :: Dup(:,:), Sup(:,:), S(:,:), Ddn(:,:), Sdn(:,:), d
  real(rkind), allocatable :: SIA(:,:)

  allocate(SIA(size(Dup,1), size(Dup,2)))
  SIA = (Dup * (Sup - S) / d - Ddn * (S - Sdn) / d) / d
end function SIA

function midpt(h1, h2)
  implicit none
  real(rkind), intent(in) :: h1(:,:), h2(:,:)
  real(rkind), allocatable :: midpt(:,:)

  allocate(midpt(size(h1,1), size(h1,2)))
  midpt = 0.5_rkind * (h1 + h2)
end function midpt

function pluss(Hm, H, Hp)
  implicit none
  real(rkind), intent(in) :: Hm(:,:), H(:,:), Hp(:,:)
  real(rkind), allocatable :: pluss(:,:)
  logical, allocatable :: mask(:,:)
  real(rkind), allocatable :: divisor(:,:), ones(:,:)

  allocate(pluss(size(Hm,1), size(Hm,2)))
  allocate(mask(size(Hm,1), size(Hm,2)),divisor(size(Hm,1), size(Hm,2)),ones(size(Hm,1), size(Hm,2)))
  ones = 1.0_rkind
  mask = (Hp /= H) .and. (Hp /= Hm) .and. (H /= Hm)
  divisor = merge(Hp - H, ones, mask)
  pluss = merge(H - 0.5_rkind * superbee(abs((H - Hm) / divisor)) * (Hp - H), H, mask)
  deallocate(mask,divisor,ones)
end function pluss

function minus(Hm, H, Hp)
  implicit none
  real(rkind), intent(in) :: Hm(:,:), H(:,:), Hp(:,:)
  real(rkind), allocatable :: minus(:,:)
  logical, allocatable :: mask(:,:)
  real(rkind), allocatable :: divisor(:,:), ones(:,:)

  allocate(minus(size(Hm,1), size(Hm,2)))
  allocate(mask(size(Hm,1), size(Hm,2)),divisor(size(Hm,1), size(Hm,2)),ones(size(Hm,1), size(Hm,2)))
  ones = 1.0_rkind
  mask = (Hp /= H) .and. (Hp /= Hm) .and. (H /= Hm)
  divisor = merge(H - Hm, ones, mask)
  minus = merge(H + 0.5_rkind * superbee(abs((Hp - H) / divisor)) * (H - Hm), H, mask)
  deallocate(mask,divisor,ones)
end function minus


! ************************************************************************************************
! public subroutine time_updateGlacArea: find date to update glacier domain area, elevation, and layering 
! ************************************************************************************************
subroutine time_updateGlacArea(&
                         ! input
                         now_iyyy, now_im, now_id, now_ih, now_imin, & ! intent(in): current time
                         latitude,                   & ! intent(in): latitude of HRU (degrees, +N, -S)
                         ! output
                         updateJulDay,               & ! intent(inout): julian day of last glacier area update (fraction of day)
                         updateJulDayNext,           & ! intent(inout): julian day of next glacier area update (fraction of day)
                         updateGlacArea,             & ! intent(inout): flag to update glacier area this time step
                         sec_since_last_update,      & ! intent(out):   seconds since last glacier area update
                         ! error control
                         err, message)                 ! intent(out):   error control
  ! ----- define downstream subroutines -----------------------------------------------------------------------------------
  USE time_utils_module,only:compjulday                       ! convert calendar date to julian day
  USE time_utils_module,only:compcalday                       ! convert julian day to calendar date
  ! ----- define dummy variables ------------------------------------------------------------------------------------------
  implicit none
  integer(i4b), intent(in)        :: now_iyyy, now_im, now_id, now_ih, now_imin ! current time
  real(rkind), intent(in)         :: latitude                 ! latitude of HRU (degrees, +N, -S)
  real(rkind), intent(inout)      :: updateJulDay             ! julian day of last area update (fraction of day)
  real(rkind), intent(inout)      :: updateJulDayNext         ! julian day of next glacier area update (fraction of day)
  logical, intent(inout)          :: updateGlacArea           ! flag to update glacier area this time step
  real(rkind), intent(out)        :: sec_since_last_update    ! seconds since last glacier area update
  integer(i4b),intent(out)        :: err                      ! error code
  character(*),intent(out)        :: message                  ! error message 
  ! ----- local variables -------------------------------------------------------------------------------------------------
  real(rkind)                     :: currentJulDay            ! current julian day
  integer(i4b)                    :: lowMonth                 ! lowest mass balance month for updating glacier area
  integer(i4b)                    :: nextLowMonthYr           ! next year of lowMonth
  real(rkind)                     :: dsec                     ! seconds fraction of day
  integer(i4b)                    :: iyyy, im, id, ih, imin   ! year, month, day, hour, minute
  integer(i4b),parameter          :: nYears=1                 ! number of years in between glacier area updates (on first of lowMonth)
  character(LEN=256)              :: cmessage                 ! error message of downwind routine
  ! -----------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message='time_updateGlacArea/'

  call compjulday(now_iyyy, now_im, now_id, now_ih, now_imin,0._rkind,  & ! input  = year, month, day, hour, minute, second 
                  currentJulDay,err,cmessage)                            ! output = julian day (fraction of day) + error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  currentJulDay = currentJulDay + data_step/secprday ! julian day at end of time step

  ! determine lowMonth based on hemisphere, based on data from Dussaillant et al 2024
  if(latitude > 25._rkind)then
    lowMonth = 10_i4b ! October for Northern Hemisphere
  elseif(latitude < -25._rkind)then
    lowMonth = 4_i4b  ! April for Southern Hemisphere
  else
    lowMonth = 1_i4b  ! January for low latitudes (between 25N and 25S)
  end if

  ! get julian day of next update  (on first of lowMonth)
  if(updateJulDay==realMissing)then ! will only be true at the start of the simulation
    updateJulDay = dJulianStart
    if(now_im>=lowMonth)then ! start of simulation so now is dJulianStart
      nextLowMonthYr = now_iyyy + nYears
    else
      nextLowMonthYr = now_iyyy + nYears - 1_i4b
    endif
    call compjulday(nextLowMonthYr, lowMonth, 1_i4b, 0_i4b, 0_i4b, 0._rkind, updateJulDayNext,err,cmessage) 
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  elseif(updateJuldayNext==realMissing)then ! will only be true at the start of the simulation
    call compcalday(updateJulDay, iyyy, im, id, ih, imin, dsec, err, cmessage) ! get calendar date of last update
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    iyyy = iyyy + nYears
    call compjulday(iyyy, im, id, ih, imin, dsec, updateJuldayNext, err, cmessage) ! next update is nYears after last update
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  ! determine if glacier area needs to be updated this time step
  sec_since_last_update = (currentJulDay - updateJulDay)*secprday ! seconds since last update
  if(currentJulDay>=updateJuldayNext)then ! update glacier area, reset updateJulDay
    updateGlacArea = .true.
    updateJulDay = updateJulDayNext
    nextLowMonthYr = now_iyyy + nYears ! now_im is always lowMonth here
    call compjulday(nextLowMonthYr, lowMonth, 1_i4b, 0_i4b, 0_i4b, 0._rkind, updateJuldayNext, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

end subroutine time_updateGlacArea


! ************************************************************************************************
! public subroutine updateGlacDomain: update glacier domain area, elevation, and layering 
! ************************************************************************************************
subroutine updateGlacDomain(&
                  ! input
                  iglac,               & ! intent(inout): glacier domain index
                  glac_elev,           & ! intent(in):    elevation of each glacier domain (m) per HRU
                  glac_area,           & ! intent(in):    area of each glacier domain (m2)
                  glac_tan_slope,      & ! intent(in):    tan local ground surface slope of the domain (m/m)
                  glac_aspect,         & ! intent(in):    azimuth in degrees East of North of the domain (degrees)
                  glac_contourLength,  & ! intent(in):    length of contour at downslope edge of the domain (m)
                  glac_ablFrac,        & ! intent(in):    fraction of glacier domain area that is ablating
                  glac_debris_thick,   & ! intent(in):    debris thickness of each glacier domain (m) per HRU
                  dom_type,            & ! intent(in):    domain type
                  nSnow,               & ! intent(in):    number of snow layers
                  nLake,               & ! intent(in):    number of lake layers
                  nSoil,               & ! intent(in):    number of soil layers
                  nGlce,               & ! intent(in):    number of glacier ice layers
                  ! data structures
                  mpar_data,           & ! intent(in):    model parameters
                  indx_data,           & ! intent(in):    model indices
                  prog_data,           & ! intent(inout): model prognostic variables
                  diag_data,           & ! intent(inout): model diagnostic variables
                  flux_data,           & ! intent(inout): model fluxes
                  ! error control
                  err, message)         ! intent(out):   error control
  ! ----- define downstream subroutines -----------------------------------------------------------------------------------
  USE var_derive_module,only:rootDensty  ! module to calculate the vertical distribution of roots
  USE var_derive_module,only:satHydCond  ! module to calculate the saturated hydraulic conductivity in each soil layer
  ! ----- define dummy variables ------------------------------------------------------------------------------------------
  implicit none
  integer(i4b), intent(inout)     :: iglac                   ! glacier domain index
  real(rkind), intent(in)         :: glac_elev(:)             ! elevation of each glacier domain (m) per HRU
  real(rkind), intent(in)         :: glac_area(:)             ! area of each glacier domain (m2)
  real(rkind), intent(in)         :: glac_tan_slope(:)        ! tan local ground surface slope of the domain (m/m)
  real(rkind), intent(in)         :: glac_aspect(:)           ! azimuth in degrees East of North of the domain (degrees)
  real(rkind), intent(in)         :: glac_contourLength(:)    ! length of contour at downslope edge of the domain (m)
  real(rkind), intent(in)         :: glac_ablFrac(:)          ! fraction of glacier domain area that is ablating
  real(rkind), intent(in)         :: glac_debris_thick(:)     ! debris thickness of each glacier domain (m) per HRU
  integer(i4b), intent(in)        :: dom_type                 ! domain type
  integer(i4b), intent(in)        :: nSnow                    ! number of snow layers
  integer(i4b), intent(in)        :: nLake                    ! number of lake layers　(should be 0)
  integer(i4b), intent(in)        :: nSoil                    ! number of soil layers
  integer(i4b), intent(in)        :: nGlce                    ! number of glacier ice layers
  type(var_dlength),intent(in)    :: mpar_data                ! model parameters
  type(var_ilength),intent(in)    :: indx_data                ! model indices
  type(var_dlength),intent(inout) :: prog_data                ! model prognostic variables
  type(var_dlength),intent(inout) :: diag_data                ! model diagnostic variables
  type(var_dlength),intent(inout) :: flux_data                ! model fluxes
  integer(i4b),intent(out)        :: err                      ! error code
  character(*),intent(out)        :: message                  ! error message 
   ! ----- define local variables ------------------------------------------------------------------------------------------
  integer(i4b)                    :: i                        ! loop index
  integer(i4b)                    :: nLayers                  ! total number of layers
  real(rkind)                     :: layers_thick             ! depth of layers modifying
  real(rkind)                     :: thick_ratio              ! ratio of new layers thickness to previous thickness
  character(len=256)              :: cmessage                 ! error message
  ! ----------------------------------------------------------------------------------
  ! associate variables in data structure
  associate(&
   ! coordinate variables
   mLayerDepth          => prog_data%var(iLookPROG%mLayerDepth)%dat,            & ! thickness of each layer (m)       
   mLayerHeight         => prog_data%var(iLookPROG%mLayerHeight)%dat,           & ! height at the mid-point of each layer (m)
   iLayerHeight         => prog_data%var(iLookPROG%iLayerHeight)%dat,           & ! height at the interface of each layer (m)
    ! glacier domain variables
   glacMass4AreaChange  => prog_data%var(iLookPROG%glacMass4AreaChange)%dat(1), & ! mass change (kg m-2)
   DOMarea              => prog_data%var(iLookPROG%DOMarea)%dat(1),             & ! area of each domain (m2)
   DOMelev              => prog_data%var(iLookPROG%DOMelev)%dat(1),             & ! elevation of each domain (m)
   DOMtan_slope         => prog_data%var(iLookPROG%DOMtan_slope)%dat(1),        & ! tan local ground surface slope of the domain (m/m)
   DOMaspect            => prog_data%var(iLookPROG%DOMaspect)%dat(1),           & ! azimuth in degrees East of North of the domain (degrees)
   DOMcontourLength     => prog_data%var(iLookPROG%DOMcontourLength)%dat(1),    & ! length of contour at downslope edge of the domain (m)
   scalarAblFrac        => prog_data%var(iLookPROG%scalarAblFrac)%dat(1)        & ! fraction of glacier domain area that is ablating
   ) ! end associate
   ! ----------------------------------------------------------------------------------------------
   ! initialize
   err=0; message='updateGlacDomain/'
   nLayers = nSnow + nLake + nSoil + nGlce

   ! update glacier domain elevation, area, and ablating fraction and reset mass change
   DOMelev = glac_elev(iglac) ! realMissing if no area
   DOMarea = glac_area(iglac) ! may be 0
   DOMtan_slope = glac_tan_slope(iglac) ! realMissing if no area
   DOMaspect = glac_aspect(iglac) ! realMissing if no area
   DOMcontourLength = glac_contourLength(iglac) ! may be 0
   scalarAblFrac = glac_ablFrac(iglac) ! may be 0
   glacMass4AreaChange = 0._rkind ! reset

   ! update glacier layering
   if(dom_type==glacDbr .and. DOMarea>0._rkind)then ! thickness of average debris cover in HRU changes with debris advection
     ! for zero area, glacDbr needs to keep previous non-zero soil layering so can adjust if debris-covered glacier re-appears
     ! scale soil layer thickness with debris thickness change, keep the same number of layers 
     layers_thick = sum(mLayerDepth(nSnow+nLake+1:nSnow+nLake+nSoil))
     thick_ratio = glac_debris_thick(iglac)/layers_thick ! glac_debris_thick will be >0 since there is a positive area of debris-covered glacier
     mLayerDepth(nSnow+nLake+1:nSnow+nLake+nSoil)  = mLayerDepth(nSnow+nLake+1:nSnow+nLake+nSoil)*thick_ratio
     mLayerHeight(nSnow+nLake+1:nSnow+nLake+nSoil) = mLayerHeight(nSnow+nLake+1:nSnow+nLake+nSoil)*thick_ratio
     iLayerHeight(nSnow+nLake+1:nSnow+nLake+nSoil) = iLayerHeight(nSnow+nLake+1:nSnow+nLake+nSoil)*thick_ratio            

     ! recalculate the layer heights below soil
     do i=nSnow+nLake+nSoil+1,nLayers
       mLayerHeight(i) = mLayerHeight(i) + glac_debris_thick(iglac) - layers_thick
       iLayerHeight(i) = iLayerHeight(i) + glac_debris_thick(iglac) - layers_thick
     enddo

     ! recalculate vertical distribution of root density
     call rootDensty(mpar_data,     & ! intent(in):    model parameters
                     indx_data,     & ! intent(in):    model indices
                     prog_data,     & ! intent(in):    model prognostic variables
                     diag_data,     & ! intent(inout): model diagnostic variables
                     err,cmessage)    ! error control
     if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

     ! recalculate saturated hydraulic conductivity in each soil layer
     call satHydCond(mpar_data,    & ! intent(in):    model parameters
                     indx_data,    & ! intent(in):    model indices
                     prog_data,    & ! intent(in):    model prognostic variables
                     flux_data,    & ! intent(inout): model fluxes
                     err,cmessage)   ! error control
     if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   elseif(DOMarea==0._rkind .and. nSnow>0)then ! no glacier area but still has snow layers, so make snow very thin
     ! scale snow layers so comes back with tiny snow (would prefer 0 thickness but then would have to relayer here)
     layers_thick = sum(mLayerDepth(1:nSnow))
     if(layers_thick>verySmall)then
       thick_ratio = verySmall/layers_thick
       mLayerDepth(1:nSnow) = mLayerDepth(1:nSnow)*thick_ratio
       mLayerHeight(1:nSnow) = mLayerHeight(1:nSnow)*thick_ratio
       iLayerHeight(1:nSnow) = iLayerHeight(1:nSnow)*thick_ratio          

       ! recalculate the layer heights below snow
       do i=nSnow+1,nLayers
         mLayerHeight(i) = mLayerHeight(i) + verySmall - layers_thick
         iLayerHeight(i) = iLayerHeight(i) + verySmall - layers_thick
       enddo
     endif

   endif ! ( if changing soil layers or snow layers)
  
 end associate
end subroutine updateGlacDomain


end module glacAreaChange_module
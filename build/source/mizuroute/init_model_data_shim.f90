!-----------------------------------------------------------------------
! mizuRoute initialization compatibility shim
!
! This module provides the subset of the original mizuRoute model
! initialization interface required by FUSE. The routines contained
! here are retained unchanged from mizuRoute so that the routing-method
! objects can be initialized without introducing a dependency on the
! full mizuRoute model initialization infrastructure.
!-----------------------------------------------------------------------
MODULE init_model_data_shim

! data types
USE nrtype,    ONLY: i4b,dp,lgt,strLen
USE dataTypes, ONLY: var_ilength         ! integer type:          var(:)%dat
USE dataTypes, ONLY: var_clength         ! integer type:          var(:)%dat
USE dataTypes, ONLY: var_dlength         ! double precision type: var(:)%dat, or dat

USE var_lookup, ONLY: ixNTOPO            ! index of variables for the network topology
USE var_lookup, ONLY: ixHRU2SEG          ! index of variables for data structure
USE var_lookup, ONLY: ixPFAF             ! index of variables for the pfafstetter code

! Shared data
USE public_var, ONLY: iulog
USE public_var, ONLY: charMissing

implicit none

private
!public :: get_mpi_omp
!public :: init_model
!public :: init_ntopo_data
!public :: init_state_data
!public :: init_qmod
!public :: update_time
!public :: init_pio
!
public :: init_ntopo
public :: init_route_method

CONTAINS

 ! *********************************************************************
 ! public subroutine: initialize river network data
 ! *********************************************************************
 SUBROUTINE init_ntopo(instance_rank,                                                & ! input:  model instance rank
                       nHRU_out, nRch_out,                                           & ! output: number of HRU and Reaches
                       structHRU, structSEG, structHRU2SEG, structNTOPO, structPFAF, & ! output: data structure for river data
                       ierr, message)                                                  ! output: error controls
  ! Shared data
  USE public_var, ONLY: ancil_dir                ! name of the ancillary directory
  USE public_var, ONLY: fname_ntopOld            ! name of the old network topology file
  USE public_var, ONLY: fname_ntopNew            ! name of the new network topology file
  USE public_var, ONLY: dname_nhru               ! dimension name for HRUs
  USE public_var, ONLY: dname_sseg               ! dimension name for stream segments
  USE public_var, ONLY: maxPfafLen               ! maximum digit of pfafstetter code (default 32)
  ! options
  USE public_var, ONLY: ntopAugmentMode          ! River network augmentation mode
  USE public_var, ONLY: idSegOut                 ! River network subset mode (idSegOut > 0)
  ! global data
  USE globalData, ONLY: meta_PFAF                ! meta for pfafstetter code
  ! external subroutines
  USE read_streamSeg,       ONLY: getData                  ! get the ancillary data
  USE write_streamSeg,      ONLY: writeData                ! write the ancillary data
  USE process_ntopo,        ONLY: check_river_properties   ! check if river network data is physically valid
  USE ncio_utils,           ONLY: get_var_dims
  USE process_ntopo,        ONLY: augment_ntopo            ! compute all the additional network topology (only compute option = on)

  implicit none
  ! dummy variables
  integer(i4b)                  , intent(in)  :: instance_rank            ! rank of model instance
  integer(i4b)                  , intent(out) :: nHRU_out                 ! number of HRUs
  integer(i4b)                  , intent(out) :: nRch_out                 ! number of reaches
  type(var_dlength), allocatable, intent(out) :: structHRU(:)             ! HRU properties
  type(var_dlength), allocatable, intent(out) :: structSeg(:)             ! stream segment properties
  type(var_ilength), allocatable, intent(out) :: structHRU2SEG(:)         ! HRU-to-segment mapping
  type(var_ilength), allocatable, intent(out) :: structNTOPO(:)           ! network topology
  type(var_clength), allocatable, intent(out) :: structPFAF(:)            ! pfafstetter code
  integer(i4b)      , intent(out)             :: ierr                     ! error code
  character(*)      , intent(out)             :: message                  ! error message
  ! Local variables
  integer(i4b)                                :: tot_upstream             ! total number of all of the upstream stream segments for all stream segments
  integer(i4b)                                :: tot_upseg                ! total number of immediate upstream segments for all  stream segments
  integer(i4b)                                :: tot_hru                  ! total number of all the upstream hrus for all stream segments
  integer(i4b)                                :: tot_uh                   ! total number of unit hydrograph from all the stream segments
  integer(i4b),      allocatable              :: ixHRU_desired(:)         ! indices of desired hrus
  integer(i4b),      allocatable              :: ixSeg_desired(:)         ! indices of desired reaches
  integer(i4b)                                :: dummy(2)                 ! dummy variable to hold dimension length for 2D variables in netCDF
  integer(i4b)   , parameter                  :: maxUpstreamFile=90000000 ! 90 million: maximum number of upstream reaches to enable writing
  character(len=strLen)                       :: cmessage                 ! error message of downwind routine

  ierr=0; message='init_ntopo/'

  ! get the variable dimensions
  ! NOTE: need to update maxPfafLen to the exact character size for pfaf code in netCDF
  if (meta_PFAF(ixPFAF%code)%varFile) then
    call get_var_dims(trim(ancil_dir)//trim(fname_ntopOld), & ! input: file name
                      trim(meta_PFAF(ixPFAF%code)%varName), & ! input: pfaf code variable name in netcdf
                      ierr, cmessage,                       & ! output: error control
                      dlen=dummy)                             ! output optional: dimension length
    if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif
    maxPfafLen = dummy(1)
  end if

  call getData(trim(ancil_dir)//trim(fname_ntopOld), & ! input: file name
               dname_nhru,    & ! input: dimension name of the HRUs
               dname_sseg,    & ! input: dimension name of the stream segments
               nHRU_out,      & ! output: number of HRUs
               nRch_out,      & ! output: number of stream segments
               structHRU,     & ! output: ancillary data for HRUs
               structSeg,     & ! output: ancillary data for stream segments
               structHRU2seg, & ! output: ancillary data for mapping hru2basin
               structNTOPO,   & ! output: ancillary data for network topology
               structPFAF,    & ! output: ancillary data for pfafstetter code
               ierr,cmessage)   ! output: error control
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  call check_river_properties(structNTOPO, structHRU, structSEG, ierr, cmessage) ! input: data structure for physical river network data
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  call augment_ntopo(nHRU_out,                         & ! number of HRUs
                     nRch_out,                         & ! number of stream segments
                     structHRU,                        & ! ancillary data for HRUs
                     structSeg,                        & ! ancillary data for stream segments
                     structHRU2seg,                    & ! ancillary data for mapping hru2basin
                     structNTOPO,                      & ! ancillary data for network toopology
                     ierr, cmessage,                   & ! error control
                     tot_hru       = tot_hru,          & ! total number of all the upstream hrus for all stream segments
                     tot_upseg     = tot_upseg,        & ! total number of all the immediate upstream segments for all stream segments
                     tot_upstream  = tot_upstream,     & ! total number of all the upstream segments for all stream segments
                     tot_uh        = tot_uh,           & ! total number of unit hydrograph for all stream segments
                     ixHRU_desired = ixHRU_desired,    & ! indices of desired hrus
                     ixSeg_desired = ixSeg_desired)      ! indices of desired reaches
  if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! write network topology (if augment mode or subset mode)
  if(ntopAugmentMode .or. idSegOut>0)then

    ! disable the dimension containing all upstream reaches
    ! NOTE: For the CONUS this is 1,872,516,819 reaches !!
    !        --> it will always be quicker to recompute than read+write
    !        --> users can modify the hard-coded parameter "maxUpstreamFile" if desired
    if(tot_upstream > maxUpstreamFile) tot_upstream=0

    ! only write if model instance rank is zero
    if (instance_rank == 0)then

      call writeData(trim(ancil_dir)//trim(fname_ntopNew), & ! input: file name
                     tot_hru,       & ! input: total number of all the upstream hrus for all stream segments
                     tot_upseg,     & ! input: total number of immediate upstream segments for all  stream segments
                     tot_upstream,  & ! input: total number of all of the upstream stream segments for all stream segments
                     tot_uh,        & ! input: total number of unit hydrograph for all stream segments
                     ixHRU_desired, & ! input: indices of desired hrus
                     ixSeg_desired, & ! input: indices of desired reaches
                     structHRU,     & ! input: ancillary data for HRUs
                     structSeg,     & ! input: ancillary data for stream segments
                     structHRU2seg, & ! input: ancillary data for mapping hru2basin
                     structNTOPO,   & ! input: ancillary data for network topology
                     structPFAF,    & ! input: ancillary data for pfafstetter code
                     ierr,cmessage) ! output: error control
      if(ierr/=0)then; message=trim(message)//trim(cmessage); return; endif

    endif ! (if model instance rank is zero)

    if (idSegOut>0) write(iulog,'(a)') 'Running in river network subset mode'
    if (ntopAugmentMode) write(iulog,'(a)') 'Running in river network augmentation mode'
    write(iulog,'(a)') 'Created a new network topology file '//trim(fname_ntopNew)
    write(iulog,'(a)') ' --> Run again using the new network topology file '
    return
  endif

 END SUBROUTINE init_ntopo

 ! *********************************************************************
 ! public subroutine: initialize routing method object
 ! *********************************************************************
 SUBROUTINE init_route_method(ierr, message)
   !
   ! DESCRIPTION:
   ! Instantiate a collection of routing method objects
   !
   USE globalData,         ONLY: rch_routes            ! routing methods instantiated
   USE globalData,         ONLY: routeMethods          ! Active routing method
   USE public_var,         ONLY: accumRunoff           ! routing method ID
   USE public_var,         ONLY: impulseResponseFunc   ! routing method ID
   USE public_var,         ONLY: kinematicWaveTracking ! routing method ID
   USE public_var,         ONLY: kinematicWave         ! routing method ID
   USE public_var,         ONLY: muskingumCunge        ! routing method ID
   USE public_var,         ONLY: diffusiveWave         ! routing method ID
   USE accum_runoff_module,ONLY: accum_runoff_rch      ! routing routine: accumulation instantaneous runoff
   USE irf_route_module,   ONLY: irf_route_rch         ! routing routine: Impulse response function
   USE kwt_route_module,   ONLY: kwt_route_rch         ! routing routine: Lagrangian kinematic
   USE kw_route_module,    ONLY: kwe_route_rch         ! routing routine: kinematic
   USE mc_route_module,    ONLY: mc_route_rch          ! routing routine: muskingum
   USE dfw_route_module,   ONLY: dfw_route_rch         ! routing routine: diffusive

   implicit none
   ! Argument variables:
   integer(i4b),          intent(out)   :: ierr        ! error code
   character(*),          intent(out)   :: message     ! error message
   ! Local variables:
   character(len=strLen)                :: cmessage     ! error message from subroutine
   integer(i4b)                         :: ix

   ierr=0; message='init_route_method/'

   allocate(rch_routes(size(routeMethods)), stat=ierr, errmsg=cmessage)
   if(ierr/=0)then; message=trim(message)//trim(cmessage)//' [rch_routes]'; return; endif

   do ix=1, size(routeMethods)
     select case (routeMethods(ix))
       case (accumRunoff)
         allocate(accum_runoff_rch :: rch_routes(ix)%rch_route)
       case (impulseResponseFunc)
         allocate(irf_route_rch :: rch_routes(ix)%rch_route)
       case (kinematicWaveTracking)
         allocate(kwt_route_rch :: rch_routes(ix)%rch_route)
       case (kinematicWave)
         allocate(kwe_route_rch :: rch_routes(ix)%rch_route)
       case (muskingumCunge)
         allocate(mc_route_rch  :: rch_routes(ix)%rch_route)
       case (diffusiveWave)
         allocate(dfw_route_rch :: rch_routes(ix)%rch_route)
       case default
         ierr=20; message=trim(message)//'no valid routing method'; return
     end select
   end do

 END SUBROUTINE init_route_method

END MODULE init_model_data_shim

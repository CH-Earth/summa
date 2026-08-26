module mizuroute_types

  use nrtype, only: wp, i4b, lgt

  use dataTypes, only: mizu_var_dlength => var_dlength
  use dataTypes, only: mizu_var_ilength => var_ilength
  use dataTypes, only: mizu_var_clength => var_clength
  use dataTypes, only: mizu_remap       => remap
  use dataTypes, only: mizu_runoff      => runoff

  use dataTypes, only: mizu_RCHPRP      => RCHPRP
  use dataTypes, only: mizu_RCHTOPO     => RCHTOPO

  use dataTypes, ONLY: mizu_STRFLX      => STRFLX
  use dataTypes, ONLY: mizu_STRSTA      => STRSTA

  implicit none
  private

  public :: rout_info
  public :: topo_info
  public :: remap_info
  public :: mizuroute_info

  public :: routing_time_data
  public :: mizuroute_topology
  public :: river_network_data
  public :: spatial_remap_data
  public :: mizuroute_domain

  !---------------------------------------------------------------------
  ! mizuRoute configuration
  !---------------------------------------------------------------------
  
  !---------------------------------------------------------------------

  type :: rout_info

   character(len=:), allocatable :: namelist_path       ! namelist path
   character(len=:), allocatable :: namelist_file       ! namelist file

   character(len=:), allocatable :: methods             ! string of integers defining methods

   real(wp)                      :: dt = 3600._wp       ! routing time step (s)

  end type rout_info

  !---------------------------------------------------------------------

  type :: topo_info

   ! Hydrofabric path/filenames
   character(len=:), allocatable :: hfabric_path       ! hydrofabric path
   character(len=:), allocatable :: hfabric_file       ! hydrofabric file
   character(len=:), allocatable :: hfabric_newfile    ! hydrofabric file (new)

   ! NetCDF dimensions
   character(len=:), allocatable :: dname_hru          ! dimension name: hru ID
   character(len=:), allocatable :: dname_seg          ! dimension name: segment

   ! NetCDF variable names
   character(len=:), allocatable :: varname_HRUid      ! variable name: HRU ID
   character(len=:), allocatable :: varname_segId      ! variable name: segment ID
   character(len=:), allocatable :: varname_hruSegId   ! variable name: ID of segment in HRU
   character(len=:), allocatable :: varname_downSegId  ! variable name: downstream segment ID

   character(len=:), allocatable :: varname_area       ! variable name: HRU area
   character(len=:), allocatable :: varname_slope      ! variable name: segment slope
   character(len=:), allocatable :: varname_length     ! variable name: segment length

   ! Network topology
   integer(i4b)                  :: idSegOut = -9999   ! ID of outlet segment
   integer(i4b)                  :: ixSegOut = -9999   ! index of outlet segment

  end type topo_info

  !---------------------------------------------------------------------
  
  type :: remap_info

   ! Remapping filename
   character(len=:), allocatable :: remap_file         ! remapping file

   ! NetCDF dimensions
   character(len=:), allocatable :: dname_hru          ! name of dimension of river network HRU ID
   character(len=:), allocatable :: dname_data         ! name of dimension of runoff HRU overlapping with river network HRU

   ! NetCDF variable names
   character(len=:), allocatable :: vname_hruid        ! name of variable containing ID of river network HRU
   character(len=:), allocatable :: vname_weight       ! name of variable contating areal weights of runoff HRUs within each river network HRU
   character(len=:), allocatable :: vname_num_qhru     ! name of variable containing numbers of runoff HRUs within each river network HRU
   character(len=:), allocatable :: vname_i_index      ! name of variable containing index of xlon dimension in runoff grid (if runoff file is grid)
   character(len=:), allocatable :: vname_j_index      ! name of variable containing index of ylat dimension in runoff grid (if runoff file is grid)

  end type remap_info

  !---------------------------------------------------------------------
  
  type :: mizuroute_info
  
   type(rout_info)  :: mrout
   type(topo_info)  :: ntopo
   type(remap_info) :: remap
  
   integer(i4b) :: n_hru
   integer(i4b) :: n_seg
  
   real(wp)     :: dt_landmodel
  
   logical(lgt) :: is_gridded

   logical(lgt) :: do_mizuroute = .false.
   logical(lgt) :: do_remapping = .false.
   logical(lgt) :: is_print     = .false.

  end type mizuroute_info

  !---------------------------------------------------------------------
  ! Time information
  !---------------------------------------------------------------------
  type :: routing_time_data
  
    integer(i4b) :: n_sub  = 1
    real(wp)     :: dt_sub = 0._wp
  
  end type routing_time_data

  !---------------------------------------------------------------------
  ! River-network topology and attributes
  !---------------------------------------------------------------------
  type :: mizuroute_topology

    integer(i4b) :: n_hru = 0
    integer(i4b) :: n_seg = 0

    type(mizu_var_dlength), allocatable :: hru(:)
    type(mizu_var_dlength), allocatable :: seg(:)
    type(mizu_var_ilength), allocatable :: hru2seg(:)
    type(mizu_var_ilength), allocatable :: ntopo(:)
    type(mizu_var_clength), allocatable :: pfaf(:)

    logical(lgt) :: is_initialized = .false.

  end type mizuroute_topology

  !---------------------------------------------------------------------
  ! FUSE-step mean streamflow for each reach and time step [m3/s] 
  !---------------------------------------------------------------------
  type :: fusestep_mean
  
    real(wp),               allocatable :: streamflow(:,:) 
  
  end type fusestep_mean

  !---------------------------------------------------------------------
  ! Information to couple with a host land model
  !---------------------------------------------------------------------
  type :: reach_data

    ! IDs for HRU and stream segments
    integer(i4b), allocatable :: hru_id(:)
    integer(i4b), allocatable :: seg_id(:)

    ! Reach properties needed for FUSE coupling
    real(wp)    , allocatable :: totArea(:)

  end type reach_data

  !---------------------------------------------------------------------
  ! Persistent river-network data
  !---------------------------------------------------------------------
  type :: river_network_data

    type(mizuroute_topology)         :: topology  ! static network topology and attributes
    type(mizu_runoff)                :: runoff    ! FUSE runoff in mizuRoute structures

    ! mizuRoute routing: reach properties and network topology
    type(mizu_RCHPRP),   allocatable :: param(:)  ! reach properties
    type(mizu_RCHTOPO),  allocatable :: ntopo(:)  ! network topology

    ! mizuRoute routing state and fluxes
    type(mizu_STRSTA),   allocatable :: state(:)  ! model states
    type(mizu_STRFLX),   allocatable :: flux(:)   ! model fluxes

    ! coupling workspace
    real(wp),            allocatable :: reach_inflow(:)  ! lateral inflow to each reach [m3/s]

    ! outputs for each routing method
    type(fusestep_mean), allocatable :: method(:)

    ! time data for routing substeps
    type(routing_time_data)          :: time

  end type river_network_data

  !---------------------------------------------------------------------
  ! Spatial mappings between model discretizations
  !---------------------------------------------------------------------
  type :: spatial_remap_data

    type(mizu_remap) :: forcing
    type(mizu_remap) :: routing

  end type spatial_remap_data

  !---------------------------------------------------------------------
  ! Combined domain structures for mizuroute
  !---------------------------------------------------------------------

  type :: mizuroute_domain

    type(river_network_data) :: river_network
    type(reach_data)         :: reach
    type(spatial_remap_data) :: remap

  end type mizuroute_domain


end module mizuroute_types

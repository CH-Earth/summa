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
MODULE summa_type

! used to define the top-level summa data structure

! *****************************************************************************
! * higher-level derived data types
! *****************************************************************************
USE nr_type                                ! variable types, etc.
USE iso_fortran_env, only: output_unit     ! output unit (normally=6)

! general summa data types
USE data_types,  only : &
                    ! no spatial dimension
                    var_i,                 & ! x%var(:)            (i4b)
                    var_d,                 & ! x%var(:)            (dp)
                    var_ilength,           & ! x%var(:)%dat        (i4b)
                    var_dlength,           & ! x%var(:)%dat        (dp)
                    ! no variable dimension
                    hru_i,                 & ! x%hru(:)            (i4b)
                    hru_d,                 & ! x%hru(:)            (dp)
                    gru_hru_i,             & ! x%gru(:)%hru(:)     (i4b)
                    gru_hru_d,             & ! x%gru(:)%hru(:)     (dp)
                    gru_hru_dom_d,         & ! x%gru(:)%hru(:)%dom(:) (dp)
                    ! gru dimension
                    gru_int,               & ! x%gru(:)%var(:)     (i4b)
                    gru_double,            & ! x%gru(:)%var(:)     (dp)
                    gru_intVec,            & ! x%gru(:)%var(:)%dat (i4b)
                    gru_doubleVec,         & ! x%gru(:)%var(:)%dat (dp)
                    ! gru+hru dimension
                    gru_hru_int,           & ! x%gru(:)%hru(:)%var(:)     (i4b)
                    gru_hru_int8,          & ! x%gru(:)%hru(:)%var(:)     (i8b)
                    gru_hru_double,        & ! x%gru(:)%hru(:)%var(:)     (dp)
                    gru_hru_intVec,        & ! x%gru(:)%hru(:)%var(:)%dat (i4b)
                    gru_hru_doubleVec,     & ! x%gru(:)%hru(:)%var(:)%dat (dp)
                    ! gru+hru+dom dimension
                    gru_hru_dom_int,       & ! x%gru(:)%hru(:)%dom(:)%var(:)     (i4b)
                    gru_hru_dom_int8,      & ! x%gru(:)%hru(:)%dom(:)%var(:)     (i8b)
                    gru_hru_dom_double,    & ! x%gru(:)%hru(:)%dom(:)%var(:)     (dp)
                    gru_hru_dom_intVec,    & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat (i4b)
                    gru_hru_dom_doubleVec, & ! x%gru(:)%hru(:)%dom(:)%var(:)%dat (dp)
                    ! gru+hru+dom+z dimension
                    gru_hru_dom_z_vLookup, & ! x%gru(:)%hru(:)%dom(:)%z(:)%var(:)%lookup(:)  (dp)
                    ! gru+grid dimension
                    gru_grid_double,       & ! x%gru(:)%grid(:)%var(:)%dat2(:,:) (dp)
                    ! mapping between the GRUs and HRUs
                    gru2hru_map,           & ! x(iGRU)%hruinfo(iHRU)%y
                    hru2gru_map              ! x(iHRU)%y
USE data_types,      only: q_coupling      ! x(:)%id, x(:)%qsim

! access missing values
USE globalData,only:integerMissing      ! missing integer

! objective function
USE data_types,      only: obs_fileinfo    ! information on the observation file
USE data_types,      only: calib_info      ! calibration configuration

! mizuRoute coupling
#ifdef MIZUROUTE_ACTIVE
USE mizuroute_types, only: mizuroute_info
USE mizuroute_types, only: mizuroute_domain
#endif
implicit none

private

! ***********************************************************************************************************
! Configuration information shared across SUMMA simulations.
!
! Contains settings that are established during initial model configuration
! and can be reused when initializing individual SUMMA model instances.
! ***********************************************************************************************************
type, public :: config_info
  ! logging
  integer(i4b)                   :: iulog_summa = output_unit ! output unit for log files
  ! paths for a specific experiment
  character(len=1024)            :: cwd                     ! Current working directory
  character(len=:),  allocatable :: home_path_override      ! Command-line override for home_path
  character(len=:),  allocatable :: persistent_output       ! Directory for persistent output
  ! configuration flags
  logical(lgt)                   :: read_cli = .true.       ! .true. = read command-line interface
  logical(lgt)                   :: read_config = .true.    ! .true. = read configuration files
  ! Multi-case configuration
  logical(lgt)                   :: is_case_root = .true.   ! .true. for rank 0 within a model case
  integer(i4b)                   :: cases_per_node = 1      ! Number of concurrent cases per node
  character(len=:),  allocatable :: manifest_file           ! Path and name of the multi-case manifest
  character(len=64), allocatable :: case_names(:)           ! Names of cases defined in the manifest
  character(len=:),  allocatable :: manifest_casename       ! Case name selected from the run manifest
  character(len=:),  allocatable :: template_path           ! Path to the SUMMA configuration template
  character(len=:),  allocatable :: template_file           ! SUMMA configuration template filename
  ! SUMMA configuration options from the CLI (-g and -h)
  integer(i4b)                   :: nGRU_user = -1          ! Number of GRUs requested by the user
  integer(i4b)                   :: nHRU_check = 1          ! HRU used for diagnostic checks
  ! Parameter overrides
  character(len=64), allocatable :: param_name(:)           ! Names of parameters to override
  real(rkind),       allocatable :: param_value(:)          ! Values of parameter overrides
  ! Simulation
  character(len=:), allocatable  :: home_path               ! Root path for user-specific files
  character(len=:), allocatable  :: basin_dir               ! Directory containing basin-specific input data
  character(len=:), allocatable  :: case_name               ! Name of the simulation case
  character(len=:), allocatable  :: work_path               ! Path for simulation output
  character(len=:), allocatable  :: start_time              ! Start time of the simulation
  character(len=:), allocatable  :: end_time                ! End time of the simulation
  character(len=:), allocatable  :: time_zone               ! Time zone for simulation times
  ! SUMMA files and paths
  character(len=:), allocatable  :: settings_path           ! Path containing SUMMA settings files
  character(len=:), allocatable  :: forcing_path            ! Path containing forcing files
  character(len=:), allocatable  :: output_path             ! Path for SUMMA output files
  character(len=:), allocatable  :: state_path              ! Path containing model state files
  character(len=:), allocatable  :: init_condition          ! Initial-condition file
  character(len=:), allocatable  :: attributes              ! Local attributes file
  character(len=:), allocatable  :: trial_params            ! Trial parameter file
  character(len=:), allocatable  :: forcing_list            ! Forcing file list
  character(len=:), allocatable  :: decisions               ! Model decisions file
  character(len=:), allocatable  :: output_control          ! Output control file
  character(len=:), allocatable  :: local_parameters        ! Local (HRU) parameter information file
  character(len=:), allocatable  :: basin_parameters        ! Basin (GRU) parameter information file
  character(len=:), allocatable  :: vegetation_table        ! Vegetation parameter table
  character(len=:), allocatable  :: soil_table              ! Soil parameter table
  character(len=:), allocatable  :: general_table           ! General parameter table
  character(len=:), allocatable  :: noahmp_table            ! Noah-MP parameter table
  ! Observations and objective function
  type(obs_fileinfo)             :: obs                     ! Observation file configuration
  type(calib_info)               :: calib                   ! Calibration configuration
  ! Configuration sources
  character(len=:), allocatable  :: control_file            ! Legacy SUMMA control file
  character(len=:), allocatable  :: config_file             ! SUMMA TOML configuration file
  ! User configuration options
  logical(lgt)                   :: use_mizuroute = .false. ! Enable coupled mizuRoute for this simulation
  logical(lgt)                   :: write_timeseries = .true.  ! Write SUMMA time-series output file
#ifdef MIZUROUTE_ACTIVE
  type(mizuroute_info)           :: mizu_info               ! mizuRoute configuration infirmation
#endif
end type config_info

! ************************************************************************
! * parallel communication context
! ************************************************************************
type, public :: parallel_context_type
  integer(I4B) :: comm = -1
  integer(I4B) :: rank = 0
  integer(I4B) :: size = 1
end type parallel_context_type

! ************************************************************************
! * top-level summa data type
! *****************************************************************************
type, public :: summa1_type_dec    
    ! summa/mizuroute configuration settings
    type(config_info)                :: config                     ! CLI, file paths, observations, calibration
    ! MPI communication contexts (x%comm, x%rank, x%size)
    type(parallel_context_type)      :: domain_parallel            ! parallelization within one model instance
    type(parallel_context_type)      :: instance_parallel          ! parallelization across model instances
    ! define the lookup tables
    type(gru_hru_dom_z_vLookup)      :: lookupStruct               ! x%gru(:)%hru(:)%dom(:)%z(:)%var(:)%lookup(:) -- lookup tables
    ! define the statistics structures
    type(gru_hru_doubleVec)          :: forcStat                   ! x%gru(:)%hru(:)%var(:)%dat        -- model forcing data, does not need %dat but use so can use same structure as other data
    type(gru_hru_dom_doubleVec)      :: progStat                   ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model prognostic (state) variables
    type(gru_hru_dom_doubleVec)      :: diagStat                   ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model diagnostic variables
    type(gru_hru_dom_doubleVec)      :: fluxStat                   ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model fluxes
    type(gru_hru_dom_doubleVec)      :: indxStat                   ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model indices
    type(gru_doubleVec)              :: bvarStat                   ! x%gru(:)%var(:)%dat               -- basin-average variables
    ! define the primary data structures (scalars)
    type(var_i)                      :: timeStruct                 ! x%var(:)               -- model time data
    type(gru_hru_double)             :: forcStruct                 ! x%gru(:)%hru(:)%var(:) -- model forcing data
    type(gru_hru_double)             :: attrStruct                 ! x%gru(:)%hru(:)%var(:) -- local attributes for each HRU
    type(gru_hru_int)                :: typeStruct                 ! x%gru(:)%hru(:)%var(:) -- local classification of soil veg etc. for each HRU
    type(gru_hru_int8)               :: idStruct                   ! x%gru(:)%hru(:)%var(:) -- local values of hru and gru IDs
    ! define the primary data structures (variable length vectors)
    type(gru_hru_dom_intVec)         :: indxStruct                 ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model indices
    type(gru_hru_dom_doubleVec)      :: mparStruct                 ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model parameters
    type(gru_hru_dom_doubleVec)      :: progStruct                 ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model prognostic (state) variables
    type(gru_hru_dom_doubleVec)      :: diagStruct                 ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model diagnostic variables
    type(gru_hru_dom_doubleVec)      :: fluxStruct                 ! x%gru(:)%hru(:)%dom(:)%var(:)%dat -- model fluxes
    ! define the basin-average structures
    type(gru_double)                 :: bparStruct                 ! x%gru(:)%var(:)                   -- basin-average parameters
    type(gru_doubleVec)              :: bvarStruct                 ! x%gru(:)%var(:)%dat               -- basin-average variables
    type(gru_grid_double)            :: gridStruct                 ! x%gru(:)%grid(:)%var(:)%dat2(:,:) -- basin grid parameters and variables
    ! define the ancillary data structures
    type(gru_hru_double)             :: dparStruct                 ! x%gru(:)%hru(:)%var(:) -- default model parameters
    ! define the run-time variables
    type(gru_hru_i)                  :: computeVegFlux             ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
    type(gru_hru_dom_d)              :: dt_init                    ! used to initialize the length of the sub-step for each HRU
    type(gru_hru_d)                  :: upArea                     ! area upslope of each HRU
    ! GRU and HRU dimensions (this rank)
    integer(i4b)                     :: nGRU_local = 0             ! number of GRUs assigned to this rank
    integer(i4b)                     :: nHRU_local = 0             ! number of HRUs assigned to this rank
    integer(i4b)                     :: nDOM                       ! number of domains from the initial conditions file (same on all ranks)
    ! gru2hru mapping structures
    type(gru2hru_map), allocatable   :: gru_struc(:)               ! gru2hru map
    type(hru2gru_map), allocatable   :: index_map(:)               ! hru2gru map
    ! global time step information
    real(rkind)                      :: data_step                  ! length of the data window (seconds)
    integer(i4b)                     :: n_write                    ! length of the output buffer
    ! generic runoff coupling data
    type(q_coupling), allocatable    :: coupling(:)                ! x(:)%id, x(:)%qsim
#ifdef MIZUROUTE_ACTIVE
    type(mizuroute_domain)           :: mizu_domain                ! mizuroute domain data
#endif
    ! file managers
    character(len=256)               :: summaFileManagerFile       ! path/name of file defining directories and files
    character(len=256)               :: summaConfigFile = ''       ! path/name of the TOML configuration file
end type summa1_type_dec

END MODULE summa_type

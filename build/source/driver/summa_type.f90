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
                    var_i,               & ! x%var(:)            (i4b)
                    var_d,               & ! x%var(:)            (dp)
                    var_ilength,         & ! x%var(:)%dat        (i4b)
                    var_dlength,         & ! x%var(:)%dat        (dp)
                    ! no variable dimension
                    hru_i,               & ! x%hru(:)            (i4b)
                    hru_d,               & ! x%hru(:)            (dp)
                    gru_i,               & ! x%gru(:)%hru(:)     (i4b)
                    gru_d,               & ! x%gru(:)%hru(:)     (dp)
                    ! gru dimension
                    gru_int,             & ! x%gru(:)%var(:)     (i4b)
                    gru_double,          & ! x%gru(:)%var(:)     (dp)
                    gru_intVec,          & ! x%gru(:)%var(:)%dat (i4b)
                    gru_doubleVec,       & ! x%gru(:)%var(:)%dat (dp)
                    ! gru+hru dimension
                    gru_hru_int,         & ! x%gru(:)%hru(:)%var(:)     (i4b)
                    gru_hru_int8,        & ! x%gru(:)%hru(:)%var(:)     (i8b)
                    gru_hru_double,      & ! x%gru(:)%hru(:)%var(:)     (dp)
                    gru_hru_intVec,      & ! x%gru(:)%hru(:)%var(:)%dat (i4b)
                    gru_hru_doubleVec,   & ! x%gru(:)%hru(:)%var(:)%dat (dp)
                    ! gru+hru+z dimension
                    gru_hru_z_vLookup,   & ! x%gru(:)%hru(:)%z(:)%var(:)%lookup(:)  (dp)
                    ! mapping between the GRUs and HRUs
                    gru2hru_map,         & ! x(iGRU)%hruinfo(iHRU)%y 
                    hru2gru_map            ! x(iHRU)%y

USE data_types,      only: q_coupling      ! x(:)%id, x(:)%qsim

USE data_types,      only: obs_fileinfo    ! information on the observation file
USE data_types,      only: calib_info      ! calibration configuration 

! mizuRoute coupling
#ifdef MIZUROUTE_ACTIVE
use mizuroute_types, only: mizuroute_info
use mizuroute_types, only: mizuroute_domain
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

  ! SUMMA configuration options from the CLI (-g and -h)
  integer(i4b)                   :: nGRU_user = -1          ! Number of GRUs requested by the user
  integer(i4b)                   :: nHRU_check = 1          ! HRU used for diagnostic checks

  ! Parameter overrides
  character(len=64), allocatable :: param_name(:)           ! Names of parameters to override
  real(rkind),       allocatable :: param_value(:)          ! Values of parameter overrides

  ! Simulation
  character(len=:), allocatable  :: case_name               ! Name of the simulation case
  character(len=:), allocatable  :: start_time              ! Simulation start time
  character(len=:), allocatable  :: end_time                ! Simulation end time
  character(len=:), allocatable  :: time_zone               ! Time zone used for simulation times

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

! summa/mizuroute information
type(config_info)                :: config                     ! summa/mizuroute configuration settings

! MPI communication context (x%comm, x%rank, x%size)
type(parallel_context_type)      :: domain_parallel            ! parallelization within one model instance
type(parallel_context_type)      :: instance_parallel          ! parallelization across model instances

! the lookup tables
type(gru_hru_z_vLookup)          :: lookupStruct               ! x%gru(:)%hru(:)%z(:)%var(:)%lookup(:) -- lookup tables

! the statistics structures
type(gru_hru_doubleVec)          :: forcStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model forcing data
type(gru_hru_doubleVec)          :: progStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
type(gru_hru_doubleVec)          :: diagStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
type(gru_hru_doubleVec)          :: fluxStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes
type(gru_hru_doubleVec)          :: indxStat                   ! x%gru(:)%hru(:)%var(:)%dat -- model indices
type(gru_doubleVec)              :: bvarStat                   ! x%gru(:)%var(:)%dat        -- basin-average variable

! the primary data structures (scalars)
type(var_i)                      :: timeStruct                 ! x%var(:)                   -- model time data
type(gru_hru_double)             :: forcStruct                 ! x%gru(:)%hru(:)%var(:)     -- model forcing data
type(gru_hru_double)             :: attrStruct                 ! x%gru(:)%hru(:)%var(:)     -- local attributes for each HRU
type(gru_hru_int)                :: typeStruct                 ! x%gru(:)%hru(:)%var(:)     -- local classification of soil veg etc. for each HRU
type(gru_hru_int8)               :: idStruct                   ! x%gru(:)%hru(:)%var(:)     -- local values of hru and gru IDs

! the primary data structures (variable length vectors)
type(gru_hru_intVec)             :: indxStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model indices
type(gru_hru_doubleVec)          :: mparStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model parameters
type(gru_hru_doubleVec)          :: progStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model prognostic (state) variables
type(gru_hru_doubleVec)          :: diagStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model diagnostic variables
type(gru_hru_doubleVec)          :: fluxStruct                 ! x%gru(:)%hru(:)%var(:)%dat -- model fluxes

! the basin-average structures
type(gru_double)                 :: bparStruct                 ! x%gru(:)%var(:)            -- basin-average parameters
type(gru_doubleVec)              :: bvarStruct                 ! x%gru(:)%var(:)%dat        -- basin-average variables

! the ancillary data structures
type(gru_hru_double)             :: dparStruct                 ! x%gru(:)%hru(:)%var(:)     -- default model parameters

! the run-time variables
type(gru_i)                      :: computeVegFlux             ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
type(gru_d)                      :: dt_init                    ! used to initialize the length of the sub-step for each HRU
type(gru_d)                      :: upArea                     ! area upslope of each HRU

! GRU and HRU dimensions
integer(i4b)                     :: nGRU_local = 0             ! number of GRUs assigned to this rank
integer(i4b)                     :: nHRU_local = 0             ! number of HRUs assigned to this rank

! gru2hru mapping structures
type(gru2hru_map), allocatable   :: gru_struc(:)               ! gru2hru map
type(hru2gru_map), allocatable   :: index_map(:)               ! hru2gru map

! global time step information
real(dp)                         :: data_step                  ! length of the data window (seconds)
integer(i4b)                     :: n_write                    ! length of the output buffer

! generic runoff coupling data
type(q_coupling), allocatable    :: coupling(:)                ! x(:)%id, x(:)%qsim

#ifdef MIZUROUTE_ACTIVE
type(mizuroute_domain)           :: mizu_domain                ! mizuroute domain data
#endif

end type summa1_type_dec

END MODULE summa_type

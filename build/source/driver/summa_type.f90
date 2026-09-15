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

! used to define master summa data structure

! *****************************************************************************
! * higher-level derived data types
! *****************************************************************************

USE nr_type         ! variable types, etc.

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

! generic runoff coupling structure
USE data_types,      only: q_coupling      ! x(:)%id, x(:)%qsim

! access missing values
USE globalData,only:integerMissing      ! missing integer

! mizuRoute coupling
#ifdef MIZUROUTE_ACTIVE
USE mizuroute_types, only: mizuroute_info
USE mizuroute_types, only: mizuroute_domain
#endif

implicit none

private

! ************************************************************************
! * parallel communication context
! *****************************************************************************

type, public :: parallel_context_type
  integer(I4B) :: comm = -1
  integer(I4B) :: rank = 0
  integer(I4B) :: size = 1
end type parallel_context_type

! ************************************************************************
! * master summa data type
! *****************************************************************************
type, public :: summa1_type_dec    
    ! MPI communication context
    type(parallel_context_type)      :: parallel                   ! x%comm, x%rank, x%size
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
    ! GRU and HRU dimensions
    integer(i4b)                     :: nGRU_user = integerMissing ! number of GRUs requested with CLI -g
    integer(i4b)                     :: nHRU_check = 1             ! number of HRUs requested with CLI -h
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
    type(mizuroute_info)             :: mizu_info                  ! mizuroute information structure
    type(mizuroute_domain)           :: mizu_domain                ! mizuroute domain data
#endif
    ! file managers
    character(len=256)               :: summaFileManagerFile       ! path/name of file defining directories and files
    character(len=256)               :: summaConfigFile = ''       ! path/name of the TOML configuration file
end type summa1_type_dec

END MODULE summa_type

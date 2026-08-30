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
                    gru_hru_z_vLookup      ! x%gru(:)%hru(:)%z(:)%var(:)%lookup(:)  (dp)

! mizuRoute coupling
#ifdef MIZUROUTE_ACTIVE
use data_types,      only: q_coupling
use mizuroute_types, only: mizuroute_info
use mizuroute_types, only: mizuroute_domain
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

! runoff exchanged between SUMMA and mizuRoute
type(q_coupling), allocatable    :: coupling(:)                ! x(:)%id, x(:)%qsim

! GRU and HRU dimensions
integer(i4b)                     :: nGRU_user = -1             ! number of GRUs requested with CLI -g
integer(i4b)                     :: nHRU_check = 1             ! number of HRUs requested with CLI -h
integer(i4b)                     :: nGRU_local = 0             ! number of GRUs assigned to this rank
integer(i4b)                     :: nHRU_local = 0             ! number of HRUs assigned to this rank

! global time step information
real(dp)                         :: data_step                  ! length of the data window (seconds)
integer(i4b)                     :: n_write                    ! length of the output buffer

! file managers
character(len=256)               :: summaFileManagerFile       ! path/name of file defining directories and files
character(len=256)               :: summaConfigFile =''        ! path/name of the TOML configuration file 

#ifdef MIZUROUTE_ACTIVE
type(mizuroute_info)   :: mizu_info
type(mizuroute_domain) :: mizu_domain
#endif

end type summa1_type_dec

END MODULE summa_type

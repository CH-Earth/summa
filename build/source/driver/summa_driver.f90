! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
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

program summa_driver
! driver program for summa simulations
! *****************************************************************************
! * use desired modules
! *****************************************************************************
! provide access to data types
USE nrtype                                                  ! variable types, etc.
USE summa_type, only:summa1_type_dec                        ! master summa data type
! provide access to subroutines and functions
USE summa_init, only:summa_initialize                       ! used to allocate/initialize summa data structures
USE summa_util, only:stop_program                           ! used to stop the summa program (with errors)
USE summa_util, only:handle_err                             ! used to process errors
implicit none

! *****************************************************************************
! * variable definitions
! *****************************************************************************
! define the master summa data structure
type(summa1_type_dec), allocatable :: summa1_struc(:)
! define parameters for the model simulation
integer(i4b), parameter            :: nInstantiation=1           ! number of instantiations
integer(i4b)                       :: n                          ! index of model instantiation
! error control
integer(i4b)                       :: err=0                      ! error code
character(len=1024)                :: message=''                 ! error message

! *****************************************************************************
! * main driver
! *****************************************************************************

! allocate space for the master summa structure
allocate(summa1_struc(nInstantiation), stat=err)
if(err/=0) call stop_program(1, 'problem allocating master summa structure')

! declare and allocate summa data structures and initialize model state to known values
do n = 1, nInstantiation
 call summa_initialize(summa1_struc(n), err, message)
 call handle_err(err, message)
end do ! loop through model instantiations

! successful end
call stop_program(0, 'finished simulation successfully.')
end program summa_driver

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

program summa_driver

  USE nr_type, only: i4b, rkind
  USE summa_simulation, only: evaluate_objective
  USE summa_util, only: handle_err, stop_program

  implicit none

  ! parallel dummy variables
  integer(i4b), parameter :: comm=0, rank=0, nproc=1

  ! parameter overrides
  character(len=64), allocatable :: param_name(:)
  real(rkind),       allocatable :: param_value(:)

  ! objective function
  real(rkind)                    :: objective

  ! error control
  integer(i4b)        :: err=0
  character(len=1024) :: message=''

  ! no externally supplied parameter overrides
  allocate(param_name(0))
  allocate(param_value(0))

  ! run SUMMA and evaluate the objective function
  call evaluate_objective(comm, rank, nproc,        &
                          param_name, param_value,  &
                          objective,                &
                          err, message)
  call handle_err(err,message)

  call stop_program(0,'finished simulation successfully.')

end program summa_driver

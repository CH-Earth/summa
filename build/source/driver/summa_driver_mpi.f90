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

program summa_driver_mpi
  
  ! **** Driver program for SUMMA simulations ****

  ! * module access *
  ! data types
  USE nr_type, only: i4b, rkind                               ! variable types, etc.
  USE mpi, only : MPI_COMM_WORLD, MPI_SUCCESS                 ! MPI constants
  
  ! subroutines and functions: MPI 
  USE mpi, only : MPI_Init, MPI_Finalize                      ! MPI subroutine interfaces
  USE mpi_context, only : set_mpi_context                     ! define MPI communication context
  USE error_utils, only : check_mpi, abort_mpi                ! check MPI errors
 
  ! subroutines and functions: SUMMA
  USE summa_simulation, only: run_simulation                  ! run a model simulation
  USE summa_util,       only: handle_err, stop_program        ! error handling

  implicit none

  ! * driver variables *

  ! MPI
  integer(i4b) :: rank = 0
  integer(i4b) :: size = 1

  ! parameters
  character(len=64), allocatable :: param_name(:)
  real(rkind),       allocatable :: param_value(:)

  ! flow
  real(rkind),       allocatable :: timeSim(:)
  real(rkind),       allocatable :: flowSim(:)

  ! units
  character(len=:),  allocatable :: timeUnits
  character(len=:),  allocatable :: flowUnits

  ! error codes
  integer(i4b) :: err = 0
  integer(i4b) :: mpi_err = 0

  ! error messages
  character(len=1024) :: message = ''
  character(len=256)  :: mpi_message = ''

  ! -------------------------------------------------------------------

  call MPI_Init(mpi_err)
  call check_mpi(-1, mpi_err, 'MPI_Init failed')

  call set_mpi_context(MPI_COMM_WORLD, rank, size, mpi_err, mpi_message)
  if (mpi_err /= MPI_SUCCESS) call abort_mpi(rank, trim(mpi_message)) 

  call run_simulation(MPI_COMM_WORLD, rank, size,              &
                      timeSim, flowSim, timeUnits, flowUnits,  &
                      param_name, param_value,  &
                      err, message)
  call handle_err(err, message)

  call MPI_Finalize(mpi_err)
  if (mpi_err /= MPI_SUCCESS)then
    write(message,'(A,I0,A)') 'ERROR [rank ', rank, ']: MPI_Finalize failed'
    call handle_err(mpi_err, message)
  endif

  if (rank == 0) then
    call stop_program(0, 'finished simulation successfully.')
  end if

end program summa_driver_mpi

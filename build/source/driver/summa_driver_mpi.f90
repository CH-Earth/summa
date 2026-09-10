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
  USE summa_type, only: config_info                           ! summa configuration
  USE summa_type, only: parallel_context_type                 ! parallel context
  USE mpi, only : MPI_COMM_WORLD, MPI_COMM_SELF, MPI_SUCCESS  ! MPI constants
  
  ! subroutines and functions: MPI 
  USE mpi, only : MPI_Init, MPI_Finalize                      ! MPI subroutine interfaces
  USE mpi_context, only : set_mpi_context                     ! define MPI communication context
  USE error_utils, only : check_mpi, abort_mpi                ! check MPI errors
 
  ! subroutines and functions: SUMMA
  USE summa_simulation, only: run_simulation                  ! run a model simulation
  USE summa_util,       only: handle_err, stop_program        ! error handling

  implicit none

  ! * driver variables *

  ! MPI contexts for domain and model-instance parallelism
  type(parallel_context_type)    :: domain_parallel
  type(parallel_context_type)    :: instance_parallel

  ! configuration info
  type(config_info)              :: config

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

  ! ---------------------------------------------------------------------------------------
  ! Initialize MPI
  ! ---------------------------------------------------------------------------------------
  
  call MPI_Init(mpi_err)
  call check_mpi(-1,mpi_err,'MPI_Init failed')
  
  ! distribute the model domain across MPI processes
  domain_parallel%comm=MPI_COMM_WORLD
  
  call set_mpi_context(domain_parallel%comm,  &
                       domain_parallel%rank,  &
                       domain_parallel%size,  &
                       mpi_err,mpi_message)
  
  if(mpi_err/=MPI_SUCCESS)then
    call abort_mpi(domain_parallel%rank,trim(mpi_message))
  endif
  
  ! run a single model instance
  instance_parallel%comm=MPI_COMM_SELF
  instance_parallel%rank=0
  instance_parallel%size=1

  ! ---------------------------------------------------------------------------------------
  ! Run a SUMMA simulation
  ! ---------------------------------------------------------------------------------------

  ! no externally supplied parameter overrides
  allocate(param_name(0))
  allocate(param_value(0))

  call run_simulation(config,                 & ! SUMMA configuration structure
                      domain_parallel,        & ! MPI context for domain parallelism
                      instance_parallel,      & ! MPI context for model-instance parallelism
                      timeSim,flowSim,        & ! simulated time and streamflow
                      timeUnits,flowUnits,    & ! time and streamflow units
                      param_name,param_value, & ! parameter names and values
                      err, message)             ! error code and message
  call handle_err(err,message)

  ! ---------------------------------------------------------------------------------------
  ! Finalize MPI
  ! ---------------------------------------------------------------------------------------

  call MPI_Finalize(mpi_err)

  if(mpi_err/=MPI_SUCCESS)then
    write(message,'(A,I0,A)') 'ERROR [rank ',domain_parallel%rank,']: MPI_Finalize failed'
    call handle_err(mpi_err,message)
  endif
  
  if(domain_parallel%rank==0)then
    call stop_program(0,'finished simulation successfully.')
  endif

end program summa_driver_mpi

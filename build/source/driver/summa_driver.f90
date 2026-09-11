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
  USE summa_type, only: config_info
  USE summa_type, only: parallel_context_type

  USE summa_simulation, only: evaluate_objective
  USE summa_util, only: handle_err, stop_program

  implicit none

  ! configuration info
  type(config_info)              :: config

  ! MPI contexts (initialize both as serial)
  type(parallel_context_type)    :: domain_parallel
  type(parallel_context_type)    :: instance_parallel

  ! parameter overrides
  character(len=64), allocatable :: param_name(:)
  real(rkind),       allocatable :: param_value(:)

  ! objective function
  integer(i4b), parameter        :: sample_id = 0
  real(rkind)                    :: objective

  ! error control
  integer(i4b)        :: err=0
  character(len=1024) :: message=''

  ! serial domain execution
  domain_parallel%comm=0
  domain_parallel%rank=0
  domain_parallel%size=1
  
  ! single model instance
  instance_parallel%comm=0
  instance_parallel%rank=0
  instance_parallel%size=1

  ! no externally supplied parameter overrides
  allocate(param_name(0))
  allocate(param_value(0))

  ! run SUMMA and evaluate the objective function
  call evaluate_objective(config,                            & ! SUMMA configuration structure
                          domain_parallel,                   & ! MPI context for domain parallelism
                          instance_parallel,                 & ! MPI context for model-instance parallelism
                          sample_id,param_name,param_value,  & ! sample ID + parameter names and values
                          objective,                         & ! objective function value
                          err, message)                        ! error code and message
  call handle_err(err,message)

  call stop_program(0,'finished simulation successfully.')

end program summa_driver

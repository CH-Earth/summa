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
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program. If not, see <http://www.gnu.org/licenses/>.

! **************************************************************************************************
! SUMMA optimization driver
!
! Performs a common cold-start spinup and dynamically evaluates parameter
! samples across independent MPI model instances.
! **************************************************************************************************
program summa_driver_opt

  ! data types
  USE nr_type,    only: i4b
  USE summa_type, only: config_info,parallel_context_type
  USE summa_util, only: stop_program

  ! logging
  USE iso_fortran_env, only: output_unit

  ! SUMMA-specific calibration
  USE summa_calibration_case, only: initialize_summa_calibration_config
  USE summa_calibration_case, only: summa_case_context,execute_summa_case
  USE summa_calibration_case, only: run_case

  ! calibration MPI infrastructure
  USE calibration_parallel, only: init_calibration_parallel

  ! dynamic scheduling of independent calibration cases
  USE dynamic_case_dispatch, only: dispatch_calibration_cases

  ! MPI
  USE mpi, only: MPI_Init,MPI_Finalize,MPI_Barrier
  USE mpi, only: MPI_Bcast
  USE mpi, only: MPI_Comm_split_type,MPI_Comm_split
  USE mpi, only: MPI_COMM_TYPE_SHARED,MPI_INFO_NULL,MPI_UNDEFINED
  USE mpi, only: MPI_INTEGER,MPI_COMM_WORLD,MPI_COMM_SELF,MPI_SUCCESS
  USE mpi_context, only: set_mpi_context
  USE error_utils, only: check_mpi,abort_mpi

  implicit none

  type(config_info), target   :: config                    ! SUMMA configuration information
  type(parallel_context_type) :: world_parallel            ! all MPI ranks
  type(parallel_context_type) :: node_parallel             ! all MPI ranks on the same physical node
  type(parallel_context_type) :: leader_parallel           ! one leader rank from each physical node
  type(parallel_context_type) :: instance_parallel         ! ranks assigned to one calibration case
    
  type(summa_case_context)    :: dispatch_context          ! module-level callback for a calibration case

  integer(i4b) :: node_index                               ! Index of the physical node
  integer(i4b) :: nNodes                                   ! Number of physical nodes
  integer(i4b) :: case_group                               ! Case group on the local node
  integer(i4b) :: global_case_group                        ! Case group across all nodes
  integer(i4b) :: nCaseGroups                              ! Total number of case groups
  integer(i4b) :: ranks_per_case                           ! MPI ranks assigned to each case
  integer(i4b) :: leader_color                             ! Color used to construct node-leader communicator
  integer(i4b) :: requested_cases                          ! cases_per_node as configured, before any reduction
  integer(i4b) :: nCases                                   ! Total number of cases available for execution
  integer(i4b)        :: err=0                             ! SUMMA error code
  integer(i4b)        :: mpi_err=0                         ! MPI error code
  character(len=1024) :: message=''                        ! SUMMA error message
  character(len=256)  :: mpi_message=''                    ! MPI error message
  ! ---------------------------------------------------------------------------------------
  ! Initialize MPI
  ! ---------------------------------------------------------------------------------------

  call MPI_Init(mpi_err)
  call check_mpi(-1,mpi_err,'MPI_Init failed')

  ! establish the global MPI context
  world_parallel%comm=MPI_COMM_WORLD

  call set_mpi_context(world_parallel%comm,  & ! global MPI communicator
                       world_parallel%rank,  & ! rank in global communicator
                       world_parallel%size,  & ! number of global MPI ranks
                       mpi_err,mpi_message)    ! MPI error information

  if(mpi_err/=MPI_SUCCESS) call abort_mpi(world_parallel%rank,trim(mpi_message))

  ! ---------------------------------------------------------------------------------------
  ! Read run configuration
  ! ---------------------------------------------------------------------------------------

  call initialize_summa_calibration_config(config,err,message)
  if(err/=0) call abort_mpi(world_parallel%rank,trim(message))

  ! ---------------------------------------------------------------------------------------
  ! Build the MPI communicator hierarchy
  ! ---------------------------------------------------------------------------------------

  call init_calibration_parallel(allocated(config%manifest_file),  & ! Single- or multi-case execution
                                 config%cases_per_node,            & ! Concurrent cases per node
                                 world_parallel,                   & ! Global MPI context
                                 node_parallel,instance_parallel,  & ! Node and instance MPI contexts
                                 node_index,case_group,            & ! Local node and case indices
                                 global_case_group,nCaseGroups,    & ! Global case-group information
                                 err,message)                        ! Error code and message

  if(err/=0) call abort_mpi(world_parallel%rank,trim(message))

  ! Log Processor layout: report the MPI rank and calibration-group assignment
  if(allocated(config%manifest_file))then
    write(*,'(A,I4,A,I3,A,I3,A,I3,A,I3,A,I3,A,I3)') &
      '  rank = '          , world_parallel%rank,    &
      ', node = '          , node_index,             &
      ', node_rank = '     , node_parallel%rank,     &
      ', case_group = '    , case_group,             &
      ', instance_rank = ' , instance_parallel%rank, &
      ', global_group = '  , global_case_group,      &
      ', nCaseGroups = '   , nCaseGroups
  endif

  ! ---------------------------------------------------------------------------------------
  ! Run SUMMA cases
  ! ---------------------------------------------------------------------------------------

  if(allocated(config%manifest_file))then

    ! number of independent calibration cases
    nCases=size(config%case_names)*config%n_runs

    ! supply the module-level callback with this rank's case configuration.
    dispatch_context%config            => config
    dispatch_context%node_index        =  node_index
    dispatch_context%instance_parallel =  instance_parallel
    dispatch_context%global_case_group =  global_case_group

    ! dynamically schedule cases across independent MPI instance groups.
    call dispatch_calibration_cases(nCases,world_parallel,instance_parallel, &
                                    execute_summa_case,dispatch_context,err,message)
    if(err/=0) call abort_mpi(world_parallel%rank,trim(message))

  else

    ! execute a single calibration case
    call run_case(config,instance_parallel,node_index,global_case_group, &
                  err=err,message=message)
    if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))

  endif

  ! ---------------------------------------------------------------------------------------
  ! Finalize MPI
  ! ---------------------------------------------------------------------------------------

  call MPI_Barrier(world_parallel%comm,mpi_err)
  call check_mpi(world_parallel%rank,mpi_err,'MPI_Barrier failed before MPI_Finalize')

  call MPI_Finalize(mpi_err)
  call check_mpi(world_parallel%rank,mpi_err,'MPI_Finalize failed')

  if(world_parallel%rank==0) call stop_program(0,'finished parallel parameter evaluation successfully.')
  
end program summa_driver_opt

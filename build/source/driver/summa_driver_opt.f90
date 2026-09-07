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

program summa_driver_opt

  ! **** Driver program for parallel SUMMA parameter evaluation ****

  ! data types
  USE nr_type,    only: i4b, rkind
  USE summa_type, only: config_info

  ! SUMMA global data
  USE globalData, only: iulog 
  USE globalData, only: initConfig 
  USE globalData, only: output_fileSuffix       ! modify based on context/rank

  USE globalData, only: ixRestart               ! define frequency to write restart files
  USE globalData, only: ixRestart_end           ! named variable to write a re-start file at the end of a run
  USE globalData, only: ixRestart_never         ! named variable to never write a re-start file
  USE globalData, only: restart_filename        ! name of the restart file

  ! model control
  USE summaFileManager, only: SIM_START_TM      ! simulation start time
  USE summaFileManager, only: SIM_END_TM        ! simulation end time
  USE summaFileManager, only: MODEL_INITCOND    ! name of the model initial conditions file

  ! MPI constants
  USE mpi, only: MPI_COMM_WORLD
  USE mpi, only: MPI_COMM_SELF
  USE mpi, only: MPI_SUCCESS

  ! MPI routines
  USE mpi, only: MPI_Init
  USE mpi, only: MPI_Finalize

  USE mpi_context, only: set_mpi_context
  USE error_utils, only: check_mpi, abort_mpi

  ! SUMMA
  USE summa_init,             only: init_config
  USE summa_simulation,       only: run_simulation
  USE summa_simulation,       only: evaluate_objective
  USE summa_util,             only: handle_err, stop_program

  ! SUMMA parameter sampling
  USE summa_parameter_search, only: parameter_search_info
  USE summa_parameter_search, only: initialize_parameter_search
  USE summa_parameter_search, only: sample_parameters

  implicit none

  ! MPI context for the ensemble
  integer(i4b) :: rank  = 0
  integer(i4b) :: nproc = 1

  ! looping
  integer(i4b) :: i,j

  ! configuration info
  type(config_info)              :: config
  character(len=4)               :: rankString
  character(len=:),  allocatable :: outputFileSuffix_orig
  
  character(len=:),  allocatable :: simStartOriginal
  character(len=:),  allocatable :: simEndOriginal
  character(len=:),  allocatable :: spinStart
  integer(i4b)                   :: iyear

  ! parameters
  
  character(len=64), allocatable :: param_name(:)
  real(rkind),       allocatable :: param_value(:)

  type(parameter_search_info)    :: search

  integer(i4b)                   :: nseed
  integer(i4b), allocatable      :: seed(:)

  ! flow
  real(rkind),       allocatable :: timeSim(:)
  real(rkind),       allocatable :: flowSim(:)

  ! units
  character(len=:),  allocatable :: timeUnits
  character(len=:),  allocatable :: flowUnits

  ! objective function
  real(rkind) :: objective

  ! errors
  integer(i4b) :: err     = 0
  integer(i4b) :: mpi_err = 0

  character(len=1024) :: message     = ''
  character(len=256)  :: mpi_message = ''

  ! -------------------------------------------------------------------

  ! initialize MPI
  call MPI_Init(mpi_err)
  call check_mpi(-1,mpi_err,'MPI_Init failed')

  ! ensemble-level MPI context
  call set_mpi_context(MPI_COMM_WORLD,rank,nproc,mpi_err,mpi_message)
  if(mpi_err/=MPI_SUCCESS) call abort_mpi(rank,trim(mpi_message))

  ! get the rank string for use in output files
  write(rankString,'(I4.4)') rank

  allocate(param_name(0))
  allocate(param_value(0))
  
  ! ---------------------------------------------------------------------------------------
  ! Initialize SUMMA and perform a one-year spinup from a cold start.
  ! 
  ! This initial run populates the model data structures and parameter setup, and generates
  ! a restart state that provides the initial conditions for subsequent simulations.
  ! ---------------------------------------------------------------------------------------

  ! initialize SUMMA configuration once
  call init_config(config,err,message)
  call handle_err(err,message)

  initConfig = .false.

  ! save the original start and end times
  simStartOriginal = trim(SIM_START_TM)
  simEndOriginal   = trim(SIM_END_TM)

  ! specify an initial spinup from a cold start for one year before the normal simulation start
  spinStart = simStartOriginal
  read(spinStart(1:4),*) iyear
  write(spinStart(1:4),'(I4.4)') (iyear-1)

  SIM_START_TM = spinStart
  SIM_END_TM   = simStartOriginal

  ! define a new output file suffix
  outputFileSuffix_orig = trim(output_fileSuffix)
  output_fileSuffix     = trim(outputFileSuffix_orig)//'_spinup_rank'//rankString

  ! force writing the restart file at the end of the warm-up period
  ! (ixRestart shared in global data)
  ixRestart = ixRestart_end  

  ! Run baseline SUMMA independently on every MPI process
  ! (this initial run populates the model data structures and parameter setup)
  call run_simulation(config,                         &
                      MPI_COMM_SELF,0,1,              &
                      timeSim,flowSim,                &
                      timeUnits,flowUnits,            &
                      param_name,param_value,         &
                      err,message)
  call handle_err(err,message)

  ! set initial conditions filename to the state generated during spinup
  if(allocated(restart_filename))then
    MODEL_INITCOND = trim(restart_filename)
  else
    call handle_err(20, 'restart filename not defined')
  endif

  ! restore original settings
  SIM_START_TM      = simStartOriginal
  SIM_END_TM        = simEndOriginal
  output_fileSuffix = outputFileSuffix_orig 
  ixRestart         = ixRestart_never

 
  ! initialize parameter search
  call initialize_parameter_search(config,search,err,message)
  call handle_err(err,message)

  deallocate(param_name,param_value)

  allocate(param_name(size(search%param_names)))
  allocate(param_value(size(search%param_names)))

  param_name = search%param_names

  write(iulog,'(/,A)') 'Calibration parameter bounds:'
  write(iulog,'(A)')   '  Parameter                         Lower              Upper'

  do i=1,size(search%param_names)
     write(iulog,'(2X,A30,2X,ES16.8,2X,ES16.8)') &
      trim(search%param_names(i)), search%lower(i), search%upper(i)
  enddo

  ! initialize random-number generator
  call random_seed(size=nseed)
  allocate(seed(nseed))

  seed = rank + 42
  call random_seed(put=seed)


  do j=1,10

    ! Sample a feasible parameter vector
    call sample_parameters(search,param_value,err,message)
    call handle_err(err,message)

    ! Evaluate SUMMA for this parameter vector
    call evaluate_objective(config,                     &
                            MPI_COMM_SELF,0,1,           &
                            search%param_names,          &
                            param_value,                 &
                            objective,                   &
                            err,message)
    call handle_err(err,message)
   
  ! Print sample, parameters, and objective
  write(iulog,'(I6,*(1X,ES16.8))') j,param_value,objective

  enddo

  ! finalize MPI
  call MPI_Finalize(mpi_err)

  if(mpi_err/=MPI_SUCCESS)then
    write(message,'(A,I0,A)') 'ERROR [rank ',rank,']: MPI_Finalize failed'
    call handle_err(mpi_err,message)
  endif

  if(rank==0)then
    call stop_program(0,'finished parallel parameter evaluation successfully.')
  endif

end program summa_driver_opt

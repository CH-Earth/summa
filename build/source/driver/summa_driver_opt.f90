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


! *********************************************************************
! SUMMA optimization driver
!
! Initializes SUMMA and the parameter search infrastructure, generates
! trial parameter vectors, applies the corresponding SUMMA parameter
! overrides, and evaluates the calibration objective for each trial.
!
! The parameter search algorithms are implemented in parameter_search;
! SUMMA-specific parameter metadata and overrides are handled by
! summa_parameter_spec.
! *********************************************************************

program summa_driver_opt

  ! MPI
  USE mpi, only: MPI_COMM_WORLD,MPI_COMM_SELF,MPI_SUCCESS,MPI_Init,MPI_Finalize
  USE mpi_context, only: set_mpi_context
  USE error_utils, only: check_mpi,abort_mpi

  ! logging
  USE globalData,       only: iulog
  USE iso_fortran_env,  only: error_unit
  
  ! data types
  USE nr_type,    only: i4b, rkind
  USE summa_type, only: config_info
  USE summa_type, only: parallel_context_type

  ! error handling
  USE summa_util, only: handle_err,stop_program

  implicit none

  ! SUMMA configuration
  type(config_info) :: config

  ! MPI contexts for domain and model-instance parallelism
  type(parallel_context_type) :: domain_parallel
  type(parallel_context_type) :: instance_parallel

  ! number of parameter samples
  integer(i4b), parameter :: nSamples=1000   ! total number of trials across all model instances
  integer(i4b)            :: nLocal          ! number of trials assigned to this model instance

  ! error control
  integer(i4b) :: err=0,mpi_err=0

  character(len=1024) :: message=''
  character(len=256)  :: mpi_message=''


  ! ---------------------------------------------------------------------------------------
  ! Initialize MPI
  ! ---------------------------------------------------------------------------------------

  call MPI_Init(mpi_err)
  call check_mpi(-1,mpi_err,'MPI_Init failed')

  ! get the MPI context for the optimization ensemble
  instance_parallel%comm=MPI_COMM_WORLD

  call set_mpi_context(instance_parallel%comm,  &
                       instance_parallel%rank,  &
                       instance_parallel%size,  &
                       mpi_err,mpi_message)
  if(mpi_err/=MPI_SUCCESS) &
    call abort_mpi(instance_parallel%rank,trim(mpi_message))

  ! each model instance runs independently without domain parallelism
  domain_parallel%comm=MPI_COMM_SELF
  domain_parallel%rank=0
  domain_parallel%size=1

  ! divide the total number of parameter trials across model instances
  nLocal=nSamples/instance_parallel%size

  if(instance_parallel%rank < mod(nSamples,instance_parallel%size))then
    nLocal=nLocal+1
  endif

  ! ---------------------------------------------------------------------------------------
  ! Initialize SUMMA
  ! ---------------------------------------------------------------------------------------

  ! use stderr for logging during initialization
  iulog = error_unit

  call spinup_from_cold(config,             & ! SUMMA configuration structure
                        domain_parallel,    & ! MPI context for domain parallelism
                        instance_parallel,  & ! MPI context for model-instance parallelism
                        err, message)         ! error code and error message
  call handle_err(err,message)

  ! subsequent output is written to the rank-specific log

  ! ---------------------------------------------------------------------------------------
  ! Sample parameters and evaluate SUMMA
  ! ---------------------------------------------------------------------------------------

  ! identify case and MPI rank in the rank-specific log
  write(iulog,'(A)') repeat('-',80)
  write(iulog,'(A,A,A,I0)') 'INFO: running case ',trim(config%case_name), &
                            ' on rank ',instance_parallel%rank
  write(iulog,'(A)') repeat('-',80)

  call evaluate_parameter_samples(config,                 & ! SUMMA configuration structure
                                  domain_parallel,        & ! MPI context for domain parallelism
                                  instance_parallel,      & ! MPI context for model-instance parallelism
                                  nLocal,                 & ! number of parameter trials assigned to this instance
                                  err, message)             ! error code and message
  call handle_err(err,message)


  ! ---------------------------------------------------------------------------------------
  ! Finalize MPI
  ! ---------------------------------------------------------------------------------------

  ! close the logging file
  close(iulog)

  call MPI_Finalize(mpi_err)

  if(mpi_err/=MPI_SUCCESS)then
    write(message,'(A,I0,A)') 'ERROR [rank ',instance_parallel%rank, &
                              ']: MPI_Finalize failed'
    call handle_err(mpi_err,message)
  endif

  if(instance_parallel%rank==0) call stop_program(0,'finished parallel parameter evaluation successfully.')


contains


  ! **************************************************************************************************
  ! Spinup SUMMA from a cold state.
  !
  ! Initializes the SUMMA configuration, performs a one-year cold-start spinup, writes the resulting
  ! restart state, and configures that restart file as the initial condition for subsequent model
  ! evaluations. Additional warmup will be required for individual parameters after this point.
  !
  ! This initial warmup also provides access to the parameter data structures.
  ! **************************************************************************************************

  subroutine spinup_from_cold(config,             & ! SUMMA configuration structure
                              domain_parallel,    & ! MPI context for domain parallelism
                              instance_parallel,  & ! MPI context for model-instance parallelism
                              err, message)         ! error code and error message
   
    ! SUMMA global configuration
    USE globalData, only: initConfig,ixRestart
    USE globalData, only: ixRestart_end,ixRestart_never,restart_filename

    ! filenames
    USE summaFileManager, only: SIM_START_TM,SIM_END_TM,MODEL_INITCOND
    USE summaFileManager, only: OUTPUT_PATH
    USE globalData, only: output_fileSuffix

    ! SUMMA initialization and simulation
    USE summa_init,       only: init_config
    USE summa_simulation, only: run_simulation

    implicit none
   
    ! dummy variables
    type(config_info),           intent(out) :: config             ! SUMMA configuration structure
    type(parallel_context_type), intent(in)  :: domain_parallel    ! MPI context for domain parallelism
    type(parallel_context_type), intent(in)  :: instance_parallel  ! MPI context for model-instance parallelism
    integer(i4b),                intent(out) :: err                ! error code
    character(*),                intent(out) :: message            ! error message

    ! local variables
    character(len=4)   :: rankString
    character(len=256) :: log_file
    character(len=256) :: cmessage

    character(len=:), allocatable :: outputFileSuffix_orig
    character(len=:), allocatable :: simStartOriginal,simEndOriginal,spinStart
    character(len=:), allocatable :: timeUnits,flowUnits

    character(len=64), allocatable :: param_name(:)

    real(rkind), allocatable :: param_value(:),timeSim(:),flowSim(:)

    integer(i4b) :: iyear

    err=0
    message='spinup_from_cold/'

    ! zero-length parameter vectors for baseline initialization
    allocate(param_name(0),param_value(0),stat=err)
    if(err/=0)then
      message=trim(message)//'unable to allocate baseline parameter vectors'
      return
    endif


    ! -----------------------------------------------------------------------------------------------
    ! Initialize SUMMA configuration
    ! -----------------------------------------------------------------------------------------------

    call init_config(config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    initConfig=.false.

    ! ---------------------------------------------------------------------------------------
    ! Configure logging files
    ! ---------------------------------------------------------------------------------------

    ! set the logging unit and log files
    iulog = 99
    config%iulog_summa = iulog
 
    write(rankString,'(I4.4)') instance_parallel%rank
    log_file=trim(OUTPUT_PATH)//'logs/'//trim(config%case_name)//'_rank'//rankString//'.log'
    call execute_command_line('mkdir -p "'//trim(OUTPUT_PATH)//'logs"')

    open(unit=iulog,file=trim(log_file),status='replace',action='write')


    ! -----------------------------------------------------------------------------------------------
    ! Perform a one-year cold-start spinup prior to the start of the simulation period
    ! -----------------------------------------------------------------------------------------------

    simStartOriginal=trim(SIM_START_TM)
    simEndOriginal=trim(SIM_END_TM)

    spinStart=simStartOriginal

    read(spinStart(1:4),*) iyear
    write(spinStart(1:4),'(I4.4)') iyear-1

    SIM_START_TM=spinStart
    SIM_END_TM=simStartOriginal

    outputFileSuffix_orig = trim(output_fileSuffix)
    output_fileSuffix     = trim(outputFileSuffix_orig)//'_'//trim(config%case_name)// &
                            '_spinup_rank'//rankString

    ! force writing a restart file at the end of the spinup
    ixRestart=ixRestart_end

    ! run SUMMA for one year following the cold start
    call run_simulation(config,                 & ! SUMMA configuration structure
                        domain_parallel,        & ! MPI context for domain parallelism
                        instance_parallel,      & ! MPI context for model-instance parallelism
                        timeSim,flowSim,        & ! simulated time and streamflow
                        timeUnits,flowUnits,    & ! time and streamflow units
                        param_name,param_value, & ! parameter names and values
                        err,cmessage)             ! error code and message
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


    ! -----------------------------------------------------------------------------------------------
    ! Prepare files/settings for subsequent warm(er) starts. Additonal spinup will still be required.
    ! -----------------------------------------------------------------------------------------------

    if(allocated(restart_filename))then
      MODEL_INITCOND=trim(restart_filename)
    else
      message=trim(message)//'restart filename not defined'
      err=20; return
    endif

    ! restore original simulation settings
    SIM_START_TM=simStartOriginal
    SIM_END_TM=simEndOriginal
    output_fileSuffix=outputFileSuffix_orig
    ixRestart=ixRestart_never

  end subroutine spinup_from_cold


  ! **************************************************************************************************
  ! Sample parameter vectors and evaluate SUMMA.
  !
  ! Initializes the SUMMA-specific parameter specification and model-agnostic search information,
  ! seeds the random number generator, samples feasible parameter vectors, constructs the complete
  ! spatially uniform SUMMA override vector, and evaluates the objective function for each sample.
  ! **************************************************************************************************

  subroutine evaluate_parameter_samples(config,                 & ! SUMMA configuration structure
                                        domain_parallel,        & ! MPI context for domain parallelism
                                        instance_parallel,      & ! MPI context for model-instance parallelism
                                        nLocal,                 & ! parameter trials assigned to this instance
                                        err, message)             ! error code and error message
    
    ! model-agnostic parameter search
    USE parameter_search, only: parameter_spec,parameter_search_info
    USE parameter_search, only: initialize_parameter_search,sample_parameters

    ! SUMMA-specific parameter information and overrides
    USE summa_parameter_spec, only: get_summa_parameter_spec
    USE summa_parameter_spec, only: build_summa_parameter_overrides

    ! objective-function evaluation
    USE summa_simulation, only: evaluate_objective

    ! file paths/names
    USE summaFileManager, only: OUTPUT_PATH
    USE globalData,       only: output_fileSuffix

    ! calibration output
    USE calibration_output_module, only: create_calibration_output
    USE calibration_output_module, only: write_calibration_output
    USE calibration_output_module, only: close_calibration_output

    implicit none

    ! dummy variables
    type(config_info),           intent(inout) :: config             ! SUMMA configuration structure
    type(parallel_context_type), intent(in)    :: domain_parallel    ! MPI context for domain parallelism
    type(parallel_context_type), intent(in)    :: instance_parallel  ! MPI context for model-instance parallelism
    integer(i4b),                intent(in)    :: nLocal             ! parameter trials assigned to this instance
    integer(i4b),                intent(out)   :: err                ! error code
    character(*),                intent(out)   :: message            ! error message

    ! ncid/name for calibration output
    integer(i4b)       :: ncid_calib
    character(len=256) :: calib_file
    
    ! strings to create unique file names
    character(len=4)              :: rankString
    character(len=6)              :: sampleString
    character(len=:), allocatable :: outputFileSuffix_orig

    ! parameter-search information
    type(parameter_spec)        :: param_spec
    type(parameter_search_info) :: search

    ! sampled parameter values
    real(rkind), allocatable :: param_value(:)

    ! complete SUMMA parameter overrides
    character(len=64), allocatable :: param_name(:)
    real(rkind),       allocatable :: param_override(:)

    ! random number generator seed
    integer(i4b)              :: nseed
    integer(i4b), allocatable :: seed(:)

    ! local variables
    integer(i4b) :: i,j
    real(rkind)  :: objective

    integer(i4b) :: startModelRun(8),endModelRun(8)

    character(len=256) :: cmessage

    err=0
    message='evaluate_parameter_samples/'

    ! save the original output file suffix to build unique filenames
    outputFileSuffix_orig=trim(output_fileSuffix)

    ! -----------------------------------------------------------------------------------------------
    ! Initialize parameter search
    ! -----------------------------------------------------------------------------------------------

    call get_summa_parameter_spec(config,param_spec,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    call initialize_parameter_search(param_spec,search,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


    ! -----------------------------------------------------------------------------------------------
    ! Allocate sampled parameter vector
    ! -----------------------------------------------------------------------------------------------

    allocate(param_value(size(search%param_names)),stat=err)

    if(err/=0)then
      message=trim(message)//'unable to allocate sampled parameter vector'
      return
    endif


    ! -----------------------------------------------------------------------------------------------
    ! Print parameter bounds
    ! -----------------------------------------------------------------------------------------------

    write(iulog,'(/,A)') 'Calibration parameter bounds:'
    write(iulog,'(A)')   '  Parameter                         Lower              Upper'

    do i=1,size(search%param_names)
      write(iulog,'(2X,A30,2X,ES16.8,2X,ES16.8)') &
        trim(search%param_names(i)),search%lower(i),search%upper(i)
    enddo


    ! -----------------------------------------------------------------------------------------------
    ! Initialize random number generator
    ! -----------------------------------------------------------------------------------------------

    call random_seed(size=nseed)

    allocate(seed(nseed),stat=err)

    if(err/=0)then
      message=trim(message)//'unable to allocate random-number seed'
      return
    endif

    seed = instance_parallel%rank + 42
    call random_seed(put=seed)


    ! -----------------------------------------------------------------------------------------------
    ! Initialize calibration output
    ! -----------------------------------------------------------------------------------------------
   
    ! define a rank-specific calibration output file 
    write(rankString,'(I4.4)') instance_parallel%rank
    calib_file=trim(OUTPUT_PATH)//trim(config%case_name)//'_calibration_rank'//rankString//'.nc'
    
    call create_calibration_output(calib_file,                  &
                                   param_spec,                  &
                                   instance_parallel%rank,      &
                                   config%case_name,            &
                                   config%calib%metric,         &
                                   config%calib%obs_transform,  &
                                   ncid_calib,                  &
                                   err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


    ! -----------------------------------------------------------------------------------------------
    ! Sample and evaluate parameter vectors
    ! -----------------------------------------------------------------------------------------------

    do j=1,nLocal

      ! record start time for this parameter trial
      call date_and_time(values=startModelRun)

      ! define a unique SUMMA output suffix for this parameter trial
      ! include the case, instance rank, and local sample index to prevent file collisions
      write(sampleString,'(I6.6)') j
      output_fileSuffix=trim(outputFileSuffix_orig)//'_'//trim(config%case_name)// &
                        '_rank'//rankString//'_sample'//sampleString

      ! sample a feasible parameter vector
      call sample_parameters(search,param_value,err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! ensure all parameters in ordered constraints are spatially uniform by
      ! overriding non-sampled parameters with their default values
      call build_summa_parameter_overrides(param_spec,         &
                                           search%param_names, &
                                           param_value,        &
                                           param_name,         &
                                           param_override,     &
                                           err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! evaluate SUMMA for the complete parameter override vector
      call evaluate_objective(config,                 & ! SUMMA configuration structure
                              domain_parallel,        & ! MPI context for domain parallelism
                              instance_parallel,      & ! MPI context for model-instance parallelism
                              param_name,             & ! parameter names
                              param_override,         & ! parameter values
                              objective,              & ! objective function value
                              err,cmessage)             ! error code and message
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! record end time for this parameter trial
      call date_and_time(values=endModelRun)
      
      call write_calibration_output(ncid_calib,      &
                                    j,               &
                                    param_name,      &
                                    param_override,  &
                                    objective,       &
                                    startModelRun,   &
                                    endModelRun,     &
                                    err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      
      ! print sample, sampled parameters, and objective
      write(iulog,'(I6,*(1X,ES16.8))') j,param_value,objective

    enddo

    ! close calibration output
    call close_calibration_output(ncid_calib,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  end subroutine evaluate_parameter_samples


end program summa_driver_opt

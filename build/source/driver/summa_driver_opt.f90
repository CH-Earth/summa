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

  ! data types
  USE nr_type,    only: i4b
  USE summa_type, only: config_info

  ! error handling
  USE summa_util, only: handle_err,stop_program

  implicit none

  ! MPI context
  integer(i4b) :: rank=0,nproc=1

  ! SUMMA configuration
  type(config_info) :: config

  ! number of parameter samples
  integer(i4b), parameter :: nSamples=3

  ! error control
  integer(i4b) :: err=0,mpi_err=0

  character(len=1024) :: message=''
  character(len=256)  :: mpi_message=''


  ! ---------------------------------------------------------------------------------------
  ! Initialize MPI
  ! ---------------------------------------------------------------------------------------

  call MPI_Init(mpi_err)
  call check_mpi(-1,mpi_err,'MPI_Init failed')

  call set_mpi_context(MPI_COMM_WORLD,rank,nproc,mpi_err,mpi_message)
  if(mpi_err/=MPI_SUCCESS) call abort_mpi(rank,trim(mpi_message))


  ! ---------------------------------------------------------------------------------------
  ! Initialize SUMMA
  ! ---------------------------------------------------------------------------------------

  call initialize_summa(config,MPI_COMM_SELF,rank,err,message)
  call handle_err(err,message)


  ! ---------------------------------------------------------------------------------------
  ! Sample parameters and evaluate SUMMA
  ! ---------------------------------------------------------------------------------------

  call evaluate_parameter_samples(config,MPI_COMM_SELF,rank,nSamples,err,message)
  call handle_err(err,message)


  ! ---------------------------------------------------------------------------------------
  ! Finalize MPI
  ! ---------------------------------------------------------------------------------------

  call MPI_Finalize(mpi_err)

  if(mpi_err/=MPI_SUCCESS)then
    write(message,'(A,I0,A)') 'ERROR [rank ',rank,']: MPI_Finalize failed'
    call handle_err(mpi_err,message)
  endif

  if(rank==0) call stop_program(0,'finished parallel parameter evaluation successfully.')


contains


  ! **************************************************************************************************
  ! Initialize SUMMA.
  !
  ! Initializes the SUMMA configuration, performs a one-year cold-start spinup, writes the resulting
  ! restart state, and configures that restart file as the initial condition for subsequent model
  ! evaluations.
  ! **************************************************************************************************

  subroutine initialize_summa(config,comm,rank,err,message)

    ! data types
    USE nr_type,    only: i4b,rkind
    USE summa_type, only: config_info

    ! SUMMA global configuration
    USE globalData, only: initConfig,output_fileSuffix,ixRestart
    USE globalData, only: ixRestart_end,ixRestart_never,restart_filename

    ! model control
    USE summaFileManager, only: SIM_START_TM,SIM_END_TM,MODEL_INITCOND

    ! SUMMA initialization and simulation
    USE summa_init,       only: init_config
    USE summa_simulation, only: run_simulation

    implicit none

    ! dummy variables
    type(config_info), intent(out) :: config

    integer(i4b), intent(in)  :: comm,rank
    integer(i4b), intent(out) :: err

    character(*), intent(out) :: message

    ! local variables
    character(len=4)   :: rankString
    character(len=256) :: cmessage

    character(len=:), allocatable :: outputFileSuffix_orig
    character(len=:), allocatable :: simStartOriginal,simEndOriginal,spinStart
    character(len=:), allocatable :: timeUnits,flowUnits

    character(len=64), allocatable :: param_name(:)

    real(rkind), allocatable :: param_value(:),timeSim(:),flowSim(:)

    integer(i4b) :: iyear

    err=0
    message='initialize_summa/'

    write(rankString,'(I4.4)') rank

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


    ! -----------------------------------------------------------------------------------------------
    ! Perform a one-year cold-start spinup
    ! -----------------------------------------------------------------------------------------------

    simStartOriginal=trim(SIM_START_TM)
    simEndOriginal=trim(SIM_END_TM)

    spinStart=simStartOriginal

    read(spinStart(1:4),*) iyear
    write(spinStart(1:4),'(I4.4)') iyear-1

    SIM_START_TM=spinStart
    SIM_END_TM=simStartOriginal

    outputFileSuffix_orig=trim(output_fileSuffix)
    output_fileSuffix=trim(outputFileSuffix_orig)//'_spinup_rank'//rankString

    ! force writing a restart file at the end of the spinup
    ixRestart=ixRestart_end

    call run_simulation(config,                 &
                        comm,0,1,               &
                        timeSim,flowSim,        &
                        timeUnits,flowUnits,    &
                        param_name,param_value, &
                        err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


    ! -----------------------------------------------------------------------------------------------
    ! Configure spinup restart state
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

  end subroutine initialize_summa


  ! **************************************************************************************************
  ! Sample parameter vectors and evaluate SUMMA.
  !
  ! Initializes the SUMMA-specific parameter specification and model-agnostic search information,
  ! seeds the random number generator, samples feasible parameter vectors, constructs the complete
  ! spatially uniform SUMMA override vector, and evaluates the objective function for each sample.
  ! **************************************************************************************************

  subroutine evaluate_parameter_samples(config,comm,rank,nSample,err,message)

    ! data types
    USE nr_type,    only: i4b,rkind
    USE summa_type, only: config_info

    ! logging
    USE globalData, only: iulog

    ! model-agnostic parameter search
    USE parameter_search, only: parameter_spec,parameter_search_info
    USE parameter_search, only: initialize_parameter_search,sample_parameters

    ! SUMMA-specific parameter information and overrides
    USE summa_parameter_spec, only: get_summa_parameter_spec
    USE summa_parameter_spec, only: build_summa_parameter_overrides

    ! objective-function evaluation
    USE summa_simulation, only: evaluate_objective

    ! calibration output
    USE summaFileManager,          only: OUTPUT_PATH
    USE calibration_output_module, only: create_calibration_output
    USE calibration_output_module, only: write_calibration_output
    USE calibration_output_module, only: close_calibration_output

    implicit none

    ! dummy variables
    type(config_info), intent(inout) :: config
    integer(i4b), intent(in)  :: comm,rank,nSample
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message

    ! calibration output
    integer(i4b)       :: ncid_calib
    character(len=4)   :: rankString
    character(len=256) :: calib_file
    
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

    seed=rank+42
    call random_seed(put=seed)


    ! -----------------------------------------------------------------------------------------------
    ! Initialize calibration output
    ! -----------------------------------------------------------------------------------------------
    
    write(rankString,'(I4.4)') rank
    calib_file=trim(OUTPUT_PATH)//trim(config%case_name)//'_calibration_rank'//rankString//'.nc'
    
    call create_calibration_output(calib_file,                  &
                                   param_spec,                  &
                                   rank,                        &
                                   config%case_name,            &
                                   config%calib%metric,         &
                                   config%calib%obs_transform,  &
                                   ncid_calib,                  &
                                   err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


    ! -----------------------------------------------------------------------------------------------
    ! Sample and evaluate parameter vectors
    ! -----------------------------------------------------------------------------------------------

    do j=1,nSample

      ! record start time for this parameter trial
      call date_and_time(values=startModelRun)

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
      call evaluate_objective(config,                   &
                              comm,0,1,                 &
                              param_name,param_override, &
                              objective,                 &
                              err,cmessage)
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

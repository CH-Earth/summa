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
  USE nr_type,    only: i4b,rkind,lgt
  USE summa_type, only: config_info
  USE summa_type, only: parallel_context_type

  ! logging
  USE globalData,      only: iulog
  USE iso_fortran_env, only: error_unit

  ! MPI
  USE mpi, only: MPI_Init,MPI_Finalize
  USE mpi, only: MPI_Send,MPI_Recv
  USE mpi, only: MPI_Bcast

  USE mpi, only: MPI_INTEGER,MPI_DOUBLE_PRECISION
  USE mpi, only: MPI_CHARACTER
  
  USE mpi, only: MPI_COMM_WORLD,MPI_COMM_SELF
  USE mpi, only: MPI_ANY_TAG,MPI_ANY_SOURCE,MPI_STATUS_SIZE
  USE mpi, only: MPI_SOURCE,MPI_TAG
  USE mpi, only: MPI_SUCCESS
  
  USE mpi_context, only: set_mpi_context
  USE error_utils, only: check_mpi,abort_mpi

  ! SUMMA globals
  USE globalData, only: initConfig
  USE globalData, only: restart_filename

  ! SUMMA paths/filenames
  USE summaFileManager, only: OUTPUT_PATH
  USE summaFileManager, only: MODEL_INITCOND

  ! SUMMA subroutines/functions
  USE summa_init, only: init_config
  USE summa_util, only: handle_err,stop_program
  USE summa_parameter_evaluation, only: spinup_from_cold

  implicit none

  ! SUMMA configuration
  type(config_info) :: config

  ! MPI contexts for domain and model-instance parallelism
  type(parallel_context_type) :: domain_parallel
  type(parallel_context_type) :: instance_parallel

  ! MPI message tags
  integer(i4b), parameter :: tag_work = 1
  integer(i4b), parameter :: tag_done = 2
  integer(i4b), parameter :: tag_stop = 3

  ! rank-specific logging
  character(len=4)    :: rankString
  character(len=256)  :: log_file

  ! number of parameter samples
  integer(i4b), parameter :: nSamples=1000

  ! error control
  integer(i4b)        :: err=0
  integer(i4b)        :: mpi_err=0
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

  ! dynamic parameter evaluation requires one dispatcher and at least one worker
  if(instance_parallel%size < 2)then
    call abort_mpi(instance_parallel%rank, &
                   'parallel parameter evaluation requires at least two MPI ranks')
  endif

  ! each model instance runs independently without domain parallelism
  domain_parallel%comm=MPI_COMM_SELF
  domain_parallel%rank=0
  domain_parallel%size=1

  ! ---------------------------------------------------------------------------------------
  ! Configure SUMMA
  ! ---------------------------------------------------------------------------------------
  
  ! use stderr until the configuration and output paths are known
  iulog=error_unit
 
  ! read the configuration files to establish file paths, simulation settings, and calibration options
  call init_config(config,err,message)
  call handle_err(err,message)
  
  initConfig=.false.
  
  ! configure rank-specific logging
  iulog=99
  config%iulog_summa=iulog
  
  write(rankString,'(I4.4)') instance_parallel%rank
  log_file=trim(OUTPUT_PATH)//'logs/'//trim(config%case_name)// &
           '_rank'//rankString//'.log'
  
  call execute_command_line('mkdir -p "'//trim(OUTPUT_PATH)//'logs"')
  open(unit=iulog,file=trim(log_file),status='replace',action='write')

  ! ---------------------------------------------------------------------------------------
  ! Spin up SUMMA from a cold state
  ! ---------------------------------------------------------------------------------------

  ! perform a one-year cold-start spinup to establish the initial model state
  call spinup_from_cold(config,             & ! SUMMA configuration structure
                        domain_parallel,    & ! MPI context for domain parallelism
                        instance_parallel,  & ! MPI context for model-instance parallelism
                        err,message)          ! error code and message
  call handle_err(err,message)

  ! broadcast the rank-0 restart filename to all model instances
  call MPI_Bcast(restart_filename,len(restart_filename),MPI_CHARACTER,0, &
                 instance_parallel%comm,mpi_err)
  call check_mpi(instance_parallel%rank,mpi_err, &
                 'unable to broadcast restart filename')
  
  ! use the common cold-start restart state as the initial conditions
  MODEL_INITCOND=trim(restart_filename)

  ! ---------------------------------------------------------------------------------------
  ! Evaluate parameter samples
  ! ---------------------------------------------------------------------------------------

  call evaluate_parameter_samples(config,                 & ! SUMMA configuration structure
                                  domain_parallel,        & ! MPI context for domain parallelism
                                  instance_parallel,      & ! MPI context for model-instance parallelism
                                  nSamples,               & ! total number of parameter samples
                                  err,message)             ! error code and message
  call handle_err(err,message)

  ! ---------------------------------------------------------------------------------------
  ! Finalize MPI
  ! ---------------------------------------------------------------------------------------
  
  ! close the rank-specific logging file
  close(iulog);  iulog=error_unit

  call MPI_Finalize(mpi_err)
  call check_mpi(instance_parallel%rank,mpi_err,'MPI_Finalize failed')

  if(instance_parallel%rank==0) &
    call stop_program(0,'finished parallel parameter evaluation successfully.')

contains

  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------

  ! **************************************************************************************************
  ! Dynamically distribute and evaluate parameter samples.
  !
  ! Rank 0 generates parameter samples and assigns work to available workers. Each worker evaluates
  ! one sample at a time and returns the resulting objective value. Workers that finish early are
  ! immediately assigned additional samples, reducing load imbalance from variable model runtimes.
  ! **************************************************************************************************

  subroutine evaluate_parameter_samples(config,                &
                                        domain_parallel,       &
                                        instance_parallel,     &
                                        nSamples,              &
                                        err,message)

    ! parameter search types
    USE parameter_search, only: parameter_spec,parameter_search_info
    USE parameter_search, only: sample_parameters
    
    ! SUMMA parameter evaluation subroutines
    USE summa_parameter_evaluation, only: initialize_parameter_evaluation
    USE summa_parameter_evaluation, only: evaluate_parameter_sample

    ! calibration output
    USE calibration_output_module, only: close_calibration_output

    implicit none

    ! dummy variables
    type(config_info),           intent(inout) :: config             ! SUMMA configuration structure
    type(parallel_context_type), intent(in)    :: domain_parallel    ! MPI context for domain parallelism
    type(parallel_context_type), intent(in)    :: instance_parallel  ! MPI context for model-instance parallelism
    integer(i4b),                intent(in)    :: nSamples           ! total number of parameter trials
    integer(i4b),                intent(out)   :: err                ! error code
    character(*),                intent(out)   :: message            ! error message

    ! parameter-search information
    type(parameter_spec)        :: param_spec
    type(parameter_search_info) :: search

    ! sampled parameter values
    real(rkind), allocatable :: param_value(:)

    ! MPI work-queue state
    integer(i4b) :: worker
    integer(i4b) :: sample_id
    integer(i4b) :: next_sample
    integer(i4b) :: nComplete
    logical(lgt) :: stop_worker

    ! worker-local calibration output
    integer(i4b) :: ncid_calib
    integer(i4b) :: local_sample

    ! objective value
    real(rkind) :: objective

    ! error control
    integer(i4b)        :: mpi_err
    character(len=256)  :: cmessage

    err=0
    message='evaluate_parameter_samples/'

    ! -----------------------------------------------------------------------------------------------
    ! Initialize parameter evaluation
    ! -----------------------------------------------------------------------------------------------

    ! build the parameter specification and search information, allocate the sampled parameter vector,
    ! initialize the dispatcher random-number generator, and create worker calibration output files
    call initialize_parameter_evaluation(config,                 & ! SUMMA configuration structure
                                         instance_parallel,      & ! MPI context for model-instance parallelism
                                         param_spec,search,      & ! parameter specification and search information
                                         param_value,ncid_calib, & ! sampled parameter vector and calibration output NetCDF ID
                                         err,cmessage)             ! error code and message

    if(err/=0)then
      message=trim(message)//trim(cmessage)
      call abort_mpi(instance_parallel%rank,trim(message))
    endif

    ! -----------------------------------------------------------------------------------------------
    ! Dispatcher
    ! -----------------------------------------------------------------------------------------------

    if(instance_parallel%rank == 0)then

      next_sample=1
      nComplete=0

      ! assign one initial parameter sample to each available worker
      do worker=1,instance_parallel%size-1

        if(next_sample <= nSamples)then

          ! generate a feasible parameter vector
          call sample_parameters(search,param_value,err,cmessage)
          
          if(err/=0)then
            message=trim(message)//trim(cmessage)
            call abort_mpi(instance_parallel%rank,trim(message))
          endif

          ! send the sample index and parameter vector to this worker
          call send_sample(worker,next_sample,param_value, &
                           instance_parallel%comm,mpi_err)
          call check_mpi(instance_parallel%rank,mpi_err, &
                         'unable to send parameter sample')

          next_sample=next_sample+1

        else

          ! no work is available for this worker
          call send_stop(worker,instance_parallel%comm,mpi_err)
          call check_mpi(instance_parallel%rank,mpi_err, &
                         'unable to send stop message')

        endif

      enddo


      ! wait for completed trials and immediately refill available workers
      do while(nComplete < nSamples)

        ! receive the objective value from whichever worker finishes next
        call receive_objective(objective,worker, &
                               instance_parallel%comm,mpi_err)
        call check_mpi(instance_parallel%rank,mpi_err, &
                       'unable to receive objective value')

        nComplete=nComplete+1

        ! immediately give the completed worker another sample if work remains
        if(next_sample <= nSamples)then

          ! generate the next feasible parameter vector
          call sample_parameters(search,param_value,err,cmessage)
          
          if(err/=0)then
            message=trim(message)//trim(cmessage)
            call abort_mpi(instance_parallel%rank,trim(message))
          endif

          ! send the next sample to the worker that just became available
          call send_sample(worker,next_sample,param_value, &
                           instance_parallel%comm,mpi_err)
          call check_mpi(instance_parallel%rank,mpi_err, &
                         'unable to send parameter sample')

          next_sample=next_sample+1

        else

          ! all samples have been dispatched, so this worker is finished
          call send_stop(worker,instance_parallel%comm,mpi_err)
          call check_mpi(instance_parallel%rank,mpi_err, &
                         'unable to send stop message')

        endif

      enddo


    ! -----------------------------------------------------------------------------------------------
    ! Workers
    ! -----------------------------------------------------------------------------------------------

    else

      local_sample=0

      do

        ! wait for either another parameter sample or a stop instruction
        call receive_sample(sample_id,param_value,stop_worker, &
                            instance_parallel%comm,mpi_err)
        call check_mpi(instance_parallel%rank,mpi_err, &
                       'unable to receive parameter sample')

        if(stop_worker) exit

        ! advance the local record index in this worker's calibration file
        local_sample=local_sample+1

        ! evaluate the supplied parameter sample
        call evaluate_parameter_sample(config,                 &
                                       domain_parallel,        &
                                       instance_parallel,      &
                                       param_spec,search,      &
                                       sample_id,local_sample, &
                                       param_value,ncid_calib, &
                                       objective,              &
                                       err, cmessage)
          
        if(err/=0)then
          message=trim(message)//trim(cmessage)
          call abort_mpi(instance_parallel%rank,trim(message))
        endif

        ! return the objective value and become available for additional work
        call send_objective(objective,instance_parallel%comm,mpi_err)
        call check_mpi(instance_parallel%rank,mpi_err, &
                       'unable to send objective value')

      enddo

      ! close rank-specific calibration output
      call close_calibration_output(ncid_calib,err,cmessage)

      if(err/=0)then
        message=trim(message)//trim(cmessage)
        call abort_mpi(instance_parallel%rank,trim(message))
      endif

    endif

  end subroutine evaluate_parameter_samples

  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------

  ! **************************************************************************************************
  ! Send one parameter sample to a worker.
  ! **************************************************************************************************

  subroutine send_sample(worker,sample_id,param_value,comm,mpi_err)

    integer(i4b), intent(in)  :: worker
    integer(i4b), intent(in)  :: sample_id
    real(rkind),  intent(in)  :: param_value(:)
    integer(i4b), intent(in)  :: comm
    integer(i4b), intent(out) :: mpi_err

    ! first send the work instruction and global sample index
    call MPI_Send(sample_id,1,MPI_INTEGER,worker,tag_work,comm,mpi_err)
    if(mpi_err/=MPI_SUCCESS) return

    ! then send the corresponding parameter vector
    call MPI_Send(param_value,size(param_value),MPI_DOUBLE_PRECISION, &
                  worker,tag_work,comm,mpi_err)

  end subroutine send_sample


  ! **************************************************************************************************
  ! Tell a worker that no additional parameter samples remain.
  ! **************************************************************************************************

  subroutine send_stop(worker,comm,mpi_err)

    integer(i4b), intent(in)  :: worker
    integer(i4b), intent(in)  :: comm
    integer(i4b), intent(out) :: mpi_err

    integer(i4b) :: dummy

    dummy=0

    call MPI_Send(dummy,1,MPI_INTEGER,worker,tag_stop,comm,mpi_err)

  end subroutine send_stop


  ! **************************************************************************************************
  ! Receive either a parameter sample or a stop instruction from rank 0.
  ! **************************************************************************************************

  subroutine receive_sample(sample_id,param_value,stop_worker,comm,mpi_err)

    integer(i4b), intent(out) :: sample_id
    real(rkind),  intent(out) :: param_value(:)
    logical(lgt), intent(out) :: stop_worker
    integer(i4b), intent(in)  :: comm
    integer(i4b), intent(out) :: mpi_err

    integer(i4b) :: status(MPI_STATUS_SIZE)

    stop_worker=.false.

    ! receive either a work or stop instruction
    call MPI_Recv(sample_id,1,MPI_INTEGER,0,MPI_ANY_TAG,comm,status,mpi_err)
    if(mpi_err/=MPI_SUCCESS) return

    select case(status(MPI_TAG))

      case (tag_stop)
        stop_worker=.true.
        return
    
      case (tag_work)
        call MPI_Recv(param_value,size(param_value),MPI_DOUBLE_PRECISION, &
                      0,tag_work,comm,status,mpi_err)
    
      case default
        ! all message tags are defined internally, so an unknown tag is fatal
        call abort_mpi(instance_parallel%rank,'unknown MPI tag')
    
    end select

  end subroutine receive_sample


  ! **************************************************************************************************
  ! Return a completed objective-function value to rank 0.
  ! **************************************************************************************************

  subroutine send_objective(objective,comm,mpi_err)

    real(rkind),  intent(in)  :: objective
    integer(i4b), intent(in)  :: comm
    integer(i4b), intent(out) :: mpi_err

    call MPI_Send(objective,1,MPI_DOUBLE_PRECISION,0,tag_done,comm,mpi_err)

  end subroutine send_objective


  ! **************************************************************************************************
  ! Receive an objective-function value from whichever worker finishes first.
  ! **************************************************************************************************

  subroutine receive_objective(objective,worker,comm,mpi_err)

    real(rkind),  intent(out) :: objective
    integer(i4b), intent(out) :: worker
    integer(i4b), intent(in)  :: comm
    integer(i4b), intent(out) :: mpi_err

    integer(i4b) :: status(MPI_STATUS_SIZE)

    ! wait for the next completed parameter trial
    call MPI_Recv(objective,1,MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,tag_done, &
                  comm,status,mpi_err)
    if(mpi_err/=MPI_SUCCESS) return

    ! identify the worker that is now available for additional work
    worker=status(MPI_SOURCE)

  end subroutine receive_objective

end program summa_driver_opt

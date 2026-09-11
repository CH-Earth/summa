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
  USE iso_fortran_env, only: output_unit

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

  ! SUMMA parameter information
  USE parameter_search, only: parameter_spec,parameter_search_info

  ! SUMMA subroutines/functions
  USE summa_init,   only: init_config
  USE summa_spinup, only: spinup_from_cold
  USE summa_util,   only: stop_program

  USE calibration_output_module, only: create_calibration_output
  USE calibration_output_module, only: close_calibration_output

  implicit none

  type(config_info) :: config                              ! SUMMA configuration information

  type(parallel_context_type) :: domain_parallel           ! MPI context for domain parallelism
  type(parallel_context_type) :: instance_parallel         ! MPI context for model-instance parallelism

  integer(i4b), parameter :: tag_work=1                    ! MPI tag for work messages
  integer(i4b), parameter :: tag_done=2                    ! MPI tag for completed work
  integer(i4b), parameter :: tag_stop=3                    ! MPI tag for stop messages

  character(len=4)   :: rankString                         ! MPI rank formatted as a string
  character(len=256) :: log_file                           ! rank-specific log file

  type(parameter_spec)                   :: param_spec     ! complete SUMMA parameter specification
  type(parameter_search_info)            :: search         ! sampled parameter search information
  character(len=64), allocatable         :: param_name(:)  ! complete SUMMA parameter names

  character(len=*), parameter :: sampling_method='dds'     ! parameter sampling method
  integer(i4b),     parameter :: nSamples=5000             ! total number of parameter samples

  real(rkind), allocatable :: x_best(:)                    ! best parameter vector found so far
  real(rkind)              :: F_best                       ! objective value associated with x_best
  integer(i4b)             :: sample_best                  ! sample index associated with x_best

  integer(i4b)       :: ncid_calib                         ! NetCDF identifier for calibration output
  character(len=256) :: calib_file                         ! calibration output file

  integer(i4b)        :: err=0                             ! SUMMA error code
  integer(i4b)        :: mpi_err=0                         ! MPI error code
  character(len=1024) :: message=''                        ! SUMMA error message
  character(len=256)  :: mpi_message=''                    ! MPI error message

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
  if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))
  
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
  if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))

  ! broadcast the rank-0 restart filename to all model instances
  call MPI_Bcast(restart_filename,len(restart_filename),MPI_CHARACTER,0, &
                 instance_parallel%comm,mpi_err)
  call check_mpi(instance_parallel%rank,mpi_err, &
                 'unable to broadcast restart filename')
  
  ! use the common cold-start restart state as the initial conditions
  MODEL_INITCOND=trim(restart_filename)

  ! ---------------------------------------------------------------------------------------
  ! Initialize parameter sampling
  ! ---------------------------------------------------------------------------------------

  ! Initialize parameter-evaluation state that is invariant across samples and shared by all ranks
  call initialize_parameter_evaluation(config,            & ! SUMMA configuration
                                       instance_parallel, & ! MPI instance-parallel context
                                       param_spec,search, & ! parameter specification and search information
                                       param_name,        & ! complete SUMMA parameter-name vector
                                       x_best,F_best,     & ! current best solution and objective value
                                       sample_best,       & ! sample index associated with current best solution
                                       err,message)         ! error code and message
  if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))

  ! ---------------------------------------------------------------------------------------
  ! Create calibration output file
  ! ---------------------------------------------------------------------------------------
  
  if(instance_parallel%rank == 0)then
    calib_file=trim(OUTPUT_PATH)//trim(config%case_name)//'_calibration.nc'
    call create_calibration_output(calib_file,param_spec,nSamples,instance_parallel%size-1,           &
                                   config%case_name,config%calib%metric,config%calib%obs_transform,   &
                                   ncid_calib,err,message)
    if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))
  else
    ncid_calib=-1
  endif

  ! ---------------------------------------------------------------------------------------
  ! Evaluate parameter samples
  ! ---------------------------------------------------------------------------------------

  call evaluate_parameter_samples(config,                     & ! SUMMA configuration structure
                                  domain_parallel,            & ! MPI context for domain parallelism
                                  instance_parallel,          & ! MPI context for model-instance parallelism
                                  param_spec,search,          & ! parameter specification and search information
                                  param_name,ncid_calib,      & ! parameter names and calibration output
                                  x_best,F_best,sample_best,  & ! current best solution and objective value; sample index
                                  nSamples,err,message)         ! total number of parameter samples; error code and message
  if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))

  ! ---------------------------------------------------------------------------------------
  ! Close calibration output
  ! ---------------------------------------------------------------------------------------

  if(instance_parallel%rank == 0)then
    call close_calibration_output(ncid_calib,err,message)
    if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))
  endif

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
  ! Initialize parameter evaluation.
  !
  ! Constructs the SUMMA parameter specification and parameter-search information used for all
  ! parameter evaluations. Builds the complete parameter-name vector, including sampled parameters
  ! and non-sampled parameters required by calibration constraints. On the dispatcher rank, also
  ! initializes the random-number generator used for parameter sampling.
  !
  ! Initializes parameter information shared across evaluations and establishes
  ! dispatcher-owned state used to track the current best solution.!
  ! **************************************************************************************************
  
  subroutine initialize_parameter_evaluation(config,instance_parallel, &
                                             param_spec,search,param_name, &
                                             x_best,F_best,sample_best, &
                                             err,message)

    USE parameter_search, only: initialize_parameter_search
    USE summa_parameter_spec, only: get_summa_parameter_spec
  
    implicit none

    type(config_info),           intent(in)  :: config             ! SUMMA configuration information
    type(parallel_context_type), intent(in)  :: instance_parallel  ! MPI context for model-instance parallelism

    type(parameter_spec),        intent(out) :: param_spec         ! complete SUMMA parameter specification
    type(parameter_search_info), intent(out) :: search             ! sampled parameter search information
    character(len=64), allocatable, intent(out) :: param_name(:)   ! complete SUMMA parameter names

    real(rkind), allocatable, intent(out) :: x_best(:)             ! best decision-variable vector
    real(rkind),              intent(out) :: F_best                ! best objective value
    integer(i4b),             intent(out) :: sample_best           ! sample index associated with best objective

    integer(i4b), intent(out) :: err                               ! error code
    character(*), intent(out) :: message                           ! error message

    integer(i4b)              :: i                                 ! parameter index
    integer(i4b)              :: nseed                             ! random-number seed vector size
    integer(i4b), allocatable :: seed(:)                           ! random-number seed vector

    character(len=256) :: cmessage                                 ! message returned by called routines
  
    err=0
    message='initialize_parameter_evaluation/'
  
    ! build the SUMMA parameter specification from the calibration configuration
    call get_summa_parameter_spec(config,param_spec,err,cmessage)
    if(err/=0)then
      message=trim(message)//trim(cmessage)
      return
    endif
  
    ! construct and validate the model-agnostic parameter-search information
    call initialize_parameter_search(param_spec,search,err,cmessage)
    if(err/=0)then
      message=trim(message)//trim(cmessage)
      return
    endif
 
    ! construct the complete invariant SUMMA parameter-name vector
    allocate(param_name(size(param_spec%params)),stat=err)
  
    if(err/=0)then
      message=trim(message)//'unable to allocate parameter-name vector'
      return
    endif
  
    do i=1,size(param_spec%params)
      param_name(i)=param_spec%params(i)%name
    enddo
  
    if(instance_parallel%rank == 0)then

      ! initialize random-number generator on the dispatcher
      
      call random_seed(size=nseed)

      allocate(seed(nseed),stat=err)

      if(err/=0)then
        message=trim(message)//'unable to allocate random-number seed'
        return
      endif

      seed=42
      call random_seed(put=seed)

      ! initialize best parameter vector on the dispatcher
      
      allocate(x_best(size(search%param_names)),stat=err)

      if(err/=0)then
        message=trim(message)//'unable to allocate best parameter vector'
        return
      endif

      x_best=0._rkind
      F_best=-huge(1._rkind)
      sample_best=0

    endif

  end subroutine initialize_parameter_evaluation


  ! **************************************************************************************************
  ! Dynamically distribute and evaluate parameter samples.
  !
  ! Rank 0 generates parameter samples and assigns work to available workers. Each worker evaluates
  ! one sample at a time and returns the resulting objective value. Workers that finish early are
  ! immediately assigned additional samples, reducing load imbalance from variable model runtimes.
  ! **************************************************************************************************

  subroutine evaluate_parameter_samples(config,                    &
                                        domain_parallel,           &
                                        instance_parallel,         &
                                        param_spec,search,         &
                                        param_name,ncid_calib,     &
                                        x_best,F_best,sample_best, &
                                        nSamples,err,message)

    ! parameter search
    USE parameter_search, only: parameter_spec,parameter_search_info
    
    ! objective-function evaluation
    USE summa_parameter_sampling, only: generate_dds_sample
    USE summa_parameter_sampling, only: generate_parameter_sample
    USE summa_simulation,         only: evaluate_objective
    
    ! calibration output
    USE calibration_output_module, only: write_calibration_output

    implicit none

    ! dummy variables
    type(config_info),              intent(inout) :: config             ! SUMMA configuration structure

    type(parallel_context_type),    intent(in)    :: domain_parallel    ! MPI context for domain parallelism
    type(parallel_context_type),    intent(in)    :: instance_parallel  ! MPI context for model-instance parallelism
    
    type(parameter_spec),           intent(in)    :: param_spec         ! complete SUMMA parameter specification
    type(parameter_search_info),    intent(in)    :: search             ! parameter-search configuration and metadata
    character(len=64),              intent(in)    :: param_name(:)      ! complete SUMMA parameter-name vector
    integer(i4b),                   intent(in)    :: ncid_calib         ! calibration output NetCDF file ID
    
    real(rkind), allocatable,       intent(inout) :: x_best(:)          ! current best decision-variable vector
    real(rkind),                    intent(inout) :: F_best             ! objective value associated with x_best
    integer(i4b),                   intent(inout) :: sample_best        ! sample index associated with F_best
    
    integer(i4b),                   intent(in)    :: nSamples           ! total number of parameter trials
    
    integer(i4b),                   intent(out)   :: err                ! error code
    character(*),                   intent(out)   :: message            ! error message

    ! sampled parameter values
    real(rkind),       allocatable   :: param_value(:)
    real(rkind),       allocatable   :: param_override(:)

    ! complete parameter samples and parameter overrides retained by rank 0
    real(rkind),       allocatable   :: param_samples(:,:)
    real(rkind),       allocatable   :: param_overrides(:,:)

    ! parameter-evaluation timing
    integer(i4b),      allocatable  :: startModelRun(:,:)
    integer(i4b),      allocatable  :: endModelRun(:,:)

    ! MPI work-queue state
    integer(i4b) :: worker
    integer(i4b) :: worker_sample(instance_parallel%size-1)
    integer(i4b) :: sample_id
    integer(i4b) :: next_sample
    integer(i4b) :: nComplete
    logical(lgt) :: stop_worker

    ! objective value
    real(rkind) :: objective

    ! error control
    integer(i4b)        :: mpi_err
    character(len=256)  :: cmessage

    err=0
    message='evaluate_parameter_samples/'

    ! -----------------------------------------------------------------------------------------------
    ! Allocate local arrays
    ! -----------------------------------------------------------------------------------------------

    ! available on all ranks
    allocate(param_override(size(param_spec%params)),stat=err)
    if(err/=0)then
      message=trim(message)//'unable to allocate parameter override vector'
      return
    endif

    ! rank 0 responsible for parameter sampling and dispatch
    if(instance_parallel%rank == 0)then

      allocate(param_value(size(search%param_names)),stat=err)
      if(err/=0)then
        message=trim(message)//'unable to allocate sampled parameter vector'
        return
      endif

      allocate(param_samples(size(search%param_names),nSamples),stat=err)
      if(err/=0)then
        message=trim(message)//'unable to allocate sampled parameter storage'
        return
      endif

      allocate(param_overrides(size(param_spec%params),nSamples),stat=err)
      if(err/=0)then
        message=trim(message)//'unable to allocate parameter override storage'
        return
      endif

      allocate(startModelRun(8,nSamples), stat=err)
      if(err/=0)then
        message=trim(message)//'unable to allocate parameter start-time storage'
        return
      endif

      allocate(endModelRun(8,nSamples), stat=err)
      if(err/=0)then
        message=trim(message)//'unable to allocate parameter end-time storage'
        return
      endif
    
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

          ! generate the next parameter sample and complete SUMMA override vector
          call generate_parameter_sample(param_spec,search,       &
                                         param_value,param_override, &
                                         err,cmessage)
          if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

          ! save parameter samples (rank 0)
          param_samples(:,next_sample)=param_value       ! retain sampled decision-variable vector
          param_overrides(:,next_sample)=param_override  ! retain complete SUMMA override vector
          call date_and_time(values=startModelRun(:,next_sample))

          ! send the sample index and parameter vector to this worker
          call send_sample(worker,next_sample,param_override, &
                           instance_parallel%comm,mpi_err)
          call check_mpi(instance_parallel%rank,mpi_err, &
                         'unable to send parameter sample')

          worker_sample(worker)=next_sample
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

        sample_id=worker_sample(worker)
        call date_and_time(values=endModelRun(:,sample_id))

        ! update best parameter sample
        if(objective > F_best)then
          F_best=objective                     ! update best objective value
          x_best=param_samples(:,sample_id)    ! update best decision-variable vector
          sample_best=sample_id                ! record sample associated with current best

          write(output_unit,'(A,I0,A,F14.6)') 'new DDS best: sample=',sample_best,', objective=',F_best

        endif

        call write_calibration_output(ncid_calib,                   &
                                      sample_id, worker,            &
                                      param_name,                   &
                                      param_overrides(:,sample_id), &
                                      objective,                    &
                                      startModelRun(:,sample_id),   &
                                      endModelRun(:,sample_id),     &
                                      err,cmessage)
        
        if(err/=0)then
          message=trim(message)//trim(cmessage)
          call abort_mpi(instance_parallel%rank,trim(message))
        endif


        nComplete=nComplete+1

        ! immediately give the completed worker another sample if work remains
        if(next_sample <= nSamples)then

          ! generate the next parameter sample and complete SUMMA override vector
          select case(trim(sampling_method))

            case ('dds')

              call generate_dds_sample(param_spec,search,             &
                                       x_best,                         &
                                       next_sample,nSamples,           &
                                       param_value,param_override,     &
                                       err,cmessage)
              if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
         
            case ('random')
         
              call generate_parameter_sample(param_spec,search,       &
                                             param_value,param_override, &
                                             err,cmessage)
              if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
         
            case default
         
              message=trim(message)//'unknown parameter sampling method: '//trim(sampling_method)
              err=20; return
       
          end select


          ! save parameter samples (rank 0)
          param_samples(:,next_sample)=param_value       ! retain sampled decision-variable vector
          param_overrides(:,next_sample)=param_override  ! retain complete SUMMA override vector
          call date_and_time(values=startModelRun(:,next_sample))

          ! send the next sample to the worker that just became available
          call send_sample(worker,next_sample,param_override, &
                           instance_parallel%comm,mpi_err)
          call check_mpi(instance_parallel%rank,mpi_err, &
                         'unable to send parameter sample')

          worker_sample(worker)=next_sample
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

      do

        ! wait for either another parameter sample or a stop instruction
        call receive_sample(sample_id,param_override,stop_worker, &
                            instance_parallel%comm,mpi_err)
        call check_mpi(instance_parallel%rank,mpi_err, &
                       'unable to receive parameter sample')

        if(stop_worker) exit

        ! run SUMMA and evaluate the objective function
        call evaluate_objective(config,                              & ! SUMMA configuration structure
                                domain_parallel,instance_parallel,   & ! MPI context for model domain and model-instance parallelism
                                sample_id,param_name,param_override, & ! complete parameter overrides
                                objective,err,cmessage)                ! objective function value and error control
        if(err/=0)then
          message=trim(message)//trim(cmessage)
          call abort_mpi(instance_parallel%rank,trim(message))
        endif

        ! return the objective value and become available for additional work
        call send_objective(objective, &
                            instance_parallel%comm,mpi_err)
        call check_mpi(instance_parallel%rank,mpi_err, &
                       'unable to send objective value')

      enddo

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

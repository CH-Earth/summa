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
! SUMMA parameter sampling
!
! Provides routines for initializing, generating, distributing, and evaluating SUMMA parameter
! samples. Parameter sampling and the asynchronous MPI work queue are encapsulated here so that
! the optimization driver is responsible only for case-level orchestration.
! **************************************************************************************************
module summa_parameter_sampling

  ! data types
  USE nr_type,    only: i4b,rkind,lgt
  USE data_types, only: var_d
  USE data_types, only: var_dlength
  USE summa_type, only: config_info
  USE summa_type, only: parallel_context_type

  ! logging
  USE iso_fortran_env, only: output_unit

  ! MPI
  USE mpi, only: MPI_Send,MPI_Recv
  USE mpi, only: MPI_INTEGER,MPI_DOUBLE_PRECISION
  USE mpi, only: MPI_ANY_TAG,MPI_ANY_SOURCE
  USE mpi, only: MPI_STATUS_SIZE
  USE mpi, only: MPI_SOURCE,MPI_TAG
  USE mpi, only: MPI_SUCCESS

  USE error_utils, only: check_mpi,abort_mpi

  ! parameter-search information
  USE parameter_search, only: parameter_spec
  USE parameter_search, only: parameter_search_info
  USE parameter_search, only: search_state_type

  implicit none
  private

  ! MPI tags for asynchronous parameter evaluation
  integer(i4b), parameter :: tag_work=1
  integer(i4b), parameter :: tag_done=2
  integer(i4b), parameter :: tag_stop=3
  integer(i4b), parameter :: tag_fail=4

  ! parameter sampling method
  character(len=*), parameter :: sampling_method='dds'

  ! every trial of a calibration, kept by rank 0
  type :: trial_ledger
    real(rkind),  allocatable :: value(:,:)       ! sampled decision vector (nSampled,nSamples)
    real(rkind),  allocatable :: override(:,:)    ! complete SUMMA override vector (nParam,nSamples)
    real(rkind),  allocatable :: objective(:,:)   ! one value per calibration target (nTarget,nSamples)
    logical(lgt), allocatable :: failed(:)        ! the model refused or failed on the trial (nSamples)
    integer(i4b), allocatable :: start_time(:,:)  ! dispatch time (8,nSamples)
    integer(i4b), allocatable :: end_time(:,:)    ! completion time (8,nSamples)
    integer(i4b)              :: nComplete = 0    ! trials returned so far
  end type trial_ledger

  ! public parameter-sampling interface
  public :: check_search_settings
  public :: initialize_parameter_evaluation
  public :: dispatch_parameter_samples

contains

  ! **************************************************************************************************
  ! Initialize parameter evaluation.
  !
  ! Constructs the SUMMA parameter specification and parameter-search information used for all
  ! parameter evaluations. Builds the complete parameter-name vector, including sampled parameters
  ! and non-sampled parameters required by calibration constraints. On the dispatcher rank, also
  ! initializes the random-number generator used for parameter sampling.
  !
  ! Initializes parameter information shared across evaluations and establishes the dispatcher-owned
  ! state used by the parameter-search algorithm.
  ! **************************************************************************************************
  subroutine initialize_parameter_evaluation(config,instance_parallel,param_spec,search,param_name, &
                                             search_state,sample_best,err,message,mpar_init,bpar_init)
    USE parameter_search, only: initialize_parameter_search
    USE summa_parameter_spec, only: get_summa_parameter_spec
    implicit none
    type(config_info),           intent(in)  :: config             ! SUMMA configuration information
    type(parallel_context_type), intent(in)  :: instance_parallel  ! MPI context for model-instance parallelism
    type(var_dlength), optional, intent(in)  :: mpar_init          ! initialized HRU parameters, single-domain runs only
    type(var_d),       optional, intent(in)  :: bpar_init          ! initialized basin parameters, single-domain runs only
    type(parameter_spec),        intent(out) :: param_spec         ! complete SUMMA parameter specification
    type(parameter_search_info), intent(out) :: search             ! sampled parameter search information
    character(len=64), allocatable, intent(out) :: param_name(:)   ! complete SUMMA parameter names
    type(search_state_type),     intent(out) :: search_state       ! information from previously evaluated parameter sets
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
    call get_summa_parameter_spec(config, param_spec, err, cmessage, mpar_init, bpar_init)
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
      seed=42 + config%run_index - 1
      call random_seed(put=seed)

      ! initialize parameter-search state on the dispatcher
      search_state%elite_fitness=-huge(1._rkind)
      sample_best=0

    endif

  end subroutine initialize_parameter_evaluation

  ! **************************************************************************************************
  ! Evaluate the calibration's parameter samples.
  !
  ! Workers evaluate whatever trial rank 0 sends until told to stop. Rank 0 runs the configured search:
  ! DDS over the whole budget at once, or NSGA-II one generation at a time. Both hand their trials to
  ! the same asynchronous dispatcher, so the worker pool persists throughout.
  ! **************************************************************************************************
  subroutine dispatch_parameter_samples(config,                    &
                                        domain_parallel,           &
                                        instance_parallel,         &
                                        param_spec,search,         &
                                        param_name,ncid_calib,     &
                                        search_state,sample_best,  &
                                        nSamples,err,message)
    USE summa_simulation, only: n_calibration_targets
    implicit none
    type(config_info),              intent(inout) :: config             ! SUMMA configuration structure
    type(parallel_context_type),    intent(in)    :: domain_parallel    ! MPI context for domain parallelism
    type(parallel_context_type),    intent(in)    :: instance_parallel  ! MPI context for model-instance parallelism
    type(parameter_spec),           intent(in)    :: param_spec         ! complete SUMMA parameter specification
    type(parameter_search_info),    intent(in)    :: search             ! parameter-search configuration and metadata
    character(len=64),              intent(in)    :: param_name(:)      ! complete SUMMA parameter-name vector
    integer(i4b),                   intent(in)    :: ncid_calib         ! calibration output NetCDF file ID
    type(search_state_type),        intent(inout) :: search_state       ! information from previously evaluated parameter sets
    integer(i4b),                   intent(inout) :: sample_best        ! sample index associated with elite fitness
    integer(i4b),                   intent(in)    :: nSamples           ! total number of parameter trials
    integer(i4b),                   intent(out)   :: err                ! error code
    character(*),                   intent(out)   :: message            ! error message
    type(trial_ledger)  :: ledger                                       ! every trial, kept by rank 0
    integer(i4b)        :: nTarget                                      ! objectives per trial
    character(len=256)  :: cmessage

    err=0
    message='dispatch_parameter_samples/'

    ! every rank reads the same configuration, so the objective messages agree in length on both sides
    nTarget=max(n_calibration_targets(config),1)

    if(instance_parallel%rank /= 0)then
      call serve_trials(config,domain_parallel,instance_parallel,param_name, &
                        size(param_spec%params),nTarget,err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      return
    endif

    allocate(ledger%value(size(search%param_names),nSamples),  &
             ledger%override(size(param_spec%params),nSamples), &
             ledger%objective(nTarget,nSamples),                 &
             ledger%failed(nSamples),                            &
             ledger%start_time(8,nSamples),                      &
             ledger%end_time(8,nSamples),stat=err)
    if(err/=0)then
      message=trim(message)//'unable to allocate the trial record'
      return
    endif
    ledger%failed=.false.

    select case(trim(config%calib%algorithm))
      case('dds')
        call dispatch_trials(config,instance_parallel,param_spec,search,param_name,ncid_calib, &
                             1,nSamples,nSamples,ledger,search_state,sample_best,err,cmessage)
      case('nsga2')
        call run_nsga2(config,instance_parallel,param_spec,search,param_name,ncid_calib, &
                       nSamples,ledger,search_state,sample_best,err,cmessage)
    end select
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    call release_workers(instance_parallel)

  end subroutine dispatch_parameter_samples

  ! **************************************************************************************************
  ! Dispatch trials first..last to the workers and wait for all of them (rank 0 only).
  !
  ! Each worker gets one trial at a time and is refilled the moment it reports, which absorbs the
  ! spread in model runtimes. DDS proposes each trial as a worker frees up, from the best so far;
  ! NSGA-II has already written the whole generation into the ledger.
  ! **************************************************************************************************
  subroutine dispatch_trials(config,instance_parallel,param_spec,search,param_name,ncid_calib, &
                             first,last,nSamples,ledger,search_state,sample_best,err,message)
    USE calibration_output_module, only: write_calibration_output
    implicit none
    type(config_info),           intent(in)    :: config
    type(parallel_context_type), intent(in)    :: instance_parallel
    type(parameter_spec),        intent(in)    :: param_spec
    type(parameter_search_info), intent(in)    :: search
    character(len=64),           intent(in)    :: param_name(:)
    integer(i4b),                intent(in)    :: ncid_calib
    integer(i4b),                intent(in)    :: first,last        ! trials to evaluate
    integer(i4b),                intent(in)    :: nSamples          ! total sample budget
    type(trial_ledger),          intent(inout) :: ledger
    type(search_state_type),     intent(inout) :: search_state
    integer(i4b),                intent(inout) :: sample_best
    integer(i4b),                intent(out)   :: err
    character(*),                intent(out)   :: message
    integer(i4b) :: worker
    integer(i4b) :: worker_sample(instance_parallel%size-1)
    integer(i4b) :: sample_id
    integer(i4b) :: next_sample
    integer(i4b) :: nActive
    integer(i4b) :: mpi_err
    real(rkind)  :: objective(size(ledger%objective,1))
    logical(lgt) :: failed
    character(len=256) :: cmessage

    err=0
    message='dispatch_trials/'
    next_sample=first
    nActive=0

    ! one trial to each worker
    do worker=1,instance_parallel%size-1
      if(next_sample > last) exit
      call start_trial(worker,next_sample)
      if(err/=0) return
      next_sample=next_sample+1
    enddo

    ! wait for completed trials and immediately refill available workers
    do while(nActive > 0)
      call receive_objective(objective,failed,worker,instance_parallel%comm,mpi_err)
      call check_mpi(instance_parallel%rank,mpi_err,'unable to receive objective value')
      nActive=nActive-1
      ledger%nComplete=ledger%nComplete+1
      sample_id=worker_sample(worker)
      call date_and_time(values=ledger%end_time(:,sample_id))
      ledger%objective(:,sample_id)=objective
      ledger%failed(sample_id)=failed
      if(failed) write(output_unit,'(A,I0,A,I0,A)') 'trial ',sample_id,' failed on worker ',worker, &
                                                    '; it ranks last on every target'

      call record_trial(config,sample_id,ledger,search_state,sample_best,err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      call write_calibration_output(ncid_calib,                    &
                                    sample_id,worker,              &
                                    param_name,                    &
                                    ledger%override(:,sample_id),  &
                                    ledger%objective(:,sample_id), &
                                    ledger%failed(sample_id),      &
                                    ledger%start_time(:,sample_id),&
                                    ledger%end_time(:,sample_id),  &
                                    err,cmessage)
      if(err/=0)then
        message=trim(message)//trim(cmessage)
        call abort_mpi(instance_parallel%rank,trim(message))
      endif

      if(next_sample <= last)then
        call start_trial(worker,next_sample)
        if(err/=0) return
        next_sample=next_sample+1
      endif
    enddo

  contains

    ! propose trial i if the search has not already, and send it to a worker
    subroutine start_trial(worker,i)
      integer(i4b), intent(in) :: worker,i
      call propose_trial(config,param_spec,search,i,nSamples,search_state,ledger,err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      call date_and_time(values=ledger%start_time(:,i))
      call send_sample(worker,i,ledger%override(:,i),instance_parallel%comm,mpi_err)
      call check_mpi(instance_parallel%rank,mpi_err,'unable to send parameter sample')
      worker_sample(worker)=i
      nActive=nActive+1
    end subroutine start_trial

  end subroutine dispatch_trials

  ! **************************************************************************************************
  ! Propose trial i into the ledger. DDS samples at random while it has no elite to perturb: the first
  ! round of trials, one per worker, and after that until a trial succeeds; its progress counts the
  ! trials returned so far. NSGA-II proposes a generation before dispatching it.
  ! **************************************************************************************************
  subroutine propose_trial(config,param_spec,search,i,nSamples,search_state,ledger,err,message)
    implicit none
    type(config_info),           intent(in)    :: config
    type(parameter_spec),        intent(in)    :: param_spec
    type(parameter_search_info), intent(in)    :: search
    integer(i4b),                intent(in)    :: i                 ! trial to propose
    integer(i4b),                intent(in)    :: nSamples
    type(search_state_type),     intent(in)    :: search_state
    type(trial_ledger),          intent(inout) :: ledger
    integer(i4b),                intent(out)   :: err
    character(*),                intent(out)   :: message
    character(len=256) :: cmessage

    err=0
    message='propose_trial/'
    if(trim(config%calib%algorithm) /= 'dds') return

    ! NOTE: before any trial returns, nComplete=1 so that p = 1 - log(1) / log(nSamples) = 1
    call generate_parameter_sample(sampling_method,param_spec,search,search_state,max(ledger%nComplete,1),nSamples, &
                                   ledger%value(:,i),ledger%override(:,i),err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  end subroutine propose_trial

  ! **************************************************************************************************
  ! Take in the objectives of trial i. DDS searches on one number, so it keeps the trial with the best
  ! weighted, consistently oriented sum of its targets; with one target of weight one, the metric.
  ! **************************************************************************************************
  subroutine record_trial(config,i,ledger,search_state,sample_best,err,message)
    USE summa_simulation, only: scalarize_objectives
    implicit none
    type(config_info),       intent(in)    :: config
    integer(i4b),            intent(in)    :: i
    type(trial_ledger),      intent(in)    :: ledger
    type(search_state_type), intent(inout) :: search_state
    integer(i4b),            intent(inout) :: sample_best
    integer(i4b),            intent(out)   :: err
    character(*),            intent(out)   :: message
    real(rkind) :: F_sample

    err=0
    message='record_trial/'
    if(trim(config%calib%algorithm) /= 'dds' .or. ledger%failed(i)) return
    F_sample=scalarize_objectives(config,ledger%objective(:,i))
    if(F_sample > search_state%elite_fitness)then
      search_state%elite_fitness=F_sample
      if(.not.allocated(search_state%elite_individual))then
        allocate(search_state%elite_individual(size(ledger%value,1)),stat=err)
        if(err/=0)then
          message=trim(message)//'unable to allocate elite individual'
          return
        endif
      endif
      search_state%elite_individual=ledger%value(:,i)
      sample_best=i
      write(output_unit,'(A,I0,A,F14.6)') 'new elite: sample=',sample_best,', objective=',search_state%elite_fitness
    endif

  end subroutine record_trial

  ! **************************************************************************************************
  ! Worker loop: evaluate each trial rank 0 sends, return its objectives, until told to stop.
  !
  ! A trial the model refuses or fails on is saved as a repro case and reported to rank 0 as failed,
  ! and the worker goes on to the next; any other error ends the calibration.
  ! **************************************************************************************************
  subroutine serve_trials(config,domain_parallel,instance_parallel,param_name,nParam,nTarget,err,message)
    USE globalData,                only: iulog
    USE summa_simulation,          only: evaluate_objectives
    USE calibration_output_module, only: write_failed_trial
    implicit none
    type(config_info),           intent(inout) :: config
    type(parallel_context_type), intent(in)    :: domain_parallel
    type(parallel_context_type), intent(in)    :: instance_parallel
    character(len=64),           intent(in)    :: param_name(:)
    integer(i4b),                intent(in)    :: nParam            ! length of the override vector
    integer(i4b),                intent(in)    :: nTarget           ! objectives per trial
    integer(i4b),                intent(out)   :: err
    character(*),                intent(out)   :: message
    real(rkind)  :: param_override(nParam)
    real(rkind)  :: objective(nTarget)
    integer(i4b) :: sample_id
    integer(i4b) :: mpi_err
    logical(lgt) :: stop_worker
    logical(lgt) :: failed
    character(len=:), allocatable :: dir
    character(len=1024) :: cmessage

    err=0
    message='serve_trials/'
    do
      call receive_sample(sample_id,param_override,stop_worker, &
                          instance_parallel%comm,instance_parallel%rank,mpi_err)
      call check_mpi(instance_parallel%rank,mpi_err,'unable to receive parameter sample')
      if(stop_worker) exit

      call evaluate_objectives(config,                              & ! SUMMA configuration structure
                               domain_parallel,instance_parallel,   & ! MPI context for model domain and model-instance parallelism
                               sample_id,param_name,param_override, & ! complete parameter overrides
                               objective,err,cmessage,failed)         ! objective value per target and error control
      if(err/=0 .and. .not.failed)then
        message=trim(message)//trim(cmessage)
        call abort_mpi(instance_parallel%rank,trim(message))
      endif

      if(failed)then
        call write_failed_trial(config,sample_id,instance_parallel%rank,param_name,param_override,cmessage,dir)
        write(output_unit,'(A,I0,A)') 'trial ',sample_id,' failed; saved in '//dir//': '//trim(cmessage)
        write(iulog,'(A,I0,A)') 'trial ',sample_id,' failed; saved in '//dir//': '//trim(cmessage)
        err=0
      endif

      call send_objective(objective,failed,instance_parallel%comm,mpi_err)
      call check_mpi(instance_parallel%rank,mpi_err,'unable to send objective value')
    enddo

  end subroutine serve_trials

  ! **************************************************************************************************
  ! Tell every worker the calibration is over (rank 0 only).
  ! **************************************************************************************************
  subroutine release_workers(instance_parallel)
    implicit none
    type(parallel_context_type), intent(in) :: instance_parallel
    integer(i4b) :: worker
    integer(i4b) :: mpi_err

    do worker=1,instance_parallel%size-1
      call send_stop(worker,instance_parallel%comm,mpi_err)
      call check_mpi(instance_parallel%rank,mpi_err,'unable to send stop message')
    enddo

  end subroutine release_workers

  ! **************************************************************************************************
  ! NSGA-II over the calibration targets (rank 0 only).
  !
  ! The first generation is a random population; each later one breeds as many offspring, evaluates
  ! them, and keeps the best population_size of parents and offspring together. Objectives are
  ! oriented so smaller is better; weights play no part. Writes each generation's population and,
  ! at the end, which of all the trials are non-dominated.
  ! **************************************************************************************************
  subroutine run_nsga2(config,instance_parallel,param_spec,search,param_name,ncid_calib, &
                       nSamples,ledger,search_state,sample_best,err,message)
    USE nsga2,                     only: nondominated_sort,nondominated_set
    USE nsga2,                     only: crowding_distance,select_survivors,make_offspring
    USE summa_simulation,          only: oriented_objectives
    USE summa_parameter_spec,      only: build_summa_parameter_overrides
    USE calibration_output_module, only: write_generation_output,write_pareto_front
    implicit none
    type(config_info),           intent(in)    :: config
    type(parallel_context_type), intent(in)    :: instance_parallel
    type(parameter_spec),        intent(in)    :: param_spec
    type(parameter_search_info), intent(in)    :: search
    character(len=64),           intent(in)    :: param_name(:)
    integer(i4b),                intent(in)    :: ncid_calib
    integer(i4b),                intent(in)    :: nSamples
    type(trial_ledger),          intent(inout) :: ledger
    type(search_state_type),     intent(inout) :: search_state
    integer(i4b),                intent(inout) :: sample_best
    integer(i4b),                intent(out)   :: err
    character(*),                intent(out)   :: message
    real(rkind),  allocatable :: f(:,:)          ! oriented objectives of every trial
    integer(i4b), allocatable :: member(:)       ! sample index of each member of the population
    integer(i4b), allocatable :: rank(:)         ! their fronts
    real(rkind),  allocatable :: distance(:)     ! their crowding distances
    integer(i4b), allocatable :: combined(:)     ! parents then offspring
    integer(i4b), allocatable :: keep(:)         ! survivors' positions in combined
    integer(i4b) :: nPop,nGen,iGen,first,last,i
    real(rkind)  :: pm
    character(len=256) :: cmessage

    err=0
    message='run_nsga2/'
    associate(settings => config%calib%nsga2)
    nPop=settings%population_size
    nGen=nSamples/nPop
    pm=settings%mutation_probability
    if(pm < 0._rkind) pm=1._rkind/real(size(search%param_names),rkind)
    allocate(f(size(ledger%objective,1),nSamples),member(nPop),rank(nPop),distance(nPop), &
             combined(2*nPop),keep(nPop),stat=err)
    if(err/=0)then
      message=trim(message)//'unable to allocate the population'
      return
    endif

    do iGen=1,nGen
      first=(iGen-1)*nPop+1
      last=iGen*nPop

      ! propose the generation: a random population, then offspring of the one before
      if(iGen == 1)then
        do i=first,last
          call generate_parameter_sample('random',param_spec,search,search_state,1,nSamples, &
                                         ledger%value(:,i),ledger%override(:,i),err,cmessage)
          if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
        enddo
      else
        call make_offspring(search,ledger%value(:,member),rank,distance,                         &
                            settings%crossover_probability,settings%crossover_eta,pm,settings%mutation_eta, &
                            ledger%value(:,first:last),err,cmessage)
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
        do i=first,last
          call build_summa_parameter_overrides(param_spec,search%param_names,ledger%value(:,i), &
                                               ledger%override(:,i),err,cmessage)
          if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
        enddo
      endif

      call dispatch_trials(config,instance_parallel,param_spec,search,param_name,ncid_calib, &
                           first,last,nSamples,ledger,search_state,sample_best,err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      do i=first,last
        if(ledger%failed(i))then
          f(:,i)=0.1_rkind*huge(1._rkind)
        else
          f(:,i)=oriented_objectives(config,ledger%objective(:,i))
        endif
      enddo

      ! select the population the next generation breeds from
      if(iGen == 1)then
        member=[(i,i=first,last)]
        call nondominated_sort(f(:,member),rank)
        call crowding_distance(f(:,member),rank,distance)
      else
        combined=[member,(i,i=first,last)]
        call select_survivors(f(:,combined),nPop,keep,rank,distance)
        member=combined(keep)
      endif

      call write_generation_output(ncid_calib,iGen,first,last,member,rank,distance,err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      call report_generation(config,iGen,nGen,ledger%objective(:,member),rank)
    enddo

    call write_pareto_front(ncid_calib,nondominated_set(f),err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    end associate

  end subroutine run_nsga2

  ! **************************************************************************************************
  ! One line per generation: the size of the non-dominated front and each target's best in it.
  ! **************************************************************************************************
  subroutine report_generation(config,iGen,nGen,objective,rank)
    USE metrics, only: metric_is_maximized
    implicit none
    type(config_info), intent(in) :: config
    integer(i4b),      intent(in) :: iGen,nGen
    real(rkind),       intent(in) :: objective(:,:)  ! population objectives (nTarget,nPop)
    integer(i4b),      intent(in) :: rank(:)
    character(len=1024) :: line
    character(len=32)   :: val
    real(rkind)  :: best
    integer(i4b) :: iTarget

    write(line,'(A,I0,A,I0,A,I0,A)') 'NSGA-II generation ',iGen,' of ',nGen,': ',count(rank==1),' non-dominated'
    do iTarget=1,min(size(objective,1),size(config%calib%targets))
      if(metric_is_maximized(config%calib%targets(iTarget)%metric))then
        best=maxval(objective(iTarget,:),mask=rank==1)
      else
        best=minval(objective(iTarget,:),mask=rank==1)
      endif
      write(val,'(F14.6)') best
      line=trim(line)//'; best '//trim(config%calib%targets(iTarget)%name)//' '// &
           trim(config%calib%targets(iTarget)%metric)//' = '//trim(adjustl(val))
    enddo
    write(output_unit,'(A)') trim(line)

  end subroutine report_generation

  ! **************************************************************************************************
  ! Check the search algorithm and, for NSGA-II, that the sample budget is a whole number of generations.
  ! **************************************************************************************************
  subroutine check_search_settings(config,err,message)
    USE summa_simulation, only: n_calibration_targets
    implicit none
    type(config_info), intent(in)  :: config
    integer(i4b),      intent(out) :: err
    character(*),      intent(out) :: message

    err=0
    message='check_search_settings/'
    select case(trim(config%calib%algorithm))
      case('dds')
      case('nsga2')
        associate(nsga2 => config%calib%nsga2, nSamples => config%calib%n_samples)
        if(n_calibration_targets(config) < 1)then
          message=trim(message)//'NSGA-II needs at least one [[calibration.target]] or an [observations] file'
          err=20; return
        endif
        if(nsga2%population_size < 2)then
          message=trim(message)//'calibration.nsga2.population_size must be at least 2'
          err=20; return
        endif
        if(mod(nSamples,nsga2%population_size) /= 0 .or. nSamples < 2*nsga2%population_size)then
          write(message,'(A,I0,A,I0,A,I0)') trim(message)//'calibration.n_samples = ',nSamples, &
            ' must be a whole number of generations of population_size = ',nsga2%population_size, &
            ', at least two: a multiple of it no smaller than ',2*nsga2%population_size
          err=20; return
        endif
        if(nsga2%crossover_probability < 0._rkind .or. nsga2%crossover_probability > 1._rkind .or. &
           nsga2%mutation_probability > 1._rkind)then
          message=trim(message)//'calibration.nsga2 probabilities must lie between 0 and 1'
          err=20; return
        endif
        if(nsga2%crossover_eta < 0._rkind .or. nsga2%mutation_eta < 0._rkind)then
          message=trim(message)//'calibration.nsga2 distribution indices must not be negative'
          err=20; return
        endif
        end associate
      case default
        message=trim(message)//'unknown calibration.algorithm "'//trim(config%calib%algorithm)//'"; use dds or nsga2'
        err=20; return
    end select

  end subroutine check_search_settings


  ! **************************************************************************************************
  ! Generate a SUMMA parameter sample.
  !
  ! Generates one feasible parameter vector using the configured parameter-search strategy and
  ! constructs the complete SUMMA parameter override vector for model evaluation. The sampled
  ! parameter vector contains only parameters included in the search, whereas the override vector
  ! also includes non-sampled parameters required by calibration constraints.
  ! **************************************************************************************************
  subroutine generate_parameter_sample(method,param_spec,search,search_state,nComplete,nSamples, &
                                       param_value,param_override,err,message)
    USE parameter_search, only: parameter_spec
    USE parameter_search, only: parameter_search_info
    USE parameter_search, only: search_state_type
    USE parameter_search, only: generate_search_sample
    USE summa_parameter_spec, only: build_summa_parameter_overrides

    implicit none

    character(*),                intent(in)  :: method
    type(parameter_spec),        intent(in)  :: param_spec
    type(parameter_search_info), intent(in)  :: search
    type(search_state_type),     intent(in)  :: search_state
    integer(i4b),                intent(in)  :: nComplete
    integer(i4b),                intent(in)  :: nSamples
    real(rkind),                 intent(out) :: param_value(:)
    real(rkind),                 intent(out) :: param_override(:)
    integer(i4b),                intent(out) :: err
    character(*),                intent(out) :: message

    character(len=256) :: cmessage

    err=0
    message='generate_parameter_sample/'

    ! generate the next sampled parameter vector
    call generate_search_sample(method,search,search_state,nComplete, &
                                nSamples,param_value,err,cmessage)
    if(err/=0)then
      message=trim(message)//trim(cmessage)
      return
    endif

    ! construct the complete SUMMA parameter override vector
    call build_summa_parameter_overrides(param_spec,search%param_names,param_value, &
                                         param_override,err,cmessage)
    if(err/=0)then
      message=trim(message)//trim(cmessage)
      return
    endif

  end subroutine generate_parameter_sample

  ! --------------------------------------------------------------------------------------------------
  ! --- PRIVATE HELPER ROUTINES ----------------------------------------------------------------------
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
    call MPI_Send(param_value,size(param_value),MPI_DOUBLE_PRECISION, worker,tag_work,comm,mpi_err)

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
  subroutine receive_sample(sample_id,param_value,stop_worker,comm,rank,mpi_err)
    integer(i4b), intent(out) :: sample_id
    real(rkind),  intent(out) :: param_value(:)
    logical(lgt), intent(out) :: stop_worker
    integer(i4b), intent(in)  :: comm
    integer(i4b), intent(in)  :: rank
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
        call MPI_Recv(param_value,size(param_value),MPI_DOUBLE_PRECISION, 0,tag_work,comm,status,mpi_err)
      case default
        ! all message tags are defined internally, so an unknown tag is fatal
        call abort_mpi(rank,'unknown MPI tag')
    end select

  end subroutine receive_sample

  ! **************************************************************************************************
  ! Return the objective value of every calibration target to rank 0, tagged as failed if the trial was.
  ! **************************************************************************************************
  subroutine send_objective(objective,failed,comm,mpi_err)
    real(rkind),  intent(in)  :: objective(:)
    logical(lgt), intent(in)  :: failed
    integer(i4b), intent(in)  :: comm
    integer(i4b), intent(out) :: mpi_err

    call MPI_Send(objective,size(objective),MPI_DOUBLE_PRECISION,0,merge(tag_fail,tag_done,failed),comm,mpi_err)

  end subroutine send_objective

  ! **************************************************************************************************
  ! Receive the objective values from whichever worker finishes first.
  !
  ! Every rank builds its objective vector from the same configuration, so the message is the length
  ! the dispatcher expects.
  ! **************************************************************************************************
  subroutine receive_objective(objective,failed,worker,comm,mpi_err)
    real(rkind),  intent(out) :: objective(:)
    logical(lgt), intent(out) :: failed
    integer(i4b), intent(out) :: worker
    integer(i4b), intent(in)  :: comm
    integer(i4b), intent(out) :: mpi_err
    integer(i4b) :: status(MPI_STATUS_SIZE)

    ! wait for the next completed parameter trial
    call MPI_Recv(objective,size(objective),MPI_DOUBLE_PRECISION,MPI_ANY_SOURCE,MPI_ANY_TAG, comm,status,mpi_err)
    if(mpi_err/=MPI_SUCCESS) return
    failed=(status(MPI_TAG) == tag_fail)

    ! identify the worker that is now available for additional work
    worker=status(MPI_SOURCE)

  end subroutine receive_objective

end module summa_parameter_sampling

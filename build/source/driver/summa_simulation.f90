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
module summa_simulation

USE nr_type,    only: i4b, i8b, rkind, lgt
USE data_types, only: gru_double
USE data_types, only: gru_hru_dom_doublevec
USE summa_type, only: config_info
USE summa_type, only: summa1_type_dec
USE summa_type, only: parallel_context_type

USE summa_init, only: summa_initialize
USE summa_setup, only: summa_paramSetup
USE summa_restart, only: summa_readRestart
USE summa_forcing, only: summa_readForcing
USE summa_modelRun, only: summa_runPhysics
USE summa_writeOutput, only: summa_writeOutputFiles
USE globalData, only: integerMissing
USE globalData, only: realMissing
USE globalData, only: iulog

USE build_options, only: mizuroute_active
USE build_options, only: openwq_active

#ifdef MIZUROUTE_ACTIVE
USE mizuroute_coupling,        only: get_mizuroute_streamflow
USE mizuroute_coupling,        only: get_mizuroute_reach_index
USE finalize_mizuroute_module, only: finalize_mizuroute
#endif

#ifdef OPENWQ_ACTIVE
USE summa_openwq, only: openwq_init
USE summa_openwq, only: openwq_run_time_start
USE summa_openwq, only: openwq_run_space_step
USE summa_openwq, only: openwq_run_time_end
USE summa_openwq, only: openwq_finalize
#endif

! module-level data structure to share configurations
implicit none
private

public :: run_simulation
public :: evaluate_objective
public :: evaluate_objectives
public :: n_calibration_targets
public :: scalarize_objectives
public :: oriented_objectives

contains

  ! **************************************************************************************************
  ! Run a complete SUMMA simulation and return the simulated streamflow time series.
  ! The interface is model agnostic: model-specific initialization, parameter updates,
  ! simulation, and finalization are handled internally.
  ! **************************************************************************************************
  subroutine run_simulation(config,                 & ! SUMMA configuration structure
                            domain_parallel,        & ! MPI context for domain parallelism
                            instance_parallel,      & ! MPI context for model-instance parallelism
                            timeSim,flowSim,        & ! simulated time and streamflow
                            timeUnits,flowUnits,    & ! time and streamflow units
                            param_name,param_value, & ! parameter names and values
                            mpar_out,bpar_out,      & ! model parameter structures
                            err, message)             ! error code and message
    ! dummy arguments
    type(config_info)                    , intent(inout) :: config
    type(parallel_context_type)          , intent(in)    :: domain_parallel
    type(parallel_context_type)          , intent(in)    :: instance_parallel
    real(rkind)             , allocatable, intent(out)   :: timeSim(:)
    real(rkind)             , allocatable, intent(out)   :: flowSim(:)
    character(len=:)        , allocatable, intent(out)   :: timeUnits
    character(len=:)        , allocatable, intent(out)   :: flowUnits
    character(*)                         , intent(in)    :: param_name(:)
    real(rkind)                          , intent(in)    :: param_value(:)
    type(gru_hru_dom_doublevec), optional, intent(out)   :: mpar_out
    type(gru_double)           , optional, intent(out)   :: bpar_out
    integer(i4b)                         , intent(out)   :: err
    character(*)                         , intent(out)   :: message
    ! locals
    type(summa1_type_dec), allocatable :: summa1_struc(:)
    integer(i4b), parameter            :: n=1
    character(len=1024)                :: cmessage

    err=0
    message='run_simulation/'
    allocate(summa1_struc(n),stat=err)
    if(err/=0)then
      message=trim(message)//'problem allocating top-level summa structure'
      return
    endif

    ! populate domain and model-instance parallel contexts
    summa1_struc(n)%domain_parallel=domain_parallel
    summa1_struc(n)%instance_parallel=instance_parallel

    call initialize_summa(config, summa1_struc(n), param_name, param_value, &
                          mpar_out=mpar_out, bpar_out=bpar_out, err=err, message=cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    call run_summa(summa1_struc(n), timeSim,flowSim, timeUnits,flowUnits, err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    call finalize_summa(summa1_struc(n),err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  end subroutine run_simulation

  ! **************************************************************************************************
  ! Report the number of calibration targets, which is the number of objectives per parameter trial.
  !
  ! A configuration that names no target is a single-objective calibration of the one observed series
  ! in [observations], which the configuration reader has already expressed as one target.  A run with
  ! nothing to calibrate against reports zero.
  ! **************************************************************************************************
  pure function n_calibration_targets(config) result(nTarget)
    type(config_info), intent(in) :: config
    integer(i4b)                  :: nTarget

    nTarget = 0
    if(allocated(config%calib%targets)) nTarget = size(config%calib%targets)

  end function n_calibration_targets

  ! **************************************************************************************************
  ! Collapse the per-target objectives into the single value a single-objective search maximizes.
  !
  ! Targets are oriented before they are combined: efficiencies count as they stand and error metrics
  ! count negatively, so a larger scalar is a better fit whatever mix of metrics the targets use.  A
  ! single target of weight one therefore returns exactly the metric itself, which is what DDS
  ! maximized before calibration grew a target list.
  ! **************************************************************************************************
  pure function scalarize_objectives(config,objective) result(F)
    use metrics, only: metric_is_maximized
    type(config_info), intent(in) :: config
    real(rkind),       intent(in) :: objective(:)   ! one value per calibration target
    real(rkind)                   :: F              ! scalarized objective
    integer(i4b) :: iTarget

    F = 0._rkind
    if(.not.allocated(config%calib%targets)) return

    do iTarget=1,min(size(config%calib%targets),size(objective))
      if(metric_is_maximized(config%calib%targets(iTarget)%metric))then
        F = F + config%calib%targets(iTarget)%weight*objective(iTarget)
      else
        F = F - config%calib%targets(iTarget)%weight*objective(iTarget)
      endif
    enddo

  end function scalarize_objectives

  ! **************************************************************************************************
  ! Orient each target's objective so that smaller is better, as a Pareto search compares them:
  ! efficiencies are negated, error metrics kept. A value that is not finite is the worst possible.
  ! **************************************************************************************************
  pure function oriented_objectives(config,objective) result(f)
    use metrics, only: metric_is_maximized
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    type(config_info), intent(in) :: config
    real(rkind),       intent(in) :: objective(:)   ! one value per calibration target
    real(rkind)                   :: f(size(objective))
    integer(i4b) :: iTarget

    f = objective
    if(.not.allocated(config%calib%targets)) return

    do iTarget=1,min(size(config%calib%targets),size(objective))
      if(metric_is_maximized(config%calib%targets(iTarget)%metric)) f(iTarget) = -objective(iTarget)
      if(.not.ieee_is_finite(f(iTarget))) f(iTarget) = 0.1_rkind*huge(1._rkind)
    enddo

  end function oriented_objectives

  ! **************************************************************************************************
  ! Evaluate the objective function for a specified parameter vector.
  !
  ! Retained for callers that want one number: it evaluates every calibration target and returns the
  ! scalarized objective.  With a single target, which is what a configuration naming no target has,
  ! that number is the metric itself.
  ! **************************************************************************************************
  subroutine evaluate_objective(config,                            & ! SUMMA configuration structure
                                domain_parallel,                   & ! MPI context for domain parallelism
                                instance_parallel,                 & ! MPI context for model-instance parallelism
                                sample_id,param_name,param_value,  & ! sample ID + parameter names and values
                                metric,                            & ! objective function value
                                err, message)                        ! error code and message
    ! dummy arguments
    type(config_info),           intent(inout) :: config
    type(parallel_context_type), intent(in)    :: domain_parallel
    type(parallel_context_type), intent(in)    :: instance_parallel
    integer(i4b), intent(in)  :: sample_id
    character(*), intent(in)  :: param_name(:)
    real(rkind),  intent(in)  :: param_value(:)
    real(rkind),  intent(out) :: metric
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    ! locals
    real(rkind), allocatable :: objective(:)   ! one value per calibration target
    integer(i4b)             :: nTarget        ! number of calibration targets
    character(len=1024)      :: cmessage       ! error message of downwind routine

    err=0
    message='evaluate_objective/'
    metric=realMissing

    nTarget=max(n_calibration_targets(config),1)
    allocate(objective(nTarget),stat=err)
    if(err/=0)then
      message=trim(message)//'problem allocating the objective vector'
      return
    endif

    call evaluate_objectives(config,                              &
                             domain_parallel,instance_parallel,   &
                             sample_id, param_name,param_value,   &
                             objective,                           &
                             err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    metric=scalarize_objectives(config,objective)

  end subroutine evaluate_objective

  ! **************************************************************************************************
  ! Return the position of a variable over a spatial unit in a list, or integerMissing when it is not
  ! there.
  !
  ! Two targets scoring the same variable over the same unit - one storage series against two
  ! products, say - collect it once.  The same variable over two units is two series.
  ! **************************************************************************************************
  pure function find_series(names,units,unitIds,name,unitKind,unitId) result(iName)
    use simulated_series, only: ix_unit_domain
    character(*), intent(in) :: names(:)
    integer(i4b), intent(in) :: units(:)
    integer(i8b), intent(in) :: unitIds(:)
    character(*), intent(in) :: name
    integer(i4b), intent(in) :: unitKind
    integer(i8b), intent(in) :: unitId
    integer(i4b)             :: iName
    integer(i4b) :: i

    iName=integerMissing
    do i=1,size(names)
      if(trim(names(i))/=trim(name)) cycle
      if(units(i)/=unitKind) cycle
      if(unitKind/=ix_unit_domain .and. unitIds(i)/=unitId) cycle
      iName=i
      return
    enddo

  end function find_series

  ! **************************************************************************************************
  ! Evaluate every calibration target for a specified parameter vector.
  !
  ! The routine initializes SUMMA, runs the model once, and then scores that one simulation against
  ! each configured target in turn, returning one objective value per target.  Running the model once
  ! for all targets is the point: the targets are different views of the same simulation, so they must
  ! come from the same run to be comparable.
  !
  ! Objectives are returned as the metrics themselves, in the order the targets were configured.  It
  ! is the caller that decides what to do with them - maximize a weighted sum, or rank a population by
  ! Pareto dominance - so nothing here assumes which way any of them points.
  ! **************************************************************************************************
  subroutine evaluate_objectives(config,                            & ! SUMMA configuration structure
                                 domain_parallel,                   & ! MPI context for domain parallelism
                                 instance_parallel,                 & ! MPI context for model-instance parallelism
                                 sample_id, param_name,param_value, & ! sample ID + parameter names and values
                                 objective,                         & ! objective value for each target
                                 err, message,                      & ! error code and message
                                 trial_failed)                        ! .true. if the error came from the parameter set's own run
    use iso_fortran_env, only: output_unit, error_unit
    use globalData, only: ncid
    USE globalData, only: output_fileSuffix
    use var_lookup, only: iLookFREQ
    use globalData,              only: numtim
    use read_flowobs_module,     only: read_observations
    use timeseries_alignment,    only: align_timeseries
    use metrics,                 only: compute_metric
    use write_evaluation_module, only: write_evaluation
    use series_transform,        only: accumulate_series
    use series_transform,        only: remove_baseline_mean
    use series_transform,        only: balance_between_extremes
    use simulated_series,        only: sim_series_type
    use simulated_series,        only: init_simulated_series
    use simulated_series,        only: find_simulated_series
    use simulated_series,        only: is_routed_streamflow
    use simulated_series,        only: spatial_unit_index
    use simulated_series,        only: ix_unit_domain, ix_unit_gru, ix_unit_reach
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    ! dummy arguments
    type(config_info),           intent(inout) :: config
    type(parallel_context_type), intent(in)    :: domain_parallel
    type(parallel_context_type), intent(in)    :: instance_parallel
    integer(i4b), intent(in)  :: sample_id
    character(*), intent(in)  :: param_name(:)
    real(rkind),  intent(in)  :: param_value(:)
    real(rkind),  intent(out) :: objective(:)
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    logical(lgt), optional, intent(out) :: trial_failed
    ! locals
    logical(lgt)                       :: failed             ! the model refused or failed on this parameter set
    logical(lgt)                       :: unscored           ! no calibration target, so nothing to score
    integer(i4b)                       :: errFinal           ! error code of finalize_summa
    type(summa1_type_dec), allocatable :: summa1_struc(:)    ! top-level SUMMA data structure
    integer(i4b), parameter            :: n=1                ! number of SUMMA data structures
    integer(i4b)                       :: i                  ! looping
    integer(i4b)                       :: iTarget            ! calibration target index
    integer(i4b)                       :: nTarget            ! number of calibration targets
    character(len=4)                   :: runString          ! include run in the output filename
    character(len=4)                   :: rankString         ! include rank in the output filename
    character(len=6)                   :: sampleString       ! include sample index in the output filename
    character(len=:), allocatable      :: outputFileSuffix_orig  ! orig suffix (to restore output suffix)
    real(rkind), allocatable           :: timeSim(:)         ! simulated time
    real(rkind), allocatable           :: flowSim(:)         ! simulated streamflow
    real(rkind), allocatable           :: timeObs(:)         ! observed time
    real(rkind), allocatable           :: timeBounds(:,:)    ! observed time bounds
    real(rkind), allocatable           :: flowObs(:)         ! observed streamflow
    character(len=:), allocatable      :: timeSimUnits       ! simulated time units
    character(len=:), allocatable      :: flowSimUnits       ! simulated flow units
    character(len=:), allocatable      :: timeObsUnits       ! observed time units
    character(len=:), allocatable      :: flowObsUnits       ! observed flow units
    real(rkind), allocatable           :: timeAligned(:)     ! common time vector
    real(rkind), allocatable           :: flowSimAligned(:)  ! flow simulations aligned to the common time period
    real(rkind), allocatable           :: flowObsAligned(:)  ! flow observations aligned to the common time period
    type(sim_series_type), allocatable :: series(:)          ! simulated series the targets asked for
    character(len=64),     allocatable :: seriesName(:)      ! names of their variables
    integer(i4b),          allocatable :: seriesUnit(:)      ! spatial unit each is reduced over
    integer(i8b),          allocatable :: seriesUnitId(:)    ! id of that unit
    integer(i4b),          allocatable :: targetUnit(:)      ! spatial unit of each target
    integer(i4b)                       :: nSeries            ! number of them
    integer(i4b)                       :: iSeries            ! index of one of them
    real(rkind), allocatable           :: valSim(:)          ! the simulated series for the current target
    character(len=:), allocatable      :: valSimUnits        ! its units
    real(rkind), allocatable           :: valAccum(:)        ! it, integrated from a rate
    real(rkind), allocatable           :: valAnom(:)         ! a series as departures from its baseline mean
    real(rkind), allocatable           :: timeVal(:)         ! time coordinate of the simulated series for the current target
    character(len=:), allocatable      :: timeValUnits       ! its units
    character(len=1024)                :: cmessage           ! error message of downwind routine

    err=0
    message='evaluate_objectives/'
    objective=realMissing
    failed=.false.; unscored=.false.
    if(present(trial_failed)) trial_failed=.false.

    nTarget=n_calibration_targets(config)
    if(nTarget > size(objective))then
      message=trim(message)//'the objective vector is too short for the configured calibration targets'
      err=20; return
    endif

    ! allocate top-level SUMMA structure
    allocate(summa1_struc(n),stat=err)
    if(err/=0)then
      message=trim(message)//'problem allocating top-level summa structure'
      return
    endif

    ! populate domain and model-instance parallel contexts
    summa1_struc(n)%domain_parallel=domain_parallel
    summa1_struc(n)%instance_parallel=instance_parallel

    ! define unique output filenames for each rank and sample
    outputFileSuffix_orig=trim(output_fileSuffix)
    if(instance_parallel%size > 1)then
      write(rankString,'(I4.4)') instance_parallel%rank
      output_fileSuffix=trim(output_fileSuffix)//'_rank'//rankString
    endif
    if(sample_id > 0)then
      write(sampleString,'(I6.6)') sample_id
      output_fileSuffix=trim(output_fileSuffix)//'_sample'//sampleString
    endif
    if(summa1_struc(n)%config%n_runs > 1)then
      write(runString,'(I4.4)') summa1_struc(n)%config%run_index
      output_fileSuffix=trim(output_fileSuffix)//'_run'//runString
    endif

    ! the model is released after this block however it is left, so a worker can take the next trial
    trial: block

    ! initialize SUMMA; a parameter set the model refuses is a failed trial
    call initialize_summa(config, summa1_struc(n), param_name, param_value, &
                          err=err, message=cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); failed=.true.; exit trial; endif

    ! an objective function needs observations to compare against; without them this is an
    ! ordinary SUMMA run, so run the model and return rather than failing
    if(n_calibration_targets(summa1_struc(n)%config) == 0)then
      unscored=.true.
      call run_summa(summa1_struc(n), timeSim,flowSim, timeSimUnits,flowSimUnits, err,cmessage, physics_failed=failed)
      if(err/=0) message=trim(message)//trim(cmessage)
      exit trial
    endif

    ! calibration run: send model chatter to stderr so stdout carries only the metric, unless a caller owns iulog
    if(iulog == output_unit) iulog = error_unit

    ! the simulated series the targets ask for, resolved before the run so an unknown variable or
    ! spatial unit is refused here rather than after a simulation has been paid for.  Routed streamflow
    ! at the network outlet is not among them: run_summa returns it for every run.
    nSeries=0
    allocate(seriesName(nTarget),seriesUnit(nTarget),seriesUnitId(nTarget),targetUnit(nTarget),stat=err)
    if(err/=0)then
      message=trim(message)//'problem allocating the simulated-series list'
      exit trial
    endif
    do iTarget=1,nTarget
      associate(calTarget => summa1_struc(n)%config%calib%targets(iTarget))
      targetUnit(iTarget)=spatial_unit_index(calTarget%spatial_unit)
      if(targetUnit(iTarget) == integerMissing)then
        message=trim(message)//'calibration target "'//trim(calTarget%name)//'" asks for spatial unit "'// &
                trim(calTarget%spatial_unit)//'"; use gru, hru or reach'
        err=20; exit trial
      endif
      ! with mizuRoute running, streamflow comes from the network, so a GRU's outflow is named by its reach
      if(is_routed_streamflow(calTarget%variable) .and. targetUnit(iTarget) == ix_unit_gru .and. &
         mizuroute_active .and. summa1_struc(n)%config%use_mizuroute)then
        message=trim(message)//'calibration target "'//trim(calTarget%name)//'" asks for streamflow out of a gru '// &
                'while mizuRoute routes it; name the reach instead'
        err=20; exit trial
      endif
      if(is_routed_streamflow(calTarget%variable) .and. targetUnit(iTarget) == ix_unit_domain) cycle
      if(find_series(seriesName(1:nSeries),seriesUnit(1:nSeries),seriesUnitId(1:nSeries), &
                     calTarget%variable,targetUnit(iTarget),calTarget%spatial_id) == integerMissing)then
        nSeries=nSeries+1
        seriesName(nSeries)  =calTarget%variable
        seriesUnit(nSeries)  =targetUnit(iTarget)
        seriesUnitId(nSeries)=calTarget%spatial_id
      endif
      end associate
    enddo

    call init_simulated_series(seriesName(1:nSeries),seriesUnit(1:nSeries),seriesUnitId(1:nSeries), &
                               numtim,series,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); exit trial; endif

    ! a reach belongs to the river network, so it is resolved against the mizuRoute topology
    call resolve_reach_series(summa1_struc(n),series,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); exit trial; endif

    ! run SUMMA once; every target scores this same simulation
    if(nSeries > 0)then
      call run_summa(summa1_struc(n), timeSim,flowSim, timeSimUnits,flowSimUnits, err,cmessage, series=series, physics_failed=failed)
    else
      call run_summa(summa1_struc(n), timeSim,flowSim, timeSimUnits,flowSimUnits, err,cmessage, physics_failed=failed)
    endif
    if(err/=0)then; message=trim(message)//trim(cmessage); exit trial; endif

    ! score the simulation against each calibration target
    do iTarget=1,nTarget
      associate(calTarget => summa1_struc(n)%config%calib%targets(iTarget))

        ! the simulated series this target is compared against: routed streamflow at the network
        ! outlet, or the series collected for its variable over its spatial unit
        if(allocated(valSim)) deallocate(valSim)
        timeVal=timeSim
        timeValUnits=timeSimUnits
        if(is_routed_streamflow(calTarget%variable) .and. targetUnit(iTarget) == ix_unit_domain)then
          valSim=flowSim
          valSimUnits=flowSimUnits
        else
          iSeries=find_simulated_series(series,calTarget%variable,targetUnit(iTarget),calTarget%spatial_id)
          if(iSeries == integerMissing)then
            message=trim(message)//'calibration target "'//trim(calTarget%name)// &
                    '" asks for simulated variable "'//trim(calTarget%variable)//'", which was not collected'
            err=20; exit trial
          endif
          valSim=series(iSeries)%values
          valSimUnits=trim(series(iSeries)%units)
        endif

        ! read this target's observations
        if(allocated(timeObs)) deallocate(timeObs)
        if(allocated(flowObs)) deallocate(flowObs)
        if(allocated(timeBounds)) deallocate(timeBounds)
        call read_observations(trim(calTarget%obs_path),trim(calTarget%obs_file),trim(calTarget%vname_obs), &
                               timeObs,timeBounds,flowObs,timeObsUnits,flowObsUnits,err,cmessage)
        if(err/=0)then
          message=trim(message)//'calibration target "'//trim(calTarget%name)//'": '//trim(cmessage)
          exit trial
        endif

        ! units the target states override the file's
        if(len_trim(calTarget%obs_units) > 0) flowObsUnits=trim(calTarget%obs_units)

        ! -------------------------------------------------------------------------------------------
        ! Put the two series on the same footing, as far as this target asks for
        ! -------------------------------------------------------------------------------------------

        ! integrate a simulated rate into the quantity the observations report
        if(calTarget%accumulate)then
          call accumulate_series(timeVal,valSim,timeValUnits,valAccum,err,cmessage)
          if(err/=0)then
            message=trim(message)//'calibration target "'//trim(calTarget%name)//'": '//trim(cmessage)
            exit trial
          endif
          call move_alloc(valAccum,valSim)
          ! integrating kg m-2 s-1 over seconds leaves kg m-2, which is millimetres of water
          if(trim(valSimUnits)=='kg m-2 s-1') valSimUnits='mm'
        endif

        ! express both sides as departures from the same baseline, which is what an anomaly product
        ! reports and what removes the constant of integration an accumulated series carries
        if(allocated(calTarget%baseline_start) .and. allocated(calTarget%baseline_end))then
          call remove_baseline_mean(timeVal,valSim,timeValUnits,                              &
                                    calTarget%baseline_start,calTarget%baseline_end,          &
                                    valAnom,err,cmessage)
          if(err/=0)then
            message=trim(message)//'calibration target "'//trim(calTarget%name)//'" (simulated): '//trim(cmessage)
            exit trial
          endif
          call move_alloc(valAnom,valSim)

          call remove_baseline_mean(timeObs,flowObs,timeObsUnits,                             &
                                    calTarget%baseline_start,calTarget%baseline_end,          &
                                    valAnom,err,cmessage)
          if(err/=0)then
            message=trim(message)//'calibration target "'//trim(calTarget%name)//'" (observed): '//trim(cmessage)
            exit trial
          endif
          call move_alloc(valAnom,flowObs)
        endif

        ! a balance is the simulated change between seasonal extremes, one value at each observed time
        if(len_trim(calTarget%balance) > 0)then
          call balance_between_extremes(timeVal,valSim,timeValUnits,timeObs,timeObsUnits,            &
                                        calTarget%balance,calTarget%balance_window,valAnom,err,cmessage)
          if(err/=0)then
            message=trim(message)//'calibration target "'//trim(calTarget%name)//'": '//trim(cmessage)
            exit trial
          endif
          call move_alloc(valAnom,valSim)
          timeVal=timeObs
          timeValUnits=timeObsUnits
        endif

        ! align simulated and observed series
        if(allocated(timeAligned))    deallocate(timeAligned)
        if(allocated(flowSimAligned)) deallocate(flowSimAligned)
        if(allocated(flowObsAligned)) deallocate(flowObsAligned)
        call align_timeseries(timeVal,valSim,timeValUnits,valSimUnits,   &
                              timeObs,flowObs,timeObsUnits,flowObsUnits, &
                              summa1_struc(n)%config%calib%start_date,   &
                              summa1_struc(n)%config%calib%end_date,     &
                              timeBounds,                                &
                              timeAligned,flowSimAligned,flowObsAligned, &
                              err,cmessage)
        if(err/=0)then
          message=trim(message)//'calibration target "'//trim(calTarget%name)//'": '//trim(cmessage)
          exit trial
        endif

      ! compute this target's metric
      call compute_metric(flowObsAligned,flowSimAligned,        &
                          calTarget%metric,calTarget%obs_transform,   &
                          objective(iTarget),err,cmessage)
      if(err/=0)then
        message=trim(message)//'calibration target "'//trim(calTarget%name)//'": '//trim(cmessage)
        ! observations in the scored period that this simulation gives nothing to pair with: the parameter set failed
        failed = any(ieee_is_finite(flowObsAligned)) .and. .not.any(ieee_is_finite(flowSimAligned))
        exit trial
      endif

      ! write the aligned evaluation time series and objective value.  The SUMMA output file holds one
      ! such series, so it holds the first target's; the trials file carries every target's value.
      if(summa1_struc(n)%config%write_timeseries .and. iTarget == 1)then
        call write_evaluation(ncid(iLookFREQ%timestep),                    &
                              summa1_struc(n)%config%calib%write_aligned,  &
                              timeAligned,                                 &
                              flowObsAligned,flowSimAligned,               &
                              timeObsUnits,flowObsUnits,                   &
                              calTarget%metric,calTarget%obs_transform,          &
                              objective(iTarget),                          &
                              err,cmessage)
        if(err/=0)then; message=trim(message)//trim(cmessage); exit trial; endif
      endif

      end associate
    enddo

    end block trial
    if(present(trial_failed)) trial_failed=failed

    ! finalize SUMMA and release model resources, keeping the first error
    call finalize_summa(summa1_struc(n),errFinal,cmessage)
    if(err==0 .and. errFinal/=0)then; err=errFinal; message=trim(message)//trim(cmessage); endif

    ! release top-level SUMMA data structure
    ! NOTE: Deallocate here because finalize_summa operates on a single array element
    if(allocated(summa1_struc)) deallocate(summa1_struc)

    ! restore output file suffix
    output_fileSuffix=outputFileSuffix_orig
    if(err/=0 .or. unscored) return

    ! write objective function to standard output
    if(instance_parallel%size == 1)then
      write(output_unit,'(ES24.16)') objective(1)
    else
      write(output_unit,'(A,A,A,I0,A,*(1X,F12.9))') &
           'case=',trim(config%case_name),', rank=',instance_parallel%rank,', objective=', &
           (objective(iTarget), iTarget=1,nTarget)
    endif

  end subroutine evaluate_objectives

  ! ---- PRIVATE SUBROUTINES --------------------------------------------------------------------------

  ! **************************************************************************************************
  ! initialize SUMMA
  ! **************************************************************************************************
  subroutine initialize_summa(config, summa_struct, param_name, param_value, &
                              mpar_out, bpar_out, err, message)
    type(config_info)                    , intent(inout)    :: config
    type(summa1_type_dec)                , intent(inout)    :: summa_struct
    character(*)                         , intent(in)       :: param_name(:)
    real(rkind)                          , intent(in)       :: param_value(:)
    type(gru_hru_dom_doublevec), optional, intent(out)      :: mpar_out
    type(gru_double)           , optional, intent(out)      :: bpar_out
    integer(i4b)                         , intent(out)      :: err
    character(*)                         , intent(out)      :: message
    character(len=1024) :: cmessage

    err = 0
    message = 'initialize_summa/'

    ! declare and allocate SUMMA data structures and initialize model state
    call summa_initialize(config, summa_struct, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! initialize parameter data structures
    call summa_paramSetup(summa_struct, param_name, param_value, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! read restart data and reset model state
    call summa_readRestart(summa_struct, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! initialize OpenWQ
    if(openwq_active)then
      call openwq_init(summa_struct, err)
      if(err/=0)then; message=trim(message)//'problem initializing OpenWQ'; return; endif
    endif

    ! Retain initialized parameter structures
    if(present(mpar_out)) mpar_out = summa_struct%mparStruct
    if(present(bpar_out)) bpar_out = summa_struct%bparStruct

  end subroutine initialize_summa

  ! **************************************************************************************************
  ! Resolve the reach of every series scored on one, against the mizuRoute river network.
  !
  ! A reach id the network does not hold is refused here, before the run, as an unknown SUMMA variable
  ! or GRU is; a reach asked for in a run that does not route is refused the same way.
  ! **************************************************************************************************
  subroutine resolve_reach_series(summa_struct, series, err, message)
    use simulated_series, only: sim_series_type
    use simulated_series, only: ix_unit_reach
    type(summa1_type_dec), intent(in)    :: summa_struct
    type(sim_series_type), intent(inout) :: series(:)
    integer(i4b),          intent(out)   :: err
    character(*),          intent(out)   :: message
    integer(i4b)                         :: iSeries
    character(len=256)                   :: cmessage

    err=0
    message='resolve_reach_series/'

    do iSeries=1,size(series)
      if(series(iSeries)%ix_unit /= ix_unit_reach) cycle
      if(.not.(mizuroute_active .and. summa_struct%config%use_mizuroute))then
        message=trim(message)//'a calibration target names a reach, which needs mizuRoute routing'
        err=20; return
      endif
      if(mizuroute_active)then ! build-time capability
        call get_mizuroute_reach_index(summa_struct, series(iSeries)%unit_id, series(iSeries)%ix_seg, err, cmessage)
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      endif
    enddo

  end subroutine resolve_reach_series

  ! **************************************************************************************************
  ! run SUMMA
  ! **************************************************************************************************
  subroutine run_summa(summa_struct, timeSim,flowSim, timeUnits,flowUnits, err,message, series, physics_failed)
    USE var_lookup, only: iLookFORCE
    USE globalData, only: forc_meta
    USE globalData, only: numtim
    USE simulated_series, only: sim_series_type
    USE simulated_series, only: collect_simulated_series
    USE simulated_series, only: ix_unit_reach
    ! dummy arguments
    type(summa1_type_dec), intent(inout)       :: summa_struct  ! top-level SUMMA data structure
    real(rkind), allocatable, intent(out)      :: timeSim(:)    ! simulation time
    real(rkind), allocatable, intent(out)      :: flowSim(:)    ! simulated streamflow
    character(len=:), allocatable, intent(out) :: timeUnits     ! units and reference time for simulation time
    character(len=:), allocatable, intent(out) :: flowUnits     ! units for simulated streamflow
    integer(i4b), intent(out)                  :: err           ! error code
    character(*), intent(out)                  :: message       ! error message
    ! SUMMA variables a calibration target asked for, recorded as the run proceeds
    type(sim_series_type), optional, intent(inout) :: series(:)
    ! .true. if the error came from the model physics rather than from reading or writing
    logical(lgt), optional, intent(out)            :: physics_failed
    ! locals
    integer(i4b)                               :: modelTimeStep ! index of model time step
    integer(i4b)                               :: iSeries       ! index of a simulated series
    character(len=1024)                        :: cmessage      ! error message of downwind routine

    err=0
    message='run_summa/'
    if(present(physics_failed)) physics_failed=.false.

    ! define units for time and flow
    timeUnits = trim(forc_meta(iLookFORCE%time)%varunit) ! time since reference (varies)
    flowUnits = 'm3/s' ! always in mizuroute

    ! routed streamflow time series
    allocate(timeSim(numtim), source=realMissing)
    allocate(flowSim(numtim), source=realMissing)

    ! loop through time
    do modelTimeStep=1,numtim

      ! read model forcing data
      call summa_readForcing(modelTimeStep, summa_struct, err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); exit; endif

      ! initialize OpenWQ time step
      if(openwq_active) call openwq_run_time_start(summa_struct)

      ! run SUMMA physics and mizuRoute
      call summa_runPhysics(modelTimeStep, summa_struct, err, cmessage)
      if(err/=0)then
        message=trim(message)//trim(cmessage)
        if(present(physics_failed)) physics_failed=.true.
        exit
      endif

      ! transfer SUMMA fluxes to OpenWQ
      if(openwq_active) call openwq_run_space_step(summa_struct)

      ! save streamflow time series at the outlet: routed by mizuRoute, or else SUMMA's routed runoff
      call get_streamflow(modelTimeStep, summa_struct, &
                          timeSim(modelTimeStep),flowSim(modelTimeStep), &
                          err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); exit; endif

      ! flow at any reach a target named (unavailable when mizuRoute is not active)
#ifdef MIZUROUTE_ACTIVE
      if(mizuroute_active)then ! build-time capability
       if(summa_struct%config%use_mizuroute)then
        if(present(series))then
          do iSeries=1,size(series)
            if(series(iSeries)%ix_unit /= ix_unit_reach) cycle
            call get_mizuroute_streamflow(modelTimeStep, summa_struct, series(iSeries)%values(modelTimeStep), &
                                          series(iSeries)%ix_seg)
          enddo
        endif
       endif
      endif
#endif

      ! record the SUMMA variables any calibration target asked for
      if(present(series))then
        call collect_simulated_series(modelTimeStep, summa_struct, series, err, cmessage)
        if(err/=0)then; message=trim(message)//trim(cmessage); exit; endif
      endif

      ! write the model output
      call summa_writeOutputFiles(modelTimeStep, summa_struct, err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); exit; endif

      ! finalize OpenWQ time step
      if(openwq_active) call openwq_run_time_end(summa_struct)
    enddo

  end subroutine run_summa

  ! **************************************************************************************************
  ! finalize SUMMA
  ! **************************************************************************************************
  subroutine finalize_summa(summa_struct, err, message)
    ! SUMMA global data
    use globalData, only: forcNcid                ! netcdf id for current netcdf forcing file
    use globalData, only: ncid                    ! vector of file ids of netcdf output files
    ! SUMMA buffered output structures
    use globalData, only: fullIndxSave
    use globalData, only: fullForcSave
    use globalData, only: fullProgSave
    use globalData, only: fullDiagSave
    use globalData, only: fullFluxSave
    use globalData, only: fullBvarSave
    use netcdf_util_module, only: nc_file_close   ! module to handle netcdf stuff for inputs and outputs
    type(summa1_type_dec), intent(inout) :: summa_struct
    integer(i4b),          intent(out)   :: err
    character(*),          intent(out)   :: message
    integer(i4b)                         :: iFreq
    character(len=1024)                  :: cmessage

    err = 0
    message = 'finalize_summa/'

    ! deallocate SUMMA buffered output structures
    if(allocated(fullIndxSave)) deallocate(fullIndxSave)
    if(allocated(fullForcSave)) deallocate(fullForcSave)
    if(allocated(fullProgSave)) deallocate(fullProgSave)
    if(allocated(fullDiagSave)) deallocate(fullDiagSave)
    if(allocated(fullFluxSave)) deallocate(fullFluxSave)
    if(allocated(fullBvarSave)) deallocate(fullBvarSave)

    ! close NetCDF forcing file
    if(forcNcid/=integerMissing)then
      call nc_file_close(forcNcid, err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      forcNcid = integerMissing
    endif

    ! close SUMMA NetCDF output files
    do iFreq=1,size(ncid)
      if(ncid(iFreq)/=integerMissing)then
        call nc_file_close(ncid(iFreq), err, cmessage)
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
        ncid(iFreq) = integerMissing
      endif
    enddo

    ! report the water check of the OpenWQ coupling
    if(openwq_active) call openwq_finalize()

    ! deallocate mizuroute structures
    if(mizuroute_active)then ! build-time capability
     if(summa_struct%config%use_mizuroute) then
      call finalize_mizuroute(err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
     endif
    endif

    ! more cleanup operations can be added here as required

    ! Allow output libraries to complete file closure
    call sleep(2)

  end subroutine finalize_summa

  ! ---------------------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------------------
  ! ---- HELPER SUBROUTINES ---------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------------------

  subroutine get_streamflow(modelTimeStep,summa_struct,timeSim,flowSim,err,message)

  USE var_lookup, only: iLookFORCE
  USE var_lookup, only: iLookBVAR

  implicit none

  integer(i4b),          intent(in)    :: modelTimeStep
  type(summa1_type_dec), intent(inout) :: summa_struct
  real(rkind),           intent(out)   :: timeSim
  real(rkind),           intent(out)   :: flowSim
  integer(i4b),          intent(out)   :: err
  character(*),          intent(out)   :: message

  logical(lgt) :: get_routed
  integer(i4b) :: iGRU

  err=0
  message='get_streamflow/'

  ! if streamflow is routed using mizuroute
  get_routed = mizuroute_active .and. summa_struct%config%use_mizuroute

  ! simulation time
  timeSim = summa_struct%forcStruct%gru(1)%hru(1)%var(iLookFORCE%time)

  ! get streamflow from mizuroute
  if(get_routed)then
#ifdef MIZUROUTE_ACTIVE
    call get_mizuroute_streamflow(modelTimeStep,summa_struct,flowSim)
#endif

  ! get streamflow from SUMMA routed runoff
  else

    flowSim = 0._rkind

    do iGRU=1,size(summa_struct%bvarStruct%gru)
      flowSim = flowSim + &
        summa_struct%bvarStruct%gru(iGRU)%var(iLookBVAR%averageRoutedRunoff)%dat(1) * &
        summa_struct%bvarStruct%gru(iGRU)%var(iLookBVAR%basin__totalArea)%dat(1)
    enddo

  endif

end subroutine get_streamflow


end module summa_simulation

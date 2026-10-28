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

    call run_summa_streamflow(summa1_struc(n), timeSim,flowSim, timeUnits,flowUnits, err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    call finalize_summa(summa1_struc(n),err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  end subroutine run_simulation


  ! **************************************************************************************************
  ! Run SUMMA and return streamflow at the basin outlet.
  ! **************************************************************************************************
  subroutine run_summa_streamflow(summa_struct, timeSim,flowSim, timeUnits,flowUnits, err,message)

    use globalData, only: numtim
    use simulated_series, only: sim_series_type, init_simulated_series, ix_unit_domain

    implicit none

    type(summa1_type_dec), intent(inout) :: summa_struct
    real(rkind), allocatable, intent(out) :: timeSim(:),flowSim(:)
    character(len=:), allocatable, intent(out) :: timeUnits,flowUnits
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message

    type(sim_series_type), allocatable :: series(:)
    logical(lgt) :: physics_failed
    character(len=1024) :: cmessage

    err=0
    message='run_summa_streamflow/'

    ! Request streamflow at the basin outlet
    call init_simulated_series([character(len=64) :: 'streamflow'], &
                               [ix_unit_domain], [0_i8b], numtim, &
                               series, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! Run SUMMA
    call run_summa(summa_struct, timeSim, timeUnits, &
                   series, physics_failed, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! Return streamflow
    flowSim=series(1)%values
    flowUnits=trim(series(1)%units)

  end subroutine run_summa_streamflow


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
  ! Identify and initialize the distinct simulated series requested by calibration targets.
  ! Multiple targets may use the same simulated series.
  ! **************************************************************************************************
  subroutine prepare_target_series(config, targetUnit, series, err, message)

    use globalData, only: numtim
    use simulated_series, only: sim_series_type, init_simulated_series
    use simulated_series, only: spatial_unit_index, is_routed_streamflow
    use simulated_series, only: ix_unit_gru

    type(config_info), intent(in) :: config
    integer(i4b), allocatable, intent(out) :: targetUnit(:)
    type(sim_series_type), allocatable, intent(out) :: series(:)
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message

    character(len=64), allocatable :: seriesName(:)
    integer(i4b), allocatable :: seriesUnit(:)
    integer(i8b), allocatable :: seriesUnitId(:)
    integer(i4b) :: nTarget, nSeries, iTarget
    character(len=1024) :: cmessage

    err=0
    message='prepare_target_series/'

    nTarget=size(config%calib%targets)
    nSeries=0

    allocate(seriesName(nTarget),seriesUnit(nTarget),seriesUnitId(nTarget),targetUnit(nTarget),stat=err)
    if(err/=0)then
      message=trim(message)//'problem allocating the simulated-series list'
      return
    endif

    do iTarget=1,nTarget
      associate(calTarget => config%calib%targets(iTarget))

        ! Identify the spatial unit requested by this target
        targetUnit(iTarget)=spatial_unit_index(calTarget%spatial_unit)

        if(targetUnit(iTarget) == integerMissing)then
          message=trim(message)//'calibration target "'//trim(calTarget%name)//'" asks for spatial unit "'// &
                  trim(calTarget%spatial_unit)//'"; use domain, gru, hru or reach'
          err=20
          return
        endif

        ! GRU streamflow is not available when mizuRoute is active
        if(is_routed_streamflow(calTarget%variable) .and. targetUnit(iTarget) == ix_unit_gru .and. &
           mizuroute_active .and. config%use_mizuroute)then
          message=trim(message)//'calibration target "'//trim(calTarget%name)//'" asks for streamflow out of a gru '// &
                  'while mizuRoute routes it; name the reach instead'
          err=20
          return
        endif

        ! Collect each distinct simulated series only once
        if(find_series(seriesName(1:nSeries),seriesUnit(1:nSeries),seriesUnitId(1:nSeries), &
                       calTarget%variable,targetUnit(iTarget),calTarget%spatial_id) == integerMissing)then

          nSeries=nSeries+1
          seriesName(nSeries)  =calTarget%variable
          seriesUnit(nSeries)  =targetUnit(iTarget)
          seriesUnitId(nSeries)=calTarget%spatial_id

        endif

      end associate
    enddo

    ! Resolve variables and spatial units, then allocate the time series
    call init_simulated_series(seriesName(1:nSeries),seriesUnit(1:nSeries),seriesUnitId(1:nSeries), &
                               numtim,series,err,cmessage)
    if(err/=0) message=trim(message)//trim(cmessage)

  end subroutine prepare_target_series

  ! **************************************************************************************************
  ! Score one calibration target against the simulated time series.
  ! **************************************************************************************************
  subroutine score_target(config, iTarget, targetUnit, series, timeSim, timeSimUnits, &
                          objective, failed, err, message)

    use globalData, only: ncid
    use var_lookup, only: iLookFREQ
    use simulated_series, only: sim_series_type, find_simulated_series
    use read_flowobs_module, only: read_observations
    use timeseries_alignment, only: align_timeseries
    use series_transform, only: accumulate_series, remove_baseline_mean, balance_between_extremes
    use metrics, only: compute_metric
    use write_evaluation_module, only: write_evaluation
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    type(config_info), intent(in) :: config
    integer(i4b), intent(in) :: iTarget, targetUnit
    type(sim_series_type), intent(in) :: series(:)
    real(rkind), intent(in) :: timeSim(:)
    character(*), intent(in) :: timeSimUnits
    real(rkind), intent(out) :: objective
    logical(lgt), intent(inout) :: failed
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message

    integer(i4b) :: iSeries
    real(rkind), allocatable :: timeObs(:),timeBounds(:,:),flowObs(:)
    real(rkind), allocatable :: timeVal(:),valSim(:),valAccum(:),valAnom(:)
    real(rkind), allocatable :: timeAligned(:),flowSimAligned(:),flowObsAligned(:)
    character(len=:), allocatable :: timeObsUnits,flowObsUnits
    character(len=:), allocatable :: timeValUnits,valSimUnits
    character(len=1024) :: cmessage

    err=0
    message=''
    objective=realMissing

    associate(calTarget => config%calib%targets(iTarget))

      ! Get the simulated series requested by this target
      iSeries=find_simulated_series(series,calTarget%variable,targetUnit,calTarget%spatial_id)
      if(iSeries == integerMissing)then
        message='calibration target "'//trim(calTarget%name)// &
                '" asks for simulated variable "'//trim(calTarget%variable)//'", which was not collected'
        err=20; return
      endif

      timeVal=timeSim
      timeValUnits=timeSimUnits
      valSim=series(iSeries)%values
      valSimUnits=trim(series(iSeries)%units)

      ! Read observations
      call read_observations(trim(calTarget%obs_path),trim(calTarget%obs_file),trim(calTarget%vname_obs), &
                             timeObs,timeBounds,flowObs,timeObsUnits,flowObsUnits,err,cmessage)
      if(err/=0)then
        message='calibration target "'//trim(calTarget%name)//'": '//trim(cmessage)
        return
      endif

      if(len_trim(calTarget%obs_units) > 0) flowObsUnits=trim(calTarget%obs_units)

      ! Accumulate a simulated rate when requested
      if(calTarget%accumulate)then
        call accumulate_series(timeVal,valSim,timeValUnits,valAccum,err,cmessage)
        if(err/=0)then
          message='calibration target "'//trim(calTarget%name)//'": '//trim(cmessage)
          return
        endif
        call move_alloc(valAccum,valSim)
        if(trim(valSimUnits)=='kg m-2 s-1') valSimUnits='mm'
      endif

      ! Remove the baseline mean from simulated and observed series
      if(allocated(calTarget%baseline_start) .and. allocated(calTarget%baseline_end))then
        call remove_baseline_mean(timeVal,valSim,timeValUnits,                     &
                                  calTarget%baseline_start,calTarget%baseline_end, &
                                  valAnom,err,cmessage)
        if(err/=0)then
          message='calibration target "'//trim(calTarget%name)//'" (simulated): '//trim(cmessage)
          return
        endif
        call move_alloc(valAnom,valSim)

        call remove_baseline_mean(timeObs,flowObs,timeObsUnits,                    &
                                  calTarget%baseline_start,calTarget%baseline_end, &
                                  valAnom,err,cmessage)
        if(err/=0)then
          message='calibration target "'//trim(calTarget%name)//'" (observed): '//trim(cmessage)
          return
        endif
        call move_alloc(valAnom,flowObs)
      endif

      ! Calculate the change between seasonal extremes
      if(len_trim(calTarget%balance) > 0)then
        call balance_between_extremes(timeVal,valSim,timeValUnits,timeObs,timeObsUnits, &
                                      calTarget%balance,calTarget%balance_window,      &
                                      valAnom,err,cmessage)
        if(err/=0)then
          message='calibration target "'//trim(calTarget%name)//'": '//trim(cmessage)
          return
        endif
        call move_alloc(valAnom,valSim)
        timeVal=timeObs
        timeValUnits=timeObsUnits
      endif

      ! Align simulated and observed time series
      call align_timeseries(timeVal,valSim,timeValUnits,valSimUnits,   &
                            timeObs,flowObs,timeObsUnits,flowObsUnits, &
                            config%calib%start_date,                  &
                            config%calib%end_date,                    &
                            timeBounds,                               &
                            timeAligned,flowSimAligned,flowObsAligned, &
                            err,cmessage)
      if(err/=0)then
        message='calibration target "'//trim(calTarget%name)//'": '//trim(cmessage)
        return
      endif

      ! Compute the calibration metric
      call compute_metric(flowObsAligned,flowSimAligned,             &
                          calTarget%metric,calTarget%obs_transform, &
                          objective,err,cmessage)
      if(err/=0)then
        message='calibration target "'//trim(calTarget%name)//'": '//trim(cmessage)
        failed=any(ieee_is_finite(flowObsAligned)) .and. .not.any(ieee_is_finite(flowSimAligned))
        return
      endif

      ! Write evaluation time series for the first target
      if(config%write_timeseries .and. iTarget == 1)then
        call write_evaluation(ncid(iLookFREQ%timestep),   &
                              config%calib%write_aligned, &
                              timeAligned,                &
                              flowObsAligned,flowSimAligned, &
                              timeObsUnits,flowObsUnits,  &
                              calTarget%metric,calTarget%obs_transform, &
                              objective,err,cmessage)
        if(err/=0)then
          message=trim(cmessage)
          return
        endif
      endif

    end associate

  end subroutine score_target

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
    use globalData, only: output_fileSuffix
    use simulated_series, only: sim_series_type
    ! dummy arguments
    type(config_info),           intent(inout) :: config
    type(parallel_context_type), intent(in)    :: domain_parallel
    type(parallel_context_type), intent(in)    :: instance_parallel
    integer(i4b)               , intent(in)    :: sample_id
    character(*)               , intent(in)    :: param_name(:)
    real(rkind)                , intent(in)    :: param_value(:)
    real(rkind), allocatable   , intent(out)   :: objective(:)
    integer(i4b)               , intent(out)   :: err
    character(*)               , intent(out)   :: message
    logical(lgt), optional, intent(out) :: trial_failed
    ! locals
    logical(lgt)                       :: failed                 ! true if the parameter set causes a model or scoring failure
    logical(lgt)                       :: unscored               ! true if no calibration targets are configured
    integer(i4b)                       :: errFinal               ! error code from model finalization
    type(summa1_type_dec), allocatable :: summa1_struc(:)        ! top-level SUMMA data structures
    integer(i4b), parameter            :: n=1                    ! index of the SUMMA instance
    integer(i4b)                       :: iTarget                ! calibration target index
    integer(i4b)                       :: nTarget                ! number of calibration targets
    character(len=4)                   :: runString              ! run index for output filenames
    character(len=4)                   :: rankString             ! MPI rank for output filenames
    character(len=6)                   :: sampleString           ! sample index for output filenames
    character(len=:), allocatable      :: outputFileSuffix_orig  ! original output filename suffix
    real(rkind), allocatable           :: timeSim(:)             ! simulated time series
    character(len=:), allocatable      :: timeSimUnits           ! simulation time units and reference date
    type(sim_series_type), allocatable :: series(:)              ! simulated variables requested by calibration targets
    integer(i4b), allocatable          :: targetUnit(:)          ! spatial unit index for each calibration target
    character(len=1024)                :: cmessage               ! error message from called routines


    err=0
    message='evaluate_objectives/'
    failed=.false.; unscored=.false.
    if(present(trial_failed)) trial_failed=.false.
    
    ! Allocate objective values
    nTarget = 0
    if(allocated(config%calib%targets)) nTarget = size(config%calib%targets)
    
    allocate(objective(max(1,nTarget)),stat=err)
    if(err/=0)then
      message=trim(message)//'problem allocating objective vector'
      return
    endif
    objective=realMissing

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
    if(config%n_runs > 1)then
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

      allocate(series(0))

      call run_summa(summa1_struc(n), timeSim, timeSimUnits, &
                     series, failed, err, cmessage)
      if(err/=0) message=trim(message)//trim(cmessage)

      exit trial
    endif

    ! calibration run: send model chatter to stderr so stdout carries only the metric, unless a caller owns iulog
    if(iulog == output_unit) iulog = error_unit

    ! Prepare the simulated time series required by calibration
    call prepare_target_series(summa1_struc(n)%config, targetUnit, series, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); exit trial; endif

    ! a reach belongs to the river network, so it is resolved against the mizuRoute topology
    call resolve_reach_series(summa1_struc(n),series,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); exit trial; endif

    ! Run SUMMA once and collect every requested simulated series
    call run_summa(summa1_struc(n), timeSim, timeSimUnits, &
                   series, failed, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); exit trial; endif

    ! Score the simulation against each calibration target
    do iTarget=1,nTarget
      call score_target(summa1_struc(n)%config, iTarget, targetUnit(iTarget), &
                        series, timeSim, timeSimUnits, objective(iTarget), &
                        failed, err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); exit trial; endif
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
  subroutine run_summa(summa_struct, timeSim, timeUnits, series, physics_failed, err, message)
    USE var_lookup, only: iLookFORCE
    USE globalData, only: forc_meta
    USE globalData, only: numtim
    USE simulated_series, only: sim_series_type
    USE simulated_series, only: collect_simulated_series
    ! dummy arguments
    type(summa1_type_dec), intent(inout)       :: summa_struct   ! top-level SUMMA data structure
    real(rkind), allocatable, intent(out)      :: timeSim(:)     ! simulation time
    character(len=:), allocatable, intent(out) :: timeUnits      ! units and reference time for simulation time
    type(sim_series_type), intent(inout)       :: series(:)      ! simulated time series
    logical(lgt), intent(out)                  :: physics_failed ! flag if physics fails
    integer(i4b), intent(out)                  :: err            ! error code
    character(*), intent(out)                  :: message        ! error message
    ! locals
    integer(i4b)                               :: modelTimeStep ! index of model time step
    character(len=1024)                        :: cmessage      ! error message of downwind routine

    err=0
    message='run_summa/'
    physics_failed=.false.

    ! Define simulation time
    timeUnits=trim(forc_meta(iLookFORCE%time)%varunit)
    allocate(timeSim(numtim),source=realMissing)

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
        physics_failed=.true.
        exit
      endif

      ! transfer SUMMA fluxes to OpenWQ
      if(openwq_active) call openwq_run_space_step(summa_struct)

      ! Save simulation time
      timeSim(modelTimeStep) = summa_struct%forcStruct%gru(1)%hru(1)%var(iLookFORCE%time)

      ! Collect all requested simulated variables, including streamflow
      call collect_simulated_series(modelTimeStep, summa_struct, series, err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); exit; endif

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

end module summa_simulation

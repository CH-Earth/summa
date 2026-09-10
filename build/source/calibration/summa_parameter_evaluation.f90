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
! SUMMA parameter evaluation
!
! Provides initialization and model-evaluation routines used for parameter
! estimation. MPI scheduling of parameter trials is handled separately.
! **************************************************************************************************

module summa_parameter_evaluation

  ! data types
  USE nr_type,    only: i4b,rkind
  USE summa_type, only: config_info,parallel_context_type

  ! parameter-search data types
  USE parameter_search, only: parameter_spec,parameter_search_info

  implicit none
  private

  public :: spinup_from_cold
  public :: initialize_parameter_evaluation
  public :: evaluate_parameter_sample

contains


  ! **************************************************************************************************
  ! Spin up SUMMA from a cold state.
  !
  ! Performs a one-year cold-start spinup.
  ! Additional parameter-specific warmup will be required after this point.
  !
  ! This initial warmup also provides access to the parameter data structures.
  ! **************************************************************************************************

  subroutine spinup_from_cold(config,             & ! SUMMA configuration structure
                              domain_parallel,    & ! MPI context for domain parallelism
                              instance_parallel,  & ! MPI context for model-instance parallelism
                              err,message)          ! error code and error message

    ! SUMMA global configuration
    USE globalData, only: ixRestart
    USE globalData, only: ixRestart_end,ixRestart_never
    USE globalData, only: output_fileSuffix
    
    ! filenames
    USE summaFileManager, only: SIM_START_TM,SIM_END_TM
    
    ! SUMMA simulation
    USE summa_simulation, only: run_simulation

    implicit none

    ! dummy variables
    type(config_info),           intent(inout) :: config             ! SUMMA configuration structure
    type(parallel_context_type), intent(in)    :: domain_parallel    ! MPI context for domain parallelism
    type(parallel_context_type), intent(in)    :: instance_parallel  ! MPI context for model-instance parallelism
    integer(i4b),                intent(out)   :: err                ! error code
    character(*),                intent(out)   :: message            ! error message

    ! strings
    character(len=4)   :: rankString
    character(len=256) :: cmessage

    character(len=:), allocatable :: outputFileSuffix_orig
    character(len=:), allocatable :: simStartOriginal,simEndOriginal,spinStart
    character(len=:), allocatable :: timeUnits,flowUnits

    ! baseline parameter vectors
    character(len=64), allocatable :: param_name(:)
    real(rkind),       allocatable :: param_value(:)

    ! simulated time and flow
    real(rkind), allocatable :: timeSim(:),flowSim(:)

    ! local variables
    integer(i4b) :: iyear

    err=0
    message='spinup_from_cold/'

    ! zero-length parameter vectors for baseline initialization
    allocate(param_name(0),param_value(0),stat=err)
    if(err/=0)then
      message=trim(message)//'unable to allocate baseline parameter vectors'
      return
    endif

    ! Modify dates to define a one-year cold-start spinup
    ! prior to the start of the simulation period

    simStartOriginal=trim(SIM_START_TM)
    simEndOriginal=trim(SIM_END_TM)

    spinStart=simStartOriginal

    read(spinStart(1:4),*) iyear
    write(spinStart(1:4),'(I4.4)') iyear-1

    SIM_START_TM=spinStart
    SIM_END_TM=simStartOriginal

    ! use rank-specific output filenames during the common spinup

    write(rankString,'(I4.4)') instance_parallel%rank

    outputFileSuffix_orig=trim(output_fileSuffix)
    output_fileSuffix=trim(outputFileSuffix_orig)//'_spinup_rank'//rankString

    ! write the common restart state on rank 0 only
    ixRestart=merge(ixRestart_end, ixRestart_never, instance_parallel%rank == 0)

    ! run SUMMA for one year following the cold start
    call run_simulation(config,                 & ! SUMMA configuration structure
                        domain_parallel,        & ! MPI context for domain parallelism
                        instance_parallel,      & ! MPI context for model-instance parallelism
                        timeSim,flowSim,        & ! simulated time and streamflow
                        timeUnits,flowUnits,    & ! time and streamflow units
                        param_name,param_value, & ! parameter names and values
                        err,cmessage)             ! error code and message
    
    ! restore original global simulation settings
    ! NOTE: Do this before processing the error code
    SIM_START_TM=simStartOriginal
    SIM_END_TM=simEndOriginal
    output_fileSuffix=outputFileSuffix_orig
    ixRestart=ixRestart_never

    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  end subroutine spinup_from_cold


  ! **************************************************************************************************
  ! Initialize parameter evaluation.
  !
  ! Defines the parameter-search information on all ranks, allocates the sampled parameter vector,
  ! initializes the random-number generator on rank 0, and creates rank-specific calibration output
  ! files on worker ranks.
  ! **************************************************************************************************

  subroutine initialize_parameter_evaluation(config,                 &
                                             instance_parallel,      &
                                             param_spec,search,      &
                                             param_value,ncid_calib, &
                                             err,message)

    ! logging
    USE globalData, only: iulog

    ! file paths
    USE summaFileManager, only: OUTPUT_PATH

    ! parameter search
    USE parameter_search, only: initialize_parameter_search

    ! SUMMA parameter information
    USE summa_parameter_spec, only: get_summa_parameter_spec

    ! calibration output
    USE calibration_output_module, only: create_calibration_output

    implicit none

    ! dummy variables
    type(config_info),           intent(in)    :: config             ! SUMMA configuration structure
    type(parallel_context_type), intent(in)    :: instance_parallel  ! MPI context for model-instance parallelism
    type(parameter_spec),        intent(out)   :: param_spec         ! SUMMA parameter specifications
    type(parameter_search_info), intent(out)   :: search             ! parameter-search information
    real(rkind), allocatable,    intent(out)   :: param_value(:)     ! parameter values for the current trial
    integer(i4b),                intent(out)   :: ncid_calib         ! calibration output NetCDF ID
    integer(i4b),                intent(out)   :: err                ! error code
    character(*),                intent(out)   :: message            ! error message

    ! random-number generator
    integer(i4b)              :: nseed
    integer(i4b), allocatable :: seed(:)

    ! calibration output
    character(len=4)   :: rankString
    character(len=256) :: calib_file

    ! local variables
    integer(i4b)       :: i
    character(len=256) :: cmessage

    err=0
    message='initialize_parameter_evaluation/'

    ! -----------------------------------------------------------------------------------------------
    ! Initialize parameter search on all ranks
    ! -----------------------------------------------------------------------------------------------

    ! build the SUMMA parameter specification from the calibration configuration:
    ! identify sampled and constraint-only parameters, retrieve SUMMA defaults and bounds, apply
    ! parameter transformations, and map parameter names in ordered constraints to indices
    call get_summa_parameter_spec(config,param_spec,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! construct the parameter-search information and validate the parameter specification,
    ! including the sampled parameter names and bounds and the ordered parameter constraints
    call initialize_parameter_search(param_spec,search,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! -----------------------------------------------------------------------------------------------
    ! Allocate sampled parameter vector on all ranks
    ! -----------------------------------------------------------------------------------------------

    allocate(param_value(size(search%param_names)),stat=err)

    if(err/=0)then
      message=trim(message)//'unable to allocate sampled parameter vector'
      return
    endif

    ! -----------------------------------------------------------------------------------------------
    ! Initialize random-number generator (rank 0) 
    ! -----------------------------------------------------------------------------------------------

    if(instance_parallel%rank == 0)then

      ! print parameter bounds once
      write(iulog,'(/,A)') 'Calibration parameter bounds:'
      write(iulog,'(A)')   '  Parameter                         Lower              Upper'

      do i=1,size(search%param_names)
        write(iulog,'(2X,A30,2X,ES16.8,2X,ES16.8)') &
          trim(search%param_names(i)),search%lower(i),search%upper(i)
      enddo

      ! initialize random-number generator on the dispatcher
      call random_seed(size=nseed)

      allocate(seed(nseed),stat=err)

      if(err/=0)then
        message=trim(message)//'unable to allocate random-number seed'
        return
      endif

      seed=42
      call random_seed(put=seed)

      ! dispatcher does not write a calibration output file
      ncid_calib=-1

    endif


    ! -----------------------------------------------------------------------------------------------
    ! Create calibration output file (workers)
    ! -----------------------------------------------------------------------------------------------

    if(instance_parallel%rank /= 0)then

      ! define a rank-specific calibration output file
      write(rankString,'(I4.4)') instance_parallel%rank

      calib_file=trim(OUTPUT_PATH)//trim(config%case_name)// &
                 '_calibration_rank'//rankString//'.nc'

      call create_calibration_output(calib_file,                  &
                                     param_spec,                  &
                                     instance_parallel%rank,      &
                                     config%case_name,            &
                                     config%calib%metric,         &
                                     config%calib%obs_transform,  &
                                     ncid_calib,                  &
                                     err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    endif

  end subroutine initialize_parameter_evaluation


  ! **************************************************************************************************
  ! Evaluate one parameter sample.
  !
  ! Constructs the complete SUMMA parameter override vector, evaluates the objective function,
  ! writes the calibration result, and records the parameter trial in the rank-specific log.
  ! **************************************************************************************************

  subroutine evaluate_parameter_sample(config,                 &
                                       domain_parallel,        &
                                       instance_parallel,      &
                                       param_spec,search,      &
                                       sample_id,local_sample, &
                                       param_value,ncid_calib, &
                                       objective,              &
                                       err,message)

    ! SUMMA-specific parameter information and overrides
    USE summa_parameter_spec, only: build_summa_parameter_overrides

    ! objective-function evaluation
    USE summa_simulation, only: evaluate_objective

    ! output filename suffix and logging
    USE globalData, only: output_fileSuffix,iulog

    ! calibration output
    USE calibration_output_module, only: write_calibration_output

    implicit none

    ! dummy variables
    type(config_info),           intent(inout) :: config             ! SUMMA configuration structure
    type(parallel_context_type), intent(in)    :: domain_parallel    ! MPI context for domain parallelism
    type(parallel_context_type), intent(in)    :: instance_parallel  ! MPI context for model-instance parallelism

    type(parameter_spec),        intent(in)    :: param_spec         ! SUMMA parameter specification
    type(parameter_search_info), intent(in)    :: search             ! parameter-search information

    integer(i4b),                intent(in)    :: sample_id          ! global parameter trial index
    integer(i4b),                intent(in)    :: local_sample       ! local calibration output index
    real(rkind),                 intent(in)    :: param_value(:)     ! sampled parameter values
    integer(i4b),                intent(in)    :: ncid_calib         ! calibration output NetCDF ID

    real(rkind),                 intent(out)   :: objective          ! objective function value
    integer(i4b),                intent(out)   :: err                ! error code
    character(*),                intent(out)   :: message            ! error message

    ! complete SUMMA parameter overrides
    character(len=64), allocatable :: param_name(:)
    real(rkind),       allocatable :: param_override(:)

    ! strings to create unique output filenames
    character(len=4)              :: rankString
    character(len=6)              :: sampleString
    character(len=:), allocatable :: outputFileSuffix_orig

    ! timing
    integer(i4b) :: startModelRun(8),endModelRun(8)

    ! local variables
    character(len=256) :: cmessage

    err=0
    message='evaluate_parameter_sample/'

    ! record start time for this parameter trial
    call date_and_time(values=startModelRun)

    ! construct the complete spatially uniform SUMMA parameter override vector, using
    ! sampled values where available and defaults for required non-sampled parameters
    call build_summa_parameter_overrides(param_spec,         &
                                         search%param_names, &
                                         param_value,        &
                                         param_name,         &
                                         param_override,     &
                                         err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! define a unique SUMMA output suffix for this parameter trial
    write(rankString,'(I4.4)') instance_parallel%rank
    write(sampleString,'(I6.6)') sample_id

    outputFileSuffix_orig=trim(output_fileSuffix)
    output_fileSuffix=trim(outputFileSuffix_orig)// &
                      '_rank'//rankString//'_sample'//sampleString

    ! evaluate SUMMA for the complete parameter override vector
    call evaluate_objective(config,                 & ! SUMMA configuration structure
                            domain_parallel,        & ! MPI context for domain parallelism
                            instance_parallel,      & ! MPI context for model-instance parallelism
                            param_name,             & ! parameter names
                            param_override,         & ! parameter values
                            objective,              & ! objective function value
                            err,cmessage)             ! error code and message

    ! restore the output suffix before returning to the caller
    output_fileSuffix=outputFileSuffix_orig

    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! record end time for this parameter trial
    call date_and_time(values=endModelRun)

    ! write parameter values, objective value, and timing information
    call write_calibration_output(ncid_calib,      &
                                  local_sample,    &
                                  param_name,      &
                                  param_override,  &
                                  objective,       &
                                  startModelRun,   &
                                  endModelRun,     &
                                  err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! print global sample number, sampled parameters, and objective
    write(iulog,'(I6,*(1X,ES16.8))') sample_id,param_value,objective

  end subroutine evaluate_parameter_sample


end module summa_parameter_evaluation

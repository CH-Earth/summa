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

USE nr_type, only: i4b, rkind
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
USE finalize_mizuroute_module, only: finalize_mizuroute
#endif

#ifdef OPENWQ_ACTIVE
USE summa_openwq, only: openwq_init
USE summa_openwq, only: openwq_run_time_start
USE summa_openwq, only: openwq_run_space_step
USE summa_openwq, only: openwq_run_time_end
#endif

! module-level data structure to share configurations

implicit none
private

public :: run_simulation
public :: evaluate_objective

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
                            err, message)             ! error code and message
  
    ! dummy arguments
  
    type(config_info),           intent(inout) :: config
    type(parallel_context_type), intent(in)    :: domain_parallel
    type(parallel_context_type), intent(in)    :: instance_parallel
  
    real(rkind), allocatable, intent(out) :: timeSim(:)
    real(rkind), allocatable, intent(out) :: flowSim(:)
  
    character(len=:), allocatable, intent(out) :: timeUnits
    character(len=:), allocatable, intent(out) :: flowUnits
  
    character(*), intent(in) :: param_name(:)
    real(rkind),  intent(in) :: param_value(:)
  
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
  
    ! locals
  
    type(summa1_type_dec), allocatable :: summa1_struc(:)
    integer(i4b), parameter            :: n=1
    character(len=256)                 :: cmessage
  
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
  
    call initialize_summa(config,                 &
                          summa1_struc(n),        &
                          param_name,param_value, &
                          err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  
    call run_summa(summa1_struc(n),      &
                   timeSim,flowSim,      &
                   timeUnits,flowUnits,  &
                   err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  
    call finalize_summa(summa1_struc(n),err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  
  end subroutine run_simulation
  
  ! **************************************************************************************************
  ! Evaluate the objective function for a specified parameter vector.
  ! The routine initializes SUMMA, reads the observed streamflow time series,
  ! runs the model, computes the objective function, and finalizes the simulation.
  ! **************************************************************************************************
 
  subroutine evaluate_objective(config,                 & ! SUMMA configuration structure
                                domain_parallel,        & ! MPI context for domain parallelism
                                instance_parallel,      & ! MPI context for model-instance parallelism
                                param_name,param_value, & ! parameter names and values
                                metric,                 & ! objective function value
                                err, message)             ! error code and message
 
    use iso_fortran_env, only: output_unit

    use globalData, only: ncid
    use var_lookup, only: iLookFREQ
    
    use read_flowobs_module,     only: read_flow_observations
    use timeseries_alignment,    only: align_timeseries
    use metrics,                 only: compute_metric

    use write_evaluation_module, only: write_evaluation

    ! dummy arguments
    
    type(config_info),           intent(inout) :: config
    type(parallel_context_type), intent(in)    :: domain_parallel
    type(parallel_context_type), intent(in)    :: instance_parallel
   
    character(*), intent(in)  :: param_name(:)
    real(rkind),  intent(in)  :: param_value(:)
   
    real(rkind),  intent(out) :: metric
   
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message

    ! locals

    type(summa1_type_dec), allocatable :: summa1_struc(:)    ! top-level SUMMA data structure
    integer(i4b), parameter            :: n=1                ! number of SUMMA data structures
  
    integer(i4b)                       :: i                  ! looping

    real(rkind), allocatable           :: timeSim(:)         ! simulated time
    real(rkind), allocatable           :: flowSim(:)         ! simulated streamflow
    real(rkind), allocatable           :: timeObs(:)         ! observed time
    real(rkind), allocatable           :: flowObs(:)         ! observed streamflow
  
    character(len=:), allocatable      :: timeSimUnits       ! simulated time units
    character(len=:), allocatable      :: flowSimUnits       ! simulated flow units
    character(len=:), allocatable      :: timeObsUnits       ! observed time units
    character(len=:), allocatable      :: flowObsUnits       ! observed flow units
  
    real(rkind), allocatable           :: timeAligned(:)     ! common time vector
    real(rkind), allocatable           :: flowSimAligned(:)  ! flow simulations aligned to the common time period 
    real(rkind), allocatable           :: flowObsAligned(:)  ! flow observations aligned to the common time period

    character(len=256)                 :: cmessage           ! error message of downwind routine
  
    err=0
    message='evaluate_objective/'
 
    ! allocate top-level SUMMA structure
    allocate(summa1_struc(n),stat=err)
    if(err/=0)then
      message=trim(message)//'problem allocating top-level summa structure'
      return
    endif
  
    ! populate domain and model-instance parallel contexts
    summa1_struc(n)%domain_parallel=domain_parallel
    summa1_struc(n)%instance_parallel=instance_parallel
    
    ! initialize SUMMA
    call initialize_summa(config,                &
                          summa1_struc(n),       &
                          param_name,param_value,&
                          err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! read observed streamflow
    call read_flow_observations(summa1_struc(n),              &
                                timeObs,flowObs,              &
                                timeObsUnits,flowObsUnits,    &
                                err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    
    ! run SUMMA
    call run_summa(summa1_struc(n),           &
                   timeSim,flowSim,           &
                   timeSimUnits,flowSimUnits, &
                   err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! align simulated and observed streamflow
    call align_timeseries(timeSim,flowSim,timeSimUnits,flowSimUnits, &
                          timeObs,flowObs,timeObsUnits,flowObsUnits, &
                          summa1_struc(n)%config%calib%start_date,   &
                          summa1_struc(n)%config%calib%end_date,     &
                          timeAligned,flowSimAligned,flowObsAligned, &
                          err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   
    ! compute objective function
    call compute_metric(flowObsAligned,flowSimAligned,                 &
                        summa1_struc(n)%config%calib%metric,           &
                        summa1_struc(n)%config%calib%obs_transform,    &
                        metric,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! write aligned evaluation time series and objective value
    call write_evaluation(ncid(iLookFREQ%timestep),                    &
                          summa1_struc(n)%config%calib%write_aligned,  &
                          timeAligned,                                 &
                          flowObsAligned,flowSimAligned,               &
                          timeObsUnits,flowObsUnits,                   &
                          summa1_struc(n)%config%calib%metric,         &
                          summa1_struc(n)%config%calib%obs_transform,  &
                          metric,                                      &
                          err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    call finalize_summa(summa1_struc(n),err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! release top-level SUMMA data structure
    ! NOTE: Deallocate here because finalize_summa operates on a single array element
    if(allocated(summa1_struc)) deallocate(summa1_struc)

    ! write objective function to standard output
    if(instance_parallel%size == 1)then
      write(output_unit,'(ES24.16)') metric
    else
      write(output_unit,'(A,A,A,I0,A,F12.9)') &
           'case=',trim(config%case_name),', rank=',instance_parallel%rank,', objective=',metric
    endif

  end subroutine evaluate_objective

  ! ---------------------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------------------
  ! ---- PRIVATE SUBROUTINES --------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------------------

  ! **************************************************************************************************
  ! initialize SUMMA
  ! **************************************************************************************************
  subroutine initialize_summa(config, summa_struct, param_name, param_value, err, message)

    type(config_info),       intent(inout)    :: config
    type(summa1_type_dec)  , intent(inout)    :: summa_struct
    character(*)           , intent(in)       :: param_name(:)
    real(rkind)            , intent(in)       :: param_value(:)
    integer(i4b)           , intent(out)      :: err
    character(*)           , intent(out)      :: message

    character(len=256) :: cmessage

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
      call openwq_init(err)
      if(err/=0)then; message=trim(message)//'problem initializing OpenWQ'; return; endif
    endif

  end subroutine initialize_summa


  ! **************************************************************************************************
  ! run SUMMA
  ! **************************************************************************************************
  subroutine run_summa(summa_struct,      &
                       timeSim,flowSim,   &
                       timeUnits,flowUnits, &
                       err,message)

   
    USE var_lookup, only: iLookFORCE
    USE globalData, only: forc_meta
    USE globalData, only: numtim
    
    ! dummy arguments
    type(summa1_type_dec), intent(inout)       :: summa_struct  ! top-level SUMMA data structure
    real(rkind), allocatable, intent(out)      :: timeSim(:)    ! simulation time
    real(rkind), allocatable, intent(out)      :: flowSim(:)    ! simulated streamflow
    character(len=:), allocatable, intent(out) :: timeUnits     ! units and reference time for simulation time
    character(len=:), allocatable, intent(out) :: flowUnits     ! units for simulated streamflow
    integer(i4b), intent(out)                  :: err           ! error code
    character(*), intent(out)                  :: message       ! error message
   
    ! locals
    integer(i4b)                               :: modelTimeStep ! index of model time step
    character(len=256)                         :: cmessage      ! error message of downwind routine
   
    err=0
    message='run_summa/'

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
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! initialize OpenWQ time step
      if(openwq_active) call openwq_run_time_start(summa_struct)

      ! run SUMMA physics and mizuRoute
      call summa_runPhysics(modelTimeStep, summa_struct, err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! transfer SUMMA fluxes to OpenWQ
      if(openwq_active) call openwq_run_space_step(summa_struct)

      ! save streamflow time series (unavailable when mizuRoute is not active)
      if(mizuroute_active)then ! build-time capability
       if(summa_struct%config%use_mizuroute)then
        timeSim(modelTimeStep) = summa_struct%forcStruct%gru(1)%hru(1)%var(iLookFORCE%time)
        call get_mizuroute_streamflow(modelTimeStep, summa_struct, flowSim(modelTimeStep))
       endif
      endif

      ! write the model output
      call summa_writeOutputFiles(modelTimeStep, summa_struct, err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

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
    character(len=256)                   :: cmessage

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

end module summa_simulation

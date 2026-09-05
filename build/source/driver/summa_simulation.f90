module summa_simulation

USE nr_type, only: i4b, rkind
USE summa_type, only: summa1_type_dec

USE summa_init, only: summa_initialize
USE summa_setup, only: summa_paramSetup
USE summa_restart, only: summa_readRestart
USE summa_forcing, only: summa_readForcing
USE summa_modelRun, only: summa_runPhysics
USE summa_writeOutput, only: summa_writeOutputFiles

USE globalData, only: realMissing
USE globalData, only: iulog

USE build_options, only: mizuroute_active
USE build_options, only: openwq_active

USE read_flowobs_module,  only: read_flow_observations
USE timeseries_alignment, only: align_timeseries
USE metrics,              only: compute_metric

#ifdef MIZUROUTE_ACTIVE
USE mizuroute_coupling,  only: get_mizuroute_streamflow
#endif

#ifdef OPENWQ_ACTIVE
USE summa_openwq, only: openwq_init
USE summa_openwq, only: openwq_run_time_start
USE summa_openwq, only: openwq_run_space_step
USE summa_openwq, only: openwq_run_time_end
#endif

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

  subroutine run_simulation(comm, rank, nproc,         &
                            timeSim, flowSim,         &
                            timeUnits, flowUnits,     &
                            param_name, param_value,  &
                            err, message)
    ! dummy arguments
   
    integer(i4b), intent(in)                   :: comm            ! MPI communicator
    integer(i4b), intent(in)                   :: rank            ! MPI rank
    integer(i4b), intent(in)                   :: nproc           ! number of MPI processes
   
    real(rkind), allocatable, intent(out)      :: timeSim(:)      ! simulation time
    real(rkind), allocatable, intent(out)      :: flowSim(:)      ! simulated streamflow
   
    character(len=:), allocatable, intent(out) :: timeUnits       ! units and reference time for simulation time
    character(len=:), allocatable, intent(out) :: flowUnits       ! units for simulated streamflow
   
    character(*), intent(in)                   :: param_name(:)   ! parameter names
    real(rkind),  intent(in)                   :: param_value(:)  ! parameter values
   
    integer(i4b), intent(out)                  :: err             ! error code
    character(*), intent(out)                  :: message         ! error message
   
    ! locals
   
    type(summa1_type_dec), allocatable         :: summa1_struc(:)  ! master SUMMA data structure
    integer(i4b), parameter                    :: n=1              ! n copies of the SUMMA data structure
    character(len=256)                         :: cmessage         ! error message of downwind routine
   
    err=0
    message='run_simulation/'
   
    allocate(summa1_struc(n), stat=err)
    if (err/=0) then
      message=trim(message)//'problem allocating master summa structure'
      return
    endif
   
    summa1_struc(n)%parallel%comm = comm
    summa1_struc(n)%parallel%rank = rank
    summa1_struc(n)%parallel%size = nproc
   
    call initialize_summa(summa1_struc(n), param_name, param_value, err, cmessage)
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
  
  subroutine evaluate_objective(comm, rank, nproc,        &
                                param_name, param_value,  &
                                metric,                   &
                                err, message)
 
    use iso_fortran_env, only: error_unit
    use iso_fortran_env, only: output_unit

    ! dummy arguments
  
    integer(i4b), intent(in)           :: comm               ! MPI communicator
    integer(i4b), intent(in)           :: rank               ! MPI rank
    integer(i4b), intent(in)           :: nproc              ! number of MPI processes
  
    character(*), intent(in)           :: param_name(:)      ! parameter names
    real(rkind),  intent(in)           :: param_value(:)     ! parameter values
  
    real(rkind), intent(out)           :: metric             ! objective-function value
  
    integer(i4b), intent(out)          :: err                ! error code
    character(*), intent(out)          :: message            ! error message
  
    ! locals
  
    type(summa1_type_dec), allocatable :: summa1_struc(:)    ! master SUMMA data structure
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
 
    ! send log information to stderr
    iulog = error_unit

    ! allocate master SUMMA structure
    allocate(summa1_struc(n),stat=err)
    if(err/=0)then
      message=trim(message)//'problem allocating master summa structure'
      return
    endif
  
    ! define parallel context
    summa1_struc(n)%parallel%comm = comm
    summa1_struc(n)%parallel%rank = rank
    summa1_struc(n)%parallel%size = nproc
  
    ! initialize SUMMA
    call initialize_summa(summa1_struc(n),param_name,param_value,err,cmessage)
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
                          summa1_struc(n)%obj%start_date,            &
                          summa1_struc(n)%obj%end_date,              &
                          timeAligned,flowSimAligned,flowObsAligned, &
                          err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   
    ! compute objective function
    call compute_metric(flowObsAligned,flowSimAligned,           &
                        summa1_struc(n)%obj%metric,              &
                        summa1_struc(n)%obj%transformation,      &
                        metric,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    call finalize_summa(summa1_struc(n),err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! write error metric to standard output
    write(output_unit,'(ES24.16)') metric

  end subroutine evaluate_objective

  ! ---------------------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------------------
  ! ---- PRIVATE SUBROUTINES --------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------------------
  ! ---------------------------------------------------------------------------------------------------

  ! **************************************************************************************************
  ! initialize SUMMA
  ! **************************************************************************************************
  subroutine initialize_summa(summa_struct, param_name, param_value, err, message)

    type(summa1_type_dec)  , intent(inout)    :: summa_struct
    character(*)           , intent(in)       :: param_name(:)
    real(rkind)            , intent(in)       :: param_value(:)
    integer(i4b)           , intent(out)      :: err
    character(*)           , intent(out)      :: message

    character(len=256) :: cmessage

    err = 0
    message = 'initialize_summa/'

    ! declare and allocate SUMMA data structures and initialize model state
    call summa_initialize(summa_struct, err, cmessage)
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
    type(summa1_type_dec), intent(inout)       :: summa_struct  ! master SUMMA data structure
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
      if(mizuroute_active)then
        timeSim(modelTimeStep) = summa_struct%forcStruct%gru(1)%hru(1)%var(iLookFORCE%time)
        call get_mizuroute_streamflow(modelTimeStep, summa_struct, flowSim(modelTimeStep))
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

    type(summa1_type_dec), intent(inout) :: summa_struct
    integer(i4b),          intent(out)   :: err
    character(*),          intent(out)   :: message

    err = 0
    message = 'finalize_summa/'

    ! cleanup operations can be added here as required
 
    ! Allow output libraries to complete file closure
    call sleep(2)

  end subroutine finalize_summa

end module summa_simulation

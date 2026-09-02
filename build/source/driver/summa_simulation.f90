module summa_simulation

USE nr_type, only: i4b, rkind
USE summa_type, only: summa1_type_dec

USE summa_init, only: summa_initialize
USE summa_setup, only: summa_paramSetup
USE summa_restart, only: summa_readRestart
USE summa_forcing, only: summa_readForcing
USE summa_modelRun, only: summa_runPhysics
USE summa_writeOutput, only: summa_writeOutputFiles

USE globalData, only: numtim
USE globalData, only: realMissing
USE globalData, only: print_step_freq
USE globalData, only: isPrint

USE build_options, only: mizuroute_active

#ifdef MIZUROUTE_ACTIVE
USE mizuroute_coupling, only: get_mizuroute_streamflow
#endif

implicit none
private

public :: initialize_simulation
public :: run_simulation
public :: finalize_simulation

contains

  ! **************************************************************************************************
  ! initialize SUMMA
  ! **************************************************************************************************
  subroutine initialize_simulation(summa_struct, err, message)

    type(summa1_type_dec), intent(inout) :: summa_struct
    integer(i4b),          intent(out)   :: err
    character(*),          intent(out)   :: message

    character(len=256) :: cmessage

    err = 0
    message = 'initialize_simulation/'

    ! declare and allocate SUMMA data structures and initialize model state
    call summa_initialize(summa_struct, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! initialize parameter data structures
    call summa_paramSetup(summa_struct, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! read restart data and reset model state
    call summa_readRestart(summa_struct, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  end subroutine initialize_simulation


  ! **************************************************************************************************
  ! run SUMMA
  ! **************************************************************************************************
  subroutine run_simulation(summa_struct, simTime, simFlow, err, message)

    USE var_lookup, only: iLookFORCE

    type(summa1_type_dec),    intent(inout)      :: summa_struct
    real(rkind), allocatable, intent(out)        :: simTime(:)
    real(rkind), allocatable, intent(out)        :: simFlow(:)
    integer(i4b),             intent(out)        :: err
    character(*),             intent(out)        :: message

    integer(i4b)       :: modelTimeStep
    character(len=256) :: cmessage

    err = 0
    message = 'run_simulation/'

    ! routed streamflow time series
    allocate(simTime(numtim), source=realMissing)
    allocate(simFlow(numtim), source=realMissing)

    ! loop through time
    do modelTimeStep=1,numtim

      ! read model forcing data
      call summa_readForcing(modelTimeStep, summa_struct, err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

      if(isPrint .and. mod(modelTimeStep,print_step_freq)==0)then
        print *, 'step ---> ', modelTimeStep
      endif

      ! run SUMMA physics and mizuRoute
      call summa_runPhysics(modelTimeStep, summa_struct, err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! save streamflow time series (unavailable when mizuRoute is not active)
      if(mizuroute_active)then
        simTime(modelTimeStep) = summa_struct%forcStruct%gru(1)%hru(1)%var(iLookFORCE%time)
        call get_mizuroute_streamflow(modelTimeStep, summa_struct, simFlow(modelTimeStep))
      endif

      ! write the model output
      call summa_writeOutputFiles(modelTimeStep, summa_struct, err, message)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    enddo

  end subroutine run_simulation


  ! **************************************************************************************************
  ! finalize SUMMA
  ! **************************************************************************************************
  subroutine finalize_simulation(summa_struct, err, message)

    type(summa1_type_dec), intent(inout) :: summa_struct
    integer(i4b),          intent(out)   :: err
    character(*),          intent(out)   :: message

    err = 0
    message = 'finalize_simulation/'

    ! cleanup operations can be added here as required
 
    ! Allow output libraries to complete file closure
    call sleep(2)

  end subroutine finalize_simulation

end module summa_simulation

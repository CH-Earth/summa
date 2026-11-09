! *************************************************************************************************
! dynamic_case_dispatch.f90
! Utilities for dynamically assigning calibration cases to model-instance groups.
! *************************************************************************************************

module dynamic_case_dispatch

  USE mpi
  USE nr_type,   only: i4b

  USE parallel_context, only: parallel_context_type
  USE error_utils, only: check_mpi,abort_mpi

  implicit none
  private

  public :: init_case_counter
  public :: get_next_case
  public :: finalize_case_counter

  public :: dispatch_calibration_cases

  ! Interface for executing an assigned calibration case
  abstract interface
    subroutine execute_case_interface(context,case_index,err,message)
      import :: i4b
      class(*),         intent(inout) :: context
      integer(i4b),     intent(in)    :: case_index
      integer(i4b),     intent(out)   :: err
      character(len=*), intent(out)   :: message
    end subroutine execute_case_interface
  end interface

contains


  ! **************************************************************************************************
  ! Dispatch calibration cases dynamically across independent MPI instance groups.
  !
  ! Each instance-group leader claims cases from a shared MPI counter and broadcasts the
  ! assignment to its workers. A model-specific callback executes each assigned case.
  ! All MPI ranks must call this routine collectively.
  ! **************************************************************************************************
  subroutine dispatch_calibration_cases(nCases,world_parallel,instance_parallel, &
                                        execute_case,context,err,message)

    ! ---------------------------------------------------------------------------------------
    ! Dummy arguments
    ! ---------------------------------------------------------------------------------------
    integer(i4b),                intent(in)    :: nCases             ! Total number of calibration cases
    type(parallel_context_type), intent(in)    :: world_parallel     ! Global MPI context
    type(parallel_context_type), intent(in)    :: instance_parallel  ! Calibration instance MPI context
    procedure(execute_case_interface)          :: execute_case       ! Model-specific case execution
    class(*),                    intent(inout) :: context            ! Calibration context
    integer(i4b),                intent(out)   :: err                ! Error code
    character(len=*),            intent(out)   :: message            ! Error message

    ! ---------------------------------------------------------------------------------------
    ! Local variables
    ! ---------------------------------------------------------------------------------------
    integer(i4b) :: next_case   ! Global index of the next unassigned case
    integer(i4b) :: case_index  ! Index of the assigned calibration case
    integer(i4b) :: iDispatch   ! Number of case-dispatch attempts
    integer      :: case_win    ! MPI window exposing the global case counter
    integer      :: mpi_err     ! MPI error code

    err=0
    message='dispatch_calibration_cases/'

    ! initialize the shared case counter
    call init_case_counter(world_parallel%comm,world_parallel%rank, &
                           next_case,case_win,err,message)
    if(err/=0) call abort_mpi(world_parallel%rank,trim(message))

    ! dynamically assign cases to independent instance groups
    do iDispatch=1,nCases

      ! instance leader claims the next available case
      if(instance_parallel%rank==0)then
        call get_next_case(case_win,nCases,case_index,err,message)
        if(err/=0) call abort_mpi(world_parallel%rank,trim(message))
      endif

      ! broadcast the assigned case to all instance-group ranks
      call MPI_Bcast(case_index,1,MPI_INTEGER,0,instance_parallel%comm,mpi_err)
      call check_mpi(world_parallel%rank,mpi_err,'unable to broadcast case index')

      ! zero indicates that all cases have been assigned
      if(case_index==0) exit

      ! execute the assigned case using the model-specific callback
      call execute_case(context,case_index,err,message)
      if(err/=0) call abort_mpi(instance_parallel%rank,trim(message))

    enddo

    ! all ranks collectively release the global case counter
    call finalize_case_counter(case_win,err,message)
    if(err/=0) call abort_mpi(world_parallel%rank,trim(message))

  end subroutine dispatch_calibration_cases

  ! ===============================================================================================
  ! Initialize the global case counter
  ! ===============================================================================================
  subroutine init_case_counter(world_comm,world_rank,next_case,case_win,err,message)
    ! dummy
    integer(i4b),                intent(in)    :: world_comm  ! MPI communicator containing all ranks
    integer(i4b),                intent(in)    :: world_rank  ! Rank in the global MPI communicator
    integer(i4b),                intent(out)   :: next_case   ! Next case index stored in the global counter
    integer,                     intent(out)   :: case_win    ! MPI window exposing the global case counter
    integer(i4b),                intent(out)   :: err         ! Error code
    character(*),                intent(out)   :: message     ! Error message
    ! local
    integer                                    :: mpi_err    ! MPI error code
    integer                                    :: int_size   ! Size of the case counter in bytes
    integer(MPI_ADDRESS_KIND)                  :: win_size   ! Size of the MPI window in bytes

    err     = 0
    message = 'init_case_counter/'

    ! size of the case counter in bytes
    int_size = storage_size(next_case)/8

    ! world rank 0 owns the global counter
    if(world_rank == 0)then
      next_case = 1_i4b
      win_size  = int(int_size,MPI_ADDRESS_KIND)
    else
      next_case = 0_i4b
      win_size  = 0_MPI_ADDRESS_KIND
    endif

    ! expose the counter through an MPI window
    call MPI_Win_create(next_case, win_size,int_size, MPI_INFO_NULL, &
                        world_comm,case_win,mpi_err)
    if(mpi_err /= MPI_SUCCESS)then
      err = mpi_err
      message = trim(message)//'MPI_Win_create failed'
      return
    endif

    ! open the passive-target access epoch for dynamic case dispatch
    call MPI_Win_lock_all(0,case_win,mpi_err)
    if(mpi_err /= MPI_SUCCESS)then
      err = mpi_err
      message = trim(message)//'MPI_Win_lock_all failed'
      return
    endif

  end subroutine init_case_counter


  ! ===============================================================================================
  ! Atomically claim the next calibration case
  ! NOTE: This routine should only be called by instance_parallel rank 0.
  !       iCase=0 indicates that no cases remain.
  ! ===============================================================================================
  subroutine get_next_case(case_win,nCases,iCase,err,message)
    integer,      intent(in)    :: case_win
    integer(i4b), intent(in)    :: nCases
    integer(i4b), intent(out)   :: iCase
    integer(i4b), intent(out)   :: err
    character(*), intent(out)   :: message
    integer(i4b), parameter     :: one =  1_i4b
    integer                     :: mpi_err

    err     = 0
    message = 'get_next_case/'

    ! Atomically:
    !   iCase    = next_case
    !   next_case = next_case + 1
    ! The counter is stored on world rank 0 at displacement zero.
    call MPI_Fetch_and_op(one, iCase, MPI_INTEGER, 0, 0_MPI_ADDRESS_KIND, &
                          MPI_SUM, case_win, mpi_err)
    if(mpi_err /= MPI_SUCCESS)then
      message = trim(message)//'MPI_Fetch_and_op failed'
      err = mpi_err; return
    endif

    ! ensure the remote update and local fetch are complete
    call MPI_Win_flush(0,case_win,mpi_err)
    if(mpi_err /= MPI_SUCCESS)then
      message = trim(message)//'MPI_Win_flush failed'
      err = mpi_err; return
    endif

    ! zero signals that all cases have already been assigned
    if(iCase > nCases) iCase = 0_i4b

  end subroutine get_next_case


  ! ===============================================================================================
  ! Finalize the global case counter
  ! ===============================================================================================
  subroutine finalize_case_counter(case_win, err, message)
    integer,      intent(inout) :: case_win
    integer(i4b), intent(out)   :: err
    character(*), intent(out)   :: message
    integer                     :: mpi_err

    err     = 0
    message = 'finalize_case_counter/'

    call MPI_Win_unlock_all(case_win,mpi_err)
    if(mpi_err /= MPI_SUCCESS)then
      message = trim(message)//'MPI_Win_unlock_all failed'
      err = mpi_err; return
    endif

    call MPI_Win_free(case_win,mpi_err)
    if(mpi_err /= MPI_SUCCESS)then
      message = trim(message)//'MPI_Win_free failed'
      err = mpi_err; return
    endif

  end subroutine finalize_case_counter

end module dynamic_case_dispatch

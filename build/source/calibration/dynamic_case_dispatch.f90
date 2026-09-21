! *************************************************************************************************
! dynamic_case_dispatch.f90
!
! Utilities for dynamically assigning calibration cases to model-instance groups.
! *************************************************************************************************

module dynamic_case_dispatch

  USE mpi
  USE nr_type,   only: i4b
  USE summa_type,only: parallel_context_type

  implicit none
  private

  public :: init_case_counter
  public :: get_next_case
  public :: finalize_case_counter

contains

  ! ===============================================================================================
  ! Initialize the global case counter
  ! ===============================================================================================
  subroutine init_case_counter(world_parallel, next_case, case_win, err, message)

    type(parallel_context_type), intent(in)    :: world_parallel
    integer(i4b),                intent(out)   :: next_case
    integer,                     intent(out)   :: case_win
    integer(i4b),                intent(out)   :: err
    character(*),                intent(out)   :: message

    integer                                   :: mpi_err
    integer                                   :: int_size
    integer(MPI_ADDRESS_KIND)                 :: win_size

    err     = 0
    message = 'init_case_counter/'

    ! size of the case counter in bytes
    int_size = storage_size(next_case)/8

    ! world rank 0 owns the global counter
    if(world_parallel%rank == 0)then
      next_case = 1_i4b
      win_size  = int(int_size,MPI_ADDRESS_KIND)
    else
      next_case = 0_i4b
      win_size  = 0_MPI_ADDRESS_KIND
    endif

    ! expose the counter through an MPI window
    call MPI_Win_create(next_case,          &
                        win_size,           &
                        int_size,           &
                        MPI_INFO_NULL,       &
                        world_parallel%comm, &
                        case_win,            &
                        mpi_err)

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
  !
  ! NOTE: This routine should only be called by instance_parallel rank 0.
  !       iCase=0 indicates that no cases remain.
  ! ===============================================================================================
  subroutine get_next_case(case_win, nCases, iCase, err, message)

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
    !
    !   iCase    = next_case
    !   next_case = next_case + 1
    !
    ! The counter is stored on world rank 0 at displacement zero.
    call MPI_Fetch_and_op(one,                 &
                          iCase,               &
                          MPI_INTEGER,         &
                          0,                   &
                          0_MPI_ADDRESS_KIND,  &
                          MPI_SUM,             &
                          case_win,            &
                          mpi_err)

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

module mpi_context

  USE nr_type, only: i4b

  use mpi

  implicit none
  private

  public :: set_mpi_context
  public :: init_mpi

contains
  ! start MPI for threads that never call it; threads default to one per rank, not every core per rank
  subroutine init_mpi(ierr)
    !$ use omp_lib, only: omp_set_num_threads
    integer(I4B), intent(out) :: ierr
    integer(I4B)              :: provided, status

    call MPI_Init_thread(MPI_THREAD_FUNNELED, provided, ierr)
    !$ if (ierr == MPI_SUCCESS .and. provided < MPI_THREAD_FUNNELED) ierr = MPI_ERR_OTHER
    !$ call get_environment_variable('OMP_NUM_THREADS', status=status)
    !$ if (status /= 0) call omp_set_num_threads(1)

  end subroutine init_mpi

  subroutine set_mpi_context(comm, rank, size, ierr, message)
    integer(I4B),     intent(in)  :: comm
    integer(I4B),     intent(out) :: rank
    integer(I4B),     intent(out) :: size
    integer(I4B),     intent(out) :: ierr
    character(len=*), intent(out) :: message

    ierr    = MPI_SUCCESS
    message = ''
    call MPI_Comm_rank(comm, rank, ierr)
    if (ierr /= MPI_SUCCESS) then
      message = 'set_mpi_context: MPI_Comm_rank failed'
      return
    end if
    call MPI_Comm_size(comm, size, ierr)
    if (ierr /= MPI_SUCCESS) then
      message = 'set_mpi_context: MPI_Comm_size failed'
      return
    end if

  end subroutine set_mpi_context

end module mpi_context

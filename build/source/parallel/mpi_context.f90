module mpi_context

  USE nr_type, only: i4b

  use mpi

  implicit none
  private

  public :: set_mpi_context

contains

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

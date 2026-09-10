module error_utils

  USE nr_type, only: i4b

  use mpi

  implicit none

  private

  public :: check_mpi
  public :: abort_mpi

contains

  subroutine check_mpi(rank, ierr, message)

    integer(i4b),     intent(in) :: rank
    integer(i4b),     intent(in) :: ierr
    character(len=*), intent(in) :: message

    character(len=MPI_MAX_ERROR_STRING) :: mpi_message
    integer(i4b)  :: message_length
    integer(i4b)  :: ierr_string

    if (ierr /= MPI_SUCCESS) then

       call MPI_Error_string(ierr, mpi_message, message_length, ierr_string)

       if (ierr_string == MPI_SUCCESS) then
          call abort_mpi(rank, trim(message) // ': ' // &
               trim(mpi_message(1:message_length)))
       else
          call abort_mpi(rank, trim(message))
       end if

    end if

  end subroutine check_mpi


  subroutine abort_mpi(rank, message)

    integer(i4b),     intent(in) :: rank
    character(len=*), intent(in) :: message

    integer(i4b) :: ierr

    write(*,'(A,I0,A,A)') &
         'ERROR [rank ', rank, ']: ', trim(message)

    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)

  end subroutine abort_mpi

end module error_utils

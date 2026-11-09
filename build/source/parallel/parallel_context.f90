module parallel_context

  USE nr_type, only: i4b

  implicit none
  private

  ! **************************************************************************************************
  ! MPI communication context.
  !
  ! Stores the communicator, rank, and number of ranks for a parallel execution group.
  ! **************************************************************************************************
  type, public :: parallel_context_type
    integer(i4b) :: comm = -1  ! MPI communicator
    integer(i4b) :: rank = 0   ! Rank within the communicator
    integer(i4b) :: size = 1   ! Number of ranks in the communicator
  end type parallel_context_type

end module parallel_context

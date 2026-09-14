module summa_work_balance

  use nr_type, only : I4B

  implicit none

  private
  public :: balance_even

contains

  subroutine balance_even(first_work, nwork, rank, size, &
                          first_work_local, nwork_local, ierr, message)

    implicit none

    integer(I4B),           intent(in)  :: first_work
    integer(I4B),           intent(in)  :: nwork
    integer(I4B),           intent(in)  :: rank
    integer(I4B),           intent(in)  :: size

    integer(I4B),           intent(out) :: first_work_local
    integer(I4B),           intent(out) :: nwork_local
    integer(I4B),           intent(out) :: ierr
    character(len=*),       intent(out) :: message

    integer(I4B) :: nwork_per_rank
    integer(I4B) :: n_larger_blocks

    ! Initialize outputs
    ierr             = 0
    message          = ''
    first_work_local = 0
    nwork_local      = 0

    ! Check inputs
    if (nwork < 0) then
      ierr = 1
      message = 'balance_even: nwork must be non-negative'
      return
    end if

    if (size <= 0) then
      ierr = 2
      message = 'balance_even: communicator size must be positive'
      return
    end if

    if (rank < 0 .or. rank >= size) then
      ierr = 3
      message = 'balance_even: invalid communicator rank'
      return
    end if

    ! Nothing to distribute
    if (nwork == 0) return

    ! Ceiling of nwork / number of ranks, using integer arithmetic
    nwork_per_rank = (nwork + size - 1) / size

    ! Number of ranks receiving nwork_per_rank work units.
    ! Remaining ranks receive nwork_per_rank - 1 work units.
    n_larger_blocks = nwork - (nwork_per_rank - 1) * size

    if (rank < n_larger_blocks) then

      nwork_local      = nwork_per_rank
      first_work_local = first_work + rank*nwork_per_rank

    else

      nwork_local = nwork_per_rank - 1

      first_work_local = first_work + &
        n_larger_blocks * nwork_per_rank + &
        (rank - n_larger_blocks) * &
        (nwork_per_rank - 1)

    end if

  end subroutine balance_even

end module summa_work_balance

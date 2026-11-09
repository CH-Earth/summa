module test_model

  USE nr_type, only: rkind

  implicit none
  private

  public :: initialize_model
  public :: evaluate_model

contains

  subroutine initialize_model()
  end subroutine initialize_model

  subroutine evaluate_model(theta,simulated)

    real(rkind), intent(in) :: theta(:)
    real(rkind), allocatable, intent(out) :: simulated(:)

    simulated = [theta(1),2._rkind*theta(1),3._rkind*theta(1)]

  end subroutine evaluate_model

end module test_model

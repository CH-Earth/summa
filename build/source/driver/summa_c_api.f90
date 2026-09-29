! ===============================================================
!  summa_c_api.f90  --  C-bindable wrapper around evaluate_objective
! ===============================================================
!  Exposes SUMMA's objective function to C and, through ctypes, to
!  Python:
!
!    summa_evaluate(param_values, n, objective, err)  bind(C)
!
!  Parameter values are paired with PARAM_NAMES below and applied
!  through SUMMA's existing override path (param_override.f90), the
!  same route used by the --param command-line flag. An unrecognised
!  name returns err=20 rather than being silently ignored.
!
!  The objective is KGE (get_kge in objfunc/metrics.f90): higher is
!  better, 1 is perfect.
!
!  Build and test: see calibration/README.md ("make libsumma").
!
!  Limitation: PARAM_NAMES is fixed at compile time, so changing the
!  calibration parameters requires a rebuild. evaluate_objective
!  already accepts a name array, so passing names from Python is the
!  natural next step.
! ===============================================================

module summa_c_api

  use iso_c_binding
  use nr_type,          only: i4b, rkind
  use summa_simulation, only: evaluate_objective

  implicit none

  ! Calibration parameters. Names must match SUMMA's parameter tables
  ! (localParamInfo.txt / basinParamInfo.txt); see calibration/README.md.
  integer, parameter :: N_PARAMS = 15
  character(len=64), parameter :: PARAM_NAMES(N_PARAMS) = &
      [character(len=64) :: "k_soil", "theta_sat", "aquiferBaseflowExp", &
       "aquiferBaseflowRate", "qSurfScale", "summerLAI", "frozenPrecipMultip", &
       "Fcapil", "tempCritRain", "heightCanopyTop", "heightCanopyBottom", &
       "windReductionParam", "vGn_n", "routingGammaScale", "routingGammaShape"]

contains

  ! ─────────────────────────────────────────────────────────────
  ! summa_evaluate: takes a flat array of parameter values (in the
  ! order defined by PARAM_NAMES above), runs SUMMA, and returns the
  ! objective function value.
  ! ─────────────────────────────────────────────────────────────
  subroutine summa_evaluate(param_values, n, objective, err) &
      bind(C, name="summa_evaluate")

    integer(c_int),  intent(in), value :: n
    real(c_double),  intent(in)        :: param_values(n)
    real(c_double),  intent(out)       :: objective
    integer(c_int),  intent(out)       :: err

    ! parallel dummy variables (matches summa_driver.f90 — serial only)
    integer(i4b), parameter :: comm=0, rank=0, nproc=1

    character(len=64), allocatable :: param_name(:)
    real(rkind),       allocatable :: param_value(:)
    real(rkind)                    :: obj_internal
    integer(i4b)                   :: err_internal
    character(len=1024)            :: message

    integer :: i

    if (n /= N_PARAMS) then
      err = 1
      objective = -huge(1.0_c_double)
      return
    end if

    allocate(param_name(n))
    allocate(param_value(n))

    param_name = PARAM_NAMES
    do i = 1, n
      param_value(i) = real(param_values(i), rkind)
    end do

    call evaluate_objective(comm, rank, nproc,       &
                            param_name, param_value, &
                            obj_internal,             &
                            err_internal, message)

    objective = real(obj_internal, c_double)
    err       = int(err_internal, c_int)

    if (err_internal /= 0) then
      ! don't call stop_program here — let Python decide what to do
      ! with a failed evaluation (e.g. return a penalty value)
      print *, "[summa_c_api] evaluate_objective error: ", trim(message)
    end if

    deallocate(param_name, param_value)

  end subroutine summa_evaluate

end module summa_c_api
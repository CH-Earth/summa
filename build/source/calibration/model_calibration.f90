module model_calibration
  use nr_type, only: rkind,i4b,lgt
  use parameter_search, only: parameter_spec,parameter_search_info
  use parameter_fields, only: parameter_field,parameter_field_plan,initialize_parameter_fields
  use parameter_constraints, only: constraint_plan,initialize_parameter_constraints,check_parameter_constraints
  implicit none
  private

  ! One context per case/concurrent checker. After initialization, only effective
  ! field values and the selected scalar values may change, not the specification,
  ! bounds, initial values or spatial/local ordering. Reinitialize after such edits.
  type, public :: model_calibration_context
    type(parameter_spec) :: spec
    type(parameter_search_info) :: search
    type(parameter_field), allocatable :: fields(:)
    real(rkind), allocatable :: values(:)
    type(parameter_field_plan) :: parameterization
    type(constraint_plan) :: constraints
  end type model_calibration_context

  public :: initialize_model_constraints,is_feasible
contains
  subroutine initialize_model_constraints(context,err,message)
    type(model_calibration_context), intent(inout) :: context
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    ! Clear both readiness states even if a later setup operation fails.
    context%parameterization=parameter_field_plan()
    context%constraints=constraint_plan()
    if(.not.allocated(context%fields))then
      err=20; message='calibration fields are not initialized'; return
    endif
    call initialize_parameter_fields(context%spec,context%search,context%fields, &
                                     context%parameterization,err,message)
    if(err/=0) return
    call initialize_parameter_constraints(context%spec,context%fields,context%parameterization, &
                                          context%constraints,err,message)
  end subroutine initialize_model_constraints

  ! Common feasibility callback for random/DDS sampling and NSGA-II offspring.
  ! Ordinary rejection: feasible=.false., err=0. Checking failure: err/=0.
  subroutine is_feasible(context,x,feasible,err,message)
    class(*), intent(inout) :: context
    real(rkind), intent(in) :: x(:)
    logical(lgt), intent(out) :: feasible
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    feasible=.false.
    select type(ctx=>context)
    type is(model_calibration_context)
      if(.not.allocated(ctx%fields))then
        err=20; message='calibration fields are not initialized'; return
      endif
      call check_parameter_constraints(ctx%spec,ctx%search,ctx%fields,ctx%parameterization, &
                                       ctx%constraints,x,feasible,err,message)
    class default
      err=20; message='incorrect feasibility context type'
    end select
  end subroutine is_feasible
end module model_calibration

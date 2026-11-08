module test_support
  use nr_type, only: i4b,rkind,lgt
  use parameter_search, only: parameter_info,parameter_spec,parameter_search_info,initialize_parameter_search
  use parameter_fields, only: parameter_field
  use model_calibration, only: model_calibration_context,initialize_model_constraints
  implicit none
  integer :: assertions=0
contains
  subroutine require(ok,label)
    logical, intent(in) :: ok
    character(*), intent(in) :: label
    assertions=assertions+1
    if(.not.ok)then
      print '(A)', 'FAIL: '//trim(label)
      stop 1
    endif
  end subroutine require

  subroutine set_param(p,name,sampled,lo,hi,kind,reference)
    type(parameter_info), intent(out) :: p
    character(*), intent(in) :: name
    logical(lgt), intent(in) :: sampled
    real(rkind), intent(in) :: lo,hi
    integer(i4b), intent(in) :: kind
    character(*), optional, intent(in) :: reference
    p%name=name
    p%sampled=sampled
    p%lower=lo
    p%upper=hi
    p%trial_value=lo/2._rkind+hi/2._rkind
    p%param_type=kind
    if(present(reference)) p%reference_name=reference
  end subroutine set_param

  subroutine set_field(f,first,second)
    type(parameter_field), intent(out) :: f
    real(rkind), intent(in) :: first(:)
    real(rkind), optional, intent(in) :: second(:)
    integer :: n
    n=1
    if(present(second)) n=2
    allocate(f%initial%hspatial(n))
    f%initial%hspatial(1)%dat=first
    if(present(second)) f%initial%hspatial(2)%dat=second
  end subroutine set_field

  ! The production search initializer now preserves configured bounds and
  ! delegates all ordering checks to the generic feasibility checker.
  subroutine build_search(spec,search,err,message)
    type(parameter_spec), intent(in) :: spec
    type(parameter_search_info), intent(out) :: search
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    integer :: n
    if(count(spec%params%sampled)>0)then
      call initialize_parameter_search(spec,search,err,message)
    else
      ! The generic predicate also supports an entirely fixed specification.
      n=size(spec%params)
      search%params=spec%params
      allocate(search%search_index(n),search%param_names(0),search%lower(0),search%upper(0), &
               search%transformation(0),search%search_lower(0),search%search_upper(0))
      search%search_index=0
      err=0
      message=''
    endif
  end subroutine build_search

  subroutine prepare(ctx,err,message)
    type(model_calibration_context), intent(inout) :: ctx
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    call build_search(ctx%spec,ctx%search,err,message)
    if(err/=0) return
    call initialize_model_constraints(ctx,err,message)
  end subroutine prepare
end module test_support

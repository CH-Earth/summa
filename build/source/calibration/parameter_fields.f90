! Generic physical parameter fields: horizontal locations containing local vectors.
module parameter_fields
  use nr_type, only: i4b,rkind,lgt
  use parameter_search, only: parameter_spec,parameter_search_info, &
                             param_physical,param_spatial_multiplier,param_ratio
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  type, public :: parameter_values
    real(rkind), allocatable :: dat(:)
  end type parameter_values

  type, public :: spatial_parameter_field
    type(parameter_values), allocatable :: hspatial(:)
  end type spatial_parameter_field

  type, public :: parameter_field
    type(spatial_parameter_field) :: initial
    type(spatial_parameter_field) :: effective
  end type parameter_field

  ! Resolved once, rather than looking up names/sorting ratios for every proposal.
  type, public :: parameter_field_plan
    integer(i4b), allocatable :: reference(:),order(:)
    logical(lgt) :: ready=.false.
  end type parameter_field_plan

  public :: initialize_parameter_fields,apply_parameter_sample,reset_parameter_fields
  public :: field_is_finite,field_extrema
contains

  ! The adapter establishes physical correspondence. This routine verifies shapes,
  ! search mapping and ratio dependencies, then allocates reusable effective fields.
  ! The specification, search and initial fields must remain immutable after setup.
  subroutine initialize_parameter_fields(spec,search,field,plan,err,message)
    type(parameter_spec), intent(in) :: spec
    type(parameter_search_info), intent(in) :: search
    type(parameter_field), intent(inout) :: field(:)
    type(parameter_field_plan), intent(out) :: plan
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    integer(i4b) :: i,j,h,k,n,nx,nDone,pass
    logical(lgt) :: changed
    logical(lgt), allocatable :: complete(:),seen(:)

    err=20
    message='initialize_parameter_fields/'
    if(.not.allocated(spec%params))then
      message=trim(message)//'missing parameter specification'; return
    endif
    n=size(spec%params)
    if(n<1 .or. size(field)/=n)then
      message=trim(message)//'incorrect parameter field count'; return
    endif
    if(.not.allocated(search%param_names) .or. .not.allocated(search%search_index) .or. &
       .not.allocated(search%lower) .or. .not.allocated(search%upper) .or. &
       .not.allocated(search%transformation))then
      message=trim(message)//'search mapping/bounds not initialized'; return
    endif
    nx=size(search%param_names)
    if(size(search%search_index)/=n .or. size(search%lower)/=nx .or. &
       size(search%upper)/=nx .or. size(search%transformation)/=nx)then
      message=trim(message)//'inconsistent search sizes'; return
    endif
    allocate(plan%reference(n),plan%order(nx),complete(n),seen(nx),stat=err)
    if(err/=0)then
      message=trim(message)//'cannot allocate parameterization plan'; return
    endif
    err=20
    plan%reference=0
    seen=.false.

    do i=1,n
      if(len_trim(spec%params(i)%name)==0)then
        message=trim(message)//'empty parameter name'; return
      endif
      do j=1,i-1
        if(spec%params(j)%name==spec%params(i)%name)then
          message=trim(message)//'duplicate parameter: '//trim(spec%params(i)%name); return
        endif
      enddo
      if(.not.allocated(field(i)%initial%hspatial))then
        message=trim(message)//'missing initial field: '//trim(spec%params(i)%name); return
      endif
      if(size(field(i)%initial%hspatial)==0)then
        message=trim(message)//'empty spatial field: '//trim(spec%params(i)%name); return
      endif
      do h=1,size(field(i)%initial%hspatial)
        if(.not.allocated(field(i)%initial%hspatial(h)%dat))then
          message=trim(message)//'missing local data: '//trim(spec%params(i)%name); return
        endif
        if(size(field(i)%initial%hspatial(h)%dat)==0)then
          message=trim(message)//'empty local data: '//trim(spec%params(i)%name); return
        endif
      enddo
      if(.not.field_is_finite(field(i)%initial))then
        message=trim(message)//'non-finite initial field: '//trim(spec%params(i)%name); return
      endif

      ! Reject invalid metadata even for callers that bypass the TOML parser.
      if(trim(spec%params(i)%reference_reduction)/='none')then
        if(.not.spec%params(i)%sampled .or. spec%params(i)%param_type/=param_ratio)then
          message=trim(message)//'reference reduction requires a sampled ratio: '//trim(spec%params(i)%name); return
        endif
        if(trim(spec%params(i)%reference_reduction)/='min')then
          message=trim(message)//'unsupported reference reduction: '//trim(spec%params(i)%name); return
        endif
      endif

      k=search%search_index(i)
      if(.not.spec%params(i)%sampled)then
        if(k/=0)then
          message=trim(message)//'fixed parameter has a search index'; return
        endif
        cycle
      endif
      if(k<1 .or. k>nx)then
        message=trim(message)//'invalid search index: '//trim(spec%params(i)%name); return
      endif
      if(seen(k) .or. search%param_names(k)/=spec%params(i)%name)then
        message=trim(message)//'inconsistent search mapping: '//trim(spec%params(i)%name); return
      endif
      seen(k)=.true.
      if(.not.ieee_is_finite(spec%params(i)%lower) .or. .not.ieee_is_finite(spec%params(i)%upper))then
        message=trim(message)//'non-finite configured bounds'; return
      endif
      if(spec%params(i)%lower>spec%params(i)%upper)then
        message=trim(message)//'inverted configured bounds'; return
      endif
      ! Exact equality is intentional: do not silently accept old bound tightening.
      if(search%lower(k)/=spec%params(i)%lower .or. search%upper(k)/=spec%params(i)%upper)then
        message=trim(message)//'search bounds differ from configured bounds'; return
      endif
      if(search%transformation(k)/=spec%params(i)%transformation)then
        message=trim(message)//'inconsistent search transformation'; return
      endif
      select case(trim(spec%params(i)%transformation))
      case('none')
      case('log')
        if(spec%params(i)%lower<=0._rkind)then
          message=trim(message)//'log bounds must be positive'; return
        endif
      case default
        message=trim(message)//'unknown transformation'; return
      end select
      select case(spec%params(i)%param_type)
      case(param_physical,param_spatial_multiplier,param_ratio)
      case default
        message=trim(message)//'unknown parameterization'; return
      end select
    enddo
    if(.not.all(seen))then
      message=trim(message)//'unmapped search variable'; return
    endif

    do i=1,n
      if(.not.spec%params(i)%sampled .or. spec%params(i)%param_type/=param_ratio) cycle
      do j=1,n
        if(spec%params(j)%name==spec%params(i)%reference_name) exit
      enddo
      if(j>n .or. j==i)then
        message=trim(message)//'missing or self-referencing ratio: '//trim(spec%params(i)%name); return
      endif
      plan%reference(i)=j
      if(size(field(i)%initial%hspatial)/=size(field(j)%initial%hspatial))then
        message=trim(message)//'ratio reference has a different spatial count'; return
      endif
      do h=1,size(field(i)%initial%hspatial)
        select case(trim(spec%params(i)%reference_reduction))
        case('none')
          ! Preserve the existing elementwise contract; no implicit broadcasting.
          if(size(field(i)%initial%hspatial(h)%dat)/=size(field(j)%initial%hspatial(h)%dat))then
            message=trim(message)//'ratio reference has different local lengths'; return
          endif
        case('min')
          ! Reduce only the local reference vector; never combine spatial locations.
          if(size(field(i)%initial%hspatial(h)%dat)/=1)then
            message=trim(message)//'min reference reduction requires a local scalar target: '// &
                    trim(spec%params(i)%name); return
          endif
        end select
      enddo
    enddo

    complete=.not.spec%params%sampled
    nDone=0
    do pass=1,n
      changed=.false.
      do i=1,n
        if(complete(i)) cycle
        j=plan%reference(i)
        if(j>0)then
          if(.not.complete(j)) cycle
        endif
        nDone=nDone+1
        plan%order(nDone)=i
        complete(i)=.true.
        changed=.true.
      enddo
      if(all(complete)) exit
      if(.not.changed)then
        message=trim(message)//'cyclic parameter-ratio dependencies'; return
      endif
    enddo

    ! Deep allocation happens only at setup. Candidate construction uses slices.
    do i=1,n
      field(i)%effective=field(i)%initial
    enddo
    plan%ready=.true.
    err=0
    message=''
  end subroutine initialize_parameter_fields

  subroutine reset_parameter_fields(field)
    type(parameter_field), intent(inout) :: field(:)
    integer(i4b) :: i,h
    ! Requires successful initialization and unchanged shapes.
    do i=1,size(field)
      do h=1,size(field(i)%initial%hspatial)
        field(i)%effective%hspatial(h)%dat(:)=field(i)%initial%hspatial(h)%dat(:)
      enddo
    enddo
  end subroutine reset_parameter_fields

  ! Construct native-sized fields only. Ordering/bounds are checked separately.
  ! x contains untransformed calibration variables, NOT log coordinates.
  subroutine apply_parameter_sample(spec,search,x,field,plan,err,message)
    type(parameter_spec), intent(in) :: spec
    type(parameter_search_info), intent(in) :: search
    real(rkind), intent(in) :: x(:)
    type(parameter_field), intent(inout) :: field(:)
    type(parameter_field_plan), intent(in) :: plan
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    integer(i4b) :: k,i,j,h
    real(rkind) :: value
    err=20
    message='apply_parameter_sample/'
    if(.not.plan%ready)then
      message=trim(message)//'parameterization is not initialized'; return
    endif
    if(size(x)/=size(search%param_names) .or. size(field)/=size(spec%params))then
      message=trim(message)//'incorrect vector or field size'; return
    endif
    if(any(.not.ieee_is_finite(x)))then
      message=trim(message)//'non-finite calibration variable'; return
    endif
    call reset_parameter_fields(field)
    do k=1,size(plan%order)
      i=plan%order(k)
      j=plan%reference(i)
      value=x(search%search_index(i))
      do h=1,size(field(i)%effective%hspatial)
        select case(spec%params(i)%param_type)
        case(param_physical)
          field(i)%effective%hspatial(h)%dat(:)=value
        case(param_spatial_multiplier)
          field(i)%effective%hspatial(h)%dat(:)=value*field(i)%initial%hspatial(h)%dat(:)
        case(param_ratio)
          if(trim(spec%params(i)%reference_reduction)=='min')then
            ! The reference already includes its multiplier or upstream ratio.
            field(i)%effective%hspatial(h)%dat(1)= &
              value*minval(field(j)%effective%hspatial(h)%dat(:))
          else
            field(i)%effective%hspatial(h)%dat(:)=value*field(j)%effective%hspatial(h)%dat(:)
          endif
        end select
      enddo
    enddo
    err=0
    message=''
  end subroutine apply_parameter_sample

  logical(lgt) function field_is_finite(values)
    type(spatial_parameter_field), intent(in) :: values
    integer(i4b) :: h
    field_is_finite=.false.
    do h=1,size(values%hspatial)
      if(any(.not.ieee_is_finite(values%hspatial(h)%dat))) return
    enddo
    field_is_finite=.true.
  end function field_is_finite

  ! Whole-field extrema: do not change the gap to a per-location definition.
  subroutine field_extrema(values,lo,hi)
    type(spatial_parameter_field), intent(in) :: values
    real(rkind), intent(out) :: lo,hi
    integer(i4b) :: h
    lo=huge(1._rkind)
    hi=-huge(1._rkind)
    do h=1,size(values%hspatial)
      lo=min(lo,minval(values%hspatial(h)%dat))
      hi=max(hi,maxval(values%hspatial(h)%dat))
    enddo
  end subroutine field_extrema
end module parameter_fields

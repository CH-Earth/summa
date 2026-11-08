! Model-independent constraints on complete, untransformed calibration vectors.
! No sampling, bound refinement, SUMMA storage, or temporary broadcast arrays.
module parameter_constraints
  use nr_type, only: i4b,rkind,lgt
  use parameter_search, only: parameter_spec,parameter_search_info,param_physical
  use parameter_fields, only: parameter_field,parameter_field_plan,apply_parameter_sample, &
                             field_is_finite,field_extrema
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  type, public :: constraint_plan
    private
    logical(lgt), allocatable :: fixed(:),dynamic_gap(:)
    real(rkind), allocatable :: gap(:)
    logical(lgt) :: ready=.false.
  end type constraint_plan

  public :: initialize_parameter_constraints,check_parameter_constraints
contains

  ! Physical correspondence is an adapter contract. Equal lengths alone cannot
  ! distinguish different locations, soil layers, spectral bands, etc.
  subroutine initialize_parameter_constraints(spec,field,parameterization,plan,err,message)
    type(parameter_spec), intent(in) :: spec
    type(parameter_field), intent(in) :: field(:)
    type(parameter_field_plan), intent(in) :: parameterization
    type(constraint_plan), intent(out) :: plan
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    integer(i4b) :: c,j,i,a,b,h,n,nc,nMember,na,nb
    real(rkind) :: gap
    logical(lgt) :: valid

    err=20
    message='initialize_parameter_constraints/'
    if(.not.parameterization%ready)then
      message=trim(message)//'parameter fields are not initialized'; return
    endif
    n=size(spec%params)
    if(size(field)/=n)then
      message=trim(message)//'incorrect parameter field count'; return
    endif
    nc=0
    if(allocated(spec%ordered)) nc=size(spec%ordered)
    allocate(plan%fixed(nc),plan%dynamic_gap(nc),plan%gap(nc),stat=err)
    if(err/=0)then
      message=trim(message)//'cannot allocate constraint plan'; return
    endif
    err=20
    do c=1,nc
      if(.not.allocated(spec%ordered(c)%param_index))then
        message=trim(message)//'ordered chain has no parameter indices'; return
      endif
      nMember=size(spec%ordered(c)%param_index)
      if(nMember<2)then
        message=trim(message)//'ordered chain needs at least two parameters'; return
      endif
      if(.not.ieee_is_finite(spec%ordered(c)%gap_fraction))then
        message=trim(message)//'non-finite gap_fraction'; return
      endif
      if(spec%ordered(c)%gap_fraction<0._rkind)then
        message=trim(message)//'negative gap_fraction'; return
      endif
      plan%fixed(c)=.true.
      do j=1,nMember
        i=spec%ordered(c)%param_index(j)
        if(i<1 .or. i>n)then
          message=trim(message)//'ordered parameter index out of range'; return
        endif
        if(any(spec%ordered(c)%param_index(:j-1)==i))then
          message=trim(message)//'duplicate ordered parameter: '//trim(spec%params(i)%name); return
        endif
        if(spec%params(i)%sampled) plan%fixed(c)=.false.
      enddo
      do j=2,nMember
        a=spec%ordered(c)%param_index(j-1)
        b=spec%ordered(c)%param_index(j)
        if(size(field(a)%initial%hspatial)/=size(field(b)%initial%hspatial))then
          message=trim(message)//'ordered fields have different spatial counts'; return
        endif
        do h=1,size(field(a)%initial%hspatial)
          na=size(field(a)%initial%hspatial(h)%dat)
          nb=size(field(b)%initial%hspatial(h)%dat)
          if(na/=nb .and. na/=1 .and. nb/=1)then
            message=trim(message)//'ordered fields have incompatible local lengths'; return
          endif
        enddo
      enddo
      a=spec%ordered(c)%param_index(1)
      b=spec%ordered(c)%param_index(nMember)
      plan%dynamic_gap(c)=parameterization%reference(a)/=0 .or. parameterization%reference(b)/=0
      plan%gap(c)=0._rkind
      if(.not.plan%dynamic_gap(c))then
        call constraint_gap(spec,field,parameterization,c,gap,valid)
        if(.not.valid)then
          message=trim(message)//'invalid physical endpoint range or gap'; return
        endif
        plan%gap(c)=gap
      endif
      if(plan%fixed(c))then
        call check_chain(spec,field,c,plan%gap(c),valid,message)
        if(.not.valid)then
          message='fixed ordered constraint is infeasible: '//trim(message); return
        endif
      endif
    enddo
    plan%ready=.true.
    err=0
    message=''
  end subroutine initialize_parameter_constraints

  ! One candidate decision. A normal rejection has err=0. Configuration/API
  ! failures have err/=0 and must not be treated as resampling opportunities.
  ! Only effective fields are mutable scratch; no model state is accessed.
  subroutine check_parameter_constraints(spec,search,field,parameterization,plan,x,feasible,err,message)
    type(parameter_spec), intent(in) :: spec
    type(parameter_search_info), intent(in) :: search
    type(parameter_field), intent(inout) :: field(:)
    type(parameter_field_plan), intent(in) :: parameterization
    type(constraint_plan), intent(in) :: plan
    real(rkind), intent(in) :: x(:)
    logical(lgt), intent(out) :: feasible
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    integer(i4b) :: i,c
    real(rkind) :: gap,tol
    logical(lgt) :: valid

    feasible=.false.
    err=0
    message=''
    if(.not.plan%ready .or. .not.parameterization%ready)then
      err=20; message='feasibility context is not initialized'; return
    endif
    if(size(x)/=size(search%param_names))then
      err=20; message='incorrect calibration vector length'; return
    endif
    if(any(.not.ieee_is_finite(x)))then
      message='non-finite calibration variable'; return
    endif
    do i=1,size(x)
      tol=128._rkind*epsilon(1._rkind)*max(tiny(1._rkind),abs(x(i)), &
                                         abs(search%lower(i)),abs(search%upper(i)))
      if(x(i)<search%lower(i)-tol .or. x(i)>search%upper(i)+tol)then
        message='calibration bounds violated: '//trim(search%param_names(i)); return
      endif
      if(trim(search%transformation(i))=='log' .and. x(i)<=0._rkind)then
        message='non-positive log variable: '//trim(search%param_names(i)); return
      endif
    enddo
    call apply_parameter_sample(spec,search,x,field,parameterization,err,message)
    if(err/=0) return
    do i=1,size(field)
      if(.not.field_is_finite(field(i)%effective))then
        message='non-finite effective field: '//trim(spec%params(i)%name); return
      endif
    enddo
    do c=1,size(plan%fixed)
      if(plan%fixed(c)) cycle
      gap=plan%gap(c)
      if(plan%dynamic_gap(c))then
        call constraint_gap(spec,field,parameterization,c,gap,valid)
        if(.not.valid)then
          message='candidate gives an invalid physical endpoint range or gap'; return
        endif
      endif
      call check_chain(spec,field,c,gap,valid,message)
      if(.not.valid) return
    enddo
    feasible=.true.
    message=''
  end subroutine check_parameter_constraints

  ! Same whole-field envelope as the earlier implementation. Ratio endpoints use
  ! EFFECTIVE reference fields. These limits do not refine the search bounds.
  subroutine parameter_limits(spec,field,parameterization,i,lo,hi)
    type(parameter_spec), intent(in) :: spec
    type(parameter_field), intent(in) :: field(:)
    type(parameter_field_plan), intent(in) :: parameterization
    integer(i4b), intent(in) :: i
    real(rkind), intent(out) :: lo,hi
    integer(i4b) :: j,h
    real(rkind) :: smin,smax,local_reference,v(4)
    if(.not.spec%params(i)%sampled)then
      call field_extrema(field(i)%initial,lo,hi)
    else if(spec%params(i)%param_type==param_physical)then
      lo=spec%params(i)%lower
      hi=spec%params(i)%upper
    else
      j=parameterization%reference(i)
      if(j==0)then
        call field_extrema(field(i)%initial,smin,smax)
      else if(trim(spec%params(i)%reference_reduction)=='min')then
        ! Match construction: extrema ACROSS local minima, not all reference values.
        smin=huge(1._rkind)
        smax=-huge(1._rkind)
        do h=1,size(field(j)%effective%hspatial)
          local_reference=minval(field(j)%effective%hspatial(h)%dat(:))
          smin=min(smin,local_reference)
          smax=max(smax,local_reference)
        enddo
      else
        call field_extrema(field(j)%effective,smin,smax)
      endif
      v=[spec%params(i)%lower*smin,spec%params(i)%lower*smax, &
         spec%params(i)%upper*smin,spec%params(i)%upper*smax]
      lo=minval(v)
      hi=maxval(v)
    endif
  end subroutine parameter_limits

  subroutine constraint_gap(spec,field,parameterization,c,gap,valid)
    type(parameter_spec), intent(in) :: spec
    type(parameter_field), intent(in) :: field(:)
    type(parameter_field_plan), intent(in) :: parameterization
    integer(i4b), intent(in) :: c
    real(rkind), intent(out) :: gap
    logical(lgt), intent(out) :: valid
    integer(i4b) :: a,b,n
    real(rkind) :: lo,hi,unused
    n=size(spec%ordered(c)%param_index)
    a=spec%ordered(c)%param_index(1)
    b=spec%ordered(c)%param_index(n)
    call parameter_limits(spec,field,parameterization,a,lo,unused)
    call parameter_limits(spec,field,parameterization,b,unused,hi)
    gap=0._rkind
    valid=ieee_is_finite(lo).and.ieee_is_finite(hi).and.hi>=lo
    if(.not.valid) return
    if(spec%ordered(c)%gap_fraction>0._rkind) gap=spec%ordered(c)%gap_fraction*(hi-lo)
    valid=ieee_is_finite(gap)
  end subroutine constraint_gap

  ! Local singleton broadcasting only. Horizontal correspondence is established
  ! at setup; there is no automatic broadcasting between horizontal locations.
  subroutine check_chain(spec,field,c,gap,feasible,message)
    type(parameter_spec), intent(in) :: spec
    type(parameter_field), intent(in) :: field(:)
    integer(i4b), intent(in) :: c
    real(rkind), intent(in) :: gap
    logical(lgt), intent(out) :: feasible
    character(*), intent(out) :: message
    integer(i4b) :: j,a,b,h,k,na,nb
    real(rkind) :: pa,pb,tol
    character(len=256) :: detail
    feasible=.false.
    message=''
    do j=2,size(spec%ordered(c)%param_index)
      a=spec%ordered(c)%param_index(j-1)
      b=spec%ordered(c)%param_index(j)
      do h=1,size(field(a)%effective%hspatial)
        associate(va=>field(a)%effective%hspatial(h)%dat, vb=>field(b)%effective%hspatial(h)%dat)
          na=size(va)
          nb=size(vb)
          do k=1,max(na,nb)
            pa=va(min(k,na))
            pb=vb(min(k,nb))
            tol=128._rkind*epsilon(1._rkind)*max(tiny(1._rkind),abs(pa),abs(pb),abs(gap))
            if(pb-pa<gap-tol)then
              write(detail,'(A,I0,A,I0,A,I0,3(A,ES12.5))') &
                ' chain=',c,', spatial=',h,', local=',k, &
                ', lower=',pa,', upper=',pb,', gap=',gap
              message=trim(spec%params(a)%name)//' < '//trim(spec%params(b)%name)//trim(detail)
              return
            endif
          enddo
        end associate
      enddo
    enddo
    feasible=.true.
  end subroutine check_chain
end module parameter_constraints

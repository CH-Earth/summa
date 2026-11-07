! Same specifications/candidates against the earlier block prototype and Step 1.
! This file is a regression fixture, not part of the production modules.
program parity_fixture
  use nr_type, only: rkind,i4b,lgt
  use parameter_search
  use parameter_fields, only: parameter_field
  use model_calibration
  implicit none
  type(model_calibration_context) :: ctx
  integer :: i,j,k

  allocate(ctx%spec%params(4),ctx%spec%ordered(2),ctx%fields(4))
  call param(1,'theta_sat',.true.,.8_rkind,1.2_rkind,param_spatial_multiplier)
  call param(2,'residual',.false.,0._rkind,1._rkind,param_physical)
  call param(3,'wilting',.false.,0._rkind,1._rkind,param_physical)
  call param(4,'transpiration',.false.,0._rkind,1._rkind,param_physical)
  call set(ctx%fields(1),[.4_rkind,.42_rkind,.44_rkind],[.38_rkind,.41_rkind])
  call set(ctx%fields(2),[.04_rkind,.05_rkind,.06_rkind],[.03_rkind,.07_rkind])
  call set(ctx%fields(3),[.1_rkind],[.12_rkind])
  call set(ctx%fields(4),[.3_rkind],[.35_rkind])
  ctx%spec%ordered(1)%param_index=[2,3]; ctx%spec%ordered(1)%gap_fraction=.01_rkind
  ctx%spec%ordered(2)%param_index=[3,4,1]; ctx%spec%ordered(2)%gap_fraction=.025_rkind
  call prepare()
  do i=0,2000
    call report([.8_rkind+.4_rkind*i/2000._rkind])
  enddo

  ctx=model_calibration_context()
  allocate(ctx%spec%params(3),ctx%spec%ordered(1),ctx%fields(3))
  call param(1,'mixed_scales',.true.,-2._rkind,2._rkind,param_spatial_multiplier)
  call param(2,'low',.false.,-1._rkind,0._rkind,param_physical)
  call param(3,'high',.false.,0._rkind,1._rkind,param_physical)
  call set(ctx%fields(1),[-1._rkind,0._rkind,1._rkind])
  call set(ctx%fields(2),[-.5_rkind])
  call set(ctx%fields(3),[.5_rkind])
  ctx%spec%ordered(1)%param_index=[2,1,3]; ctx%spec%ordered(1)%gap_fraction=.1_rkind
  call prepare()
  do i=0,2000
    call report([-2._rkind+4._rkind*i/2000._rkind])
  enddo

  ctx=model_calibration_context()
  allocate(ctx%spec%params(3),ctx%spec%ordered(1),ctx%fields(3))
  call param(1,'reference',.true.,1._rkind,4._rkind,param_spatial_multiplier)
  call param(2,'ratio',.true.,.05_rkind,1._rkind,param_ratio)
  ctx%spec%params(2)%reference_name='reference'
  call param(3,'zero',.false.,0._rkind,1._rkind,param_physical)
  call set(ctx%fields(1),[2._rkind],[3._rkind])
  call set(ctx%fields(2),[.1_rkind],[.1_rkind])
  call set(ctx%fields(3),[0._rkind],[0._rkind])
  ctx%spec%ordered(1)%param_index=[3,2]; ctx%spec%ordered(1)%gap_fraction=.25_rkind
  call prepare()
  do i=0,40
    do j=0,40
      call report([1._rkind+3._rkind*i/40._rkind,.05_rkind+.95_rkind*j/40._rkind])
    enddo
  enddo
contains
  subroutine param(n,name,sampled,lo,hi,kind)
    integer, intent(in) :: n,kind
    character(*), intent(in) :: name
    logical, intent(in) :: sampled
    real(rkind), intent(in) :: lo,hi
    ctx%spec%params(n)%name=name; ctx%spec%params(n)%sampled=sampled
    ctx%spec%params(n)%lower=lo; ctx%spec%params(n)%upper=hi
    ctx%spec%params(n)%trial_value=lo/2+hi/2; ctx%spec%params(n)%param_type=kind
  end subroutine param
  subroutine set(f,a,b)
    type(parameter_field), intent(out) :: f
    real(rkind), intent(in) :: a(:)
    real(rkind), optional, intent(in) :: b(:)
#ifdef NESTED
    integer :: n
    n=1
    if(present(b)) n=2
    allocate(f%initial%hspatial(n))
    f%initial%hspatial(1)%dat=a
    if(present(b)) f%initial%hspatial(2)%dat=b
#else
    f%initial=a
    f%block_id=[1]
    f%offset=[1,size(a)+1]
    if(present(b))then
      f%initial=[a,b]
      f%block_id=[1,2]
      f%offset=[1,size(a)+1,size(a)+size(b)+1]
    endif
#endif
  end subroutine set
  subroutine prepare()
    type(parameter_spec) :: copy
    integer(i4b) :: err
    character(len=1024) :: message
    copy=ctx%spec
    if(allocated(copy%ordered)) deallocate(copy%ordered)
    call initialize_parameter_search(copy,ctx%search,err,message)
    if(err/=0)then; print *,message; stop 1; endif
    call initialize_model_constraints(ctx,err,message)
    if(err/=0)then; print *,message; stop 2; endif
  end subroutine prepare
  subroutine report(x)
    real(rkind), intent(in) :: x(:)
    integer(i4b) :: err
    character(len=1024) :: message
    logical(lgt) :: feasible
    real(rkind) :: total
#ifdef NESTED
    integer :: h
#endif
    call is_feasible(ctx,x,feasible,err,message)
    if(err/=0)then; print *,message; stop 3; endif
    total=0._rkind
    do k=1,size(ctx%fields)
#ifdef NESTED
      do h=1,size(ctx%fields(k)%effective%hspatial)
        total=total+sum(ctx%fields(k)%effective%hspatial(h)%dat)
      enddo
#else
      total=total+sum(ctx%fields(k)%effective)
#endif
    enddo
    write(*,'(I1,1X,ES25.17)') merge(1,0,feasible),total
  end subroutine report
end program parity_fixture

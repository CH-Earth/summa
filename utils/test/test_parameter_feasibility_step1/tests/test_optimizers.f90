! Integration tests: generic candidate generators use the SAME checker and
! perform rejection inside their existing loops (not in dispatcher wrappers).
program test_optimizers
  use test_support
  use parameter_search, only: param_physical,param_spatial_multiplier
  use parameter_search, only: search_state_type,generate_search_sample,sample_parameters
  use parameter_search, only: perturb_parameters_dds,perturb_parameters
  use nsga2, only: make_offspring
  use model_calibration, only: is_feasible
  implicit none
  type(model_calibration_context) :: ctx
  type(search_state_type) :: state
  integer(i4b) :: err,i,ncalls
  logical(lgt) :: feasible
  real(rkind) :: x(1),children(1,16)
  real(rkind) :: parents(1,4),distance(4)
  integer(i4b) :: rank(4)
  integer, allocatable :: seed(:)
  character(len=1024) :: message

  call setup(ctx,err,message)
  call require(err==0,'optimizer fixture setup: '//trim(message))
  call require(abs(ctx%search%lower(1)-0.8_rkind)<1.e-14_rkind,'configured lower bound not tightened')

  call random_seed(size=i)
  allocate(seed(i))
  seed=221
  call random_seed(put=seed)

  ncalls=0
  do i=1,100
    call generate_search_sample('random',ctx%search,state,1,250,x,counted_check,ctx,err,message)
    call require(err==0,'random candidate generation: '//trim(message))
    call is_feasible(ctx,x,feasible,err,message)
    call require(err==0.and.feasible,'random sampler returns only feasible vectors')
  enddo
  call require(ncalls>100,'random sampler retries rejected candidates internally')

  call generate_search_sample('dds',ctx%search,state,1,250,x,counted_check,ctx,err,message)
  call require(err==0,'DDS initializes from a feasible random vector')
  call is_feasible(ctx,x,feasible,err,message)
  call require(err==0.and.feasible,'DDS no-elite initialization is feasible')

  state%elite_individual=[1.10_rkind]
  ncalls=0
  do i=1,100
    call generate_search_sample('dds',ctx%search,state,i,250,x,counted_check,ctx,err,message)
    call require(err==0,'DDS proposal: '//trim(message))
    call is_feasible(ctx,x,feasible,err,message)
    call require(err==0.and.feasible,'DDS returns only feasible vectors')
  enddo
  call require(ncalls>=100,'DDS sampler invokes feasibility predicate')

  ncalls=0
  call perturb_parameters_dds(ctx%search,[1.10_rkind],10,250,0.2_rkind, &
                              x,counted_check,ctx,err,message)
  call require(err==0.and.ncalls>=1,'direct DDS sampler uses feasibility predicate')

  ncalls=0
  call perturb_parameters(ctx%search,[1.1_rkind],0.2_rkind, &
                          x,counted_check,ctx,err,message)
  call require(err==0.and.ncalls>=1,'generic neighborhood perturbation uses feasibility predicate')
  call is_feasible(ctx,x,feasible,err,message)
  call require(err==0.and.feasible,'generic perturbation returns feasible candidate')

  parents(1,:)=[1.02_rkind,1.08_rkind,1.14_rkind,1.19_rkind]
  rank=[1_i4b,1_i4b,2_i4b,2_i4b]
  distance=[0.1_rkind,0.2_rkind,0.3_rkind,0.4_rkind]
  ncalls=0
  call make_offspring(ctx%search,parents,rank,distance,0.9_rkind,10._rkind, &
                      1._rkind,15._rkind,children,counted_check,ctx,err,message)
  call require(err==0,'NSGA-II offspring generation: '//trim(message))
  call require(ncalls>=size(children,2),'NSGA-II checks every child')
  do i=1,size(children,2)
    call is_feasible(ctx,children(:,i),feasible,err,message)
    call require(err==0.and.feasible,'NSGA-II offspring is feasible')
  enddo

  ! Fatal checker errors are not silently converted into ordinary rejections.
  call sample_parameters(ctx%search,x,broken_check,ctx,err,message)
  call require(err==1234,'fatal feasibility callback error stops random sampler')
  call perturb_parameters_dds(ctx%search,[1.10_rkind],10,250,0.2_rkind, &
                              x,broken_check,ctx,err,message)
  call require(err==1234,'fatal feasibility callback error stops DDS')
  call perturb_parameters(ctx%search,[1.1_rkind],0.2_rkind, &
                          x,broken_check,ctx,err,message)
  call require(err==1234,'fatal feasibility callback error stops generic perturbation')
  call make_offspring(ctx%search,parents,rank,distance,0.9_rkind,10._rkind, &
                      1._rkind,15._rkind,children,broken_check,ctx,err,message)
  call require(err==1234,'fatal feasibility callback error stops NSGA-II offspring')

  print '(A,I0,A)', 'PASS: production DDS/random/NSGA-II feasibility integration (',assertions,' assertions).'
contains
  subroutine setup(context,err,message)
    type(model_calibration_context), intent(out) :: context
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    allocate(context%spec%params(2),context%spec%ordered(1),context%fields(2))
    call set_param(context%spec%params(1),'threshold',.false.,0._rkind,1._rkind,param_physical)
    call set_param(context%spec%params(2),'multiplier',.true.,0.8_rkind,1.2_rkind,param_spatial_multiplier)
    call set_field(context%fields(1),[0.78_rkind],[0.85_rkind])
    call set_field(context%fields(2),[0.80_rkind],[0.90_rkind])
    context%spec%ordered(1)%param_index=[1,2]
    context%spec%ordered(1)%gap_fraction=0.05_rkind
    call prepare(context,err,message)
  end subroutine setup

  subroutine counted_check(context,values,valid,code,detail)
    class(*), intent(inout) :: context
    real(rkind), intent(in) :: values(:)
    logical(lgt), intent(out) :: valid
    integer(i4b), intent(out) :: code
    character(*), intent(out) :: detail
    ncalls=ncalls+1
    call is_feasible(context,values,valid,code,detail)
  end subroutine counted_check

  subroutine broken_check(context,values,valid,code,detail)
    class(*), intent(inout) :: context
    real(rkind), intent(in) :: values(:)
    logical(lgt), intent(out) :: valid
    integer(i4b), intent(out) :: code
    character(*), intent(out) :: detail
    valid=.false.
    code=1234
    detail='deliberate fatal error'
  end subroutine broken_check
end program test_optimizers

program test_core
  use test_support
  use parameter_search, only: param_physical,param_spatial_multiplier,param_ratio
  use parameter_fields, only: field_extrema,apply_parameter_sample
  use model_calibration, only: is_feasible
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan,ieee_positive_inf
  implicit none
  call soil_tests()
  call scalar_and_local_tests()
  call ratio_tests()
  call signed_tests()
  call invalid_setup_tests()
  call numerical_tests()
  print '(A,I0,A)', 'PASS: nested feasibility core (',assertions,' assertions).'
contains
  subroutine make_soil(ctx,err,message)
    type(model_calibration_context), intent(out) :: ctx
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    allocate(ctx%spec%params(4),ctx%spec%ordered(2),ctx%fields(4))
    call set_param(ctx%spec%params(1),'theta_sat',.true.,.8_rkind,1.2_rkind,param_spatial_multiplier)
    call set_param(ctx%spec%params(2),'theta_res',.false.,0._rkind,1._rkind,param_physical)
    call set_param(ctx%spec%params(3),'wilting',.false.,0._rkind,1._rkind,param_physical)
    call set_param(ctx%spec%params(4),'transpiration',.false.,0._rkind,1._rkind,param_physical)
    call set_field(ctx%fields(1),[.40_rkind,.42_rkind,.44_rkind],[.38_rkind,.41_rkind])
    call set_field(ctx%fields(2),[.04_rkind,.05_rkind,.06_rkind],[.03_rkind,.07_rkind])
    call set_field(ctx%fields(3),[.10_rkind],[.12_rkind])
    call set_field(ctx%fields(4),[.30_rkind],[.35_rkind])
    ctx%spec%ordered(1)%param_index=[2,3]
    ctx%spec%ordered(1)%gap_fraction=.01_rkind
    ctx%spec%ordered(2)%param_index=[3,4,1]
    ctx%spec%ordered(2)%gap_fraction=.025_rkind
    call prepare(ctx,err,message)
  end subroutine make_soil

  subroutine soil_tests()
    type(model_calibration_context) :: ctx,bad
    integer(i4b) :: err,i
    character(len=1024) :: message
    logical(lgt) :: feasible
    real(rkind) :: x,limit,gap,lo,hi
    call make_soil(ctx,err,message)
    call require(err==0,'soil setup: '//trim(message))
    call require(size(ctx%fields(1)%initial%hspatial)==2,'two horizontal locations')
    call require(size(ctx%fields(1)%initial%hspatial(1)%dat)==3,'three local elements at first location')
    call require(size(ctx%fields(1)%initial%hspatial(2)%dat)==2,'two local elements at second location')
    call require(size(ctx%fields(3)%effective%hspatial(1)%dat)==1,'local scalar remains length one')
    call field_extrema(ctx%fields(1)%initial,lo,hi)
    call require(abs(lo-.38_rkind)<1.e-14_rkind.and.abs(hi-.44_rkind)<1.e-14_rkind,'whole-field extrema')
    call is_feasible(ctx,[.8_rkind],feasible,err,message)
    call require(err==0.and..not.feasible,'low multiplier rejected: '//trim(message))
    call is_feasible(ctx,[.94_rkind],feasible,err,message)
    call require(err==0.and..not.feasible,'second location constrains multiplier')
    call require(index(message,'spatial=2')>0.and.index(message,'local=1')>0,'spatial/local diagnostics')
    call is_feasible(ctx,[1.1_rkind],feasible,err,message)
    call require(err==0.and.feasible,'feasible multiplier')
    call require(all(abs(ctx%fields(1)%effective%hspatial(2)%dat-[.418_rkind,.451_rkind])<1.e-14_rkind), &
                 'multiplication at second location')
    call is_feasible(ctx,[1._rkind],feasible,err,message)
    call require(err==0.and.feasible,'next accepted sample')
    call require(all(abs(ctx%fields(1)%effective%hspatial(1)%dat-ctx%fields(1)%initial%hspatial(1)%dat) &
                     <1.e-14_rkind),'multipliers never accumulate')
    call require(abs(ctx%search%lower(1)-.8_rkind)<1.e-14_rkind,'no lower bound tightening')
    call require(abs(ctx%search%upper(1)-1.2_rkind)<1.e-14_rkind,'no upper bound tightening')

    ! Analytic feasible boundary; these comparisons are not bound-refinement code.
    gap=.025_rkind*(1.2_rkind*.44_rkind-.10_rkind)
    limit=(.35_rkind+gap)/.38_rkind
    do i=0,2000
      x=.8_rkind+.4_rkind*real(i,rkind)/2000._rkind
      call is_feasible(ctx,[x],feasible,err,message)
      call require(err==0.and.(feasible.eqv.(x>=limit)),'analytical soil feasibility grid')
    enddo
    call is_feasible(ctx,[limit],feasible,err,message)
    call require(err==0.and.feasible,'boundary candidate accepted')
    call is_feasible(ctx,[limit-1.e-9_rkind],feasible,err,message)
    call require(err==0.and..not.feasible,'outside feasible boundary rejected')
    bad=ctx
    bad%fields(2)%initial%hspatial(1)%dat(3)=.11_rkind
    call initialize_model_constraints(bad,err,message)
    call require(err/=0.and.index(message,'fixed ordered')>0,'fixed constraint checked before proposals')
    call is_feasible(bad,[1.1_rkind],feasible,err,message)
    call require(err/=0.and..not.feasible,'failed reinitialization cannot leave a ready checker')
    bad=ctx
    bad%fields(1)%initial%hspatial(2)%dat=[.50_rkind,.34_rkind]
    call initialize_model_constraints(bad,err,message)
    call require(err==0,'new layered case initializes without probing a lower-bound candidate')
    call is_feasible(bad,[1._rkind],feasible,err,message)
    call require(err==0.and..not.feasible.and.index(message,'local=2')>0,'all local depths checked')
    call require(all(abs(ctx%fields(2)%initial%hspatial(1)%dat-[.04_rkind,.05_rkind,.06_rkind]) &
                     <1.e-14_rkind),'initialized field preserved')
  end subroutine soil_tests

  subroutine scalar_and_local_tests()
    type(model_calibration_context) :: ctx
    integer(i4b) :: err,i
    character(len=1024) :: message
    logical(lgt) :: feasible
    allocate(ctx%spec%params(2),ctx%spec%ordered(1),ctx%fields(2))
    call set_param(ctx%spec%params(1),'a',.true.,0._rkind,1._rkind,param_physical)
    call set_param(ctx%spec%params(2),'b',.true.,0._rkind,1._rkind,param_physical)
    call set_field(ctx%fields(1),[.5_rkind])
    call set_field(ctx%fields(2),[.5_rkind])
    ctx%spec%ordered(1)%param_index=[1,2]
    ctx%spec%ordered(1)%gap_fraction=.1_rkind
    call prepare(ctx,err,message)
    call require(err==0,'scalar-only setup uses same checker')
    call is_feasible(ctx,[.2_rkind,.3_rkind],feasible,err,message)
    call require(err==0.and.feasible,'scalar boundary with floating tolerance')
    call is_feasible(ctx,[.2_rkind,.29_rkind],feasible,err,message)
    call require(err==0.and..not.feasible,'scalar gap violation')
    call require(abs(ctx%search%upper(1)-1._rkind)<1.e-14_rkind,'scalar bounds not propagated')
    ctx%spec%ordered(1)%gap_fraction=0._rkind
    call initialize_model_constraints(ctx,err,message)
    call is_feasible(ctx,[.3_rkind,.3_rkind],feasible,err,message)
    call require(err==0.and.feasible,'zero gap permits equality')

    ! The same local vector storage can describe spectral bands, without soil code.
    ctx=model_calibration_context()
    allocate(ctx%spec%params(2),ctx%spec%ordered(1),ctx%fields(2))
    call set_param(ctx%spec%params(1),'common_limit',.true.,0._rkind,1._rkind,param_physical)
    call set_param(ctx%spec%params(2),'band_values',.false.,0._rkind,1._rkind,param_physical)
    call set_field(ctx%fields(1),[0._rkind],[0._rkind])
    call set_field(ctx%fields(2),[.2_rkind,.6_rkind,.8_rkind],[.1_rkind,.4_rkind])
    ctx%spec%ordered(1)%param_index=[1,2]
    call prepare(ctx,err,message)
    call require(err==0,'spectral example setup')
    do i=0,100
      call is_feasible(ctx,[real(i,rkind)/100._rkind],feasible,err,message)
      call require(err==0.and.(feasible.eqv.(i<=10)),'spectral local scalar broadcasting')
    enddo
    call require(size(ctx%fields(1)%effective%hspatial(2)%dat)==1,'spectral scalar unexpanded')
    ctx%spec%params%sampled=.false.
    ctx%fields(1)%initial%hspatial(1)%dat=.05_rkind
    ctx%fields(1)%initial%hspatial(2)%dat=.05_rkind
    call prepare(ctx,err,message)
    call require(err==0,'entirely fixed context')
    call is_feasible(ctx,[real(rkind)::],feasible,err,message)
    call require(err==0.and.feasible,'empty calibration vector with fixed fields')
  end subroutine scalar_and_local_tests

  subroutine ratio_tests()
    type(model_calibration_context) :: ctx,bad
    integer(i4b) :: err
    character(len=1024) :: message
    logical(lgt) :: feasible
    allocate(ctx%spec%params(3),ctx%fields(3),ctx%spec%ordered(1))
    call set_param(ctx%spec%params(1),'bottom',.true.,.01_rkind,.5_rkind,param_ratio,'top')
    call set_param(ctx%spec%params(2),'top',.true.,.5_rkind,1.5_rkind,param_spatial_multiplier)
    call set_param(ctx%spec%params(3),'middle',.true.,1._rkind,1.5_rkind,param_ratio,'bottom')
    call set_field(ctx%fields(1),[.1_rkind],[10._rkind])
    call set_field(ctx%fields(2),[13._rkind],[16._rkind])
    call set_field(ctx%fields(3),[.2_rkind],[12._rkind])
    ctx%spec%ordered(1)%param_index=[1,3,2]
    call prepare(ctx,err,message)
    call require(err==0,'out-of-registry-order ratio references')
    call require(all(ctx%parameterization%order==[2,1,3]),'cached reference evaluation order')
    call is_feasible(ctx,[.25_rkind,1.2_rkind,1.5_rkind],feasible,err,message)
    call require(err==0.and.feasible,'nested ratio candidate')
    call require(abs(ctx%fields(2)%effective%hspatial(2)%dat(1)-19.2_rkind)<1.e-14_rkind,'top multiplier')
    call require(abs(ctx%fields(1)%effective%hspatial(2)%dat(1)-4.8_rkind)<1.e-14_rkind,'effective reference ratio')
    call require(abs(ctx%fields(3)%effective%hspatial(2)%dat(1)-7.2_rkind)<1.e-14_rkind,'ratio of ratio')
    call is_feasible(ctx,[.1_rkind,.8_rkind,1._rkind],feasible,err,message)
    call require(err==0.and.feasible,'ratio reset')
    call require(abs(ctx%fields(1)%effective%hspatial(1)%dat(1)-1.04_rkind)<1.e-14_rkind,'no ratio accumulation')
    bad=ctx
    bad%spec%params(2)%param_type=param_ratio
    bad%spec%params(2)%reference_name='middle'
    call initialize_model_constraints(bad,err,message)
    call require(err/=0.and.index(message,'cyclic')>0,'ratio cycle rejected')
    bad=ctx
    bad%spec%params(1)%reference_name='missing'
    call initialize_model_constraints(bad,err,message)
    call require(err/=0,'missing ratio reference')
    bad=ctx
    bad%spec%params(1)%reference_name='bottom'
    call initialize_model_constraints(bad,err,message)
    call require(err/=0,'self ratio reference')
    bad=ctx
    bad%fields(1)%initial%hspatial(1)%dat=[.1_rkind,.2_rkind]
    call initialize_model_constraints(bad,err,message)
    call require(err/=0,'ratio construction does not silently broadcast native arrays')

    ! Same effective value, different reference => different physical-envelope gap.
    ctx=model_calibration_context()
    allocate(ctx%spec%params(3),ctx%fields(3),ctx%spec%ordered(1))
    call set_param(ctx%spec%params(1),'reference',.true.,1._rkind,4._rkind,param_spatial_multiplier)
    call set_param(ctx%spec%params(2),'ratio',.true.,.05_rkind,1._rkind,param_ratio,'reference')
    call set_param(ctx%spec%params(3),'fixed_zero',.false.,0._rkind,1._rkind,param_physical)
    call set_field(ctx%fields(1),[2._rkind])
    call set_field(ctx%fields(2),[.1_rkind])
    call set_field(ctx%fields(3),[0._rkind])
    ctx%spec%ordered(1)%param_index=[3,2]
    ctx%spec%ordered(1)%gap_fraction=.25_rkind
    call prepare(ctx,err,message)
    call require(err==0,'do not reject setup for an infeasible lower-bound sample')
    call is_feasible(ctx,[1._rkind,.5_rkind],feasible,err,message)
    call require(err==0.and.feasible,'ratio gap 0.5 with effective ratio value 1')
    call is_feasible(ctx,[4._rkind,.125_rkind],feasible,err,message)
    call require(err==0.and..not.feasible,'ratio gap updates to 2 with same effective ratio value 1')
    call require(abs(ctx%fields(2)%effective%hspatial(1)%dat(1)-1._rkind)<1.e-14_rkind,'same physical value')
    ctx%spec%params(1)%sampled=.false.
    call prepare(ctx,err,message)
    call require(err==0,'ratio with an unsampled reference')
    call is_feasible(ctx,[.5_rkind],feasible,err,message)
    call require(err==0.and.feasible,'unsampled reference retains initialized field')
  end subroutine ratio_tests

  subroutine signed_tests()
    type(model_calibration_context) :: ctx
    integer(i4b) :: err,i
    character(len=1024) :: message
    logical(lgt) :: feasible
    real(rkind) :: x
    allocate(ctx%spec%params(3),ctx%spec%ordered(1),ctx%fields(3))
    call set_param(ctx%spec%params(1),'scaled',.true.,-2._rkind,2._rkind,param_spatial_multiplier)
    call set_param(ctx%spec%params(2),'low',.false.,-1._rkind,0._rkind,param_physical)
    call set_param(ctx%spec%params(3),'high',.false.,0._rkind,1._rkind,param_physical)
    call set_field(ctx%fields(1),[-1._rkind,0._rkind,1._rkind])
    call set_field(ctx%fields(2),[-.5_rkind])
    call set_field(ctx%fields(3),[.5_rkind])
    ctx%spec%ordered(1)%param_index=[2,1,3]
    call prepare(ctx,err,message)
    call require(err==0,'mixed-sign initial scales')
    do i=0,2000
      x=-2._rkind+4._rkind*real(i,rkind)/2000._rkind
      call is_feasible(ctx,[x],feasible,err,message)
      call require(err==0.and.(feasible.eqv.(abs(x)<=.5_rkind)),'signed/zero scale analytic feasibility')
    enddo
    ctx%spec%ordered(1)%gap_fraction=.1_rkind
    call initialize_model_constraints(ctx,err,message)
    call is_feasible(ctx,[.4_rkind],feasible,err,message)
    call require(err==0.and.feasible,'signed multiplier exact gap boundary')
    call is_feasible(ctx,[.41_rkind],feasible,err,message)
    call require(err==0.and..not.feasible,'signed multiplier positive gap')
  end subroutine signed_tests

  subroutine invalid_setup_tests()
    type(model_calibration_context) :: base,ctx
    integer(i4b) :: err
    character(len=1024) :: message
    logical(lgt) :: feasible
    integer :: other_context
    call make_soil(base,err,message)
    call require(err==0,'invalid setup test baseline')
    ctx=base
    deallocate(ctx%fields(1)%initial%hspatial)
    call initialize_model_constraints(ctx,err,message)
    call require(err/=0,'unallocated horizontal array')
    ctx=base
    deallocate(ctx%fields(1)%initial%hspatial)
    allocate(ctx%fields(1)%initial%hspatial(0))
    call initialize_model_constraints(ctx,err,message)
    call require(err/=0,'empty horizontal array')
    ctx=base
    deallocate(ctx%fields(1)%initial%hspatial(2)%dat)
    call initialize_model_constraints(ctx,err,message)
    call require(err/=0,'unallocated local array')
    ctx=base
    ctx%fields(1)%initial%hspatial(2)%dat=[real(rkind)::]
    call initialize_model_constraints(ctx,err,message)
    call require(err/=0,'empty local array')
    ctx=base
    call set_field(ctx%fields(3),[.1_rkind])
    call initialize_model_constraints(ctx,err,message)
    call require(err/=0,'no horizontal singleton broadcasting')
    ctx=base
    ctx%fields(3)%initial%hspatial(1)%dat=[.1_rkind,.1_rkind]
    call initialize_model_constraints(ctx,err,message)
    call require(err/=0,'2-to-3 local correspondence is invalid')
    ctx=base
    ctx%spec%ordered(2)%param_index=[3,4,4,1]
    call initialize_model_constraints(ctx,err,message)
    call require(err/=0,'duplicate ordered member')
    ctx=base
    ctx%spec%ordered(2)%param_index=[3,4,0]
    call initialize_model_constraints(ctx,err,message)
    call require(err/=0,'invalid parameter index')
    ctx=base
    deallocate(ctx%spec%ordered(1)%param_index)
    call initialize_model_constraints(ctx,err,message)
    call require(err/=0,'unallocated ordered members')
    ctx=base
    ctx%spec%ordered(1)%param_index=[1]
    call initialize_model_constraints(ctx,err,message)
    call require(err/=0,'one-member ordered chain')
    ctx=base
    ctx%spec%ordered(1)%gap_fraction=-.1_rkind
    call initialize_model_constraints(ctx,err,message)
    call require(err/=0,'negative minimum gap fraction')
    ctx=base
    ctx%search%lower(1)=.9_rkind
    call initialize_model_constraints(ctx,err,message)
    call require(err/=0,'silently tightened search bounds rejected')
    ctx=base
    ctx%search%search_index(2)=1
    call initialize_model_constraints(ctx,err,message)
    call require(err/=0,'fixed parameter cannot occupy sampled coordinate')
    ctx=base
    ctx%search%param_names(1)='wrong'
    call initialize_model_constraints(ctx,err,message)
    call require(err/=0,'search name mismatch')
    ctx=base
    ctx%spec%params(1)%param_type=99
    call initialize_model_constraints(ctx,err,message)
    call require(err/=0,'unknown parameterization')
    ctx=base
    ctx%spec%params(2)%name='theta_sat'
    call initialize_model_constraints(ctx,err,message)
    call require(err/=0,'duplicate registry name')
    ctx=base
    ctx%spec%params(1)%reference_name='unused'
    call initialize_model_constraints(ctx,err,message)
    call require(err==0,'irrelevant reference name does not change multiplier')
    ctx=model_calibration_context()
    call is_feasible(ctx,[1._rkind],feasible,err,message)
    call require(err/=0.and..not.feasible,'uninitialized context')
    other_context=1
    call is_feasible(other_context,[1._rkind],feasible,err,message)
    call require(err/=0.and..not.feasible,'incorrect opaque context type')
  end subroutine invalid_setup_tests

  subroutine numerical_tests()
    type(model_calibration_context) :: ctx,bad
    integer(i4b) :: err
    character(len=1024) :: message
    logical(lgt) :: feasible
    real(rkind) :: nanv,infv
    nanv=ieee_value(0._rkind,ieee_quiet_nan)
    infv=ieee_value(0._rkind,ieee_positive_inf)
    call make_soil(ctx,err,message)
    call is_feasible(ctx,[nanv],feasible,err,message)
    call require(err==0.and..not.feasible,'NaN candidate is ordinary rejection')
    call is_feasible(ctx,[infv],feasible,err,message)
    call require(err==0.and..not.feasible,'infinite candidate is ordinary rejection')
    call is_feasible(ctx,[1.3_rkind],feasible,err,message)
    call require(err==0.and..not.feasible,'configured bounds enforced')
    call is_feasible(ctx,[1._rkind,1._rkind],feasible,err,message)
    call require(err/=0.and..not.feasible,'wrong vector length is API error')
    bad=ctx
    bad%fields(1)%initial%hspatial(2)%dat(2)=nanv
    call initialize_model_constraints(bad,err,message)
    call require(err/=0,'non-finite initial value is setup error')
    bad=ctx
    bad%spec%ordered(1)%gap_fraction=nanv
    call initialize_model_constraints(bad,err,message)
    call require(err/=0,'non-finite gap is setup error')
    bad=ctx
    bad%spec%params(1)%lower=nanv
    call initialize_model_constraints(bad,err,message)
    call require(err/=0,'non-finite configured bounds are setup error')

    ctx=model_calibration_context()
    allocate(ctx%spec%params(1),ctx%fields(1))
    call set_param(ctx%spec%params(1),'conductivity',.true.,1.e-300_rkind,1._rkind,param_spatial_multiplier)
    ctx%spec%params(1)%transformation='log'
    call set_field(ctx%fields(1),[2._rkind])
    call prepare(ctx,err,message)
    call require(err==0,'positive log variable setup')
    call is_feasible(ctx,[0._rkind],feasible,err,message)
    call require(err==0.and..not.feasible,'log zero rejected even near tolerance')
    call is_feasible(ctx,[.1_rkind],feasible,err,message)
    call require(err==0.and.feasible,'checker receives untransformed log variable')
    call require(abs(ctx%fields(1)%effective%hspatial(1)%dat(1)-.2_rkind)<1.e-14_rkind,'no second exponentiation')
    call is_feasible(ctx,[-1.e-310_rkind],feasible,err,message)
    call require(err==0.and..not.feasible,'negative near-zero log variable rejected')

    ctx%spec%params(1)%transformation='none'
    ctx%spec%params(1)%lower=.01_rkind
    ctx%spec%params(1)%upper=3._rkind
    ctx%fields(1)%initial%hspatial(1)%dat=huge(1._rkind)/2._rkind
    call prepare(ctx,err,message)
    call require(err==0,'large finite initial value')
    call is_feasible(ctx,[3._rkind],feasible,err,message)
    call require(err==0.and..not.feasible,'overflow in constructed field rejected')
    call is_feasible(ctx,[.1_rkind],feasible,err,message)
    call require(err==0.and.feasible,'valid candidate after overflow rejection resets scratch')
  end subroutine numerical_tests
end program test_core

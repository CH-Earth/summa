program test_ratio_reduction
  use nr_type, only: i4b,rkind,lgt
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use parameter_search, only: parameter_spec,parameter_info,parameter_search_info, &
                             initialize_parameter_search,param_physical,param_spatial_multiplier,param_ratio
  use parameter_fields, only: parameter_field
  use model_calibration, only: model_calibration_context,initialize_model_constraints,is_feasible
  use summa_type, only: config_info
  use data_types, only: gru_hru_dom_doubleVec,gru_double
  use globalData, only: localParFallback,mpar_meta,iulog
  use var_lookup, only: iLookVarType,maxvarMpar
  use get_ixname_module, only: get_ixParam
  use summa_parameter_spec, only: get_summa_parameter_spec
  use summa_parameter_fields, only: summa_parameter_map,summa_calibration_adapter, &
                                   initialize_summa_calibration_adapter,apply_summa_calibration_adapter
  implicit none
  type(model_calibration_context) :: ctx
  type(config_info) :: config
  type(parameter_spec) :: spec
  type(parameter_search_info) :: search
  type(summa_parameter_map),allocatable :: mapping(:)
  type(summa_calibration_adapter) :: adapter
  type(gru_hru_dom_doubleVec) :: mpar,initial
  type(gru_double) :: bpar
  integer :: err,i,h,k,checks=0,it,iw,is,ir
  logical(lgt) :: feasible
  real(rkind) :: expected_s,expected_t,expected_w,scale
  character(len=1024) :: message

  ! Reverse registry order deliberately: W depends on T, which depends on S.
  call chain_fixture(ctx)
  call prepare(ctx)
  call require(err==0,'initialize chained scalar/vector ratios: '//trim(message))
  call require(all(ctx%parameterization%order==[3,2,1]),'dependency order is S,T,W')
  do i=1,5
    scale=0.8_rkind+0.1_rkind*real(i-1,rkind)
    call is_feasible(ctx,[0.5_rkind,0.8_rkind,scale],feasible,err,message)
    call require(err==0.and.feasible,'positive chain passes ordered constraints')
    do h=1,2
      expected_s=scale*minval(ctx%fields(3)%initial%hspatial(h)%dat)
      call close(ctx%fields(2)%effective%hspatial(h)%dat(1),0.8_rkind*expected_s,'local transpiration minimum')
      call close(ctx%fields(1)%effective%hspatial(h)%dat(1),0.4_rkind*expected_s,'downstream effective ratio')
      call require(size(ctx%fields(1)%effective%hspatial(h)%dat)==1,'wilting stays scalar')
      call require(size(ctx%fields(2)%effective%hspatial(h)%dat)==1,'transpiration stays scalar')
      call require(all(abs(ctx%fields(3)%effective%hspatial(h)%dat- &
                          scale*ctx%fields(3)%initial%hspatial(h)%dat)<1.e-13_rkind), &
                   'porosity stays layered and scales from initial values')
    enddo
  enddo

  ! Legacy, equal-size elementwise ratios keep their meaning. A missing reduction
  ! must not silently turn a previously invalid vector-to-scalar relation into a minimum.
  call chain_fixture(ctx)
  ctx%spec%params(2)%reference_reduction='none'
  call prepare(ctx)
  call require(err/=0,'mismatched ratio without explicit reduction is rejected')
  call chain_fixture(ctx)
  ctx%spec%params(2)%reference_reduction='mean'
  call prepare(ctx)
  call require(err/=0,'unsupported reduction is rejected')
  call chain_fixture(ctx)
  ctx%spec%params(3)%reference_reduction='min'
  call prepare(ctx)
  call require(err/=0,'reduction on a multiplier is rejected')
  call chain_fixture(ctx)
  ctx%fields(2)%initial%hspatial(1)%dat=[0.2_rkind,0.3_rkind]
  call prepare(ctx)
  call require(err/=0,'min reduction requires scalar target at every location')
  call chain_fixture(ctx)
  deallocate(ctx%fields(2)%initial%hspatial)
  allocate(ctx%fields(2)%initial%hspatial(1))
  ctx%fields(2)%initial%hspatial(1)%dat=[0.2_rkind]
  call prepare(ctx)
  call require(err/=0,'no horizontal broadcasting')
  call chain_fixture(ctx)
  ctx%spec%params(2)%reference_name='wilting'
  call prepare(ctx)
  call require(err/=0.and.index(message,'cyclic')>0,'cycle rejected before evaluation')
  call chain_fixture(ctx)
  ctx%spec%params(1)%reference_name='wilting'
  call prepare(ctx)
  call require(err/=0,'self reference rejected')

  ! Ratio endpoint gaps must use the extrema of local minima: min([2,100])=2,
  ! min([4,50])=4. With ratio bounds [0,1], the required gap is 0.1*4=0.4.
  ctx=model_calibration_context()
  allocate(ctx%spec%params(3),ctx%fields(3),ctx%spec%ordered(1))
  call param(ctx%spec%params(1),'reference',param_spatial_multiplier,0.8_rkind,1.2_rkind)
  call param(ctx%spec%params(2),'ratio',param_ratio,0._rkind,1._rkind,'reference','min')
  call param(ctx%spec%params(3),'zero',param_physical,0._rkind,1._rkind)
  ctx%spec%params(3)%sampled=.false.
  call field(ctx%fields(1),[2._rkind,100._rkind],[4._rkind,50._rkind])
  call field(ctx%fields(2),[0.5_rkind],[0.5_rkind])
  call field(ctx%fields(3),[0._rkind],[0._rkind])
  ctx%spec%ordered(1)%param_index=[3,2]
  ctx%spec%ordered(1)%gap_fraction=0.1_rkind
  call prepare(ctx)
  call require(err==0,'initialize reduced-reference gap fixture')
  call is_feasible(ctx,[1._rkind,0.25_rkind],feasible,err,message)
  call require(err==0.and.feasible,'gap uses local minima, not all layer values')
  call is_feasible(ctx,[1._rkind,0.15_rkind],feasible,err,message)
  call require(err==0.and..not.feasible,'gap uses maximum of local minima, not global minimum')
  call is_feasible(ctx,[0.8_rkind,0.25_rkind],feasible,err,message)
  call require(err==0.and.feasible,'gap follows current effective reference')

  ! Signed ratios multiply AFTER taking the minimum; positivity is not assumed by
  ! the generic numerical rule (physical ordering guarantees require positive fields).
  deallocate(ctx%spec%ordered)
  ctx%spec%params(2)%lower=-1._rkind
  ctx%spec%params(2)%upper=1._rkind
  ctx%fields(1)%initial%hspatial(1)%dat=[-4._rkind,-2._rkind]
  ctx%fields(1)%initial%hspatial(2)%dat=[-3._rkind,-1._rkind]
  call prepare(ctx)
  call require(err==0,'initialize signed reference fixture')
  call is_feasible(ctx,[1._rkind,-0.5_rkind],feasible,err,message)
  call require(err==0.and.feasible,'signed ratio construction is finite')
  call close(ctx%fields(2)%effective%hspatial(1)%dat(1),2._rkind,'reduction precedes multiplication')
  call is_feasible(ctx,[1._rkind,ieee_value(0._rkind,ieee_quiet_nan)],feasible,err,message)
  call require(err==0.and..not.feasible,'nonfinite ratio is ordinary rejection')

  ! A fixed layered reference is valid too.
  ctx%spec%params(1)%sampled=.false.
  call prepare(ctx)
  call require(err==0,'fixed reference is supported')
  call is_feasible(ctx,[-0.5_rkind],feasible,err,message)
  call require(err==0.and.feasible,'fixed reference candidate accepted')
  call close(ctx%fields(2)%effective%hspatial(2)%dat(1),1.5_rkind,'fixed local reference minimum')

  ! Exercise the REAL SUMMA specification builder and adapter with actual SUMMA
  ! data/lookup types. Only configuration/global-state holders are substituted.
  localParFallback%default_val=0.4_rkind
  localParFallback%lower_limit=0._rkind
  localParFallback%upper_limit=1._rkind
  mpar_meta%varUnit='-'
  mpar_meta%varDesc='synthetic fixture'
  mpar_meta%varType=iLookVarType%scalarv
  iw=get_ixParam('critSoilWilting'); it=get_ixParam('critSoilTranspire')
  is=get_ixParam('theta_sat'); ir=get_ixParam('theta_res')
  mpar_meta(is)%varType=iLookVarType%parSoil
  mpar_meta(ir)%varType=iLookVarType%parSoil
  allocate(config%calib%param_list(3),config%calib%multiplier(2),config%calib%ordered(2))
  config%calib%param_list%has_lower=.true.
  config%calib%param_list%has_upper=.true.
  config%calib%param_list(1)%name='critSoilWilting'
  config%calib%param_list(1)%param_type=param_ratio
  config%calib%param_list(1)%reference_name='critSoilTranspire'
  config%calib%param_list(1)%lower_bound=0.2_rkind
  config%calib%param_list(1)%upper_bound=0.9_rkind
  config%calib%param_list(2)%name='critSoilTranspire'
  config%calib%param_list(2)%param_type=param_ratio
  config%calib%param_list(2)%reference_name='theta_sat'
  config%calib%param_list(2)%reference_reduction='min'
  config%calib%param_list(2)%lower_bound=0.5_rkind
  config%calib%param_list(2)%upper_bound=0.95_rkind
  config%calib%param_list(3)%name='theta_sat'
  config%calib%param_list(3)%param_type=param_spatial_multiplier
  config%calib%param_list(3)%lower_bound=0.8_rkind
  config%calib%param_list(3)%upper_bound=1.2_rkind
  config%calib%multiplier(1)%parameters=[character(len=64)::'critSoilTranspire','critSoilWilting']
  config%calib%multiplier(2)%parameters=[character(len=64)::'theta_sat','critSoilTranspire']
  config%calib%ordered(1)%parameters=[character(len=64)::'theta_res','critSoilWilting']
  config%calib%ordered(1)%gap_fraction=0.01_rkind
  config%calib%ordered(2)%parameters=[character(len=64)::'critSoilWilting','critSoilTranspire','theta_sat']
  config%calib%ordered(2)%gap_fraction=0.025_rkind
  call get_summa_parameter_spec(config,spec,mapping,err=err,message=message)
  call require(err==0,'SUMMA spec accepts chained ratios in ordered constraints: '//trim(message))
  call require(spec%params(2)%reference_reduction=='min','SUMMA spec copies reduction metadata')
  call initialize_parameter_search(spec,search,err,message)
  call require(err==0,'SUMMA parameter search setup')
  allocate(mpar%gru(2))
  do h=1,2
    allocate(mpar%gru(h)%hru(1))
    allocate(mpar%gru(h)%hru(1)%dom(1))
    allocate(mpar%gru(h)%hru(1)%dom(1)%var(maxvarMpar))
    mpar%gru(h)%hru(1)%dom(1)%var(iw)%dat=[0.1_rkind]
    mpar%gru(h)%hru(1)%dom(1)%var(it)%dat=[0.2_rkind]
  enddo
  mpar%gru(1)%hru(1)%dom(1)%var(is)%dat=[0.4_rkind,0.5_rkind,0.45_rkind]
  mpar%gru(2)%hru(1)%dom(1)%var(is)%dat=[0.35_rkind,0.6_rkind]
  mpar%gru(1)%hru(1)%dom(1)%var(ir)%dat=[0.02_rkind,0.03_rkind,0.04_rkind]
  mpar%gru(2)%hru(1)%dom(1)%var(ir)%dat=[0.01_rkind,0.05_rkind]
  initial=mpar
  call initialize_summa_calibration_adapter(spec,search,mapping,mpar,bpar,adapter,err,message)
  call require(err==0,'SUMMA scalar/layer adapter initializes: '//trim(message))
  do k=1,3
    scale=1.2_rkind-0.1_rkind*real(k-1,rkind)
    adapter%calibration%values=[0.5_rkind,0.8_rkind,scale]
    call is_feasible(adapter%calibration,adapter%calibration%values,feasible,err,message)
    call require(err==0.and.feasible,'SUMMA candidate accepted: '//trim(message))
    call apply_summa_calibration_adapter(adapter,mpar,bpar,err,message)
    call require(err==0,'SUMMA write-back succeeds')
    do h=1,2
      expected_t=0.8_rkind*scale*minval(initial%gru(h)%hru(1)%dom(1)%var(is)%dat)
      expected_w=0.5_rkind*expected_t
      call close(mpar%gru(h)%hru(1)%dom(1)%var(it)%dat(1),expected_t,'native scalar transpiration')
      call close(mpar%gru(h)%hru(1)%dom(1)%var(iw)%dat(1),expected_w,'native scalar wilting')
      call require(size(mpar%gru(h)%hru(1)%dom(1)%var(it)%dat)==1,'native scalar allocation unchanged')
      call require(all(mpar%gru(h)%hru(1)%dom(1)%var(ir)%dat==initial%gru(h)%hru(1)%dom(1)%var(ir)%dat), &
                   'fixed residual field unchanged')
      call require(all(abs(mpar%gru(h)%hru(1)%dom(1)%var(is)%dat - &
                          scale*initial%gru(h)%hru(1)%dom(1)%var(is)%dat)<1.e-13_rkind), &
                   'no accumulation during repeated native write-back')
    enddo
  enddo
  ! Same lengths are insufficient when native metadata says the target is layered.
  mpar_meta(it)%varType=iLookVarType%parSoil
  call initialize_summa_calibration_adapter(spec,search,mapping,initial,bpar,adapter,err,message)
  call require(err/=0.and.index(message,'scalarv')>0,'reduced target must be a true native scalar')
  mpar_meta(it)%varType=iLookVarType%scalarv
  ! Old scalar/vector mismatch remains a useful error without opt-in reduction.
  spec%params(2)%reference_reduction='none'
  call initialize_summa_calibration_adapter(spec,search,mapping,initial,bpar,adapter,err,message)
  call require(err/=0,'SUMMA default rejects vector-to-scalar mismatch')
  spec%params(2)%reference_reduction='min'
  spec%params(2)%reference_name='critSoilWilting'
  call initialize_summa_calibration_adapter(spec,search,mapping,initial,bpar,adapter,err,message)
  call require(err/=0.and.index(message,'cyclic')>0,'SUMMA cycle still rejected during setup')

  print '(A,I0,A)', 'PASS: explicit local-minimum ratios and SUMMA chain integration (',checks,' assertions).'
contains
  subroutine require(ok,label)
    logical,intent(in) :: ok
    character(*),intent(in) :: label
    checks=checks+1
    if(.not.ok)then
      print '(A)', 'FAIL: '//trim(label)
      error stop 1
    endif
  end subroutine
  subroutine close(actual,expected,label)
    real(rkind),intent(in) :: actual,expected
    character(*),intent(in) :: label
    call require(abs(actual-expected)<1.e-13_rkind,label)
  end subroutine
  subroutine param(p,name,kind,lo,hi,reference,reduction)
    type(parameter_info),intent(out) :: p
    character(*),intent(in) :: name
    integer,intent(in) :: kind
    real(rkind),intent(in) :: lo,hi
    character(*),optional,intent(in) :: reference,reduction
    p%name=name;p%param_type=kind;p%sampled=.true.
    p%lower=lo;p%upper=hi;p%trial_value=(lo+hi)/2._rkind
    if(present(reference))p%reference_name=reference
    if(present(reduction))p%reference_reduction=reduction
  end subroutine
  subroutine field(f,a,b)
    type(parameter_field),intent(out) :: f
    real(rkind),intent(in) :: a(:),b(:)
    allocate(f%initial%hspatial(2))
    f%initial%hspatial(1)%dat=a;f%initial%hspatial(2)%dat=b
  end subroutine
  subroutine chain_fixture(c)
    type(model_calibration_context),intent(out) :: c
    allocate(c%spec%params(3),c%fields(3),c%spec%ordered(1))
    call param(c%spec%params(1),'wilting',param_ratio,0.2_rkind,0.9_rkind,'transpire')
    call param(c%spec%params(2),'transpire',param_ratio,0.5_rkind,0.95_rkind,'porosity','min')
    call param(c%spec%params(3),'porosity',param_spatial_multiplier,0.8_rkind,1.2_rkind)
    call field(c%fields(1),[0.1_rkind],[0.1_rkind])
    call field(c%fields(2),[0.2_rkind],[0.2_rkind])
    call field(c%fields(3),[0.4_rkind,0.5_rkind,0.45_rkind],[0.35_rkind,0.6_rkind])
    c%spec%ordered(1)%param_index=[1,2,3]
    c%spec%ordered(1)%gap_fraction=0.025_rkind
  end subroutine
  subroutine prepare(c)
    type(model_calibration_context),intent(inout) :: c
    call initialize_parameter_search(c%spec,c%search,err,message)
    if(err/=0)return
    call initialize_model_constraints(c,err,message)
  end subroutine
end program test_ratio_reduction

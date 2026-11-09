program test_adapter
  use test_support
  use nr_type, only: nr_integerMissing
  use parameter_search, only: param_physical,param_spatial_multiplier,param_ratio
  use data_types, only: gru_hru_dom_doubleVec,gru_double
  use globalData, only: gru_struc,mpar_meta
  use var_lookup, only: iLookVarType
  use summa_parameter_fields, only: summa_parameter_map,summa_calibration_adapter, &
                                   initialize_summa_calibration_adapter,apply_summa_calibration_adapter
  use model_calibration, only: is_feasible
  implicit none
  type(gru_hru_dom_doubleVec) :: mpar,before
  type(gru_double) :: bpar
  type(parameter_spec) :: spec
  type(parameter_search_info) :: search
  type(summa_parameter_map), allocatable :: map(:)
  type(summa_calibration_adapter) :: adapter
  integer(i4b) :: g,d,i,n,err
  character(len=1024) :: message
  logical(lgt) :: feasible

  call setup_model()
  allocate(spec%params(4),spec%ordered(2),map(4))
  call set_param(spec%params(1),'theta_sat',.true.,.8_rkind,1.2_rkind,param_spatial_multiplier)
  call set_param(spec%params(2),'theta_res',.false.,0._rkind,1._rkind,param_physical)
  call set_param(spec%params(3),'wilting',.false.,0._rkind,1._rkind,param_physical)
  call set_param(spec%params(4),'transpiration',.false.,0._rkind,1._rkind,param_physical)
  do i=1,4
    map(i)%ixM=i
  enddo
  spec%ordered(1)%param_index=[2,3]
  spec%ordered(1)%gap_fraction=.01_rkind
  spec%ordered(2)%param_index=[3,4,1]
  spec%ordered(2)%gap_fraction=.025_rkind
  call build_search(spec,search,err,message)
  call require(err==0,'adapter search setup')
  before=mpar
  call initialize_summa_calibration_adapter(spec,search,map,mpar,bpar,adapter,err,message)
  call require(err==0,'adapter setup: '//trim(message))

  ! Production lookups return a negative missing index for the unused family.
  map%ixB=nr_integerMissing
  call initialize_summa_calibration_adapter(spec,search,map,mpar,bpar,adapter,err,message)
  call require(err==0,'HRU mapping accepts missing basin index: '//trim(message))
  call require(size(adapter%calibration%fields(1)%initial%hspatial)==4,'one spatial entry per HRU/domain')
  call require(size(adapter%calibration%fields(1)%initial%hspatial(1)%dat)==3,'first domain has three layers')
  call require(size(adapter%calibration%fields(1)%initial%hspatial(2)%dat)==1,'second domain has one layer')
  call require(size(adapter%calibration%fields(1)%initial%hspatial(3)%dat)==2,'third domain has two layers')
  call require(size(adapter%calibration%fields(1)%initial%hspatial(4)%dat)==4,'fourth domain has four layers')
  call is_feasible(adapter%calibration,[.8_rkind],feasible,err,message)
  call require(err==0.and..not.feasible,'adapter rejects infeasible porosity multiplier')
  call require(same_native(mpar,before),'rejected candidate does not write SUMMA')
  call is_feasible(adapter%calibration,[1.1_rkind],feasible,err,message)
  call require(err==0.and.feasible,'adapter accepts valid porosity multiplier')
  call require(same_native(mpar,before),'accepted feasibility check also does not write SUMMA')
  adapter%calibration%values=[1.1_rkind]
  call apply_summa_calibration_adapter(adapter,mpar,bpar,err,message)
  call require(err==0,'write accepted fields: '//trim(message))
  do g=1,2
    do d=1,2
      call require(all(abs(mpar%gru(g)%hru(1)%dom(d)%var(1)%dat-.44_rkind)<1.e-14_rkind), &
                   'porosity written in every selected layer')
      call require(size(mpar%gru(g)%hru(1)%dom(d)%var(3)%dat)==1,'native threshold remains scalar')
      call require(all(abs(mpar%gru(g)%hru(1)%dom(d)%var(2)%dat-.05_rkind)<1.e-14_rkind),'fixed residual unchanged')
    enddo
  enddo
  adapter%calibration%values=[1.2_rkind]
  call is_feasible(adapter%calibration,adapter%calibration%values,feasible,err,message)
  call require(err==0.and.feasible,'second accepted candidate')
  call apply_summa_calibration_adapter(adapter,mpar,bpar,err,message)
  call require(err==0,'second write-back')
  call require(all(abs(mpar%gru(2)%hru(1)%dom(2)%var(1)%dat-.48_rkind)<1.e-14_rkind),'no native multiplier compounding')

  ! Domain-qualified fields stay nested. Non-selected domains are untouched.
  mpar=before
  map%domain_type=1
  call initialize_summa_calibration_adapter(spec,search,map,mpar,bpar,adapter,err,message)
  call require(err==0,'domain-filtered setup')
  call require(size(adapter%calibration%fields(1)%initial%hspatial)==2,'filtered spatial count')
  call is_feasible(adapter%calibration,[1.1_rkind],feasible,err,message)
  call require(err==0.and.feasible,'filtered feasibility')
  adapter%calibration%values=[1.1_rkind]
  call apply_summa_calibration_adapter(adapter,mpar,bpar,err,message)
  call require(err==0,'filtered write-back')
  do g=1,2
    call require(all(abs(mpar%gru(g)%hru(1)%dom(1)%var(1)%dat-.44_rkind)<1.e-14_rkind),'selected domain updated')
    call require(all(abs(mpar%gru(g)%hru(1)%dom(2)%var(1)%dat-.4_rkind)<1.e-14_rkind),'other domain unchanged')
  enddo

  ! Both domain selections have two locations: count equality is not enough.
  map(3)%domain_type=2
  call initialize_summa_calibration_adapter(spec,search,map,mpar,bpar,adapter,err,message)
  call require(err/=0.and.index(message,'different SUMMA')>0,'same count/different domain rejected')
  map(3)%domain_type=0
  map(3)%ixM=0
  map(3)%ixB=1
  call initialize_summa_calibration_adapter(spec,search,map,mpar,bpar,adapter,err,message)
  call require(err/=0.and.index(message,'different SUMMA')>0,'same count GRU/HRU mismatch rejected')
  map(3)%ixM=3
  map(3)%ixB=0
  map%domain_type=99
  call initialize_summa_calibration_adapter(spec,search,map,mpar,bpar,adapter,err,message)
  call require(err/=0,'empty selected domain rejected')
  map%domain_type=0
  map(1)%ixM=999
  call initialize_summa_calibration_adapter(spec,search,map,mpar,bpar,adapter,err,message)
  call require(err/=0,'out-of-range native parameter rejected without indexing it')
  map(1)%ixM=1
  map(1)%ixB=1
  call initialize_summa_calibration_adapter(spec,search,map,mpar,bpar,adapter,err,message)
  call require(err/=0,'ambiguous native parameter mapping rejected')
  ! Neither family may be missing: check both sentinel conventions and mixtures.
  map(1)%ixM=0
  map(1)%ixB=0
  call initialize_summa_calibration_adapter(spec,search,map,mpar,bpar,adapter,err,message)
  call require(err/=0,'both zero native indices rejected')
  map(1)%ixB=nr_integerMissing
  call initialize_summa_calibration_adapter(spec,search,map,mpar,bpar,adapter,err,message)
  call require(err/=0,'zero model and missing basin indices rejected')
  map(1)%ixM=nr_integerMissing
  call initialize_summa_calibration_adapter(spec,search,map,mpar,bpar,adapter,err,message)
  call require(err/=0,'both missing native indices rejected')
  map(1)%ixB=0
  call initialize_summa_calibration_adapter(spec,search,map,mpar,bpar,adapter,err,message)
  call require(err/=0,'missing model and zero basin indices rejected')
  map(1)%ixM=1
  map(1)%ixB=nr_integerMissing
  call initialize_summa_calibration_adapter(spec,search,map,mpar,bpar,adapter,err,message)
  call require(err==0,'valid native mapping restored')
  ! Equal-length nonscalar fields on different local axes must not compare.
  mpar_meta(2)%varType=iLookVarType%wLength
  spec%ordered(1)%param_index=[2,1]
  call initialize_summa_calibration_adapter(spec,search,map,mpar,bpar,adapter,err,message)
  call require(err/=0.and.index(message,'local axes')>0,'soil versus spectral axis rejected despite equal lengths')
  spec%ordered(1)%param_index=[2,3]
  mpar_meta(2)%varType=iLookVarType%parSoil

  ! Check all write destinations before writing any sampled parameter.
  mpar=before
  spec%params(3)%sampled=.true.
  spec%params(3)%lower=.10_rkind
  spec%params(3)%upper=.20_rkind
  call build_search(spec,search,err,message)
  call initialize_summa_calibration_adapter(spec,search,map,mpar,bpar,adapter,err,message)
  call require(err==0,'two sampled parameters for atomic write test')
  call is_feasible(adapter%calibration,[1.1_rkind,.15_rkind],feasible,err,message)
  call require(err==0.and.feasible,'two-parameter candidate accepted')
  adapter%calibration%values=[1.1_rkind,.15_rkind]
  mpar%gru(2)%hru(1)%dom(2)%var(3)%dat=[.12_rkind,.12_rkind]
  call apply_summa_calibration_adapter(adapter,mpar,bpar,err,message)
  call require(err/=0,'changed last destination shape detected')
  do g=1,2
    do d=1,2
      call require(all(abs(mpar%gru(g)%hru(1)%dom(d)%var(1)%dat-.4_rkind)<1.e-14_rkind), &
                   'earlier destinations not partly written')
    enddo
  enddo

  ! Basin parameters are also nested spatial fields with singleton local arrays.
  mpar=before
  spec=parameter_spec()
  deallocate(map)
  allocate(spec%params(2),spec%ordered(1),map(2))
  call set_param(spec%params(1),'basin_scaled',.true.,.8_rkind,1.2_rkind,param_spatial_multiplier)
  call set_param(spec%params(2),'basin_upper',.false.,0._rkind,1._rkind,param_physical)
  map(1)%ixB=1
  map(2)%ixB=2
  spec%ordered(1)%param_index=[1,2]
  spec%ordered(1)%gap_fraction=.025_rkind
  call build_search(spec,search,err,message)
  call initialize_summa_calibration_adapter(spec,search,map,mpar,bpar,adapter,err,message)
  call require(err==0,'basin-only adapter setup')
  map%ixM=nr_integerMissing
  call initialize_summa_calibration_adapter(spec,search,map,mpar,bpar,adapter,err,message)
  call require(err==0,'basin mapping accepts missing HRU index: '//trim(message))
  call require(size(adapter%calibration%fields(1)%initial%hspatial)==2,'one spatial entry per GRU')
  call is_feasible(adapter%calibration,[1.2_rkind],feasible,err,message)
  call require(err==0.and..not.feasible,'basin constraint enforced')
  call is_feasible(adapter%calibration,[1._rkind],feasible,err,message)
  call require(err==0.and.feasible,'basin feasible candidate')
  adapter%calibration%values=[1._rkind]
  call apply_summa_calibration_adapter(adapter,mpar,bpar,err,message)
  call require(err==0.and.same_native(mpar,before),'basin write does not alter HRU values')

  ! Joint multiplier/ratio construction using the same native adapter.
  spec=parameter_spec()
  allocate(spec%params(2),spec%ordered(1))
  call set_param(spec%params(1),'top',.true.,.5_rkind,1.5_rkind,param_spatial_multiplier)
  call set_param(spec%params(2),'bottom',.true.,.01_rkind,.5_rkind,param_ratio,'top')
  map(1)=summa_parameter_map(ixM=5)
  map(2)=summa_parameter_map(ixM=6)
  spec%ordered(1)%param_index=[2,1]
  call build_search(spec,search,err,message)
  call initialize_summa_calibration_adapter(spec,search,map,mpar,bpar,adapter,err,message)
  call require(err==0,'canopy ratio adapter setup')
  call is_feasible(adapter%calibration,[1.2_rkind,.25_rkind],feasible,err,message)
  call require(err==0.and.feasible,'joint multiplier/ratio feasible')
  adapter%calibration%values=[1.2_rkind,.25_rkind]
  call apply_summa_calibration_adapter(adapter,mpar,bpar,err,message)
  call require(err==0,'write joint multiplier/ratio')
  do g=1,2
    do d=1,2
      call require(abs(mpar%gru(g)%hru(1)%dom(d)%var(6)%dat(1)/ &
                       mpar%gru(g)%hru(1)%dom(d)%var(5)%dat(1)-.25_rkind)<1.e-14_rkind,'canopy ratio matches')
    enddo
  enddo
  map(1)%domain_type=1
  map(2)%domain_type=2
  call initialize_summa_calibration_adapter(spec,search,map,mpar,bpar,adapter,err,message)
  call require(err/=0.and.index(message,'ratio parameters')>0,'ratio wrong domain rejected')
  print '(A,I0,A)', 'PASS: SUMMA nested adapter (',assertions,' assertions).'
contains
  subroutine setup_model()
    mpar_meta(1:2)%varType=iLookVarType%parSoil
    mpar_meta(3:6)%varType=iLookVarType%scalarv
    allocate(mpar%gru(2),bpar%gru(2),gru_struc(2))
    do g=1,2
      allocate(mpar%gru(g)%hru(1),bpar%gru(g)%var(2),gru_struc(g)%hruInfo(1))
      allocate(mpar%gru(g)%hru(1)%dom(2),gru_struc(g)%hruInfo(1)%domInfo(2))
      do d=1,2
        gru_struc(g)%hruInfo(1)%domInfo(d)%dom_type=d
        allocate(mpar%gru(g)%hru(1)%dom(d)%var(6))
        if(d==1)then
          n=4-g
        else
          n=3*g-2
        endif
        do i=1,2
          allocate(mpar%gru(g)%hru(1)%dom(d)%var(i)%dat(n))
        enddo
        mpar%gru(g)%hru(1)%dom(d)%var(1)%dat=.4_rkind
        mpar%gru(g)%hru(1)%dom(d)%var(2)%dat=.05_rkind
        mpar%gru(g)%hru(1)%dom(d)%var(3)%dat=[.12_rkind]
        mpar%gru(g)%hru(1)%dom(d)%var(4)%dat=[.35_rkind]
        mpar%gru(g)%hru(1)%dom(d)%var(5)%dat=[real(10+g+d,rkind)]
        mpar%gru(g)%hru(1)%dom(d)%var(6)%dat=[.1_rkind]
      enddo
      bpar%gru(g)%var=[.25_rkind+.25_rkind*g,.8_rkind]
    enddo
  end subroutine setup_model

  logical function same_native(a,b)
    type(gru_hru_dom_doubleVec), intent(in) :: a,b
    integer :: gg,dd,ii
    same_native=.false.
    do gg=1,2
      do dd=1,2
        do ii=1,6
          if(size(a%gru(gg)%hru(1)%dom(dd)%var(ii)%dat)/=size(b%gru(gg)%hru(1)%dom(dd)%var(ii)%dat)) return
          if(any(abs(a%gru(gg)%hru(1)%dom(dd)%var(ii)%dat-b%gru(gg)%hru(1)%dom(dd)%var(ii)%dat) &
                 >1.e-14_rkind)) return
        enddo
      enddo
    enddo
    same_native=.true.
  end function same_native
end program test_adapter

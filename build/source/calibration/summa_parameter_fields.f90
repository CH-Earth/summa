module summa_parameter_fields
  use nr_type, only: i4b,rkind,lgt
  use data_types, only: gru_hru_dom_doubleVec,gru_double
  use globalData, only: gru_struc,mpar_meta
  use var_lookup, only: iLookVarType
  use parameter_search, only: parameter_spec,parameter_search_info,param_ratio
  use parameter_fields, only: parameter_field,spatial_parameter_field,field_is_finite,apply_parameter_sample
  use model_calibration, only: model_calibration_context,initialize_model_constraints
  implicit none
  private

  type, public :: summa_parameter_map
    integer(i4b) :: ixM=0,ixB=0,domain_type=0
  end type summa_parameter_map

  type, public :: summa_calibration_adapter
    type(model_calibration_context) :: calibration
    type(summa_parameter_map), allocatable :: map(:)
  end type summa_calibration_adapter

  public :: initialize_summa_calibration_adapter,apply_summa_calibration_adapter
  public :: initialize_summa_parameter_fields
contains

  subroutine initialize_summa_calibration_adapter(spec,search,param_map,mpar_init,bpar_init,adapter,err,message)
    type(parameter_spec), intent(in) :: spec
    type(parameter_search_info), intent(in) :: search
    type(summa_parameter_map), intent(in) :: param_map(:)
    type(gru_hru_dom_doubleVec), intent(in) :: mpar_init
    type(gru_double), intent(in) :: bpar_init
    type(summa_calibration_adapter), intent(out) :: adapter
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message

    adapter%calibration%spec=spec
    adapter%calibration%search=search
    adapter%map=param_map
    call initialize_summa_parameter_fields(spec,param_map,mpar_init,bpar_init, &
                                           adapter%calibration%fields,err,message)
    if(err/=0) return
    call initialize_model_constraints(adapter%calibration,err,message)
    if(err/=0) return
    allocate(adapter%calibration%values(size(search%param_names)),stat=err)
    if(err/=0)then
      message='cannot allocate calibration variables'; return
    endif
    ! This placeholder is not an accepted candidate; callers must set/check values.
    adapter%calibration%values=search%lower
    message=''
  end subroutine initialize_summa_calibration_adapter

  ! The ordering of selected HRU/domain entries is the horizontal correspondence.
  ! All paired parameters must select the same model space and domain qualifier.
  ! Local component meaning/order must also agree for nonscalar pairs: for example,
  ! corresponding soil layers, not a soil vector compared with spectral bands.
  subroutine initialize_summa_parameter_fields(spec,param_map,mpar_init,bpar_init,field,err,message)
    type(parameter_spec), intent(in) :: spec
    type(summa_parameter_map), intent(in) :: param_map(:)
    type(gru_hru_dom_doubleVec), intent(in) :: mpar_init
    type(gru_double), intent(in) :: bpar_init
    type(parameter_field), allocatable, intent(out) :: field(:)
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    integer(i4b) :: i,j,c,a,b,n

    err=20
    message='initialize_summa_parameter_fields/'
    if(.not.allocated(spec%params))then
      message=trim(message)//'missing parameter specification'; return
    endif
    n=size(spec%params)
    if(size(param_map)/=n)then
      message=trim(message)//'incorrect parameter map size'; return
    endif
    allocate(field(n),stat=err)
    if(err/=0)then
      message=trim(message)//'cannot allocate fields'; return
    endif
    do i=1,n
      call get_summa_parameter_field(param_map(i),mpar_init,bpar_init,field(i)%initial,err,message)
      if(err/=0) return
    enddo
    err=20
    do i=1,n
      if(.not.spec%params(i)%sampled .or. spec%params(i)%param_type/=param_ratio) cycle
      do j=1,n
        if(spec%params(j)%name==spec%params(i)%reference_name) exit
      enddo
      if(j>n)then
        message='ratio reference not found: '//trim(spec%params(i)%name); return
      endif
      if(.not.corresponding_maps(param_map(i),param_map(j)))then
        message='ratio parameters select different SUMMA domains or local axes'; return
      endif
      if(trim(spec%params(i)%reference_reduction)=='min')then
        if(param_map(i)%ixM>0)then
          ! A one-layer vector is not a native scalar parameter.
          if(mpar_meta(param_map(i)%ixM)%varType/=iLookVarType%scalarv)then
            message='min ratio reduction requires a SUMMA scalarv target'; return
          endif
        endif
      endif
    enddo
    if(allocated(spec%ordered))then
      do c=1,size(spec%ordered)
        if(.not.allocated(spec%ordered(c)%param_index))then
          message='ordered constraint has no parameter indices'; return
        endif
        if(any(spec%ordered(c)%param_index<1) .or. any(spec%ordered(c)%param_index>n))then
          message='ordered parameter index out of range'; return
        endif
        do j=2,size(spec%ordered(c)%param_index)
          a=spec%ordered(c)%param_index(j-1)
          b=spec%ordered(c)%param_index(j)
          if(.not.corresponding_maps(param_map(a),param_map(b)))then
            message='ordered parameters select different SUMMA domains or local axes'; return
          endif
        enddo
      enddo
    endif
    err=0
    message=''
  end subroutine initialize_summa_parameter_fields

  logical(lgt) function corresponding_maps(a,b)
    type(summa_parameter_map), intent(in) :: a,b
    integer(i4b) :: ta,tb
    corresponding_maps=.false.
    if((a%ixB>0).neqv.(b%ixB>0)) return
    if(a%domain_type/=b%domain_type) return
    if(a%ixM>0)then
      if(a%ixM>size(mpar_meta) .or. b%ixM>size(mpar_meta)) return
      ta=mpar_meta(a%ixM)%varType
      tb=mpar_meta(b%ixM)%varType
      if(ta<1 .or. tb<1 .or. ta==iLookVarType%unknown .or. tb==iLookVarType%unknown) return
      ! Equal lengths do not make a soil vector and a spectral vector comparable.
      if(ta/=tb .and. ta/=iLookVarType%scalarv .and. tb/=iLookVarType%scalarv) return
    endif
    corresponding_maps=.true.
  end function corresponding_maps

  logical(lgt) function selected_domain(map,g,h,d)
    type(summa_parameter_map), intent(in) :: map
    integer(i4b), intent(in) :: g,h,d
    selected_domain=.true.
    if(map%domain_type/=0) selected_domain=map%domain_type==gru_struc(g)%hruInfo(h)%domInfo(d)%dom_type
  end function selected_domain

  ! Validate the model-specific boundary, then copy native arrays directly.
  ! No flattened values, offsets or persistent location IDs are created.
  subroutine get_summa_parameter_field(map,mpar,bpar,values,err,message)
    type(summa_parameter_map), intent(in) :: map
    type(gru_hru_dom_doubleVec), intent(in) :: mpar
    type(gru_double), intent(in) :: bpar
    type(spatial_parameter_field), intent(out) :: values
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    integer(i4b) :: g,h,d,q,n

    err=20
    message='get_summa_parameter_field/invalid SUMMA storage or mapping'
    ! The unused family may be zero or SUMMA's negative missing-value sentinel.
    ! Exactly one family must have a positive index.
    if((map%ixM>0).eqv.(map%ixB>0)) return
    if(map%ixB>0)then
      if(map%domain_type/=0 .or. .not.allocated(bpar%gru)) return
      n=size(bpar%gru)
      if(n<1) return
      do g=1,n
        if(.not.allocated(bpar%gru(g)%var)) return
        if(map%ixB>size(bpar%gru(g)%var)) return
      enddo
      allocate(values%hspatial(n),stat=err)
      if(err/=0) return
      do g=1,n
        values%hspatial(g)%dat=[bpar%gru(g)%var(map%ixB)]
      enddo
    else
      if(.not.allocated(mpar%gru)) return
      if(map%ixM>size(mpar_meta)) return
      if(map%domain_type/=0)then
        if(.not.allocated(gru_struc)) return
        if(size(gru_struc)/=size(mpar%gru)) return
      endif
      n=0
      do g=1,size(mpar%gru)
        if(.not.allocated(mpar%gru(g)%hru)) return
        if(map%domain_type/=0)then
          if(.not.allocated(gru_struc(g)%hruInfo)) return
          if(size(gru_struc(g)%hruInfo)/=size(mpar%gru(g)%hru)) return
        endif
        do h=1,size(mpar%gru(g)%hru)
          if(.not.allocated(mpar%gru(g)%hru(h)%dom)) return
          if(map%domain_type/=0)then
            if(.not.allocated(gru_struc(g)%hruInfo(h)%domInfo)) return
            if(size(gru_struc(g)%hruInfo(h)%domInfo)/=size(mpar%gru(g)%hru(h)%dom)) return
          endif
          do d=1,size(mpar%gru(g)%hru(h)%dom)
            if(.not.selected_domain(map,g,h,d)) cycle
            if(.not.allocated(mpar%gru(g)%hru(h)%dom(d)%var)) return
            if(map%ixM>size(mpar%gru(g)%hru(h)%dom(d)%var)) return
            if(.not.allocated(mpar%gru(g)%hru(h)%dom(d)%var(map%ixM)%dat)) return
            if(size(mpar%gru(g)%hru(h)%dom(d)%var(map%ixM)%dat)<1) return
            if(mpar_meta(map%ixM)%varType==iLookVarType%scalarv)then
              if(size(mpar%gru(g)%hru(h)%dom(d)%var(map%ixM)%dat)/=1) return
            endif
            n=n+1
          enddo
        enddo
      enddo
      if(n<1)then
        message='parameter selects no SUMMA locations'; return
      endif
      allocate(values%hspatial(n),stat=err)
      if(err/=0) return
      q=0
      do g=1,size(mpar%gru)
        do h=1,size(mpar%gru(g)%hru)
          do d=1,size(mpar%gru(g)%hru(h)%dom)
            if(.not.selected_domain(map,g,h,d)) cycle
            q=q+1
            values%hspatial(q)%dat=mpar%gru(g)%hru(h)%dom(d)%var(map%ixM)%dat
          enddo
        enddo
      enddo
    endif
    err=0
    message=''
  end subroutine get_summa_parameter_field

  ! Called only after all sampled destinations have been verified. Slice
  ! assignment preserves SUMMA's native dimensions and allocation.
  subroutine put_summa_parameter_field(map,values,mpar,bpar)
    type(summa_parameter_map), intent(in) :: map
    type(spatial_parameter_field), intent(in) :: values
    type(gru_hru_dom_doubleVec), intent(inout) :: mpar
    type(gru_double), intent(inout) :: bpar
    integer(i4b) :: g,h,d,q
    if(map%ixB>0)then
      do g=1,size(bpar%gru)
        bpar%gru(g)%var(map%ixB)=values%hspatial(g)%dat(1)
      enddo
      return
    endif
    q=0
    do g=1,size(mpar%gru)
      do h=1,size(mpar%gru(g)%hru)
        do d=1,size(mpar%gru(g)%hru(h)%dom)
          if(.not.selected_domain(map,g,h,d)) cycle
          q=q+1
          mpar%gru(g)%hru(h)%dom(d)%var(map%ixM)%dat(:)=values%hspatial(q)%dat(:)
        enddo
      enddo
    enddo
  end subroutine put_summa_parameter_field

  ! Apply an already accepted vector. It is the caller's responsibility to call
  ! is_feasible first. This routine retains structural/finite-value safety checks,
  ! but does not re-evaluate calibration ordered inequalities.
  subroutine apply_summa_calibration_adapter(adapter,mpar,bpar,err,message)
    type(summa_calibration_adapter), intent(inout) :: adapter
    type(gru_hru_dom_doubleVec), intent(inout) :: mpar
    type(gru_double), intent(inout) :: bpar
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    type(spatial_parameter_field) :: destination
    integer(i4b) :: i,h

    if(.not.allocated(adapter%map) .or. .not.allocated(adapter%calibration%values))then
      err=20; message='SUMMA adapter is not initialized'; return
    endif
    if(.not.allocated(adapter%calibration%fields))then
      err=20; message='SUMMA calibration fields are missing'; return
    endif
    if(size(adapter%map)/=size(adapter%calibration%fields))then
      err=20; message='SUMMA parameter mapping size changed'; return
    endif
    call apply_parameter_sample(adapter%calibration%spec,adapter%calibration%search, &
                                adapter%calibration%values,adapter%calibration%fields, &
                                adapter%calibration%parameterization,err,message)
    if(err/=0) return
    ! Validate every destination before the first model write.
    do i=1,size(adapter%map)
      if(.not.adapter%calibration%spec%params(i)%sampled) cycle
      if(.not.field_is_finite(adapter%calibration%fields(i)%effective))then
        err=20; message='non-finite effective SUMMA field'; return
      endif
      call get_summa_parameter_field(adapter%map(i),mpar,bpar,destination,err,message)
      if(err/=0) return
      err=20
      message='SUMMA destination shape changed'
      if(size(destination%hspatial)/=size(adapter%calibration%fields(i)%effective%hspatial)) return
      do h=1,size(destination%hspatial)
        if(size(destination%hspatial(h)%dat)/=size(adapter%calibration%fields(i)%effective%hspatial(h)%dat)) return
      enddo
    enddo
    do i=1,size(adapter%map)
      if(.not.adapter%calibration%spec%params(i)%sampled) cycle
      call put_summa_parameter_field(adapter%map(i),adapter%calibration%fields(i)%effective,mpar,bpar)
    enddo
    err=0
    message=''
  end subroutine apply_summa_calibration_adapter
end module summa_parameter_fields

module summa_parameter_fields

  USE nr_type,    only: i4b,rkind
  USE data_types, only: gru_hru_dom_doubleVec,gru_double
  USE globalData, only: gru_struc

  USE parameter_search,  only: parameter_spec,parameter_search_info,param_ratio
  USE parameter_fields,  only: parameter_field
  USE parameter_search,  only: param_physical
  USE parameter_fields,  only: construct_effective_parameter_fields => apply_parameter_sample

  USE model_calibration, only: model_calibration_context

  implicit none
  private

  ! **************************************************************************************************
  ! SUMMA-specific adapter for the model-independent calibration context.
  !
  ! The generic context holds the calibration specification, search variables,
  ! and physical parameter fields. The SUMMA map identifies where those fields
  ! are stored in the model.
  ! **************************************************************************************************
  
  ! Location of a parameter within the SUMMA data structures.
  type, public :: summa_parameter_map
    integer(i4b) :: ixM = 0           ! Model parameter index, or zero
    integer(i4b) :: ixB = 0           ! Basin parameter index, or zero
    integer(i4b) :: domain_type = 0   ! Domain qualifier, or zero for all domains
  end type summa_parameter_map

  ! calibration adaptor 
  type, public :: summa_calibration_adapter
    type(model_calibration_context) :: calibration
    type(summa_parameter_map), allocatable :: map(:)
  end type summa_calibration_adapter

  public :: initialize_summa_parameter_fields
  public :: get_summa_parameter_vector
  public :: put_summa_parameter_vector
  public :: apply_summa_parameter_sample

  public :: initialize_summa_calibration_adapter
  public :: apply_summa_calibration_adapter

contains

  
  ! **************************************************************************************************
  ! Initialize the SUMMA calibration adapter.
  !
  ! Store the generic parameter specification and search information, retain
  ! the SUMMA parameter mapping, and extract the initialized physical fields.
  !
  ! Called once per calibration case on each model worker, after the common
  ! cold-start spinup. Subsequent trials reuse these initialized fields.
  ! **************************************************************************************************
  subroutine initialize_summa_calibration_adapter(spec,search,param_map,        &
                                                  mpar_init,bpar_init,adapter,  &
                                                  err,message)
    implicit none
    type(parameter_spec)            , intent(in)  :: spec          ! Generic parameter specification
    type(parameter_search_info)     , intent(in)  :: search        ! Parameter search information
    type(summa_parameter_map)       , intent(in)  :: param_map(:)  ! Mapping to SUMMA parameter locations
    type(gru_hru_dom_doubleVec)     , intent(in)  :: mpar_init     ! Initialized model parameter fields
    type(gru_double)                , intent(in)  :: bpar_init     ! Initialized basin parameter fields
    type(summa_calibration_adapter) , intent(out) :: adapter       ! Initialized SUMMA calibration adapter
    integer(i4b)                    , intent(out) :: err           ! Error code
    character(*)                    , intent(out) :: message       ! Error message

    character(len=1024) :: cmessage

    err=0
    message='initialize_summa_calibration_adapter/'

    ! Retain generic calibration information and SUMMA-specific parameter locations.
    adapter%calibration%spec=spec
    adapter%calibration%search=search
    adapter%map=param_map

    ! Extract initialized fields, including non-sampled reference parameters.
    call initialize_summa_parameter_fields(spec,param_map,mpar_init,bpar_init, &
                                           adapter%calibration%fields,err,cmessage)
    if(err/=0)then
      message=trim(message)//trim(cmessage)
      return
    endif

    ! Allocate the scalar calibration variables received from the dispatcher.
    allocate(adapter%calibration%values(size(search%param_names)),stat=err)
    if(err/=0)then
      message=trim(message)//'unable to allocate calibration variables'
      return
    endif

  end subroutine initialize_summa_calibration_adapter

  ! **************************************************************************************************
  ! Construct and apply the calibrated physical parameter fields to SUMMA.
  ! **************************************************************************************************
  subroutine apply_summa_calibration_adapter(adapter,mpar,bpar,err,message)
    implicit none
    ! dummy
    type(summa_calibration_adapter), intent(inout) :: adapter  ! Calibration context and parameter mapping
    type(gru_hru_dom_doubleVec)    , intent(inout) :: mpar     ! SUMMA model parameter fields
    type(gru_double)               , intent(inout) :: bpar     ! SUMMA basin parameter fields
    integer(i4b)                   , intent(out)   :: err      ! Error code
    character(*)                   , intent(out)   :: message  ! Error message

    err=0
    message='apply_summa_calibration_adapter/'

    if(.not.allocated(adapter%map))then
      message=trim(message)//'parameter mapping is not initialized'
      err=20; return
    endif

    if(.not.allocated(adapter%calibration%fields) .or. &
       .not.allocated(adapter%calibration%values))then
      message=trim(message)//'calibration fields or values are not initialized'
      err=20; return
    endif

    call apply_summa_parameter_sample(                 &
      adapter%calibration%spec,                        &
      adapter%calibration%search,                      &
      adapter%calibration%values,                      &
      adapter%map,                                     &
      adapter%calibration%fields,                      &
      mpar,bpar,err,message)

  end subroutine apply_summa_calibration_adapter


  ! **************************************************************************************************
  ! Extract and validate the initialized SUMMA parameter fields.
  !
  ! All parameters in spec%params are included, even when not sampled, because
  ! fixed parameters may serve as references for calibrated parameter ratios.
  !
  ! Parameter ratios require corresponding physical values at each spatial
  ! element. This routine verifies that reference and dependent fields have
  ! compatible domain selections and array dimensions.
  !
  ! The initialized fields are retained separately from the effective fields
  ! so that subsequent calibration trials do not accumulate multipliers.
  ! **************************************************************************************************
  subroutine initialize_summa_parameter_fields(spec,param_map,mpar_init,bpar_init,field,err,message)

    implicit none

    ! Dummy arguments
    type(parameter_spec),               intent(in)  :: spec          ! Generic calibration parameter specification
    type(summa_parameter_map),          intent(in)  :: param_map(:)  ! Mapping between calibration and SUMMA parameters
    type(gru_hru_dom_doubleVec),         intent(in)  :: mpar_init     ! Initialized SUMMA model parameter fields
    type(gru_double),                   intent(in)  :: bpar_init     ! Initialized SUMMA basin parameter fields
    type(parameter_field), allocatable, intent(out) :: field(:)      ! Flattened initial and effective parameter fields
    integer(i4b),                       intent(out) :: err           ! Error code
    character(*),                       intent(out) :: message       ! Error message

    ! Local variables
    integer(i4b) :: i,j           ! indices
    integer(i4b) :: g,h,d         ! GRU index, HRU index, Domain index
    integer(i4b) :: ixReference   ! Reference parameter index for parameter ratios
    integer(i4b) :: iConstraint   ! Ordered constraint index
    integer(i4b) :: iMember       ! Parameter position within an ordered constraint
    integer(i4b) :: ix            ! Parameter registry index
    real(rkind)  :: tolerance     ! Tolerance for comparing initialized scalar values

    err=0
    message='initialize_summa_parameter_fields/'

    if(.not.allocated(spec%params))then
      message=trim(message)//'parameter specification is not allocated'
      err=20; return
    endif

    if(size(param_map)/=size(spec%params))then
      message=trim(message)//'parameter map size does not match parameter specification'
      err=20; return
    endif
   
    ! Allocate the initialized parameter fields.
    allocate(field(size(spec%params)),stat=err)
    if(err/=0)then
      message=trim(message)//'unable to allocate SUMMA parameter fields'
      return
    endif

    ! -----------------------------------------------------------------------------------------------
    ! Extract all initialized parameter fields, including fixed references.
    ! -----------------------------------------------------------------------------------------------
    do i=1,size(field)

      if((param_map(i)%ixM>0) .eqv. (param_map(i)%ixB>0))then
        message=trim(message)//'invalid SUMMA parameter map for '//trim(spec%params(i)%name)
        err=20; return
      endif

      if(param_map(i)%ixB>0 .and. param_map(i)%domain_type/=0)then
        message=trim(message)//'basin parameter cannot have a domain qualifier: '//trim(spec%params(i)%name)
        err=20; return
      endif

      call get_summa_parameter_vector(param_map(i),mpar_init,bpar_init,field(i)%initial)

      if(size(field(i)%initial)==0)then
        message=trim(message)//'no initialized spatial elements for '//trim(spec%params(i)%name)
        err=20; return
      endif

    enddo

    ! -----------------------------------------------------------------------------------------------
    ! Support ordered constraints involving uniform physical parameters.
    !
    ! The generic parameter search already enforces these constraints using
    ! scalar physical values. This is valid only when non-sampled parameters
    ! have initialized fields consistent with their scalar trial values.
    !
    ! Constraints involving spatial multipliers, parameter ratios, or
    ! spatially varying fixed parameters require field-based handling.
    ! -----------------------------------------------------------------------------------------------
    if(allocated(spec%ordered))then

      do iConstraint=1,size(spec%ordered)

        do iMember=1,size(spec%ordered(iConstraint)%param_index)

          ix=spec%ordered(iConstraint)%param_index(iMember)

          ! Only uniform physical parameterizations can use scalar constraints.
          if(spec%params(ix)%param_type/=param_physical)then
            message=trim(message)//'field-based ordered constraints are not yet implemented for '// &
                    trim(spec%params(ix)%name)
            err=20; return
          endif

          ! Sampled physical parameters are assigned uniformly during each trial.
          if(spec%params(ix)%sampled) cycle

          ! Fixed parameters must match the scalar values used by the search.
          tolerance=1000._rkind*epsilon(1._rkind)* &
                    max(1._rkind,abs(spec%params(ix)%trial_value))

          if(any(abs(field(ix)%initial-spec%params(ix)%trial_value)>tolerance))then
            message=trim(message)//'spatially varying fixed parameter requires field-based '// &
                    'ordered constraints: '//trim(spec%params(ix)%name)
            err=20; return
          endif

        enddo

      enddo

    endif
    
    ! -----------------------------------------------------------------------------------------------
    ! Verify the spatial alignment of every parameter ratio.
    ! -----------------------------------------------------------------------------------------------
    do i=1,size(field)

      if(.not.spec%params(i)%sampled) cycle
      if(spec%params(i)%param_type/=param_ratio) cycle

      ixReference=0

      do j=1,size(spec%params)
        if(trim(spec%params(j)%name)==trim(spec%params(i)%reference_name))then
          ixReference=j
          exit
        endif
      enddo

      if(ixReference==0)then
        message=trim(message)//'unknown ratio reference for '//trim(spec%params(i)%name)
        err=20; return
      endif

      if(size(field(i)%initial)/=size(field(ixReference)%initial))then
        message=trim(message)//'different field lengths in ratio for '//trim(spec%params(i)%name)
        err=20; return
      endif

      ! Basin and model parameters have different spatial indexing.
      if((param_map(i)%ixB>0) .neqv. (param_map(ixReference)%ixB>0))then
        message=trim(message)//'ratio cannot mix basin and model parameters: '//trim(spec%params(i)%name)
        err=20; return
      endif

      ! Both basin fields contain one value per GRU.
      if(param_map(i)%ixB>0) cycle

      ! Model fields must select the same domains.
      if(param_map(i)%domain_type/=param_map(ixReference)%domain_type)then
        message=trim(message)//'ratio fields have different domain qualifiers: '//trim(spec%params(i)%name)
        err=20; return
      endif

      ! Verify that array elements align within every selected domain.
      do g=1,size(mpar_init%gru)
        do h=1,size(mpar_init%gru(g)%hru)
          do d=1,size(mpar_init%gru(g)%hru(h)%dom)

            if(param_map(i)%domain_type/=0 .and. &
               param_map(i)%domain_type/=gru_struc(g)%hruInfo(h)%domInfo(d)%dom_type) cycle

            if(size(mpar_init%gru(g)%hru(h)%dom(d)%var(param_map(i)%ixM)%dat)/= &
               size(mpar_init%gru(g)%hru(h)%dom(d)%var(param_map(ixReference)%ixM)%dat))then
              message=trim(message)//'ratio fields are not aligned for '//trim(spec%params(i)%name)
              err=20; return
            endif

          enddo
        enddo
      enddo

    enddo

  end subroutine initialize_summa_parameter_fields


  ! **************************************************************************************************
  ! Extract one SUMMA parameter field into a one-dimensional vector.
  !
  ! Model parameters are traversed in GRU, HRU, domain, and array-element
  ! order. Basin parameters contain one value per GRU.
  !
  ! The parameter map identifies the SUMMA parameter and any domain
  ! restriction. The same traversal order is used when writing values back.
  ! **************************************************************************************************
  subroutine get_summa_parameter_vector(map,mpar,bpar,values)
    implicit none

    type(summa_parameter_map), intent(in) :: map          ! Location of the parameter in SUMMA
    type(gru_hru_dom_doubleVec), intent(in) :: mpar        ! Model parameter fields
    type(gru_double), intent(in) :: bpar                 ! Basin parameter fields
    real(rkind), allocatable, intent(out) :: values(:)   ! Extracted parameter values

    integer(i4b) :: g                                    ! GRU index
    integer(i4b) :: h                                    ! HRU index
    integer(i4b) :: d                                    ! Domain index
    integer(i4b) :: n                                    ! Total number of values
    integer(i4b) :: k                                    ! Current vector offset
    integer(i4b) :: nLocal                               ! Number of values in current domain

    ! Basin parameters: one value per GRU.
    if(map%ixB > 0)then
      allocate(values(size(bpar%gru)))

      do g=1,size(bpar%gru)
        values(g)=bpar%gru(g)%var(map%ixB)
      enddo

      return
    endif

    ! Model parameters: first count all applicable elements.
    n=0

    do g=1,size(mpar%gru)
      do h=1,size(mpar%gru(g)%hru)
        do d=1,size(mpar%gru(g)%hru(h)%dom)

          ! Skip domains outside the parameter's domain qualifier.
          if(map%domain_type/=0 .and. &
             map%domain_type/=gru_struc(g)%hruInfo(h)%domInfo(d)%dom_type) cycle

          n=n+size(mpar%gru(g)%hru(h)%dom(d)%var(map%ixM)%dat)

        enddo
      enddo
    enddo

    allocate(values(n))

    ! Extract the values in the same order used to count them.
    k=0

    do g=1,size(mpar%gru)
      do h=1,size(mpar%gru(g)%hru)
        do d=1,size(mpar%gru(g)%hru(h)%dom)

          if(map%domain_type/=0 .and. &
             map%domain_type/=gru_struc(g)%hruInfo(h)%domInfo(d)%dom_type) cycle

          nLocal=size(mpar%gru(g)%hru(h)%dom(d)%var(map%ixM)%dat)

          values(k+1:k+nLocal)=mpar%gru(g)%hru(h)%dom(d)%var(map%ixM)%dat

          k=k+nLocal

        enddo
      enddo
    enddo

  end subroutine get_summa_parameter_vector


  ! **************************************************************************************************
  ! Assign a one-dimensional parameter vector to the SUMMA parameter fields.
  !
  ! Uses the same GRU, HRU, domain, and array-element ordering as the
  ! extraction routine. Only the selected parameter is modified.
  !
  ! The caller must supply a vector with the correct number of elements.
  ! **************************************************************************************************
  subroutine put_summa_parameter_vector(map,values,mpar,bpar)
    implicit none

    type(summa_parameter_map), intent(in) :: map          ! Location of the parameter in SUMMA
    real(rkind), intent(in) :: values(:)                 ! Physical parameter values to assign
    type(gru_hru_dom_doubleVec), intent(inout) :: mpar    ! Model parameter fields to update
    type(gru_double), intent(inout) :: bpar               ! Basin parameter fields to update

    integer(i4b) :: g                                    ! GRU index
    integer(i4b) :: h                                    ! HRU index
    integer(i4b) :: d                                    ! Domain index
    integer(i4b) :: k                                    ! Current vector offset
    integer(i4b) :: nLocal                               ! Number of values in current domain

    ! Basin parameters: assign one value to each GRU.
    if(map%ixB > 0)then

      do g=1,size(bpar%gru)
        bpar%gru(g)%var(map%ixB)=values(g)
      enddo

      return
    endif

    ! Model parameters: reconstruct the nested spatial fields.
    k=0

    do g=1,size(mpar%gru)
      do h=1,size(mpar%gru(g)%hru)
        do d=1,size(mpar%gru(g)%hru(h)%dom)

          if(map%domain_type/=0 .and. &
             map%domain_type/=gru_struc(g)%hruInfo(h)%domInfo(d)%dom_type) cycle

          nLocal=size(mpar%gru(g)%hru(h)%dom(d)%var(map%ixM)%dat)

          mpar%gru(g)%hru(h)%dom(d)%var(map%ixM)%dat=values(k+1:k+nLocal)

          k=k+nLocal

        enddo
      enddo
    enddo

  end subroutine put_summa_parameter_vector

  
  ! **************************************************************************************************
  ! Apply one calibration sample to the SUMMA parameter fields.
  !
  ! The scalar calibration values have already been transformed from search
  ! coordinates. The generic parameter_fields module constructs the effective
  ! physical fields using the appropriate parameterization.
  !
  ! Only sampled parameters are written to SUMMA. Non-sampled parameters retain
  ! their initialized physical fields, but may serve as ratio references.
  !
  ! The destination parameter structures must have the same spatial layout
  ! as the initialized structures used to construct field(:).
  ! **************************************************************************************************
  subroutine apply_summa_parameter_sample(spec,search,sampled_values,param_map,field,mpar,bpar,err,message)
    USE, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    type(parameter_spec), intent(in) :: spec
    type(parameter_search_info), intent(in) :: search
    real(rkind), intent(in) :: sampled_values(:)
    type(summa_parameter_map), intent(in) :: param_map(:)
    type(parameter_field), intent(inout) :: field(:)
    type(gru_hru_dom_doubleVec), intent(inout) :: mpar
    type(gru_double), intent(inout) :: bpar
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message

    integer(i4b) :: i
    real(rkind), allocatable :: destination(:)
    character(len=1024) :: cmessage

    err=0
    message='apply_summa_parameter_sample/'

    if(.not.allocated(spec%params))then
      message=trim(message)//'parameter specification is not allocated'
      err=20; return
    endif

    if(size(param_map)/=size(spec%params) .or. size(field)/=size(spec%params))then
      message=trim(message)//'parameter map or fields have the wrong size'
      err=20; return
    endif

    ! Construct the effective physical fields using the generic parameterization routines.
    call construct_effective_parameter_fields(spec,search,sampled_values,field,err,cmessage)

    if(err/=0)then
      message=trim(message)//trim(cmessage)
      return
    endif

    ! -----------------------------------------------------------------------------------------------
    ! Validate the physical fields and destination sizes before changing SUMMA.
    ! -----------------------------------------------------------------------------------------------
    do i=1,size(field)

      if(.not.spec%params(i)%sampled) cycle

      if(.not.all(ieee_is_finite(field(i)%effective)))then
        message=trim(message)//'non-finite effective field for '//trim(spec%params(i)%name)
        err=20; return
      endif

      call get_summa_parameter_vector(param_map(i),mpar,bpar,destination)

      if(size(destination)/=size(field(i)%effective))then
        message=trim(message)//'SUMMA destination size changed for '//trim(spec%params(i)%name)
        err=20; return
      endif

    enddo

    ! -----------------------------------------------------------------------------------------------
    ! Write sampled physical parameter fields back into SUMMA.
    ! -----------------------------------------------------------------------------------------------
    do i=1,size(field)

      if(.not.spec%params(i)%sampled) cycle

      call put_summa_parameter_vector(param_map(i),field(i)%effective,mpar,bpar)

    enddo

  end subroutine apply_summa_parameter_sample

end module summa_parameter_fields

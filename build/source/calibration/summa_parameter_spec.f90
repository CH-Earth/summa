module summa_parameter_spec

  USE nr_type,          only: i4b, rkind, lgt
  USE data_types,       only: var_d
  USE data_types,       only: var_dlength
  USE summa_type,       only: config_info
  USE globalData,       only: iulog

  USE parameter_search, only: parameter_spec
  USE parameter_search, only: param_physical,param_spatial_multiplier,param_ratio

  implicit none
  private

  public :: get_summa_parameter_spec
  public :: build_summa_parameter_overrides

contains

  ! **************************************************************************************************
  ! Construct the SUMMA parameter specification.
  !
  ! Populates the model-agnostic parameter specification used by the parameter-search routines.
  !
  ! Parameters explicitly included in the calibration list are marked as sampled. Additional
  ! parameters required by ordered dependencies are included as non-sampled parameters, except in a
  ! chain searched by its gaps, whose every member is sampled.
  !
  ! All parameters represented by this interface are currently assumed to be spatially uniform.
  ! Initial trial values are obtained from the initialized SUMMA parameter structures rather than
  ! the fallback defaults when the run has one GRU, one HRU and one domain, and from the fallback
  ! defaults otherwise. For sampled parameters, these values are replaced during calibration.
  ! For non-sampled parameters, the initialized values are retained when evaluating constraints
  ! and constructing parameter overrides.
  !
  ! Spatially varying calibration parameters and spatial regularization are not currently supported.
  ! **************************************************************************************************
  subroutine get_summa_parameter_spec(config,spec,param_map,mpar_init,bpar_init,err,message)
    USE summa_parameter_fields, only: summa_parameter_map
    USE get_ixname_module, only: get_ixParam,get_ixBpar
    USE parameter_overrides_module, only: split_domain_target
    implicit none
    ! dummy
    type(config_info)                     , intent(in)  :: config
    type(parameter_spec)                  , intent(out) :: spec
    type(var_dlength),            optional, intent(in)  :: mpar_init     ! initialized HRU parameters, single-domain runs only
    type(var_d),                  optional, intent(in)  :: bpar_init     ! initialized basin parameters, single-domain runs only
    type(summa_parameter_map), allocatable, intent(out) :: param_map(:)  ! Resolved SUMMA parameter locations
    integer(i4b)                          , intent(out) :: err
    character(*)                          , intent(out) :: message
    ! local
    character(len=64), allocatable :: param_names(:)
    logical(lgt),      allocatable :: sampled(:)
    integer(i4b) :: i
    integer(i4b) :: j
    integer(i4b) :: ix
    integer(i4b) :: iAt  ! position of @ in a domain-qualified parameter name
    integer(i4b) :: nParam
    integer(i4b) :: nMax
    integer(i4b) :: nConstraint
    integer(i4b) :: nSumma
    integer(i4b) :: nMultiplier
    integer(i4b) :: ixReference, ixDependent
    real(rkind)  :: lower,upper
    real(rkind)  :: original_value
    logical(lgt) :: isMultiplier
    logical(lgt) :: hasMultiplierBounds
    logical(lgt) :: isBasin                                 ! True if the ordered constraint contains basin parameters
    character(len=64) :: base_name                          ! Parameter name without its domain qualifier
    character(len=1024) :: cmessage

    err = 0
    message = 'get_summa_parameter_spec/'

    ! calibration parameter list is required
    if(.not.allocated(config%calib%param_list))then
      message=trim(message)//'calibration parameter list is not defined'
      err=20; return
    endif
    nSumma = size(config%calib%param_list)
    if(nSumma == 0)then
      message=trim(message)//'calibration parameter list is empty'
      err=20; return
    endif

    ! -----------------------------------------------------------------------------------------------
    ! Determine maximum possible number of unique parameters.
    ! -----------------------------------------------------------------------------------------------
    nMax = nSumma
    if(allocated(config%calib%ordered))then
      do i=1,size(config%calib%ordered)
        if(.not.allocated(config%calib%ordered(i)%parameters))then
          write(message,'(A,I0)') trim(message)// 'parameter list not defined for ordered constraint = ',i
          err=20; return
        endif
        nMax = nMax + size(config%calib%ordered(i)%parameters)
      enddo
    endif
    ! Include parameters referenced by parameter ratios.
    if(allocated(config%calib%multiplier))then
      nMax = nMax + size(config%calib%multiplier)
    endif
    allocate(param_names(nMax),sampled(nMax),stat=err)
    if(err/=0)then
      message=trim(message)//'unable to allocate SUMMA parameter registry'
      return
    endif
    param_names = ''
    sampled     = .false.
    nParam      = 0

    ! -----------------------------------------------------------------------------------------------
    ! Add sampled calibration parameters.
    ! -----------------------------------------------------------------------------------------------
    do i=1,nSumma
      ix = find_parameter(trim(config%calib%param_list(i)%name), param_names(1:nParam))
      if(ix > 0)then
        message=trim(message)//'duplicate calibration parameter: '// trim(config%calib%param_list(i)%name)
        err=20; return
      endif
      nParam = nParam + 1
      param_names(nParam) = trim(config%calib%param_list(i)%name)
      sampled(nParam)     = .true.
    enddo
    ! -----------------------------------------------------------------------------------------------
    ! Add parameters required only by ordered dependencies.
    ! -----------------------------------------------------------------------------------------------
    if(allocated(config%calib%ordered))then
      do i=1,size(config%calib%ordered)
        do j=1,size(config%calib%ordered(i)%parameters)
          ix = find_parameter(trim(config%calib%ordered(i)%parameters(j)), param_names(1:nParam))
          if(ix == 0)then
            nParam = nParam + 1
            param_names(nParam) = trim(config%calib%ordered(i)%parameters(j))
            sampled(nParam)     = .false.
          endif
        enddo
      enddo
    endif
    ! -----------------------------------------------------------------------------------------------
    ! Add parameters required only as references for parameter ratios.
    ! -----------------------------------------------------------------------------------------------
    if(allocated(config%calib%multiplier))then
      do i=1,size(config%calib%multiplier)
    
        ix=find_parameter(trim(config%calib%multiplier(i)%parameters(1)),param_names(1:nParam))
    
        if(ix==0)then
          nParam=nParam+1
          param_names(nParam)=trim(config%calib%multiplier(i)%parameters(1))
          sampled(nParam)=.false.
        endif
    
      enddo
    endif
    ! -----------------------------------------------------------------------------------------------
    ! Reject combinations of global and domain-specific calibration parameters.
    !
    ! An unqualified parameter applies to all domains, so combining it with a domain-specific
    ! version would make the result depend on the order in which parameter overrides are applied.
    ! -----------------------------------------------------------------------------------------------
    do i=1,nParam
      iAt=index(trim(param_names(i)),'@')
      if(iAt == 0) cycle

      ix=find_parameter(param_names(i)(1:iAt-1),param_names(1:nParam))
      if(ix > 0)then
        message=trim(message)//'cannot combine global and domain-specific parameters: "'// &
                trim(param_names(ix))//'" and "'//trim(param_names(i))//'"'
        err=20; return
      endif
    enddo

    ! -----------------------------------------------------------------------------------------------
    ! Populate master parameter registry from SUMMA metadata and initialized parameter values.
    !
    ! trial_value is initialized from the actual SUMMA parameter values. For sampled parameters,
    ! this value is replaced during calibration. For non-sampled constraint parameters, the
    ! initialized value is retained and applied across the model domain for every trial.
    ! -----------------------------------------------------------------------------------------------
    allocate(spec%params(nParam),stat=err)
    if(err/=0)then
      message=trim(message)//'unable to allocate parameter specification'
      return
    endif
    do i=1,nParam
      spec%params(i)%name    = trim(param_names(i))
      spec%params(i)%sampled = sampled(i)
      spec%params(i)%transformation = 'none'
      call get_summa_parameter_info(trim(param_names(i)),     &
                                    spec%params(i)%trial_value, &
                                    spec%params(i)%lower,       &
                                    spec%params(i)%upper,       &
                                    spec%params(i)%units,       &
                                    spec%params(i)%long_name,   &
                                    mpar_init,bpar_init,err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    enddo

    ! Copy the parameterization metadata for the sampled calibration parameters.
    do i=1,nSumma
      spec%params(i)%param_type=config%calib%param_list(i)%param_type
      spec%params(i)%reference_name=config%calib%param_list(i)%reference_name
      spec%params(i)%reference_reduction=config%calib%param_list(i)%reference_reduction
    enddo

    ! -----------------------------------------------------------------------------------------------
    ! Resolve calibration parameters to SUMMA model and basin parameter indices.
    ! -----------------------------------------------------------------------------------------------
    allocate(param_map(nParam))

    do i=1,nParam
      call split_domain_target(spec%params(i)%name,base_name,param_map(i)%domain_type)
      param_map(i)%ixM=get_ixParam(trim(base_name))
      param_map(i)%ixB=get_ixBpar(trim(base_name))
    enddo

    ! -----------------------------------------------------------------------------------------------
    ! Replace the parameter files' bounds with configured ones; a default outside them moves inside
    ! -----------------------------------------------------------------------------------------------
    if(allocated(config%calib%param_bounds))then
      do i=1,size(config%calib%param_bounds)
        ix = find_parameter(trim(config%calib%param_bounds(i)%name), param_names(1:nParam))
        if(ix == 0)then
          message=trim(message)//'bounds given for a parameter that is not calibrated: '// &
                  trim(config%calib%param_bounds(i)%name)
          err=20; return
        endif
        if(config%calib%param_bounds(i)%lower > config%calib%param_bounds(i)%upper)then
          message=trim(message)//'invalid bounds for parameter: '//trim(config%calib%param_bounds(i)%name)
          err=20; return
        endif
        ! Apply configured bounds and clamp the initial parameter value if necessary
        original_value = spec%params(ix)%trial_value
        spec%params(ix)%lower = config%calib%param_bounds(i)%lower
        spec%params(ix)%upper = config%calib%param_bounds(i)%upper
        spec%params(ix)%trial_value = min(max(original_value, spec%params(ix)%lower), spec%params(ix)%upper)
        ! Warn if clamping changes a parameter that is not sampled
        if(.not.spec%params(ix)%sampled)then
          if(original_value < spec%params(ix)%lower .or. &
             original_value > spec%params(ix)%upper)then
            if(config%is_case_root)then
              write(iulog,'(3A,ES12.4,A,ES12.4)') &
                'WARNING: non-sampled parameter "',trim(spec%params(ix)%name), &
                '" clamped from ',original_value,' to ',spec%params(ix)%trial_value
            endif
          endif
        endif
      enddo
    endif

    ! -----------------------------------------------------------------------------------------------
    ! Apply configured calibration parameter bounds and transformations.
    !
    ! Calibration bounds specified in the TOML replace the native SUMMA search bounds
    ! and may extend beyond them. Emit a warning when the configured range exceeds
    ! the native SUMMA bounds.
    ! -----------------------------------------------------------------------------------------------
    if(allocated(config%calib%param_list))then

      do i=1,size(config%calib%param_list)

        ix = find_parameter(trim(config%calib%param_list(i)%name), &
                            param_names(1:nParam))

        if(ix == 0)then
          message=trim(message)//'calibration parameter not found: '// &
                  trim(config%calib%param_list(i)%name)
          err=20
          return
        endif

        ! Check whether this parameter is sampled as a multiplier
        isMultiplier=(spec%params(ix)%param_type /= param_physical)

        ! A multiplier searches a ratio, not the native physical parameter.
        ! Require explicit search bounds, either inline or in parameter_bounds.
        if(isMultiplier)then
          if(.not.config%calib%param_list(i)%has_lower .or. &
             .not.config%calib%param_list(i)%has_upper)then
            message=trim(message)//'spatial multiplier or parameter ratio requires explicit bounds: '// &
                    trim(config%calib%param_list(i)%name)
            err=20; return
          endif
        endif


        ! Bounds not given in the configuration retain the existing bounds
        lower = spec%params(ix)%lower
        upper = spec%params(ix)%upper

        if(config%calib%param_list(i)%has_lower) &
          lower = config%calib%param_list(i)%lower_bound

        if(config%calib%param_list(i)%has_upper) &
          upper = config%calib%param_list(i)%upper_bound

        if(lower >= upper)then
          write(message,'(A,ES12.4,1X,ES12.4)') trim(message)// &
            'lower bound must be less than upper bound for calibration parameter '// &
            trim(config%calib%param_list(i)%name)//': ',lower,upper
          err=20
          return
        endif

        ! Compare physical bounds with native SUMMA bounds.
        ! Multiplier bounds describe search coordinates, not physical values.
        if(.not.isMultiplier)then
          if(lower < spec%params(ix)%lower .or. upper > spec%params(ix)%upper)then

            write(iulog,'(A,A,A,2(ES12.4,1X),A,2(ES12.4,1X))') &
              'WARNING: calibration bounds exceed native SUMMA bounds for ', &
              trim(config%calib%param_list(i)%name),                         &
              ': calibration = ',                                            &
              lower,upper,                                                   &
              ', native = ',                                                 &
              spec%params(ix)%lower,                                         &
              spec%params(ix)%upper

          endif
        endif

        ! apply calibration-specific bounds
        spec%params(ix)%lower = lower
        spec%params(ix)%upper = upper

        ! apply calibration-specific transformation
        spec%params(ix)%transformation = &
          trim(config%calib%param_list(i)%transformation)

      enddo

    endif

    ! -----------------------------------------------------------------------------------------------
    ! Convert named SUMMA constraints to master-registry indices.
    ! -----------------------------------------------------------------------------------------------
    nConstraint = 0
    if(allocated(config%calib%ordered)) nConstraint = size(config%calib%ordered)
    
    allocate(spec%ordered(nConstraint),stat=err)
    if(err/=0)then
      message=trim(message)//'unable to allocate ordered parameter specifications'
      return
    endif
    do i=1,nConstraint
      allocate(spec%ordered(i)%param_index(size(config%calib%ordered(i)%parameters)),stat=err)
      if(err/=0)then
        message=trim(message)//'unable to allocate ordered parameter indices'
        return
      endif
      do j=1,size(config%calib%ordered(i)%parameters)
        ix = find_parameter(trim(config%calib%ordered(i)%parameters(j)), param_names(1:nParam))
        if(ix == 0)then
          message=trim(message)// &
                  'unable to resolve parameter in ordered constraint: '// &
                  trim(config%calib%ordered(i)%parameters(j))
          err=20; return
        endif
        spec%ordered(i)%param_index(j) = ix
      enddo
      spec%ordered(i)%gap_fraction = config%calib%ordered(i)%gap_fraction
    enddo

    ! -----------------------------------------------------------------------------------------------
    ! Ensure ordered constraints do not mix model and basin parameters.
    ! -----------------------------------------------------------------------------------------------
    do i=1,size(spec%ordered)
   
      ix=spec%ordered(i)%param_index(1)
      isBasin=(param_map(ix)%ixB > 0)
   
      do j=2,size(spec%ordered(i)%param_index)
   
        ix=spec%ordered(i)%param_index(j)
   
        if((param_map(ix)%ixB > 0) .neqv. isBasin)then
          message=trim(message)//'ordered constraint cannot mix mpar and bpar parameters'
          err=20; return
        endif
   
      enddo
    enddo

    ! -----------------------------------------------------------------------------------------------
    ! Convert named multiplier dependencies to master-registry indices.
    ! -----------------------------------------------------------------------------------------------

    nMultiplier = 0
    if(allocated(config%calib%multiplier)) nMultiplier = size(config%calib%multiplier)

    allocate(spec%multiplier(nMultiplier),stat=err)
    if(err/=0)then
      message=trim(message)//'unable to allocate multiplier parameter specifications'
      return
    endif

    do i=1,nMultiplier

      ! Resolve the reference and dependent parameters.
      ixReference = find_parameter( &
        trim(config%calib%multiplier(i)%parameters(1)), param_names(1:nParam))

      ixDependent = find_parameter( &
        trim(config%calib%multiplier(i)%parameters(2)), param_names(1:nParam))

      if(ixReference == 0 .or. ixDependent == 0)then
        message=trim(message)//'multiplier parameter not found in parameter registry'
        err=20; return
      endif

      ! A dependent parameter can appear in only one multiplier pair.
      do j=1,i-1
        if(spec%multiplier(j)%dependent_index == ixDependent)then
          message=trim(message)//'duplicate multiplier dependent: '// &
                  trim(param_names(ixDependent))
          err=20; return
        endif
      enddo

      ! Ratio dependents may participate in ordered constraints. The generic
      ! checker tests the effective physical fields after ratio construction.

      ! Store the parameter indices.
      spec%multiplier(i)%reference_index = ixReference
      spec%multiplier(i)%dependent_index = ixDependent

    enddo

    ! Chained ratios are resolved in dependency order by initialize_parameter_fields.
    ! Its cycle check runs during adapter setup, before any trial is dispatched.

    ! -----------------------------------------------------------------------------------------------
    ! Log multiplier dependencies 
    ! -----------------------------------------------------------------------------------------------
    write(iulog,'(A,I0)') 'Number of multiplier dependencies: ',nMultiplier

    do i=1,nMultiplier
      ixReference = spec%multiplier(i)%reference_index
      ixDependent = spec%multiplier(i)%dependent_index

      write(iulog,'(A,A,A,A)') &
        '  reference = ',trim(spec%params(ixReference)%name), &
        ', dependent = ',trim(spec%params(ixDependent)%name)
    enddo

  end subroutine get_summa_parameter_spec

  ! **************************************************************************************************
  ! Build the complete scalar SUMMA parameter override vector for one trial.
  !
  ! Sampled parameters receive their values from sampled_values. Non-sampled parameters represented
  ! in the parameter specification receive their scalar trial_value. A chain searched by its gaps is
  ! then turned from fractions into values.
  !
  ! The ordering of param_values matches spec%params and therefore the invariant parameter-name
  ! vector constructed from spec%params(:)%name during parameter-evaluation initialization.
  !
  ! Because summa_paramSetup applies caller-supplied overrides after reading the trial parameter
  ! file, these values ensure that all parameters participating in calibration constraints are
  ! spatially uniform for the current implementation.
  ! **************************************************************************************************
  subroutine build_summa_parameter_overrides(spec, sampled_names, sampled_values, param_values, err, message)
    implicit none
    type(parameter_spec), intent(in) :: spec
    character(*), intent(in)  :: sampled_names(:)
    real(rkind),  intent(in)  :: sampled_values(:)
    real(rkind),  intent(out) :: param_values(:)
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    integer(i4b) :: i
    integer(i4b) :: ixSample
    integer(i4b) :: ixReference
    integer(i4b) :: ixDependent

    err = 0
    message = 'build_summa_parameter_overrides/'

    ! sampled name-value vectors must be consistent
    if(size(sampled_names) /= size(sampled_values))then
      message=trim(message)// 'sampled parameter name and value vectors have different sizes'
      err=20; return
    endif

    ! override vector must be consistent with the param specs
    if(size(param_values) /= size(spec%params))then
      message=trim(message)//'parameter override vector has incorrect size'
      err=20; return
    endif
    do i=1,size(spec%params)
      if(spec%params(i)%sampled)then
        ixSample = find_parameter(trim(spec%params(i)%name),sampled_names)
        if(ixSample == 0)then
          message=trim(message)// 'sampled value not provided for parameter: '// trim(spec%params(i)%name)
          err=20; return
        endif
        param_values(i) = sampled_values(ixSample)

      else

        ! Constraint-only parameters are reset to the scalar SUMMA default
        ! stored as trial_value in the model parameter specification.
        param_values(i) = spec%params(i)%trial_value
      endif
    enddo

    ! Apply multiplicative parameter dependencies.
    do i=1,size(spec%multiplier)
    
      ixReference = spec%multiplier(i)%reference_index
      ixDependent = spec%multiplier(i)%dependent_index
    
      param_values(ixDependent) = &
        param_values(ixReference) * param_values(ixDependent)
    
    enddo

  end subroutine build_summa_parameter_overrides

  ! **************************************************************************************************
  ! Get SUMMA parameter metadata.
  !
  ! Returns the SUMMA default value and native lower and upper bounds for a local or basin parameter.
  !
  ! The returned default value is used as the scalar trial value for parameters that participate in
  ! calibration constraints but are not sampled. The scalar override is subsequently applied across
  ! all spatial elements by summa_paramSetup.
  ! **************************************************************************************************
  subroutine get_summa_parameter_info(param_name,trial_value,lower,upper,units,long_name, &
                                      mpar_init,bpar_init,err,message)
    USE get_ixname_module, only: get_ixParam,get_ixBpar
    USE globalData, only: localParFallback,basinParFallback
    USE globalData, only: mpar_meta,bpar_meta
    USE globalData, only: integerMissing
    USE parameter_overrides_module, only: split_domain_target
    implicit none
    ! dummy
    character(*)               , intent(in)  :: param_name
    real(rkind)                , intent(out) :: trial_value
    real(rkind)                , intent(out) :: lower
    real(rkind)                , intent(out) :: upper
    character(*)               , intent(out) :: units
    character(*)               , intent(out) :: long_name
    type(var_dlength), optional, intent(in)  :: mpar_init  ! initialized HRU parameters, single-domain runs only
    type(var_d),       optional, intent(in)  :: bpar_init  ! initialized basin parameters, single-domain runs only
    integer(i4b)               , intent(out) :: err
    character(*)               , intent(out) :: message
    ! local
    integer(i4b) :: ixParam,ixBasin
    character(len=64) :: name     ! parameter name without its domain type
    integer(i4b)      :: domType  ! domain type a name@domainType targets, 0 for every domain
  
    err=0
    message='get_summa_parameter_info/'
    ! a name@domainType parameter takes the default and bounds of the parameter it names
    call split_domain_target(param_name,name,domType)
    if(domType==integerMissing)then
      message=trim(message)//'unknown domain type in "'//trim(param_name)//'" (use upland, glacCln1, glacCln2, glacDbr, wetland or stream)'
      err=20; return
    endif
  
    ! local parameter
    ixParam=get_ixParam(trim(name))
    if(ixParam > 0)then
      if(present(mpar_init))then
        trial_value=mpar_init%var(ixParam)%dat(1)
      else
        trial_value=localParFallback(ixParam)%default_val
      endif
      lower      =localParFallback(ixParam)%lower_limit
      upper      =localParFallback(ixParam)%upper_limit
      units      =trim(mpar_meta(ixParam)%varUnit)
      long_name  =trim(mpar_meta(ixParam)%varDesc)
      if(lower > upper)then
        message=trim(message)//'invalid bounds for parameter: '//trim(param_name)
        err=20; return
      endif
      return
    endif

    ! basin parameter
    ixBasin=get_ixBpar(trim(name))
    if(ixBasin > 0)then
      if(domType>0)then
        message=trim(message)//'basin parameter "'//trim(name)//'" has no domain type to target'
        err=20; return
      endif
      if(present(bpar_init))then
        trial_value=bpar_init%var(ixBasin)
      else
        trial_value=basinParFallback(ixBasin)%default_val
      endif 
      lower      =basinParFallback(ixBasin)%lower_limit
      upper      =basinParFallback(ixBasin)%upper_limit
      units      =trim(bpar_meta(ixBasin)%varUnit)
      long_name  =trim(bpar_meta(ixBasin)%varDesc)
      if(lower > upper)then
        message=trim(message)//'invalid bounds for parameter: '//trim(param_name)
        err=20; return
      endif
      return
    endif
    message=trim(message)//'parameter not found: '//trim(param_name)
    err=20; return

  end subroutine get_summa_parameter_info

  ! **************************************************************************************************
  ! Find a parameter in a parameter-name vector.
  !
  ! Returns the one-based parameter index, or zero if the parameter is not present.
  ! **************************************************************************************************
  integer(i4b) function find_parameter(param_name,param_names)
    implicit none
    character(*), intent(in) :: param_name
    character(*), intent(in) :: param_names(:)
    integer(i4b) :: i

    find_parameter = 0
    do i=1,size(param_names)
      if(trim(param_names(i)) == trim(param_name))then
        find_parameter = i
        return
      endif
    enddo

  end function find_parameter

end module summa_parameter_spec

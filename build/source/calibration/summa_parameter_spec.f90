module summa_parameter_spec

  USE nr_type,          only: i4b, rkind, lgt
  USE summa_type,       only: config_info

  USE parameter_search, only: parameter_spec

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
  ! parameters required by ordered dependencies are included as non-sampled parameters.
  !
  ! All parameters represented by this interface are currently assumed to be spatially uniform.
  ! For non-sampled parameters participating in constraints, trial_value is set to the SUMMA default
  ! value. build_summa_parameter_overrides() subsequently ensures that this scalar value overrides
  ! any spatially varying values read from a trial parameter file.
  !
  ! Spatially varying calibration parameters and spatial regularization are not currently supported.
  ! **************************************************************************************************

  subroutine get_summa_parameter_spec(config, spec, err, message)

    implicit none

    type(config_info),    intent(in)  :: config
    type(parameter_spec), intent(out) :: spec
    integer(i4b),         intent(out) :: err
    character(*),         intent(out) :: message

    character(len=64), allocatable :: param_names(:)
    logical(lgt),      allocatable :: sampled(:)

    integer(i4b) :: i
    integer(i4b) :: j
    integer(i4b) :: ix
    integer(i4b) :: nParam
    integer(i4b) :: nMax
    integer(i4b) :: nConstraint

    character(len=256) :: cmessage

    err = 0
    message = 'get_summa_parameter_spec/'

    ! calibration parameter list is required
    if(.not.allocated(config%calib%param_list))then
      message=trim(message)//'calibration parameter list is not defined'
      err=20; return
    endif

    if(size(config%calib%param_list) == 0)then
      message=trim(message)//'calibration parameter list is empty'
      err=20; return
    endif


    ! -----------------------------------------------------------------------------------------------
    ! Determine maximum possible number of unique parameters.
    ! -----------------------------------------------------------------------------------------------

    nMax = size(config%calib%param_list)

    if(allocated(config%calib%ordered))then

      do i=1,size(config%calib%ordered)

        if(.not.allocated(config%calib%ordered(i)%parameters))then
          write(message,'(A,I0)') trim(message)// &
            'parameter list not defined for ordered constraint = ',i
          err=20; return
        endif

        nMax = nMax + size(config%calib%ordered(i)%parameters)

      enddo

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

    do i=1,size(config%calib%param_list)

      ix = find_parameter(trim(config%calib%param_list(i)), &
                          param_names(1:nParam))

      if(ix > 0)then
        message=trim(message)//'duplicate calibration parameter: '// &
                trim(config%calib%param_list(i))
        err=20; return
      endif

      nParam = nParam + 1

      param_names(nParam) = trim(config%calib%param_list(i))
      sampled(nParam)     = .true.

    enddo


    ! -----------------------------------------------------------------------------------------------
    ! Add parameters required only by ordered dependencies.
    ! -----------------------------------------------------------------------------------------------

    if(allocated(config%calib%ordered))then

      do i=1,size(config%calib%ordered)

        do j=1,size(config%calib%ordered(i)%parameters)

          ix = find_parameter(trim(config%calib%ordered(i)%parameters(j)), &
                              param_names(1:nParam))

          if(ix == 0)then

            nParam = nParam + 1

            param_names(nParam) = trim(config%calib%ordered(i)%parameters(j))
            sampled(nParam)     = .false.

          endif

        enddo

      enddo

    endif


    ! -----------------------------------------------------------------------------------------------
    ! Populate master parameter registry from SUMMA metadata.
    !
    ! trial_value is initialized from the SUMMA default. For sampled parameters this value is not
    ! used during sampling. For non-sampled constraint parameters it becomes the scalar value that
    ! will be applied across the model domain for every trial.
    ! -----------------------------------------------------------------------------------------------

    allocate(spec%params(nParam),stat=err)

    if(err/=0)then
      message=trim(message)//'unable to allocate parameter specification'
      return
    endif

    do i=1,nParam

      spec%params(i)%name    = trim(param_names(i))
      spec%params(i)%sampled = sampled(i)

      call get_summa_parameter_info(trim(param_names(i)),                &
                                    spec%params(i)%trial_value,          &
                                    spec%params(i)%lower,                &
                                    spec%params(i)%upper,                &
                                    spec%params(i)%units,                &
                                    spec%params(i)%long_name,            &
                                    err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! default parameter transformation
      spec%params(i)%transformation = 'none'

    enddo


    ! -----------------------------------------------------------------------------------------------
    ! Apply configured parameter transformations.
    !
    ! Transformations are meaningful only for sampled parameters.
    ! -----------------------------------------------------------------------------------------------

    if(allocated(config%calib%param_transform))then

      do i=1,size(config%calib%param_transform)

        ix = find_parameter(                               &
               trim(config%calib%param_transform(i)%name), &
               param_names(1:nParam))

        if(ix == 0)then
          message=trim(message)// &
                  'parameter transformation defined for unknown parameter: '// &
                  trim(config%calib%param_transform(i)%name)
          err=20; return
        endif

        if(.not.spec%params(ix)%sampled)then
          message=trim(message)// &
                  'parameter transformation defined for non-sampled parameter: '// &
                  trim(config%calib%param_transform(i)%name)
          err=20; return
        endif

        spec%params(ix)%transformation = &
          trim(config%calib%param_transform(i)%transformation)

      enddo

    endif


    ! -----------------------------------------------------------------------------------------------
    ! Convert named SUMMA constraints to master-registry indices.
    ! -----------------------------------------------------------------------------------------------

    if(.not.allocated(config%calib%ordered))then
      allocate(spec%ordered(0))
      return
    endif

    nConstraint = size(config%calib%ordered)

    allocate(spec%ordered(nConstraint),stat=err)

    if(err/=0)then
      message=trim(message)//'unable to allocate ordered parameter specifications'
      return
    endif


    do i=1,nConstraint

      allocate(spec%ordered(i)%param_index( &
                 size(config%calib%ordered(i)%parameters)),stat=err)

      if(err/=0)then
        message=trim(message)//'unable to allocate ordered parameter indices'
        return
      endif

      do j=1,size(config%calib%ordered(i)%parameters)

        ix = find_parameter(trim(config%calib%ordered(i)%parameters(j)), &
                            param_names(1:nParam))

        if(ix == 0)then
          message=trim(message)// &
                  'unable to resolve parameter in ordered constraint: '// &
                  trim(config%calib%ordered(i)%parameters(j))
          err=20; return
        endif

        spec%ordered(i)%param_index(j) = ix

      enddo

      spec%ordered(i)%gap_fraction = &
        config%calib%ordered(i)%gap_fraction

    enddo

  end subroutine get_summa_parameter_spec


  ! **************************************************************************************************
  ! Build the complete scalar SUMMA parameter override vector for one trial.
  !
  ! Sampled parameters receive their values from sampled_values. Non-sampled parameters represented
  ! in the parameter specification receive their scalar trial_value.
  !
  ! Because summa_paramSetup applies caller-supplied overrides after reading the trial parameter
  ! file, these overrides ensure that all parameters participating in calibration constraints are
  ! spatially uniform for the current implementation.
  !
  ! The returned param_names and param_values vectors should be passed directly to run_simulation()
  ! or evaluate_objective().
  ! **************************************************************************************************

  subroutine build_summa_parameter_overrides(spec, sampled_names, sampled_values, &
                                             param_names, param_values, err, message)

    implicit none

    type(parameter_spec), intent(in) :: spec

    character(*), intent(in) :: sampled_names(:)
    real(rkind),  intent(in) :: sampled_values(:)

    character(len=64), allocatable, intent(out) :: param_names(:)
    real(rkind),       allocatable, intent(out) :: param_values(:)

    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message

    integer(i4b) :: i
    integer(i4b) :: ixSample

    err = 0
    message = 'build_summa_parameter_overrides/'


    ! sampled name-value vectors must be consistent
    if(size(sampled_names) /= size(sampled_values))then
      message=trim(message)// &
              'sampled parameter name and value vectors have different sizes'
      err=20; return
    endif


    allocate(param_names(size(spec%params)), &
             param_values(size(spec%params)), &
             stat=err)

    if(err/=0)then
      message=trim(message)//'unable to allocate SUMMA parameter overrides'
      return
    endif


    do i=1,size(spec%params)

      param_names(i) = spec%params(i)%name

      if(spec%params(i)%sampled)then

        ixSample = find_parameter(trim(spec%params(i)%name),sampled_names)

        if(ixSample == 0)then
          message=trim(message)// &
                  'sampled value not provided for parameter: '// &
                  trim(spec%params(i)%name)
          err=20; return
        endif

        param_values(i) = sampled_values(ixSample)

      else

        ! Constraint-only parameters are reset to the scalar SUMMA default
        ! stored as trial_value in the model parameter specification.
        param_values(i) = spec%params(i)%trial_value

      endif

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

  subroutine get_summa_parameter_info(param_name,trial_value,lower,upper, &
                                      units,long_name,err,message)
  
    USE get_ixname_module, only: get_ixParam,get_ixBpar
  
    USE globalData, only: localParFallback,basinParFallback
    USE globalData, only: mpar_meta,bpar_meta
  
    implicit none
  
    character(*), intent(in)  :: param_name
    real(rkind),  intent(out) :: trial_value
    real(rkind),  intent(out) :: lower
    real(rkind),  intent(out) :: upper
    character(*), intent(out) :: units
    character(*), intent(out) :: long_name
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
  
    integer(i4b) :: ixParam,ixBasin
  
    err=0
    message='get_summa_parameter_info/'
  
    ! local parameter
    ixParam=get_ixParam(trim(param_name))
  
    if(ixParam > 0)then
  
      trial_value=localParFallback(ixParam)%default_val
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
    ixBasin=get_ixBpar(trim(param_name))
  
    if(ixBasin > 0)then
  
      trial_value=basinParFallback(ixBasin)%default_val
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

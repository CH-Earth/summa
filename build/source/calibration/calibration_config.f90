module calibration_config

  USE, INTRINSIC :: ieee_arithmetic, only: ieee_is_finite
  USE nr_type,           only: i4b,i8b,rkind,lgt
  USE parameter_search,  only: param_physical,param_spatial_multiplier,param_ratio
  USE calibration_types, only: calib_info,calib_param_info,target_info,nsga2_info
  
  USE tomlf_all, only: toml_table,toml_key,get_value

  implicit none
  private

  public :: parse_calibration_settings

  public :: parse_word_list
  public :: parse_nsga2_settings
  public :: parse_calibration_parameters
  public :: parse_parameter_dependencies
  public :: resolve_parameter_ratio_references

  public :: parse_target_spatial_unit
  public :: parse_calibration_targets

  public :: finalize_calibration_settings

contains

  ! --------------------------------------------------------------------------------------------------
  ! ---- PARSERS -------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------

  ! **************************************************************************************************
  ! Parse general calibration settings.
  !
  ! Reads scalar settings from [calibration]. Nested parameter specifications,
  ! optimization settings, and target definitions are parsed separately.
  ! **************************************************************************************************
  subroutine parse_calibration_settings(calib_table,calib,ierr,message)

    implicit none

    ! Dummy arguments
    type(toml_table), pointer, intent(in)    :: calib_table  ! Calibration TOML table
    type(calib_info),          intent(inout) :: calib        ! Calibration configuration
    integer(i4b),              intent(out)   :: ierr         ! Error code
    character(len=*),          intent(out)   :: message      ! Error message

    ! Local variables
    type(toml_key), allocatable      :: keys(:)  ! Calibration setting names
    character(len=:), allocatable   :: cvalue   ! String setting value
    integer(i4b)                    :: i        ! Setting index
    integer(i4b)                    :: istat    ! TOML read status

    ierr=0
    message='parse_calibration_settings/'

    call calib_table%get_keys(keys)
    if(.not.allocated(keys)) return

    do i=1,size(keys)

      istat=0

      select case(trim(keys(i)%key))

        ! Objective function settings
        case('metric')
          call get_value(calib_table,trim(keys(i)%key),calib%metric,stat=istat)

        case('obs_transform')
          call get_value(calib_table,trim(keys(i)%key),calib%obs_transform,stat=istat)

        case('boxcox_lambda')
          call get_value(calib_table,trim(keys(i)%key),calib%boxcox_lambda,stat=istat)
          if(istat==0)then
            if(.not.ieee_is_finite(calib%boxcox_lambda))then
              message=trim(message)//'calibration.boxcox_lambda must be finite'
              ierr=20; return
            endif
          endif

        ! Calibration period
        case('start_date')
          call get_value(calib_table,trim(keys(i)%key),calib%start_date,stat=istat)

        case('end_date')
          call get_value(calib_table,trim(keys(i)%key),calib%end_date,stat=istat)

        ! Optimization settings
        case('n_samples')
          call get_value(calib_table,trim(keys(i)%key),calib%n_samples,stat=istat)

        case('algorithm')
          call get_value(calib_table,trim(keys(i)%key),cvalue,stat=istat)
          if(istat/=0 .or. .not.allocated(cvalue))then
            message=trim(message)//'unable to read calibration.algorithm'
            ierr=20; return
          endif
          if(len_trim(cvalue)>len(calib%algorithm))then
            message=trim(message)//'calibration.algorithm exceeds maximum length'
            ierr=20; return
          endif
          calib%algorithm=trim(cvalue)

        ! Output settings
        case('write_aligned')
          call get_value(calib_table,trim(keys(i)%key),calib%write_aligned,stat=istat)

        ! Nested calibration settings are parsed separately.
        case('parameters','spatial_multipliers','parameter_ratios','nsga2','target')
          cycle

        case default
          message=trim(message)//'unexpected calibration setting: '//trim(keys(i)%key)
          ierr=20; return

      end select

      if(istat/=0)then
        message=trim(message)//'unable to read calibration.'//trim(keys(i)%key)
        ierr=20; return
      endif

    enddo

  end subroutine parse_calibration_settings

  ! **************************************************************************************************
  ! Parse the NSGA-II settings in [calibration.nsga2].  Keys left out keep their defaults.
  ! **************************************************************************************************
  subroutine parse_nsga2_settings(calib_table, calib, ierr, message)
    use tomlf_all, only: toml_table, toml_key, get_value
    implicit none

    ! Dummy arguments
    type(toml_table), pointer, intent(in)    :: calib_table  ! Calibration TOML table
    type(calib_info),          intent(inout) :: calib        ! Calibration configuration
    integer(i4b),              intent(out)   :: ierr         ! Error code
    character(len=*),          intent(out)   :: message      ! Error message

    ! Local variables
    type(toml_table), pointer    :: nsga2_table  ! NSGA-II settings table
    type(toml_key), allocatable  :: keys(:)      ! Keys in the NSGA-II settings table
    integer(i4b)                 :: i            ! Index of the current setting
    integer(i4b)                 :: istat        ! TOML read status

    ierr = 0
    message = 'parse_nsga2_settings/'

    call get_value(calib_table, 'nsga2', nsga2_table, stat=istat)
    if(istat/=0 .or. .not.associated(nsga2_table))then
      message=trim(message)//'unable to read the calibration.nsga2 table'
      ierr=20; return
    endif

    ! a multi-case run reuses one config structure, so start from the defaults
    calib%nsga2 = nsga2_info()

    call nsga2_table%get_keys(keys)
    if(.not.allocated(keys)) return
    do i=1,size(keys)
      associate(nsga2 => calib%nsga2)
      select case(trim(keys(i)%key))
        case('population_size');       call get_value(nsga2_table, trim(keys(i)%key), nsga2%population_size,       stat=istat)
        case('crossover_probability'); call get_value(nsga2_table, trim(keys(i)%key), nsga2%crossover_probability, stat=istat)
        case('crossover_eta');         call get_value(nsga2_table, trim(keys(i)%key), nsga2%crossover_eta,         stat=istat)
        case('mutation_probability');  call get_value(nsga2_table, trim(keys(i)%key), nsga2%mutation_probability,  stat=istat)
        case('mutation_eta');          call get_value(nsga2_table, trim(keys(i)%key), nsga2%mutation_eta,          stat=istat)
        case default
          message=trim(message)//'unexpected entry: calibration.nsga2.'//trim(keys(i)%key)
          ierr=20; return
      end select
      end associate
      if(istat/=0)then
        message=trim(message)//'unable to read calibration.nsga2.'//trim(keys(i)%key)
        ierr=20; return
      endif
    enddo

  end subroutine parse_nsga2_settings

  ! **************************************************************************************************
  ! Parse calibration parameter specifications.
  !
  ! Three parameterizations are supported:
  !
  !   [calibration.parameters]          Physical parameter values
  !   [calibration.spatial_multipliers] Spatial multipliers
  !   [calibration.parameter_ratios]    Parameter ratios
  !
  ! Each table contains parameter names and optional search bounds and
  ! transformations. All parameters are stored in a single calibration
  ! parameter list, with param_type identifying their parameterization.
  !
  ! Parameter-ratio references are resolved separately from the
  ! [parameter_dependencies] section after loading the complete TOML file.
  ! **************************************************************************************************
  
  subroutine parse_calibration_parameters(calib_table,calib,ierr,message)

    USE tomlf_all, only: toml_table,toml_key,get_value
    USE parameter_search, only: param_physical,param_spatial_multiplier,param_ratio
    implicit none

    type(toml_table), pointer, intent(in)    :: calib_table  ! Calibration TOML table
    type(calib_info),          intent(inout) :: calib        ! Calibration configuration
    integer(i4b),              intent(out)   :: ierr         ! Error code
    character(len=*),          intent(out)   :: message      ! Error message

    ! Supported calibration parameter tables.
    character(len=24), parameter :: table_names(3) = [character(len=24) :: &
      'parameters', 'spatial_multipliers', 'parameter_ratios']

    integer(i4b), parameter :: table_types(3) = [ &
      param_physical,param_spatial_multiplier,param_ratio]

    type(toml_table), pointer :: param_table             ! Current parameterization table
    type(toml_table), pointer :: parameter               ! Inline parameter specification
    type(toml_key), allocatable :: keys(:)               ! Parameter names in current table

    type(calib_param_info), allocatable :: entries(:)   ! Temporary parameter registry
    character(len=:), allocatable :: transformation     ! Parsed transformation
    character(len=:), allocatable :: reduction          ! Optional local ratio-reference reduction

    integer(i4b) :: iType                               ! Parameterization table index
    integer(i4b) :: i                                   ! Parameter index within table
    integer(i4b) :: j                                   ! Registry search index
    integer(i4b) :: nMax                                ! Maximum number of parameters
    integer(i4b) :: nParam                              ! Number of registered parameters
    integer(i4b) :: istat                               ! TOML read status

    ierr=0
    message='parse_calibration_parameters/'

    ! Multi-case runs reuse the configuration structure.
    if(allocated(calib%param_list)) deallocate(calib%param_list)

    ! -----------------------------------------------------------------------------------------------
    ! Reject obsolete calibration configuration keys.
    ! -----------------------------------------------------------------------------------------------
    if(calib_table%has_key('param_list') .or. &
       calib_table%has_key('parameter_transformations'))then
      message=trim(message)// &
        'param_list and parameter_transformations are no longer supported'
      ierr=20; return
    endif

    ! -----------------------------------------------------------------------------------------------
    ! Determine the maximum number of calibration parameters.
    !
    ! Names appearing in multiple tables are counted more than once here.
    ! Duplicate parameter names are rejected when constructing the registry.
    ! -----------------------------------------------------------------------------------------------
    nMax=0

    do iType=1,size(table_names)

      nullify(param_table)
      call get_value(calib_table,trim(table_names(iType)),param_table,requested=.false.,stat=istat)

      if(.not.associated(param_table))then
        if(calib_table%has_key(trim(table_names(iType))))then
          message=trim(message)//'invalid calibration parameter table: '//trim(table_names(iType))
          ierr=20; return
        endif
        cycle
      endif

      if(allocated(keys)) deallocate(keys)
      call param_table%get_keys(keys)

      if(allocated(keys)) nMax=nMax+size(keys)

    enddo

    if(nMax==0)then
      message=trim(message)//'no calibration parameters defined'
      ierr=20
      return
    endif

    allocate(entries(nMax),stat=ierr)
    if(ierr/=0)then
      message=trim(message)//'unable to allocate calibration parameter registry'
      return
    endif

    nParam=0

    ! -----------------------------------------------------------------------------------------------
    ! Read the three calibration parameter tables.
    !
    ! All three use the same bounds and transformation syntax. The table
    ! determines the parameterization assigned to each calibration variable.
    ! -----------------------------------------------------------------------------------------------
    do iType=1,size(table_names)

      nullify(param_table)
      call get_value(calib_table,trim(table_names(iType)),param_table,requested=.false.,stat=istat)
      if(.not.associated(param_table)) cycle

      if(allocated(keys)) deallocate(keys)
      call param_table%get_keys(keys)
      if(.not.allocated(keys)) cycle

      do i=1,size(keys)

        ! Parameters cannot appear in multiple calibration tables.
        do j=1,nParam
          if(trim(entries(j)%name)==trim(keys(i)%key))then
            message=trim(message)//'parameter appears in multiple calibration tables: '// &
                    trim(keys(i)%key)
            ierr=20; return
          endif
        enddo

        nParam=nParam+1

        if(len_trim(keys(i)%key)>len(entries(nParam)%name))then
          message=trim(message)//'calibration parameter name is too long: '//trim(keys(i)%key)
          ierr=20; return
        endif

        entries(nParam)%name=trim(keys(i)%key)
        entries(nParam)%param_type=table_types(iType)

        ! Read the inline parameter specification.
        nullify(parameter)
        call get_value(param_table,trim(keys(i)%key),parameter,stat=istat)

        if(istat/=0 .or. .not.associated(parameter))then
          message=trim(message)//'unable to read specification for parameter: '// &
                  trim(keys(i)%key)
          ierr=20; return
        endif

        ! Optional lower bound.
        if(parameter%has_key('lower'))then

          call get_value(parameter,'lower',entries(nParam)%lower_bound,stat=istat)

          if(istat/=0)then
            message=trim(message)//'unable to read lower bound for parameter: '//trim(keys(i)%key)
            ierr=20; return
          endif

          entries(nParam)%has_lower=.true.

        endif

        ! Optional upper bound.
        if(parameter%has_key('upper'))then

          call get_value(parameter,'upper',entries(nParam)%upper_bound,stat=istat)

          if(istat/=0)then
            message=trim(message)//'unable to read upper bound for parameter: '//trim(keys(i)%key)
            ierr=20; return
          endif

          entries(nParam)%has_upper=.true.

        endif

        ! Validate bounds when both are explicitly supplied.
        if(entries(nParam)%has_lower .and. entries(nParam)%has_upper)then

          if(entries(nParam)%lower_bound>=entries(nParam)%upper_bound)then
            message=trim(message)//'lower bound must be less than upper bound for parameter: '// &
                    trim(keys(i)%key)
            ierr=20; return
          endif

        endif

        ! Optional search transformation.
        if(parameter%has_key('transformation'))then

          if(allocated(transformation)) deallocate(transformation)

          call get_value(parameter,'transformation',transformation,stat=istat)

          if(istat/=0 .or. .not.allocated(transformation))then
            message=trim(message)//'unable to read transformation for parameter: '//trim(keys(i)%key)
            ierr=20; return
          endif

          if(len_trim(transformation)>len(entries(nParam)%transformation))then
            message=trim(message)//'transformation name too long for parameter: '//trim(keys(i)%key)
            ierr=20; return
          endif

          entries(nParam)%transformation=trim(transformation)

        endif

        ! Ratios remain elementwise unless a local reduction is explicitly requested.
        if(parameter%has_key('reference_reduction'))then
          if(entries(nParam)%param_type/=param_ratio)then
            message=trim(message)//'reference_reduction is only valid for parameter ratios: '//trim(keys(i)%key)
            ierr=20; return
          endif
          if(allocated(reduction)) deallocate(reduction)
          call get_value(parameter,'reference_reduction',reduction,stat=istat)
          if(istat/=0 .or. .not.allocated(reduction))then
            message=trim(message)//'unable to read reference_reduction for parameter: '//trim(keys(i)%key)
            ierr=20; return
          endif
          select case(trim(reduction))
          case('none','min')
            entries(nParam)%reference_reduction=trim(reduction)
          case default
            message=trim(message)//'reference_reduction must be none or min: '//trim(keys(i)%key)
            ierr=20; return
          end select
        endif

      enddo
    enddo

    ! -----------------------------------------------------------------------------------------------
    ! Store the calibration parameter registry.
    ! -----------------------------------------------------------------------------------------------
    allocate(calib%param_list(nParam),stat=ierr)

    if(ierr/=0)then
      message=trim(message)//'unable to allocate calibration parameter list'
      return
    endif

    calib%param_list=entries(1:nParam)

  end subroutine parse_calibration_parameters
  
  ! **************************************************************************************************
  ! Parse calibration parameter dependencies from TOML.
  !
  ! Supported dependencies:
  !
  !   [[parameter_dependencies.ordered]]
  !   parameters = ["parameter1", "parameter2", ...]
  !   gap_fraction = 0.025
  !
  !   [[parameter_dependencies.parameter_ratio]]
  !   reference = "k_soil"
  !   dependent = "k_macropore"
  !
  ! For parameter ratios:
  !   physical(dependent) = sampled(dependent) * physical(reference)
  ! **************************************************************************************************
  subroutine parse_parameter_dependencies(subtable,calib,ierr,message)

    USE tomlf_all, only: toml_table,toml_array,get_value,len

    implicit none

    ! ---------------------------------------------------------------------------------------
    ! Dummy arguments
    ! ---------------------------------------------------------------------------------------
    type(toml_table), pointer, intent(in)    :: subtable  ! Parameter dependencies TOML table
    type(calib_info),          intent(inout) :: calib     ! Calibration configuration
    integer(i4b),              intent(out)   :: ierr      ! Error code
    character(len=*),          intent(out)   :: message   ! Error message

    ! ---------------------------------------------------------------------------------------
    ! Local variables
    ! ---------------------------------------------------------------------------------------
    type(toml_array), pointer :: ordered                 ! Ordered parameter constraints
    type(toml_array), pointer :: ratios                  ! Parameter-ratio dependencies
    type(toml_array), pointer :: param_list              ! Parameter names in an ordered constraint
    type(toml_table), pointer :: constraint              ! Individual dependency specification

    character(len=:), allocatable :: reference           ! Reference parameter name
    character(len=:), allocatable :: dependent           ! Dependent parameter name
    character(len=1024)           :: cmessage            ! Error message from helper routines

    integer(i4b) :: i                                    ! Constraint index
    integer(i4b) :: j                                    ! Index for duplicate checks
    integer(i4b) :: istat                                ! TOML read status
    integer(i4b) :: nconstraints                         ! Number of constraints

    ierr=0
    message='parse_parameter_dependencies/'

    nullify(ordered,ratios,constraint,param_list)

    ! Reset parameter dependencies for each calibration case.
    if(allocated(calib%ordered)) deallocate(calib%ordered)
    if(allocated(calib%multiplier)) deallocate(calib%multiplier)

    ! ---------------------------------------------------------------------------------------
    ! Ordered parameter dependencies
    ! ---------------------------------------------------------------------------------------
    call get_value(subtable,'ordered',ordered,requested=.false.,stat=istat)

    if(associated(ordered))then

      nconstraints=len(ordered)
      allocate(calib%ordered(nconstraints),stat=ierr)
      if(ierr/=0)then
        message=trim(message)//'unable to allocate ordered parameter constraints'
        return
      endif

      do i=1,nconstraints

        nullify(constraint,param_list)

        call get_value(ordered,i,constraint,stat=istat)
        if(istat/=0 .or. .not.associated(constraint))then
          write(cmessage,'(A,I0)') 'unable to read ordered constraint, i = ',i
          message=trim(message)//trim(cmessage)
          ierr=20; return
        endif

        ! Reject the obsolete gap-sampling option.
        if(constraint%has_key('sample_gaps'))then
          message=trim(message)//'sample_gaps is no longer supported'
          ierr=20; return
        endif

        ! Parameter list
        call get_value(constraint,'parameters',param_list,stat=istat)
        if(istat/=0 .or. .not.associated(param_list))then
          write(cmessage,'(A,I0)') 'parameter list not defined for ordered constraint, i = ',i
          message=trim(message)//trim(cmessage)
          ierr=20; return
        endif

        call parse_word_list(param_list,calib%ordered(i)%parameters,ierr,cmessage)
        if(ierr/=0)then
          message=trim(message)//trim(cmessage)
          return
        endif

        ! Minimum gap between adjacent parameters
        call get_value(constraint,'gap_fraction',calib%ordered(i)%gap_fraction,stat=istat)
        if(istat/=0)then
          write(cmessage,'(A,I0)') 'gap_fraction not defined for ordered constraint, i = ',i
          message=trim(message)//trim(cmessage)
          ierr=20; return
        endif
    
      enddo

    endif

    ! ---------------------------------------------------------------------------------------
    ! Parameter-ratio dependencies
    ! ---------------------------------------------------------------------------------------

    ! Reject the obsolete multiplier syntax.
    if(subtable%has_key('multiplier'))then
      message=trim(message)//'multiplier dependencies are obsolete; use parameter_ratio'
      ierr=20; return
    endif

    call get_value(subtable,'parameter_ratio',ratios,requested=.false.,stat=istat)

    if(subtable%has_key('parameter_ratio') .and. .not.associated(ratios))then
      message=trim(message)//'unable to read parameter_ratio dependencies'
      ierr=20; return
    endif

    if(associated(ratios))then

      nconstraints=len(ratios)
      allocate(calib%multiplier(nconstraints),stat=ierr)
      if(ierr/=0)then
        message=trim(message)//'unable to allocate parameter-ratio dependencies'
        return
      endif

      do i=1,nconstraints

        nullify(constraint)

        call get_value(ratios,i,constraint,stat=istat)
        if(istat/=0 .or. .not.associated(constraint))then
          write(cmessage,'(A,I0)') 'unable to read parameter-ratio dependency, i = ',i
          message=trim(message)//trim(cmessage)
          ierr=20; return
        endif

        ! Read the reference parameter.
        if(allocated(reference)) deallocate(reference)
        call get_value(constraint,'reference',reference,stat=istat)

        if(istat/=0 .or. .not.allocated(reference))then
          write(cmessage,'(A,I0)') 'reference not defined for parameter ratio, i = ',i
          message=trim(message)//trim(cmessage)
          ierr=20; return
        endif

        ! Read the dependent parameter.
        if(allocated(dependent)) deallocate(dependent)
        call get_value(constraint,'dependent',dependent,stat=istat)

        if(istat/=0 .or. .not.allocated(dependent))then
          write(cmessage,'(A,I0)') 'dependent not defined for parameter ratio, i = ',i
          message=trim(message)//trim(cmessage)
          ierr=20; return
        endif

        ! Validate parameter names.
        if(len_trim(reference)==0 .or. len_trim(dependent)==0)then
          write(cmessage,'(A,I0)') 'empty parameter name in parameter ratio, i = ',i
          message=trim(message)//trim(cmessage)
          ierr=20; return
        endif

        if(trim(reference)==trim(dependent))then
          write(cmessage,'(A,I0)') 'reference and dependent cannot be identical, i = ',i
          message=trim(message)//trim(cmessage)
          ierr=20; return
        endif

        ! Store the dependency using the existing configuration structure.
        allocate(calib%multiplier(i)%parameters(2),stat=ierr)
        if(ierr/=0)then
          message=trim(message)//'unable to allocate parameter-ratio names'
          return
        endif

        if(len_trim(reference)>len(calib%multiplier(i)%parameters(1)) .or. &
           len_trim(dependent)>len(calib%multiplier(i)%parameters(2)))then
          message=trim(message)//'parameter name too long in parameter ratio'
          ierr=20; return
        endif

        ! A dependent can have only one ratio reference.
        do j=1,i-1
          if(trim(calib%multiplier(j)%parameters(2))==trim(dependent))then
            message=trim(message)//'duplicate parameter-ratio dependent: '//trim(dependent)
            ierr=20; return
          endif
        enddo

        calib%multiplier(i)%parameters(1)=trim(reference)
        calib%multiplier(i)%parameters(2)=trim(dependent)

      enddo

    endif

  end subroutine parse_parameter_dependencies

  ! **************************************************************************************************
  ! Parse the spatial selection for a calibration target.
  !
  ! A target may identify a spatial element using a model-defined spatial unit
  ! and identifier. If no spatial unit is specified, the target applies to the
  ! entire model domain. Interpretation of spatial units is model-specific.
  ! **************************************************************************************************
  subroutine parse_target_spatial_unit(target_table,calTarget,ierr,message)

    implicit none

    ! Dummy arguments
    type(toml_table), pointer, intent(in)    :: target_table  ! Calibration target TOML table
    type(target_info),         intent(inout) :: calTarget     ! Calibration target
    integer(i4b),              intent(out)   :: ierr          ! Error code
    character(len=*),          intent(inout) :: message       ! Error message

    ! Local variables
    character(len=:), allocatable :: unit_name  ! Model-defined spatial unit
    logical                       :: has_unit   ! Whether spatial_unit is specified
    logical                       :: has_id     ! Whether spatial_id is specified
    integer(i4b)                  :: istat      ! TOML read status

    ierr=0

    has_unit=target_table%has_key('spatial_unit')
    has_id=target_table%has_key('spatial_id')

    ! No spatial selection means the entire model domain.
    if(.not.has_unit)then

      if(has_id)then
        message=trim(message)//'spatial_id requires spatial_unit for target: '// &
                trim(calTarget%name)
        ierr=20; return
      endif

      return

    endif

    ! Read the model-defined spatial unit.
    call get_value(target_table,'spatial_unit',unit_name,stat=istat)

    if(istat/=0 .or. .not.allocated(unit_name))then
      message=trim(message)//'unable to read spatial_unit for target: '//trim(calTarget%name)
      ierr=20; return
    endif

    if(len_trim(unit_name)==0 .or. len_trim(unit_name)>len(calTarget%spatial_unit))then
      message=trim(message)//'invalid spatial_unit for target: '//trim(calTarget%name)
      ierr=20; return
    endif

    calTarget%spatial_unit=trim(unit_name)

    ! The entire domain has no individual spatial identifier.
    if(trim(unit_name)=='domain')then

      if(has_id)then
        message=trim(message)//'domain target must not specify spatial_id: '// &
                trim(calTarget%name)
        ierr=20; return
      endif

      return

    endif

    ! All other spatial units require an identifier.
    if(.not.has_id)then
      message=trim(message)//'spatial_id not defined for target: '//trim(calTarget%name)
      ierr=20; return
    endif

    call get_value(target_table,'spatial_id',calTarget%spatial_id,stat=istat)

    if(istat/=0)then
      message=trim(message)//'unable to read spatial_id for target: '//trim(calTarget%name)
      ierr=20; return
    endif

  end subroutine parse_target_spatial_unit


  ! **************************************************************************************************
  ! Parse calibration targets.
  !
  ! Reads the [[calibration.target]] array of tables. Each target defines an
  ! observed series, the simulated variable it is compared against, and the
  ! metric used to evaluate their agreement.
  !
  ! Observation paths, variable names, and objective settings may inherit
  ! defaults from the enclosing configuration.
  !
  ! Spatial selections are specified using a model-defined spatial unit
  ! and identifier. Interpretation of these values is model-specific.
  ! **************************************************************************************************
  subroutine parse_calibration_targets(subtable,calib,default_obs_path,default_vname_obs,ierr,message)

    USE tomlf_all, only: toml_table,toml_array,toml_key,get_value,len
    USE, INTRINSIC :: iso_fortran_env, only: error_unit

    implicit none

    ! Dummy arguments
    type(toml_table), pointer, intent(in)      :: subtable          ! Calibration TOML table
    type(calib_info),          intent(inout)   :: calib             ! Calibration configuration
    character(len=:), allocatable, intent(in) :: default_obs_path  ! Default observation path
    character(len=:), allocatable, intent(in) :: default_vname_obs ! Default observation variable
    integer(i4b),              intent(out)     :: ierr              ! Error code
    character(len=*),          intent(out)     :: message           ! Error message

    ! Local variables
    type(toml_array), pointer       :: targets         ! Array of calibration targets
    type(toml_table), pointer       :: target_table    ! Current target table
    type(toml_key), allocatable     :: target_keys(:)  ! Keys in current target
    character(len=:), allocatable  :: cvalue          ! String setting value
    integer(i4b)                   :: i               ! Target index
    integer(i4b)                   :: j               ! Target key index
    integer(i4b)                   :: istat           ! TOML read status
    integer(i4b)                   :: ntargets        ! Number of calibration targets

    ierr=0
    message='parse_calibration_targets/'

    ! -----------------------------------------------------------------------------------------------
    ! Get the array of calibration targets.
    ! -----------------------------------------------------------------------------------------------
    nullify(targets)
    call get_value(subtable,'target',targets,requested=.false.,stat=istat)

    if(.not.associated(targets)) return

    ntargets=len(targets)

    if(ntargets<1)then
      message=trim(message)//'[[calibration.target]] is present but names no targets'
      ierr=20; return
    endif

    ! A multi-case run may reuse the same calibration configuration.
    if(allocated(calib%targets)) deallocate(calib%targets)

    allocate(calib%targets(ntargets),stat=ierr)
    if(ierr/=0)then
      message=trim(message)//'unable to allocate calibration targets'
      return
    endif

    ! -----------------------------------------------------------------------------------------------
    ! Parse each calibration target.
    ! -----------------------------------------------------------------------------------------------
    do i=1,ntargets

      nullify(target_table)
      call get_value(targets,i,target_table,stat=istat)

      if(istat/=0 .or. .not.associated(target_table))then
        write(message,'(A,I0)') trim(message)//'unable to read calibration target, i = ',i
        ierr=20; return
      endif

      ! ---------------------------------------------------------------------------------------------
      ! Check for unknown or obsolete target settings.
      ! ---------------------------------------------------------------------------------------------
      if(allocated(target_keys)) deallocate(target_keys)
      call target_table%get_keys(target_keys)

      if(allocated(target_keys))then

        do j=1,size(target_keys)

          select case(trim(target_keys(j)%key))

            case('name','variable','obs_path','obs_file','vname_obs', &
                 'metric','obs_transform','boxcox_lambda','weight','obs_units', &
                 'spatial_unit','spatial_id','accumulate',           &
                 'baseline_start','baseline_end','balance',         &
                 'balance_window','low_flow','flow_duration')
              cycle

            case('gru','hru','reach')
              message=trim(message)//'obsolete spatial selection key: '// &
                      trim(target_keys(j)%key)//'; use spatial_unit and spatial_id'
              ierr=20; return

            case default
              write(error_unit,'(3A)') &
                'WARNING: ignoring unknown calibration target setting "', &
                trim(target_keys(j)%key), '"'

          end select

        enddo

      endif

      ! ---------------------------------------------------------------------------------------------
      ! Observation file.
      ! ---------------------------------------------------------------------------------------------
      if(allocated(cvalue)) deallocate(cvalue)
      call get_value(target_table,'obs_file',cvalue,stat=istat)

      if(istat/=0 .or. .not.allocated(cvalue))then
        write(message,'(A,I0)') trim(message)//'obs_file not defined for calibration target, i = ',i
        ierr=20; return
      endif

      calib%targets(i)%obs_file=trim(cvalue)

      ! ---------------------------------------------------------------------------------------------
      ! Target name.
      ! ---------------------------------------------------------------------------------------------
      if(allocated(cvalue)) deallocate(cvalue)
      call get_value(target_table,'name',cvalue,stat=istat)

      if(istat==0 .and. allocated(cvalue))then
        calib%targets(i)%name=trim(cvalue)
      else
        write(calib%targets(i)%name,'(A,I0)') 'target',i
      endif

      ! ---------------------------------------------------------------------------------------------
      ! Simulated variable.
      ! ---------------------------------------------------------------------------------------------
      if(allocated(cvalue)) deallocate(cvalue)
      call get_value(target_table,'variable',cvalue,stat=istat)

      if(istat==0 .and. allocated(cvalue))then
        calib%targets(i)%variable=trim(cvalue)
      endif

      ! ---------------------------------------------------------------------------------------------
      ! Observation path.
      ! ---------------------------------------------------------------------------------------------
      if(allocated(cvalue)) deallocate(cvalue)
      call get_value(target_table,'obs_path',cvalue,stat=istat)

      if(istat==0 .and. allocated(cvalue))then
        calib%targets(i)%obs_path=trim(cvalue)
      else if(allocated(default_obs_path))then
        calib%targets(i)%obs_path=trim(default_obs_path)
      else
        calib%targets(i)%obs_path=''
      endif

      ! ---------------------------------------------------------------------------------------------
      ! Observed variable name.
      ! ---------------------------------------------------------------------------------------------
      if(allocated(cvalue)) deallocate(cvalue)
      call get_value(target_table,'vname_obs',cvalue,stat=istat)

      if(istat==0 .and. allocated(cvalue))then
        calib%targets(i)%vname_obs=trim(cvalue)
      else if(allocated(default_vname_obs))then
        calib%targets(i)%vname_obs=trim(default_vname_obs)
      else
        calib%targets(i)%vname_obs='q_obs'
      endif

      ! ---------------------------------------------------------------------------------------------
      ! Objective metric.
      ! ---------------------------------------------------------------------------------------------
      if(allocated(cvalue)) deallocate(cvalue)
      call get_value(target_table,'metric',cvalue,stat=istat)

      if(istat==0 .and. allocated(cvalue))then
        calib%targets(i)%metric=trim(cvalue)
      else if(allocated(calib%metric))then
        calib%targets(i)%metric=trim(calib%metric)
      endif

      ! ---------------------------------------------------------------------------------------------
      ! Observation transformation.
      ! ---------------------------------------------------------------------------------------------
      if(allocated(cvalue)) deallocate(cvalue)
      call get_value(target_table,'obs_transform',cvalue,stat=istat)

      if(istat==0 .and. allocated(cvalue))then
        calib%targets(i)%obs_transform=trim(cvalue)
      else if(allocated(calib%obs_transform))then
        calib%targets(i)%obs_transform=trim(calib%obs_transform)
      endif

      ! ---------------------------------------------------------------------------------------------
      ! Box-Cox exponent: inherit the calibration default, then override locally.
      ! It is not a sampled model parameter, and is ignored by other transforms.
      ! ---------------------------------------------------------------------------------------------
      calib%targets(i)%boxcox_lambda=calib%boxcox_lambda
      if(target_table%has_key('boxcox_lambda'))then
        call get_value(target_table,'boxcox_lambda',calib%targets(i)%boxcox_lambda,stat=istat)
        if(istat/=0)then
          message=trim(message)//'invalid boxcox_lambda for calibration target "'// &
                  trim(calib%targets(i)%name)//'"; expected a real number'
          ierr=20; return
        endif
      endif
      if(.not.ieee_is_finite(calib%targets(i)%boxcox_lambda))then
        message=trim(message)//'boxcox_lambda must be finite for calibration target "'// &
                trim(calib%targets(i)%name)//'"'
        ierr=20; return
      endif

      ! ---------------------------------------------------------------------------------------------
      ! Objective weight.
      ! ---------------------------------------------------------------------------------------------
      if(target_table%has_key('weight'))then
        call get_value(target_table,'weight',calib%targets(i)%weight,stat=istat)
        if(istat/=0)then
          message=trim(message)//'invalid weight for calibration target "'// &
                  trim(calib%targets(i)%name)//'"'
          ierr=20; return
        endif
      endif

      ! ---------------------------------------------------------------------------------------------
      ! Observation units.
      ! ---------------------------------------------------------------------------------------------
      if(target_table%has_key('obs_units'))then
        if(allocated(cvalue)) deallocate(cvalue)
        call get_value(target_table,'obs_units',cvalue,stat=istat)
        if(istat/=0 .or. .not.allocated(cvalue))then
          message=trim(message)//'invalid obs_units for calibration target "'// &
                  trim(calib%targets(i)%name)//'"'
          ierr=20; return
        endif
        calib%targets(i)%obs_units=trim(cvalue)
      endif

      ! ---------------------------------------------------------------------------------------------
      ! Spatial selection.
      ! ---------------------------------------------------------------------------------------------
      call parse_target_spatial_unit(target_table,calib%targets(i),ierr,message)
      if(ierr/=0) return

      ! ---------------------------------------------------------------------------------------------
      ! Accumulation of simulated values.
      ! ---------------------------------------------------------------------------------------------
      if(target_table%has_key('accumulate'))then
        call get_value(target_table,'accumulate',calib%targets(i)%accumulate,stat=istat)
        if(istat/=0)then
          message=trim(message)//'invalid accumulate for calibration target "'// &
                  trim(calib%targets(i)%name)//'"; expected a logical value'
          ierr=20; return
        endif
      endif

      ! ---------------------------------------------------------------------------------------------
      ! Anomaly baseline.
      ! ---------------------------------------------------------------------------------------------
      if(target_table%has_key('baseline_start'))then
        if(allocated(cvalue)) deallocate(cvalue)
        call get_value(target_table,'baseline_start',cvalue,stat=istat)
        if(istat/=0 .or. .not.allocated(cvalue))then
          message=trim(message)//'invalid baseline_start for calibration target "'// &
                  trim(calib%targets(i)%name)//'"'
          ierr=20; return
        endif
        calib%targets(i)%baseline_start=trim(cvalue)
      endif

      if(target_table%has_key('baseline_end'))then
        if(allocated(cvalue)) deallocate(cvalue)
        call get_value(target_table,'baseline_end',cvalue,stat=istat)
        if(istat/=0 .or. .not.allocated(cvalue))then
          message=trim(message)//'invalid baseline_end for calibration target "'// &
                  trim(calib%targets(i)%name)//'"'
          ierr=20; return
        endif
        calib%targets(i)%baseline_end=trim(cvalue)
      endif

      ! ---------------------------------------------------------------------------------------------
      ! Seasonal or annual balance.
      ! ---------------------------------------------------------------------------------------------
      if(target_table%has_key('balance'))then
        if(allocated(cvalue)) deallocate(cvalue)
        call get_value(target_table,'balance',cvalue,stat=istat)
        if(istat/=0 .or. .not.allocated(cvalue))then
          message=trim(message)//'invalid balance for calibration target "'// &
                  trim(calib%targets(i)%name)//'"'
          ierr=20; return
        endif

        select case(trim(cvalue))
          case('seasonal','annual')
            calib%targets(i)%balance=trim(cvalue)
          case default
            message=trim(message)//'calibration target "'//trim(calib%targets(i)%name)// &
                    '" asks for balance "'//trim(cvalue)//'"; use seasonal or annual'
            ierr=20; return
        end select
      endif

      if(target_table%has_key('balance_window'))then
        call get_value(target_table,'balance_window',calib%targets(i)%balance_window,stat=istat)
        if(istat/=0)then
          message=trim(message)//'invalid balance_window for calibration target "'// &
                  trim(calib%targets(i)%name)//'"'
          ierr=20; return
        endif
      endif

      ! a low-flow target: the lowest fraction of the observations, paired in time or along the flow-duration curves
      if(target_table%has_key('low_flow'))then
        call get_value(target_table,'low_flow',calib%targets(i)%low_flow,stat=istat)
        if(istat/=0 .or. calib%targets(i)%low_flow <= 0._rkind .or. calib%targets(i)%low_flow > 1._rkind)then
          message=trim(message)//'calibration target "'//trim(calib%targets(i)%name)// &
                  '" needs a low_flow fraction in (0,1]'
          ierr=20; return
        endif
      endif
      if(target_table%has_key('flow_duration'))then
        call get_value(target_table,'flow_duration',calib%targets(i)%flow_duration,stat=istat)
        if(istat/=0)then
          message=trim(message)//'invalid flow_duration for calibration target "'// &
                  trim(calib%targets(i)%name)//'"'
          ierr=20; return
        endif
      endif

    enddo

  end subroutine parse_calibration_targets

  ! --------------------------------------------------------------------------------------------------
  ! ---- HELPERS -------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------

  ! **************************************************************************************************
  ! Resolve parameter-ratio references after parsing the TOML configuration.
  !
  ! Calibration parameters are defined in [calibration.parameter_ratios], while their
  ! dependencies are defined separately in [[parameter_dependencies.parameter_ratio]].
  ! This routine connects these definitions by storing the reference parameter name
  ! in the corresponding calibration parameter entry.
  !
  ! Reference parameters do not need to be sampled; they may use initialized
  ! parameter fields. Their existence and spatial compatibility are checked later
  ! when constructing the model parameter specification and fields.
  ! **************************************************************************************************
  subroutine resolve_parameter_ratio_references(calib,err,message)

    implicit none

    ! ---------------------------------------------------------------------------------------
    ! Dummy arguments
    ! ---------------------------------------------------------------------------------------
    type(calib_info),  intent(inout) :: calib    ! Calibration configuration
    integer(i4b),     intent(out)   :: err      ! Error code
    character(len=*), intent(out)   :: message  ! Error message

    ! ---------------------------------------------------------------------------------------
    ! Local variables
    ! ---------------------------------------------------------------------------------------
    integer(i4b)      :: i             ! Dependency or parameter index
    integer(i4b)      :: j             ! Parameter registry search index
    integer(i4b)      :: ixDependent   ! Index of the dependent parameter
    character(len=64) :: reference     ! Reference parameter name
    character(len=64) :: dependent     ! Dependent parameter name

    err=0
    message='resolve_parameter_ratio_references/'

    ! Calibration parameters must be defined before resolving dependencies.
    if(.not.allocated(calib%param_list))then
      if(allocated(calib%multiplier))then
        message=trim(message)//'calibration parameters are not defined'
        err=20; return
      endif
      return
    endif

    ! ---------------------------------------------------------------------------------------
    ! Resolve parameter-ratio dependencies.
    ! ---------------------------------------------------------------------------------------
    if(allocated(calib%multiplier))then

      do i=1,size(calib%multiplier)

        if(.not.allocated(calib%multiplier(i)%parameters))then
          message=trim(message)//'parameter-ratio dependency is not defined'
          err=20; return
        endif

        if(size(calib%multiplier(i)%parameters)/=2)then
          message=trim(message)//'parameter-ratio dependency must contain reference and dependent'
          err=20; return
        endif

        reference=trim(calib%multiplier(i)%parameters(1))
        dependent=trim(calib%multiplier(i)%parameters(2))

        ! Find the dependent parameter in the calibration registry.
        ixDependent=0

        do j=1,size(calib%param_list)
          if(trim(calib%param_list(j)%name)==trim(dependent))then
            ixDependent=j
            exit
          endif
        enddo

        if(ixDependent==0)then
          message=trim(message)//'parameter-ratio dependent is not calibrated: '//trim(dependent)
          err=20; return
        endif

        ! The dependent must be configured as a parameter ratio.
        if(calib%param_list(ixDependent)%param_type/=param_ratio)then
          message=trim(message)//'dependent is not configured as a parameter ratio: '//trim(dependent)
          err=20; return
        endif

        ! Each dependent parameter can have only one reference.
        if(len_trim(calib%param_list(ixDependent)%reference_name)>0)then
          message=trim(message)//'duplicate parameter-ratio dependency: '//trim(dependent)
          err=20; return
        endif

        ! Store the reference parameter name.
        calib%param_list(ixDependent)%reference_name=reference

      enddo

    endif

    ! ---------------------------------------------------------------------------------------
    ! Verify that every calibrated parameter ratio has a reference.
    ! ---------------------------------------------------------------------------------------
    do i=1,size(calib%param_list)

      if(calib%param_list(i)%param_type/=param_ratio) cycle

      if(len_trim(calib%param_list(i)%reference_name)==0)then
        message=trim(message)//'missing reference for parameter ratio: '// &
                trim(calib%param_list(i)%name)
        err=20; return
      endif

    enddo

  end subroutine resolve_parameter_ratio_references


  ! **************************************************************************************************
  ! Finalize calibration settings.
  !
  ! Apply default settings and validate the objective-function metric.
  ! Observation transformations are initialized only when observations are
  ! available, since a simulation may be run without calibration.
  ! **************************************************************************************************
  subroutine finalize_calibration_settings(calib,has_observations,err,message)

    implicit none

    ! Dummy arguments
    type(calib_info), intent(inout) :: calib             ! Calibration configuration
    logical(lgt),     intent(in)    :: has_observations  ! Whether observations are available
    integer(i4b),     intent(out)   :: err               ! Error code
    character(len=*), intent(out)   :: message           ! Error message

    err=0
    message='finalize_calibration_settings/'

    ! Default objective-function metric.
    if(.not.allocated(calib%metric)) calib%metric='kge'

    ! Default observation transformation.
    if(has_observations)then
      if(.not.allocated(calib%obs_transform)) calib%obs_transform='none'
    endif

    ! Validate metric and determine the default objective-function value.
    select case(trim(calib%metric))

      case('kge','kgep','nse')
        calib%default_metric=-9999._rkind

      case('mae','rmse')
        calib%default_metric=9999._rkind

      case default
        message=trim(message)//'unsupported calibration metric "'//trim(calib%metric)// &
                '"; expected one of: kge, kgep, nse, mae, rmse (case sensitive)'
        err=20
        return

    end select

  end subroutine finalize_calibration_settings

  ! **************************************************************************************************
  ! Parse a TOML array containing a list of words.
  ! **************************************************************************************************
  subroutine parse_word_list(word_list, words, ierr, message)
    use tomlf_all, only: toml_array, get_value, len
    implicit none

    type(toml_array), pointer, intent(in)          :: word_list
    character(len=*), allocatable, intent(out)     :: words(:)
    integer(i4b), intent(out)                      :: ierr
    character(*), intent(out)                      :: message
    integer(i4b)                  :: i
    integer(i4b)                  :: nwords
    character(len=:), allocatable :: word

    ierr = 0
    message = 'parse_word_list/'

    ! check that the TOML array exists
    if(.not.associated(word_list))then
      if(allocated(words)) deallocate(words)
      allocate(words(0))
      return
    endif
    nwords = len(word_list)

    ! NOTE: as elsewhere, a multi-case run reuses one config, so the caller's array may
    !       still hold the previous case's words
    if(allocated(words)) deallocate(words)

    ! allocate output array
    allocate(words(nwords), stat=ierr)
    if(ierr/=0)then
      message=trim(message)//'unable to allocate word list'
      return
    endif

    ! populate output array
    do i=1,nwords
      call get_value(word_list, i, word, stat=ierr)
      if(ierr/=0)then
        write(message,'(A,I0)') trim(message)//'unable to read word, i = ',i
        return
      endif

      ! prevent silent truncation
      if(len_trim(word) > len(words))then
        write(message,'(A,I0,A,I0)') trim(message)// &
          'word exceeds maximum character length, i = ', i, ', maximum length = ', len(words)
        ierr=20; return
      endif
      words(i) = trim(word)

    enddo

  end subroutine parse_word_list

end module calibration_config


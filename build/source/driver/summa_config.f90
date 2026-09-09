module summa_config

use build_options, only: mizuroute_active

USE nr_type
USE summa_type, only: config_info       ! summa configuation info

USE globalData, only: iulog             ! I/O unit for logging messages

USE globalData, only: iRunMode, iRunModeFull

#ifdef MIZUROUTE_ACTIVE
USE mizuroute_config, ONLY: parse_mizuroute_config
#endif

implicit none
private

public :: read_summa_config 

contains

  ! **************************************************************************************************
  ! Read SUMMA TOML configuration.
  !
  ! The TOML configuration may contain settings that overlap with values previously read
  ! from the legacy SUMMA file manager. When such values are provided in the TOML file,
  ! they overwrite the corresponding legacy file-manager values.
  ! **************************************************************************************************
  
  subroutine read_summa_config(config_file, config, err, message)
  
    implicit none
  
    character(*),      intent(in)    :: config_file
    type(config_info), intent(inout) :: config
    integer(i4b),      intent(out)   :: err
    character(*),      intent(out)   :: message
  
    character(len=256) :: cmessage
  
    err = 0
    message = 'read_summa_config/'
  
    ! load configuration values from the TOML file
    call load_summa_config(config_file, config, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  
    ! apply TOML values that supersede legacy file-manager settings
    call apply_summa_config(config, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  
  end subroutine read_summa_config

  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! ---- PRIVATE SUBROUTINES -------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  !
  ! **************************************************************************************************
  ! Load SUMMA configuration from a TOML file.
  !
  ! Loads the TOML configuration file and delegates parsing to component-specific readers.
  ! TOML configuration may include general SUMMA settings as well as configuration for optional
  ! components such as mizuRoute. Settings provided in the TOML file take precedence over
  ! corresponding values previously read from the legacy SUMMA file manager.
  ! **************************************************************************************************

  subroutine load_summa_config(config_file, config, err, message)


  use tomlf_all, only: toml_table, toml_array, toml_error, toml_key, toml_value ! data types
  use tomlf_all, only: toml_load, get_value, len                                ! procedures

  implicit none

  character(*),            intent(in)    :: config_file
  type(config_info),       intent(inout) :: config
  integer,                 intent(out)   :: err
  character(*),            intent(out)   :: message

  ! TOML table
  type(toml_table),        allocatable   :: table       ! root TOML table
  type(toml_table),        pointer       :: subtable    ! sub-table for a given section
  type(toml_key),          allocatable   :: sections(:) ! top-level sections
  type(toml_key),          allocatable   :: keys(:)     ! sub-table keys
  type(toml_error),        allocatable   :: error

  ! locals
  integer(i4b)       :: i,j,k
  character(len=256) :: cmessage
  logical(lgt)       :: mizuroute_config_present = .false.

  err = 0
  message = 'load_summa_config/'

  ! ----- initial checks with early return -----

  ! No configuration file is required unless mizuRoute is enabled
  if (len_trim(config_file) == 0) then
    if (mizuroute_active) then
      message = trim(message)//'mizuRoute is enabled but no TOML configuration file was specified; use -c <config_file>'
      err = 20
    endif
    return
  endif

  ! ----- load the root TOML table -----
  call toml_load(table, trim(config_file), error=error)

  if (allocated(error)) then
    message = "problem loading TOML file ['"//trim(config_file)//"']: "//trim(error%message)
    err = 10; return
  endif

  ! ----- get the top-level sections -----
  call table%get_keys(sections)
  if(.not.allocated(sections)) then
    message = trim(message)//"problem loading toml sections['"//trim(config%config_file)//"']"
    err=10; return
  endif

  ! ----- loop through sections -----
  do i = 1, size(sections)

    ! ----- load the TOML sub-table for the current section -----
    call get_value(table, trim(sections(i)%key), subtable, requested=.false.)
    if(.not.associated(subtable)) then
      message = trim(message)//"problem loading toml sub-sections['"//trim(config%config_file)//"']:"//trim(sections(i)%key)
      err=10; return
    endif

    ! ----- parameter dependencies are parsed as a complete section -----
    if(trim(sections(i)%key) == "parameter_dependencies")then

      call parse_parameter_dependencies(subtable, config, err, cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

      cycle

    endif

    ! ----- get keys for a given section (sub-table) -----
    call subtable%get_keys(keys)

    ! ----- loop through the sub-table -----
    do j = 1, size(keys)


      ! ----- parameter transformations are parsed as a complete sub-table -----
      if(trim(sections(i)%key) == "calibration" .and. &
         trim(keys(j)%key)     == "parameter_transformations")then
     
        call parse_parameter_transformations(subtable, config, err, cmessage)
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

        cycle
     
      endif

      ! select section
      select case (trim(sections(i)%key))

        ! ----- parse the summa sections of the TOML table -----
        case ("simulation", "summa_files", "observations", "calibration")

          call parse_summa_config(subtable,              &
                                  trim(sections(i)%key), &
                                  trim(keys(j)%key),     &
                                  config, err, cmessage)
          if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

        ! ----- parse the mizuRoute sections of the TOML table -----
        case ("mizuRoute", "hydrofabric", "remapping")

          mizuroute_config_present = .true.

          if (mizuroute_active) then
            call parse_mizuroute_config(subtable,              &
                                        trim(sections(i)%key), &
                                        trim(keys(j)%key),     &
                                        config, err, cmessage)
            if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
          endif

        ! ----- no other sections implemented -----
        case default
          message=trim(message)//'section ['//trim(sections(i)%key)//'] not implemented: remove section from TOML file'
          err=10; return

      end select   ! (cases within a desired section)

    end do  ! (looping through sub-sections)
  end do  ! (looping through sections)

  ! ----- checks -----

  ! Coupled mizuRoute requires the complete set of SUMMA GRUs because runoff from upstream
  ! GRUs may contribute to river reaches outside the selected SUMMA subdomain.
  if (config%use_mizuroute .and. iRunMode /= iRunModeFull) then
    message=trim(message)//'The -g subdomain option cannot be used with coupled mizuRoute because '// &
                           'the selected GRUs may not contain the complete upstream river network.'
    err=20; return
  endif

  ! Coupled mizuRoute requires an executable built with mizuRoute support because the
  ! required routing functionality is only available when mizuRoute is enabled at compile time.
  if (config%use_mizuroute .and. .not.mizuroute_active) then
    message=trim(message)//'mizuRoute was requested for this simulation, but this executable '// &
                           'was not built with mizuRoute support.'
    err=20; return
  endif

  ! mizuRoute configuration is ignored unless coupled mizuRoute is explicitly enabled
  ! for the simulation. Warn the user because the supplied configuration may indicate
  ! that they intended to run mizuRoute.
  if (mizuroute_config_present .and. .not.config%use_mizuroute) then
    write(iulog,*) 'WARNING: mizuRoute configuration was provided, but use_mizuroute is false.'
    write(iulog,*) '         Set simulation.use_mizuroute = true to run coupled mizuRoute, or remove the '
    write(iulog,*) '         mizuRoute, hydrofabric, and remapping sections if routing is not required.'
  endif

  ! ----- set default objective function settings -----

  ! set default objective-function metric
  if(.not.allocated(config%calib%metric))then
    config%calib%metric = 'kge'
    write(iulog,*) 'WARNING: objective metric not specified; using kge'
  endif
  
  ! set default obs transformation
  if(.not.allocated(config%calib%obs_transform))then
    config%calib%obs_transform = 'none'
    write(iulog,*) 'WARNING: observation transformation not specified; using none'
  endif

  end subroutine load_summa_config

  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  
  ! **************************************************************************************************
  ! Parse summa configuration.
  ! **************************************************************************************************
  
  subroutine parse_summa_config(subtable, section, key, config, ierr, message)
  
  use tomlf_all, only: toml_table, toml_array, toml_error, toml_key, toml_value ! data types
  use tomlf_all, only: toml_load, get_value, len                                ! procedures
  
  type(toml_table), pointer, intent(in)    :: subtable
  character(*),              intent(in)    :: section
  character(*),              intent(in)    :: key
  type(config_info),         intent(inout) :: config
  integer,                   intent(out)   :: ierr
  character(*),              intent(out)   :: message
 
  type(toml_array), pointer     :: param_list  ! sub-table for the list of parameters to vary
  character(len=256)            :: cmessage    ! error message from downwind routine
  integer(i4b)                  :: istat       ! error code

  associate(obs   => config%obs, &
            calib => config%calib)
  
  ierr    = 0
  message = 'parse_summa_config/'
  
  ! extract configuration values and populate structures
  
  select case(trim(section)//'.'//trim(key))
  
    ! ---- simulation: times ----
    case ("simulation.start_time"        ); call get_value(subtable, trim(key), config%start_time       , stat=istat)
    case ("simulation.end_time"          ); call get_value(subtable, trim(key), config%end_time         , stat=istat)
    case ("simulation.time_zone"         ); call get_value(subtable, trim(key), config%time_zone        , stat=istat)
 
    ! ---- simulation: settings ----
    case ("simulation.case_name"         ); call get_value(subtable, trim(key), config%case_name        , stat=istat)
    case ("simulation.use_mizuroute"     ); call get_value(subtable, trim(key), config%use_mizuroute    , stat=istat)
    case ("simulation.write_timeseries"  ); call get_value(subtable, trim(key), config%write_timeseries , stat=istat)

    ! ---- SUMMA files: paths ----
    case ("summa_files.settings_path"    ); call get_value(subtable, trim(key), config%settings_path    , stat=istat)
    case ("summa_files.forcing_path"     ); call get_value(subtable, trim(key), config%forcing_path     , stat=istat)
    case ("summa_files.output_path"      ); call get_value(subtable, trim(key), config%output_path      , stat=istat)
    case ("summa_files.state_path"       ); call get_value(subtable, trim(key), config%state_path       , stat=istat)
  
    ! ---- SUMMA files: model input files ----
    case ("summa_files.init_condition"   ); call get_value(subtable, trim(key), config%init_condition   , stat=istat)
    case ("summa_files.attributes"       ); call get_value(subtable, trim(key), config%attributes       , stat=istat)
    case ("summa_files.trial_params"     ); call get_value(subtable, trim(key), config%trial_params     , stat=istat)
    case ("summa_files.forcing_list"     ); call get_value(subtable, trim(key), config%forcing_list     , stat=istat)
    case ("summa_files.decisions"        ); call get_value(subtable, trim(key), config%decisions        , stat=istat)
    case ("summa_files.output_control"   ); call get_value(subtable, trim(key), config%output_control   , stat=istat)
  
    ! ---- SUMMA files: parameter files ----
    case ("summa_files.local_parameters" ); call get_value(subtable, trim(key), config%local_parameters , stat=istat)
    case ("summa_files.basin_parameters" ); call get_value(subtable, trim(key), config%basin_parameters , stat=istat)
  
    ! ---- SUMMA files: parameter tables ----
    case ("summa_files.vegetation_table" ); call get_value(subtable, trim(key), config%vegetation_table , stat=istat)
    case ("summa_files.soil_table"       ); call get_value(subtable, trim(key), config%soil_table       , stat=istat)
    case ("summa_files.general_table"    ); call get_value(subtable, trim(key), config%general_table    , stat=istat)
    case ("summa_files.noahmp_table"     ); call get_value(subtable, trim(key), config%noahmp_table     , stat=istat)
  
    ! ---- observations: filename ----
    case ("observations.obs_path"        ); call get_value(subtable, trim(key), obs%obs_path            , stat=istat)
    case ("observations.obs_file"        ); call get_value(subtable, trim(key), obs%obs_file            , stat=istat)
  
    ! ---- observations: variable names ----
    case ("observations.vname_obsflow"   ); call get_value(subtable, trim(key), obs%vname_obsflow       , stat=istat)
  
    ! ---- objective function: metrics  ----
    case ("calibration.metric"           ); call get_value(subtable, trim(key), calib%metric            , stat=istat)
    case ("calibration.obs_transform"    ); call get_value(subtable, trim(key), calib%obs_transform     , stat=istat)
  
    ! ---- objective function: calibration period  ----
    case ("calibration.start_date"       ); call get_value(subtable, trim(key), calib%start_date        , stat=istat)
    case ("calibration.end_date"         ); call get_value(subtable, trim(key), calib%end_date          , stat=istat)
  
    ! ---- objective function: list of parameters to modify  ----
    case ("calibration.param_list"       ); call get_value(subtable, trim(key), param_list              , stat=istat)
 
      if(istat == 0)then
        call parse_word_list(param_list, calib%param_list, ierr, cmessage)
        if(ierr/=0) then; message=trim(message)//trim(cmessage); return; endif
      endif   

    ! ---- objective function: flag to write aligned sim/obs time series  ----
    case ("calibration.write_aligned"    ); call get_value(subtable, trim(key), calib%write_aligned     , stat=istat)
    
    ! ---- default case (something in the table that is not specified above) -----
    case default
      message = trim(message)// "unexpected entry: section = "//trim(section)//"; sub-section = "//trim(key)
      ierr=20; return
  
  end select ! (select key/value pair based on lookup)
  
  ! ---- error checking -----
  if(istat /= 0)then
    message=trim(message)// "get_value error: section = "//trim(section)//"; sub-section = "//trim(key)
    ierr=20; return
  endif
  
  end associate
  
  end subroutine parse_summa_config
  
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------
  
  ! **************************************************************************************************
  ! Apply SUMMA configuration values parsed from TOML.
  !
  ! Values supplied in the TOML configuration overwrite values previously read
  ! from the legacy SUMMA file manager. Unspecified TOML values leave the legacy
  ! configuration unchanged.
  ! **************************************************************************************************
  
  subroutine apply_summa_config(config, err, message)
  
    USE summaFileManager, only: SIM_START_TM
    USE summaFileManager, only: SIM_END_TM
    USE summaFileManager, only: NC_TIME_ZONE
  
    USE summaFileManager, only: SETTINGS_PATH
    USE summaFileManager, only: FORCING_PATH
    USE summaFileManager, only: OUTPUT_PATH
    USE summaFileManager, only: STATE_PATH
  
    USE summaFileManager, only: M_DECISIONS
    USE summaFileManager, only: OUTPUT_CONTROL
    USE summaFileManager, only: LOCAL_ATTRIBUTES
    USE summaFileManager, only: LOCALPARAM_INFO
    USE summaFileManager, only: BASINPARAM_INFO
    USE summaFileManager, only: VEGPARM
    USE summaFileManager, only: SOILPARM
    USE summaFileManager, only: GENPARM
    USE summaFileManager, only: MPTABLE
    USE summaFileManager, only: FORCING_FILELIST
    USE summaFileManager, only: MODEL_INITCOND
    USE summaFileManager, only: PARAMETER_TRIAL
    USE summaFileManager, only: OUTPUT_PREFIX
  
    implicit none
  
    type(config_info), intent(in)  :: config
    integer(i4b),      intent(out) :: err
    character(*),      intent(out) :: message
  
    err = 0
    message = 'apply_summa_config/'
  
    ! ---- simulation ----
    if(allocated(config%case_name))  OUTPUT_PREFIX = trim(config%case_name)
    if(allocated(config%start_time)) SIM_START_TM  = trim(config%start_time)
    if(allocated(config%end_time))   SIM_END_TM    = trim(config%end_time)
    if(allocated(config%time_zone))  NC_TIME_ZONE  = trim(config%time_zone)
  
    ! ---- SUMMA files: paths ----
    if(allocated(config%settings_path)) SETTINGS_PATH = trim(config%settings_path)
    if(allocated(config%forcing_path))  FORCING_PATH  = trim(config%forcing_path)
    if(allocated(config%output_path))   OUTPUT_PATH   = trim(config%output_path)
    if(allocated(config%state_path))    STATE_PATH    = trim(config%state_path)
  
    ! ---- SUMMA files: model input files ----
    if(allocated(config%init_condition)) MODEL_INITCOND   = trim(config%init_condition)
    if(allocated(config%attributes))     LOCAL_ATTRIBUTES = trim(config%attributes)
    if(allocated(config%trial_params))   PARAMETER_TRIAL  = trim(config%trial_params)
    if(allocated(config%forcing_list))   FORCING_FILELIST = trim(config%forcing_list)
    if(allocated(config%decisions))      M_DECISIONS      = trim(config%decisions)
    if(allocated(config%output_control)) OUTPUT_CONTROL   = trim(config%output_control)
  
    ! ---- SUMMA files: parameter files ----
    if(allocated(config%local_parameters)) LOCALPARAM_INFO = trim(config%local_parameters)
    if(allocated(config%basin_parameters)) BASINPARAM_INFO = trim(config%basin_parameters)
  
    ! ---- SUMMA files: parameter tables ----
    if(allocated(config%vegetation_table)) VEGPARM = trim(config%vegetation_table)
    if(allocated(config%soil_table))       SOILPARM = trim(config%soil_table)
    if(allocated(config%general_table))    GENPARM  = trim(config%general_table)
    if(allocated(config%noahmp_table))     MPTABLE  = trim(config%noahmp_table)
  
  end subroutine apply_summa_config

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
      allocate(words(0))
      return
    endif
  
    nwords = len(word_list)
  
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


  ! **************************************************************************************************
  ! Parse parameter transformations.
  !
  ! Reads parameter transformations from a TOML key-value table where each key is a parameter name
  ! and each value defines the transformation used for that parameter during parameter search.
  ! **************************************************************************************************
  
  subroutine parse_parameter_transformations(calib_table, config, ierr, message)
  
    use tomlf_all, only: toml_table, toml_key, get_value
  
    implicit none
  
    type(toml_table), pointer, intent(in)    :: calib_table
    type(config_info),         intent(inout) :: config
    integer(i4b),              intent(out)   :: ierr
    character(*),              intent(out)   :: message
  
    type(toml_table), pointer    :: transform_table
    type(toml_key), allocatable  :: keys(:)
  
    integer(i4b) :: i
    integer(i4b) :: istat
  
    character(len=:), allocatable :: transform
  
    ierr = 0
    message = 'parse_parameter_transformations/'
  
    ! get parameter transformation sub-table
    call get_value(calib_table, 'parameter_transformations', &
                   transform_table, stat=istat)
  
    if(istat/=0 .or. .not.associated(transform_table))then
      message=trim(message)//'unable to read parameter_transformations table'
      ierr=20; return
    endif
  
    ! get parameter names from table keys
    call transform_table%get_keys(keys)
  
    if(.not.allocated(keys))then
      allocate(config%calib%param_transform(0))
      return
    endif
  
    ! allocate transformation information
    allocate(config%calib%param_transform(size(keys)),stat=ierr)
    if(ierr/=0)then
      message=trim(message)//'unable to allocate parameter transformations'
      return
    endif
  
    ! read parameter -> transformation mappings
    do i=1,size(keys)
  
      ! check parameter-name length
      if(len_trim(keys(i)%key) > &
         len(config%calib%param_transform(i)%name))then
  
        write(message,'(A,I0)') trim(message)// &
          'parameter name exceeds maximum character length, i = ',i
        ierr=20; return
      endif
  
      config%calib%param_transform(i)%name = trim(keys(i)%key)
  
      call get_value(transform_table,trim(keys(i)%key),transform,stat=istat)
  
      if(istat/=0)then
        message=trim(message)//'unable to read transformation for parameter: '// &
                trim(keys(i)%key)
        ierr=20; return
      endif
  
      ! check transformation-name length
      if(len_trim(transform) > &
         len(config%calib%param_transform(i)%transformation))then
  
        message=trim(message)//'transformation name exceeds maximum character length for parameter: '// &
                trim(keys(i)%key)
        ierr=20; return
      endif
  
      config%calib%param_transform(i)%transformation = trim(transform)
  
    enddo
  
  end subroutine parse_parameter_transformations


  ! **************************************************************************************************
  ! Parse parameter dependency configuration.
  !
  ! Reads ordered parameter constraints from the TOML configuration. Each constraint defines
  ! an ordered list of parameters and the minimum gap between adjacent parameters as a fraction
  ! of the total parameter range.
  ! **************************************************************************************************
  
  subroutine parse_parameter_dependencies(subtable, config, ierr, message)
  
    use tomlf_all, only: toml_table, toml_array, get_value, len
  
    implicit none
  
    type(toml_table), pointer, intent(in)    :: subtable
    type(config_info),         intent(inout) :: config
    integer(i4b),              intent(out)   :: ierr
    character(*),              intent(out)   :: message
  
    type(toml_array), pointer :: ordered
    type(toml_table), pointer :: constraint
    type(toml_array), pointer :: param_list
  
    integer(i4b) :: i
    integer(i4b) :: istat
    integer(i4b) :: nconstraints
  
    character(len=256) :: cmessage
  
    ierr = 0
    message = 'parse_parameter_dependencies/'
  
    ! get array of ordered constraints
    call get_value(subtable, 'ordered', ordered, requested=.false., stat=istat)
  
    if(.not.associated(ordered)) return
  
    nconstraints = len(ordered)
  
    ! allocate constraint structures
    allocate(config%calib%ordered(nconstraints), stat=ierr)
    if(ierr/=0)then
      message=trim(message)//'unable to allocate ordered parameter constraints'
      return
    endif
  
    ! parse each ordered constraint
    do i=1,nconstraints
  
      call get_value(ordered, i, constraint, stat=istat)
      if(istat/=0 .or. .not.associated(constraint))then
        write(message,'(A,I0)') trim(message)// &
          'unable to read ordered constraint, i = ',i
        ierr=20; return
      endif
  
      ! parameter list
      call get_value(constraint, 'parameters', param_list, stat=istat)
      if(istat/=0 .or. .not.associated(param_list))then
        write(message,'(A,I0)') trim(message)// &
          'parameter list not defined for ordered constraint, i = ',i
        ierr=20; return
      endif
  
      call parse_word_list(param_list,                         &
                           config%calib%ordered(i)%parameters, &
                           ierr,cmessage)
      if(ierr/=0)then
        message=trim(message)//trim(cmessage)
        return
      endif
  
      ! gap fraction
      call get_value(constraint, 'gap_fraction', &
                     config%calib%ordered(i)%gap_fraction, stat=istat)
  
      if(istat/=0)then
        write(message,'(A,I0)') trim(message)// &
          'gap_fraction not defined for ordered constraint, i = ',i
        ierr=20; return
      endif
  
    enddo
  
  end subroutine parse_parameter_dependencies


end module summa_config

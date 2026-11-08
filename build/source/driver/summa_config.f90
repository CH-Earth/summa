module summa_config

use build_options, only: mizuroute_active

USE nr_type
USE summa_type, only: config_info       ! summa configuation info
USE data_types, only: calib_info        ! calibration info
USE data_types, only: obs_fileinfo      ! obs file info

USE globalData, only: iulog             ! I/O unit for logging messages
USE globalData, only: startGRU          ! index of the first GRU of the run domain
USE globalData, only: integerMissing    ! missing integer

USE calibration_config, only: parse_word_list
USE calibration_config, only: parse_nsga2_settings
USE calibration_config, only: parse_calibration_targets
USE calibration_config, only: parse_calibration_parameters
USE calibration_config, only: parse_parameter_dependencies
USE calibration_config, only: parse_calibration_settings

USE calibration_config, only: resolve_parameter_ratio_references
USE calibration_config, only: finalize_calibration_settings

USE parameter_search,   only: param_spatial_multiplier,param_ratio

#ifdef MIZUROUTE_ACTIVE
USE mizuroute_config, ONLY: parse_mizuroute_config
#endif

implicit none
private

public :: read_manifest
public :: read_summa_config 

contains

  ! **************************************************************************************************
  ! Read a multi-case SUMMA run manifest.
  ! Loads the TOML manifest, extracts the [multi_case] table, and parses the run-level configuration
  ! used to define and distribute independent SUMMA cases.
  ! **************************************************************************************************
  subroutine read_manifest(config,err,message)
    USE tomlf_all, only: toml_table,toml_error
    USE tomlf_all, only: toml_load,get_value
    implicit none
  
    type(config_info), intent(inout) :: config
    integer(i4b),      intent(out)   :: err
    character(*),      intent(out)   :: message
    type(toml_table), allocatable :: table
    type(toml_table), pointer     :: subtable
    type(toml_error), allocatable :: toml_err
    character(len=:), allocatable :: manifest_temp
    integer(i4b)         :: istat
    character(len=1024)  :: cmessage
  
    err=0
    message='read_manifest/'

    ! construct the full path to the manifest
    if(config%manifest_file(1:1) == '/')then
     manifest_temp = trim(config%manifest_file)
    else
     manifest_temp = trim(config%cwd)//'/'//trim(config%manifest_file)
    endif
    config%manifest_file = trim(manifest_temp)

    ! load the TOML manifest
    call toml_load(table,trim(config%manifest_file),error=toml_err)
    if(allocated(toml_err))then
      message=trim(message)//"problem loading manifest ['"// &
              trim(config%manifest_file)//"']: "//trim(toml_err%message)
      err=20
      return
    endif
  
    ! extract the multi-case configuration table
    call get_value(table,'multi_case',subtable,stat=istat)
    if(istat/=0 .or. .not.associated(subtable))then
      message=trim(message)//'manifest does not contain [multi_case]'
      err=20
      return
    endif
  
    ! parse the multi-case configuration
    call parse_manifest(subtable,config,err,cmessage)
    if(err/=0)then
      message=trim(message)//trim(cmessage)
      return
    endif
  
  end subroutine read_manifest

  ! **************************************************************************************************
  ! Read SUMMA TOML configuration.
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
    character(len=1024) :: cmessage
  
    err = 0
    message = 'read_summa_config/'

    ! load configuration values from the TOML file
    call load_summa_config(config_file, config, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! apply command-line overrides
    if(allocated(config%home_path_override)) config%home_path = config%home_path_override

    ! expand case-specific placeholders
    call expand_summa_config(config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! apply TOML values that supersede legacy file-manager settings
    call apply_summa_config(config, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  
  end subroutine read_summa_config

  ! --------------------------------------------------------------------------------------------------
  ! ---- PRIVATE SUBROUTINES -------------------------------------------------------------------------
  ! **************************************************************************************************
  ! Load SUMMA configuration from a TOML file.
  ! Loads the TOML configuration file and delegates parsing to component-specific readers.
  ! TOML configuration may include general SUMMA settings as well as configuration for optional
  ! components such as mizuRoute. Settings provided in the TOML file take precedence over
  ! corresponding values previously read from the legacy SUMMA file manager.
  ! **************************************************************************************************
  subroutine load_summa_config(config_file, config, err, message)
    USE tomlf_all, only: toml_table, toml_array, toml_error, toml_key, toml_value ! data types
    USE tomlf_all, only: toml_load, get_value, len                                ! procedures
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
    integer(i4b)        :: i,j
    character(len=1024) :: cmessage
    logical(lgt)        :: mizuroute_config_present
    logical(lgt)        :: hasObs                      ! .true. if streamflow observations are configured

    err = 0
    message = 'load_summa_config/'
    mizuroute_config_present=.false.

    ! Reset parameter dependencies when reading a new case
    config%calib = calib_info()
    config%obs   = obs_fileinfo()

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
      
      nullify(subtable)
      
      ! ----- load the TOML sub-table for the current section -----
      call get_value(table, trim(sections(i)%key), subtable, requested=.false.)
      if(.not.associated(subtable)) then
        message = trim(message)//"problem loading toml sub-sections['"//trim(config%config_file)//"']:"//trim(sections(i)%key)
        err=10; return
      endif

      ! ----- parameter dependencies are parsed as a complete section -----
      if(trim(sections(i)%key) == "parameter_dependencies")then
        call parse_parameter_dependencies(subtable,config%calib,err,cmessage)
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
        cycle
      endif

      ! ----- parse calibration configuration -----
      if(trim(sections(i)%key) == 'calibration')then

        ! General calibration settings
        call parse_calibration_settings(subtable,config%calib,err,cmessage)
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

        ! Calibration parameter specifications
        if(subtable%has_key('parameters')          .or. &
           subtable%has_key('spatial_multipliers') .or. &
           subtable%has_key('parameter_ratios'))then
          call parse_calibration_parameters(subtable,config%calib,err,cmessage)
          if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
        endif

        ! NSGA-II optimization settings
        if(subtable%has_key('nsga2'))then
          call parse_nsga2_settings(subtable,config%calib,err,cmessage)
          if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
        endif

        ! Calibration targets are parsed after all sections have been read.
        cycle

      endif

      ! ----- get keys for a given section (sub-table) -----
      call subtable%get_keys(keys)
      if(.not.allocated(keys)) cycle

      ! ----- loop through the sub-table -----
      do j = 1, size(keys)

        ! select section
        select case (trim(sections(i)%key))

          ! ----- parse the summa sections of the TOML table -----
          case ("simulation", "summa_files", "observations")
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

    ! Resolve calibration parameter dependencies.
    call resolve_parameter_ratio_references(config%calib,err,cmessage)
    if(err/=0)then
      message=trim(message)//trim(cmessage)
      return
    endif

    ! -----------------------------------------------------------------------------------------------
    ! Parse calibration targets after the remaining configuration settings.
    ! This ensures that observation and calibration defaults are available
    ! regardless of the order of sections and keys in the TOML file.
    ! -----------------------------------------------------------------------------------------------
    nullify(subtable)
    call get_value(table,'calibration',subtable,requested=.false.)

    if(associated(subtable))then
      if(subtable%has_key('target'))then
        call parse_calibration_targets(subtable,config%calib,          &
                                       config%obs%obs_path,            &
                                       config%obs%vname_obsflow,err,cmessage)
        if(err/=0)then
          message=trim(message)//trim(cmessage)
          return
        endif
      endif
    endif

    ! Print the completed calibration parameter registry.
    call print_calibration_parameters(config)
  
    ! ----- checks -----
    ! Coupled mizuRoute requires the complete set of SUMMA GRUs because runoff from upstream
    ! GRUs may contribute to river reaches outside the selected SUMMA subdomain.  It is the GRU
    ! range that says whether the run is the whole domain, not the run mode: a host owning the
    ! command line runs every GRU of the attributes file, but under iRunModeGRU for the file suffix.
    if (config%use_mizuroute .and. .not.(startGRU == 1 .and. config%nGRU_user == integerMissing)) then
      message=trim(message)//'A GRU subdomain cannot be used with coupled mizuRoute because '// &
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

    ! ----- set calibration defaults and derived metric settings -----

    ! NOTE: warn only when observations are configured. A configuration file is also used for a
    !       plain run or for mizuRoute alone, and neither has an objective function to evaluate.
    hasObs = .false.
    if(allocated(config%obs%obs_file))then
      if(len_trim(config%obs%obs_file) > 0) hasObs = .true.
    endif
    if(allocated(config%calib%targets))then
      if(size(config%calib%targets) > 0) hasObs = .true.
    endif

    ! ----- finalize generic calibration settings -----

    ! Warn when observations are configured but calibration defaults are used.
    if(hasObs .and. config%is_case_root)then

      if(.not.allocated(config%calib%metric))then
        write(iulog,*) 'WARNING: objective metric not specified; using kge'
      endif

      if(.not.allocated(config%calib%obs_transform))then
        write(iulog,*) 'WARNING: observation transformation not specified; using none'
      endif

    endif

    call finalize_calibration_settings(config%calib,hasObs,err,cmessage)
    if(err/=0)then
      message=trim(message)//trim(cmessage)
      return
    endif

    ! ----- provide backward compatibility with legacy SUMMA streamflow observations -----

    if(hasObs)then

      call set_default_calibration_target(config,err,cmessage)
      if(err/=0)then
        message=trim(message)//trim(cmessage)
        return
      endif

    endif

  end subroutine load_summa_config
  
  ! **************************************************************************************************
  ! Print calibration parameter information.
  ! **************************************************************************************************
  subroutine print_calibration_parameters(config)
    USE parameter_search, only: param_physical,param_spatial_multiplier,param_ratio
    implicit none

    type(config_info), intent(in) :: config
    integer(i4b) :: i
    character(len=16) :: lower_text,upper_text
    character(len=18) :: type_text

    if(.not.config%is_case_root) return
    if(.not.allocated(config%calib%param_list)) return

    write(iulog,'(/,A)') 'Calibration parameters:'
    write(iulog,'(5X,A24,2X,A18,2X,A16,2X,A16,2X,A16,2X,A)') &
      'Parameter','Type','Lower bound','Upper bound','Transformation','Reference'

    do i=1,size(config%calib%param_list)

      associate(param => config%calib%param_list(i))

        ! Parameterization
        select case(param%param_type)
        case (param_physical)
          type_text='physical'
        case (param_spatial_multiplier)
          type_text='spatial multiplier'
        case (param_ratio)
          type_text='parameter ratio'
        case default
          type_text='unknown'
        end select

        ! Calibration bounds
        lower_text='native'
        upper_text='native'

        if(param%has_lower) write(lower_text,'(ES16.8)') param%lower_bound
        if(param%has_upper) write(upper_text,'(ES16.8)') param%upper_bound

        ! Parameter information
        write(iulog,'(I3,2X,A24,2X,A18,2X,A16,2X,A16,2X,A16,2X,A)') &
          i,trim(param%name),trim(type_text),lower_text,upper_text, &
          trim(param%transformation),trim(param%reference_name)

      end associate

    enddo

  end subroutine print_calibration_parameters


  ! **************************************************************************************************
  ! Express a configuration that names no calibration target as a single-target calibration.
  !
  ! The legacy configuration carries one observed streamflow series in [observations] and one metric
  ! in [calibration].  That is a calibration with exactly one target, so it is stored as one, and the
  ! objective code needs only the target list.
  ! **************************************************************************************************
  subroutine set_default_calibration_target(config, ierr, message)
    implicit none
    type(config_info), intent(inout) :: config
    integer(i4b),      intent(out)   :: ierr
    character(*),      intent(out)   :: message

    ierr = 0
    message = 'set_default_calibration_target/'

    ! an explicit target list wins
    if(allocated(config%calib%targets))then
      if(size(config%calib%targets) > 0) return
      deallocate(config%calib%targets)
    endif

    ! without an observation file there is nothing to build a target from
    if(.not.allocated(config%obs%obs_file)) return
    if(len_trim(config%obs%obs_file) == 0) return

    allocate(config%calib%targets(1), stat=ierr)
    if(ierr/=0)then
      message=trim(message)//'unable to allocate the default calibration target'
      return
    endif

    config%calib%targets(1)%name     = 'streamflow'
    config%calib%targets(1)%variable = 'streamflow'
    if(allocated(config%obs%obs_path)) config%calib%targets(1)%obs_path = trim(config%obs%obs_path)
    config%calib%targets(1)%obs_file = trim(config%obs%obs_file)
    if(allocated(config%obs%vname_obsflow))then
      config%calib%targets(1)%vname_obs = trim(config%obs%vname_obsflow)
    else
      config%calib%targets(1)%vname_obs = 'q_obs'
    endif
    config%calib%targets(1)%metric        = trim(config%calib%metric)
    config%calib%targets(1)%obs_transform = trim(config%calib%obs_transform)
    config%calib%targets(1)%boxcox_lambda = config%calib%boxcox_lambda
    config%calib%targets(1)%weight        = 1.0_rkind

  end subroutine set_default_calibration_target

  ! **************************************************************************************************
  ! Parse summa configuration.
  ! **************************************************************************************************
  subroutine parse_summa_config(subtable, section, key, config, ierr, message)
    USE tomlf_all, only: toml_table, toml_array, toml_error, toml_key, toml_value ! data types
    USE tomlf_all, only: toml_load, get_value, len                                ! procedures
    
    type(toml_table), pointer, intent(in)    :: subtable
    character(*),              intent(in)    :: section
    character(*),              intent(in)    :: key
    type(config_info),         intent(inout) :: config
    integer,                   intent(out)   :: ierr
    character(*),              intent(out)   :: message
    type(toml_array), pointer     :: param_list  ! sub-table for the list of parameters to vary
    character(len=:), allocatable :: cvalue      ! a string value read into a fixed-length field
    character(len=1024)           :: cmessage    ! error message from downwind routine
    integer(i4b)                  :: istat       ! error code
    integer(i4b)                  :: i

    associate(obs   => config%obs)
    
    ierr    = 0
    message = 'parse_summa_config/'
    
    ! extract configuration values and populate structures
    
    select case(trim(section)//'.'//trim(key))
    
      ! ---- simulation: times ----
      case ("simulation.start_time"        ); call get_value(subtable, trim(key), config%start_time       , stat=istat)
      case ("simulation.end_time"          ); call get_value(subtable, trim(key), config%end_time         , stat=istat)
      case ("simulation.time_zone"         ); call get_value(subtable, trim(key), config%time_zone        , stat=istat)
    
      ! ---- simulation: settings ----
      case ("simulation.home_path"         ); call get_value(subtable, trim(key), config%home_path        , stat=istat)
      case ("simulation.basin_dir"         ); call get_value(subtable, trim(key), config%basin_dir        , stat=istat)
      case ("simulation.work_path"         ); call get_value(subtable, trim(key), config%work_path        , stat=istat)
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
  
  ! **************************************************************************************************
  ! Parse the multi-case configuration from a SUMMA run manifest.
  ! Extracts settings from the [multi_case] TOML table, including the configuration template,
  ! case names, and the number of concurrent cases assigned to each compute node.
  ! **************************************************************************************************
  subroutine parse_manifest(subtable,config,err,message)
    USE tomlf_all, only: toml_table,toml_key,toml_array,get_value
    implicit none
  
    type(toml_table), pointer, intent(in)    :: subtable
    type(config_info),         intent(inout) :: config
    integer(i4b),              intent(out)   :: err
    character(*),              intent(out)   :: message
    type(toml_key), allocatable :: keys(:)
    type(toml_array), pointer   :: case_names
    integer(i4b)        :: i,idx,istat
    character(len=1024) :: key,cmessage
    logical(lgt)        :: have_template_path
  
    err=0
    message='parse_manifest/'
    have_template_path = .false.
    call subtable%get_keys(keys)
  
    do i=1,size(keys)
      istat=0 
      key='multi_case.'//trim(keys(i)%key)
  
      select case(trim(key))
  
        ! ---- multi-case configuration ----
        case ("multi_case.cases_per_node"   ); call get_value(subtable,trim(keys(i)%key),config%cases_per_node    , stat=istat)
        case ("multi_case.n_runs"           ); call get_value(subtable,trim(keys(i)%key),config%n_runs            , stat=istat)
        case ("multi_case.experiment"       ); call get_value(subtable,trim(keys(i)%key),config%experiment        , stat=istat)

        ! ----- path/name of toml template (template_path is optional: default is the manifest directory) -----
        case ("multi_case.template_path"    ); call get_value(subtable,trim(keys(i)%key),config%template_path     , stat=istat)
                                               have_template_path = (istat==0)
        case ("multi_case.template_file"    ); call get_value(subtable,trim(keys(i)%key),config%template_file     , stat=istat)

        ! ---- case names ----
        case ("multi_case.case_names")
          call get_value(subtable,trim(keys(i)%key),case_names,stat=istat)
          if(istat==0) then
            call parse_word_list(case_names,config%case_names,err,cmessage)
            if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
          endif

        case default
          message=trim(message)//"unknown manifest option '"//trim(key)//"'"
          err=20; return
  
      end select
  
      if(istat/=0)then
        message=trim(message)//"unable to read manifest option '"//trim(key)//"'"
        err=20; return
      endif
  
    enddo

    ! ----- set the template path (template file assumed to be in the same directory as the manifest) -----

    idx = max( scan(trim(config%manifest_file), '/',       back=.true.), & ! mac/linux
               scan(trim(config%manifest_file), achar(92), back=.true.)  ) ! windows

    if(.not.have_template_path)then
      if(idx > 0)then
        config%template_path = config%manifest_file(:idx)
      else
        config%template_path = ''
      endif
    endif

    ! ----- validate required manifest settings -----

    if(.not.allocated(config%template_file))then
      message=trim(message)//'template_file is not defined in the multi-case manifest'
      err=20; return
    endif

    if(.not.allocated(config%case_names))then
      message=trim(message)//'case_names are not defined in the multi-case manifest'
      err=20; return
    endif

    if(size(config%case_names)==0)then
      message=trim(message)//'case_names empty in the multi-case manifest'
      err=20; return
    endif

    if(config%cases_per_node < 1)then
      message=trim(message)//'cases_per_node must be greater than zero'
      err=20; return
    endif

    if(config%n_runs < 1)then
      message=trim(message)//'n_runs must be greater than zero'
      err=20; return
    endif

  end subroutine parse_manifest

  ! **************************************************************************************************
  ! Expand case-specific placeholders in the SUMMA configuration.
  ! Resolves the case-path template using the current case name, then expands supported placeholders
  ! in configuration strings before the values are applied to SUMMA data structures.
  ! **************************************************************************************************
  subroutine expand_summa_config(config,err,message)
    implicit none
  
    type(config_info), intent(inout) :: config              ! SUMMA configuration information
    integer(i4b),      intent(out)   :: err                 ! error code
    character(*),      intent(out)   :: message             ! error message
    logical(lgt), parameter :: isPrint=.false.              ! temporary diagnostic output
    integer(i4b) :: iTarget                                 ! calibration target index
    character(len=256) :: cmessage                          ! message returned by called routines

    err=0
    message='expand_summa_config/'
  
    ! ----- use manifest values to populate case name -----
    if(allocated(config%manifest_file))then
      if(.not.allocated(config%manifest_casename))then
        message=trim(message)//'manifest_casename has not been assigned for the current case'
        err=20; return
      endif
      config%case_name=trim(config%manifest_casename)
    endif

    ! ----- resolve base path templates first -----
    ! home_path must be fully resolved because other paths may depend on it
    call expand_config_string(config%home_path,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    
    ! basin_dir may depend on home_path and case_name
    call expand_config_string(config%basin_dir,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    
    ! work_path may depend on home_path and basin_dir
    call expand_config_string(config%work_path,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! ----- check that template variables do not contain unresolved placeholders -----
    if(allocated(config%case_name))then
      if(index(config%case_name,'{')>0 .or. index(config%case_name,'}')>0)then
        message=trim(message)//"case_name contains an unresolved template placeholder: '"// &
                trim(config%case_name)//"'"
        err=20; return
      endif
    endif

    if(allocated(config%home_path))then
      if(index(config%home_path,'{')>0 .or. index(config%home_path,'}')>0)then
        message=trim(message)//"home_path contains an unresolved template placeholder: '"// &
                trim(config%home_path)//"'"
        err=20; return
      endif
    endif

    if(allocated(config%basin_dir))then
      if(index(config%basin_dir,'{')>0 .or. index(config%basin_dir,'}')>0)then
        message=trim(message)//"basin_dir contains an unresolved template placeholder: '"// &
                trim(config%basin_dir)//"'"
        err=20; return
      endif
    endif

    if(allocated(config%work_path))then
      if(index(config%work_path,'{')>0 .or. index(config%work_path,'}')>0)then
        message=trim(message)//"work_path contains an unresolved template placeholder: '"// &
                trim(config%work_path)//"'"
        err=20; return
      endif
    endif

    ! ---- SUMMA paths ----
    call expand_config_string(config%settings_path,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  
    call expand_config_string(config%forcing_path,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  
    call expand_config_string(config%output_path,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  
    call expand_config_string(config%state_path,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! ---- SUMMA files ----
    call expand_config_string(config%init_condition,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    
    call expand_config_string(config%attributes,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    
    call expand_config_string(config%trial_params,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    
    call expand_config_string(config%forcing_list,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! ---- Common SUMMA configuration files ----
    call expand_config_string(config%decisions,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    call expand_config_string(config%output_control,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! ---- Common SUMMA parameter files ----
    call expand_config_string(config%local_parameters,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    call expand_config_string(config%basin_parameters,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! ---- Common SUMMA parameter tables ----
    call expand_config_string(config%vegetation_table,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    call expand_config_string(config%soil_table,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    call expand_config_string(config%general_table,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    call expand_config_string(config%noahmp_table,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

#ifdef MIZUROUTE_ACTIVE
    ! ---- mizuRoute paths ----
    call expand_config_string(config%mizu_info%mrout%namelist_path,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! ---- hydrofabric paths ----
    call expand_config_string(config%mizu_info%ntopo%hfabric_path,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! ---- remapping paths ----
    call expand_config_string(config%mizu_info%remap%remap_path,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
#endif
  
    ! ---- observation paths and filenames ----
    call expand_config_string(config%obs%obs_path,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    call expand_config_string(config%obs%obs_file,config,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! ---- per-target observation paths and filenames ----
    ! A target names its own observation file, so it needs the same placeholder expansion the
    ! [observations] section gets; a multi-case run relies on it to reach each case's own data.
    if(allocated(config%calib%targets))then
      do iTarget=1,size(config%calib%targets)
        call expand_config_string(config%calib%targets(iTarget)%obs_path,config,err,cmessage)
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

        call expand_config_string(config%calib%targets(iTarget)%obs_file,config,err,cmessage)
        if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
      enddo
    endif

    ! temporary diagnostic output
    if(isPrint)then
      write(*,'(A)') 'Expanded SUMMA configuration:'
      write(*,'(A)') '  case_name      = '//trim(config%case_name)
      write(*,'(A)') '  basin_dir      = '//trim(config%basin_dir)
      write(*,'(A)') '  settings_path  = '//trim(config%settings_path)
      write(*,'(A)') '  forcing_path   = '//trim(config%forcing_path)
      write(*,'(A)') '  output_path    = '//trim(config%output_path)
      write(*,'(A)') '  state_path     = '//trim(config%state_path)
      write(*,'(A)') '  init_condition = '//trim(config%init_condition)
      write(*,'(A)') '  attributes     = '//trim(config%attributes)
      write(*,'(A)') '  trial_params   = '//trim(config%trial_params)
      write(*,'(A)') '  forcing_list   = '//trim(config%forcing_list)
      write(*,'(A)') '  decisions      = '//trim(config%decisions)
      write(*,'(A)') '  local_params   = '//trim(config%local_parameters)
#ifdef MIZUROUTE_ACTIVE
      write(*,'(A)') '  namelist_path  = '//trim(config%mizu_info%mrout%namelist_path)
      write(*,'(A)') '  hfabric_path   = '//trim(config%mizu_info%ntopo%hfabric_path)
      write(*,'(A)') '  remap_path     = '//trim(config%mizu_info%remap%remap_path)
#endif
      write(*,'(A)') '  obs_path       = '//trim(config%obs%obs_path)
      write(*,'(A)') '  obs_file       = '//trim(config%obs%obs_file)
      write(*,*)
    endif

  end subroutine expand_summa_config

  ! **************************************************************************************************
  ! Apply SUMMA configuration values parsed from TOML.
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

  ! --------------------------------------------------------------------------------------------------
  ! ---- HELPERS -------------------------------------------------------------------------------------
  ! --------------------------------------------------------------------------------------------------

  ! **************************************************************************************************
  ! Expand supported placeholders in a configuration string.
  ! Replaces case-specific template variables with their resolved values. Unallocated configuration
  ! strings are ignored.
  ! **************************************************************************************************
  subroutine expand_config_string(value,config,err,message)
    implicit none
  
    character(len=:), allocatable, intent(inout)  :: value   ! configuration string to expand
    type(config_info),              intent(in)    :: config  ! SUMMA configuration information
    integer(i4b),                   intent(out)   :: err     ! error code
    character(*),                   intent(out)   :: message ! error message
    character(len=1024)                           :: cmessage

    err=0
    message='expand_config_string/'
  
    ! nothing to expand when this configuration value was not provided
    if(.not.allocated(value)) return

    ! expand the resolved home path
    if(index(value,'{home}')>0)then
      if(.not.allocated(config%home_path))then
        message=trim(message)//"placeholder '{home}' used but home_path is not defined"
        err=20; return
      endif
      call replace_string(value,'{home}',trim(config%home_path),err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    endif

    ! expand the resolved case name
    if(index(value,'{case_name}')>0)then
      if(.not.allocated(config%case_name))then
        message=trim(message)//"placeholder '{case_name}' used but case_name is not defined"
        err=20; return
      endif
      call replace_string(value,'{case_name}',trim(config%case_name),err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    endif

    ! expand the resolved basin directory
    if(index(value,'{basin_dir}')>0)then
      if(.not.allocated(config%basin_dir))then
        message=trim(message)//"placeholder '{basin_dir}' used but basin_dir is not defined"
        err=20; return
      endif
      call replace_string(value,'{basin_dir}',trim(config%basin_dir),err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    endif

    ! expand the resolved work path
    if(index(value,'{work_path}')>0)then
      if(.not.allocated(config%work_path))then
        message=trim(message)//"placeholder '{work_path}' used but work_path is not defined"
        err=20; return
      endif
      call replace_string(value,'{work_path}',trim(config%work_path),err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    endif
 
    ! expand the experiment identifier.
    if(index(value,'{experiment}')>0)then
      if(.not.allocated(config%experiment))then
        message=trim(message)//"placeholder '{experiment}' used but experiment is not defined"
        err=20; return
      endif
      call replace_string(value,'{experiment}',trim(config%experiment),err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    endif

  end subroutine expand_config_string

  ! **************************************************************************************************
  ! Replace all occurrences of a substring within a string.
  ! **************************************************************************************************
  subroutine replace_string(string,pattern,replacement,err,message)
    implicit none
  
    character(len=:), allocatable, intent(inout) :: string
    character(*),                  intent(in)    :: pattern
    character(*),                  intent(in)    :: replacement
    integer(i4b),                  intent(out)   :: err
    character(*),                  intent(out)   :: message
    integer(i4b), parameter :: maxTry=100
    integer(i4b) :: iTry
    integer(i4b) :: ipos
  
    err=0
    message='replace_string/'
  
    ! ignore empty search patterns
    if(len(pattern)==0) return
 
    ! prevent substitutions that reproduce the search pattern
    if(index(replacement,pattern)>0)then
      message=trim(message)//"invalid template expansion: placeholder '"//trim(pattern)// &
              "' resolves to '"//trim(replacement)//"' while expanding '"//trim(string)//"'"
      err=20; return
    endif
  
    ! replace all occurrences of the search pattern
    do iTry=1,maxTry
      ipos=index(string,pattern)
      if(ipos==0) return
      string=string(:ipos-1)//trim(replacement)// &
             string(ipos+len(pattern):)
    enddo
  
    ! maximum number of substitutions exceeded
    write(message,'(A,A,A,I0,A,A,A)') &
      "replace_string/maximum number of replacements for pattern '", &
      trim(pattern),"' exceeded (",maxTry,") in string '",trim(string),"'"
    err=20
  
  end subroutine replace_string

end module summa_config

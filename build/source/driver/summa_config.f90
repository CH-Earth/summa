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

public :: load_summa_config

! NOTE:
! This module provides the top-level interface for TOML configuration.
! It loads the TOML file and delegates parsing to component-specific readers.
!
! At present, TOML is used only for configuration of optional components
! such as mizuRoute. SUMMA configuration continues to use the existing
! ASCII file-manager infrastructure. Support for reading general SUMMA
! configuration from TOML may be added in the future.

! -------------------------------------------------------------------------------------
! -------------------------------------------------------------------------------------

contains

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
  integer(i4b)       :: i, j
  character(len=256) :: cmessage

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
    message = trim(message)//"problem loading toml sections['"//trim(config%summaConfigFile)//"']"
    err=10; return
  endif

  ! ----- loop through sections -----
  do i = 1, size(sections)

    ! ----- load the TOML sub-table for the current section -----
    call get_value(table, trim(sections(i)%key), subtable, requested=.false.)
    if(.not.associated(subtable)) then
      message = trim(message)//"problem loading toml sub-sections['"//trim(config%summaConfigFile)//"']:"//trim(sections(i)%key)
      err=10; return
    endif

    ! ----- get keys for a given section (sub-table) -----
    call subtable%get_keys(keys)

    ! ----- loop through the sub-table -----
    do j = 1, size(keys)

      ! select section
      select case (trim(sections(i)%key))

        ! ----- parse the mizuRoute sections of the TOML table -----
        case ("mizuRoute", "hydrofabric", "remapping")

          if (mizuroute_active) then
            call parse_mizuroute_config(subtable,              &
                                        trim(sections(i)%key), &
                                        trim(keys(j)%key),     &
                                        config, err, cmessage)
            if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
          endif

        ! ----- parse the objective function sections of the TOML table -----
        case ("observations", "objective")

          call parse_objective_config(subtable,              &
                                      trim(sections(i)%key), &
                                      trim(keys(j)%key),     &
                                      config, err, cmessage)
          if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

        ! ----- no other sections implemented -----
        case default
          message=trim(message)//'section ['//trim(sections(i)%key)//'] not implemented: remove section from TOML file'
          err=10; return

      end select   ! (cases within a desired section)

    end do  ! (looping through sub-sections)
  end do  ! (looping through sections)

  ! ----- check mizuRoute execution constraints -----

  ! Coupled mizuRoute requires the complete set of SUMMA GRUs because runoff from upstream
  ! GRUs may contribute to river reaches outside the selected SUMMA subdomain.
  if (mizuroute_active .and. iRunMode /= iRunModeFull) then
    message=trim(message)//'The -g subdomain option cannot be used with coupled mizuRoute because '// &
                           'the selected GRUs may not contain the complete upstream river network.'
    err=20; return
  endif

  ! ----- set default objective function settings -----

  ! set default objective-function metric
  if(.not.allocated(config%obj%metric))then
    config%obj%metric = 'kge'
    write(iulog,*) 'WARNING: objective metric not specified; using kge'
  endif
  
  ! set default objective-function transformation
  if(.not.allocated(config%obj%transformation))then
    config%obj%transformation = 'none'
    write(iulog,*) 'WARNING: objective transformation not specified; using none'
  endif

  ! check start_date and end_date are defined
  if(.not.allocated(config%obj%start_date) .or. .not.allocated(config%obj%end_date) )then
    message=trim(message)//'Objective function start_date or end_date are not defined'
    err=20; return
  endif

  end subroutine load_summa_config

 ! --------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------

 ! **************************************************************************************************
 ! Parse observation and objective-function configuration.
 ! **************************************************************************************************

  subroutine parse_objective_config(subtable, section, key, config, ierr, message)

  use tomlf_all, only: toml_table, toml_array, toml_error, toml_key, toml_value ! data types
  use tomlf_all, only: toml_load, get_value, len                                ! procedures

  type(toml_table), pointer, intent(in)    :: subtable
  character(*),              intent(in)    :: section
  character(*),              intent(in)    :: key
  type(config_info),         intent(inout) :: config
  integer,                   intent(out)   :: ierr
  character(*),              intent(out)   :: message

  integer(i4b)       :: istat

  associate(obs => config%obs, &
            obj => config%obj)

  ierr    = 0
  message = 'parse_obs_config/'

  ! extract configuration values and populate structures

  select case(trim(section)//'.'//trim(key))

    ! ---- observations: filename ----
    case ("observations.obs_path"        ); call get_value(subtable, trim(key), obs%obs_path           , stat=istat)
    case ("observations.obs_file"        ); call get_value(subtable, trim(key), obs%obs_file           , stat=istat)

    ! ---- observations: variable names ----
    case ("observations.vname_obsflow"   ); call get_value(subtable, trim(key), obs%vname_obsflow      , stat=istat)

    ! ---- objective function: metrics  ----
    case ("objective.metric"             ); call get_value(subtable, trim(key), obj%metric             , stat=istat)
    case ("objective.transformation"     ); call get_value(subtable, trim(key), obj%transformation     , stat=istat)

    ! ---- objective function: calibration period  ----
    case ("objective.start_date"         ); call get_value(subtable, trim(key), obj%start_date         , stat=istat)
    case ("objective.end_date"           ); call get_value(subtable, trim(key), obj%end_date           , stat=istat)

    ! ---- objective function: flag to write aligned sim/obs time series  ----
    case ("objective.write_aligned"      ); call get_value(subtable, trim(key), obj%write_aligned      , stat=istat)
    
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

  end subroutine parse_objective_config





end module summa_config

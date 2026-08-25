module summa_config

use build_options, only: mizuroute_active

USE nr_type
USE summa_type, only:summa1_type_dec

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

  subroutine load_summa_config(config_file, summaStruct, err, message)


  use tomlf_all, only: toml_table, toml_array, toml_error, toml_key, toml_value ! data types
  use tomlf_all, only: toml_load, get_value, len                                ! procedures

  implicit none

  character(*),          intent(in)    :: config_file
  type(summa1_type_dec), intent(inout) :: summaStruct
  integer,               intent(out)   :: err
  character(*),          intent(out)   :: message

  ! TOML table
  type(toml_table), allocatable :: table       ! root TOML table
  type(toml_table), pointer     :: subtable    ! sub-table for a given section
  type(toml_key),   allocatable :: sections(:) ! top-level sections
  type(toml_key),   allocatable :: keys(:)     ! sub-table keys
  type(toml_error), allocatable :: error

  ! locals
  integer(i4b)       :: i, j
  character(len=256) :: cmessage

  err = 0
  message = 'load_summa_config/'
  print*, trim(message), mizuroute_active

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
    message = trim(message)//"problem loading toml sections['"//trim(summaStruct%summaConfigFile)//"']"
    err=10; return
  endif

  ! ----- loop through sections -----
  do i = 1, size(sections)

    ! ----- load the TOML sub-table for the current section -----
    call get_value(table, trim(sections(i)%key), subtable, requested=.false.)
    if(.not.associated(subtable)) then
      message = trim(message)//"problem loading toml sub-sections['"//trim(summaStruct%summaConfigFile)//"']:"//trim(sections(i)%key)
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
                                        summaStruct, err, cmessage)
            if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
          endif

        ! ----- no other sections implemented -----
        case default
          message=trim(message)//'section ['//trim(sections(i)%key)//'] not implemented: remove section from TOML file'
          err=10; return

      end select   ! (cases within a desired section)

    end do  ! (looping through sub-sections)
  end do  ! (looping through sections)

  end subroutine load_summa_config

end module summa_config

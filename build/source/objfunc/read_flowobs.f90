module read_flowobs_module

  USE netcdf
  USE nr_type, only: i4b, i8b, rkind
  USE summa_type, only: summa1_type_dec

  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

  implicit none
  private

  public :: read_flow_observations

contains

  ! **************************************************************************************************
  ! Read observed streamflow and time coordinates from a NetCDF file
  ! **************************************************************************************************
  subroutine read_flow_observations(summaStruc,            &
                                    timeObs, flowObs,      &
                                    timeUnits, flowUnits,  &
                                    err, message)

    type(summa1_type_dec), intent(in)           :: summaStruc    ! master summa data structure

    real(rkind), allocatable, intent(out)       :: timeObs(:)    ! observation time coordinate
    real(rkind), allocatable, intent(out)       :: flowObs(:)    ! observed streamflow

    character(len=:), allocatable, intent(out)  :: timeUnits     ! units and reference time
    character(len=:), allocatable, intent(out)  :: flowUnits     ! streamflow units

    integer(i4b), intent(out)                   :: err           ! error code
    character(*), intent(out)                   :: message       ! error message

    integer(i4b) :: ncid
    integer(i4b) :: dimid
    integer(i4b) :: varid_time
    integer(i4b) :: varid_flow
    integer(i4b) :: nTime
    integer(i4b) :: attLen
    integer(i4b) :: err_close

    integer(i8b), allocatable :: timeInt(:)

    character(len=:), allocatable :: units
    character(len=:), allocatable :: vname_obsflow

    logical :: file_open

    character(len=256) :: cmessage

    err = 0
    message = 'read_flow_observations/'

    associate(obs => summaStruc%config%obs)

    ! check that the observation file information is defined
    if(.not.allocated(obs%obs_path) .or. &
       .not.allocated(obs%obs_file))then
       message=trim(message)//'observation file path or filename is not defined'
       err=20; return
    endif

    ! check that the variable name is defined
    if(allocated(obs%vname_obsflow))then
      vname_obsflow = trim(obs%vname_obsflow)
    else
      vname_obsflow = 'q_obs'
    endif

    file_open = .false.

    netcdf_block: block

      ! open observation file
      err = nf90_open(trim(obs%obs_path)// &
                      trim(obs%obs_file), NF90_NOWRITE, ncid)
      if(err/=nf90_noerr) exit netcdf_block
      file_open = .true.

      ! get time dimension
      err = nf90_inq_dimid(ncid, 'time', dimid)
      if(err/=nf90_noerr) exit netcdf_block

      err = nf90_inquire_dimension(ncid, dimid, len=nTime)
      if(err/=nf90_noerr) exit netcdf_block

      ! get variable IDs
      err = nf90_inq_varid(ncid, 'time', varid_time)
      if(err/=nf90_noerr) exit netcdf_block

      err = nf90_inq_varid(ncid, trim(vname_obsflow), varid_flow)
      if(err/=nf90_noerr) exit netcdf_block

      ! allocate time series
      allocate(timeInt(nTime), timeObs(nTime), flowObs(nTime), stat=err)
      if(err/=0)then; message=trim(message)//'problem allocating'; return; endif

      ! read time
      err = nf90_get_var(ncid, varid_time, timeInt)
      if(err/=nf90_noerr) exit netcdf_block

      timeObs = real(timeInt, rkind)

      ! read streamflow
      err = nf90_get_var(ncid, varid_flow, flowObs)
      if(err/=nf90_noerr) exit netcdf_block

      ! normalize missing values
      call normalize_missing_values(ncid,varid_flow,flowObs,err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! read time units
      err = nf90_inquire_attribute(ncid, varid_time, 'units', len=attLen)
      if(err/=nf90_noerr) exit netcdf_block

      allocate(character(len=attLen) :: units)

      err = nf90_get_att(ncid, varid_time, 'units', units)
      if(err/=nf90_noerr) exit netcdf_block

      timeUnits = trim(units)
      deallocate(units)

      ! read flow units
      err = nf90_inquire_attribute(ncid, varid_flow, 'units', len=attLen)
      if(err/=nf90_noerr) exit netcdf_block

      allocate(character(len=attLen) :: units)

      err = nf90_get_att(ncid, varid_flow, 'units', units)
      if(err/=nf90_noerr) exit netcdf_block

      flowUnits = trim(units)
      deallocate(units)

      ! close observation file
      err = nf90_close(ncid)
      if(err/=nf90_noerr) exit netcdf_block
      file_open = .false.

    end block netcdf_block

    ! process NetCDF errors
    if(err/=nf90_noerr)then
      message = trim(message)//trim(nf90_strerror(err))
      if(file_open) err_close = nf90_close(ncid)
      return
    endif

    err = 0

    end associate

  end subroutine read_flow_observations

  ! **************************************************************************************************
  ! Normalize missing values in a NetCDF variable.
  !
  ! Missing values are converted to IEEE quiet NaN so that downstream routines
  ! can handle missing data independently of the conventions used in the input
  ! file. The routine first checks for an explicit _FillValue attribute, then
  ! for a missing_value attribute, and otherwise uses the default NetCDF fill
  ! value.
  !
  ! Input values that are already NaN are left unchanged.
  ! **************************************************************************************************

  subroutine normalize_missing_values(ncid, varid, values, err, message)

    USE netcdf
    USE, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    USE, intrinsic :: ieee_arithmetic, only: ieee_value
    USE, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan
   
    integer(i4b), intent(in)    :: ncid
    integer(i4b), intent(in)    :: varid
    real(rkind), intent(inout)  :: values(:)
   
    integer(i4b), intent(out)   :: err
    character(*), intent(out)   :: message
   
    real(rkind) :: fillValue
    real(rkind) :: nanValue
   
    integer(i4b) :: ierr
   
    logical :: hasFillValue
   
    err = 0
    message = 'normalize_missing_values/'
   
    hasFillValue = .false.
    nanValue = ieee_value(0._rkind,ieee_quiet_nan)
   
    ! first try _FillValue
    ierr = nf90_get_att(ncid,varid,'_FillValue',fillValue)
   
    if(ierr==nf90_noerr)then
   
      hasFillValue = .true.
   
    else if(ierr==nf90_enotatt)then
   
      ! try missing_value
      ierr = nf90_get_att(ncid,varid,'missing_value',fillValue)
   
      if(ierr==nf90_noerr)then
   
        hasFillValue = .true.
   
      else if(ierr==nf90_enotatt)then
   
        ! no explicit missing-value attribute
        ierr = nf90_noerr
   
      else
   
        message=trim(message)//trim(nf90_strerror(ierr))
        err=ierr
        return
   
      endif
   
    else
   
      message=trim(message)//trim(nf90_strerror(ierr))
      err=ierr
      return
   
    endif
   
    ! replace explicitly defined missing values with NaN
    if(hasFillValue)then
   
      if(.not.ieee_is_nan(fillValue))then
        where(values==fillValue)
          values=nanValue
        endwhere
      endif
   
    else
   
      ! replace NetCDF default double fill value
      where(values==NF90_FILL_DOUBLE)
        values=nanValue
      endwhere
   
    endif
   
  end subroutine normalize_missing_values



end module read_flowobs_module

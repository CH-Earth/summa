module read_flowobs_module

  USE netcdf
  USE nr_type, only: i4b, i8b, rkind, lgt
  USE summa_type, only: summa1_type_dec

  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

  implicit none
  private

  public :: read_flow_observations
  public :: read_observations

contains

  ! **************************************************************************************************
  ! Read observed streamflow and time coordinates from a NetCDF file
  !
  ! Thin wrapper over read_observations for the single streamflow series named in the [observations]
  ! section, which is what a single-objective calibration compares against.
  ! **************************************************************************************************
  subroutine read_flow_observations(summaStruc, timeObs, timeBounds, flowObs, timeUnits, flowUnits, err, message)
    type(summa1_type_dec), intent(in)           :: summaStruc       ! master summa data structure
    real(rkind), allocatable, intent(out)       :: timeObs(:)       ! observation time coordinate
    real(rkind), allocatable, intent(out)       :: timeBounds(:,:)  ! observation time bounds (2,nTime)
    real(rkind), allocatable, intent(out)       :: flowObs(:)       ! observed streamflow
    character(len=:), allocatable, intent(out)  :: timeUnits        ! units and reference time
    character(len=:), allocatable, intent(out)  :: flowUnits        ! streamflow units
    integer(i4b), intent(out)                   :: err              ! error code
    character(*), intent(out)                   :: message          ! error message
    character(len=:), allocatable :: vname_obsflow
    character(len=256) :: cmessage

    err = 0
    message = 'read_flow_observations/'
    associate(obs => summaStruc%config%obs, calib => summaStruc%config%calib)

    ! check calibration period is defined
    if(.not.allocated(calib%start_date) .or. .not.allocated(calib%end_date))then
      message=trim(message)//'calibration start_date or end_date are not defined'
      err=20; return
    endif

    ! check that the observation file information is defined
    if(.not.allocated(obs%obs_path) .or. .not.allocated(obs%obs_file))then
       message=trim(message)//'observation file path or filename is not defined'
       err=20; return
    endif

    ! check that the variable name is defined
    if(allocated(obs%vname_obsflow))then
      vname_obsflow = trim(obs%vname_obsflow)
    else
      vname_obsflow = 'q_obs'
    endif

    call read_observations(trim(obs%obs_path), trim(obs%obs_file), vname_obsflow, &
                           timeObs, timeBounds, flowObs, timeUnits, flowUnits, err, cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    end associate

  end subroutine read_flow_observations

  ! **************************************************************************************************
  ! Read one observed series and its time coordinate from a NetCDF file
  !
  ! The series is named explicitly rather than taken from the configuration, so a calibration with
  ! several targets can read a different variable, from a different file, for each of them.
  ! **************************************************************************************************
  subroutine read_observations(obs_path, obs_file, vname_obs, timeObs, timeBounds, valObs, timeUnits, valUnits, err, message)
    character(*), intent(in)                    :: obs_path      ! path to the observation file
    character(*), intent(in)                    :: obs_file      ! observation file
    character(*), intent(in)                    :: vname_obs     ! name of the observed variable
    real(rkind), allocatable, intent(out)       :: timeObs(:)    ! observation time coordinate
    real(rkind), allocatable, intent(out)       :: timeBounds(:,:) ! observation time bounds (2,nTime)
    real(rkind), allocatable, intent(out)       :: valObs(:)     ! observed values
    character(len=:), allocatable, intent(out)  :: timeUnits     ! units and reference time
    character(len=:), allocatable, intent(out)  :: valUnits      ! units of the observed variable
    integer(i4b), intent(out)                   :: err           ! error code
    character(*), intent(out)                   :: message       ! error message
    integer(i4b) :: ncid
    integer(i4b) :: dimid
    integer(i4b) :: varid_time
    integer(i4b) :: varid_bounds
    integer(i4b) :: varid_flow
    integer(i4b) :: nTime
    integer(i4b) :: attLen
    integer(i4b) :: err_close
    logical(lgt) :: file_exists
    logical(lgt) :: has_bounds
    integer(i8b), allocatable :: timeInt(:)
    integer(i8b), allocatable :: boundsInt(:,:)
    integer(i8b), allocatable :: boundsByTime(:,:)
    integer(i4b) :: boundDims(2)
    character(len=:), allocatable :: units
    logical :: file_open
    character(len=256) :: cmessage

    err = 0
    message = 'read_observations/'

    ! check that the observation file exists
    inquire(file=trim(obs_path)//trim(obs_file),exist=file_exists)
    if(.not.file_exists)then
      message=trim(message)//'observation file does not exist: '// trim(obs_path)//trim(obs_file)
      err=20; return
    endif

    file_open = .false.
    netcdf_block: block

      ! open observation file
      err = nf90_open(trim(obs_path)// trim(obs_file), NF90_NOWRITE, ncid)
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
      err = nf90_inq_varid(ncid, trim(vname_obs), varid_flow)
      if(err/=nf90_noerr) exit netcdf_block

      ! time bounds are optional: without them each observation covers the span since the one before it
      has_bounds = (nf90_inq_varid(ncid, 'time_bnds', varid_bounds) == nf90_noerr)
      if(err/=nf90_noerr) exit netcdf_block

      ! allocate time series
      allocate(timeInt(nTime), timeObs(nTime), timeBounds(2,nTime), valObs(nTime), stat=err)
      if(err/=0)then; message=trim(message)//'problem allocating'; return; endif

      ! read time
      err = nf90_get_var(ncid, varid_time, timeInt)
      if(err/=nf90_noerr) exit netcdf_block
      timeObs = real(timeInt, rkind)

      ! read time bounds, stored as time_bnds(time,nbnds) or time_bnds(nbnds,time), or take each period as ending at its time stamp
      if(has_bounds)then
        err = nf90_inquire_variable(ncid, varid_bounds, dimids=boundDims)
        if(err/=nf90_noerr) exit netcdf_block
        if(boundDims(1) == dimid)then ! time varies fastest: time_bnds(nbnds,time) in the file
          allocate(boundsByTime(nTime,2), stat=err)
          if(err/=0)then; message=trim(message)//'problem allocating'; return; endif
          err = nf90_get_var(ncid, varid_bounds, boundsByTime)
          if(err/=nf90_noerr) exit netcdf_block
          boundsInt = transpose(boundsByTime)
        else
          allocate(boundsInt(2,nTime), stat=err)
          if(err/=0)then; message=trim(message)//'problem allocating'; return; endif
          err = nf90_get_var(ncid, varid_bounds, boundsInt)
          if(err/=nf90_noerr) exit netcdf_block
        endif
        timeBounds = real(boundsInt, rkind)
      else
        timeBounds(2,:) = timeObs
        if(nTime > 1)then
          timeBounds(1,2:nTime) = timeObs(1:nTime-1)
          timeBounds(1,1) = timeObs(1) - (timeObs(2) - timeObs(1))
        else
          timeBounds(1,1) = timeObs(1)
        endif
      endif

      ! read the observed series
      err = nf90_get_var(ncid, varid_flow, valObs)
      if(err/=nf90_noerr) exit netcdf_block

      ! normalize missing values
      call normalize_missing_values(ncid,varid_flow,valObs,err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

      ! read time units
      err = nf90_inquire_attribute(ncid, varid_time, 'units', len=attLen)
      if(err/=nf90_noerr) exit netcdf_block
      allocate(character(len=attLen) :: units)
      err = nf90_get_att(ncid, varid_time, 'units', units)
      if(err/=nf90_noerr) exit netcdf_block
      timeUnits = trim(units)
      deallocate(units)

      ! read the units of the observed series
      err = nf90_inquire_attribute(ncid, varid_flow, 'units', len=attLen)
      if(err/=nf90_noerr) exit netcdf_block
      allocate(character(len=attLen) :: units)
      err = nf90_get_att(ncid, varid_flow, 'units', units)
      if(err/=nf90_noerr) exit netcdf_block
      valUnits = trim(units)
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

  end subroutine read_observations

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

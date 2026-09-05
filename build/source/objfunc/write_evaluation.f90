module write_evaluation_module

  use nr_type, only: i4b, rkind
  use netcdf

  implicit none
  private

  public :: write_evaluation

contains


  ! **************************************************************************************************
  ! Write objective-function evaluation results to the SUMMA NetCDF output file.
  !
  ! The objective-function value is always written. The aligned observed and
  ! simulated streamflow time series are written only when write_aligned is true.
  ! **************************************************************************************************
  subroutine write_evaluation(ncid,write_aligned,          &
                              timeEval,flowObs,flowSim,   &
                              timeUnits,flowUnits,        &
                              metric,transformation,      &
                              objective,                  &
                              err,message)

    integer(i4b), intent(in)  :: ncid

    logical, intent(in)       :: write_aligned

    real(rkind), intent(in)   :: timeEval(:)
    real(rkind), intent(in)   :: flowObs(:)
    real(rkind), intent(in)   :: flowSim(:)

    character(*), intent(in)  :: timeUnits
    character(*), intent(in)  :: flowUnits
    character(*), intent(in)  :: metric
    character(*), intent(in)  :: transformation

    real(rkind), intent(in)   :: objective

    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message

    character(len=256) :: cmessage

    err = 0
    message = 'write_evaluation/'

    ! define objective-function variable
    call define_objective(ncid,metric,transformation,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! optionally define aligned evaluation time series
    if(write_aligned)then

      call define_evaluation_series(ncid,size(timeEval), &
                                    timeUnits,flowUnits,  &
                                    err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    endif

    ! write objective-function value
    call write_objective(ncid,objective,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! optionally write aligned evaluation time series
    if(write_aligned)then

      call write_evaluation_series(ncid,                   &
                                   timeEval,flowObs,flowSim,&
                                   err,cmessage)
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    endif

  end subroutine write_evaluation


  ! **************************************************************************************************
  ! Add the objective-function variable and metadata to an existing SUMMA NetCDF file.
  ! **************************************************************************************************
  subroutine define_objective(ncid,metric,transformation,ierr,message)

    integer(i4b), intent(in)  :: ncid

    character(*), intent(in)  :: metric
    character(*), intent(in)  :: transformation

    integer(i4b), intent(out) :: ierr
    character(*), intent(out) :: message

    integer(i4b) :: varid_obj
    integer(i4b) :: ierr_enddef

    logical :: in_define

    ierr = 0
    message = 'define_objective/'

    in_define = .false.

    netcdf_block: block

      ! enter (re)-define mode
      ierr = nf90_redef(ncid)
      if(ierr/=nf90_noerr) exit netcdf_block

      in_define = .true.

      ! objective-function value
      ierr = nf90_def_var(ncid,'objective',NF90_DOUBLE,varid=varid_obj)
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr = nf90_put_att(ncid,varid_obj,'long_name', &
                          'objective-function value')
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr = nf90_put_att(ncid,varid_obj,'metric',trim(metric))
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr = nf90_put_att(ncid,varid_obj,'transformation',trim(transformation))
      if(ierr/=nf90_noerr) exit netcdf_block

      ! leave define mode
      ierr = nf90_enddef(ncid)
      if(ierr/=nf90_noerr) exit netcdf_block

      in_define = .false.

    end block netcdf_block

    ! process NetCDF errors
    if(ierr/=nf90_noerr)then
      message = trim(message)//trim(nf90_strerror(ierr))
      if(in_define) ierr_enddef = nf90_enddef(ncid)
      return
    endif

    ierr = 0

  end subroutine define_objective


  ! **************************************************************************************************
  ! Add aligned objective-function evaluation time series to an existing SUMMA NetCDF file.
  !
  ! The evaluation time series has its own time dimension because the observation
  ! timestep may differ from the native SUMMA simulation timestep.
  ! **************************************************************************************************
  subroutine define_evaluation_series(ncid,nEval,timeUnits,flowUnits, &
                                      ierr,message)

    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

    integer(i4b), intent(in)  :: ncid
    integer(i4b), intent(in)  :: nEval

    character(*), intent(in)  :: timeUnits
    character(*), intent(in)  :: flowUnits

    integer(i4b), intent(out) :: ierr
    character(*), intent(out) :: message

    integer(i4b) :: dim_eval

    integer(i4b) :: varid_time
    integer(i4b) :: varid_obs
    integer(i4b) :: varid_sim

    integer(i4b) :: ierr_enddef

    real(rkind) :: nanValue

    logical :: in_define

    ierr = 0
    message = 'define_evaluation_series/'

    in_define = .false.

    nanValue = ieee_value(0._rkind,ieee_quiet_nan)

    netcdf_block: block

      ! enter (re)-define mode
      ierr = nf90_redef(ncid)
      if(ierr/=nf90_noerr) exit netcdf_block

      in_define = .true.

      ! evaluation time dimension
      ierr = nf90_def_dim(ncid,'eval_time',nEval,dim_eval)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! evaluation time coordinate
      ierr = nf90_def_var(ncid,'eval_time',NF90_DOUBLE,(/dim_eval/),varid_time)
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr = nf90_put_att(ncid,varid_time,'long_name', &
                          'objective-function evaluation time')
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr = nf90_put_att(ncid,varid_time,'standard_name','time')
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr = nf90_put_att(ncid,varid_time,'units',trim(timeUnits))
      if(ierr/=nf90_noerr) exit netcdf_block

      ! observed streamflow
      ierr = nf90_def_var(ncid,'eval_qobs',NF90_DOUBLE,(/dim_eval/),varid_obs)
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr = nf90_put_att(ncid,varid_obs,'long_name', &
                          'observed streamflow used for objective-function evaluation')
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr = nf90_put_att(ncid,varid_obs,'units',trim(flowUnits))
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr = nf90_put_att(ncid,varid_obs,'_FillValue',nanValue)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! simulated streamflow aligned to observation periods
      ierr = nf90_def_var(ncid,'eval_qsim',NF90_DOUBLE,(/dim_eval/),varid_sim)
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr = nf90_put_att(ncid,varid_sim,'long_name', &
                          'simulated streamflow aligned to observation periods')
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr = nf90_put_att(ncid,varid_sim,'units',trim(flowUnits))
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr = nf90_put_att(ncid,varid_sim,'_FillValue',nanValue)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! leave define mode
      ierr = nf90_enddef(ncid)
      if(ierr/=nf90_noerr) exit netcdf_block

      in_define = .false.

    end block netcdf_block

    ! process NetCDF errors
    if(ierr/=nf90_noerr)then
      message = trim(message)//trim(nf90_strerror(ierr))
      if(in_define) ierr_enddef = nf90_enddef(ncid)
      return
    endif

    ierr = 0

  end subroutine define_evaluation_series


  ! **************************************************************************************************
  ! Write the objective-function value to an existing SUMMA NetCDF file.
  ! **************************************************************************************************
  subroutine write_objective(ncid,objective,ierr,message)

    integer(i4b), intent(in)  :: ncid
    real(rkind), intent(in)   :: objective

    integer(i4b), intent(out) :: ierr
    character(*), intent(out) :: message

    integer(i4b) :: varid_obj

    ierr = 0
    message = 'write_objective/'

    netcdf_block: block

      ierr = nf90_inq_varid(ncid,'objective',varid_obj)
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr = nf90_put_var(ncid,varid_obj,objective)
      if(ierr/=nf90_noerr) exit netcdf_block

    end block netcdf_block

    ! process NetCDF errors
    if(ierr/=nf90_noerr)then
      message = trim(message)//trim(nf90_strerror(ierr))
      return
    endif

    ierr = 0

  end subroutine write_objective


  ! **************************************************************************************************
  ! Write aligned observed and simulated streamflow time series to an existing SUMMA NetCDF file.
  ! **************************************************************************************************
  subroutine write_evaluation_series(ncid,timeEval,flowObs,flowSim, &
                                     ierr,message)

    integer(i4b), intent(in)  :: ncid

    real(rkind), intent(in)   :: timeEval(:)
    real(rkind), intent(in)   :: flowObs(:)
    real(rkind), intent(in)   :: flowSim(:)

    integer(i4b), intent(out) :: ierr
    character(*), intent(out) :: message

    integer(i4b) :: varid_time
    integer(i4b) :: varid_obs
    integer(i4b) :: varid_sim

    ierr = 0
    message = 'write_evaluation_series/'

    ! check dimensions
    if(size(flowObs)/=size(timeEval) .or. &
       size(flowSim)/=size(timeEval))then
      message=trim(message)//'evaluation time-series dimensions differ'
      ierr=20; return
    endif

    netcdf_block: block

      ! get variable IDs
      ierr = nf90_inq_varid(ncid,'eval_time',varid_time)
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr = nf90_inq_varid(ncid,'eval_qobs',varid_obs)
      if(ierr/=nf90_noerr) exit netcdf_block

      ierr = nf90_inq_varid(ncid,'eval_qsim',varid_sim)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! evaluation time coordinate
      ierr = nf90_put_var(ncid,varid_time,timeEval)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! observed streamflow
      ierr = nf90_put_var(ncid,varid_obs,flowObs)
      if(ierr/=nf90_noerr) exit netcdf_block

      ! simulated streamflow
      ierr = nf90_put_var(ncid,varid_sim,flowSim)
      if(ierr/=nf90_noerr) exit netcdf_block

    end block netcdf_block

    ! process NetCDF errors
    if(ierr/=nf90_noerr)then
      message = trim(message)//trim(nf90_strerror(ierr))
      return
    endif

    ierr = 0

  end subroutine write_evaluation_series


end module write_evaluation_module

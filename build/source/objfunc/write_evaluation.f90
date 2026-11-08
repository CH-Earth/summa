module write_evaluation_module

  use nr_type, only: i4b, rkind
  use netcdf

  implicit none
  private

  public :: write_evaluation

contains

  ! **************************************************************************************************
  ! Write one calibration target to its own NetCDF4 group, under /evaluation/<targetName>.
  !
  ! Each target has an independent time dimension. This is essential when, for example, daily
  ! streamflow is scored alongside irregularly spaced GRACE storage observations.
  ! The exact observation time bounds are retained, rather than inferred from adjacent timestamps.
  ! Optional inputs precede err/message; callers use keyword error arguments.
  ! **************************************************************************************************
  subroutine write_evaluation(ncid,targetName,write_aligned,                    &
                              timeEval,timeBounds,valObs,valSim,nSimSteps,     &
                              timeUnits,valUnits,metric,transformation,       &
                              objective,boxcox_lambda,err,message)
    integer(i4b), intent(in)  :: ncid
    character(*), intent(in)  :: targetName
    logical, intent(in)       :: write_aligned
    real(rkind), intent(in)   :: timeEval(:)
    real(rkind), intent(in)   :: timeBounds(:,:) ! (2,nEval), in timeUnits
    real(rkind), intent(in)   :: valObs(:)
    real(rkind), intent(in)   :: valSim(:)
    integer(i4b), intent(in) :: nSimSteps(:)
    character(*), intent(in)  :: timeUnits
    character(*), intent(in)  :: valUnits
    character(*), intent(in)  :: metric
    character(*), intent(in)  :: transformation
    real(rkind), intent(in)   :: objective
    real(rkind), optional, intent(in) :: boxcox_lambda
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    integer(i4b) :: groupTarget
    character(len=512) :: cmessage

    err = 0
    message = 'write_evaluation/'

    ! All diagnostics must refer to the same observation intervals.
    if(write_aligned)then
      if(size(valObs)/=size(timeEval) .or. size(valSim)/=size(timeEval) .or. &
         size(nSimSteps)/=size(timeEval))then
        message=trim(message)//'evaluation time-series dimensions differ'
        err=20; return
      endif
      if(size(timeBounds,1)/=2 .or. size(timeBounds,2)/=size(timeEval))then
        message=trim(message)//'evaluation time bounds must have shape (2,nEval)'
        err=20; return
      endif
    endif
    if(len_trim(targetName)==0)then
      message=trim(message)//'calibration target name is empty'
      err=20; return
    endif

    call define_evaluation_target(ncid,targetName,write_aligned,size(timeEval),timeUnits,valUnits,metric, &
                                  transformation,groupTarget,boxcox_lambda=boxcox_lambda,ierr=err, &
                                  message=cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    call write_evaluation_values(groupTarget,write_aligned,timeEval,timeBounds, &
                                 valObs,valSim,nSimSteps,objective,err,cmessage)
    if(err/=0) message=trim(message)//trim(cmessage)

  end subroutine write_evaluation

  ! **************************************************************************************************
  ! Define one target group and all its variables in the existing SUMMA NetCDF4 output file.
  ! Repeated target names are an error: each target must own a unique group.
  ! **************************************************************************************************
  subroutine define_evaluation_target(ncid,targetName,write_aligned,nEval,            &
                                      timeUnits,valUnits,metric,transformation,      &
                                      groupTarget,boxcox_lambda,ierr,message)
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
    integer(i4b), intent(in)  :: ncid,nEval
    character(*), intent(in)  :: targetName,timeUnits,valUnits,metric,transformation
    logical, intent(in)       :: write_aligned
    integer(i4b), intent(out) :: groupTarget
    real(rkind), optional, intent(in) :: boxcox_lambda
    integer(i4b), intent(out) :: ierr
    character(*), intent(out) :: message
    integer(i4b) :: groupEval
    integer(i4b) :: dimTime,dimNv
    integer(i4b) :: varObjective,varTime,varBounds,varObs,varSim,varSteps
    integer(i4b) :: ierr_enddef
    real(rkind) :: nanValue,lambda
    logical :: in_define

    ierr=0
    message='define_evaluation_target/'
    in_define=.false.
    nanValue=ieee_value(0._rkind,ieee_quiet_nan)
    lambda=0.25_rkind
    if(present(boxcox_lambda)) lambda=boxcox_lambda
    if(trim(transformation)=='boxcox')then
      if(.not.ieee_is_finite(lambda))then
        message=trim(message)//'Box-Cox exponent must be finite'
        ierr=20; return
      endif
    endif

    netcdf_block: block
      ! SUMMA uses NetCDF4, so groups and independent time dimensions are supported.
      ierr=nf90_redef(ncid)
      if(ierr/=nf90_noerr) exit netcdf_block
      in_define=.true.

      ierr=nf90_inq_ncid(ncid,'evaluation',groupEval)
      if(ierr/=nf90_noerr)then
        ierr=nf90_def_grp(ncid,'evaluation',groupEval)
        if(ierr/=nf90_noerr) exit netcdf_block
      endif

      ierr=nf90_def_grp(groupEval,trim(targetName),groupTarget)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(groupTarget,NF90_GLOBAL,'target_name',trim(targetName))
      if(ierr/=nf90_noerr) exit netcdf_block

      ! A scalar objective is always present, even if diagnostics are disabled.
      ierr=nf90_def_var(groupTarget,'objective',NF90_DOUBLE,varid=varObjective)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(groupTarget,varObjective,'long_name','objective-function value')
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(groupTarget,varObjective,'metric',trim(metric))
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_att(groupTarget,varObjective,'transformation',trim(transformation))
      if(ierr/=nf90_noerr) exit netcdf_block
      if(trim(transformation)=='boxcox')then
        ierr=nf90_put_att(groupTarget,varObjective,'boxcox_lambda',lambda)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(groupTarget,varObjective,'boxcox_definition', &
                          '(Q**lambda-1)/lambda; lambda=0: log(Q); no offset')
        if(ierr/=nf90_noerr) exit netcdf_block
      endif

      if(write_aligned)then
        ! Target-local dimensions: different targets may have different numbers of observations.
        ierr=nf90_def_dim(groupTarget,'time',nEval,dimTime)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_def_dim(groupTarget,'nv',2,dimNv)
        if(ierr/=nf90_noerr) exit netcdf_block

        ierr=nf90_def_var(groupTarget,'time',NF90_DOUBLE,(/dimTime/),varTime)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(groupTarget,varTime,'standard_name','time')
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(groupTarget,varTime,'units',trim(timeUnits))
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(groupTarget,varTime,'bounds','time_bnds')
        if(ierr/=nf90_noerr) exit netcdf_block

        ! Fortran dimensions (nv,time) appear as (time,nv) in ncdump / CF files.
        ierr=nf90_def_var(groupTarget,'time_bnds',NF90_DOUBLE,(/dimNv,dimTime/),varBounds)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(groupTarget,varBounds,'units',trim(timeUnits))
        if(ierr/=nf90_noerr) exit netcdf_block

        ierr=nf90_def_var(groupTarget,'obs',NF90_DOUBLE,(/dimTime/),varObs)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(groupTarget,varObs,'long_name','observed values used for evaluation')
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(groupTarget,varObs,'units',trim(valUnits))
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(groupTarget,varObs,'_FillValue',nanValue)
        if(ierr/=nf90_noerr) exit netcdf_block

        ierr=nf90_def_var(groupTarget,'sim',NF90_DOUBLE,(/dimTime/),varSim)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(groupTarget,varSim,'long_name','simulated values aligned to observation periods')
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(groupTarget,varSim,'units',trim(valUnits))
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(groupTarget,varSim,'_FillValue',nanValue)
        if(ierr/=nf90_noerr) exit netcdf_block

        ierr=nf90_def_var(groupTarget,'n_sim_steps',NF90_INT,(/dimTime/),varSteps)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_att(groupTarget,varSteps,'long_name', &
                          'number of simulated time steps in observation period')
        if(ierr/=nf90_noerr) exit netcdf_block
      endif

      ierr=nf90_enddef(ncid)
      if(ierr/=nf90_noerr) exit netcdf_block
      in_define=.false.
    end block netcdf_block

    if(ierr/=nf90_noerr)then
      message=trim(message)//'target "'//trim(targetName)//'": '//trim(nf90_strerror(ierr))
      if(in_define) ierr_enddef=nf90_enddef(ncid)
      return
    endif
    ierr=0

  end subroutine define_evaluation_target

  ! **************************************************************************************************
  ! Write scalar objective and optional target-specific aligned diagnostics.
  ! **************************************************************************************************
  subroutine write_evaluation_values(groupTarget,write_aligned,       &
                                     timeEval,timeBounds,valObs,valSim, &
                                     nSimSteps,objective,ierr,message)
    integer(i4b), intent(in)  :: groupTarget
    logical, intent(in)       :: write_aligned
    real(rkind), intent(in)   :: timeEval(:),timeBounds(:,:),valObs(:),valSim(:)
    integer(i4b), intent(in) :: nSimSteps(:)
    real(rkind), intent(in)   :: objective
    integer(i4b), intent(out) :: ierr
    character(*), intent(out) :: message
    integer(i4b) :: varObjective,varTime,varBounds,varObs,varSim,varSteps

    ierr=0
    message='write_evaluation_values/'

    netcdf_block: block
      ierr=nf90_inq_varid(groupTarget,'objective',varObjective)
      if(ierr/=nf90_noerr) exit netcdf_block
      ierr=nf90_put_var(groupTarget,varObjective,objective)
      if(ierr/=nf90_noerr) exit netcdf_block

      if(write_aligned)then
        ierr=nf90_inq_varid(groupTarget,'time',varTime)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_inq_varid(groupTarget,'time_bnds',varBounds)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_inq_varid(groupTarget,'obs',varObs)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_inq_varid(groupTarget,'sim',varSim)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_inq_varid(groupTarget,'n_sim_steps',varSteps)
        if(ierr/=nf90_noerr) exit netcdf_block

        ierr=nf90_put_var(groupTarget,varTime,timeEval)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_var(groupTarget,varBounds,timeBounds)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_var(groupTarget,varObs,valObs)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_var(groupTarget,varSim,valSim)
        if(ierr/=nf90_noerr) exit netcdf_block
        ierr=nf90_put_var(groupTarget,varSteps,nSimSteps)
        if(ierr/=nf90_noerr) exit netcdf_block
      endif
    end block netcdf_block

    if(ierr/=nf90_noerr)then
      message=trim(message)//trim(nf90_strerror(ierr))
      return
    endif
    ierr=0

  end subroutine write_evaluation_values

end module write_evaluation_module

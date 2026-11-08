module metrics

  use nr_type, only: i4b, rkind
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan

  implicit none
  private

  public :: compute_metric
  public :: metric_is_maximized

  ! A mathematically undefined score is different from a configuration error (20).
  ! Callers may retain an undefined diagnostic without failing the simulation.
  integer(i4b), parameter, public :: metric_undefined = 21

contains

  ! **************************************************************************************************
  ! Report whether a larger value of a metric is the better fit.
  !
  ! Efficiencies (KGE, KGE', NSE) improve as they increase; error metrics (MAE, RMSE) improve as they
  ! decrease.  A search combining several targets has to orient them the same way before it can compare
  ! or combine them, and it is the metric, not the target, that decides which way each one points.
  !
  ! An unrecognized metric is reported as maximized: compute_metric rejects it, so the caller never
  ! gets far enough to use the orientation.
  ! **************************************************************************************************
  pure function metric_is_maximized(metric) result(maximize)
    character(*), intent(in) :: metric       ! objective-function metric
    logical                  :: maximize     ! .true. if larger is better

    select case(trim(metric))
      case ('mae','rmse'); maximize = .false.
      case default;        maximize = .true.
    end select

  end function metric_is_maximized

  ! **************************************************************************************************
  ! Compute the selected objective-function metric.
  !
  ! Supported metrics:
  !   kge   : Kling-Gupta efficiency
  !   kgep  : modified Kling-Gupta efficiency
  !   nse   : Nash-Sutcliffe efficiency
  !   mae   : mean absolute error
  !   rmse  : root mean square error
  !
  ! Missing values are removed and the selected streamflow transformation is
  ! applied before evaluating the metric.
  ! **************************************************************************************************
  subroutine compute_metric(obs,sim,metric,transfo,objective,err,message)
    ! dummies
    real(rkind), intent(in)             :: obs(:)       ! observed streamflow
    real(rkind), intent(in)             :: sim(:)       ! simulated streamflow
    character(*), intent(in)            :: metric       ! objective-function metric
    character(*), intent(in)            :: transfo      ! streamflow transformation
    real(rkind), intent(out)            :: objective    ! objective-function value
    integer(i4b), intent(out)           :: err          ! error code
    character(*), intent(out)           :: message      ! error message
    ! locals
    real(rkind), allocatable :: obsUse(:)
    real(rkind), allocatable :: simUse(:)
    character(len=256) :: cmessage

    err = 0
    message = 'compute_metric/'
    objective = ieee_value(0._rkind,ieee_quiet_nan)

    ! Validate the metric even when the series contain no usable observations.
    select case(trim(metric))
      case ('kge','kgep','nse','mae','rmse')
      case default
        message=trim(message)//'unknown objective-function metric "'//trim(metric)//'"'
        err=20; return
    end select

    ! remove missing values and conduct any transformations
    call prepare_series(obs,sim,obsUse,simUse,transfo,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! Check denominators before evaluating efficiencies. Do not turn undefined
    ! scores into finite penalties: those would look like valid diagnostics.
    select case(trim(metric))
      case ('kge','kgep','nse')
        if(size(obsUse)<2)then
          message=trim(message)//'efficiency requires at least two valid pairs'
          err=metric_undefined; return
        endif
        if(all(obsUse==obsUse(1)) .or. sum((obsUse-sum(obsUse)/size(obsUse))**2)<=0._rkind)then
          message=trim(message)//'observed series has zero variance'
          err=metric_undefined; return
        endif
    end select
    if(trim(metric)=='kge' .or. trim(metric)=='kgep')then
      if(all(simUse==simUse(1)) .or. sum((simUse-sum(simUse)/size(simUse))**2)<=0._rkind)then
        message=trim(message)//'simulated series has zero variance for KGE'
        err=metric_undefined; return
      endif
      if(sum(obsUse)/size(obsUse)==0._rkind)then
        message=trim(message)//'observed series has zero mean for KGE'
        err=metric_undefined; return
      endif
      if(trim(metric)=='kgep')then
        if(sum(simUse)/size(simUse)==0._rkind)then
          message=trim(message)//'simulated series has zero mean for modified KGE'
          err=metric_undefined; return
        endif
      endif
    endif

    ! compute selected metric (unchanged formulas for valid inputs)
    select case(trim(metric))
      case ('kge');  objective = get_kge( obsUse,simUse)
      case ('kgep'); objective = get_kgep(obsUse,simUse)
      case ('nse');  objective = get_nse( obsUse,simUse)
      case ('mae');  objective = get_mae( obsUse,simUse)
      case ('rmse'); objective = get_rmse(obsUse,simUse)
      case default
        message=trim(message)//'unknown objective-function metric "'//trim(metric)//'"'
        err=20; return
    end select

    if(.not.ieee_is_finite(objective))then
      objective=ieee_value(0._rkind,ieee_quiet_nan)
      message=trim(message)//'metric calculation produced a non-finite value'
      err=metric_undefined
    endif

  end subroutine compute_metric

  ! **************************************************************************************************
  ! Prepare observed and simulated streamflow for metric calculation.
  !
  ! Retain pairs where both observed and simulated streamflow are finite, then
  ! apply the selected transformation to both series.
  ! **************************************************************************************************
  subroutine prepare_series(obs,sim,obsUse,simUse,transfo,err,message)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    real(rkind), intent(in)               :: obs(:)
    real(rkind), intent(in)               :: sim(:)
    real(rkind), allocatable, intent(out) :: obsUse(:)
    real(rkind), allocatable, intent(out) :: simUse(:)
    character(*), intent(in)              :: transfo
    integer(i4b), intent(out)             :: err
    character(*), intent(out)             :: message
    logical, allocatable :: valid(:)
    integer(i4b) :: n

    err = 0
    message = 'prepare_series/'

    ! check dimensions
    if(size(obs)/=size(sim))then
      message=trim(message)//'observed and simulated streamflow dimensions differ'
      err=20; return
    endif

    ! identify finite observation-simulation pairs
    valid = ieee_is_finite(obs) .and. ieee_is_finite(sim)
    n = count(valid)
    ! retain only valid pairs
    allocate(obsUse(n),simUse(n))
    obsUse = pack(obs,valid)
    simUse = pack(sim,valid)

    ! Validate the transformation even for an empty series. Invalid configuration
    ! must never be hidden by the diagnostic-only handling of undefined metrics.
    call apply_transformation(obsUse,simUse,transfo,err,message)
    if(err/=0) return
    if(n==0)then
      message='prepare_series/no valid observation-simulation pairs'
      err=metric_undefined; return
    endif

  end subroutine prepare_series

  ! **************************************************************************************************
  ! Apply a transformation to observed and simulated streamflow before calculating a metric.
  !
  ! Supported transformations:
  !   log     : logarithmic transformation
  !   boxcox  : Box-Cox transformation with exponent 0.25
  !   numeric : power transformation, where the character string contains the exponent
  ! **************************************************************************************************
  subroutine apply_transformation(obs,sim,transfo,err,message)
    real(rkind), intent(inout) :: obs(:)
    real(rkind), intent(inout) :: sim(:)
    character(*), intent(in)   :: transfo
    integer(i4b), intent(out)  :: err
    character(*), intent(out) :: message
    real(rkind) :: eps,transfoVal
    integer(i4b) :: i,ioStat
    logical :: integerPower

    err=0
    message='apply_transformation/'
    transfoVal=0._rkind

    select case(trim(transfo))
      case ('none','log','boxcox')
      case default
        read(transfo,*,iostat=ioStat) transfoVal
        if(ioStat/=0)then
          message=trim(message)//'unknown transformation "'//trim(transfo)//'"'
          err=20; return
        endif
        if(.not.ieee_is_finite(transfoVal))then
          message=trim(message)//'power exponent must be finite'
          err=20; return
        endif
    end select
    if(size(obs)==0 .or. trim(transfo)=='none') return

    if(trim(transfo)=='log')then
      eps=sum(obs)/(100*size(obs))
      if(any(eps+obs<=0._rkind) .or. any(eps+sim<=0._rkind))then
        message=trim(message)//'log transformation requires positive Q + epsilon'
        err=metric_undefined; return
      endif
      do i=1,size(obs)
        obs(i)=log(eps+obs(i))
        sim(i)=log(eps+sim(i))
      enddo

    else if(trim(transfo)=='boxcox')then
      if(any(obs<0._rkind) .or. any(sim<0._rkind))then
        message=trim(message)//'Box-Cox transformation requires nonnegative values'
        err=metric_undefined; return
      endif
      do i=1,size(obs)
        obs(i)=(obs(i)**0.25_rkind-1._rkind)/0.25_rkind
        sim(i)=(sim(i)**0.25_rkind-1._rkind)/0.25_rkind
      enddo

    else
      eps=0._rkind
      if(transfoVal<0._rkind) eps=sum(obs)/(100*size(obs))
      integerPower=(transfoVal==anint(transfoVal))
      if(.not.integerPower)then
        if(any(eps+obs<0._rkind) .or. any(eps+sim<0._rkind))then
          message=trim(message)//'fractional power requires nonnegative arguments'
          err=metric_undefined; return
        endif
      endif
      if(transfoVal<0._rkind)then
        if(any(eps+obs==0._rkind) .or. any(eps+sim==0._rkind))then
          message=trim(message)//'negative power is undefined at zero'
          err=metric_undefined; return
        endif
      endif
      do i=1,size(obs)
        obs(i)=(eps+obs(i))**transfoVal
        sim(i)=(eps+sim(i))**transfoVal
      enddo
    endif

    if(.not.all(ieee_is_finite(obs)) .or. .not.all(ieee_is_finite(sim)))then
      message=trim(message)//'transformation produced a non-finite value'
      err=metric_undefined
    endif

  end subroutine apply_transformation

  ! **************************************************************************************************
  ! Compute Kling-Gupta efficiency.
  ! **************************************************************************************************
  function get_kge(obs,sim) result(kge)
    real(rkind), intent(in) :: obs(:)
    real(rkind), intent(in) :: sim(:)
    real(rkind) :: kge
    real(rkind) :: r
    real(rkind) :: alpha
    real(rkind) :: beta
    real(rkind) :: meanObs
    real(rkind) :: meanSim
    real(rkind) :: sdObs
    real(rkind) :: sdSim

    ! compute means and standard deviations
    meanObs = sum(obs)/size(obs)
    meanSim = sum(sim)/size(sim)
    sdObs = standard_deviation(obs)
    sdSim = standard_deviation(sim)

    ! compute KGE components
    r     = correlation(sim,obs)
    alpha = sdSim/sdObs
    beta  = meanSim/meanObs

    ! compute KGE
    kge = 1._rkind - sqrt((r-1._rkind)**2 + (alpha-1._rkind)**2 + (beta-1._rkind)**2)

  end function get_kge

  ! **************************************************************************************************
  ! Compute modified Kling-Gupta efficiency.
  ! **************************************************************************************************
  function get_kgep(obs,sim) result(kgep)
    real(rkind), intent(in) :: obs(:)
    real(rkind), intent(in) :: sim(:)
    real(rkind) :: kgep
    real(rkind) :: r
    real(rkind) :: alphaP
    real(rkind) :: beta
    real(rkind) :: meanObs
    real(rkind) :: meanSim
    real(rkind) :: sdObs
    real(rkind) :: sdSim

    ! compute means and standard deviations
    meanObs = sum(obs)/size(obs)
    meanSim = sum(sim)/size(sim)
    sdObs = standard_deviation(obs)
    sdSim = standard_deviation(sim)

    ! compute modified KGE components
    r      = correlation(sim,obs)
    alphaP = (sdSim/meanSim)/(sdObs/meanObs)
    beta   = meanSim/meanObs

    ! compute modified KGE
    kgep = 1._rkind - sqrt((r-1._rkind)**2 + (alphaP-1._rkind)**2 + (beta-1._rkind)**2)

  end function get_kgep

  ! **************************************************************************************************
  ! Compute Nash-Sutcliffe efficiency.
  ! **************************************************************************************************
  function get_nse(obs,sim) result(nse)
    real(rkind), intent(in) :: obs(:)
    real(rkind), intent(in) :: sim(:)
    real(rkind) :: nse

    nse = 1._rkind - sum((obs-sim)**2) / sum((obs-sum(obs)/size(obs))**2)

  end function get_nse

  ! **************************************************************************************************
  ! Compute mean absolute error.
  ! **************************************************************************************************
  function get_mae(obs,sim) result(mae)
    real(rkind), intent(in) :: obs(:)
    real(rkind), intent(in) :: sim(:)
    real(rkind) :: mae

    mae = sum(abs(obs-sim))/size(obs)

  end function get_mae

  ! **************************************************************************************************
  ! Compute root mean square error.
  ! **************************************************************************************************
  function get_rmse(obs,sim) result(rmse)
    real(rkind), intent(in) :: obs(:)
    real(rkind), intent(in) :: sim(:)
    real(rkind) :: rmse

    rmse = sqrt(sum((obs-sim)**2)/size(obs))

  end function get_rmse



  ! **************************************************************************************************
  ! Compute the sample standard deviation.
  ! **************************************************************************************************
  function standard_deviation(x) result(sd)
    real(rkind), intent(in) :: x(:)
    real(rkind) :: sd
    real(rkind) :: meanX
    integer(i4b) :: n

    n = size(x)
    meanX = sum(x)/n
    sd = sqrt(sum((x-meanX)**2)/(n-1))

  end function standard_deviation

  ! **************************************************************************************************
  ! Compute the correlation coefficient.
  ! **************************************************************************************************
  function correlation(x,y) result(r)
    real(rkind), intent(in) :: x(:)
    real(rkind), intent(in) :: y(:)
    real(rkind) :: r
    real(rkind) :: meanX
    real(rkind) :: meanY
    real(rkind) :: sdX
    real(rkind) :: sdY
    integer(i4b) :: n

    n = size(x)
    meanX = sum(x)/n
    meanY = sum(y)/n
    sdX = standard_deviation(x)
    sdY = standard_deviation(y)
    r = sum((x-meanX)*(y-meanY))/((n-1)*sdX*sdY)

  end function correlation



end module metrics

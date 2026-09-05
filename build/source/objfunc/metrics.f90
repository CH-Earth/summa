module metrics

  use nr_type, only: i4b, rkind

  implicit none
  private

  public :: compute_metric

contains

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

    ! remove missing values and conduct any transformations
    call prepare_series(obs,sim,obsUse,simUse,transfo,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! compute selected metric
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

    if(n==0)then
      message=trim(message)//'no valid observation-simulation pairs'
      err=20; return
    endif

    ! retain only valid pairs
    allocate(obsUse(n),simUse(n))

    obsUse = pack(obs,valid)
    simUse = pack(sim,valid)

    ! apply streamflow transformation
    if(trim(transfo)/='none')then
      call apply_transformation(obsUse,simUse,transfo)
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
  subroutine apply_transformation(obs,sim,transfo)

    real(rkind), intent(inout) :: obs(:)
    real(rkind), intent(inout) :: sim(:)

    character(*), intent(in)   :: transfo

    real(rkind) :: eps
    real(rkind) :: transfoVal

    integer(i4b) :: i

    ! log transformation
    if(trim(transfo)=='log')then

      eps = sum(obs)/(100*size(obs))

      do i=1,size(obs)
        obs(i) = log(eps+obs(i))
        sim(i) = log(eps+sim(i))
      enddo

    ! box-cox transformation
    else if(trim(transfo)=='boxcox')then

      do i=1,size(obs)
        obs(i) = (obs(i)**0.25_rkind-1._rkind)/0.25_rkind
        sim(i) = (sim(i)**0.25_rkind-1._rkind)/0.25_rkind
      enddo

    ! power transformation
    else

      transfoVal = char_to_float(transfo)

      if(transfoVal<0._rkind)then

        eps = sum(obs)/(100*size(obs))

        do i=1,size(obs)
          obs(i) = (eps+obs(i))**transfoVal
          sim(i) = (eps+sim(i))**transfoVal
        enddo

      else

        do i=1,size(obs)
          obs(i) = obs(i)**transfoVal
          sim(i) = sim(i)**transfoVal
        enddo

      endif

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

    meanObs = sum(obs)/size(obs)
    meanSim = sum(sim)/size(sim)

    sdObs = standard_deviation(obs)
    sdSim = standard_deviation(sim)

    r     = correlation(sim,obs)
    alpha = sdSim/sdObs
    beta  = meanSim/meanObs

    kge = 1._rkind - sqrt((r-1._rkind)**2 + &
                          (alpha-1._rkind)**2 + &
                          (beta-1._rkind)**2)

    if(is_nan(kge)) kge = -1.e6_rkind

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

    meanObs = sum(obs)/size(obs)
    meanSim = sum(sim)/size(sim)

    sdObs = standard_deviation(obs)
    sdSim = standard_deviation(sim)

    r      = correlation(sim,obs)
    alphaP = (sdSim/meanSim)/(sdObs/meanObs)
    beta   = meanSim/meanObs

    kgep = 1._rkind - sqrt((r-1._rkind)**2 + &
                           (alphaP-1._rkind)**2 + &
                           (beta-1._rkind)**2)

    if(is_nan(kgep)) kgep = -1.e6_rkind

  end function get_kgep


  ! **************************************************************************************************
  ! Compute Nash-Sutcliffe efficiency.
  ! **************************************************************************************************
  function get_nse(obs,sim) result(nse)

    real(rkind), intent(in) :: obs(:)
    real(rkind), intent(in) :: sim(:)

    real(rkind) :: nse

    nse = 1._rkind - &
          sum((obs-sim)**2) / &
          sum((obs-sum(obs)/size(obs))**2)

    if(is_nan(nse)) nse = -1.e6_rkind

  end function get_nse


  ! **************************************************************************************************
  ! Compute mean absolute error.
  ! **************************************************************************************************
  function get_mae(obs,sim) result(mae)

    real(rkind), intent(in) :: obs(:)
    real(rkind), intent(in) :: sim(:)

    real(rkind) :: mae

    mae = sum(abs(obs-sim))/size(obs)

    if(is_nan(mae)) mae = 1.e6_rkind

  end function get_mae


  ! **************************************************************************************************
  ! Compute root mean square error.
  ! **************************************************************************************************
  function get_rmse(obs,sim) result(rmse)

    real(rkind), intent(in) :: obs(:)
    real(rkind), intent(in) :: sim(:)

    real(rkind) :: rmse

    rmse = sqrt(sum((obs-sim)**2)/size(obs))

    if(is_nan(rmse)) rmse = 1.e6_rkind

  end function get_rmse


  ! **************************************************************************************************
  ! Convert a character value to a floating-point value.
  ! **************************************************************************************************
  function char_to_float(charVal) result(floatVal)

    character(*), intent(in) :: charVal

    real(rkind) :: floatVal

    integer(i4b) :: ioStat

    read(charVal,*,iostat=ioStat) floatVal

    if(ioStat/=0)then
      print *, "Error: Unable to convert '",trim(charVal),"' to float"
      floatVal = 1._rkind
    endif

  end function char_to_float


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


  ! **************************************************************************************************
  ! Return true when a floating-point value is NaN.
  ! **************************************************************************************************
  elemental function is_nan(x) result(isNaN)

    real(rkind), intent(in) :: x

    logical :: isNaN

    isNaN = (x/=x)

  end function is_nan

end module metrics

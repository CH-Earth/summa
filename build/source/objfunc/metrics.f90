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
  ! Optional streamflow transformations are applied before evaluating the metric.
  ! **************************************************************************************************
  subroutine compute_metric(obs,sim,metric,transfo,objective,err,message)

    real(rkind), intent(in)             :: obs(:)       ! observed streamflow
    real(rkind), intent(in)             :: sim(:)       ! simulated streamflow

    character(*), intent(in)            :: metric       ! objective-function metric
    character(*), intent(in)            :: transfo      ! streamflow transformation

    real(rkind), intent(out)            :: objective    ! objective-function value

    integer(i4b), intent(out)           :: err          ! error code
    character(*), intent(out)           :: message      ! error message

    err = 0
    message = 'compute_metric/'

    select case(trim(metric))

      case ('kge');  objective = get_kge(obs,sim,transfo)
      case ('kgep'); objective = get_kgep(obs,sim,transfo)
      case ('nse');  objective = get_nse(obs,sim,transfo)
      case ('mae');  objective = get_mae(obs,sim,transfo)
      case ('rmse'); objective = get_rmse(obs,sim,transfo)

      case default
        message=trim(message)//'unknown objective-function metric "'//trim(metric)//'"'
        err=20; return

    end select

  end subroutine compute_metric


  ! **************************************************************************************************
  ! Prepare observed and simulated streamflow for metric calculation.
  !
  ! Missing values are removed and the selected transformation is applied to
  ! both observed and simulated streamflow.
  ! **************************************************************************************************
  subroutine prepare_series(obs,sim,obsUse,simUse,transfo)

    real(rkind), intent(in)               :: obs(:)
    real(rkind), intent(in)               :: sim(:)

    real(rkind), allocatable, intent(out) :: obsUse(:)
    real(rkind), allocatable, intent(out) :: simUse(:)

    character(*), intent(in)              :: transfo

    logical, allocatable :: valid(:)

    integer(i4b) :: n

    valid = .not.(is_nan(obs) .or. is_nan(sim))

    n = count(valid)

    allocate(obsUse(n),simUse(n))

    obsUse = pack(obs,valid)
    simUse = pack(sim,valid)

    if(trim(transfo)/='none') &
      call apply_transformation(obsUse,simUse,transfo)

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
  function get_kge(obs,sim,transfo) result(kge)

    real(rkind), intent(in) :: obs(:)
    real(rkind), intent(in) :: sim(:)

    character(*), intent(in) :: transfo

    real(rkind) :: kge
    real(rkind) :: r
    real(rkind) :: alpha
    real(rkind) :: beta
    real(rkind) :: meanObs
    real(rkind) :: meanSim
    real(rkind) :: sdObs
    real(rkind) :: sdSim

    real(rkind), allocatable :: obsUse(:)
    real(rkind), allocatable :: simUse(:)

    call prepare_series(obs,sim,obsUse,simUse,transfo)

    meanObs = sum(obsUse)/size(obsUse)
    meanSim = sum(simUse)/size(simUse)

    sdObs = standard_deviation(obsUse)
    sdSim = standard_deviation(simUse)

    r     = correlation(simUse,obsUse)
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
  function get_kgep(obs,sim,transfo) result(kgep)

    real(rkind), intent(in) :: obs(:)
    real(rkind), intent(in) :: sim(:)

    character(*), intent(in) :: transfo

    real(rkind) :: kgep
    real(rkind) :: r
    real(rkind) :: alphaP
    real(rkind) :: beta
    real(rkind) :: meanObs
    real(rkind) :: meanSim
    real(rkind) :: sdObs
    real(rkind) :: sdSim

    real(rkind), allocatable :: obsUse(:)
    real(rkind), allocatable :: simUse(:)

    call prepare_series(obs,sim,obsUse,simUse,transfo)

    meanObs = sum(obsUse)/size(obsUse)
    meanSim = sum(simUse)/size(simUse)

    sdObs = standard_deviation(obsUse)
    sdSim = standard_deviation(simUse)

    r      = correlation(simUse,obsUse)
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
  function get_nse(obs,sim,transfo) result(nse)

    real(rkind), intent(in) :: obs(:)
    real(rkind), intent(in) :: sim(:)

    character(*), intent(in) :: transfo

    real(rkind) :: nse

    real(rkind), allocatable :: obsUse(:)
    real(rkind), allocatable :: simUse(:)

    call prepare_series(obs,sim,obsUse,simUse,transfo)

    nse = 1._rkind - &
          sum((obsUse-simUse)**2) / &
          sum((obsUse-sum(obsUse)/size(obsUse))**2)

    if(is_nan(nse)) nse = -1.e6_rkind

  end function get_nse


  ! **************************************************************************************************
  ! Compute mean absolute error.
  ! **************************************************************************************************
  function get_mae(obs,sim,transfo) result(mae)

    real(rkind), intent(in) :: obs(:)
    real(rkind), intent(in) :: sim(:)

    character(*), intent(in) :: transfo

    real(rkind) :: mae

    real(rkind), allocatable :: obsUse(:)
    real(rkind), allocatable :: simUse(:)

    call prepare_series(obs,sim,obsUse,simUse,transfo)

    mae = sum(abs(obsUse-simUse))/size(obsUse)

    if(is_nan(mae)) mae = 1.e6_rkind

  end function get_mae


  ! **************************************************************************************************
  ! Compute root mean square error.
  ! **************************************************************************************************
  function get_rmse(obs,sim,transfo) result(rmse)

    real(rkind), intent(in) :: obs(:)
    real(rkind), intent(in) :: sim(:)

    character(*), intent(in) :: transfo

    real(rkind) :: rmse

    real(rkind), allocatable :: obsUse(:)
    real(rkind), allocatable :: simUse(:)

    call prepare_series(obs,sim,obsUse,simUse,transfo)

    rmse = sqrt(sum((obsUse-simUse)**2)/size(obsUse))

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

    r = sum((x-meanX)*(y-meanY))/(n*sdX*sdY)

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

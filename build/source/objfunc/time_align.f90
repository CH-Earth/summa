module timeseries_alignment

  USE nr_type, only: i4b, rkind

  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan

  implicit none
  private

  public :: align_timeseries
  public :: select_aligned_window

contains


  ! **************************************************************************************************
  ! Align simulated and observed time series.
  !
  !  The routine:
  !    - converts simulation and observation times to a common absolute time axis
  !    - averages simulated values over each observation period
  !    - preserves the original observation time bounds
  !    - returns the number of contributing simulation steps
  !    - retains missing values for subsequent objective-function masking
  ! **************************************************************************************************
  subroutine align_timeseries(timeSim,valSim,timeSimUnits,valSimUnits,   &
                              timeObs,valObs,timeObsUnits,valObsUnits,   &
                              startDate,endDate,timeBounds,              &
                              timeAligned,timeBoundsAligned,nSimAligned, &
                              valSimAligned,valObsAligned,err,message)
    ! dummy arguments
    real(rkind)   , intent(in)                  :: timeSim(:)          ! simulation time
    real(rkind)   , intent(in)                  :: valSim(:)           ! simulated values
    character(*)  , intent(in)                  :: timeSimUnits        ! simulation time units
    character(*)  , intent(in)                  :: valSimUnits         ! simulation units
    real(rkind)   , intent(in)                  :: timeObs(:)          ! observation time
    real(rkind)   , intent(in)                  :: valObs(:)           ! observation
    character(*)  , intent(in)                  :: timeObsUnits        ! observation time units
    character(*)  , intent(in)                  :: valObsUnits         ! observation units
    character(*)  , intent(in)                  :: startDate           ! start date (YYYY-MM-DD)
    character(*)  , intent(in)                  :: endDate             ! end date (YYYY-MM-DD)
    real(rkind)   , intent(in)                  :: timeBounds(:,:)     ! observation time bounds (2,nObs)
    real(rkind)   , intent(out) , allocatable   :: timeAligned(:)      ! aligned observation time
    real(rkind)   , intent(out) , allocatable   :: timeBoundsAligned(:,:) ! original observation interval boundaries 
    integer(i4b)  , intent(out) , allocatable   :: nSimAligned(:)         ! number of model time steps selected in each interval 
    real(rkind)   , intent(out) , allocatable   :: valSimAligned(:)    ! aligned simulation
    real(rkind)   , intent(out) , allocatable   :: valObsAligned(:)    ! aligned observation
    integer(i4b)  , intent(out)                 :: err                 ! error code
    character(*)  , intent(out)                 :: message             ! error message
    ! locals
    real(rkind), allocatable :: timeSimSec(:)   ! simulation time on common time axis
    real(rkind), allocatable :: timeObsSec(:)   ! observation time on common time axis
    real(rkind), allocatable :: timeBegSec(:)   ! observation time bounds on common time axis
    real(rkind), allocatable :: timeEndSec(:)   ! observation time bounds on common time axis
    real(rkind), allocatable :: simMean(:)      ! simulations averaged to observation periods
    integer(i4b), allocatable :: nSim(:)        ! number of simulation values in each period
    real(rkind) :: tStart                       ! start of observation period
    real(rkind) :: tEnd                         ! end of observation period
    real(rkind) :: tol                          ! time comparison tolerance
    integer(i4b) :: year, month, day, ios       ! check startDate, endDate
    real(rkind) :: evalStart                    ! start of evaluation period
    real(rkind) :: evalEnd                      ! end of evaluation period
    integer(i4b) :: iObs
    integer(i4b) :: iSim
    integer(i4b) :: nMatch
    real(rkind)  :: nanValue
    character(len=256) :: cmessage

    err = 0
    message = 'align_timeseries/'

    ! set NaN
    nanValue = ieee_value(0._rkind, ieee_quiet_nan)

    ! check array dimensions
    if(size(timeSim)/=size(valSim))then
      message=trim(message)//'simulation time and val dimensions differ'
      err=20; return
    endif
    if(size(timeObs)/=size(valObs))then
      message=trim(message)//'observation time and val dimensions differ'
      err=20; return
    endif
    if(size(timeObs)<2)then
      message=trim(message)//'at least two observations are required'
      err=20; return
    endif
    if(size(timeBounds,1)/=2 .or. size(timeBounds,2)/=size(timeObs))then
      message=trim(message)//'observation time bounds must be (2, number of observations)'
      err=20; return
    endif

    ! check start/end dates
    read(startDate,'(i4,1x,i2,1x,i2)',iostat=ios) year,month,day
    if(ios/=0)then
      message=trim(message)//'invalid start date; expected YYYY-MM-DD'
      err=20; return
    endif
    read(endDate,'(i4,1x,i2,1x,i2)',iostat=ios) year,month,day
    if(ios/=0)then
      message=trim(message)//'invalid end date; expected YYYY-MM-DD'
      err=20; return
    endif

    ! Check consistency of simulated and observed units
    if(.not.val_units_equivalent(valSimUnits,valObsUnits))then
      message=trim(message)//'simulation and observation units differ: "'// &
              trim(valSimUnits)//'" and "'//trim(valObsUnits)//'"'
      err=20; return
    endif

    ! convert times to seconds on a common absolute time axis
    call convert_time_to_seconds(timeSim,timeSimUnits,timeSimSec,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    call convert_time_to_seconds(timeObs,timeObsUnits,timeObsSec,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    call convert_time_to_seconds(timeBounds(1,:),timeObsUnits,timeBegSec,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    call convert_time_to_seconds(timeBounds(2,:),timeObsUnits,timeEndSec,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! convert evaluation dates to the same absolute time axis
    call reference_time_seconds(trim(startDate),evalStart,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    call reference_time_seconds(trim(endDate),evalEnd,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! allow small floating-point differences in time coordinates
    tol = max(1.e-6_rkind, &
              1.e-8_rkind*abs(timeEndSec(1)-timeBegSec(1)))

    ! observation times may be irregular (a monthly product's months are 28 to 31 days long), but must increase
    do iObs=2,size(timeObsSec)
      if(timeObsSec(iObs)-timeObsSec(iObs-1) <= 0._rkind)then
        message=trim(message)//'observation times are not increasing'
        err=20; return
      endif
    enddo

    ! allocate simulation values aggregated to observation periods
    allocate(simMean(size(timeObs)))
    allocate(nSim(size(timeObs)))
    simMean = 0._rkind
    nSim    = 0

    ! aggregate simulated val to observation periods
    do iObs=1,size(timeObsSec)

      ! only process observations within the evaluation period
      if(timeObsSec(iObs) >= evalStart .and. timeObsSec(iObs) <  evalEnd)then
     
        ! use actual observation time bounds
        tStart = timeBegSec(iObs)
        tEnd   = timeEndSec(iObs)

        do iSim=1,size(timeSimSec)
          if(timeSimSec(iSim)>tStart+tol .and. timeSimSec(iSim)<=tEnd+tol)then
            simMean(iObs) = simMean(iObs)+valSim(iSim)
            nSim(iObs)    = nSim(iObs)+1
          endif
        enddo
     
        ! compute mean simulated val over the observation period
        if(nSim(iObs)>0)then
          simMean(iObs) = simMean(iObs)/real(nSim(iObs),rkind)
        else
          simMean(iObs) = nanValue
        endif
      endif
    enddo

    ! count valid overlapping periods
    nMatch = 0
    do iObs=1,size(timeObs)
      if(timeObsSec(iObs) >= evalStart .and. timeObsSec(iObs) <  evalEnd)then
        nMatch=nMatch+1
      endif
    enddo
    if(nMatch==0)then
      message=trim(message)//'no overlapping simulation and observation periods'
      err=20; return
    endif

    ! Allocate aligned arrays
    allocate(timeAligned(nMatch))
    allocate(timeBoundsAligned(2,nMatch))
    allocate(nSimAligned(nMatch))
    allocate(valSimAligned(nMatch))
    allocate(valObsAligned(nMatch))

   ! Populate aligned arrays
   nMatch = 0
   
   do iObs=1,size(timeObs)
   
     if(timeObsSec(iObs)>=evalStart .and. &
        timeObsSec(iObs)<evalEnd)then
   
       nMatch = nMatch+1
   
       ! Preserve original observation coordinates and intervals
       timeAligned(nMatch)         = timeObs(iObs)
       timeBoundsAligned(:,nMatch) = timeBounds(:,iObs)
   
       ! Number of simulated values used in this interval
       nSimAligned(nMatch) = nSim(iObs)
   
       ! Aligned observed and simulated values
       valSimAligned(nMatch) = simMean(iObs)
       valObsAligned(nMatch) = valObs(iObs)
   
     endif
   
   enddo

  end subroutine align_timeseries

  ! **************************************************************************************************
  ! Retain only aligned observations within [startDate,endDate).
  ! Allows alignment over the union of baseline and evaluation periods, while ensuring
  ! the objective and diagnostic NetCDF still contain only evaluation observations.
  ! **************************************************************************************************
  subroutine select_aligned_window(time,bounds,nSteps,valSim,valObs,timeUnits,startDate,endDate, &
                                   err,message)
    real(rkind), allocatable, intent(inout) :: time(:),bounds(:,:),valSim(:),valObs(:)
    integer(i4b), allocatable, intent(inout) :: nSteps(:)
    character(*), intent(in) :: timeUnits,startDate,endDate
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message

    real(rkind), allocatable :: timeSec(:)
    real(rkind), allocatable :: timeKeep(:),boundsKeep(:,:),simKeep(:),obsKeep(:)
    integer(i4b), allocatable :: stepsKeep(:)
    real(rkind) :: evalStart,evalEnd
    integer(i4b) :: i,j,nKeep
    character(len=256) :: cmessage

    err=0
    message='select_aligned_window/'

    ! Check dimensions
    if(size(bounds,1)/=2 .or. size(bounds,2)/=size(time) .or. &
       size(nSteps)/=size(time) .or. size(valSim)/=size(time) .or. &
       size(valObs)/=size(time))then
      message=trim(message)//'aligned series lengths differ'
      err=20; return
    endif

    ! Convert observation times to the common absolute time axis
    call convert_time_to_seconds(time,timeUnits,timeSec,err,cmessage)
    if(err/=0)then
      message=trim(message)//trim(cmessage)
      return
    endif

    ! Convert evaluation start date
    call reference_time_seconds(trim(startDate),evalStart,err,cmessage)
    if(err/=0)then
      message=trim(message)//trim(cmessage)
      return
    endif

    ! Convert evaluation end date
    call reference_time_seconds(trim(endDate),evalEnd,err,cmessage)
    if(err/=0)then
      message=trim(message)//trim(cmessage)
      return
    endif

    if(evalEnd<=evalStart)then
      message=trim(message)//'evaluation end must be later than start'
      err=20; return
    endif

    ! Count observations inside the evaluation period
    nKeep=count(timeSec>=evalStart .and. timeSec<evalEnd)

    if(nKeep==0)then
      message=trim(message)//'no observations within evaluation period'
      err=20; return
    endif

    ! Allocate evaluation-only arrays
    allocate(timeKeep(nKeep),boundsKeep(2,nKeep),stepsKeep(nKeep), &
             simKeep(nKeep),obsKeep(nKeep))

    ! Copy evaluation observations
    j=0
    do i=1,size(time)

      if(timeSec(i)<evalStart .or. timeSec(i)>=evalEnd) cycle

      j=j+1

      timeKeep(j)=time(i)
      boundsKeep(:,j)=bounds(:,i)
      stepsKeep(j)=nSteps(i)
      simKeep(j)=valSim(i)
      obsKeep(j)=valObs(i)

    enddo

    ! Replace the original arrays
    call move_alloc(timeKeep,time)
    call move_alloc(boundsKeep,bounds)
    call move_alloc(stepsKeep,nSteps)
    call move_alloc(simKeep,valSim)
    call move_alloc(obsKeep,valObs)

  end subroutine select_aligned_window


  ! **************************************************************************************************
  ! Convert a CF-style time coordinate to seconds on a common absolute time axis.
  !
  ! Supported units:
  !   seconds since YYYY-MM-DD [...]
  !   minutes since YYYY-MM-DD [...]
  !   hours   since YYYY-MM-DD [...]
  !   days    since YYYY-MM-DD [...]
  ! **************************************************************************************************
  subroutine convert_time_to_seconds(time,timeUnits,timeSec,err,message)
    real(rkind), intent(in)               :: time(:)
    character(*), intent(in)              :: timeUnits
    real(rkind), allocatable, intent(out) :: timeSec(:)
    integer(i4b), intent(out)             :: err
    character(*), intent(out)             :: message
    real(rkind) :: scale
    real(rkind) :: refSeconds
    character(len=512) :: unitsLower
    character(len=512) :: reference
    integer(i4b) :: iSince

    err = 0
    message = 'convert_time_to_seconds/'
    unitsLower = lower_case(trim(timeUnits))

    ! identify reference-time delimiter
    iSince = index(unitsLower,'since')
    if(iSince==0)then
      message=trim(message)//'unable to parse time units "'//trim(timeUnits)//'"'
      err=20; return
    endif

    ! identify time-unit scaling
    select case(trim(adjustl(unitsLower(:iSince-1))))
      case ('second','seconds'); scale = 1._rkind
      case ('minute','minutes'); scale = 60._rkind
      case ('hour','hours');     scale = 3600._rkind
      case ('day','days');       scale = 86400._rkind
      case default; message=trim(message)//'unsupported time units "'//trim(timeUnits)//'"'; err=20; return
    end select

    ! extract reference date/time
    reference = adjustl(unitsLower(iSince+5:))
    call reference_time_seconds(trim(reference),refSeconds,err,message)
    if(err/=0) return

    ! convert to common absolute seconds
    allocate(timeSec(size(time)))
    timeSec = refSeconds + time*scale

  end subroutine convert_time_to_seconds

  ! **************************************************************************************************
  ! Convert a reference date/time to absolute seconds.
  !
  ! The absolute origin is arbitrary; only consistency between the two time
  ! coordinates is required.
  ! **************************************************************************************************
  subroutine reference_time_seconds(reference,seconds,err,message)
    character(*), intent(in)  :: reference
    real(rkind), intent(out)  :: seconds
    integer(i4b), intent(out) :: err
    character(*), intent(out) :: message
    character(len=512) :: work
    integer(i4b) :: year
    integer(i4b) :: month
    integer(i4b) :: day
    integer(i4b) :: hour
    integer(i4b) :: minute
    integer(i4b) :: ios
    real(rkind) :: second
    integer(i4b) :: nDay

    err = 0
    message = 'reference_time_seconds/'
    work = trim(reference)

    ! replace date/time delimiters with spaces
    call replace_character(work,'-',' ')
    call replace_character(work,':',' ')
    call replace_character(work,'T',' ')
    call replace_character(work,'t',' ')

    ! defaults for date-only reference strings
    hour   = 0
    minute = 0
    second = 0._rkind

    ! first attempt: date and time
    read(work,*,iostat=ios) year,month,day,hour,minute,second

    ! second attempt: date only
    if(ios/=0)then
      hour   = 0
      minute = 0
      second = 0._rkind
      read(work,*,iostat=ios) year,month,day
    endif
    if(ios/=0)then
      message=trim(message)//'unable to parse reference time "'//trim(reference)//'"'
      err=20; return
    endif

    ! convert calendar date to an absolute day number
    nDay = absolute_day(year,month,day)
    seconds = real(nDay,rkind)*86400._rkind + &
              real(hour,rkind)*3600._rkind   + &
              real(minute,rkind)*60._rkind  + second

  end subroutine reference_time_seconds

  ! **************************************************************************************************
  ! Return the number of days preceding the specified Gregorian calendar date.
  ! **************************************************************************************************
  integer(i4b) function absolute_day(year,month,day)
    integer(i4b), intent(in) :: year
    integer(i4b), intent(in) :: month
    integer(i4b), intent(in) :: day
    integer(i4b), parameter :: monthDays(12) = [31,28,31,30,31,30,31,31,30,31,30,31]
    integer(i4b) :: iYear
    integer(i4b) :: iMonth

    absolute_day = 0

    ! complete years
    do iYear=1,year-1
      absolute_day = absolute_day + 365
      if(is_leap_year(iYear)) absolute_day = absolute_day + 1
    enddo

    ! complete months in current year
    do iMonth=1,month-1
      absolute_day = absolute_day + monthDays(iMonth)
      if(iMonth==2 .and. is_leap_year(year)) absolute_day = absolute_day + 1
    enddo

    ! completed days in current month
    absolute_day = absolute_day + day-1

  end function absolute_day

  ! **************************************************************************************************
  ! Determine whether a year is a Gregorian leap year.
  ! **************************************************************************************************
  logical function is_leap_year(year)
    integer(i4b), intent(in) :: year

    is_leap_year = mod(year,4)==0 .and. (mod(year,100)/=0 .or. mod(year,400)==0)

  end function is_leap_year

  ! **************************************************************************************************
  ! Check whether streamflow unit strings describe cubic metres per second.
  ! **************************************************************************************************
  logical function val_units_equivalent(units1,units2)
    character(*), intent(in) :: units1
    character(*), intent(in) :: units2
    character(len=128) :: u1
    character(len=128) :: u2

    u1 = canonical_value_units(units1)
    u2 = canonical_value_units(units2)
    val_units_equivalent = trim(u1)==trim(u2)

  end function val_units_equivalent

  ! **************************************************************************************************
  ! Convert common streamflow-unit spellings to a canonical representation.
  ! **************************************************************************************************
  function canonical_value_units(units) result(canonical)
    character(*), intent(in) :: units
    character(len=128) :: canonical
    character(len=128) :: work

    work = lower_case(trim(units))
    call remove_character(work,' ')
    call remove_character(work,'^')
    select case(trim(work))
      case ('m3/s','m3s-1','m3s^-1'); canonical = 'm3/s'
      case default; canonical = trim(work)
    end select

  end function canonical_value_units

  ! **************************************************************************************************
  ! Convert a character string to lower case.
  ! **************************************************************************************************
  function lower_case(string) result(lower)
    character(*), intent(in) :: string
    character(len=len(string)) :: lower
    integer(i4b) :: i
    integer(i4b) :: ia

    lower = string
    do i=1,len(string)
      ia = iachar(lower(i:i))
      if(ia>=iachar('A') .and. ia<=iachar('Z')) lower(i:i)=achar(ia+32)
    enddo

  end function lower_case

  ! **************************************************************************************************
  ! Replace one character with another in a string.
  ! **************************************************************************************************
  subroutine replace_character(string,old,new)
    character(*), intent(inout) :: string
    character, intent(in)       :: old
    character, intent(in)       :: new
    integer(i4b) :: i

    do i=1,len_trim(string)
      if(string(i:i)==old) string(i:i)=new
    enddo

  end subroutine replace_character

  ! **************************************************************************************************
  ! Remove a character from a string.
  ! **************************************************************************************************
  subroutine remove_character(string,target)
    character(*), intent(inout) :: string
    character, intent(in)       :: target
    character(len=len(string)) :: work
    integer(i4b) :: i
    integer(i4b) :: j

    work = ''
    j = 0
    do i=1,len_trim(string)
      if(string(i:i)/=target)then
        j=j+1
        work(j:j)=string(i:i)
      endif
    enddo
    string = work

  end subroutine remove_character

end module timeseries_alignment

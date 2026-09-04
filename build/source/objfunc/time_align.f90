module timeseries_alignment

  USE nr_type, only: i4b, rkind

  USE, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none
  private

  public :: align_timeseries

contains

  ! **************************************************************************************************
  ! Align simulated and observed streamflow time series.
  !
  ! The routine:
  !   - converts simulation and observation times to a common absolute time axis
  !   - assumes observation timestamps represent the end of each averaging period
  !   - averages higher-frequency simulated flow over each observation period
  !   - removes periods with missing observed or simulated values
  !
  ! Observation time steps are currently assumed to be regular.
  ! **************************************************************************************************
  subroutine align_timeseries(timeSim,flowSim,timeSimUnits,flowSimUnits, &
                              timeObs,flowObs,timeObsUnits,flowObsUnits, &
                              startDate,endDate,                         &
                              timeAligned,flowSimAligned,flowObsAligned, &
                              err,message)

    ! dummy arguments

    real(rkind), intent(in)                   :: timeSim(:)          ! simulation time
    real(rkind), intent(in)                   :: flowSim(:)          ! simulated streamflow
    character(*), intent(in)                  :: timeSimUnits        ! simulation time units
    character(*), intent(in)                  :: flowSimUnits        ! simulation flow units

    real(rkind), intent(in)                   :: timeObs(:)          ! observation time
    real(rkind), intent(in)                   :: flowObs(:)          ! observed streamflow
    character(*), intent(in)                  :: timeObsUnits        ! observation time units
    character(*), intent(in)                  :: flowObsUnits        ! observation flow units

    character(*), intent(in)                  :: startDate           ! start date (YYYY-MM-DD)
    character(*), intent(in)                  :: endDate             ! end date (YYYY-MM-DD)

    real(rkind), allocatable, intent(out)     :: timeAligned(:)      ! aligned observation time
    real(rkind), allocatable, intent(out)     :: flowSimAligned(:)   ! aligned simulated streamflow
    real(rkind), allocatable, intent(out)     :: flowObsAligned(:)   ! aligned observed streamflow

    integer(i4b), intent(out)                 :: err                 ! error code
    character(*), intent(out)                 :: message             ! error message

    ! locals

    real(rkind), allocatable :: timeSimSec(:)   ! simulation time on common time axis
    real(rkind), allocatable :: timeObsSec(:)   ! observation time on common time axis
    real(rkind), allocatable :: simMean(:)      ! simulated flow averaged to observation periods

    integer(i4b), allocatable :: nSim(:)        ! number of simulation values in each period

    real(rkind) :: obsStep                      ! observation time step in seconds
    real(rkind) :: tStart                       ! start of observation period
    real(rkind) :: tEnd                         ! end of observation period
    real(rkind) :: tol                          ! time comparison tolerance

    integer(i4b) :: year, month, day, ios       ! check startDate, endDate
    real(rkind) :: evalStart                    ! start of evaluation period
    real(rkind) :: evalEnd                      ! end of evaluation period

    integer(i4b) :: iObs
    integer(i4b) :: iSim
    integer(i4b) :: nMatch

    character(len=256) :: cmessage

    err = 0
    message = 'align_timeseries/'

    ! check array dimensions

    if(size(timeSim)/=size(flowSim))then
      message=trim(message)//'simulation time and flow dimensions differ'
      err=20; return
    endif

    if(size(timeObs)/=size(flowObs))then
      message=trim(message)//'observation time and flow dimensions differ'
      err=20; return
    endif

    if(size(timeObs)<2)then
      message=trim(message)//'at least two observations are required'
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

    ! check flow units

    if(.not.flow_units_equivalent(flowSimUnits,flowObsUnits))then
      message=trim(message)//'simulation and observation flow units differ: "'// &
              trim(flowSimUnits)//'" and "'//trim(flowObsUnits)//'"'
      err=20; return
    endif

    ! convert times to seconds on a common absolute time axis

    call convert_time_to_seconds(timeSim,timeSimUnits,timeSimSec,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    call convert_time_to_seconds(timeObs,timeObsUnits,timeObsSec,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! convert evaluation dates to the same absolute time axis

    call reference_time_seconds(trim(startDate),evalStart,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    call reference_time_seconds(trim(endDate),evalEnd,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! determine observation time step

    obsStep = timeObsSec(2)-timeObsSec(1)

    if(obsStep<=0._rkind)then
      message=trim(message)//'observation times are not increasing'
      err=20; return
    endif

    ! allow small floating-point differences in time coordinates

    tol = max(1.e-6_rkind,1.e-8_rkind*obsStep)

    ! check that observation time step is regular

    do iObs=2,size(timeObsSec)-1
      if(abs((timeObsSec(iObs+1)-timeObsSec(iObs))-obsStep)>tol)then
        message=trim(message)//'observation time step is not regular'
        err=20; return
      endif
    enddo

    ! allocate simulation values aggregated to observation periods

    allocate(simMean(size(timeObs)))
    allocate(nSim(size(timeObs)))

    simMean = 0._rkind
    nSim    = 0

    ! aggregate simulated flow to observation periods
    !
    ! timestamps are assumed to represent period ending, so observation i
    ! represents the interval:
    !
    !        timeObs(i)-obsStep < t <= timeObs(i)

    do iObs=1,size(timeObsSec)

      ! only process observations within the evaluation period
      if(timeObsSec(iObs) >= evalStart .and. &
         timeObsSec(iObs) <  evalEnd)then
     
        ! observation timestamps are assumed to be period ending
        tStart = timeObsSec(iObs)-obsStep
        tEnd   = timeObsSec(iObs)
     
        do iSim=1,size(timeSimSec)
     
          if(timeSimSec(iSim)>tStart+tol .and. &
             timeSimSec(iSim)<=tEnd+tol)then
     
            if(ieee_is_finite(flowSim(iSim)))then
              simMean(iObs) = simMean(iObs)+flowSim(iSim)
              nSim(iObs)    = nSim(iObs)+1
            endif
     
          endif
     
        enddo
     
        ! compute mean simulated flow over the observation period
        if(nSim(iObs)>0) &
          simMean(iObs)=simMean(iObs)/real(nSim(iObs),rkind)
     
      endif

    enddo

    ! count valid overlapping periods

    nMatch = 0

    do iObs=1,size(timeObs)

      if(timeObsSec(iObs) >= evalStart .and. &
         timeObsSec(iObs) <  evalEnd   .and. &
         nSim(iObs)>0 .and. ieee_is_finite(flowObs(iObs)))then
        nMatch=nMatch+1
      endif

    enddo

    if(nMatch==0)then
      message=trim(message)//'no overlapping simulation and observation periods'
      err=20; return
    endif

    ! allocate aligned arrays

    allocate(timeAligned(nMatch))
    allocate(flowSimAligned(nMatch))
    allocate(flowObsAligned(nMatch))

    ! populate aligned arrays

    nMatch = 0

    do iObs=1,size(timeObs)

      if(nSim(iObs)>0 .and. ieee_is_finite(flowObs(iObs)))then

        nMatch=nMatch+1

        ! retain the original observation time coordinate
        timeAligned(nMatch)    = timeObs(iObs)
        flowSimAligned(nMatch) = simMean(iObs)
        flowObsAligned(nMatch) = flowObs(iObs)

      endif

    enddo

  end subroutine align_timeseries


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

      case ('second','seconds')
        scale = 1._rkind

      case ('minute','minutes')
        scale = 60._rkind

      case ('hour','hours')
        scale = 3600._rkind

      case ('day','days')
        scale = 86400._rkind

      case default
        message=trim(message)//'unsupported time units "'//trim(timeUnits)//'"'
        err=20; return

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

    integer(i4b), parameter :: monthDays(12) = &
      [31,28,31,30,31,30,31,31,30,31,30,31]

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

      if(iMonth==2 .and. is_leap_year(year)) &
        absolute_day = absolute_day + 1

    enddo

    ! completed days in current month

    absolute_day = absolute_day + day-1

  end function absolute_day


  ! **************************************************************************************************
  ! Determine whether a year is a Gregorian leap year.
  ! **************************************************************************************************
  logical function is_leap_year(year)

    integer(i4b), intent(in) :: year

    is_leap_year = mod(year,4)==0 .and. &
                   (mod(year,100)/=0 .or. mod(year,400)==0)

  end function is_leap_year


  ! **************************************************************************************************
  ! Check whether streamflow unit strings describe cubic metres per second.
  ! **************************************************************************************************
  logical function flow_units_equivalent(units1,units2)

    character(*), intent(in) :: units1
    character(*), intent(in) :: units2

    character(len=128) :: u1
    character(len=128) :: u2

    u1 = canonical_flow_units(units1)
    u2 = canonical_flow_units(units2)

    flow_units_equivalent = trim(u1)==trim(u2)

  end function flow_units_equivalent


  ! **************************************************************************************************
  ! Convert common streamflow-unit spellings to a canonical representation.
  ! **************************************************************************************************
  function canonical_flow_units(units) result(canonical)

    character(*), intent(in) :: units

    character(len=128) :: canonical
    character(len=128) :: work

    work = lower_case(trim(units))

    call remove_character(work,' ')
    call remove_character(work,'^')

    select case(trim(work))

      case ('m3/s','m3s-1','m3s^-1')
        canonical = 'm3/s'

      case default
        canonical = trim(work)

    end select

  end function canonical_flow_units


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

      if(ia>=iachar('A') .and. ia<=iachar('Z')) &
        lower(i:i)=achar(ia+32)

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

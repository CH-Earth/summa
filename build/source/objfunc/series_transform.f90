! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
!
! This file is part of SUMMA
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! **************************************************************************************************
! Putting a simulated series and an observed series on the same footing
!
! Not every observation can be compared against a model variable as it stands.  GRACE reports the
! terrestrial water storage of a basin once a month, as a departure in millimetres from a multi-year
! mean; SUMMA carries the rate that storage is changing at, every time step.  Comparing them means
! integrating the rate into a storage and expressing both as departures from the same baseline.
! Averaging to the month is the alignment's job, which averages the simulation over each observation.
!
! Each of those is a separate step here, and each is asked for by the target that needs it, so a
! target compares what it means to compare and nothing happens to a series that did not ask for it.
! **************************************************************************************************
module series_transform

  USE nr_type, only: i4b,rkind,lgt

  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan

  implicit none
  private

  public :: accumulate_series
  public :: remove_baseline_mean
  public :: balance_between_extremes

contains

  ! **************************************************************************************************
  ! Integrate a rate into the quantity it is the rate of.
  !
  ! basin__StorageChange is a storage flux, in kg m-2 s-1, which is millimetres of water per second;
  ! GRACE reports the storage itself.  Multiplying by the step length and running the sum turns one
  ! into the other, in millimetres, counted from zero at the start of the run.
  !
  ! The series starts at zero rather than at the first increment, because what follows expresses it
  ! as a departure from a baseline mean and the constant of integration drops out there.
  ! **************************************************************************************************
  subroutine accumulate_series(time,values,timeUnits,accumulated,err,message)
    real(rkind),              intent(in)  :: time(:)        ! time coordinate
    real(rkind),              intent(in)  :: values(:)      ! the rate to integrate
    character(*),             intent(in)  :: timeUnits      ! units of the time coordinate
    real(rkind), allocatable, intent(out) :: accumulated(:) ! the integral
    integer(i4b),             intent(out) :: err            ! error code
    character(*),             intent(out) :: message        ! error message
    real(rkind) :: scale       ! seconds per unit of the time coordinate
    real(rkind) :: dt          ! length of one step, in seconds
    real(rkind) :: total       ! running integral
    integer(i4b) :: i

    err=0
    message='accumulate_series/'
    if(size(time) /= size(values))then
      message=trim(message)//'the time and value vectors have different lengths'
      err=20; return
    endif

    call time_scale(timeUnits,scale,err,message)
    if(err/=0) return

    allocate(accumulated(size(values)),stat=err)
    if(err/=0)then
      message=trim(message)//'problem allocating the accumulated series'
      return
    endif

    total=0._rkind
    do i=1,size(values)
      if(i > 1)then
        dt=(time(i)-time(i-1))*scale
      else
        dt=0._rkind
      endif
      ! a step with no value adds nothing, so a gap holds the total rather than destroying it
      if(ieee_is_finite(values(i)) .and. dt > 0._rkind) total=total+values(i)*dt
      accumulated(i)=total
    enddo

  end subroutine accumulate_series

  ! **************************************************************************************************
  ! Express a series as a departure from its own mean over a baseline period.
  !
  ! A GRACE anomaly is a departure from a multi-year mean, so the model has to be expressed the same
  ! way, over the same years, before the two can be compared.  Doing it to both sides removes the
  ! constant of integration the accumulated series carries, which is what makes an integrated rate
  ! comparable to a storage anomaly at all.
  !
  ! An empty baseline period means the whole series.
  ! **************************************************************************************************
  subroutine remove_baseline_mean(time,values,timeUnits,baselineStart,baselineEnd,anomaly,err,message)
    real(rkind),              intent(in)  :: time(:)         ! time coordinate
    real(rkind),              intent(in)  :: values(:)       ! values to reference
    character(*),             intent(in)  :: timeUnits       ! units of the time coordinate
    character(*),             intent(in)  :: baselineStart   ! YYYY-MM-DD, or empty for the whole series
    character(*),             intent(in)  :: baselineEnd     ! YYYY-MM-DD, or empty for the whole series
    real(rkind), allocatable, intent(out) :: anomaly(:)      ! the departures
    integer(i4b),             intent(out) :: err             ! error code
    character(*),             intent(out) :: message         ! error message
    real(rkind)  :: scale
    real(rkind)  :: dayStart,dayEnd,dayHere
    real(rkind)  :: total
    integer(i4b) :: refYear,refMonth,refDay
    integer(i4b) :: year,month,day
    integer(i4b) :: i,count
    logical(lgt) :: bounded

    err=0
    message='remove_baseline_mean/'
    if(size(time) /= size(values))then
      message=trim(message)//'the time and value vectors have different lengths'
      err=20; return
    endif

    call time_scale(timeUnits,scale,err,message)
    if(err/=0) return
    call reference_date(timeUnits,refYear,refMonth,refDay,err,message)
    if(err/=0) return

    bounded = len_trim(baselineStart) > 0 .and. len_trim(baselineEnd) > 0
    dayStart=0._rkind
    dayEnd=0._rkind
    if(bounded)then
      call date_to_days(baselineStart,refYear,refMonth,refDay,dayStart,err,message)
      if(err/=0) return
      call date_to_days(baselineEnd,refYear,refMonth,refDay,dayEnd,err,message)
      if(err/=0) return
    endif

    ! the mean over the baseline
    total=0._rkind
    count=0
    do i=1,size(values)
      if(.not.ieee_is_finite(values(i))) cycle
      if(bounded)then
        dayHere=time(i)*scale/86400._rkind
        if(dayHere < dayStart .or. dayHere > dayEnd) cycle
      endif
      total=total+values(i)
      count=count+1
    enddo
    if(count < 1)then
      message=trim(message)//'the baseline period holds no finite values'
      err=20; return
    endif

    allocate(anomaly(size(values)),stat=err)
    if(err/=0)then
      message=trim(message)//'problem allocating the anomaly series'
      return
    endif
    anomaly=values-total/real(count,rkind)

  end subroutine remove_baseline_mean

  ! **************************************************************************************************
  ! The stratigraphic balance of a storage series, sampled at the times of the observed balances.
  !
  ! A glacier's balances are measured between its own extremes: winter from the autumn minimum to the
  ! spring maximum, summer from that maximum to the next minimum, annual from minimum to minimum.  The
  ! simulated extremes are found where the simulation puts them, so a melt season that runs early or
  ! late is scored on its mass and not on its dates:
  !   seasonal  the extreme within windowDays of the observed date that is a turning point - the
  !             maximum in spring, the minimum in autumn - less the opposite extreme in the year before it
  !   annual    the minimum within windowDays of the observed date, less the minimum within windowDays
  !             of a year before that
  ! A balance whose search reaches outside the simulation, or whose turning point is not clear, is NaN.
  ! **************************************************************************************************
  subroutine balance_between_extremes(timeSim,storage,timeSimUnits,timeObs,timeObsUnits, &
                                      balance,windowDays,simAtObs,err,message)
    real(rkind),              intent(in)  :: timeSim(:)     ! simulated time coordinate
    real(rkind),              intent(in)  :: storage(:)     ! simulated storage
    character(*),             intent(in)  :: timeSimUnits   ! units of the simulated time coordinate
    real(rkind),              intent(in)  :: timeObs(:)     ! times of the observed balances
    character(*),             intent(in)  :: timeObsUnits   ! units of the observed time coordinate
    character(*),             intent(in)  :: balance        ! seasonal or annual
    real(rkind),              intent(in)  :: windowDays     ! half-width of the search for an extreme (days)
    real(rkind), allocatable, intent(out) :: simAtObs(:)    ! simulated balance at each observed time
    integer(i4b),             intent(out) :: err            ! error code
    character(*),             intent(out) :: message        ! error message
    real(rkind), parameter    :: yearDays=365._rkind
    real(rkind), allocatable  :: daySim(:),dayObs(:)
    real(rkind)  :: nan
    integer(i4b) :: i,iEnd,iStart,iMax,iMin
    logical(lgt) :: maxTurns,minTurns

    err=0
    message='balance_between_extremes/'
    if(size(timeSim) /= size(storage))then
      message=trim(message)//'the time and value vectors have different lengths'
      err=20; return
    endif
    if(trim(balance)/='seasonal' .and. trim(balance)/='annual')then
      message=trim(message)//'balance "'//trim(balance)//'" is not seasonal or annual'
      err=20; return
    endif

    call absolute_days(timeSim,timeSimUnits,daySim,err,message); if(err/=0) return
    call absolute_days(timeObs,timeObsUnits,dayObs,err,message); if(err/=0) return

    nan=ieee_value(nan,ieee_quiet_nan)
    allocate(simAtObs(size(timeObs)),stat=err)
    if(err/=0)then
      message=trim(message)//'problem allocating the simulated balances'
      return
    endif

    do i=1,size(dayObs)
      simAtObs(i)=nan
      if(trim(balance)=='annual')then
        iEnd=extreme_in(dayObs(i)-windowDays,dayObs(i)+windowDays,.false.)
        if(iEnd==0) cycle
        iStart=extreme_in(daySim(iEnd)-yearDays-windowDays,daySim(iEnd)-yearDays+windowDays,.false.)
      else
        ! the balance ends at whichever extreme turns inside the window; the other lies on its edge
        iMax=extreme_in(dayObs(i)-windowDays,dayObs(i)+windowDays,.true.)
        iMin=extreme_in(dayObs(i)-windowDays,dayObs(i)+windowDays,.false.)
        if(iMax==0 .or. iMin==0) cycle
        maxTurns=turns(iMax,dayObs(i)-windowDays,dayObs(i)+windowDays)
        minTurns=turns(iMin,dayObs(i)-windowDays,dayObs(i)+windowDays)
        if(maxTurns .eqv. minTurns) cycle
        iEnd=merge(iMax,iMin,maxTurns)
        iStart=extreme_in(daySim(iEnd)-yearDays,daySim(iEnd),.not.maxTurns)
      endif
      if(iStart==0) cycle
      simAtObs(i)=storage(iEnd)-storage(iStart)
    enddo

  contains

    ! index of the largest or smallest finite storage in [dayLo,dayHi], or 0 if the simulation does not cover it
    function extreme_in(dayLo,dayHi,isMax) result(ix)
      real(rkind),  intent(in) :: dayLo,dayHi
      logical(lgt), intent(in) :: isMax
      integer(i4b)             :: ix
      real(rkind), parameter   :: tol=1.e-6_rkind
      integer(i4b) :: j

      ix=0
      if(dayLo < daySim(1)-tol .or. dayHi > daySim(size(daySim))+tol) return
      do j=1,size(daySim)
        if(daySim(j) < dayLo-tol .or. daySim(j) > dayHi+tol) cycle
        if(.not.ieee_is_finite(storage(j))) cycle
        if(ix==0)then
          ix=j
        else if(isMax .and. storage(j) > storage(ix))then
          ix=j
        else if(.not.isMax .and. storage(j) < storage(ix))then
          ix=j
        endif
      enddo
    end function extreme_in

    ! whether an extreme lies inside [dayLo,dayHi] by more than a day, rather than on its edge
    pure function turns(ix,dayLo,dayHi) result(inside)
      integer(i4b), intent(in) :: ix
      real(rkind),  intent(in) :: dayLo,dayHi
      logical(lgt)             :: inside
      inside = daySim(ix) > dayLo+1._rkind .and. daySim(ix) < dayHi-1._rkind
    end function turns

  end subroutine balance_between_extremes

  ! ---- PRIVATE HELPERS -----------------------------------------------------------------------------

  ! **************************************************************************************************
  ! Days of a "<unit> since <reference>" time coordinate on the fixed epoch of day_number.
  ! **************************************************************************************************
  subroutine absolute_days(time,timeUnits,days,err,message)
    real(rkind),              intent(in)    :: time(:)
    character(*),             intent(in)    :: timeUnits
    real(rkind), allocatable, intent(out)   :: days(:)
    integer(i4b),             intent(out)   :: err
    character(*),             intent(inout) :: message
    real(rkind)  :: scale
    integer(i4b) :: refYear,refMonth,refDay

    call time_scale(timeUnits,scale,err,message); if(err/=0) return
    call reference_date(timeUnits,refYear,refMonth,refDay,err,message); if(err/=0) return
    days=real(day_number(refYear,refMonth,refDay),rkind)+time*scale/86400._rkind

  end subroutine absolute_days

  ! **************************************************************************************************
  ! Seconds per unit of a "<unit> since <reference>" time coordinate.
  ! **************************************************************************************************
  subroutine time_scale(timeUnits,scale,err,message)
    character(*), intent(in)    :: timeUnits
    real(rkind),  intent(out)   :: scale
    integer(i4b), intent(out)   :: err
    character(*), intent(inout) :: message
    character(len=512) :: work
    integer(i4b) :: iSince

    err=0
    scale=1._rkind
    work=lower_case(trim(timeUnits))
    iSince=index(work,'since')
    if(iSince==0)then
      message=trim(message)//'unable to parse time units "'//trim(timeUnits)//'"'
      err=20; return
    endif

    select case(trim(adjustl(work(:iSince-1))))
      case ('second','seconds'); scale=1._rkind
      case ('minute','minutes'); scale=60._rkind
      case ('hour','hours');     scale=3600._rkind
      case ('day','days');       scale=86400._rkind
      case default
        message=trim(message)//'unsupported time units "'//trim(timeUnits)//'"'
        err=20; return
    end select

  end subroutine time_scale

  ! **************************************************************************************************
  ! The reference date of a "<unit> since <reference>" time coordinate.
  ! **************************************************************************************************
  subroutine reference_date(timeUnits,year,month,day,err,message)
    character(*), intent(in)    :: timeUnits
    integer(i4b), intent(out)   :: year,month,day
    integer(i4b), intent(out)   :: err
    character(*), intent(inout) :: message
    character(len=512) :: work
    integer(i4b) :: iSince
    integer(i4b) :: iostat

    err=0
    year=1900; month=1; day=1
    work=lower_case(trim(timeUnits))
    iSince=index(work,'since')
    if(iSince==0)then
      message=trim(message)//'unable to parse time units "'//trim(timeUnits)//'"'
      err=20; return
    endif
    work=adjustl(work(iSince+5:))

    ! reference dates are written every which way - "1990-01-01", "1990-1-1 0:0:0.0 -0:00" - so the
    ! delimiters become blanks and the leading three numbers are the date
    call blank_delimiters(work)
    read(work,*,iostat=iostat) year,month,day
    if(iostat/=0)then
      message=trim(message)//'unable to parse the reference date in "'//trim(timeUnits)//'"'
      err=20; return
    endif

  end subroutine reference_date

  ! **************************************************************************************************
  ! Replace date and time delimiters with blanks, so a list-directed read can take the numbers.
  ! **************************************************************************************************
  pure subroutine blank_delimiters(text)
    character(*), intent(inout) :: text
    integer(i4b) :: i

    do i=1,len(text)
      select case(text(i:i))
        case ('-',':','T','t','/'); text(i:i)=' '
      end select
    enddo

  end subroutine blank_delimiters

  ! **************************************************************************************************
  ! Days from a reference date to a YYYY-MM-DD date.
  ! **************************************************************************************************
  subroutine date_to_days(text,refYear,refMonth,refDay,days,err,message)
    character(*), intent(in)    :: text
    integer(i4b), intent(in)    :: refYear,refMonth,refDay
    real(rkind),  intent(out)   :: days
    integer(i4b), intent(out)   :: err
    character(*), intent(inout) :: message
    integer(i4b) :: year,month,day
    integer(i4b) :: iostat
    character(len=64) :: work

    err=0
    days=0._rkind
    work=adjustl(text)
    call blank_delimiters(work)
    read(work,*,iostat=iostat) year,month,day
    if(iostat/=0)then
      message=trim(message)//'"'//trim(text)//'" is not a YYYY-MM-DD date'
      err=20; return
    endif

    days=real(day_number(year,month,day)-day_number(refYear,refMonth,refDay),rkind)

  end subroutine date_to_days

  ! **************************************************************************************************
  ! Day number of a proleptic Gregorian date, counted from a fixed epoch.
  ! **************************************************************************************************
  pure function day_number(year,month,day) result(n)
    integer(i4b), intent(in) :: year,month,day
    integer(i4b)             :: n
    integer(i4b) :: y,m,era,yoe,doy,doe

    y=year
    m=month
    if(m <= 2) y=y-1
    if(y >= 0)then
      era=y/400
    else
      era=(y-399)/400
    endif
    yoe=y-era*400
    if(m > 2)then
      doy=(153*(m-3)+2)/5+day-1
    else
      doy=(153*(m+9)+2)/5+day-1
    endif
    doe=yoe*365+yoe/4-yoe/100+doy
    n=era*146097+doe

  end function day_number

  ! **************************************************************************************************
  ! Lower-case a string.
  ! **************************************************************************************************
  pure function lower_case(string) result(lower)
    character(*), intent(in)   :: string
    character(len=len(string)) :: lower
    integer(i4b) :: i,ic

    lower=string
    do i=1,len(string)
      ic=iachar(string(i:i))
      if(ic >= iachar('A') .and. ic <= iachar('Z')) lower(i:i)=achar(ic+32)
    enddo

  end function lower_case

end module series_transform

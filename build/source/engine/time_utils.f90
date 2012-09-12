module time_utils_module
USE nrtype
implicit none
private
public::extractTime
public::compjulday
contains 

 ! *********************************************************************************
 ! new subroutine: extract year/month/day/hour/minute/second from units string
 ! *********************************************************************************
 subroutine extractTime(refdate,iyyy,im,id,ih,imin,dsec,err,message)
 implicit none
 ! dummy variables
 character(*),intent(in)    :: refdate             ! units string (time since...)
 integer(i4b),intent(out)   :: iyyy,im,id,ih,imin  ! time (year/month/day/hour/minute)
 real(dp),intent(out)       :: dsec                ! seconds
 integer(i4b),intent(out)   :: err                 ! error code
 character(*),intent(out)   :: message             ! error message
 ! local variables
 integer(i4b)               :: n                   ! length of the string
 integer(i4b)               :: istart,iend         ! position in string
 ! iniitalize error control
 err=0; message="extractTime/"

 ! get the length of the string
 n      = len_trim(refdate)
 ! move to a position in string past the time units (seconds since , days since , hours since )
 istart = index(refdate,'since')  ! get the index at the beginning of the word "since"
 if (istart>0) then ! if the word "since" exists
  iend   = index(refdate(istart:n)," ")
  istart = istart+iend
 else
  istart=1
 endif

 ! get the year
 call extract(refdate(istart:n),"-",iend,iyyy,err,message); if (err/=0) return
 ! get the month
 istart=istart+iend
 call extract(refdate(istart:n),"-",iend,im,err,message);   if (err/=0) return
 ! get the day
 istart=istart+iend
 call extract(refdate(istart:n)," ",iend,id,err,message);   if (err/=0) return
 ! check if we are at the end of the string
 if (istart+(iend-2)==n) then
  ih=0; imin=0; dsec=0._dp; return
 endif

 ! get the hour (":" at end of hour)
 istart = istart+iend
 call extract(refdate(istart:n),":",iend,ih,err,message);   if (err/=0) return
 ! get the minute (":" at end of hour)
 istart = istart+iend
 call extract(refdate(istart:n),":",iend,imin,err,message); if (err/=0) return
 ! get the second
 istart = istart+iend
 iend   = index(refdate(istart:n)," ")
 read(refdate(istart:n),*) dsec

 contains
  subroutine extract(substring,cdelim,iend,itemp,err,message)
  implicit none
  ! input
  character(*),intent(in)     :: substring  ! sub-string to process
  character(len=1),intent(in) :: cdelim     ! string delimiter
  ! output
  integer(i4b),intent(out)    :: iend       ! index at the end of desired string
  integer(i4b),intent(out)    :: itemp      ! output date
  integer(i4b),intent(out)    :: err        ! error code
  character(*),intent(out)    :: message    ! error message
  ! initialize error code and message
  err=0; message="extract"
  ! identify end-point of string
  iend = index(substring,cdelim)
  ! if sub-string does not exist, assume end is at end of string
  if (iend==0) iend=len_trim(substring)+1
  ! convert string to integer
  read(substring(1:iend-1),*,iostat=err) itemp
  ! read error
  if (err/=0) then
   err=20; message=trim(message)//"unexpectedCharacters/[string='"//trim(substring)//"']"; return
  endif
  end subroutine extract

 end subroutine extractTime

 ! *********************************************************************************
 ! new subroutine: convert date to julian day (units of days)
 ! *********************************************************************************
 subroutine compjulday(iyyy,mm,id,ih,imin,dsec,&  ! input
                       juldayss,err,message)      ! output
 USE multiconst,only:secprday,secprhour,secprmin  ! seconds in an (day, hour, minute)
 implicit none
 ! input variables
 integer(i4b),intent(in)   :: iyyy,mm,id   ! year, month, day
 integer(i4b),intent(in)   :: ih,imin      ! hour, minute
 real(dp),intent(in)       :: dsec         ! seconds
 ! output
 real(dp),intent(out)      :: juldayss
  integer(i4b),intent(out) :: err          ! error code
  character(*),intent(out) :: message      ! error message
 ! local variables
 integer(i4b)              :: julday       ! julian day
 integer(i4b),parameter    :: igreg=15+31*(10+12*1582)  !IGREG = 588829
 integer(i4b)              :: ja,jm,jy
 real(dp)                  :: jfrac        ! fraction of julian day

 ! initialize errors
 err=0; message="f-juldayss"

 ! compute julian day
 jy=iyyy
 if (jy.eq.0) then; err=10; message=trim(message)//"noYearZero/"; return; endif
 if (jy.lt.0) jy=jy+1
 if (mm.gt.2) then
  jm=mm+1
 else
  jy=jy-1
  jm=mm+13
 endif
 julday=int(365.25*jy)+int(30.6001*jm)+id+1720995
 if (id+31*(mm+12*iyyy).ge.IGREG) then
  ja=int(0.01*jy)
  julday=julday+2-ja+int(0.25*ja)
 endif

 ! compute fraction of the day
 jfrac = (real(ih,kind(dp))*secprhour + real(imin,kind(dp))*secprmin + dsec) / secprday

 ! and return the julian day, expressed in fraction of a day
 juldayss = real(julday,kind(dp)) + jfrac

 end subroutine compjulday



end module time_utils_module

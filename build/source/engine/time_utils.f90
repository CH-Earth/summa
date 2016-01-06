! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
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

module time_utils_module
USE nrtype
implicit none
private
public::extractTime
public::compjulday
contains


 ! ******************************************************************************************
 ! public subroutine extractTime: extract year/month/day/hour/minute/second from units string
 ! ******************************************************************************************
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
 if(iyyy < 1900)then; err=20; message=trim(message)//'year < 1900'; return; endif
 if(iyyy > 2100)then; err=20; message=trim(message)//'year > 2100'; return; endif
 ! get the month
 istart=istart+iend
 call extract(refdate(istart:n),"-",iend,im,err,message);   if (err/=0) return
 if(im <  1)then; err=20; message=trim(message)//'month < 1'; return; endif
 if(im > 12)then; err=20; message=trim(message)//'month > 12'; return; endif
 ! get the day
 istart=istart+iend
 call extract(refdate(istart:n)," ",iend,id,err,message);   if (err/=0) return
 if(id <  1)then; err=20; message=trim(message)//'day < 1'; return; endif
 if(id > 31)then; err=20; message=trim(message)//'day > 31'; return; endif
 ! check if we are at the end of the string
 if (istart+(iend-2)==n) then
  ih=0; imin=0; dsec=0._dp; return
 endif

 ! get the hour (":" at end of hour)
 istart = istart+iend
 if(istart > len_trim(refdate))then; err=20; message=trim(message)//'string does not include hours'; return; endif
 call extract(refdate(istart:n),":",iend,ih,err,message);   if (err/=0) return
 if(ih <  0)then; err=20; message=trim(message)//'hour < 0'; return; endif
 if(ih > 24)then; err=20; message=trim(message)//'hour > 24'; return; endif
 ! get the minute (":" at end of minute)
 istart = istart+iend
 if(istart > len_trim(refdate))then; err=20; message=trim(message)//'string does not include minutes'; return; endif
 call extract(refdate(istart:n),":",iend,imin,err,message); if (err/=0) return
 if(imin <  0)then; err=20; message=trim(message)//'minute < 0'; return; endif
 if(imin > 60)then; err=20; message=trim(message)//'minute > 60'; return; endif
 ! get the second
 istart = istart+iend
 if(istart > len_trim(refdate)) return
 iend   = index(refdate(istart:n)," ")
 read(refdate(istart:n),*) dsec

 contains


  ! ******************************************************************************************
  ! internal subroutine extract: extract substring
  ! ******************************************************************************************
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
  err=0; message="extract/"
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


 ! ***************************************************************************************
 ! public subroutine compjulday: convert date to julian day (units of days)
 ! ***************************************************************************************
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
 err=0; message="juldayss"

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

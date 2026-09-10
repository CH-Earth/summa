! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
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

module read_force_module

! data types
USE nr_type                                   ! variable types, etc.
USE netcdf                                    ! netcdf routines

! derived data types
USE data_types,only:var_ilength               ! x%var(:)%dat(:)            (i4b)
USE data_types,only:gru_hru_double            ! x%gru(:)%hru(:)%var(:)     (rkind)
USE data_types,only:model_options             ! defines the model decisions

! constants
USE multiconst,only:secprday                  ! number of seconds in a day

! access missing values
USE globalData,only:realMissing               ! real missing value
USE globalData,only:integerMissing            ! integer missing value

! access the mapping betweeen GRUs and HRUs
USE globalData,only:gru_struc                 ! gru-hru mapping structures

! global data on the forcing file
USE globalData,only:numtim                    ! number time steps
USE globalData,only:data_step                 ! length of the data step (s)
USE globalData,only:forcFileInfo              ! forcing file info
USE globalData,only:dJulianStart              ! julian day of start time of simulation
USE globalData,only:refJulDay                 ! reference time (fractional julian days)
USE globalData,only:refJulDay_data            ! reference time for data files (fractional julian days)
USE globalData,only:fracJulDay                ! fractional julian days since the start of year
USE globalData,only:yearLength                ! number of days in the current year
USE globalData,only:nHRUfile                  ! number of days in the data file

! global data holding the buffered read
USE globalData,only:ixStartRead               ! start index of the data read
USE globalData,only:fulltimeVec               ! full time vector in an input file (nRead)
USE globalData,only:fullforcingStruct         ! full forcing data structure

! global metadata
USE globalData,only:time_meta,forc_meta       ! metadata structures
USE var_lookup,only:iLookTIME,iLookFORCE      ! named variables to define structure elements
USE var_lookup,only:iLookDECISIONS            ! named variables for elements of the decision structure

! look-up values for option to read input data
USE mDecisions_module,only: &
  readPerStep             , &  ! read forcing data per time step (default)
  readFullSeries               ! read full forcing series

! file paths
USE summaFileManager,only:FORCING_PATH        ! path of the forcing data file

! privacy
implicit none
private
public::read_force

! global parameters
real(rkind),parameter  :: timeDiffTol=1.e-3_rkind   ! time difference tolerance (units=days) to check if the time is correct
real(rkind),parameter  :: smallOffset=1.e-8_rkind   ! small offset (units=days) to force ih=0 at the start of the day

contains

 ! ************************************************************************************************
 ! public subroutine read_force: read in forcing data
 ! ************************************************************************************************
 subroutine read_force(iStep,model_decisions,iFile,iRead,ncid,time_data,forcStruct,err,message)
 ! provide access to subroutines
 USE time_utils_module,only:compJulDay                   ! convert calendar date to julian day
 USE time_utils_module,only:compcalday                   ! convert julian day to calendar date
 USE time_utils_module,only:elapsedSec                   ! calculate the elapsed time
 implicit none
 ! define input variables
 integer(i4b),intent(in)           :: iStep              ! time index AFTER the start index
 type(model_options),intent(in)    :: model_decisions(:) ! model decisions
 ! define input-output variables
 integer(i4b),intent(inout)        :: iFile              ! index of current forcing file in forcing file list
 integer(i4b),intent(inout)        :: iRead              ! index of read position in time dimension in current netcdf file
 integer(i4b),intent(inout)        :: ncid               ! netcdf file identifier
 ! define output variables
 integer(i4b),intent(out)          :: time_data(:)       ! vector of time data for a given time step
 type(gru_hru_double)              :: forcStruct         ! x%gru(:)%hru(:)%var(:)     -- model forcing data
 integer(i4b),intent(out)          :: err                ! error code
 character(*),intent(out)          :: message            ! error message
 ! define local variables
 integer(i4b)                      :: nRemain            ! number of steps remaining in the simulation
 integer(i4b)                      :: nData              ! number of steps remaining in the file
 integer(i4b)                      :: nRead              ! number of steps in the data read
 integer(i4b)                      :: iline              ! loop through lines in the file
 integer(i4b)                      :: jRead              ! index of time in data subset
 integer(i4b)                      :: iGRU,iHRU          ! index of GRU and HRU
 character(len=256),save           :: infile             ! filename
 real(rkind)                       :: dsec               ! double precision seconds (not used)
 real(rkind)                       :: dataJulDay         ! julian day of current forcing data step being read
 real(rkind)                       :: startJulDay        ! julian day at the start of the year
 real(rkind)                       :: currentJulDay      ! Julian day of current time step
 logical(lgt)                      :: isNewFile          ! .true. if reading a new forcing file
 logical(lgt)                      :: isRead             ! .true. if reading data
 logical(lgt),parameter            :: checkTime=.false.  ! flag to check the time
 ! ensure accurate mapping between data structures and HRUs in the output files
 integer(i4b)                      :: iHRU_file_min      ! first required HRU in forcing file
 integer(i4b)                      :: iHRU_file_max      ! last required HRU in forcing file
 integer(i4b)                      :: nHRU_read          ! size of contiguous forcing-file read
 ! error control
 integer(i4b)                      :: ierr               ! local error code
 character(len=256)                :: cmessage           ! error message for downwind routine
 ! Start procedure here
 err=0; message="read_force/"

 ! initialize new file
 isNewFile = .false.

 ! determine forcing-file HRU range needed by this rank
 ! NOTE: forcing-file indices are defined by hru_nc (LocalAttributes ordering)
 !       and are independent of the ordering used in the initial-conditions file.

 iHRU_file_min=huge(1_i4b)
 iHRU_file_max=0

 do iGRU=1,size(gru_struc)
   do iHRU=1,gru_struc(iGRU)%hruCount
     iHRU_file_min=min(iHRU_file_min,gru_struc(iGRU)%hruInfo(iHRU)%hru_nc)
     iHRU_file_max=max(iHRU_file_max,gru_struc(iGRU)%hruInfo(iHRU)%hru_nc)
   enddo
 enddo

 nHRU_read=iHRU_file_max-iHRU_file_min+1

 ! determine the julDay of current model step (iStep) we need to read
 if(iStep==1)then
  currentJulDay = dJulianStart
 else
  currentJulDay = dJulianStart + (data_step*real(iStep-1,dp))/secprday
 end if

#ifdef NGEN_FORCING_ACTIVE
 ! **********************************************************************************************
 ! ***** part 0-1: if using NGEN forcing will be using forcing read with BMI and only need time
 ! **********************************************************************************************
 ! get forcing time data
  call createForcingTimeData(currentJulDay,time_data,err,message)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

#else
 ! **********************************************************************************************
 ! ***** part 0: if initial step, then open first file and find initial model time step
 ! *****         loop through as many forcing files as necessary to find the initial model step
 ! **********************************************************************************************
 ! check if file is open
 if(ncid==integerMissing)then ! file is closed if ncid==integerMissing

  ! identify the first time step
  call getFirstTimestep(currentJulDay,iFile,iRead,ncid,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

  ! flag new file
  isNewFile=.true.

 end if  ! if the file is not yet open

 ! **********************************************************************************************
 ! ***** part 1a: if file open, check to see if we've reached the end of the file, if so close it,
 ! *****          and open new file
 ! **********************************************************************************************
 if(ncid>0)then

  ! check to see if we've passed end of netcdf file
  if(iRead>forcFileInfo(iFile)%nTimeSteps)then

   ! close the NetCDF file
   err = nf90_close(ncid)
   if(err/=nf90_noerr)then; message=trim(message)//'problem closing file ['//trim(infile)//']'; return; endif

   ! increment iFile so we open next forcing file
   iFile = iFile+1
   if(size(forcFileInfo)<iFile)then
    message=trim(message)//'files in list do not include desired data'
    err=20; return
   endif

   ! define new forcing filename
   infile=trim(FORCING_PATH)//trim(forcFileInfo(iFile)%filenmData)

   ! open up the forcing file
   call openForcingFile(iFile,trim(infile),ncid,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

   ! flag new file
   isNewFile=.true.

   ! reset iRead since we opened a new file
   iRead=1

  end if  ! if we've passed the end of the NetCDF file

 ! check that the file was in fact open
 else
  message=trim(message)//'expect the file to be open'
  err=20; return
 end if  ! end ncid open check

 ! **********************************************************************************************
 ! ***** part 1b: allocate space for data structures and read the data
 ! **********************************************************************************************

 ! re-define data length and arrays if new file
 if(isNewFile)then

  ! deallocate space
  if(allocated(fulltimeVec))       deallocate(fulltimeVec)
  if(allocated(fullforcingStruct)) deallocate(fullforcingStruct)

  ! get the number of time steps to read
  select case(model_decisions(iLookDECISIONS%read_force)%iDecision)
   case(readPerStep);    nRead=1     ! ** read forcing data per time step (default)
   case(readFullSeries)              ! ** read full forcing series
     nRemain = (numtim - iStep) + 1                         ! number of remaining time steps in simulation
     nData   = (forcFileInfo(iFile)%nTimeSteps - iRead) + 1 ! number of remaining data steps in file (starting with iRead)
     nRead   = min(nRemain,nData)                           ! number of data steps to read
   case default; err=10; message=trim(message)//'unable to identify option to read forcing data'; return
  end select

  ! allocate space
  allocate(fulltimeVec(nRead), stat=ierr)
  if(ierr/=0)then; err=20; message=trim(message)//'problem allocating space for fulltimeVec'; return; endif
  allocate(fullforcingStruct(nRead), source=forcStruct, stat=ierr)
  if(ierr/=0)then; err=20; message=trim(message)//'problem allocating space for fullforcingStruct'; return; endif

  ! define starting index for the data read
  ixStartRead = iRead

 endif  ! if new file

 ! update ixStartRead and nRead for readPerStep
 if(model_decisions(iLookDECISIONS%read_force)%iDecision == readPerStep)then
   ixStartRead = iRead
   nRead       = 1
  endif

 ! check if we need to read the data
 isRead = (isNewFile .or. model_decisions(iLookDECISIONS%read_force)%iDecision == readPerStep)

 ! read forcing data
 ! NOTE: reads data into global variables fulltimeVec and fullforcingStruct
 if(isRead)then

   call readForcingData(ncid,iFile,ixStartRead,nRead, &
                        iHRU_file_min,nHRU_read,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
 
  endif  ! of reading the forcing data

 ! get the index in the data structures
 jRead = (iRead - ixStartRead) + 1

 ! get forcing structure for the desired time
 forcStruct = fullforcingStruct(jRead)

#endif

 ! **********************************************************************************************
 ! ***** part 2: compute time
 ! **********************************************************************************************

 ! check that the computed julian day matches the time information in the NetCDF file
 dataJulDay = fulltimeVec(jRead)/forcFileInfo(iFile)%convTime2Days + refJulDay_data
 if(abs(currentJulDay - dataJulDay) > timeDiffTol)then
  write(message,'(a,f18.8,a,f18.8)') trim(message)//'date for time step: ',dataJulDay,' differs from the expected date: ',currentJulDay
  err=40; return
 end if

 ! convert julian day to time vector
 ! NOTE: use small offset to force ih=0 at the start of the day
 call compcalday(dataJulDay+smallOffset,         & ! input  = julian day
                 time_data(iLookTIME%iyyy),      & ! output = year
                 time_data(iLookTIME%im),        & ! output = month
                 time_data(iLookTIME%id),        & ! output = day
                 time_data(iLookTIME%ih),        & ! output = hour
                 time_data(iLookTIME%imin),dsec, & ! output = minute/second
                 err,cmessage)                     ! output = error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! check to see if any of the time data is missing -- note that it is OK if ih_tz or imin_tz are missing
 if((time_data(iLookTIME%iyyy)==integerMissing) .or. &
    (time_data(iLookTIME%im)  ==integerMissing) .or. &
    (time_data(iLookTIME%id)  ==integerMissing) .or. &
    (time_data(iLookTIME%ih)  ==integerMissing) .or. &
    (time_data(iLookTIME%imin)==integerMissing) )then
  do iline=1,size(time_data)
   if(time_data(iline)==integerMissing)then; err=40; message=trim(message)//"variableMissing[var='"//trim(time_meta(iline)%varName)//"']"; return; end if
  end do
 end if ! if time data is missing

 ! compute the julian day at the start of the year
 call compjulday(time_data(iLookTIME%iyyy),           & ! input  = year
                 1, 1, 1, 1, 0._rkind,                & ! input  = month, day, hour, minute, second
                 startJulDay,err,cmessage)              ! output = julian day (fraction of day) + error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! compute the fractional julian day for the current time step
 call compjulday(time_data(iLookTIME%iyyy),           & ! input  = year
                 time_data(iLookTIME%im),             & ! input  = month
                 time_data(iLookTIME%id),             & ! input  = day
                 time_data(iLookTIME%ih),             & ! input  = hour
                 time_data(iLookTIME%imin),0._rkind,     & ! input  = minute/second
                 currentJulDay,err,cmessage)            ! output = julian day (fraction of day) + error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if
 ! compute the time since the start of the year (in fractional days)
 fracJulDay = currentJulDay - startJulDay
 ! set timing of current forcing vector (in seconds since reference day)
 ! NOTE: It is a bit silly to have time information for each HRU and GRU
 do iGRU=1,size(gru_struc)
  do iHRU=1,gru_struc(iGRU)%hruCount
   forcStruct%gru(iGRU)%hru(iHRU)%var(iLookFORCE%time) = (currentJulDay-refJulDay)*secprday
  end do  ! looping through HRUs
 end do  ! looping through GRUs

 ! compute the number of days in the current year
 yearLength = 365
 if(mod(time_data(iLookTIME%iyyy),4) == 0)then
  yearLength = 366
  if(mod(time_data(iLookTIME%iyyy),100) == 0)then
   yearLength = 365
   if(mod(time_data(iLookTIME%iyyy),400) == 0)then
    yearLength = 366
   end if
  end if
 end if

 ! test
 if(checkTime)then
  write(*,'(i4,1x,4(i2,1x),f9.3,1x,i4)') time_data(iLookTIME%iyyy),           & ! year
                                         time_data(iLookTIME%im),             & ! month
                                         time_data(iLookTIME%id),             & ! day
                                         time_data(iLookTIME%ih),             & ! hour
                                         time_data(iLookTIME%imin),           & ! minute
                                         fracJulDay,                          & ! fractional julian day for the current time step
                                         yearLength                             ! number of days in the current year
  !pause ' checking time'
 end if

 end subroutine read_force


 ! *************************************************************************
 ! * private subroutine: find first timestep in any of the forcing files...
 ! *************************************************************************
 subroutine getFirstTimestep(currentJulDay,iFile,iRead,ncid,err,message)
 USE nr_utils_module,only:arth                         ! use to build vectors with regular increments
 implicit none
 ! define input
 real(rkind),intent(in)            :: currentJulDay    ! Julian day of current time step
 ! define input-output variables
 integer(i4b),intent(inout)        :: iFile            ! index of current forcing file in forcing file list
 integer(i4b),intent(inout)        :: iRead            ! index of read position in time dimension in current netcdf file
 integer(i4b),intent(inout)        :: ncid             ! netcdf file identifier
 ! define output variables
 integer(i4b),intent(out)          :: err              ! error code
 character(*),intent(out)          :: message          ! error message
 ! ------------------------------------------------------------------------------------------------------------------
 ! netcdf related
 integer(i4b)                      :: varId            ! variable identifier
 integer(i4b)                      :: dimId            ! dimension identifier
 integer(i4b)                      :: dimLen           ! dimension length
 ! other local variables
 character(len=256),save           :: infile           ! filename
 character(len=256)                :: cmessage         ! error message for downwind routine
 integer(i4b)                      :: nFiles           ! number of forcing files
 real(rkind)                          :: timeVal(1)       ! single time value (restrict time read)
 real(rkind),allocatable              :: fileTime(:)      ! array of time from netcdf file
 real(rkind),allocatable              :: diffTime(:)      ! array of time differences
 ! Start procedure here
 err=0; message="getFirstTimestep/"

 ! get the number of forcing files
 nFiles=size(forcFileInfo)  ! number of forcing files

 ! keep going until we find the file containing the first time step
 do iFile=1,nFiles

  ! define new forcing filename
  infile=trim(FORCING_PATH)//trim(forcFileInfo(iFile)%filenmData)

  ! open netCDF file
  call openForcingFile(iFile,trim(infile),ncid,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

  ! how many time steps in current file?
  err = nf90_inq_dimid(ncid,'time',dimId);             if(err/=nf90_noerr)then; message=trim(message)//'trouble finding time dimension/'//trim(nf90_strerror(err)); return; endif
  err = nf90_inquire_dimension(ncid,dimId,len=dimLen); if(err/=nf90_noerr)then; message=trim(message)//'trouble reading time dimension size/'//trim(nf90_strerror(err)); return; endif

  ! how many HRUs in current file?
  err = nf90_inq_dimid(ncid,'hru',dimId);                if(err/=nf90_noerr)then; message=trim(message)//'trouble finding hru dimension/'//trim(nf90_strerror(err)); return; endif
  err = nf90_inquire_dimension(ncid,dimId,len=nHRUfile); if(err/=nf90_noerr)then; message=trim(message)//'trouble reading hru dimension size/'//trim(nf90_strerror(err)); return; endif

  ! allocate space for time vectors
  if(allocated(fileTime)) deallocate(fileTime)
  if(allocated(diffTime)) deallocate(diffTime)
  allocate(fileTime(dimLen),diffTime(dimLen),stat=err)
  if(err/=0)then; message=trim(message)//'problem allocating time vectors'; return; end if

  ! read time vector from current file
  ! NOTE: This could be faster by checking just the start and the end times

  ! (get variable ID)
  err = nf90_inq_varid(ncid,'time',varId)
  if(err/=nf90_noerr)then; message=trim(message)//'trouble finding time variable/'//trim(nf90_strerror(err)); return; endif

  ! (get single time data value)
  err = nf90_get_var(ncid,varId,timeVal,start=(/1/),count=(/1/))
  if(err/=nf90_noerr)then; message=trim(message)//'trouble reading time vector/'//trim(nf90_strerror(err)); return; endif

  ! get time vector & convert units based on offset and data step
  fileTime = arth(0,1,dimLen) * data_step/secprday + refJulDay_data &
             + timeVal(1)/forcFileInfo(iFile)%convTime2Days

  ! find difference of fileTime from currentJulDay
  diffTime=abs(fileTime-currentJulDay)

  ! start time is in the current file
  if(any(diffTime < timeDiffTol))then

   iRead=minloc(diffTime,1)
   exit

  ! time step is not in current file
  else
   ! close file
   err = nf90_close(ncid)
   if(err/=nf90_noerr)then; message=trim(message)//'trouble closing file '//trim(infile); return; endif

   ! check that it is not the last file
   if(iFile==nFiles)then; err=99; message=trim(message)//'first requested simulation timestep not in any forcing file'; return; end if

  end if  ! first time step is not in any forcing files

 end do ! end of search for model first time step in forcing files

 end subroutine getFirstTimestep

 ! *************************************************************************
 ! * open the NetCDF forcing file and get the time information
 ! *************************************************************************
 subroutine openForcingFile(iFile,infile,ncid,err,message)
 USE netcdf_util_module,only:nc_file_open                ! open netcdf file
 USE time_utils_module,only:fracDay                      ! compute fractional day
 USE time_utils_module,only:extractTime                  ! extract time info from units string
 USE time_utils_module,only:compJulDay                   ! convert calendar date to julian day
 USE globalData,only:tmZoneOffsetFracDay                 ! time zone offset in fractional days
 USE globalData,only:ncTime                              ! time zone information from NetCDF file (timeOffset = longitude/15. - ncTimeOffset)
 USE globalData,only:utcTime                             ! all times in UTC (timeOffset = longitude/15. hours)
 USE globalData,only:localTime                           ! all times local (timeOffset = 0)
 USE summafilemanager,only:NC_TIME_ZONE
 ! dummy variables
 integer(i4b),intent(in)           :: iFile              ! index of current forcing file in forcing file list
 character(*) ,intent(in)          :: infile             ! input file
 integer(i4b) ,intent(out)         :: ncid               ! NetCDF ID
 integer(i4b) ,intent(out)         :: err                ! error code
 character(*) ,intent(out)         :: message            ! error message
 ! local variables
 character(len=256)                :: cmessage           ! error message for downwind routine
 integer(i4b)                      :: iyyy,im,id,ih,imin ! date
 integer(i4b)                      :: ih_tz,imin_tz      ! time zone information
 real(rkind)                       :: dsec,dsec_tz       ! seconds
 integer(i4b)                      :: varId              ! variable identifier
 integer(i4b)                      :: mode               ! netcdf file mode
 integer(i4b)                      :: attLen             ! attribute length
 character(len=256)                :: refTimeString      ! reference time string

 ! initialize error control
 err=0; message='openForcingFile/'

 ! open file
 mode=nf90_NoWrite
 call nc_file_open(trim(infile),mode,ncid,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! get definition of time data
 err = nf90_inq_varid(ncid,'time',varId);                       if(err/=nf90_noerr)then; message=trim(message)//'cannot find time variable/'//trim(nf90_strerror(err)); return; endif
 err = nf90_inquire_attribute(ncid,varId,'units',len = attLen); if(err/=nf90_noerr)then; message=trim(message)//'cannot find time units/'//trim(nf90_strerror(err));    return; endif
 err = nf90_get_att(ncid,varid,'units',refTimeString);          if(err/=nf90_noerr)then; message=trim(message)//'cannot read time units/'//trim(nf90_strerror(err));    return; endif

 ! define the reference time for the model simulation
 call extractTime(refTimeString,                         & ! input  = units string for time data
                  iyyy,im,id,ih,imin,dsec,               & ! output = year, month, day, hour, minute, second
                  ih_tz, imin_tz, dsec_tz,               & ! output = time zone information (hour, minute, second)
                  err,cmessage)                            ! output = error code and error message
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 select case(trim(NC_TIME_ZONE))
  case('ncTime'); tmZoneOffsetFracDay = sign(1, ih_tz) * fracDay(ih_tz,   & ! time zone hour
                                                               imin_tz, & ! time zone minute
                                                               dsec_tz)                        ! time zone second
  case('utcTime');   tmZoneOffsetFracDay = 0._rkind
  case('localTime'); tmZoneOffsetFracDay = 0._rkind
  case default; err=20; message=trim(message)//'unable to identify time zone info option'; return
 end select ! (option time zone option)


 ! convert the reference time to days since the beginning of time
 call compjulday(iyyy,im,id,ih,imin,dsec,                & ! input = year, month, day, hour, minute, second
                 refJulDay_data,err,cmessage)              ! output = julian day (fraction of day) + error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 ! get the time multiplier needed to convert time to units of days
 select case( trim( refTimeString(1:index(refTimeString,' ')) ) )
  case('seconds'); forcFileInfo(iFile)%convTime2Days=86400._rkind
  case('minutes'); forcFileInfo(iFile)%convTime2Days=1440._rkind
  case('hours');   forcFileInfo(iFile)%convTime2Days=24._rkind
  case('days');    forcFileInfo(iFile)%convTime2Days=1._rkind
  case default;    message=trim(message)//'unable to identify time units'; err=20; return
 end select

 end subroutine openForcingFile

 ! *************************************************************************
 ! * read the NetCDF forcing data
 ! *************************************************************************
 subroutine readForcingData(ncid,iFile,ixStartRead,nRead, &
                            iHRU_file_min,nHRU_read,err,message)
 ! dummy variables
 integer(i4b) ,intent(in)                :: ncid               ! NetCDF ID
 integer(i4b) ,intent(in)                :: iFile              ! index of forcing file
 integer(i4b) ,intent(in)                :: ixStartRead        ! starting index in data file
 integer(i4b) ,intent(in)                :: nRead              ! number of time steps for the local data read
 integer(i4b) ,intent(in)                :: iHRU_file_min      ! first HRU index in forcing-file read
 integer(i4b) ,intent(in)                :: nHRU_read          ! number of HRU positions in forcing-file read
 integer(i4b) ,intent(out)               :: err                ! error code
 character(*) ,intent(out)               :: message            ! error message
 ! local variables
 integer(i4b)                            :: varId              ! variable identifier
 character(len = nf90_max_name)          :: varName            ! dimenison name
 ! other local variables
 integer(i4b)                            :: iTime              ! time index
 integer(i4b)                            :: iGRU,iHRU          ! index of GRU and HRU
 integer(i4b)                            :: iHRU_file          ! index of HRU in the NetCDF file
 integer(i4b)                            :: iHRU_read          ! index of HRU in the data subset
 integer(i4b)                            :: iline              ! loop through lines in the file
 integer(i4b)                            :: iNC                ! loop through variables in forcing file
 integer(i4b)                            :: iVar               ! index of forcing variable in forcing data vector
 real(rkind),dimension(nHRU_read,nRead)  :: dataMatrix         ! vector of data
 real(rkind),parameter                   :: dataMin=-1._rkind  ! minimum allowable data value (all forcing variables should be positive)
 logical(lgt),dimension(size(forc_meta)) :: checkForce         ! flags to check forcing data variables exist
 ! Start procedure here
 err=0; message="readForcingData/"

 ! read time data starting from iRead location in netcdf file
 err = nf90_inq_varid(ncid,'time',varId);                                          if(err/=nf90_noerr)then; message=trim(message)//'trouble finding time variable/'//trim(nf90_strerror(err)); return; endif
 err = nf90_get_var(ncid,varId,fulltimeVec,start=(/ixStartRead/),count=(/nRead/)); if(err/=nf90_noerr)then; message=trim(message)//'trouble reading time variable/'//trim(nf90_strerror(err)); return; endif

 ! initialize flags for forcing data
 checkForce(:) = .false.
 checkForce(iLookFORCE%time) = .true.  ! time is handled separately

 ! loop through forcing data variables
 do iNC=1,forcFileInfo(iFile)%nVars

  ! check variable is desired
  if(forcFileInfo(iFile)%var_ix(iNC)==integerMissing) cycle

  ! get index in forcing structure
  iVar = forcFileInfo(iFile)%var_ix(iNC)
  checkForce(iVar) = .true.

  ! get variable name for error reporting
  err=nf90_inquire_variable(ncid,iNC,name=varName)
  if(err/=nf90_noerr)then; message=trim(message)//'problem reading forcing variable name from netCDF: '//trim(nf90_strerror(err)); return; endif

  ! read forcing data for the HRU range needed by this rank
  
  err=nf90_get_var(ncid,forcFileInfo(iFile)%data_id(iVar),dataMatrix, &
                   start=[iHRU_file_min,ixStartRead], &
                   count=[nHRU_read,nRead])
  
  if(err/=nf90_noerr)then
    message=trim(message)//'problem reading forcing data: '// &
            trim(varName)//'/'//trim(nf90_strerror(err))
    return
  endif

  ! loop through time
  do iTime=1,nRead

   ! loop through GRUs and HRUs
   do iGRU=1,size(gru_struc)
    do iHRU=1,gru_struc(iGRU)%hruCount
 
     ! map forcing-file HRU index to the local data read
     iHRU_file = gru_struc(iGRU)%hruInfo(iHRU)%hru_nc
     iHRU_read = iHRU_file-iHRU_file_min+1
  
     ! check the HRU indices
 
     if(iHRU_file < 1 .or. iHRU_file > nHRUfile)then
       message=trim(message)//'HRU index is outside the forcing-file HRU dimension'
       err=20; return
     endif
 
     if(iHRU_read < 1 .or. iHRU_read > nHRU_read)then
       message=trim(message)//'HRU index is outside the forcing-data read'
       err=20; return
     endif

     ! check individual data value
     if(dataMatrix(iHRU_read,iTime) < dataMin)then
      write(message,'(a,f13.5)') trim(message)//'forcing data for variable '//trim(varName)//' is less than minimum allowable value ', dataMin
      err=20; return
     endif
  
     ! put the data into structures
     fullforcingStruct(iTime)%gru(iGRU)%hru(iHRU)%var(iVar) = dataMatrix(iHRU_read,iTime)
  
    end do  ! looping through HRUs within a given GRU
   end do  ! looping through GRUs

  end do  ! looping through time

 end do  ! loop through forcing variables

 ! check if any forcing data is missing
 if(count(checkForce)<size(forc_meta))then
  do iline=1,size(forc_meta)
   if(.not.checkForce(iline))then
    message=trim(message)//"variableMissing[var='"//trim(forc_meta(iline)%varName)//"']"
    err=20; return
   endif    ! if variable is missing
  end do   ! looping through variables
 end if   ! if any variables are missing

 end subroutine readForcingData

 ! *************************************************************************
 ! * create forcing time data if using NGEN forcing
 ! *************************************************************************
 subroutine createForcingTimeData(currentJulDay,time_data,err,message)
 USE time_utils_module,only:compcalday                 ! convert julian day to calendar date
 ! dummy variables
 real(rkind), intent(in)           :: currentJulDay    ! Julian day of current time step
 integer(i4b), intent(out)         :: time_data(:)     ! vector of time data for a given time step
 integer(i4b) ,intent(out)         :: err              ! error code
 character(*) ,intent(out)         :: message          ! error message
 ! local variables
 character(len=256)                :: cmessage         ! error message for downwind routine
 ! other local variables
 real(rkind)                       :: dsec             ! double precision seconds (not used)

 ! Start procedure here
 err=0; message="createForcingTimeData/"

 ! initialize time and forcing data structures
 time_data(:) = integerMissing

 ! convert julian day to time vector
 ! NOTE: use small offset to force ih=0 at the start of the day
 call compcalday(currentJulDay+smallOffset,      & ! input  = julian day
                 time_data(iLookTIME%iyyy),      & ! output = year
                 time_data(iLookTIME%im),        & ! output = month
                 time_data(iLookTIME%id),        & ! output = day
                 time_data(iLookTIME%ih),        & ! output = hour
                 time_data(iLookTIME%imin),dsec, & ! output = minute/second
                 err,cmessage)                     ! output = error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; end if

 end subroutine createForcingTimeData


end module read_force_module

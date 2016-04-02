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

module read_force_module
implicit none
private
public::read_force
contains


 ! ************************************************************************************************
 ! public subroutine read_force: read in forcing data
 ! ************************************************************************************************
 subroutine read_force(istep,iGRU,iHRU_local,iHRU_global,iFile,iRead,ncid,time_data,forc_data,err,message)
 ! provide access to subroutines
 USE nrtype                                            ! variable types, etc.
 USE netcdf                                            ! netcdf capability
 USE netcdf_util_module,only:nc_file_open              ! open netcdf file
 USE summaFileManager,only:INPUT_PATH                  ! path of the forcing data file
 USE time_utils_module,only:extractTime                ! extract time info from units string
 USE time_utils_module,only:compJulday                 ! convert calendar date to julian day
 USE time_utils_module,only:compcalday                 ! convert julian day to calendar date
 USE multiconst,only:secprday                          ! number of seconds in a day
 USE globalData,only:forcFileInfo                      ! forcing file info
 USE globalData,only:data_step                         ! length of the data step (s)
 USE globalData,only:dJulianStart                      ! julian day of start time of simulation
 USE globalData,only:refTime,refJulday                 ! reference time
 USE globalData,only:fracJulDay                        ! fractional julian days since the start of year
 USE globalData,only:yearLength                        ! number of days in the current year
 USE globalData,only:time_meta,forc_meta               ! metadata structures
 USE globalData,only:gru_struc                         ! gru-hru mapping structure
 USE var_lookup,only:iLookTIME,iLookFORCE              ! named variables to define structure elements
 USE get_ixname_module,only:get_ixforce                ! identify index of named variable
 USE multiconst,only:integerMissing                    ! integer missing value
 implicit none
 ! define input variables
 integer(i4b),intent(in)           :: istep            ! time index AFTER the start index
 integer(i4b),intent(in)           :: iGRU             ! index of grouped response unit
 integer(i4b),intent(in)           :: iHRU_local       ! index of local hydrologic response unit
 integer(i4b),intent(in)           :: iHRU_global      ! index of global hydrologic response unit
 ! define input-output variables
 integer(i4b),intent(inout)        :: iFile            ! index of current forcing file in forcing file list
 integer(i4b),intent(inout)        :: iRead            ! index of read position in time dimension in current netcdf file
 integer(i4b),intent(inout)        :: ncid             ! netcdf file identifier
 ! define output variables
 integer(i4b),intent(out)          :: time_data(:)     ! vector of time data for a given time step
 real(dp),    intent(out)          :: forc_data(:)     ! vector of forcing data for a given time step
 integer(i4b),intent(out)          :: err              ! error code
 character(*),intent(out)          :: message          ! error message
 ! define local variables
 ! netcdf related
 integer(i4b)                      :: varId            ! variable identifier
 integer(i4b)                      :: dimId            ! dimension identifier
 integer(i4b)                      :: mode             ! netcdf file mode
 integer(i4b)                      :: dimLen           ! dimension length
 integer(i4b)                      :: attLen           ! attribute length
 character(len = nf90_max_name)    :: varName          ! dimenison name
 integer(i4b)                      :: ncStart(2)       ! start array for reading hru forcing
 ! rest
 real(dp),parameter                :: amiss= -1.d+30   ! missing real
 real(dp),parameter                :: verySmall=1e-9   ! tiny number
 character(len=256)                :: infile           ! filename
 character(len=256)                :: cmessage         ! error message for downwind routine
 character(len=256)                :: refTimeString    ! reference time string
 logical(lgt)                      :: xist             ! .TRUE. if the file exists
 integer(i4b),parameter            :: baseUnit=28      ! DK: need to either define units globally, or use getSpareUnit
 integer(i4b)                      :: iline            ! loop through lines in the file
 integer(i4b)                      :: iNC              ! loop through variables in forcing file
 integer(i4b)                      :: iVar             ! index of forcing variable in forcing data vector
 integer(i4b)                      :: iFFile           ! forcing file counter
 integer(i4b)                      :: hruId            ! unique hru id
 real(dp)                          :: dsec             ! double precision seconds (not used)
 real(dp)                          :: startJulDay      ! julian day at the start of the year
 real(dp)                          :: currentJulday    ! Julian day of current time step
 logical(lgt),parameter            :: checkTime=.false.  ! flag to check the time
 real(dp)                          :: dataJulDay       ! julian day of current forcing data step being read
 integer(i4b)                      :: nFiles           ! number of forcing files
 real(dp),allocatable              :: fileTime(:)      ! array of time from netcdf file
 real(dp),allocatable              :: diffTime(:)      ! array of time differences

 ! Start procedure here
 err=0; message="read_force/"

 ! determine the julDay of current model step (istep) we need to read
 if(istep==1)then
  currentJulDay = dJulianStart
 else
  currentJulDay = dJulianStart + (data_step*real(iStep-1,dp))/secprday
 endif

 ! get the number of forcing files
 nFiles=size(forcFileInfo)  ! number of forcing files

 ! **********************************************************************************************
 ! ***** part 0: if initial step, then open first file and find initial model time step
 ! *****         loop through as many forcing files as necessary to find the initial model step
 ! **********************************************************************************************

 ! check if file is open
 if(ncid==integerMissing)then ! file is closed if ncid==integerMissing

  ! ***
  ! * get the time information from the NetCDF file...
  ! **************************************************

  ! define filename
  infile=trim(INPUT_PATH)//trim(forcFileInfo(iFile)%filenmData)

  ! check if the forcing file exists
  inquire(file=trim(infile),exist=xist)
  if(.not.xist)then
   message=trim(message)//"FileNotFound[file='"//trim(infile)//"']"
   err=10; return
  endif

  ! open forcing data file
  mode=nf90_NoWrite
  call nc_file_open(trim(infile),mode,ncid,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! get definition of time data
  err = nf90_inq_varid(ncid,'time',varId);                       if(err/=0)then; message=trim(message)//'cannot find time variable'; return; endif
  err = nf90_inquire_attribute(ncid,varId,'units',len = attLen); if(err/=0)then; message=trim(message)//'cannot find time units';    return; endif
  err = nf90_get_att(ncid,varid,'units',refTimeString);          if(err/=0)then; message=trim(message)//'cannot read time units';    return; endif
  
  ! define the reference time for the model simulation
  call extractTime(refTimeString,                         & ! input  = units string for time data
                   refTime%var(iLookTIME%iyyy),           & ! output = year
                   refTime%var(iLookTIME%im),             & ! output = month
                   refTime%var(iLookTIME%id),             & ! output = day
                   refTime%var(iLookTIME%ih),             & ! output = hour
                   refTime%var(iLookTIME%imin),dsec,      & ! output = minute/second
                   err,cmessage)                            ! output = error code and error message
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! convert the reference time to days since the beginning of time
  call compjulday(refTime%var(iLookTIME%iyyy),            & ! input  = year
                  refTime%var(iLookTIME%im),              & ! input  = month
                  refTime%var(iLookTIME%id),              & ! input  = day
                  refTime%var(iLookTIME%ih),              & ! input  = hour
                  refTime%var(iLookTIME%imin),dsec,       & ! input  = minute/second
                  refJulday,err,cmessage)                   ! output = julian day (fraction of day) + error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! close netCDF file
  err = nf90_close(ncid)
  if(err/=0)then; message=trim(message)//'trouble closing file '//trim(infile); return; endif

  ! ***
  ! * find first timestep in any of the forcing files...
  ! ****************************************************

  ! keep going until we find the file containing the first time step
  do iFFile=1,nFiles

   ! create netCDF file name
   infile=trim(INPUT_PATH)//trim(forcFileInfo(iFFile)%filenmData)
   inquire(file=trim(infile),exist=xist)
   if(.not.xist)then
    message=trim(message)//"FileNotFound[file='"//trim(infile)//"']"
    err=10; return
   endif

   ! open file
   mode=nf90_NoWrite
   call nc_file_open(trim(infile),mode,ncid,err,cmessage)

   ! how many time steps in current file?
   err = nf90_inq_dimid(ncid,'time',dimId);             if(err/=0)then; message=trim(message)//'trouble finding time dimension'; return; endif
   err = nf90_inquire_dimension(ncid,dimId,len=dimLen); if(err/=0)then; message=trim(message)//'trouble reading time dimension size'; return; endif

   ! allocate space for time vectors
   if(allocated(fileTime)) deallocate(fileTime)
   if(allocated(diffTime)) deallocate(diffTime)
   allocate(fileTime(dimLen),diffTime(dimLen),stat=err)
   if(err/=0)then; message=trim(message)//'problem allocating time vectors'; return; endif

   ! read time vector from current file
   ! NOTE: This could be faster by checking just the start and the end times
   err = nf90_get_var(ncid,varId,fileTime,start=(/1/),count=(/dimLen/))
   if(err/=0)then; message=trim(message)//'trouble reading time vector'; return; endif
   fileTime=fileTime+refJulday ! add reference julian day

   ! find difference of fileTime from currentJulday
   diffTime=abs(fileTime-currentJulday)

   ! start time is in the current file
   if(any(diffTime < verySmall))then
    iRead=minloc(diffTime,1)
    iFile=iFFile
    exit

   else ! time step is not in current file

    ! close file
    err = nf90_close(ncid)
    if(err/=0)then; message=trim(message)//'trouble closing file '//trim(infile); return; endif

    ! check that it is not the last file
    if(iFFile==nFiles)then; err=99; message=trim(message)//'first requested simulation timestep not in any forcing file'; return; endif

   endif  ! first time step is not in any forcing files

  end do ! end of search for model first time step in forcing files
   
 endif  ! if the file is not yet open

 ! **********************************************************************************************
 ! ***** part 1: if file open, check to see if we've reached the end of the file, if so close it, 
 ! *****         and open new file
 ! *****         Then read the data
 ! **********************************************************************************************
 if(ncid>0)then

  ! check to see if we've passed end of netcdf file
  if(iRead>forcFileInfo(iFile)%nTimeSteps)then

   ! close the NetCDF file
   err = nf90_close(ncid)
   if(err/=0)then; message=trim(message)//'problem closing file ['//trim(infile)//']'; return; endif

   ! increment iFile so we open next forcing file
   iFile = iFile+1

   ! define new filename, and check that it exists
   infile=trim(INPUT_PATH)//trim(forcFileInfo(iFile)%filenmData)
   inquire(file=trim(infile),exist=xist)
   if(.not.xist)then
    message=trim(message)//"FileNotFound[file='"//trim(infile)//"']"
    err=10; return
   endif

   ! open forcing data file
   mode=nf90_NoWrite
   call nc_file_open(trim(infile),mode,ncid,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! reset iRead since we opened a new file
   iRead=1

  endif  ! if we've passed the end of the NetCDF file

  ! **********************************************************************************************
  ! ***** part 1b: read data
  ! **********************************************************************************************

  ! initialize time and forcing data structures
  time_data(:) = integerMissing
  forc_data(:) = amiss

  ! read time data from iRead location in netcdf file
  err = nf90_inq_varid(ncid,'time',varId);                   if(err/=0)then; message=trim(message)//'trouble finding time variable'; return; endif
  err = nf90_get_var(ncid,varId,dataJulDay,start=(/iRead/)); if(err/=0)then; message=trim(message)//'trouble reading time variable'; return; endif

  ! check that the compted julian day matches the time information in the NetCDF file
  dataJulDay = dataJulDay + refJulday
  if(abs(currentJulday - dataJulDay) > verySmall)then
   write(message,'(a,i0,f18.8,a,f18.8,a)') trim(message)//'date for time step: ',iStep,dataJulDay,' differs from the expected date: ',currentJulDay,' in file: '//trim(infile)
   err=40; return
  endif

  ! convert julian day to time vector
  call compcalday(dataJulDay,                     & ! input  = julian day
                  time_data(iLookTIME%iyyy),      & ! output = year
                  time_data(iLookTIME%im),        & ! output = month
                  time_data(iLookTIME%id),        & ! output = day
                  time_data(iLookTIME%ih),        & ! output = hour
                  time_data(iLookTIME%imin),dsec, & ! output = minute/second
                  err,cmessage)                     ! output = error control

  ! check to see if any of the time data is missing
  if(any(time_data(:)==integerMissing))then
   do iline=1,size(time_data)
    if(time_data(iline)==integerMissing)then; err=40; message=trim(message)//"variableMissing[var='"//trim(time_meta(iline)%varname)//"']"; return; endif
   enddo
  endif

  ! setup count,start arrays
  ncStart = (/iHRU_global,iRead/)

  ! get hruId
  err = nf90_inq_varid(ncid,'hruId',varId);                   if(err/=0)then; message=trim(message)//'trouble finding HRUid variable'; return; endif
  err = nf90_get_var(ncid,varId,hruId,start=(/iHRU_global/)); if(err/=0)then; message=trim(message)//'trouble reading HRUid variable'; return; endif

  ! check that the HRUid is what we expect
  ! NOTE: we enforce that the HRU order in the forcing files is the same as in the zLocalAttributes files (too slow otherwise)
  if(gru_struc(iGRU)%hruInfo(iHRU_local)%hru_id /= hruId)then
   write(message,'(a,i0,i0,a,i0,a)') trim(message)//'hruId for global HRU: ', iHRU_global, hruId, 'differs from the expected:',     &
                                                     gru_struc(iGRU)%hruInfo(iHRU_local)%hru_id, 'in file', trim(infile)
   write(message,'(a)') trim(message)//'order of hruId in forcing file needs to match order in zLocalAttributes.nc'
   err=40; return
  endif

  ! read data into forcing structure
  do iNC=1,forcFileInfo(iFile)%nVars

   ! inqure about current variable name
   err = nf90_inquire_variable(ncid,iNC,name=varName)
   if(err/=0)then; message=trim(message)//'problem inquiring variable: '//trim(varName); return; endif

   ! make sure the variable name is one desired
   select case(trim(varname))
    case('time','pptrate','SWRadAtm','LWRadAtm','airtemp','windspd','airpres','spechum')
    case default; cycle
   end select

   ! get index of forcing variable in forcing data vector
   ivar = get_ixforce(trim(varname))
   if(ivar < 0)then;                                 err=40; message=trim(message)//"variableNotFound [var="//trim(varname)//"]"; return; endif
   if(ivar > size(forcFileInfo(iFile)%data_id))then; err=40; message=trim(message)//"indexOutOfRange  [var="//trim(varname)//"]"; return; endif

   ! get forcing data
   err=nf90_get_var(ncid,forcFileInfo(iFile)%data_id(ivar),forc_data(ivar),start=ncStart)

   ! for time, convert days since reference to seconds since reference
   if(trim(varName)=='time') forc_data(ivar) = forc_data(ivar)*secprday

  end do  ! loop through forcing variables

 ! check that the file was in fact open
 else
  message=trim(message)//'expect the file to be open'
  err=20; return
 endif  ! end ncid open check

 ! **********************************************************************************************
 ! ***** part 2: compute time
 ! **********************************************************************************************

 ! compute the julian day at the start of the year
 call compjulday(time_data(iLookTIME%iyyy),          & ! input  = year
                 1, 1, 1, 1, 0._dp,                  & ! input  = month, day, hour, minute, second
                 startJulDay,err,cmessage)                 ! output = julian day (fraction of day) + error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! compute the fractional julian day for the current time step
 call compjulday(time_data(iLookTIME%iyyy),           & ! input  = year
                 time_data(iLookTIME%im),             & ! input  = month
                 time_data(iLookTIME%id),             & ! input  = day
                 time_data(iLookTIME%ih),             & ! input  = hour
                 time_data(iLookTIME%imin),0._dp,     & ! input  = minute/second
                 currentJulday,err,cmessage)            ! output = julian day (fraction of day) + error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! compute the time since the start of the year (in fractional days)
 fracJulday = currentJulday - startJulDay

 ! compute the number of days in the current year
 yearLength = 365
 if(mod(time_data(iLookTIME%iyyy),4) == 0)then
  yearLength = 366
  if(mod(time_data(iLookTIME%iyyy),100) == 0)then
   yearLength = 365
   if(mod(time_data(iLookTIME%iyyy),400) == 0)then
    yearLength = 366
   endif
  endif
 endif

 ! check to see if any of the forcing data is missing
 ! NOTE: The 0.99 multiplier is used to avoid precision issues when comparing real numbers.
 if(any(forc_data(:)<amiss*0.99_dp))then
  do iline=1,size(forc_meta)
   if(forc_data(iline)<amiss*0.99_dp)then; err=40; message=trim(message)//"variableMissing[var='"//trim(forc_meta(iline)%varname)//"']"; return; endif
  end do
 endif

 ! test
 if(checkTime)then
  write(*,'(i4,1x,4(i2,1x),f9.3,1x,i4)') time_data(iLookTIME%iyyy),           & ! year
                                         time_data(iLookTIME%im),             & ! month
                                         time_data(iLookTIME%id),             & ! day
                                         time_data(iLookTIME%ih),             & ! hour
                                         time_data(iLookTIME%imin),           & ! minute
                                         fracJulday,                          & ! fractional julian day for the current time step
                                         yearLength                             ! number of days in the current year
  !pause ' checking time'
 endif

 end subroutine read_force

end module read_force_module

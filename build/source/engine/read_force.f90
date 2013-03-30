module read_force_module
implicit none
private
public::read_force
contains

 subroutine read_force(istep,err,message)
 ! used to read in initial conditions
 USE nrtype                                            ! variable types, etc.
 USE snow_fileManager,only:INPUT_PATH                  ! path of the forcing data file
 USE allocspace_module,only:alloc_time,alloc_forc      ! allocate space for forcing data
 USE time_utils_module,only:extractTime,compJulday     ! extract time info from units string
 USE multiconst,only:secprday                          ! number of seconds in a day
 USE data_struc,only:forcFileInfo                      ! forcing file info
 USE data_struc,only:refTime,refJulday                 ! reference time
 USE data_struc,only:fracJulDay                        ! fractional julian days since the start of year
 USE data_struc,only:yearLength                        ! number of days in the current year
 USE data_struc,only:time_meta,forc_meta               ! metadata structures
 USE data_struc,only:time_data,forc_data               ! data structures
 USE var_lookup,only:iLookTIME,iLookFORCE              ! named variables to define structure elements 
 implicit none
 ! define dummy variables
 integer(i4b),intent(in)           :: istep            ! time index AFTER the start index
 integer(i4b),intent(out)          :: err              ! error code
 character(*),intent(out)          :: message          ! error message
 ! define local variables
 integer(i4b),parameter            :: imiss= -9999     ! missing integer
 real(dp),parameter                :: amiss= -1.d+30   ! missing real
 character(len=256)                :: infile           ! filename
 character(len=256)                :: cmessage         ! error message for downwind routine
 logical(lgt)                      :: xist             ! .TRUE. if the file exists
 logical(lgt)                      :: xopn             ! .TRUE. if the file is open
 integer(i4b),parameter            :: unt=28           ! DK: need to either define units globally, or use getSpareUnit
 integer(i4b)                      :: iline            ! loop through lines in the file
 character(len=1024),allocatable   :: cline(:)         ! a line of data
 real(dp)                          :: dsec             ! double precision seconds (not used)
 real(dp)                          :: startJulDay      ! julian day at the start of the year
 real(dp)                          :: currentJulday    ! Julian day of current time step
 logical(lgt),parameter            :: checkTime=.false.  ! flag to check the time
 ! local pointers to data structures
 integer(i4b),pointer              :: ncols            ! number of columns in the forcing data file
 integer(i4b),pointer              :: istart           ! start index of the simulation
 integer(i4b),pointer              :: time_ix(:)       ! column index for time
 integer(i4b),pointer              :: data_ix(:)       ! column index for forcing data
 ! Start procedure here
 err=0; message="fuseReadForce/"
 ! define local pointers to data structures
 ncols   => forcFileInfo%ncols   ! number of columns in the forcing data file
 istart  => forcFileInfo%istart  ! start index in the forcing data file
 time_ix => forcFileInfo%time_ix ! column index for time
 data_ix => forcFileInfo%data_ix ! column index for forcing data
 ! define file
 infile=trim(INPUT_PATH)//trim(forcFileInfo%filenm)
 ! check if the forcing info file exists
 inquire(file=trim(infile),exist=xist)  ! Check for existence of forcing datafile
 if(.not.xist)then
   message=trim(message)//"FileNotFound[file='"//trim(infile)//"']"
   err=10; return
 endif
 ! check if the file is open
 inquire(file=trim(infile),opened=xopn) ! Check if the file is open

 ! **********************************************************************************************
 ! ***** part 1: if file not open, then open file and get to the appropriate position in the file
 ! **********************************************************************************************
 if(.not.xopn)then
  ! open forcing data file
  open(unt,file=trim(infile),status="old",action="read",iostat=err)
  if(err/=0)then
   message=trim(message)//"OpenError['"//trim(infile)//"']"
   err=20; return
  endif
  ! allocate space for the time structure
  call alloc_time(time_data,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! allocate space for the forcing structure
  call alloc_forc(err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! allocate space for the reference time
  call alloc_time(refTime,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! define the reference time for the model simulation
  call extractTime(forc_meta(iLookFORCE%time)%varunit,    & ! input  = units string for time data
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
  ! read until just before start index
  do iline=1,istart-1; read(unt,'(a)'); end do
  ! handle situation where istep>1
  if (istep>1) then
   ! read until just before the time step index
   do iline=1,istep-1; read(unt,'(a)'); end do
   ! set a warning message
   err=-20; message="w-"//trim(message)//"UnexpectedFileOpen"
  endif
 endif  ! if the file is not yet open

 ! **********************************************************************************************
 ! ***** part 2: read data
 ! **********************************************************************************************
 
 ! initialize time and forcing data structures
 time_data%var(:) = imiss
 forc_data%var(:) = amiss
 ! allocate space for the character vector
 allocate(cline(ncols),stat=err)
 if (err/=0) then; err=10; message=trim(message)//"problemAllocate"; return; endif
 ! read data using the appropriate format code
 read(unt,*,iostat=err) cline
 if(err/=0)then; err=20; write(message,'(a,i0,a)')trim(message)//"ProblemLineRead[iStep=",istep,"]"; return; endif
 do iline=1,ncols
  !print*,trim(cline(iline))
 end do 
 ! put data in time structure
 do iline=1,size(time_ix)
  if (time_ix(iline)<1 .or. time_ix(iline)>ncols) cycle
  read(cline(time_ix(iline)),*,iostat=err) time_data%var(iline)
  if(err/=0)then; err=30; message=trim(message)//"ProblemTimeRead[var='"//trim(time_meta(iline)%varname)//"']"; return; endif
  !print*,trim(time_meta(iline)%varname),time_data%var(iline)
 end do
 ! check to see if any of the time data is missing
 if(any(time_data%var(:)==imiss))then
  do iline=1,size(time_ix)
   if(time_data%var(iline)==imiss)then; err=40; message=trim(message)//"variableMissing[var='"//trim(time_meta(iline)%varname)//"']"; return; endif
  end do
 endif
 ! put data in forcing structure
 do iline=1,size(data_ix)
  !print*,data_ix(iline)
  if (data_ix(iline)<1 .or. data_ix(iline)>ncols) cycle
  read(cline(data_ix(iline)),*,iostat=err) forc_data%var(iline)
  if(err/=0)then; err=30; message=trim(message)//"ProblemDataRead[var='"//trim(forc_meta(iline)%varname)//"']"; return; endif
  print*,trim(forc_meta(iline)%varname),forc_data%var(iline)
 end do
 ! compute the julian day at the start of the year
 call compjulday(time_data%var(iLookTIME%iyyy),          & ! input  = year
                 1, 1, 1, 1, 0._dp,                      & ! input  = month, day, hour, minute, second
                 startJulDay,err,cmessage)                 ! output = julian day (fraction of day) + error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! compute the fractional julian day for the current time step
 call compjulday(time_data%var(iLookTIME%iyyy),           & ! input  = year
                 time_data%var(iLookTIME%im),             & ! input  = month
                 time_data%var(iLookTIME%id),             & ! input  = day
                 time_data%var(iLookTIME%ih),             & ! input  = hour
                 time_data%var(iLookTIME%imin),0._dp,     & ! input  = minute/second
                 currentJulday,err,cmessage)                ! output = julian day (fraction of day) + error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! compute the time since the start of the year (in fractional days)
 fracJulday = currentJulday - startJulDay
 ! compute time since the reference time (in seconds)
 forc_data%var(iLookFORCE%time) = (currentJulday-refJulday)*secprday
 ! compute the number of days in the current year
 yearLength = 365
 if(mod(time_data%var(iLookTIME%iyyy),4) == 0)then
  yearLength = 366
  if(mod(time_data%var(iLookTIME%iyyy),100) == 0)then
   yearLength = 365
   if(mod(time_data%var(iLookTIME%iyyy),400) == 0)then
    yearLength = 366
   endif
  endif
 endif
 ! check to see if any of the forcing data is missing
 if(any(forc_data%var(:)<amiss*0.99_dp))then
  do iline=1,size(data_ix)
   if(forc_data%var(iline)<amiss*0.99_dp)then; err=40; message=trim(message)//"variableMissing[var='"//trim(forc_meta(iline)%varname)//"']"; return; endif
  end do
 endif
 ! test
 if(checkTime)then
  write(*,'(i4,1x,4(i2,1x),f9.3,1x,i4)') time_data%var(iLookTIME%iyyy),           & ! year
                                         time_data%var(iLookTIME%im),             & ! month
                                         time_data%var(iLookTIME%id),             & ! day
                                         time_data%var(iLookTIME%ih),             & ! hour
                                         time_data%var(iLookTIME%imin),           & ! minute
                                         fracJulday,                              & ! fractional julian day for the current time step
                                         yearLength                                 ! number of days in the current year
  pause ' checking time'
 endif
 ! deallocate cline
 deallocate(cline,stat=err)
 if (err/=0) then; err=10; message=trim(message)//"problemDeallocate"; return; endif
 end subroutine read_force

end module read_force_module

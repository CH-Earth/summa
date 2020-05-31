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

!******************************************************************
! (C) Copyright 2009-2010  ---  Dmitri Kavetski and Martyn Clark ---  All rights reserved
! 2020 substantially rewritten by A. Wood
!******************************************************************
MODULE summafilemanager
use nrtype
!USE var_lookup, only: maxvarControl  ! maximum number of control file entries
implicit none
public
! summa-wide pathlength
integer(i4b),parameter::summaPathLen=4096
! defines the path for data files (and default values)
CHARACTER(LEN=summaPathLen)  :: SETNGS_PATH      = 'settings/'         ! SETNGS_PATH
CHARACTER(LEN=summaPathLen)  :: INPUT_PATH       = 'input/default/'    ! INPUT_PATH
CHARACTER(LEN=summaPathLen)  :: OUTPUT_PATH      = 'output/default/'   ! OUTPUT_PATH
! define name of control files    (and default values)
CHARACTER(LEN=summaPathLen)  :: M_DECISIONS      = 'summa_zDecisions.txt'           ! definition of model decisions
CHARACTER(LEN=summaPathLen)  :: META_TIME        = 'summa_zTimeMeta.txt'            ! metadata for time
CHARACTER(LEN=summaPathLen)  :: META_ATTR        = 'summa_zLocalAttributeMeta.txt'  ! metadata for local attributes
CHARACTER(LEN=summaPathLen)  :: META_TYPE        = 'summa_zCategoryMeta.txt'        ! metadata for local classification of veg, soil, etc.
CHARACTER(LEN=summaPathLen)  :: META_FORCE       = 'summa_zForceMeta.txt'           ! metadata for model forcing variables
CHARACTER(LEN=summaPathLen)  :: META_LOCALPARAM  = 'summa_zLocalParamMeta.txt'      ! metadata for model parameters
CHARACTER(LEN=summaPathLen)  :: OUTPUT_CONTROL   = 'summa_zLocalModelVarMeta.txt'   ! metadata for model variables
CHARACTER(LEN=summaPathLen)  :: META_LOCALINDEX  = 'summa_zLocalModelIndexMeta.txt' ! metadata for model indices
CHARACTER(LEN=summaPathLen)  :: META_BASINPARAM  = 'summa_zBasinParamMeta.txt'      ! metadata for model parameters
CHARACTER(LEN=summaPathLen)  :: META_BASINMVAR   = 'summa_zBasinModelVarMeta.txt'   ! metadata for model variables
CHARACTER(LEN=summaPathLen)  :: LOCAL_ATTRIBUTES = 'summa_zLocalAttributes.txt'     ! local attributes
CHARACTER(LEN=summaPathLen)  :: LOCALPARAM_INFO  = 'summa_zLocalParamInfo.txt'      ! default values and constraints for local model parameters
CHARACTER(LEN=summaPathLen)  :: BASINPARAM_INFO  = 'summa_zBasinParamInfo.txt'      ! default values and constraints for basin model parameters
CHARACTER(LEN=summaPathLen)  :: FORCING_FILELIST = 'summa_zForcingFileList.txt'     ! list of focing files for each HRU
CHARACTER(LEN=summaPathLen)  :: MODEL_INITCOND   = 'summa_zInitialCond.txt'         ! model initial conditions
CHARACTER(LEN=summaPathLen)  :: PARAMETER_TRIAL  = 'summa_zParamTrial.txt'          ! trial values for model parameters
CHARACTER(LEN=summaPathLen)  :: OUTPUT_PREFIX    = 'summa_output_'                  ! prefix for the output file

contains

 ! **************************************************************************************************
 ! public subroutine summa_SetTimesDirsUndPhiles: Sets times, directories and filenames for summa run
 ! **************************************************************************************************
 subroutine summa_SetTimesDirsUndPhiles(summaFileManagerIn,err,message)
 ! Purpose: Sets run times, directories and philenames for summa.
 ! ---
 ! Programmer: Dmitri Kavetski and Martyn Clark
 ! Last modified:  May 30 2020 AW Wood rewrote to move times into file manager, create control structures

 USE globaldata,only:model_control_info     ! model decision structure
 USE multiconst,only:secprday               ! number of seconds in a day
 USE var_lookup,only:iLookTIME              ! named variables that identify indices in the time structures
 USE globalData,only:refTime,refJulday      ! reference time
 USE globalData,only:oldTime                ! time from the previous time step
 USE globalData,only:startTime,finshTime    ! start/end time of simulation
 USE globalData,only:dJulianStart           ! julian day of start time of simulation
 USE globalData,only:dJulianFinsh           ! julian day of end time of simulation
 USE globalData,only:data_step              ! length of data step (s)
 USE globalData,only:numtim                 ! number of time steps in the simulation
 ! forcing metadata
 USE globalData,only:forc_meta              ! metadata structures
 USE var_lookup,only:iLookFORCE             ! named variables to define structure elements
 ! time utility programs
 USE time_utils_module,only:extractTime     ! extract time info from units string
 USE time_utils_module,only:compjulday      ! compute the julian day
 USE time_utils_module,only:fracDay         ! compute fractional day

 implicit none

 ! input/output vars
 character(*),intent(in)    ::  summaFileManagerIn
 integer(i4b),intent(out)   ::  err
 character(*),intent(out)   ::  message
 ! locals
 logical(lgt)               ::  xist
 character(*),parameter     ::  summaFileManagerHeader="SUMMA_FILE_MANAGER_V2.0"
 character(LEN=100)         ::  temp
 character(LEN=64)          ::  varOption             ! option given in file manager
 integer(i4b)               ::  ierr                  ! temporary error code
 integer(i4b),parameter     ::  runinfo_fileunit=67   ! file unit for run time information
 character(len=8)           ::  cdate
 character(len=10)          ::  ctime
 character(len=256)         ::  cmessage       ! error message for downwind routine
 real(dp)                   ::  dsec,dsec_tz   ! second

 ! Start procedure here
 err=0; message="summa_SetTimezDirsUndPhiles/"
 
 ! read information from model control file, and populate model control structure
 ! populates global control information structure
 call read_control_info(summaFileManagerIn, err, cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 ! -------------------------------------------------------------------------------------------------

 ! loop over control file entries, assign to public varnames, add checking if needed
 do iControl=1,nControl
  select case(trim(model_control_info(iControl)%cOption))
   case('controlVersion' ); if(trim(model_control_info(iControl)%cControl)/=summaFileManagerHeader)then; &&
     message=trim(message)//"unknown control file version in '"//trim(summaFileManagerIn)//; err=20; return; end if
   case('settingsPath'   ); SETNGS_PATH = trim(model_control_info(iControl)%cControl)       ! settings directory
   case('forcingPath'    ); INPUT_PATH = trim(model_control_info(iControl)%cControl)        ! input forcing directory
   case('outputPath'     ); OUTPUT_PATH = trim(model_control_info(iControl)%cControl)       ! output directory
   case('decisionsFile'  ); M_DECISIONS = trim(model_control_info(iControl)%cControl)       ! model decisions file
   case('outputDefFile'  ); OUTPUT_CONTROL = trim(model_control_info(iControl)%cControl)    ! output control file
   case('hruParamFile'   ); LOCALPARAM_INFO = trim(model_control_info(iControl)%cControl)   ! default hru-level param file
   case('gruParamFile'   ); BASINPARAM_INFO = trim(model_control_info(iControl)%cControl)   ! default gru-level param file
   case('attributeFile'  ); LOCAL_ATTRIBUTES = trim(model_control_info(iControl)%cControl)  ! attribute file
   case('trialParamFile' ); PARAMETER_TRIAL = trim(model_control_info(iControl)%cControl)   ! trial parameters file
   case('forcingList'    ); FORCING_FILELIST = trim(model_control_info(iControl)%cControl)  ! file listing forcing filenames
   case('initCondFile'   ); MODEL_INITCOND = trim(model_control_info(iControl)%cControl)    ! initial conditions file (cold State)
   case('outFilePrefix'  ); OUTPUT_PREFIX = trim(model_control_info(iControl)%cControl)     ! filename root for output files
   case('metaTime'       ); META_TIME = trim(model_control_info(iControl)%cControl)         ! 
   case('metaAttr'       ); META_ATTR = trim(model_control_info(iControl)%cControl)         ! 
   case('metaType'       ); META_TYPE = trim(model_control_info(iControl)%cControl)         ! 
   case('metaForc'       ); META_FORCE = trim(model_control_info(iControl)%cControl)        ! 
   case('metaLocParam'   ); META_LOCALPARAM = trim(model_control_info(iControl)%cControl)   ! 
   case('metaBasParam'   ); META_BASINPARAM = trim(model_control_info(iControl)%cControl)   ! 
   case('metaBasMvar'    ); META_BASINMVAR = trim(model_control_info(iControl)%cControl)    ! 
   case('metaLocIndex'   ); META_LOCALINDEX = trim(model_control_info(iControl)%cControl)   ! 
   ! get to here if cannot find the variable
   case default
     err=10; message=trim(message)//"unknown control file option"//trim(model_control_info(iControl)%cControl); return
  end select
 end do

! leftover, should not be needed
 ! identify the choice of the time zone option
! select case(trim(model_control_info(iLookCONTROL%tmZoneInfo)%cControl))
!  case('ncTime'   ); model_control_info(iLookCONTROL%tmZoneInfo)%iControl = ncTime       ! time zone information from NetCDF file
!  case('utcTime'  ); model_control_info(iLookCONTROL%tmZoneInfo)%iControl = utcTime      ! all times in UTC
!  case('localTime'); model_control_info(iLookCONTROL%tmZoneInfo)%iControl = localTime    ! all times local
!  case default
!   err=10; message=trim(message)//"unknown time zone info option [option="//trim(model_control(iLookCONTROL%tmZoneInfo)%cControl)//"]"; return
! end select

 ! process time information 

 ! put reference time information into the time structures
 call extractTime(forc_meta(iLookFORCE%time)%varunit,                    & ! date-time string
                  refTime%var(iLookTIME%iyyy),                           & ! year
                  refTime%var(iLookTIME%im),                             & ! month
                  refTime%var(iLookTIME%id),                             & ! day
                  refTime%var(iLookTIME%ih),                             & ! hour
                  refTime%var(iLookTIME%imin),                           & ! minute
                  dsec,                                                  & ! second
                  refTime%var(iLookTIME%ih_tz),                          & ! time zone hour
                  refTime%var(iLookTIME%imin_tz),                        & ! time zone minute
                  dsec_tz,                                               & ! time zone seconds
                  err,cmessage)                                            ! error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if


 ! compute the julian date (fraction of day) for the reference time
 call compjulday(&
                 refTime%var(iLookTIME%iyyy),                           & ! year
                 refTime%var(iLookTIME%im),                             & ! month
                 refTime%var(iLookTIME%id),                             & ! day
                 refTime%var(iLookTIME%ih),                             & ! hour
                 refTime%var(iLookTIME%imin),                           & ! minute
                 0._dp,                                                 & ! second
                 refJulday,                                             & ! julian date for the start of the simulation
                 err, cmessage)                                           ! error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 ! put simulation start time information into the time structures
 call extractTime(model_control_info(iLookCONTROL%simStartTime)%cControl,& ! date-time string
                  startTime%var(iLookTIME%iyyy),                         & ! year
                  startTime%var(iLookTIME%im),                           & ! month
                  startTime%var(iLookTIME%id),                           & ! day
                  startTime%var(iLookTIME%ih),                           & ! hour
                  startTime%var(iLookTIME%imin),                         & ! minute
                  dsec,                                                  & ! second
                  startTime%var(iLookTIME%ih_tz),                        & ! time zone hour
                  startTime%var(iLookTIME%imin_tz),                      & ! time zone minnute
                  dsec_tz,                                               & ! time zone seconds
                  err,cmessage)                                            ! error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 ! compute the julian date (fraction of day) for the start of the simulation
 call compjulday(&
                 startTime%var(iLookTIME%iyyy),                         & ! year
                 startTime%var(iLookTIME%im),                           & ! month
                 startTime%var(iLookTIME%id),                           & ! day
                 startTime%var(iLookTIME%ih),                           & ! hour
                 startTime%var(iLookTIME%imin),                         & ! minute
                 0._dp,                                                 & ! second
                 dJulianStart,                                          & ! julian date for the start of the simulation
                 err, cmessage)                                           ! error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 ! put simulation end time information into the time structures
 call extractTime(model_control_info(iLookCONTROL%simEndTime)%cControl,  & ! date-time string
                  finshTime%var(iLookTIME%iyyy),                         & ! year
                  finshTime%var(iLookTIME%im),                           & ! month
                  finshTime%var(iLookTIME%id),                           & ! day
                  finshTime%var(iLookTIME%ih),                           & ! hour
                  finshTime%var(iLookTIME%imin),                         & ! minute
                  dsec,                                                  & ! second
                  finshTime%var(iLookTIME%ih_tz),                        & ! time zone hour
                  finshTime%var(iLookTIME%imin_tz),                      & ! time zone minnute
                  dsec_tz,                                               & ! time zone seconds
                  err,cmessage)                                            ! error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 ! compute the julian date (fraction of day) for the end of the simulation
 call compjulday(&
                 finshTime%var(iLookTIME%iyyy),                         & ! year
                 finshTime%var(iLookTIME%im),                           & ! month
                 finshTime%var(iLookTIME%id),                           & ! day
                 finshTime%var(iLookTIME%ih),                           & ! hour
                 finshTime%var(iLookTIME%imin),                         & ! minute
                 0._dp,                                                 & ! second
                 dJulianFinsh,                                          & ! julian date for the end of the simulation
                 err, cmessage)                                           ! error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; end if

 ! check start and finish time
 write(*,'(a,i4,1x,4(i2,1x))') 'startTime: iyyy, im, id, ih, imin = ', startTime%var(1:5)
 write(*,'(a,i4,1x,4(i2,1x))') 'finshTime: iyyy, im, id, ih, imin = ', finshTime%var(1:5)

 ! check that simulation end time is > start time
 if(dJulianFinsh < dJulianStart)then; err=20; message=trim(message)//'end time of simulation occurs before start time'; return; end if

 ! initialize the old time vector (time from the previous time step)
 oldTime%var(:) = startTime%var(:)

 ! compute the number of time steps
 numtim = nint( (dJulianFinsh - dJulianStart)*secprday/data_step ) + 1
 write(*,'(a,1x,i10)') 'number of time steps = ', numtim

 ! -------------------------------------------------------------------------------------------------
 
 ! before embarking on a run, check that the output directory is writable; write system date and time to a log file there
 open(runinfo_fileunit,file=trim(OUTPUT_PATH)//"runinfo.txt",iostat=err)
 if(err/=0)then; err=10; message=trim(message)//"cannot write to output directory '"//trim(OUTPUT_PATH)//"'"; return; end if
 call date_and_time(cdate,ctime)
 write(runinfo_fileunit,*) 'Run start time on system:  ccyy='//cdate(1:4)//' - mm='//cdate(5:6)//' - dd='//cdate(7:8), &
                         ' - hh='//ctime(1:2)//' - mi='//ctime(3:4)//' - ss='//ctime(5:10)
 close(runinfo_fileunit)

 ! End procedure here
 end subroutine summa_SetTimesDirsUndPhiles

 ! ************************************************************************************************
 ! private subroutine read_control_info: read information from filemanager / control file
 ! ************************************************************************************************
 subroutine read_control_info(summaFileManagerIn, err, message)
 ! used to read information from model control file
 USE ascii_util_module,only:file_open       ! open file
 USE ascii_util_module,only:linewidth       ! max character number for one line
 USE ascii_util_module,only:get_vlines      ! get a vector of non-comment lines
 USE get_ixname_module,only:get_ixcontrol   ! identify index of named variable
 USE globalData,only:model_control_info     ! model control file info structure
 implicit none

 ! input
 character(*),intent(in)             ::  summaFileManagerIn  ! model control file

 ! output
 integer(i4b),intent(out)             :: err            ! error code
 character(*),intent(out)             :: message        ! error message

 ! local variables
 character(len=256)                   :: cmessage       ! error message for downwind routine
 character(LEN=2048)                  :: infile         ! input filename
 integer(i4b)                         :: unt            ! file unit (free unit output from file_open)
 character(LEN=linewidth),allocatable :: charline(:)    ! vector of character strings
 integer(i4b)                         :: nControl       ! number of model info
 integer(i4b)                         :: iControl       ! index of model info
 character(len=summaPathLen)          :: entry          ! name of model info
 character(len=32)                    :: option         ! option for model info
 integer(i4b)                         :: iVar           ! index of the info in the data structure
 
 ! Start routine here
 err=0; message='read_control_info/'
 
 ! open file, read non-comment lines, close file
 infile=trim(summaFileManagerIn)
 write(*,'(2(a,1x))') 'model control (file manager) file = ', infile)
 call file_open(trim(infile),unt,err,cmessage)
 if(err/=0) then; message=trim(message)//trim(cmessage)//"/Failed to open control file [''"//infile//"']"; err=-10; return; end if
 call get_vlines(unt,charline,err,cmessage)  ! 'charline' is a list of strings from non-comment lines
 if(err/=0) then; message=trim(message)//trim(cmessage)//"/Control file read issue in get_vlines()"; return; end if
 close(unt)
  
 ! get the number of model control file entries
 nControl = size(charline)

 ! populate the model control info structure
 do iControl=1,nControl
  ! extract name of decision and the decision selected
  read(charline(iControl),*,iostat=err) option, entry
  if (err/=0) then; err=30; message=trim(message)//"error reading charline array"; return; end if
  ! get the index of the control file entry in the data structure
  iVar = get_ixControl(trim(option))
  write(*,'(i4,1x,a)') iControl, trim(option)//': '//trim(entry)
  if(iVar<=0)then; err=40; message=trim(message)//"cannot recognize control file option [name='"//trim(option)//"']"; return; end if
  ! populate the model control info structure
  model_control_info(iVar)%cOption  = trim(option)
  model_control_info(iVar)%cControl = trim(entry)
 end do
 
 end subroutine read_model_info

END MODULE summafilemanager

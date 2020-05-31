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
! 2020 substantially rewritten by A. Wood to use keywords and read time info
!******************************************************************
MODULE summaFileManager
use nrtype
implicit none
public
! summa-wide pathlength
integer(i4b),parameter       :: summaPathLen=4096
! defines the time of the run
CHARACTER(LEN=summaPathLen)  :: CONTROL_VRS      = 'SUMMA_FILEMANAGER_V2.0'         ! control version
CHARACTER(LEN=summaPathLen)  :: SIM_START_TM     = '2000-01-01 00:00'               ! simulation start time
CHARACTER(LEN=summaPathLen)  :: SIM_END_TM       = '2000-01-01 00:00'               ! simulation end time
CHARACTER(LEN=summaPathLen)  :: NC_TIME_ZONE     = 'utcTime'                        ! time zone info
! defines the path for data files (and default values)
CHARACTER(LEN=summaPathLen)  :: SETNGS_PATH      = 'settings/'                      ! settings dir path
CHARACTER(LEN=summaPathLen)  :: INPUT_PATH       = 'input/default/'                 ! input_dir_path
CHARACTER(LEN=summaPathLen)  :: OUTPUT_PATH      = 'output/default/'                ! output_dir_path
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
 USE ascii_util_module,only:file_open       ! function to open file
 USE ascii_util_module,only:linewidth       ! max character number for one line
 USE ascii_util_module,only:get_vlines      ! function to get a vector of non-comment lines

 implicit none

 ! input/output vars
 character(*),intent(in)              :: summaFileManagerIn
 integer(i4b),intent(out)             :: err
 character(*),intent(out)             :: message
 ! local vars
 character(*),parameter               :: summaFileManagerHeader="SUMMA_FILE_MANAGER_V2.0"
 integer(i4b),parameter               :: runinfo_fileunit=67   ! file unit for run time information
 character(len=8)                     :: cdate
 character(len=10)                    :: ctime
 character(len=256)                   :: cmessage       ! error message for downwind routine
 integer(i4b)                         :: unt            ! file unit (free unit output from file_open)
 character(LEN=linewidth),allocatable :: charline(:)    ! vector of character strings
 integer(i4b)                         :: iControl, nControl       ! number of model info
 character(len=summaPathLen)          :: varEntry          ! name of model info
 character(len=32)                    :: option         ! option for model info

 ! Start procedure here
 err=0; message="summa_SetTimesDirsUndPhiles/"
 
 ! read information from model control file, and populate model control structure
 ! populates global control information structure
 
 ! open file, read non-comment lines, close file
 write(*,'(2(a,1x))') 'model control (file manager) file = ', trim(summaFileManagerIn)
 call file_open(trim(summaFileManagerIn),unt,err,cmessage)
 if(err/=0) then; message=trim(message)//trim(cmessage)//"/Failed to open control file [''"//trim(summaFileManagerIn)//"']"; err=-10; return; end if
 call get_vlines(unt,charline,err,cmessage)  ! 'charline' is a list of strings from non-comment lines
 if(err/=0) then; message=trim(message)//trim(cmessage)//"/Control file read issue in get_vlines()"; return; end if
 close(unt)
  
 ! get the number of model control file entries
 nControl = size(charline)
 print *, 'Read ',nControl,' entries in control file (fileManager)'

 ! populate the model control info structure
 do iControl=1,nControl
  ! extract name of decision and the decision selected
  read(charline(iControl),*,iostat=err) option, varEntry
  if (err/=0) then; err=30; message=trim(message)//"error reading charline array"; return; end if
  ! get the index of the control file entry in the data structure
  write(*,'(i4,1x,a)') iControl, trim(option)//': '//trim(varEntry)
 
  ! assign entries from control file to module public variables; add checking as needed
  select case(trim(option))
   case('controlVersion' ); CONTROL_VRS = trim(varEntry); if(trim(varEntry)/=summaFileManagerHeader)then;&
     message=trim(message)//"unknown control file version in '"//trim(summaFileManagerIn); err=20; return; end if
   case('simStartTime'   ); SIM_START_TM = trim(varEntry)                  ! start simulation time
   case('simEndTime'     ); SIM_END_TM = trim(varEntry)                    ! end simulation time
   case('tmZoneInfo'     ); NC_TIME_ZONE = trim(varEntry)                  ! time zone info
   case('settingsPath'   ); SETNGS_PATH = trim(varEntry)                   ! settings directory
   case('forcingPath'    ); INPUT_PATH = trim(varEntry)                    ! input forcing directory
   case('outputPath'     ); OUTPUT_PATH = trim(varEntry)                   ! output directory
   case('decisionsFile'  ); M_DECISIONS = trim(varEntry)                   ! model decisions file
   case('outputDefFile'  ); OUTPUT_CONTROL = trim(varEntry)                ! output control file
   case('hruParamFile'   ); LOCALPARAM_INFO = trim(varEntry)               ! default hru-level param file
   case('gruParamFile'   ); BASINPARAM_INFO = trim(varEntry)               ! default gru-level param file
   case('attributeFile'  ); LOCAL_ATTRIBUTES = trim(varEntry)              ! attribute file
   case('trialParamFile' ); PARAMETER_TRIAL = trim(varEntry)               ! trial parameters file
   case('forcingList'    ); FORCING_FILELIST = trim(varEntry)              ! file listing forcing filenames
   case('initCondFile'   ); MODEL_INITCOND = trim(varEntry)                ! initial conditions file (cold State)
   case('outFilePrefix'  ); OUTPUT_PREFIX = trim(varEntry)                 ! filename root for output files
   case('metaTime'       ); META_TIME = trim(varEntry)                     ! 
   case('metaAttr'       ); META_ATTR = trim(varEntry)                     !
   case('metaType'       ); META_TYPE = trim(varEntry)                     ! 
   case('metaForc'       ); META_FORCE = trim(varEntry)                    ! 
   case('metaLocParam'   ); META_LOCALPARAM = trim(varEntry)               ! 
   case('metaBasParam'   ); META_BASINPARAM = trim(varEntry)               ! 
   case('metaBasMvar'    ); META_BASINMVAR = trim(varEntry)                ! 
   case('metaLocIndex'   ); META_LOCALINDEX = trim(varEntry)               ! 
   ! get to here if cannot find the variable
   case default
     err=10; message=trim(message)//"unknown control file option: "//trim(option); return
  end select
 end do

 ! before embarking on a run, check that the output directory is writable; write system date and time to a log file there
 open(runinfo_fileunit,file=trim(OUTPUT_PATH)//"runinfo.txt",iostat=err)
 if(err/=0)then; err=10; message=trim(message)//"cannot write to output directory '"//trim(OUTPUT_PATH)//"'"; return; end if
 call date_and_time(cdate,ctime)
 write(runinfo_fileunit,*) 'Run start time on system:  ccyy='//cdate(1:4)//' - mm='//cdate(5:6)//' - dd='//cdate(7:8), &
                         ' - hh='//ctime(1:2)//' - mi='//ctime(3:4)//' - ss='//ctime(5:10)
 close(runinfo_fileunit)

 ! End procedure here
 end subroutine summa_SetTimesDirsUndPhiles


END MODULE summaFileManager

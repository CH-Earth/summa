!******************************************************************
! (C) Copyright 2009-2010  ---  Dmitri Kavetski and Martyn Clark ---  All rights reserved
!******************************************************************
MODULE snow_filemanager
use nrtype
implicit none
public
! FUSE-wide pathlength
integer(i4b),parameter::fusePathLen=256
! defines the path for data files (and default values)
CHARACTER(LEN=fusePathLen)  :: SETNGS_PATH='/d1/mclark/FUSE_SNOW/settings/'         ! dinoris
CHARACTER(LEN=fusePathLen)  :: INPUT_PATH ='/d1/mclark/FUSE_SNOW/input/default/'    ! dinoris
CHARACTER(LEN=fusePathLen)  :: OUTPUT_PATH='/d1/mclark/FUSE_SNOW/output/default/'   ! dinoris
! define name of control files    (and default values)
CHARACTER(LEN=fusePathLen)  :: M_DECISIONS    ='snow_zDecisions.txt'      ! definition of model decisions
CHARACTER(LEN=fusePathLen)  :: META_TIME      ='snow_zTimeMeta.txt'       ! metadata for time
CHARACTER(LEN=fusePathLen)  :: META_FORCE     ='snow_zForceMeta.txt'      ! metadata for model forcing variables
CHARACTER(LEN=fusePathLen)  :: META_PARAM     ='snow_zParamMeta.txt'      ! metadata for model parameters
CHARACTER(LEN=fusePathLen)  :: META_MVAR      ='snow_zModelVarMeta.txt'   ! metadata for model variables
CHARACTER(LEN=fusePathLen)  :: META_INDEX     ='snow_zModelIndexMeta.txt' ! metadata for model indices
CHARACTER(LEN=fusePathLen)  :: PARAMETER_INFO ='snow_zParamInfo.txt'      ! default values and constraints for model parameters
CHARACTER(LEN=fusePathLen)  :: FORCEFILE_DESC ='snow_zforcingInfo.txt'    ! description of forcing data file
CHARACTER(LEN=fusePathLen)  :: MODEL_INITCOND ='snow_zInitialCond.txt'    ! model initial conditions
CHARACTER(LEN=fusePathLen)  :: PARAMETER_TRIAL='snow_zParamTrial.txt'     ! trial values for model parameters
CHARACTER(LEN=fusePathLen)  :: OUTPUT_PREFIX  ='snow_zOutputPrefix.txt'   ! prefix for the output file
!----------------------------------------------------
contains
!----------------------------------------------------
subroutine fuse_SetDirsUndPhiles(fuseFileManagerIn,err,message)
! Purpose: Sets direcotries and philenames for FUSE.
! ---
! Programmer: Dmitri Kavetski and Martyn Clark
! Last modified: NCAR, 20110408
! ---
! Usage
! fuseFileManagerIn     = global names/path file
implicit none
! dummies
character(*),intent(in) ::fuseFileManagerIn
integer(i4b),intent(out)::err
character(*),intent(out)::message
! locals
logical(lgt)::xist
integer(i4b),parameter::unt=99 !DK: need to either define units globally, or use getSpareUnit
character(*),parameter::fuseFileManagerHeader="SNOW_FILEMANAGER_V1.0"
character(LEN=100)::temp
integer(i4b)::ierr ! temporary error code
integer(i4b),parameter :: runinfo_fileunit=67 ! file unit for run time information
character(len=8)  :: cdate
character(len=10) :: ctime

! Start procedure here
err=0; message="fuseSetDirsUndPhiles/"
! check if the file manager file exists
inquire(file=fuseFileManagerIn,exist=xist) ! Check for existence of masterfile
if(.not.xist)then
  message="f-fuseSetDirsUndPhiles/fuseFileManager/FileNotFound['"//trim(fuseFileManagerIn)//"']"&
              //'/ProceedingWithDefaults'
  err=-10; return
endif
! open file manager file
open(unt,file=fuseFileManagerIn,status="old",action="read",iostat=err)
if(err/=0)then
  message="f-fuseSetDirsUndPhiles/fileManagerOpenError['"//trim(fuseFileManagerIn)//"']"
  err=10; return
endif
! check the header matches the code
read(unt,*)temp
if(trim(temp)/=fuseFileManagerHeader)then
  message="f-fuseSetDirsUndPhiles/unknownHeader&[file='"//trim(fuseFileManagerIn)//"']&&
    &[header="//trim(temp)//"]"
  err=20; return
endif
! read information from file
ierr=0  ! initialize errors
read(unt,'(a)')temp
read(unt,'(a)')temp
read(unt,*)SETNGS_PATH    ; call checkLineRead(SETNGS_PATH,    err,message); if(err/=0)return
read(unt,*)INPUT_PATH     ; call checkLineRead(INPUT_PATH,     err,message); if(err/=0)return
read(unt,*)OUTPUT_PATH    ; call checkLineRead(OUTPUT_PATH,    err,message); if(err/=0)return
read(unt,'(a)')temp
read(unt,*)M_DECISIONS    ; call checkLineRead(M_DECISIONS,    err,message); if(err/=0)return
read(unt,*)META_TIME      ; call checkLineRead(META_TIME,      err,message); if(err/=0)return
read(unt,*)META_FORCE     ; call checkLineRead(META_FORCE,     err,message); if(err/=0)return
read(unt,*)META_PARAM     ; call checkLineRead(META_PARAM,     err,message); if(err/=0)return
read(unt,*)META_MVAR      ; call checkLineRead(META_MVAR,      err,message); if(err/=0)return
read(unt,*)META_INDEX     ; call checkLineRead(META_INDEX,     err,message); if(err/=0)return
read(unt,*)PARAMETER_INFO ; call checkLineRead(PARAMETER_INFO, err,message); if(err/=0)return
read(unt,*)FORCEFILE_DESC ; call checkLineRead(FORCEFILE_DESC, err,message); if(err/=0)return
read(unt,*)MODEL_INITCOND ; call checkLineRead(MODEL_INITCOND, err,message); if(err/=0)return
read(unt,*)PARAMETER_TRIAL; call checkLineRead(PARAMETER_TRIAL,err,message); if(err/=0)return
read(unt,*)OUTPUT_PREFIX  ; call checkLineRead(OUTPUT_PREFIX,  err,message); if(err/=0)return
close(unt)
! check that the output directory exists and write the date and time to a log file
open(runinfo_fileunit,file=trim(OUTPUT_PATH)//"runinfo.txt",iostat=err)
if(err/=0)then; err=10; message=trim(message)//"cannot write to directory '"//trim(OUTPUT_PATH)//"'"; return; endif
call date_and_time(cdate,ctime)
write(runinfo_fileunit,*) 'ccyy='//cdate(1:4)//' - mm='//cdate(5:6)//' - dd='//cdate(7:8), &
                         ' - hh='//ctime(1:2)//' - mi='//ctime(3:4)//' - ss='//ctime(5:10)
close(runinfo_fileunit)
! End procedure here
endsubroutine fuse_SetDirsUndPhiles
!----------------------------------------------------
! check if there is a space in the character string
subroutine checkLineRead(stringInput,err,message)
implicit none
character(*),intent(in)   :: stringInput
integer(i4b),intent(inout):: err
character(*),intent(inout):: message
if(index(trim(stringInput),' ')/=0) then
 err=30; message="f-fuseSetDirsUndPhiles/spaceInString[string="//trim(stringInput)//"]"
endif
end subroutine checkLineRead
!----------------------------------------------------
END MODULE snow_filemanager

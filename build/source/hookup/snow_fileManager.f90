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
CHARACTER(LEN=fusePathLen)  :: M_DECISIONS      ='snow_zDecisions.txt'           ! definition of model decisions
CHARACTER(LEN=fusePathLen)  :: META_TIME        ='snow_zTimeMeta.txt'            ! metadata for time
CHARACTER(LEN=fusePathLen)  :: META_ATTR        ='snow_zLocalAttributeMeta.txt'  ! metadata for local attributes
CHARACTER(LEN=fusePathLen)  :: META_TYPE        ='snow_zCatergoryMeta.txt'       ! metadata for local classification of veg, soil, etc.
CHARACTER(LEN=fusePathLen)  :: META_FORCE       ='snow_zForceMeta.txt'           ! metadata for model forcing variables
CHARACTER(LEN=fusePathLen)  :: META_LOCALPARAM  ='snow_zLocalParamMeta.txt'      ! metadata for model parameters
CHARACTER(LEN=fusePathLen)  :: META_LOCALMVAR   ='snow_zLocalModelVarMeta.txt'   ! metadata for model variables
CHARACTER(LEN=fusePathLen)  :: META_LOCALINDEX  ='snow_zLocalModelIndexMeta.txt' ! metadata for model indices
CHARACTER(LEN=fusePathLen)  :: META_BASINPARAM  ='snow_zBasinParamMeta.txt'      ! metadata for model parameters
CHARACTER(LEN=fusePathLen)  :: META_BASINMVAR   ='snow_zBasinModelVarMeta.txt'   ! metadata for model variables
CHARACTER(LEN=fusePathLen)  :: LOCAL_ATTRIBUTES ='snow_zLocalAttributes.txt'     ! local attributes
CHARACTER(LEN=fusePathLen)  :: LOCALPARAM_INFO  ='snow_zLocalParamInfo.txt'      ! default values and constraints for local model parameters
CHARACTER(LEN=fusePathLen)  :: BASINPARAM_INFO  ='snow_zBasinParamInfo.txt'      ! default values and constraints for basin model parameters
CHARACTER(LEN=fusePathLen)  :: FORCING_FILELIST ='snow_zForcingFileList.txt'     ! list of focing files for each HRU
CHARACTER(LEN=fusePathLen)  :: MODEL_INITCOND   ='snow_zInitialCond.txt'         ! model initial conditions
CHARACTER(LEN=fusePathLen)  :: PARAMETER_TRIAL  ='snow_zParamTrial.txt'          ! trial values for model parameters
CHARACTER(LEN=fusePathLen)  :: OUTPUT_PREFIX    ='xx'                            ! prefix for the output file
contains


! *************************************************************************************************
! public subroutine fuse_SetDirsUndPhiles: Sets directories and filenames for FUSE
! *************************************************************************************************
subroutine fuse_SetDirsUndPhiles(fuseFileManagerIn,err,message)
! Purpose: Sets directories and philenames for FUSE.
! ---
! Programmer: Dmitri Kavetski and Martyn Clark
! Last modified: Vienna, 14 April 2013
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
character(*),parameter::fuseFileManagerHeader="SNOW_FILEMANAGER_V1.4"
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
  message=trim(message)//"FileNotFound['"//trim(fuseFileManagerIn)//"']"&
                       //'/ProceedingWithDefaults'
  err=-10; return
endif
! open file manager file
open(unt,file=fuseFileManagerIn,status="old",action="read",iostat=err)
if(err/=0)then
  message=trim(message)//"fileManagerOpenError['"//trim(fuseFileManagerIn)//"']"
  err=10; return
endif
! check the header matches the code
read(unt,*)temp
if(trim(temp)/=fuseFileManagerHeader)then
  message=trim(message)//"unknownHeader&[file='"//trim(fuseFileManagerIn)//"']&&
    &[header="//trim(temp)//"]"
  err=20; return
endif
! read information from file
ierr=0  ! initialize errors
read(unt,'(a)')temp
read(unt,'(a)')temp
read(unt,*)SETNGS_PATH     ; call checkLineRead(SETNGS_PATH,      err,message); if(err/=0)return
read(unt,*)INPUT_PATH      ; call checkLineRead(INPUT_PATH,       err,message); if(err/=0)return
read(unt,*)OUTPUT_PATH     ; call checkLineRead(OUTPUT_PATH,      err,message); if(err/=0)return
read(unt,'(a)')temp
read(unt,*)M_DECISIONS     ; call checkLineRead(M_DECISIONS,      err,message); if(err/=0)return
read(unt,*)META_TIME       ; call checkLineRead(META_TIME,        err,message); if(err/=0)return
read(unt,*)META_ATTR       ; call checkLineRead(META_ATTR,        err,message); if(err/=0)return
read(unt,*)META_TYPE       ; call checkLineRead(META_TYPE,        err,message); if(err/=0)return
read(unt,*)META_FORCE      ; call checkLineRead(META_FORCE,       err,message); if(err/=0)return
read(unt,*)META_LOCALPARAM ; call checkLineRead(META_LOCALPARAM,  err,message); if(err/=0)return
read(unt,*)META_LOCALMVAR  ; call checkLineRead(META_LOCALMVAR,   err,message); if(err/=0)return
read(unt,*)META_LOCALINDEX ; call checkLineRead(META_LOCALINDEX,  err,message); if(err/=0)return
read(unt,*)META_BASINPARAM ; call checkLineRead(META_BASINPARAM,  err,message); if(err/=0)return
read(unt,*)META_BASINMVAR  ; call checkLineRead(META_BASINMVAR,   err,message); if(err/=0)return
read(unt,*)LOCAL_ATTRIBUTES; call checkLineRead(LOCAL_ATTRIBUTES, err,message); if(err/=0)return
read(unt,*)LOCALPARAM_INFO ; call checkLineRead(LOCALPARAM_INFO,  err,message); if(err/=0)return
read(unt,*)BASINPARAM_INFO ; call checkLineRead(BASINPARAM_INFO,  err,message); if(err/=0)return
read(unt,*)FORCING_FILELIST; call checkLineRead(FORCING_FILELIST, err,message); if(err/=0)return
read(unt,*)MODEL_INITCOND  ; call checkLineRead(MODEL_INITCOND,   err,message); if(err/=0)return
read(unt,*)PARAMETER_TRIAL ; call checkLineRead(PARAMETER_TRIAL,  err,message); if(err/=0)return
read(unt,*)OUTPUT_PREFIX   ; call checkLineRead(OUTPUT_PREFIX,    err,message); if(err/=0)return
close(unt)
! check that the output directory exists and write the date and time to a log file
open(runinfo_fileunit,file=trim(OUTPUT_PATH)//"runinfo.txt",iostat=err)
if(err/=0)then; err=10; message=trim(message)//"cannot write to directory '"//trim(OUTPUT_PATH)//"'"; return; endif
call date_and_time(cdate,ctime)
write(runinfo_fileunit,*) 'ccyy='//cdate(1:4)//' - mm='//cdate(5:6)//' - dd='//cdate(7:8), &
                         ' - hh='//ctime(1:2)//' - mi='//ctime(3:4)//' - ss='//ctime(5:10)
close(runinfo_fileunit)
! End procedure here
end subroutine fuse_SetDirsUndPhiles


! *************************************************************************************************
! public subroutine checkLineRead: check if there is a space in the character string
! *************************************************************************************************
subroutine checkLineRead(stringInput,err,message)
implicit none
character(*),intent(in)   :: stringInput
integer(i4b),intent(inout):: err
character(*),intent(inout):: message
if(index(trim(stringInput),' ')/=0) then
 err=30; message="f-fuseSetDirsUndPhiles/spaceInString[string="//trim(stringInput)//"]"
endif
end subroutine checkLineRead


END MODULE snow_filemanager

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
MODULE summafilemanager
use nrtype
implicit none
public
! summa-wide pathlength
integer(i4b),parameter::summaPathLen=256
! defines the path for data files (and default values)
CHARACTER(LEN=summaPathLen)  :: SETNGS_PATH='settings/'         ! SETNGS_PATH
CHARACTER(LEN=summaPathLen)  :: INPUT_PATH ='input/default/'    ! INPUT_PATH
CHARACTER(LEN=summaPathLen)  :: OUTPUT_PATH='output/default/'   ! OUTPUT_PATH
! define name of control files    (and default values)
CHARACTER(LEN=summaPathLen)  :: M_DECISIONS      ='summa_zDecisions.txt'           ! definition of model decisions
CHARACTER(LEN=summaPathLen)  :: META_TIME        ='summa_zTimeMeta.txt'            ! metadata for time
CHARACTER(LEN=summaPathLen)  :: META_ATTR        ='summa_zLocalAttributeMeta.txt'  ! metadata for local attributes
CHARACTER(LEN=summaPathLen)  :: META_TYPE        ='summa_zCatergoryMeta.txt'       ! metadata for local classification of veg, soil, etc.
CHARACTER(LEN=summaPathLen)  :: META_FORCE       ='summa_zForceMeta.txt'           ! metadata for model forcing variables
CHARACTER(LEN=summaPathLen)  :: META_LOCALPARAM  ='summa_zLocalParamMeta.txt'      ! metadata for model parameters
CHARACTER(LEN=summaPathLen)  :: META_LOCALMVAR   ='summa_zLocalModelVarMeta.txt'   ! metadata for model variables
CHARACTER(LEN=summaPathLen)  :: META_LOCALINDEX  ='summa_zLocalModelIndexMeta.txt' ! metadata for model indices
CHARACTER(LEN=summaPathLen)  :: META_BASINPARAM  ='summa_zBasinParamMeta.txt'      ! metadata for model parameters
CHARACTER(LEN=summaPathLen)  :: META_BASINMVAR   ='summa_zBasinModelVarMeta.txt'   ! metadata for model variables
CHARACTER(LEN=summaPathLen)  :: LOCAL_ATTRIBUTES ='summa_zLocalAttributes.txt'     ! local attributes
CHARACTER(LEN=summaPathLen)  :: LOCALPARAM_INFO  ='summa_zLocalParamInfo.txt'      ! default values and constraints for local model parameters
CHARACTER(LEN=summaPathLen)  :: BASINPARAM_INFO  ='summa_zBasinParamInfo.txt'      ! default values and constraints for basin model parameters
CHARACTER(LEN=summaPathLen)  :: FORCING_FILELIST ='summa_zForcingFileList.txt'     ! list of focing files for each HRU
CHARACTER(LEN=summaPathLen)  :: MODEL_INITCOND   ='summa_zInitialCond.txt'         ! model initial conditions
CHARACTER(LEN=summaPathLen)  :: PARAMETER_TRIAL  ='summa_zParamTrial.txt'          ! trial values for model parameters
CHARACTER(LEN=summaPathLen)  :: OUTPUT_PREFIX    ='summa_output_'                  ! prefix for the output file
contains


! *************************************************************************************************
! public subroutine summa_SetDirsUndPhiles: Sets directories and filenames for summa
! *************************************************************************************************
subroutine summa_SetDirsUndPhiles(summaFileManagerIn,err,message)
! Purpose: Sets directories and philenames for summa.
! ---
! Programmer: Dmitri Kavetski and Martyn Clark
! Last modified: Vienna, 14 April 2013
! ---
! Usage
! summaFileManagerIn     = global names/path file
implicit none
! dummies
character(*),intent(in) ::summaFileManagerIn
integer(i4b),intent(out)::err
character(*),intent(out)::message
! locals
logical(lgt)::xist
integer(i4b),parameter::unt=99 !DK: need to either define units globally, or use getSpareUnit
character(*),parameter::summaFileManagerHeader="SUMMA_FILE_MANAGER_V1.0"
character(LEN=100)::temp
integer(i4b)::ierr ! temporary error code
integer(i4b),parameter :: runinfo_fileunit=67 ! file unit for run time information
character(len=8)  :: cdate
character(len=10) :: ctime

! Start procedure here
err=0; message="summaSetDirsUndPhiles/"
! check if the file manager file exists
inquire(file=summaFileManagerIn,exist=xist) ! Check for existence of masterfile
if(.not.xist)then
  message=trim(message)//"FileNotFound['"//trim(summaFileManagerIn)//"']"&
                       //'/ProceedingWithDefaults'
  err=-10; return
endif
! open file manager file
open(unt,file=summaFileManagerIn,status="old",action="read",iostat=err)
if(err/=0)then
  message=trim(message)//"fileManagerOpenError['"//trim(summaFileManagerIn)//"']"
  err=10; return
endif
! check the header matches the code
read(unt,*)temp
if(trim(temp)/=summaFileManagerHeader)then
  message=trim(message)//"unknownHeader&[file='"//trim(summaFileManagerIn)//"']&&
    &[header="//trim(temp)//"]"
  err=20; return
endif
! read information from file
ierr=0 ! initialize errors

call readLine(unt,SETNGS_PATH,      err,message)
call readLine(unt,INPUT_PATH,       err,message)
call readLine(unt,OUTPUT_PATH,      err,message)

call readLine(unt,M_DECISIONS ,     err,message)
call readLine(unt,META_TIME ,       err,message)
call readLine(unt,META_ATTR,        err,message)
call readLine(unt,META_TYPE,        err,message)
call readLine(unt,META_FORCE,       err,message)
call readLine(unt,META_LOCALPARAM,  err,message)
call readLine(unt,META_LOCALMVAR ,  err,message)
call readLine(unt,META_LOCALINDEX,  err,message)
call readLine(unt,META_BASINPARAM,  err,message)
call readLine(unt,META_BASINMVAR,   err,message)
call readLine(unt,LOCAL_ATTRIBUTES, err,message)
call readLine(unt,LOCALPARAM_INFO,  err,message)
call readLine(unt,BASINPARAM_INFO,  err,message)
call readLine(unt,FORCING_FILELIST, err,message)
call readLine(unt,MODEL_INITCOND,   err,message)
call readLine(unt,PARAMETER_TRIAL,  err,message)
call readLine(unt,OUTPUT_PREFIX,    err,message)
close(unt)
! check that the output directory exists and write the date and time to a log file
open(runinfo_fileunit,file=trim(OUTPUT_PATH)//"runinfo.txt",iostat=err)
if(err/=0)then; err=10; message=trim(message)//"cannot write to directory '"//trim(OUTPUT_PATH)//"'"; return; endif
call date_and_time(cdate,ctime)
write(runinfo_fileunit,*) 'ccyy='//cdate(1:4)//' - mm='//cdate(5:6)//' - dd='//cdate(7:8), &
                         ' - hh='//ctime(1:2)//' - mi='//ctime(3:4)//' - ss='//ctime(5:10)
close(runinfo_fileunit)
! End procedure here
end subroutine summa_SetDirsUndPhiles



! *************************************************************************************************
! public subroutine readLine: read the first string to a string variable
! *************************************************************************************************
subroutine readLine(funt,str,err,message)
implicit none
integer(i4b),intent(in)   :: funt
character(*),intent(inout):: str
integer(i4b),intent(inout):: err
character(*),intent(inout):: message




do
 ! read line that is not comment
 read(funt,*) str
 if (str(1:1) /= '!') exit
end do

! check if there is a space in the character string
if(index(trim(str),' ')/=0) then
 err=30; message="f-summaSetDirsUndPhiles/spaceInString[string="//trim(str)//"]"
endif
end subroutine readLine



END MODULE summafilemanager


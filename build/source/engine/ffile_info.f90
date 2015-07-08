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

module ffile_info_module
USE nrtype
implicit none
private
public::ffile_info
contains


 ! ************************************************************************************************
 ! public subroutine ffile_info: read information on model forcing files
 ! ************************************************************************************************
 subroutine ffile_info(nHRU,err,message)
 ! used to read metadata on the forcing data file
 USE ascii_util_module,only:file_open
 USE summaFileManager,only:SETNGS_PATH       ! path for metadata files
 USE summaFileManager,only:FORCING_FILELIST  ! list of model forcing files
 USE data_struc,only:time_meta,forc_meta     ! model forcing metadata
 USE data_struc,only:forcFileInfo,data_step  ! info on model forcing file
 USE data_struc,only:type_hru                ! data structure for categorical data
 USE var_lookup,only:iLookTYPE               ! named variables to index elements of the data vectors
 USE get_ixname_module,only:get_ixtime,get_ixforce  ! identify index of named variable
 USE ascii_util_module,only:get_vlines      ! get a vector of non-comment lines
 USE ascii_util_module,only:split_line      ! split a line into words
 implicit none
 ! define output
 integer(i4b),intent(in)              :: nHRU           ! number of hydrologic response units
 integer(i4b),intent(out)             :: err            ! error code
 character(*),intent(out)             :: message        ! error message
 ! define local variables
 character(LEN=1024),allocatable      :: dataLines(:)   ! vector of lines of information (non-comment lines)
 integer(i4b),parameter               :: imiss = -999   ! missing data
 character(len=256)                   :: cmessage       ! error message for downwind routine
 character(LEN=256)                   :: infile         ! input filename
 integer(i4b),parameter               :: unt=99         ! DK: need to either define units globally, or use getSpareUnit
 integer(i4b)                         :: iline          ! loop through lines in the file
 integer(i4b),parameter               :: maxLines=1000  ! maximum lines in the file
 character(LEN=256)                   :: filenameDesc   ! name of file that describes the forcing datafile
 character(LEN=256)                   :: temp='uninitialized'  ! single lime of information
 integer(i4b)                         :: iend           ! check for the end of the file
 character(LEN=256)                   :: ffmt           ! file format
 character(LEN=32)                    :: varname        ! name of variable
 character(LEN=64)                    :: vardata        ! data on variable
 character(len=2)                     :: dLim           ! column delimiter
 integer(i4b)                         :: ivar           ! index of model variable
 integer(i4b)                         :: iHRU,jHRU,kHRU ! index of HRUs (position in vector)
 integer(i4b)                         :: hruIndex       ! identifier of each HRU
 real(dp)                             :: dataStep_iHRU  ! data step for a given forcing data file
 ! Start procedure here
 err=0; message="ffile_info/"
 ! ------------------------------------------------------------------------------------------------------------------
 ! (1) read in the list of forcing files
 ! ------------------------------------------------------------------------------------------------------------------
 ! allocate space for forcing information
 if(associated(forcFileInfo)) deallocate(forcFileInfo)
 allocate(forcFileInfo(nHRU), stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for forcFileInfo'; return; endif
 ! build filename
 infile = trim(SETNGS_PATH)//trim(FORCING_FILELIST)
 ! open file
 call file_open(trim(infile),unt,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! get a list of character strings from non-comment lines
 call get_vlines(unt,dataLines,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
 ! check that we have the correct number of HRUs
 if(size(dataLines) /= nHRU)then; err=20; message=trim(message)//'incorrect number of HRUs in file ['//trim(infile)//']'; return; endif
 ! loop through list of forcing descriptor files and put in the appropriate place in the data structure
 do iHRU=1,nHRU
  ! split the line into "words" (expect two words: the HRU index, and the file describing forcing data for that index)
  read(dataLines(iHRU),*,iostat=err) hruIndex, filenameDesc
  if(err/=0)then; message=trim(message)//'problem reading a line of data from file ['//trim(infile)//']'; return; endif
  ! identify the HRU index
  do jHRU=1,nHRU
   if(hruIndex == type_hru(jHRU)%var(iLookTYPE%hruIndex))then
    kHRU=jHRU
    exit
   endif
   if(jHRU == nHRU)then ! we get to here if we have tested the last HRU and have not exited the loop
    write(message,'(a,i0,a)') trim(message)//'unable to identify HRU in forcing file description [index = ',hruIndex,'; file='//trim(infile)//']'
    err=20; return
   endif
  end do
  ! put the filename in the structure
  forcFileInfo(kHRU)%filenmDesc = trim(filenameDesc)
  write(*,'(2(a,1x),2(i6,1x))') 'filenameDesc, hruIndex, kHRU = ', trim(filenameDesc), hruIndex, kHRU
 end do  ! (looping through files)
 close(unt)
 ! ------------------------------------------------------------------------------------------------------------------
 ! (2) read in the information that describes each forcing file
 ! ------------------------------------------------------------------------------------------------------------------
 ! check that the time metadata is already populated
 if(.not.associated(time_meta))then; err=30; message=trim(message)//"TimeMetadataNonexistent"; return; endif
 ! check that the forcing metadata is already populated
 if(.not.associated(forc_meta))then; err=30; message=trim(message)//"ForcingMetadataNonexistent"; return; endif
 ! read description of file that is used in each HRU
 do iHRU=1,nHRU
  ! allocate space for the column indices
  if(associated(forcFileInfo(iHRU)%time_ix)) deallocate(forcFileInfo(iHRU)%time_ix)
  if(associated(forcFileInfo(iHRU)%data_ix)) deallocate(forcFileInfo(iHRU)%data_ix)
  allocate(forcFileInfo(iHRU)%time_ix(size(time_meta)),&
           forcFileInfo(iHRU)%data_ix(size(forc_meta)),stat=err)
  if(err/=0)then; err=40; message=trim(message)//"problemAllocateStructureElement"; return; endif
  ! initialize column indices to missing
  forcFileInfo(iHRU)%time_ix(:) = imiss
  forcFileInfo(iHRU)%data_ix(:) = imiss
  ! build filename
  infile = trim(SETNGS_PATH)//trim(forcFileInfo(iHRU)%filenmDesc)
  !print*, 'infile = ', trim(infile)
  ! open file
  call file_open(trim(infile),unt,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! get to the start of the variable descriptions
  do iline=1,maxLines
   read(unt,'(a)',iostat=iend)temp; if (iend/=0)exit    ! read line of data
   if (temp(1:1)/='!') exit  ! assume first line not comment is format code
  end do ! looping through file to find the format code
  ! read in format string
  read(temp,*)ffmt
  ! loop through the lines in the file
  do iline=1,maxLines
   ! read a line of data and exit if an error code (character read, so only possible error is end of file)
   read(unt,'(a)',iostat=iend)temp; if (iend/=0)exit
   ! check that the line is not a comment
   if (temp(1:1)=='!')cycle
   ! save data into a temporary variables
   read(temp,trim(ffmt),iostat=err) varname, dLim, vardata
   if (err/=0) then; err=30; message=trim(message)//"errorReadLine[file="//trim(infile)//"; line="//trim(temp)//"]"; return; endif
   ! check the delimiter
   if(dLim(1:1)/='|')then; err=30; message=trim(message)//"incorrectFormat"//trim(infile); return; endif
   !print*, 'varname = ', trim(varname)
   !print*, 'vardata = ', trim(vardata)
   ! put data into data structure
   select case(trim(varname))
    case('filenmData'); read(vardata,*) forcFileInfo(iHRU)%filenmData
    case('ncols'     ); read(vardata,*) forcFileInfo(iHRU)%ncols
    ! process the data step
    case('data_step' )
     read(vardata,*) dataStep_iHRU
     if(iHRU == 1)then
      data_step = dataStep_iHRU
     else
      if(abs(dataStep_iHRU - data_step) > epsilon(dataStep_iHRU))then
       write(message,'(a,i0,a)') trim(message)//'data step for HRU ',iHRU,'differs from the datastep of the first HRU'
       err=20; return
      endif
     endif
    ! ***** identify the index of the time data variable
    case('iyyy','im','id','ih','imin')
     ivar = get_ixtime(trim(varname))
     if(ivar < 0)then; err=40; message=trim(message)//"variableNotFound[var="//trim(varname)//"]"; return; endif
     if(ivar>size(forcFileInfo(iHRU)%time_ix))then
      err=35; message=trim(message)//"indexOutOfRange[var="//trim(varname)//"]"; return
     endif
     ! put column index in the structure
     read(vardata,*) forcFileInfo(iHRU)%time_ix(ivar)
    ! ***** identity index for the forcing data variable
    case('pptrate','SWRadAtm','LWRadAtm','airtemp','windspd','airpres','spechum')
     ivar = get_ixforce(trim(varname))
     if(ivar < 0)then; err=40; message=trim(message)//"variableNotFound[var="//trim(varname)//"]"; return; endif
     if(ivar>size(forcFileInfo(iHRU)%data_ix))then
      err=35; message=trim(message)//"indexOutOfRange[var="//trim(varname)//"]"; return
     endif
     ! put column index in the structure
     read(vardata,*) forcFileInfo(iHRU)%data_ix(ivar)
    ! ***** error check
    case default
     message=trim(message)//'variableNotFound[var='//trim(varname)//'; file='//trim(infile)//']'
     err=20; return
   endselect
  enddo ! (loop through lines in the file)
  ! close file unit
  close(unt)
 end do  ! (looping through files describing each HRU)
 ! identify the first HRU to use a given data file
 do iHRU=1,nHRU
  do jHRU=1,iHRU-1
   if(trim(forcFileInfo(iHRU)%filenmData) == trim(forcFileInfo(jHRU)%filenmData))then
    forcFileInfo(iHRU)%ixFirstHRU = jHRU  ! index of first HRU to share the same data
   else
    forcFileInfo(iHRU)%ixFirstHRU = 0
   endif
  end do
 end do
 end subroutine ffile_info


end module ffile_info_module

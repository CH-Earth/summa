module ffile_info_module
USE nrtype
implicit none
private
public::ffile_info
contains

 ! ************************************************************************************************
 ! new subroutine: read information on model forcing fils
 ! ************************************************************************************************
 subroutine ffile_info(err,message)
 ! used to read metadata on the forcing data file
 USE ascii_util_module,only:file_open
 USE snow_fileManager,only:SETNGS_PATH     ! path for metadata files
 USE snow_fileManager,only:FORCEFILE_DESC  ! description of model forcing datafile
 USE data_struc,only:time_meta,forc_meta   ! model forcing metadata
 USE data_struc,only:forcFileInfo          ! info on model forcing file
 USE get_ixname_module,only:get_ixtime,get_ixforce    ! identify index of named variable
 implicit none
 ! define output
 integer(i4b),intent(out)             :: err            ! error code
 character(*),intent(out)             :: message        ! error message
 ! define local variables
 integer(i4b),parameter               :: imiss = -999   ! missing data
 character(len=256)                   :: cmessage       ! error message for downwind routine
 character(LEN=256)                   :: infile         ! input filename
 integer(i4b),parameter               :: unt=99         ! DK: need to either define units globally, or use getSpareUnit
 integer(i4b)                         :: iline          ! loop through lines in the file 
 integer(i4b),parameter               :: maxLines=1000  ! maximum lines in the file 
 character(LEN=256)                   :: temp           ! single lime of information
 integer(i4b)                         :: iend           ! check for the end of the file
 character(LEN=256)                   :: ffmt           ! file format
 character(LEN=32)                    :: varname        ! name of variable
 character(LEN=32)                    :: vardata        ! data on variable
 character(len=2)                     :: dLim           ! column delimiter
 integer(i4b)                         :: ivar           ! index of model variable
 ! Start procedure here
 err=0; message="f-fuse/ffile_info/"
 ! build filename
 infile = trim(SETNGS_PATH)//trim(FORCEFILE_DESC)
 ! open file
 call file_open(trim(infile),unt,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! check that the time metadata is already populated
 if(.not.associated(time_meta))then; err=30; message=trim(message)//"TimeMetadataNonexistent"; return; endif
 ! check that the forcing metadata is already populated
 if(.not.associated(forc_meta))then; err=30; message=trim(message)//"ForcingMetadataNonexistent"; return; endif
 ! allocate space for the forcing file structure
 if (associated(forcFileInfo)) deallocate(forcFileInfo)
 allocate(forcFileInfo,stat=err)
 if(err/=0)then; err=40; message=trim(message)//"problemAllocateStructure"; return; endif
 ! allocate space for the column indices
 if(associated(forcFileInfo%time_ix)) deallocate(forcFileInfo%time_ix)
 if(associated(forcFileInfo%data_ix)) deallocate(forcFileInfo%data_ix)
 allocate(forcFileInfo%time_ix(size(time_meta)),&
          forcFileInfo%data_ix(size(forc_meta)),stat=err)
 if(err/=0)then; err=40; message=trim(message)//"problemAllocateStructureElement"; return; endif
 ! initialize column indices to missing
 forcFileInfo%time_ix(:) = imiss
 forcFileInfo%data_ix(:) = imiss
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
  if (err/=0) then; err=30; message=trim(message)//"errorReadLine"; return; endif
  ! put data into data structure
  select case(trim(varname))
   case('filenm')   ; read(vardata,*) forcFileInfo%filenm
   case('ncols')    ; read(vardata,*) forcFileInfo%ncols
   case('data_step'); read(vardata,*) forcFileInfo%data_step
   case('istart')   ; read(vardata,*) forcFileInfo%istart
   case('numtim')   ; read(vardata,*) forcFileInfo%numtim
   case default
    ! ***** identity index for the forcing data variable
    ivar = get_ixforce(trim(varname))
    if(ivar>0)then
     if(ivar>size(forcFileInfo%data_ix))then
      err=35; message=trim(message)//"indexOutOfRange[var="//trim(varname)//"]"; return
     endif
     ! put column index in the structure
     read(vardata,*) forcFileInfo%data_ix(ivar)
    else
     ! ***** identify the index of the time data variable
     ivar = get_ixtime(trim(varname))
     if(ivar>0)then
      if(ivar>size(forcFileInfo%time_ix))then
       err=35; message=trim(message)//"indexOutOfRange[var="//trim(varname)//"]"; return
      endif
      ! put column index in the structure
      read(vardata,*) forcFileInfo%time_ix(ivar)
     else
      err=40; message="f-fuse/ffile_info/variableNotFound[var="//trim(varname)//"]"; return
     endif
    endif
  endselect
 enddo ! (loop through lines in the file)
 ! close file unit
 close(unt)
 end subroutine ffile_info

end module ffile_info_module

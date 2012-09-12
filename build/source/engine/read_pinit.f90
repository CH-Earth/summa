module read_pinit_module
USE nrtype
implicit none
private
public::read_pinit
contains

 ! ************************************************************************************************
 ! (1) new subroutine: read information on model forcing fils
 ! ************************************************************************************************
 subroutine read_pinit(err,message)
 ! used to read metadata on the forcing data file
 USE ascii_util_module,only:file_open      ! open ascii file
 USE snow_fileManager,only:SETNGS_PATH     ! path for metadata files
 USE snow_fileManager,only:PARAMETER_INFO  ! default values and constraints for model parameters
 USE data_struc,only:mpar_meta             ! parameter metadata
 USE data_struc,only:par_info,parFallback  ! default values and constraints for model parameters
 USE get_ixname_module,only:get_ixparam    ! identify index of named variable
 implicit none
 ! define output
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! define local variables
 real(dp),parameter                   :: amiss=1.d+30   ! missing data
 character(len=256)                   :: cmessage       ! error message for downwind routine
 character(LEN=256)                   :: infile         ! input filename
 integer(i4b),parameter               :: unt=99         ! DK: need to either define units globally, or use getSpareUnit
 integer(i4b)                         :: iline          ! loop through lines in the file 
 integer(i4b),parameter               :: maxLines=1000  ! maximum lines in the file 
 character(LEN=256)                   :: temp           ! single lime of information
 integer(i4b)                         :: iend           ! check for the end of the file
 character(LEN=256)                   :: ffmt           ! file format
 character(LEN=32)                    :: varname        ! name of variable
 type(par_info)                       :: parTemp        ! temporary parameter structure
 character(len=2)                     :: dLim           ! column delimiter
 integer(i4b)                         :: ivar           ! index of model variable
 ! Start procedure here
 err=0; message="f-fuse/read_pinit/"
 ! build filename
 infile = trim(SETNGS_PATH)//trim(PARAMETER_INFO)
 ! open file
 call file_open(trim(infile),unt,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! check that the parameter metadata is already populated
 if(.not.associated(mpar_meta))then; err=30; message=trim(message)//"ParameterMetadataNonexistent"; return; endif
 ! allocate space for the parameter structure
 if (associated(parFallback)) deallocate(parFallback)
 allocate(parFallback(size(mpar_meta)),stat=err)
 if(err/=0)then; err=40; message=trim(message)//"problemAllocateStructure"; return; endif
 ! fill parameter vector with missing data
 parFallback(:)%default_val = amiss
 parFallback(:)%lower_limit = amiss
 parFallback(:)%upper_limit = amiss
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
  read(temp,trim(ffmt),iostat=err) varname, dLim, parTemp%default_val, dLim, parTemp%lower_limit, dLim, parTemp%upper_limit
  if (err/=0) then; err=30; message=trim(message)//"errorReadLine"; return; endif
  ! identify the index of the variable in the data structure
  ivar = get_ixparam(trim(varname))
  if(ivar>0)then
   if(ivar>size(parFallback))then
    err=35; message=trim(message)//"indexOutOfRange[var="//trim(varname)//"]"; return
   endif
   ! put data in the structure
   parFallback(ivar)=parTemp
  else
   err=40; message=trim(message)//"variableNotFound[var="//trim(varname)//"]"; return
  endif
 end do  ! looping through variables in the file
 ! check we have populated all variables
 if(any(parFallback(:)%default_val > 0.99_dp*amiss))then
  do ivar=1,size(parFallback)
   if(parFallback(ivar)%default_val > 0.99_dp*amiss)then
    err=40; message=trim(message)//"variableNonexistent[var="//trim(mpar_meta(ivar)%varname)//"]"; return
   endif
  end do
 endif
 ! close file unit
 close(unt)
 end subroutine read_pinit

end module read_pinit_module

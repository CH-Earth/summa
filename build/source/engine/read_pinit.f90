module read_pinit_module
USE nrtype
implicit none
private
public::read_pinit
contains

 ! ************************************************************************************************
 ! (1) new subroutine: read default model parameter values and constraints
 ! ************************************************************************************************
 subroutine read_pinit(err,message)
 ! used to read metadata on the forcing data file
 USE snow_fileManager,only:SETNGS_PATH     ! path for metadata files
 USE snow_fileManager,only:PARAMETER_INFO  ! file containing site characteristix and fefault values and constraints for model parameters
 USE ascii_util_module,only:file_open      ! open ascii file
 USE ascii_util_module,only:split_line     ! extract the list of variable names from the character string
 USE data_struc,only:mpar_meta             ! parameter metadata
 USE data_struc,only:par_info,parFallback  ! data structures for default values and constraints for model parameters
 USE get_ixname_module,only:get_ixParam    ! identify index of named variable for model parameters
 implicit none
 ! define output
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! define general variables
 real(dp),parameter                   :: amiss=1.d+30   ! missing data
 character(len=256)                   :: cmessage       ! error message for downwind routine
 character(LEN=256)                   :: infile         ! input filename
 integer(i4b),parameter               :: unt=99         ! DK: need to either define units globally, or use getSpareUnit
 integer(i4b)                         :: iline          ! loop through lines in the file 
 integer(i4b),parameter               :: maxLines=1000  ! maximum lines in the file 
 character(LEN=256)                   :: temp           ! single lime of information
 ! define local variables for the default model parameters
 integer(i4b)                         :: iend           ! check for the end of the file
 character(LEN=256)                   :: ffmt           ! file format
 character(LEN=32)                    :: varname        ! name of variable
 type(par_info)                       :: parTemp        ! temporary parameter structure
 character(len=2)                     :: dLim           ! column delimiter
 integer(i4b)                         :: ivar           ! index of model variable
 ! Start procedure here
 err=0; message="f-fuse/read_pinit/"
 ! **********************************************************************************************
 ! (1) open files, etc.
 ! **********************************************************************************************
 ! build filename
 infile = trim(SETNGS_PATH)//trim(PARAMETER_INFO)
 ! open file
 call file_open(trim(infile),unt,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! **********************************************************************************************
 ! (2) read default model parameter values and constraints
 ! **********************************************************************************************
 ! check that the parameter metadata is already populated
 if(.not.associated(mpar_meta))then; err=30; message=trim(message)//"Parameter metadata is non-existent"; return; endif
 ! allocate space for the parameter structure
 if (associated(parFallback)) deallocate(parFallback)
 allocate(parFallback(size(mpar_meta)),stat=err)
 if(err/=0)then; err=40; message=trim(message)//"problemAllocateStructure"; return; endif
 ! fill parameter vector with missing data
 parFallback(:)%default_val = amiss
 parFallback(:)%lower_limit = amiss
 parFallback(:)%upper_limit = amiss
 ! ---------------------------------------------------------------------------------------------
 ! read format code
 ! ---------------------------------------------------------------------------------------------
 do iline=1,maxLines
  ! (read through comment lines)
  read(unt,'(a)',iostat=iend) temp  ! read a line of data
  if(iend/=0)then; err=20; message=trim(message)//'got to end of file before found the format code'; return; endif
  if (temp(1:1)=='!')cycle 
  ! (read in format string -- assume that the first non-comment line is the format code)
  read(temp,*)ffmt  ! read in format string
  exit
  if(iLine==maxLines)then; err=20; message=trim(message)//'problem finding format code -- no non-comment line after start of parameter definitions'; return; endif 
 end do ! looping through lines
 ! ---------------------------------------------------------------------------------------------
 ! read in default values of model parameters, and parameter constraints
 ! ---------------------------------------------------------------------------------------------
 do iline=1,maxLines
  ! (read through comment lines)
  read(unt,'(a)',iostat=iend) temp  ! read a line of data
  if(iend/=0)exit !end of file
  if (temp(1:1)=='!')cycle
  ! (save data into a temporary variables)
  read(temp,trim(ffmt),iostat=err) varname, dLim, parTemp%default_val, dLim, parTemp%lower_limit, dLim, parTemp%upper_limit
  if (err/=0) then; err=30; message=trim(message)//"errorReadLine"; return; endif
  ! (identify the index of the variable in the data structure)
  ivar = get_ixParam(trim(varname))
  if(ivar>0)then
   if(ivar>size(parFallback))then
    err=35; message=trim(message)//"indexOutOfRange[var="//trim(varname)//"]"; return
   endif
   ! (put data in the structure)
   parFallback(ivar)=parTemp
  else
   err=40; message=trim(message)//"variableNotFound[var="//trim(varname)//"]"; return
  endif
 end do  ! (looping through lines in the file)
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

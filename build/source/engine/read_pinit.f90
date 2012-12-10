module read_pinit_module
USE nrtype
implicit none
private
public::read_pinit
contains

 ! ************************************************************************************************
 ! (1) new subroutine: read information on model site characteristix and parameter info
 ! ************************************************************************************************
 subroutine read_pinit(err,message)
 ! used to read metadata on the forcing data file
 USE snow_fileManager,only:SETNGS_PATH     ! path for metadata files
 USE snow_fileManager,only:PARAMETER_INFO  ! file containing site characteristix and fefault values and constraints for model parameters
 USE ascii_util_module,only:file_open      ! open ascii file
 USE ascii_util_module,only:split_line     ! extract the list of variable names from the character string
 USE allocspace_module,only:alloc_site     ! module to allocate space for the site characteristix
 USE data_struc,only:site_meta,mpar_meta   ! site characteristix and parameter metadata
 USE data_struc,only:site_data             ! data structures for site characteristix
 USE data_struc,only:par_info,parFallback  ! data structures for default values and constraints for model parameters
 USE get_ixname_module,only:get_ixSite     ! identify index of named variable for site characteristix
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
 ! define local variables for site characteristics
 character(len=256),parameter         :: site_tag='siteCharacteristix'
 logical(lgt)                         :: site_flag      ! determines if in the characteristix part of the file
 integer(i4b)                         :: iWord          ! index of a "word" in a character string
 integer(i4b)                         :: nElements      ! number of elements for a given variable
 character(LEN=256),allocatable       :: chardata(:)    ! vector of character data
 ! define local variables for the default model parameters
 character(len=256),parameter         :: param_tag='modelParameters'
 logical(lgt)                         :: param_flag     ! determines if in the parameter part of the file
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
 ! (2) read site characteristics
 ! **********************************************************************************************
 ! allocate space for the site characteristix structure
 call alloc_site(err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! initialize flag to identify the site characteristics part of the file
 site_flag=.false.
 ! loop through file until reach the site characteristics tag
 do iline=1,maxLines
  ! ---------------------------------------------------------------------------------------------
  ! identify site characteristics in the file, and read a given characteristic
  ! ---------------------------------------------------------------------------------------------
  ! (read through comment lines)
  read(unt,'(a)',iostat=iend)temp; if(iend/=0)then; rewind(unt); exit; endif    ! read line of data, and exit if reach the end of file
  if (temp(1:1)=='!')cycle
  ! (check if the reached the start of the site characteristix definitions)
  if (trim(temp)=='<start_'//trim(site_tag)//'>')then
   site_flag=.true.
   cycle  ! (don't do anything with the line that contains the tag)
  endif
  ! (check if the reached the start of the site characteristix definitions)
  if (trim(temp)=='<end_'//trim(site_tag)//'>')then
   site_flag=.false.   ! set site flag to false
   rewind(unt); exit   ! exit the do loop
  endif
  ! ---------------------------------------------------------------------------------------------
  ! put the site characteristic in the data structures
  ! ---------------------------------------------------------------------------------------------
  ! (if .true. here if there is a line that we want)
  if(site_flag)then
   ! (split the line into an array of words)
   call split_line(temp,chardata,err,cmessage)
   ! (identify the variable in the data structure)
   ivar = get_ixSite(trim(charData(1)))
   ! (check that the variable can in fact be identified)
   if(ivar>0)then
    ! (check that the variable index is within range of the data vector)
    if(ivar > size(site_data%var(:)) )then
     message=trim(message)//"index is out of range of the site data vector [var="//trim(varname)//"]"
     err=15; return
    endif
    ! (get the number of data elements associated with the current variable)
    nElements = size(charData)-1
    if(nElements /= size(site_data%var(iVar)%dat(:)) )then
     message=trim(message)//'unexpected number of elements in variable [var='//trim(chardata(1))//'] -- check spaces between each data element'
     err=20; return
    endif
    ! (put data in the structure)
    do iWord=1,nElements
     read(charData(iWord+1),*,iostat=err) site_data%var(iVar)%dat(iWord)
    end do
    !write(*,'(i4,1x,a30,1x,12(f9.3,1x))') ivar, ' --> '//trim(charData(1)), site_data%var(iVar)%dat
   ! (case where the variable is not found)
   else
    err=40; message=trim(message)//"variableNotFound[var="//trim(varname)//"]"; return
   endif
  endif  ! if the line is part of the site characteristics
 end do  ! looping through maximum number of lines in the file
 ! check we have populated all variables
 do ivar=1,size(site_data%var)
  if(site_data%var(ivar)%dat(1) > 0.99_dp*amiss)then
   err=40; message=trim(message)//"variableNonexistent[var="//trim(site_meta(ivar)%varname)//"]"; return
  endif
 end do


 ! **********************************************************************************************
 ! (3) read default model parameter values and constraints
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
 ! identify position of model parameters in the file
 ! ---------------------------------------------------------------------------------------------
 param_flag=.false. ! initialize flag to identify the parameter portion of the file
 ! loop through file until reach the parameter tag
 do iline=1,maxLines
  ! (read through comment lines)
  read(unt,'(a)',iostat=iend) temp  ! read a line of data
  if(iend/=0)then; err=20; message=trim(message)//'got to end of file before found the parameter block'; return; endif
  if (temp(1:1)=='!')cycle
  ! (check if the reached the start of the parameter definitions
  if(trim(temp) == '<start_'//trim(param_tag)//'>')then
   param_flag=.true.
   exit
  endif
 end do  ! looping through until identify the start of the parameter definitions
 if(.not.param_flag)then; err=20; message=trim(message)//'cannot find the start of the parameter definitions'; return; endif
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
  if(iend/=0)then; err=20; message=trim(message)//'got to end of file before found end of parameter data'; return; endif
  if (temp(1:1)=='!')cycle
  ! (check if have already reached the end of parameter descriptions -- if so, rewind and exit)
  if (trim(temp)=='<end_'//trim(param_tag)//'>')then
   site_flag=.false.   ! set site flag to false
   rewind(unt); exit   ! exit the do loop
  endif
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

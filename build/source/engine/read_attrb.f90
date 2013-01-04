module read_attrb_module
USE nrtype
implicit none
private
public::read_attrb
contains

 ! ************************************************************************************************
 ! (1) new subroutine: read information on local attributes
 ! ************************************************************************************************
 subroutine read_attrb(err,message)
 ! used to read metadata on the forcing data file
 USE snow_fileManager,only:SETNGS_PATH             ! path for metadata files
 USE snow_fileManager,only:LOCAL_ATTRIBUTES        ! file containing information on local attributes
 USE ascii_util_module,only:file_open              ! open ascii file
 USE ascii_util_module,only:split_line             ! extract the list of variable names from the character string
 USE allocspace_module,only:alloc_attr,alloc_type  ! allocate space for the local attributes
 USE data_struc,only:attr_meta,type_meta           ! metadata structures
 USE data_struc,only:attr_data,type_data           ! data structures
 USE var_lookup,only:iLookATTR,iLookTYPE           ! named variables for elements of the data structures
 USE get_ixname_module,only:get_ixAttr,get_ixType  ! access function to find index of elements in structure

 implicit none
 ! define output
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! define general variables
 real(dp),parameter                   :: amiss=-9999._dp ! missing data
 real(dp),parameter                   :: imiss=-9999     ! missing data
 character(len=256)                   :: cmessage        ! error message for downwind routine
 character(LEN=256)                   :: infile          ! input filename
 integer(i4b),parameter               :: unt=99          ! DK: need to either define units globally, or use getSpareUnit
 integer(i4b)                         :: iline           ! loop through lines in the file 
 integer(i4b),parameter               :: maxLines=1000   ! maximum lines in the file 
 character(LEN=256)                   :: temp            ! single lime of information
 ! define local variables
 integer(i4b)                         :: iend            ! check for the end of the file
 character(LEN=256)                   :: ffmt            ! file format
 character(LEN=32)                    :: varName         ! name of variable
 character(LEN=32)                    :: varValue        ! value for variable
 character(LEN=32)                    :: varType         ! type of variable
 character(len=1)                     :: dLim1,dLim2     ! column delimiters
 integer(i4b)                         :: ivar            ! index of model variable
 ! Start procedure here
 err=0; message="f-fuse/read_attrb/"
 ! **********************************************************************************************
 ! (1) open files, etc.
 ! **********************************************************************************************
 ! build filename
 infile = trim(SETNGS_PATH)//trim(LOCAL_ATTRIBUTES)
 ! open file
 call file_open(trim(infile),unt,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! allocate space
  call alloc_attr(err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  call alloc_type(err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! initialize data as missing
 attr_data%var(:) = amiss
 type_data%var(:) = imiss

 ! **********************************************************************************************
 ! (2) read local attributes
 ! **********************************************************************************************
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
  if(iLine==maxLines)then; err=20; message=trim(message)//'problem finding format code'; return; endif 
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
  read(temp,trim(ffmt),iostat=err) varName, dLim1, varValue, dLim2, varType
  if(err/=0)then; err=30; message=trim(message)//"errorInternalRead(entireLine)"; return; endif
  ! (check the deliminators are where we think they are)
  if(dLim1/='|' .or. dLim2/='|')then; err=20; message='incorrectDataFormat[line='//trim(temp)//']'; return; endif
  ! (identify the type of variable)
  select case(trim(varType))
   ! (put categorical data into the structure)
   case('categorical')
    iVar=get_ixType(trim(varName))
    if(iVar<=0)then; err=40; message=trim(message)//"variableNotFound[var="//trim(varname)//"]"; return; endif
    if(iVar>size(type_data%var(:)))then; err=40; message=trim(message)//"indexOutOfRange[var="//trim(varname)//"]"; return; endif
    read(varValue,*,iostat=err) type_data%var(iVar)
    if(err/=0)then; err=30; message=trim(message)//"errorInternalRead(varValue)[var="//trim(varname)//"]"; return; endif
   ! (put numerical data into the structure)
   case('numericData')
    iVar=get_ixAttr(trim(varName))
    if(iVar<=0)then; err=40; message=trim(message)//"variableNotFound[var="//trim(varname)//"]"; return; endif
    if(iVar>size(attr_data%var(:)))then; err=40; message=trim(message)//"indexOutOfRange[var="//trim(varname)//"]"; return; endif
    read(varValue,*,iostat=err) attr_data%var(iVar)
    if(err/=0)then; err=30; message=trim(message)//"errorInternalRead(varValue)[var="//trim(varname)//"]"; return; endif
   ! (check that we actually identify the correct data structure)
   case default; err=20; message=trim(message)//'unknown data structure [varType='//trim(varType)//', line='//trim(temp)//']'; return
  end select
 end do ! looping through lines in the file
 ! close file unit
 close(unt)
 ! check we have populated all real variables
 if(any(attr_data%var(:) < 0.99_dp*amiss))then
  do ivar=1,size(attr_data%var(:))
   if(attr_data%var(iVar) < 0.99_dp*amiss)then
    err=40; message=trim(message)//"variableNonexistent[var="//trim(attr_meta(ivar)%varname)//"]"; return
   endif
  end do
 endif
 ! check we have populated all integer variables
 if(any(type_data%var(:) == imiss))then
  do ivar=1,size(type_data%var(:))
   if(type_data%var(iVar) == imiss)then
    err=40; message=trim(message)//"variableNonexistent[var="//trim(type_meta(ivar)%varname)//"]"; return
   endif
  end do
 endif
 end subroutine read_attrb

end module read_attrb_module

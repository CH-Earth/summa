module read_attrb_module
USE nrtype
implicit none
private
public::read_attrb
contains

 ! ************************************************************************************************
 ! (1) new subroutine: read information on local attributes
 ! ************************************************************************************************
 subroutine read_attrb(nHRU,err,message)
 ! provide access to subroutines
 USE ascii_util_module,only:file_open              ! open ascii file
 USE ascii_util_module,only:split_line             ! extract the list of variable names from the character string
 USE ascii_util_module,only:get_vlines             ! read a vector of non-comment lines from an ASCII file
 USE allocspace_module,only:alloc_attr             ! module to allocate space for local attributes 
 USE allocspace_module,only:alloc_type             ! module to allocate space for categorical data
 ! provide access to data
 USE snow_fileManager,only:SETNGS_PATH             ! path for metadata files
 USE snow_fileManager,only:LOCAL_ATTRIBUTES        ! file containing information on local attributes
 USE data_struc,only:attr_meta,type_meta           ! metadata structures
 USE data_struc,only:attr_hru,type_hru             ! data structures
 USE var_lookup,only:iLookATTR,iLookTYPE           ! named variables for elements of the data structures
 USE get_ixname_module,only:get_ixAttr,get_ixType  ! access function to find index of elements in structure
 implicit none
 ! define output
 integer(i4b),intent(out)             :: nHRU        ! number of hydrologic response units
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! define general variables
 real(dp),parameter                   :: missingDouble=-9999._dp  ! missing data
 integer(i4b),parameter               :: missingInteger=-9999     ! missing data
 character(len=256)                   :: cmessage        ! error message for downwind routine
 character(LEN=256)                   :: infile          ! input filename
 integer(i4b),parameter               :: unt=99          ! DK: need to either define units globally, or use getSpareUnit
 integer(i4b)                         :: iline           ! loop through lines in the file 
 integer(i4b),parameter               :: maxLines=1000   ! maximum lines in the file 
 character(LEN=256)                   :: temp            ! single lime of information
 ! define local variables
 integer(i4b)                         :: iend            ! check for the end of the file
 character(LEN=512)                   :: nameString      ! string containing the list of attribute names
 character(LEN=256),allocatable       :: attNames(:)     ! vector of attribute names
 character(LEN=256),allocatable       :: attData(:)      ! vector of attribute data for a given HRU
 character(LEN=256),allocatable       :: dataLines(:)    ! vector of character strings from non-comment lines
 integer(i4b),parameter               :: categorical=101 ! named variable to denote categorical data
 integer(i4b),parameter               :: numerical=102   ! named variable to denote numerical data
 integer(i4b),allocatable             :: varType(:)      ! type of variable (categorical or numerical)
 integer(i4b),allocatable             :: varIndx(:)      ! index of variable within its data structure
 integer(i4b)                         :: iAtt            ! index of an attribute name 
 integer(i4b)                         :: iHRU            ! index of an HRU
 integer(i4b)                         :: nAtt            ! number of model attributes
 integer(i4b)                         :: nVar_attr       ! number of variables in the model attribute structure
 integer(i4b)                         :: nVar_type       ! number of variables in the model category structure
 logical(lgt),allocatable             :: checkType(:)    ! vector to check if we have all desired categorical values
 logical(lgt),allocatable             :: checkAttr(:)    ! vector to check if we have all desired local attributes
 ! Start procedure here
 err=0; message="read_attrb/"

 ! **********************************************************************************************
 ! (0) get number of variables in each data structure
 ! **********************************************************************************************
 ! check that metadata structures are initialized
 if(.not.associated(attr_meta) .or. .not.associated(type_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 nVar_attr = size(attr_meta)
 nVar_type = size(type_meta)
 ! allocate space for the check vectors
 allocate(checkType(nVar_type),checkAttr(nVar_attr),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for variable check vectors'; return; endif
 checkType(:) = .false.
 checkAttr(:) = .false.


 ! **********************************************************************************************
 ! (1) open files, etc.
 ! **********************************************************************************************
 ! build filename
 infile = trim(SETNGS_PATH)//trim(LOCAL_ATTRIBUTES)
 ! open file
 call file_open(trim(infile),unt,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


 ! **********************************************************************************************
 ! (2) read local attributes
 ! **********************************************************************************************
 ! ---------------------------------------------------------------------------------------------
 ! read attribute names
 ! ---------------------------------------------------------------------------------------------
 do iline=1,maxLines
  ! (read through comment lines)
  read(unt,'(a)',iostat=iend) temp  ! read a line of data
  if(iend/=0)then; err=20; message=trim(message)//'got to end of file before found the format code'; return; endif
  if (temp(1:1)=='!')cycle
  ! (read in format string -- assume that the first non-comment line is the list of attribute names)
  read(temp,'(a)')nameString  ! read in list of attribute names
  exit
  if(iLine==maxLines)then; err=20; message=trim(message)//'problem finding list of attribute names'; return; endif
 end do ! looping through lines
 ! ---------------------------------------------------------------------------------------------
 ! identify the type of each attribute
 ! ---------------------------------------------------------------------------------------------
 ! split the line into an array of words
 call split_line(nameString,attNames,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! identify the number of attributes
 nAtt = size(attNames)
 ! allocate space for the variable type and index
 allocate(varType(nAtt),varIndx(nAtt), stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the variable type and index'; return; endif
 ! initialize variables as missing
 varType(:) = missingInteger
 varIndx(:) = missingInteger
 ! loop through the attribute names
 do iAtt=1,nAtt
  ! find attribute name
  select case(trim(attNames(iAtt)))
   ! categorical data
   case('hruIndex','vegTypeIndex','soilTypeIndex','slopeTypeIndex','downHRUindex')
    varType(iAtt) = categorical
    varIndx(iAtt) = get_ixType(attNames(iAtt))
    checkType(varIndx(iAtt)) = .true.
   ! numerical data
   case('latitude','longitude','elevation','tan_slope','contourLength','HRUarea','mHeight')
    varType(iAtt) = numerical
    varIndx(iAtt) = get_ixAttr(attNames(iAtt))
    checkAttr(varIndx(iAtt)) = .true.
   ! check that variables are what we expect
   case default
    message=trim(message)//'unknown variable ['//trim(attNames(iAtt))//'] in local attributes file'
    err=20; return
  end select
  ! check that the variable could be identified in the data structure
  if(varIndx(iAtt) < 1)then; err=20; message=trim(message)//'unable to find variable ['//trim(attNames(iAtt))//'] in data structure'; return; endif
  ! print progress
  !print*, (varType(iAtt)==categorical), varIndx(iAtt), trim(attNames(iAtt))
 end do  ! (looping through attribute names)
 ! check that we have all desired categorical variables
 if(any(.not.checkType))then
  do iAtt=1,nVar_type
   if(.not.checkType(iAtt))then; err=20; message=trim(message)//'missing variable ['//trim(type_meta(iAtt)%varname)//'] in local attributes file'; return; endif
  end do
 endif
 ! check that we have all desired local attributes
 if(any(.not.checkAttr))then
  do iAtt=1,nVar_attr
   if(.not.checkAttr(iAtt))then; err=20; message=trim(message)//'missing variable ['//trim(attr_meta(iAtt)%varname)//'] in local attributes file'; return; endif
  end do
 endif


 ! **********************************************************************************************
 ! (3) read attributes for each HRU, and allocate space
 ! **********************************************************************************************
 ! get a list of character strings from non-comment lines
 call get_vlines(unt,dataLines,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! get the number of HRUs
 nHRU = size(dataLines)
 ! allocate space
 call alloc_attr(nHRU,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 call alloc_type(nHRU,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


 ! **********************************************************************************************
 ! (4) put data in the structures
 ! **********************************************************************************************
 ! loop through HRUs
 do iHRU=1,nHRU
  ! split the line into an array of words
  call split_line(dataLines(iHRU),attData,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  if(size(attData) /= nAtt)then; err=20; message=trim(message)//'number of attributes does not match expected number of attributes'; return; endif
  ! put attributes in the appropriate structures
  do iAtt=1,nAtt
   select case(varType(iAtt))
    case(numerical);   read(attData(iAtt),*,iostat=err) attr_hru(iHRU)%var(varIndx(iAtt))
    case(categorical); read(attData(iAtt),*,iostat=err) type_hru(iHRU)%var(varIndx(iAtt))
    case default; err=20; message=trim(message)//'unable to find type of attribute (categorical or numerical)'; return
   end select
   if(err/=0)then; err=20; message=trim(message)//'problem with internal read of attribute data'; return; endif
  end do  ! (looping through model attributes)
 end do  ! (looping through HRUs)
 
 ! **********************************************************************************************
 ! (5) deallocate space
 ! **********************************************************************************************
 deallocate(attNames,attData,dataLines,varType,varIndx,checkType,checkAttr, stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem deallocating space'; return; endif

 ! test
 !do iHRU=1,nHRU
 ! print*, '*****'
 ! print*, 'hruIndex       = ', type_hru(iHRU)%var(iLookTYPE%hruIndex)
 ! print*, 'latitude       = ', attr_hru(iHRU)%var(iLookATTR%latitude)
 ! print*, 'longitude      = ', attr_hru(iHRU)%var(iLookATTR%longitude)
 ! print*, 'elevation      = ', attr_hru(iHRU)%var(iLookATTR%elevation)
 ! print*, 'vegTypeIndex   = ', type_hru(iHRU)%var(iLookTYPE%vegTypeIndex)
 ! print*, 'soilTypeIndex  = ', type_hru(iHRU)%var(iLookTYPE%soilTypeIndex)
 ! print*, 'slopeTypeIndex = ', type_hru(iHRU)%var(iLookTYPE%slopeTypeIndex)
 !end do ! (looping through HRUs)
 !pause


 end subroutine read_attrb

end module read_attrb_module

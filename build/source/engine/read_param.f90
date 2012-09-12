module read_param_module
USE nrtype
implicit none
private
public::read_param
contains

 ! ************************************************************************************************
 ! (1) new subroutine: read trial model parameter values
 ! ************************************************************************************************
 subroutine read_param(nParSets,err,message)
 ! used to read model initial conditions
 USE snow_fileManager,only:SETNGS_PATH               ! path for metadata files
 USE snow_fileManager,only:PARAMETER_TRIAL           ! file with parameter trial values
 USE ascii_util_module,only:file_open                ! open file
 USE ascii_util_module,only:split_line               ! extract the list of variable names from the character string
 USE ascii_util_module,only:get_vlines               ! get a list of character strings from non-comment lines
 USE allocspace_module,only:alloc_mpar               ! allocate space for model parameters
 USE data_struc,only:mpar_data,mpar_sets,parFallback ! data for model parameter sets
 USE get_ixname_module,only:get_ixparam              ! access function to find index of elements in structure
 USE var_lookup,only:iLookPARAM                      ! named variables to index elements of the data vector
 implicit none
 ! define output
 integer(i4b),intent(out)       :: nParSets          ! number of parameter sets
 integer(i4b),intent(out)       :: err               ! error code
 character(*),intent(out)       :: message           ! error message
 ! define local variables
 character(len=256)             :: cmessage          ! error message for downwind routine
 character(LEN=256)             :: infile            ! input filename
 integer(i4b),parameter         :: unt=99            ! DK: need to either define units globally, or use getSpareUnit
 integer(i4b)                   :: iline             ! loop through lines in the file 
 integer(i4b),parameter         :: maxLines=1000     ! maximum lines in the file 
 character(LEN=256)             :: temp              ! single line of information
 integer(i4b)                   :: iend              ! check for the end of the file
 character(LEN=256),allocatable :: varnames(:)       ! vector of variable names
 character(LEN=256),allocatable :: charline(:)       ! vector of character strings
 character(LEN=256),allocatable :: chardata(:)       ! vector of character data
 integer(i4b)                   :: ipar,jpar         ! index of model parameter
 integer(i4b)                   :: nPars             ! number of model parameters
 integer(i4b)                   :: iParSet           ! index of parameter sets
 ! Start procedure here
 err=0; message="f-fuse/read_param/"
 ! **********************************************************************************************
 ! (1) open files, etc.
 ! **********************************************************************************************
 ! build filename
 infile = trim(SETNGS_PATH)//trim(PARAMETER_TRIAL)
 ! open file
 call file_open(trim(infile),unt,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! **********************************************************************************************
 ! (2) read the parameter names
 ! **********************************************************************************************
 ! loop through file until reach the first non-comment line (list of variable names)
 do iline=1,maxLines
  read(unt,'(a)',iostat=iend)temp; if(iend/=0)exit    ! read line of data
  if (temp(1:1)=='!')cycle
  ! extract the list of variable names from the character string
  call split_line(temp,varnames,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  exit
 end do
 ! save the number of parameters
 nPars = size(varnames)
 ! **********************************************************************************************
 ! (3) read the initial conditions data (continue reading from previous point in the file)
 ! **********************************************************************************************
 ! get a list of character strings from non-comment lines
 call get_vlines(unt,charline,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! **********************************************************************************************
 ! (4) allocate space for the model variable and model index vectors
 ! **********************************************************************************************
 nParSets = size(charline)
 call alloc_mpar(nParSets,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! **********************************************************************************************
 ! (5) populate the model parameter vectors
 ! **********************************************************************************************
 ! check that the default parameter structures exist
 if(.not.associated(parFallback))then
  err=20;message=trim(message)//"parFallbackUninitialized"; return
 endif
 ! allocate space for the character data
 allocate(chardata(nPars),stat=err)
 if(err/=0)then;err=30;message=trim(message)//"problemAllocateChardata"; return; endif
 ! loop through parameter sets
 do iParSet=1,nParSets
  ! assign mpar_data to the given parameter set
  mpar_data => mpar_sets(iParSet)
  ! ***** populate parameter set with default model parameters *****
  mpar_data%var(:) = parFallback(:)%default_val
  ! get the vector of parameters for a given layer
  read(charline(iParSet),*,iostat=err) chardata
  if(err/=0)then;err=40;message=trim(message)//"problemInternalRead[data='"//trim(charline(iParSet))//"']"; return; endif
  ! loop through the model parameters
  do ipar=1,nPars
   ! get the variable index
   jpar = get_ixparam(trim(varnames(ipar)))
   if(jpar<=0)then; err=40; message=trim(message)//"cannotFindVariableIndex[name='"//trim(varnames(ipar))//"']"; return; endif
   ! populate the appropriate element of the parameter vector
   read(chardata(ipar),*,iostat=err) mpar_data%var(jpar)
   if(err/=0)then;err=40;message=trim(message)//"problemInternalRead[data='"//trim(chardata(ipar))//"']"; return; endif
  end do    ! (looping through model parameters)
 end do    ! (looping through model parameter sets)
 ! **********************************************************************************************
 deallocate(varnames,charline,chardata,stat=err)
 if(err/=0)then;err=30;message=trim(message)//"problemDeallocate"; return; endif
 ! **********************************************************************************************
 end subroutine read_param

end module read_param_module

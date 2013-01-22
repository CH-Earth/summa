module allocspace_module
USE nrtype
implicit none
private
public::init_metad
public::alloc_time
public::alloc_forc
public::alloc_attr
public::alloc_type
public::alloc_mpar
public::alloc_mvar
public::alloc_indx
! define missing values
real(dp),parameter ::missingInteger=-9999
real(dp),parameter ::missingDouble=-9999._dp
contains

 ! ************************************************************************************************
 ! new subroutine: initialize metadata structures
 ! ************************************************************************************************
 subroutine init_metad(err,message)
 ! used to initialize the metadata structures
 USE var_lookup,only:maxvarTime,maxvarForc,maxvarAttr,maxvarType,maxvarMpar,maxvarMvar,maxvarIndx   ! maximum number variables in each data structure
 USE data_struc,only:time_meta,forc_meta,attr_meta,type_meta,mpar_meta,mvar_meta,indx_meta          ! metadata structures
 implicit none
 ! declare variables
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! initialize errors
 err=0; message="f-init_model/"
 ! ensure metadata structures are deallocated
 if (associated(time_meta)) deallocate(time_meta)
 if (associated(forc_meta)) deallocate(forc_meta)
 if (associated(attr_meta)) deallocate(attr_meta)
 if (associated(type_meta)) deallocate(type_meta)
 if (associated(mpar_meta)) deallocate(mpar_meta)
 if (associated(mvar_meta)) deallocate(mvar_meta)
 if (associated(indx_meta)) deallocate(indx_meta)
 ! allocate metadata structures
 allocate(time_meta(maxvarTime),forc_meta(maxvarForc),attr_meta(maxvarAttr),type_meta(maxvarType),&
          mpar_meta(maxvarMpar),mvar_meta(maxvarMvar),indx_meta(maxvarIndx),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateMetadata"; return; endif
 end subroutine init_metad

 ! ************************************************************************************************
 ! new subroutine: initialize data structures for time structures
 ! ************************************************************************************************
 subroutine alloc_time(datastr,err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:var_i,time_meta                 ! data structures
 implicit none
 ! dummy variables
 type(var_i),intent(out),pointer      :: datastr     ! data structure to allocate
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! initialize errors
 err=0; message="f-alloc_time/"
 ! check that the metadata structure is allocated
 if(.not.associated(time_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(datastr)) deallocate(datastr)
 allocate(datastr,stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! initialize second level data structure
 allocate(datastr%var(size(time_meta)),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
 ! set values to missing
 datastr%var(:) = missingInteger
 end subroutine alloc_time

 
 ! ************************************************************************************************
 ! new subroutine: initialize data structures for model forcing data
 ! ************************************************************************************************
 subroutine alloc_forc(err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:forc_data,forc_meta             ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! initialize errors
 err=0; message="f-alloc_forc/"
 ! check that the metadata structure is allocated
 if(.not.associated(forc_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(forc_data)) deallocate(forc_data)
 allocate(forc_data,stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! initialize second level data structure
 allocate(forc_data%var(size(forc_meta)),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
 ! set values to missing
 forc_data%var(:) = missingDouble
 end subroutine alloc_forc


 ! ************************************************************************************************
 ! new subroutine: initialize data structures for local attributes
 ! ************************************************************************************************
 subroutine alloc_attr(err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:attr_data,attr_meta             ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: ivar        ! loop throough variables
 integer(i4b)                         :: nVar        ! number of variables
 ! initialize errors
 err=0; message="f-alloc_attr/"
 ! check that the metadata structure is allocated
 if(.not.associated(attr_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(attr_data)) deallocate(attr_data)
 allocate(attr_data,stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! initialize second level data structure
 nVar = size(attr_meta)
 allocate(attr_data%var(nVar),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
 ! fill data with missing values
 attr_data%var(:) = missingDouble
 end subroutine alloc_attr


 ! ************************************************************************************************
 ! new subroutine: initialize data structures for local classification of veg, soil, etc.
 ! ************************************************************************************************
 subroutine alloc_type(err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:type_data,type_meta             ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: ivar        ! loop throough variables
 integer(i4b)                         :: nVar        ! number of variables
 ! initialize errors
 err=0; message="f-alloc_type/"
 ! check that the metadata structure is allocated
 if(.not.associated(type_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(type_data)) deallocate(type_data)
 allocate(type_data,stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! initialize second level data structure
 nVar = size(type_meta)
 allocate(type_data%var(nVar),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
 ! fill data with missing values
 type_data%var(:) = missingInteger
 end subroutine alloc_type


 ! ************************************************************************************************
 ! new subroutine: initialize data structures for model parameters
 ! ************************************************************************************************
 subroutine alloc_mpar(nParSets,err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:mpar_sets,mpar_meta             ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(in)              :: nParSets    ! number of parameter sets
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: iparset     ! loop through parameter sets
 integer(i4b)                         :: nPar        ! number of parameters
 ! initialize errors
 err=0; message="f-alloc_mpar/"
 ! check that the metadata structure is allocated
 if(.not.associated(mpar_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(mpar_sets)) deallocate(mpar_sets)
 allocate(mpar_sets(nParSets),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! get the number of parameters
 nPar = size(mpar_meta)
 ! loop through parameter sets
 do iparset=1,nParSets
  ! initialize second level data structure
  allocate(mpar_sets(iparset)%var(nPar),stat=err)
  if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
  ! set values to missing
  mpar_sets(iparset)%var(:) = missingDouble
 end do  ! looping through parameter sets
 end subroutine alloc_mpar

 ! ************************************************************************************************
 ! new subroutine: initialize data structures for model variables
 ! ************************************************************************************************
 subroutine alloc_mvar(nSnow,nSoil,err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:mvar_data,mvar_meta             ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(in)              :: nSnow       ! number of snow layers
 integer(i4b),intent(in)              :: nSoil       ! number of soil layers
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: ivar        ! loop throough variables
 integer(i4b)                         :: nVar        ! number of variables
 integer(i4b)                         :: nLayers     ! total number of layers
 integer(i4b),parameter               :: nBand=2     ! number of spectral bands
 integer(i4b),parameter               :: nTimeDelay=1000  ! number of elements in the time delay histogram
 ! initialize errors
 err=0; message="f-alloc_mvar/"
 ! check that the metadata structure is allocated
 if(.not.associated(mvar_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(mvar_data)) deallocate(mvar_data)
 allocate(mvar_data,stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! initialize second level data structure
 nVar = size(mvar_meta)
 allocate(mvar_data%var(nVar),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
 ! compute the number of layers
 nLayers = nSnow + nSoil
 ! loop through variables
 do ivar=1,nVar
  print*,mvar_meta(ivar)%varname
  if (associated(mvar_data%var(ivar)%dat)) deallocate(mvar_data%var(ivar)%dat)
  select case(mvar_meta(ivar)%vartype)
   case('scalarv'); allocate(mvar_data%var(ivar)%dat(1),stat=err)
   case('wLength'); allocate(mvar_data%var(ivar)%dat(nBand),stat=err)
   case('midSnow'); allocate(mvar_data%var(ivar)%dat(nSnow),stat=err)
   case('midSoil'); allocate(mvar_data%var(ivar)%dat(nSoil),stat=err)
   case('midToto'); allocate(mvar_data%var(ivar)%dat(nLayers),stat=err)
   case('ifcSnow'); allocate(mvar_data%var(ivar)%dat(0:nSnow),stat=err)
   case('ifcSoil'); allocate(mvar_data%var(ivar)%dat(0:nSoil),stat=err)
   case('ifcToto'); allocate(mvar_data%var(ivar)%dat(0:nLayers),stat=err)
   case('routing'); allocate(mvar_data%var(ivar)%dat(nTimeDelay),stat=err)
   case default
    err=40; message=trim(message)//"unknownVariableType[name='"//trim(mvar_meta(ivar)%varname)//"'; &
                                   &type='"//trim(mvar_meta(ivar)%vartype)//"']"; return
  endselect
  if(err/=0)then;err=30;message=trim(message)//"problemAllocate[var='"//trim(mvar_meta(ivar)%varname)//"']"; return; endif
  ! fill data with missing values
  mvar_data%var(ivar)%dat(:) = missingDouble
 enddo ! (looping through variables in the data structure)
 end subroutine alloc_mvar

 ! ************************************************************************************************
 ! new subroutine: initialize structure components for model indices
 ! ************************************************************************************************
 subroutine alloc_indx(nLayers,err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:indx_data,indx_meta             ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(in)              :: nLayers     ! number of layers
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: ivar        ! loop throough variables
 integer(i4b)                         :: nVar        ! number of variables
 ! initialize errors
 err=0; message="f-alloc_indx/"
 ! check that the metadata structure is allocated
 if(.not.associated(indx_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(indx_data)) deallocate(indx_data)
 allocate(indx_data,stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! initialize second level data structure
 nVar = size(indx_meta)
 allocate(indx_data%var(nVar),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
 ! loop through variables
 do ivar=1,nVar
  if (associated(indx_data%var(ivar)%dat)) deallocate(indx_data%var(ivar)%dat)
  select case(indx_meta(ivar)%vartype)
   case('scalarv'); allocate(indx_data%var(ivar)%dat(1),stat=err)
   case('midToto'); allocate(indx_data%var(ivar)%dat(nLayers),stat=err)
   case default
    err=40; message=trim(message)//"unknownVariableType[name='"//trim(indx_meta(ivar)%varname)//"'; &
                                   &type='"//trim(indx_meta(ivar)%vartype)//"']"; return
  endselect
  if(err/=0)then;err=30;message=trim(message)//"problemAllocate[var='"//trim(indx_meta(ivar)%varname)//"']"; return; endif
  ! fill data with missing values
  indx_data%var(ivar)%dat(:) = missingInteger
 enddo ! (looping through variables in the data structure)
 end subroutine alloc_indx


end module allocspace_module

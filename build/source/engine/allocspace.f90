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

module allocspace_module
! data types
USE nrtype
! provide access to the derived types to define the data structures
USE data_types,only:&
                    ! no spatial dimension
                    var_i,               & ! x%var(:)            (i4b)
                    var_d,               & ! x%var(:)            (dp)
                    var_ilength,         & ! x%var(:)%dat        (i4b)
                    var_dlength,         & ! x%var(:)%dat        (dp)
                    ! gru dimension
                    gru_int,             & ! x%gru(:)%var(:)     (i4b)
                    gru_double,          & ! x%gru(:)%var(:)     (dp)
                    gru_intVec,          & ! x%gru(:)%var(:)%dat (i4b)
                    gru_doubleVec,       & ! x%gru(:)%var(:)%dat (dp)
                    ! gru+hru dimension
                    gru_hru_int,         & ! x%gru(:)%hru(:)%var(:)     (i4b)
                    gru_hru_double,      & ! x%gru(:)%hru(:)%var(:)     (dp)
                    gru_hru_intVec,      & ! x%gru(:)%hru(:)%var(:)%dat (i4b)
                    gru_hru_doubleVec      ! x%gru(:)%hru(:)%var(:)%dat (dp)
! metadata structure
USE data_types,only:var_info               ! data type for metadata
implicit none
private
public::allocGlobal
public::allocLocal
public::allocate_gru_struc
! define missing values
integer(i4b),parameter :: missingInteger=-9999
real(dp),parameter     :: missingDouble=-9999._dp
! define fixed dimensions
integer(i4b),parameter :: nBand=2         ! number of spectral bands
integer(i4b),parameter :: nTimeDelay=2000 ! number of elements in the time delay histogram
! -----------------------------------------------------------------------------------------------------------------------------------
contains

 ! ************************************************************************************************
 ! public subroutine allocate_gru_struc: allocate space for GRU-HRU mapping structures
 ! ************************************************************************************************
 subroutine allocate_gru_struc(nGRU,startGRU,nHRU,err,message)
 USE netcdf
 USE nr_utility_module,only:arth
 ! provide access to subroutines
 USE netcdf_util_module,only:nc_file_open          ! open netCDF file
 ! provide access to data
 USE summaFileManager,only:SETNGS_PATH             ! path for metadata files
 USE summaFileManager,only:LOCAL_ATTRIBUTES        ! file containing information on local attributes
 USE globalData,only: gru_struc                    ! gru-hru mapping structures
 use multiconst,only: integerMissing               ! value for missing integer
 
 implicit none
 ! define output
 integer(i4b),intent(inout)           :: nGRU               ! number of grouped response units
 integer(i4b),intent(inout)           :: nHRU               ! number of hydrologic response units
 integer(i4b),intent(in)              :: startGRU           ! index of the starting GRU for parallelization run
 integer(i4b),intent(out)             :: err                ! error code
 character(*),intent(out)             :: message            ! error message
 ! define general variables
 character(len=256)                   :: cmessage           ! error message for downwind routine
 character(LEN=256)                   :: infile             ! input filename
 
 ! define variables for NetCDF file operation
 integer(i4b)                         :: mode               ! netCDF file open mode
 integer(i4b)                         :: ncid               ! integer variables for NetCDF IDs
 integer(i4b)                         :: varid              ! variable id from netcdf file
 integer(i4b)                         :: maxHRU,maxGRU      ! total numbers of HRUs and GRUs in the netCDF file
 integer(i4b)                         :: hruDimID           ! integer variables for NetCDF IDs
 integer(i4b)                         :: gruDimID           ! integer variables for NetCDF IDs
 integer(i4b),allocatable             :: gru_id(:)          ! unique ids of GRUs stored in the netCDF file
 integer(i4b),allocatable             :: hru2gru_id(:)      ! unique GRU ids at each HRU
      
 ! define local variables
 integer(i4b)                         :: hruCount           ! number of hrus in a gru
 integer(i4b)                         :: iGRU               ! loop index of GRU

 ! Start procedure here
 err=0; message="allocate_gru_hru_map/"

 ! check that gru_struc structure is initialized
 if(allocated(gru_struc)) deallocate(gru_struc)

 ! build filename
 infile = trim(SETNGS_PATH)//trim(LOCAL_ATTRIBUTES)
 ! open file
 mode=nf90_noWrite
 call nc_file_open(trim(infile), mode, ncid, err, cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 
 ! get gru_ix dimension length
 err = nf90_inq_dimid(ncid, "gru", gruDimID);             if(err/=0)then; message=trim(message)//'problem finding gru dimension'; return; endif
 err = nf90_inquire_dimension(ncid, gruDimID, len = maxGRU); if(err/=0)then; message=trim(message)//'problem reading gru dimension'; return; endif
  
 if (nGRU == integerMissing) then
  ! this is a full run
  nGRU = maxGRU
 else
  ! check startGRU and nGRU
  if (startGRU+nGRU-1>maxGRU) then 
   ! try to reduce nGRU
   if (startGRU<=maxGRU) then 
    nGRU=maxGRU-startGRU+1
   else
    err=1; message=trim(message)//'startGRU is larger than then GRU dimension'; return
   end if
  end if 
 end if 
 
 ! get hru_dim dimension length (ie., the number of hrus in entire domain)
 err = nf90_inq_dimid(ncid, "hru", hruDimID);             if(err/=0)then; message=trim(message)//'problem finding hru dimension'; return; endif
 err = nf90_inquire_dimension(ncid, hruDimID, len = maxHRU); if(err/=0)then; message=trim(message)//'problem reading hru dimension'; return; endif  
 
 ! read the HRU to GRU mapping
 allocate(hru2gru_id(maxHRU),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for zLocalAttributes gru-hru correspondence vectors'; return; endif
 err = nf90_inq_varid(ncid, "hru2gruId", varid); if(err/=0)then; message=trim(message)//'problem finding hru2gruId variable'; return; endif
 err = nf90_get_var(ncid,varid,hru2gru_id);      if(err/=0)then; message=trim(message)//'problem reading hru2gruId variable'; return; endif
 
 ! allocate mapping array 
 allocate(gru_struc(nGRU))  
 if (nHRU == integerMissing) then
  ! get GRU Id
  allocate(gru_id(nGRU))  
  err = nf90_inq_varid(ncid, "gruId", varid);     if(err/=0)then; message=trim(message)//'problem finding gruId variable'; return; endif
  err = nf90_get_var(ncid,varid,gru_id,start=(/startGRU/),count=(/nGRU/));          if(err/=0)then; message=trim(message)//'problem reading gruId variable'; return; endif  
  ! allocate HRUs for each GRU 
  do iGRU=1,nGRU
   hruCount = count(hru2gru_id==gru_id(iGRU))
   gru_struc(iGRU)%hruCount = hruCount
   allocate(gru_struc(iGRU)%hruInfo(hruCount),stat=err)
   if(err/=0)then; err=20; message=trim(message)//'problem allocating space for hru in gru_struc'; return; endif
  end do
  nHRU = sum(gru_struc%hruCount)
 else  
  ! allocate space for single-HRU run
  gru_struc(1)%hruCount=1
  allocate(gru_struc(1)%hruInfo(1))
 endif
 
 ! close the HRU_ATTRIBUTES netCDF file
 err = nf90_close(ncid)
 if(err/=0)then; err=20; message=trim(message)//'error closing zLocalAttributes file'; return; endif

 ! check allocation was successful
 if(.not.allocated(gru_struc))then
  message=trim(message)//'gru_struc is not allocated'
  err=20; return
 endif

 end subroutine allocate_gru_struc

 ! ************************************************************************************************
 ! public subroutine allocGlobal: allocate space for global data structures 
 ! ************************************************************************************************
 subroutine allocGlobal(metaStruct,dataStruct,err,message)
 USE globalData,only: gru_struc                    ! gru-hru mapping structures
 implicit none
 ! input
 type(var_info),intent(in)       :: metaStruct(:)  ! metadata structure
 ! output
 class(*),intent(out)            :: dataStruct     ! data structure
 integer(i4b),intent(out)        :: err            ! error code
 character(*),intent(out)        :: message        ! error message
 ! local variables
 logical(lgt)                    :: check          ! .true. if structure is already allocated
 integer(i4b)                    :: iHRU           ! loop through HRUs
 integer(i4b)                    :: iGRU           ! loop through GRUs
 integer(i4b)                    :: nGRU           ! number of GRUs
 logical(lgt)                    :: spatial        ! spatial flag
 character(len=256)              :: cmessage       ! error message of the downwind routine
 ! initialize error control
 err=0; message='allocGlobal/'
 
 ! initialize allocation check
 check=.false.
 
 ! get the number of GRUs
 nGRU = size(gru_struc)

 ! * allocate GRU dimension
 select type(dataStruct)
  ! gru dimension only
  type is (gru_int);           if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); endif
  type is (gru_intVec);        if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); endif
  type is (gru_double);        if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); endif
  type is (gru_doubleVec);     if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); endif
  ! gru+hru dimensions
  type is (gru_hru_int);       if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); endif
  type is (gru_hru_intVec);    if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); endif
  type is (gru_hru_double);    if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); endif
  type is (gru_hru_doubleVec); if(allocated(dataStruct%gru))then; check=.true.; else; allocate(dataStruct%gru(nGRU),stat=err); endif
 end select

 ! check errors
 if(check) then; err=20; message=trim(message)//'GRU structure was unexpectedly allocated already'; return; endif
 if(err/=0)then; err=20; message=trim(message)//'problem allocating GRU dimension'; return; endif

 ! * allocate HRU dimension
 do iGRU=1,nGRU
  ! allocate the HRU dimension
  select type(dataStruct)
   type is (gru_hru_int);       if(allocated(dataStruct%gru(iGRU)%hru))then; check=.true.; else; allocate(dataStruct%gru(iGRU)%hru(gru_struc(iGRU)%hruCount),stat=err); endif
   type is (gru_hru_intVec);    if(allocated(dataStruct%gru(iGRU)%hru))then; check=.true.; else; allocate(dataStruct%gru(iGRU)%hru(gru_struc(iGRU)%hruCount),stat=err); endif
   type is (gru_hru_double);    if(allocated(dataStruct%gru(iGRU)%hru))then; check=.true.; else; allocate(dataStruct%gru(iGRU)%hru(gru_struc(iGRU)%hruCount),stat=err); endif
   type is (gru_hru_doubleVec); if(allocated(dataStruct%gru(iGRU)%hru))then; check=.true.; else; allocate(dataStruct%gru(iGRU)%hru(gru_struc(iGRU)%hruCount),stat=err); endif
   class default  ! do nothing: It is acceptable to not be any of these specified cases
  end select
  ! check errors
  if(check) then; err=20; message=trim(message)//'HRU structure was unexpectedly allocated already'; return; endif
  if(err/=0)then; err=20; message=trim(message)//'problem allocating HRU dimension'; return; endif
 end do

 ! * allocate local data structures where there is a spatial dimension
 gruLoop: do iGRU=1,nGRU

  ! initialize the spatial flag
  spatial=.false.

  ! loop through HRUs
  hruLoop: do iHRU=1,gru_struc(iGRU)%hruCount

   ! get the number of snow and soil layers
   associate(&
   nSnow => gru_struc(iGRU)%hruInfo(iHRU)%nSnow, & ! number of snow layers for each HRU
   nSoil => gru_struc(iGRU)%hruInfo(iHRU)%nSoil  ) ! number of soil layers for each HRU

   ! allocate space for structures WITH an HRU dimension
   select type(dataStruct)
    type is (gru_hru_int);       call allocLocal(metaStruct,dataStruct%gru(iGRU)%hru(iHRU),nSnow,nSoil,err,cmessage); spatial=.true.
    type is (gru_hru_intVec);    call allocLocal(metaStruct,dataStruct%gru(iGRU)%hru(iHRU),nSnow,nSoil,err,cmessage); spatial=.true.
    type is (gru_hru_double);    call allocLocal(metaStruct,dataStruct%gru(iGRU)%hru(iHRU),nSnow,nSoil,err,cmessage); spatial=.true.
    type is (gru_hru_doubleVec); call allocLocal(metaStruct,dataStruct%gru(iGRU)%hru(iHRU),nSnow,nSoil,err,cmessage); spatial=.true.
    class default; exit hruLoop
   end select

   ! error check
   if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

   ! end association to info in data structures
   end associate

  end do hruLoop ! loop through HRUs

  ! allocate space for structures *WITHOUT* an HRU dimension
  select type(dataStruct)
   type is (gru_double);    call allocLocal(metaStruct,dataStruct%gru(iGRU),nSnow=0,nSoil=0,err=err,message=cmessage); spatial=.true.
   type is (gru_doubleVec); call allocLocal(metaStruct,dataStruct%gru(iGRU),nSnow=0,nSoil=0,err=err,message=cmessage); spatial=.true.
   class default
    if(.not.spatial) exit gruLoop  ! no need to allocate spatial dimensions if none exist for a given variable
    cycle gruLoop  ! can have an HRU dimension if we get to here
  end select

  ! error check
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 end do gruLoop ! loop through GRUs

 ! * allocate local data structures where there is no spatial dimension
 select type(dataStruct)
  type is (var_i);         call allocLocal(metaStruct,dataStruct,err=err,message=cmessage) 
  type is (var_d);         call allocLocal(metaStruct,dataStruct,err=err,message=cmessage)
  type is (var_ilength);   call allocLocal(metaStruct,dataStruct,err=err,message=cmessage)
  type is (var_dlength);   call allocLocal(metaStruct,dataStruct,err=err,message=cmessage)
  ! check identified the data type
  class default; if(.not.spatial)then; err=20; message=trim(message)//'unable to identify derived data type'; return; endif
 end select
 ! error check
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 end subroutine allocGlobal

 ! ************************************************************************************************
 ! public subroutine allocLocal: allocate space for local data structures
 ! ************************************************************************************************
 subroutine allocLocal(metaStruct,dataStruct,nSnow,nSoil,err,message)
 implicit none
 ! input-output
 type(var_info),intent(in)        :: metaStruct(:)  ! metadata structure
 class(*),intent(inout)           :: dataStruct     ! data structure
 ! optional input
 integer(i4b),intent(in),optional :: nSnow          ! number of snow layers
 integer(i4b),intent(in),optional :: nSoil          ! number of soil layers
 ! output
 integer(i4b),intent(out)         :: err            ! error code
 character(*),intent(out)         :: message        ! error message
 ! local
 logical(lgt)                     :: check          ! .true. if the variables are allocated
 integer(i4b)                     :: nVars          ! number of variables in the metadata structure
 integer(i4b)                     :: nLayers        ! total number of layers
 character(len=256)               :: cmessage       ! error message of the downwind routine
 ! initialize error control
 err=0; message='allocLocal/'

 ! get the number of variables in the metadata structure
 nVars = size(metaStruct)

 ! check if nSnow and nSoil are present
 if(present(nSnow) .or. present(nSoil))then
  ! check both are present
  if(.not.present(nSoil))then; err=20; message=trim(message)//'expect nSoil to be present when nSnow is present'; return; endif
  if(.not.present(nSnow))then; err=20; message=trim(message)//'expect nSnow to be present when nSoil is present'; return; endif
  nLayers = nSnow+nSoil

 ! It is possible that nSnow and nSoil are actually needed here, so we return an error if the optional arguments are missing when needed
 else
  select type(dataStruct)
   type is (var_ilength); err=20
   type is (var_dlength); err=20
  end select
  if(err/=0)then; message=trim(message)//'expect nSnow and nSoil to be present for variable-length data structures'; return; endif
 endif

 ! initialize allocation check
 check=.false.

 ! allocate the dimension for model variables
 select type(dataStruct)
  type is (var_i);       if(allocated(dataStruct%var))then; check=.true.; else; allocate(dataStruct%var(nVars),stat=err); endif; return
  type is (var_d);       if(allocated(dataStruct%var))then; check=.true.; else; allocate(dataStruct%var(nVars),stat=err); endif; return
  type is (var_ilength); if(allocated(dataStruct%var))then; check=.true.; else; allocate(dataStruct%var(nVars),stat=err); endif
  type is (var_dlength); if(allocated(dataStruct%var))then; check=.true.; else; allocate(dataStruct%var(nVars),stat=err); endif
  class default; err=20; message=trim(message)//'unable to identify derived data type for the variable dimension'; return
 end select
 ! check errors
 if(check) then; err=20; message=trim(message)//'structure was unexpectedly allocated already'; return; endif
 if(err/=0)then; err=20; message=trim(message)//'problem allocating'; return; endif

 ! allocate the dimension for model data
 select type(dataStruct)
  type is (var_ilength); call allocateDat_int(metaStruct,nSnow,nSoil,nLayers,dataStruct,err,cmessage) 
  type is (var_dlength); call allocateDat_dp( metaStruct,nSnow,nSoil,nLayers,dataStruct,err,cmessage) 
  class default; err=20; message=trim(message)//'unable to identify derived data type for the data dimension'; return
 end select
 ! check errors
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine allocLocal

 ! ************************************************************************************************
 ! private subroutine allocateDat_dp: initialize data dimension of the data structures
 ! ************************************************************************************************
 subroutine allocateDat_dp(metadata,nSnow,nSoil,nLayers, & ! input
                           varData,err,message)            ! output
 implicit none
 ! input variables
 type(var_info),intent(in)         :: metadata(:) ! metadata structure
 integer(i4b),intent(in)           :: nSnow       ! number of snow layers
 integer(i4b),intent(in)           :: nSoil       ! number of soil layers
 integer(i4b),intent(in)           :: nLayers     ! total number of soil layers in the snow+soil domian (nSnow+nSoil)
 ! output variables
 type(var_dlength),intent(inout)   :: varData     ! model variables for a local HRU
 integer(i4b),intent(out)          :: err         ! error code
 character(*),intent(out)          :: message     ! error message
 ! local variables
 integer(i4b)                      :: iVar        ! variable index
 ! initialize error control
 err=0; message='allocateDat_dp/'

 ! loop through variables in the data structure
 do iVar=1,size(metadata)

  ! check allocated
  if(allocated(varData%var(iVar)%dat))then
   message=trim(message)//'variable '//trim(metadata(iVar)%varname)//' is unexpectedly allocated'
   err=20; return

  ! allocate structures
  else
   select case(trim(metadata(iVar)%vartype))
    case('scalarv'); allocate(varData%var(iVar)%dat(1),stat=err)
    case('wLength'); allocate(varData%var(iVar)%dat(nBand),stat=err)
    case('midSnow'); allocate(varData%var(iVar)%dat(nSnow),stat=err)
    case('midSoil'); allocate(varData%var(iVar)%dat(nSoil),stat=err)
    case('midToto'); allocate(varData%var(iVar)%dat(nLayers),stat=err)
    case('ifcSnow'); allocate(varData%var(iVar)%dat(0:nSnow),stat=err)
    case('ifcSoil'); allocate(varData%var(iVar)%dat(0:nSoil),stat=err)
    case('ifcToto'); allocate(varData%var(iVar)%dat(0:nLayers),stat=err)
    case('routing'); allocate(varData%var(iVar)%dat(nTimeDelay),stat=err)
    case('unknown'); allocate(varData%var(iVar)%dat(0),stat=err)  ! unknown=special (and valid) case that is allocated later (initialize with zero-length vector)
    case default; err=40; message=trim(message)//"unknownVariableType[name='"//trim(metadata(iVar)%varname)//"'; type='"//trim(metadata(iVar)%vartype)//"']"; return
   endselect
   ! check error
   if(err/=0)then; err=20; message=trim(message)//'problem allocating variable '//trim(metadata(iVar)%varname); return; endif
   ! set to missing
   varData%var(iVar)%dat(:) = missingDouble
  endif  ! if not allocated

 end do  ! looping through variables

 end subroutine allocateDat_dp

 ! ************************************************************************************************
 ! private subroutine allocateDat_int: initialize data dimension of the data structures
 ! ************************************************************************************************
 subroutine allocateDat_int(metadata,nSnow,nSoil,nLayers, & ! input
                            varData,err,message)            ! output
 implicit none
 ! input variables
 type(var_info),intent(in)         :: metadata(:) ! metadata structure
 integer(i4b),intent(in)           :: nSnow       ! number of snow layers
 integer(i4b),intent(in)           :: nSoil       ! number of soil layers
 integer(i4b),intent(in)           :: nLayers     ! total number of soil layers in the snow+soil domian (nSnow+nSoil)
 ! output variables
 type(var_ilength),intent(inout)   :: varData     ! model variables for a local HRU
 integer(i4b),intent(out)          :: err         ! error code
 character(*),intent(out)          :: message     ! error message
 ! local variables
 integer(i4b)                      :: iVar        ! variable index
 character(len=1024)               :: cmessage    ! error message
 ! initialize error control
 err=0; message='allocateDat_int/'
 cmessage=''
 ! loop through variables in the data structure
 do iVar=1,size(metadata)

  ! check allocated
  if(allocated(varData%var(iVar)%dat))then
   message=trim(message)//'variable '//trim(metadata(iVar)%varname)//' is unexpectedly allocated'
   err=20; return

  ! allocate structures
  else
   select case(trim(metadata(iVar)%vartype))
    case('scalarv'); allocate(varData%var(iVar)%dat(1),stat=err,errmsg=cmessage)
    case('wLength'); allocate(varData%var(iVar)%dat(nBand),stat=err,errmsg=cmessage)
    case('midSnow'); allocate(varData%var(iVar)%dat(nSnow),stat=err,errmsg=cmessage)
    case('midSoil'); allocate(varData%var(iVar)%dat(nSoil),stat=err,errmsg=cmessage)
    case('midToto'); allocate(varData%var(iVar)%dat(nLayers),stat=err,errmsg=cmessage)
    case('ifcSnow'); allocate(varData%var(iVar)%dat(0:nSnow),stat=err,errmsg=cmessage)
    case('ifcSoil'); allocate(varData%var(iVar)%dat(0:nSoil),stat=err,errmsg=cmessage)
    case('ifcToto'); allocate(varData%var(iVar)%dat(0:nLayers),stat=err,errmsg=cmessage)
    case('routing'); allocate(varData%var(iVar)%dat(nTimeDelay),stat=err,errmsg=cmessage)
    case('unknown'); allocate(varData%var(iVar)%dat(0),stat=err,errmsg=cmessage)  ! unknown=special (and valid) case that is allocated later (initialize with zero-length vector)
    case default; err=40; message=trim(message)//"unknownVariableType[name='"//trim(metadata(iVar)%varname)//"'; type='"//trim(metadata(iVar)%vartype)//"']"; return
   endselect
   ! check error
   if(err/=0)then; err=20; message=trim(message)//'problem allocating variable '//trim(metadata(iVar)%varname)//NEW_LINE('A')//trim(cmessage); return; endif
   ! set to missing
   varData%var(iVar)%dat(:) = missingInteger
  endif  ! if not allocated

 end do  ! looping through variables

 end subroutine allocateDat_int

end module allocspace_module

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
USE globalData,only:gru_struc              ! gru-hru mapping structures
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
 subroutine allocate_gru_struc(nGRU,nHRU,err,message)
 USE netcdf
 USE nr_utility_module,only:arth
 ! provide access to subroutines
 USE netcdf_util_module,only:nc_file_open          ! open netCDF file
 ! provide access to data
 USE summaFileManager,only:SETNGS_PATH             ! path for metadata files
 USE summaFileManager,only:LOCAL_ATTRIBUTES        ! file containing information on local attributes
 USE globalData,only: index_map                    ! relating different indexing system
 USE globalData,only: gru_struc                    ! gru-hru mapping structures
 
 implicit none
 ! define output
 integer(i4b),intent(out)             :: nGRU               ! number of grouped response units
 integer(i4b),intent(out)             :: nHRU               ! number of hydrologic response units
 integer(i4b),intent(out)             :: err                ! error code
 character(*),intent(out)             :: message            ! error message
 ! define general variables
 character(len=256)                   :: cmessage           ! error message for downwind routine
 character(LEN=256)                   :: infile             ! input filename
 
 ! define variables for NetCDF file operation
 integer(i4b)                         :: mode               ! netCDF file open mode
 integer(i4b)                         :: ncid               ! integer variables for NetCDF IDs
 integer(i4b)                         :: varid              ! variable id from netcdf file
 integer(i4b)                         :: hruDimID           ! integer variables for NetCDF IDs
 integer(i4b)                         :: gruDimID           ! integer variables for NetCDF IDs
 integer(i4b),allocatable             :: hru_id(:)          ! unique id of hru over entire domain
 integer(i4b),allocatable             :: gru_id(:)          ! unique ids of GRUs
 integer(i4b),allocatable             :: hru2gru_id(:)      ! unique GRU ids at each HRU
      
 ! define local variables
 integer(i4b)                         :: hruCount           ! number of hrus in a gru
 integer(i4b)                         :: iGRU               ! index of a GRU and HRU

 ! Start procedure here
 err=0; message="allocate_gru_hru_map/"

 ! check that gru_struc structure is initialized
 if(allocated(gru_struc)) deallocate(gru_struc)
 if(allocated(index_map)) deallocate(index_map)

 ! build filename
 infile = trim(SETNGS_PATH)//trim(LOCAL_ATTRIBUTES)
 ! open file
 mode=nf90_noWrite
 call nc_file_open(trim(infile), mode, ncid, err, cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get gru_ix dimension length
 err = nf90_inq_dimid(ncid, "ngru", gruDimID);             if(err/=0)then; message=trim(message)//'problem finding nGRU dimension'; return; endif
 err = nf90_inquire_dimension(ncid, gruDimID, len = nGRU); if(err/=0)then; message=trim(message)//'problem reading nGRU dimension'; return; endif

 ! get hru_dim dimension length (ie., the number of hrus in entire domain)
 err = nf90_inq_dimid(ncid, "nhru", hruDimID);             if(err/=0)then; message=trim(message)//'problem finding nHRU dimension'; return; endif
 err = nf90_inquire_dimension(ncid, hruDimID, len = nHRU); if(err/=0)then; message=trim(message)//'problem reading nHRU dimension'; return; endif
 
 allocate(gru_struc(nGRU), index_map(nHRU), stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for mapping structures'; return; endif
 
 allocate(gru_id(nGRU),hru_id(nHRU),hru2gru_id(nHRU),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for zLocalAttributes gru-hru correspondence vectors'; return; endif

 ! read gru_id from netcdf file
 err = nf90_inq_varid(ncid, "gruId", varid);     if(err/=0)then; message=trim(message)//'problem finding gruId variable'; return; endif
 err = nf90_get_var(ncid,varid,gru_id);          if(err/=0)then; message=trim(message)//'problem reading gruId variable'; return; endif

 ! read the HRU to GRU mapping
 err = nf90_inq_varid(ncid, "hru2gruId", varid); if(err/=0)then; message=trim(message)//'problem finding hru2gruId variable'; return; endif
 err = nf90_get_var(ncid,varid,hru2gru_id);      if(err/=0)then; message=trim(message)//'problem reading hru2gruId variable'; return; endif

 ! allocate mapping array
 do iGRU=1,nGRU
  hruCount = count(hru2gru_id==gru_id(iGRU))
  gru_struc(iGRU)%hruCount = hruCount
  allocate(gru_struc(iGRU)%hruInfo(hruCount),stat=err)
  if(err/=0)then; err=20; message=trim(message)//'problem allocating space for hru in gru_struc'; return; endif
 enddo

 ! close the HRU_ATTRIBUTES netCDF file
 err = nf90_close(ncid)
 if(err/=0)then; err=20; message=trim(message)//'error closing zLocalAttributes file'; return; endif
 end subroutine allocate_gru_struc

 ! ************************************************************************************************
 ! public subroutine allocGlobal: allocate space for global data structures 
 ! ************************************************************************************************
 subroutine allocGlobal(metaStruct,dataStruct,err,message)
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
  end select
  ! check errors
  if(check) then; err=20; message=trim(message)//'HRU structure was unexpectedly allocated already'; return; endif
  if(err/=0)then; err=20; message=trim(message)//'problem allocating HRU dimension'; return; endif
 end do

 ! * allocate local data structures where there is a spatial dimension
 do iGRU=1,nGRU
  ! initialize the spatial flag
  spatial=.true.
  ! get the number of snow and soil layers
  associate(&
  nHRU  => gru_struc(iGRU)%hruCount,         & ! number of HRUs
  nSnow => gru_struc(iGRU)%hruInfo(:)%nSnow, & ! number of snow layers for each HRU
  nSoil => gru_struc(iGRU)%hruInfo(:)%nSoil  ) ! number of soil layers for each HRU
  ! allocate space
  select type(dataStruct)
   ! structures with an HRU dimension
   type is (gru_hru_int);       do iHRU=1,nHRU; call allocLocal(metaStruct,dataStruct%gru(iGRU)%hru(iHRU),nSnow(iHRU),nSoil(iHRU),err,cmessage); if(err/=0)exit; end do
   type is (gru_hru_intVec);    do iHRU=1,nHRU; call allocLocal(metaStruct,dataStruct%gru(iGRU)%hru(iHRU),nSnow(iHRU),nSoil(iHRU),err,cmessage); if(err/=0)exit; end do
   type is (gru_hru_double);    do iHRU=1,nHRU; call allocLocal(metaStruct,dataStruct%gru(iGRU)%hru(iHRU),nSnow(iHRU),nSoil(iHRU),err,cmessage); if(err/=0)exit; end do
   type is (gru_hru_doubleVec); do iHRU=1,nHRU; call allocLocal(metaStruct,dataStruct%gru(iGRU)%hru(iHRU),nSnow(iHRU),nSoil(iHRU),err,cmessage); if(err/=0)exit; end do
   ! structures without an HRU dimension
   type is (gru_double);    call allocLocal(metaStruct,dataStruct%gru(iGRU),err=err,message=cmessage)
   type is (gru_doubleVec); call allocLocal(metaStruct,dataStruct%gru(iGRU),err=err,message=cmessage)
   class default; spatial=.false.
  end select
  ! end association to info in data structures
  end associate
  ! error check
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
  ! check that we found the structure
  if(.not.spatial)exit
 end do  ! loop through GRUs

 ! * allocate local data structures where there is a spatial dimension
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

 ! check that nSnow and nSoil are not needed
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
    case('unknown'); allocate(varData%var(iVar)%dat(0),stat=err)  ! unknown=initialize with zero-length vector
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
 ! initialize error control
 err=0; message='allocateDat_int/'

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
    case('unknown'); allocate(varData%var(iVar)%dat(0),stat=err)  ! unknown=initialize with zero-length vector
    case default; err=40; message=trim(message)//"unknownVariableType[name='"//trim(metadata(iVar)%varname)//"'; type='"//trim(metadata(iVar)%vartype)//"']"; return
   endselect
   ! check error
   if(err/=0)then; err=20; message=trim(message)//'problem allocating variable '//trim(metadata(iVar)%varname); return; endif
   ! set to missing
   varData%var(iVar)%dat(:) = missingInteger
  endif  ! if not allocated

 end do  ! looping through variables

 end subroutine allocateDat_int

end module allocspace_module

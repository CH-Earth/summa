! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2020 NCAR/RAL; University of Saskatchewan; University of Washington
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

module read_attrb_module
USE nr_type
! provide access to global data
USE globalData,only:gru_struc                              ! gru->hru mapping structure
USE globalData,only:index_map                              ! hru->gru mapping structure
USE globalData,only:attr_meta,type_meta,id_meta,grid_meta  ! metadata structures
USE globalData,only:int8Missing                            ! missing long integer

implicit none
private
public::read_dimension
public::read_dimensionGrid
public::read_attrb

contains

 ! ************************************************************************************************
 ! public subroutine read_dimension: read HRU and GRU dimension information on local attributes
 ! ************************************************************************************************
 subroutine read_dimension(attrFile,fileGRU,fileHRU,nGRU,nHRU,err,message,startGRU,checkHRU)
 USE netcdf
 USE netcdf_util_module,only:nc_file_open                   ! open netcdf file
 USE netcdf_util_module,only:nc_file_close                  ! close netcdf file
 USE nr_utils_module ,only:arth                             ! use to build vectors with regular increments
 implicit none

 character(*),intent(in)              :: attrFile           ! name of attributed file
 integer(i4b),intent(out)             :: fileGRU            ! number of GRUs in the input file
 integer(i4b),intent(out)             :: fileHRU            ! number of HRUs in the input file
 integer(i4b),intent(inout)           :: nGRU               ! number of GRUs in the run space
 integer(i4b),intent(inout)           :: nHRU               ! number of HRUs in the run space
 integer(i4b),intent(out)             :: err                ! error code
 character(*),intent(out)             :: message            ! error message
 integer(i4b),intent(in),optional     :: startGRU           ! index of the starting GRU for parallelization run
 integer(i4b),intent(in),optional     :: checkHRU           ! index of the HRU for a single HRU run
 ! locals
 integer(i4b)                         :: sGRU               ! starting GRU
 integer(i4b)                         :: iHRU               ! HRU counting index
 integer(i4b)                         :: iGRU               ! GRU loop index
 integer(i8b),allocatable             :: gru_id(:),hru_id(:)! read gru/hru IDs in from attributes file
 integer(i8b),allocatable             :: hru2gru_id(:)      ! read hru->gru mapping in from attributes file
 integer(i4b),allocatable             :: hru_ix(:)          ! hru index for search
 ! define variables for NetCDF file operation
 integer(i4b)                         :: ncid               ! NetCDF file ID
 integer(i4b)                         :: varID              ! NetCDF variable ID
 integer(i4b)                         :: dimID              ! netcdf file dimension id
 integer(i4b),allocatable             :: nGlac(:)           ! number of glaciers in gru
 integer(i4b),allocatable             :: nWtld(:)           ! number of wetlands/lakes in gru
 character(len=256)                   :: cmessage           ! error message for downwind routine

 ! Start procedure here
 err=0; message="read_dimension/"

 ! check that we do not have conflicting flags
 if(present(startGRU).and.present(checkHRU))then; message=trim(message)//'startGRU and checkHRU both exist, which is not supported'; return; end if

 ! open nc file
 call nc_file_open(trim(attrFile),nf90_noWrite,ncid,err,cmessage)
 if(err/=nf90_noerr)then; message=trim(message)//trim(cmessage); return; end if

 ! *********************************************************************************************
 ! read and set GRU dimensions
 ! **********************************************************************************************
 ! get gru dimension of whole file
 err = nf90_inq_dimid(ncid,"gru",dimID);                   if(err/=nf90_noerr)then; message=trim(message)//'problem finding gru dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncid, dimID, len = fileGRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading gru dimension/'//trim(nf90_strerror(err)); return; end if

 ! get hru dimension of whole file
 err = nf90_inq_dimid(ncid,"hru",dimID);                   if(err/=nf90_noerr)then; message=trim(message)//'problem finding hru dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncid, dimID, len = fileHRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading hru dimension/'//trim(nf90_strerror(err)); return; end if

 ! get runtime GRU dimensions
 if (present(startGRU)) then
   if (nGRU < 1) then; err=20; message=trim(message)//'nGRU < 1 for a startGRU run'; return; end if
   sGRU = startGRU
 elseif (present(checkHRU)) then
   nGRU = 1
 else
   sGRU = 1
   nGRU = fileGRU
 endif

 ! check dimensions
 if (present(startGRU)) then
   if(startGRU + nGRU - 1  > fileGRU) then; err=20; message=trim(message)//'startGRU + nGRU is larger than then the GRU dimension'; return; end if
 end if
 if (present(checkHRU)) then
   if(checkHRU > fileHRU) then; err=20; message=trim(message)//'checkHRU is larger than then the HRU dimension'; return; end if
 end if

 ! *********************************************************************************************
 ! read mapping vectors and populate mapping structures
 ! **********************************************************************************************
 ! allocate space for indices and types
 allocate(gru_id(fileGRU),nGlac(fileGRU),nWtld(fileGRU),hru_ix(fileHRU),hru_id(fileHRU),hru2gru_id(fileHRU))

 ! read gru_id from netcdf file
 err = nf90_inq_varid(ncid,"gruId",varID);     if (err/=nf90_noerr) then; message=trim(message)//'problem finding gruId'; return; end if
 err = nf90_get_var(ncid,varID,gru_id);        if (err/=nf90_noerr) then; message=trim(message)//'problem reading gruId'; return; end if

 ! read hru_id from netcdf file
 err = nf90_inq_varid(ncid,"hruId",varID);     if (err/=nf90_noerr) then; message=trim(message)//'problem finding hruId'; return; end if
 err = nf90_get_var(ncid,varID,hru_id);        if (err/=nf90_noerr) then; message=trim(message)//'problem reading hruId'; return; end if

 ! read hru2gru_id from netcdf file
 err = nf90_inq_varid(ncid,"hru2gruId",varID); if (err/=nf90_noerr) then; message=trim(message)//'problem finding hru2gruId'; return; end if
 err = nf90_get_var(ncid,varID,hru2gru_id);    if (err/=nf90_noerr) then; message=trim(message)//'problem reading hru2gruId'; return; end if

 ! read domain information from netcdf file
 err = nf90_inq_varid(ncid,"nGlac",varID)
 if (err/=nf90_noerr) then
   nGlac = 0 ! backwards compatibility
 else
   err = nf90_get_var(ncid,varID,nGlac);       if (err/=nf90_noerr) then; message=trim(message)//'problem reading nGlac'; return; end if
 end if
 err = nf90_inq_varid(ncid,"nWtld",varID) 
 if (err/=nf90_noerr) then
   nWtld = 0 ! backwards compatibility
 else
   err = nf90_get_var(ncid,varID,nWtld);       if (err/=nf90_noerr) then; message=trim(message)//'problem reading nWtld'; return; end if
 end if

 ! array from 1 to total # of HRUs in attributes file
 hru_ix=arth(1,1,fileHRU)

! check that the mappings are not alreaday allocated
#ifdef NGEN_ACTIVE
if (allocated(gru_struc)) deallocate(gru_struc) ! free existing mapping possibly from previous GRU run
if (allocated(index_map)) deallocate(index_map) ! free existing mapping possibly from previous GRU run
#else
if (allocated(gru_struc)) then; err=20; message=trim(message)//'gru_struc is unexpectedly allocated'; return; end if
if (allocated(index_map)) then; err=20; message=trim(message)//'index_map is unexpectedly allocated'; return; end if
#endif

 ! allocate first level of gru to hru mapping
 allocate(gru_struc(nGRU))

 ! set gru to hru mapping
 if (present(checkHRU)) then                                  ! allocate space for single-HRU run
   ! gru to hru mapping
   iGRU = 1
   gru_struc(iGRU)%hruCount             = 1                    ! number of HRUs in each GRU
   gru_struc(iGRU)%gru_id               = hru2gru_id(checkHRU) ! set gru id
   gru_struc(iGRU)%gru_nc               = sGRU                 ! set gru index within the netcdf file
   allocate(gru_struc(iGRU)%hruInfo(gru_struc(iGRU)%hruCount)) ! allocate second level of gru to hru map
   gru_struc(iGRU)%hruInfo(iGRU)%hru_nc = checkHRU             ! set hru id in attributes netcdf file
   gru_struc(iGRU)%hruInfo(iGRU)%hru_ix = 1                    ! set index of hru in run space
   gru_struc(iGRU)%hruInfo(iGRU)%hru_id = hru_id(checkHRU)     ! set id of hru
   gru_struc(iGRU)%nGlac                = nGlac(sGRU)          ! set number of glaciers in the gru
   gru_struc(iGRU)%nWtld                = nWtld(sGRU)          ! set number of wetlands in the gru

 else ! allocate space for anything except a single HRU run
   iHRU = 1
   do iGRU = 1,nGRU
     if (count(hru2gru_Id == gru_id(iGRU+sGRU-1)) < 1) then; err=20; message=trim(message)//'problem finding HRUs belonging to GRU'; return; end if
     gru_struc(iGRU)%hruCount = count(hru2gru_Id == gru_id(iGRU+sGRU-1))          ! number of HRUs in each GRU
#ifdef NGEN_ACTIVE
     if (gru_struc(iGRU)%hruCount > 1) then; err=20; message=trim(message)//'NGEN currently only supports single-HRU per GRU'; return; end if
     print *, 'GRU id is ', gru_id(iGRU+sGRU-1)
#endif
     gru_struc(iGRU)%gru_id            = gru_id(iGRU+sGRU-1)                               ! set gru id
     gru_struc(iGRU)%gru_nc            = iGRU+sGRU-1                                       ! set gru index in the netcdf file
     allocate(gru_struc(iGRU)%hruInfo(gru_struc(iGRU)%hruCount))                           ! allocate second level of gru to hru map
     gru_struc(iGRU)%hruInfo(:)%hru_nc = pack(hru_ix,hru2gru_id == gru_struc(iGRU)%gru_id) ! set hru id in attributes netcdf file
     gru_struc(iGRU)%hruInfo(:)%hru_ix = arth(iHRU,1,gru_struc(iGRU)%hruCount)             ! set index of hru in run space
     gru_struc(iGRU)%hruInfo(:)%hru_id = hru_id(gru_struc(iGRU)%hruInfo(:)%hru_nc)         ! set id of hru
     gru_struc(iGRU)%nGlac             = nGlac(iGRU+sGRU-1)                                ! set number of glaciers in the gru
     gru_struc(iGRU)%nWtld             = nWtld(iGRU+sGRU-1)                                ! set number of wetlands in the gru
     iHRU = iHRU + gru_struc(iGRU)%hruCount
   enddo ! iGRU = 1,nGRU
 end if ! not checkHRU

 ! set hru to gru mapping
 nHRU = sum(gru_struc%hruCount)  
 allocate(index_map(nHRU))                                                                      ! allocate first level of hru to gru mapping

 if (present(checkHRU)) then                                                                    ! allocate space for single-HRU run
   if (nHRU/=1) then; err=-20; message=trim(message)//'wrong # of HRUs for checkHRU run'; return; end if
   iGRU = 1;
   index_map(1)%gru_ix      = iGRU                                                               ! index of gru in run space to which the hru belongs
   index_map(1)%localHRU_ix = hru_ix(1)                                                          ! index of hru within the gru

 else ! anything other than a single HRU run
   do iGRU = 1,nGRU
     index_map(gru_struc(iGRU)%hruInfo(:)%hru_ix)%gru_ix      = iGRU                              ! index of gru in run space to which the hru belongs
     index_map(gru_struc(iGRU)%hruInfo(:)%hru_ix)%localHRU_ix = hru_ix(1:gru_struc(iGRU)%hruCount)! index of hru within the gru
   enddo ! iGRU = 1,nGRU

 end if ! not checkHRU

 ! if necessary, do grid dimensions
 if(sum(gru_struc(1:nGRU)%nGlac)>0)then
   call read_dimensionGrid(ncid,fileGRU,nGRU,err,message); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 else ! no grids
   do iGRU = 1, nGRU
     gru_struc(iGRU)%nGrid = 0
     allocate(gru_struc(iGRU)%gridInfo(0))
   enddo
 endif

 deallocate(gru_id, nGlac, nWtld, hru_ix, hru_id, hru2gru_id)
 ! close netcdf file
 call nc_file_close(ncid,err,cmessage)
 if (err/=nf90_noerr) then; message=trim(message)//trim(cmessage); return; end if

end subroutine read_dimension


! ************************************************************************************************
! private subroutine read_dimensionGrid: read grid dimension information on local attributes
! ************************************************************************************************
subroutine read_dimensionGrid(ncid,fileGRU,nGRU,err,message)
  USE netcdf
  USE netcdf_util_module,only:netcdf_err                     ! netcdf error handling function
  USE var_lookup,only:iLookGRID                              ! named variables for the glacier grid information

  ! io vars
  integer(i4b),intent(in)              :: ncid                    ! netcdf file ID
  integer(i4b),intent(in)              :: fileGRU                 ! number of GRUs in the input file
  integer(i4b),intent(in)              :: nGRU                    ! number of grouped response units
  integer(i4b),intent(out)             :: err                     ! error code
  character(*),intent(out)             :: message                 ! error message
  ! define local variables
  integer(i4b)                         :: iGRU                    ! loop indices
  integer(i4b)                         :: varID                   ! NetCDF variable ID
  integer(i4b)                         :: dimID                   ! netcdf file dimension id
  integer(i4b)                         :: filegrid                ! number of grids in the input file
  integer(i4b)                         :: filexgrid               ! number of x grid points in the input file
  integer(i4b)                         :: fileygrid               ! number of y grid points in the input file
  integer(i4b),allocatable             :: grid_id(:,:)            ! read glacier grid IDs in from attributes file
  real(rkind),allocatable              :: dx(:,:),dy(:,:)         ! glacier grid information (spacing)
  integer(i4b),allocatable             :: nx(:,:),ny(:,:)         ! glacier grid information (size)
  integer(i4b),allocatable             :: nGrid(:)                ! number of grids in a GRU
 
  ! Start procedure here
  err=0; message="read_dimensionGrid/"
 
  ! *********************************************************************************************
  ! read and set grid dimensions
  ! **********************************************************************************************
 
  ! get grid dimension of whole file
  err = nf90_inq_dimid(ncid,"grid",dimID);                   if(err/=nf90_noerr)then; message=trim(message)//'problem finding grid dimension/'//trim(nf90_strerror(err)); return; end if
  err = nf90_inquire_dimension(ncid, dimID, len = filegrid); if(err/=nf90_noerr)then; message=trim(message)//'problem reading grid dimension/'//trim(nf90_strerror(err)); return; end if

  ! get ygrid dimension of whole file
  err = nf90_inq_dimid(ncid,"ygrid",dimID);                   if(err/=nf90_noerr)then; message=trim(message)//'problem finding ygrid dimension/'//trim(nf90_strerror(err)); return; end if
  err = nf90_inquire_dimension(ncid, dimID, len = fileygrid); if(err/=nf90_noerr)then; message=trim(message)//'problem reading ygrid dimension/'//trim(nf90_strerror(err)); return; end if

  ! get xgrid dimension of whole file
  err = nf90_inq_dimid(ncid,"xgrid",dimID);                   if(err/=nf90_noerr)then; message=trim(message)//'problem finding xgrid dimension/'//trim(nf90_strerror(err)); return; end if
  err = nf90_inquire_dimension(ncid, dimID, len = filexgrid); if(err/=nf90_noerr)then; message=trim(message)//'problem reading xgrid dimension/'//trim(nf90_strerror(err)); return; end if
  
  ! *********************************************************************************************
  ! read glacier grid information and populate structures
  ! **********************************************************************************************
  ! allocate space for indices
  allocate(grid_id(fileGRU,filegrid),nGrid(fileGRU),dx(fileGRU,filegrid),dy(fileGRU,filegrid),nx(fileGRU,filegrid),ny(fileGRU,filegrid))
 
  ! read grid_id from netcdf file
  err = nf90_inq_varid(ncid,"gridId",varID); if (err/=nf90_noerr) then; message=trim(message)//'problem finding gridId'; return; end if
  err = nf90_get_var(ncid,varID,grid_id);    if (err/=nf90_noerr) then; message=trim(message)//'problem reading gridId'; return; end if

  ! read nGrid from netcdf file
  err = nf90_inq_varid(ncid,"nGrid",varID);  if (err/=nf90_noerr) then; message=trim(message)//'problem finding nGrid'; return; end if
  err = nf90_get_var(ncid,varID,nGrid);      if (err/=nf90_noerr) then; message=trim(message)//'problem reading nGrid'; return; end if
 
  ! read dx from netcdf file
  err = nf90_inq_varid(ncid,"dx",varID);     if (err/=nf90_noerr) then; message=trim(message)//'problem finding dx'; return; end if
  err = nf90_get_var(ncid,varID,dx);         if (err/=nf90_noerr) then; message=trim(message)//'problem reading dx'; return; end if
 
  ! read dy from netcdf file
  err = nf90_inq_varid(ncid,"dy",varID);     if (err/=nf90_noerr) then; message=trim(message)//'problem finding dy'; return; end if
  err = nf90_get_var(ncid,varID,dy);         if (err/=nf90_noerr) then; message=trim(message)//'problem reading dy'; return; end if
 
  ! read nx from netcdf file
  err = nf90_inq_varid(ncid,"nx",varID);     if (err/=nf90_noerr) then; message=trim(message)//'problem finding nx'; return; end if
  err = nf90_get_var(ncid,varID,nx);         if (err/=nf90_noerr) then; message=trim(message)//'problem reading nx'; return; end if
 
  ! read ny from netcdf file
  err = nf90_inq_varid(ncid,"ny",varID);     if (err/=nf90_noerr) then; message=trim(message)//'problem finding ny'; return; end if
  err = nf90_get_var(ncid,varID,ny);         if (err/=nf90_noerr) then; message=trim(message)//'problem reading ny'; return; end if

  ! Main loop to set grid information
  do iGRU = 1, nGRU
    gru_struc(iGRU)%nGrid = nGrid(gru_struc(iGRU)%gru_nc)
    ! consistency check, require there is a grid for every glacier
    if (gru_struc(iGRU)%nGlac>gru_struc(iGRU)%nGrid) then; err=20; message=trim(message)//'needs to be as many grids as glaciers'; return; end if
    allocate(gru_struc(iGRU)%gridInfo(gru_struc(iGRU)%nGrid))
    if (gru_struc(iGRU)%nGrid > 0) then
      gru_struc(iGRU)%gridInfo(:)%grid_id = grid_id(gru_struc(iGRU)%gru_nc,1:gru_struc(iGRU)%nGrid)  ! set grid information (id)
      gru_struc(iGRU)%gridInfo(:)%dx      = dx(gru_struc(iGRU)%gru_nc,1:gru_struc(iGRU)%nGrid)       ! set grid information (spacing)
      gru_struc(iGRU)%gridInfo(:)%dy      = dy(gru_struc(iGRU)%gru_nc,1:gru_struc(iGRU)%nGrid)       ! set grid information (spacing)
      gru_struc(iGRU)%gridInfo(:)%nx      = nx(gru_struc(iGRU)%gru_nc,1:gru_struc(iGRU)%nGrid)       ! set grid information (size)
      gru_struc(iGRU)%gridInfo(:)%ny      = ny(gru_struc(iGRU)%gru_nc,1:gru_struc(iGRU)%nGrid)       ! set grid information (size)
    endif ! if has grid
  end do ! gru loop

  deallocate(grid_id,nGrid,dx,dy,nx,ny)

end subroutine read_dimensionGrid
  

! ************************************************************************************************
! public subroutine read_attrb: read information on local attributes
! ************************************************************************************************
subroutine read_attrb(attrFile,nGRU,attrStruct,typeStruct,idStruct,gridStruct,err,message)
 ! provide access to subroutines
 USE netcdf
 USE netcdf_util_module,only:nc_file_open                   ! open netcdf file
 USE netcdf_util_module,only:nc_file_close                  ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                     ! netcdf error handling function
 ! provide access to derived data types
 USE data_types,only:gru_hru_int                            ! x%gru(:)%hru(:)%var(:)     (i4b)
 USE data_types,only:gru_hru_int8                           ! x%gru(:)%hru(:)%var(:)     (i8b)
 USE data_types,only:gru_hru_double                         ! x%gru(:)%hru(:)%var(:)     (rkind)
 USE data_types,only:gru_grid_double                        ! x%gru(:)%grid(:)%var(:)%dat2(:,:)     (rkind)
 USE get_ixname_module,only:get_ixAttr,get_ixType,get_ixId  ! access function to find index of elements in structure
 implicit none

 ! io vars
 character(*)                         :: attrFile           ! input filename
 integer(i4b),intent(in)              :: nGRU               ! number of grouped response units
 type(gru_hru_double),intent(inout)   :: attrStruct         ! local attributes for each HRU
 type(gru_hru_int),intent(inout)      :: typeStruct         ! local classification of soil veg etc. for each HRU
 type(gru_hru_int8),intent(inout)     :: idStruct           ! local values of hru and gru IDs
 type(gru_grid_double),intent(inout)  :: gridStruct         ! grid information
 integer(i4b),intent(out)             :: err                ! error code
 character(*),intent(out)             :: message            ! error message
 ! define local variables
 character(len=256)                   :: cmessage           ! error message for downwind routine
 integer(i4b)                         :: iVar               ! loop through varibles in the netcdf file
 integer(i4b)                         :: iHRU               ! index of an HRU within a GRU
 integer(i4b)                         :: iGRU               ! index of an GRU
 integer(i4b)                         :: varType            ! type of variable (categorica, numerical, idrelated)
 integer(i4b)                         :: varIndx            ! index of variable within its data structure
 ! check structures
 integer(i4b)                         :: iCheck             ! index of an attribute name
 logical(lgt),allocatable             :: checkType(:)       ! vector to check if we have all desired categorical values
 logical(lgt),allocatable             :: checkId(:)         ! vector to check if we have all desired IDs
 logical(lgt),allocatable             :: checkAttr(:)       ! vector to check if we have all desired local attributes
 ! netcdf variables
 integer(i4b)                         :: ncid               ! netcdf file id
 character(LEN=nf90_max_name)         :: varName            ! character array of netcdf variable name
 integer(i4b)                         :: nVar               ! number of variables in netcdf local attribute file
 integer(i4b),parameter               :: categorical=101    ! named variable to denote categorical data
 integer(i4b),parameter               :: numerical=102      ! named variable to denote numerical data
 integer(i4b),parameter               :: idrelated=103      ! named variable to denote ID related data
 integer(i4b)                         :: categorical_var(1) ! temporary categorical variable from local attributes netcdf file
 real(rkind)                          :: numeric_var(1)     ! temporary numeric variable from local attributes netcdf file
 integer(i8b)                         :: idrelated_var(1)   ! temporary ID related variable from local attributes netcdf file

 ! Start procedure here
 err=0; message="read_attrb/"

 ! **********************************************************************************************
 ! (1) prepare check vectors
 ! **********************************************************************************************
 allocate(checkType(size(type_meta)),checkAttr(size(attr_meta)),checkId(size(id_meta)),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for variable check vectors'; return; endif
 checkType(:) = .false.
 checkAttr(:) = .false.
 checkId(:)   = .false.

 ! **********************************************************************************************
 ! (2) open netcdf file
 ! **********************************************************************************************
 ! open file
 call nc_file_open(trim(attrFile),nf90_noWrite,ncid,err,cmessage)
 if(err/=nf90_noerr)then; message=trim(message)//trim(cmessage); return; endif

 ! get number of variables total in netcdf file
 err = nf90_inquire(ncid,nvariables=nVar)
 call netcdf_err(err,message); if (err/=nf90_noerr) return

 ! **********************************************************************************************
 ! (3) read local attributes
 ! **********************************************************************************************
 ! loop through variables in netcdf file and pull out local attributes
 iCheck = 1
 do iVar = 1,nVar

   ! inqure about current variable name, type, number of dimensions
   err = nf90_inquire_variable(ncid,iVar,name=varName)
   if(err/=nf90_noerr)then; message=trim(message)//'problem inquiring variable: '//trim(varName)//'/'//trim(nf90_strerror(err)); return; endif

   ! find attribute name
   select case(trim(varName))

     ! ** categorical data
     case('vegTypeIndex','soilTypeIndex','slopeTypeIndex','downHRUindex')

      ! get the index of the variable
      varType = categorical
      varIndx = get_ixType(varName)
      checkType(varIndx) = .true.

      ! check that the variable could be identified in the data structure
      if(varIndx < 1)then; err=20; message=trim(message)//'unable to find variable ['//trim(varName)//'] in data structure'; return; endif

      ! get data from netcdf file and store in vector
      do iGRU=1,nGRU
        do iHRU = 1,gru_struc(iGRU)%hruCount
          err = nf90_get_var(ncid,iVar,categorical_var,start=(/gru_struc(iGRU)%hruInfo(iHRU)%hru_nc/),count=(/1/))
          if(err/=nf90_noerr)then; message=trim(message)//'problem reading: '//trim(varName); return; end if
          typeStruct%gru(iGRU)%hru(iHRU)%var(varIndx) = categorical_var(1)
        end do
      end do

     ! ** ID related data
     case('hruId')
      ! get the index of the variable
      varType = idrelated
      varIndx = get_ixId(varName)
      checkId(varIndx) = .true.

      ! check that the variable could be identified in the data structure
      if(varIndx < 1)then; err=20; message=trim(message)//'unable to find variable ['//trim(varName)//'] in data structure'; return; endif

      ! get data from netcdf file and store in vector
      do iGRU=1,nGRU
        do iHRU = 1,gru_struc(iGRU)%hruCount
          err = nf90_get_var(ncid,iVar,idrelated_var,start=(/gru_struc(iGRU)%hruInfo(iHRU)%hru_nc/),count=(/1/))
          if(err/=nf90_noerr)then; message=trim(message)//'problem reading: '//trim(varName); return; end if
          idStruct%gru(iGRU)%hru(iHRU)%var(varIndx) = idrelated_var(1)
        end do
      end do

     ! ** numerical data
     case('latitude','longitude','elevation','tan_slope','contourLength','HRUarea','mHeight','aspect')

      ! get the index of the variable
      varType = numerical
      varIndx = get_ixAttr(varName)
      checkAttr(varIndx) = .true.

      ! check that the variable could be identified in the data structure
      if(varIndx < 1)then; err=20; message=trim(message)//'unable to find variable ['//trim(varName)//'] in data structure'; return; endif

      ! get data from netcdf file and store in vector
      do iGRU=1,nGRU
        do iHRU = 1, gru_struc(iGRU)%hruCount
          err = nf90_get_var(ncid,iVar,numeric_var,start=(/gru_struc(iGRU)%hruInfo(iHRU)%hru_nc/),count=(/1/))
          if(err/=nf90_noerr)then; message=trim(message)//'problem reading: '//trim(varName); return; end if
          attrStruct%gru(iGRU)%hru(iHRU)%var(varIndx) = numeric_var(1)
        end do
      end do

     ! for GRU domain quantity and grid variables, do nothing (information read above in read_dimension or read_dimensionGrid)
     case('nGlac','nWtld','nGrid','gridId','nx','ny','dx','dy'); cycle

     ! for mapping variables, do nothing (information read above in read_dimension)   
     case('hru2gruId','gruId')
      ! get the index of the variable
      varType = idrelated
      varIndx = get_ixId(varName)
      checkId(varIndx) = .true.

     ! for glacier grid variables, will do later in read_attrbGlac call
     case('cell2hruId','bed_elev','glacierMask'); cycle

     ! check that variables are what we expect
     case default; message=trim(message)//'unknown variable ['//trim(varName)//'] in local attributes file'; err=20; return

   end select ! select variable

 end do ! (looping through netcdf local attribute file)
 
 ! ** now handle the optional aspect variable if it's missing
 varIndx = get_ixAttr('aspect')
 ! check that the variable was not found in the attribute file
 if(.not. checkAttr(varIndx)) then
   write(*,*) NEW_LINE('A')//'INFO: aspect not found in the input attribute file, continuing ...'//NEW_LINE('A')

   do iGRU=1,nGRU
     do iHRU = 1, gru_struc(iGRU)%hruCount
       attrStruct%gru(iGRU)%hru(iHRU)%var(varIndx) = nr_realMissing      ! populate variable with out-of-range value, used later
     end do
   end do
   checkAttr(varIndx) = .true.
 endif

 ! **********************************************************************************************
 ! (4) check that we have all the desired variables
 ! **********************************************************************************************
 ! check that we have all desired categorical variables
 if(any(.not.checkType))then
   do iCheck = 1,size(type_meta)
     if(.not.checkType(iCheck))then; err=20; message=trim(message)//'missing variable ['//trim(type_meta(iCheck)%varName)//'] in local attributes file'; return; endif
   end do
 endif

 ! check that we have all desired ID variables
 if(any(.not.checkId))then
   do iCheck = 1,size(id_meta)
     if(.not.checkId(iCheck))then; err=20; message=trim(message)//'missing variable ['//trim(id_meta(iCheck)%varName)//'] in local attributes file'; return; endif
   end do
 endif

 ! check that we have all desired local attributes
 if(any(.not.checkAttr))then
   do iCheck = 1,size(attr_meta)
     if(.not.checkAttr(iCheck))then; err=20; message=trim(message)//'missing variable ['//trim(attr_meta(iCheck)%varName)//'] in local attributes file'; return; endif
   end do
 endif

 ! **********************************************************************************************
 ! (5) get glacier grid variables if needed
 ! **********************************************************************************************
 if(sum(gru_struc(1:nGRU)%nGlac)>0)then
   call read_attrbGlac(ncid,nGRU,gridStruct,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 endif

 ! **********************************************************************************************
 ! (6) close netcdf file
 ! **********************************************************************************************
! free memory
 deallocate(checkType)
 deallocate(checkId)
 deallocate(checkAttr)

 call nc_file_close(ncid,err,cmessage)
 if (err/=nf90_noerr)then; message=trim(message)//trim(cmessage); return; end if

 end subroutine read_attrb

 ! ************************************************************************************************
 ! public subroutine read_attrbGlac: read information on glacier attributes
 ! ************************************************************************************************
 subroutine read_attrbGlac(ncid,nGRU,gridStruct,err,message)
 ! provide access to subroutines
 USE netcdf
 USE netcdf_util_module,only:netcdf_err                     ! netcdf error handling function
 ! provide access to derived data types
 USE var_lookup,only:iLookGRID                              ! named variables for the glacier grid information
 USE data_types,only:gru_grid_double                        ! x%gru(:)%grid(:)%var(:)%dat2(:,:)     (rkind)
 implicit none

 ! io vars
 integer(i4b),intent(in)              :: ncid                    ! netcdf file ID
 integer(i4b),intent(in)              :: nGRU                    ! number of grouped response units in simulation domain
 type(gru_grid_double),intent(inout)  :: gridStruct              ! grid information
 integer(i4b),intent(out)             :: err                     ! error code
 character(*),intent(out)             :: message                 ! error message
 ! define local variables
 integer(i4b)                         :: iGRU,iGrid,iHRU,i       ! loop indices
 integer(i4b)                         :: iVar                    ! variable index
 integer(i4b)                         :: varID                   ! NetCDF variable ID
 integer(i4b)                         :: nGrid                   ! number of grids in a GRU
 integer(i4b)                         :: nx,ny                   ! number of grid points in a glacier
 integer(i4b),dimension(3)            :: ngdx                    ! intermediate array of loop indices for grid variables
 real(rkind),allocatable              :: grid_var(:,:,:,:)       ! grid variable read from netcdf file
 real(rkind),allocatable              :: iHRUgrid(:,:)           ! index of HRU as a grid
 character(LEN=128)                   :: varName                 ! variable name

 ! Start procedure here
 err=0; message="read_attrbGlac/"

 ! vars to read
 ngdx = (/iLookGRID%bed_elev,iLookGRID%glacierMask,iLookGRID%cell2hru/)

 ! fill structure with data
 do iGRU = 1, nGRU
  nGrid = size(gru_struc(iGRU)%gridInfo(:)%grid_id)
  do iGrid = 1, nGrid
    nx = gru_struc(iGRU)%gridInfo(iGrid)%nx
    ny = gru_struc(iGRU)%gridInfo(iGrid)%ny
    do i = 1, size(ngdx)
      iVar = ngdx(i)
      varName = trim(grid_meta(iVar)%varName)
      if(iVar==iLookGRID%cell2hru) varName = "cell2hruId" ! cell2hru is the HRU index, but the parameter is an hruId
      ! get data from netcdf file
      allocate(grid_var(1,1,nx,ny))
      err = nf90_inq_varid(ncid,varName,varID) 
      if (err/=nf90_noerr) then; message=trim(message)//'problem finding '//varName; return; end if
      err = nf90_get_var(ncid,varID,grid_var,start=(/gru_struc(iGRU)%gru_nc,iGrid,1,1/),count=(/1,1,nx,ny/))
      if (err/=nf90_noerr) then; message=trim(message)//'problem reading '//varName; return; end if
      ! put in structure
      if(iVar==iLookGRID%cell2hru)then
        allocate(iHRUgrid(nx,ny))
        gridStruct%gru(iGRU)%grid(iGrid)%var(iVar)%dat2(1:nx,1:ny) = -1.0_rkind ! initialize with an invalid index
        do iHRU = 1, gru_struc(iGRU)%hruCount
          iHRUgrid = real(iHRU, rkind)
          gridStruct%gru(iGRU)%grid(iGrid)%var(iVar)%dat2(1:nx,1:ny) = merge(iHRUgrid, gridStruct%gru(iGRU)%grid(iGrid)%var(iVar)%dat2(1:nx,1:ny), &
                                                                              grid_var(1,1,1:nx,1:ny) == gru_struc(iGRU)%hruInfo(iHRU)%hru_id)
        end do
        deallocate(iHRUgrid)
      else
        gridStruct%gru(iGRU)%grid(iGrid)%var(iVar)%dat2(1:nx,1:ny)  = grid_var(1,1,1:nx,1:ny)
      endif
      deallocate(grid_var)
    end do ! variable loop
  end do ! grid loop
 end do ! gru loop

end subroutine read_attrbGlac


end module read_attrb_module

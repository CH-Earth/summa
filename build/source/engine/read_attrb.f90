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

USE netcdf
USE netcdf_util_module,only:nc_file_open                   ! open netcdf file
USE netcdf_util_module,only:nc_file_close                  ! close netcdf file
USE nr_utils_module ,only:arth                             ! use to build vectors with regular increments

USE build_options, only : ngen_active

USE globalData, only: integerMissing

implicit none

private

public::read_dimension
public::read_mapping_vectors
public::read_attrb

contains

 ! ************************************************************************************************
 ! public subroutine read_dimension:
 !   read GRU/HRU dimensions from the LocalAttributes file and
 !   determine the GRU range for the current run domain
 ! ************************************************************************************************
 subroutine read_dimension(attrFile,                       & ! LocalAttributes filename
                           nGRU_file, nHRU_file,           & ! dimensions of the complete file
                           startGRU_domain, nGRU_domain,   & ! GRU range for the current run domain
                           err, message,                   & ! error control
                           startGRU_user,                  & ! optional: first GRU specified with CLI -g
                           checkHRU)                         ! optional: HRU file index specified with CLI -h

 implicit none

 ! input
 character(*),intent(in)              :: attrFile           ! LocalAttributes filename

 ! file dimensions
 integer(i4b),intent(out)             :: nGRU_file          ! number of GRUs in the complete file
 integer(i4b),intent(out)             :: nHRU_file          ! number of HRUs in the complete file

 ! run-domain dimensions
 integer(i4b),intent(inout)           :: startGRU_domain    ! file index of first GRU in the run domain
 integer(i4b),intent(inout)           :: nGRU_domain        ! number of GRUs in the run domain

 ! error control
 integer(i4b),intent(out)             :: err                ! error code
 character(*),intent(out)             :: message            ! error message

 ! optional user selections
 integer(i4b),intent(in),optional     :: startGRU_user      ! file index of first GRU specified with CLI -g
 integer(i4b),intent(in),optional     :: checkHRU           ! file index of HRU specified with CLI -h

 ! define variables for NetCDF file operations
 integer(i4b)       :: ncid      ! NetCDF file ID
 integer(i4b)       :: dimID     ! NetCDF dimension ID
 character(len=256) :: cmessage  ! error message from downwind routine

 ! initialize error control
 err=0; message='read_dimension/'

 ! check for conflicting run selections
 if(present(startGRU_user).and.present(checkHRU))then
   message=trim(message)//'startGRU_user and checkHRU both exist, which is not supported'
   err=10; return
 endif

 ! open LocalAttributes file
 call nc_file_open(trim(attrFile),nf90_noWrite,ncid,err,cmessage)
 if(err/=nf90_noerr)then; message=trim(message)//trim(cmessage); return; endif

 ! get GRU dimension of complete file
 err=nf90_inq_dimid(ncid,'gru',dimID); if(err/=nf90_noerr)then; message=trim(message)//'problem finding GRU dimension/'//trim(nf90_strerror(err)); return; endif
 err=nf90_inquire_dimension(ncid,dimID,len=nGRU_file); if(err/=nf90_noerr)then; message=trim(message)//'problem reading GRU dimension/'//trim(nf90_strerror(err)); return; endif

 ! get HRU dimension of complete file
 err=nf90_inq_dimid(ncid,'hru',dimID); if(err/=nf90_noerr)then; message=trim(message)//'problem finding HRU dimension/'//trim(nf90_strerror(err)); return; endif
 err=nf90_inquire_dimension(ncid,dimID,len=nHRU_file); if(err/=nf90_noerr)then; message=trim(message)//'problem reading HRU dimension/'//trim(nf90_strerror(err)); return; endif

 ! define GRU range for the run domain
 if(present(startGRU_user))then
   if(nGRU_domain<1)then; err=20; message=trim(message)//'nGRU_domain < 1 for a -g run'; return; endif
   startGRU_domain=startGRU_user
 elseif(present(checkHRU))then
   startGRU_domain=integerMissing  ! resolved later from the HRU-to-GRU mapping
   nGRU_domain=1
 else
   startGRU_domain=1
   nGRU_domain=nGRU_file
 endif

 ! check requested domain against file dimensions
 
 ! check bounds
 if(present(checkHRU))then
   if(checkHRU<1 .or. checkHRU>nHRU_file)then
     message=trim(message)//'checkHRU is outside the HRU dimension'
     err=20; return
   endif
 endif

 if(present(startGRU_user))then
   if(startGRU_domain+nGRU_domain-1>nGRU_file)then
     err=20; message=trim(message)//'GRU run domain extends beyond the GRU dimension'; return
   endif
 endif

 ! close LocalAttributes file
 call nc_file_close(ncid,err,cmessage)
 if(err/=nf90_noerr)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine read_dimension

 ! -----------------------------------------------------------------------------------------------------------------------
 ! -----------------------------------------------------------------------------------------------------------------------
 ! -----------------------------------------------------------------------------------------------------------------------
 ! -----------------------------------------------------------------------------------------------------------------------

 subroutine read_mapping_vectors(attrFile, nGRU_file, nHRU_file,         &
                                 startGRU_local, nGRU_local, nHRU_local, &
                                 checkHRU, err, message)

 ! provide access to global mapping structures
 USE globalData, only : gru_struc   ! local GRU-to-HRU mapping
 USE globalData, only : index_map   ! local HRU-to-GRU mapping

 implicit none

 character(*), intent(in)   :: attrFile               ! LocalAttributes filename
 integer(i4b), intent(in)   :: nGRU_file              ! number of GRUs in the complete input file
 integer(i4b), intent(in)   :: nHRU_file              ! number of HRUs in the complete input file
 integer(i4b), intent(in)   :: startGRU_local         ! file index of first GRU assigned to this rank
 integer(i4b), intent(in)   :: nGRU_local             ! number of GRUs assigned to this rank
 integer(i4b), intent(out)  :: nHRU_local             ! number of HRUs assigned to this rank
 integer(i4b), intent(in)   :: checkHRU               ! file index of HRU for single-HRU run; missing otherwise
 integer(i4b), intent(out)  :: err                    ! error code
 character(*), intent(out)  :: message                ! error message

 ! local indices and mapping vectors
 integer(i4b)               :: iHRU                   ! HRU counting index in the local run domain
 integer(i4b)               :: iGRU                   ! GRU index in the local run domain
 integer(i4b)               :: iGRU_file              ! GRU index in the LocalAttributes file
 integer(i8b)               :: checkGRU_id            ! GRU ID containing checkHRU
 
 integer(i8b)               :: gru_id(nGRU_local)     ! GRU IDs assigned to this rank
 integer(i8b)               :: gru_id_file(nGRU_file) ! GRU IDs in the complete LocalAttributes file
 integer(i8b)               :: hru_id(nHRU_file)      ! HRU IDs in the complete LocalAttributes file
 integer(i8b)               :: hru2gru_id(nHRU_file)  ! GRU ID associated with each HRU in the file
 integer(i4b)               :: hru_ix(nHRU_file)      ! HRU file indices used for mapping searches

 ! define variables for NetCDF file operation
 integer(i4b)               :: ncid                   ! NetCDF file ID
 integer(i4b)               :: varID_gruID            ! NetCDF variable ID for gruId
 integer(i4b)               :: varID_hruID            ! NetCDF variable ID for hruId
 integer(i4b)               :: varID_hru2gruID        ! NetCDF variable ID for hru2gruId
 character(len=256)         :: cmessage               ! error message from downwind routine

 err=0
 message='read_mapping_vectors/'

 ! *********************************************************************************************
 ! read mapping vectors and populate mapping structures
 ! **********************************************************************************************

 !  NOTE: GRUs assigned to each rank are assumed to form a contiguous block in
 !  the input files. This supports efficient NetCDF hyperslab reads using
 !  startGRU_local and nGRU_local. More general load balancing should reorder
 !  the input files during preprocessing so that each rank retains contiguous I/O.

 ! open LocalAttributes file
 call nc_file_open(trim(attrFile),nf90_noWrite,ncid,err,cmessage)
 if(err/=nf90_noerr)then; message=trim(message)//trim(cmessage); return; endif

 ! get variable IDs
 err=nf90_inq_varid(ncid,'gruId',varID_gruID);        if(err/=nf90_noerr)then; message=trim(message)//'problem finding gruId'; return; endif
 err=nf90_inq_varid(ncid,'hruId',varID_hruID);        if(err/=nf90_noerr)then; message=trim(message)//'problem finding hruId'; return; endif
 err=nf90_inq_varid(ncid,'hru2gruId',varID_hru2gruID);if(err/=nf90_noerr)then; message=trim(message)//'problem finding hru2gruId'; return; endif

 ! read GRU IDs for the local GRU block when its file range is known
 if(startGRU_local /= integerMissing)then ! skip case where checkHRU/=integerMissing
   err=nf90_get_var(ncid,varID_gruID,gru_id,start=[startGRU_local],count=[nGRU_local])
   if(err/=nf90_noerr)then; message=trim(message)//'problem reading gruId'; return; endif
 endif

 ! read complete HRU mapping vectors
 err=nf90_get_var(ncid,varID_hruID,hru_id);        if(err/=nf90_noerr)then; message=trim(message)//'problem reading hruId'; return; endif
 err=nf90_get_var(ncid,varID_hru2gruID,hru2gru_id);if(err/=nf90_noerr)then; message=trim(message)//'problem reading hru2gruId'; return; endif
 
 ! HRU indices in the LocalAttributes file
 hru_ix=arth(1,1,nHRU_file)

 ! check that mapping structures are not already allocated
 if(NGEN_ACTIVE)then
   if(allocated(gru_struc)) deallocate(gru_struc)
   if(allocated(index_map)) deallocate(index_map)
 else
   if(allocated(gru_struc))then; err=20; message=trim(message)//'gru_struc is unexpectedly allocated'; return; endif
   if(allocated(index_map))then; err=20; message=trim(message)//'index_map is unexpectedly allocated'; return; endif
 endif

 ! allocate local GRU-to-HRU mapping
 allocate(gru_struc(nGRU_local))

 ! set local GRU-to-HRU mapping
 if(checkHRU /= integerMissing)then
 
   ! single-HRU run: identify the HRU and containing GRU in the input file
   
   err=nf90_get_var(ncid,varID_gruID,gru_id_file)
   if(err/=nf90_noerr)then; message=trim(message)//'problem reading gruId'; return; endif
   
   checkGRU_id = hru2gru_id(checkHRU)
   
   iGRU_file = findloc(gru_id_file,checkGRU_id,dim=1)
   if(iGRU_file < 1)then
     message=trim(message)//'problem finding GRU containing checkHRU'
     err=20; return
   endif
   
   ! construct local mapping
   iGRU = 1
   gru_struc(iGRU)%hruCount = 1
   gru_struc(iGRU)%gru_id   = checkGRU_id
   gru_struc(iGRU)%gru_nc   = iGRU_file
   
   allocate(gru_struc(iGRU)%hruInfo(1))
   
   gru_struc(iGRU)%hruInfo(1)%hru_nc = checkHRU
   gru_struc(iGRU)%hruInfo(1)%hru_ix = 1
   gru_struc(iGRU)%hruInfo(1)%hru_id = hru_id(checkHRU)

 else
 
   ! GRU-domain run
   iHRU = 1
 
   do iGRU = 1,nGRU_local
 
     gru_struc(iGRU)%hruCount = count(hru2gru_id == gru_id(iGRU))
     if(gru_struc(iGRU)%hruCount < 1)then
       err=20; message=trim(message)//'problem finding HRUs belonging to GRU'; return
     endif
 
     if(NGEN_ACTIVE)then
       if(gru_struc(iGRU)%hruCount > 1)then
         err=20; message=trim(message)//'NGEN currently only supports single-HRU per GRU'; return
       endif
       print *, 'GRU id is ', gru_id(iGRU)
     endif
 
     gru_struc(iGRU)%gru_id = gru_id(iGRU)
     gru_struc(iGRU)%gru_nc = startGRU_local + iGRU - 1
 
     allocate(gru_struc(iGRU)%hruInfo(gru_struc(iGRU)%hruCount))
 
     gru_struc(iGRU)%hruInfo(:)%hru_nc = pack(hru_ix,hru2gru_id == gru_struc(iGRU)%gru_id)
     gru_struc(iGRU)%hruInfo(:)%hru_ix = arth(iHRU,1,gru_struc(iGRU)%hruCount)
     gru_struc(iGRU)%hruInfo(:)%hru_id = hru_id(gru_struc(iGRU)%hruInfo(:)%hru_nc)
 
     iHRU = iHRU + gru_struc(iGRU)%hruCount
 
   end do
 
 endif

 ! set local HRU-to-GRU mapping
 nHRU_local = sum(gru_struc%hruCount)
 allocate(index_map(nHRU_local))
 
 if(checkHRU /= integerMissing)then
 
   ! single-HRU run
   if(nHRU_local /= 1)then
     err=-20; message=trim(message)//'wrong # of HRUs for checkHRU run'; return
   endif
 
   index_map(1)%gru_ix      = 1  ! local GRU index
   index_map(1)%localHRU_ix = 1  ! HRU index within the GRU
 
 else
 
   ! full or GRU-subset run
   do iGRU = 1,nGRU_local
 
     index_map(gru_struc(iGRU)%hruInfo(:)%hru_ix)%gru_ix = iGRU
 
     index_map(gru_struc(iGRU)%hruInfo(:)%hru_ix)%localHRU_ix = &
       arth(1,1,gru_struc(iGRU)%hruCount)
 
   enddo
 
 endif

 ! close netcdf file
 call nc_file_close(ncid,err,cmessage)
 if (err/=nf90_noerr) then; message=trim(message)//trim(cmessage); return; end if

 end subroutine read_mapping_vectors

 ! ************************************************************************************************
 ! public subroutine read_attrb: read information on local attributes
 ! ************************************************************************************************
 subroutine read_attrb(attrFile,nGRU_local,attrStruct,typeStruct,idStruct,upArea,err,message)
 ! subroutines
 USE netcdf
 USE netcdf_util_module,only:nc_file_open                   ! open netcdf file
 USE netcdf_util_module,only:nc_file_close                  ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                     ! netcdf error handling function
 ! derived data types
 USE data_types,only:gru_hru_int                            ! x%gru(:)%hru(:)%var(:)     (i4b)
 USE data_types,only:gru_hru_int8                           ! x%gru(:)%hru(:)%var(:)     (i8b)
 USE data_types,only:gru_hru_double                         ! x%gru(:)%hru(:)%var(:)     (rkind)
 USE data_types,only:gru_d                                  ! x%gru(:)%hru(:)            (rkind)
 ! named variables
 USE var_lookup,only:iLookID                                ! look-up values for local column model ids
 USE var_lookup,only:iLookTYPE                              ! look-up values for classification of veg, soils etc.
 USE var_lookup,only:iLookATTR                              ! look-up values for local attributes
 ! global data
 USE globalData,only:gru_struc                              ! gru-hru mapping structure
 USE globalData,only:attr_meta,type_meta,id_meta            ! metadata structures
 USE get_ixname_module,only:get_ixAttr,get_ixType,get_ixId  ! access function to find index of elements in structure
 implicit none

 ! io vars
 character(*)                         :: attrFile           ! input filename
 integer(i4b),intent(in)              :: nGRU_local         ! number of grouped response units in local rank
 type(gru_hru_double),intent(inout)   :: attrStruct         ! local attributes for each HRU
 type(gru_hru_int),intent(inout)      :: typeStruct         ! local classification of soil veg etc. for each HRU
 type(gru_hru_int8),intent(inout)     :: idStruct           ! local values of hru and gru IDs
 type(gru_d),intent(inout)            :: upArea             ! area upslope of each HRU
 integer(i4b),intent(out)             :: err                ! error code
 character(*),intent(out)             :: message            ! error message
 ! define local variables
 character(len=256)                   :: cmessage           ! error message for downwind routine
 integer(i4b)                         :: iVar               ! loop through varibles in the netcdf file
 integer(i4b)                         :: iHRU,jHRU,kHRU     ! index of an HRU within a GRU
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

 ! define mapping variables

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

      ! check that the variable could be identified in the data structure
      if(varIndx < 1)then; err=20; message=trim(message)//'unable to find variable ['//trim(varName)//'] in data structure'; return; endif
      checkType(varIndx) = .true.

      ! get data from netcdf file and store in vector
      do iGRU=1,nGRU_local
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

      ! check that the variable could be identified in the data structure
      if(varIndx < 1)then; err=20; message=trim(message)//'unable to find variable ['//trim(varName)//'] in data structure'; return; endif
      checkId(varIndx) = .true.

      ! get data from netcdf file and store in vector
      do iGRU=1,nGRU_local
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

      ! check that the variable could be identified in the data structure
      if(varIndx < 1)then; err=20; message=trim(message)//'unable to find variable ['//trim(varName)//'] in data structure'; return; endif
      checkAttr(varIndx) = .true.

      ! get data from netcdf file and store in vector
      do iGRU=1,nGRU_local
        do iHRU = 1, gru_struc(iGRU)%hruCount
          err = nf90_get_var(ncid,iVar,numeric_var,start=(/gru_struc(iGRU)%hruInfo(iHRU)%hru_nc/),count=(/1/))
          if(err/=nf90_noerr)then; message=trim(message)//'problem reading: '//trim(varName); return; end if
          attrStruct%gru(iGRU)%hru(iHRU)%var(varIndx) = numeric_var(1)
        end do
      end do

     ! mapping variables were already read in read_mapping_vectors
     case('hru2gruId','gruId')
      ! get the index of the variable
      varType = idrelated
      varIndx = get_ixId(varName)
      checkId(varIndx) = .true.

     ! check that variables are what we expect
     case default; message=trim(message)//'unknown variable ['//trim(varName)//'] in local attributes file'; err=20; return

   end select ! select variable

 end do ! (looping through netcdf local attribute file)
 
 ! ** now handle the optional aspect variable if it's missing
 varIndx = get_ixAttr('aspect')
 ! check that the variable was not found in the attribute file
 if(.not. checkAttr(varIndx)) then
   write(*,*) NEW_LINE('A')//'INFO: aspect not found in the input attribute file, continuing ...'//NEW_LINE('A')

   do iGRU=1,nGRU_local
     do iHRU = 1, gru_struc(iGRU)%hruCount
       attrStruct%gru(iGRU)%hru(iHRU)%var(varIndx) = nr_realMissing      ! populate variable with out-of-range value, used later
     end do
   end do
   checkAttr(varIndx) = .true.
 endif

 ! **********************************************************************************************
 ! (4) check that we have all the desired varaibles
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
 ! (5) close netcdf file
 ! **********************************************************************************************
! free memory
 deallocate(checkType)
 deallocate(checkId)
 deallocate(checkAttr)

 call nc_file_close(ncid,err,cmessage)
 if (err/=nf90_noerr)then; message=trim(message)//trim(cmessage); return; end if

 ! *****************************************************************************
 ! (6) validate HRU connectivity and compute directly contributing area
 ! *****************************************************************************

 do iGRU=1,nGRU_local

  do iHRU=1,gru_struc(iGRU)%hruCount

    kHRU = 0
    upArea%gru(iGRU)%hru(iHRU) = 0._rkind

    do jHRU=1,gru_struc(iGRU)%hruCount

      ! check whether iHRU drains to jHRU
      if(typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%downHRUindex) == &
         idStruct%gru(iGRU)%hru(jHRU)%var(iLookID%hruId))then

        if(kHRU==0)then
          kHRU = jHRU
        else
          message=trim(message)//'downslope HRU identifier is not unique'
          err=20; return
        endif

      endif

      ! check whether jHRU drains directly to iHRU
      if(typeStruct%gru(iGRU)%hru(jHRU)%var(iLookTYPE%downHRUindex) == &
         idStruct%gru(iGRU)%hru(iHRU)%var(iLookID%hruId))then

        upArea%gru(iGRU)%hru(iHRU) = upArea%gru(iGRU)%hru(iHRU) + &
          attrStruct%gru(iGRU)%hru(jHRU)%var(iLookATTR%HRUarea)

      endif

    enddo ! jHRU

  enddo ! iHRU

 enddo ! iGRU

 end subroutine read_attrb

end module read_attrb_module

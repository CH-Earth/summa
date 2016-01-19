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

module read_attrb_module
USE nrtype
implicit none
private
public::read_attrb
contains

 ! ************************************************************************************************
 ! public subroutine read_attrb: read information on local attributes
 ! ************************************************************************************************
 subroutine read_attrb(err,message)
 ! provide access to subroutines
 USE netcdf
 USE netcdf_util_module,only:nc_file_open          ! open netcdf file
 USE nr_utility_module,only:arth
 USE ascii_util_module,only:file_open              ! open ascii file
 USE ascii_util_module,only:split_line             ! extract the list of variable names from the character string
 USE ascii_util_module,only:get_vlines             ! read a vector of non-comment lines from an ASCII file
 USE allocspace_module,only:alloc_attr             ! module to allocate space for local attributes
 USE allocspace_module,only:alloc_type             ! module to allocate space for categorical data
 ! provide access to data
 USE summaFileManager,only:SETNGS_PATH             ! path for metadata files
 USE summaFileManager,only:LOCAL_ATTRIBUTES        ! file containing information on local attributes
 USE data_struc,only: gru_struc,nGRU,nHRU          ! gru-hru structure
 USE data_struc,only:index_map                     ! index mapping 
 USE data_struc,only:attr_meta,type_meta           ! metadata structures
 USE data_struc,only:attr_gru,type_gru             ! data structures
 USE var_lookup,only:iLookATTR,iLookTYPE           ! named variables for elements of the data structures
 USE get_ixname_module,only:get_ixAttr,get_ixType  ! access function to find index of elements in structure
 implicit none
 ! define output
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! define general variables
 real(dp),parameter                   :: missingDouble=-9999._dp  ! missing data
 integer(i4b),parameter               :: missingInteger=-9999     ! missing data
 character(len=256)                   :: cmessage        ! error message for downwind routine
 character(LEN=256)                   :: infile          ! input filename
 integer(i4b)                         :: ivars            ! loop through varibles in the netcdf file
 ! define local variables
 integer(i4b)                         :: ncid            ! netcdf file id
 integer(i4b)                         :: numDims         ! number of dimensions for netcdf variable
 integer(i4b)                         :: mode            ! netCDF file open mode
 integer(i4b)                         :: var_type        ! type of netcdf variable (e.g. integer, double, float)
 integer(i4b)                         :: varid           ! netcdf variable id
 character(LEN=nf90_max_name)         :: var_name        ! character array of netcdf variable name
 integer(i4b),parameter               :: categorical=101 ! named variable to denote categorical data
 integer(i4b),parameter               :: numerical=102   ! named variable to denote numerical data
 integer(i4b),allocatable             :: varType(:)      ! type of variable (categorical or numerical)
 integer(i4b),allocatable             :: varIndx(:)      ! index of variable within its data structure
 integer(i4b)                         :: nvar            ! number of variables in netcdf local attribute file
 integer(i4b)                         :: iAtt            ! index of an attribute name
 integer(i4b)                         :: iHRU            ! index of an HRU
 integer(i4b)                         :: iGRU            ! index of an GRU
 integer(i4b)                         :: nAtt            ! number of model attributes
 integer(i4b)                         :: nVar_attr       ! number of variables in the model attribute structure
 integer(i4b)                         :: nVar_type       ! number of variables in the model category structure
 logical(lgt),allocatable             :: checkType(:)    ! vector to check if we have all desired categorical values
 logical(lgt),allocatable             :: checkAttr(:)    ! vector to check if we have all desired local attributes
 integer(i4b)                         :: ix_gru             ! index of current gru
 integer(i4b)                         :: ix_hru             ! index of current hru
 integer(i4b)                         :: ix                 ! index of current GRU
 integer(i4b),allocatable             :: gru_ix(:)          ! sequential index of GRUs
 integer(i4b),allocatable             :: hru_ix(:)          ! sequential index of hru over entire domain
 integer(i4b),allocatable             :: hru2gru_ixGRU(:)   ! sequential index of GRU for each HRU
 integer(i4b),allocatable             :: hru2gru_ixHRU(:)   ! squential index of HRUs in each GRUs
 integer(i4b),allocatable             :: hru_id(:)          ! unique id of hru over entire domain
 integer(i4b),allocatable             :: gru_id(:)          ! unique ids of GRUs
 integer(i4b),allocatable             :: hru2gru_id(:)      ! unique GRU ids at each HRUare
 integer(i4b),allocatable             :: categorical_var(:) ! temporary categorical variable from local attributes netcdf file
 real(dp),allocatable                 :: numeric_var(:)     ! temporary numeric variable from local attributes netcdf file
 ! define indicators for different data structures
 integer(i4b),parameter               :: ixType=1           ! indicator for type_gru structure
 integer(i4b),parameter               :: ixAttr=2           ! indicator for attr_gru structure 
 integer(i4b), parameter              :: imiss = -999       ! missing value for variable id

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

 !total number of local attributes
 nAtt = nVar_attr + nVar_type

 ! **********************************************************************************************
 ! (1) open files, etc.
 ! **********************************************************************************************
 ! build filename
 infile = trim(SETNGS_PATH)//trim(LOCAL_ATTRIBUTES)
 ! open file
 mode=nf90_NoWrite
 call nc_file_open(trim(infile), mode, ncid, err, cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for local attributes, both numerical and categorical
 call alloc_attr(nGRU,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 call alloc_type(nGRU,err,cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! allocate space for user input mapping vectors
 allocate(gru_id(nGRU),hru_id(nHRU),hru2gru_id(nHRU),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for user input gru-hru mapping vectors'; return; endif

 !allocate space for local mapping vectors
 allocate(hru2gru_ixGRU(nHRU),hru2gru_ixHRU(nHRU),gru_ix(nGRU),hru_ix(nHRU),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for local mapping vectors'; return; endif


 ! **********************************************************************************************
 ! (2) read mapping vectors and populate mapping structures
 ! **********************************************************************************************

 ! read gru_id from netcdf file
 err = nf90_inq_varid(ncid, "gruId", varid)
 if(err/=0)then; message=trim(message)//'problem inquiring gruId'; return; endif
 err = nf90_get_var(ncid,varid,gru_id)
 if(err/=0)then; message=trim(message)//'problem reading gruId'; return; endif
 ! read hruIndex from netcdf file
 err = nf90_inq_varid(ncid, "hruId", varid)
 if(err/=0)then; message=trim(message)//'problem inquiring hruId'; return; endif
 err = nf90_get_var(ncid,varid,hru_id)
 if(err/=0)then; message=trim(message)//'problem reading hruId'; return; endif
 ! read hru2gru_id from netcdf file
 err = nf90_inq_varid(ncid, "hru2gruId", varid)
 if(err/=0)then; message=trim(message)//'problem inquiring hru2gruId'; return; endif
 err = nf90_get_var(ncid,varid,hru2gru_id)
 if(err/=0)then; message=trim(message)//'problem reading hru2gruId'; return; endif

 ! defining the indexes for the GRUs and HRUs
 gru_ix = arth(1,1,nGRU)
 hru_ix = arth(1,1,nHRU)

 ! create local mapping vectors
 do iGRU=1,nGRU
  ix=gru_ix(iGRU)
  where(hru2gru_id == gru_id(iGRU))
    hru2gru_ixGRU = ix
    hru2gru_ixHRU = arth(1,1,gru_struc(iGRU)%hruCount)
  endwhere
 enddo

 ! fill global mapping structures
 !gru_struc contains the HRU index for HRUs in GRUs and the unique HRU id for each HRU in each GRU
 !index_map contains the GRU index associated with each HRU and the global HRU index (1-nHRU)
 if(associated(gru_struc) .and. associated (index_map) )then
  do iHRU=1, nHRU
   gru_struc(hru2gru_ixGRU(iHRU))%hru(hru2gru_ixHRU(iHRU))%hru_ix = hru_ix(iHRU)   ! HRU index
   gru_struc(hru2gru_ixGRU(iHRU))%hru(hru2gru_ixHRU(iHRU))%hru_id = hru_id(iHRU)   ! HRU id

   index_map(hru2gru_ixHRU(iHRU))%gru_ix = hru2gru_ixGRU(iHRU)  !GRU index for each HRU
   index_map(hru2gru_ixHRU(iHRU))%ihru   = hru2gru_ixHRU(iHRU)  !HRU index for each HRU
  enddo
 else
  err = 40
  if(err/=0)then; message=trim(message)//'gru_struc and/or index_map not associated'; return; endif
 end if

 ! **********************************************************************************************
 ! (3) read local attributes
 ! **********************************************************************************************

 !run through netcdf file
 !how many variables are there?
 err = nf90_inquire(ncid, nvariables=nvar)
 call netcdf_err(err,message); if (err/=0) return

 allocate(varType(nAtt),varIndx(nAtt), stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the variable type and index'; return; endif
 ! initialize variables as missing
 varType(:) = missingInteger
 varIndx(:) = missingInteger

 !loop through variables in netcdf file and pull out local attributes
 iAtt = 1
 do ivars=1,nvar
  !inqure about current variable name, type, number of dimensions
  err = nf90_inquire_variable(ncid,ivars,name=var_name,xtype=var_type,ndims=numDims)
  if(err/=0)then; message=trim(message)//'problem inquiring variable: '//trim(var_name); return; endif

  ! find attribute name
  select case(trim(var_name))
   ! categorical data
   case('hruId','vegTypeIndex','soilTypeIndex','slopeTypeIndex','downHRUindex')
    varType(iAtt) = categorical
    if(trim(var_name) == "hruId") then
     varIndx(iAtt) = get_ixType("hruIndex")
    else
     varIndx(iAtt) = get_ixType(var_name)
    end if
    checkType(varIndx(iAtt)) = .true.
    ! check that the variable could be identified in the data structure
    if(varIndx(iAtt) < 1)then; err=20; message=trim(message)//'unable to find variable ['//trim(var_name)//'] in data structure'; return; endif
    !allocate space for attribute
    if(allocated(categorical_var))then
      deallocate(categorical_var)
    end if
    allocate(categorical_var(nHRU),stat=err)
    if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for netcdf variable ['//trim(var_name)//']'; return; endif
    !grab variable from netcdf file
    err = nf90_get_var(ncid,ivars,categorical_var)
    if(err/=0)then; message=trim(message)//'problem reading: '//trim(var_name); return; endif
    ! set attribute to GRU & HRU
    do iHRU=1,nHRU
     ix_gru = index_map(iHRU)%gru_ix
     ix_hru = index_map(iHRU)%ihru
     type_gru(ix_gru)%hru(ix_hru)%var(varIndx(iAtt)) = categorical_var(iHRU)
    end do
    iAtt = iAtt + 1

   ! numerical data
   case('latitude','longitude','elevation','tan_slope','contourLength','HRUarea','mHeight')
    varType(iAtt) = numerical
    varIndx(iAtt) = get_ixAttr(var_name)
    checkAttr(varIndx(iAtt)) = .true.
    ! check that the variable could be identified in the data structure
    if(varIndx(iAtt) < 1)then; err=20; message=trim(message)//'unable to find variable ['//trim(var_name)//'] in data structure'; return; endif
    !allocate space for attribute
    if(allocated(numeric_var))then
      deallocate(numeric_var)
    end if
    allocate(numeric_var(nHRU),stat=err)
    if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for netcdf variable ['//trim(var_name)//']'; return; endif
    !grab variable from netcdf file
    err = nf90_get_var(ncid,ivars,numeric_var)
    if(err/=0)then; message=trim(message)//'problem reading: '//trim(var_name); return; endif
    ! set attribute to GRU & HRU
    do iHRU=1,nHRU
     ix_gru = index_map(iHRU)%gru_ix
     ix_hru = index_map(iHRU)%ihru
     attr_gru(ix_gru)%hru(ix_hru)%var(varIndx(iAtt)) = numeric_var(iHRU)
    end do
    iAtt = iAtt + 1
   !for mapping varibles, do nothing
   case('hru2gruId','gruId')
   ! check that variables are what we expect
   case default
    message=trim(message)//'unknown variable ['//trim(var_name)//'] in local attributes file'
    err=20; return
  end select
 end do! (looping through netcdf local attribute file)

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

! test
! do iGRU=1,nGRU
! do iHRU=1,nHRU
!  print*, '*****'
!  print*, 'hruIndex       = ', type_gru(iGRU)%hru(iHRU)%var(iLookTYPE%hruIndex)
!  print*, 'latitude       = ', attr_gru(iGRU)%hru(iHRU)%var(iLookATTR%latitude)
!  print*, 'longitude      = ', attr_gru(iGRU)%hru(iHRU)%var(iLookATTR%longitude)
!  print*, 'elevation      = ', attr_gru(iGRU)%hru(iHRU)%var(iLookATTR%elevation)
!  print*, 'mHeight        = ', attr_gru(iGRU)%hru(iHRU)%var(iLookATTR%mHeight)
!  print*, 'vegTypeIndex   = ', type_gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex)
!  print*, 'soilTypeIndex  = ', type_gru(iGRU)%hru(iHRU)%var(iLookTYPE%soilTypeIndex)
!  print*, 'slopeTypeIndex = ', type_gru(iGRU)%hru(iHRU)%var(iLookTYPE%slopeTypeIndex)
! end do ! (looping through HRUs)
! end do ! (looping through GRUs)
 !pause

 ! **********************************************************************************************
 ! (5) deallocate space
 ! **********************************************************************************************
 deallocate(varType,varIndx,checkType,checkAttr, stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem deallocating space'; return; endif

 end subroutine read_attrb


 ! **********************************************************************************************************
 ! private subroutine netcdf_err: error control
 ! **********************************************************************************************************
 subroutine netcdf_err(err,message)
 ! used to handle errors for NetCDF calls
 use netcdf
 implicit none
 ! declare dummies
 integer(i4b), intent(inout)   :: err
 character(*), intent(inout)   :: message
 ! start procedure here
 if (err/=nf90_noerr) then
  message=trim(message)//"["//trim(nf90_strerror(err))//"]"
  err=200
 endif
 end subroutine netcdf_err


end module read_attrb_module

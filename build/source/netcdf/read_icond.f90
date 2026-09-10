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

module read_icond_module
USE nr_type
USE netcdf
USE globalData,only:ixHRUfile_min,ixHRUfile_max ! first and last HRUs in the forcing file
USE globalData,only:nTimeDelay        ! number of timesteps in the time delay histogram
USE globalData,only:nSpecBand         ! number of spectral bands
USE globalData,only:nMeltingIceLayers ! number of glacier ice layers that can have a change in total water content
USE globalData,only:thick4area        ! an arbitrary small threshold for glacier thickness to be considered as glacier area

! access domain types
USE globalData,only:upland             ! horizontal domain type for upland areas
USE globalData,only:glacCln1           ! first horizontal domain type for glacier clean areas
USE globalData,only:glacCln2           ! second horizontal domain type for glacier clean areas
USE globalData,only:glacDbr            ! horizontal domain type for glacier debris areas
USE globalData,only:wetland            ! horizontal domain type for wetland areas

implicit none
private
public::read_icond
public::read_icond_nlayers

contains

 ! ************************************************************************************************
 ! public subroutine read_icond_nlayers: read model initial conditions file for number of snow/soil layers
 ! ************************************************************************************************
 subroutine read_icond_nlayers(iconFile,nGRU,nDOM,indx_meta,err,message)
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nr_type
 USE var_lookup,only:iLookINDEX                        ! variable lookup structure
 USE globalData,only:gru_struc                         ! gru-hru mapping structures
 USE globalData,only:startGRU                          ! index of first gru for parallel runs
 USE netcdf_util_module,only:nc_file_close             ! close netcdf file
 USE netcdf_util_module,only:nc_file_open              ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                ! netcdf error handling
 USE data_types,only:var_info                          ! metadata
 implicit none

 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 character(*)  ,intent(in)   :: iconFile            ! name of input (restart) file
 integer(i4b)  ,intent(in)   :: nGRU                ! total # of GRUs in run space
 integer(i4b)  ,intent(inout):: nDOM                ! max number of domains in any HRU
 type(var_info),intent(in)   :: indx_meta(:)        ! metadata
 integer(i4b)  ,intent(out)  :: err                 ! error code
 character(*)  ,intent(out)  :: message             ! returned error message
 ! locals
 integer(i4b)                :: i,j                 ! loop indices
 integer(i4b)                :: ncid                ! netcdf file id
 integer(i4b)                :: dimID               ! netcdf file dimension id
 integer(i4b)                :: ncVarID             ! netcdf variable id
 integer(i4b)                :: fileGRU             ! number of GRUs in netcdf file
 integer(i4b)                :: fileHRU             ! number of HRUs in netcdf file
 integer(i4b)                :: fileDOM             ! number of domains in netcdf file
 integer(i4b)                :: snowID, soilID      ! netcdf variable ids
 integer(i4b)                :: glceID, lakeID      ! netcdf variable ids
 integer(i4b)                :: iGRU, iHRU, iDOM    ! loop indexes
 integer(i4b)                :: iHRU_global         ! index of HRU in the netcdf file
 logical(lgt)                :: no_glceData         ! flag that no ice data in icond
 logical(lgt)                :: no_lakeData         ! flag that no lake data in icond
 logical(lgt)                :: no_dom              ! flag that no domain variable in file
 integer(i4b),allocatable    :: snowData1(:)        ! number of snow layers in all HRUs
 integer(i4b),allocatable    :: soilData1(:)        ! number of soil layers in all HRUs
 integer(i4b),allocatable    :: glceData1(:)        ! number of glacier ice layers in all HRUs
 integer(i4b),allocatable    :: lakeData1(:)        ! number of lake layers in all HRUs
 integer(i4b),allocatable    :: snowData2(:,:)      ! number of snow layers in all HRUs
 integer(i4b),allocatable    :: soilData2(:,:)      ! number of soil layers in all HRUs
 integer(i4b),allocatable    :: glceData2(:,:)      ! number of glacier ice layers in all HRUs
 integer(i4b),allocatable    :: lakeData2(:,:)      ! number of lake layers in all HRUs
 integer(i4b),allocatable    :: dom_type(:,:)       ! read domain type in from netcdf file
 character(len=256)          :: cmessage            ! downstream error message
 integer(i8b),allocatable    :: gru_id(:)           ! GRU id
 integer(i8b),allocatable    :: hru_id(:)           ! HRU id
 integer(i4b),allocatable    :: index_to_gruid(:)   ! mapping from index to gru_id in gru_struc
 integer(i4b),allocatable    :: index_to_hrunc(:,:) ! mapping from index to hru_nc in gru_struc
 logical(lgt)                :: has_gru_id          ! flag for whether the file has gru_id
 logical(lgt)                :: has_hru_id          ! flag for whether the file has hru_id
 ! --------------------------------------------------------------------------------------------------------
 ! initialize error message
 err=0
 message = 'read_icond_nlayers/'
 no_glceData = .false.
 no_lakeData = .false.
 no_dom = .false.
 has_gru_id = .true.
 has_hru_id = .true.

 ! open netcdf file
 call nc_file_open(iconFile,nf90_nowrite,ncid,err,cmessage);
 if (err/=nf90_noerr) then; message=trim(message)//trim(cmessage); return; end if

 ! get number of HRUs in file
 err = nf90_inq_dimid(ncid,"hru",dimId);               if(err/=nf90_noerr)then; message=trim(message)//'problem finding hru dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncid,dimId,len=fileHRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading hru dimension/'//trim(nf90_strerror(err)); return; end if

 ! get number of domains with type in file, if present
 err = nf90_inq_dimid(ncid,"dom",dimId)               
 if(err/=nf90_noerr)then
  fileDOM = 1 ! backwards compatible, just upland domain
  no_dom = .true.
  allocate(dom_type(1,fileHRU))
  dom_type = upland
  err=nf90_noerr    ! reset this err
 else
  err = nf90_inquire_dimension(ncid,dimId,len=fileDOM); if(err/=nf90_noerr)then; message=trim(message)//'problem reading dom dimension/'//trim(nf90_strerror(err)); return; end if
  ! read dom_type from netcdf file
  allocate(dom_type(fileDOM,fileHRU))
  err = nf90_inq_varid(ncid,"domType",ncVarID);  if (err/=nf90_noerr) then; message=trim(message)//'problem finding domType'; return; end if
  err = nf90_get_var(ncid,ncVarID,dom_type);     if (err/=nf90_noerr) then; message=trim(message)//'problem reading domType'; return; end if
 end if
 nDOM = fileDOM

 ! check if the file has the hruId variable
 allocate(hru_id(fileHRU))
 err = nf90_inq_varid(ncid,"hruId",ncVarID)
 if (err/=nf90_noerr)then
   write(*,*) 'WARNING: hruId is not in the initial conditions file ... assuming HRUs in attribute order'
   has_hru_id = .false.
   err=nf90_noerr    ! reset this err
 else
  ! read hru_id from netcdf file
   err = nf90_get_var(ncid,ncVarID,hru_id); if (err/=nf90_noerr) then; message=trim(message)//'problem reading hruId'; return; end if
 end if

 ! check if the file has the GRU dimension
 err = nf90_inq_dimid(ncid,"gru",dimID);    
 if(err/=nf90_noerr)then         
   write(*,*) 'WARNING: GRU is not in the initial conditions file ... assuming GRUs in attribute order'
   has_gru_id = .false.
   err=nf90_noerr    ! reset this err
 else
   err = nf90_inquire_dimension(ncid,dimID,len=fileGRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading gru dimension/'//trim(nf90_strerror(err)); return; end if
   ! read gru_id from netcdf file
   allocate(gru_id(fileGRU))
   err = nf90_inq_varid(ncid,"gruId",ncVarID)
   if (err/=nf90_noerr) then
     write(*,*) 'WARNING: gruId is not in the initial conditions file ... assuming GRUs in attribute order'
     has_gru_id = .false.
     err=nf90_noerr    ! reset this err
   else
     ! read gru_id from netcdf file
     err = nf90_get_var(ncid,ncVarID,gru_id); if (err/=nf90_noerr) then; message=trim(message)//'problem reading gruId'; return; end if
   endif
 end if

 ! Allocate the mapping arrays
 allocate(index_to_gruid(nGRU), index_to_hrunc(nGRU,maxval(gru_struc(:)%hruCount)))

 ! Populate the mapping arrays
 if(has_gru_id .and. has_hru_id)then
   ! if the file has both gru_id and hru_id, use these to populate the mapping arrays
   do iGRU = 1, nGRU
     index_to_gruid(iGRU) = -1  ! Initialize with an invalid index
     do i = 1, fileGRU
       if (gru_struc(iGRU)%gru_id == gru_id(i)) then
         index_to_gruid(iGRU) = i
         do iHRU = 1, gru_struc(iGRU)%hruCount
           index_to_hrunc(iGRU,iHRU) = -1 
           do j = 1, fileHRU
             if (gru_struc(iGRU)%hruInfo(iHRU)%hru_id == hru_id(j)) then
               index_to_hrunc(iGRU,iHRU) = j
               exit
             endif
           end do
         end do ! HRU id loop
         exit
       endif
     end do
   end do ! GRU id loop
 else
   ! assume that the order of the HRUs in the file matches the order of the HRUs in the model attributes
   do iGRU = 1, nGRU
     index_to_gruid(iGRU) = iGRU + startGRU - 1
     do iHRU = 1, gru_struc(iGRU)%hruCount
       index_to_hrunc(iGRU,iHRU) = gru_struc(iGRU)%hruInfo(iHRU)%hru_nc
     end do
   end do
 endif

 ! allocate storage for reading from file (allocate entire file size, even when doing subdomain run)
 allocate(snowData1(fileHRU),snowData2(fileDOM,fileHRU))
 allocate(soilData1(fileHRU),soilData2(fileDOM,fileHRU))
 allocate(glceData1(fileHRU),glceData2(fileDOM,fileHRU))
 allocate(lakeData1(fileHRU),lakeData2(fileDOM,fileHRU))
 snowData1 = 0
 soilData1 = 0
 glceData1 = 0
 lakeData1 = 0
 snowData2 = 0
 soilData2 = 0
 glceData2 = 0
 lakeData2 = 0

 ! count domains and set domain type, DOM not in attributes so don't have to worry about order of dom_type
 ! NOTE: dom_type 0 will be no domain
 do iGRU = 1,nGRU
   do iHRU = 1,gru_struc(iGRU)%hruCount
     iHRU_global = gru_struc(iGRU)%hruInfo(iHRU)%hru_nc
     gru_struc(iGRU)%hruInfo(iHRU)%domCount = 1                                              ! upland domain always present, for changing size glaciers and lakes
     if (any(dom_type(1:fileDOM,iHRU_global)==glacCln1)) &
       gru_struc(iGRU)%hruInfo(iHRU)%domCount = gru_struc(iGRU)%hruInfo(iHRU)%domCount + 1   ! glacier clean domain 1 possible
     if (any(dom_type(1:fileDOM,iHRU_global)==glacCln2)) &
       gru_struc(iGRU)%hruInfo(iHRU)%domCount = gru_struc(iGRU)%hruInfo(iHRU)%domCount + 1   ! glacier clean domain 2 possible
      if (any(dom_type(1:fileDOM,iHRU_global)==glacDbr)) &
       gru_struc(iGRU)%hruInfo(iHRU)%domCount = gru_struc(iGRU)%hruInfo(iHRU)%domCount + 1   ! glacier debris domain possible
     if (any(dom_type(1:fileDOM,iHRU_global)==wetland)) &
       gru_struc(iGRU)%hruInfo(iHRU)%domCount = gru_struc(iGRU)%hruInfo(iHRU)%domCount + 1   ! wetland domain possible
     allocate(gru_struc(iGRU)%hruInfo(iHRU)%domInfo(gru_struc(iGRU)%hruInfo(iHRU)%domCount)) ! allocate third level of gru to hru map
     gru_struc(iGRU)%hruInfo(iHRU)%domInfo(:)%dom_type = dom_type(1:gru_struc(iGRU)%hruInfo(iHRU)%domCount,iHRU_global)
   enddo
 enddo

 ! get netcdf ids for the variables holding number of layers in each domain or hru
 err = nf90_inq_varid(ncid,trim(indx_meta(iLookINDEX%nSnow)%varName),snowID); call netcdf_err(err,message)
 err = nf90_inq_varid(ncid,trim(indx_meta(iLookINDEX%nLake)%varName),lakeID)
 if(err/=nf90_noerr ) no_lakeData = .true.
 err = nf90_inq_varid(ncid,trim(indx_meta(iLookINDEX%nSoil)%varName),soilID); call netcdf_err(err,message)
 err = nf90_inq_varid(ncid,trim(indx_meta(iLookINDEX%nGlce)%varName), glceID)
 if(err/=nf90_noerr) no_glceData = .true.

 ! get nLayer data (reads entire state file)
 if(no_dom)then
   err = nf90_get_var(ncid,snowID,snowData1); call netcdf_err(err,message)
   err = nf90_get_var(ncid,soilID,soilData1); call netcdf_err(err,message)
   if (.not. no_glceData) err = nf90_get_var(ncid,glceID,glceData1); call netcdf_err(err,message)
   if (.not. no_lakeData) err = nf90_get_var(ncid,lakeID,lakeData1); call netcdf_err(err,message)
 else
   err = nf90_get_var(ncid,snowID,snowData2); call netcdf_err(err,message)
   err = nf90_get_var(ncid,soilID,soilData2); call netcdf_err(err,message)
   if (.not. no_glceData) err = nf90_get_var(ncid,glceID,glceData2); call netcdf_err(err,message)
   if (.not. no_lakeData) err = nf90_get_var(ncid,lakeID,lakeData2); call netcdf_err(err,message)
 endif

 ! find the min and max hru indices in the state file
 ixHRUfile_min=huge(1)
 ixHRUfile_max=0
 do iGRU = 1,nGRU
   do iHRU = 1,gru_struc(iGRU)%hruCount
     iHRU_global = gru_struc(iGRU)%hruInfo(iHRU)%hru_nc
     if(iHRU_global < ixHRUfile_min) ixHRUfile_min = iHRU_global
     if(iHRU_global > ixHRUfile_max) ixHRUfile_max = iHRU_global
   end do
 end do

 ! loop over grus in current run to update layer information
 do iGRU = 1,nGRU
   do iHRU = 1,gru_struc(iGRU)%hruCount
    iHRU_global = index_to_hrunc(iGRU,iHRU) ! index of HRU in the netcdf file
    do iDOM = 1, gru_struc(iGRU)%hruInfo(iHRU)%domCount
     if(no_dom)then
       gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSnow = snowData1(iHRU_global)
       gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nLake = lakeData1(iHRU_global)
       gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSoil = soilData1(iHRU_global)
       gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nGlce = glceData1(iHRU_global)
     else
       gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSnow = snowData2(iDOM,iHRU_global)
       gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nLake = lakeData2(iDOM,iHRU_global)
       gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSoil = soilData2(iDOM,iHRU_global)
       gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nGlce = glceData2(iDOM,iHRU_global)
     endif
    end do
  end do
 end do

 ! close file
 call nc_file_close(ncid,err,cmessage)
 if(err/=nf90_noerr)then;message=trim(message)//trim(cmessage);return;end if

 ! cleanup
 deallocate(snowData1,lakeData1,soilData1,glceData1,snowData2,lakeData2,soilData2,glceData2,dom_type)
 deallocate(hru_id,index_to_gruid,index_to_hrunc)
 if(allocated(gru_id)) deallocate(gru_id)  ! not allocated when the file has no gru dimension

 end subroutine read_icond_nlayers


 ! ************************************************************************************************
 ! public subroutine read_icond: read model initial conditions
 ! ************************************************************************************************
 subroutine read_icond(iconFile,                      & ! intent(in):    name of initial conditions file
                       nGRU,                          & ! intent(in):    number of GRUs
                       mparData,                      & ! intent(in):    model parameters
                       progData,                      & ! intent(inout): model prognostic variables
                       bvarData,                      & ! intent(inout): model basin (GRU) variables
                       indxData,                      & ! intent(inout): model indices
                       gridData,                      & ! intent(inout): basin grid parameters and variables
                       no_dom_vars,                   & ! intent(out):   flag that domain variables are not in initial conditions
                       no_ice_vars,                   & ! intent(out):   flag that glacier ice variables are not in initial conditions
                       no_ablfrac,                    & ! intent(out):   flag that glacier ablation fraction variable is not in initial conditions
                       no_icond_enth,                 & ! intent(out):   flag that enthalpy variables are not in the initial conditions
                       err,message)                     ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nr_type
 USE var_lookup,only:iLookVarType                       ! variable lookup structure
 USE var_lookup,only:iLookPROG                          ! variable lookup structure
 USE var_lookup,only:iLookPARAM                         ! variable lookup structure
 USE var_lookup,only:iLookBVAR                          ! variable lookup structure
 USE var_lookup,only:iLookINDEX                         ! variable lookup structure
 USE globalData,only:prog_meta                          ! metadata for prognostic variables
 USE globalData,only:bvar_meta                          ! metadata for basin (GRU) variables
 USE globalData,only:gru_struc                          ! gru-hru mapping structures
 USE globalData,only:startGRU                           ! index of first gru for parallel runs
 USE globalData,only:iname_soil,iname_snow,iname_glce,iname_lake ! named variables to describe the type of 
 USE globalData,only:maxGlaciers                        ! maximum number of glaciers in a GRU
 USE netcdf_util_module,only:nc_file_open               ! open netcdf file
 USE netcdf_util_module,only:nc_file_close              ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                 ! netcdf error handling
 USE data_types,only:gru_hru_dom_doubleVec              ! full double precision structure
 USE data_types,only:gru_hru_dom_intVec                 ! full integer structure
 USE data_types,only:gru_doubleVec                      ! gru-length double precision structure (basin variables)
 USE data_types,only:gru_grid_double                    ! full grid double precision structure
 USE get_ixName_module,only:get_varTypeName             ! to access type strings for error messages
 USE updatState_module,only:updatSoil                   ! update soil states

 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 character(*)               ,intent(in)    :: iconFile                      ! name of netcdf file containing the initial conditions
 integer(i4b)               ,intent(in)    :: nGRU                          ! number of grouped response units in simulation domain
 type(gru_hru_dom_doubleVec),intent(in)    :: mparData                      ! model parameters
 type(gru_hru_dom_doubleVec),intent(inout) :: progData                      ! model prognostic variables
 type(gru_doubleVec)        ,intent(inout) :: bvarData                      ! model basin (GRU) variables
 type(gru_hru_dom_intVec)   ,intent(inout) :: indxData                      ! model indices
 type(gru_grid_double)      ,intent(inout) :: gridData                      ! grid data structure
 logical                    ,intent(out)   :: no_dom_vars                   ! flag that domain variables are not in initial conditions
 logical                    ,intent(out)   :: no_ice_vars                   ! flag that glacier ice variables are not in initial conditions
 logical                    ,intent(out)   :: no_ablfrac                    ! flag that glacier ablation fraction variable is not in initial conditions
 logical                    ,intent(out)   :: no_icond_enth                 ! flag that enthalpy variables are not in initial conditions
 integer(i4b)               ,intent(out)   :: err                           ! error code
 character(*)               ,intent(out)   :: message                       ! returned error message
 ! locals
 character(len=256)                        :: cmessage                      ! downstream error message
 integer(i4b)                              :: fileHRU                       ! number of HRUs in file
 integer(i4b)                              :: fileGRU                       ! number of GRUs in file
 integer(i4b)                              :: fileDOM                       ! number of domains in netcdf file
 integer(i4b)                              :: iVar,i,j                      ! loop indices
 integer(i4b),dimension(1)                 :: nrdx                          ! intermediate array of loop indices for basin variables
 integer(i4b),dimension(7)                 :: ngdx                          ! intermediate array of loop indices for glacier variables
 integer(i4b)                              :: iGRU,iHRU,iDOM,iGlac,iGrid    ! loop indices
 integer(i4b)                              :: dimID                         ! varible dimension ids
 integer(i4b)                              :: ncVarID                       ! variable ID in netcdf file
 character(256)                            :: dimName                       ! not used except as a placeholder in call to inq_dim function
 integer(i4b)                              :: dimLen                        ! data dimensions
 integer(i4b)                              :: ncid                          ! netcdf file ID
 integer(i4b)                              :: iGRU_global                   ! index of GRU in the netcdf file
 integer(i4b)                              :: iHRU_global                   ! index of HRU in the netcdf file
 real(rkind),allocatable                   :: varData2(:,:)                 ! variable data storage
 real(rkind),allocatable                   :: varData3(:,:,:)               ! variable data storage
 integer(i4b)                              :: nSnow,nLake,nSoil,nGlce,nToto ! # layers
 integer(i4b)                              :: noThetaChange                 ! number of layers with no change in total water content (bottom layers)
 integer(i4b)                              :: nTDH                          ! number of points in time-delay 
 integer(i4b)                              :: nGlac                         ! number of glaciers in basin
 integer(i4b)                              :: fileglac                      ! max number of glaciers in any GRU
 integer(i4b)                              :: iLayer,jLayer                 ! layer indices
 logical(lgt)                              :: has_glacier                   ! flag for glacier presence in at least one GRU
 logical(lgt)                              :: has_wetland                   ! flag for wetlandpresence in at least one GRU
 logical(lgt)                              :: no_dom                        ! flag for no domain variable in file
 ! currently only writing restart for progressive variables with these dimensions
 character(len=32),parameter               :: scalDimName   ='scalarv'      ! dimension name for scalar data
 character(len=32),parameter               :: midSoilDimName='midSoil'      ! dimension name for soil-only layers
 character(len=32),parameter               :: midTotoDimName='midToto'      ! dimension name for layered varaiables
 character(len=32),parameter               :: ifcTotoDimName='ifcToto'      ! dimension name for layered varaiables
 character(len=32),parameter               :: tdhDimName    ='tdh'          ! dimension name for time-delay basin variables
 character(len=32),parameter               :: glacierDimName='glac'         ! dimension name for glacier variables
 integer(i8b),allocatable                  :: gru_id(:)                     ! GRU id
 integer(i8b),allocatable                  :: hru_id(:)                     ! HRU id
 integer(i8b),allocatable                  :: glac_id(:,:)                  ! glac id
 integer(i4b),allocatable                  :: index_to_gruid(:)             ! mapping from index to gru_id in gru_struc
 integer(i4b),allocatable                  :: index_to_hrunc(:,:)           ! mapping from index to hru_nc in gru_struc
 integer(i4b),allocatable                  :: index_to_glacid(:,:)          ! mapping from index to glac_id in gru_struc
 logical(lgt)                              :: has_gru_id                    ! flag for whether the file has gru_id
 logical(lgt)                              :: has_hru_id                    ! flag for whether the file has hru_id
 ! --------------------------------------------------------------------------------------------------------
 ! Start procedure here
 err=0; message="read_icond/"
 no_dom = .false.
 has_gru_id = .true.
 has_hru_id = .true.

 ! --------------------------------------------------------------------------------------------------------
 ! (1) read the file
 ! --------------------------------------------------------------------------------------------------------
 ! open netcdf file
 call nc_file_open(iconFile,nf90_nowrite,ncid,err,cmessage)
 if (err/=nf90_noerr) then; message=trim(message)//trim(cmessage); return; end if

 ! get number of HRUs in file
 err = nf90_inq_dimid(ncid,"hru",dimID);               if(err/=nf90_noerr)then; message=trim(message)//'problem finding hru dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncid,dimID,len=fileHRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading hru dimension/'//trim(nf90_strerror(err)); return; end if

 ! get max number of DOMs any HRU in file, if present
 err = nf90_inq_dimid(ncid,"dom",dimId)               
 if(err/=nf90_noerr)then
  fileDOM = 1
  no_dom = .true.
  err=nf90_noerr    ! reset this err
 else
  err = nf90_inquire_dimension(ncid,dimId,len=fileDOM); if(err/=nf90_noerr)then; message=trim(message)//'problem reading dom dimension/'//trim(nf90_strerror(err)); return; end if
 end if

 ! check if the file has the hruId variable
 allocate(hru_id(fileHRU))
 err = nf90_inq_varid(ncid,"hruId",ncVarID)
 if (err/=nf90_noerr)then
   has_hru_id = .false.
   err=nf90_noerr    ! reset this err
 else
  ! read hru_id from netcdf file
   err = nf90_get_var(ncid,ncVarID,hru_id); if (err/=nf90_noerr) then; message=trim(message)//'problem reading hruId'; return; end if
 end if

 ! check if the file has the GRU dimension
 err = nf90_inq_dimid(ncid,"gru",dimID)    
 if(err/=nf90_noerr)then         
   has_gru_id = .false.
   err=nf90_noerr    ! reset this err
 else
   err = nf90_inquire_dimension(ncid,dimID,len=fileGRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading gru dimension/'//trim(nf90_strerror(err)); return; end if
   ! read gru_id from netcdf file
   allocate(gru_id(fileGRU))
   err = nf90_inq_varid(ncid,"gruId",ncVarID)
   if (err/=nf90_noerr) then
     has_gru_id = .false.
     err=nf90_noerr    ! reset this err
   else
     ! read gru_id from netcdf file
     err = nf90_get_var(ncid,ncVarID,gru_id); if (err/=nf90_noerr) then; message=trim(message)//'problem reading gruId'; return; end if
   endif
 end if

 ! Allocate the mapping arrays
 allocate(index_to_gruid(nGRU), index_to_hrunc(nGRU,maxval(gru_struc(:)%hruCount)), index_to_glacid(nGRU,maxGlaciers))

 ! Populate the mapping arrays
 if(has_gru_id .and. has_hru_id)then
   ! if the file has both gru_id and hru_id, use these to populate the mapping arrays
   do iGRU = 1, nGRU
     index_to_gruid(iGRU) = -1  ! Initialize with an invalid index
     do i = 1, fileGRU
       if (gru_struc(iGRU)%gru_id == gru_id(i)) then
         index_to_gruid(iGRU) = i
         do iHRU = 1, gru_struc(iGRU)%hruCount
           index_to_hrunc(iGRU,iHRU) = -1 
           do j = 1, fileHRU
             if (gru_struc(iGRU)%hruInfo(iHRU)%hru_id == hru_id(j)) then
               index_to_hrunc(iGRU,iHRU) = j
               exit
             endif
           end do
         end do ! HRU id loop
         exit
       endif
     end do
   end do ! GRU id loop
 else
   ! assume that the order of the HRUs in the file matches the order of the HRUs in the model attributes
   do iGRU = 1, nGRU
     index_to_gruid(iGRU) = iGRU + startGRU - 1
     do iHRU = 1, gru_struc(iGRU)%hruCount
       index_to_hrunc(iGRU,iHRU) = gru_struc(iGRU)%hruInfo(iHRU)%hru_nc
     end do
   end do
 endif

! get dimension of basin glac variables from initial conditions file 
err = nf90_inq_dimid(ncid,"glac",dimID) ! max number of glaciers in any GRU
if(err/=nf90_noerr)then
  do iGRU = 1,nGRU
    nGlac = 0 ! no glaciers in this GRU
    gru_struc(iGRU)%nGlac = nGlac
  end do
  has_glacier = .false.
  err=nf90_noerr    ! reset this err
else
  has_glacier = .true.
  err = nf90_inquire_dimension(ncid,dimID,len=fileglac);
  if(err/=nf90_noerr)then; message=trim(message)//'problem reading glac dimension from initial condition file/'//trim(nf90_strerror(err)); return; end if

  ! read glac_id from netcdf file
  allocate(glac_id(filegru,fileglac))
  err = nf90_inq_varid(ncid,"glacId ",ncVarID);   if (err/=nf90_noerr) then; message=trim(message)//'problem finding glacId '; return; end if
  err = nf90_get_var(ncid,ncVarID,glac_id);       if (err/=nf90_noerr) then; message=trim(message)//'problem reading glacId '; return; end if

  ! set glacier information (id)
   do iGRU = 1,nGRU
    iGRU_global = index_to_gruid(iGRU) ! index of GRU in the netcdf file
    nGlac = gru_struc(iGRU)%nGlac ! get dimension of basin glacier variables from attribute file, per GRU
    if(nGlac > 0)then
      if (.not. allocated(gru_struc(iGRU)%glacInfo)) then
        allocate(gru_struc(iGRU)%glacInfo(nGlac))
      endif
      gru_struc(iGRU)%glacInfo(1:nGlac)%glac_id = gru_struc(iGRU)%gridInfo(1:nGlac)%grid_id ! glaciers are first grids

      ! Populate the mapping arrays
      do iGlac = 1, nGlac
        index_to_glacid(iGRU,iGlac) = -1  ! Initialize with an invalid index
        do iGrid = 1, gru_struc(iGRU)%nGrid
          if (gru_struc(iGRU)%glacInfo(iGlac)%glac_id == glac_id(iGRU_global,iGrid)) then
            index_to_glacid(iGRU,iGlac) = iGrid
            exit
          endif
        end do
        if (index_to_glacid(iGRU,iGlac) == -1) then
          message=trim(message)//'glacier glacId needs to match a gridId'; err=20; return
        endif
      enddo ! glac loop
    endif
  end do ! GRU loop

  ! set glacier variables to read
  ngdx = (/iLookBVAR%basin__GlacierStorage,iLookBVAR%updateJulDay,iLookBVAR%glacierAblArea,iLookBVAR%glacierAccArea,iLookBVAR%glacIceRunoffFuture,iLookBVAR%glacSnowRunoffFuture,iLookBVAR%glacFirnRunoffFuture/)
  deallocate(glac_id)
 end if

 ! should be looking for wetland domain as well, but no basin variables yet
  has_wetland = .false.

 ! --------------------------------------------------------------------------------------------------------
 ! (2) read the prognostic variables
 ! --------------------------------------------------------------------------------------------------------
 ! loop through prognostic variables
 no_dom_vars=.false.
 no_ice_vars=.false.
 no_ablfrac=.false.
 no_icond_enth=.false.
 do iVar = 1,size(prog_meta)
  ! skip variables that are computed later
  if(prog_meta(iVar)%varName=='scalarCanopyWat'           .or. &
     prog_meta(iVar)%varName=='spectralSnowAlbedoDiffuse' .or. &
     prog_meta(iVar)%varName=='scalarSurfaceTemp'         .or. &
     prog_meta(iVar)%varName=='mLayerVolFracWat'          .or. &
     prog_meta(iVar)%varName=='mLayerHeight'                   ) cycle
  ! get variable id
  err = nf90_inq_varid(ncid,trim(prog_meta(iVar)%varName),ncVarID)
  if(err/=nf90_noerr)then
   if(prog_meta(iVar)%varName=='DOMarea'              .or. &
      prog_meta(iVar)%varName=='DOMelev'              .or. &
      prog_meta(iVar)%varName=='DOMtan_slope '        .or. &
      prog_meta(iVar)%varName=='DOMaspect'            .or. &
      prog_meta(iVar)%varName=='DOMcontourLength'          )then; err=nf90_noerr; no_dom_vars=.true.;cycle; endif ! backwards compatible, may be missing, correct in check_icond
   if(prog_meta(iVar)%varName=='scalarGlceWE'         .or. &
      prog_meta(iVar)%varName=='glacMass4AreaChange'       )then; err=nf90_noerr; no_ice_vars=.true.; cycle; endif ! backwards compatible, may be missing, correct in check_icond
   if(prog_meta(iVar)%varName=='scalarAblFrac'             )then; err=nf90_noerr; no_ablfrac=.true.; cycle; endif ! backwards compatible, may be missing, correct in check_icond
   if(prog_meta(iVar)%varName=='scalarCanairEnthalpy' .or. &
      prog_meta(iVar)%varName=='scalarCanopyEnthalpy' .or. &  
      prog_meta(iVar)%varName=='mLayerEnthalpy'            )then; err=nf90_noerr; no_icond_enth=.true.; cycle; endif ! skip enthalpy variables if not in file
   call netcdf_err(err,message)
   message=trim(message)//': problem with getting variable id, var='//trim(prog_meta(iVar)%varName)
   return
  endif

  ! get variable dimension IDs
  select case (prog_meta(iVar)%varType)
   case (iLookVarType%scalarv); err = nf90_inq_dimid(ncid,trim(scalDimName)   ,dimID); call netcdf_err(err,message)
   case (iLookVarType%midSoil); err = nf90_inq_dimid(ncid,trim(midSoilDimName),dimID); call netcdf_err(err,message)
   case (iLookVarType%midToto); err = nf90_inq_dimid(ncid,trim(midTotoDimName),dimID); call netcdf_err(err,message)
   case (iLookVarType%ifcToto); err = nf90_inq_dimid(ncid,trim(ifcTotoDimName),dimID); call netcdf_err(err,message)
   case default
    message=trim(message)//"unexpectedVariableType[name='"//trim(prog_meta(iVar)%varName)//"';type='"//trim(get_varTypeName(prog_meta(iVar)%varType))//"']"
    err=20; return
  end select
  if(err/=nf90_noerr)then; message=trim(message)//': problem with dimension ids, var='//trim(prog_meta(iVar)%varName); return; endif

  ! get the dimension length
  err = nf90_inquire_dimension(ncid,dimID,dimName,dimLen); call netcdf_err(err,message)
  if(err/=nf90_noerr)then; message=trim(message)//': problem getting the dimension length'; return; endif

  ! initialize the variable data
  allocate(varData3(fileDOM,fileHRU,dimLen),stat=err)
  allocate(varData2(fileHRU,dimLen),stat=err)
  if(err/=0)then; message=trim(message)//'problem allocating HRU variable data'; return; endif

  ! get data
  if(no_dom)then
   err = nf90_get_var(ncid,ncVarID,varData2); call netcdf_err(err,message)
  else
   err = nf90_get_var(ncid,ncVarID,varData3); call netcdf_err(err,message)
  endif
  if(err/=nf90_noerr)then; message=trim(message)//': problem getting the data for variable '//trim(prog_meta(iVar)%varName); return; endif

  ! store data in prognostics structure
  do iGRU = 1,nGRU
   do iHRU = 1,gru_struc(iGRU)%hruCount
    iHRU_global = index_to_hrunc(iGRU,iHRU) ! index of HRU in the netcdf file
    do iDOM = 1, gru_struc(iGRU)%hruInfo(iHRU)%domCount
     ! get the number of layers
     nSnow = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSnow
     nLake = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nLake
     nSoil = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSoil
     nGlce = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nGlce
     nToto = nSnow + nLake + nSoil + nGlce

     ! put the data into data structures and check that none of the values are set to nf90_fill_double
     if(no_dom)then
      select case (prog_meta(iVar)%varType)
       case (iLookVarType%scalarv)
        progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1) = varData2(iHRU_global,1)
        if(abs(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1) - nf90_fill_double) < epsilon(varData2))then; err=20; endif
       case (iLookVarType%midSoil)
        progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1:nSoil) = varData2(iHRU_global,1:nSoil)
        if(any(abs(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1:nSoil) - nf90_fill_double) < epsilon(varData2)))then; err=20; endif
       case (iLookVarType%midToto)
        progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1:nToto) = varData2(iHRU_global,1:nToto)
        if(any(abs(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1:nToto) - nf90_fill_double) < epsilon(varData2)))then; err=20; endif
       case (iLookVarType%ifcToto)
        progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(0:nToto) = varData2(iHRU_global,1:nToto+1)
        if(any(abs(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(0:nToto) - nf90_fill_double) < epsilon(varData2)))then; err=20; endif
       case default
        message=trim(message)//"unexpectedVariableType[name='"//trim(prog_meta(iVar)%varName)//"';type='"//trim(get_varTypeName(prog_meta(iVar)%varType))//"']"
        err=20; return
      end select
     else
      select case (prog_meta(iVar)%varType)
       case (iLookVarType%scalarv)
        progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1)       = varData3(iDOM,iHRU_global,1)
        if(abs(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1) - nf90_fill_double) < epsilon(varData3))then; err=20; endif
       case (iLookVarType%midSoil)
        progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1:nSoil) = varData3(iDOM,iHRU_global,1:nSoil)
        if(any(abs(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1:nSoil) - nf90_fill_double) < epsilon(varData3)))then; err=20; endif
       case (iLookVarType%midToto)
        progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1:nToto) = varData3(iDOM,iHRU_global,1:nToto)
        if(any(abs(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(1:nToto) - nf90_fill_double) < epsilon(varData3)))then; err=20; endif
       case (iLookVarType%ifcToto)
        progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(0:nToto) = varData3(iDOM,iHRU_global,1:nToto+1)
        if(any(abs(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iVar)%dat(0:nToto) - nf90_fill_double) < epsilon(varData3)))then; err=20; endif
       case default
        message=trim(message)//"unexpectedVariableType[name='"//trim(prog_meta(iVar)%varName)//"';type='"//trim(get_varTypeName(prog_meta(iVar)%varType))//"']"
        err=20; return
      end select
     endif
     if(err==20)then; message=trim(message)//"data set to the fill value (name='"//trim(prog_meta(iVar)%varName)//"')"; return; endif

     if(prog_meta(iVar)%varName=='iLayerHeight')then ! last variable in the loop, so we can correct prognostic variables if had legacy starting values
      ! make sure snow albedo is not negative
      if(progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSnowAlbedo)%dat(1) < 0._rkind)then
       progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSnowAlbedo)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%albedoMax)%dat(1)
      endif

      ! initialize the spectral albedo
      progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1:nSpecBand) = progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%scalarSnowAlbedo)%dat(1)
     endif ! (if last variable in the loop)

    end do ! iDOM
   end do ! iHRU
  end do ! iGRU

  ! deallocate storage vector for next variable
  deallocate(varData2,varData3,stat=err)
  if(err/=0)then; message=trim(message)//'problem deallocating HRU variable data'; return; endif

 end do ! end looping through prognostic variables (iVar)

 ! --------------------------------------------------------------------------------------------------------
 ! (3) set number of layers
 ! --------------------------------------------------------------------------------------------------------
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   do iDOM = 1, gru_struc(iGRU)%hruInfo(iHRU)%domCount

    ! save the number of layers
    nSnow = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSnow
    nLake = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nLake
    nSoil = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nSoil
    nGlce = gru_struc(iGRU)%hruInfo(iHRU)%domInfo(iDOM)%nGlce
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nSnow)%dat(1)   = nSnow
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nLake)%dat(1)   = nLake
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nSoil)%dat(1)   = nSoil
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nGlce)%dat(1)   = nGlce
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nLayers)%dat(1) = nSnow + nLake + nSoil + nGlce

    ! define layers that will not have a change in total water content
    noThetaChange = 0
    if(nGlce>0)then
      noThetaChange = nGlce - nMeltingIceLayers
      ! need at least one glacier top layer with a theta change
      if(noThetaChange>=nGlce)then; err=20; message=trim(message)//'number of glacier ice layers without a change in total water content is not less than the number of glacier ice layers'; return; endif
    endif
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%noThetaChange)%dat(1) = noThetaChange

    ! set layer type
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%layerType)%dat(1:nSnow) = iname_snow
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%layerType)%dat((nSnow+1):(nSnow+nLake)) = iname_lake
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%layerType)%dat((nSnow+nLake+1):(nSnow+nLake+nSoil)) = iname_soil
    indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%layerType)%dat((nSnow+nLake+nSoil+1):(nSnow+nLake+nSoil+nGlce)) = iname_glce
   end do
  end do
 end do

 ! --------------------------------------------------------------------------------------------------------
 ! (4) update soil layers (diagnostic variables)
 ! --------------------------------------------------------------------------------------------------------
 ! loop through GRUs and HRUs
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   do iDOM = 1, gru_struc(iGRU)%hruInfo(iHRU)%domCount

    ! loop through soil layers
    do iLayer = 1,indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nSoil)%dat(1)

     ! get layer in the total vector
     jLayer = iLayer+indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nSnow)%dat(1)+indxData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookINDEX%nLake)%dat(1)

     ! update soil layers
     call updatSoil(&
                    ! input
                    progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerTemp          )%dat(jLayer),& ! intent(in): temperature vector (K)
                    progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerMatricHead    )%dat(iLayer),& ! intent(in): matric head (m)
                    mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%vGn_alpha          )%dat(iLayer),& ! intent(in): van Genutchen "alpha" parameter
                    mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%vGn_n              )%dat(iLayer),& ! intent(in): van Genutchen "n" parameter
                    mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%theta_sat          )%dat(iLayer),& ! intent(in): soil porosity (-)
                    mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%theta_res          )%dat(iLayer),& ! intent(in): soil residual volumetric water content (-)
                    1._rkind - 1._rkind/mparData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPARAM%vGn_n)%dat(iLayer),& ! intent(in): van Genutchen "m" parameter (-)
                    ! output
                    progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerVolFracWat    )%dat(jLayer),& ! intent(out): volumetric fraction of total water (-)
                    progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerVolFracLiq    )%dat(jLayer),& ! intent(out): volumetric fraction of liquid water (-)
                    progData%gru(iGRU)%hru(iHRU)%dom(iDOM)%var(iLookPROG%mLayerVolFracIce    )%dat(jLayer),& ! intent(out): volumetric fraction of ice (-)
                    err,cmessage)                                                                            ! intent(out): error control
     if (err/=0) then; message=trim(message)//trim(cmessage); return; end if

    end do  ! looping through soil layers
   end do  ! looping through DOMs
  end do  ! looping through HRUs
 end do  ! looping through GRUs

 ! --------------------------------------------------------------------------------------------------------
 ! (5) get the basin variable(s)
 ! --------------------------------------------------------------------------------------------------------
 ! get dimension of time delay histogram (TDH) from initial conditions file
 err = nf90_inq_dimid(ncid,"tdh",dimID)
 if(err/=nf90_noerr)then
  write(*,*) 'WARNING: routingRunoffFuture is not in the initial conditions file ... using zeros'  ! previously created in var_derive.f90
  err=nf90_noerr    ! reset this err

 else
  ! the state file *does* have the basin variable(s), so process them
  err = nf90_inquire_dimension(ncid,dimID,len=nTDH)
  if(err/=nf90_noerr)then; message=trim(message)//'problem reading tdh dimension from initial condition file/'//trim(nf90_strerror(err)); return; end if

  ! check vs hardwired value set in globalData.f90
  if(nTDH /= nTimeDelay)then
   write(*,*) 'tdh=',nTDH,' nTimeDelay=',nTimeDelay
   message=trim(message)//': state file time delay dimension tdh does not match summa expectation of nTimeDelay set in globalData()'
   return
  endif

  ! loop through specific basin variables (currently 1 but loop provided to enable inclusion of others)
  nrdx = (/iLookBVAR%routingRunoffFuture/)   ! array of desired variable indices
  do i = 1,size(nrdx)
   iVar = nrdx(i)

   ! get tdh dimension Id in file (should be 'tdh')
   err = nf90_inq_dimid(ncid,trim(tdhDimName), dimID)
   if(err/=nf90_noerr)then; message=trim(message)//': problem with dimension ids for tdh vars'; return; endif

   ! get the tdh dimension length (dimName and dimLen are outputs of this call)
   err = nf90_inquire_dimension(ncid,dimID,dimName,dimLen); call netcdf_err(err,message)
   if(err/=nf90_noerr)then; message=trim(message)//': problem getting the dimension length for tdh vars'; return; endif

   ! get tdh-based variable id
   err = nf90_inq_varid(ncid,trim(bvar_meta(iVar)%varName),ncVarID); call netcdf_err(err,message)
   if(err/=nf90_noerr)then; message=trim(message)//': problem with getting basin variable id, var='//trim(bvar_meta(iVar)%varName); return; endif

   ! initialize the tdh variable data
   allocate(varData2(fileGRU,dimLen),stat=err)
   if(err/=0)then; message=trim(message)//'problem allocating GRU variable data'; return; endif

   ! get data
   err = nf90_get_var(ncid,ncVarID,varData2); call netcdf_err(err,message)
   if(err/=nf90_noerr)then; message=trim(message)//': problem getting the data'; return; endif

   ! store data in basin var (bvar) structure
   do iGRU = 1,nGRU
    iGRU_global = index_to_gruid(iGRU) ! index of GRU in the netcdf file
    ! put the data into data structures
    bvarData%gru(iGRU)%var(iVar)%dat(1:nTDH) = varData2(iGRU_global,1:nTDH)
    ! check whether the first values is set to nf90_fill_double
    if(any(abs(bvarData%gru(iGRU)%var(iVar)%dat(1:nTDH) - nf90_fill_double) < epsilon(varData2)))then; err=20; endif
    if(err==20)then; message=trim(message)//"data set to the fill value (name='"//trim(bvar_meta(iVar)%varName)//"')"; return; endif
   end do ! end iGRU loop

   ! deallocate temporary data array for next variable
   deallocate(varData2)
   if(err/=nf90_noerr)then; message=trim(message)//'problem deallocating GRU variable data'; return; endif

  end do ! end looping through basin variables
 endif  ! end if case for tdh variables being in init. cond. file

 ! --------------------------------------------------------------------------------------------------------
 ! (6) get the glacier variables
 ! --------------------------------------------------------------------------------------------------------
 if (has_glacier)then

   do i = 1,size(ngdx) ! loop through specific basin variables
    iVar = ngdx(i)

    ! get glac-based variable id
    err = nf90_inq_varid(ncid,trim(bvar_meta(iVar)%varName),ncVarID)
    if(err/=nf90_noerr)then
      if(iVar == iLookBVAR%updateJulDay)then
        write(*,*) 'WARNING: updateJulDay for last time glacier geometry updated is not in the initial conditions file ... assuming start of simulation'  ! previously created in var_derive.f90
        err=nf90_noerr    ! reset this err
        cycle
      elseif(iVar == iLookBVAR%glacIceRunoffFuture)then ! either all glacier runoff variables are in the file or none
        write(*,*) 'WARNING: glac(Ice,Snow,Firn)RunoffFuture is not in the initial conditions file ... using zeros'  ! previously created in var_derive.f90
        err=nf90_noerr    ! reset this err
        exit ! exit the loop, don't need to check the other glacier runoff variables
      else
        message=trim(message)//': problem with getting basin variable id, var='//trim(bvar_meta(iVar)%varName); return
      endif
    endif

    ! get variable dimension IDs
    select case (bvar_meta(iVar)%varType)
     case (iLookVarType%scalarv); err = nf90_inq_dimid(ncid,trim(scalDimName)   ,dimID); call netcdf_err(err,message)
     case (iLookVarType%glacier); err = nf90_inq_dimid(ncid,trim(glacierDimName),dimID); call netcdf_err(err,message)
     case default
      message=trim(message)//"unexpectedVariableType[name='"//trim(bvar_meta(iVar)%varName)//"';type='"//trim(get_varTypeName(bvar_meta(iVar)%varType))//"']"
      err=20; return
    end select
    if(err/=nf90_noerr)then; message=trim(message)//': problem with dimension ids, var='//trim(bvar_meta(iVar)%varName); return; endif

    ! get the dimension length
    err = nf90_inquire_dimension(ncid,dimID,dimName,dimLen); call netcdf_err(err,message)
    if(err/=nf90_noerr)then; message=trim(message)//': problem getting the dimension length'; return; endif

    ! initialize the variable data
    allocate(varData2(filegru,dimLen),stat=err)
    if(err/=nf90_noerr)then; message=trim(message)//'problem allocating HRU variable data'; return; endif

    ! get data
    err = nf90_get_var(ncid,ncVarID,varData2); call netcdf_err(err,message)
    if(err/=nf90_noerr)then; message=trim(message)//': problem getting the data for variable '//trim(bvar_meta(iVar)%varName); return; endif

    ! store data in basin var structure
    do iGRU = 1,nGRU
      iGRU_global = index_to_gruid(iGRU) ! index of GRU in the netcdf file
      nGlac = gru_struc(iGRU)%nGlac ! get dimension of basin glacier variables from attribute file, per GRU
      select case (bvar_meta(iVar)%varType)
       case (iLookVarType%scalarv)
         bvarData%gru(iGRU)%var(iVar)%dat(1) = varData2(iGRU_global,1)
         if(abs(bvarData%gru(iGRU)%var(iVar)%dat(1) - nf90_fill_double) < epsilon(varData2))then; err=20; endif
       case (iLookVarType%glacier)
         do iGlac = 1, nGlac
           j = index_to_glacid(iGRU,iGlac) ! index of grid in the netcdf file
           bvarData%gru(iGRU)%var(iVar)%dat(iGlac) = varData2(iGRU_global,j)
         enddo
         if(any(abs(bvarData%gru(iGRU)%var(iVar)%dat(1:nGlac) - nf90_fill_double) < epsilon(varData2)))then; err=20; endif
       case default
         message=trim(message)//"unexpectedVariableType[name='"//trim(bvar_meta(iVar)%varName)//"';type='"//trim(get_varTypeName(bvar_meta(iVar)%varType))//"']"
         err=20; return
      end select
      if(err==20)then; message=trim(message)//"data set to the fill value (name='"//trim(bvar_meta(iVar)%varName)//"')"; return; endif
    end do ! iGRU loop

    ! deallocate storage vector for next variable
    deallocate(varData2,stat=err)
    if(err/=0)then; message=trim(message)//'problem deallocating glacier variable data'; return; endif
 
   end do ! end looping through basin variables

   call read_icondGlac(ncid, nGRU, fileGRU, index_to_gruid, bvarData, gridData, err, cmessage); if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 endif  ! end if case for glac variables being in init. cond. file

 call nc_file_close(ncid,err,cmessage)
 if(err/=nf90_noerr)then;message=trim(message)//trim(cmessage);return;end if

 ! cleanup
 deallocate(hru_id,index_to_gruid,index_to_hrunc,index_to_glacid)
 if(allocated(gru_id)) deallocate(gru_id)  ! not allocated when the file has no gru dimension

 end subroutine read_icond

 ! ************************************************************************************************
 ! private subroutine read_icondGlac: read model initial grid conditions for glaciers
 ! ************************************************************************************************
 subroutine read_icondGlac(ncid,                          & ! intent(in):    netcdf file ID
                           nGRU,                          & ! intent(in):    number of GRUs
                           fileGRU,                       & ! intent(in):    number of GRUs in file
                           index_to_gruid,                & ! intent(in):    mapping from gru_id to index in gridData
                           bvarData,                      & ! intent(inout): basin variable data structure
                           gridData,                      & ! intent(inout): grid data structure
                           err,message)                     ! intent(out):   error control
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nr_type
 USE var_lookup,only:iLookGRID                          ! variable glacier grid 
 USE var_lookup,only:iLookBVAR                          ! variable basin variables
 USE globalData,only:grid_meta                          ! metadata for grid variables
 USE globalData,only:gru_struc                          ! gru-hru mapping structures
 USE globalData,only:maxGrid                            ! maximum number of grids in a GRU
 USE netcdf_util_module,only:netcdf_err                 ! netcdf error handling
 USE data_types,only:gru_grid_double                    ! full grid double precision structure
 USE data_types,only:gru_doubleVec                      ! gru-length double precision structure (basin variables)

 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 integer(i4b)           ,intent(in)        :: ncid                     ! netcdf file ID
 integer(i4b)           ,intent(in)        :: nGRU                     ! number of GRUs
 integer(i4b)           ,intent(in)        :: fileGRU                  ! number of GRUs in file
 integer(i4b)           ,intent(in)        :: index_to_gruid(:)        ! mapping from gru_id to index in gridData
 type(gru_doubleVec)    ,intent(in)        :: bvarData                 ! basin variable data structure
 type(gru_grid_double)  ,intent(inout)     :: gridData                 ! grid data structure
 integer(i4b)           ,intent(out)       :: err                      ! error code
 character(*)           ,intent(out)       :: message                  ! returned error message
 ! locals
 integer(i4b)                              :: iGRU_global              ! index of GRU in the netcdf file
 integer(i4b)                              :: filegrid                 ! max number of glacier grids in any GRU
 integer(i4b)                              :: filexgrid                ! number of xgrid points in file
 integer(i4b)                              :: fileygrid                ! number of ygrid points in file
 integer(i4b)                              :: ny,nx                    ! number of grid points for a glacier 
 integer(i4b)                              :: iVar,i,j                 ! loop indices
 integer(i4b),dimension(2)                 :: ngdx                     ! intermediate array of loop indices
 integer(i4b)                              :: iGRU,iGrid               ! loop index
 integer(i4b)                              :: dimID                    ! varible dimension ids
 integer(i4b)                              :: ncVarID                  ! variable ID in netcdf file
 real(rkind),allocatable                   :: varData(:,:,:,:)         ! variable data storage
 integer(i8b),allocatable                  :: grid_id(:,:)             ! grid id
 integer(i4b),allocatable                  :: index_to_gridid(:,:)     ! mapping from index to grid_id in gridData
 real(rkind),parameter                     :: areaTol=1.e-2_rkind      ! tolerance to address precision issues in glacier area summation
 real(rkind)                               :: area,area_grid           ! glacier area for a single glacier (m2)
 ! --------------------------------------------------------------------------------------------------------
 ! Start procedure here
 err=0; message="read_icondGlac/"

 ! --------------------------------------------------------------------------------------------------------
 ! (1) read the file
 ! --------------------------------------------------------------------------------------------------------
 ! get dimension of basin grid variables from initial conditions file
 err = nf90_inq_dimid(ncid,"grid",dimID) ! max number of grids in any GRU
 err = nf90_inquire_dimension(ncid,dimID,len=filegrid);
 if(err/=nf90_noerr)then; message=trim(message)//'problem reading grid dimension from initial condition file/'//trim(nf90_strerror(err)); return; end if

 ! read grid_id from netcdf file
 allocate(grid_id(fileGRU,filegrid))
 err = nf90_inq_varid(ncid,"gridId",ncVarID); if (err/=nf90_noerr) then; message=trim(message)//'problem finding gridId'; return; end if
 err = nf90_get_var(ncid,ncVarID,grid_id);    if (err/=nf90_noerr) then; message=trim(message)//'problem reading gridId'; return; end if

 ! get ygrid dimension of whole file
 err = nf90_inq_dimid(ncid,"ygrid",dimId);                   if(err/=nf90_noerr)then; message=trim(message)//'problem finding ygrid dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncid, dimId, len = fileygrid); if(err/=nf90_noerr)then; message=trim(message)//'problem reading ygrid dimension/'//trim(nf90_strerror(err)); return; end if
 ! get xgrid dimension of whole file
 err = nf90_inq_dimid(ncid,"xgrid",dimId);                   if(err/=nf90_noerr)then; message=trim(message)//'problem finding xgrid dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncid, dimId, len = filexgrid); if(err/=nf90_noerr)then; message=trim(message)//'problem reading xgrid dimension/'//trim(nf90_strerror(err)); return; end if

 ! Allocate the mapping array
 allocate(index_to_gridid(nGRU,maxGrid))

 ! Populate the mapping arrays
 do iGRU = 1, nGRU
   do iGrid = 1, gru_struc(iGRU)%nGrid
     index_to_gridid(iGRU,iGrid) = -1
     do j = 1, gru_struc(iGRU)%nGrid
       if (gru_struc(iGRU)%gridInfo(iGrid)%grid_id == grid_id(index_to_gruid(iGRU),j)) then
         index_to_gridid(iGRU,iGrid) = j
         exit
       endif
     end do  ! grid id loop
   end do ! grid id loop
 end do ! GRU id loop

 ! loop through specific basin grid variables
 ngdx = (/iLookGRID%surface_elev, iLookGRID%debris_thick/)   ! array of desired variable indices
 do i = 1,size(ngdx)
  iVar = ngdx(i)

  ! get grid-based variable id
  err = nf90_inq_varid(ncid,trim(grid_meta(iVar)%varName),ncVarID)
  if(err/=nf90_noerr)then; message=trim(message)//': problem with getting grid variable id, var= '//trim(grid_meta(iVar)%varName); return; endif

  ! initialize the grid variable data
  allocate(varData(fileGRU,filegrid,filexgrid,fileygrid),stat=err)
  if(err/=0)then; message=trim(message)//'problem allocating GRU variable data'; return; endif

  ! get data
  err = nf90_get_var(ncid,ncVarID,varData); call netcdf_err(err,message)
  if(err/=nf90_noerr)then; message=trim(message)//': problem getting the data'; return; endif

  ! store data in grid structure
   do iGRU = 1,nGRU
    iGRU_global = index_to_gruid(iGRU) ! index of GRU in the netcdf file
    do iGrid = 1, gru_struc(iGRU)%nGrid
      j = index_to_gridid(iGRU,iGrid) ! index of grid in the netcdf file
      ny = gru_struc(iGRU)%gridInfo(iGrid)%ny
      nx = gru_struc(iGRU)%gridInfo(iGrid)%nx
      ! put the data into data structures
      gridData%gru(iGRU)%grid(iGrid)%var(iVar)%dat2(1:nx,1:ny) = varData(iGRU_global,j,1:nx,1:ny) 

      ! correct for negative glacier thickness and note area differences with grid approximation
      ! Note, glacier grids are first grids in order
      if (iVar==iLookGRID%surface_elev .and. iGrid<=gru_struc(iGRU)%nGlac) then
        associate(& 
         area_abl     => bvarData%gru(iGRU)%var(iLookBVAR%glacierAblArea)%dat(iGrid)                ,&
         area_acc     => bvarData%gru(iGRU)%var(iLookBVAR%glacierAccArea)%dat(iGrid)                ,&
         bed_elev     => gridData%gru(iGRU)%grid(iGrid)%var(iLookGRID%bed_elev)%dat2(1:nx,1:ny)     ,&
         surface_elev => gridData%gru(iGRU)%grid(iGrid)%var(iLookGRID%surface_elev)%dat2(1:nx,1:ny) ,&
         dx           => gru_struc(iGRU)%gridInfo(iGrid)%dx                                         ,&
         dy           => gru_struc(iGRU)%gridInfo(iGrid)%dy                                          &
         )
         area = area_abl + area_acc
         ! correct if glacier thickness is negative
         surface_elev = merge(bed_elev, surface_elev, surface_elev - bed_elev <= thick4area)
         area_grid = sum(merge(0._rkind, dx*dy, surface_elev - bed_elev <= 0._rkind))
         ! check if area is significantly different from grid approximation
         ! this is okay for the first year as it is considered spin up, subsequent years will use the grid data to compute the area
         if (area_grid>0._rkind) then 
           if (abs(area/area_grid-1._rkind) >areaTol) then
             write(*,*) 'WARNING: Area of glacier ',iGrid,' in GRU ',iGRU,' starts at ', area/area_grid, ' times the grid approximation but will be calculated from the grid data after a year.'
             write(*,*) 'If this is not expected, check the glacier grid data in the initial grid conditions file (perhaps the grid should be refined or the glacier area should be adjusted).'
           endif      
         endif
        end associate
      endif

      ! check whether the first values is set to nf90_fill_double
      if(any(abs(gridData%gru(iGRU)%grid(iGrid)%var(iVar)%dat2(1:nx,1:ny) - nf90_fill_double) < epsilon(varData)))then; err=20; endif
      if(err==20)then; message=trim(message)//"data set to the fill value (name='"//trim(grid_meta(iVar)%varName)//"')"; return; endif
    enddo
  end do ! end iGRU loop

  ! deallocate temporary data array for next variable
  deallocate(varData, stat=err)
  if(err/=0)then; message=trim(message)//'problem deallocating GRU variable data'; return; endif
 enddo ! end looping through basin variables

 ! cleanup
 deallocate(grid_id,index_to_gridid)

 end subroutine read_icondGlac


end module read_icond_module

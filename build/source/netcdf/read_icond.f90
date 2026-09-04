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
USE globalData,only:nTimeDelay    ! number of timesteps in the time delay histogram
USE globalData,only:nSpecBand     ! number of spectral bands

implicit none
private
public::read_icond
public::read_icond_nlayers

contains

 ! ************************************************************************************************
 ! public subroutine read_icond_nlayers: read model initial conditions file for number of snow/soil layers
 ! ************************************************************************************************
 subroutine read_icond_nlayers(iconFile,nGRU,indx_meta,err,message)
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nr_type
 USE var_lookup,only:iLookINDEX                        ! variable lookup structure
 USE globalData,only:gru_struc                         ! gru-hru mapping structures
 USE netcdf_util_module,only:nc_file_close             ! close netcdf file
 USE netcdf_util_module,only:nc_file_open              ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                ! netcdf error handling
 USE data_types,only:gru_hru_intVec                    ! actual data
 USE data_types,only:var_info                          ! metadata
 implicit none

 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 character(*)  ,intent(in)   :: iconFile            ! name of input (restart) file
 integer(i4b)  ,intent(in)   :: nGRU                ! total # of GRUs in run domain
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
 integer(i4b)                :: snowID, soilID      ! netcdf variable ids
 integer(i4b)                :: iGRU, iHRU          ! loop indexes
 integer(i4b)                :: iHRU_global         ! index of HRU in the netcdf file
 integer(i4b),allocatable    :: snowData(:)         ! number of snow layers in all HRUs
 integer(i4b),allocatable    :: soilData(:)         ! number of soil layers in all HRUs
 character(len=256)          :: cmessage            ! downstream error message
 integer(i8b),allocatable    :: gru_id(:)           ! GRU id
 integer(i8b),allocatable    :: hru_id(:)           ! HRU id
 integer(i4b),allocatable    :: index_to_hrunc(:,:) ! mapping from index to hru_nc in gru_struc
 logical(lgt)                :: has_gru_id          ! flag for whether the file has gru_id
 logical(lgt)                :: has_hru_id          ! flag for whether the file has hru_id
 ! --------------------------------------------------------------------------------------------------------
 ! initialize error message
 err=0
 message = 'read_icond_nlayers/'
 has_gru_id = .true.
 has_hru_id = .true.

 ! open netcdf file
 call nc_file_open(iconFile,nf90_nowrite,ncid,err,cmessage);
 if (err/=nf90_noerr) then; message=trim(message)//trim(cmessage); return; end if

 ! get number of HRUs in file
 err = nf90_inq_dimid(ncid,"hru",dimId);               if(err/=nf90_noerr)then; message=trim(message)//'problem finding hru dimension/'//trim(nf90_strerror(err)); return; end if
 err = nf90_inquire_dimension(ncid,dimId,len=fileHRU); if(err/=nf90_noerr)then; message=trim(message)//'problem reading hru dimension/'//trim(nf90_strerror(err)); return; end if

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
 err = nf90_inq_dimid(ncid,"gru",dimID)    
 if(err/=nf90_noerr)then         
   write(*,*) 'WARNING: GRU is not in the initial conditions file ... assuming GRUs in attribute order'
   has_gru_id = .false.
   allocate(gru_id(1)) ! just allocate something to avoid problems with the deallocation at the end
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
 allocate(index_to_hrunc(nGRU,maxval(gru_struc(:)%hruCount)))

 ! Populate the mapping arrays
 if(has_gru_id .and. has_hru_id)then
   ! if the file has both gru_id and hru_id, use these to populate the mapping arrays
   do iGRU = 1, nGRU
     do i = 1, fileGRU
       if (gru_struc(iGRU)%gru_id == gru_id(i)) then
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
     do iHRU = 1, gru_struc(iGRU)%hruCount
       index_to_hrunc(iGRU,iHRU) = gru_struc(iGRU)%hruInfo(iHRU)%hru_nc
     end do
   end do
 endif

 ! allocate storage for reading from file (allocate entire file size, even when doing subdomain run)
 allocate(snowData(fileHRU))
 allocate(soilData(fileHRU))
 snowData = 0
 soilData = 0

 ! get netcdf ids for the variables holding number of snow and soil layers in each hru
 err = nf90_inq_varid(ncid,trim(indx_meta(iLookINDEX%nSnow)%varName),snowID); call netcdf_err(err,message)
 err = nf90_inq_varid(ncid,trim(indx_meta(iLookINDEX%nSoil)%varName),soilID); call netcdf_err(err,message)

 ! get nSnow and nSoil data (reads entire state file)
 err = nf90_get_var(ncid,snowID,snowData); call netcdf_err(err,message)
 err = nf90_get_var(ncid,soilID,soilData); call netcdf_err(err,message)

 ! find the min and max hru indices in the state file
 ixHRUfile_min=huge(1)
 ixHRUfile_max=0
 do iGRU = 1,nGRU
   do iHRU = 1,gru_struc(iGRU)%hruCount
     if(gru_struc(iGRU)%hruInfo(iHRU)%hru_nc < ixHRUfile_min) ixHRUfile_min = gru_struc(iGRU)%hruInfo(iHRU)%hru_nc
     if(gru_struc(iGRU)%hruInfo(iHRU)%hru_nc > ixHRUfile_max) ixHRUfile_max = gru_struc(iGRU)%hruInfo(iHRU)%hru_nc
   end do
 end do

 ! loop over grus in current run to update snow/soil layer information
 do iGRU = 1,nGRU
   do iHRU = 1,gru_struc(iGRU)%hruCount
     iHRU_global = index_to_hrunc(iGRU,iHRU) ! index of HRU in the netcdf file
     gru_struc(iGRU)%hruInfo(iHRU)%nSnow = snowData(iHRU_global)
     gru_struc(iGRU)%hruInfo(iHRU)%nSoil = soilData(iHRU_global)
   end do
 end do

 ! close file
 call nc_file_close(ncid,err,cmessage)
 if(err/=nf90_noerr)then;message=trim(message)//trim(cmessage);return;end if

 ! cleanup
 deallocate(snowData,soilData)
 deallocate(gru_id,hru_id,index_to_hrunc)

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
                       no_icond_enth,                 & ! intent(out):   flag that enthalpy variables are not in the file
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
 USE globalData,only:iname_soil,iname_snow              ! named variables to describe the type of layer
 USE netcdf_util_module,only:nc_file_open               ! open netcdf file
 USE netcdf_util_module,only:nc_file_close              ! close netcdf file
 USE netcdf_util_module,only:netcdf_err                 ! netcdf error handling
 USE data_types,only:gru_hru_doubleVec                  ! full double precision structure
 USE data_types,only:gru_hru_intVec                     ! full integer structure
 USE data_types,only:gru_doubleVec                      ! gru-length double precision structure (basin variables)
 USE get_ixName_module,only:get_varTypeName             ! to access type strings for error messages
 USE updatState_module,only:updatSoil                   ! update soil states

 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 character(*)           ,intent(in)     :: iconFile                 ! name of netcdf file containing the initial conditions
 integer(i4b)           ,intent(in)     :: nGRU                     ! number of grouped response units in simulation domain
 type(gru_hru_doubleVec),intent(in)     :: mparData                 ! model parameters
 type(gru_hru_doubleVec),intent(inout)  :: progData                 ! model prognostic variables
 type(gru_doubleVec)    ,intent(inout)  :: bvarData                 ! model basin (GRU) variables
 type(gru_hru_intVec)   ,intent(inout)  :: indxData                 ! model indices
 logical                ,intent(out)    :: no_icond_enth            ! flag that enthalpy variables are not in the file
 integer(i4b)           ,intent(out)    :: err                      ! error code
 character(*)           ,intent(out)    :: message                  ! returned error message
 ! locals
 character(len=256)                     :: cmessage                 ! downstream error message
 integer(i4b)                           :: fileHRU                  ! number of HRUs in file
 integer(i4b)                           :: fileGRU                  ! number of GRUs in file
 integer(i4b)                           :: iVar,i,j                 ! loop indices
 integer(i4b),dimension(1)              :: ndx                      ! intermediate array of loop indices
 integer(i4b)                           :: iGRU                     ! loop index
 integer(i4b)                           :: iHRU                     ! loop index
 integer(i4b)                           :: dimID                    ! varible dimension ids
 integer(i4b)                           :: ncVarID                  ! variable ID in netcdf file
 character(256)                         :: dimName                  ! not used except as a placeholder in call to inq_dim function
 integer(i4b)                           :: dimLen                   ! data dimensions
 integer(i4b)                           :: ncid                     ! netcdf file ID
 integer(i4b)                           :: iGRU_global              ! index of GRU in the netcdf file
 integer(i4b)                           :: iHRU_global              ! index of HRU in the netcdf file
 real(rkind),allocatable                :: varData(:,:)             ! variable data storage
 integer(i4b)                           :: nSoil, nSnow, nToto      ! # layers
 integer(i4b)                           :: nTDH                     ! number of points in time-delay histogram
 integer(i4b)                           :: iLayer,jLayer            ! layer indices
 character(len=32),parameter            :: scalDimName   ='scalarv' ! dimension name for scalar data
 character(len=32),parameter            :: midSoilDimName='midSoil' ! dimension name for soil-only layers
 character(len=32),parameter            :: midTotoDimName='midToto' ! dimension name for layered varaiables
 character(len=32),parameter            :: ifcTotoDimName='ifcToto' ! dimension name for layered varaiables
 character(len=32),parameter            :: tdhDimName    ='tdh'     ! dimension name for time-delay basin variables
 integer(i8b),allocatable               :: gru_id(:)                ! GRU id
 integer(i8b),allocatable               :: hru_id(:)                ! HRU id
 integer(i4b),allocatable               :: index_to_gruid(:)        ! mapping from index to gru_id in gru_struc
 integer(i4b),allocatable               :: index_to_hrunc(:,:)      ! mapping from index to hru_nc in gru_struc
 logical(lgt)                           :: has_gru_id               ! flag for whether the file has gru_id
 logical(lgt)                           :: has_hru_id               ! flag for whether the file has hru_id
 ! --------------------------------------------------------------------------------------------------------
 ! Start procedure here
 err=0; message="read_icond/"
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
   allocate(gru_id(1)) ! just allocate something to avoid problems with the deallocation at the end
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

 ! --------------------------------------------------------------------------------------------------------
 ! (2) read the prognostic variables
 ! -------------------------------------------------------------------------------------------------------- 
 ! loop through prognostic variables
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
   if(prog_meta(iVar)%varName=='scalarCanairEnthalpy'     .or. &
      prog_meta(iVar)%varName=='scalarCanopyEnthalpy'     .or. &  
      prog_meta(iVar)%varName=='mLayerEnthalpy'                )then; err=nf90_noerr; no_icond_enth=.true.; cycle; endif ! skip enthalpy variables if not in file
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
  allocate(varData(fileHRU,dimLen),stat=err)
  if(err/=0)then; message=trim(message)//'problem allocating HRU variable data'; return; endif

  ! get data
  err = nf90_get_var(ncid,ncVarID,varData); call netcdf_err(err,message)
  if(err/=nf90_noerr)then; message=trim(message)//': problem getting the data for variable '//trim(prog_meta(iVar)%varName); return; endif

  ! store data in prognostics structure
  do iGRU = 1,nGRU
   do iHRU = 1,gru_struc(iGRU)%hruCount
    iHRU_global = index_to_hrunc(iGRU,iHRU) ! index of HRU in the netcdf file
    ! get the number of layers
    nSnow = gru_struc(iGRU)%hruInfo(iHRU)%nSnow
    nSoil = gru_struc(iGRU)%hruInfo(iHRU)%nSoil
    nToto = nSnow + nSoil

    ! put the data into data structures and check that none of the values are set to nf90_fill_double
    select case (prog_meta(iVar)%varType)
     case (iLookVarType%scalarv)
      progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat(1)       = varData(iHRU_global,1)
      if(abs(progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat(1) - nf90_fill_double) < epsilon(varData))then; err=20; endif
     case (iLookVarType%midSoil)
      progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat(1:nSoil) = varData(iHRU_global,1:nSoil)
      if(any(abs(progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat(1:nSoil) - nf90_fill_double) < epsilon(varData)))then; err=20; endif
     case (iLookVarType%midToto)
      progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat(1:nToto) = varData(iHRU_global,1:nToto)
      if(any(abs(progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat(1:nToto) - nf90_fill_double) < epsilon(varData)))then; err=20; endif
     case (iLookVarType%ifcToto)
      progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat(0:nToto) = varData(iHRU_global,1:nToto+1)
      if(any(abs(progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat(0:nToto) - nf90_fill_double) < epsilon(varData)))then; err=20; endif
     case default
      message=trim(message)//"unexpectedVariableType[name='"//trim(prog_meta(iVar)%varName)//"';type='"//trim(get_varTypeName(prog_meta(iVar)%varType))//"']"
      err=20; return
    end select
    if(err==20)then; message=trim(message)//"data set to the fill value (name='"//trim(prog_meta(iVar)%varName)//"')"; return; endif

    if(prog_meta(iVar)%varName=='iLayerHeight')then ! last variable in the loop, so we can correct prognostic variables if had legacy starting values
     ! make sure snow albedo is not negative
     if(progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSnowAlbedo)%dat(1) < 0._rkind)then
      progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSnowAlbedo)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMax)%dat(1)
     endif

     ! initialize the spectral albedo
     progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1:nSpecBand) = progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSnowAlbedo)%dat(1)
    endif ! (if last variable in the loop)

   end do ! iHRU
  end do ! iGRU

  ! deallocate storage vector for next variable
  deallocate(varData, stat=err)
  if(err/=0)then; message=trim(message)//'problem deallocating HRU variable data'; return; endif

 end do ! end looping through prognostic variables (iVar)

 ! --------------------------------------------------------------------------------------------------------
 ! (3) set number of layers
 ! --------------------------------------------------------------------------------------------------------
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount

   ! save the number of layers
   nSnow = gru_struc(iGRU)%hruInfo(iHRU)%nSnow
   nSoil = gru_struc(iGRU)%hruInfo(iHRU)%nSoil
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSnow)%dat(1)   = nSnow
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSoil)%dat(1)   = nSoil
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nLayers)%dat(1) = nSnow + nSoil

   ! set layer type
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%layerType)%dat(1:nSnow) = iname_snow
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%layerType)%dat((nSnow+1):(nSnow+nSoil)) = iname_soil
  end do
 end do

 ! --------------------------------------------------------------------------------------------------------
 ! (4) update soil layers (diagnostic variables)
 ! --------------------------------------------------------------------------------------------------------
 ! loop through GRUs and HRUs
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount

   ! loop through soil layers
   do iLayer = 1,indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSoil)%dat(1)

    ! get layer in the total vector
    jLayer = iLayer+indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSnow)%dat(1)

    ! update soil layers
    call updatSoil(&
                    ! input
                    progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerTemp          )%dat(jLayer),& ! intent(in): temperature vector (K)
                    progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerMatricHead    )%dat(iLayer),& ! intent(in): matric head (m)
                    mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%vGn_alpha          )%dat(iLayer),& ! intent(in): van Genutchen "alpha" parameter
                    mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%vGn_n              )%dat(iLayer),& ! intent(in): van Genutchen "n" parameter
                    mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%theta_sat          )%dat(iLayer),& ! intent(in): soil porosity (-)
                    mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%theta_res          )%dat(iLayer),& ! intent(in): soil residual volumetric water content (-)
                    1._rkind - 1._rkind/mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%vGn_n)%dat(iLayer),& ! intent(in): van Genutchen "m" parameter (-)
                    ! output
                    progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracWat    )%dat(jLayer),& ! intent(out): volumetric fraction of total water (-)
                    progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracLiq    )%dat(jLayer),& ! intent(out): volumetric fraction of liquid water (-)
                    progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracIce    )%dat(jLayer),& ! intent(out): volumetric fraction of ice (-)
                    err,cmessage)                                                                   ! intent(out): error control
    if (err/=0) then; message=trim(message)//trim(cmessage); return; end if

   end do  ! looping through soil layers
  end do  ! looping through HRUs
 end do  ! looping through GRUs

 ! --------------------------------------------------------------------------------------------------------
 ! (5) get the basin variable(s)
 ! --------------------------------------------------------------------------------------------------------

 ! get dimension of time delay histogram (TDH) from initial conditions file
 err = nf90_inq_dimid(ncid,"tdh",dimID);
 if(err/=nf90_noerr)then
  write(*,*) 'WARNING: routingRunoffFuture is not in the initial conditions file ... using zeros'  ! previously created in var_derive.f90
  err=nf90_noerr    ! reset this err

 else
  ! the state file *does* have the basin variable(s), so process them
  err = nf90_inquire_dimension(ncid,dimID,len=nTDH);
  if(err/=nf90_noerr)then; message=trim(message)//'problem reading tdh dimension from initial condition file/'//trim(nf90_strerror(err)); return; end if

  ! check vs hardwired value set in globalData.f90
  if(nTDH /= nTimeDelay)then
   write(*,*) 'tdh=',nTDH,' nTimeDelay=',nTimeDelay
   message=trim(message)//': state file time delay dimension tdh does not match summa expectation of nTimeDelay set in globalData()'
   return
  endif

  ! loop through specific basin variables (currently 1 but loop provided to enable inclusion of others)
  ndx = (/iLookBVAR%routingRunoffFuture/)   ! array of desired variable indices
  do i = 1,size(ndx)
   iVar = ndx(i)

   ! get tdh dimension Id in file (should be 'tdh')
   err = nf90_inq_dimid(ncid,trim(tdhDimName), dimID);
   if(err/=nf90_noerr)then; message=trim(message)//': problem with dimension ids for tdh vars'; return; endif

   ! get the tdh dimension length (dimName and dimLen are outputs of this call)
   err = nf90_inquire_dimension(ncid,dimID,dimName,dimLen); call netcdf_err(err,message)
   if(err/=nf90_noerr)then; message=trim(message)//': problem getting the dimension length for tdh vars'; return; endif

   ! get tdh-based variable id
   err = nf90_inq_varid(ncid,trim(bvar_meta(iVar)%varName),ncVarID); call netcdf_err(err,message)
   if(err/=nf90_noerr)then; message=trim(message)//': problem with getting basin variable id, var='//trim(bvar_meta(iVar)%varName); return; endif

   ! initialize the tdh variable data
   allocate(varData(fileGRU,dimLen),stat=err)
   if(err/=0)then; message=trim(message)//'problem allocating GRU variable data'; return; endif

   ! get data
   err = nf90_get_var(ncid,ncVarID,varData); call netcdf_err(err,message)
   if(err/=nf90_noerr)then; message=trim(message)//': problem getting the data'; return; endif

   ! store data in basin var (bvar) structure
   do iGRU = 1,nGRU
    iGRU_global = index_to_gruid(iGRU) ! index of GRU in the netcdf file
    ! put the data into data structures
    bvarData%gru(iGRU)%var(iVar)%dat(1:nTDH) = varData(iGRU_global,1:nTDH)
    ! check whether the first values is set to nf90_fill_double
    if(any(abs(bvarData%gru(iGRU)%var(iVar)%dat(1:nTDH) - nf90_fill_double) < epsilon(varData)))then; err=20; endif
    if(err==20)then; message=trim(message)//"data set to the fill value (name='"//trim(bvar_meta(iVar)%varName)//"')"; return; endif
   end do ! end iGRU loop

   ! deallocate temporary data array for next variable
   deallocate(varData, stat=err)
   if(err/=nf90_noerr)then; message=trim(message)//'problem deallocating GRU variable data'; return; endif

  end do ! end looping through basin variables
 endif  ! end if case for tdh variables being in init. cond. file
 
 call nc_file_close(ncid,err,cmessage)
 if(err/=0)then;message=trim(message)//trim(cmessage);return;end if

 ! cleanup
 deallocate(gru_id,hru_id,index_to_gruid,index_to_hrunc)

 end subroutine read_icond

end module read_icond_module

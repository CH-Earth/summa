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

module summa_read_param_module

! missing values
USE globalData,only:integerMissing  ! missing integer
USE globalData,only:realMissing     ! missing real number

! input sizes
USE globalData,only:maxSoilLayers          ! maximum number of soil layers

! common modules
USE nr_type
USE netcdf
USE netcdf_util_module,only:nc_file_close  ! close netcdf file
USE netcdf_util_module,only:nc_file_open   ! open netcdf file
USE netcdf_util_module,only:netcdf_err     ! netcdf error handling function

! data types
USE data_types,only:gru_double             ! spatial double data type:  x%gru(:)%var(:)
USE data_types,only:gru_hru_int8           ! spatial integer data type: x%gru(:)%hru(:)%var(:)
USE data_types,only:gru_hru_doubleVec      ! spatial double data type:  x%gru(:)%hru(:)%var(:)%dat(:)

implicit none
private
public::read_param
contains


 ! ************************************************************************************************
 ! public subroutine read_param: read trial model parameter values
 ! ************************************************************************************************
 subroutine read_param(nGRU_local, nHRU_local, &
                      idStruct, mparStruct, bparStruct, err, message)

 USE summaFileManager,only:SETTINGS_PATH                     ! path for metadata files
 USE summaFileManager,only:PARAMETER_TRIAL                   ! file with parameter trial values
 USE get_ixname_module,only:get_ixParam,get_ixBpar           ! access function to find index of elements in structure
 USE globalData,only:index_map,gru_struc                     ! mapping from global HRUs to the elements in the data structures
 USE var_lookup,only:iLookPARAM,iLookTYPE,iLookID            ! named variables to index elements of the data vectors
 
 implicit none
 
 integer(i4b)            , intent(in)    :: nGRU_local       ! number of GRUs assigned to this rank
 integer(i4b)            , intent(in)    :: nHRU_local       ! number of HRUs assigned to this rank
 type(gru_hru_int8)      , intent(in)    :: idStruct         ! GRU/HRU identifiers for the local domain
 type(gru_hru_doubleVec) , intent(inout) :: mparStruct       ! model parameters for each local HRU
 type(gru_double)        , intent(inout) :: bparStruct       ! basin parameters for each local GRU
 integer(i4b)            , intent(out)   :: err              ! error code
 character(*)            , intent(out)   :: message          ! error message
 
 ! local variables
 character(len=1024)                     :: cmessage         ! error message for downwind routine
 character(LEN=1024)                     :: infile           ! input filename
 integer(i4b)                            :: ixParam          ! index of the model parameter in the data structure
 logical(lgt)                            :: found_hru_id     ! flag if HRU ID exists

 ! indices/metadata in the NetCDF file   
 integer(i4b)                            :: ncid             ! netcdf id
 integer(i4b)                            :: nDims            ! number of dimensions
 integer(i4b)                            :: nVars            ! number of variables
 integer(i4b)                            :: iDimId           ! dimension index
 integer(i4b)                            :: iVarId           ! variable index
 character(LEN=64)                       :: dimName          ! dimension name
 character(LEN=64)                       :: parName          ! parameter name
 integer(i4b)                            :: dimLength        ! dimension length
 integer(i4b)                            :: nHRU_file        ! number of HRUs in the parafile
 integer(i4b)                            :: nGRU_file        ! number of GRUs in the parafile
 integer(i4b)                            :: nSoil_file       ! number of soil layers in the file
 integer(i4b)                            :: idim_list(2)     ! list of dimension ids
 
 ! data in the netcdf file
 integer(i4b)                            :: parLength        ! length of the parameter data
 real(rkind),allocatable                 :: parVector(:)     ! model parameter vector
 logical                                 :: fexist           ! inquire whether the paramTrial file exists

 ! mapping arrays
 integer(i4b),allocatable                :: index_to_hrunc(:) ! local HRU -> parameter-file HRU index
 integer(i4b),allocatable                :: index_to_grunc(:) ! local GRU -> parameter-file GRU index
 
 ! identifiers in the parameter file
 integer(i8b),allocatable                :: hruId(:)           ! HRU identifiers in the parameter file
 integer(i8b),allocatable                :: gruId(:)           ! GRU identifiers in the parameter file
 
 ! mapping indices
 integer(i4b)                            :: iHRU               ! local HRU index
 integer(i4b)                            :: iGRU               ! local GRU index
 integer(i4b)                            :: localHRU_ix        ! HRU index within local GRU
 integer(i4b)                            :: fHRU               ! HRU index in parameter file
 
 ! file metadata
 logical(lgt)                            :: has_gru_id         ! .true. if gruId exists in parameter file

 ! Start procedure here
 err=0; message="read_param/"

 ! **********************************************************************************************
 ! * open files, etc.
 ! **********************************************************************************************

 ! build filename
 infile = trim(SETTINGS_PATH)//trim(PARAMETER_TRIAL)

 ! check whether the user-specified parameter file exists
 inquire(file=trim(infile),exist=fexist)

 if (.not.fexist) then
   write(*,'(/,A,/)') 'WARNING: trial parameter file not found; using default parameters. '// &
                      'Check the file manager path if this is not the intended behavior.'
   return
 endif

 ! open trial parameters file if it exists
 call nc_file_open(trim(infile),nf90_nowrite,ncid,err,cmessage)
 if(err/=nf90_noerr)then; message=trim(message)//trim(cmessage); return; end if

 ! get the number of variables in the parameter file
 err=nf90_inquire(ncid, nDimensions=nDims, nVariables=nVars)
 call netcdf_err(err,message); if (err/=nf90_noerr) then; err=20; return; end if

 ! initialize the number of HRUs
 nHRU_file=integerMissing
 nGRU_file=integerMissing

 ! get the length of the dimensions
 do iDimId=1,nDims
  ! get the dimension name and length
  err=nf90_inquire_dimension(ncid, iDimId, name=dimName, len=dimLength)
  if(err/=nf90_noerr)then; message=trim(message)//trim(cmessage); return; end if
  ! get the number of HRUs
  if(trim(dimName)=='hru') nHRU_file=dimLength
  if(trim(dimName)=='gru') nGRU_file=dimLength
 end do

 ! check HRU dimension exists
 if(nHRU_file==integerMissing)then
  message=trim(message)//'unable to identify HRU dimension in file '//trim(infile)
  err=20; return
 endif

 ! allocate hruID vector
 allocate(hruId(nHRU_file))

 ! **********************************************************************************************
 ! * read the GRU index and build local-to-file mapping
 ! **********************************************************************************************

 has_gru_id=.false.
 
 if(nGRU_file/=integerMissing)then
 
   allocate(gruId(nGRU_file))
 
   err=nf90_inq_varid(ncid,'gruId',iVarId)
 
   if(err==nf90_noerr)then
 
     has_gru_id=.true.
 
     err=nf90_get_var(ncid,iVarId,gruId)
     if(err/=nf90_noerr)then
       message=trim(message)//'problem reading gruId'
       return
     endif
 
     allocate(index_to_grunc(nGRU_local))
     index_to_grunc=-1
 
     do iGRU=1,nGRU_local
 
       index_to_grunc(iGRU)=findloc(gruId,gru_struc(iGRU)%gru_id,dim=1)
 
       if(index_to_grunc(iGRU)<1)then
         message=trim(message)//'problem finding GRU in parameter file'
         err=20; return
       endif
 
     enddo
 
   else
     err=nf90_noerr
   endif
 
 endif

 ! **********************************************************************************************
 ! * read the HRU index and build local-to-file mapping
 ! **********************************************************************************************

 found_hru_id = .false.

 ! loop through the parameters in the NetCDF file
 do iVarId=1,nVars

   ! get the parameter name
   err=nf90_inquire_variable(ncid,iVarId,name=parName)
   call netcdf_err(err,message); if(err/=nf90_noerr)then; err=20; return; endif
  
   ! special case of the HRU id
   if(trim(parName)=='hruIndex' .or. trim(parName)=='hruId')then

     found_hru_id = .true.

     ! read HRU IDs from the parameter file
     err=nf90_get_var(ncid,iVarId,hruId)
     if(err/=nf90_noerr)then; message=trim(message)//trim(cmessage); return; endif
  
     ! build mapping from local HRUs to parameter-file HRU indices
     allocate(index_to_hrunc(nHRU_local))
     index_to_hrunc=-1
  
     do iHRU=1,nHRU_local
  
       iGRU        = index_map(iHRU)%gru_ix
       localHRU_ix = index_map(iHRU)%localHRU_ix
  
       index_to_hrunc(iHRU)=findloc(hruId, &
         idStruct%gru(iGRU)%hru(localHRU_ix)%var(iLookID%hruId),dim=1)
  
       if(index_to_hrunc(iHRU)<1)then
         err=20; message=trim(message)//'problem finding HRU in parameter file'; return
       endif
  
     enddo
  
     exit  ! read the hruID successfully
 
   endif
 
 enddo

 if(.not.found_hru_id)then
   message=trim(message)//'parameter file does not contain hruId or hruIndex'
   err=20; return
 endif

 ! **********************************************************************************************
 ! * read the local parameters and the basin parameters
 ! **********************************************************************************************

 ! loop through the parameters in the NetCDF file
 do iVarId=1,nVars

  ! get the parameter name
  err=nf90_inquire_variable(ncid, iVarId, name=parName)
  call netcdf_err(err,message); if (err/=nf90_noerr) then; err=20; return; end if

  ! get the local parameters
  ixParam = get_ixParam( trim(parName) )
  if(ixParam/=integerMissing)then

   ! **********************************************************************************************
   ! * read the local parameters
   ! **********************************************************************************************

   ! get the variable shape
   err=nf90_inquire_variable(ncid, iVarId, nDims=nDims, dimids=idim_list)
   if(err/=nf90_noerr)then; message=trim(message)//trim(cmessage); return; end if

   select case(nDims)
  
     case (1); parLength=1

     case (2)

       ! get the information on the 2nd dimension for 2-d variables
       err=nf90_inquire_dimension(ncid, idim_list(2), dimName, nSoil_file)
       if(err/=nf90_noerr)then; message=trim(message)//trim(cmessage); return; end if

       ! check that it is the depth dimension
       if(trim(dimName)/='depth')then
        message=trim(message)//'expect 2nd dimension of 2-d variable to be depth (dimension name = '//trim(dimName)//')'
        err=20; return
       endif

       ! check that the dimension length is correct (maxSoilLayers is the maximum number of soil layers in the model)
       if(maxSoilLayers /= nSoil_file)then
        message=trim(message)//'unexpected number of soil layers in parameter file'
        err=20; return
       endif

       ! define parameter length
       parLength = nSoil_file

     case default
       message=trim(message)//'unexpected number of dimensions for parameter '//trim(parName)
       err=20; return

   end select

   ! allocate space for model parameters
   allocate(parVector(parLength),stat=err)
   if(err/=0)then
    message=trim(message)//'problem allocating space for parameter vector'
    err=20; return
   endif

   ! loop through HRUs
   do iHRU=1,nHRU_local

    ! map to the GRUs and HRUs
    iGRU        = index_map(iHRU)%gru_ix
    localHRU_ix = index_map(iHRU)%localHRU_ix
    fHRU        = index_to_hrunc(iHRU)

    ! read parameter data
    select case(nDims)
     case(1); err=nf90_get_var(ncid, iVarId, parVector, start=(/fHRU/), count=(/1/) )
     case(2); err=nf90_get_var(ncid, iVarId, parVector, start=(/fHRU,1/), count=(/1,nSoil_file/) )
     case default; err=20; message=trim(message)//'unexpected number of dimensions for parameter '//trim(parName)
    end select

    ! error check for the parameter read
    if(err/=nf90_noerr)then; message=trim(message)//trim(cmessage); return; end if

    ! populate parameter structures
    select case(nDims)
     case(1); mparStruct%gru(iGRU)%hru(localHRU_ix)%var(ixParam)%dat(:) = parVector(1)  ! also distributes scalar across depth dimension
     case(2); mparStruct%gru(iGRU)%hru(localHRU_ix)%var(ixParam)%dat(:) = parVector(:)
     case default; err=20; message=trim(message)//'unexpected number of dimensions for parameter '//trim(parName)
    end select

   end do  ! looping through HRUs

   ! deallocate space for model parameters
   deallocate(parVector,stat=err)
   if(err/=0)then
    message=trim(message)//'problem deallocating space for parameter vector'
    err=20; return
   endif

  ! **********************************************************************************************
  ! * read the basin parameters
  ! **********************************************************************************************

  ! get the basin parameters
  else

    ! get the parameter index
    ixParam = get_ixBpar(trim(parName))

    ! allow extra variables in the file that are not used
    if(ixParam==integerMissing) cycle

    ! check that basin parameters can be indexed by GRU
    if(nGRU_file==integerMissing)then
      message=trim(message)//'basin parameter requires GRU dimension in parameter file'
      err=20; return
    endif

    ! allocate space for basin parameter data
    allocate(parVector(nGRU_file),stat=err)
    if(err/=0)then
      message=trim(message)//'problem allocating space for parameter vector'
      err=20; return
    endif

    ! read parameter data
    err=nf90_get_var(ncid,iVarId,parVector)
    if(err/=nf90_noerr)then; message=trim(message)//trim(cmessage); return; endif

    ! populate basin parameters for GRUs assigned to this rank
    do iGRU=1,nGRU_local
    
      if(has_gru_id)then
        bparStruct%gru(iGRU)%var(ixParam) = parVector(index_to_grunc(iGRU))
      else
        if(gru_struc(iGRU)%gru_nc<1 .or. gru_struc(iGRU)%gru_nc>nGRU_file)then
          message=trim(message)//'GRU index is inconsistent with parameter file'
          err=20; return
        endif
        bparStruct%gru(iGRU)%var(ixParam) = parVector(gru_struc(iGRU)%gru_nc)
      endif
    
    enddo

    ! deallocate space for model parameters
    deallocate(parVector,stat=err)
    if(err/=0)then
     message=trim(message)//'problem deallocating space for parameter vector'
     err=20; return
    endif

  endif  ! reading the basin parameters

 end do ! (looping through the parameters in the NetCDF file)

 ! **********************************************************************************************
 ! * finalize 
 ! **********************************************************************************************

 ! deallocate temporary arrays
 if(allocated(hruId))           deallocate(hruId)
 if(allocated(gruId))           deallocate(gruId)
 if(allocated(index_to_hrunc))  deallocate(index_to_hrunc)
 if(allocated(index_to_grunc))  deallocate(index_to_grunc)

 ! close the NetCDF file
 call nc_file_close(ncid,err,cmessage)
 if(err/=nf90_noerr)then
  message=trim(message)//'problem closing parameter file '//trim(infile)//': '//trim(cmessage)
  err=20; return
 endif

 end subroutine read_param

end module summa_read_param_module

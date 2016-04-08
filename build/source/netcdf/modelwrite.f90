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

module modelwrite_module
USE nrtype
USE netcdf
USE netcdf_util_module,only:netcdf_err                    ! netcdf error handling function
implicit none
private
public::writeForce
!public::writeAttrb
public::writeParm
public::writeModel
public::writeBasin
! define dimension lengths
integer(i4b),parameter      :: maxSpectral=2              ! maximum number of spectral bands
contains

 ! **********************************************************************************************************
 ! public subroutine writeParm: write model parameters
 ! **********************************************************************************************************
 subroutine writeParm(iHRU,dat,meta,err,message)
 USE data_types,only:var_info,var_d,var_i        ! metadata structure type
 USE var_lookup,only:iLookStat                   ! to index into write flag
 USE multiconst,only:integerMissing
 USE globalData,only:ncid                        ! netcdf file ids
 implicit none

 ! declare input variables
 integer(i4b)  ,intent(in)   :: iHRU             ! hydrologic response unit
 class(*)      ,intent(in)   :: dat(:)           ! local attributes
 type(var_info),intent(in)   :: meta(:)
 integer(i4b)  ,intent(out)  :: err              ! error code
 character(*)  ,intent(out)  :: message          ! error message
 ! local variables
 integer(i4b)                :: iVar             ! loop through variables
 integer(i4b)  ,parameter    :: modelTime=1      ! these particular data are only output in the timestep file

 ! initialize error control
 err=0;message="f-writeParm/"

 ! loop through local column model parameters
 do iVar = 1,size(meta)

  ! check that the variable is desired
  if (.not.meta(iVar)%statFlag(iLookStat%inst)) cycle

  ! initialize message
  message=trim(message)//trim(meta(iVar)%varName)//'/'
print*,'here',iHRU

  ! write data
  if (iHRU.ne.integerMissing) then
   selecttype (dat)
    typeis (integer)
print*, 'i,HRU',(meta(iVar)%varName),dat(iVar)
     err = nf90_put_var(ncid(modelTime),meta(iVar)%ncVarID(iLookStat%inst),(/dat(iVar)/),start=(/iHRU/),count=(/1/))
    typeis (real(dp))
print*, 'd,HRU',(meta(iVar)%varName),dat(iVar)
     err = nf90_put_var(ncid(modelTime),meta(iVar)%ncVarID(iLookStat%inst),(/dat(iVar)/),start=(/iHRU/),count=(/1/))
   endselect
 else
   selecttype (dat)
    typeis (real(dp))
print*, 'd,noHRU',(meta(iVar)%varName),dat(iVar)
     err = nf90_put_var(ncid(modelTime),meta(iVar)%ncVarID(iLookStat%inst),(/dat(iVar)/),start=(/1/),count=(/1/))
   endselect
  endif
  call netcdf_err(err,message); if (err/=0) return

  ! re-initialize message
  message="f-writeParm/"
 end do  ! looping through local column model parameters

 end subroutine writeParm

 ! **********************************************************************************************************
 ! public subroutine writeParam: write model forcing data
 ! **********************************************************************************************************
 subroutine writeForce(forc_data,iHRU,istep,err,message)
 USE data_types,only:var_d                       ! data structure: x%var(:)    (dp)
 USE globalData,only:forc_meta                   ! forcing metadata
 USE var_lookup,only:iLookFORCE                  ! identifies element of the forcing structure
 USE globalData,only:ncid                        ! id of netcdf output file
 implicit none
 ! declare dummy variables
 type(var_d),  intent(in)    :: forc_data        ! forcing data structure
 integer(i4b), intent(in)    :: iHRU             ! hydrologic response unit
 integer(i4b), intent(in)    :: istep            ! model time step
 integer(i4b),intent(out)    :: err              ! error code
 character(*),intent(out)    :: message          ! error message
 ! local variables
 integer(i4b)                :: iforce           ! loop through model forcing variables
 integer(i4b)                :: iVarId           ! variable ID
 integer(i4b)                 :: nncid
 ! initialize error control
 err=0;message="f-writeForce/"
 nncid = ncid(1)

 ! write the time coordinate variable
 if(iHRU == 1)then
  message=trim(message)//'writeTime/'
  !print*, 'iHRU, istep, forc_data%var(iLookFORCE%time) = ', iHRU, istep, forc_data%var(iLookFORCE%time)
  err = nf90_inq_varid(nncid,'time',iVarId); call netcdf_err(err,message); if (err/=0) return
  err = nf90_put_var(nncid,iVarId,(/forc_data%var(iLookFORCE%time)/),start=(/istep/),count=(/1/))
  call netcdf_err(err,message); if (err/=0) return
  message="f-writeForce/"
 endif

 ! loop through model forcing variables
 do iforce=1,size(forc_meta)
  ! ignore the time variable (used as a coordinate variable above)
  if(forc_meta(iforce)%varname == 'time') cycle
  ! check that the variable is desired
  if (.not.forc_meta(iforce)%v_write) cycle
  ! initialize message
  message=trim(message)//trim(forc_meta(iforce)%varname)//'/'
  ! get variable ID
  err = nf90_inq_varid(nncid,trim(forc_meta(iforce)%varname),iVarId)
  call netcdf_err(err,message); if (err/=0) return
  ! write data
  err = nf90_put_var(nncid,iVarId,(/forc_data%var(iforce)/),start=(/iHRU,istep/),count=(/1,1/))
  call netcdf_err(err,message); if (err/=0) return
 end do  ! looping through forcing data variables

 end subroutine writeForce

 ! **********************************************************************************************************
 ! public subroutine writeModel: write local column model variables
 ! **********************************************************************************************************
 subroutine writeModel(indx_data,metaStruct,dataStruct,iHRU,istep,err,message)
 USE var_lookup,only:iLookVarType                ! look up structure for variable typed
 USE get_ixName_module,only:get_varTypeName      ! to access type strings for error messages
 USE var_lookup,only:iLookINDEX                  ! identifies element of the index structure
 USE data_types,only:var_ilength                 ! data structure: x%var(:)%dat (i4b)
 USE data_types,only:var_dlength                 ! data structure: x%var(:)%dat (dp)
 USE data_types,only:var_info                    ! metadata structure
 USE globalData,only:ncid                        ! id of netcdf output file
 implicit none
 ! input variables
 type(var_ilength),intent(in)  :: indx_data      ! model indices
 type(var_info),intent(in)     :: metaStruct(:)  ! metadata structure
 class(*),intent(in)           :: dataStruct     ! data structure
 integer(i4b), intent(in)      :: iHRU           ! hydrologic response unit
 integer(i4b), intent(in)      :: istep          ! model time step
 ! output variables
 integer(i4b),intent(out)      :: err            ! error code
 character(*),intent(out)      :: message        ! error message
 ! local variables
 integer(i4b)                  :: ivar           ! variable index
 integer(i4b)                  :: iVarId         ! variable ID
 integer(i4b)                 :: nncid
 ! initialize error control
 err=0;message="writeModel/"
 nncid = ncid(1)

 ! associate local variables with information in the data structures
 associate(&
 ! model layers
 nSoil             => indx_data%var(iLookINDEX%nSoil)%dat(1)             ,&          
 nSnow             => indx_data%var(iLookINDEX%nSnow)%dat(1)             ,&  
 nLayers           => indx_data%var(iLookINDEX%nLayers)%dat(1)           ,&   
 ! model indices
 midSnowStartIndex => indx_data%var(iLookINDEX%midSnowStartIndex)%dat(1) ,&             
 midSoilStartIndex => indx_data%var(iLookINDEX%midSoilStartIndex)%dat(1) ,&             
 midTotoStartIndex => indx_data%var(iLookINDEX%midTotoStartIndex)%dat(1) ,&             
 ifcSnowStartIndex => indx_data%var(iLookINDEX%ifcSnowStartIndex)%dat(1) ,&             
 ifcSoilStartIndex => indx_data%var(iLookINDEX%ifcSoilStartIndex)%dat(1) ,&             
 ifcTotoStartIndex => indx_data%var(iLookINDEX%ifcTotoStartIndex)%dat(1)  &
 )  ! (associating local variables with information in the data structures)

 ! loop through model variables
 do ivar=1,size(metaStruct)

  ! check that the variable is desired
  if (.not.metaStruct(ivar)%v_write .or. trim(metaStruct(ivar)%varname)=='unknown') cycle

  ! get variable ID
  err = nf90_inq_varid(nncid,trim(metaStruct(ivar)%varname),iVarId)
  call netcdf_err(err,message); if (err/=0) return

  ! select data type
  select type(dataStruct)

   ! write integer data
   type is (var_ilength)
    err = nf90_put_var(nncid,iVarId,dataStruct%var(ivar)%dat,start=(/iHRU,istep/),count=(/1,1/))

   ! write double precision data
   type is (var_dlength)

    ! write model data 
    select case(metaStruct(ivar)%vartype)
     case(iLookVarType%scalarv); err = nf90_put_var(nncid,iVarId,dataStruct%var(ivar)%dat,start=(/iHRU,istep/),count=(/1,1/))
     case(iLookVarType%wLength); err = nf90_put_var(nncid,iVarId,dataStruct%var(ivar)%dat,start=(/iHRU,1,istep/),count=(/1,maxSpectral,1/))
     case(iLookVarType%midSnow); err = nf90_put_var(nncid,iVarId,dataStruct%var(ivar)%dat,start=(/iHRU,midSnowStartIndex/),count=(/1,nSnow/))
     case(iLookVarType%midSoil); err = nf90_put_var(nncid,iVarId,dataStruct%var(ivar)%dat,start=(/iHRU,midSoilStartIndex/),count=(/1,nSoil/))
     case(iLookVarType%midToto); err = nf90_put_var(nncid,iVarId,dataStruct%var(ivar)%dat,start=(/iHRU,midTotoStartIndex/),count=(/1,nLayers/))
     case(iLookVarType%ifcSnow); err = nf90_put_var(nncid,iVarId,dataStruct%var(ivar)%dat,start=(/iHRU,ifcSnowStartIndex/),count=(/1,nSnow+1/))
     case(iLookVarType%ifcSoil); err = nf90_put_var(nncid,iVarId,dataStruct%var(ivar)%dat,start=(/iHRU,ifcSoilStartIndex/),count=(/1,nSoil+1/))
     case(iLookVarType%ifcToto); err = nf90_put_var(nncid,iVarId,dataStruct%var(ivar)%dat,start=(/iHRU,ifcTotoStartIndex/),count=(/1,nLayers+1/))
     case default
      err=40; message=trim(message)//"unknownVariableType[name='"//trim(metaStruct(ivar)%varname)//"'; &
                                     &type='"//trim(get_varTypeName(metaStruct(ivar)%vartype))//"']"; return
    endselect ! selecting the variable type

   ! check that we found the data type
   class default; err=20; message=trim(message)//'unable to identify the data structure'; return

  end select  ! selecting the data structure

  ! process error code
  call netcdf_err(err,message); if (err/=0) return

 end do  ! looping through model variables

 ! end associating local variables with information in the data structures
 end associate

 end subroutine writeModel


 ! **********************************************************************************************************
 ! public subroutine writeBasin: write basin-average variables
 ! **********************************************************************************************************
 subroutine writeBasin(bvar_data,istep,err,message)
 USE var_lookup,only:iLookVarType              ! look up structure for variable typed
 USE get_ixName_module,only:get_varTypeName    ! to access type strings for error messages
 USE data_types,only:var_dlength               ! data structure: x%var(:)%dat (dp)
 USE globalData,only:bvar_meta                 ! metadata structures
 USE var_lookup,only:iLookINDEX                ! identifies element of the index structure
 USE globalData,only:ncid                      ! id of netcdf output file
 implicit none
 ! declare dummy variables
 type(var_dlength),intent(in) :: bvar_data     ! model variables for the local basin
 integer(i4b), intent(in)     :: istep         ! model time step
 integer(i4b),intent(out)     :: err           ! error code
 character(*),intent(out)     :: message       ! error message
 ! local variables
 integer(i4b)                 :: imodel        ! loop through model variables
 integer(i4b)                 :: iVarId        ! variable ID
 integer(i4b)                 :: nncid
 ! initialize error control
 err=0;message="f-writeModel/"
 nncid = ncid(1)

 ! loop through model variables
 ! ----------------------------
 do imodel=1,size(bvar_meta)
  ! check that the variable is desired
  if (.not.bvar_meta(imodel)%v_write) cycle
  ! initialize message
  message=trim(message)//trim(bvar_meta(imodel)%varname)//'/'
  ! get variable ID
  err = nf90_inq_varid(nncid,trim(bvar_meta(imodel)%varname),iVarId)
  call netcdf_err(err,message); if (err/=0) return
  ! write data
  select case(bvar_meta(imodel)%vartype)
   case(iLookVarType%scalarv); err = nf90_put_var(nncid,iVarId,bvar_data%var(imodel)%dat,start=(/istep/),count=(/1/))
   case(iLookVarType%routing)
    if(istep==1)then
     err = nf90_put_var(nncid,iVarId,bvar_data%var(imodel)%dat,start=(/1/),count=(/1000/))
    endif
   case default
    err=40; message=trim(message)//"unknownVariableType[name='"//trim(bvar_meta(imodel)%varname)//"'; &
                                   &type='"//trim(get_varTypeName(bvar_meta(imodel)%vartype))//"']"; return
  endselect
  call netcdf_err(err,message); if (err/=0) return
  message="f-writeBasin/"
 end do  ! looping through model variables

 end subroutine writeBasin

end module modelwrite_module

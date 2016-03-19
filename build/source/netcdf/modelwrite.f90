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
implicit none
private
public::writeForce
public::writeAttrb
public::writeParam
public::writeModel
public::writeBasin
! define dimension lengths
integer(i4b),parameter      :: maxSpectral=2              ! maximum number of spectral bands
contains


 ! **********************************************************************************************************
 ! public subroutine writeAttrb: write local attributes
 ! **********************************************************************************************************
 subroutine writeAttrb(fileout,iHRU,attr_data,type_data,err,message)
 USE globalData,only:attr_meta                   ! metadata for local attributes
 USE globalData,only:type_meta                   ! metadata for local classification of veg, soil, etc.
 implicit none
 ! declare input variables
 character(*), intent(in)    :: fileout          ! output file
 integer(i4b), intent(in)    :: iHRU             ! hydrologic response unit
 real(dp),     intent(in)    :: attr_data(:)     ! local attributes
 integer(i4b), intent(in)    :: type_data(:)     ! local classification of veg, soil, etc.
 ! declare output variables
 integer(i4b),intent(out)    :: err              ! error code
 character(*),intent(out)    :: message          ! error message
 ! local variables
 integer(i4b)                :: ncid             ! NetCDF file ID
 integer(i4b)                :: iVar             ! loop through variables
 integer(i4b)                :: iVarId           ! variable ID
 ! initialize error control
 err=0;message="f-writeAttrb/"

 ! open NetCDF file
 err = nf90_open(trim(fileout),nf90_write,ncid)
 call netcdf_err(err,message); if (err/=0) return

 ! loop through local attributes
 do iVar=1,size(attr_meta)
  ! check that the variable is desired
  if (.not.attr_meta(iVar)%v_write) cycle
  ! initialize message
  message=trim(message)//trim(attr_meta(iVar)%varname)//'/'
  ! get variable ID
  err = nf90_inq_varid(ncid,trim(attr_meta(iVar)%varname),iVarId)
  call netcdf_err(err,message); if (err/=0) return
  ! write data
  err = nf90_put_var(ncid,iVarId,(/attr_data(iVar)/),start=(/iHRU/),count=(/1/))
  call netcdf_err(err,message); if (err/=0) return
  ! re-initialize message
  message="f-writeAttrb/"
 end do  ! looping through local attributes

 ! loop through local classification of veg, soil, etc.
 do iVar=1,size(type_meta)
  ! check that the variable is desired
  if (.not.type_meta(iVar)%v_write) cycle
  ! initialize message
  message=trim(message)//trim(type_meta(iVar)%varname)//'/'
  ! get variable ID
  err = nf90_inq_varid(ncid,trim(type_meta(iVar)%varname),iVarId)
  call netcdf_err(err,message); if (err/=0) return
  ! write data
  err = nf90_put_var(ncid,iVarId,(/type_data(iVar)/),start=(/iHRU/),count=(/1/))
  call netcdf_err(err,message); if (err/=0) return
  ! re-initialize message
  message="f-writeAttrb/"
 end do  ! looping through local classification of veg, soil, etc.

 ! close output file
 err = nf90_close(ncid); call netcdf_err(err,message); if (err/=0) return
 end subroutine writeAttrb


 ! **********************************************************************************************************
 ! public subroutine writeParam: write model parameters
 ! **********************************************************************************************************
 subroutine writeParam(fileout,iHRU,mpar_data,bpar_data,err,message)
 USE globalData,only:mpar_meta                   ! metadata for local-column model parameter structures
 USE globalData,only:bpar_meta                   ! metadata for basin-average model parameter structures
 implicit none
 ! declare input variables
 character(*), intent(in)    :: fileout          ! output file
 integer(i4b), intent(in)    :: iHRU             ! hydrologic response unit
 real(dp),     intent(in)    :: mpar_data(:)     ! vector of local-column model parameters
 real(dp),     intent(in)    :: bpar_data(:)     ! vector of basin-aaverage model parameters
 ! declare output variables
 integer(i4b),intent(out)    :: err              ! error code
 character(*),intent(out)    :: message          ! error message
 ! local variables
 integer(i4b)                :: ncid             ! NetCDF file ID
 integer(i4b)                :: ipar             ! loop through model parameters
 integer(i4b)                :: iVarId           ! variable ID
 ! initialize error control
 err=0;message="f-writeParam/"

 ! open NetCDF file
 err = nf90_open(trim(fileout),nf90_write,ncid)
 call netcdf_err(err,message); if (err/=0) return

 ! loop through local column model parameters
 do ipar=1,size(mpar_meta)
  ! check that the variable is desired
  if (.not.mpar_meta(ipar)%v_write) cycle
  ! initialize message
  message=trim(message)//trim(mpar_meta(ipar)%varname)//'/'
  ! get variable ID
  err = nf90_inq_varid(ncid,trim(mpar_meta(ipar)%varname),iVarId)
  call netcdf_err(err,message); if (err/=0) return
  ! write data
  err = nf90_put_var(ncid,iVarId,(/mpar_data(ipar)/),start=(/iHRU/),count=(/1/))
  call netcdf_err(err,message); if (err/=0) return
  ! re-initialize message
  message="f-writeParam/"
 end do  ! looping through local column model parameters

 ! loop through basin-average model parameters
 do ipar=1,size(bpar_meta)
  ! check that the variable is desired
  if (.not.bpar_meta(ipar)%v_write) cycle
  ! initialize message
  message=trim(message)//trim(bpar_meta(ipar)%varname)//'/'
  ! get variable ID
  err = nf90_inq_varid(ncid,trim(bpar_meta(ipar)%varname),iVarId)
  call netcdf_err(err,message); if (err/=0) return
  ! write data
  err = nf90_put_var(ncid,iVarId,(/bpar_data(ipar)/),start=(/1/),count=(/1/))
  call netcdf_err(err,message); if (err/=0) return
  ! re-initialize message
  message="f-writeParam/"
 end do  ! looping through basin-average model parameters

 ! close output file
 err = nf90_close(ncid); call netcdf_err(err,message); if (err/=0) return
 end subroutine writeParam


 ! **********************************************************************************************************
 ! public subroutine writeParam: write model forcing data
 ! **********************************************************************************************************
 subroutine writeForce(fileout,forc_data,iHRU,istep,err,message)
 USE data_types,only:var_d                       ! data structure: x%var(:)    (dp)
 USE globalData,only:forc_meta                   ! forcing metadata
 USE var_lookup,only:iLookFORCE                  ! identifies element of the forcing structure
 implicit none
 ! declare dummy variables
 character(*), intent(in)    :: fileout          ! output file
 type(var_d),  intent(in)    :: forc_data        ! forcing data structure
 integer(i4b), intent(in)    :: iHRU             ! hydrologic response unit
 integer(i4b), intent(in)    :: istep            ! model time step
 integer(i4b),intent(out)    :: err              ! error code
 character(*),intent(out)    :: message          ! error message
 ! local variables
 integer(i4b)                :: ncid             ! NetCDF file ID
 integer(i4b)                :: iforce           ! loop through model forcing variables
 integer(i4b)                :: iVarId           ! variable ID
 ! initialize error control
 err=0;message="f-writeForce/"

 ! open NetCDF file
 err = nf90_open(trim(fileout),nf90_write,ncid)
 call netcdf_err(err,message); if (err/=0) return

 ! write the time coordinate variable
 if(iHRU == 1)then
  message=trim(message)//'writeTime/'
  !print*, 'iHRU, istep, forc_data%var(iLookFORCE%time) = ', iHRU, istep, forc_data%var(iLookFORCE%time)
  err = nf90_inq_varid(ncid,'time',iVarId); call netcdf_err(err,message); if (err/=0) return
  err = nf90_put_var(ncid,iVarId,(/forc_data%var(iLookFORCE%time)/),start=(/istep/),count=(/1/))
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
  err = nf90_inq_varid(ncid,trim(forc_meta(iforce)%varname),iVarId)
  call netcdf_err(err,message); if (err/=0) return
  ! write data
  err = nf90_put_var(ncid,iVarId,(/forc_data%var(iforce)/),start=(/iHRU,istep/),count=(/1,1/))
  call netcdf_err(err,message); if (err/=0) return
 end do  ! looping through forcing data variables

 ! close output file
 message="f-writeForce/"
 err = nf90_close(ncid); call netcdf_err(err,message); if (err/=0) return
 end subroutine writeForce

 ! **********************************************************************************************************
 ! public subroutine writeModel: write local column model variables
 ! **********************************************************************************************************
 subroutine writeModel(fileout,indx_data,metaStruct,dataStruct,iHRU,istep,err,message)
 USE var_lookup,only:iLookINDEX                  ! identifies element of the index structure
 USE data_types,only:var_ilength                 ! data structure: x%var(:)%dat (i4b)
 USE data_types,only:var_dlength                 ! data structure: x%var(:)%dat (dp)
 USE data_types,only:var_info                    ! metadata structure
 implicit none
 ! input variables
 character(*), intent(in)      :: fileout        ! output file
 type(var_ilength),intent(in)  :: indx_data      ! model indices
 type(var_info),intent(in)     :: metaStruct(:)  ! metadata structure
 class(*),intent(in)           :: dataStruct     ! data structure
 integer(i4b), intent(in)      :: iHRU           ! hydrologic response unit
 integer(i4b), intent(in)      :: istep          ! model time step
 ! output variables
 integer(i4b),intent(out)      :: err            ! error code
 character(*),intent(out)      :: message        ! error message
 ! local variables
 integer(i4b)                  :: ncid           ! NetCDF file ID
 integer(i4b)                  :: ivar           ! variable index
 integer(i4b)                  :: iVarId         ! variable ID
 ! initialize error control
 err=0;message="writeModel/"

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

 ! open NetCDF file
 err = nf90_open(trim(fileout),nf90_write,ncid)
 call netcdf_err(err,message); if (err/=0) return

 ! loop through model variables
 do ivar=1,size(metaStruct)

  ! check that the variable is desired
  if (.not.metaStruct(ivar)%v_write .or. trim(metaStruct(ivar)%varname)=='unknown') cycle

  ! get variable ID
  err = nf90_inq_varid(ncid,trim(metaStruct(ivar)%varname),iVarId)
  call netcdf_err(err,message); if (err/=0) return

  ! select data type
  select type(dataStruct)

   ! write integer data
   type is (var_ilength)
    err = nf90_put_var(ncid,iVarId,dataStruct%var(ivar)%dat,start=(/iHRU,istep/),count=(/1,1/))

   ! write double precision data
   type is (var_dlength)

    ! write model data 
    select case(trim(metaStruct(ivar)%vartype))
     case('scalarv'); err = nf90_put_var(ncid,iVarId,dataStruct%var(ivar)%dat,start=(/iHRU,istep/),count=(/1,1/))
     case('wLength'); err = nf90_put_var(ncid,iVarId,dataStruct%var(ivar)%dat,start=(/iHRU,1,istep/),count=(/1,maxSpectral,1/))
     case('midSnow'); err = nf90_put_var(ncid,iVarId,dataStruct%var(ivar)%dat,start=(/iHRU,midSnowStartIndex/),count=(/1,nSnow/))
     case('midSoil'); err = nf90_put_var(ncid,iVarId,dataStruct%var(ivar)%dat,start=(/iHRU,midSoilStartIndex/),count=(/1,nSoil/))
     case('midToto'); err = nf90_put_var(ncid,iVarId,dataStruct%var(ivar)%dat,start=(/iHRU,midTotoStartIndex/),count=(/1,nLayers/))
     case('ifcSnow'); err = nf90_put_var(ncid,iVarId,dataStruct%var(ivar)%dat,start=(/iHRU,ifcSnowStartIndex/),count=(/1,nSnow+1/))
     case('ifcSoil'); err = nf90_put_var(ncid,iVarId,dataStruct%var(ivar)%dat,start=(/iHRU,ifcSoilStartIndex/),count=(/1,nSoil+1/))
     case('ifcToto'); err = nf90_put_var(ncid,iVarId,dataStruct%var(ivar)%dat,start=(/iHRU,ifcTotoStartIndex/),count=(/1,nLayers+1/))
     case default
      err=40; message=trim(message)//"unknownVariableType[name='"//trim(metaStruct(ivar)%varname)//"'; &
                                     &type='"//trim(metaStruct(ivar)%vartype)//"']"; return
    endselect ! selecting the variable type

   ! check that we found the data type
   class default; err=20; message=trim(message)//'unable to identify the data structure'

  end select  ! selecting the data structure

  ! process error code
  call netcdf_err(err,message); if (err/=0) return

 end do  ! looping through model variables

 ! close output file
 err = nf90_close(ncid); call netcdf_err(err,message); if (err/=0) return

 ! end associating local variables with information in the data structures
 end associate

 end subroutine writeModel


 ! **********************************************************************************************************
 ! public subroutine writeBasin: write basin-average variables
 ! **********************************************************************************************************
 subroutine writeBasin(fileout,bvar_data,istep,err,message)
 USE data_types,only:var_dlength               ! data structure: x%var(:)%dat (dp)
 USE globalData,only:bvar_meta                 ! metadata structures
 USE var_lookup,only:iLookINDEX                ! identifies element of the index structure
 implicit none
 ! declare dummy variables
 character(*), intent(in)     :: fileout       ! output file
 type(var_dlength),intent(in) :: bvar_data     ! model variables for the local basin
 integer(i4b), intent(in)     :: istep         ! model time step
 integer(i4b),intent(out)     :: err           ! error code
 character(*),intent(out)     :: message       ! error message
 ! local variables
 integer(i4b)                 :: ncid          ! NetCDF file ID
 integer(i4b)                 :: imodel        ! loop through model variables
 integer(i4b)                 :: iVarId        ! variable ID
 ! initialize error control
 err=0;message="f-writeModel/"

 ! open NetCDF file
 err = nf90_open(trim(fileout),nf90_write,ncid)
 call netcdf_err(err,message); if (err/=0) return

 ! loop through model variables
 ! ----------------------------
 do imodel=1,size(bvar_meta)
  ! check that the variable is desired
  if (.not.bvar_meta(imodel)%v_write) cycle
  ! initialize message
  message=trim(message)//trim(bvar_meta(imodel)%varname)//'/'
  ! get variable ID
  err = nf90_inq_varid(ncid,trim(bvar_meta(imodel)%varname),iVarId)
  call netcdf_err(err,message); if (err/=0) return
  ! write data
  select case(trim(bvar_meta(imodel)%vartype))
   case('scalarv'); err = nf90_put_var(ncid,iVarId,bvar_data%var(imodel)%dat,start=(/istep/),count=(/1/))
   case('routing')
    if(istep==1)then
     err = nf90_put_var(ncid,iVarId,bvar_data%var(imodel)%dat,start=(/1/),count=(/1000/))
    endif
   case default
    err=40; message=trim(message)//"unknownVariableType[name='"//trim(bvar_meta(imodel)%varname)//"'; &
                                   &type='"//trim(bvar_meta(imodel)%vartype)//"']"; return
  endselect
  call netcdf_err(err,message); if (err/=0) return
  message="f-writeBasin/"
 end do  ! looping through model variables

 ! close output file
 err = nf90_close(ncid); call netcdf_err(err,message); if (err/=0) return
 end subroutine writeBasin


 ! **********************************************************************************************************
 ! private subroutine netcdf_err: error control
 ! **********************************************************************************************************
 subroutine netcdf_err(err,message)
 ! used to handle errors for NetCDF calls
 implicit none
 ! declare dummies
 integer(i4b), intent(inout)   :: err
 character(*), intent(inout)   :: message
 ! start procedure here
 if (err/=nf90_noerr) then
  message=trim(message)//"["//trim(nf90_strerror(err))//"]"
  err=200
 else
  err=0
 endif
 end subroutine netcdf_err


end module modelwrite_module

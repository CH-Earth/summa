module modelwrite_module
USE nrtype
USE netcdf
implicit none
private
public::writeForce
public::writeAttrb
public::writeParam
public::writeModel
! define dimension names
character(len=32),parameter :: parSet_DimName='parSet'                 ! dimension name for the parameter sets
character(len=32),parameter :: timestep_DimName='time'                 ! dimension name for the time step (unlimited)
character(len=32),parameter :: midSnowAndTime_DimName='midSnowAndTime' ! dimension name for midSnow-time (unlimited)
character(len=32),parameter :: midSoilAndTime_DimName='midSoilAndTime' ! dimension name for midSoil-time (unlimited)
character(len=32),parameter :: midTotoAndTime_DimName='midTotoAndTime' ! dimension name for midToto-time (unlimited)
character(len=32),parameter :: ifcSnowAndTime_DimName='ifcSnowAndTime' ! dimension name for ifcSnow-time (unlimited)
character(len=32),parameter :: ifcSoilAndTime_DimName='ifcSoilAndTime' ! dimension name for ifcSoil-time (unlimited)
character(len=32),parameter :: ifcTotoAndTime_DimName='ifcTotoAndTime' ! dimension name for ifcToto-time (unlimited)
contains

 ! **********************************************************************************************************
 ! new subroutine: write local attributes
 ! **********************************************************************************************************
 subroutine writeAttrb(fileout,err,message)
 USE data_struc,only:attr_data,attr_meta                   ! local attributes
 USE data_struc,only:type_data,type_meta                   ! local classification of veg, soil, etc.
 implicit none
 ! declare dummy variables
 character(*), intent(in)    :: fileout                    ! output file
 integer(i4b),intent(out)    :: err                        ! error code
 character(*),intent(out)    :: message                    ! error message
 ! local variables
 integer(i4b)                :: ncid                       ! NetCDF file ID
 integer(i4b)                :: iVar                       ! loop through variables
 integer(i4b)                :: iVarId                     ! variable ID
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
  err = nf90_put_var(ncid,iVarId,(/attr_data%var(iVar)/),start=(/1/),count=(/1/))
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
  err = nf90_put_var(ncid,iVarId,(/type_data%var(iVar)/),start=(/1/),count=(/1/))
  call netcdf_err(err,message); if (err/=0) return
  ! re-initialize message
  message="f-writeAttrb/"
 end do  ! looping through local classification of veg, soil, etc.

 ! close output file
 err = nf90_close(ncid); call netcdf_err(err,message); if (err/=0) return
 end subroutine writeAttrb



 ! **********************************************************************************************************
 ! new subroutine: write model parameters
 ! **********************************************************************************************************
 subroutine writeParam(fileout,iParSet,err,message)
 USE data_struc,only:mpar_data,mpar_meta                   ! model parameter structures
 implicit none
 ! declare dummy variables
 character(*), intent(in)    :: fileout                    ! output file
 integer(i4b), intent(in)    :: iParSet                    ! model parameter set
 integer(i4b),intent(out)    :: err                        ! error code
 character(*),intent(out)    :: message                    ! error message
 ! local variables
 integer(i4b)                :: ncid                       ! NetCDF file ID
 integer(i4b)                :: ipar                       ! loop through model parameters
 integer(i4b)                :: iVarId                     ! variable ID
 ! initialize error control
 err=0;message="f-writeParam/"

 ! open NetCDF file
 err = nf90_open(trim(fileout),nf90_write,ncid)
 call netcdf_err(err,message); if (err/=0) return

 ! loop through model parameters
 do ipar=1,size(mpar_meta)
  ! check that the variable is desired
  if (.not.mpar_meta(ipar)%v_write) cycle
  ! initialize message
  message=trim(message)//trim(mpar_meta(ipar)%varname)//'/'
  ! get variable ID
  err = nf90_inq_varid(ncid,trim(mpar_meta(ipar)%varname),iVarId)
  call netcdf_err(err,message); if (err/=0) return
  ! write data
  err = nf90_put_var(ncid,iVarId,(/mpar_data%var(ipar)/),start=(/iParSet/),count=(/1/))
  call netcdf_err(err,message); if (err/=0) return
  ! re-initialize message
  message="f-writeParam/"
 end do  ! looping through model parameters

 ! close output file
 err = nf90_close(ncid); call netcdf_err(err,message); if (err/=0) return
 end subroutine writeParam


 ! **********************************************************************************************************
 ! new subroutine: write model forcing data
 ! **********************************************************************************************************
 subroutine writeForce(fileout,istep,err,message)
 USE data_struc,only:forc_data,forc_meta                   ! forcing data structures
 USE var_lookup,only:iLookFORCE                            ! identifies element of the forcing structure
 implicit none
 ! declare dummy variables
 character(*), intent(in)    :: fileout                    ! output file
 integer(i4b), intent(in)    :: istep                      ! model time step
 integer(i4b),intent(out)    :: err                        ! error code
 character(*),intent(out)    :: message                    ! error message
 ! pointers to time variables
 real(dp),pointer            :: dtime                      ! time since reference time (seconds)
 ! local variables
 integer(i4b)                :: ncid                       ! NetCDF file ID
 integer(i4b)                :: iforce                     ! loop through model forcing variables
 integer(i4b)                :: iVarId                     ! variable ID
 ! initialize error control
 err=0;message="f-writeForce/"

 ! assign pointers
 dtime => forc_data%var(iLookFORCE%time)

 ! open NetCDF file
 err = nf90_open(trim(fileout),nf90_write,ncid)
 call netcdf_err(err,message); if (err/=0) return

 ! write the time coordinate variable
 message=trim(message)//'writeTime/'
 err = nf90_inq_varid(ncid,'time',iVarId); call netcdf_err(err,message); if (err/=0) return
 err = nf90_put_var(ncid,iVarId,(/dtime/),start=(/istep/),count=(/1/))
 call netcdf_err(err,message); if (err/=0) return

 ! loop through model forcing variables
 do iforce=1,size(forc_meta)
  ! check that the variable is desired
  if (.not.forc_meta(iforce)%v_write) cycle
  ! initialize message
  message=trim(message)//trim(forc_meta(iforce)%varname)//'/'
  ! get variable ID
  err = nf90_inq_varid(ncid,trim(forc_meta(iforce)%varname),iVarId)
  call netcdf_err(err,message); if (err/=0) return
  ! write data
  err = nf90_put_var(ncid,iVarId,(/forc_data%var(iforce)/),start=(/istep/),count=(/1/))
  call netcdf_err(err,message); if (err/=0) return
 end do  ! looping through forcing data variables

 ! close output file
 message="f-writeForce/"
 err = nf90_close(ncid); call netcdf_err(err,message); if (err/=0) return
 end subroutine writeForce

 ! **********************************************************************************************************
 ! new subroutine: write model variables
 ! **********************************************************************************************************
 subroutine writeModel(fileout,ipar,istep,err,message)
 USE data_struc,only:indx_data,indx_meta                   ! index data structures
 USE data_struc,only:mvar_data,mvar_meta                   ! model data structures
 USE var_lookup,only:iLookINDEX                            ! identifies element of the index structure
 implicit none
 ! declare dummy variables
 character(*), intent(in)    :: fileout                    ! output file
 integer(i4b), intent(in)    :: ipar                       ! model parameter index
 integer(i4b), intent(in)    :: istep                      ! model time step
 integer(i4b),intent(out)    :: err                        ! error code
 character(*),intent(out)    :: message                    ! error message
 ! declare pointers to the model index variables (used to identify the appropriate write position)
 integer(i4b),pointer       :: nSnow                          ! number of snow layers
 integer(i4b),pointer       :: nSoil                          ! number of soil layers
 integer(i4b),pointer       :: nLayers                        ! total number of layers
 integer(i4b),pointer       :: midSnowStartIndex              ! start index of the midSnow vector for a given timestep
 integer(i4b),pointer       :: midSoilStartIndex              ! start index of the midSoil vector for a given timestep
 integer(i4b),pointer       :: midTotoStartIndex              ! start index of the midToto vector for a given timestep
 integer(i4b),pointer       :: ifcSnowStartIndex              ! start index of the ifcSnow vector for a given timestep
 integer(i4b),pointer       :: ifcSoilStartIndex              ! start index of the ifcSoil vector for a given timestep
 integer(i4b),pointer       :: ifcTotoStartIndex              ! start index of the ifcToto vector for a given timestep
 ! local variables
 integer(i4b)                :: ncid                       ! NetCDF file ID
 integer(i4b)                :: iindex                     ! loop through model index variables
 integer(i4b)                :: imodel                     ! loop through model variables
 integer(i4b)                :: iVarId                     ! variable ID
 ! initialize error control
 err=0;message="f-writeModel/"

 ! assign pointers to model layers
 nSoil   => indx_data%var(iLookINDEX%nSoil)%dat(1)
 nSnow   => indx_data%var(iLookINDEX%nSnow)%dat(1)
 nLayers => indx_data%var(iLookINDEX%nLayers)%dat(1)

 ! assign pointers to model indices
 midSnowStartIndex => indx_data%var(iLookINDEX%midSnowStartIndex)%dat(1)
 midSoilStartIndex => indx_data%var(iLookINDEX%midSoilStartIndex)%dat(1)
 midTotoStartIndex => indx_data%var(iLookINDEX%midTotoStartIndex)%dat(1)
 ifcSnowStartIndex => indx_data%var(iLookINDEX%ifcSnowStartIndex)%dat(1)
 ifcSoilStartIndex => indx_data%var(iLookINDEX%ifcSoilStartIndex)%dat(1)
 ifcTotoStartIndex => indx_data%var(iLookINDEX%ifcTotoStartIndex)%dat(1)

 ! open NetCDF file
 err = nf90_open(trim(fileout),nf90_write,ncid)
 call netcdf_err(err,message); if (err/=0) return

 ! loop through model index variables
 ! ----------------------------------
 do iindex=1,size(indx_meta)
  ! check that the variable is desired
  if (.not.indx_meta(iindex)%v_write) cycle
  ! initialize message
  message=trim(message)//trim(indx_meta(iindex)%varname)//'/'
  ! get variable ID
  err = nf90_inq_varid(ncid,trim(indx_meta(iindex)%varname),iVarId)
  call netcdf_err(err,message); if (err/=0) return
  ! write data
  err = nf90_put_var(ncid,iVarId,indx_data%var(iindex)%dat,start=(/ipar,istep/),count=(/1,1/))
  call netcdf_err(err,message); if (err/=0) return
  message="f-writeModel/"
 end do

 ! loop through model variables
 ! ----------------------------
 do imodel=1,size(mvar_meta)
  ! check that the variable is desired
  if (.not.mvar_meta(imodel)%v_write) cycle
  ! initialize message
  message=trim(message)//trim(mvar_meta(imodel)%varname)//'/'
  ! get variable ID
  err = nf90_inq_varid(ncid,trim(mvar_meta(imodel)%varname),iVarId)
  call netcdf_err(err,message); if (err/=0) return
  ! write data
  select case(trim(mvar_meta(imodel)%vartype))
   case('scalarv'); err = nf90_put_var(ncid,iVarId,mvar_data%var(imodel)%dat,start=(/ipar,istep/),count=(/1,1/))
   case('midSnow'); err = nf90_put_var(ncid,iVarId,mvar_data%var(imodel)%dat,start=(/ipar,midSnowStartIndex/),count=(/1,nSnow/))
   case('midSoil'); err = nf90_put_var(ncid,iVarId,mvar_data%var(imodel)%dat,start=(/ipar,midSoilStartIndex/),count=(/1,nSoil/))
   case('midToto'); err = nf90_put_var(ncid,iVarId,mvar_data%var(imodel)%dat,start=(/ipar,midTotoStartIndex/),count=(/1,nLayers/))
   case('ifcSnow'); err = nf90_put_var(ncid,iVarId,mvar_data%var(imodel)%dat,start=(/ipar,ifcSnowStartIndex/),count=(/1,nSnow+1/))
   case('ifcSoil'); err = nf90_put_var(ncid,iVarId,mvar_data%var(imodel)%dat,start=(/ipar,ifcSoilStartIndex/),count=(/1,nSoil+1/))
   case('ifcToto'); err = nf90_put_var(ncid,iVarId,mvar_data%var(imodel)%dat,start=(/ipar,ifcTotoStartIndex/),count=(/1,nLayers+1/))
   case('routing')
    if(istep==1)then
     err = nf90_put_var(ncid,iVarId,mvar_data%var(imodel)%dat,start=(/ipar,1/),count=(/1,1000/))
    endif
   case default
    err=40; message=trim(message)//"unknownVariableType[name='"//trim(mvar_meta(imodel)%varname)//"'; &
                                   &type='"//trim(mvar_meta(imodel)%vartype)//"']"; return
  endselect
  call netcdf_err(err,message); if (err/=0) return
  message="f-writeModel/"
 end do  ! looping through model variables

 ! close output file
 err = nf90_close(ncid); call netcdf_err(err,message); if (err/=0) return
 end subroutine writeModel

 ! **********************************************************************************************************
 ! subroutine X: error control
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

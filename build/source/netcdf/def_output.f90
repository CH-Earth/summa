module def_output_module
USE nrtype
USE netcdf
implicit none
private
public :: def_output
! define dimension names
character(len=32),parameter :: parSet_DimName='parSet'                 ! dimension name for the parameter sets (unlimited)
character(len=32),parameter :: timestep_DimName='time'                 ! dimension name for the time step (unlimited)
character(len=32),parameter :: midSnowAndTime_DimName='midSnowAndTime' ! dimension name for midSnow-time (unlimited)
character(len=32),parameter :: midSoilAndTime_DimName='midSoilAndTime' ! dimension name for midSoil-time (unlimited)
character(len=32),parameter :: midTotoAndTime_DimName='midTotoAndTime' ! dimension name for midToto-time (unlimited)
character(len=32),parameter :: ifcSnowAndTime_DimName='ifcSnowAndTime' ! dimension name for ifcSnow-time (unlimited)
character(len=32),parameter :: ifcSoilAndTime_DimName='ifcSoilAndTime' ! dimension name for ifcSoil-time (unlimited)
character(len=32),parameter :: ifcTotoAndTime_DimName='ifcTotoAndTime' ! dimension name for ifcToto-time (unlimited)
contains

 ! **********************************************************************************************************
 ! new subroutine: define model output file
 ! **********************************************************************************************************
 subroutine def_output(infile,err,message)
 USE data_struc,only:forc_meta,mpar_meta,mvar_meta,indx_meta  ! metadata structures
 USE data_struc,only:model_decisions
 ! declare dummy variables
 character(*), intent(in)    :: infile                       ! file suffix
 integer(i4b),intent(out)    :: err                          ! error code
 character(*),intent(out)    :: message                      ! error message
 ! local variables
 integer(i4b)                :: ivar                         ! loop through model variables
 character(len=256)          :: cmessage                     ! temporary error message
 ! initialize errors
 err=0; message="f-fuse/def_output/"
 ! **********************************************************************************************************
 ! ***** create initial file
 ! **********************************************************************************************************
 call ini_create(trim(infile),err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! **********************************************************************************************************
 ! ***** define model decisions
 ! **********************************************************************************************************
 do ivar=1,size(model_decisions)
  call put_attrib(trim(infile),model_decisions(ivar)%cOption,model_decisions(ivar)%cDecision,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do
 ! **********************************************************************************************************
 ! ***** define model parameters
 ! **********************************************************************************************************
 do ivar=1,size(mpar_meta)
  if (.not.mpar_meta(ivar)%v_write) cycle
  call def_variab(trim(infile),(/parSet_DimName/),mpar_meta(ivar),nf90_double,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do  ! looping through model parameters
 ! **********************************************************************************************************
 ! ***** define model forcing data
 ! **********************************************************************************************************
 do ivar=1,size(forc_meta)
  if (.not.forc_meta(ivar)%v_write) cycle
  call def_variab(trim(infile),(/Timestep_DimName/),forc_meta(ivar),nf90_double,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do ! looping through forcing variables
 ! **********************************************************************************************************
 ! ***** define model variables -- dimensions depend on the variable type
 ! **********************************************************************************************************
 do ivar=1,size(mvar_meta)
  if (.not.mvar_meta(ivar)%v_write) cycle
  select case(trim(mvar_meta(ivar)%vartype))
   case('scalarv'); call def_variab(trim(infile),(/parSet_DimName,Timestep_DimName/),mvar_meta(ivar),nf90_double,err,cmessage)
   case('midSnow'); call def_variab(trim(infile),(/parSet_DimName,midSnowAndTime_DimName/),mvar_meta(ivar),nf90_double,err,cmessage)
   case('midSoil'); call def_variab(trim(infile),(/parSet_DimName,midSoilAndTime_DimName/),mvar_meta(ivar),nf90_double,err,cmessage)
   case('midToto'); call def_variab(trim(infile),(/parSet_DimName,midTotoAndTime_DimName/),mvar_meta(ivar),nf90_double,err,cmessage)
   case('ifcSnow'); call def_variab(trim(infile),(/parSet_DimName,ifcSnowAndTime_DimName/),mvar_meta(ivar),nf90_double,err,cmessage)
   case('ifcSoil'); call def_variab(trim(infile),(/parSet_DimName,ifcSoilAndTime_DimName/),mvar_meta(ivar),nf90_double,err,cmessage)
   case('ifcToto'); call def_variab(trim(infile),(/parSet_DimName,ifcTotoAndTime_DimName/),mvar_meta(ivar),nf90_double,err,cmessage)
   case default; err=35; message="f-fuse/def_output/varTypeNotFound"; return
  endselect
  ! check variable definition was OK
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do ! loop through model variables
 ! **********************************************************************************************************
 ! ***** define model indices -- dimensions depend on the variable type
 ! **********************************************************************************************************
 do ivar=1,size(indx_meta)
  if (.not.indx_meta(ivar)%v_write) cycle
  select case(trim(indx_meta(ivar)%vartype))
   case('scalarv'); call def_variab(trim(infile),(/parSet_DimName,Timestep_DimName/),indx_meta(ivar),nf90_int,err,cmessage)
   case('midToto'); call def_variab(trim(infile),(/parSet_DimName,midTotoAndTime_DimName/),indx_meta(ivar),nf90_int,err,cmessage)
   case default; err=35; message="f-fuse/def_output/varTypeNotFound"; return
  endselect
  ! check variable definition was OK
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do ! loop through model variable

 end subroutine def_output

 ! **********************************************************************************************************
 ! new subroutine: initial create
 ! **********************************************************************************************************
 subroutine ini_create(infile,err,message)
 implicit none
 ! declare dummy variables
 character(*),intent(in)     :: infile                     ! filename
 integer(i4b),intent(out)    :: err                        ! error code
 character(*),intent(out)    :: message                    ! error message
 ! define local variables
 integer(i4b)                :: ncid                       ! NetCDF file ID
 integer(i4b)                :: dimID
 integer(i4b),parameter      :: maxLength=25000         ! maximum length of the variable vector
 integer(i4b),parameter      :: maxParSets=1               ! maximum number of parameter sets
 ! initialize error control
 err=0;message="f-iniCreate/"
 ! create output file
 err = nf90_create(trim(infile),nf90_classic_model,ncid)
 call netcdf_err(err,message); if (err/=0) return
 ! create parameter dimension (unlimited)
 err = nf90_def_dim(ncid, trim(parSet_DimName), maxParSets, dimId)
 call netcdf_err(err,message); if (err/=0) return
 ! create time dimension (unlimited)
 err = nf90_def_dim(ncid, trim(timestep_DimName), nf90_unlimited, dimId)
 call netcdf_err(err,message); if (err/=0) return
 ! create dimension for midSnow+time (unlimited)
 err = nf90_def_dim(ncid, trim(midSnowAndTime_DimName), maxLength, dimId)
 call netcdf_err(err,message); if (err/=0) return
 ! create dimension for midSoil+time (unlimited)
 err = nf90_def_dim(ncid, trim(midSoilAndTime_DimName), maxLength, dimId)
 call netcdf_err(err,message); if (err/=0) return
 ! create dimension for midToto+time (unlimited)
 err = nf90_def_dim(ncid, trim(midTotoAndTime_DimName), maxLength, dimId)
 call netcdf_err(err,message); if (err/=0) return
 ! create dimension for ifcSnow+time (unlimited)
 err = nf90_def_dim(ncid, trim(ifcSnowAndTime_DimName), maxLength, dimId)
 call netcdf_err(err,message); if (err/=0) return
 ! create dimension for ifcSoil+time (unlimited)
 err = nf90_def_dim(ncid, trim(ifcSoilAndTime_DimName), maxLength, dimId)
 call netcdf_err(err,message); if (err/=0) return
 ! create dimension for ifcToto+time (unlimited)
 err = nf90_def_dim(ncid, trim(ifcTotoAndTime_DimName), maxLength, dimId)
 call netcdf_err(err,message); if (err/=0) return
 ! close NetCDF file
 err = nf90_enddef(ncid); call netcdf_err(err,message); if (err/=0) return
 err = nf90_close(ncid); call netcdf_err(err,message); if (err/=0) return
 end subroutine ini_create


 ! **********************************************************************************************************
 ! new subroutine: put global attributes as character string
 ! **********************************************************************************************************
 subroutine put_attrib(infile,attname,attvalue,err,message)
 USE data_struc,only:var_info                              ! derived type for metadata
 implicit none
 ! declare dummy variables
 character(*), intent(in)   :: infile      ! filename
 character(*), intent(in)   :: attname     ! attribute name
 character(*), intent(in)   :: attvalue    ! attribute vaue
 integer(i4b),intent(out)   :: err         ! error code
 character(*),intent(out)   :: message     ! error message
 ! local variables
 integer(i4b)               :: ncid        ! NetCDF file ID
 ! initialize error control
 err=0;message="f-defAttrib/"//trim(attname)//"/"//trim(attvalue)//"/"
 ! open NetCDF file
 err = nf90_open(infile,nf90_write,ncid)
 call netcdf_err(err,message); if (err/=0) return
 ! allow re-definition of variables
 err = nf90_redef(ncid); call netcdf_err(err,message); if (err/=0) return
 ! put the attribute
 err = nf90_put_att(ncid,nf90_global,trim(attname),trim(attvalue))
 call netcdf_err(err,message); if (err/=0) return
 ! close output file
 err = nf90_enddef(ncid); call netcdf_err(err,message); if (err/=0) return
 err = nf90_close(ncid); call netcdf_err(err,message); if (err/=0) return
 end subroutine put_attrib


 ! **********************************************************************************************************
 ! new subroutine: define variables
 ! **********************************************************************************************************
 subroutine def_variab(infile,dimNames,metadata,ivtype,err,message)
 USE data_struc,only:var_info                              ! derived type for metadata
 implicit none
 ! declare dummy variables
 character(*), intent(in)   :: infile      ! filename
 character(*), intent(in)   :: dimNames(:) ! dimension namess
 type(var_info),intent(in)  :: metadata    ! metadata structure for a given variable
 integer(i4b),intent(in)    :: ivtype      ! variable type
 integer(i4b),intent(out)   :: err         ! error code
 character(*),intent(out)   :: message     ! error message
 ! local variables
 integer(i4b)               :: id          ! loop through dimensions
 integer(i4b)               :: dimIDs(size(dimNames))
 integer(i4b)               :: ncid        ! NetCDF file ID
 integer(i4b)               :: iVarId      ! variable ID
 ! initialize error control
 err=0;message="f-defVariab/"//trim(metadata%varname)//"/"

 ! open NetCDF file
 err = nf90_open(infile,nf90_write,ncid)
 call netcdf_err(err,message); if (err/=0) return
 ! allow re-definition of variables
 err = nf90_redef(ncid); call netcdf_err(err,message); if (err/=0) return

 ! define dimension IDs
 do id=1,size(dimNames)
  err=nf90_inq_dimid(ncid,trim(dimNames(id)),dimIDs(id)); call netcdf_err(err,message); if (err/=0) return
 end do

 ! define variable
 err = nf90_def_var(ncid,trim(metadata%varname),ivtype,dimIds,iVarId)
 call netcdf_err(err,message); if (err/=0) return
 ! add parameter description
 err = nf90_put_att(ncid,iVarId,'long_name',trim(metadata%vardesc))
 call netcdf_err(err,message); if (err/=0) return
 ! add parameter units
 err = nf90_put_att(ncid,iVarId,'units',trim(metadata%varunit))
 call netcdf_err(err,message); if (err/=0) return

 ! close output file
 err = nf90_enddef(ncid); call netcdf_err(err,message); if (err/=0) return
 err = nf90_close(ncid); call netcdf_err(err,message); if (err/=0) return
 end subroutine def_variab

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
  message=trim(message)//"a-OK"
  err=0
 endif
 end subroutine netcdf_err

end module def_output_module

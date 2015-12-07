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

module def_output_module
USE nrtype
USE netcdf
implicit none
private
public :: def_output
! define dimension names
character(len=32),parameter :: hru_DimName='hru'                       ! dimension name for the HRUs
character(len=32),parameter :: scalar_DimName='scalar'                 ! dimension name for scalar variables
character(len=32),parameter :: wLength_dimName='spectral_bands'        ! dimension name for the number of spectral bands
character(len=32),parameter :: timestep_DimName='time'                 ! dimension name for the time step
character(len=32),parameter :: routing_DimName='timeDelayRouting'      ! dimension name for thetime delay routing vectors
character(len=32),parameter :: midSnowAndTime_DimName='midSnowAndTime' ! dimension name for midSnow-time
character(len=32),parameter :: midSoilAndTime_DimName='midSoilAndTime' ! dimension name for midSoil-time
character(len=32),parameter :: midTotoAndTime_DimName='midTotoAndTime' ! dimension name for midToto-time
character(len=32),parameter :: ifcSnowAndTime_DimName='ifcSnowAndTime' ! dimension name for ifcSnow-time
character(len=32),parameter :: ifcSoilAndTime_DimName='ifcSoilAndTime' ! dimension name for ifcSoil-time
character(len=32),parameter :: ifcTotoAndTime_DimName='ifcTotoAndTime' ! dimension name for ifcToto-time
contains


 ! **********************************************************************************************************
 ! public subroutine def_output: define model output file
 ! **********************************************************************************************************
 subroutine def_output(nHRU,infile,err,message)
 USE data_struc,only:forc_meta,attr_meta,type_meta  ! metadata structures
 USE data_struc,only:mpar_meta,mvar_meta,indx_meta  ! metadata structures
 USE data_struc,only:bpar_meta,bvar_meta            ! metadata structures
 USE data_struc,only:model_decisions
 ! declare dummy variables
 integer(i4b), intent(in)    :: nHRU                         ! number of HRUs
 character(*), intent(in)    :: infile                       ! file suffix
 integer(i4b),intent(out)    :: err                          ! error code
 character(*),intent(out)    :: message                      ! error message
 ! local variables
 integer(i4b)                :: ivar                         ! loop through model variables
 character(len=256)          :: cmessage                     ! temporary error message
 ! initialize errors
 err=0; message="def_output/"
 ! **********************************************************************************************************
 ! ***** create initial file
 ! **********************************************************************************************************
 call ini_create(nHRU,trim(infile),err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! **********************************************************************************************************
 ! ***** define model decisions
 ! **********************************************************************************************************
 do ivar=1,size(model_decisions)
  call put_attrib(trim(infile),model_decisions(ivar)%cOption,model_decisions(ivar)%cDecision,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do
 ! **********************************************************************************************************
 ! ***** define model forcing data
 ! **********************************************************************************************************
 do ivar=1,size(forc_meta)
  if(.not.forc_meta(ivar)%v_write) cycle
  if(forc_meta(ivar)%varname == 'time')then
   call def_variab(trim(infile),(/Timestep_DimName/),forc_meta(ivar),nf90_double,err,cmessage)
  else
   call def_variab(trim(infile),(/hru_DimName,Timestep_DimName/),forc_meta(ivar),nf90_double,err,cmessage)
  endif
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do ! looping through forcing variables
 ! **********************************************************************************************************
 ! ***** define local attributes
 ! **********************************************************************************************************
 do ivar=1,size(attr_meta)
  if (.not.attr_meta(ivar)%v_write) cycle
  call def_variab(trim(infile),(/hru_DimName/),attr_meta(ivar),nf90_double,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do  ! looping through local attributes
 ! **********************************************************************************************************
 ! ***** define local classification of veg, soil, etc.
 ! **********************************************************************************************************
 do ivar=1,size(type_meta)
  if (.not.type_meta(ivar)%v_write) cycle
  call def_variab(trim(infile),(/hru_DimName/),type_meta(ivar),nf90_int,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do  ! looping through local classification of veg, soil, etc.
 ! **********************************************************************************************************
 ! ***** define local column model parameters
 ! **********************************************************************************************************
 do ivar=1,size(mpar_meta)
  if (.not.mpar_meta(ivar)%v_write) cycle
  call def_variab(trim(infile),(/hru_DimName/),mpar_meta(ivar),nf90_double,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do  ! looping through model parameters
 ! **********************************************************************************************************
 ! ***** define basin-average model parameters
 ! **********************************************************************************************************
 do ivar=1,size(bpar_meta)
  if (.not.bpar_meta(ivar)%v_write) cycle
  call def_variab(trim(infile),(/scalar_DimName/),bpar_meta(ivar),nf90_double,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do  ! looping through model parameters
 ! **********************************************************************************************************
 ! ***** define local column model variables -- dimensions depend on the variable type
 ! **********************************************************************************************************
 do ivar=1,size(mvar_meta)
  if (.not.mvar_meta(ivar)%v_write) cycle
  select case(trim(mvar_meta(ivar)%vartype))
   case('scalarv'); call def_variab(trim(infile),(/hru_DimName,Timestep_DimName/),mvar_meta(ivar),nf90_double,err,cmessage)
   case('wLength'); call def_variab(trim(infile),(/hru_DimName,wLength_DimName,Timestep_DimName/),mvar_meta(ivar),nf90_double,err,cmessage)
   case('midSnow'); call def_variab(trim(infile),(/hru_DimName,midSnowAndTime_DimName/),mvar_meta(ivar),nf90_double,err,cmessage)
   case('midSoil'); call def_variab(trim(infile),(/hru_DimName,midSoilAndTime_DimName/),mvar_meta(ivar),nf90_double,err,cmessage)
   case('midToto'); call def_variab(trim(infile),(/hru_DimName,midTotoAndTime_DimName/),mvar_meta(ivar),nf90_double,err,cmessage)
   case('ifcSnow'); call def_variab(trim(infile),(/hru_DimName,ifcSnowAndTime_DimName/),mvar_meta(ivar),nf90_double,err,cmessage)
   case('ifcSoil'); call def_variab(trim(infile),(/hru_DimName,ifcSoilAndTime_DimName/),mvar_meta(ivar),nf90_double,err,cmessage)
   case('ifcToto'); call def_variab(trim(infile),(/hru_DimName,ifcTotoAndTime_DimName/),mvar_meta(ivar),nf90_double,err,cmessage)
   case default; err=35; message=trim(message)//"varTypeNotFound"; return
  endselect
  ! check variable definition was OK
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do ! loop through model variables
 ! **********************************************************************************************************
 ! ***** define local column model indices -- dimensions depend on the variable type
 ! **********************************************************************************************************
 do ivar=1,size(indx_meta)
  if (.not.indx_meta(ivar)%v_write) cycle
  select case(trim(indx_meta(ivar)%vartype))
   case('scalarv'); call def_variab(trim(infile),(/hru_DimName,Timestep_DimName/),indx_meta(ivar),nf90_int,err,cmessage)
   case('midToto'); call def_variab(trim(infile),(/hru_DimName,midTotoAndTime_DimName/),indx_meta(ivar),nf90_int,err,cmessage)
   case default; err=35; message=trim(message)//"varTypeNotFound"; return
  endselect
  ! check variable definition was OK
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do ! loop through model variable
 ! **********************************************************************************************************
 ! ***** define local column model variables -- dimensions depend on the variable type
 ! **********************************************************************************************************
 do ivar=1,size(bvar_meta)
  if (.not.bvar_meta(ivar)%v_write) cycle
  select case(trim(bvar_meta(ivar)%vartype))
   case('scalarv'); call def_variab(trim(infile),(/Timestep_DimName/),bvar_meta(ivar),nf90_double,err,cmessage)
   case('routing'); call def_variab(trim(infile),(/routing_DimName/), bvar_meta(ivar),nf90_double,err,cmessage)
   case default; err=35; message=trim(message)//"varTypeNotFound"; return
  endselect
  ! check variable definition was OK
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do ! loop through model variables

 end subroutine def_output


 ! **********************************************************************************************************
 ! private subroutine ini_create: initial create
 ! **********************************************************************************************************
 subroutine ini_create(nHRU,infile,err,message)
 ! variables to define number of steps per file (total number of time steps, step length, etc.)
 USE multiconst,only:secprday           ! number of seconds per day
 USE data_struc,only:data_step          ! time step of model forcing data (s)
 USE data_struc,only:numtim             ! number of time steps
 ! model model index structures
 USE data_struc,only:indx_data          ! data structures
 USE data_struc,only:ix_soil            ! named variable to identify a soil layer
 USE var_lookup,only:iLookINDEX         ! named variables for structure elements
 ! model decisions
 USE data_struc,only:model_decisions    ! model decision structure
 USE var_lookup,only:iLookDECISIONS     ! named variables for elements of the decision structure
 USE mDecisions_module,only:&
  sameRulesAllLayers, & ! SNTHERM option: same combination/sub-dividion rules applied to all layers
  rulesDependLayerIndex ! CLM option: combination/sub-dividion rules depend on layer index
 implicit none
 ! declare dummy variables
 integer(i4b),intent(in)     :: nHRU                       ! number of HRUs
 character(*),intent(in)     :: infile                     ! filename
 integer(i4b),intent(out)    :: err                        ! error code
 character(*),intent(out)    :: message                    ! error message
 ! define local variables
 integer(i4b)                :: ncid                       ! NetCDF file ID
 integer(i4b)                :: dimID
 integer(i4b)                :: maxRouting=1000            ! maximum length of routing vector
 integer(i4b),parameter      :: maxSpectral=2              ! maximum number of spectral bands
 integer(i4b),parameter      :: scalarLength=1             ! length of scalar variable
 integer(i4b)                :: meanSnowLayersPerStep      ! mean number of snow layers per time step
 integer(i4b)                :: maxStepsPerFile            ! maximum number of time steps to be stored in each file
 integer(i4b)                :: maxLength                  ! maximum length of the variable vector
 integer(i4b)                :: nSoil                      ! number of soil layers
 ! initialize error control
 err=0;message="f-iniCreate/"
 ! define number of soil layers
 nSoil = count(indx_data%var(iLookINDEX%layerType)%dat == ix_soil)  ! number of soil layers
 ! identify length of the variable vector
 maxStepsPerFile = min(numtim, nint(366._dp * secprday/data_step) )
 select case(model_decisions(iLookDECISIONS%snowLayers)%iDecision)
  case(sameRulesAllLayers);    meanSnowLayersPerStep = 1000
  case(rulesDependLayerIndex); meanSnowLayersPerStep = 5
  case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
 end select ! (option to combine/sub-divide snow layers)
 maxLength = maxStepsPerFile*(nSoil+1 + meanSnowLayersPerStep)
 print*, 'maxStepsPerFile, maxLength = ', maxStepsPerFile, maxLength
 ! create output file
 err = nf90_create(trim(infile),nf90_classic_model,ncid)
 message='iCreate[create]'; call netcdf_err(err,message); if (err/=0) return
 ! create time dimension (unlimited)
 err = nf90_def_dim(ncid, trim(timestep_DimName), nf90_unlimited, dimId)
 message='iCreate[time]'; call netcdf_err(err,message); if (err/=0) return
 ! create scalar dimension
 err = nf90_def_dim(ncid, trim(scalar_DimName), scalarLength, dimId)
 message='iCreate[scalar]'; call netcdf_err(err,message); if (err/=0) return
 ! create HRU dimension
 err = nf90_def_dim(ncid, trim(hru_DimName), nHRU, dimId)
 message='iCreate[HRU]'; call netcdf_err(err,message); if (err/=0) return
 ! create spectral band dimension
 err = nf90_def_dim(ncid, trim(wLength_DimName), maxSpectral, dimId)
 message='iCreate[spectral]'; call netcdf_err(err,message); if (err/=0) return
 ! create dimension for the time-delay routing variables
 err = nf90_def_dim(ncid, trim(routing_DimName), maxRouting, dimId)
 message='iCreate[routing]'; call netcdf_err(err,message); if (err/=0) return
 ! create dimension for midSnow+time
 err = nf90_def_dim(ncid, trim(midSnowAndTime_DimName), maxLength, dimId)
 message='iCreate[midSnow]'; call netcdf_err(err,message); if (err/=0) return
 ! create dimension for midSoil+time
 err = nf90_def_dim(ncid, trim(midSoilAndTime_DimName), maxLength, dimId)
 message='iCreate[midSoil]'; call netcdf_err(err,message); if (err/=0) return
 ! create dimension for midToto+time
 err = nf90_def_dim(ncid, trim(midTotoAndTime_DimName), maxLength, dimId)
 message='iCreate[minToto]'; call netcdf_err(err,message); if (err/=0) return
 ! create dimension for ifcSnow+time
 err = nf90_def_dim(ncid, trim(ifcSnowAndTime_DimName), maxLength, dimId)
 message='iCreate[ifcSnow]'; call netcdf_err(err,message); if (err/=0) return
 ! create dimension for ifcSoil+time
 err = nf90_def_dim(ncid, trim(ifcSoilAndTime_DimName), maxLength, dimId)
 message='iCreate[ifcSoil]'; call netcdf_err(err,message); if (err/=0) return
 ! create dimension for ifcToto+time
 err = nf90_def_dim(ncid, trim(ifcTotoAndTime_DimName), maxLength, dimId)
 message='iCreate[ifcToto]'; call netcdf_err(err,message); if (err/=0) return
 ! close NetCDF file
 err = nf90_enddef(ncid); call netcdf_err(err,message); if (err/=0) return
 err = nf90_close(ncid); call netcdf_err(err,message); if (err/=0) return
 end subroutine ini_create


 ! **********************************************************************************************************
 ! private subroutine put_attrib: put global attributes as character string
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
 ! private subroutine def_variab: define variables
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
 endif
 end subroutine netcdf_err


end module def_output_module

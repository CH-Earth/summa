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

module read_icond_module
USE nrtype
USE netcdf
! define modeling decisions
USE mDecisions_module,only:  &
 moisture,                   & ! moisture-based form of Richards' equation
 mixdform                      ! mixed form of Richards' equation
implicit none
private
public::read_icond
public::read_icond_layers
contains

 ! ************************************************************************************************
 ! public subroutine read_icond_layers: read model initial conditions file for number of snow/soil layers
 ! ************************************************************************************************
 subroutine read_icond_layers(infile,nGRU,nHRU,indx_meta,err,message)
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nrtype
 USE var_lookup,only:iLookIndex             ! variable lookup structure
 USE globalData,only:gru_struc              ! gru-hru mapping structures
 USE netcdf_util_module,only:nc_file_close  ! close netcdf file
 USE netcdf_util_module,only:nc_file_open   ! close netcdf file
 USE netcdf_util_module,only:netcdf_err     ! netcdf error handling
 USE data_types,only:gru_hru_intVec         ! actual data
 USE data_types,only:var_info               ! metadata 
 implicit none

 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 character(*)        ,intent(in)     :: infile       ! name of input (restart) file 
 integer(i4b)        ,intent(in)     :: nGRU, nHRU   ! total # of GRUs and HRUs 
 type(var_info)      ,intent(in)     :: indx_meta(:) ! metadata 
 integer(i4b)        ,intent(out)    :: err          ! error code
 character(*)        ,intent(out)    :: message      ! returned error message

 ! locals
 integer(i4b)            :: ncID                     ! netcdf file id
 integer(i4b)            :: snowID, soilID           ! netcdf variable ids
 integer(i4b)            :: iGRU, iHRU, cHRU         ! loop indexes
 integer(i4b)            :: snowData(nHRU)           ! number of snow layers in all HRUs
 integer(i4b)            :: soilData(nHRU)           ! number of soil layers in all HRUs
 character(len=256)      :: cmessage                 ! downstream error message 

 ! --------------------------------------------------------------------------------------------------------
 ! initialize error message
 err=0
 message = 'read_icond_layers/'

 ! open netcdf file
 call nc_file_open(infile,nf90_nowrite,ncid,err,cmessage);
 if (err/=0) then; message=trim(message)//trim(cmessage); return; endif

 ! get variable ids
 err = nf90_inq_varid(ncid,trim(indx_meta(iLookIndex%nSnow)%varName),snowid); call netcdf_err(err,message)
 err = nf90_inq_varid(ncid,trim(indx_meta(iLookIndex%nSoil)%varName),soilid); call netcdf_err(err,message)

 ! get data
 err = nf90_get_var(ncid,snowid,snowData); call netcdf_err(err,message) 
 err = nf90_get_var(ncid,soilid,soilData); call netcdf_err(err,message) 

 ! assign to index structure - gru by hru
 cHRU = 0
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   cHRU = cHRU + 1
   gru_struc(iGRU)%hruInfo(iHRU)%nSnow = snowData(cHRU)
   gru_struc(iGRU)%hruInfo(iHRU)%nSoil = soilData(cHRU)
  enddo
 enddo

 ! close file
 call nc_file_close(ncid,err,cmessage)
 if(err/=0)then;message=trim(message)//trim(cmessage);return;endif

 end subroutine read_icond_layers


 ! ************************************************************************************************
 ! public subroutine read_icond: read model initial conditions
 ! ************************************************************************************************
 subroutine read_icond(filename,                      & ! name of initial conditions file
                       nGRU,nHRU,                     & ! number of GRUs and HRUs
                       prog_meta,                     & ! metadata
                       progData,                      & ! model prognostic (state) variables
                       mparData,                      & ! model parameters
                       indxData,                      & ! layer index data
                       err,message)                     ! error control
 ! --------------------------------------------------------------------------------------------------------
 ! modules
 USE nrtype
 USE var_lookup,only:iLookVarType           ! variable lookup structure
 USE var_lookup,only:iLookParam             ! variable lookup structure
 USE var_lookup,only:iLookProg              ! variable lookup structure
 USE var_lookup,only:iLookIndex             ! variable lookup structure
 USE globalData,only:gru_struc              ! gru-hru mapping structures
 USE netcdf_util_module,only:nc_file_close  ! close netcdf file
 USE netcdf_util_module,only:nc_file_open   ! close netcdf file
 USE netcdf_util_module,only:netcdf_err     ! netcdf error handling
 USE data_types,only:gru_hru_doubleVec      ! actual data
 USE data_types,only:gru_hru_intVec         ! actual data
 USE data_types,only:gru_hru_double         ! actual data
 USE data_types,only:var_info               ! metadata 
 USE globaldata,only:ix_soil,ix_snow        ! named variables to describe the type of layer
 USE multiconst,only:&
                       LH_fus,    &  ! latent heat of fusion                (J kg-1)
                       iden_ice,  &  ! intrinsic density of ice             (kg m-3)
                       iden_water,&  ! intrinsic density of liquid water    (kg m-3)
                       gravity,   &  ! gravitational acceleration           (m s-2)
                       Tfreeze       ! freezing point of pure water         (K)
 USE snow_utils_module,only:fracliquid             ! compute volumetric fraction of liquid water in snow based on temperature
 USE updatState_module,only:updateSnow             ! update snow states
 USE updatState_module,only:updateSoil             ! update soil states
 implicit none

 ! --------------------------------------------------------------------------------------------------------
 ! variable declarations
 ! dummies
 character(*)           ,intent(in)     :: filename     ! name of netcdf file containing the initial conditions
 integer(i4b)           ,intent(in)     :: nGRU, nHRU   ! number of response units
 type(var_info)         ,intent(in)     :: prog_meta(:) ! prognostic metadata
 type(gru_hru_doubleVec),intent(inout)  :: progData     ! prognostic vars 
 type(gru_hru_double)   ,intent(in)     :: mparData     ! parameters 
 type(gru_hru_intVec)   ,intent(inout)  :: indxData     ! layer indexes 
 integer(i4b)           ,intent(out)    :: err          ! error code
 character(*)           ,intent(out)    :: message      ! returned error message

 ! locals
 character(len=256)                     :: cmessage     ! downstream error message
 logical(lgt)                           :: snowExists   ! query whether to read snow layers
 integer(i4b)                           :: iVar         ! loop index 
 integer(i4b)                           :: iGRU         ! loop index 
 integer(i4b)                           :: iHRU         ! loop index 
 integer(i4b)                           :: cHRU         ! loop index 
 integer(i4b)                           :: dimID        ! varible dimension ids
 integer(i4b)                           :: ncVarID      ! variable ID in netcdf file
 character(256)                         :: dimName      ! not used except as a placeholder in call to inq_dim function
 integer(i4b)                           :: dimLen       ! data dimensions
 integer(i4b)                           :: ncID         ! netcdf file ID
 real(dp)    ,allocatable               :: varData(:,:) ! variable data storage        

 character(len=32),parameter            :: hruDimName    ='hru'      ! dimension name for HRUs
 character(len=32),parameter            :: specDimName   ='spectral' ! dimension name for spectral bands
 character(len=32),parameter            :: midSnowDimName='midSnow'  ! dimension name for snow-only layers
 character(len=32),parameter            :: midSoilDimName='midSoil'  ! dimension name for soil-only layers
 character(len=32),parameter            :: midTotoDimName='midToto'  ! dimension name for layered varaiables
 character(len=32),parameter            :: ifcSnowDimName='ifcSnow'  ! dimension name for snow-only layers
 character(len=32),parameter            :: ifcSoilDimName='ifcSoil'  ! dimension name for soil-only layers
 character(len=32),parameter            :: ifcTotoDimName='ifcToto'  ! dimension name for layered varaiables

 ! temporary variables for realism checks
 integer(i4b)                      :: iLayer              ! layer index
 real(dp)                          :: fLiq                ! fraction of liquid water on the vegetation canopy (-)
 real(dp)                          :: vGn_m               ! van Genutchen "m" parameter (-)
 real(dp)                          :: tWat                ! total water on the vegetation canopy (kg m-2)
 real(dp)                          :: scalarTheta         ! liquid water equivalent of total water [liquid water + ice] (-)
 real(dp)                          :: h1,h2               ! used to check depth and height are consistent
 integer(i4b)                      :: nLayers         ! total number of layers
 real(dp)                          :: kappa               ! constant in the freezing curve function (m K-1)
 real(dp)                          :: maxVolFracLiq       ! maximum volumetric fraction of liquid water (used in moisture-based form of Richards' equation)
 integer(i4b)                      :: nSnow           ! number of snow layers
 integer(i4b)                      :: nSoil           ! number of soil layers

 ! --------------------------------------------------------------------------------------------------------

 ! Start procedure here
 err=0; message="read_icond/"

 ! --------------------------------------------------------------------------------------------------------
 ! (1) read the file
 ! --------------------------------------------------------------------------------------------------------

 ! open netcdf file
 call nc_file_open(filename,nf90_nowrite,ncID,err,cmessage)
 if (err/=0) then; message=trim(message)//trim(cmessage); return; endif

 ! is there any snow in initial conditions file?
 snowExists = .true.
 err = nf90_inq_dimid(ncID,trim(midSnowDimName),dimID)
 if (err/=0) snowExists = .false.

 ! loop through prognostic variables
 do iVar = 1,size(prog_meta)

  ! if there is no snow, do not look for snow variables
  if ((.not.snowExists).and.((prog_meta(iVar)%varType==iLookVarType%midSnow).or.(prog_meta(iVar)%varType==iLookVarType%ifcSnow))) then
   do iGRU = 1,nGRU
    do iHRU = 1,gru_struc(iGRU)%hruCount
     progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat = 0.
    enddo
   enddo
   cycle
  endif

  ! get variable id
  err = nf90_inq_varid(ncID,trim(prog_meta(iVar)%varName),ncVarID); call netcdf_err(err,message)

  ! get variable dimension IDs
  selectcase (prog_meta(iVar)%varType)
   case (iLookVarType%scalarv); err = nf90_inq_dimid(ncID,trim(hruDimName)    ,dimID); call netcdf_err(err,message)
   case (iLookVarType%wlength); err = nf90_inq_dimid(ncID,trim(specDimName)   ,dimID); call netcdf_err(err,message)
   case (iLookVarType%midSoil); err = nf90_inq_dimid(ncID,trim(midSoilDimName),dimID); call netcdf_err(err,message)
   case (iLookVarType%midToto); err = nf90_inq_dimid(ncID,trim(midTotoDimName),dimID); call netcdf_err(err,message)
   case (iLookVarType%ifcSoil); err = nf90_inq_dimid(ncID,trim(ifcSoilDimName),dimID); call netcdf_err(err,message)
   case (iLookVarType%ifcToto); err = nf90_inq_dimid(ncID,trim(ifcTotoDimName),dimID); call netcdf_err(err,message)
   case (iLookVarType%midSnow); err = nf90_inq_dimid(ncID,trim(midSnowDimName),dimID); call netcdf_err(err,message)
   case (iLookVarType%ifcSnow); err = nf90_inq_dimid(ncID,trim(ifcSnowDimName),dimID); call netcdf_err(err,message)
  endselect

  ! get variable dimensions
  err = nf90_inquire_dimension(ncID,dimID,dimName,dimLen); call netcdf_err(err,message)

  ! iniitialize the varialbe data
  allocate(varData(nHRU,dimLen))

  ! get data
  err = nf90_get_var(ncID,ncVarID,varData); call netcdf_err(err,message) 
 
  ! store data in prognostics structure 
  ! loop through GRUs
  cHRU = 0
  do iGRU = 1,nGRU
   do iHRU = 1,gru_struc(iGRU)%hruCount
    cHRU = cHRU + 1
    progData%gru(iGRU)%hru(iHRU)%var(iVar)%dat = varData(cHRU,:) 
   enddo ! iHRU
  enddo ! iGRU

  ! deallocate storage vector for next variable
  deallocate(varData)

 enddo ! iVar 

 ! --------------------------------------------------------------------------------------------------------
 ! (2) set number of layers 
 ! --------------------------------------------------------------------------------------------------------
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount

   ! save the number of layers
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSnow)%dat(1)   = gru_struc(iGRU)%hruInfo(iHRU)%nSnow 
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSoil)%dat(1)   = gru_struc(iGRU)%hruInfo(iHRU)%nSoil
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nLayers)%dat(1) = gru_struc(iGRU)%hruInfo(iHRU)%nSnow + gru_struc(iGRU)%hruInfo(iHRU)%nSoil

   ! initalize the indices for midSnow, midSoil, midToto, and ifcToto
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%midSnowStartIndex)%dat(1) = 1
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%midSoilStartIndex)%dat(1) = 1
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%midTotoStartIndex)%dat(1) = 1
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%ifcSnowStartIndex)%dat(1) = 1
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%ifcSoilStartIndex)%dat(1) = 1
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%ifcTotoStartIndex)%dat(1) = 1

   ! set layer type
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%layerType)%dat(1:gru_struc(iGRU)%hruInfo(iHRU)%nSnow) = ix_snow
   indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%layerType)%dat(1:gru_struc(iGRU)%hruInfo(iHRU)%nSoil) = ix_soil

  enddo
 enddo

 ! --------------------------------------------------------------------------------------------------------
 ! (3) check that the initial conditions do not conflict with parameters, structure, etc.
 ! --------------------------------------------------------------------------------------------------------
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
   ! ensure the spectral average albedo is realistic
   if(progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSnowAlbedo)%dat(1) > mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMax)) &
      progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSnowAlbedo)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMax)
   if(progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSnowAlbedo)%dat(1) < mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMinWinter)) &
      progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSnowAlbedo)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMinWinter)
   ! ensure the visible albedo is realistic
   if(progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1) > mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMaxVisible)) &
      progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMaxVisible)
   if(progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1) < mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMinVisible)) &
      progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1) = mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMinVisible)
   ! ensure the nearIR albedo is realistic
   if(progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(2) > mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMaxNearIR)) &
      progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(2) = mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMaxNearIR)
   if(progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(2) < mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMinNearIR)) &
      progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(2) = mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%albedoMinNearIR)
  enddo
 enddo

 
 ! ensure the initial conditions are consistent with the constitutive functions
 do iGRU = 1,nGRU
  do iHRU = 1,gru_struc(iGRU)%hruCount
 
   ! associate local variables with variables in the data structures
   associate(&
   ! state variables in the vegetation canopy
   scalarCanopyTemp  => progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanopyTemp)%dat(1)   , & ! canopy temperature
   scalarCanopyIce   => progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanopyIce)%dat(1)    , & ! mass of ice on the vegetation canopy (kg m-2)
   scalarCanopyLiq   => progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanopyLiq)%dat(1)    , & ! mass of liquid water on the vegetation canopy (kg m-2)
   ! state variables in the snow+soil domain
   mLayerTemp        => progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerTemp)%dat            , & ! temperature (K)
   mLayerVolFracLiq  => progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracLiq)%dat      , & ! volumetric fraction of liquid water in each snow layer (-)
   mLayerVolFracIce  => progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracIce)%dat      , & ! volumetric fraction of ice in each snow layer (-)
   mLayerMatricHead  => progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerMatricHead)%dat      , & ! matric head (m)
   mLayerLayerType   => indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%layerType)%dat            , & ! type of layer (ix_soil or ix_snow)
   ! model parameters
   vGn_alpha         => mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%vGn_alpha)                , & ! van Genutchen "alpha" parameter (m-1)
   vGn_n             => mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%vGn_n)                    , & ! van Genutchen "n" parameter (-)
   theta_sat         => mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%theta_sat)                , & ! soil porosity (-)
   theta_res         => mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%theta_res)                , & ! soil residual volumetric water content (-)
   snowfrz_scale     => mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%snowfrz_scale)            , & ! scaling parameter for the snow freezing curve (K-1)
   FCapil            => mparData%gru(iGRU)%hru(iHRU)%var(iLookPARAM%FCapil)                     & ! fraction of pore space in tension storage (-)
   )  ! (associate local variables with model parameters)

   ! compute the maximum volumetric fraction of liquid water -- used to avoid problems of super-saturation in the moisture-based form of Richards' equation
   maxVolFracLiq = theta_sat - 1.e-4_dp

   ! compute the van Genutchen "m" parameter (-)
   vGn_m = 1._dp - 1._dp/vGn_n

   ! compute the constant in the freezing curve function (m K-1)
   kappa  = (iden_ice/iden_water)*(LH_fus/(gravity*Tfreeze))  ! NOTE: J = kg m2 s-2

   ! modify the liquid water and ice in the canopy
   if(scalarCanopyIce > 0._dp .and. scalarCanopyTemp > Tfreeze)then
    message=trim(message)//'canopy ice > 0 when canopy temperature > Tfreeze'
    err=20; return
   endif
   fLiq = fracliquid(scalarCanopyTemp,snowfrz_scale)  ! fraction of liquid water (-)
   tWat = scalarCanopyLiq + scalarCanopyIce           ! total water (kg m-2)
   scalarCanopyLiq = fLiq*tWat                        ! mass of liquid water on the canopy (kg m-2)
   scalarCanopyIce = (1._dp - fLiq)*tWat              ! mass of ice on the canopy (kg m-2)

   ! number of layers
   nLayers = gru_struc(iGRU)%hruInfo(iHRU)%nSnow + gru_struc(iGRU)%hruInfo(iHRU)%nSoil
   nSnow   = gru_struc(iGRU)%hruInfo(iHRU)%nSnow

   ! loop through all layers
   do iLayer=1,nLayers

    ! compute liquid water equivalent of total water (liquid plus ice)
    if (iLayer>nSnow) then ! soil layer = no volume expansion
     scalarTheta = mLayerVolFracIce(iLayer) + mLayerVolFracLiq(iLayer)
    else ! snow layer = volume expansion allowed
     scalarTheta = mLayerVolFracIce(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiq(iLayer)
    endif

    ! check that the initial volumetric fraction of liquid water and ice is reasonable
    select case(mlayerLayerType(iLayer))

     ! ***** snow
     case(ix_snow)
      ! (check liquid water)
      if(mLayerVolFracLiq(iLayer) < 0._dp .or. mLayerVolFracLiq(iLayer) > 1._dp)then
       write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of liquid water < 0 or > 1: layer = ',iLayer
       err=20; return
      endif
      ! (check ice)
      if(mLayerVolFracIce(iLayer) < 0.05_dp .or. mLayerVolFracIce(iLayer) > 0.80_dp)then
       write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of ice < 0.05 or > 0.70: layer = ',iLayer
       err=20; return
      endif
      ! check total water
      if(scalarTheta < 0.05_dp .or. scalarTheta > 0.80_dp)then
       write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with theta (total water fraction [liquid + ice]) < 0.05 or > 0.70: layer = ',iLayer
       err=20; return
      endif

     ! ***** soil
     case(ix_soil)
      ! (check liquid water)
      if(mLayerVolFracLiq(iLayer) < theta_res .or. mLayerVolFracLiq(iLayer) > theta_sat)then
       write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of liquid water < theta_res or > theta_sat: layer = ',iLayer
       err=20; return
      endif
      ! (check ice)
      if(mLayerVolFracIce(iLayer) < 0._dp .or. mLayerVolFracIce(iLayer) > theta_sat)then
       write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with volumetric fraction of ice < 0 or > theta_sat: layer = ',iLayer
       err=20; return
      endif
      ! check total water
      if(scalarTheta < theta_res .or. scalarTheta > theta_sat)then
       write(message,'(a,1x,i0)') trim(message)//'cannot initialize the model with theta (total water fraction [liquid + ice]) < theta_res or > theta_sat: layer = ',iLayer
       err=20; return
      endif

     case default
print*,'debug 1: ',ix_snow,ix_soil
print*,'debug 2: ',ilayer,nlayers,mlayerLayerType(iLayer)
print*,indxData%gru(iGRU)%hru(iHRU)%var(iLookINDEX%layerType)%dat  
      err=20; message=trim(message)//'cannot identify layer type'; return
    end select

    ! process snow and soil separately
    select case(mLayerLayerType(iLayer))
  
     ! ** snow
     case(ix_snow)
  
      ! check that snow temperature is less than freezing
      if(mLayerTemp(iLayer) > Tfreeze)then
       message=trim(message)//'initial snow temperature is greater than freezing'
       err=20; return
      endif
  
      ! ensure consistency among state variables
      call updateSnow(&
                      ! input
                      mLayerTemp(iLayer),                        & ! intent(in): temperature (K)
                      scalarTheta,                               & ! intent(in): mass fraction of total water (-)
                      snowfrz_scale,                             & ! intent(in): scaling parameter for the snow freezing curve (K-1)
                      ! output
                      mLayerVolFracLiq(iLayer),                  & ! intent(out): volumetric fraction of liquid water (-)
                      mLayerVolFracIce(iLayer),                  & ! intent(out): volumetric fraction of ice (-)
                      fLiq,                                      & ! intent(out): fraction of liquid water (-)
                      err,cmessage)                                ! intent(out): error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
  
     ! ** soil
     case(ix_soil)
  
      ! ensure consistency among state variables
      call updateSoil(&
                      ! input
                      mLayerTemp(iLayer),                        & ! intent(in): layer temperature (K)
                      mLayerMatricHead(iLayer-nSnow),            & ! intent(in): matric head (m)
                      vGn_alpha,vGn_n,theta_sat,theta_res,vGn_m, & ! intent(in): van Genutchen soil parameters
                      ! output
                      scalarTheta,                               & ! intent(out): volumetric fraction of total water (-)
                      mLayerVolFracLiq(iLayer),                  & ! intent(out): volumetric fraction of liquid water (-)
                      mLayerVolFracIce(iLayer),                  & ! intent(out): volumetric fraction of ice (-)
                      err,cmessage)                                ! intent(out): error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
  
     case default; err=10; message=trim(message)//'unknown case for model layer'; return
    endselect
  
   end do  ! (looping through layers)
  
   ! end association to variables in the data structures
   end associate
  
   ! if snow layers exist, compute snow depth and SWE
   if(nSnow > 0)then
    progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSnowDepth)%dat(1) =      -progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%iLayerHeight)%dat(0)
    progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSWE)%dat(1)       = sum( (progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                                               progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracIce)%dat(1:nSnow)*iden_ice)        * &
                                                                               progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerDepth)%dat(1:nSnow) )
   endif  ! if snow layers exist
 
   ! check that the layering is consistent
   do iLayer=1,nLayers
    h1 = sum(progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerDepth)%dat(1:iLayer)) ! sum of the depths up to the current layer
    h2 = progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%iLayerHeight)%dat(iLayer) - progData%gru(iGRU)%hru(iHRU)%var(iLookPROG%iLayerHeight)%dat(0)  ! difference between snow-atm interface and bottom of layer
    if(abs(h1 - h2) > 1.e-12_dp)then
     write(message,'(a,1x,i0)') trim(message)//'mis-match between layer depth and layer height [suggest round numbers in initial conditions file]; layer = ', iLayer
     err=20; return
    endif
   end do

  enddo ! iHRU
 enddo ! iGRU

 end subroutine read_icond

end module read_icond_module

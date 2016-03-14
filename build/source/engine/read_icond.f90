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
! define modeling decisions
USE mDecisions_module,only:  &
 moisture,                   & ! moisture-based form of Richards' equation
 mixdform                      ! mixed form of Richards' equation
implicit none
private
public::read_icond
contains


 ! ************************************************************************************************
 ! public subroutine read_icond: read model initial conditions
 ! ************************************************************************************************
 subroutine read_icond(nSnow,      &   ! number of snow layers
                       nSoil,      &   ! number of soil layers
                       parData  ,  &   ! vector of model parameters
                       indx_data,  &   ! data structure of model indices
                       prog_data,  &   ! model prognostic (state) variables
                       err,message)  ! error control
 USE multiconst, only:&
                       LH_fus,    &  ! latent heat of fusion                (J kg-1)
                       iden_ice,  &  ! intrinsic density of ice             (kg m-3)
                       iden_water,&  ! intrinsic density of liquid water    (kg m-3)
                       gravity,   &  ! gravitational acceleration           (m s-2)
                       Tfreeze       ! freezing point of pure water         (K)
 ! modules
 USE summaFileManager,only:SETNGS_PATH             ! path for metadata files
 USE summaFileManager,only:MODEL_INITCOND          ! model initial conditions file
 USE snow_utils_module,only:fracliquid             ! compute volumetric fraction of liquid water in snow based on temperature
 USE snow_utils_module,only:templiquid             ! compute temperature of snow based on volumetric fraction of liquid water
 USE soil_utils_module,only:volFracLiq             ! compute volumetric fraction of liquid water based on matric head
 USE soil_utils_module,only:matricHead             ! compute matric head based on volumetric fraction of liquid water
 USE soil_utils_module,only:crit_soilT             ! compute temperature above which all water is unfrozen
 USE updatState_module,only:updateSnow             ! update snow states
 USE updatState_module,only:updateSoil             ! update soil states
 USE ascii_util_module,only:file_open              ! open file
 USE ascii_util_module,only:split_line             ! extract the list of variable names from the character string
 USE ascii_util_module,only:get_vlines             ! get a list of character strings from non-comment lines
 USE get_ixname_module,only:get_ixindex            ! access function to find index of elements in structure
 USE get_ixname_module,only:get_ixprog             ! access function to find index of elements in structure
 ! named variabes for model decisions
 USE var_lookup,only:iLookDECISIONS                ! named variables for elements of the decision structure
 USE var_lookup,only:iLookPARAM,iLookINDEX,iLookPROG ! named variables to describe structure elements
 ! metadata
 USE data_struc,only:prog_meta                     ! metadata for model prognostic (state) variables
 USE data_struc,only:ix_soil,ix_snow               ! named variables to describe the type of layer
 ! data types
 USE data_struc,only:var_dlength    ! x%var(:)%dat (dp)
 USE data_struc,only:var_ilength    ! x%var(:)%dat (i4b)
 implicit none
 ! define input
 integer(i4b),intent(in)           :: nSnow           ! number of snow layers
 integer(i4b),intent(in)           :: nSoil           ! number of soil layers
 real(dp),intent(in)               :: parData(:)      ! vector of model parameters
 type(var_ilength),intent(inout)   :: indx_data       ! data structure of model indices for a local HRU
 type(var_dlength),intent(inout)   :: prog_data       ! data structure of model prognostic (state) variables for a local HRU
 ! define output
 integer(i4b),intent(out)          :: err             ! error code
 character(*),intent(out)          :: message         ! error message
 ! define local variables
 integer(i4b)                      :: nLayers         ! total number of layers
 integer(i4b),parameter            :: missingInteger=-9999     ! missing value for integers
 real(dp),parameter                :: missingDouble=-9999._dp  ! missing value for double
 character(len=256)                :: cmessage        ! error message for downwind routine
 character(LEN=256)                :: infile          ! input filename
 integer(i4b),parameter            :: nBand=2         ! number of spectral bands
 integer(i4b),parameter            :: ix_miss=-999    ! index for missing data
 integer(i4b)                      :: unt            ! file unit (free unit output from file_open)
 integer(i4b)                      :: iline           ! loop through lines in the file
 integer(i4b)                      :: iword           ! loop through words in a line
 integer(i4b),parameter            :: maxLines=10000  ! maximum lines in the file
 character(LEN=256)                :: temp            ! single line of information
 integer(i4b)                      :: iend            ! check for the end of the file
 character(LEN=256)                :: namesScalarDesired(10) ! names of desired scalar variables
 logical(lgt),allocatable          :: checkGotVars(:) ! used to check if we have got desired variables
 character(LEN=256),allocatable    :: varnames(:)     ! vector of variable names
 character(LEN=256),allocatable    :: chardata(:)     ! vector of character data
 integer(i4b)                      :: ivar,jvar       ! index of model variable
 integer(i4b)                      :: layerType       ! ix_snow or ix_soil
 integer(i4b)                      :: nVars           ! number of model variables
 integer(i4b)                      :: iSnow           ! index of snow model layers
 integer(i4b)                      :: iSoil           ! index of soil model layers
 integer(i4b)                      :: iToto           ! index of model layers
 character(len=256),parameter      :: scalar_tag='scalar_icond' ! tag for the scalar initial conditions
 character(len=256),parameter      :: layer_tag='layer_icond'   ! tag for the layer initial conditions
 logical(lgt)                      :: scalar_flag=.false. ! flag determines if in the scalar portion of the file
 logical(lgt)                      :: layer_flag=.false.  ! flag determines if in the layer portion of the file
 logical(lgt)                      :: first_flag=.false.  ! flag determines if reading the variable names
 ! (ensure the initial conditions are consistent with the constitutive functions)
 integer(i4b)                      :: iLayer              ! layer index
 real(dp)                          :: scalarTheta         ! liquid water equivalent of total water [liquid water + ice] (-)
 real(dp)                          :: vGn_m               ! van Genutchen "m" parameter (-)
 real(dp)                          :: kappa               ! constant in the freezing curve function (m K-1)
 real(dp)                          :: maxVolFracLiq       ! maximum volumetric fraction of liquid water (used in moisture-based form of Richards' equation)
 real(dp)                          :: h1,h2               ! used to check depth and height are consistent
 real(dp)                          :: fLiq                ! fraction of liquid water on the vegetation canopy (-)
 real(dp)                          :: tWat                ! total water on the vegetation canopy (kg m-2)
 logical(lgt),parameter            :: doPrintStates=.false.  ! flag to print states
 ! Start procedure here
 err=0; message="read_icond/"

 ! check the missing data flag is OK
 if(ix_miss==ix_snow .or. ix_miss==ix_soil)then; err=20; message=trim(message)//&
  'missing value index is the same as ix_snow or ix_soil'; return; endif

 ! allocate space for the variable check vector
 allocate(checkGotVars(size(prog_meta)),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'allocating logical check vector'; return; endif
 checkGotVars(:) = .false.  ! initialize vector

 ! define desired scalar variables
 if(size(namesScalarDesired)/=10)then
  err=20; message=trim(message)//'expect 10 variables in namesScalarDesired'; return
 endif
 namesScalarDesired( 1) = 'dt_init'
 namesScalarDesired( 2) = 'scalarCanopyIce'
 namesScalarDesired( 3) = 'scalarCanopyLiq'
 namesScalarDesired( 4) = 'scalarCanairTemp'
 namesScalarDesired( 5) = 'scalarCanopyTemp'
 namesScalarDesired( 6) = 'scalarSnowAlbedo'
 namesScalarDesired( 7) = 'scalarSWE'
 namesScalarDesired( 8) = 'scalarSnowDepth'
 namesScalarDesired( 9) = 'scalarSfcMeltPond'
 namesScalarDesired(10) = 'scalarAquiferStorage'

 ! save the number of layers
 nLayers = nSnow+nSoil
 indx_data%var(iLookINDEX%nSnow)%dat(1)   = nSnow
 indx_data%var(iLookINDEX%nSoil)%dat(1)   = nSoil
 indx_data%var(iLookINDEX%nLayers)%dat(1) = nLayers

 ! initalize the indices for midSnow, midSoil, midToto, and ifcToto
 indx_data%var(iLookINDEX%midSnowStartIndex)%dat(1) = 1
 indx_data%var(iLookINDEX%midSoilStartIndex)%dat(1) = 1
 indx_data%var(iLookINDEX%midTotoStartIndex)%dat(1) = 1
 indx_data%var(iLookINDEX%ifcSnowStartIndex)%dat(1) = 1
 indx_data%var(iLookINDEX%ifcSoilStartIndex)%dat(1) = 1
 indx_data%var(iLookINDEX%ifcTotoStartIndex)%dat(1) = 1

 ! **********************************************************************************************
 ! (1) open files, etc.
 ! **********************************************************************************************
 ! build filename
 infile = trim(SETNGS_PATH)//trim(MODEL_INITCOND)
 ! open file
 call file_open(trim(infile),unt,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! **********************************************************************************************
 ! (4) read the scalar initial conditions
 ! **********************************************************************************************
 scalar_flag=.false. ! initialize scalar flag
 ! loop through file until reach the scalar_tag
 do iline=1,maxLines
  read(unt,'(a)',iostat=iend)temp; if(iend/=0)then; rewind(unt); exit; endif    ! read line of data, and exit if reach the end of file
  if (temp(1:1)=='!')cycle
  ! check if reached the end of the scalar definitions
  if (trim(temp)=='<end_'//trim(scalar_tag)//'>')then; rewind(unt); exit; endif
  ! check if in the scalar portion of the file
  if(scalar_flag)then
   ! split the line -- variable name followed by variable value
   call split_line(temp,chardata,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   ! find the variable -- variable name is first
   jvar = get_ixprog(trim(chardata(1)))
   if(jvar<=0)then; err=30; message=trim(message)//'variableNotFound[var='//trim(chardata(1))//']'; return; endif
   checkGotVars(jvar)=.true.
   ! read the data -- value is second
   read(chardata(2),*,iostat=err) prog_data%var(jvar)%dat(1)
   if(err/=0)then; err=40; message=trim(message)//"problemInternalRead[data='"//trim(chardata(2))//"']"; return; endif
   deallocate(chardata)
   !print*, jVar, trim(prog_meta(jvar)%vardesc), prog_data%var(jvar)%dat(1)
  endif    ! if we are in the scalar part of the file
  ! check if reached the start of the scalar definitions
  if (trim(temp)=='<start_'//trim(scalar_tag)//'>') scalar_flag=.true.
  ! check if reached the end of the file
  if (iline==maxLines) rewind(unt)
 end do  ! looping through lines
 ! check if we got the desired scalar variables
 do ivar=1,size(namesScalarDesired)
  jvar=get_ixprog(trim(namesScalarDesired(ivar)))
  if(.not.checkGotVars(jvar))then
   message=trim(message)//'initial condion undefined for variable '//trim(namesScalarDesired(ivar))
   err=20; return
  endif
 end do
 ! initialize the spectral albedo
 prog_data%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1:nBand) = prog_data%var(iLookPROG%scalarSnowAlbedo)%dat(1)

 ! **********************************************************************************************
 ! (5) read the layer initial conditions
 ! **********************************************************************************************
 iSnow=0            ! initialize the index of the snow vector
 iSoil=0            ! initialize the index of the soil vector
 iToto=0            ! initialize the index of the toto vector
 first_flag=.true.  ! flag to define first non-comment line, which defines the layer variables
 layer_flag=.false. ! initialize layer flag
 ! loop through file until reach the layer_tag
 do iline=1,maxLines

  read(unt,'(a)',iostat=iend)temp; if(iend/=0)exit    ! read line of data, and exit if reach the end of file
  if (temp(1:1)=='!')cycle
  ! check if reached the end of the layer definitions
  if(trim(temp)=='<end_'//trim(layer_tag)//'>')then; rewind(unt); exit; endif

  ! read layer data
  if(layer_flag) then

   ! ***** process the layer names and allocate space for the character data
   if(first_flag)then
    ! split the line into an array of words
    call split_line(temp,varnames,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    ! check that required data are present
    if(count(varnames=='layerType'       )==0)then; err=50; message=trim(message)//"layerType=missing"; return; endif
    if(count(varnames=='mLayerDepth'     )==0)then; err=50; message=trim(message)//"mLayerDepth=missing"; return; endif
    if(count(varnames=='iLayerHeight'    )==0)then; err=50; message=trim(message)//"iLayerHeight=missing"; return; endif
    if(count(varnames=='mLayerTemp'      )==0)then; err=50; message=trim(message)//"mLayerTemp=missing"; return; endif
    if(count(varnames=='mLayerVolFracIce')==0)then; err=50; message=trim(message)//"mLayerVolFracIce=missing"; return; endif
    if(count(varnames=='mLayerVolFracLiq')==0)then; err=50; message=trim(message)//"mLayerVolFracLiq=missing"; return; endif
    if(count(varnames=='mLayerMatricHead')==0)then; err=50; message=trim(message)//"mLayerMatricHead=missing"; return; endif
    ! allocate space for character data
    nVars = size(varnames)
    allocate(chardata(nVars),stat=err)
    if(err/=0)then;err=30;message=trim(message)//"problemAllocateChardata"; return; endif
    ! set flag to .false. -- now read data
    first_flag=.false.
    cycle
   endif ! (if reading the variable names

   ! ***** get the vector of data for a given layer
   read(temp,*,iostat=err) chardata
   if(err/=0)then;err=40;message=trim(message)//"problemInternalRead[data='"//trim(temp)//"']"; return; endif
   ! identify the layer type (snow or soil)
   layerType=ix_miss
   do iword=1,size(chardata)
    if(chardata(iword)=='snow') layerType = ix_snow
    if(chardata(iword)=='soil') layerType = ix_soil
    if(chardata(iword)=='snow' .or. chardata(iword)=='soil') exit ! exit once read the layer type
   end do
   if(layerType==ix_miss)then; err=40; message=trim(message)//"cannot identify the layer type"; return; endif
   ! increment the index of the snow or soil Layer
   if(layerType==ix_soil) iSoil = iSoil+1
   if(layerType==ix_snow) iSnow = iSnow+1
   ! increment the index of the concatanated vector
   iToto = iToto+1
   ! loop through initial conditions variables
   do ivar=1,nVars
    ! check if it is the layerType variable (special case)
    if(trim(varnames(ivar))=='layerType')then
     indx_data%var(iLookINDEX%layerType)%dat(iToto) = layerType
     cycle
    endif
    ! get the variable index
    jvar = get_ixprog(trim(varnames(ivar)))
    if(jvar<=0)then; err=40; message=trim(message)//"cannotFindVariableIndex[name='"//trim(varnames(ivar))//"']"; return; endif
    ! ***** populate the data variable *****
    select case(trim(prog_meta(jvar)%vartype))
     case('midSoil'); if(layerType==ix_soil) read(chardata(ivar),*,iostat=err) prog_data%var(jvar)%dat(iSoil)
     case('midSnow'); if(layerType==ix_snow) read(chardata(ivar),*,iostat=err) prog_data%var(jvar)%dat(iSnow)
     case('midToto');                        read(chardata(ivar),*,iostat=err) prog_data%var(jvar)%dat(iToto)
     case('ifcSnow');                        read(chardata(ivar),*,iostat=err) prog_data%var(jvar)%dat(iSnow-1)  ! IC = top interface
     case('ifcSoil');                        read(chardata(ivar),*,iostat=err) prog_data%var(jvar)%dat(iSoil-1)  ! IC = top interface
     case('ifcToto');                        read(chardata(ivar),*,iostat=err) prog_data%var(jvar)%dat(iToto-1)  ! IC = top interface
     case default
     err=40; message=trim(message)//"unknownInitCondType[name='"//trim(prog_meta(jvar)%varname)//"']"; return
    endselect
    if(err/=0)then;err=40;message=trim(message)//"problemInternalRead[data='"//trim(chardata(ivar))//"']"; return; endif
   end do    ! (looping through initial conditions variables)
  endif   ! (if layer flag)
  ! check if reached the start of the layer definitions
  if (trim(temp)=='<start_'//trim(layer_tag)//'>') layer_flag=.true.
 end do  ! looping through lines in the file
 ! close file
 close(unt)
 ! set iLayerHeight for the bottom layer
 prog_data%var(iLookPROG%iLayerHeight)%dat(nLayers) = &
 prog_data%var(iLookPROG%iLayerHeight)%dat(nLayers-1) + prog_data%var(iLookPROG%mLayerDepth)%dat(nLayers)
 ! check matric head is read correctly
 !print*,'mLayerMatricHead ', prog_data%var(iLookPROG%mLayerMatricHead)%dat(:)
 ! ***************************************************************************************
 ! ***************************************************************************************
 ! ensure the snow albedo is realistic
 ! ***************************************************************************************
 ! ***************************************************************************************
 ! ensure the spectral average albedo is realistic
 if(prog_data%var(iLookPROG%scalarSnowAlbedo)%dat(1) > parData(iLookPARAM%albedoMax)) &
    prog_data%var(iLookPROG%scalarSnowAlbedo)%dat(1) = parData(iLookPARAM%albedoMax)
 if(prog_data%var(iLookPROG%scalarSnowAlbedo)%dat(1) < parData(iLookPARAM%albedoMinWinter)) &
    prog_data%var(iLookPROG%scalarSnowAlbedo)%dat(1) = parData(iLookPARAM%albedoMinWinter)
 ! ensure the visible albedo is realistic
 if(prog_data%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1) > parData(iLookPARAM%albedoMaxVisible)) &
    prog_data%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1) = parData(iLookPARAM%albedoMaxVisible)
 if(prog_data%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1) < parData(iLookPARAM%albedoMinVisible)) &
    prog_data%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(1) = parData(iLookPARAM%albedoMinVisible)
 ! ensure the nearIR albedo is realistic
 if(prog_data%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(2) > parData(iLookPARAM%albedoMaxNearIR)) &
    prog_data%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(2) = parData(iLookPARAM%albedoMaxNearIR)
 if(prog_data%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(2) < parData(iLookPARAM%albedoMinNearIR)) &
    prog_data%var(iLookPROG%spectralSnowAlbedoDiffuse)%dat(2) = parData(iLookPARAM%albedoMinNearIR)
 ! ***************************************************************************************
 ! ***************************************************************************************
 ! ensure the initial conditions are consistent with the constitutive functions
 ! ***************************************************************************************
 ! ***************************************************************************************
 ! associate local variables with variables in the data structures
 associate(&
 ! state variables in the vegetation canopy
 scalarCanopyTemp  => prog_data%var(iLookPROG%scalarCanopyTemp)%dat(1) , & ! canopy temperature
 scalarCanopyIce   => prog_data%var(iLookPROG%scalarCanopyIce)%dat(1)  , & ! mass of ice on the vegetation canopy (kg m-2)
 scalarCanopyLiq   => prog_data%var(iLookPROG%scalarCanopyLiq)%dat(1)  , & ! mass of liquid water on the vegetation canopy (kg m-2)
 ! state variables in the snow+soil domain
 mLayerTemp        => prog_data%var(iLookPROG%mLayerTemp)%dat          , & ! temperature (K)
 mLayerVolFracLiq  => prog_data%var(iLookPROG%mLayerVolFracLiq)%dat    , & ! volumetric fraction of liquid water in each snow layer (-)
 mLayerVolFracIce  => prog_data%var(iLookPROG%mLayerVolFracIce)%dat    , & ! volumetric fraction of ice in each snow layer (-)
 mLayerMatricHead  => prog_data%var(iLookPROG%mLayerMatricHead)%dat    , & ! matric head (m)
 mLayerLayerType   => indx_data%var(iLookINDEX%layerType)%dat          , & ! type of layer (ix_soil or ix_snow)
 ! model parameters
 vGn_alpha         => parData(iLookPARAM%vGn_alpha)                    , & ! van Genutchen "alpha" parameter (m-1)
 vGn_n             => parData(iLookPARAM%vGn_n)                        , & ! van Genutchen "n" parameter (-)
 theta_sat         => parData(iLookPARAM%theta_sat)                    , & ! soil porosity (-)
 theta_res         => parData(iLookPARAM%theta_res)                    , & ! soil residual volumetric water content (-)
 snowfrz_scale     => parData(iLookPARAM%snowfrz_scale)                , & ! scaling parameter for the snow freezing curve (K-1)
 FCapil            => parData(iLookPARAM%FCapil)                         & ! fraction of pore space in tension storage (-)
 )  ! (associate local variables with model parameters)
 ! ***************************************************************************************

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

 ! loop through all layers
 do iLayer=1,nLayers

  ! compute liquid water equivalent of total water (liquid plus ice)
  scalarTheta = mLayerVolFracIce(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiq(iLayer)

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

   case default; err=20; message=trim(message)//'cannot identify layer type'; return

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
  prog_data%var(iLookPROG%scalarSnowDepth)%dat(1) =      -prog_data%var(iLookPROG%iLayerHeight)%dat(0)
  prog_data%var(iLookPROG%scalarSWE)%dat(1)       = sum( (prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(1:nSnow)*iden_water + &
                                                          prog_data%var(iLookPROG%mLayerVolFracIce)%dat(1:nSnow)*iden_ice) &
                                                        * prog_data%var(iLookPROG%mLayerDepth)%dat(1:nSnow) )
 endif  ! (if snow layers exist

 ! check that the layering is consistent
 do iLayer=1,nLayers
  h1 = sum(prog_data%var(iLookPROG%mLayerDepth)%dat(1:iLayer)) ! sum of the depths up to the current layer
  h2 = prog_data%var(iLookPROG%iLayerHeight)%dat(iLayer) - prog_data%var(iLookPROG%iLayerHeight)%dat(0)  ! difference between snow-atm interface and bottom of layer
  !write(*,'(a,1x,10(e20.10,1x))') 'h1, h2, (h1 - h2) = ', h1, h2, (h1 - h2)
  if(abs(h1 - h2) > 1.e-12_dp)then
   write(message,'(a,1x,i0)') trim(message)//'mis-match between layer depth and layer height [suggest round numbers in initial conditions file]; layer = ', iLayer
   err=20; return
  endif
 end do

 ! **********************************************************************************************
 ! deallocate variable names vector
 deallocate(varnames,chardata,stat=err)
 if(err/=0)then;err=30;message=trim(message)//'deallocating variable names vector'; return; endif
 ! deallocate variable check vector
 deallocate(checkGotVars,stat=err)
 if(err/=0)then; err=20; message=trim(message)//'deallocating logical check vector'; return; endif
 ! print states
 if(doPrintStates)then
  print*,'****************************************************************************************'
  print*, 'mLayerDepth      ', prog_data%var(iLookPROG%mLayerDepth)%dat(:)
  print*, 'iLayerHeight     ', prog_data%var(iLookPROG%iLayerHeight)%dat(:)
  print*, 'mLayerTemp       ', prog_data%var(iLookPROG%mLayerTemp)%dat(:)
  print*, 'mLayerVolFracIce ', prog_data%var(iLookPROG%mLayerVolFracIce)%dat(:)
  print*, 'mLayerVolFracLiq ', prog_data%var(iLookPROG%mLayerVolFracLiq)%dat(:)
  print*, 'mLayerMatricHead ', prog_data%var(iLookPROG%mLayerMatricHead)%dat(:)
  print*, 'scalarCanopyIce  ', prog_data%var(iLookPROG%scalarCanopyIce)%dat(:)
  print*, 'scalarCanopyLiq  ', prog_data%var(iLookPROG%scalarCanopyLiq)%dat(:)
  print*, 'scalarCanairTemp ', prog_data%var(iLookPROG%scalarCanairTemp)%dat(:)
  print*, 'scalarCanopyTemp ', prog_data%var(iLookPROG%scalarCanopyTemp)%dat(:)
  print*, 'scalarSnowAlbedo ', prog_data%var(iLookPROG%scalarSnowAlbedo)%dat(:)
  print*, 'scalarSnowDepth  ', prog_data%var(iLookPROG%scalarSnowDepth)%dat(:)
  print*, 'scalarSWE        ', prog_data%var(iLookPROG%scalarSWE)%dat(:)
  print*, 'layerType        ', indx_data%var(iLookINDEX%layerType)%dat(:)
  print*,'****************************************************************************************'
  !pause
 endif  ! (if printing states)
 end subroutine read_icond


end module read_icond_module

module read_icond_module
USE nrtype
USE mDecisions_module,only:  &
 moisture,                   & ! moisture-based form of Richards' equation
 mixdform                      ! mixed form of Richards' equation
implicit none
private
public::read_icond
contains

 ! ************************************************************************************************
 ! (1) new subroutine: read model initial conditions
 ! ************************************************************************************************
 subroutine read_icond(err,message)
 ! used to read model initial conditions
 USE multiconst, only:&
                       LH_fus,    &  ! latent heat of fusion                (J kg-1)
                       iden_ice,  &  ! intrinsic density of ice             (kg m-3)
                       iden_water,&  ! intrinsic density of liquid water    (kg m-3)
                       gravity,   &  ! gravitational acceleration           (m s-2)
                       Tfreeze       ! freezing point of pure water         (K)
 ! modules
 USE snow_utils_module,only:fracliquid             ! compute volumetric fraction of liquid water
 USE soil_utils_module,only:volFracLiq             ! compute volumetric fraction of liquid water based on matric head
 USE soil_utils_module,only:matricHead             ! compute matric head based on volumetric fraction of liquid water
 USE soil_utils_module,only:crit_soilT             ! compute temperature above which all water is unfrozen
 USE snow_fileManager,only:SETNGS_PATH             ! path for metadata files
 USE snow_fileManager,only:MODEL_INITCOND          ! model initial conditions file
 USE ascii_util_module,only:file_open              ! open file
 USE ascii_util_module,only:split_line             ! extract the list of variable names from the character string
 USE ascii_util_module,only:get_vlines             ! get a list of character strings from non-comment lines
 USE allocspace_module,only:alloc_mvar             ! allocate space for model variables
 USE allocspace_module,only:alloc_indx             ! allocate space for model variables
 ! data structures
 USE data_struc,only:model_decisions               ! model decision structure
 USE var_lookup,only:iLookDECISIONS                ! named variables for elements of the decision structure
 USE data_struc,only:mpar_data                     ! data for model parameetrs
 USE data_struc,only:mvar_data,mvar_meta           ! data/metadata for model variables
 USE data_struc,only:indx_data,indx_meta           ! data/metadata for model indices
 USE data_struc,only:ix_soil,ix_snow               ! named variables to describe the type of layer 
 USE var_lookup,only:iLookMVAR,iLookPARAM,iLookINDEX ! named variables to describe structure elements
 USE get_ixname_module,only:get_ixmvar,get_ixindex ! access function to find index of elements in structure
 implicit none
 ! define output
 integer(i4b),intent(out)       :: err             ! error code
 character(*),intent(out)       :: message         ! error message
 ! define local variables
 character(len=256)             :: cmessage        ! error message for downwind routine
 character(LEN=256)             :: infile          ! input filename
 integer(i4b),parameter         :: ix_miss=-999    ! index for missing data
 integer(i4b),parameter         :: unt=99          ! DK: need to either define units globally, or use getSpareUnit
 integer(i4b)                   :: iline           ! loop through lines in the file 
 integer(i4b)                   :: iword           ! loop through words in a line
 integer(i4b),parameter         :: maxLines=10000  ! maximum lines in the file 
 character(LEN=256)             :: temp            ! single line of information
 integer(i4b)                   :: iend            ! check for the end of the file
 character(LEN=256)             :: namesScalarDesired(7) ! names of desired scalar variables
 logical(lgt),allocatable       :: checkGotVars(:) ! used to check if we have got desired variables
 character(LEN=256),allocatable :: varnames(:)     ! vector of variable names
 character(LEN=256),allocatable :: chardata(:)     ! vector of character data
 integer(i4b)                   :: ivar,jvar       ! index of model variable
 integer(i4b)                   :: layerType       ! ix_snow or ix_soil
 integer(i4b)                   :: nVars           ! number of model variables
 integer(i4b)                   :: nSnow           ! number of snow layers
 integer(i4b)                   :: nSoil           ! number of soil layers
 integer(i4b)                   :: iSnow           ! index of snow model layers
 integer(i4b)                   :: iSoil           ! index of soil model layers
 integer(i4b)                   :: iToto           ! index of model layers
 integer(i4b)                   :: nLayers         ! number of model layers
 character(len=256),parameter   :: scalar_tag='scalar_icond' ! tag for the scalar initial conditions
 character(len=256),parameter   :: layer_tag='layer_icond'   ! tag for the layer initial conditions
 logical(lgt)                   :: scalar_flag=.false. ! flag determines if in the scalar portion of the file
 logical(lgt)                   :: layer_flag=.false.  ! flag determines if in the layer portion of the file
 logical(lgt)                   :: first_flag=.false.  ! flag determines if reading the variable names
 ! (ensure the initial conditions are consistent with the constitutive functions)
 integer(i4b)                   :: iLayer              ! layer index
 real(dp),pointer               :: scalarTemp          ! temperature (K)
 real(dp)                       :: scalarTheta         ! liquid water equivalent of total water [liquid water + ice] (-)
 integer(i4b),pointer           :: scalarLayerType     ! layer type
 real(dp),pointer               :: scalarVolFracIce    ! volumetric fraction of ice (-)
 real(dp),pointer               :: scalarVolFracLiq    ! volumetric fraction of liquid water (-)
 real(dp),pointer               :: scalarMatricHead    ! matric head (m)
 real(dp),pointer               :: snowfrz_scale       ! scaling parameter for the snow freezing curve (K-1)
 real(dp),pointer               :: vGn_alpha           ! van Genutchen "alpha" parameter
 real(dp),pointer               :: vGn_n               ! van Genutchen "n" parameter
 real(dp),pointer               :: theta_sat           ! soil porosity (-)
 real(dp),pointer               :: theta_res           ! soil residual volumetric water content (-)
 real(dp)                       :: Tcrit               ! temperature above which all water is unfrozen (K)
 real(dp)                       :: vGn_m               ! van Genutchen "m" parameter (-)
 real(dp)                       :: kappa               ! constant in the freezing curve function (m K-1) 
 real(dp)                       :: maxVolFracLiq       ! maximum volumetric fraction of liquid water (used in moisture-based form of Richards' equation)
 ! Start procedure here
 err=0; message="read_icond/"
 ! check the missing data flag is OK
 if(ix_miss==ix_snow .or. ix_miss==ix_soil)then; err=20; message=trim(message)//&
  'missing value index is the same as ix_snow or ix_soil'; return; endif
 ! allocate space for the variable check vector
 allocate(checkGotVars(size(mvar_meta)),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'allocating logical check vector'; return; endif
 checkGotVars(:) = .false.  ! initialize vector
 ! define desired scalar variables
 if(size(namesScalarDesired)/=7)then
  err=20; message=trim(message)//'expect 7 variables in namesScalarDesired'; return
 endif
 namesScalarDesired(1) = 'scalarVegetationTemp'
 namesScalarDesired(2) = 'scalarAlbedo'
 namesScalarDesired(3) = 'scalarSWE'
 namesScalarDesired(4) = 'scalarSnowDepth'
 namesScalarDesired(5) = 'scalarSfcMeltPond'
 namesScalarDesired(6) = 'scalarAquiferStorage'
 namesScalarDesired(7) = 'scalarWaterTableDepth'

 ! **********************************************************************************************
 ! (1) open files, etc.
 ! **********************************************************************************************
 ! build filename
 infile = trim(SETNGS_PATH)//trim(MODEL_INITCOND)
 ! open file
 call file_open(trim(infile),unt,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! **********************************************************************************************
 ! (2) identify the number of layers and allocate space
 ! **********************************************************************************************
 nSnow=0           ! initialize the number of snow layers
 nSoil=0           ! initialize the number of soil layers
 layer_flag=.false. ! initialize layer flag
 ! loop through file until reach the layer_tag
 do iline=1,maxLines
  read(unt,'(a)',iostat=iend)temp; if(iend/=0)then; rewind(unt); exit; endif    ! read line of data, and exit if reach the end of file
  if (temp(1:1)=='!')cycle
  ! check if reached the end of the layer definitions
  if(trim(temp)=='<end_'//trim(layer_tag)//'>')then; rewind(unt); exit; endif
  ! read layer data
  if(layer_flag) then
   ! split the line into an array of words
   call split_line(temp,chardata,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   ! check if the line contains initial conditions data (contains the word "snow" or "soil")
   do iword=1,size(chardata)
    if(chardata(iword)=='snow') nSnow = nSnow+1
    if(chardata(iword)=='soil') nSoil = nSoil+1
    if(chardata(iword)=='snow' .or. chardata(iword)=='soil') exit ! exit once read the layer type
   end do
   deallocate(chardata)
  endif  ! if in the layer section of the file
  ! check if reached the start of the layer definitions
  if (trim(temp)=='<start_'//trim(layer_tag)//'>') layer_flag=.true.
  ! check if reached the end of the file
  if (iline==maxLines) rewind(unt)
 end do  ! looping through lines
 nLayers = nSnow + nSoil
 print *, 'nLayers = ', nLayers
 ! allocate space
 call alloc_mvar(nSnow,nSoil,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 call alloc_indx(nLayers,err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 ! save the number of layers
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
 ! (2) read the scalar initial conditions
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
   jvar = get_ixmvar(trim(chardata(1)))
   if(jvar<=0)then; err=30; message=trim(message)//'variableNotFound[var='//trim(chardata(1))//']'; return; endif
   checkGotVars(jvar)=.true.
   ! read the data -- value is second
   read(chardata(2),*,iostat=err) mvar_data%var(jvar)%dat(1)
   if(err/=0)then; err=40; message=trim(message)//"problemInternalRead[data='"//trim(chardata(2))//"']"; return; endif
   deallocate(chardata)
  endif    ! if we are in the scalar part of the file
  ! check if reached the start of the scalar definitions
  if (trim(temp)=='<start_'//trim(scalar_tag)//'>') scalar_flag=.true.
  ! check if reached the end of the file
  if (iline==maxLines) rewind(unt)
 end do  ! looping through lines
 ! check if we got the desired scalar variables
 do ivar=1,size(namesScalarDesired)
  jvar=get_ixmvar(trim(namesScalarDesired(ivar)))
  if(.not.checkGotVars(jvar))then
   message=trim(message)//'initial condion undefined for variable '//trim(namesScalarDesired(ivar))
   err=20; return
  endif
 end do
 ! **********************************************************************************************
 ! (3) read the layer initial conditions
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
    jvar = get_ixmvar(trim(varnames(ivar)))
    if(jvar<=0)then; err=40; message=trim(message)//"cannotFindVariableIndex[name='"//trim(varnames(ivar))//"']"; return; endif
    ! ***** populate the data variable *****
    select case(trim(mvar_meta(jvar)%vartype))
     case('midSoil'); if(layerType==ix_soil) read(chardata(ivar),*,iostat=err) mvar_data%var(jvar)%dat(iSoil)
     case('midSnow'); if(layerType==ix_snow) read(chardata(ivar),*,iostat=err) mvar_data%var(jvar)%dat(iSnow)
     case('midToto'); read(chardata(ivar),*,iostat=err) mvar_data%var(jvar)%dat(iToto)
     case('ifcSnow'); read(chardata(ivar),*,iostat=err) mvar_data%var(jvar)%dat(iSnow-1)  ! IC = top interface
     case('ifcSoil'); read(chardata(ivar),*,iostat=err) mvar_data%var(jvar)%dat(iSoil-1)  ! IC = top interface
     case('ifcToto'); read(chardata(ivar),*,iostat=err) mvar_data%var(jvar)%dat(iToto-1)  ! IC = top interface
     case default
     err=40; message=trim(message)//"unknownInitCondType[name='"//trim(mvar_meta(jvar)%varname)//"']"; return
    endselect
    if(err/=0)then;err=40;message=trim(message)//"problemInternalRead[data='"//trim(chardata(ivar))//"']"; return; endif
   end do    ! (looping through initial conditions variables)
  endif     ! (if in the "layer" portion of the file)
  ! check if reached the start of the layer definitions
  if (trim(temp)=='<start_'//trim(layer_tag)//'>') layer_flag=.true.
 end do  ! looping through lines in the file
 ! set iLayerHeight for the bottom layer
 mvar_data%var(iLookMVAR%iLayerHeight)%dat(nLayers) = &
 mvar_data%var(iLookMVAR%iLayerHeight)%dat(nLayers-1) + mvar_data%var(iLookMVAR%mLayerDepth)%dat(nLayers)
 ! ***************************************************************************************
 ! ***************************************************************************************
 ! ensure the initial conditions are consistent with the constitutive functions
 ! ***************************************************************************************
 ! ***************************************************************************************
 ! assign pointers to model parameters
 snowfrz_scale       => mpar_data%var(iLookPARAM%snowfrz_scale)              ! scaling parameter for the snow freezing curve (K-1)
 vGn_alpha           => mpar_data%var(iLookPARAM%vGn_alpha)                  ! van Genutchen "alpha" parameter (m-1)
 vGn_n               => mpar_data%var(iLookPARAM%vGn_n)                      ! van Genutchen "n" parameter (-)
 theta_sat           => mpar_data%var(iLookPARAM%theta_sat)                  ! soil porosity (-)
 theta_res           => mpar_data%var(iLookPARAM%theta_res)                  ! soil residual volumetric water content (-)
 ! compute the maximum volumetric fraction of liquid water -- used to avoid problems of super-saturation in the moisture-based form of Richards' equation
 maxVolFracLiq = theta_sat - 1.e-4_dp
 ! compute the van Genutchen "m" parameter (-)
 vGn_m = 1._dp - 1._dp/vGn_n
 ! compute the constant in the freezing curve function (m K-1)
 kappa  = (iden_ice/iden_water)*(LH_fus/(gravity*Tfreeze))  ! NOTE: J = kg m2 s-2
 ! loop through all layers
 do iLayer=1,nLayers
  ! define short-cuts
  scalarTemp       => mvar_data%var(iLookMVAR%mLayerTemp)%dat(iLayer)          ! temperature (K) 
  scalarVolFracLiq => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(iLayer)    ! volumetric fraction of liquid water in each snow layer (-)
  scalarVolFracIce => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(iLayer)    ! volumetric fraction of ice in each snow layer (-)
  scalarLayerType  => indx_data%var(iLookINDEX%layerType)%dat(iLayer)          ! type of layer (ix_soil or ix_snow)
  ! compute liquid water equivalent of total water (liquid plus ice)
  scalarTheta = scalarVolFracIce*(iden_ice/iden_water) + scalarVolFracLiq
  ! process snow and soil separately
  select case(scalarLayerType)
   ! ** snow
   case(ix_snow)
    scalarVolFracLiq = fracliquid(scalarTemp,snowfrz_scale)*scalarTheta        ! volumetric fraction of liquid water
    scalarVolFracIce = (scalarTheta - scalarVolFracLiq)*(iden_water/iden_ice)  ! volumetric fraction of ice
   ! ** soil
   case(ix_soil)
    ! assign pointer to the matric head
    scalarMatricHead => mvar_data%var(iLookMVAR%mLayerMatricHead)%dat(iLayer-nSnow)
    ! compute the critical temperature above which all water is unfrozen (K)
    Tcrit = crit_soilT(scalarTheta,theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m)
    ! compute the matric head (m) and volumetric fractio of liquid water and ice (-) for cases of frozen soil
    if(scalarTemp < Tcrit) then
     scalarMatricHead = kappa*(scalarTemp - Tfreeze)
     scalarVolFracLiq = volFracLiq(scalarMatricHead,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
     scalarVolFracIce = (scalarTheta - scalarVolFracLiq)*(iden_water/iden_ice)
    ! compute volumetric fraction of liquid water and ice (-) for cases of unfrozen soil
    else
     select case(model_decisions(iLookDECISIONS%f_Richards)%iDecision)
      case(moisture)
       scalarVolFracLiq = min(scalarVolFracLiq,maxVolFracLiq)
       scalarMatricHead = matricHead(scalarVolFracLiq,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
      case(mixdform)
       scalarVolFracLiq = volFracLiq(scalarMatricHead,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
      case default; err=20; message=trim(message)//"unknown form of Richards' equation"; return
     end select
     scalarVolFracIce = 0._dp
    endif
   case default; err=10; message=trim(message)//'unknown case for model layer'; return
  endselect
 end do  ! (looping through layers)
 ! **********************************************************************************************
 ! deallocate variable names vector
 deallocate(varnames,chardata,stat=err)
 if(err/=0)then;err=30;message=trim(message)//'deallocating variable names vector'; return; endif
 ! deallocate variable check vector
 deallocate(checkGotVars,stat=err)
 if(err/=0)then; err=20; message=trim(message)//'deallocating logical check vector'; return; endif
 print*,'****************************************************************************************'
 print*,'mLayerDepth      ', mvar_data%var(iLookMVAR%mLayerDepth)%dat(:)
 print*,'iLayerHeight     ', mvar_data%var(iLookMVAR%iLayerHeight)%dat(:)
 print*,'mLayerTemp       ', mvar_data%var(iLookMVAR%mLayerTemp)%dat(:)
 print*,'mLayerVolFracIce ', mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(:)
 print*,'mLayerVolFracLiq ', mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(:)
 print*,'mLayerMatricHead ', mvar_data%var(iLookMVAR%mLayerMatricHead)%dat(:)
 print*,'scalarAlbedo     ', mvar_data%var(iLookMVAR%scalarAlbedo)%dat(:)
 print*,'scalarSnowDepth  ', mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(:)
 print*,'layerType        ', indx_data%var(iLookINDEX%layerType)%dat(:)
 print*,'****************************************************************************************'
 end subroutine read_icond


end module read_icond_module

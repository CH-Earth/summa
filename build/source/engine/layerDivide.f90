module layerDivide_module

! variable types
USE nrtype

! physical constants
USE multiconst,only:&
                    iden_ice,       & ! intrinsic density of ice             (kg m-3)
                    iden_water        ! intrinsic density of liquid water    (kg m-3)

! access the number of snow and soil layers
USE data_struc,only:&
                    nSnow,   & ! number of snow layers  
                    nSoil,   & ! number of soil layers  
                    nLayers    ! total number of layers

! access named variables for snow and soil
USE data_struc,only:ix_soil,ix_snow            ! named variables for snow and soil

! define look-up values for the choice of method to combine and sub-divide snow layers
USE mDecisions_module,only:&
 sameRulesAllLayers,       & ! SNTHERM option: same combination/sub-dividion rules applied to all layers
 rulesDependLayerIndex       ! CLM option: combination/sub-dividion rules depend on layer index

! define look-up values for the choice of canopy shortwave radiation method
USE mDecisions_module,only:&
 noah_mp,                  & ! full Noah-MP implementation (including albedo)
 CLM_2stream,              & ! CLM 2-stream model (see CLM documentation)
 UEB_2stream,              & ! UEB 2-stream model (Mahat and Tarboton, WRR 2011)
 NL_scatter,               & ! Simplified method Nijssen and Lettenmaier (JGR 1999)
 BeersLaw                    ! Beer's Law (as implemented in VIC)

! define look-up values for the choice of albedo method
USE mDecisions_module,only:& ! identify model options for snow albedo
 constantDecay,            & ! constant decay in snow albedo (e.g., VIC, CLASS)
 variableDecay               ! variable decay in snow albedo (e.g., BATS approach, with destructive metamorphism + soot content)

implicit none
private
public::layerDivide
interface addOneLayer
 module procedure AddOneLayer_rv, AddOneLayer_iv
end interface AddOneLayer

contains

 ! ************************************************************************************************
 ! new subroutine: add new snowfall to the system, and increase number of snow layers if needed
 ! ************************************************************************************************
 subroutine layerDivide(&
                        ! input/output: model data structures
                        model_decisions,             & ! intent(in):    model decisions
                        mpar_data,                   & ! intent(in):    model parameters
                        indx_data,                   & ! intent(inout): type of each layer
                        mvar_data,                   & ! intent(inout): model variables for a local HRU
                        ! output: error control
                        err,message)                   ! intent(out): error control
 ! --------------------------------------------------------------------------------------------------------
 ! --------------------------------------------------------------------------------------------------------
 ! access the derived types to define the data structures
 USE data_struc,only:&
                     var_d,            & ! data vector (dp)
                     var_ilength,      & ! data vector with variable length dimension (i4b)
                     var_dlength,      & ! data vector with variable length dimension (dp)
                     model_options       ! defines the model decisions
 ! access named variables defining elements in the data structures
 USE var_lookup,only:iLookTIME,iLookTYPE,iLookATTR,iLookFORCE,iLookPARAM,iLookMVAR,iLookBVAR,iLookINDEX  ! named variables for structure elements
 USE var_lookup,only:iLookDECISIONS                               ! named variables for elements of the decision structure
 ! computational modules
 USE snow_utils_module,only:fracliquid,templiquid               ! functions to compute temperature/liquid water
 implicit none
 ! --------------------------------------------------------------------------------------------------------
 ! input/output: model data structures
 type(model_options),intent(in)  :: model_decisions(:)  ! model decisions
 type(var_d),intent(in)          :: mpar_data           ! model parameters
 type(var_ilength),intent(inout) :: indx_data           ! type of each layer
 type(var_dlength),intent(inout) :: mvar_data           ! model variables for a local HRU
 ! output: error control
 integer(i4b),intent(out)        :: err                 ! error code
 character(*),intent(out)        :: message             ! error message
 ! --------------------------------------------------------------------------------------------------------
 ! variables in the data structures
 ! model decisions
 integer(i4b)                    :: ix_snowLayers       ! decision for snow combination
 ! model parameters (new snow density)
 real(dp)                        :: newSnowDenMin       ! minimum new snow density (kg m-3)
 real(dp)                        :: newSnowDenMult      ! multiplier for new snow density (kg m-3)
 real(dp)                        :: newSnowDenScal      ! scaling factor for new snow density (K)
 ! model parameters (control on the depth of snow layers)
 real(dp)                        :: zmax                ! maximum layer depth (m)
 real(dp)                        :: zmaxLayer1_lower    ! maximum layer depth for the 1st (top) layer when only 1 layer (m) 
 real(dp)                        :: zmaxLayer2_lower    ! maximum layer depth for the 2nd layer when only 2 layers (m) 
 real(dp)                        :: zmaxLayer3_lower    ! maximum layer depth for the 3rd layer when only 3 layers (m) 
 real(dp)                        :: zmaxLayer4_lower    ! maximum layer depth for the 4th layer when only 4 layers (m) 
 real(dp)                        :: zmaxLayer1_upper    ! maximum layer depth for the 1st (top) layer when > 1 layer (m) 
 real(dp)                        :: zmaxLayer2_upper    ! maximum layer depth for the 2nd layer when > 2 layers (m) 
 real(dp)                        :: zmaxLayer3_upper    ! maximum layer depth for the 3rd layer when > 3 layers (m) 
 real(dp)                        :: zmaxLayer4_upper    ! maximum layer depth for the 4th layer when > 4 layers (m) 
 ! model parameters (compute layer temperature)
 real(dp)                        :: fc_param            ! freeezing curve parameter for snow (K-1)
 ! diagnostic scalar variables
 real(dp)                        :: scalarSnowDepth     ! total snow depth (m)
 real(dp)                        :: scalarSWE           ! SWE (kg m-2)
 real(dp)                        :: scalarSnowfall      ! snowfall flux (kg m-2 s-1)
 real(dp)                        :: scalarSnowfallTemp  ! computed temperature of fresh snow (K) 
 ! model state variables (all layers)
 ! NOTE: use pointers because dimension length changes
 real(dp),pointer                :: mLayerTemp(:)       ! temperature of each layer (K)
 real(dp),pointer                :: mLayerVolFracIce(:) ! volumetric fraction of ice in each layer (-)
 real(dp),pointer                :: mLayerVolFracLiq(:) ! volumetric fraction of liquid water in each layer (-)
 ! model coordinate variables
 ! NOTE: use pointers because dimension length changes
 real(dp),pointer                :: mLayerDepth(:)      ! depth of the layer (m)
 real(dp),pointer                :: mLayerHeight(:)     ! height of the layer mid-point (m)
 real(dp),pointer                :: iLayerHeight(:)     ! height of the layer interface (m)
 ! model index variables
 ! NOTE: use pointers because dimension length changes
 integer(i4b),pointer            :: layerType(:)        ! type of the layer (ix_soil or ix_snow)
 ! --------------------------------------------------------------------------------------------------------
 ! define local variables
 character(LEN=256)              :: cmessage            ! error message of downwind routine
 integer(i4b)                    :: iLayer              ! layer index
 integer(i4b)                    :: jLayer              ! layer index
 integer(i4b)                    :: kLayer              ! layer index
 integer(i4b)                    :: ivar                ! variable index
 real(dp),dimension(4)           :: zmax_lower          ! lower value of maximum layer depth 
 real(dp),dimension(4)           :: zmax_upper          ! upper value of maximum layer depth 
 real(dp)                        :: zmaxCheck           ! value of zmax for a given snow layer
 integer(i4b)                    :: nCheck              ! number of layers to check to divide
 logical(lgt)                    :: createLayer         ! flag to indicate we are creating a new snow layer 
 real(dp)                        :: surfaceLayerSoilTemp  ! temperature of the top soil layer (K)
 real(dp)                        :: maxFrozenSnowTemp   ! maximum temperature when effectively all water is frozen (K)
 real(dp),parameter              :: unfrozenLiq=0.01_dp ! unfrozen liquid water used to compute maxFrozenSnowTemp (-)
 real(dp)                        :: volFracWater        ! volumetric fraction of total water, liquid and ice (-)
 real(dp)                        :: fracLiq             ! fraction of liquid water (-)
 integer(i4b),parameter          :: ixVisible=1         ! named variable to define index in array of visible part of the spectrum
 integer(i4b),parameter          :: ixNearIR=2          ! named variable to define index in array of near IR part of the spectrum
 ! --------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="layerDivide/"
 ! --------------------------------------------------------------------------------------------------------
 ! associate variables in the data structures
 associate(&
 ! model decisions
 ix_snowLayers          => model_decisions(iLookDECISIONS%snowLayers)%iDecision, & ! decision for snow combination
 ! model parameters (new snow density)
 newSnowDenMin          => mpar_data%var(iLookPARAM%newSnowDenMin),            & ! minimum new snow density (kg m-3)
 newSnowDenMult         => mpar_data%var(iLookPARAM%newSnowDenMult),           & ! multiplier for new snow density (kg m-3)
 newSnowDenScal         => mpar_data%var(iLookPARAM%newSnowDenScal),           & ! scaling factor for new snow density (K)
 ! model parameters (control the depth of snow layers)
 zmax                   => mpar_data%var(iLookPARAM%zmax),                     & ! maximum layer depth (m)
 zmaxLayer1_lower       => mpar_data%var(iLookPARAM%zmaxLayer1_lower),         & ! maximum layer depth for the 1st (top) layer when only 1 layer (m) 
 zmaxLayer2_lower       => mpar_data%var(iLookPARAM%zmaxLayer2_lower),         & ! maximum layer depth for the 2nd layer when only 2 layers (m) 
 zmaxLayer3_lower       => mpar_data%var(iLookPARAM%zmaxLayer3_lower),         & ! maximum layer depth for the 3rd layer when only 3 layers (m) 
 zmaxLayer4_lower       => mpar_data%var(iLookPARAM%zmaxLayer4_lower),         & ! maximum layer depth for the 4th layer when only 4 layers (m) 
 zmaxLayer1_upper       => mpar_data%var(iLookPARAM%zmaxLayer1_upper),         & ! maximum layer depth for the 1st (top) layer when > 1 layer (m) 
 zmaxLayer2_upper       => mpar_data%var(iLookPARAM%zmaxLayer2_upper),         & ! maximum layer depth for the 2nd layer when > 2 layers (m) 
 zmaxLayer3_upper       => mpar_data%var(iLookPARAM%zmaxLayer3_upper),         & ! maximum layer depth for the 3rd layer when > 3 layers (m) 
 zmaxLayer4_upper       => mpar_data%var(iLookPARAM%zmaxLayer4_upper),         & ! maximum layer depth for the 4th layer when > 4 layers (m) 
 ! model parameters (compute layer temperature)
 fc_param               => mpar_data%var(iLookPARAM%snowfrz_scale),            & ! freezing curve parameter for snow (K-1)
 ! diagnostic scalar variables
 scalarSnowfall         => mvar_data%var(iLookMVAR%scalarSnowfall)%dat(1),     & ! snowfall flux (kg m-2 s-1)
 scalarSnowfallTemp     => mvar_data%var(iLookMVAR%scalarSnowfallTemp)%dat(1), & ! computed temperature of fresh snow (K)
 scalarSnowDepth        => mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1),    & ! total snow depth (m)
 scalarSWE              => mvar_data%var(iLookMVAR%scalarSWE)%dat(1)           & ! SWE (kg m-2)
 )  ! end associate statement

 ! assign pointers to model state variables
 mLayerDepth            => mvar_data%var(iLookMVAR%mLayerDepth)%dat           ! depth of the layer (m)
 mLayerTemp             => mvar_data%var(iLookMVAR%mLayerTemp)%dat            ! temperature of each layer (K)
 mLayerVolFracIce       => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat      ! volumetric fraction of ice in each layer (-)
 mLayerVolFracLiq       => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat      ! volumetric fraction of liquid water in each layer (-)
 ! assign local pointers to the model index structures
 layerType              => indx_data%var(iLookINDEX%layerType)%dat            ! layer type (ix_soil or ix_snow)

 ! --------------------------------------------------------------------------------------------------------

 ! identify algorithmic control parameters to syb-divide and combine snow layers
 zmax_lower = (/zmaxLayer1_lower, zmaxLayer2_lower, zmaxLayer3_lower, zmaxLayer4_lower/)
 zmax_upper = (/zmaxLayer1_upper, zmaxLayer2_upper, zmaxLayer3_upper, zmaxLayer4_upper/) 

 ! ***** special case of no snow layers
 if(nSnow==0)then

  ! check if create the first snow layer
  select case(ix_snowLayers)
   case(sameRulesAllLayers);    createLayer = (scalarSnowDepth > zmax)
   case(rulesDependLayerIndex); createLayer = (scalarSnowDepth > zmaxLayer1_lower)
   case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
  end select ! (option to combine/sub-divide snow layers)

  ! ** create a new snow layer
  if(createLayer)then

   ! add a layer to all model variables
   iLayer=0 ! (layer to divide: 0 is the special case of "snow without a layer")
   call addModelLayer(mvar_data,indx_data,iLayer,err,cmessage) 
   if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif

   ! re-assign pointers to the model state variables
   ! NOTE: need to do this here, since state vectors have just been modified
   mLayerTemp       => mvar_data%var(iLookMVAR%mLayerTemp)%dat           ! temperature of each layer (K)
   mLayerVolFracIce => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat     ! volumetric fraction of ice in each layer (-)
   mLayerVolFracLiq => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat     ! volumetric fraction of liquid water in each layer (-)

   ! compute surface layer temperature
   surfaceLayerSoilTemp = mLayerTemp(2)    ! temperature of the top soil layer (K)
   maxFrozenSnowTemp    = templiquid(unfrozenLiq,fc_param)               ! snow temperature at fraction "unfrozenLiq" (K)
   mLayerTemp(1)        = min(maxFrozenSnowTemp,surfaceLayerSoilTemp)    ! snow temperature  (K)

   ! compute the fraction of liquid water associated with the layer temperature
   fracLiq      = fracliquid(mLayerTemp(1),fc_param)

   ! compute volumeteric fraction of liquid water and ice
   volFracWater = (scalarSWE/scalarSnowDepth)/iden_water  ! volumetric fraction of total water (liquid and ice)
   mLayerVolFracIce(1) = (1._dp - fracLiq)*volFracWater*(iden_water/iden_ice)   ! volumetric fraction of ice (-)
   mLayerVolFracLiq(1) =          fracLiq *volFracWater                         ! volumetric fraction of liquid water (-)

   ! initialize albedo
   ! NOTE: albedo is computed within the Noah-MP radiation routine
   if(model_decisions(iLookDECISIONS%canopySrad)%iDecision /= noah_mp)then
    select case(model_decisions(iLookDECISIONS%alb_method)%iDecision)
     ! (constant decay rate -- albedo the same for all spectral bands)
     case(constantDecay)
      mvar_data%var(iLookMVAR%scalarSnowAlbedo)%dat(1)          = mpar_data%var(iLookPARAM%albedoMax)
      mvar_data%var(iLookMVAR%spectralSnowAlbedoDiffuse)%dat(:) = mpar_data%var(iLookPARAM%albedoMax)
     ! (variable decay rate)
     case(variableDecay)
      mvar_data%var(iLookMVAR%spectralSnowAlbedoDiffuse)%dat(ixVisible) = mpar_data%var(iLookPARAM%albedoMaxVisible)
      mvar_data%var(iLookMVAR%spectralSnowAlbedoDiffuse)%dat(ixNearIR)  = mpar_data%var(iLookPARAM%albedoMaxNearIR)
      mvar_data%var(iLookMVAR%scalarSnowAlbedo)%dat(1)                  = (        mpar_data%var(iLookPARAM%Frad_vis))*mpar_data%var(iLookPARAM%albedoMaxVisible) + &
                                                                          (1._dp - mpar_data%var(iLookPARAM%Frad_vis))*mpar_data%var(iLookPARAM%albedoMaxNearIR)
     case default; err=20; message=trim(message)//'unable to identify option for snow albedo'; return
    end select  ! identify option for snow albedo
    ! set direct albedo to diffuse albedo
    mvar_data%var(iLookMVAR%spectralSnowAlbedoDirect)%dat(:) = mvar_data%var(iLookMVAR%spectralSnowAlbedoDiffuse)%dat(:)
   endif  ! (if NOT using the Noah-MP radiation routine)

   ! check
   print*, trim(message)
   do kLayer=1,nLayers
    write(*,'(i4,1x,4(f9.3,1x))') layerType(kLayer), mLayerDepth(kLayer), mLayerTemp(kLayer), mLayerVolFracIce(kLayer), mLayerVolFracLiq(kLayer)
   end do
   print*, 'created a new layer, nSnow = ', count(indx_data%var(iLookINDEX%layerType)%dat==ix_snow)
   print*, 'snow albedo = ', mvar_data%var(iLookMVAR%scalarSnowAlbedo)%dat(1)

   !pause ' check layer sub-division'

  endif  ! if creating a new layer
  return
 endif

 ! end special case of nSnow=0
 ! ********************************************************************************************************************
 ! ********************************************************************************************************************

 ! check
 !print*, 'before sub-division'
 !do kLayer=1,nLayers
 ! write(*,'(i4,1x,4(f9.3,1x))') layerType(kLayer), mLayerDepth(kLayer), mLayerTemp(kLayer), mLayerVolFracIce(kLayer), mLayerVolFracLiq(kLayer)
 !end do
 !if(scalarSnowDepth > 0.5_dp) pause ' deep snow'

 ! ***** sub-divide snow layers, if necessary

 ! identify the number of layers to check for need for sub-division
 select case(ix_snowLayers)
  case(sameRulesAllLayers);    nCheck = nSnow
  case(rulesDependLayerIndex); nCheck = min(nSnow,4)  ! the depth of the 5th layer, if it exists, does not have a maximum value
  case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
 end select ! (option to combine/sub-divide snow layers)

 ! loop through all layers, and sub-divide a given layer, if necessary
 do iLayer=1,nCheck

  ! identify the maximum depth of the layer
  select case(ix_snowLayers)
   case(sameRulesAllLayers);    zmaxCheck = zmax
   case(rulesDependLayerIndex)
    if(iLayer == nSnow)then
     zmaxCheck = zmax_lower(iLayer)
    else
     zmaxCheck = zmax_upper(iLayer)
    endif
   case default; err=20; message=trim(message)//'unable to identify option to combine/sub-divide snow layers'; return
  end select ! (option to combine/sub-divide snow layers)

  ! check the need to sub-divide
  if(mvar_data%var(iLookMVAR%mLayerDepth)%dat(iLayer) > zmaxCheck)then

   ! add a layer to all model variables
   call addModelLayer(mvar_data,indx_data,iLayer,err,cmessage)  ! adds model layer to the index BELOW the layer that is too thick
   if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif

   ! identify the number of snow and soil layers, and check all is a-OK
   nSnow   = count(layerType==ix_snow)
   nSoil   = count(layerType==ix_soil)
   nLayers = nSnow + nSoil

   ! check
   mLayerTemp       => mvar_data%var(iLookMVAR%mLayerTemp)%dat           ! temperature of each layer (K)
   mLayerDepth      => mvar_data%var(iLookMVAR%mLayerDepth)%dat          ! depth of each layer (m)
   mLayerVolFracIce => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat     ! volumetric fraction of ice in each layer (-)
   mLayerVolFracLiq => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat     ! volumetric fraction of liquid water in each layer (-)
   !print*, 'after sub-division'
   !do kLayer=1,nLayers
   ! write(*,'(i4,1x,4(f9.3,1x))') layerType(kLayer), mLayerDepth(kLayer), mLayerTemp(kLayer), mLayerVolFracIce(kLayer), mLayerVolFracLiq(kLayer)
   !end do
   print*, 'created a new layer, nSnow = ', count(indx_data%var(iLookINDEX%layerType)%dat==ix_snow)
   !pause ' check layer sub-division'

   exit  ! NOTE: only sub-divide one layer per substep

  endif   ! (if sub-dividing layer)

 end do  ! (looping through layers)

 ! end associate variables in data structure
 end associate
 
 end subroutine layerDivide


 ! *********************************************************************************************
 ! new subroutine: add an additional layer to all model vectors
 ! *********************************************************************************************
 subroutine addModelLayer(mvar_data,indx_data,ix_divide,err,message)
 ! provide access to variables in the data structures
 USE data_struc,only:mvar_meta,indx_meta      ! metadata
 USE data_struc,only:var_ilength,var_dlength  ! data vectors with variable length dimension
 USE var_lookup,only:iLookPARAM,iLookMVAR,iLookINDEX  ! named variables for structure elements
 implicit none
 ! ---------------------------------------------------------------------------------------------
 ! input/output: data structures
 type(var_dlength),intent(inout) :: mvar_data ! model variables for a local HRU
 type(var_ilength),intent(inout) :: indx_data ! type of model layer
 ! input: snow layer indices
 integer(i4b),intent(in)         :: ix_divide ! index of the layer to divide
 ! output: error control
 integer(i4b),intent(out)        :: err       ! error code
 character(*),intent(out)        :: message   ! error message
 ! ---------------------------------------------------------------------------------------------
 ! variables in the data structures
 ! diagnostic variables
 real(dp)                        :: scalarSnowDepth     ! total snow depth (m)
 ! model coordinate variables
 ! NOTE: use pointers because dimension length changes
 real(dp),pointer                :: mLayerDepth(:)      ! depth of the layer (m)
 real(dp),pointer                :: mLayerHeight(:)     ! height of the layer mid-point (m)
 real(dp),pointer                :: iLayerHeight(:)     ! height of the layer interface (m)
 ! model index variables
 ! NOTE: use pointers because dimension length changes
 integer(i4b),pointer            :: layerType(:)        ! type of the layer (ix_soil or ix_snow)
 ! ---------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                    :: ivar                ! index of model variable
 integer(i4b)                    :: jLayer              ! index of model layer
 integer(i4b)                    :: ix_lower            ! lower bound of the vector
 integer(i4b)                    :: ix_upper            ! upper bound of the vector
 logical(lgt)                    :: stateVariable       ! .true. if variable is a state variable
 real(dp)                        :: depthOriginal       ! original layer depth before sub-division (m)
 real(dp),parameter              :: fracTop=0.5_dp      ! fraction of old layer used for the top layer
 character(LEN=256)              :: cmessage            ! error message of downwind routine
 ! ---------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='addModelLayer/'

 ! associate variables in data structure
 associate(&
 scalarSnowDepth        => mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1)    & ! total snow depth (m)
 )  ! associate

 ! ***** add a layer to each model variable
 do ivar=1,size(mvar_data%var)

  ! define bounds
  select case(trim(mvar_meta(ivar)%vartype))
   case('midSnow'); ix_lower=1; ix_upper=nSnow
   case('midToto'); ix_lower=1; ix_upper=nLayers
   case('ifcSnow'); ix_lower=0; ix_upper=nSnow
   case('ifcToto'); ix_lower=0; ix_upper=nLayers
   case default; cycle
  end select

  ! identify whether it is a state variable
  select case(trim(mvar_meta(ivar)%varname))
   case('mLayerDepth','mLayerTemp','mLayerVolFracIce','mLayerVolFracLiq'); stateVariable=.true.
   case default; stateVariable=.false.
  end select

  ! add an additional layer -- only get to here if snow in the layer
  call AddOneLayer(mvar_data%var(ivar)%dat,ix_lower,ix_upper,ix_divide,stateVariable,err,cmessage)
  if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif

 end do  ! looping through variables


 ! ***** modify the layer indices
 call AddOneLayer(indx_data%var(iLookINDEX%layerType)%dat,1,nLayers,ix_divide,.false.,err,cmessage)
 if(err/=0)then; err=10; message=trim(message)//trim(cmessage); return; endif
 indx_data%var(iLookINDEX%layerType)%dat(1:nSnow+1)         = ix_snow
 indx_data%var(iLookINDEX%layerType)%dat(nSnow+2:nLayers+1) = ix_soil
 nLayers = nLayers + 1

 ! assign pointers to model coordinate variables
 mLayerDepth      => mvar_data%var(iLookMVAR%mLayerDepth)%dat          ! depth of the layer (m)
 mLayerHeight     => mvar_data%var(iLookMVAR%mLayerHeight)%dat         ! height of the layer mid-point (m)
 iLayerHeight     => mvar_data%var(iLookMVAR%iLayerHeight)%dat         ! height of the layer interface (m)
 layerType        => indx_data%var(iLookINDEX%layerType)%dat           ! type of each layer (ix_snow or ix_soil)

 ! ***** modify the layer depth
 if(ix_divide==0)then ! no layers exist currently
  mLayerDepth(1) = scalarSnowDepth
 else ! layers already exist
  depthOriginal = mLayerDepth(ix_divide)
  mLayerDepth(ix_divide)   = fracTop*depthOriginal
  mLayerDepth(ix_divide+1) = (1._dp - fracTop)*depthOriginal
 endif

 ! check
 if(scalarSnowDepth - sum(mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow)) < epsilon(scalarSnowDepth))then
  message=trim(message)//'problem sub-dividing snow layer'
  err=20; return
 endif

 ! ***** re-set coordinate variables
 iLayerHeight(0) = -scalarSnowDepth
 do jLayer=1,nLayers
  iLayerHeight(jLayer) = iLayerHeight(jLayer-1) + mLayerDepth(jLayer)
  mLayerHeight(jLayer) = (iLayerHeight(jLayer-1) + iLayerHeight(jLayer))/2._dp
 end do

 ! end associate variables in data structure
 end associate

 end subroutine addModelLayer



 ! ************************************************************************************************
 ! new subroutine: add an additional snow layer
 ! ************************************************************************************************
 subroutine AddOneLayer_rv(datavec,ix_lower,ix_upper,ix_divide,stateVariable,err,message)
 ! Returns a new vector which has one more element than the input vector
 !  -- optionally copies data from the original vector to the new vector for elements (2:n)->(3:n+1),
 !      and copies element 1 into elements 1:2, and copies element 0 into element 0
 implicit none
 ! dummies
 real(dp),pointer,intent(inout)     :: datavec(:)    ! the original and the new vector
 integer(i4b),intent(in)            :: ix_lower      ! lower bound of the old vector
 integer(i4b),intent(in)            :: ix_upper      ! upper bound of the old vector
 integer(i4b),intent(in)            :: ix_divide     ! index of the layer to divide
 logical(lgt),intent(in)            :: stateVariable ! .true. if a state variable
 integer(i4b),intent(out)           :: err           ! error code
 character(*),intent(out)           :: message       ! error message
 ! locals
 real(dp)                           :: tempvec(ix_lower:ix_upper)  ! temporary vector
 real(dp),parameter                 :: missingReal=-9999._dp
 ! initialize error control
 err=0; message='AddOneLayer_rv/'
 ! check the data vector is associated
 if(.not.associated(datavec))then; err=20; message='data vector is not associated'; return; endif 
 ! assign the data vector to the temporary vector
 tempvec=datavec
 ! reallocate space for the new vector
 deallocate(datavec,stat=err)
 if(err/=0)then; err=20; message='problem in attempt to deallocate memory for data vector'; return; endif
 allocate(datavec(ix_lower:ix_upper+1),stat=err)
 if(err/=0)then; err=20; message='problem in attempt to reallocate memory for data vector'; return; endif
 if(stateVariable)then
  if(ix_upper > 0)then  ! (only copy data if the vector exists -- can be a variable for snow, with no layers)
   if(ix_divide > 0)then
    datavec(1:ix_divide) = tempvec(1:ix_divide)  ! copy data
    datavec(ix_divide+1) = tempvec(ix_divide)    ! repeat data for the sub-divided layer
   endif
   if(ix_upper > ix_divide) datavec(ix_divide+2:ix_upper+1) = tempvec(ix_divide+1:ix_upper)
  endif
 else
  datavec = missingReal
 endif
 end subroutine AddOneLayer_rv

 subroutine AddOneLayer_iv(datavec,ix_lower,ix_upper,ix_divide,stateVariable,err,message)
 ! Returns a new vector which has one more element than the input vector
 !  -- optionally copies data from the original vector to the new vector for elements (2:n)->(3:n+1),
 !      and copies element 1 into elements 1:2, and copies element 0 into element 0
 implicit none
 ! dummies
 integer(i4b),pointer,intent(inout) :: datavec(:)    ! the original and the new vector
 integer(i4b),intent(in)            :: ix_lower      ! lower bound of the old vector
 integer(i4b),intent(in)            :: ix_upper      ! upper bound of the old vector
 integer(i4b),intent(in)            :: ix_divide     ! index of the layer to divide
 logical(lgt),intent(in)            :: stateVariable ! .true. if a state variable
 integer(i4b),intent(out)           :: err           ! error code
 character(*),intent(out)           :: message       ! error message
 ! locals
 integer(i4b)                       :: tempvec(ix_lower:ix_upper)  ! temporary vector
 integer(i4b),parameter             :: missingInteger=-9999
 ! initialize error control
 err=0; message='AddOneLayer_iv/'
 ! check the data vector is associated
 if(.not.associated(datavec))then; err=20; message='data vector is not associated'; return; endif
 ! assign the data vector to the temporary vector
 tempvec=datavec
 ! reallocate space for the new vector
 deallocate(datavec,stat=err)
 if(err/=0)then; err=20; message='problem in attempt to deallocate memory for data vector'; return; endif
 allocate(datavec(ix_lower:ix_upper+1),stat=err)
 if(err/=0)then; err=20; message='problem in attempt to reallocate memory for data vector'; return; endif
 if(stateVariable)then
  if(ix_upper > 0)then  ! (only copy data if the vector exists -- can be a variable for snow, with no layers)
   if(ix_divide > 0)then
    datavec(1:ix_divide) = tempvec(1:ix_divide)  ! copy data
    datavec(ix_divide+1) = tempvec(ix_divide)    ! repeat data for the sub-divided layer
   endif
   if(ix_upper > ix_divide) datavec(ix_divide+2:ix_upper+1) = tempvec(ix_divide+1:ix_upper)
  endif
 else
  datavec = missingInteger
 endif
 end subroutine AddOneLayer_iv

end module layerDivide_module

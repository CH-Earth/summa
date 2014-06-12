module systemSolv_module
! data types
USE nrtype
! layer types
USE data_struc,only:ix_soil,ix_snow ! named variables for snow and soil
! access the number of snow and soil layers
USE data_struc,only:&
                    nSnow,        & ! number of snow layers  
                    nSoil,        & ! number of soil layers  
                    nLayers         ! total number of layers
! constants
USE multiconst,only:&
                    gravity,      & ! acceleration of gravity              (m s-2)
                    Tfreeze,      & ! temperature at freezing              (K)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)
! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:  &
 qbaseTopmodel,              & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                  & ! a big bucket (lumped aquifer model)
 noExplicit                    ! no explicit groundwater parameterization
! look-up values for the form of Richards' equation
USE mDecisions_module,only:  &
 moisture,                   & ! moisture-based form of Richards' equation
 mixdform                      ! mixed form of Richards' equation
! look-up values for the choice of boundary conditions for hydrology
USE mDecisions_module,only:  &
 prescribedHead,             & ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
 funcBottomHead,             & ! function of matric head in the lower-most layer
 freeDrainage,               & ! free drainage
 liquidFlux,                 & ! liquid water flux
 zeroFlux                      ! zero flux
implicit none
private
public::systemSolv
! number of soil levels
integer(i4b)        :: nLevels      ! number of soil layers to use in the soil hydrology routine
! control parameters
real(dp),parameter  :: valueMissing=-9999._dp     ! missing value
real(dp),parameter  :: verySmall=1.e-6_dp         ! a very small number
contains

 ! ************************************************************************************************
 ! new subroutine: run the coupled energy-mass model for one timestep
 ! ************************************************************************************************
 subroutine systemSolv(&
                       ! input
                       dt,            & ! time step (s)
                       maxiter,       & ! maximum number of iterations
                       firstSubstep,  & ! flag to denote first sub-step
                       computeVegFlux,& ! flag to denote if computing energy flux over vegetation
                       ! output
                       niter,         & ! number of iterations taken
                       err,message)     ! error control
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:attr_data,type_data,mpar_data,forc_data,mvar_data,indx_data    ! data structures
 USE var_lookup,only:iLookATTR,iLookTYPE,iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX ! named variables for structure elements
 ! model decision structures
 USE data_struc,only:model_decisions    ! model decision structure
 USE var_lookup,only:iLookDECISIONS     ! named variables for elements of the decision structure
 ! snow densification
 USE snwDensify_module,only:snwDensify
 implicit none
 ! input
 real(dp),intent(in)                  :: dt                       ! time step (seconds)
 integer(i4b),intent(in)              :: maxiter                  ! maximum number of iterations
 logical(lgt),intent(in)              :: firstSubStep             ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(in)              :: computeVegFlux           ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! output
 integer(i4b),intent(out)             :: niter                    ! number of iterations
 integer(i4b),intent(out)             :: err                      ! error code
 character(*),intent(out)             :: message                  ! error message
 ! local
 character(LEN=256)                   :: cmessage                 ! error message of downwind routine
 real(dp)                             :: scalarCanairTempNew      ! temperature of the canopy air space at the end of the sub-step (K)
 real(dp)                             :: scalarCanopyTempNew      ! temperature of the vegetation canopy at the end of the sub-step (K)
 real(dp)                             :: scalarCanopyIceNew       ! mass of ice on the vegetation canopy at the end of the sub-step (kg m-2)
 real(dp)                             :: scalarCanopyLiqNew       ! mass of liquid water on the vegetation canopy at the end of the sub-step (kg m-2)
 real(dp),allocatable                 :: mLayerTempNew(:)         ! temperature of each snow/soil layer at the end of the sub-step (K)
 real(dp),allocatable                 :: mLayerVolFracIceNew(:)   ! volumetric fraction of ice at the end of the sub-step (-)
 real(dp),allocatable                 :: mLayerVolFracLiqNew(:)   ! volumetric fraction of liquid water at the end of the sub-step (-)
 real(dp),allocatable                 :: mLayerMatricHeadNew(:)   ! matric head at the end of the sub-step (m)
 real(dp)                             :: scalarAquiferStorageNew  ! aquifer storage at the end of the sub-step (m)
 real(dp)                             :: scalarCanopyWater        ! total storage of wate rin the canopy; liquidwater plus ice (kg m-2)
 real(dp)                             :: volSub                   ! volumetric sublimation (kg m-3)
 real(dp)                             :: testSWE                  ! test SWE (kg m-2) 

 ! initialize error control
 err=0; message="systemSolv/"

 ! call wrapper routine for the system solver
 call systemSolv_muster(&
                        ! input: model control
                        dt,                                                   & ! time step (s)
                        maxiter,                                              & ! maximum number of iterations
                        firstSubstep,                                         & ! flag to denote first sub-step
                        computeVegFlux,                                       & ! flag to denote if computing energy flux over vegetation
                        ! input/output: model state variables (vegetation canopy)
                        mvar_data%var(iLookMVAR%scalarCanairTemp)%dat(1),     & ! intent(inout): temperature of the canopy air space (K)
                        mvar_data%var(iLookMVAR%scalarCanopyTemp)%dat(1),     & ! intent(inout): temperature of the vegetation canopy (K)
                        mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1),      & ! intent(inout): mass of ice on the vegetation canopy (kg m-2)
                        mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1),      & ! intent(inout): mass of liquid water on the vegetation canopy (kg m-2)
                        ! input/output: model state variables (snow and soil domains)
                        mvar_data%var(iLookMVAR%mLayerTemp)%dat,              & ! intent(inout): temperature of each snow/soil layer (K)
                        mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat,        & ! intent(inout): volumetric fraction of ice (-)
                        mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat,        & ! intent(inout): volumetric fraction of liquid water (-)
                        mvar_data%var(iLookMVAR%mLayerMatricHead)%dat,        & ! intent(inout): matric head (m)
                        mvar_data%var(iLookMVAR%scalarAquiferStorage)%dat(1), & ! intent(inout): aquifer storage (m)
                        ! output: model control
                        niter,                                                & ! intent(out): number of iterations taken
                        err,cmessage)                                           ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

 end subroutine systemSolv

 ! ****************************************************************************************************************************************
 ! ****************************************************************************************************************************************
 ! ****************************************************************************************************************************************
 ! ****************************************************************************************************************************************
 ! ****************************************************************************************************************************************

 ! ************************************************************************************************
 ! wrapper subroutine: solve the model equations
 ! ************************************************************************************************
 subroutine systemSolv_muster(&
                              ! input: model control
                              dt,                   & ! intent(in): time step (s)
                              maxiter,              & ! intent(in): maximum number of iterations
                              firstSubstep,         & ! intent(in): flag to denote first sub-step
                              computeVegFlux,       & ! intent(in): flag to denote if computing energy flux over vegetation
                              ! input/output: model state variables (vegetation canopy)
                              scalarCanairTemp,     & ! intent(inout): temperature of the canopy air space (K)
                              scalarCanopyTemp,     & ! intent(inout): temperature of the vegetation canopy (K)
                              scalarCanopyIce,      & ! intent(inout): mass of ice on the vegetation canopy (kg m-2)
                              scalarCanopyLiq,      & ! intent(inout): mass of liquid water on the vegetation canopy (kg m-2)
                              ! input/output: model state variables (snow and soil domains)
                              mLayerTemp,           & ! intent(inout): temperature of each snow/soil layer (K)
                              mLayerVolFracIce,     & ! intent(inout): volumetric fraction of ice (-)
                              mLayerVolFracLiq,     & ! intent(inout): volumetric fraction of liquid water (-)
                              mLayerMatricHead,     & ! intent(inout): matric head (m)
                              scalarAquiferStorage, & ! intent(inout): aquifer storage (m)
                              ! output: model control
                              niter,                & ! intent(out): number of iterations taken
                              err,cmessage)           ! intent(out): error control
 implicit none
 ! input: model control
 real(dp),intent(in)            :: dt                           ! time step (seconds)
 integer(i4b),intent(in)        :: maxiter                      ! maximum number of iterations
 logical(lgt),intent(in)        :: firstSubStep                 ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(in)        :: computeVegFlux               ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! input/output: model state variables (vegetation canopy)
 real(dp),intent(inout)         :: scalarCanairTemp             ! temperature of the canopy air space (K)
 real(dp),intent(inout)         :: scalarCanopyTemp             ! temperature of the vegetation canopy (K)
 real(dp),intent(inout)         :: scalarCanopyIce              ! mass of ice on the vegetation canopy (kg m-2)
 real(dp),intent(inout)         :: scalarCanopyLiq              ! mass of liquid water on the vegetation canopy (kg m-2)
 ! input/output: model state variables (snow and soil domains)
 real(dp),intent(inout)         :: mLayerTemp(:)                ! temperature of each snow/soil layer (K)
 real(dp),intent(inout)         :: mLayerVolFracIce(:)          ! volumetric fraction of ice (-)
 real(dp),intent(inout)         :: mLayerVolFracLiq(:)          ! volumetric fraction of liquid water (-)
 real(dp),intent(inout)         :: mLayerMatricHead(:)          ! matric head (m)
 real(dp),intent(inout)         :: scalarAquiferStorage         ! aquifer storage (m)
 ! output: model control
 integer(i4b),intent(out)       :: niter                        ! number of iterations
 integer(i4b),intent(out)       :: err                          ! error code
 character(*),intent(out)       :: message                      ! error message
 ! ------------------------------------------------------------------------------------------------------
 ! ------------------------------------------------------------------------------------------------------
 ! ------------------------------------------------------------------------------------------------------
 ! general local variables
 integer(i4b)                   :: iter                         ! iteration index
 ! ------------------------------------------------------------------------------------------------------
 ! * trial state variables
 ! ------------------------------------------------------------------------------------------------------
 ! trial state variables (vegetation canopy)
 real(dp)                       :: scalarCanairTempTrial        ! trial value for temperature of the canopy air space (K)
 real(dp)                       :: scalarCanopyTempTrial        ! trial value for temperature of the vegetation canopy (K)
 real(dp)                       :: scalarCanopyIceTrial         ! trial value for mass of ice on the vegetation canopy (kg m-2)
 real(dp)                       :: scalarCanopyLiqTrial         ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
 ! trial state variables (snow and soil domains)
 real(dp),dimension(nLayers)    :: mLayerTempTrial(:)           ! trial value for temperature of each snow/soil layer (K)
 real(dp),dimension(nLayers)    :: mLayerVolFracIceTrial(:)     ! trial value for volumetric fraction of ice (-)
 real(dp),dimension(nLayers)    :: mLayerVolFracLiqTrial(:)     ! trial value for volumetric fraction of liquid water (-)
 real(dp),dimension(nSoil)      :: mLayerMatricHeadTrial(:)     ! trial value for matric head (m)
 real(dp)                       :: scalarAquiferStorageTrial    ! trial value for aquifer storage (m)
 ! ------------------------------------------------------------------------------------------------------
 ! * model fluxes
 ! ------------------------------------------------------------------------------------------------------
 ! energy fluxes and derivatives for the vegetation domain
 real(dp)                       :: canairNetFlux                ! net energy flux for the canopy air space (W m-2)
 real(dp)                       :: canopyNetFlux                ! net energy flux for the vegetation canopy (W m-2)
 real(dp)                       :: groundNetFlux                ! net energy flux for the ground surface (W m-2)
 real(dp)                       :: dCanairNetFlux_dCanairTemp   ! derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp)                       :: dCanairNetFlux_dCanopyTemp   ! derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp)                       :: dCanairNetFlux_dGroundTemp   ! derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
 real(dp)                       :: dCanopyNetFlux_dCanairTemp   ! derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp)                       :: dCanopyNetFlux_dCanopyTemp   ! derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp)                       :: dCanopyNetFlux_dGroundTemp   ! derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
 real(dp)                       :: dGroundNetFlux_dCanairTemp   ! derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp)                       :: dGroundNetFlux_dCanopyTemp   ! derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp)                       :: dGroundNetFlux_dGroundTemp   ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
 ! liquid water fluxes associated with transpiration
 real(dp)                       :: scalarCanopyTranspiration    ! canopy transpiration (kg m-2 s-1)
 real(dp)                       :: scalarCanopyEvaporation      ! canopy evaporation/condensation (kg m-2 s-1)
 real(dp)                       :: scalarGroundEvaporation      ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
 ! energy fluxes and derivatives for the snow and soil domains
 real(dp),dimension(0:nLayers)  :: iLayerNrgFlux                ! energy flux at the layer interfaces (W m-2)
 real(dp),dimension(0:nLayers)  :: dNrgFlux_dTempAbove          ! derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
 real(dp),dimension(0:nLayers)  :: dNrgFlux_dTempBelow          ! derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
 ! liquid water fluxes and derivatives for the vegetation domain
 real(dp)                       :: scalarThroughfallRain        ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp)                       :: scalarCanopyLiqDrainage      ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
 real(dp)                       :: scalarCanopyLiqDrainageDeriv ! derivative in canopy drainage w.r.t. canopy liquid water (s-1)
 ! liquid water fluxes and derivatives for the snow domain 
 real(dp),dimension(0:nSnow)    :: iLayerLiqFluxSnow            ! vertical liquid water flux at layer interfaces (m s-1)
 real(dp),dimension(0:nSnow)    :: iLayerLiqFluxSnowDeriv       ! derivative in vertical liquid water flux at layer interfaces (m s-1)
 ! liquid water fluxes and derivatives for the soil domain
 real(dp),dimension(0:nSoil)    :: iLayerLiqFluxSoil            ! liquid flux at soil layer interfaces (m s-1)
 real(dp),dimension(nSoil)      :: mLayerTranspire              ! transpiration loss from each soil layer (m s-1)
 real(dp),dimension(nSoil)      :: mLayerBaseflow               ! baseflow from each soil layer -- only compute at the start of the step (m s-1)
 real(dp),dimension(0:nSoil)    :: dq_dStateAbove               ! change in the flux in layer interfaces w.r.t. state variables in the layer above
 real(dp),dimension(0:nSoil)    :: dq_dStateBelow               ! change in the flux in layer interfaces w.r.t. state variables in the layer below
 ! liquid water fluxes and derivatives for the aquifer
 real(dp)                       :: scalarAquiferTranspire       ! transpiration loss from the aquifer at the start-of-step (m s-1)
 real(dp)                       :: scalarAquiferRecharge        ! recharge to the aquifer (m s-1)
 real(dp)                       :: scalarAquiferBaseflow        ! total baseflow from the aquifer (m s-1)
 real(dp)                       :: scalarAquiferRechargeDeriv   ! derivative in recharge to the aquifer w.r.t. aquifer storage (s-1)
 real(dp)                       :: scalarAquiferBaseflowDeriv   ! derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
 ! ------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='systemSolv_muster/'

 ! *****
 ! (0) PRELIMINARIES...
 ! ********************

 ! get an initial canopy temperature if veg just starts protruding through snow on the ground
 if(computeVegFlux)then
  ! (NOTE: if canopy temperature is below absolute zero then canopy was previously buried by snow)
  if(scalarCanopyTemp < 0._dp .or. scalarCanairTemp < 0._dp)then
   ! check there is snow (there really has to be)
   if(nSnow == 0)then
    message=trim(message)//'no snow when canopy temperature or canopy air temperature is undefined -- canopy temps can only be undefined when buried with snow'
    err=20; return
   endif
   ! set canopy temperature to the temperature of the top snow layer + small offset to check derivative calculations
   scalarCanairTemp = mLayerTemp(1) + 0.1_dp
   scalarCanopyTemp = mLayerTemp(1) + 0.1_dp
  endif  ! (if canopy temperature undefined -- means canopy previously buried with snow)
 endif  ! (if computing vegetation fluxes -- canopy exposed)

 ! initialize state variables for the vegetation canopy
 scalarCanairTempTrial     = scalarCanairTemp     ! trial value of the canopy air space temperature (K)
 scalarCanopyTempTrial     = scalarCanopyTemp     ! trial value of canopy temperature (K)
 scalarCanopyIceTrial      = scalarCanopyIce      ! trial value of mass of ice on the vegetation canopy (kg m-2)
 scalarCanopyLiqTrial      = scalarCanopyLiq      ! trial value of mass of liquid water on the vegetation canopy (kg m-2)

 ! initialize state variables for the snow and soil domains
 mLayerTempTrial           = mLayerTemp            ! trial value of temperature (K)
 mLayerMatricHeadTrial     = mLayerMatricHead      ! trial value of matric head (m)
 mLayerVolFracLiqTrial     = mLayerVolFracLiq      ! trial value of volumetric fraction of liquid water (-)
 mLayerVolFracIceTrial     = mLayerVolFracIce      ! trial value of volumetric fraction of ice (-)
 scalarAquiferStorageTrial = scalarAquiferStorage  ! trial value of aquifer storage (m)


 ! *****
 ! (1) MAIN ITERATION LOOP...
 ! **************************

 ! iterate
 do iter=1,maxiter


  ! *****
  ! (A) COMPUTE MODEL FLUXES AND DERIVATIVES...
  ! *******************************************
  call computFlux(&
                  ! input: model control
                  iter,                                   & ! intent(in): iteration index
                  firstSubStep,                           & ! intent(in): flag to indicate if we are processing the first sub-step
                  computeVegFlux,                         & ! intent(in): flag to indicate if we need to compute fluxes over vegetation
                  ! input: model state variables for the vegetation canopy
                  scalarCanairTempTrial,                  & ! intent(in): trial value of the canopy air space temperature (K)
                  scalarCanopyTempTrial,                  & ! intent(in): trial value of canopy temperature (K)
                  scalarCanopyIceTrial,                   & ! intent(in): trial value of mass of ice on the vegetation canopy (kg m-2)
                  scalarCanopyLiqTrial,                   & ! intent(in): trial value of mass of liquid water on the vegetation canopy (kg m-2)
                  ! input: model state variables for the snow and soil domain
                  mLayerTempTrial,                        & ! intent(in): trial value of temperature (K)
                  mLayerMatricHeadTrial,                  & ! intent(in): trial value of matric head (m)
                  mLayerVolFracLiqTrial,                  & ! intent(in): trial value of volumetric fraction of liquid water (-)
                  mLayerVolFracIceTrial,                  & ! intent(in): trial value of volumetric fraction of ice (-)
                  scalarAquiferStorageTrial,              & ! intent(in): trial value of aquifer storage (m)
                  ! output: energy fluxes and derivatives for the vegetation domain
                  canairNetFlux,                          & ! intent(out): net energy flux for the canopy air space (W m-2)
                  canopyNetFlux,                          & ! intent(out): net energy flux for the vegetation canopy (W m-2)
                  groundNetFlux,                          & ! intent(out): net energy flux for the ground surface (W m-2)
                  dCanairNetFlux_dCanairTemp,             & ! intent(out): derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
                  dCanairNetFlux_dCanopyTemp,             & ! intent(out): derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
                  dCanairNetFlux_dGroundTemp,             & ! intent(out): derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
                  dCanopyNetFlux_dCanairTemp,             & ! intent(out): derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
                  dCanopyNetFlux_dCanopyTemp,             & ! intent(out): derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
                  dCanopyNetFlux_dGroundTemp,             & ! intent(out): derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
                  dGroundNetFlux_dCanairTemp,             & ! intent(out): derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
                  dGroundNetFlux_dCanopyTemp,             & ! intent(out): derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
                  dGroundNetFlux_dGroundTemp,             & ! intent(out): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
                  ! output: liquid water fluxes associated with transpiration
                  scalarCanopyTranspiration,              & ! intent(out): canopy transpiration (kg m-2 s-1)
                  scalarCanopyEvaporation,                & ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
                  scalarGroundEvaporation,                & ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
                  ! output: energy fluxes and derivatives for the snow and soil domains
                  iLayerNrgFlux,                          & ! intent(out): energy flux at the layer interfaces (W m-2)
                  dFlux_dTempAbove,                       & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
                  dFlux_dTempBelow,                       & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
                  ! output: liquid water fluxes and derivatives for the vegetation domain
                  scalarThroughfallRain,                  & ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
                  scalarCanopyLiqDrainage,                & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
                  scalarCanopyLiqDrainageDeriv,           & ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
                  ! output: liquid water fluxes and derivatives for the snow domain
                  iLayerLiqFluxSnow,                      & ! intent(out): vertical liquid water flux at layer interfaces (m s-1)
                  iLayerLiqFluxSnowDeriv,                 & ! intent(out): derivative in vertical liquid water flux at layer interfaces (m s-1)
                  ! output: liquid water fluxes and derivatives for the soil domain
                  iLayerLiqFluxSoil,                      & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                  mLayerTranspire,                        & ! intent(out): transpiration loss from each soil layer (m s-1)
                  mLayerBaseflow,                         & ! intent(inout): baseflow from each soil layer -- only compute at the start of the step (m s-1)
                  dq_dStateAbove,                         & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer above (m s-1)
                  dq_dStateBelow,                         & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer below (m s-1)
                  ! output: liquid water fluxes and derivatives for the aquifer storage
                  scalarAquiferTranspire,                 & ! intent(out): transpiration loss from the aquifer at the start-of-step (m s-1)
                  scalarAquiferRecharge,                  & ! intent(out): recharge to the aquifer (m s-1)
                  scalarAquiferBaseflow,                  & ! intent(out): total baseflow from the aquifer (m s-1)
                  scalarAquiferRechargeDeriv,             & ! intent(out): derivative in recharge to the aquifer w.r.t. aquifer storage (s-1)
                  scalarAquiferBaseflowDeriv,             & ! intent(out): derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
                  ! output: error control
                  err,cmessage)                             ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


 end do  ! iterating


 end subroutine systemSolv



 ! *********************************************************************************************************************************************************************
 ! *********************************************************************************************************************************************************************
 ! ***** PRIVATE SUBROUTINES *******************************************************************************************************************************************
 ! *********************************************************************************************************************************************************************
 ! *********************************************************************************************************************************************************************
 
 ! ************************************************************************************************
 ! private subroutine: compute model fluxes
 ! ************************************************************************************************
 subroutine computFlux(&
                       iter,                                   & ! intent(in): iteration index
                       firstSubStep,                           & ! intent(in): flag to indicate if we are processing the first sub-step
                       computeVegFlux,                         & ! intent(in): flag to indicate if we need to compute fluxes over vegetation
                       ! input: model state variables for the vegetation canopy
                       scalarCanairTempTrial,                  & ! intent(in): trial value of the canopy air space temperature (K)
                       scalarCanopyTempTrial,                  & ! intent(in): trial value of canopy temperature (K)
                       scalarCanopyIceTrial,                   & ! intent(in): trial value of mass of ice on the vegetation canopy (kg m-2)
                       scalarCanopyLiqTrial,                   & ! intent(in): trial value of mass of liquid water on the vegetation canopy (kg m-2)
                       ! input: model state variables for the snow and soil domain
                       mLayerTempTrial,                        & ! intent(in): trial value of temperature (K)
                       mLayerMatricHeadTrial,                  & ! intent(in): trial value of matric head (m)
                       mLayerVolFracLiqTrial,                  & ! intent(in): trial value of volumetric fraction of liquid water (-)
                       mLayerVolFracIceTrial,                  & ! intent(in): trial value of volumetric fraction of ice (-)
                       scalarAquiferStorageTrial,              & ! intent(in): trial value of aquifer storage (m)
                       ! output: energy fluxes and derivatives for the vegetation domain
                       canairNetFlux,                          & ! intent(out): net energy flux for the canopy air space (W m-2)
                       canopyNetFlux,                          & ! intent(out): net energy flux for the vegetation canopy (W m-2)
                       groundNetFlux,                          & ! intent(out): net energy flux for the ground surface (W m-2)
                       dCanairNetFlux_dCanairTemp,             & ! intent(out): derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
                       dCanairNetFlux_dCanopyTemp,             & ! intent(out): derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
                       dCanairNetFlux_dGroundTemp,             & ! intent(out): derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
                       dCanopyNetFlux_dCanairTemp,             & ! intent(out): derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
                       dCanopyNetFlux_dCanopyTemp,             & ! intent(out): derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
                       dCanopyNetFlux_dGroundTemp,             & ! intent(out): derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
                       dGroundNetFlux_dCanairTemp,             & ! intent(out): derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
                       dGroundNetFlux_dCanopyTemp,             & ! intent(out): derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
                       dGroundNetFlux_dGroundTemp,             & ! intent(out): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
                       ! output: liquid water fluxes associated with transpiration
                       scalarCanopyTranspiration,              & ! intent(out): canopy transpiration (kg m-2 s-1)
                       scalarCanopyEvaporation,                & ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
                       scalarGroundEvaporation,                & ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
                       ! output: energy fluxes and derivatives for the snow and soil domains
                       iLayerNrgFlux,                          & ! intent(out): energy flux at the layer interfaces (W m-2)
                       dNrgFlux_dTempAbove,                    & ! intent(out): derivatives in the energy flux w.r.t. temperature in the layer above (W m-2 K-1)
                       dNrgFlux_dTempBelow,                    & ! intent(out): derivatives in the energy flux w.r.t. temperature in the layer below (W m-2 K-1)
                       ! output: liquid water fluxes and derivatives for the vegetation domain
                       scalarThroughfallRain,                  & ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
                       scalarCanopyLiqDrainage,                & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
                       scalarCanopyLiqDrainageDeriv,           & ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
                       ! output: liquid water fluxes and derivatives for the snow domain
                       iLayerLiqFluxSnow,                      & ! intent(inout): vertical liquid water flux at layer interfaces (m s-1)
                       iLayerLiqFluxSnowDeriv,                 & ! intent(out): derivative in vertical liquid water flux at layer interfaces (m s-1)
                       ! output: liquid water fluxes and derivatives for the soil domain
                       iLayerLiqFluxSoil,                      & ! intent(inout): liquid fluxes at layer interfaces (m s-1)
                       mLayerTranspire,                        & ! intent(out): transpiration loss from each soil layer (m s-1)
                       mLayerBaseflow,                         & ! intent(inout): baseflow from each soil layer -- only compute at the start of the step (m s-1)
                       dq_dStateAbove,                         & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer above (m s-1)
                       dq_dStateBelow,                         & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer below (m s-1)
                       ! output: liquid water fluxes and derivatives for the aquifer storage
                       scalarAquiferTranspire,                 & ! intent(out): transpiration loss from the aquifer at the start-of-step (m s-1)
                       scalarAquiferRecharge,                  & ! intent(out): recharge to the aquifer (m s-1)
                       scalarAquiferBaseflow,                  & ! intent(out): total baseflow from the aquifer (m s-1)
                       scalarAquiferRechargeDeriv,             & ! intent(out): derivative in recharge to the aquifer w.r.t. aquifer storage (s-1)
                       scalarAquiferBaseflowDeriv,             & ! intent(out): derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
                       ! output: error control
                       err,message)                              ! intent(out): error control
 implicit none
 ! input: model control
 integer(i4b),intent(in)        :: iter                          ! iteration index
 logical(lgt),intent(in)        :: firstSubStep                  ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(in)        :: computeVegFlux                ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! input: model state variables for the vegetation canopy 
 real(dp),intent(in)            :: scalarCanairTempTrial         ! trial value of canopy air space temperature (K)
 real(dp),intent(in)            :: scalarCanopyTempTrial         ! trial value of canopy temperature (K)
 real(dp),intent(in)            :: scalarCanopyIceTrial          ! trial value of mass of ice on the vegetation canopy (kg m-2)
 real(dp),intent(in)            :: scalarCanopyLiqTrial          ! trial value of mass of liquid water on the vegetation canopy (kg m-2)
 ! input: model state variables for the snow and soil domain
 real(dp),intent(in)            :: mLayerTempTrial(:)            ! trial value of temperature of each snow/soil layer (K)
 real(dp),intent(in)            :: mLayerMatricHeadTrial(:)      ! trial value of matric head in each layer (m)
 real(dp),intent(in)            :: mLayerVolFracLiqTrial(:)      ! trial value of volumetric fraction of liquid water (-)
 real(dp),intent(in)            :: mLayerVolFracIceTrial(:)      ! trial value of volumetric fraction of ice (-)
 real(dp),intent(in)            :: scalarAquiferStorageTrial     ! aquifer storage (m)
 ! output: energy fluxes and derivatives for the vegetation domain
 real(dp),intent(out)           :: canairNetFlux                 ! net energy flux for the canopy air space (W m-2)
 real(dp),intent(out)           :: canopyNetFlux                 ! net energy flux for the vegetation canopy (W m-2)
 real(dp),intent(out)           :: groundNetFlux                 ! net energy flux for the ground surface (W m-2)
 real(dp),intent(out)           :: dCanairNetFlux_dCanairTemp    ! derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp),intent(out)           :: dCanairNetFlux_dCanopyTemp    ! derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp),intent(out)           :: dCanairNetFlux_dGroundTemp    ! derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
 real(dp),intent(out)           :: dCanopyNetFlux_dCanairTemp    ! derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp),intent(out)           :: dCanopyNetFlux_dCanopyTemp    ! derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp),intent(out)           :: dCanopyNetFlux_dGroundTemp    ! derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
 real(dp),intent(out)           :: dGroundNetFlux_dCanairTemp    ! derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp),intent(out)           :: dGroundNetFlux_dCanopyTemp    ! derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp),intent(out)           :: dGroundNetFlux_dGroundTemp    ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
 ! output: liquid water fluxes associated with transpiration
 real(dp),intent(out)           :: scalarCanopyTranspiration     ! canopy transpiration (kg m-2 s-1)
 real(dp),intent(out)           :: scalarCanopyEvaporation       ! canopy evaporation/condensation (kg m-2 s-1)
 real(dp),intent(out)           :: scalarGroundEvaporation       ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
 ! output: energy fluxes and derivatives for the snow and soil domains
 real(dp),intent(out)           :: iLayerNrgFlux(0:)             ! energy flux at the layer interfaces (W m-2)
 real(dp),intent(out)           :: dNrgFlux_dTempAbove(0:)       ! derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
 real(dp),intent(out)           :: dNrgFlux_dTempBelow(0:)       ! derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
 ! liquid water fluxes and derivatives for the vegetation domain
 real(dp),intent(out)           :: scalarThroughfallRain         ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp),intent(out)           :: scalarCanopyLiqDrainage       ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
 real(dp),intent(out)           :: scalarCanopyLiqDrainageDeriv  ! derivative in canopy drainage w.r.t. canopy liquid water (s-1)
 ! output: liquid water fluxes and derivatives for the snow domain
 real(dp),intent(inout)         :: iLayerLiqFluxSnow(0:)         ! vertical liquid water flux at layer interfaces (m s-1)
 real(dp),intent(out)           :: iLayerLiqFluxSnowDeriv(0:)    ! derivative in vertical liquid water flux at layer interfaces (m s-1)
 ! output: liquid water fluxes and derivatives for the soil domain
 real(dp),intent(inout)         :: iLayerLiqFluxSoil(0:)       ! liquid flux at soil layer interfaces (m s-1)
 real(dp),intent(out)           :: mLayerTranspire(:)            ! transpiration loss from each soil layer (m s-1)
 real(dp),intent(inout)         :: mLayerBaseflow(:)             ! baseflow from each soil layer -- only compute at the start of the step (m s-1)
 real(dp),intent(out)           :: dq_dStateAbove(0:)            ! change in the flux in layer interfaces w.r.t. state variables in the layer above
 real(dp),intent(out)           :: dq_dStateBelow(0:)            ! change in the flux in layer interfaces w.r.t. state variables in the layer below
 ! output: liquid water fluxes and derivatives for the aquifer storage
 real(dp),intent(out)           :: scalarAquiferTranspire        ! transpiration loss from the aquifer at the start-of-step (m s-1)
 real(dp),intent(out)           :: scalarAquiferRecharge         ! recharge to the aquifer (m s-1)
 real(dp),intent(out)           :: scalarAquiferBaseflow         ! total baseflow from the aquifer (m s-1)
 real(dp),intent(out)           :: scalarAquiferRechargeDeriv    ! derivative in recharge to the aquifer w.r.t. aquifer storage (s-1)
 real(dp),intent(out)           :: scalarAquiferBaseflowDeriv    ! derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
 ! output: error control
 integer(i4b),intent(out)       :: err                           ! error code
 character(*),intent(out)       :: message                       ! error message
 ! ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='computFlux/'

 ! *****
 ! (0) PRELIMINARIES...
 ! ********************

 ! initialize liquid water fluxes throughout the snow and soil domains
 ! NOTE: used in the energy routines, which is called before the hydrology routines
 if(iter==1)then
  if(nSnow > 0)&
  iLayerLiqFluxSnow(0:nSnow) = 0._dp
  iLayerLiqFluxSoil(0:nSoil) = 0._dp
 endif

 ! *****
 ! (1) CALCULATE ENERGY FLUXES OVER VEGETATION...
 ! **********************************************
 call vegNrgFlux(&
                 ! input: model control
                 iter,                                   & ! intent(in): iteration index
                 firstSubStep,                           & ! intent(in): flag to indicate if we are processing the first sub-step
                 computeVegFlux,                         & ! intent(in): flag to indicate if we need to compute fluxes over vegetation
                 ! input: model state variables
                 scalarCanairTempTrial,                  & ! intent(in): trial value of the canopy air space temperature (K)
                 scalarCanopyTempTrial,                  & ! intent(in): trial value of canopy temperature (K)
                 mLayerTempTrial(1),,                    & ! intent(in): trial value of ground temperature (K)
                 scalarCanopyIceTrial,                   & ! intent(in): trial value of mass of ice on the vegetation canopy (kg m-2)
                 scalarCanopyLiqTrial,                   & ! intent(in): trial value of mass of liquid water on the vegetation canopy (kg m-2)
                 upperBoundTemp,                         & ! intent(in): temperature of the upper boundary (K) --> NOTE: use air temperature
                 ! output: liquid water fluxes associated with evaporation/transpiration
                 scalarCanopyTranspiration,              & ! intent(out): canopy transpiration (kg m-2 s-1)
                 scalarCanopyEvaporation,                & ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
                 scalarGroundEvaporation,                & ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
                 ! output: fluxes
                 canairNetFlux,                          & ! intent(out): net energy flux for the canopy air space (W m-2)
                 canopyNetFlux,                          & ! intent(out): net energy flux for the vegetation canopy (W m-2)
                 groundNetFlux,                          & ! intent(out): net energy flux for the ground surface (W m-2)
                 ! output: flux derivatives
                 dCanairNetFlux_dCanairTemp,             & ! intent(out): derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
                 dCanairNetFlux_dCanopyTemp,             & ! intent(out): derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
                 dCanairNetFlux_dGroundTemp,             & ! intent(out): derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
                 dCanopyNetFlux_dCanairTemp,             & ! intent(out): derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
                 dCanopyNetFlux_dCanopyTemp,             & ! intent(out): derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
                 dCanopyNetFlux_dGroundTemp,             & ! intent(out): derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
                 dGroundNetFlux_dCanairTemp,             & ! intent(out): derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
                 dGroundNetFlux_dCanopyTemp,             & ! intent(out): derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
                 dGroundNetFlux_dGroundTemp,             & ! intent(out): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
                 ! output: error control
                 err,cmessage)                             ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

 ! *****
 ! (2) CALCULATE ENERGY FLUXES THROUGH THE SNOW-SOIL DOMAIN...
 ! ***********************************************************
 call ssdNrgFlux(&
                 ! input: fluxes and derivatives at the upper boundary
                 groundNetFlux,                          & ! intent(in): total flux at the ground surface (W m-2)
                 dGroundNetFlux_dGroundTemp,             & ! intent(in): derivative in total ground surface flux w.r.t. ground temperature (W m-2 K-1)
                 ! input: liquid water fluxes throughout the snow and soil domains
                 iLayerLiqFluxSnow,                      & ! intent(in): liquid flux at the interface of each snow layer (m s-1)
                 iLayerLiqFluxSoil,                      & ! intent(in): liquid flux at the interface of each soil layer (m s-1)
                 ! input: trial value of model state variabes
                 mLayerTempTrial,                        & ! intent(in): trial temperature at the current iteration (K)
                 ! output: fluxes and derivatives at all layer interfaces
                 iLayerNrgFlux,                          & ! intent(out): energy flux at the layer interfaces (W m-2)
                 dFlux_dTempAbove,                       & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
                 dFlux_dTempBelow,                       & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
                 ! output: error control
                 err,cmessage)                             ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! *****
 ! (3) CALCULATE THE LIQUID FLUX THROUGH VEGETATION...
 ! ***************************************************
 call vegLiqFlux(&
                 ! input
                 computeVegFlux,                         & ! intent(in): flag to denote if computing energy flux over vegetation
                 scalarCanopyLiqTrial,                   & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                 scalarRainfall,                         & ! intent(in): rainfall (kg m-2 s-1)
                 scalarCanopyLiqMax,                     & ! intent(in): maximum storage before canopy drainage begins (kg m-2 s-1)
                 scalarCanopyDrainageCoeff,              & ! intent(in): canopy drainage coefficient (s-1)
                 ! output
                 scalarThroughfallRain,                  & ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
                 scalarCanopyLiqDrainage,                & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
                 scalarCanopyLiqDrainageDeriv,           & ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
                 err,cmessage)                             ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! *****
 ! (4) CALCULATE THE LIQUID FLUX THROUGH SNOW...
 ! *********************************************
 if(nSnow > 0)then
  call snowLiqFlx(&
                  ! input: model control
                  iter,                                  & ! intent(in): iteration index
                  ! input: forcing for the snow domain
                  scalarThroughfallRain,                 & ! intent(in): rain that reaches the snow surface without ever touching vegetation (kg m-2 s-1)
                  scalarCanopyLiqDrainage,               & ! intent(in): liquid drainage from the vegetation canopy (kg m-2 s-1)
                  ! input: model state vector
                  mLayerVolFracLiqTrial(1:nSnow),        & ! intent(in): trial value of volumetric fraction of liquid water at the current iteration (-)
                  ! output: fluxes and derivatives
                  iLayerLiqFluxSnow(0:nSnow),            & ! intent(out): vertical liquid water flux at layer interfaces (m s-1)
                  iLayerLiqFluxSnowDeriv(0:nSnow),       & ! intent(out): derivative in vertical liquid water flux at layer interfaces (m s-1)
                  ! output: error control
                  err,cmessage)                            ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 endif

 ! *****
 ! (5) CALCULATE THE LIQUID FLUX THROUGH SOIL...
 ! *********************************************
 call soilLiqFlx(&
                 ! input: model control
                 iter,                                   & ! intent(in): iteration index
                 deriv_desired,                          & ! intent(in): flag indicating if derivatives are desired
                 ! input: trial state variables
                 mLayerMatricHeadTrial,                  & ! intent(in): matric head (m)
                 mLayerVolFracLiqTrial(nSnow+1:nLayers), & ! intent(in): volumetric fraction of liquid water (-)
                 mLayerVolFracIceTrial(nSnow+1:nLayers), & ! intent(in): volumetric fraction of ice (-)
                 scalarAquiferStorageTrial,              & ! intent(in): aquifer storage at the start of the step (m)
                 ! output: fluxes at layer interfaces
                 iLayerLiqFluxSoil,                      & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                 mLayerTranspire,                        & ! intent(out): transpiration loss from each soil layer (m s-1)
                 mLayerBaseflow,                         & ! intent(inout): baseflow from each soil layer -- only compute at the start of the step (m s-1)
                 ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                 dq_dStateAbove,                         & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer above (m s-1)
                 dq_dStateBelow,                         & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer below (m s-1)
                 ! output: fluxes and derivatives for the aquifer storage
                 scalarAquiferTranspire,                 & ! intent(out): transpiration loss from the aquifer at the start-of-step (m s-1)
                 scalarAquiferRecharge,                  & ! intent(out): recharge to the aquifer (m s-1)
                 scalarAquiferBaseflow,                  & ! intent(out): total baseflow from the aquifer (m s-1)
                 scalarAquiferRechargeDeriv,             & ! intent(out): derivative in recharge to the aquifer w.r.t. aquifer storage (s-1)
                 scalarAquiferBaseflowDeriv,             & ! intent(out): derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
                 ! output: error control
                 err,cmessage)                             ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine computFlux





end module systemSolv_module

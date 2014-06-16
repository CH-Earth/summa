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
                    LH_vap,       & ! latent heat of vaporization          (J kg-1)
                    LH_sub,       & ! latent heat of sublimation           (J kg-1)
                    Cp_air,       & ! specific heat of air                 (J kg-1 K-1)
                    iden_air,     & ! intrinsic density of air             (kg m-3)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)
! look-up values for the choice of groundwater representation (local-column, or single-basin)
USE mDecisions_module,only:  &
 localColumn,                & ! separate groundwater representation in each local soil column
 singleBasin                   ! single groundwater store over the entire basin
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
! control parameters
real(dp),parameter  :: valueMissing=-9999._dp     ! missing value
real(dp),parameter  :: verySmall=1.e-6_dp         ! a very small number
real(dp),parameter  :: dx = 1.e-8_dp              ! finite difference increment
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
 real(dp)                             :: canopyDepth              ! depth of the vegetation canopy (m)
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

 ! define canopy depth (m)
 canopyDepth = mpar_data%var(iLookPARAM%heightCanopyTop) - mpar_data%var(iLookPARAM%heightCanopyBottom)

 print*, 'nLayers = ', nLayers

 ! call wrapper routine for the system solver
 call systemSolv_muster(&
                        ! input: model control
                        dt,                                                      & ! intent(in): time step (s)
                        maxiter,                                                 & ! intent(in): maximum number of iterations
                        firstSubstep,                                            & ! intent(in): flag to denote first sub-step
                        computeVegFlux,                                          & ! intent(in): flag to denote if computing energy flux over vegetation
                        ! input: model decisions
                        model_decisions(iLookDECISIONS%groundwatr)%iDecision,    & ! intent(in): groundwater parameterization
                        model_decisions(iLookDECISIONS%spatial_gw)%iDecision,    & ! intent(in): spatial representation of groundwater (local-column or single-basin)
                        ! input: domain boundary conditions
                        forc_data%var(iLookFORCE%airtemp),                       & ! intent(in): temperature of the upper boundary of the snow and soil domains (K)
                        mvar_data%var(iLookMVAR%scalarRainfall)%dat(1),          & ! intent(in): rainfall rate (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarSfcMeltPond)%dat(1),       & ! intent(in): ponded water caused by melt of the "snow without a layer" (kg m-2)
                        ! input: diagnostic variables
                        canopyDepth,                                             & ! intent(in): depth of the vegetation canopy (m)
                        mvar_data%var(iLookMVAR%mLayerDepth)%dat,                & ! intent(in): depth of each layer in the snow-soil sub-domain (m)
                        mvar_data%var(iLookMVAR%scalarBulkVolHeatCapVeg)%dat(1), & ! intent(in): bulk volumetric heat capacity of vegetation (J m-3 K-1)
                        mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat,         & ! intent(in): bulk volumetric heat capacity in each snow and soil layer (J m-3 K-1)
                        ! input/output: model state variables (vegetation canopy)
                        mvar_data%var(iLookMVAR%scalarCanairTemp)%dat(1),        & ! intent(inout): temperature of the canopy air space (K)
                        mvar_data%var(iLookMVAR%scalarCanopyTemp)%dat(1),        & ! intent(inout): temperature of the vegetation canopy (K)
                        mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1),         & ! intent(inout): mass of ice on the vegetation canopy (kg m-2)
                        mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1),         & ! intent(inout): mass of liquid water on the vegetation canopy (kg m-2)
                        ! input/output: model state variables (snow and soil domains)
                        mvar_data%var(iLookMVAR%mLayerTemp)%dat,                 & ! intent(inout): temperature of each snow/soil layer (K)
                        mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat,           & ! intent(inout): volumetric fraction of ice (-)
                        mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat,           & ! intent(inout): volumetric fraction of liquid water (-)
                        mvar_data%var(iLookMVAR%mLayerMatricHead)%dat,           & ! intent(inout): matric head (m)
                        mvar_data%var(iLookMVAR%scalarAquiferStorage)%dat(1),    & ! intent(inout): aquifer storage (m)
                        ! output: model control
                        niter,                                                   & ! intent(out): number of iterations taken
                        err,cmessage)                                              ! intent(out): error control
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
                              dt,                      & ! intent(in): time step (s)
                              maxiter,                 & ! intent(in): maximum number of iterations
                              firstSubstep,            & ! intent(in): flag to denote first sub-step
                              computeVegFlux,          & ! intent(in): flag to denote if computing energy flux over vegetation
                              ! input: model decisions
                              ixGroundwater,           & ! intent(in): choice of groundwater parameterization
                              ixSpatialGroundwater,    & ! intent(in): spatial representation of groundwater (local-column or single-basin)
                              ! input: domain boundary conditions
                              upperBoundTemp,          & ! intent(in): temperature of the upper boundary of the snow and soil domains (K)
                              scalarRainfall,          & ! intent(in): rainfall (kg m-2 s-1)
                              scalarSfcMeltPond,       & ! intent(in): ponded water caused by melt of the "snow without a layer" (kg m-2)
                              ! input: diagnostic variables
                              canopyDepth,             & ! intent(in): depth of the vegetation canopy (m)
                              mLayerDepth,             & ! intent(in): depth of each layer in the snow-soil sub-domain (m)
                              scalarBulkVolHeatCapVeg, & ! intent(in): bulk volumetric heat capacity of vegetation (J m-3 K-1)
                              mLayerVolHtCapBulk,      & ! intent(in): bulk volumetric heat capacity in each snow and soil layer (J m-3 K-1)
                              ! input/output: model state variables (vegetation canopy)
                              scalarCanairTemp,        & ! intent(inout): temperature of the canopy air space (K)
                              scalarCanopyTemp,        & ! intent(inout): temperature of the vegetation canopy (K)
                              scalarCanopyIce,         & ! intent(inout): mass of ice on the vegetation canopy (kg m-2)
                              scalarCanopyLiq,         & ! intent(inout): mass of liquid water on the vegetation canopy (kg m-2)
                              ! input/output: model state variables (snow and soil domains)
                              mLayerTemp,              & ! intent(inout): temperature of each snow/soil layer (K)
                              mLayerVolFracIce,        & ! intent(inout): volumetric fraction of ice (-)
                              mLayerVolFracLiq,        & ! intent(inout): volumetric fraction of liquid water (-)
                              mLayerMatricHead,        & ! intent(inout): matric head (m)
                              scalarAquiferStorage,    & ! intent(inout): aquifer storage (m)
                              ! output: model control
                              niter,                   & ! intent(out): number of iterations taken
                              err,message)               ! intent(out): error control
 ! provide access to the flux modules
 USE vegnrgflux_module,only:vegnrgflux                ! compute energy fluxes over vegetation
 USE ssdnrgflux_module,only:ssdnrgflux                ! compute energy fluxes throughout the snow and soil subdomains
 USE vegliqflux_module,only:vegliqflux                ! compute liquid water fluxes through vegetation
 USE snowliqflx_module,only:snowliqflx                ! compute liquid water fluxes through snow 
 USE soilliqflx_module,only:soilliqflx                ! compute liquid water fluxes through soil
 USE groundwatr_module,only:satstorage                ! compute saturated storage in the soil profile
 USE groundwatr_module,only:soilbsflow                ! compute total baseflow from the soil profile
 USE groundwatr_module,only:disaggflow                ! diaggregate total baseflow to individual soil layers
 USE groundwatr_module,only:aquifrflux                ! compute liquid water fluxes in the aquifer
 ! provide access to the solver modules
 USE matrixSolv_module,only:matrixSolv                ! solve full matrix
 implicit none
 ! input: model control
 real(dp),intent(in)            :: dt                           ! time step (seconds)
 integer(i4b),intent(in)        :: maxiter                      ! maximum number of iterations
 logical(lgt),intent(in)        :: firstSubStep                 ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(in)        :: computeVegFlux               ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! input: model decisions
 integer(i4b),intent(in)        :: ixGroundwater                ! choice of groundwater parameterization
 integer(i4b),intent(in)        :: ixSpatialGroundwater         ! spatial representation of groundwater (local-column or single-basin)
 ! input: domain boundary conditions
 real(dp),intent(in)            :: upperBoundTemp               ! temperature of the upper boundary of the snow and soil domains (K)
 real(dp),intent(in)            :: scalarRainfall               ! rainfall (kg m-2 s-1)
 real(dp),intent(in)            :: scalarSfcMeltPond            ! ponded water caused by melt of the "snow without a layer" (kg m-2)
 ! input: diagnostic variables
 real(dp),intent(in)            :: canopyDepth                  ! depth of the vegetation canopy (m)
 real(dp),intent(in)            :: mLayerDepth(:)               ! depth of each layer in the snow-soil sub-domain (m)
 real(dp),intent(in)            :: scalarBulkVolHeatCapVeg      ! bulk volumetric heat capacity of vegetation (J m-3 K-1)
 real(dp),intent(in)            :: mLayerVolHtCapBulk(:)        ! bulk volumetric heat capacity in each snow and soil layer (J m-3 K-1)
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
 integer(i4b)                   :: iLayer                       ! index of model layer
 integer(i4b)                   :: jLayer                       ! index of model layer
 character(LEN=256)             :: cmessage                     ! error message of downwind routine
 integer(i4b)                   :: local_ixGroundwater          ! local index for groundwater representation
 ! ------------------------------------------------------------------------------------------------------
 ! * trial state variables
 ! ------------------------------------------------------------------------------------------------------
 ! trial state variables (vegetation canopy)
 real(dp)                       :: scalarCanairTempTrial        ! trial value for temperature of the canopy air space (K)
 real(dp)                       :: scalarCanopyTempTrial        ! trial value for temperature of the vegetation canopy (K)
 real(dp)                       :: scalarCanopyIceTrial         ! trial value for mass of ice on the vegetation canopy (kg m-2)
 real(dp)                       :: scalarCanopyLiqTrial         ! trial value for mass of liquid water on the vegetation canopy (kg m-2)
 ! trial state variables (snow and soil domains)
 real(dp),dimension(nLayers)    :: mLayerTempTrial              ! trial value for temperature of each snow/soil layer (K)
 real(dp),dimension(nLayers)    :: mLayerVolFracIceTrial        ! trial value for volumetric fraction of ice (-)
 real(dp),dimension(nLayers)    :: mLayerVolFracLiqTrial        ! trial value for volumetric fraction of liquid water (-)
 real(dp),dimension(nSoil)      :: mLayerMatricHeadTrial        ! trial value for matric head (m)
 real(dp)                       :: scalarAquiferStorageTrial    ! trial value for aquifer storage (m)
 ! ------------------------------------------------------------------------------------------------------
 ! * model fluxes and derivatives
 ! ------------------------------------------------------------------------------------------------------
 ! energy fluxes and derivatives for the vegetation domain
 real(dp)                       :: canairNetNrgFlux             ! net energy flux for the canopy air space (W m-2)
 real(dp)                       :: canopyNetNrgFlux             ! net energy flux for the vegetation canopy (W m-2)
 real(dp)                       :: groundNetNrgFlux             ! net energy flux for the ground surface (W m-2)
 real(dp)                       :: dCanairNetFlux_dCanairTemp   ! derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp)                       :: dCanairNetFlux_dCanopyTemp   ! derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp)                       :: dCanairNetFlux_dGroundTemp   ! derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
 real(dp)                       :: dCanopyNetFlux_dCanairTemp   ! derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp)                       :: dCanopyNetFlux_dCanopyTemp   ! derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp)                       :: dCanopyNetFlux_dGroundTemp   ! derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
 real(dp)                       :: dGroundNetFlux_dCanairTemp   ! derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp)                       :: dGroundNetFlux_dCanopyTemp   ! derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp)                       :: dGroundNetFlux_dGroundTemp   ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
 real(dp)                       :: dCanopyNetFlux_dCanLiq       ! derivative in net canopy fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
 real(dp)                       :: dGroundNetFlux_dCanLiq       ! derivative in net ground fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
 ! liquid water fluxes and derivatives associated with transpiration
 real(dp)                       :: scalarCanopyTranspiration    ! canopy transpiration (kg m-2 s-1)
 real(dp)                       :: scalarCanopyEvaporation      ! canopy evaporation/condensation (kg m-2 s-1)
 real(dp)                       :: scalarGroundEvaporation      ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
 real(dp)                       :: dCanopyEvaporation_dCanLiq   ! derivative in canopy evaporation w.r.t. canopy liquid water content (s-1)
 ! energy fluxes and derivatives for the snow and soil domains
 real(dp),dimension(nLayers)    :: ssdNetNrgFlux                ! net energy flux for each layer (J m-3 s-1)
 real(dp),dimension(0:nLayers)  :: iLayerNrgFlux                ! energy flux at the layer interfaces (W m-2)
 real(dp),dimension(0:nLayers)  :: dNrgFlux_dTempAbove          ! derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
 real(dp),dimension(0:nLayers)  :: dNrgFlux_dTempBelow          ! derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
 ! liquid water fluxes and derivatives for the vegetation domain
 real(dp)                       :: canopyNetLiqFlux             ! net liquid water flux for the vegetation canopy (kg m-2 s-1)
 real(dp)                       :: scalarThroughfallRain        ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp)                       :: scalarCanopyLiqDrainage      ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
 real(dp)                       :: scalarCanopyLiqDrainageDeriv ! derivative in canopy drainage w.r.t. canopy liquid water (s-1)
 real(dp)                       :: dCanopyEvaporation_dTCanair  ! derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
 real(dp)                       :: dCanopyEvaporation_dTCanopy  ! derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
 real(dp)                       :: dCanopyEvaporation_dTGround  ! derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
 ! liquid water fluxes and derivatives for the snow domain 
 real(dp),dimension(0:nSnow)    :: iLayerLiqFluxSnow            ! vertical liquid water flux at layer interfaces (m s-1)
 real(dp),dimension(0:nSnow)    :: iLayerLiqFluxSnowDeriv       ! derivative in vertical liquid water flux at layer interfaces (m s-1)
 real(dp)                       :: scalarRainPlusMelt           ! surface water input to the soil zone (m s-1)
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
 real(dp)                       :: scalarAquiferBaseflowDeriv   ! derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
 ! ------------------------------------------------------------------------------------------------------
 ! * model solver
 ! ------------------------------------------------------------------------------------------------------
 logical(lgt),parameter         :: numericalJacobian=.false.    ! flag to compute the Jacobian matrix
 integer(i4b),parameter         :: nVegNrg=2                    ! number of energy state variables for vegetation
 integer(i4b),parameter         :: nVegLiq=1                    ! number of hydrology state variables for vegetation
 integer(i4b)                   :: nVegState                    ! number of vegetation state variables (defines position of snow-soil states in the state vector)
 integer(i4b)                   :: nState                       ! total number of model state variables
 integer(i4b),parameter         :: ixCasNrg=1                   ! index of the canopy air space state variable
 integer(i4b),parameter         :: ixVegNrg=2                   ! index of the canopy energy state variable
 integer(i4b),parameter         :: ixVegLiq=3                   ! index of the canopy liquid water state variable
 integer(i4b)                   :: ixSfcNrg                     ! index of the upper-most energy state variable in the snow-soil subdomain
 real(dp),allocatable           :: stateVecInit(:)              ! initial state vector at the start of the time step (mixed units)
 real(dp),allocatable           :: stateVecTrial(:)             ! trial state vector (mixed units)
 real(dp),allocatable           :: fluxVec0(:)                  ! flux vector (mixed units)
 real(dp),allocatable           :: fluxVec1(:)                  ! flux vector used in the numerical Jacobian calculations (mixed units)
 real(dp),allocatable           :: aJac(:,:)                    ! analytical Jacobian matrix
 real(dp),allocatable           :: nJac(:,:)                    ! numerical Jacobian matrix
 real(dp),allocatable           :: dMat(:)                      ! constant variables on the left hand side of the state equations
 real(dp),allocatable           :: rVec(:)                      ! residual vector
 real(dp),allocatable           :: xInc(:)                      ! iteration increment
 ! ------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='systemSolv_muster/'

 ! *****
 ! (0) PRELIMINARIES...
 ! ********************

 ! modify the groundwater representation for this single-column implementation
 select case(ixSpatialGroundwater)
  case(singleBasin); local_ixGroundwater = noExplicit    ! force no explicit representation of groundwater at the local scale
  case(localColumn); local_ixGroundwater = ixGroundwater ! go with the specified decision
  case default; err=20; message=trim(message)//'unable to identify spatial representation of groundwater'; return
 end select ! (modify the groundwater representation for this single-column implementation) 

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

 ! define the number of vegetation state variables (defines position of snow-soil states in the state vector)
 nVegState = nVegNrg + nVegLiq

 ! define the number of model state variables
 nState = nVegState + nSnow + nSoil

 ! define the index of the surface layer
 ixSfcNrg = nVegState + 1

 ! allocate space for the flux vector and Jacobian matrix
 allocate(stateVecInit(nState),stateVecTrial(nState),dMat(nState),fluxVec0(nState),aJac(nState,nState),rVec(nState),xInc(nState),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the flux vector and analytical Jacobian matrix'; return; endif

 ! define variables to calculate the numerical Jacobian matrix
 if(numericalJacobian)then
  ! (allocate space for the flux vector and Jacobian matrix
  allocate(fluxVec1(nState),nJac(nState,nState),stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the flux vector and numerical Jacobian matrix'; return; endif
 endif  ! if calculating the numerical approximation of the Jacobian matrix

 ! define the constant variables on the left-hand-side of the state equations
 dMat(:) = 1._dp
 dMat(ixCasNrg) = Cp_air*iden_air          ! volumetric heat capacity of air (J m-3 K-1)
 dMat(ixVegNrg) = scalarBulkVolHeatCapVeg  ! volumetric heat capacity of the vegetation (J m-3 K-1)
 dMat(ixVegLiq) = 1._dp                    ! nothing else on the left hand side
 do iLayer=1,nLayers
  jLayer = nVegState+iLayer
  dMat(jLayer) = mLayerVolHtCapBulk(iLayer)  ! volumetric heat capacity of each snow and soil layer (J m-3 K-1)
 end do

 ! initialize the analytical Jacobian matrix
 aJac(:,:) = 0._dp

 ! build the state vector for the vegetation canopy
 stateVecInit(ixCasNrg) = scalarCanairTemp
 stateVecInit(ixVegNrg) = scalarCanopyTemp
 stateVecInit(ixVegLiq) = scalarCanopyLiq

 ! build the state vector for the snow and soil domain
 do iLayer=1,nLayers
  stateVecInit(nVegState+iLayer) = mLayerTemp(iLayer)
 end do

 ! initialize the trial state vector
 stateVecTrial = stateVecInit

 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================

 ! (1) MAIN ITERATION LOOP...
 ! **************************

 ! iterate
 do iter=1,maxiter

  ! keep track of the number of iterations
  niter = iter

  ! compute model flux for a given state vector
  call computFlux(stateVecTrial,fluxVec0,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  ! compute the analytical Jacobian matrix
  call analJacob(err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)

  ! compute the residual vector
  ! NOTE: ignore phase change for now
  rVec = dMat(:)*stateVecTrial(:) - (dMat(:)*stateVecInit(:) + fluxVec0(:)*dt)

  ! compute the iteration increment
  call matrixSolv(aJac,-rVec,xInc,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  write(*,'(a,1x,i4,1x,100(e16.8,1x))') 'iter, dMat(:)*stateVecTrial(:) = ', iter, dMat(:)*stateVecTrial(:)
  write(*,'(a,1x,i4,1x,100(e16.8,1x))') 'iter, dMat(:)*stateVecInit(:)  = ', iter, dMat(:)*stateVecInit(:)
  write(*,'(a,1x,i4,1x,100(f16.8,1x))') 'iter, stateVec                 = ', iter, stateVecTrial
  write(*,'(a,1x,i4,1x,100(f16.8,1x))') 'iter, fluxVec0                 = ', iter, fluxVec0
  write(*,'(a,1x,i4,1x,100(e16.8,1x))') 'iter, xInc                     = ', iter, xInc
  write(*,'(a,1x,i4,1x,100(e16.8,1x))') 'iter, rVec                     = ', iter, rVec
  print*, '**'
  pause

  ! update the state vector
  stateVecTrial = stateVecTrial + xInc

  ! ----------------------------------------------------------------------------------------------------------------------------------------------------------------
  ! * testing: compute the numerical approximation of the Jacobian matrix
  if(numericalJacobian)then
   call numlJacob(stateVecTrial,fluxVec0,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
  endif  ! if computing the numerical Jacobian matrix

 end do  ! iterating

 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================
 ! ==========================================================================================================================================

 ! deallocate space for the flux vector and Jacobian matrix
 deallocate(stateVecInit,stateVecTrial,dMat,fluxVec0,aJac,rVec,xInc,stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to deallocate space for the flux vector and analytical Jacobian matrix'; return; endif

 ! deallocate space for the variables used to create the numerical Jacobian matrix
 if(numericalJacobian)then
  deallocate(fluxVec1,nJac,stat=err)
  if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the flux vector and numerical Jacobian matrix'; return; endif
 endif


 contains

  ! ************************************************************************************************
  ! ************************************************************************************************
  ! ************************************************************************************************
  ! *** INTERNAL SUBROUTINES ***********************************************************************
  ! ************************************************************************************************
  ! ************************************************************************************************
  ! ************************************************************************************************

  ! ************************************************************************************************
  ! internal subroutine: compute model fluxes
  ! ************************************************************************************************
  subroutine computFlux(stateVec,fluxVec,err,message)
  implicit none
  ! dummy variables
  real(dp),intent(in)            :: stateVec(:)             ! model state vector (mixed units)
  real(dp),intent(out)           :: fluxVec(:)              ! model flux vector (mixed units)
  integer(i4b),intent(out)       :: err                     ! error code
  character(*),intent(out)       :: message                 ! error message
  ! --------------------------------------------------------------
  ! initialize error control
  err=0; message='computFlux/'

  ! *****
  ! (0) PRELIMINARIES...
  ! ********************

  ! extract the vegetation states from the state vector
  scalarCanairTempTrial = stateVec(ixCasNrg)
  scalarCanopyTempTrial = stateVec(ixVegNrg)
  scalarCanopyLiqTrial  = stateVec(ixVegLiq)

  ! extract state varibles for the snow and soil domain
  do iLayer=1,nLayers
   mLayerTempTrial(iLayer) = stateVec(nVegState+iLayer)
  end do

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
                  upperBoundTemp,                         & ! intent(in): temperature of the upper boundary (K) --> NOTE: use air temperature
                  scalarCanairTempTrial,                  & ! intent(in): trial value of the canopy air space temperature (K)
                  scalarCanopyTempTrial,                  & ! intent(in): trial value of canopy temperature (K)
                  mLayerTempTrial(1),                     & ! intent(in): trial value of ground temperature (K)
                  scalarCanopyIceTrial,                   & ! intent(in): trial value of mass of ice on the vegetation canopy (kg m-2)
                  scalarCanopyLiqTrial,                   & ! intent(in): trial value of mass of liquid water on the vegetation canopy (kg m-2)
                  ! output: liquid water fluxes associated with evaporation/transpiration
                  scalarCanopyTranspiration,              & ! intent(out): canopy transpiration (kg m-2 s-1)
                  scalarCanopyEvaporation,                & ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
                  scalarGroundEvaporation,                & ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
                  ! output: fluxes
                  canairNetNrgFlux,                       & ! intent(out): net energy flux for the canopy air space (W m-2)
                  canopyNetNrgFlux,                       & ! intent(out): net energy flux for the vegetation canopy (W m-2)
                  groundNetNrgFlux,                       & ! intent(out): net energy flux for the ground surface (W m-2)
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
                  ! output: liquid water flux derivarives
                  dCanopyEvaporation_dCanLiq,             & ! intent(out): derivative in canopy evaporation w.r.t. canopy liquid water content (s-1)
                  dCanopyEvaporation_dTCanair,            & ! intent(out): derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
                  dCanopyEvaporation_dTCanopy,            & ! intent(out): derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
                  dCanopyEvaporation_dTGround,            & ! intent(out): derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
                  ! output: cross derivative terms
                  dCanopyNetFlux_dCanLiq,                 & ! intent(out): derivative in net canopy fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
                  dGroundNetFlux_dCanLiq,                 & ! intent(out): derivative in net ground fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
                  ! output: error control
                  err,cmessage)                             ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
  print*, 'canairNetNrgFlux = ', canairNetNrgFlux
  print*, 'canopyNetNrgFlux = ', canopyNetNrgFlux
  print*, 'groundNetNrgFlux = ', groundNetNrgFlux
 
  !print*, 'in systemSolv: scalarCanopyEvaporation = ', scalarCanopyEvaporation
  !print*, 'in systemSolv: dCanopyEvaporation_dCanLiq = ', dCanopyEvaporation_dCanLiq


  ! *****
  ! (2) CALCULATE ENERGY FLUXES THROUGH THE SNOW-SOIL DOMAIN...
  ! ***********************************************************
  ! calculate energy fluxes at layer interfaces through the snow and soil domain
  call ssdNrgFlux(&
                  ! input: fluxes and derivatives at the upper boundary
                  groundNetNrgFlux,                       & ! intent(in): total flux at the ground surface (W m-2)
                  dGroundNetFlux_dGroundTemp,             & ! intent(in): derivative in total ground surface flux w.r.t. ground temperature (W m-2 K-1)
                  ! input: liquid water fluxes throughout the snow and soil domains
                  iLayerLiqFluxSnow,                      & ! intent(in): liquid flux at the interface of each snow layer (m s-1)
                  iLayerLiqFluxSoil,                      & ! intent(in): liquid flux at the interface of each soil layer (m s-1)
                  ! input: trial value of model state variabes
                  mLayerTempTrial,                        & ! intent(in): trial temperature at the current iteration (K)
                  ! output: fluxes and derivatives at all layer interfaces
                  iLayerNrgFlux,                          & ! intent(out): energy flux at the layer interfaces (W m-2)
                  dNrgFlux_dTempAbove,                    & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
                  dNrgFlux_dTempBelow,                    & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
                  ! output: error control
                  err,cmessage)                             ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! calculate net energy fluxes for each snow and soil layer (J m-3 s-1)
  do iLayer=1,nLayers
   ssdNetNrgFlux(iLayer) = -(iLayerNrgFlux(iLayer) - iLayerNrgFlux(iLayer-1))/mLayerDepth(iLayer)
  end do

  ! *****
  ! (3) CALCULATE THE LIQUID FLUX THROUGH VEGETATION...
  ! ***************************************************
  call vegLiqFlux(&
                  ! input
                  computeVegFlux,                         & ! intent(in): flag to denote if computing energy flux over vegetation
                  scalarCanopyLiqTrial,                   & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                  scalarRainfall,                         & ! intent(in): rainfall rate (kg m-2 s-1)
                  ! output
                  scalarThroughfallRain,                  & ! intent(out): rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
                  scalarCanopyLiqDrainage,                & ! intent(out): drainage of liquid water from the vegetation canopy (kg m-2 s-1)
                  scalarCanopyLiqDrainageDeriv,           & ! intent(out): derivative in canopy drainage w.r.t. canopy liquid water (s-1)
                  err,cmessage)                             ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! calculate the net liquid water flux for the vegetation canopy
  canopyNetLiqFlux = scalarRainfall + scalarCanopyEvaporation - scalarThroughfallRain - scalarCanopyLiqDrainage
  !print*, 'scalarRainfall = ', scalarRainfall
  !print*, 'scalarCanopyEvaporation = ', scalarCanopyEvaporation
  !print*, 'scalarThroughfallRain = ', scalarThroughfallRain
  !print*, 'scalarCanopyLiqDrainage = ', scalarCanopyLiqDrainage

  ! *****
  ! (X) WRAP UP...
  ! **************

  ! define model flux vector for the vegetation sub-domain
  fluxVec(ixCasNrg) = canairNetNrgFlux/canopyDepth
  fluxVec(ixVegNrg) = canopyNetNrgFlux/canopyDepth
  fluxVec(ixVegLiq) = canopyNetLiqFlux

  ! define the model flux vector for the snow and soil sub-domains
  do iLayer=1,nLayers
   fluxVec(nVegState+iLayer) = ssdNetNrgFlux(iLayer)
  end do

  ! print progress
  !print*, '**'
  !write(*,'(a,1x,100(f15.9,1x))') 'stateVec(:) = ', stateVec(:)
  !write(*,'(a,1x,100(f15.9,1x))') 'fluxVec(:)  = ', fluxVec(:)

  end subroutine computFlux



  ! ************************************************************************************************
  ! internal subroutine: compute the Jacobian matrix (analytical)
  ! ************************************************************************************************
  subroutine analJacob(err,message)
  implicit none
  integer(i4b),intent(out)       :: err                     ! error code
  character(*),intent(out)       :: message                 ! error message
  ! --------------------------------------------------------------
  ! initialize error control
  err=0; message='analJacob/'

  ! liquid water fluxes for vegetation canopy (-)
  aJac(ixVegLiq,ixVegLiq) = -(dCanopyEvaporation_dCanLiq - scalarCanopyLiqDrainageDeriv)*dt + 1._dp

  ! cross-derivative terms w.r.t. system temperatures (kg m-2 K-1)
  aJac(ixVegLiq,ixCasNrg) = -dCanopyEvaporation_dTCanair*dt
  aJac(ixVegLiq,ixVegNrg) = -dCanopyEvaporation_dTCanopy*dt
  aJac(ixVegLiq,ixSfcNrg) = -dCanopyEvaporation_dTGround*dt

  ! cross-derivative terms w.r.t. canopy liquid water (J m-1 kg-1)
  aJac(ixVegNrg,ixVegLiq) = (dt/canopyDepth)   *(-dCanopyNetFlux_dCanLiq)
  aJac(ixSfcNrg,ixVegLiq) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanLiq)

  ! energy fluxes with the canopy air space (J m-3 K-1)
  aJac(ixCasNrg,ixCasNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dCanairTemp) + Cp_air*iden_air
  aJac(ixCasNrg,ixVegNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dCanopyTemp)
  aJac(ixCasNrg,ixSfcNrg) = (dt/canopyDepth)*(-dCanairNetFlux_dGroundTemp)

  ! energy fluxes with the vegetation canopy (J m-3 K-1)
  aJac(ixVegNrg,ixCasNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dCanairTemp)
  aJac(ixVegNrg,ixVegNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dCanopyTemp) + scalarBulkVolHeatCapVeg !+ dTheta_dTkCanopy*LH_fus*iden_water
  aJac(ixVegNrg,ixSfcNrg) = (dt/canopyDepth)*(-dCanopyNetFlux_dGroundTemp)

  ! energy fluxes with for the surface (J m-3 K-1)
  aJac(ixSfcNrg,ixCasNrg) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanairTemp)
  aJac(ixSfcNrg,ixVegNrg) = (dt/mLayerDepth(1))*(-dGroundNetFlux_dCanopyTemp)

  ! energy fluxes for the snow-soil system
  do iLayer=1,nLayers  ! loop through layers in the snow-soil system
   jLayer = iLayer+nVegState
   aJac(jLayer,jLayer)   = (dt/mLayerDepth(iLayer))*(-dNrgFlux_dTempBelow(iLayer-1) + dNrgFlux_dTempAbove(iLayer)) + mLayerVolHtCapBulk(iLayer) !+ mLayerdTheta_dTk(iLayer)*LH_fus*iden_water
   if(iLayer > 1)       aJac(jLayer-1,jLayer) = (dt/mLayerDepth(iLayer-1))*( dNrgFlux_dTempBelow(iLayer-1) )
   if(iLayer < nLayers) aJac(jLayer+1,jLayer) = (dt/mLayerDepth(iLayer+1))*(-dNrgFlux_dTempAbove(iLayer  ) )
  end do  ! (looping through layers in the snow-soil system)

  ! print the Jacobian
  !print*, '** analytical Jacobian:'; do iLayer=1,nState; write(*,'(i4,1x,100(e15.5,1x))') iLayer, aJac(:,iLayer); end do

  end subroutine analJacob



  ! ************************************************************************************************
  ! internal subroutine: compute the Jacobian matrix (numerical)
  ! ************************************************************************************************
  subroutine numlJacob(stateVec,fluxVec,err,message)
  implicit none
  ! dummy
  real(dp),intent(in)            :: stateVec(:)             ! model state vector (mixed units)
  real(dp),intent(in)            :: fluxVec(:)              ! model flux vector (mixed units)
  integer(i4b),intent(out)       :: err                     ! error code
  character(*),intent(out)       :: message                 ! error message
  ! local
  real(dp),dimension(nState)     :: stateVecPerturbed       ! perturbed state vector
  real(dp),dimension(nState)     :: fluxVecJac              ! flux vector
  integer(i4b)                   :: iJac                    ! index of row of the Jacobian matrix
  integer(i4b),parameter         :: iTry=3                  ! index of trial model state variable (used for testing)
  ! --------------------------------------------------------------
  ! initialize error control
  err=0; message='numlJacob/'

  ! get a copy of the state vector to perturb
  stateVecPerturbed(:) = stateVec(:)

  ! loop through state variables
  do iJac=1,nState
   ! (perturb state vector)
   stateVecPerturbed(iJac) = stateVec(iJac) + dx
   ! (compute fluxes)
   call computFlux(stateVecPerturbed,fluxVecJac,err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif  ! (check for errors)
   ! (compute the row of the Jacobian matrix)
   nJac(:,iJac) = -dt*(fluxVecJac(:) - fluxVec(:))/dx
   ! (add in the diagonal matrix)
   nJac(iJac,iJac) = nJac(iJac,iJac) + dMat(iJac)
   ! (set the state back to the input value)
   stateVecPerturbed(iJac) = stateVec(iJac)
   ! (test)
   !if(iJac==iTry) write(*,'(a,1x,i4,1x,10(f20.8,1x))'), 'iTry, -dt*(fluxVecJac(iTry) - fluxVec(iTry))/dx = ', iTry, -dt*(fluxVecJac(iTry) - fluxVec(iTry))/dx
  end do  ! (looping through state variables)

  ! print the Jacobian
  print*, '** numerical Jacobian:'; do iJac=1,nState; write(*,'(i4,1x,100(e15.5,1x))') iJac, nJac(:,iJac); end do
  pause 'testing Jacobian'

  end subroutine numlJacob








  ! ================================================================================================
  ! ================================================================================================


 end subroutine systemSolv_muster



 ! *********************************************************************************************************************************************************************
 ! *********************************************************************************************************************************************************************
 ! ***** PRIVATE SUBROUTINES *******************************************************************************************************************************************
 ! *********************************************************************************************************************************************************************
 ! *********************************************************************************************************************************************************************
 
 ! ************************************************************************************************
 ! private subroutine: compute model fluxes
 ! ************************************************************************************************
 subroutine computFluxTemp(&
                       dt,                                     & ! intent(in): time step (s)
                       iter,                                   & ! intent(in): iteration index
                       firstSubStep,                           & ! intent(in): flag to indicate if we are processing the first sub-step
                       computeVegFlux,                         & ! intent(in): flag to indicate if we need to compute fluxes over vegetation
                       local_ixGroundwater,                    & ! intent(in): index defining the local representation of groundwater 
                       ! input: domain boundary conditions
                       upperBoundTemp,                         & ! intent(in): temperature of the upper boundary of the snow and soil domains (K)
                       scalarRainfall,                         & ! intent(in): rainfall rate (kg m-2 s-1)
                       scalarSfcMeltPond,                      & ! intent(in): ponded water caused by melt of the "snow without a layer" (kg m-2)  
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
                       dCanopyNetFlux_dCanLiq,                 & ! intent(out): derivative in net canopy fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
                       dGroundNetFlux_dCanLiq,                 & ! intent(out): derivative in net ground fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
                       ! output: liquid water fluxes associated with transpiration
                       scalarCanopyTranspiration,              & ! intent(out): canopy transpiration (kg m-2 s-1)
                       scalarCanopyEvaporation,                & ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
                       scalarGroundEvaporation,                & ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
                       dCanopyEvaporation_dCanLiq,             & ! intent(out): derivative in canopy evaporation w.r.t. canopy liquid water content (s-1)
                       dCanopyEvaporation_dTCanair,            & ! intent(out): derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
                       dCanopyEvaporation_dTCanopy,            & ! intent(out): derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
                       dCanopyEvaporation_dTGround,            & ! intent(out): derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
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
                       scalarRainPlusMelt,                     & ! intent(out): surface water input to the soil zone (m s-1)
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
                       scalarAquiferBaseflowDeriv,             & ! intent(out): derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
                       ! output: error control
                       err,message)                              ! intent(out): error control
 ! provide access to the flux modules
 USE vegnrgflux_module,only:vegnrgflux      ! compute energy fluxes over vegetation
 USE ssdnrgflux_module,only:ssdnrgflux      ! compute energy fluxes throughout the snow and soil subdomains
 USE vegliqflux_module,only:vegliqflux      ! compute liquid water fluxes through vegetation
 USE snowliqflx_module,only:snowliqflx      ! compute liquid water fluxes through snow 
 USE soilliqflx_module,only:soilliqflx      ! compute liquid water fluxes through soil
 USE groundwatr_module,only:satstorage      ! compute saturated storage in the soil profile
 USE groundwatr_module,only:soilbsflow      ! compute total baseflow from the soil profile
 USE groundwatr_module,only:disaggflow      ! diaggregate total baseflow to individual soil layers
 USE groundwatr_module,only:aquifrflux      ! compute liquid water fluxes in the aquifer
 implicit none
 ! input: model control
 real(dp),intent(in)            :: dt                            ! time step (s)
 integer(i4b),intent(in)        :: iter                          ! iteration index
 logical(lgt),intent(in)        :: firstSubStep                  ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(in)        :: computeVegFlux                ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 integer(i4b),intent(in)        :: local_ixGroundwater           ! choice of groundwater parameterization
 ! input: domain boundary conditions
 real(dp),intent(in)            :: upperBoundTemp                ! temperature of the upper boundary of the snow and soil domains (K)
 real(dp),intent(in)            :: scalarRainfall                ! rainfall (kg m-2 s-1)
 real(dp),intent(in)            :: scalarSfcMeltPond             ! ponded water caused by melt of the "snow without a layer" (kg m-2)
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
 real(dp),intent(out)           :: dCanopyNetFlux_dCanLiq        ! derivative in net canopy fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
 real(dp),intent(out)           :: dGroundNetFlux_dCanLiq        ! derivative in net ground fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
 ! output: liquid water fluxes associated with transpiration
 real(dp),intent(out)           :: scalarCanopyTranspiration     ! canopy transpiration (kg m-2 s-1)
 real(dp),intent(out)           :: scalarCanopyEvaporation       ! canopy evaporation/condensation (kg m-2 s-1)
 real(dp),intent(out)           :: scalarGroundEvaporation       ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
 real(dp),intent(out)           :: dCanopyEvaporation_dCanLiq    ! derivative in canopy evaporation w.r.t. canopy liquid water content (s-1)
 real(dp),intent(out)           :: dCanopyEvaporation_dTCanair   ! derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
 real(dp),intent(out)           :: dCanopyEvaporation_dTCanopy   ! derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
 real(dp),intent(out)           :: dCanopyEvaporation_dTGround   ! derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
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
 real(dp),intent(out)           :: scalarRainPlusMelt            ! surface water input to the soil zone (m s-1)
 ! output: liquid water fluxes and derivatives for the soil domain
 real(dp),intent(inout)         :: iLayerLiqFluxSoil(0:)         ! liquid flux at soil layer interfaces (m s-1)
 real(dp),intent(out)           :: mLayerTranspire(:)            ! transpiration loss from each soil layer (m s-1)
 real(dp),intent(inout)         :: mLayerBaseflow(:)             ! baseflow from each soil layer -- only compute at the start of the step (m s-1)
 real(dp),intent(out)           :: dq_dStateAbove(0:)            ! change in the flux in layer interfaces w.r.t. state variables in the layer above
 real(dp),intent(out)           :: dq_dStateBelow(0:)            ! change in the flux in layer interfaces w.r.t. state variables in the layer below
 ! output: liquid water fluxes and derivatives for the aquifer storage
 real(dp),intent(out)           :: scalarAquiferTranspire        ! transpiration loss from the aquifer at the start-of-step (m s-1)
 real(dp),intent(out)           :: scalarAquiferRecharge         ! recharge to the aquifer (m s-1)
 real(dp),intent(out)           :: scalarAquiferBaseflow         ! total baseflow from the aquifer (m s-1)
 real(dp),intent(out)           :: scalarAquiferBaseflowDeriv    ! derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
 ! output: error control
 integer(i4b),intent(out)       :: err                           ! error code
 character(*),intent(out)       :: message                       ! error message
 ! ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! general local variables
 character(LEN=256)             :: cmessage                 ! error message of downwind routine
 ! local variables for the lateral flux among soil columns
 integer(i4b)                   :: ixSaturation                  ! index of the lowest saturated layer
 real(dp)                       :: subSurfaceStorage             ! sub surface storage (m)
 real(dp)                       :: maximumSoilWater              ! maximum storage (m)
 real(dp)                       :: maximumFlowRate               ! flow rate under saturated conditions (m/s)
 real(dp)                       :: totalColumnInflow             ! total column inflow (m/s)
 real(dp)                       :: totalColumnOutflow            ! total outflow from the soil column (m s-1)
 real(dp)                       :: totalOutflowDeriv             ! derivative in total outflow w.r.t. storage (s-1)
 real(dp),dimension(nLayers)    :: mLayerHydCond                 ! hydraulic conductivity in each soil layer (m s-1)
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
                 upperBoundTemp,                         & ! intent(in): temperature of the upper boundary (K) --> NOTE: use air temperature
                 scalarCanairTempTrial,                  & ! intent(in): trial value of the canopy air space temperature (K)
                 scalarCanopyTempTrial,                  & ! intent(in): trial value of canopy temperature (K)
                 mLayerTempTrial(1),                     & ! intent(in): trial value of ground temperature (K)
                 scalarCanopyIceTrial,                   & ! intent(in): trial value of mass of ice on the vegetation canopy (kg m-2)
                 scalarCanopyLiqTrial,                   & ! intent(in): trial value of mass of liquid water on the vegetation canopy (kg m-2)
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
                 ! output: liquid water flux derivarives
                 dCanopyEvaporation_dCanLiq,             & ! intent(out): derivative in canopy evaporation w.r.t. canopy liquid water content (s-1)
                 dCanopyEvaporation_dTCanair,            & ! intent(out): derivative in canopy evaporation w.r.t. canopy air temperature (kg m-2 s-1 K-1)
                 dCanopyEvaporation_dTCanopy,            & ! intent(out): derivative in canopy evaporation w.r.t. canopy temperature (kg m-2 s-1 K-1)
                 dCanopyEvaporation_dTGround,            & ! intent(out): derivative in canopy evaporation w.r.t. ground temperature (kg m-2 s-1 K-1)
                 ! output: cross derivative terms
                 dCanopyNetFlux_dCanLiq,                 & ! intent(out): derivative in net canopy fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
                 dGroundNetFlux_dCanLiq,                 & ! intent(out): derivative in net ground fluxes w.r.t. canopy liquid water content (J kg-1 s-1)
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
                 dNrgFlux_dTempAbove,                    & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
                 dNrgFlux_dTempBelow,                    & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
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
                 scalarRainfall,                         & ! intent(in): rainfall rate (kg m-2 s-1)
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
  ! compute liquid fluxes
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
  ! define forcing for the soil domain
  scalarRainPlusMelt = iLayerLiqFluxSnow(nSnow)    ! drainage from the base of the snowpack
 else
  ! define forcing for the soil domain
  scalarRainPlusMelt = (scalarThroughfallRain + scalarCanopyLiqDrainage)/iden_water + &  ! liquid flux from the canopy (m s-1)
                        + (scalarSfcMeltPond/dt)/iden_water  ! melt of the snow without a layer (m s-1)
 endif

 ! *****
 ! (5) CALCULATE THE LIQUID FLUX THROUGH SOIL...
 ! *********************************************
 call soilLiqFlx(&
                 ! input: model control
                 iter,                                   & ! intent(in): iteration index
                 .true.,                                 & ! intent(in): flag indicating if derivatives are desired
                 ! input: trial state variables
                 mLayerMatricHeadTrial,                  & ! intent(in): matric head (m)
                 mLayerVolFracLiqTrial(nSnow+1:nLayers), & ! intent(in): volumetric fraction of liquid water (-)
                 mLayerVolFracIceTrial(nSnow+1:nLayers), & ! intent(in): volumetric fraction of ice (-)
                 ! input: fluxes
                 scalarCanopyTranspiration,              & ! intent(in): canopy transpiration (kg m-2 s-1)
                 scalarGroundEvaporation,                & ! intent(in): ground evaporation (kg m-2 s-1)
                 scalarRainPlusMelt,                     & ! intent(in): rain plus melt (m s-1)
                 ! output: fluxes
                 iLayerLiqFluxSoil,                      & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                 mLayerTranspire,                        & ! intent(out): transpiration loss from each soil layer (m s-1)
                 mLayerHydCond,                          & ! intent(out): hydraulic conductivity in each layer (m s-1)
                 ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                 dq_dStateAbove,                         & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer above (m s-1)
                 dq_dStateBelow,                         & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer below (m s-1)
                 ! output: error control
                 err,cmessage)                             ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! *****
 ! (5) CALCULATE THE LATERAL FLUX IN/OUT OF THE SOIL PROFILE...
 ! ************************************************************

 ! check if the option for lateral flux is invoked
 if(local_ixGroundwater==qbaseTopmodel)then

  ! * compute saturated storage at the bottom of the soil profile (m)
  call satStorage(&
                  ! input
                  mLayerVolFracLiqTrial(nSnow+1:nLayers), & ! intent(in): volumetric fraction of liquid water (-)
                  mLayerVolFracIceTrial(nSnow+1:nLayers), & ! intent(in): volumetric fraction of ice (-)
                  ! output
                  ixSaturation,                           & ! intent(out): index of the lowest saturated layer
                  subSurfaceStorage,                      & ! intent(out): sub surface storage (m)
                  maximumSoilWater,                       & ! intent(out): maximum storage (m)
                  err,cmessage)                             ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! * compute the baseflow flux and the derivative
  call soilBsFlow(&
                  ! input: model control
                  iter,                                   & ! intent(in): iteration index
                  ! input: storage
                  subSurfaceStorage,                      & ! intent(in): sub surface storage (m)
                  maximumSoilWater,                       & ! intent(in): maximum possible sub surface storage (m)
                  ! input/output: diagnostic variables and fluxes constant over iterations
                  maximumFlowRate,                        & ! intent(inout): flow rate under saturated conditions (m/s)
                  totalColumnInflow,                      & ! intent(inout): total column inflow (m/s)
                  ! output: outflow and its derivative w.r.t. storage
                  totalColumnOutflow,                     & ! intent(out): total outflow from the soil column (m s-1)
                  totalOutflowDeriv,                      & ! intent(out): derivative in total outflow w.r.t. storage (s-1)
                  ! output: error control
                  err,cmessage)                             ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! * disaggregate total inflow and total outflow to the baseflow sink term
  call disaggFlow(&
                  ! input
                  dt,                                     & ! intent(in): time step (s) -- used to calculate maximum possible inflow rate
                  ixSaturation,                           & ! intent(in): index of the lowest saturated layer
                  totalColumnInflow,                      & ! intent(in): total column inflow (m s-1)
                  totalColumnOutflow,                     & ! intent(in): total outflow from the soil column (m s-1)
                  mLayerVolFracLiqTrial(nSnow+1:nLayers), & ! intent(in): volumetric fraction of liquid water (-)
                  mLayerVolFracIceTrial(nSnow+1:nLayers), & ! intent(in): volumetric fraction of ice (-)
                  mLayerHydCond,                          & ! intent(in): hydraulic conductivity in each soil layer (m s-1)
                  ! output
                  mLayerBaseflow,                         & ! intent(out): baseflow from each soil layer (m s-1)
                  err,cmessage)                             ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! not computing lateral flux
 else

  ! source/sink term is zero
  mLayerBaseflow(:) = 0._dp

 endif

 ! *****
 ! (6) CALCULATE THE FLUXES FOR THE AQUIFER...
 ! *******************************************
 call aquifrFlux(&
                 ! input
                 local_ixGroundwater,         & ! intent(in): index defining the local representation of groundwater
                 scalarAquiferStorageTrial,   & ! intent(in): aquifer storage (m)
                 scalarCanopyTranspiration,   & ! intent(in): canopy transpiration (kg m-2 s-1)
                 iLayerLiqFluxSoil(nSoil),    & ! intent(in): drainage from the bottom of the soil profile
                 ! output
                 scalarAquiferTranspire,      & ! intent(out): transpiration loss from the aquifer (m s-1)
                 scalarAquiferRecharge,       & ! intent(out): recharge to the aquifer (m s-1)
                 scalarAquiferBaseflow,       & ! intent(out): total baseflow from the aquifer (m s-1)
                 scalarAquiferBaseflowDeriv,  & ! intent(out): derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
                 err,cmessage)                  ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif



 ! ******************************************************************************************************************************************
 ! ******************************************************************************************************************************************

 end subroutine computFluxTemp





end module systemSolv_module

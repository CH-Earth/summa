module soilHydrol_module
! -----------------------------------------------------------------------------------------------------------
USE nrtype
! provide access to model constants
USE multiconst,only:iden_ice,iden_water            ! intrinsic density of ice and water (kg m-3)
! provide access to layer types
USE data_struc,only:ix_soil,ix_snow                ! named variables for snow and soil
! provide access to look-up values for model decisions
USE mDecisions_module,only:  &
 ! look-up values for method used to compute derivative
 numerical,                  & ! numerical solution
 analytical,                 & ! analytical solution
 ! look-up values for the form of Richards' equation
 moisture,                   & ! moisture-based form of Richards' equation
 mixdform,                   & ! mixed form of Richards' equation
 ! look-up values for the type of hydraulic conductivity profile
 constant,                   & ! constant hydraulic conductivity with depth
 exp_profile,                & ! exponential profile
 powerLaw_profile,           & ! power-law profile
 linear_profile,             & ! linear profile
 ! look-up values for the choice of groundwater parameterization
 movingBoundary,             & ! moving lower boundary
 bigBucket,                  & ! a big bucket (lumped aquifer model)
 noExplicit,                 & ! no explicit groundwater parameterization
 ! look-up values for the choice of boundary conditions for hydrology
 prescribedHead,             & ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
 funcBottomHead,             & ! function of matric head in the lower-most layer
 freeDrainage,               & ! free drainage
 liquidFlux,                 & ! liquid water flux
 groundwaterCouple             ! coupled to the groundwater sub-model (matric head=0 as a moving lower boundary)
! -----------------------------------------------------------------------------------------------------------
implicit none
private
public::soilHydrol
public::waterTablePosition
! number of layers
integer(i4b)           :: nSoil                    ! number of soil layers
integer(i4b)           :: nSnow                    ! number of snow layers
integer(i4b)           :: nLevels                  ! total number of soil layers to examine
! missing value parameter
real(dp),parameter     :: valueMissing=-1.e+20_dp  ! missing value parameter
! finite difference increment
real(dp),parameter     :: dx=1.e-8_dp              ! finite difference increment
contains


 ! ************************************************************************************************
 ! new subroutine: compute change in volumetric liquid water content (or matric head) over the time step
 ! ************************************************************************************************
 subroutine soilHydrol(&
                       ! input
                       dt,                      & ! time step (seconds)
                       iter,                    & ! iteration index
                       mLayerMatricHeadIter,    & ! matric head in each layer at the current iteration (m)
                       mLayerVolFracLiqIter,    & ! volumetric fraction of liquid water at the current iteration (-)
                       mLayerVolFracIceIter,    & ! volumetric fraction of ice at the current iteration (-)
                       scalarAquiferStorageIter,& ! aquifer storage (m)
                       ! output
                       mLayerMatricHeadNew,     & ! matric head in each layer at the next iteration (m)
                       mLayerVolFracLiqNew,     & ! volumetric fraction of liquid water at the next iteration (-)
                       scalarAquiferStorageNew, & ! aquifer storage (m)
                       err,message)
 ! model decision structures
 USE data_struc,only:model_decisions        ! model decision structure
 USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:mpar_data,forc_data,mvar_data,indx_data    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX ! named variables for structure elements
 ! compute change in temperature over the time step
 implicit none
 ! input
 real(dp),intent(in)           :: dt                       ! time step (seconds)
 integer(i4b),intent(in)       :: iter                     ! iteration index
 real(dp),intent(in)           :: mLayerMatricHeadIter(:)  ! matric head in each layer at the current iteration (m)
 real(dp),intent(in)           :: mLayerVolFracLiqIter(:)  ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),intent(in)           :: mLayerVolFracIceIter(:)  ! volumetric fraction of ice at the current iteration (-)
 real(dp),intent(in)           :: scalarAquiferStorageIter ! aquifer storage (m)
 ! output
 real(dp),intent(out)          :: mLayerMatricHeadNew(:)   ! matric head in each layer at the next iteration (m)
 real(dp),intent(out)          :: mLayerVolFracLiqNew(:)   ! volumetric fraction of liquid water at the next iteration (-)
 real(dp),intent(out)          :: scalarAquiferStorageNew  ! aquifer storage (m)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! internal
 integer(i4b)                  :: ibeg,iend                ! start and end indices of the soil layers in concatanated snow-soil vector
 character(LEN=256)            :: cmessage                 ! error message of downwind routine
 real(dp),pointer              :: waterTableDepth          ! depth to the water table (m)
 real(dp),pointer              :: volFracLiqPartiallySat   ! volumetric fraction of liquid water in the unsaturated portion of the partially saturated layer (-)
 real(dp),pointer              :: iUpperPartiallySat       ! height at the upper interface of the layer that contains the water table (m)
 real(dp),pointer              :: iLowerPartiallySat       ! height at the lower interface of the layer that contains the water table (m)
 real(dp),pointer              :: mDepthPartiallySat       ! depth of the layer that contains the water table (m)
 real(dp)                      :: satFrac                  ! fraction of partially saturated layer that is saturated (-)
 ! initialize error control
 err=0; message='soilHydrol/'

 ! identify the number of layers
 nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat == ix_soil)  ! number of soil layers
 nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat == ix_snow)  ! number of snow layers

 ! get the number of levels in the "effective" soil profile
 nLevels = size(mLayerMatricHeadIter)

 ! get indices for the data structures
 ibeg = nSnow + 1
 iend = nSnow + nLevels

 ! *****
 ! wrapper for the soil hydrology sub-routine...
 ! *********************************************
 ! used as a trick to both
 !  1) save typing; and
 !  2) "protect" variables through the use of the intent attribute
 call soilHydrol_muster(&
                        ! input variables from soilHydrol routine -- NOTE: inputs are already sized appropriately
                        dt,                                                           & ! intent(in): time step (seconds)
                        iter,                                                         & ! intent(in): current iteration count
                        mLayerMatricHeadIter,                                         & ! intent(in): matric head in each layer at the current iteration (m)
                        mLayerVolFracLiqIter,                                         & ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                        mLayerVolFracIceIter,                                         & ! intent(in): volumetric fraction of ice at the current iteration (-)
                        scalarAquiferStorageIter,                                     & ! intent(in): aquifer storage (m)
                        ! named variables for model decisions
                        model_decisions(iLookDECISIONS%fDerivMeth)%iDecision,         & ! intent(in): method used to calculate flux derivatives 
                        model_decisions(iLookDECISIONS%f_Richards)%iDecision,         & ! intent(in): form of Richards' equation
                        model_decisions(iLookDECISIONS%hc_Profile)%iDecision,         & ! intent(in): option for the hydraulic conductivity profile
                        model_decisions(iLookDECISIONS%groundwatr)%iDecision,         & ! intent(in): groundwater parameterization
                        model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision,         & ! intent(in): upper boundary conditions for soil hydrology
                        model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision,         & ! intent(in): lower boundary conditions for soil hydrology
                        ! model coordinate variables -- NOTE: use of ibeg and iend 
                        mvar_data%var(iLookMVAR%mLayerDepth)%dat(ibeg:iend),          & ! intent(in): depth of the layer (m)
                        mvar_data%var(iLookMVAR%mLayerHeight)%dat(ibeg:iend),         & ! intent(in): height of the layer mid-point (m)
                        mvar_data%var(iLookMVAR%iLayerHeight)%dat(ibeg-1:iend),       & ! intent(in): height of the layer interfaces (m)
                        ! boundary conditions for matric head
                        mpar_data%var(iLookPARAM%lowerBoundHead),                     & ! intent(in): lower boundary condition (m)
                        mpar_data%var(iLookPARAM%upperBoundHead),                     & ! intent(in): upper boundary condition (m)
                        ! boundary conditions for volumetric liquid water content
                        mpar_data%var(iLookPARAM%lowerBoundTheta),                    & ! intent(in): lower boundary condition (-)
                        mpar_data%var(iLookPARAM%upperBoundTheta),                    & ! intent(in): upper boundary condition (-)
                        ! model forcing
                        mvar_data%var(iLookMVAR%scalarRainfall)%dat(1),               & ! intent(in): computed rainfall rate (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%iLayerLiqFluxSnow)%dat(nSnow),        & ! intent(in): liquid flux from the base of the snowpack (m s-1)
                        ! general model parameters
                        mpar_data%var(iLookPARAM%wimplicit),                          & ! intent(in): weight assigned to start-of-step fluxes (-)
                        ! soil parameters
                        mpar_data%var(iLookPARAM%vGn_alpha),                          & ! intent(in): van Genutchen "alpha" parameter (m-1)
                        mpar_data%var(iLookPARAM%vGn_n),                              & ! intent(in): van Genutchen "n" parameter (-)
                        mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1),                  & ! intent(in): van Genutchen "m" parameter (-)
                        mpar_data%var(iLookPARAM%theta_sat),                          & ! intent(in): soil porosity (-)
                        mpar_data%var(iLookPARAM%theta_res),                          & ! intent(in): soil residual volumetric water content (-)
                        mpar_data%var(iLookPARAM%kAnisotropic),                       & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                        mpar_data%var(iLookPARAM%zScale_TOPMODEL),                    & ! intent(in): scale factor for TOPMODEL-ish baseflow parameterization (m)
                        mpar_data%var(iLookPARAM%bpar_VIC),                           & ! intent(in): b-parameter in the VIC surface runoff parameterization (-)
                        mpar_data%var(iLookPARAM%specificYield),                      & ! intent(in): specific yield (-)
                        mpar_data%var(iLookPARAM%specificStorage),                    & ! intent(in): specific storage coefficient (m-1)
                        mpar_data%var(iLookPARAM%f_impede),                           & ! intent(in): ice impedence factor (-)
                        ! model state variables -- NOTE: use of ibeg and iend
                        mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(ibeg:iend),     & ! intent(in): volumetric fraction of ice in each layer (-)
                        mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(ibeg:iend),     & ! intent(in): volumetric fraction of liquid water in each layer (-)
                        mvar_data%var(iLookMVAR%mLayerMatricHead)%dat(1:nLevels),     & ! intent(in): matric head in each layer (m)
                        mvar_data%var(iLookMVAR%scalarAquiferStorage)%dat(1),         & ! intent(in): (scalar) aquifer storage at the start of the time step (m)
                        mvar_data%var(iLookMVAR%scalarSfcMeltPond)%dat(1),            & ! intent(in): (scalar) ponded water caused by melt of the "snow without a layer" (kg m-2)
                        ! saturated hydraulic conductivity
                        mvar_data%var(iLookMVAR%mLayerSatHydCond)%dat,                & ! intent(in): saturated hydraulic conductivity at the mid-point of each layer (m s-1)
                        mvar_data%var(iLookMVAR%iLayerSatHydCond)%dat,                & ! intent(in): saturated hydraulic conductivity at the interface of each layer (m s-1)
                        ! transpiration (from energy routines)
                        mvar_data%var(iLookMVAR%mLayerInitTranspire)%dat,             & ! intent(in): transpiration loss from each soil layer at start-of-step (m s-1)
                        mvar_data%var(iLookMVAR%mLayerTranspire)%dat,                 & ! intent(in): transpiration loss from each soil layer (m s-1)
                        mvar_data%var(iLookMVAR%scalarInitAquiferTranspire)%dat(1),   & ! intent(out): transpiration loss from the aquifer at the start-of-step (m s-1)
                        mvar_data%var(iLookMVAR%scalarAquiferTranspire)%dat(1),       & ! intent(out): transpiration loss from the aquifer at the end-of-step (m s-1)
                        ! diagnostic scalar variables
                        mvar_data%var(iLookMVAR%scalarRainPlusMelt)%dat(1),           & ! intent(out): (scalar) rain plus melt (m s-1)
                        mvar_data%var(iLookMVAR%scalarSurfaceRunoff)%dat(1),          & ! intent(out): (scalar) surface runoff (m s-1)
                        mvar_data%var(iLookMVAR%scalarWaterTableDepth)%dat(1),        & ! intent(out): (scalar) water table depth (m)
                        mvar_data%var(iLookMVAR%scalarInitAquiferRcharge)%dat(1),     & ! intent(out): (scalar) recharge to the aquifer at the start-of-step (m s-1)
                        mvar_data%var(iLookMVAR%scalarAquiferRcharge)%dat(1),         & ! intent(out): (scalar) recharge to the aquifer at the end-of-step (m s-1)
                        mvar_data%var(iLookMVAR%scalarInitAquiferBaseflow)%dat(1),    & ! intent(out): (scalar) baseflow from the aquifer at the start-of-step (m s-1)
                        mvar_data%var(iLookMVAR%scalarAquiferBaseflow)%dat(1),        & ! intent(out): (scalar) baseflow from the aquifer at the end-of-step (m s-1)
                        ! model diagnostic variables
                        mvar_data%var(iLookMVAR%mLayerdTheta_dPsi)%dat(1:nLevels),    & ! intent(out): derivative in the soil water characteristic w.r.t. psi (m-1)
                        mvar_data%var(iLookMVAR%mLayerdPsi_dTheta)%dat(1:nLevels),    & ! intent(out): derivative in the soil water characteristic w.r.t. theta (m)
                        mvar_data%var(iLookMVAR%iLayerInitLiqFluxSoil)%dat(0:nLevels),& ! intent(out): liquid flux at layer interfaces at the start of the time step (m s-1)
                        mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat(0:nLevels),    & ! intent(out): liquid flux at layer interfaces at the end of the time step (m s-1)
                        mvar_data%var(iLookMVAR%mLayerInitEjectWater)%dat(1:nLevels), & ! intent(out): water ejected from each soil layer at the start-of-step (m s-1)
                        mvar_data%var(iLookMVAR%mLayerEjectWater)%dat(1:nLevels),     & ! intent(out): water ejected from each soil layer (m s-1)
                        ! output variables from the soilHydrol routine  -- NOTE: variables are already sized appropriately
                        mLayerMatricHeadNew,                                          & ! intent(out): matric head in each layer at the next iteration (m)
                        mLayerVolFracLiqNew,                                          & ! intent(out): volumetric fraction of liquid water at the next iteration (-)
                        scalarAquiferStorageNew,                                      & ! intent(out): aquifer storage (m)
                        err,cmessage)                                                   ! intent(out): error control
 ! check for errors
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! set un-used portion of diagnostic variables to zero
 if(nLevels < nSoil)then
  mvar_data%var(iLookMVAR%mLayerdTheta_dPsi)%dat(    nLevels+1:nSoil) = 0._dp
  mvar_data%var(iLookMVAR%mLayerdPsi_dTheta)%dat(    nLevels+1:nSoil) = 0._dp
  mvar_data%var(iLookMVAR%iLayerInitLiqFluxSoil)%dat(nLevels+1:nSoil) = 0._dp
  mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat(    nLevels+1:nSoil) = 0._dp
  mvar_data%var(iLookMVAR%mLayerInitEjectWater)%dat( nLevels+1:nSoil) = 0._dp
  mvar_data%var(iLookMVAR%mLayerEjectWater)%dat(     nLevels+1:nSoil) = 0._dp
 endif

 end subroutine soilHydrol


 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ***** PRIVATE SUBROUTINES **********************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************

 ! ************************************************************************************************
 ! private subroutine: wrapper for the soil hydrology subroutine
 ! ************************************************************************************************
 subroutine soilHydrol_muster(&
                              ! input variables from soilHydrol routine
                              dt,                         & ! intent(in): time step (seconds)
                              iter,                       & ! intent(in): current iteration count
                              mLayerMatricHeadIter,       & ! intent(in): matric head in each layer at the current iteration (m)
                              mLayerVolFracLiqIter,       & ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                              mLayerVolFracIceIter,       & ! intent(in): volumetric fraction of ice at the current iteration (-)
                              scalarAquiferStorageIter,   & ! intent(in): aquifer storage (m)
                              ! model decisions
                              ixDerivMethod,              & ! intent(in): choice of method used to compute derivative
                              ixRichards,                 & ! intent(in): choice of the form of Richards' equation
                              hc_Profile,                 & ! intent(in): index defining the option for the hydraulic conductivity profile
                              ixGroundwater,              & ! intent(in): choice of groundwater parameterization
                              ixBcUpperSoilHydrology,     & ! intent(in): choice of upper boundary condition for soil hydrology
                              ixBcLowerSoilHydrology,     & ! intent(in): choice of upper boundary condition for soil hydrology
                              ! model coordinate variables -- NOTE: use of ibeg and iend 
                              mLayerDepth,                & ! intent(in): depth of the layer (m)
                              mLayerHeight,               & ! intent(in): height of the layer mid-point (m)
                              iLayerHeight,               & ! intent(in): height of the layer interfaces (m)
                              ! boundary conditions for matric head
                              lowerBoundHead,             & ! intent(in): lower boundary condition (m)
                              upperBoundHead,             & ! intent(in): upper boundary condition (m)
                              ! boundary conditions for volumetric liqquid water content
                              lowerBoundTheta,            & ! intent(in): lower boundary condition (-)
                              upperBoundTheta,            & ! intent(in): upper boundary condition (-)
                              ! model forcing
                              scalarRainfall,             & ! intent(in): computed rainfall rate (kg m-2 s-1)
                              scalarLiqFluxSnow,          & ! intent(in): liquid flux from the base of the snowpack (m s-1)
                              ! general model parameters
                              wimplicit,                  & ! intent(in): weight assigned to start-of-step fluxes (-)
                              ! soil parameters
                              vGn_alpha,                  & ! intent(in): van Genutchen "alpha" parameter (m-1)
                              vGn_n,                      & ! intent(in): van Genutchen "n" parameter (-)
                              VGn_m,                      & ! intent(in): van Genutchen "m" parameter (-)
                              theta_sat,                  & ! intent(in): soil porosity (-)
                              theta_res,                  & ! intent(in): soil residual volumetric water content (-)
                              kAnisotropic,               & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                              zScale_TOPMODEL,            & ! intent(in): scale factor for TOPMODEL-ish baseflow parameterization (m)
                              bpar_VIC,                   & ! intent(in): b-parameter in the VIC surface runoff parameterization (-)
                              specificYield,              & ! intent(in): specific yield (-)
                              specificStorage,            & ! intent(in): specific storage coefficient (m-1)
                              f_impede,                   & ! intent(in): ice impedence factor (-)
                              ! model state variables -- NOTE: use of ibeg and iend
                              mLayerVolFracIce,           & ! intent(in): volumetric fraction of ice in each layer (-)
                              mLayerVolFracLiq,           & ! intent(in): volumetric fraction of liquid water in each layer (-)
                              mLayerMatricHead,           & ! intent(in): matric head in each layer (m)
                              scalarAquiferStorage,       & ! intent(in): (scalar)    ! aquifer storage at the start of the time step (m)
                              scalarSfcMeltPond,          & ! intent(in): (scalar)    ! ponded water caused by melt of the "snow without a layer" (kg m-2)
                              ! saturated hydraulic conductivity in each layer
                              mLayerSatHydCond,           & ! intent(in): saturated hydraulic conductivity at the mid-point of each layer (m s-1)
                              iLayerSatHydCond,           & ! intent(in): saturated hydraulic conductivity at the interface of each layer (m s-1)
                              ! transpiration (from energy routines)
                              mLayerInitTranspire,        & ! intent(in): transpiration loss from each soil layer at start-of-step (m s-1)
                              mLayerTranspire,            & ! intent(in): transpiration loss from each soil layer (m s-1)
                              scalarInitAquiferTranspire, & ! intent(in): transpiration loss from the aquifer at the start-of-step (m s-1)
                              scalarAquiferTranspire,     & ! intent(in): transpiration loss from the aquifer at the end-of-step (m s-1)
                              ! diagnostic scalar variables
                              scalarRainPlusMelt,         & ! intent(out): rain plus melt (m s-1)
                              scalarSurfaceRunoff,        & ! intent(out): surface runoff (m s-1)
                              scalarWaterTableDepth,      & ! intent(out): water table depth (m)
                              scalarInitAquiferRcharge,   & ! intent(out): recharge to the aquifer at the start-of-step (m s-1)
                              scalarAquiferRcharge,       & ! intent(out): recharge to the aquifer at the end-of-step (m s-1)
                              scalarInitAquiferBaseflow,  & ! intent(out): baseflow from the aquifer at the start-of-step (m s-1)
                              scalarAquiferBaseflow,      & ! intent(out): baseflow from the aquifer at the end-of-step (m s-1)
                              ! model diagnostic variables
                              mLayerdTheta_dPsi,          & ! intent(out): derivative in the soil water characteristic w.r.t. psi (m-1)
                              mLayerdPsi_dTheta,          & ! intent(out): derivative in the soil water characteristic w.r.t. theta (m)
                              iLayerInitLiqFluxSoil,      & ! intent(out): liquid flux at layer interfaces at the start of the time step (m s-1)
                              iLayerLiqFluxSoil,          & ! intent(out): liquid flux at layer interfaces at the end of the time step (m s-1)
                              mLayerInitEjectWater,       & ! intent(out): water ejected from each soil layer at the start-of-step (m s-1)
                              mLayerEjectWater,           & ! intent(out): water ejected from each soil layer (m s-1)
                              ! output variables from the soilHydrol routine
                              mLayerMatricHeadNew,        & ! intent(out): matric head in each layer at the next iteration (m)
                              mLayerVolFracLiqNew,        & ! intent(out): volumetric fraction of liquid water at the next iteration (-)
                              scalarAquiferStorageNew,    & ! intent(out): aquifer storage (m)
                              err,message)                  ! intent(out): error control
 ! utility modules
 USE soil_utils_module,only:volFracLiq      ! compute volumetric fraction of liquid water
 USE soil_utils_module,only:matricHead      ! compute matric head (m)
 USE soil_utils_module,only:dTheta_dPsi     ! compute derivative of the soil moisture characteristic w.r.t. psi (m-1)
 USE soil_utils_module,only:dPsi_dTheta     ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
 USE tridagSolv_module,only:tridag          ! solve tridiagonal system of equations
 implicit none
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! ***** input variables
 ! ------------------------------------------------------------------------------------------------------------------------------------------------- 
 ! input variables from the soilHydrol routine
 real(dp),intent(in)              :: dt                       ! time step (seconds)
 integer(i4b),intent(in)          :: iter                     ! iteration index
 real(dp),intent(in)              :: mLayerMatricHeadIter(:)  ! matric head in each layer at the current iteration (m)
 real(dp),intent(in)              :: mLayerVolFracLiqIter(:)  ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),intent(in)              :: mLayerVolFracIceIter(:)  ! volumetric fraction of ice at the current iteration (-)
 real(dp),intent(in)              :: scalarAquiferStorageIter ! aquifer storage (m)
 ! model decisions
 integer(i4b),intent(in)          :: ixDerivMethod            ! choice of method used to compute derivative
 integer(i4b),intent(in)          :: ixRichards               ! choice of the form of Richards' equation
 integer(i4b),intent(in)          :: hc_profile               ! choice of type of hydraulic conductivity profile
 integer(i4b),intent(in)          :: ixGroundwater            ! choice of groundwater parameterization
 integer(i4b),intent(in)          :: ixBcUpperSoilHydrology   ! choice of upper boundary condition for soil hydrology
 integer(i4b),intent(in)          :: ixBcLowerSoilHydrology   ! choice of upper boundary condition for soil hydrology
 ! model coordinate variables
 real(dp),intent(in)              :: mLayerDepth(:)           ! depth of the layer (m)
 real(dp),intent(in)              :: mLayerHeight(:)          ! height of the layer mid-point (m)
 real(dp),intent(in)              :: iLayerHeight(0:)         ! height of the layer interfaces (m)
 ! diriclet boundary conditions for matric head
 real(dp),intent(in)              :: upperBoundHead           ! upper boundary condition for matric head (m)
 real(dp),intent(in)              :: lowerBoundHead           ! lower boundary condition for matric head (m)
 ! diriclet boundary conditions for volumetric liquid water content
 real(dp),intent(in)              :: upperBoundTheta          ! upper boundary condition for volumetric liquid water content (-)
 real(dp),intent(in)              :: lowerBoundTheta          ! lower boundary condition for volumetric liquid water content (-)
 ! model forcing
 real(dp),intent(in)              :: scalarRainfall           ! computed rainfall rate (kg m-2 s-1)
 real(dp),intent(in)              :: scalarLiqFluxSnow        ! liquid flux from the base of the snowpack (m s-1)
 ! general model parameters
 real(dp),intent(in)              :: wimplicit                ! weight assigned to start-of-step fluxes (-)
 ! soil parameters
 real(dp),intent(in)              :: vGn_alpha                ! van Genutchen "alpha" parameter (m-1)
 real(dp),intent(in)              :: vGn_n                    ! van Genutchen "n" parameter (-)
 real(dp),intent(in)              :: vGn_m                    ! van Genutchen "m" parameter (-)
 real(dp),intent(in)              :: theta_sat                ! soil porosity (-)
 real(dp),intent(in)              :: theta_res                ! soil residual volumetric water content (-)
 real(dp),intent(in)              :: kAnisotropic             ! anisotropy factor for lateral hydraulic conductivity (-)
 real(dp),intent(in)              :: zScale_TOPMODEL          ! scale factor for TOPMODEL-ish baseflow parameterization (m)
 real(dp),intent(in)              :: bpar_VIC                 ! b-parameter in the VIC surface runoff parameterization (-)
 real(dp),intent(in)              :: specificYield            ! specific yield (-)
 real(dp),intent(in)              :: specificStorage          ! specific storage coefficient (m-1)
 real(dp),intent(in)              :: f_impede                 ! ice impedence factor (-)
 ! state variables
 real(dp),intent(in)              :: mLayerVolFracIce(:)      ! volumetric fraction of ice at the start of the time step (-)
 real(dp),intent(in)              :: mLayerVolFracLiq(:)      ! volumetric fraction of liquid water at the start of the time step (-)
 real(dp),intent(in)              :: mLayerMatricHead(:)      ! matric head in each layer at the start of the time step (m)
 real(dp),intent(in)              :: scalarAquiferStorage     ! aquifer storage at the start of the time step (m)
 real(dp),intent(in)              :: scalarSfcMeltPond        ! ponded water caused by melt of the "snow without a layer" (kg m-2)
 ! saturated hydraulic conductivity
 real(dp),intent(in)              :: mLayerSatHydCond(:)      ! saturated hydraulic conductivity at the mid-point of each layer (m s-1)
 real(dp),intent(in)              :: iLayerSatHydCond(0:)     ! saturated hydraulic conductivity at the interface of each layer (m s-1)
 ! transpiration (from energy routines)
 real(dp),intent(in)              :: mLayerInitTranspire(:)   ! transpiration loss from each soil layer at the start of the time step (m s-1)
 real(dp),intent(in)              :: mLayerTranspire(:)       ! transpiration loss from each soil layer (m s-1)
 real(dp),intent(in)              :: scalarInitAquiferTranspire ! transpiration loss from the aquifer at the start-of-step (m s-1)
 real(dp),intent(in)              :: scalarAquiferTranspire     ! transpiration loss from the aquifer at the start-of-step (m s-1)
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! ***** output variables
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! diagnostic scalar variables
 real(dp),intent(out)             :: scalarRainPlusMelt       ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
 real(dp),intent(out)             :: scalarSurfaceRunoff      ! surface runoff (m s-1)
 real(dp),intent(out)             :: scalarWaterTableDepth    ! water table depth (m)
 real(dp),intent(out)             :: scalarInitAquiferRcharge ! recharge to the aquifer at the start-of-step (m s-1)
 real(dp),intent(out)             :: scalarAquiferRcharge     ! recharge to the aquifer at the end-of-step (m s-1)
 real(dp),intent(out)             :: scalarInitAquiferBaseflow ! baseflow from the aquifer at the start-of-step (m s-1)
 real(dp),intent(out)             :: scalarAquiferBaseflow    ! baseflow from the aquifer at the end-of-step (m s-1)
 ! diagnostic variables for each layer
 real(dp),intent(out)             :: mLayerdTheta_dPsi(:)     ! derivative in the soil water characteristic w.r.t. psi (m-1)
 real(dp),intent(out)             :: mLayerdPsi_dTheta(:)     ! derivative in the soil water characteristic w.r.t. theta (m)
 real(dp),intent(out)             :: iLayerInitLiqFluxSoil(0:) ! liquid flux at soil layer interfaces at the start of the time step (m s-1)
 real(dp),intent(out)             :: iLayerLiqFluxSoil(0:)    ! liquid flux at soil layer interfaces at the end of the time step (m s-1)
 real(dp),intent(out)             :: mLayerInitEjectWater(:)  ! water ejected from each soil layer at the start-of-step (m s-1)
 real(dp),intent(out)             :: mLayerEjectWater(:)      ! water ejected from each soil layer (m s-1)
 ! output variables from the soilHydrol routine
 real(dp),intent(out)             :: mLayerMatricHeadNew(:)   ! matric head in each layer at the next iteration (m)
 real(dp),intent(out)             :: mLayerVolFracLiqNew(:)   ! volumetric fraction of liquid water at the next iteration (m)
 real(dp),intent(out)             :: scalarAquiferStorageNew  ! aquifer storage (m)
 integer(i4b),intent(out)         :: err                      ! error code
 character(*),intent(out)         :: message                  ! error message
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! ***** local variables
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables (general)
 character(LEN=256)               :: cmessage                 ! error message of downwind routine
 logical(lgt)                     :: printflag                ! flag to print crap to the screen
 integer(i4b)                     :: iLayer                   ! layer index
 real(dp)                         :: drainablePorosity        ! drainable porosity (-)
 logical(lgt),parameter           :: calcJacobian=.true.      ! flag to compute the Jacobian matrix
 ! check need to compute initial fluxes
 integer(i4b)                     :: nFlux                    ! index for flux calculation
 integer(i4b)                     :: ifluxInit                ! starting index for flux calculations (0 or 1)
 real(dp)                         :: maxdiffMatric(1)         ! used to check if we are starting on the first iteration
 ! flux derivatives
 real(dp),dimension(0:nLevels)    :: dq_dStateAbove           ! change in the flux in layer interfaces w.r.t. state variable in the layer above
 real(dp),dimension(0:nLevels)    :: dq_dStateBelow           ! change in the flux in layer interfaces w.r.t. state variable in the layer below
 real(dp),dimension(nLevels)      :: mLayerEjectWaterDeriv    ! derivative in water ejected from each soil layer (m s-1)
 real(dp)                         :: scalarAquiferBaseflowDeriv ! derivative in the baseflow fluw w.r.t. aquifer storage (s-1)
 ! residuals
 real(dp)                         :: scalarAquiferResidual    ! residual in aquifer storage (m) 
 ! tri-diagonal solution
 integer(i4b)                     :: nState                   ! number of state variables
 real(dp)                         :: wtim                     ! weighted time (s)
 real(dp),allocatable             :: d_m1(:)                  ! sub-diagonal elements of the tridiagonal system (m-1)
 real(dp),allocatable             :: diag(:)                  ! diagonal elements of the tridiagonal system (m-1)
 real(dp),allocatable             :: d_p1(:)                  ! super-diagonal elements of the tridiagonal system (m-1)
 real(dp),dimension(nLevels)      :: rVec                     ! right-hand-side vector (-)
 real(dp),dimension(nLevels)      :: mLayerMatricHeadDiff     ! iteration increment for matric head (m)
 real(dp),dimension(nLevels)      :: mLayerVolFracLiqDiff     ! iteration increment for volumetric fraction of liquid water (m)
 ! initialize error control
 err=0; message='soilHydrol_muster/'

 ! initilaize printflag
 printflag=.false.

 ! check the size of the input arguments
 if(any((/size(mLayerMatricHeadIter),size(mLayerVolFracIceIter),size(mLayerVolFracLiqIter),size(mLayerMatricHeadNew)/) /= nLevels)) then
  err=20; message=trim(message)//'size mis-match for the input arguments'; return
 endif

 ! check the b/c make sense
 if(nSnow>0 .and. ixBcUpperSoilHydrology==prescribedHead)then
  err=20; message=trim(message)//'using diriclet bc for the top of the soil zone when snow is present'; return
 endif

 ! only allow the moving boundary gw parameterization with the moisture-based form of Richards' equation
 if(ixGroundwater==movingBoundary)then
  if(ixRichards/=moisture)then; err=20; message=trim(message)//"moving boundary gw parameterization only allowed with the moisture-based from of Richards' equation"; return; endif
 endif

 ! define upper boundary fluxes (m s-1)
 if(ixBcUpperSoilHydrology==liquidFlux)then
  if(nSnow==0) scalarRainPlusMelt = scalarRainfall/iden_water + (scalarSfcMeltPond/dt)/iden_water  ! rainfall plus melt of the snow without a layer (convert to m s-1)
  if(nSnow>0)  scalarRainPlusMelt = scalarLiqFluxSnow                                              ! liquid water flux from the base of the snowpack (m s-1)
 endif
 if(ixBcUpperSoilHydrology==prescribedHead)then
  scalarRainPlusMelt = 0._dp
 endif

 ! get the number of state variables, and allocate space for the tri-diagonal matrix
 select case(ixBcLowerSoilHydrology)
  case(groundwaterCouple); nState = nLevels+1
  case default;            nState = nLevels
 end select
 allocate(d_m1(nState-1),diag(nState),d_p1(nState-1),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to allocate space for the tri-diagonal vectors'; return; endif

 ! initialize the need to compute inital fluxes
 ifluxInit=1  ! (don't compute initial fluxes unless instructed otherwise)

 ! *****
 ! check the need to calculate initial fluxes
 if(iter==1)then
  maxdiffMatric = maxval(abs(mLayerMatricHeadIter - mLayerMatricHead))
  if(maxdiffMatric(1) > 1.e-8_dp) ifluxInit=0
 endif


 ! *****
 ! compute the initial flux if necessary (iFluxInit=0)
 do nFlux=iFluxInit,1

  ! compute **initial** fluxes in each layer of the soil profile
  ! NOTE: use state variables at the start of the step
  if(nFlux==0)then
   call computeFlux(&
                    ! input: model control variables
                    .false.,                   & ! intent(in): flag to indicate if derivatives are desired 
                    ! input: trial state variables  (NOTE: use vectors from the start of the step)
                    mLayerMatricHead,          & ! intent(in): matric head (m)
                    mLayerVolFracLiq,          & ! intent(in): volumetric fraction of liquid water (-)
                    mLayerVolFracIce,          & ! intent(in): volumetric fraction of ice (-)
                    scalarAquiferStorage,      & ! intent(in): aquifer storage at the start of the step (m)
                    ! output: derivative in the soil water characteristic
                    mLayerdPsi_dTheta,         & ! intent(out): derivative in the soil water characteristic
                    mLayerdTheta_dPsi,         & ! intent(out): derivative in the soil water characteristic
                    ! output: diagnostic variables
                    scalarSurfaceRunoff,       & ! intent(out): surface runoff (m s-1)
                    scalarWaterTableDepth,     & ! intent(out): water table depth (m)                         
                    ! output: fluxes
                    iLayerInitLiqFluxSoil,     & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                    mLayerInitEjectWater,      & ! intent(out): water ejected because pore volume is close to capacity (m s-1)
                    scalarInitAquiferRcharge,  & ! intent(out): recharge to the aquifer (m s-1)
                    scalarInitAquiferBaseflow, & ! intent(out): baseflow flux (m s-1)
                    ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                    dq_dStateAbove,            & ! intent(out): derivatives in the flux w.r.t. state variable in the layer above (m s-1 or s-1)
                    dq_dStateBelow,            & ! intent(out): derivatives in the flux w.r.t. state variable in the layer below (m s-1 or s-1)
                    mLayerEjectWaterDeriv,     & ! intent(out): derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
                    scalarAquiferBaseflowDeriv,& ! intent(out): derivative in the baseflow fluw w.r.t. aquifer storage (s-1)
                    ! output: error control
                    err,cmessage)                ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  
  ! compute fluxes
  ! NOTE: use trial state variables
  else
   call computeFlux(&
                    ! input: model control variables
                    .true.,                    & ! intent(in): flag to indicate if derivatives are desired 
                    ! input: trial state variables  (NOTE: use vectors from the start of the step)
                    mLayerMatricHeadIter,      & ! intent(in): matric head (m)
                    mLayerVolFracLiqIter,      & ! intent(in): volumetric fraction of liquid water (-)
                    mLayerVolFracIceIter,      & ! intent(in): volumetric fraction of ice (-)
                    scalarAquiferStorageIter,  & ! intent(in): aquifer storage at the start of the step (m)
                    ! output: derivative in the soil water characteristic
                    mLayerdPsi_dTheta,         & ! intent(out): derivative in the soil water characteristic
                    mLayerdTheta_dPsi,         & ! intent(out): derivative in the soil water characteristic
                    ! output: diagnostic variables
                    scalarSurfaceRunoff,       & ! intent(out): surface runoff (m s-1)
                    scalarWaterTableDepth,     & ! intent(out): water table depth (m)                         
                    ! output: fluxes
                    iLayerLiqFluxSoil,         & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                    mLayerEjectWater,          & ! intent(out): water ejected because pore volume is close to capacity (m s-1)
                    scalarAquiferRcharge,      & ! intent(out): recharge to the aquifer (m s-1)
                    scalarAquiferBaseflow,     & ! intent(out): baseflow flux (m s-1)
                    ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                    dq_dStateAbove,            & ! intent(out): derivatives in the flux w.r.t. state variable in the layer above (m s-1 or s-1)
                    dq_dStateBelow,            & ! intent(out): derivatives in the flux w.r.t. state variable in the layer below (m s-1 or s-1)
                    mLayerEjectWaterDeriv,     & ! intent(out): derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
                    scalarAquiferBaseflowDeriv,& ! intent(out): derivative in the baseflow flux w.r.t. aquifer storage (s-1)
                    ! output: error control
                    err,cmessage)                ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! ***** assign initial fluxes, if necessary
   if(iter==1)then
    if(nFlux==ifluxInit)then   ! (handles special case where there is no separate computation for start-of-step fluxes)
     mLayerInitEjectWater      = mLayerEjectWater
     iLayerInitLiqFluxSoil     = iLayerLiqFluxSoil
     scalarInitAquiferRcharge  = scalarAquiferRcharge
     scalarInitAquiferBaseflow = scalarAquiferBaseflow
    endif  ! (if computing initial fluxes)
   endif  ! (if the first iteration)

  endif ! (if computing standard fluxes)

 end do  ! looping through initial vectors

 ! *****
 ! compute the residual vector
 call liqResidual(&
                  ! input: control variables
                  dt,                         & ! intent(in): length of the time step (s)
                  wimplicit,                  & ! intent(in): weight assigned to the start-of-step
                  ! input: coordinate variables
                  mLayerDepth,                & ! intent(in): depth of each layer (m)
                  ! input: initial flux vectors (start of the time step)
                  iLayerInitLiqFluxSoil,      & ! intent(in): initial liquid flux at layer interfaces (m s-1)
                  mLayerInitEjectWater,       & ! intent(in): initial water ejected from each layer because volumetric liquid water is approaching porosity (m s-1)
                  mLayerInitTranspire,        & ! intent(in): initial transpiration from each layer -- from energy routine (m s-1)
                  ! input: trial flux vectors
                  iLayerLiqFluxSoil,          & ! intent(in): trial liquid flux at layer interfaces (m s-1)
                  mLayerEjectWater,           & ! intent(in): trial water ejected from each layer because volumetric liquid water is approaching porosity (m s-1)
                  mLayerTranspire,            & ! intent(in): trial transpiration from each layer -- from energy routine (m s-1)
                  ! input: initial state vectors (start of the time step)
                  mLayerVolFracLiq,           & ! intent(in): initial volumetric liquid water content (-)
                  mLayerVolFracIce,           & ! intent(in): initial volumetric ice content (-)
                  ! input: trial state vectors
                  mLayerVolFracLiqIter,       & ! intent(in): trial volumetric liquid water content (-)
                  mLayerVolFracIceIter,       & ! intent(in): trial volumetric ice content (-)
                  ! output: residual vector (-)
                  rVec,                       & ! intent(out): residual vector (-)
                  ! output: error control
                  err,cmessage)                 ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! *****
 ! compute the residual for the groundwater store
 call gw_residual(&
                  ! input: control variables
                  dt,                         & ! intent(in): length of the time step (s)
                  wimplicit,                  & ! intent(in): weight assigned to the start-of-step
                  ! input: model states
                  scalarAquiferStorage,       & ! intent(in): aquifer storage at the start of the step (m)
                  scalarAquiferStorageIter,   & ! intent(in): trial value of aquifer storage (m)
                  ! input: start-of-step fluxes
                  scalarInitAquiferTranspire, & ! intent(in): transpiration from the aquifer (m s-1)
                  scalarInitAquiferRcharge,   & ! intent(in): recharge to the aquifer        (m s-1)
                  scalarInitAquiferBaseflow,  & ! intent(in): baseflow from the aquifer      (m s-1)
                  ! input: end-of-step fluxes
                  scalarAquiferTranspire,     & ! intent(in): transpiration from the aquifer (m s-1)
                  scalarAquiferRcharge,       & ! intent(in): recharge to the aquifer        (m s-1)
                  scalarAquiferBaseflow,      & ! intent(in): baseflow from the aquifer      (m s-1)
                  ! output: aquifer residual
                  scalarAquiferResidual,      & ! intent(out): aquifer residual (m)
                  ! output: error control
                  err,cmessage)                 ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! *****
 ! populate the tri-diagnonal matrices for the soil layers
 wtim = (1._dp - wimplicit)*dt  ! weighted time
 ! (get off-diagonal elements)
 d_m1(1:nLevels-1) = (wtim/mLayerDepth(1:nLevels-1))*(-dq_dStateAbove(1:nLevels-1) )
 d_p1(1:nLevels-1) = (wtim/mLayerDepth(1:nLevels-1))*( dq_dStateBelow(1:nLevels-1) )
 ! (get diagonal elements)
 select case(ixRichards)
  case(moisture); diag(1:nLevels) = (wtim/mLayerDepth(1:nLevels))*(-dq_dStateBelow(0:nLevels-1) + dq_dStateAbove(1:nLevels) + mLayerEjectWaterDeriv(1:nLevels)) + 1._dp
  case(mixdform); diag(1:nLevels) = (wtim/mLayerDepth(1:nLevels))*(-dq_dStateBelow(0:nLevels-1) + dq_dStateAbove(1:nLevels) + mLayerEjectWaterDeriv(1:nLevels)) + mLayerdTheta_dPsi(1:nLevels)
  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
 endselect

 ! *****
 ! populate the tri-diagnonal matrices for the aquifer
 if(ixBcLowerSoilHydrology == groundwaterCouple)then
  ! (error check)
  if(ixRichards == mixdform)then; err=20; message=trim(message)//"require the moisture based form of Richards' equation with a moving lower boundary"; return; endif
  ! (compute the off-diagonal elements)
  d_p1(nLevels)   = (wtim/mLayerDepth(nLevels))*( dq_dStateBelow(nLevels) )  ! change in recharge flux w.r.t change in the aquifer storage (m-1)
  d_m1(nLevels)   =  wtim                      *(-dq_dStateAbove(nLevels) )  ! change in recharge flux w.r.t. change in soil moisture in the lowest unsaturated node (m)
  ! (compute the diagonal)
  diag(nLevels+1) = wtim * (-dq_dStateBelow(nLevels) + scalarAquiferBaseflowDeriv) + 1._dp  ! flux derivatives (s-1), then diagonal here is dimensionless
 endif
 
 ! *****
 ! test the Jacobian
 if(calcJacobian)then
  call cmpJacobian(&
                   ! input: model control
                   ixRichards,                   & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                   ixBcLowerSoilHydrology,       & ! intent(in): index defining the type of boundary conditions
                   ! input: state variables
                   mLayerMatricHeadIter,         & ! intent(in): matric head in each layer (m)
                   mLayerVolFracLiqIter,         & ! intent(in): volumetric liquid water content (-)
                   mLayerVolFracIceIter,         & ! intent(in): volumetric ice content in each layer (-)
                   scalarAquiferStorageIter,     & ! intent(in): aquifer storage (m)
                   ! output: error control
                   err,cmessage)                   ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 endif


 ! deallocate space for the tri-diagonal vectors
 deallocate(d_m1,diag,d_p1,stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to deallocate space for the tri-diagonal vectors'; return; endif


 stop 'initial test'








 ! ***************************************************************************************************************************
 ! ***************************************************************************************************************************
 ! ***************************************************************************************************************************
 ! ***************************************************************************************************************************
 ! ***** END OF MAIN SUBROUTINE **********************************************************************************************
 ! ***************************************************************************************************************************
 ! ***************************************************************************************************************************
 ! ***************************************************************************************************************************

 contains

  ! ************************************************************************************************
  ! internal subroutine: compute the vector of fluxes at layer interfaces
  ! ************************************************************************************************
  subroutine computeFlux(&
                         ! input: model control
                         deriv_desired,             & ! intent(in): flag indicating if derivatives are desired
                         ! input: trial state variables (NOTE: use "trial" instead of "iter" to avoid using the "iter" vectors defined in main subroutine)
                         mLayerMatricHeadTrial,     & ! intent(in): matric head (m)
                         mLayerVolFracLiqTrial,     & ! intent(in): volumetric fraction of liquid water (-)
                         mLayerVolFracIceTrial,     & ! intent(in): volumetric fraction of ice (-)
                         scalarAquiferStorageTrial, & ! intent(in): aquifer storage at the start of the step (m)
                         ! output: derivative in the soil water characteristic
                         mLayerdPsi_dTheta,         & ! intent(out): derivative in the soil water characteristic
                         mLayerdTheta_dPsi,         & ! intent(out): derivative in the soil water characteristic
                         ! output: diagnostic variables
                         scalarSurfaceRunoff,       & ! intent(out): surface runoff (m s-1)
                         scalarWaterTableDepth,     & ! intent(out): water table depth (m)                         
                         ! output: fluxes at layer interfaces
                         iLayerLiqFluxSoil,         & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                         mLayerEjectWater,          & ! intent(out): water ejected because pore volume is close to capacity (m s-1)
                         scalarAquiferRcharge,      & ! intent(out): recharge to the aquifer (m s-1)
                         scalarAquiferBaseflow,     & ! intent(out): baseflow flux (m s-1)
                         ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                         dq_dStateAbove,            & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer above (m s-1)
                         dq_dStateBelow,            & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer below (m s-1)
                         mLayerEjectWaterDeriv,     & ! intent(out): derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
                         scalarAquiferBaseflowDeriv,& ! intent(out): derivative in the baseflow fluw w.r.t. aquifer storage (s-1)
                         ! output: error control
                         err,message)                ! intent(out): error control
  USE soil_utils_module,only:hydCond_psi  ! compute hydraulic conductivity as a function of matric head
  USE soil_utils_module,only:hydCond_liq  ! compute hydraulic conductivity as a function of volumetric liquid water content
  implicit none
  ! input: model control
  logical(lgt),intent(in)          :: deriv_desired               ! flag indicating if derivatives are desired
  ! trial model state variables
  real(dp),intent(in)              :: mLayerMatricHeadTrial(:)    ! matric head in each layer at the current iteration (m)
  real(dp),intent(in)              :: mLayerVolFracLiqTrial(:)    ! volumetric fraction of liquid water at the current iteration (-)
  real(dp),intent(in)              :: mLayerVolFracIceTrial(:)    ! volumetric fraction of ice at the current iteration (-)
  real(dp),intent(in)              :: scalarAquiferStorageTrial   ! aquifer storage at the current iteration (m)
  ! output: derivative in the soil water characteristic
  real(dp),intent(out)             :: mLayerdPsi_dTheta(:)        ! derivative in the soil water characteristic
  real(dp),intent(out)             :: mLayerdTheta_dPsi(:)        ! derivative in the soil water characteristic
  ! output: diagnostic variables
  real(dp),intent(out)             :: scalarSurfaceRunoff         ! surface runoff (m s-1)
  real(dp),intent(out)             :: scalarWaterTableDepth       ! water table depth (m)
  ! output: liquid fluxes at layer interfaces
  real(dp),intent(out)             :: iLayerLiqFluxSoil(0:)       ! liquid flux at soil layer interfaces (m s-1)
  real(dp),intent(out)             :: mLayerEjectWater(:)         ! water ejected because pore volume is close to capacity (m s-1)
  real(dp),intent(out)             :: scalarAquiferRcharge        ! recharge to the aquifer (m s-1)
  real(dp),intent(out)             :: scalarAquiferBaseflow       ! baseflow flux (m s-1)
  ! output: derivatives in fluxes w.r.t. state variables in the layer above and layer below (m s-1)
  real(dp),intent(out)             :: dq_dStateAbove(0:)          ! change in the flux in layer interfaces w.r.t. state variables in the layer above
  real(dp),intent(out)             :: dq_dStateBelow(0:)          ! change in the flux in layer interfaces w.r.t. state variables in the layer below
  real(dp),intent(out)             :: mLayerEjectWaterDeriv(:)    ! derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
  real(dp),intent(out)             :: scalarAquiferBaseflowDeriv  ! derivative in the baseflow fluw w.r.t. aquifer storage (s-1)
  ! output: error control
  integer(i4b),intent(out)         :: err                         ! error code
  character(*),intent(out)         :: message                     ! error message
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! local variables: general
  character(LEN=256)               :: cmessage                    ! error message of downwind routine
  logical(lgt)                     :: desireAnal                  ! flag to determine if analytical derivatives are desired
  integer(i4b)                     :: iSoil                       ! index of soil layer
  integer(i4b)                     :: itry                        ! index of different flux calculations
  integer(i4b)                     :: nFlux                       ! number of flux calculations required (>1 = numerical derivatives with one-sided finite differences)
  integer(i4b),parameter           :: unperturbed=1               ! named variable to identify the case of unperturbed state variables
  integer(i4b),parameter           :: perturbStateAbove=3         ! named variable to identify the case where we perturb the state above
  integer(i4b),parameter           :: perturbStateBelow=2         ! named variable to identify the case where we perturb the state below
  integer(i4b)                     :: ixPerturb                   ! index of element in 2-element vector to perturb
  integer(i4b)                     :: ixOriginal                  ! index of perturbed element in the original vector
  real(dp)                         :: scalardPsi_dTheta           ! derivative in soil water characteristix, used for perturbations when computing numerical derivatives
  ! local variables: scalar trial values
  real(dp)                         :: scalarVolFracLiqTrial       ! trial value of volumetric liquid water content (-)
  real(dp)                         :: scalarMatricHeadTrial       ! trial value of matric head (m)
  real(dp)                         :: scalarWaterTableDepthTrial  ! trial value of depth to the water table (m)
  real(dp)                         :: scalarHydCondTrial          ! trial value of hydraulic conductivity (m s-1)
  ! local variables: vector trial values (2-element vectors)
  real(dp),dimension(2)            :: vectorVolFracLiqTrial       ! trial value of volumetric liquid water content (-)
  real(dp),dimension(2)            :: vectorMatricHeadTrial       ! trial value of matric head (m)
  real(dp),dimension(2)            :: vectorHydCondTrial          ! trial value of hydraulic conductivity (m s-1)
  real(dp),dimension(2)            :: vectorDiffuseTrial          ! trial value of hydraulic diffusivity (m2 s-1)
  ! local variables: temporary variables used to store fluxes (used to compute numerical derivatives)
  real(dp)                         :: scalarFlux                  ! vertical flux (m s-1)
  real(dp)                         :: scalarFlux_dStateAbove      ! vertical flux with perturbation to the state above (m s-1)
  real(dp)                         :: scalarFlux_dStateBelow      ! vertical flux with perturbation to the state below (m s-1 [soil] or s-1 [aquifer])
  ! local variables: transmittance
  real(dp),dimension(nLevels)      :: iceImpedeFac                ! ice impedence factor at layer mid-points (-)
  real(dp),dimension(nLevels)      :: mLayerHydCond               ! hydraulic conductivity at layer mid-point (m s-1)
  real(dp),dimension(nLevels)      :: mLayerDiffuse               ! diffusivity at layer mid-point (m2 s-1)
  real(dp),dimension(0:nLevels)    :: iLayerHydCond               ! hydraulic conductivity at layer interface (m s-1)
  real(dp),dimension(0:nLevels)    :: iLayerDiffuse               ! diffusivity at layer interface (m2 s-1)
  real(dp),dimension(nLevels)      :: dHydCond_dVolLiq            ! derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
  real(dp),dimension(nLevels)      :: dDiffuse_dVolLiq            ! derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
  real(dp),dimension(nLevels)      :: dHydCond_dMatric            ! derivative in hydraulic conductivity w.r.t matric head (m s-1)
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message='computeFlux/'

  ! check the need to compute analytical derivatives
  if(deriv_desired .and. ixDerivMethod==analytical)then
   desireAnal = .true.
  else
   desireAnal = .false.
  endif

  ! check the need to compute numerical derivatives
  if(deriv_desired .and. ixDerivMethod==numerical)then
   nFlux=3  ! compute the derivatives using one-sided finite differences
  else
   nFlux=1  ! compute analytical derivatives
  endif



  ! *************************************************************************************************************************************************
  ! *************************************************************************************************************************************************

  ! compute the drainable porosity, and the depth to the water table
  if(ixGroundwater==movingBoundary)then
   call waterTablePosition(&
                           mLayerVolFracLiqTrial,      & ! intent(in):  volumetric liquid water content in each soil layer (-)
                           mLayerVolFracIceTrial,      & ! intent(in):  volumetric ice content in each soil layer (-)
                           scalarAquiferStorageTrial,  & ! intent(in):  aquifer storage (m)
                           iLayerHeight,               & ! intent(in):  height of each interface (m)
                           mLayerDepth,                & ! intent(in):  depth of each soil layer (m)
                           theta_sat,                  & ! intent(in):  soil porosity (-)
                           specificYield,              & ! intent(in):  fraction of water volume drained by gravity in an unconfined aquifer (-)
                           drainablePorosity,          & ! intent(out): drainable porosity (-)
                           scalarWaterTableDepth,      & ! intent(out): water table depth (m)
                           err,cmessage)                 ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   ! check that the water table depth is below the mid-point of the lowest soil layer
   if(scalarWaterTableDepth < mLayerHeight(nLevels))then; err=20; message=trim(message)//'water table is above the mid-point of the lowest unsaturated layer'; return; endif
  else  ! (if water table depth is not used)
   scalarWaterTableDepth = valueMissing  ! ** deliberately cause problems
  endif



  ! *************************************************************************************************************************************************
  ! *************************************************************************************************************************************************

  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! compute diagnostic variables at the nodes throughout the soil profile
  do iSoil=1,nLevels ! (loop through soil layers)
   call diagv_node(&
                   ! input: model control
                   desireAnal,                      & ! intent(in): flag indicating if derivatives are desired
                   ixRichards,                      & ! intent(in): index defining the option for Richards' equation (moisture or mixdform)
                   ! input: state variables
                   mLayerMatricHeadTrial(iSoil),    & ! intent(in):  matric head in each layer
                   mLayerVolFracLiqTrial(iSoil),    & ! intent(in):  volumetric liquid water content in each soil layer (-)
                   mLayerVolFracIceTrial(iSoil),    & ! intent(in):  volumetric ice content in each soil layer (-)
                   ! input: soil parameters
                   vGn_alpha,                       & ! intent(in): van Genutchen "alpha" parameter (m-1)
                   vGn_n,                           & ! intent(in): van Genutchen "n" parameter (-)
                   VGn_m,                           & ! intent(in): van Genutchen "m" parameter (-)
                   theta_sat,                       & ! intent(in): soil porosity (-)
                   theta_res,                       & ! intent(in): soil residual volumetric water content (-)
                   f_impede,                        & ! intent(in): ice impedence factor (-)
                   ! input: saturated hydraulic conductivity
                   mLayerSatHydCond(iSoil),         & ! intent(in): saturated hydraulic conductivity at the mid-point of each layer (m s-1)
                   ! output: derivative in the soil water characteristic
                   mLayerdPsi_dTheta(iSoil),        &  ! derivative in the soil water characteristic
                   mLayerdTheta_dPsi(iSoil),        &  ! derivative in the soil water characteristic
                   ! output: transmittance
                   mLayerHydCond(iSoil),            & ! intent(out): hydraulic conductivity at layer mid-points (m s-1)
                   mLayerDiffuse(iSoil),            & ! intent(out): diffusivity at layer mid-points (m2 s-1)
                   iceImpedeFac(iSoil),             & ! intent(out): ice impedence factor in each layer (-)
                   ! output: transmittance derivatives
                   dHydCond_dVolLiq(iSoil),         & ! intent(out): derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
                   dDiffuse_dVolLiq(iSoil),         & ! intent(out): derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
                   dHydCond_dMatric(iSoil),         & ! intent(out): derivative in hydraulic conductivity w.r.t matric head (m s-1)
                   ! output: error control
                   err,cmessage)                      ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  end do  ! (looping through soil layers)


  ! *************************************************************************************************************************************************
  ! *************************************************************************************************************************************************

  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! compute infiltraton at the surface and its derivative w.r.t. mass in the upper soil layer
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! either one or multiple flux calls, depending on if using analytical or numerical derivatives
  do itry=nFlux,1,-1  ! (work backwards to ensure all computed fluxes come from the un-perturbed case)

   ! =====
   ! get input state variables...
   ! ============================
   ! identify the type of perturbation
   select case(itry)

    ! un-perturbed case
    case(unperturbed)
     scalarVolFracLiqTrial = mLayerVolFracLiqTrial(1)
     scalarMatricHeadTrial = mLayerMatricHeadTrial(1)

    ! perturb soil state (one-sided finite differences)
    case(perturbStateBelow)
     ! (perturbation depends on the form of Richards' equation)
     select case(ixRichards)
      case(moisture)
       scalarVolFracLiqTrial = mLayerVolFracLiqTrial(1) + dx
       scalarMatricHeadTrial = mLayerMatricHeadTrial(1)
      case(mixdform)
       scalarVolFracLiqTrial = mLayerVolFracLiqTrial(1)
       scalarMatricHeadTrial = mLayerMatricHeadTrial(1) + dx
      case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
     end select ! (form of Richards' equation)

    ! cannot perturb state above (does not exist) -- so keep cycling
    case(perturbStateAbove); cycle

    ! check for an unknown perturbation 
    case default; err=10; message=trim(message)//"unknown perturbation"; return

   end select ! (type of perturbation)

   ! =====
   ! compute surface flux and its derivative...
   ! ==========================================
   call surfaceFlx(&
                   ! input: model control
                   desireAnal,                      & ! intent(in): flag indicating if derivatives are desired
                   ixRichards,                      & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                   ixBcUpperSoilHydrology,          & ! intent(in): index defining the type of boundary conditions (neumann or diriclet)
                   ! input: state variables
                   scalarMatricHeadTrial,           & ! intent(in): matric head in he upper-most soil layer (m)
                   scalarVolFracLiqTrial,           & ! intent(in): volumetric liquid water content the upper-most soil layer (-)
                   mLayerVolFracIceTrial(1),        & ! intent(in): volumetric ice content in the upper-most soil layer (-)
                   ! input: depth of upper-most soil layer (m)
                   mLayerDepth(1),                  & ! intent(in): depth of upper-most soil layer (m)
                   ! input: boundary conditions
                   upperBoundHead,                  & ! intent(in): upper boundary condition (m)
                   upperBoundTheta,                 & ! intent(in): upper boundary condition (-)
                   ! input: flux at the upper boundary
                   scalarRainPlusMelt,              & ! intent(in): rain plus melt (m s-1)
                   ! input: derivative in soil water characteristix
                   mLayerdPsi_dTheta(1),            & ! intent(in): derivative of the soil moisture characteristic w.r.t. theta (m)
                   mLayerdTheta_dPsi(1),            & ! intent(in): derivative of the soil moisture characteristic w.r.t. psi (m-1)
                   ! input: transmittance
                   iLayerSatHydCond(0),             & ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
                   iceImpedeFac(1),                 & ! intent(in): ice impedence factor in the upper-most soil layer (-)
                   ! input: soil parameters
                   vGn_alpha,                       & ! intent(in): van Genutchen "alpha" parameter (m-1)
                   vGn_n,                           & ! intent(in): van Genutchen "n" parameter (-)
                   VGn_m,                           & ! intent(in): van Genutchen "m" parameter (-)
                   theta_sat,                       & ! intent(in): soil porosity (-)
                   theta_res,                       & ! intent(in): soil residual volumetric water content (-)
                   bpar_VIC,                        & ! intent(in): b-parameter in the VIC surface runoff parameterization (-)
                   ! output: hydraulic conductivity and diffusivity at the surface
                   iLayerHydCond(0),                & ! intent(out): hydraulic conductivity at the surface (m s-1)
                   iLayerDiffuse(0),                & ! intent(out): hydraulic diffusivity at the surface (m2 s-1)
                   ! output: fluxes at layer interfaces and surface runoff
                   scalarSurfaceRunoff,             & ! intent(out): surface runoff (m s-1)
                   iLayerLiqFluxSoil(0),            & ! intent(out): surface infiltration (m s-1)
                   ! output: deriavtives in surface infiltration w.r.t. volumetric liquid water (m s-1) and matric head (s-1) in the upper-most soil layer
                   dq_dStateBelow(0),               & ! intent(out): derivative in surface infiltration w.r.t. state variable in the upper-most soil layer (m s-1 or s-1)
                   ! output: error control
                   err,cmessage)                 ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! set derivative w.r.t. state above to missing (does not exist)
   dq_dStateAbove(0) = valueMissing

   ! get copies of surface flux to compute numerical derivatives
   if(deriv_desired .and. ixDerivMethod==numerical)then
    select case(itry)
     case(unperturbed);       scalarFlux             = iLayerLiqFluxSoil(0)
     case(perturbStateBelow); scalarFlux_dStateBelow = iLayerLiqFluxSoil(0)
     case default; err=10; message=trim(message)//"unknown perturbation"; return
    end select
   endif

  end do  ! (looping through different flux calculations -- one or multiple calls depending if desire for numerical or analytical derivatives)

  ! compute numerical derivatives
  if(deriv_desired .and. ixDerivMethod==numerical)then
   dq_dStateBelow(0) = (scalarFlux_dStateBelow - scalarFlux)/dx ! change in surface flux w.r.t. change in the soil moisture in the top soil layer (m s-1)
  endif
  !print*, 'ixDerivMethod, dq_dStateBelow(0) = ', ixDerivMethod, dq_dStateBelow(0)

  ! *************************************************************************************************************************************************
  ! *************************************************************************************************************************************************

  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! * compute fluxes and derivatives at layer interfaces...
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! loop through soil layers
  do iLayer=1,nLevels-1

   ! either one or multiple flux calls, depending on if using analytical or numerical derivatives
   do itry=nFlux,1,-1  ! (work backwards to ensure all computed fluxes come from the un-perturbed case)

    ! =====
    ! determine layer to perturb
    ! ============================
    select case(itry)
     case(unperturbed);       ixPerturb = 0
     case(perturbStateAbove); ixPerturb = 1
     case(perturbStateBelow); ixPerturb = 2
     case default; err=10; message=trim(message)//"unknown perturbation"; return
    end select ! (identifying layer to of perturbation)
    ! determine the index in the original vector
    ixOriginal = iLayer + (ixPerturb-1)

    ! =====
    ! get input state variables...
    ! ============================
    ! start with the un-perturbed case
    vectorVolFracLiqTrial(1:2) = mLayerVolFracLiqTrial(iLayer:iLayer+1)
    vectorMatricHeadTrial(1:2) = mLayerMatricHeadTrial(iLayer:iLayer+1)
    ! make appropriate perturbations
    if(ixPerturb > 0)then
     select case(ixRichards)
      case(moisture); vectorVolFracLiqTrial(ixPerturb) = vectorVolFracLiqTrial(ixPerturb) + dx
      case(mixdform); vectorMatricHeadTrial(ixPerturb) = vectorMatricHeadTrial(ixPerturb) + dx
      case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
     end select ! (form of Richards' equation)
    endif

    ! =====
    ! get hydraulic conductivty...
    ! ============================
    ! start with the un-perturbed case
    vectorHydCondTrial(1:2) = mLayerHydCond(iLayer:iLayer+1)
    vectorDiffuseTrial(1:2) = mLayerDiffuse(iLayer:iLayer+1)
    ! make appropriate perturbations
    if(ixPerturb > 0)then
     select case(ixRichards)
      case(moisture)
       scalardPsi_dTheta             = dPsi_dTheta(vectorVolFracLiqTrial(ixPerturb),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
       vectorHydCondTrial(ixPerturb) = hydCond_liq(vectorVolFracLiqTrial(ixPerturb),mLayerSatHydCond(ixOriginal),theta_res,theta_sat,vGn_m) * iceImpedeFac(ixOriginal)
       vectorDiffuseTrial(ixPerturb) = scalardPsi_dTheta * vectorHydCondTrial(ixPerturb)
      case(mixdform)
       vectorHydCondTrial(ixPerturb) = hydCond_psi(vectorMatricHeadTrial(ixPerturb),mLayerSatHydCond(ixOriginal),vGn_alpha,vGn_n,vGn_m) * iceImpedeFac(ixOriginal)
       vectorDiffuseTrial(ixPerturb) = valueMissing
      case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
     end select ! (form of Richards' equation)
    endif

    ! =====
    ! compute vertical flux at layer interface and its derivative w.r.t. the state above and the state below...
    ! =========================================================================================================
    call iLayerFlux(&
                    ! input: model control
                    desireAnal,                       & ! intent(in): flag indicating if derivatives are desired
                    ixRichards,                       & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                    ! input: state variables (adjacent layers)
                    vectorMatricHeadTrial,            & ! intent(in): matric head at the soil nodes (m)
                    vectorVolFracLiqTrial,            & ! intent(in): volumetric liquid water content at the soil nodes (-)
                    ! input: model coordinate variables (adjacent layers)
                    mLayerHeight(iLayer:iLayer+1),    & ! intent(in): height of the soil nodes (m)
                    ! input: transmittance (adjacent layers)
                    vectorHydCondTrial,               & ! intent(in): hydraulic conductivity at the soil nodes (m s-1)
                    vectorDiffuseTrial,               & ! intent(in): hydraulic diffusivity at the soil nodes (m2 s-1)
                    ! input: transmittance derivatives (adjacent layers)
                    dHydCond_dVolLiq(iLayer:iLayer+1),& ! intent(in): change in hydraulic conductivity w.r.t. change in volumetric liquid water content (m s-1)
                    dDiffuse_dVolLiq(iLayer:iLayer+1),& ! intent(in): change in hydraulic diffusivity w.r.t. change in volumetric liquid water content (m2 s-1)
                    dHydCond_dMatric(iLayer:iLayer+1),& ! intent(in): change in hydraulic conductivity w.r.t. change in matric head (s-1)
                    ! output: tranmsmittance at the layer interface (scalars)
                    iLayerHydCond(iLayer),            & ! intent(out): hydraulic conductivity at the interface between layers (m s-1)
                    iLayerDiffuse(iLayer),            & ! intent(out): hydraulic diffusivity at the interface between layers (m2 s-1)
                    ! output: vertical flux at the layer interface (scalars)
                    iLayerLiqFluxSoil(iLayer),        & ! intent(out): vertical flux of liquid water at the layer interface (m s-1)
                    ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                    dq_dStateAbove(iLayer),           & ! intent(out): derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer above (m s-1 or s-1)
                    dq_dStateBelow(iLayer),           & ! intent(out): derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer below (m s-1 or s-1) 
                    ! output: error control
                    err,cmessage)             ! intent(out): error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    !print*, 'iLayer, iLayerLiqFluxSoil(iLayer) = ', iLayer, iLayerLiqFluxSoil(iLayer)

    ! get copies of drainage flux to compute derivatives
    if(deriv_desired .and. ixDerivMethod==numerical)then
     select case(itry)
      case(unperturbed);       scalarFlux             = iLayerLiqFluxSoil(iLayer)
      case(perturbStateAbove); scalarFlux_dStateAbove = iLayerLiqFluxSoil(iLayer)
      case(perturbStateBelow); scalarFlux_dStateBelow = iLayerLiqFluxSoil(iLayer)
      case default; err=10; message=trim(message)//"unknown perturbation"; return
     end select
    endif

   end do  ! (looping through different flux calculations -- one or multiple calls depending if desire for numerical or analytical derivatives)

   ! compute numerical derivatives
   if(deriv_desired .and. ixDerivMethod==numerical)then
    dq_dStateAbove(iLayer) = (scalarFlux_dStateAbove - scalarFlux)/dx    ! change in drainage flux w.r.t. change in the state in the layer below (m s-1 or s-1)
    dq_dStateBelow(iLayer) = (scalarFlux_dStateBelow - scalarFlux)/dx    ! change in drainage flux w.r.t. change in the state in the layer below (m s-1 or s-1)
   endif
   !write(*,'(a,2(i4,1x),2(e20.10,1x))') 'ixDerivMethod, iLayer, dq_dStateBelow(iLayer), dq_dStateAbove(iLayer) = ', ixDerivMethod, iLayer, dq_dStateBelow(iLayer), dq_dStateAbove(iLayer)

  end do  ! (looping through soil layers)


  ! *************************************************************************************************************************************************
  ! *************************************************************************************************************************************************

  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! * compute drainage flux and its derivative
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! either one or multiple flux calls, depending on if using analytical or numerical derivatives
  do itry=nFlux,1,-1  ! (work backwards to ensure all computed fluxes come from the un-perturbed case)

   ! =====
   ! get input state variables...
   ! ============================
   ! identify the type of perturbation
   select case(itry)

    ! un-perturbed case
    case(unperturbed)
     scalarVolFracLiqTrial      = mLayerVolFracLiqTrial(nLevels)
     scalarMatricHeadTrial      = mLayerMatricHeadTrial(nLevels)
     scalarWaterTableDepthTrial = scalarWaterTableDepth

    ! perturb soil state (one-sided finite differences)
    case(perturbStateAbove)
     ! (perturbation depends on the form of Richards' equation)
     select case(ixRichards)
      case(moisture)
       scalarVolFracLiqTrial = mLayerVolFracLiqTrial(nLevels) + dx
       scalarMatricHeadTrial = mLayerMatricHeadTrial(nLevels)
      case(mixdform)
       scalarVolFracLiqTrial = mLayerVolFracLiqTrial(nLevels)
       scalarMatricHeadTrial = mLayerMatricHeadTrial(nLevels) + dx
      case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
     end select ! (form of Richards' equation)
     ! (only perturbing soil state here -- aquifer is unperturbed)
     scalarWaterTableDepthTrial = scalarWaterTableDepth

    ! perturb aquifer (one-sided finite differences)
    case(perturbStateBelow)
     scalarVolFracLiqTrial      = mLayerVolFracLiqTrial(nLevels)
     scalarMatricHeadTrial      = mLayerMatricHeadTrial(nLevels)
     scalarWaterTableDepthTrial = scalarWaterTableDepth + dx
    case default; err=10; message=trim(message)//"unknown perturbation"; return

   end select ! (type of perturbation)

   ! =====
   ! get hydraulic conductivty...
   ! ============================
   select case(itry)

    ! compute perturbed value of hydraulic conductivity
    case(perturbStateAbove)
     select case(ixRichards)
      case(moisture); scalarHydCondTrial = hydCond_liq(scalarVolFracLiqTrial,mLayerSatHydCond(nLevels),theta_res,theta_sat,vGn_m) * iceImpedeFac(nLevels)
      case(mixdform); scalarHydCondTrial = hydCond_psi(scalarMatricHeadTrial,mLayerSatHydCond(nLevels),vGn_alpha,vGn_n,vGn_m) * iceImpedeFac(nLevels)
     end select

    ! (use un-perturbed value)
    case default
     scalarHydCondTrial = mLayerHydCond(nLevels)        ! hydraulic conductivity at the mid-point of the lowest unsaturated soil layer (m s-1)

   end select ! (re-computing hydraulic conductivity)

   ! =====
   ! compute drainage flux and its derivative...
   ! ===========================================
   call qChargeFlx(&
                   ! input: model control
                   desireAnal,                      & ! intent(in): flag indicating if derivatives are desired
                   ixRichards,                      & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                   hc_profile,                      & ! intent(in): index defining the decrease of hydraulic conductivity with depth
                   ixBcLowerSoilHydrology,          & ! intent(in): index defining the type of boundary conditions
                   ! input: state variables
                   scalarMatricHeadTrial,           & ! intent(in): matric head in the lowest unsaturated node (m)
                   scalarVolFracLiqTrial,           & ! intent(in): volumetric liquid water content the lowest unsaturated node (-)
                   mLayerVolFracIceTrial(nLevels),  & ! intent(in): volumetric ice content in the lowest unbsaturated node (-)
                   ! input: model coordinate variables
                   mLayerDepth(nLevels),            & ! intent(in): depth of the lowest unsaturated soil layer (m)
                   mLayerHeight(nLevels),           & ! intent(in): height of the lowest unsaturated soil node (m)
                   iLayerHeight(nLevels),           & ! intent(in): effective depth of the unsaturated zone (m)
                   ! input: boundary conditions
                   lowerBoundHead,                  & ! intent(in): lower boundary condition (m)
                   lowerBoundTheta,                 & ! intent(in): lower boundary condition (-)
                   ! input: water table depth
                   scalarWaterTableDepthTrial,      & ! intent(in): depth to the water table (m)
                   ! input: derivative in the soil water characteristic
                   mLayerdPsi_dTheta(nLevels),      & ! intent(in): derivative in the soil water characteristic
                   mLayerdTheta_dPsi(nLevels),      & ! intent(in): derivative in the soil water characteristic
                   ! input: transmittance
                   iLayerSatHydCond(0),             & ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
                   iLayerSatHydCond(nLevels),       & ! intent(in): saturated hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
                   scalarHydCondTrial,              & ! intent(in): hydraulic conductivity at the node itself (m s-1)
                   dHydCond_dVolLiq(nLevels),       & ! intent(in): derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
                   dHydCond_dMatric(nLevels),       & ! intent(in): derivative in hydraulic conductivity w.r.t. matric head (s-1)
                   iceImpedeFac(nLevels),           & ! intent(in): ice impedence factor in the lower-most soil layer (-)
                   ! input: soil parameters
                   vGn_alpha,                       & ! intent(in): van Genutchen "alpha" parameter (m-1)
                   vGn_n,                           & ! intent(in): van Genutchen "n" parameter (-)
                   VGn_m,                           & ! intent(in): van Genutchen "m" parameter (-)
                   theta_sat,                       & ! intent(in): soil porosity (-)
                   theta_res,                       & ! intent(in): soil residual volumetric water content (-)
                   kAnisotropic,                    & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                   zScale_TOPMODEL,                 & ! intent(in): TOPMODEL scaling factor (m)
                   drainablePorosity,               & ! intent(in): drainable porosity (-)
                   ! output: hydraulic conductivity and diffusivity at the surface
                   iLayerHydCond(nLevels),          & ! intent(out): hydraulic conductivity at the bottom of the unsatuarted zone (m s-1)
                   iLayerDiffuse(nLevels),          & ! intent(out): hydraulic diffusivity at the bottom of the unsatuarted zone (m2 s-1)
                   ! output: drainage flux
                   iLayerLiqFluxSoil(nLevels),      & ! intent(out): drainage flux (m s-1)
                   ! output: derivatives in drainage flux
                   dq_dStateAbove(nLevels),         & ! intent(out): change in drainage flux w.r.t. change in state in lowest unsaturated node (m s-1 or s-1)
                   dq_dStateBelow(nLevels),         & ! intent(out): change in drainage flux w.r.t. change in the aquifer storage (s-1)
                   ! output: error control
                   err,cmessage)                ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   !if(deriv_desired) print*, 'itry, iLayerLiqFluxSoil(nLevels) = ', itry, iLayerLiqFluxSoil(nLevels)

   ! get copies of drainage flux to compute derivatives
   if(deriv_desired .and. ixDerivMethod==numerical)then
    select case(itry)
     case(unperturbed);       scalarFlux             = iLayerLiqFluxSoil(nLevels)
     case(perturbStateAbove); scalarFlux_dStateAbove = iLayerLiqFluxSoil(nLevels)
     case(perturbStateBelow); scalarFlux_dStateBelow = iLayerLiqFluxSoil(nLevels)
     case default; err=10; message=trim(message)//"unknown perturbation"; return
    end select
   endif

  end do  ! (looping through different flux calculations -- one or multiple calls depending if desire for numerical or analytical derivatives)

  ! compute numerical derivatives
  ! NOTE: drainage derivatives w.r.t. state below are *actually* w.r.t. water table depth, so need to be corrected for aquifer storage
  !       (note also negative sign to account for inverse relationship between water table depth and aquifer storage)
  if(deriv_desired .and. ixDerivMethod==numerical)then
   dq_dStateAbove(nLevels) =  (scalarFlux_dStateAbove - scalarFlux)/dx    ! change in drainage flux w.r.t. change in state in lowest unsaturated node (m s-1 or s-1)
   dq_dStateBelow(nLevels) = -(scalarFlux_dStateBelow - scalarFlux)/dx/drainablePorosity ! change in drainage flux w.r.t. change in the aquifer storage (s-1)
  endif

  ! print progress
  !if(deriv_desired)then
  ! print*, 'ixDerivMethod, dq_dStateAbove(nLevels) = ', ixDerivMethod, dq_dStateAbove(nLevels)
  ! print*, 'ixDerivMethod, dq_dStateBelow(nLevels) = ', ixDerivMethod, dq_dStateBelow(nLevels)
  !endif

  ! *************************************************************************************************************************************************
  ! *************************************************************************************************************************************************
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! * compute baseflow flux and its derivative
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! either one or multiple flux calls, depending on if using analytical or numerical derivatives
  do itry=nFlux,1,-1  ! (work backwards to ensure all computed fluxes come from the un-perturbed case)

   ! identify the type of perturbation
   select case(itry)
    case(unperturbed);       scalarWaterTableDepthTrial = scalarWaterTableDepth
    case(perturbStateAbove); scalarWaterTableDepthTrial = scalarWaterTableDepth + dx
    case(perturbStateBelow); cycle
    case default; err=10; message=trim(message)//"unknown perturbation"; return
   end select ! (type of perturbation)

   ! compute baseflow and its derivative
   call qBflowFlux(&
                   ! input: model control
                   desireAnal,                     & ! intent(in): flag indicating if derivatives are desired
                   ! input: "states"
                   scalarWaterTableDepthTrial,     & ! intent(in): water table depth at the start/end of the time step (m)
                   ! input: diagnostic variables
                   iLayerSatHydCond(0),            & ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
                   drainablePorosity,              & ! intent(in): drainablePorosity (-)
                   ! input: parameters
                   kAnisotropic,                   & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                   zScale_TOPMODEL,                & ! intent(in): scale factor for TOPMODEL-ish baseflow parameterization (m)
                   ! output: baseflow flux and its derivative
                   scalarAquiferBaseflow,          & ! intent(out): baseflow from the aquifer (m s-1)
                   scalarAquiferBaseflowDeriv,     & ! intent(out): derivative in baseflow w.r.t. aquifer storage (s-1)
                   ! output: error control
                   err,cmessage)                     ! output: error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! get copies of baseflow flux to compute derivatives
   if(deriv_desired .and. ixDerivMethod==numerical)then
    select case(itry)
     case(unperturbed);       scalarFlux             = scalarAquiferBaseflow
     case(perturbStateAbove); scalarFlux_dStateAbove = scalarAquiferBaseflow
     case default; err=10; message=trim(message)//"unknown perturbation"; return
    end select
   endif

  end do  ! (looping through different flux calculations -- one or multiple calls depending if desire for numerical or analytical derivatives)

  ! compute numerical derivatives
  ! NOTE: baseflow derivatives w.r.t. state above are *actually* w.r.t. water table depth, so need to be corrected for aquifer storage
  !       (note also the negative sign to account for the inverse relationship between water table depth and aquifer storage)
  if(deriv_desired .and. ixDerivMethod==numerical)then
   scalarAquiferBaseflowDeriv = -(scalarFlux_dStateAbove - scalarFlux)/dx/drainablePorosity    ! change in baseflow flux w.r.t. change in state in aquifer (s-1)
  endif

 

 
  ! *************************************************************************************************************************************************
  ! *************************************************************************************************************************************************

  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! * compute water ejected because the volumetric liquid water content is close to porosity...
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! loop through soil layers
  do iLayer=1,nLevels

   ! compute the flux of water ejected because the volumetric liquid water content is close to porosity
   call ejectWater(&
                   ! input: model control
                   deriv_desired,                & ! intent(in): flag indicating if derivatives are desired
                   ixDerivMethod,                & ! intent(in): index defining the method used to compute derivatives (analytical or numerical)
                   ixRichards,                   & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                   ! input: state variables
                   mLayerMatricHeadTrial(iLayer),& ! intent(in): matric head in each layer (m)
                   mLayerVolFracLiqTrial(iLayer),& ! intent(in): volumetric liquid water content (-)
                   mLayerVolFracIceTrial(iLayer),& ! intent(in): volumetric ice content in each layer (-)
                   mLayerdTheta_dPsi(iLayer),    & ! intent(in): derivative in the soil water characteristic w.r.t. psi (m-1)
                   ! input: soil parameters
                   vGn_alpha,                    & ! intent(in): van Genutchen "alpha" parameter (m-1)
                   vGn_n,                        & ! intent(in): van Genutchen "n" parameter (-)
                   VGn_m,                        & ! intent(in): van Genutchen "m" parameter (-)
                   theta_sat,                    & ! intent(in): soil porosity (-)
                   theta_res,                    & ! intent(in): soil residual volumetric water content (-)
                   ! saturated hydraulic conductivity in each layer
                   mLayerSatHydCond(iLayer),     & ! intent(in):  saturated hydraulic conductivity in each layer (m s-1)
                   ! output: ejected water flux and derivative
                   mLayerEjectWater(iLayer),     & ! intent(out): water ejected because pore volume is filled (m s-1)
                   mLayerEjectWaterDeriv(iLayer),& ! intent(out): derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
                   ! output: error control
                   err,cmessage)                   ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  end do  ! looping through soil layers


  ! *****
  ! * compute the recharge flux
  scalarAquiferRcharge = iLayerLiqFluxSoil(nLevels) !+ sum(mLayerEjectWater(1:nLevels))

  end subroutine computeFlux



  ! ************************************************************************************************
  ! internal subroutine: compute the vector of fluxes at layer interfaces
  ! ************************************************************************************************
  subroutine cmpJacobian(&
                         ! input: model control
                         ixRichards,                   & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                         bc_lower,                     & ! intent(in): index defining the type of boundary conditions
                         ! input: state variables
                         mLayerMatricHeadInput,        & ! intent(in): matric head in each layer (m)
                         mLayerVolFracLiqInput,        & ! intent(in): volumetric liquid water content (-)
                         mLayerVolFracIceInput,        & ! intent(in): volumetric ice content in each layer (-)
                         scalarAquiferStorageInput,    & ! intent(in): aquifer storage (m)
                         ! output: error control
                         err,message)                    ! intent(out): error control
  implicit none
  ! input: model control
  integer(i4b),intent(in)          :: ixRichards                  ! index defining the option for Richards' equation (moisture or mixdform)
  integer(i4b),intent(in)          :: bc_lower                    ! index defining the type of boundary conditions
  ! input: trial model state variables
  real(dp),intent(in)              :: mLayerMatricHeadInput(:)    ! input value of matric head in each layer  (m)
  real(dp),intent(in)              :: mLayerVolFracLiqInput(:)    ! input value of volumetric fraction of liquid water in each layer (-)
  real(dp),intent(in)              :: mLayerVolFracIceInput(:)    ! input value of volumetric fraction of ice in each layer (-)
  real(dp),intent(in)              :: scalarAquiferStorageInput   ! input value of aquifer storage (m)
  ! output: error control
  integer(i4b),intent(out)         :: err                         ! error code
  character(*),intent(out)         :: message                     ! error message
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! local: perturbed vectors used to construct the Jacobian
  real(dp),dimension(nLevels)      :: mLayerMatricHeadTrial       ! perturbed value of matric head in each layer  (m)
  real(dp),dimension(nLevels)      :: mLayerVolFracLiqTrial       ! perturbed value of volumetric fraction of liquid water in each layer (-)
  real(dp)                         :: scalarAquiferStorageTrial   ! perturbed value of aquifer storage (m)
  ! local: output from computeFlux, kept local to avoid over-writing variables computed elsewhere
  real(dp),dimension(nLevels)      :: local_mLayerdPsi_dTheta     ! derivative in the soil water characteristic 
  real(dp),dimension(nLevels)      :: local_mLayerdTheta_dPsi     ! derivative in the soil water characteristic
  real(dp),dimension(0:nLevels)    :: local_dq_dStateAbove        ! derivatives in the flux w.r.t. state variable in the layer above (m s-1 or s-1)
  real(dp),dimension(0:nLevels)    :: local_dq_dStateBelow        ! derivatives in the flux w.r.t. state variable in the layer below (m s-1 or s-1)
  real(dp),dimension(nLevels)      :: local_mLayerEjectWaterDeriv ! derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
  real(dp)                         :: local_scalarSurfaceRunoff   ! surface runoff (m s-1) 
  real(dp)                         :: local_scalarWaterTableDepth ! water table depth (m)
  real(dp)                         :: local_scalarBaseflowDeriv   ! derivative in the baseflow fluw w.r.t. aquifer storage (s-1)
  ! local: fluxes used to calculate residual
  real(dp),dimension(nLevels)      :: mLayerTempEjectWater        ! water ejected because pore volume is close to capacity (m s-1)
  real(dp),dimension(0:nLevels)    :: iLayerTempLiqFluxSoil       ! liquid fluxes at layer interfaces at the start of the time step (m s-1)
  real(dp)                         :: scalarTempAquiferRcharge    ! recharge flux (m s-1)
  real(dp)                         :: scalarTempAquiferBaseflow   ! baseflow flux (m s-1)
  ! Jacobian matrix (used for testing)
  integer(i4b),parameter           :: minLayer=1                  ! minimum layer to test/print
  integer(i4b),parameter           :: maxLayer=51                 ! maximum layer to test/print
  real(dp),parameter               :: eps=1.0e-8_dp               ! finite difference increment
  integer(i4b)                     :: ijac                        ! index of columns
  real(dp),dimension(nLevels)      :: fTest                       ! function test (residual vector)
  real(dp)                         :: fAquifer                    ! function test (aquifer)
  real(dp),dimension(nState,nState):: jmat                        ! Jacobian matrix
  real, dimension(3)               :: xtry,ytry                   ! base and perturbed vectors used to compute the Jacobian
  real, dimension(3)               :: tdiag                       ! tri-diagonal matrix
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message='cmpJacobian/'

  ! get copies of the state vector to perturb
  mLayerMatricHeadTrial     = mLayerMatricHeadInput
  mLayerVolFracLiqTrial     = mLayerVolFracLiqInput
  scalarAquiferStorageTrial = scalarAquiferStorageInput

  ! loop through desired layers
  do ijac=1,nState

   ! perturb states
   if(iJac > nLevels)then  ! (can ONLY happen if aquifer storage)
    scalarAquiferStorageTrial = scalarAquiferStorageInput + eps
   else
    select case(ixRichards)
     case(moisture); mLayerVolFracLiqTrial(iJac) = mLayerVolFracLiqInput(iJac) + eps
     case(mixdform); mLayerMatricHeadTrial(iJac) = mLayerMatricHeadInput(iJac) + eps
    end select
   endif

   ! *****
   ! compute the fluxes, given trial state values
   call computeFlux(&
                    ! input: model control variables
                    .false.,                     & ! intent(in): flag to indicate if derivatives are desired 
                    ! input: trial state variables  (NOTE: use vectors from the start of the step)
                    mLayerMatricHeadTrial,       & ! intent(in): matric head (m)
                    mLayerVolFracLiqTrial,       & ! intent(in): volumetric fraction of liquid water (-)
                    mLayerVolFracIceInput,       & ! intent(in): volumetric fraction of ice (-)
                    scalarAquiferStorageTrial,   & ! intent(in): aquifer storage at the start of the step (m)
                    ! output: derivative in the soil water characteristic
                    local_mLayerdPsi_dTheta,     & ! intent(out): derivative in the soil water characteristic
                    local_mLayerdTheta_dPsi,     & ! intent(out): derivative in the soil water characteristic
                    ! output: diagnostic variables
                    local_scalarSurfaceRunoff,   & ! intent(out): surface runoff (m s-1)
                    local_scalarWaterTableDepth, & ! intent(out): water table depth (m)                         
                    ! output: fluxes
                    iLayerTempLiqFluxSoil,       & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                    mLayerTempEjectWater,        & ! intent(out): water ejected because pore volume is close to capacity (m s-1)
                    scalarTempAquiferRcharge,    & ! intent(out): recharge flux (m s-1)
                    scalarTempAquiferBaseflow,   & ! intent(out): baseflow flux (m s-1)
                    ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                    local_dq_dStateAbove,        & ! intent(out): derivatives in the flux w.r.t. state variable in the layer above (m s-1 or s-1)
                    local_dq_dStateBelow,        & ! intent(out): derivatives in the flux w.r.t. state variable in the layer below (m s-1 or s-1)
                    local_mLayerEjectWaterDeriv, & ! intent(out): derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
                    local_scalarBaseflowDeriv,   & ! intent(out): derivative in the baseflow fluw w.r.t. aquifer storage (s-1)
                    ! output: error control
                    err,cmessage)               ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! *****
   ! compute the residual vector
   call liqResidual(&
                    ! control variables
                    dt,                        & ! intent(in): length of the time step (s)
                    wimplicit,                 & ! intent(in):weight assigned to the start-of-step
                    ! coordinate variables
                    mLayerDepth,               & ! intent(in): depth of each layer (m)
                    ! initial flux vectors (start of the time step)
                    iLayerInitLiqFluxSoil,     & ! intent(in): initial liquid flux at layer interfaces (m s-1)
                    mLayerInitEjectWater,      & ! intent(in): initial water ejected from each layer because volumetric liquid water is approaching porosity (m s-1)
                    mLayerInitTranspire,       & ! intent(in): initial transpiration from each layer -- from energy routine (m s-1)
                    ! trial flux vectors
                    iLayerTempLiqFluxSoil,     & ! intent(in): trial liquid flux at layer interfaces (m s-1)
                    mLayerTempEjectWater,      & ! intent(in): trial water ejected from each layer because volumetric liquid water is approaching porosity (m s-1)
                    mLayerTranspire,           & ! intent(in): trial transpiration from each layer -- from energy routine (m s-1)
                    ! initial state vectors (start of the time step)
                    mLayerVolFracLiq,          & ! intent(in): initial volumetric liquid water content (-)
                    mLayerVolFracIce,          & ! intent(in): initial volumetric ice content (-)
                    ! trial state vectors
                    mLayerVolFracLiqTrial,     & ! intent(in): trial volumetric liquid water content (-)
                    mLayerVolFracIceInput,     & ! intent(in): trial volumetric ice content (-)
                    ! intent(out): residual vector (-)
                    fTest,                     & ! intent(out): residual vector (-)
                    ! intent(out): error control
                    err,cmessage)                ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! *****
   ! compute the residual for the groundwater store
   call gw_residual(&
                    ! input: control variables
                    dt,                         & ! intent(in): length of the time step (s)
                    wimplicit,                  & ! intent(in): weight assigned to the start-of-step
                    ! input: model states
                    scalarAquiferStorage,       & ! intent(in): aquifer storage at the start of the step (m)
                    scalarAquiferStorageTrial,  & ! intent(in): trial value of aquifer storage (m)
                    ! input: start-of-step fluxes
                    scalarInitAquiferTranspire, & ! intent(in): transpiration from the aquifer (m s-1)
                    scalarInitAquiferRcharge,   & ! intent(in): recharge to the aquifer        (m s-1)
                    scalarInitAquiferBaseflow,  & ! intent(in): baseflow from the aquifer      (m s-1)
                    ! input: end-of-step fluxes
                    scalarAquiferTranspire,     & ! intent(in): transpiration from the aquifer (m s-1)
                    scalarTempAquiferRcharge,   & ! intent(in): recharge to the aquifer        (m s-1)
                    scalarTempAquiferBaseflow,  & ! intent(in): baseflow from the aquifer      (m s-1)
                    ! output: aquifer residual
                    fAquifer,                   & ! intent(out): aquifer residual (m)
                    ! output: error control
                    err,cmessage)                 ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! compute Jacobian
   if(bc_lower==groundwaterCouple)then
    jmat(:,ijac) = ( (/ftest(:),fAquifer/) - (/rVec(:),scalarAquiferResidual/) ) / eps
   else
    jmat(:,ijac) = (ftest(:) - rVec(:) ) / eps
   endif
   
   ! set the state back to the input value
   if(iJac > nLevels)then  ! (can ONLY happen if aquifer storage)
    scalarAquiferStorageTrial = scalarAquiferStorageInput
   else
    select case(ixRichards)
     case(moisture); mLayerVolFracLiqTrial(iJac) = mLayerVolFracLiqInput(iJac)
     case(mixdform); mLayerMatricHeadTrial(iJac) = mLayerMatricHeadInput(iJac)
    end select
   endif

  end do  ! looping through soil layers

  ! print the Jacobian
  do iJac=minLayer,maxLayer
   if(iJac==1)then;          write(*,'(2(i4,1x),2(a,1x,3(e20.10,1x)))') ixDerivMethod, iJac, 'testing Jacobian', (/valueMissing,      jmat(iJac,1:iJac+1)/), '--> tri-diag = ', valueMissing, diag(iJac), d_p1(iJac)
   elseif(iJac==nState)then; write(*,'(2(i4,1x),2(a,1x,3(e20.10,1x)))') ixDerivMethod, iJac, 'testing Jacobian', (/jmat(iJac,iJac-1:nState), valueMissing/), '--> tri-diag = ', d_m1(iJac-1), diag(iJac), valueMissing
   else;                     write(*,'(2(i4,1x),2(a,1x,3(e20.10,1x)))') ixDerivMethod, iJac, 'testing Jacobian', (/jmat(iJac,iJac-1:iJac+1)              /), '--> tri-diag = ', d_m1(iJac-1), diag(iJac), d_p1(iJac) 
   endif
  end do

  end subroutine cmpJacobian

 end subroutine soilHydrol_muster


 ! ***************************************************************************************************************************
 ! ***************************************************************************************************************************
 ! ***************************************************************************************************************************
 ! ***************************************************************************************************************************
 ! ***** PUBLIC SUBROUTINES **************************************************************************************************
 ! ***************************************************************************************************************************
 ! ***************************************************************************************************************************
 ! ***************************************************************************************************************************

 ! ---------------------------------------------------------------------------------------------------------------------------
 ! public subroutine: compute the water required to fill pore space at each layer interface (m)
 ! ---------------------------------------------------------------------------------------------------------------------------
 subroutine waterTablePosition(mLayerVolFracLiqTrial, & ! intent(in):  volumetric liquid water content in each soil layer (-)
                               mLayerVolFracIceTrial, & ! intent(in):  volumetric ice content in each soil layer (-)
                               aquiferStorageTrial,   & ! intent(in):  aquifer storage (m)
                               iLayerHeight,          & ! intent(in):  height of each interface (m)
                               mLayerDepth,           & ! intent(in):  depth of each soil layer (m)
                               theta_sat,             & ! intent(in):  soil porosity (-)
                               specificYield,         & ! intent(in):  fraction of water volume drained by gravity in an unconfined aquifer (-)
                               drainablePorosity,     & ! intent(out): drainable porosity (-)
                               waterTableDepth,       & ! intent(out): water table depth (m)
                               err,message)             ! intent(out): error control
 implicit none
 ! input
 real(dp),intent(in)        :: mLayerVolFracLiqTrial(:) ! volumetric liquid water content in each soil layer (-)
 real(dp),intent(in)        :: mLayerVolFracIceTrial(:) ! volumetric ice content in each soil layer (-)
 real(dp),intent(in)        :: aquiferStorageTrial      ! aquifer storage (m)
 real(dp),intent(in)        :: iLayerHeight(0:)         ! height of each interface (m)
 real(dp),intent(in)        :: mLayerDepth(:)           ! depth of each soil layer (m)
 real(dp),intent(in)        :: theta_sat                ! soil porosity (-)
 real(dp),intent(in)        :: specificYield            ! fraction of water volume drained by gravity in an unconfined aquifer (-)
 ! output
 real(dp),intent(out)       :: drainablePorosity        ! drainable porosity (-)
 real(dp),intent(out)       :: waterTableDepth          ! water table depth (m)
 integer(i4b),intent(out)   :: err                      ! error code
 character(*),intent(out)   :: message                  ! error message
 ! local variables
 integer(i4b)               :: iSoil                    ! index of soil layers
 integer(i4b)               :: nSoil                    ! number of soil layers
 real(dp)                   :: waterReq2FillPoreBot     ! water required to fill pore space up to the bottom of a given layer (m)
 real(dp)                   :: waterReq2FillPoreTop     ! water required to fill pore space up to the top of a given layer (m)
 real(dp)                   :: waterReq2FillPoreLay     ! water required to fill pore space in a given layer (m)
 real(dp)                   :: fracFilled               ! fraction layer filled with water (-)
 ! initialize error control 
 err=0; message='waterTablePosition/'
 ! get the number of soil layers
 nSoil = size(mLayerDepth)

 ! ***** water table below the soil column
 if(aquiferStorageTrial < 0._dp)then
  drainablePorosity = specificYield
  waterTableDepth   = iLayerHeight(nSoil) - aquiferStorageTrial/specificYield
 ! ***** water table within the soil column
 else
  ! initialize the water required to fill the pore space
  waterReq2FillPoreBot = 0._dp
  ! compute the water required to fill pore space at each layer interface (m)
  do iSoil=nSoil,1,-1  ! (loop through soil layers, starting at the bottom)
   drainablePorosity    = theta_sat - mLayerVolFracLiqTrial(iSoil) - mLayerVolFracIceTrial(iSoil)
   waterReq2FillPoreLay = drainablePorosity*mLayerDepth(iSoil)
   waterReq2FillPoreTop = waterReq2FillPoreBot + waterReq2FillPoreLay
   ! return if have not filled the pore space in the current layer
   if(aquiferStorageTrial < waterReq2FillPoreTop)then
    fracFilled      = (aquiferStorageTrial - waterReq2FillPoreBot) / waterReq2FillPoreLay
    waterTableDepth = iLayerHeight(iSoil) - fracFilled*mLayerDepth(iSoil)
    return
   endif
   ! check that the water table is below the soil surface
   if(iSoil==1)then; err=20; message=trim(message)//'entire soil column is saturated'; return; endif
   ! continuing to the next layer -- set bottom storage to the top storage
   waterReq2FillPoreBot = waterReq2FillPoreTop
  end do  ! (looping through soil layers)
 endif   ! (if water table is within the soil column)
 end subroutine waterTablePosition


 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ***** PRIVATE SUBROUTINES **********************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************

 
 ! ************************************************************************************************
 ! private subroutine: compute transmittance and derivatives for model nodes
 ! ************************************************************************************************
 subroutine diagv_node(&
                       ! input: model control
                       deriv_desired,         & ! intent(in): flag indicating if derivatives are desired
                       ixRichards,            & ! intent(in): index defining the option for Richards' equation (moisture or mixdform)
                       ! input: state variables
                       scalarMatricHeadTrial, & ! intent(in):  matric head in a given
                       scalarVolFracLiqTrial, & ! intent(in):  volumetric liquid water content in a given soil layer (-)
                       scalarVolFracIceTrial, & ! intent(in):  volumetric ice content in a given soil layer (-)
                       ! input: soil parameters
                       vGn_alpha,             & ! intent(in): van Genutchen "alpha" parameter (m-1)
                       vGn_n,                 & ! intent(in): van Genutchen "n" parameter (-)
                       VGn_m,                 & ! intent(in): van Genutchen "m" parameter (-)
                       theta_sat,             & ! intent(in): soil porosity (-)
                       theta_res,             & ! intent(in): soil residual volumetric water content (-)
                       f_impede,              & ! intent(in): ice impedence factor (-)
                       ! input: saturated hydraulic conductivity
                       scalarSatHydCond,      & ! intent(in): saturated hydraulic conductivity at the mid-point of a given layer (m s-1)
                       ! output: derivative in the soil water characteristic
                       scalardPsi_dTheta,     & ! derivative in the soil water characteristic
                       scalardTheta_dPsi,     & ! derivative in the soil water characteristic
                       ! output: transmittance
                       scalarHydCond,         & ! intent(out): hydraulic conductivity at layer mid-points (m s-1)
                       scalarDiffuse,         & ! intent(out): diffusivity at layer mid-points (m2 s-1)
                       iceImpedeFac,          & ! intent(out): ice impedence factor in each layer (-)
                       ! output: transmittance derivatives
                       dHydCond_dVolLiq,      & ! intent(out): derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
                       dDiffuse_dVolLiq,      & ! intent(out): derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
                       dHydCond_dMatric,      & ! intent(out): derivative in hydraulic conductivity w.r.t matric head (m s-1)
                       ! output: error control
                       err,message)             ! intent(out): error control
 USE soil_utils_module,only:iceImpede       ! compute the ice impedence factor
 USE soil_utils_module,only:hydCond_psi     ! compute hydraulic conductivity as a function of matric head
 USE soil_utils_module,only:hydCond_liq     ! compute hydraulic conductivity as a function of volumetric liquid water content
 USE soil_utils_module,only:dTheta_dPsi     ! compute derivative of the soil moisture characteristic w.r.t. psi (m-1)
 USE soil_utils_module,only:dPsi_dTheta     ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
 USE soil_utils_module,only:dPsi_dTheta2    ! compute derivative in dPsi_dTheta (m)
 USE soil_utils_module,only:dHydCond_dLiq   ! compute derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
 USE soil_utils_module,only:dHydCond_dPsi   ! compute derivative in hydraulic conductivity w.r.t. matric head
 ! compute hydraulic transmittance and derivatives for all layers
 implicit none
 ! input: model control
 logical(lgt),intent(in)       :: deriv_desired             ! flag indicating if derivatives are desired
 integer(i4b),intent(in)       :: ixRichards                ! index defining the option for Richards' equation (moisture or mixdform)
 ! input: state and diagnostic variables
 real(dp),intent(in)           :: scalarMatricHeadTrial     ! matric head in each layer (m)
 real(dp),intent(in)           :: scalarVolFracLiqTrial     ! volumetric fraction of liquid water in a given layer (-)
 real(dp),intent(in)           :: scalarVolFracIceTrial     ! volumetric fraction of ice in a given layer (-)
 ! input: soil parameters
 real(dp),intent(in)           :: vGn_alpha                 ! van Genutchen "alpha" parameter (m-1)
 real(dp),intent(in)           :: vGn_n                     ! van Genutchen "n" parameter (-)
 real(dp),intent(in)           :: vGn_m                     ! van Genutchen "m" parameter (-)
 real(dp),intent(in)           :: theta_sat                 ! soil porosity (-)
 real(dp),intent(in)           :: theta_res                 ! soil residual volumetric water content (-)
 real(dp),intent(in)           :: f_impede                  ! ice impedence factor (-)
 ! input: saturated hydraulic conductivity
 real(dp),intent(in)           :: scalarSatHydCond          ! saturated hydraulic conductivity at the mid-point of a given layer (m s-1)
 ! output: derivative in the soil water characteristic
 real(dp),intent(out)          :: scalardPsi_dTheta         ! derivative in the soil water characteristic
 real(dp),intent(out)          :: scalardTheta_dPsi         ! derivative in the soil water characteristic
 ! output: transmittance
 real(dp),intent(out)          :: scalarHydCond             ! hydraulic conductivity at layer mid-points (m s-1)
 real(dp),intent(out)          :: scalarDiffuse             ! diffusivity at layer mid-points (m2 s-1)
 real(dp),intent(out)          :: iceImpedeFac              ! ice impedence factor in each layer (-)
 ! output: transmittance derivatives
 real(dp),intent(out)          :: dHydCond_dVolLiq          ! derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
 real(dp),intent(out)          :: dDiffuse_dVolLiq          ! derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
 real(dp),intent(out)          :: dHydCond_dMatric          ! derivative in hydraulic conductivity w.r.t matric head (m s-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                       ! error code
 character(*),intent(out)      :: message                   ! error message
 ! local variables
 real(dp)                      :: dPsi_dTheta2a             ! derivative in dPsi_dTheta (analytical)
 ! initialize error control
 err=0; message="diagv_node/"

 ! ***** 
 ! compute the derivative in the soil water characteristic
 select case(ixRichards)
  case(moisture)
   scalardPsi_dTheta = dPsi_dTheta(scalarvolFracLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   scalardTheta_dPsi = valueMissing  ! (deliberately cause problems if this is ever used)
  case(mixdform)
   scalardTheta_dPsi = dTheta_dPsi(scalarMatricHeadTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   scalardPsi_dTheta = valueMissing  ! (deliberately cause problems if this is ever used)
  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
 end select


 ! *****
 ! compute hydraulic conductivity and its derivative in each soil layer

 ! compute the ice impedence factor (-)
 iceImpedeFac = iceImpede(scalarVolFracIceTrial,theta_res,theta_sat,f_impede)

 select case(ixRichards)
  ! ***** moisture-based form of Richards' equation
  case(moisture)
   ! compute the hydraulic conductivity (m s-1) and diffusivity (m2 s-1) for a given layer
   scalarHydCond = hydCond_liq(scalarVolFracLiqTrial,scalarSatHydCond,theta_res,theta_sat,vGn_m) * iceImpedeFac
   scalarDiffuse = scalardPsi_dTheta * scalarHydCond
   ! compute derivative in hydraulic conductivity (m s-1) and hydraulic diffusivity (m2 s-1)
   if(deriv_desired)then
    dHydCond_dVolLiq = dHydCond_dLiq(scalarVolFracLiqTrial,scalarSatHydCond,theta_res,theta_sat,vGn_m,iceImpedeFac,.true.)  ! [.true. = analytical] 
    dPsi_dTheta2a    = dPsi_dTheta2(scalarVolFracLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m,.true.)   ! [.true. = analytical] compute derivative in dPsi_dTheta (m)
    dDiffuse_dVolLiq = dHydCond_dVolLiq*scalardPsi_dTheta + scalarHydCond*dPsi_dTheta2a
    dHydCond_dMatric = valueMissing ! not used, so cause problems
   endif

  ! ***** mixed form of Richards' equation -- just compute hydraulic condictivity
  case(mixdform)
   ! compute the hydraulic conductivity (m s-1) and diffusivity (m2 s-1) for a given layer
   scalarHydCond = hydCond_psi(scalarMatricHeadTrial,scalarSatHydCond,vGn_alpha,vGn_n,vGn_m) * iceImpedeFac
   scalarDiffuse = valueMissing ! not used, so cause problems
   ! compute derivative in hydraulic conductivity (m s-1)
   if(deriv_desired)then
    dHydCond_dMatric = dHydCond_dPsi(scalarMatricHeadTrial,scalarSatHydCond,vGn_alpha,vGn_n,vGn_m,iceImpedeFac,.true.)  ! analytical
    dHydCond_dVolLiq = valueMissing ! not used, so cause problems
    dDiffuse_dVolLiq = valueMissing ! not used, so cause problems
   endif

  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return

 endselect

 ! if derivatives are not desired, then set values to missing
 if(.not.deriv_desired)then
  dHydCond_dVolLiq = valueMissing ! not used, so cause problems
  dDiffuse_dVolLiq = valueMissing ! not used, so cause problems
  dHydCond_dMatric = valueMissing ! not used, so cause problems
 endif

 end subroutine diagv_node





 ! ************************************************************************************************
 ! private subroutine: compute the surface flux and its derivative
 ! ************************************************************************************************
 subroutine surfaceFlx(&
                       ! input: model control
                       deriv_desired,             & ! intent(in): flag indicating if derivatives are desired
                       ixRichards,                & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                       bc_upper,                  & ! intent(in): index defining the type of boundary conditions (neumann or diriclet)
                       ! input: state variables
                       scalarMatricHead,          & ! intent(in): matric head in he upper-most soil layer (m)
                       scalarVolFracLiq,          & ! intent(in): volumetric liquid water content the upper-most soil layer (-)
                       scalarVolFracIce,          & ! intent(in): volumetric ice content in the upper-most soil layer (-)
                       ! input: depth of upper-most soil layer (m)
                       upperLayerDepth,           & ! intent(in): depth of upper-most soil layer (m)
                       ! input: boundary conditions
                       upperBoundHead,            & ! intent(in): upper boundary condition (m)
                       upperBoundTheta,           & ! intent(in): upper boundary condition (-)
                       ! input: flux at the upper boundary
                       scalarRainPlusMelt,        & ! intent(in): rain plus melt (m s-1)
                       ! input: derivative in soil water characteristix
                       scalardPsi_dTheta,         & ! intent(in): derivative of the soil moisture characteristic w.r.t. theta (m)
                       scalardTheta_dPsi,         & ! intent(in): derivative of the soil moisture characteristic w.r.t. psi (m-1)
                       ! input: transmittance
                       surfaceSatHydCond,         & ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
                       iceImpedeFac,              & ! intent(in): ice impedence factor in the upper-most soil layer (-)
                       ! input: soil parameters
                       vGn_alpha,                 & ! intent(in): van Genutchen "alpha" parameter (m-1)
                       vGn_n,                     & ! intent(in): van Genutchen "n" parameter (-)
                       VGn_m,                     & ! intent(in): van Genutchen "m" parameter (-)
                       theta_sat,                 & ! intent(in): soil porosity (-)
                       theta_res,                 & ! intent(in): soil residual volumetric water content (-)
                       bpar_VIC,                  & ! intent(in): b-parameter in the VIC surface runoff parameterization (-)
                       ! output: hydraulic conductivity and diffusivity at the surface
                       surfaceHydCond,            & ! intent(in): hydraulic conductivity at the surface (m s-1)
                       surfaceDiffuse,            & ! intent(in): hydraulic diffusivity at the surface (m2 s-1)
                       ! output: fluxes at layer interfaces and surface runoff
                       scalarSurfaceRunoff,       & ! intent(out): surface runoff (m s-1)
                       scalarSurfaceInfiltration, & ! intent(out): surface infiltration (m s-1)
                       ! output: deriavtives in surface infiltration w.r.t. volumetric liquid water (m s-1) and matric head (s-1) in the upper-most soil layer
                       dq_dState,                 & ! intent(out): derivative in surface infiltration w.r.t. state variable in the upper-most soil layer (m s-1 or s-1)
                       ! output: error control
                       err,message)                 ! intent(out): error control
 USE soil_utils_module,only:volFracLiq      ! compute volumetric fraction of liquid water as a function of matric head (-)
 USE soil_utils_module,only:hydCond_psi     ! compute hydraulic conductivity as a function of matric head (m s-1)
 USE soil_utils_module,only:hydCond_liq     ! compute hydraulic conductivity as a function of volumetric liquid water content (m s-1)
 USE soil_utils_module,only:dPsi_dTheta     ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
 ! compute infiltraton at the surface and its derivative w.r.t. mass in the upper soil layer
 implicit none
 ! -----------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 logical(lgt),intent(in)       :: deriv_desired             ! flag to indicate if derivatives are desired
 integer(i4b),intent(in)       :: bc_upper                  ! index defining the type of boundary conditions
 integer(i4b),intent(in)       :: ixRichards                ! index defining the option for Richards' equation (moisture or mixdform)
 ! input: state and diagnostic variables
 real(dp),intent(in)           :: scalarMatricHead          ! matric head in the upper-most soil layer (m)
 real(dp),intent(in)           :: scalarVolFracLiq          ! volumetric liquid water content in the upper-most soil layer (-)
 real(dp),intent(in)           :: scalarVolFracIce          ! volumetric ice content in the upper-most soil layer (-)
 ! input: depth of upper-most soil layer (m)
 real(dp),intent(in)           :: upperLayerDepth           ! depth of upper-most soil layer (m)
 ! input: diriclet boundary conditions
 real(dp),intent(in)           :: upperBoundHead            ! upper boundary condition for matric head (m)
 real(dp),intent(in)           :: upperBoundTheta           ! upper boundary condition for volumetric liquid water content (-)
 ! input: flux at the upper boundary
 real(dp),intent(in)           :: scalarRainPlusMelt        ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
 ! input: derivative in soil water characteristix
 real(dp),intent(in)           :: scalardPsi_dTheta         ! derivative of the soil moisture characteristic w.r.t. theta (m)
 real(dp),intent(in)           :: scalardTheta_dPsi         ! derivative of the soil moisture characteristic w.r.t. psi (m-1)
 ! input: transmittance
 real(dp),intent(in)           :: iceImpedeFac              ! ice impedence factor in the upper-most soil layer (-)
 real(dp),intent(in)           :: surfaceSatHydCond         ! saturated hydraulic conductivity at the surface (m s-1)
 ! input: soil parameters
 real(dp),intent(in)           :: vGn_alpha                 ! van Genutchen "alpha" parameter (m-1)
 real(dp),intent(in)           :: vGn_n                     ! van Genutchen "n" parameter (-)
 real(dp),intent(in)           :: vGn_m                     ! van Genutchen "m" parameter (-)
 real(dp),intent(in)           :: theta_sat                 ! soil porosity (-)
 real(dp),intent(in)           :: theta_res                 ! soil residual volumetric water content (-)
 real(dp),intent(in)           :: bpar_VIC                  ! b-parameter in the VIC surface runoff parameterization (-)
 ! -----------------------------------------------------------------------------------------------------------------------------
 ! output: hydraulic conductivity and diffusivity at the surface
 real(dp),intent(out)          :: surfaceHydCond            ! hydraulic conductivity (m s-1)
 real(dp),intent(out)          :: surfaceDiffuse            ! hydraulic diffusivity at the surface (m
 ! output: surface runoff and infiltration flux (m s-1)
 real(dp),intent(out)          :: scalarSurfaceRunoff       ! surface runoff (m s-1)
 real(dp),intent(out)          :: scalarSurfaceInfiltration ! surface infiltration (m s-1)
 ! output: deriavtives in surface infiltration w.r.t. volumetric liquid water (m s-1) and matric head (s-1) in the upper-most soil layer
 real(dp),intent(out)          :: dq_dState                 ! derivative in surface infiltration w.r.t. state variable in the upper-most soil layer (m s-1 or s-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                       ! error code
 character(*),intent(out)      :: message                   ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------
 ! local variables
 real(dp)                      :: cFlux                     ! capillary flux (m s-1)
 real(dp)                      :: vFracLiq                  ! volumetric fraction of liquid water
 real(dp)                      :: fracCap                   ! fraction of pore space filled with liquid water and ice (-)
 real(dp)                      :: fInfArea                  ! area of the landscape where water infiltrates (-)
 real(dp)                      :: fInfRate                  ! infiltration rate at the surface (m s-1)
 real(dp)                      :: dInfArea                  ! derivative in the infiltrating area w.r.t. matric head (m-1)
 real(dp)                      :: dInfRate                  ! derivative in the infiltration rate w.r.t matric head (s-1)
 ! initialize error control
 err=0; message="surfaceFlx/"

 ! *****
 ! compute the surface flux and its derivative
 select case(bc_upper)

  ! *****
  ! head condition
  case(prescribedHead)

   ! surface runoff iz zero for the head condition
   scalarSurfaceRunoff = 0._dp

   ! compute transmission and the capillary flux
   select case(ixRichards)  ! (form of Richards' equation)
    case(moisture)
     ! compute the hydraulic conductivity and diffusivity at the boundary
     surfaceHydCond = hydCond_liq(upperBoundTheta,surfaceSatHydCond,theta_res,theta_sat,vGn_m) * iceImpedeFac
     surfaceDiffuse = dPsi_dTheta(upperBoundTheta,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) * surfaceHydCond
     ! compute the capillary flux
     cflux = -surfaceDiffuse*(scalarVolFracLiq - upperBoundTheta) / (upperLayerDepth*0.5_dp)
    case(mixdform)
     ! compute the hydraulic conductivity and diffusivity at the boundary
     surfaceHydCond = hydCond_psi(upperBoundHead,surfaceSatHydCond,vGn_alpha,vGn_n,vGn_m) * iceImpedeFac
     surfaceDiffuse = valueMissing
     ! compute the capillary flux
     cflux = -surfaceHydCond*(scalarMatricHead - upperBoundHead) / (upperLayerDepth*0.5_dp)
    case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
   end select  ! (form of Richards' eqn)
   ! compute the total flux
   scalarSurfaceInfiltration = cflux + surfaceHydCond
   ! compute the derivative
   if(deriv_desired)then
    select case(ixRichards)  ! (form of Richards' equation)
     case(moisture); dq_dState = -surfaceDiffuse/(upperLayerDepth/2._dp)
     case(mixdform); dq_dState = -surfaceHydCond/(upperLayerDepth/2._dp)
     case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
    end select
   endif

  ! *****
  ! flux condition
  case(liquidFlux)

   ! hydraulic conductivity and diffusivity are not used
   surfaceHydCond = valueMissing
   surfaceDiffuse = valueMissing
   ! compute volumetric fraction of liquid water (-)
   select case(ixRichards)
    case(moisture); vFracLiq = scalarVolFracLiq
    case(mixdform); vFracLiq = volFracLiq(scalarMatricHead,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
    case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
   endselect
   ! compute the infiltrating area (-)
   fracCap  = min((vFracLiq + scalarVolFracIce)/theta_sat, 1._dp - epsilon(fracCap))          ! fraction of capacity filled with liquid water and ice
   fInfArea = (1._dp - fracCap)**bpar_VIC
   ! compute the rate of infiltration over the non-saturated area (m s-1)
   fInfRate = min(surfaceSatHydCond,scalarRainPlusMelt)*iceImpedeFac  ! (m s-1)
   ! compute surface runoff (m s-1)
   scalarSurfaceRunoff = (1._dp - fInfArea)*scalarRainPlusMelt + fInfArea*(scalarRainPlusMelt - fInfRate)
   ! compute the flux at the upper boundary
   scalarSurfaceInfiltration = scalarRainPlusMelt - scalarSurfaceRunoff

   ! compute analytical derivatives (product rule)
   if(deriv_desired)then
    dInfRate = 0._dp  ! NOTE: infiltration rate does not depend on state variables in the upper-most layer
    select case(ixRichards)
     case(moisture); dInfArea = (-1._dp/theta_sat) * bpar_VIC*(1._dp - fracCap)**(bpar_VIC - 1._dp) ! derivative in the infiltrating area w.r.t. volumetric liquid water content (-)
     case(mixdform); dInfArea = (-scalardTheta_dPsi/theta_sat) * bpar_VIC*(1._dp - fracCap)**(bpar_VIC - 1._dp) ! derivative in the infiltrating area w.r.t. matric head (m-1)
     case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
    end select
    dq_dState = dInfArea*fInfRate !+ fInfArea*dInfRate ! NOTE, only first term is used currently, but "future-proofing"
   else
    dq_dState = valueMissing
   endif

  ! ***** error check
  case default; err=20; message=trim(message)//'unknown upper boundary condition for soil hydrology'; return

 endselect  ! (type of upper boundary condition)

 end subroutine surfaceFlx




 ! ************************************************************************************************
 ! private subroutine: compute the fluxes and derivatives at layer interfaces
 ! ************************************************************************************************
 subroutine iLayerFlux(&
                       ! input: model control
                       deriv_desired,             & ! intent(in): flag indicating if derivatives are desired
                       ixRichards,                & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                       ! input: state variables (adjacent layers)
                       nodeMatricHeadTrial,       & ! intent(in): matric head at the soil nodes (m)
                       nodeVolFracLiqTrial,       & ! intent(in): volumetric liquid water content at the soil nodes (-)
                       ! input: model coordinate variables (adjacent layers)
                       nodeHeight,                & ! intent(in): height of the soil nodes (m)
                       ! input: transmittance (adjacent layers)
                       nodeHydCondTrial,          & ! intent(in): hydraulic conductivity at the soil nodes (m s-1)
                       nodeDiffuseTrial,          & ! intent(in): hydraulic diffusivity at the soil nodes (m2 s-1)
                       ! input: transmittance derivatives (adjacent layers)
                       dHydCond_dVolLiq,          & ! intent(in): change in hydraulic conductivity w.r.t. change in volumetric liquid water content (m s-1)
                       dDiffuse_dVolLiq,          & ! intent(in): change in hydraulic diffusivity w.r.t. change in volumetric liquid water content (m2 s-1)
                       dHydCond_dMatric,          & ! intent(in): change in hydraulic conductivity w.r.t. change in matric head (s-1)
                       ! output: tranmsmittance at the layer interface (scalars)
                       iLayerHydCond,             & ! intent(out): hydraulic conductivity at the interface between layers (m s-1)
                       iLayerDiffuse,             & ! intent(out): hydraulic diffusivity at the interface between layers (m2 s-1)
                       ! output: vertical flux at the layer interface (scalars)
                       iLayerLiqFluxSoil,         & ! intent(out): vertical flux of liquid water at the layer interface (m s-1)
                       ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                       dq_dStateAbove,            & ! intent(out): derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer above (m s-1 or s-1)
                       dq_dStateBelow,            & ! intent(out): derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer below (m s-1 or s-1) 
                       ! output: error control
                       err,message)             ! intent(out): error control
 ! input: model control
 logical(lgt),intent(in)       :: deriv_desired               ! flag indicating if derivatives are desired
 integer(i4b),intent(in)       :: ixRichards                  ! index defining the option for Richards' equation (moisture or mixdform)
 ! input: state variables
 real(dp),intent(in)           :: nodeMatricHeadTrial(:)      ! matric head at the soil nodes (m)
 real(dp),intent(in)           :: nodeVolFracLiqTrial(:)      ! volumetric fraction of liquid water at the soil nodes (-)
 ! input: model coordinate variables
 real(dp),intent(in)           :: nodeHeight(:)               ! height at the mid-point of the lower layer (m)
 ! input: transmittance
 real(dp),intent(in)           :: nodeHydCondTrial(:)         ! hydraulic conductivity at layer mid-points (m s-1)
 real(dp),intent(in)           :: nodeDiffuseTrial(:)         ! diffusivity at layer mid-points (m2 s-1)
 ! input: transmittance derivatives
 real(dp),intent(in)           :: dHydCond_dVolLiq(:)         ! derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
 real(dp),intent(in)           :: dDiffuse_dVolLiq(:)         ! derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
 real(dp),intent(in)           :: dHydCond_dMatric(:)         ! derivative in hydraulic conductivity w.r.t matric head (m s-1)
 ! output: tranmsmittance at the layer interface (scalars)
 real(dp),intent(out)          :: iLayerHydCond               ! hydraulic conductivity at the interface between layers (m s-1)
 real(dp),intent(out)          :: iLayerDiffuse               ! hydraulic diffusivity at the interface between layers (m2 s-1)
 ! output: vertical flux at the layer interface (scalars)
 real(dp),intent(out)          :: iLayerLiqFluxSoil           ! vertical flux of liquid water at the layer interface (m s-1) 
 ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
 real(dp),intent(out)          :: dq_dStateAbove              ! derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer above (m s-1 or s-1)
 real(dp),intent(out)          :: dq_dStateBelow              ! derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer below (m s-1 or s-1) 
 ! output: error control
 integer(i4b),intent(out)      :: err                         ! error code
 character(*),intent(out)      :: message                     ! error message
 ! local variables (named variables to provide index of 2-element vectors)
 integer(i4b),parameter        :: ixUpper=1                   ! index of upper node in the 2-element vectors
 integer(i4b),parameter        :: ixLower=2                   ! index of lower node in the 2-element vectors
 ! local variables (Darcy flux)
 real(dp)                      :: dPsi                        ! spatial difference in matric head (m)
 real(dp)                      :: dLiq                        ! spatial difference in volumetric liquid water (-)
 real(dp)                      :: dz                          ! spatial difference in layer mid-points (m)
 real(dp)                      :: cflux                       ! capillary flux (m s-1)
 ! local variiables (derivative in Darcy's flux)
 real(dp)                      :: dHydCondIface_dVolLiqAbove  ! derivative in hydraulic conductivity at layer interface w.r.t. volumetric liquid water content in layer above
 real(dp)                      :: dHydCondIface_dVolLiqBelow  ! derivative in hydraulic conductivity at layer interface w.r.t. volumetric liquid water content in layer below
 real(dp)                      :: dDiffuseIface_dVolLiqAbove  ! derivative in hydraulic diffusivity at layer interface w.r.t. volumetric liquid water content in layer above
 real(dp)                      :: dDiffuseIface_dVolLiqBelow  ! derivative in hydraulic diffusivity at layer interface w.r.t. volumetric liquid water content in layer below
 real(dp)                      :: dHydCondIface_dMatricAbove  ! derivative in hydraulic conductivity at layer interface w.r.t. matric head in layer above
 real(dp)                      :: dHydCondIface_dMatricBelow  ! derivative in hydraulic conductivity at layer interface w.r.t. matric head in layer below
 ! initialize error control
 err=0; message="iLayerFlux/"

 ! ** compute the flux
 ! compute the hydraulic conductivity at the interface  -- use the geometric mean
 iLayerHydCond = (nodeHydCondTrial(ixLower) * nodeHydCondTrial(ixUpper))**0.5_dp
 ! compute the height difference between nodes
 dz = nodeHeight(ixLower) - nodeHeight(ixUpper)
 ! compute the capillary flux
 select case(ixRichards)  ! (form of Richards' equation)
  case(moisture)
   iLayerDiffuse = (nodeDiffuseTrial(ixLower) * nodeDiffuseTrial(ixUpper))**0.5_dp
   dLiq          = nodeVolFracLiqTrial(ixLower) - nodeVolFracLiqTrial(ixUpper)
   cflux         = -iLayerDiffuse * dLiq/dz
  case(mixdform)
   iLayerDiffuse = valueMissing
   dPsi          = nodeMatricHeadTrial(ixLower) - nodeMatricHeadTrial(ixUpper)
   cflux         = -iLayerHydCond * dPsi/dz
  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
 end select
 ! compute the total flux (add gravity flux, positive downwards)
 iLayerLiqFluxSoil = cflux + iLayerHydCond

 ! ** compute the derivatives
 if(deriv_desired)then
  select case(ixRichards)  ! (form of Richards' equation)
   case(moisture)
    ! derivatives in hydraulic conductivity at the layer interface (m s-1)
    dHydCondIface_dVolLiqAbove = dHydCond_dVolLiq(ixUpper)*nodeHydCondTrial(ixLower) * 0.5_dp/iLayerHydCond
    dHydCondIface_dVolLiqBelow = dHydCond_dVolLiq(ixLower)*nodeHydCondTrial(ixUpper) * 0.5_dp/iLayerHydCond
    ! derivatives in hydraulic diffusivity at the layer interface (m2 s-1)
    dDiffuseIface_dVolLiqAbove = dDiffuse_dVolLiq(ixUpper)*nodeDiffuseTrial(ixLower) * 0.5_dp/iLayerDiffuse
    dDiffuseIface_dVolLiqBelow = dDiffuse_dVolLiq(ixLower)*nodeDiffuseTrial(ixUpper) * 0.5_dp/iLayerDiffuse
    ! derivatives in the flux w.r.t. volumetric liquid water content
    dq_dStateAbove = -dDiffuseIface_dVolLiqAbove*dLiq/dz + iLayerDiffuse/dz + dHydCondIface_dVolLiqAbove
    dq_dStateBelow = -dDiffuseIface_dVolLiqBelow*dLiq/dz - iLayerDiffuse/dz + dHydCondIface_dVolLiqBelow
   case(mixdform)
    ! derivatives in hydraulic conductivity
    dHydCondIface_dMatricAbove = dHydCond_dMatric(ixUpper)*nodeHydCondTrial(ixLower) * 0.5_dp/iLayerHydCond
    dHydCondIface_dMatricBelow = dHydCond_dMatric(ixLower)*nodeHydCondTrial(ixUpper) * 0.5_dp/iLayerHydCond
    ! derivatives in the flux w.r.t. matric head
    dq_dStateAbove = -dHydCondIface_dMatricAbove*dPsi/dz + iLayerHydCond/dz + dHydCondIface_dMatricAbove
    dq_dStateBelow = -dHydCondIface_dMatricBelow*dPsi/dz - iLayerHydCond/dz + dHydCondIface_dMatricBelow
   case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
  end select
 else
  dq_dStateAbove = valueMissing
  dq_dStateBelow = valueMissing
 endif

 end subroutine iLayerFlux




 ! ************************************************************************************************
 ! private subroutine: compute the recharge flux and its derivative
 ! ************************************************************************************************
 subroutine qChargeFlx(&
                       ! input: model control
                       deriv_desired,             & ! intent(in): flag indicating if derivatives are desired
                       ixRichards,                & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                       hc_profile,                & ! intent(in): index defining the decrease of hydraulic conductivity with depth
                       bc_lower,                  & ! intent(in): index defining the type of boundary conditions
                       ! input: state variables
                       nodeMatricHead,            & ! intent(in): matric head in the lowest unsaturated node (m)
                       nodeVolFracLiq,            & ! intent(in): volumetric liquid water content the lowest unsaturated node (-)
                       nodeVolFracIce,            & ! intent(in): volumetric ice content in the lowest unsaturated node (-)
                       ! input: model coordinate variables
                       nodeDepth,                 & ! intent(in): depth of the lowest unsaturated soil layer (m)
                       nodeHeight,                & ! intent(in): height of the lowest unsaturated soil node (m)
                       unSatDepth,                & ! intent(in): effective depth of the unsaturated zone (m)
                       ! input: boundary conditions
                       lowerBoundHead,            & ! intent(in): lower boundary condition (m)
                       lowerBoundTheta,           & ! intent(in): lower boundary condition (-)
                       ! input: water table depth
                       scalarWaterTableDepth,     & ! intent(in): depth to the water table (m)
                       ! input: derivative in soil water characteristix
                       node__dPsi_dTheta,         & ! intent(in): derivative of the soil moisture characteristic w.r.t. theta (m)
                       node__dTheta_dPsi,         & ! intent(in): derivative of the soil moisture characteristic w.r.t. psi (m-1)
                       ! input: transmittance
                       surfaceSatHydCond,         & ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
                       bottomSatHydCond,          & ! intent(in): saturated hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
                       nodeHydCond,               & ! intent(in): hydraulic conductivity at the node itself (m s-1)
                       dHydCond_dVolLiq,          & ! intent(in): derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
                       dHydCond_dMatric,          & ! intent(in): derivative in hydraulic conductivity w.r.t. matric head (s-1)
                       iceImpedeFac,              & ! intent(in): ice impedence factor in the lower-most soil layer (-)
                       ! input: soil parameters
                       vGn_alpha,                 & ! intent(in): van Genutchen "alpha" parameter (m-1)
                       vGn_n,                     & ! intent(in): van Genutchen "n" parameter (-)
                       VGn_m,                     & ! intent(in): van Genutchen "m" parameter (-)
                       theta_sat,                 & ! intent(in): soil porosity (-)
                       theta_res,                 & ! intent(in): soil residual volumetric water content (-)
                       kAnisotropic,              & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                       zScale_TOPMODEL,           & ! intent(in): TOPMODEL scaling factor (m)
                       drainablePorosity,         & ! intent(in): drainable porosity (-)
                       ! output: hydraulic conductivity and diffusivity at the surface
                       bottomHydCond,             & ! intent(out): hydraulic conductivity at the bottom of the unsatuarted zone (m s-1)
                       bottomDiffuse,             & ! intent(out): hydraulic diffusivity at the bottom of the unsatuarted zone (m2 s-1)
                       ! output: recharge flux
                       scalarRecharge,            & ! intent(out): recharge flux (m s-1)
                       ! output: derivatives in recharge flux
                       dq_dStateUnsat,            & ! intent(out): change in recharge flux w.r.t. change in state variable in lowest unsaturated node (m s-1 or s-1)
                       dq_dAquifer,               & ! intent(out): change in recharge flux w.r.t. change in the aquifer storage (s-1)
                       ! output: error control
                       err,message)                 ! intent(out): error control
 USE soil_utils_module,only:volFracLiq      ! compute volumetric fraction of liquid water as a function of matric head (-)
 USE soil_utils_module,only:matricHead      ! compute matric head as a function of volumetric fraction of liquid water (m)
 USE soil_utils_module,only:hydCond_psi     ! compute hydraulic conductivity as a function of matric head (m s-1)
 USE soil_utils_module,only:hydCond_liq     ! compute hydraulic conductivity as a function of volumetric liquid water content (m s-1)
 USE soil_utils_module,only:dPsi_dTheta     ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
 ! compute infiltraton at the surface and its derivative w.r.t. mass in the upper soil layer
 implicit none
 ! -----------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 logical(lgt),intent(in)       :: deriv_desired             ! flag to indicate if derivatives are desired
 integer(i4b),intent(in)       :: ixRichards                ! index defining the option for Richards' equation (moisture or mixdform)
 integer(i4b),intent(in)       :: hc_profile                ! index defining the decrease of hydraulic conductivity with depth
 integer(i4b),intent(in)       :: bc_lower                  ! index defining the type of boundary conditions
 ! input: state and diagnostic variables
 real(dp),intent(in)           :: nodeMatricHead            ! matric head in the lowest unsaturated node (m)
 real(dp),intent(in)           :: nodeVolFracLiq            ! volumetric liquid water content in the lowest unsaturated node (-)
 real(dp),intent(in)           :: nodeVolFracIce            ! volumetric ice content in the lowest unsaturated node (-)
 ! input: model coordinate variables
 real(dp),intent(in)           :: nodeDepth                 ! depth of the lowest unsaturated soil layer (m)
 real(dp),intent(in)           :: nodeHeight                ! height of the lowest unsaturated soil node (m)
 real(dp),intent(in)           :: unSatDepth                ! effective depth of the unsaturated zone (m)
 ! input: diriclet boundary conditions
 real(dp),intent(in)           :: lowerBoundHead            ! lower boundary condition for matric head (m)
 real(dp),intent(in)           :: lowerBoundTheta           ! lower boundary condition for volumetric liquid water content (-)
 ! input: water table depth
 real(dp),intent(in)           :: scalarWaterTableDepth     ! depth to the water table (m)
 ! input: derivative in soil water characteristix
 real(dp),intent(in)           :: node__dPsi_dTheta         ! derivative of the soil moisture characteristic w.r.t. theta (m)
 real(dp),intent(in)           :: node__dTheta_dPsi         ! derivative of the soil moisture characteristic w.r.t. psi (m-1)
 ! input: transmittance
 real(dp),intent(in)           :: surfaceSatHydCond         ! saturated hydraulic conductivity at the surface (m s-1)
 real(dp),intent(in)           :: bottomSatHydCond          ! saturated hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
 real(dp),intent(in)           :: nodeHydCond               ! hydraulic conductivity at the node itself (m s-1)
 real(dp),intent(in)           :: dHydCond_dVolLiq          ! derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
 real(dp),intent(in)           :: dHydCond_dMatric          ! derivative in hydraulic conductivity w.r.t. matric head (s-1)
 real(dp),intent(in)           :: iceImpedeFac              ! ice impedence factor in the upper-most soil layer (-)
 ! input: soil parameters
 real(dp),intent(in)           :: vGn_alpha                 ! van Genutchen "alpha" parameter (m-1)
 real(dp),intent(in)           :: vGn_n                     ! van Genutchen "n" parameter (-)
 real(dp),intent(in)           :: vGn_m                     ! van Genutchen "m" parameter (-)
 real(dp),intent(in)           :: theta_sat                 ! soil porosity (-)
 real(dp),intent(in)           :: theta_res                 ! soil residual volumetric water content (-)
 real(dp),intent(in)           :: kAnisotropic              ! anisotropy factor for lateral hydraulic conductivity (-)
 real(dp),intent(in)           :: zScale_TOPMODEL           ! scale factor for TOPMODEL-ish baseflow parameterization (m)
 real(dp),intent(in)           :: drainablePorosity         ! drainable porosity (-)
 ! -----------------------------------------------------------------------------------------------------------------------------
 ! output: hydraulic conductivity at the bottom of the unsaturated zone
 real(dp),intent(out)          :: bottomHydCond             ! hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
 real(dp),intent(out)          :: bottomDiffuse             ! hydraulic diffusivity at the bottom of the unsatuarted zone (m2 s-1)
 ! output: recharge flux
 real(dp),intent(out)          :: scalarRecharge            ! recharge flux (m s-1)
 ! output: derivatives in recharge flux
 real(dp),intent(out)          :: dq_dStateUnsat            ! change in recharge flux w.r.t. change in state variable in lowest unsaturated node (m s-1 or s-1)
 real(dp),intent(out)          :: dq_dAquifer               ! change in recharge flux w.r.t. change in the aquifer storage (s-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                       ! error code
 character(*),intent(out)      :: message                   ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------
 ! local variables: compute flux when coupled to groundwater
 real(dp)                      :: zWater                    ! effective water table depth (m)
 real(dp)                      :: hc_depth                  ! hydraulic conductivity at the depth of the water table (m s-1)
 real(dp)                      :: dHeight                   ! height difference between unsaturated node and water table (m) 
 real(dp)                      :: nodePsi                   ! matric head in the lowest unsaturated node (m)
 real(dp)                      :: spGrad                    ! spatial gradient in matric head (-)
 real(dp)                      :: cflux                     ! capillary flux (m s-1)
 ! local variables: compute derivative when coupled to groundwater
 real(dp)                      :: dSpGrad_dTheta            ! change in spatial gradient w.r.t. change in theta (-) 
 real(dp)                      :: dhci_dTheta               ! change in hydraulic conductivity at the "interface" w.r.t. change in theta (m s-1
 real(dp)                      :: dCflux_dTheta             ! change in the capillary flux w.r.t. change in theta (m s-1)
 real(dp)                      :: dhc_dzwt                  ! change in hyd cond w.r.t. change in water table depth (s-1)
 real(dp)                      :: dhci_dhcz                 ! change in hydraulic conductivity at the "interface" w.r.t. change in the hydraulic conductivity at the water table depth (-)
 real(dp)                      :: dhci_dzwt                 ! change in hydraulic conductivity at the "interface" w.r.t. change in water table depth (s-1)
 real(dp)                      :: dSpGrad_dzwt              ! change in the spatial derivative w.r.t change in water table depth (m-1)
 real(dp)                      :: dCflux_dzwt               ! change in capillary flux w.r.t. change in the water table depth -- product rule (s-1)
 ! -----------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="qChargeFlx/"

 ! determine lower boundary condition
 select case(bc_lower)

  ! ---------------------------------------------------------------------------------------------
  ! * coupled to groundwater
  ! ---------------------------------------------------------------------------------------------
  case(groundwaterCouple)

   ! error check
   if(ixRichards/=moisture)then; err=20; message=trim(message)//"require moisture-based form of Richards' eqn when using coupled gw"; return; endif
   if(scalarWaterTableDepth < nodeHeight)then; err=20; message=trim(message)//'water table is above the mid-point of the lowest soil layer'; return; endif

   ! compute the hydraulic conductivity at the water table depth (m s-1), and change in hyd cond w.r.t. change in water table depth (s-1)
   select case(hc_profile)  ! (determine type of hydraulic conductivity profile)
    case(constant)
     hc_depth = surfaceSatHydCond      ! constant, so hyd cond at the depth of the water table = surfaceSatHydCond (m s-1)
     if(deriv_desired) dhc_dzwt = 0._dp  ! change in hyd cond w.r.t. change in water table depth (s-1)
    case(exp_profile)
     hc_depth = surfaceSatHydCond * exp(-scalarWaterTableDepth/zScale_TOPMODEL)  ! hyd cond at the depth of the water table (m s-1)
     if(deriv_desired) dhc_dzwt = -hc_depth/zScale_TOPMODEL ! change in hyd cond w.r.t. change in water table depth (s-1)
    case(powerLaw_profile,linear_profile)
     message=trim(message)//'hydraulic conductivity profile not implemented yet for "powerLaw" and "linear" options'; err=10; return
    case default
     message=trim(message)//'unknown hydraulic conductivity profile'; err=10; return
   end select
   ! get the hydraulic conductivity at the mid-point between (1) the lowest unsaturated node and (2) the water table depth
   bottomHydCond = (nodeHydCond * hc_depth)**0.5_dp ! NOTE: this is at height = 0.5*(nodeHeight + scalarWaterTableDepth)
   ! compute the capillary flux
   dHeight = scalarWaterTableDepth - nodeHeight
   nodePsi = matricHead(nodeVolFracLiq,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   spGrad  = -(0._dp - nodePsi) / dHeight    ! spatial gradient (-)
   cflux   = bottomHydCond*spGrad            ! capillary flux (m s-1)
   ! compute the total flux
   scalarRecharge = cflux + bottomHydCond

   ! compute derivatives
   if(deriv_desired)then

    ! ** compute change in recharge flux w.r.t. change in volumetric liquid water content in the layer above
    ! compute change in spatial gradient w.r.t. change in theta (-)
    dSpGrad_dTheta = node__dPsi_dTheta/dHeight
    ! compute the change in hydraulic conductivity at the "interface" w.r.t. change in theta (m s-1)
    dhci_dTheta    = (hc_depth/(bottomHydCond*2._dp)) * dHydCond_dVolLiq
    ! compute change in the capillary flux w.r.t. change in theta (m s-1)
    dCflux_dTheta  = dhci_dTheta*spGrad + bottomHydCond*dSpGrad_dTheta
    ! compute the final derivative (m s-1)
    dq_dStateUnsat = dCflux_dTheta + dhci_dTheta

    ! ** compute change in recharge flux w.r.t. change in the aquifer storage
    ! compute the change in hydraulic conductivity at the "interface" w.r.t. change in the hydraulic conductivity at the water table depth (-)
    dhci_dhcz    = nodeHydCond/(bottomHydCond*2._dp)
    ! compute the change in hydraulic conductivity at the "interface" w.r.t. change in water table depth -- chain rule (s-1)
    dhci_dzwt    = dhc_dzwt*dhci_dhcz
    ! compute the change in the spatial derivative w.r.t change in water table depth (m-1)
    dSpGrad_dzwt = -nodePsi/(dHeight**2._dp)
    ! define derivative in capillary flux w.r.t. water table depth -- product rule (s-1)
    dCflux_dzwt  = dhci_dzwt*spGrad + bottomHydCond*dSpGrad_dzwt
    ! compute final derivative (s-1)
    ! NOTE: negative sign used because water table depth is inversely related to aquifer storage
    dq_dAquifer = -(dCflux_dzwt + dhci_dzwt)/drainablePorosity

   else     ! (do not desire derivatives)
    dq_dStateUnsat = valueMissing
    dq_dAquifer    = valueMissing
   endif    ! (if desire derivatives)


  ! ---------------------------------------------------------------------------------------------
  ! * prescribed head
  ! ---------------------------------------------------------------------------------------------
  case(prescribedHead)

   ! compute fluxes
   select case(ixRichards)  ! (moisture-based form of Richards' equation)
    case(moisture)
     ! compute the hydraulic conductivity and diffusivity at the boundary
     bottomHydCond = hydCond_liq(lowerBoundTheta,bottomSatHydCond,theta_res,theta_sat,vGn_m) * iceImpedeFac
     bottomDiffuse = dPsi_dTheta(lowerBoundTheta,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) * bottomHydCond
     ! compute the capillary flux
     cflux = -bottomDiffuse*(lowerBoundTheta - nodeVolFracLiq) / (nodeDepth*0.5_dp)
    case(mixdform)
     ! compute the hydraulic conductivity and diffusivity at the boundary     
     bottomHydCond = hydCond_psi(lowerBoundHead,bottomSatHydCond,vGn_alpha,vGn_n,vGn_m) * iceImpedeFac
     bottomDiffuse = valueMissing
     ! compute the capillary flux
     cflux = -bottomHydCond*(lowerBoundHead  - nodeMatricHead) / (nodeDepth*0.5_dp)
    case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
   end select  ! (form of Richards' eqn)
   scalarRecharge = cflux + bottomHydCond

   ! compute derivatives
   if(deriv_desired)then
    select case(ixRichards)  ! (form of Richards' equation)
     case(moisture); dq_dStateUnsat = bottomDiffuse/(nodeDepth/2._dp)
     case(mixdform); dq_dStateUnsat = bottomHydCond/(nodeDepth/2._dp)
     case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
    end select
   else     ! (do not desire derivatives)
    dq_dStateUnsat = valueMissing
   endif

   ! no aquifer, so set aquifer derivative to missing
   dq_dAquifer = valueMissing

  ! ---------------------------------------------------------------------------------------------
  ! * function of matric head in the bottom layer
  ! ---------------------------------------------------------------------------------------------
  case(funcBottomHead)

   ! compute fluxes
   select case(ixRichards)
    case(moisture); nodePsi = matricHead(nodeVolFracLiq,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
    case(mixdform); nodePsi = nodeMatricHead
   endselect
   zWater = nodeHeight - nodePsi
   scalarRecharge = kAnisotropic*surfaceSatHydCond * exp(-zWater/zScale_TOPMODEL)

   ! compute derivatives
   if(deriv_desired)then
    select case(ixRichards)  ! (form of Richards' equation)
     case(moisture); dq_dStateUnsat = kAnisotropic*surfaceSatHydCond * node__dPsi_dTheta*exp(-zWater/zScale_TOPMODEL)/zScale_TOPMODEL
     case(mixdform); dq_dStateUnsat = kAnisotropic*surfaceSatHydCond * exp(-zWater/zScale_TOPMODEL)/zScale_TOPMODEL
     case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
    end select
   else     ! (do not desire derivatives)
    dq_dStateUnsat = valueMissing
   endif

   ! no aquifer, so set aquifer derivative to missing
   dq_dAquifer = valueMissing

  ! ---------------------------------------------------------------------------------------------
  ! * free drainage
  ! ---------------------------------------------------------------------------------------------
  case(freeDrainage)

   ! compute flux
   scalarRecharge = nodeHydCond

   ! compute derivatives
   if(deriv_desired)then
    select case(ixRichards)  ! (form of Richards' equation)
     case(moisture); dq_dStateUnsat = dHydCond_dVolLiq
     case(mixdform); dq_dStateUnsat = dHydCond_dMatric
     case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
    end select
   else     ! (do not desire derivatives)
    dq_dStateUnsat = valueMissing
   endif

   ! no aquifer, so set aquifer derivative to missing
   dq_dAquifer = valueMissing


  ! ---------------------------------------------------------------------------------------------
  ! * error check
  ! ---------------------------------------------------------------------------------------------
  case default; err=20; message=trim(message)//'unknown lower boundary condition for soil hydrology'; return

 endselect ! (type of boundary condition)

 end subroutine qChargeFlx




 ! ************************************************************************************************
 ! private subroutine: compute baseflow flux and its derivative
 ! ************************************************************************************************
 subroutine qBflowFlux(&
                       ! input: model control
                       deriv_desired,                  & ! intent(in): flag indicating if derivatives are desired
                       ! input: "states"
                       scalarWaterTableDepth,          & ! output: water table depth at the start/end of the time step (m)
                       ! input: diagnostic variables
                       k_surf,                         & ! input: hydraulic conductivity at the surface (m s-1)
                       drainablePorosity,              & ! input: drainablePorosity (-)
                       ! input: parameters
                       kAnisotropic,                   & ! input: anisotropy factor for lateral hydraulic conductivity (-)
                       zScale_TOPMODEL,                & ! input: scale factor for TOPMODEL-ish baseflow parameterization (m)
                       ! output: baseflow flux and its derivative
                       scalarAquiferBaseflow,          & ! output: baseflow from the aquifer (m s-1)
                       scalarAquiferBaseflowDeriv,     & ! output: derivative in baseflow w.r.t. aquifer storage (s-1)
                       ! output: error control
                       err,message)                      ! output: error control
 ! compute change in aquifer storage over the time step
 implicit none
 ! input: model control
 logical(lgt),intent(in)   :: deriv_desired              ! flag to indicate if derivatives are desired
 ! input: "states"
 real(dp),intent(in)       :: scalarWaterTableDepth      ! water table depth at the start/end of the time step (m)
 ! input: diagnostic variables
 real(dp),intent(in)       :: k_surf                     ! hydraulic conductivity at the surface (m s-1) 
 real(dp),intent(in)       :: drainablePorosity          ! drainable porosity (-)
 ! input: parameters
 real(dp),intent(in)       :: kAnisotropic               ! anisotropy factor for lateral hydraulic conductivity (-)
 real(dp),intent(in)       :: zScale_TOPMODEL            ! scale factor for TOPMODEL-ish baseflow parameterization (m)
 ! output: baseflow flux and its derivative
 real(dp),intent(out)      :: scalarAquiferBaseflow      ! baseflow from the aquifer (m s-1)
 real(dp),intent(out)      :: scalarAquiferBaseflowDeriv ! derivative in baseflow w.r.t. aquifer storage (s-1)
 ! output: error control
 integer(i4b),intent(out)  :: err                        ! error code
 character(*),intent(out)  :: message                    ! error message
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="qBflowFlux/"
 ! calculate the baseflow term (m s-1)
 scalarAquiferBaseflow = kAnisotropic*k_surf*exp(-scalarWaterTableDepth/zScale_TOPMODEL)
 ! compute the derivative in baseflow w.r.t. storage (s-1)
 ! NOTE: negative sign because of the inverse relationship between water table depth and aquifer storage
 if(deriv_desired)then
  scalarAquiferBaseflowDeriv = -scalarAquiferBaseflow/(zScale_TOPMODEL*drainablePorosity)
 else
  scalarAquiferBaseflowDeriv = valueMissing
 endif
 end subroutine qBflowFlux




 ! ************************************************************************************************
 ! private subroutine: compute flux of water ejected because close to exceeding porosity
 ! ************************************************************************************************
 subroutine ejectWater(&
                       ! input: model control
                       deriv_desired,         & ! intent(in): flag determining if the derivative is desired
                       dMethod,               & ! intent(in): index defining the method used to compute derivatives (analytical or numerical)
                       ixRichards,            & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                       ! input: state variable and diagnostic variables
                       scalarMatricHeadTrial, & ! intent(in): matric head in each layer (m)
                       scalarVolFracLiqTrial, & ! intent(in): volumetric liquid water content (-)
                       scalarVolFracIceTrial, & ! intent(in): volumetric ice content in each layer (-)
                       scalardTheta_dPsi,     & ! intent(in): derivative in the soil water characteristic w.r.t. psi (m-1)
                       ! input: soil parameters
                       vGn_alpha,             & ! intent(in): van Genutchen "alpha" parameter (m-1)
                       vGn_n,                 & ! intent(in): van Genutchen "n" parameter (-)
                       VGn_m,                 & ! intent(in): van Genutchen "m" parameter (-)
                       theta_sat,             & ! intent(in): soil porosity (-)
                       theta_res,             & ! intent(in): soil residual volumetric water content (-)
                       ! input: saturated hydraulic conductivity
                       scalarSatHydCond,      & ! intent(in): saturated hydraulic conductivity in each layer (m s-1)
                       ! output: ejected water flux and derivative
                       scalarEjectWater,      & ! intent(out): water ejected because pore volume is filled (m s-1)
                       scalarEjectWaterDeriv, & ! intent(out): derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
                       ! output: error control
                       err,message)             ! intent(out): error control
 USE soil_utils_module,only:volFracLiq          ! compute volumetric fraction of liquid water (-)
 ! avoid super-saturation, by computing flux of water due to displacement (m s-1)
 implicit none
 ! input: model control
 logical(lgt),intent(in)       :: deriv_desired             ! flag determining if the derivative is desired
 integer(i4b),intent(in)       :: dMethod                   ! index defining the method used to compute derivatives (analytical or numerical)
 integer(i4b),intent(in)       :: ixRichards                ! index defining the form of Richards' equation (moisture or mixdform)
 ! input: state variables
 real(dp),intent(in)           :: scalarMatricHeadTrial     ! matric head in each layer (m)
 real(dp),intent(in)           :: scalarVolFracLiqTrial     ! volumetric liquid water content in each layer (-)
 real(dp),intent(in)           :: scalarVolFracIceTrial     ! volumetric ice content in each layer (-)
 real(dp),intent(in)           :: scalardTheta_dPsi         ! derivative in the soil water characteristic w.r.t. psi (m-1)
 ! input: soil parameters
 real(dp),intent(in)           :: vGn_alpha                 ! van Genutchen "alpha" parameter (m-1)
 real(dp),intent(in)           :: vGn_n                     ! van Genutchen "n" parameter (-)
 real(dp),intent(in)           :: vGn_m                     ! van Genutchen "m" parameter (-)
 real(dp),intent(in)           :: theta_sat                 ! soil porosity (-)
 real(dp),intent(in)           :: theta_res                 ! soil residual volumetric water content (-)
 ! input: saturated hydraulic conductivity
 real(dp),intent(in)           :: scalarSatHydCond          ! saturated hydraulic conductivity in each layer (m s-1)
 ! output: ejected water flux and derivative
 real(dp),intent(out)          :: scalarEjectWater          ! water ejected because pore volume is filled (m s-1)
 real(dp),intent(out)          :: scalarEjectWaterDeriv     ! derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
 ! output: error control
 integer(i4b),intent(out)      :: err                       ! error code
 character(*),intent(out)      :: message                   ! error message
 ! local variables
 real(dp),parameter            :: supersatScale=0.001_dp    ! scaling factor for the water ejection function (-)
 real(dp),parameter            :: xMatch = 0.999_dp         ! point where x-value and function value match (-)
 real(dp),parameter            :: fSmall = epsilon(xMatch)  ! smallest possible value to test
 real(dp)                      :: supersatThresh            ! threshold in super-saturation function (-)
 real(dp)                      :: fracMin                   ! minimum fraction of pore space required to be filled in order for lateral flow to occur (-)
 real(dp)                      :: fracCap                   ! fraction of pore space filled with liquid water and ice (-)
 real(dp)                      :: expFunc                   ! exponential function used as part of the flux calculation (-)
 real(dp)                      :: expTemp                   ! exponential function used as part of the flux calculation (-)
 ! initialize error control
 err=0; message="ejectWater/"

 ! define threshold in the super-saturation function (-)
 supersatThresh = supersatScale * log(1._dp/xMatch - 1._dp) + xMatch

 ! define minimum value for calculations
 fracMin = -supersatScale*log(1._dp/fSmall - 1._dp) + supersatThresh

 ! only process if moisture-based form of RE, or if ice is present in the mixed form of RE
 if(ixRichards==moisture .or. (ixRichards==mixdform .and. scalarVolFracIceTrial>0._dp))then
  ! calculate the fraction of pore space filled with liquid water and ice (-)
  select case(ixRichards)
   case(moisture); fracCap = (scalarVolFracLiqTrial + scalarVolFracIceTrial)/theta_sat
   case(mixdform); fracCap = (volFracLiq(scalarMatricHeadTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) + scalarVolFracIceTrial)/theta_sat
   case default; err=20; message=trim(message)//"unable to identify option for Richards' equation"; return
  endselect
  ! check if the fractional capacity is greater than the minimum value
  if(fracCap > fracMin)then
   ! calculate the lateral flux (m s-1)
   expFunc = exp((supersatThresh - fracCap)/supersatScale)
   scalarEjectWater = scalarSatHydCond/(1._dp + expFunc)
   ! calculate the derivative (m s-1 [moisture form] or s-1 [mixed form])
   if(deriv_desired)then
    select case(ixRichards)
     case(moisture)
      if(dMethod==analytical)then
       scalarEjectWaterDeriv = (expFunc/(theta_sat*supersatScale)) * (1._dp + expFunc)**(-2._dp) * scalarSatHydCond
      else
       fracCap = (scalarVolFracLiqTrial+dx + scalarVolFracIceTrial)/theta_sat
       expTemp = exp((supersatThresh - fracCap)/supersatScale)
       scalarEjectWaterDeriv = (scalarSatHydCond/(1._dp + expTemp) -  scalarEjectWater)/dx
      endif
     case(mixdform)
      if(dMethod==analytical)then
       scalarEjectWaterDeriv = scalardTheta_dPsi * (expFunc/(theta_sat*supersatScale)) * (1._dp + expFunc)**(-2._dp) * scalarSatHydCond
      else
       fracCap = (volFracLiq(scalarMatricHeadTrial+dx,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) + scalarVolFracIceTrial)/theta_sat
       expTemp  = exp((supersatThresh - fracCap)/supersatScale)
       scalarEjectWaterDeriv = (scalarSatHydCond/(1._dp + expTemp) -  scalarEjectWater)/dx
      endif
     case default; err=20; message=trim(message)//"unable to identify option for Richards' equation"; return
    end select
   else
    scalarEjectWaterDeriv = valueMissing  ! (derivatives not desired)
   endif
  ! (fraction of pore space filled with liquid water and ice is less than the minimum value required for lateral flow)
  else
   scalarEjectWater      = 0._dp
   scalarEjectWaterDeriv = 0._dp
  endif
 ! (mixed form of Richards' equation when no ice is present)
 else
  scalarEjectWater      = 0._dp
  scalarEjectWaterDeriv = 0._dp
 endif
 end subroutine ejectWater




 ! ************************************************************************************************
 ! private subroutine: compute flux of water ejected because close to exceeding porosity
 ! ************************************************************************************************
 subroutine liqResidual(&
                        ! control variables
                        dt,                              & ! length of the time step (s)
                        wimplicit,                       & ! weight assigned to the start-of-step
                        ! coordinate variables
                        mLayerDepth,                     & ! depth of each layer (m)
                        ! initial flux vectors (start of the time step)
                        iLayerInitLiqFluxSoil,           & ! intent(in): initial liquid flux at layer interfaces (m s-1)
                        mLayerInitEjectWater,            & ! intent(in): initial water ejected from each layer because volumetric liquid water is approaching porosity (m s-1)
                        mLayerInitTranspire,             & ! intent(in): initial transpiration from each layer -- from energy routine (m s-1)
                        ! trial flux vectors
                        iLayerTrialLiqFluxSoil,          & ! intent(in): trial liquid flux at layer interfaces (m s-1)
                        mLayerTrialEjectWater,           & ! intent(in): trial water ejected from each layer because volumetric liquid water is approaching porosity (m s-1)
                        mLayerTrialTranspire,            & ! intent(in): trial transpiration from each layer -- from energy routine (m s-1)
                        ! initial state vectors (start of the time step)
                        mLayerInitVolFracLiq,            & ! intent(in): initial volumetric liquid water content (-)
                        mLayerInitVolFracIce,            & ! intent(in): initial volumetric ice content (-)
                        ! trial state vectors
                        mLayerTrialVolFracLiq,           & ! intent(in): trial volumetric liquid water content (-)
                        mLayerTrialVolFracIce,           & ! intent(in): trial volumetric ice content (-)
                        ! intent(out): residual vector (-)
                        residualVec,                     & ! intent(out): residual vector (-)
                        ! intent(out): error control
                        err,message)                       ! intent(out): error control
 implicit none
 ! control variables
 real(dp),intent(in)          :: dt                        ! length of the time step (s)
 real(dp),intent(in)          :: wimplicit                 ! weight assigned to the start-of-step
 ! coordinate variables
 real(dp),intent(in)          :: mLayerDepth(:)            ! depth of each layer (m)
 ! initial flux vectors (start of the time step)
 real(dp),intent(in)          :: iLayerInitLiqFluxSoil(0:) ! intent(in): initial liquid flux at layer interfaces (m s-1)
 real(dp),intent(in)          :: mLayerInitEjectWater(:)   ! intent(in): initial water ejected from each layer because volumetric liquid water is approaching porosity (m s-1)
 real(dp),intent(in)          :: mLayerInitTranspire(:)    ! intent(in): initial transpiration from each layer -- from energy routine (m s-1)
 ! trial flux vectors
 real(dp),intent(in)          :: iLayerTrialLiqFluxSoil(0:) ! intent(in): trial liquid flux at layer interfaces (m s-1)
 real(dp),intent(in)          :: mLayerTrialEjectWater(:)  ! intent(in): trial water ejected from each layer because volumetric liquid water is approaching porosity (m s-1)
 real(dp),intent(in)          :: mLayerTrialTranspire(:)   ! intent(in): trial transpiration from each layer -- from energy routine (m s-1)
 ! initial state vectors (start of the time step)
 real(dp),intent(in)          :: mLayerInitVolFracLiq(:)   ! intent(in): initial volumetric liquid water content (-)
 real(dp),intent(in)          :: mLayerInitVolFracIce(:)   ! intent(in): initial volumetric ice content (-)
 ! trial state vectors
 real(dp),intent(in)          :: mLayerTrialVolFracLiq(:)  ! intent(in): trial volumetric liquid water content (-)
 real(dp),intent(in)          :: mLayerTrialVolFracIce(:)  ! intent(in): trial volumetric ice content (-)
 ! output
 real(dp), intent(out)        :: residualVec(:)            ! residual vector
 integer(i4b),intent(out)     :: err                       ! error code
 character(*),intent(out)     :: message                   ! error message
 ! local variables
 real(dp)                     :: dt_dz                     ! time/depth (s m-1)
 real(dp)                     :: flux0,flux1               ! contribution to liquid water from fluxes at start-of-step and end of step (-) 
 real(dp)                     :: mFlux                     ! overall contribution to volumetric liquid water from fluxes (-)
 real(dp)                     :: mEvap                     ! overall contribution to volumetric liquid water from transpiration (-) 
 real(dp)                     :: mEjct                     ! overall contribution to volumetric liquid water from water ejected from the layer (-) 
 real(dp)                     :: mPhse                     ! overall contribution to volumetric liquid water from phase change (-) 
 integer(i4b)                 :: iLayer                    ! index of soil layer
 ! initialize error control
 err=0; message='liqResidual/'

 ! *****
 ! compute the residual vector (-)
 do iLayer=1,nLevels
  ! time/depth
  dt_dz = dt/mLayerDepth(iLayer)
  ! fluxes (-)
  flux0 = -(iLayerInitLiqFluxSoil(iLayer)  - iLayerInitLiqFluxSoil(iLayer-1))*dt_dz
  flux1 = -(iLayerTrialLiqFluxSoil(iLayer) - iLayerTrialLiqFluxSoil(iLayer-1))*dt_dz
  mFlux = wimplicit*flux0 + (1._dp - wimplicit)*flux1
  ! transpiration (-)
  mEvap = wimplicit*mLayerInitTranspire(iLayer)*dt_dz + (1._dp - wimplicit)*mLayerTrialTranspire(iLayer)*dt_dz
  ! ejected water (-)
  mEjct = wimplicit*mLayerInitEjectWater(iLayer)*dt_dz + (1._dp - wimplicit)*mLayerTrialEjectWater(iLayer)*dt_dz
  ! phase change (-)
  mPhse = (iden_ice/iden_water)*(mLayerTrialVolFracIce(iLayer)-mLayerInitVolFracIce(iLayer))
  ! residual (-)
  residualVec(iLayer) = mLayerTrialVolFracLiq(iLayer) - (mLayerInitVolFracLiq(iLayer) + mFlux + mEvap - mEjct - mPhse)
  ! print progress
  !if(iLayer==1)   print*, 'iLayer, residualVec(iLayer), bottom flux, mFlux, mEvap, mEjct, mPhse'
  !if(iLayer > 45) write(*,'(i4,1x,10(e20.10,1x))') iLayer, residualVec(iLayer), iLayerTrialLiqFluxSoil(iLayer)*dt_dz, mFlux, mEvap, mEjct, mPhse
 end do  ! (looping through soil layers)

 end subroutine liqResidual




 ! ************************************************************************************************
 ! private subroutine: compute flux of water ejected because close to exceeding porosity
 ! ************************************************************************************************
 subroutine gw_residual(&
                        ! input: control variables
                        dt,                        & ! intent(in): length of the time step (s)
                        wimplicit,                 & ! intent(in): weight assigned to the start-of-step
                        ! input: model states
                        scalarAquiferStorage,      & ! intent(in): aquifer storage at the start of the step (m)
                        scalarAquiferStorageTrial, & ! intent(in): trial value of aquifer storage (m)
                        ! input: start-of-step fluxes
                        scalarInitAquiferTranspire,& ! intent(in): transpiration from the aquifer (m s-1)
                        scalarInitAquiferRcharge,  & ! intent(in): recharge to the aquifer        (m s-1)
                        scalarInitAquiferBaseflow, & ! intent(in): baseflow from the aquifer      (m s-1)
                        ! input: end-of-step fluxes
                        scalarAquiferTranspire,    & ! intent(in): transpiration from the aquifer (m s-1)
                        scalarAquiferRcharge,      & ! intent(in): recharge to the aquifer        (m s-1)
                        scalarAquiferBaseflow,     & ! intent(in): baseflow from the aquifer      (m s-1)
                        ! output: aquifer residual
                        scalarAquiferResidual,     & ! intent(out): aquifer residual (m)
                        ! output: error control
                        err,message)                 ! intent(out): error control
 ! ------------------------------------------------------------------------------------------------------------------------
 implicit none
 ! input: control variables
 real(dp),intent(in)          :: dt                          ! length of the time step (s)
 real(dp),intent(in)          :: wimplicit                   ! weight assigned to the start-of-step
 ! input: state variables
 real(dp),intent(in)          :: scalarAquiferStorage        ! aquifer storage at the start of the time step (m)
 real(dp),intent(in)          :: scalarAquiferStorageTrial   ! trial value of aquifer storage (m)
 ! input: start-of-step fluxes
 real(dp),intent(in)          :: scalarInitAquiferTranspire  ! aquifer transpiration averaged over the time step (m s-1)
 real(dp),intent(in)          :: scalarInitAquiferRcharge    ! aquifer recharge averaged over the time step (m s-1)
 real(dp),intent(in)          :: scalarInitAquiferBaseflow   ! baseflow from the aquifer (m s-1) 
 ! input: end-of-step fluxes
 real(dp),intent(in)          :: scalarAquiferTranspire      ! aquifer transpiration averaged over the time step (m s-1)
 real(dp),intent(in)          :: scalarAquiferRcharge        ! aquifer recharge averaged over the time step (m s-1)
 real(dp),intent(in)          :: scalarAquiferBaseflow       ! baseflow from the aquifer (m s-1) 
 ! output: aquifer residual
 real(dp),intent(out)         :: scalarAquiferResidual       ! aquifer residual (m)
 ! output: error control
 integer(i4b),intent(out)     :: err                         ! error code
 character(*),intent(out)     :: message                     ! error message
 ! ------------------------------------------------------------------------------------------------------------------------
 ! local variables
 real(dp)                     :: netFlux0                    ! net flux at the start of the step (m s-1)
 real(dp)                     :: netFlux1                    ! net flux at the end of the step (m s-1)
 real(dp)                     :: netFlux                     ! mean net flux over the step (m s-1)
 ! ------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="gw_residual/"
 ! compute the flux at the start of the step
 netFlux0 = (scalarInitAquiferRcharge + scalarInitAquiferTranspire) - scalarInitAquiferBaseflow
 ! compute the flux at the end of the step
 netFlux1 = (scalarAquiferRcharge + scalarAquiferTranspire) - scalarAquiferBaseflow
 ! compute the mean net flux
 netFlux  = wimplicit*netFlux0 + (1._dp - wimplicit)*netFlux1
 ! compute the residual
 !print*, 'transpire = ', scalarInitAquiferTranspire, scalarAquiferTranspire
 !print*, 'recharge = ',  scalarInitAquiferRcharge, scalarAquiferRcharge
 !print*, 'baseflow = ', scalarInitAquiferBaseflow, scalarAquiferBaseflow
 !print*, 'netFlux0, netFlux1, netFlux = ', netFlux0, netFlux1, netFlux
 scalarAquiferResidual = scalarAquiferStorageTrial - (scalarAquiferStorage + netFlux*dt)
 !print*, 'scalarAquiferResidual = ', scalarAquiferResidual
 end subroutine gw_residual




end module soilHydrol_module

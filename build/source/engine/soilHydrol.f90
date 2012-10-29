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
public::waterReq2FillLayer
! number of layers
integer(i4b)           :: nSoil                    ! number of soil layers
integer(i4b)           :: nSnow                    ! number of snow layers
integer(i4b)           :: nLevels                  ! total number of soil layers to examine
! missing value parameter
real(dp),parameter     :: valueMissing=-9999._dp   ! missing value parameter
! finite difference increment
real(dp),parameter     :: dx=1.e-8_dp              ! finite difference increment
contains


 ! ************************************************************************************************
 ! new subroutine: compute change in volumetric liquid water content (or matric head) over the time step
 ! ************************************************************************************************
 subroutine soilHydrol(dt,&                   ! time step (seconds)
                       iter,&                 ! iteration index
                       mLayerMatricHeadIter,& ! matric head in each layer at the current iteration (m)
                       mLayerVolFracIceIter,& ! volumetric fraction of ice at the current iteration (-)
                       mLayerVolFracLiqIter,& ! volumetric fraction of liquid water at the current iteration (-)
                       mLayerMatricHeadNew, & ! matric head in each layer at the next iteration (m)
                       mLayerVolFracLiqNew,& ! volumetric fraction of liquid water at the next iteration (-)
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
 real(dp),intent(in)           :: mLayerVolFracIceIter(:)  ! volumetric fraction of ice at the current iteration (-)
 real(dp),intent(in)           :: mLayerVolFracLiqIter(:)  ! volumetric fraction of liquid water at the current iteration (-)
 ! output
 real(dp),intent(out)          :: mLayerMatricHeadNew(:)   ! matric head in each layer at the next iteration (m)
 real(dp),intent(out)          :: mLayerVolFracLiqNew(:)   ! volumetric fraction of liquid water at the next iteration (m)
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
                        mLayerVolFracIceIter,                                         & ! intent(in): volumetric fraction of ice at the current iteration (-)
                        mLayerVolFracLiqIter,                                         & ! intent(in): volumetric fraction of liquid water at the current iteration (-)
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
                        ! height at the bottom of the effective soil profile and the depth to the water table
                        mvar_data%var(iLookMVAR%iLayerHeight)%dat(iend),              & ! intent(in): height at the bottom of the effective soil profile (m)
                        mvar_data%var(iLookMVAR%scalarWaterTableDepth)%dat(1),        & ! intent(in): (scalar) water table depth (m)
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
                        mvar_data%var(iLookMVAR%scalarSfcMeltPond)%dat(1),            & ! intent(in): (scalar) ponded water caused by melt of the "snow without a layer" (kg m-2)
                        ! saturated hydraulic conductivity
                        mvar_data%var(iLookMVAR%mLayerSatHydCond)%dat,                & ! intent(in): saturated hydraulic conductivity at the mid-point of each layer (m s-1)
                        mvar_data%var(iLookMVAR%iLayerSatHydCond)%dat,                & ! intent(in): saturated hydraulic conductivity at the interface of each layer (m s-1)
                        ! transpiration (from energy routines)
                        mvar_data%var(iLookMVAR%mLayerInitTranspire)%dat,             & ! intent(in): transpiration loss from each soil layer at start-of-step (m s-1)
                        mvar_data%var(iLookMVAR%mLayerTranspire)%dat,                 & ! intent(in): transpiration loss from each soil layer (m s-1)
                        ! diagnostic scalar variables
                        mvar_data%var(iLookMVAR%scalarRainPlusMelt)%dat(1),           & ! intent(out): (scalar) rain plus melt (m s-1)
                        mvar_data%var(iLookMVAR%scalarSurfaceRunoff)%dat(1),          & ! intent(out): (scalar) surface runoff (m s-1)
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
                              mLayerVolFracIceIter,       & ! intent(in): volumetric fraction of ice at the current iteration (-)
                              mLayerVolFracLiqIter,       & ! intent(in): volumetric fraction of liquid water at the current iteration (-)
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
                              ! height at the bottom of the effective soil profile and the depth to the water table
                              lowerBoundaryHeight,        & ! intent(in): height at the bottom of the effective soil profile
                              scalarWaterTableDepth,      & ! intent(in): (scalar)    ! water table depth (m)
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
                              scalarSfcMeltPond,          & ! intent(in): (scalar)    ! ponded water caused by melt of the "snow without a layer" (kg m-2)
                              ! saturated hydraulic conductivity in each layer
                              mLayerSatHydCond,           & ! intent(in): saturated hydraulic conductivity at the mid-point of each layer (m s-1)
                              iLayerSatHydCond,           & ! intent(in): saturated hydraulic conductivity at the interface of each layer (m s-1)
                              ! transpiration (from energy routines)
                              mLayerInitTranspire,        & ! intent(in): transpiration loss from each soil layer at start-of-step (m s-1)
                              mLayerTranspire,            & ! intent(in): transpiration loss from each soil layer (m s-1)
                              ! diagnostic scalar variables
                              scalarRainPlusMelt,         & ! intent(out): rain plus melt (m s-1)
                              scalarSurfaceRunoff,        & ! intent(out): surface runoff (m s-1)
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
 real(dp),intent(in)              :: mLayerVolFracIceIter(:)  ! volumetric fraction of ice at the current iteration (-)
 real(dp),intent(in)              :: mLayerVolFracLiqIter(:)  ! volumetric fraction of liquid water at the current iteration (-)
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
 ! height at the bottom of the effective soil profile and the depth to the water table
 real(dp),intent(in)              :: lowerBoundaryHeight      ! height at the bottom of the effective soil profile
 real(dp),intent(in)              :: scalarWaterTableDepth    ! water table depth (m)
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
 real(dp),intent(in)              :: scalarSfcMeltPond        ! ponded water caused by melt of the "snow without a layer" (kg m-2)
 ! saturated hydraulic conductivity
 real(dp),intent(in)              :: mLayerSatHydCond(:)      ! saturated hydraulic conductivity at the mid-point of each layer (m s-1)
 real(dp),intent(in)              :: iLayerSatHydCond(0:)     ! saturated hydraulic conductivity at the interface of each layer (m s-1)
 ! transpiration (from energy routines)
 real(dp),intent(in)              :: mLayerInitTranspire(:)   ! transpiration loss from each soil layer at the start of the time step (m s-1)
 real(dp),intent(in)              :: mLayerTranspire(:)       ! transpiration loss from each soil layer (m s-1)
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! ***** output variables
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! diagnostic scalar variables
 real(dp),intent(out)             :: scalarRainPlusMelt       ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
 real(dp),intent(out)             :: scalarSurfaceRunoff      ! surface runoff (m s-1)
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
 integer(i4b),intent(out)         :: err                      ! error code
 character(*),intent(out)         :: message                  ! error message
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! ***** local variables
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables (general)
 character(LEN=256)               :: cmessage                 ! error message of downwind routine
 logical(lgt)                     :: printflag                ! flag to print crap to the screen
 integer(i4b)                     :: iLayer                   ! layer index
 real(dp)                         :: availPorosity            ! pore space available to fill
 real(dp)                         :: availLiqWater            ! liquid water available to drain (negative) 
 ! hydraulic conductivity
 real(dp),dimension(nLevels)      :: iceImpedeFac             ! ice impedence factor at layer mid-points (-)
 real(dp),dimension(nLevels)      :: mLayerHydCond            ! hydraulic conductivity at layer mid-point (m s-1)
 real(dp),dimension(nLevels)      :: mLayerDiffuse            ! diffusivity at layer mid-point (m2 s-1)
 real(dp),dimension(0:nLevels)    :: iLayerHydCond            ! hydraulic conductivity at layer interface (m s-1)
 real(dp),dimension(0:nLevels)    :: iLayerDiffuse            ! diffusivity at layer interface (m2 s-1)
 ! flux of water ejected when volumetric liquid water content is close to exceeding porosity
 real(dp),dimension(nLevels)      :: mLayerEjectWaterDeriv    ! derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
 ! check need to compute initial fluxes
 integer(i4b)                     :: nFlux                    ! index for flux calculation
 integer(i4b)                     :: ifluxInit                ! starting index for flux calculations (0 or 1)
 real(dp)                         :: maxdiffMatric(1)         ! used to check if we are starting on the first iteration
 ! flux derivatives
 real(dp),dimension(0:nLevels)    :: dq_dMatricAbove          ! change in the flux in layer interfaces w.r.t. matric head in the layer above
 real(dp),dimension(0:nLevels)    :: dq_dMatricBelow          ! change in the flux in layer interfaces w.r.t. matric head in the layer below
 real(dp),dimension(0:nLevels)    :: dq_dVolLiqAbove          ! change in the flux in layer interfaces w.r.t. volumetric liquid water content in the layer above
 real(dp),dimension(0:nLevels)    :: dq_dVolLiqBelow          ! change in the flux in layer interfaces w.r.t. volumetric liquid water content in the layer below
 ! tri-diagonal solution
 real(dp)                         :: wtim                     ! weighted time (s)
 real(dp),dimension(nLevels-1)    :: d_m1                     ! sub-diagonal elements of the tridiagonal system (m-1)
 real(dp),dimension(nLevels)      :: diag                     ! diagonal elements of the tridiagonal system (m-1)
 real(dp),dimension(nLevels-1)    :: d_p1                     ! super-diagonal elements of the tridiagonal system (m-1)
 real(dp),dimension(nLevels)      :: rvec                     ! right-hand-side vector (-)
 real(dp),dimension(nLevels)      :: mLayerMatricHeadDiff     ! iteration increment for matric head (m)
 real(dp),dimension(nLevels)      :: mLayerVolFracLiqDiff     ! iteration increment for volumetric fraction of liquid water (m)
 ! Jacobian matrix (used for testing) 
 real(dp),dimension(nLevels)         :: mLayerTempEjectWater  ! water ejected because pore volume is close to capacity (m s-1)
 real(dp),dimension(0:nLevels)       :: iLayerTempLiqFluxSoil ! liquid fluxes at layer interfaces at the start of the time step (m s-1)
 real(dp),dimension(nLevels,nLevels) :: jmat                  ! jacobian matrix
 real(dp),dimension(nLevels)         :: ftest,xsav,xtry       ! compute jacobian: function vector and dependent variables
 real(dp),dimension(nLevels)         :: xph,h                 ! compute jacobian: perturbed vector and finite difference increment
 real(dp),parameter                  :: eps=-1.0e-8_dp        ! compute jacobian: finite difference increment
 integer(i4b)                        :: ijac                  ! compute jacobian: index of columns
 logical(lgt),parameter              :: calcJacobian=.false.   ! flag to compute the Jacobian matrix
 integer(i4b),parameter              :: minLayer=45           ! minimum layer to print
 integer(i4b),parameter              :: maxLayer=50           ! minimum layer to print
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
 !print*, 'scalarRainPlusMelt, rainfall/iden_water, (scalarSfcMeltPond/dt)/iden_water, scalarSfcMeltPond = ', &
 !         scalarRainPlusMelt, rainfall/iden_water, (scalarSfcMeltPond/dt)/iden_water, scalarSfcMeltPond

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

  ! ***** process initial conditions
  if(nFlux==0)then
   call computeFlux(&
                    ! input: trial state variables  (NOTE: use vectors from the start of the step)
                    mLayerMatricHead,       & ! intent(in): matric head (m)
                    mLayerVolFracIce,       & ! intent(in): volumetric fraction of ice (-)
                    mLayerVolFracLiq,       & ! intent(in): volumetric fraction of liquid water (-)
                    ! output: fluxes (NOTE: populate initial vectors)
                    mLayerInitEjectWater,   & ! intent(out): water ejected because pore volume is close to capacity (m s-1)
                    iLayerInitLiqFluxSoil,  & ! intent(out): liquid fluxes at layer interfaces at the start of the time step(m s-1)
                    ! ----
                    ! NOTE: we have no need for the optional arguments...
                    ! -------------------------------------------------
                    ! output: error control
                    err=err,message=cmessage) ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! ***** process trial vectors
  else
   call computeFlux(&
                    ! input: trial state variables  (NOTE: use trial vectors)
                    mLayerMatricHeadIter,   & ! intent(in): matric head (m)
                    mLayerVolFracIceIter,   & ! intent(in): volumetric fraction of ice (-)
                    mLayerVolFracLiqIter,   & ! intent(in): volumetric fraction of liquid water (-)
                    ! output: fluxes (NOTE: populate standard flux vectors)
                    mLayerEjectWater,       & ! intent(out): water ejected because pore volume is close to capacity (m s-1)
                    iLayerLiqFluxSoil,      & ! intent(out): liquid fluxes at layer interfaces at the start of the time step(m s-1)
                    ! ----
                    ! NOTE: we wish to return the optional arguments...
                    ! -------------------------------------------------
                    ! output: derivative in the soil water characteristic
                    mLayerdPsi_dTheta,      & ! derivative in the soil water characteristic
                    mLayerdTheta_dPsi,      & ! derivative in the soil water characteristic
                    ! output: transmittance
                    mLayerHydCond,          & ! intent(out): hydraulic conductivity at layer mid-points (m s-1)
                    mLayerDiffuse,          & ! intent(out): diffusivity at layer mid-points (m2 s-1)
                    iLayerHydCond,          & ! intent(out): hydraulic conductivity at layer interface (m s-1)
                    iLayerDiffuse,          & ! intent(out): diffusivity at layer interface (m2 s-1)
                    iceImpedeFac,           & ! intent(out): ice impedence factor in each layer (-)
                    ! output: derivative in ejected water flux
                    mLayerEjectWaterDeriv,  & ! intent(out): derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
                    ! output: surface runoff
                    scalarSurfaceRunoff,    & ! intent(out): surface runoff (m s-1)
                    ! output: error control
                    err,cmessage)             ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  endif

 end do  ! looping through initial vectors

 ! *****
 ! compute the residual vector
 call liqResidual(&
                  ! control variables
                  dt,                        & ! length of the time step (s)
                  wimplicit,                 & ! weight assigned to the start-of-step
                  ! coordinate variables
                  mLayerDepth,               & ! depth of each layer (m)
                  ! initial flux vectors (start of the time step)
                  iLayerInitLiqFluxSoil,     & ! intent(in): initial liquid flux at layer interfaces (m s-1)
                  mLayerInitEjectWater,      & ! intent(in): initial water ejected from each layer because volumetric liquid water is approaching porosity (m s-1)
                  mLayerInitTranspire,       & ! intent(in): initial transpiration from each layer -- from energy routine (m s-1)
                  ! trial flux vectors
                  iLayerLiqFluxSoil,         & ! intent(in): trial liquid flux at layer interfaces (m s-1)
                  mLayerEjectWater,          & ! intent(in): trial water ejected from each layer because volumetric liquid water is approaching porosity (m s-1)
                  mLayerTranspire,           & ! intent(in): trial transpiration from each layer -- from energy routine (m s-1)
                  ! initial state vectors (start of the time step)
                  mLayerVolFracLiq,          & ! intent(in): initial volumetric liquid water content (-)
                  mLayerVolFracIce,          & ! intent(in): initial volumetric ice content (-)
                  ! trial state vectors
                  mLayerVolFracLiqIter,      & ! intent(in): trial volumetric liquid water content (-)
                  mLayerVolFracIceIter,      & ! intent(in): trial volumetric ice content (-)
                  ! intent(out): residual vector (-)
                  rVec,                      & ! intent(out): residual vector (-)
                  ! intent(out): error control
                  err,cmessage)                ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! call groundwater





 ! *****
 ! compute derivative in flux
 select case(ixRichards)
  case(moisture)
   ! compute the derivative in flux at layer interfaces w.r.t. volumetric liquid water content in the layer above and in the layer below
   call dq_dVLiquid(&
                    ! input: model control
                    ixDerivMethod,         & ! intent(in): method used to compute derivatives (analytical or numerical)
                    ixBcUpperSoilHydrology,& ! index defining the type of upper boundary conditions 
                    ixBcLowerSoilHydrology,& ! index defining the type of upper boundary conditions 
                    ! input: model coordinate variables
                    mLayerDepth,           & ! intent(in): depth of the layer (m)
                    mLayerHeight,          & ! intent(in): height of the layer mid-point (m)
                    ! input: state and diagnostic variables
                    mLayerVolFracLiqIter,  & ! intent(in): volumetric liquid water content (-)
                    mLayerVolFracIceIter,  & ! intent(in): volumetric fraction of ice (-)
                    mLayerdPsi_dTheta,     & ! intent(in): derivative in the soil water characteritic w.r.t. theta (m)
                    ! input: height of the lower boundary and the depth to the water table
                    lowerBoundaryHeight,   & ! intent(in): height of the lower boundary
                    scalarWaterTableDepth, & ! intent(in): water table depth (m)
                    ! input: boundary conditions for volumetric liquid water content
                    lowerBoundTheta,       & ! intent(in): lower boundary condition (-)
                    upperBoundTheta,       & ! intent(in): upper boundary condition (-)
                    ! input: flux at the upper boundary
                    scalarRainPlusMelt,    & ! intent(in): rain plus melt (m s-1)
                    ! input: transmittance
                    mLayerSatHydCond,      & ! intent(in): saturated hydraulic conductivity in each layer (m s-1)
                    iLayerSatHydCond,      & ! intent(in): saturated hydraulic conductivity at layer interfaces (m s-1)
                    mLayerHydCond,         & ! intent(in): hydraulic conductivity in each layer (m s-1)
                    mLayerDiffuse,         & ! intent(in): hydraulic diffusivity in each layer (m s-1)
                    iLayerHydCond,         & ! intent(in): hydraulic conductivity at layer interfaces (m s-1)
                    iLayerDiffuse,         & ! intent(in): hydraulic diffusivity at layer interfaces (m s-1)
                    iceImpedeFac,          & ! intent(in): ice impedence factor in each layer (-)
                    ! input: soil parameters
                    vGn_alpha,             & ! intent(in): van Genutchen "alpha" parameter (m-1)
                    vGn_n,                 & ! intent(in): van Genutchen "n" parameter (-)
                    VGn_m,                 & ! intent(in): van Genutchen "m" parameter (-)
                    theta_sat,             & ! intent(in): soil porosity (-)
                    theta_res,             & ! intent(in): soil residual volumetric water content (-)
                    kAnisotropic,          & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                    zScale_TOPMODEL,       & ! intent(in): scale factor for TOPMODEL-ish baseflow parameterization (m)
                    bpar_VIC,              & ! intent(in): b-parameter in the VIC surface runoff parameterization (-)
                    f_impede,              & ! intent(in): ice impedence factor (-)
                    ! output
                    dq_dVolLiqAbove,       & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer above (m s-1)
                    dq_dVolLiqBelow,       & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer below (m s-1)
                    err,cmessage)            ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  case(mixdform)
   ! compute the derivative in flux at layer interfaces w.r.t. matric head in the layer above and in the layer below
   call dq_dMatHead(&
                    ! input: model control
                    ixDerivMethod,         & ! intent(in): method used to compute derivatives (analytical or numerical)
                    ixBcUpperSoilHydrology,& ! index defining the type of upper boundary conditions 
                    ixBcLowerSoilHydrology,& ! index defining the type of upper boundary conditions 
                    ! input: model coordinate variables
                    mLayerDepth,           & ! intent(in): depth of the layer (m)
                    mLayerHeight,          & ! intent(in): height of the layer mid-point (m)
                    ! input: state and diagnostic variables
                    mLayerMatricHeadIter,  & ! intent(in): matric head (m)
                    mLayerVolFracIceIter,  & ! intent(in): volumetric fraction of ice (-)
                    mLayerdTheta_dPsi,     & ! intent(in): derivative in the soil water characteritic w.r.t. psi (m-1)
                    ! input: boundary conditions for matric head
                    lowerBoundHead,        & ! intent(in): lower boundary condition (m)
                    upperBoundHead,        & ! intent(in): upper boundary condition (m)
                    ! input: flux at the upper boundary
                    scalarRainPlusMelt,    & ! intent(in): rain plus melt (m s-1)
                    ! input: transmittance
                    mLayerSatHydCond,      & ! intent(in): saturated hydraulic conductivity in each layer (m s-1)
                    iLayerSatHydCond,      & ! intent(in): saturated hydraulic conductivity at layer interfaces (m s-1)
                    mLayerHydCond,         & ! intent(in): hydraulic conductivity in each layer (m s-1)
                    iLayerHydCond,         & ! intent(in): hydraulic conductivity at layer interfaces (m s-1)
                    iceImpedeFac,          & ! intent(in): ice impedence factor in each layer (-)
                    ! input: soil parameters
                    vGn_alpha,             & ! intent(in): van Genutchen "alpha" parameter (m-1)
                    vGn_n,                 & ! intent(in): van Genutchen "n" parameter (-)
                    VGn_m,                 & ! intent(in): van Genutchen "m" parameter (-)
                    theta_sat,             & ! intent(in): soil porosity (-)
                    theta_res,             & ! intent(in): soil residual volumetric water content (-)
                    kAnisotropic,          & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                    zScale_TOPMODEL,       & ! intent(in): scale factor for TOPMODEL-ish baseflow parameterization (m)
                    bpar_VIC,              & ! intent(in): b-parameter in the VIC surface runoff parameterization (-)
                    f_impede,              & ! intent(in): ice impedence factor (-)
                    ! output
                    dq_dMatricAbove,       & ! intent(out): derivatives in the flux w.r.t. matric head in the layer above (s-1)
                    dq_dMatricBelow,       & ! intent(out): derivatives in the flux w.r.t. matric head in the layer below (s-1)
                    err,cmessage)            ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
 endselect


 ! *****
 ! compute the tri-diagonal matrix, and solve the tri-diagonal system of equations
 wtim = (1._dp - wimplicit)*dt  ! weighted time
 select case(ixRichards)
  ! ***** moisture-based form of Richards' equation
  case(moisture)
   !diag = (wtim/mLayerDepth)*(-dq_dVolLiqBelow(0:nLevels-1) + dq_dVolLiqAbove(1:nLevels) ) + 1._dp
   diag = (wtim/mLayerDepth)*(-dq_dVolLiqBelow(0:nLevels-1) + dq_dVolLiqAbove(1:nLevels) + mLayerEjectWaterDeriv(1:nLevels)) + 1._dp
   d_m1 = (wtim/mLayerDepth(1:nLevels-1))*(-dq_dVolLiqAbove(1:nLevels-1) )
   d_p1 = (wtim/mLayerDepth(1:nLevels-1))*( dq_dVolLiqBelow(1:nLevels-1) )
   call tridag(d_m1,diag,d_p1,-rvec,mLayerVolFracLiqDiff,err,cmessage)
   if(err/=0)then; write(*,'(50(e20.10,1x))') mLayerVolFracLiqIter; message=trim(message)//trim(cmessage); return; endif
  ! ***** mixed form of Richards' equation
  case(mixdform)
   !diag = (wtim/mLayerDepth)*(-dq_dMatricBelow(0:nLevels-1) + dq_dMatricAbove(1:nLevels)) + mLayerdTheta_dPsi
   diag = (wtim/mLayerDepth)*(-dq_dMatricBelow(0:nLevels-1) + dq_dMatricAbove(1:nLevels) + mLayerEjectWaterDeriv(1:nLevels)) + mLayerdTheta_dPsi
   d_m1 = (wtim/mLayerDepth(1:nLevels-1))*(-dq_dMatricAbove(1:nLevels-1) )
   d_p1 = (wtim/mLayerDepth(1:nLevels-1))*( dq_dMatricBelow(1:nLevels-1) )
   call tridag(d_m1,diag,d_p1,-rvec,mLayerMatricHeadDiff,err,cmessage)
   if(err/=0)then; write(*,'(50(e20.10,1x))') mLayerMatricHeadIter; message=trim(message)//trim(cmessage); return; endif
  ! ***** unknown form of Richards' equation
  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
 endselect
 !print*, 'mLayerMatricHeadDiff(1:5) = ', mLayerMatricHeadDiff(1:5)

 ! ***** compute Jacobian (used for testing)
 if(calcJacobian)then

  ! assemble desired variables and perturbations
  xtry=mLayerVolFracLiqIter
  xsav=xtry
  h=EPS*abs(xsav)
  where (h == 0.0) h=EPS
  xph=xsav+h
  h=xph-xsav

  ! loop through desired layers
  do ijac=minLayer,maxLayer

   ! get a new value of the state vector
   xtry(ijac)=xph(ijac)

   ! *****
   ! compute the fluxes, given xtry
   call computeFlux(&
                    ! input: trial state variables  (NOTE: use vectors from the start of the step)
                    mLayerMatricHead,       & ! intent(in): matric head (m)
                    mLayerVolFracIce,       & ! intent(in): volumetric fraction of ice (-)
                    xtry,                   & ! intent(in): volumetric fraction of liquid water (-)
                    ! output: fluxes (NOTE: populate temporary vectors)
                    mLayerTempEjectWater,   & ! intent(out): water ejected because pore volume is close to capacity (m s-1)
                    iLayerTempLiqFluxSoil,  & ! intent(out): liquid fluxes at layer interfaces at the start of the time step(m s-1)
                    ! ----
                    ! NOTE: we have no need for the optional arguments...
                    ! -------------------------------------------------
                    ! output: error control
                    err=err,message=cmessage) ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! *****
   ! compute the residual vector
   call liqResidual(&
                    ! control variables
                    dt,                        & ! length of the time step (s)
                    wimplicit,                 & ! weight assigned to the start-of-step
                    ! coordinate variables
                    mLayerDepth,               & ! depth of each layer (m)
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
                    xtry,                      & ! intent(in): trial volumetric liquid water content (-)
                    mLayerVolFracIceIter,      & ! intent(in): trial volumetric ice content (-)
                    ! intent(out): residual vector (-)
                    ftest,                     & ! intent(out): residual vector (-)
                    ! intent(out): error control
                    err,cmessage)                ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! populate the jacobian matrix and update x
   jmat(:,ijac)=(ftest(:)-rVec(:))/h(ijac)
   xtry(ijac)=xsav(ijac)

  end do  ! looping through soil layers

  ! test
  print*, 'Jacobian = '
  do ijac=minLayer,maxLayer
   write(*,'(i4,1x,10(e20.10,1x))') ijac, jmat(minLayer:50,ijac)
  end do
  write(*,'(a,1x,10(e20.10,1x))') 'diag = ', diag(minLayer:maxLayer)
  write(*,'(a,1x,10(e20.10,1x))') 'd_m1 = ', d_m1(minLayer:maxLayer)
  write(*,'(a,1x,10(e20.10,1x))') 'd_p1 = ', d_p1(minLayer:maxLayer)
  pause ' computing Jacobian'

 endif  ! (if computing the Jacobian)

 ! ***** update the fluxes and state variables
 do iLayer=0,nLevels
  select case(ixRichards)
   ! ***** moisture-based form of Richards' equation
   case(moisture)
    ! (check that soil moisture does not exceed constraints)
    if(iLayer > 0)then
     ! (negative error code will force time step reduction)
     !if(mLayerVolFracLiqDiff(iLayer) > theta_sat - mLayerVolFracLiqIter(iLayer))then; err=-20; message=trim(message)//'soil moisture exceeds saturation'; return; endif
     !if(mLayerVolFracLiqDiff(iLayer) < theta_res - mLayerVolFracLiqIter(iLayer))then; err=-20; message=trim(message)//'soil moisture below residual content'; return; endif
     ! (use simplified bi-section method)
     availPorosity = theta_sat - mLayerVolFracLiqIter(iLayer) - mLayerVolFracIceIter(iLayer)*(iden_ice/iden_water)
     availLiqWater = theta_res - mLayerVolFracLiqIter(iLayer)
     if(mLayerVolFracLiqDiff(iLayer) > availPorosity) mLayerVolFracLiqDiff(iLayer) = availPorosity/2._dp
     if(mLayerVolFracLiqDiff(iLayer) < availLiqWater) mLayerVolFracLiqDiff(iLayer) = availLiqWater/2._dp
    endif
    ! (update the fluxes)
    if(iLayer==0)then;           iLayerLiqFluxSoil(iLayer) = iLayerLiqFluxSoil(iLayer) + dq_dVolLiqBelow(iLayer)*mLayerVolFracLiqDiff(iLayer+1)
    elseif(iLayer==nLevels)then; iLayerLiqFluxSoil(iLayer) = iLayerLiqFluxSoil(iLayer) + dq_dVolLiqAbove(iLayer)*mLayerVolFracLiqDiff(iLayer)
    else;                        iLayerLiqFluxSoil(iLayer) = iLayerLiqFluxSoil(iLayer) + dq_dVolLiqAbove(iLayer)*mLayerVolFracLiqDiff(iLayer) &
                                                                                       + dq_dVolLiqBelow(iLayer)*mLayerVolFracLiqDiff(iLayer+1)
    endif
    if(iLayer > 0) mLayerEjectWater(iLayer) = mLayerEjectWater(iLayer) + mLayerEjectWaterDeriv(iLayer)*mLayerVolFracLiqDiff(iLayer)
    ! (update the state variables)
    if(iLayer > 0)then
     ! (update)
     mLayerVolFracLiqNew(iLayer) = mLayerVolFracLiqIter(iLayer) + mLayerVolFracLiqDiff(iLayer)
     mLayerMatricHeadNew(iLayer) = matricHead(mLayerVolFracLiqNew(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
     ! (check)
     if(mLayerVolFracLiqNew(iLayer) + mLayerVolFracIceIter(iLayer)*(iden_ice/iden_water) > theta_sat)then
      print*, 'availPorosity, theta_sat, mLayerVolFracLiqIter(iLayer), mLayerVolFracIceIter(iLayer)*(iden_ice/iden_water), mLayerVolFracLiqDiff(iLayer), mLayerVolFracLiqNew(iLayer) = '
      print*,  availPorosity, theta_sat, mLayerVolFracLiqIter(iLayer), mLayerVolFracIceIter(iLayer)*(iden_ice/iden_water), mLayerVolFracLiqDiff(iLayer), mLayerVolFracLiqNew(iLayer)
      message=trim(message)//'volumetric (liquid + ice) content exceeds soil porosity'
      err=20; return
     endif
     !if(iLayer < 5) write(*,'(2(a,1x,3(e20.10,1x)))') 'VolFracLiq = ', mLayerVolFracLiqNew(iLayer), mLayerVolFracLiqIter(iLayer), mLayerVolFracLiqDiff(iLayer), &
     !                                                 'matricHead = ', mLayerMatricHeadNew(iLayer), mLayerMatricHeadIter(iLayer), mLayerMatricHeadNew(iLayer) - mLayerMatricHeadIter(iLayer)
    endif
   ! ***** mixed form of Richards' equation
   case(mixdform)
    ! upper boundary
    if(iLayer==0)then
     iLayerLiqFluxSoil(iLayer) = iLayerLiqFluxSoil(iLayer) + dq_dMatricBelow(iLayer)*mLayerMatricHeadDiff(iLayer+1)
    ! lower boundary
    elseif(iLayer==nLevels)then
     iLayerLiqFluxSoil(iLayer) = iLayerLiqFluxSoil(iLayer) + dq_dMatricAbove(iLayer)*mLayerMatricHeadDiff(iLayer)
    ! internal interfaces
    else
     iLayerLiqFluxSoil(iLayer) = iLayerLiqFluxSoil(iLayer) + dq_dMatricAbove(iLayer)*mLayerMatricHeadDiff(iLayer) &
                                                           + dq_dMatricBelow(iLayer)*mLayerMatricHeadDiff(iLayer+1)
    endif
    if(iLayer > 0) mLayerEjectWater(iLayer) = mLayerEjectWater(iLayer) + mLayerEjectWaterDeriv(iLayer)*mLayerMatricHeadDiff(iLayer)
    ! (update the state variables)
    if(iLayer > 0)then
     mLayerMatricHeadNew(iLayer) = mLayerMatricHeadIter(iLayer) + mLayerMatricHeadDiff(iLayer)
     mLayerVolFracLiqNew(iLayer) = volFracLiq(mLayerMatricHeadNew(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)    
     ! check for ponding on the soil surface
     if(mLayerMatricHeadNew(1) > 0._dp)then
      print*, 'mLayerMatricHeadNew(1:5) = ', mLayerMatricHeadNew(1:5) 
      print*, 'iLayerLiqFluxSoil(0:1) = ', iLayerLiqFluxSoil(0:1)
      err=20; message=trim(message)//'ponding on the soil surface'; return
     endif
     !if(iLayer < 5) write(*,'(2(a,1x,3(e20.10,1x)))') 'VolFracLiq = ', mLayerVolFracLiqNew(iLayer), mLayerVolFracLiqIter(iLayer), mLayerVolFracLiqNew(iLayer) - mLayerVolFracLiqIter(iLayer), &
     !                                                 'matricHead = ', mLayerMatricHeadNew(iLayer), mLayerMatricHeadIter(iLayer), mLayerMatricHeadDiff(iLayer)
    endif
   ! ***** unknown form of Richards' equation
   case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
  endselect
 end do  ! (loop through layers)


 contains

  ! ************************************************************************************************
  ! internal subroutine: compute the vector of fluxes at layer interfaces
  ! ************************************************************************************************
  subroutine computeFlux(&
                         ! input: trial state variables (NOTE: use "trial" instead of "iter" to avoid using the "iter" vectors defined in main subroutine)
                         mLayerMatricHeadTrial, & ! intent(in): matric head (m)
                         mLayerVolFracIceTrial, & ! intent(in): volumetric fraction of ice (-)
                         mLayerVolFracLiqTrial, & ! intent(in): volumetric fraction of liquid water (-)
                         ! output: fluxes
                         mLayerEjectWater,      & ! intent(out): water ejected because pore volume is close to capacity (m s-1)
                         iLayerLiqFluxSoil,     & ! intent(out): liquid fluxes at layer interfaces at the start of the time step(m s-1)
                         ! output: derivative in the soil water characteristic
                         mLayerdPsi_dTheta,     &  ! derivative in the soil water characteristic
                         mLayerdTheta_dPsi,     &  ! derivative in the soil water characteristic
                         ! output: transmittance
                         mLayerHydCond,         & ! intent(out): hydraulic conductivity at layer mid-points (m s-1)
                         mLayerDiffuse,         & ! intent(out): diffusivity at layer mid-points (m2 s-1)
                         iLayerHydCond,         & ! intent(out): hydraulic conductivity at layer interface (m s-1)
                         iLayerDiffuse,         & ! intent(out): diffusivity at layer interface (m2 s-1)
                         iceImpedeFac,          & ! intent(out): ice impedence factor in each layer (-)
                         ! output: derivative in ejected water flux
                         mLayerEjectWaterDeriv, & ! intent(out): derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
                         ! output: surface runoff
                         scalarSurfaceRunoff,   & ! intent(out): surface runoff (m s-1)
                         ! output: error control
                         err,message)             ! intent(out): error control
  implicit none
  ! trial model state variables
  real(dp),intent(in)              :: mLayerMatricHeadTrial(:) ! matric head in each layer at the current iteration (m)
  real(dp),intent(in)              :: mLayerVolFracIceTrial(:) ! volumetric fraction of ice at the current iteration (-)
  real(dp),intent(in)              :: mLayerVolFracLiqTrial(:) ! volumetric fraction of liquid water at the current iteration (-)
  ! output: fluxes
  real(dp),intent(out)             :: mLayerEjectWater(:)      ! water ejected from each soil layer (m s-1)
  real(dp),intent(out)             :: iLayerLiqFluxSoil(0:)    ! liquid flux at soil layer interfaces at the end of the time step (m s-1)
  ! output: derivative in the soil water characteristic
  real(dp),intent(out),optional    :: mLayerdPsi_dTheta(:)     ! derivative in the soil water characteristic
  real(dp),intent(out),optional    :: mLayerdTheta_dPsi(:)     ! derivative in the soil water characteristic
  ! output: transmittance
  real(dp),intent(out),optional    :: mLayerHydCond(:)         ! hydraulic conductivity at layer mid-point (m s-1)
  real(dp),intent(out),optional    :: mLayerDiffuse(:)         ! diffusivity at layer mid-point (m2 s-1)
  real(dp),intent(out),optional    :: iLayerHydCond(0:)        ! hydraulic conductivity at layer interface (m s-1)
  real(dp),intent(out),optional    :: iLayerDiffuse(0:)        ! diffusivity at layer interface (m2 s-1)
  real(dp),intent(out),optional    :: iceImpedeFac(:)          ! ice impedence factor at layer mid-points (-)
  ! output: flux of water ejected when volumetric liquid water content is close to exceeding porosity
  real(dp),intent(out),optional    :: mLayerEjectWaterDeriv(:) ! derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
  ! output: surface runoff
  real(dp),intent(out),optional    :: scalarSurfaceRunoff      ! surface runoff (m s-1)
  ! output: error control
  integer(i4b),intent(out)         :: err                      ! error code
  character(*),intent(out)         :: message                  ! error message
  ! local variables (general)
  character(LEN=256)               :: cmessage                 ! error message of downwind routine
  real(dp),dimension(nLevels)      :: local__mLayerdPsi_dTheta ! derivative in the soil water characteristic
  real(dp),dimension(nLevels)      :: local__mLayerdTheta_dPsi ! derivative in the soil water characteristic
  real(dp),dimension(nLevels)      :: local__mLayerHydCond     ! hydraulic conductivity at layer mid-point (m s-1)
  real(dp),dimension(nLevels)      :: local__mLayerDiffuse     ! diffusivity at layer mid-point (m2 s-1)
  real(dp),dimension(0:nLevels)    :: local__iLayerHydCond     ! hydraulic conductivity at layer interface (m s-1)
  real(dp),dimension(0:nLevels)    :: local__iLayerDiffuse     ! diffusivity at layer interface (m2 s-1)
  real(dp),dimension(nLevels)      :: local__iceImpedeFac      ! ice impedence factor at layer mid-points (-)
  real(dp),dimension(nLevels)      :: local__mLayerEjectWaterDeriv ! derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
  real(dp)                         :: local__scalarSurfaceRunoff   ! surface runoff (m s-1)
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message='computeFlux/'

  ! ***** 
  ! compute the derivative in the soil water characteristic
  do iLayer=1,nLevels
   select case(ixRichards)
    case(moisture)
     local__mLayerdPsi_dTheta(iLayer) = dPsi_dTheta(mLayervolFracLiqTrial(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
     local__mLayerdTheta_dPsi(iLayer) = valueMissing  ! (deliberately cause problems if this is ever used)
    case(mixdform)
     local__mLayerdTheta_dPsi(iLayer) = dTheta_dPsi(mLayerMatricHeadTrial(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
     local__mLayerdPsi_dTheta(iLayer) = valueMissing  ! (deliberately cause problems if this is ever used)
    case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
   end select
  end do

  ! *****
  ! compute the hydraulic conductivity for all layers
  call hydCond_all(&
                   ! input: model control
                   hc_Profile,                   & ! intent(in): index defining the option for the hydraulic conductivity profile
                   ixRichards,                   & ! intent(in): index defining the  option for Richards' equation (moisture or mixdform)
                   ixBcUpperSoilHydrology,       & ! intent(in): index defining the type of upper boundary conditions 
                   ixBcLowerSoilHydrology,       & ! intent(in): index defining the type of upper boundary conditions 
                   ! input: state and diagnostic variables
                   mLayerMatricHeadTrial,        & ! intent(in): matric head (m)
                   mLayerVolFracIceTrial,        & ! intent(in): volumetric fraction of ice (-)
                   mLayerVolFracLiqTrial,        & ! intent(in): volumetric fraction of liquid water (-)
                   local__mLayerdPsi_dTheta,     & ! intent(in): derivative in the soil water characteristic w.r.t. theta (m)
                   ! input: water table depth
                   scalarWaterTableDepth,        & ! intent(in): depth to the water table (m)
                   ! input: boundary conditions for matric head
                   lowerBoundHead,               & ! intent(in): lower boundary condition (m)
                   upperBoundHead,               & ! intent(in): upper boundary condition (m)
                   ! input: boundary conditions for volumetric liquid water content
                   lowerBoundTheta,              & ! intent(in): lower boundary condition (-)
                   upperBoundTheta,              & ! intent(in): upper boundary condition (-)
                   ! input: soil parameters
                   vGn_alpha,                    & ! intent(in): van Genutchen "alpha" parameter (m-1)
                   vGn_n,                        & ! intent(in): van Genutchen "n" parameter (-)
                   VGn_m,                        & ! intent(in): van Genutchen "m" parameter (-)
                   theta_sat,                    & ! intent(in): soil porosity (-)
                   theta_res,                    & ! intent(in): soil residual volumetric water content (-)
                   f_impede,                     & ! intent(in): ice impedence factor (-)
                   zScale_TOPMODEL,              & ! intent(in): TOPMODEL scaling factor (m)
                   ! saturated hydraulic conductivity
                   mLayerSatHydCond,             & ! intent(in):  saturated hydraulic conductivity at the mid-point of each layer (m s-1)
                   iLayerSatHydCond,             & ! intent(in):  saturated hydraulic conductivity at the interface of each layer (m s-1)
                   ! output: transmittance
                   local__mLayerHydCond,         & ! intent(out): hydraulic conductivity at layer mid-points (m s-1)
                   local__mLayerDiffuse,         & ! intent(out): diffusivity at layer mid-points (m2 s-1)
                   local__iLayerHydCond,         & ! intent(out): hydraulic conductivity at layer interface (m s-1)
                   local__iLayerDiffuse,         & ! intent(out): diffusivity at layer interface (m2 s-1)
                   local__iceImpedeFac,          & ! intent(out): ice impedence factor in each layer (-)
                   ! output: error control
                   err,cmessage)             ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


  ! *****
  ! compute the flux of water ejected because the volumetric liquid water content is close to porosity
  call ejectWater(&
                  ! input: model control
                  ixDerivMethod,                & ! intent(in): index defining the method used to compute derivatives (analytical or numerical)
                  ixRichards,                   & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                  ! input: state variables
                  mLayerMatricHeadTrial,        & ! intent(in): matric head in each layer (m)
                  mLayerVolFracIceTrial,        & ! intent(in): volumetric ice content in each layer (-)
                  mLayerVolFracLiqTrial,        & ! intent(in): volumetric liquid water content (-)
                  local__mLayerdTheta_dPsi,     & ! intent(in): derivative in the soil water characteristic w.r.t. psi (m-1)
                  ! input: soil parameters
                  vGn_alpha,                    & ! intent(in): van Genutchen "alpha" parameter (m-1)
                  vGn_n,                        & ! intent(in): van Genutchen "n" parameter (-)
                  VGn_m,                        & ! intent(in): van Genutchen "m" parameter (-)
                  theta_sat,                    & ! intent(in): soil porosity (-)
                  theta_res,                    & ! intent(in): soil residual volumetric water content (-)
                  ! saturated hydraulic conductivity in each layer
                  mLayerSatHydCond,             & ! intent(in):  saturated hydraulic conductivity in each layer (m s-1)
                  ! output: ejected water flux and derivative
                  mLayerEjectWater,             & ! intent(out): water ejected because pore volume is filled (m s-1)
                  local__mLayerEjectWaterDeriv, & ! intent(out): derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
                  ! output: error control
                  err,cmessage)           ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


  ! *****
  ! compute the fluxes at layer interfaces
  call iLayer_liq(&
                  ! input: model control
                  ixRichards,                   & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                  ixBcUpperSoilHydrology,       & ! index defining the type of upper boundary conditions 
                  ixBcLowerSoilHydrology,       & ! index defining the type of upper boundary conditions 
                  ! input: model coordinate variables
                  mLayerDepth,                  & ! intent(in): depth of the layer (m)
                  mLayerHeight,                 & ! intent(in): height of the layer mid-point (m)
                  ! input: state variables
                  mLayerMatricHeadTrial,        & ! intent(in): matric head in each layer (m)
                  mLayerVolFracIceTrial,        & ! intent(in): volumetric ice content in each layer (-)
                  mLayerVolFracLiqTrial,        & ! intent(in): volumetric liquid water content (-)
                  ! input: boundary conditions for matric head
                  lowerBoundHead,               & ! intent(in): lower boundary condition (m)
                  upperBoundHead,               & ! intent(in): upper boundary condition (m)
                  ! input: boundary conditions for volumetric liqquid water content
                  lowerBoundTheta,              & ! intent(in): lower boundary condition (-)
                  upperBoundTheta,              & ! intent(in): upper boundary condition (-)
                  ! input: flux at the upper boundary
                  scalarRainPlusMelt,           & ! intent(in): rain plus melt (m s-1)
                  ! input: height of the lower boundary and the water table
                  lowerBoundaryHeight,          & ! intent(in): height of the lower boundary (m)
                  scalarWaterTableDepth,        & ! intent(in): water table depth (m)
                  ! input: transmittance
                  iLayerSatHydCond,             & ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
                  local__iLayerHydCond,         & ! intent(in): hydraulic conductivity at layer interfaces (m s-1)
                  local__iLayerDiffuse,         & ! intent(in): hydraulic diffusivity at layer interfaces (m2 s-1)
                  local__iceImpedeFac,          & ! intent(in): ice impedence factor in each layer (-)
                  ! input: soil parameters
                  vGn_alpha,                    & ! intent(in): van Genutchen "alpha" parameter (m-1)
                  vGn_n,                        & ! intent(in): van Genutchen "n" parameter (-)
                  VGn_m,                        & ! intent(in): van Genutchen "m" parameter (-)
                  theta_sat,                    & ! intent(in): soil porosity (-)
                  theta_res,                    & ! intent(in): soil residual volumetric water content (-)
                  kAnisotropic,                 & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                  zScale_TOPMODEL,              & ! intent(in): scale factor for TOPMODEL-ish baseflow parameterization (m)
                  bpar_VIC,                     & ! intent(in): b-parameter in the VIC surface runoff parameterization (-)
                  ! output: fluxes at layer interfaces and surface runoff
                  iLayerLiqFluxSoil,            & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                  local__scalarSurfaceRunoff,   & ! intent(out): surface runoff (m s-1)
                  ! output: error control
                  err,cmessage)                   ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! *****
  ! retrive optional arguments
  if( present(mLayerdPsi_dTheta)    ) mLayerdPsi_dTheta     = local__mLayerdPsi_dTheta     ! derivative in the soil water characteristic
  if( present(mLayerdTheta_dPsi)    ) mLayerdTheta_dPsi     = local__mLayerdTheta_dPsi     ! derivative in the soil water characteristic
  if( present(mLayerHydCond)        ) mLayerHydCond         = local__mLayerHydCond         ! hydraulic conductivity at layer mid-point (m s-1)
  if( present(mLayerDiffuse)        ) mLayerDiffuse         = local__mLayerDiffuse         ! diffusivity at layer mid-point (m2 s-1)
  if( present(iLayerHydCond)        ) iLayerHydCond         = local__iLayerHydCond         ! hydraulic conductivity at layer interface (m s-1)
  if( present(iLayerDiffuse)        ) iLayerDiffuse         = local__iLayerDiffuse         ! diffusivity at layer interface (m2 s-1)
  if( present(iceImpedeFac)         ) iceImpedeFac          = local__iceImpedeFac          ! ice impedence factor at layer mid-points (-)
  if( present(mLayerEjectWaterDeriv)) mLayerEjectWaterDeriv = local__mLayerEjectWaterDeriv ! derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
  if( present(scalarSurfaceRunoff)  ) scalarSurfaceRunoff   = local__scalarSurfaceRunoff   ! surface runoff (m s-1)


  end subroutine computeFlux


 end subroutine soilHydrol_muster


 ! ---------------------------------------------------------------------------------------------------------------------------
 ! public subroutine: compute the water required to fill pore space at each layer interface (m)
 ! ---------------------------------------------------------------------------------------------------------------------------
 subroutine waterReq2FillLayer(theta_sat,            & ! intent(in):  total porosity (-)
                               mLayerDepth,          & ! intent(in):  depth of each soil layer (m)
                               mLayerVolFracLiqIter, & ! intent(in):  volumetric liquid water content in each soil layer (-)
                               mLayerVolFracIceIter, & ! intent(in):  volumetric ice content in each soil layer (-)
                               waterReq2FillPore,    & ! intent(out): water required to fill pore space at each layer interface (m)
                               err,message)            ! intent(out): error control
 implicit none
 ! input
 real(dp),intent(in)        :: theta_sat               ! total porosity (-)
 real(dp),intent(in)        :: mLayerDepth(:)          ! depth of each soil layer (m)
 real(dp),intent(in)        :: mLayerVolFracLiqIter(:) ! volumetric liquid water content in each soil layer (-)
 real(dp),intent(in)        :: mLayerVolFracIceIter(:) ! volumetric ice content in each soil layer (-)
 ! output
 real(dp),intent(out)       :: waterReq2FillPore(0:)   ! water required to fill pore space at each layer interface (m)
 integer(i4b),intent(out)   :: err                     ! error code
 character(*),intent(out)   :: message                 ! error message
 ! internal
 integer(i4b)               :: iSoil                   ! index of soil layers
 integer(i4b)               :: nSoil                   ! number of soil layers
 real(dp)                   :: drainablePorosity       ! drainable porosity (-)
 ! initialize error control
 err=0; message='waterReq2FillLayer/'

 ! get the number of soil layers
 nSoil = size(mLayerVolFracLiqIter)

 ! initialize the water required to fill pore space (m)
 waterReq2FillPore(nSoil) = 0._dp ! (no water required to fill pore space at the bottom interface)

 ! compute the water required to fill pore space at each layer interface (m)
 do iSoil=nSoil,1,-1  ! (loop through soil layers, starting at the bottom)
  drainablePorosity = theta_sat - mLayerVolFracLiqIter(iSoil) - mLayerVolFracIceIter(iSoil)
  waterReq2FillPore(iSoil-1) = waterReq2FillPore(iSoil) + drainablePorosity*mLayerDepth(iSoil)
 end do

 end subroutine waterReq2FillLayer




 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ***** PRIVATE SUBROUTINES **********************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************

 ! compute the residual vector
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
  if(iLayer==1)                     print*, 'iLayer, residualVec(iLayer), bottom flux, mFlux, mEvap, mEjct, mPhse'
  if(iLayer > 25 .and. iLayer < 30) write(*,'(i4,1x,10(e20.10,1x))') iLayer, residualVec(iLayer), iLayerTrialLiqFluxSoil(iLayer)*dt_dz, mFlux, mEvap, mEjct, mPhse
 end do  ! (looping through soil layers)

 end subroutine liqResidual


 ! ************************************************************************************************
 ! private subroutine: compute hydraulic conductivity
 ! ************************************************************************************************
 subroutine hydCond_all(&
                        ! input: model control
                        hc_Profile,            & ! intent(in): index defining the option for the hydraulic conductivity profile
                        ixRichards,            & ! intent(in): index defining the option for Richards' equation (moisture or mixdform)
                        bc_upper,bc_lower,     & ! intent(in): index defining the type of boundary conditions (neumann or diriclet)
                        ! input: state and diagnostic variables
                        mLayerMatricHeadTrial, & ! intent(in): matric head (m)
                        mLayerVolFracIceTrial, & ! intent(in): volumetric fraction of ice (-)
                        mLayerVolFracLiqTrial, & ! intent(in): volumetric fraction of liquid water (-)
                        mLayerdPsi_dTheta,     & ! intent(in): derivative in the soil water characteristic w.r.t. theta (m)
                        ! input: water table depth
                        scalarWaterTableDepth, & ! intent(in): depth to the water table (m)
                        ! boundary conditions for matric head
                        lowerBoundHead,        & ! intent(in): lower boundary condition (m)
                        upperBoundHead,        & ! intent(in): upper boundary condition (m)
                        ! boundary conditions for volumetric liquid water content
                        lowerBoundTheta,       & ! intent(in): lower boundary condition (-)
                        upperBoundTheta,       & ! intent(in): upper boundary condition (-)
                        ! input: soil parameters
                        vGn_alpha,             & ! intent(in): van Genutchen "alpha" parameter (m-1)
                        vGn_n,                 & ! intent(in): van Genutchen "n" parameter (-)
                        VGn_m,                 & ! intent(in): van Genutchen "m" parameter (-)
                        theta_sat,             & ! intent(in): soil porosity (-)
                        theta_res,             & ! intent(in): soil residual volumetric water content (-)
                        f_impede,              & ! intent(in): ice impedence factor (-)
                        zScale_TOPMODEL,       & ! intent(in): TOPMODEL scaling factor (m)
                        ! input: saturated hydraulic conductivity
                        mLayerSatHydCond,      & ! intent(in): saturated hydraulic conductivity at the mid-point of each layer (m s-1)
                        iLayerSatHydCond,      & ! intent(in): saturated hydraulic conductivity at the interface of each layer (m s-1)
                        ! output: transmittance
                        mLayerHydCond,         & ! intent(out): hydraulic conductivity at layer mid-points (m s-1)
                        mLayerDiffuse,         & ! intent(out): diffusivity at layer mid-points (m2 s-1)
                        iLayerHydCond,         & ! intent(out): hydraulic conductivity at layer interface (m s-1)
                        iLayerDiffuse,         & ! intent(out): diffusivity at layer interface (m2 s-1)
                        iceImpedeFac,          & ! intent(out): ice impedence factor in each layer (-)
                        ! output: error control
                        err,message)             ! intent(out): error control
 USE soil_utils_module,only:iceImpede       ! compute the ice impedence factor
 USE soil_utils_module,only:hydCond_psi     ! compute hydraulic conductivity as a function of matric head
 USE soil_utils_module,only:hydCond_liq     ! compute hydraulic conductivity as a function of volumetric liquid water content
 USE soil_utils_module,only:dPsi_dTheta     ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
 ! compute hydraulic conductivity for all layers
 implicit none
 ! input: model control
 integer(i4b),intent(in)       :: hc_profile                ! choice of type of hydraulic conductivity profile
 integer(i4b),intent(in)       :: ixRichards                ! index defining the option for Richards' equation (moisture or mixdform)
 integer(i4b),intent(in)       :: bc_upper,bc_lower         ! index defining the type of boundary conditions (neumann or diriclet)
 ! input: state and diagnostic variables
 real(dp),intent(in)           :: mLayerMatricHeadTrial(:)  ! matric head in each layer (m)
 real(dp),intent(in)           :: mLayerVolFracIceTrial(:)  ! volumetric fraction of ice in each layer (-)
 real(dp),intent(in)           :: mLayerVolFracLiqTrial(:)  ! volumetric fraction of liquid water in each layer (-)
 real(dp),intent(in)           :: mLayerdPsi_dTheta(:)      ! derivative in the soil water characteristic w.r.t. theta (m)
 ! water table depth
 real(dp),intent(in)           :: scalarWaterTableDepth     ! depth to the water table (m)
 ! diriclet boundary conditions for matric head
 real(dp),intent(in)           :: lowerBoundHead            ! lower boundary condition for matric head (m)
 real(dp),intent(in)           :: upperBoundHead            ! upper boundary condition for matric head (m)
 ! diriclet boundary conditions for volumetric liquid water content
 real(dp),intent(in)           :: lowerBoundTheta           ! lower boundary condition for volumetric liquid water content (-)
 real(dp),intent(in)           :: upperBoundTheta           ! upper boundary condition for volumetric liquid water content (-)
 ! input: soil parameters
 real(dp),intent(in)           :: vGn_alpha                 ! van Genutchen "alpha" parameter (m-1)
 real(dp),intent(in)           :: vGn_n                     ! van Genutchen "n" parameter (-)
 real(dp),intent(in)           :: vGn_m                     ! van Genutchen "m" parameter (-)
 real(dp),intent(in)           :: theta_sat                 ! soil porosity (-)
 real(dp),intent(in)           :: theta_res                 ! soil residual volumetric water content (-)
 real(dp),intent(in)           :: f_impede                  ! ice impedence factor (-)
 real(dp),intent(in)           :: zScale_TOPMODEL           ! TOPMODEL scaling factor (m)
 ! input: saturated hydraulic conductivity
 real(dp),intent(in)           :: mLayerSatHydCond(:)       ! saturated hydraulic conductivity at the mid-point of each layer (m s-1)
 real(dp),intent(in)           :: iLayerSatHydCond(0:)      ! saturated hydraulic conductivity at the interface of each layer (m s-1)
 ! output: transmittance
 real(dp),intent(out)          :: mLayerHydCond(:)          ! hydraulic conductivity at layer mid-points (m s-1)
 real(dp),intent(out)          :: mLayerDiffuse(:)          ! diffusivity at layer mid-points (m2 s-1)
 real(dp),intent(out)          :: iLayerHydCond(0:)         ! hydraulic conductivity at layer interface (m s-1)
 real(dp),intent(out)          :: iLayerDiffuse(0:)         ! diffusivity at layer interface (m2 s-1)
 real(dp),intent(out)          :: iceImpedeFac(:)           ! ice impedence factor in each layer (-)
 ! output: error control
 integer(i4b),intent(out)      :: err                       ! error code
 character(*),intent(out)      :: message                   ! error message
 ! local variables
 real(dp)                      :: test                      ! scalar for testing
 integer(i4b)                  :: iLayer                    ! layer index
 real(dp)                      :: hk_depth                  ! hydraulic conductivity at the water table depth (m s-1)
 ! initialize error control
 err=0; message="hydCond_all/"

 ! initialize diffusivity for the mixed form of Richards' equation (not used)
 mLayerDiffuse(:) = valueMissing
 iLayerDiffuse(:) = valueMissing

 ! loop through layers
 do iLayer=1,nLevels
  ! compute the ice impedence factor (-)
  iceImpedeFac(iLayer) = iceImpede(mLayerVolFracIceTrial(iLayer),theta_res,theta_sat,f_impede)
  ! compute the hydraulic conductivity (m s-1) and diffusivity (m2 s-1) for a given layer
  select case(ixRichards)
   ! ***** moisture-based form of Richards' equation
   case(moisture)
    mLayerHydCond(iLayer) = hydCond_liq(mLayerVolFracLiqTrial(iLayer),mLayerSatHydCond(iLayer),theta_res,theta_sat,vGn_m) * iceImpedeFac(iLayer)
    mLayerDiffuse(iLayer) = mLayerdPsi_dTheta(iLayer) * mLayerHydCond(iLayer)
   ! ***** mixed form of Richards' equation -- just compute hydraulic condictivity
   case(mixdform); mLayerHydCond(iLayer) = hydCond_psi(mLayerMatricHeadTrial(iLayer),mLayerSatHydCond(iLayer),vGn_alpha,vGn_n,vGn_m) * iceImpedeFac(iLayer)
   case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
  endselect
  ! compute the hydraulic conductivity (m s-1) and diffusivity (m2 s-1) for layer interfaces
  if(iLayer>1)then  ! *** NOTE: use the geometric mean
   iLayerHydCond(iLayer-1) = (mLayerHydCond(iLayer-1) * mLayerHydCond(iLayer))**0.5_dp
   if(ixRichards == moisture)&
   iLayerDiffuse(iLayer-1) = (mLayerDiffuse(iLayer-1) * mLayerDiffuse(iLayer))**0.5_dp
  endif
 end do

 ! compute hydraulic conductivity at the upper boundaries
 select case(bc_upper)
  ! ***** prescribed head at the upper boundary
  case(prescribedHead)
   select case(ixRichards)  ! (moisture-based form of Richards' equation)
    case(moisture)
     iLayerHydCond(0) = hydCond_liq(upperBoundTheta,iLayerSatHydCond(0),theta_res,theta_sat,vGn_m) * iceImpedeFac(1)
     iLayerDiffuse(0) = dPsi_dTheta(upperBoundTheta,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) * iLayerHydCond(0)
    case(mixdform)
     iLayerHydCond(0) = hydCond_psi(upperBoundHead,iLayerSatHydCond(0),vGn_alpha,vGn_n,vGn_m) * iceImpedeFac(1)
    case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
   end select  ! (form of Richards' eqn)
  ! ***** all other cases, set hydraulic conductivity at boundaries to the layer conductivity
  case(liquidFlux)
   iLayerHydCond(0)   = mLayerHydCond(1)
  case default; err=10; message=trim(message)//"unknown upper  boundary condition"; return
 end select

 ! compute hydraulic conductivity at the lower boundaries
 select case(bc_lower)
  ! ***** use diriclet conditions at the lower boundary
  case(prescribedHead)
   select case(ixRichards)  ! (moisture-based form of Richards' equation)
    case(moisture)
     iLayerHydCond(nLevels) = hydCond_liq(lowerBoundTheta,iLayerSatHydCond(nLevels),theta_res,theta_sat,vGn_m) * iceImpedeFac(nLevels)
     iLayerDiffuse(nLevels) = dPsi_dTheta(lowerBoundTheta,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) * iLayerHydCond(nLevels)
    case(mixdform)
     iLayerHydCond(nLevels) = hydCond_psi(lowerBoundHead,iLayerSatHydCond(nLevels),vGn_alpha,vGn_n,vGn_m) * iceImpedeFac(nLevels)
    case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
   end select  ! (form of Richards' eqn)
  ! ***** moving groundwater boundary
  case(groundwaterCouple)
   ! get the hydraulic conductivity at the water table depth
   select case(hc_profile)
    case(constant);    hk_depth = iLayerSatHydCond(0)
    case(exp_profile); hk_depth = iLayerSatHydCond(0) * exp(-scalarWaterTableDepth/zScale_TOPMODEL)
    case(powerLaw_profile,linear_profile)
     message=trim(message)//'hydraulic conductivity profile not implemented yet for "powerLaw" and "linear" options'; err=10; return
    case default
     message=trim(message)//'unknown hydraulic conductivity profile'; err=10; return
   end select
   iLayerHydCond(nLevels) = (mLayerHydCond(nLevels) * hk_depth)**0.5_dp

  ! ***** all other cases, set hydraulic conductivity at boundaries to the layer conductivity
  case(funcBottomHead,freeDrainage)
   iLayerHydCond(nLevels)   = mLayerHydCond(nLevels)
  case default; err=10; message=trim(message)//"unknown lower boundary condition"; return
 end select
 
 end subroutine hydCond_all


 ! ************************************************************************************************
 ! private subroutine: compute flux of water ejected because exceeding porosity
 ! ************************************************************************************************
 subroutine ejectWater(&
                       ! input: model control
                       dMethod,               & ! intent(in): index defining the method used to compute derivatives (analytical or numerical)
                       ixRichards,            & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                       ! input: state variable and diagnostic variables
                       mLayerMatricHeadTrial, & ! intent(in): matric head in each layer (m)
                       mLayerVolFracIceTrial, & ! intent(in): volumetric ice content in each layer (-)
                       mLayerVolFracLiqTrial, & ! intent(in): volumetric liquid water content (-)
                       mLayerdTheta_dPsi,     & ! intent(in): derivative in the soil water characteristic w.r.t. psi (m-1)
                       ! input: soil parameters
                       vGn_alpha,             & ! intent(in): van Genutchen "alpha" parameter (m-1)
                       vGn_n,                 & ! intent(in): van Genutchen "n" parameter (-)
                       VGn_m,                 & ! intent(in): van Genutchen "m" parameter (-)
                       theta_sat,             & ! intent(in): soil porosity (-)
                       theta_res,             & ! intent(in): soil residual volumetric water content (-)
                       ! input: saturated hydraulic conductivity
                       mLayerSatHydCond,      & ! intent(in): saturated hydraulic conductivity in each layer (m s-1)
                       ! output: ejected water flux and derivative
                       mLayerEjectWater,      & ! intent(out): water ejected because pore volume is filled (m s-1)
                       mLayerEjectWaterDeriv, & ! intent(out): derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
                       ! output: error control
                       err,message)             ! intent(out): error control
 USE soil_utils_module,only:volFracLiq          ! compute volumetric fraction of liquid water (-)
 ! avoid super-saturation, by computing flux of water due to displacement (m s-1)
 implicit none
 ! input: model control
 integer(i4b),intent(in)       :: dMethod                   ! index defining the method used to compute derivatives (analytical or numerical)
 integer(i4b),intent(in)       :: ixRichards                ! index defining the form of Richards' equation (moisture or mixdform)
 ! input: state variables
 real(dp),intent(in)           :: mLayerMatricHeadTrial(:)  ! matric head in each layer (m)
 real(dp),intent(in)           :: mLayerVolFracIceTrial(:)  ! volumetric ice content in each layer (-)
 real(dp),intent(in)           :: mLayerVolFracLiqTrial(:)  ! volumetric liquid water content in each layer (-)
 real(dp),intent(in)           :: mLayerdTheta_dPsi(:)      ! derivative in the soil water characteristic w.r.t. psi (m-1)
 ! input: soil parameters
 real(dp),intent(in)           :: vGn_alpha                 ! van Genutchen "alpha" parameter (m-1)
 real(dp),intent(in)           :: vGn_n                     ! van Genutchen "n" parameter (-)
 real(dp),intent(in)           :: vGn_m                     ! van Genutchen "m" parameter (-)
 real(dp),intent(in)           :: theta_sat                 ! soil porosity (-)
 real(dp),intent(in)           :: theta_res                 ! soil residual volumetric water content (-)
 ! input: saturated hydraulic conductivity
 real(dp),intent(in)           :: mLayerSatHydCond(:)       ! saturated hydraulic conductivity in each layer (m s-1)
 ! output: ejected water flux and derivative
 real(dp),intent(out)          :: mLayerEjectWater(:)       ! water ejected because pore volume is filled (m s-1)
 real(dp),intent(out)          :: mLayerEjectWaterDeriv(:)  ! derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
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
 integer(i4b)                  :: iSoil                     ! index of soil layers
 ! initialize error control
 err=0; message="ejectWater/"

 ! define threshold in the super-saturation function (-)
 supersatThresh = supersatScale * log(1._dp/xMatch - 1._dp) + xMatch

 ! define minimum value for calculations
 fracMin = -supersatScale*log(1._dp/fSmall - 1._dp) + supersatThresh

 ! loop through soil layers
 do iSoil=1,nLevels
  ! only process if moisture-based form of RE, or if ice is present in the mixed form of RE
  if(ixRichards==moisture .or. (ixRichards==mixdform .and. mLayerVolFracIceTrial(iSoil)>0._dp))then
   ! calculate the fraction of pore space filled with liquid water and ice (-)
   select case(ixRichards)
    case(moisture); fracCap = (mLayerVolFracLiqTrial(iSoil) + mLayerVolFracIceTrial(iSoil))/theta_sat
    case(mixdform); fracCap = (volFracLiq(mLayerMatricHeadTrial(iSoil),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) + mLayerVolFracIceTrial(iSoil))/theta_sat
    case default; err=20; message=trim(message)//"unable to identify option for Richards' equation"; return
   endselect
   ! check if the fractional capacity is greater than the minimum value
   if(fracCap > fracMin)then
    ! calculate the lateral flux (m s-1)
    expFunc = exp((supersatThresh - fracCap)/supersatScale)
    mLayerEjectWater(iSoil) = mLayerSatHydCond(iSoil)/(1._dp + expFunc)
    ! calculate the derivative (m s-1 [moisture form] or s-1 [mixed form])
    select case(ixRichards)
     case(moisture)
      if(dMethod==analytical)then
       mLayerEjectWaterDeriv(iSoil) = (expFunc/(theta_sat*supersatScale)) * (1._dp + expFunc)**(-2._dp) * mLayerSatHydCond(iSoil)
      else
       fracCap = (mLayerVolFracLiqTrial(iSoil)+dx + mLayerVolFracIceTrial(iSoil))/theta_sat
       expTemp = exp((supersatThresh - fracCap)/supersatScale)
       mLayerEjectWaterDeriv(iSoil) = (mLayerSatHydCond(iSoil)/(1._dp + expTemp) -  mLayerEjectWater(iSoil))/dx
      endif
     case(mixdform)
      if(dMethod==analytical)then
       mLayerEjectWaterDeriv(iSoil) = mLayerdTheta_dPsi(iSoil) * (expFunc/(theta_sat*supersatScale)) * (1._dp + expFunc)**(-2._dp) * mLayerSatHydCond(iSoil)
      else
       fracCap = (volFracLiq(mLayerMatricHeadTrial(iSoil)+dx,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) + mLayerVolFracIceTrial(iSoil))/theta_sat
       expTemp  = exp((supersatThresh - fracCap)/supersatScale)
       mLayerEjectWaterDeriv(iSoil) = (mLayerSatHydCond(iSoil)/(1._dp + expTemp) -  mLayerEjectWater(iSoil))/dx
      endif
     case default; err=20; message=trim(message)//"unable to identify option for Richards' equation"; return
    end select
   ! (fraction of pore space filled with liquid water and ice is less than the minimum value required for lateral flow)
   else
    mLayerEjectWater(iSoil)      = 0._dp
    mLayerEjectWaterDeriv(iSoil) = 0._dp
   endif
  ! (mixed form of Richards' equation when no ice is present)
  else
   mLayerEjectWater(iSoil)      = 0._dp
   mLayerEjectWaterDeriv(iSoil) = 0._dp
  endif
  !write(*,'(a,i4,1x,10(f20.9,1x))') 'iSoil, mLayerEjectWater(iSoil) = ', iSoil, mLayerEjectWater(iSoil)
 enddo  ! (loop through soil layers)
 end subroutine ejectWater


 ! ************************************************************************************************
 ! private subroutine: compute liquid fluxes at layer interfaces
 ! ************************************************************************************************
 subroutine iLayer_liq(&
                       ! input: model control
                       ixRichards,            & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                       bc_upper,bc_lower,     & ! intent(in): index defining the type of boundary conditions (neumann or diriclet)
                       ! input: model coordinate variables
                       mLayerDepth,           & ! intent(in): depth of the layer (m)
                       mLayerHeight,          & ! intent(in): height of the layer mid-point (m)
                       ! input: state variables
                       mLayerMatricHeadTrial, & ! intent(in): matric head in each layer (m)
                       mLayerVolFracIceTrial, & ! intent(in): volumetric ice content in each layer (-)
                       mLayerVolFracLiqTrial, & ! intent(in): volumetric liquid water content (-)
                       ! input: boundary conditions for matric head
                       lowerBoundHead,        & ! intent(in): lower boundary condition (m)
                       upperBoundHead,        & ! intent(in): upper boundary condition (m)
                       ! input: boundary conditions for volumetric liqquid water content
                       lowerBoundTheta,       & ! intent(in): lower boundary condition (-)
                       upperBoundTheta,       & ! intent(in): upper boundary condition (-)
                       ! input: flux at the upper boundary
                       scalarRainPlusMelt,    & ! intent(in): rain plus melt (m s-1)
                       ! input: height of the lower boundary and the water table
                       lowerBoundaryHeight,   & ! intent(in): height of the lower boundary (m)
                       scalarWaterTableDepth, & ! intent(in): water table depth (m)
                       ! input: transmittance
                       iLayerSatHydCond,      & ! intent(in): saturated hydraulic conductivity at the layer interfaces (m s-1)
                       iLayerHydCond,         & ! intent(in): hydraulic conductivity at layer interfaces (m s-1)
                       iLayerDiffuse,         & ! intent(in): hydraulic diffusivity at layer interfaces (m2 s-1)
                       iceImpedeFac,          & ! intent(in): ice impedence factor in each layer (-)
                       ! input: soil parameters
                       vGn_alpha,             & ! intent(in): van Genutchen "alpha" parameter (m-1)
                       vGn_n,                 & ! intent(in): van Genutchen "n" parameter (-)
                       VGn_m,                 & ! intent(in): van Genutchen "m" parameter (-)
                       theta_sat,             & ! intent(in): soil porosity (-)
                       theta_res,             & ! intent(in): soil residual volumetric water content (-)
                       kAnisotropic,          & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                       zScale_TOPMODEL,       & ! intent(in): scale factor for TOPMODEL-ish baseflow parameterization (m)
                       bpar_VIC,              & ! intent(in): b-parameter in the VIC surface runoff parameterization (-)
                       ! output: fluxes at layer interfaces and surface runoff
                       iLayerLiqFlux,         & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                       scalarSurfaceRunoff,   & ! intent(out): surface runoff (m s-1)
                       ! output: error control
                       err,message)             ! intent(out): error control
 USE soil_utils_module,only:matricHead      ! compute matric head (m)
 USE soil_utils_module,only:satArea_matric  ! compute saturated area as a function of matric head (-)
 USE soil_utils_module,only:satArea_liquid  ! compute saturated area as a function of volumetric liquid water content (-)
 ! compute liquid flux at layer interfaces
 implicit none
 ! input: model control
 integer(i4b),intent(in)       :: ixRichards                ! index defining the option for Richards' equation (moisture or mixdform)
 integer(i4b),intent(in)       :: bc_upper,bc_lower         ! index defining the type of boundary conditions
 ! model coordinate variables
 real(dp),intent(in)           :: mLayerDepth(:)            ! depth of the layer (m)
 real(dp),intent(in)           :: mLayerHeight(:)           ! height of the layer mid-point (m)
 ! input: state and diagnostic variables
 real(dp),intent(in)           :: mLayerMatricHeadTrial(:)  ! matric head in each layer (m)
 real(dp),intent(in)           :: mLayerVolFracIceTrial(:)  ! volumetric ice content in each layer (-)
 real(dp),intent(in)           :: mLayerVolFracLiqTrial(:)  ! volumetric liquid water content in each layer (-)
 ! input: diriclet boundary conditions for matric head
 real(dp),intent(in)           :: upperBoundHead            ! upper boundary condition for matric head (m)
 real(dp),intent(in)           :: lowerBoundHead            ! lower boundary condition for matric head (m)
 ! input: diriclet boundary conditions for volumetric liquid water content
 real(dp),intent(in)           :: upperBoundTheta           ! upper boundary condition for volumetric liquid water content (-)
 real(dp),intent(in)           :: lowerBoundTheta           ! lower boundary condition for volumetric liquid water content (-)
 ! input: flux at the upper boundary
 real(dp),intent(in)           :: scalarRainPlusMelt        ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
 ! input: height of the lower boundary and the water table
 real(dp),intent(in)           :: lowerBoundaryHeight       ! height of the lower boundary (m)
 real(dp),intent(in)           :: scalarWaterTableDepth     ! water table depth (m)
 ! input: transmittance
 real(dp),intent(in)           :: iLayerSatHydCond(0:)      ! saturated hydraulic conductivity at the layer interfaces (m s-1)
 real(dp),intent(in)           :: iLayerHydCond(0:)         ! hydraulic conductivity at layer interfaces (m s-1)
 real(dp),intent(in)           :: iLayerDiffuse(0:)         ! hydraulic diffusivity at layer interfaces (m s-1)
 real(dp),intent(in)           :: iceImpedeFac(:)           ! ice impedence factor in each layer (-)
 ! input: soil parameters
 real(dp),intent(in)           :: vGn_alpha                 ! van Genutchen "alpha" parameter (m-1)
 real(dp),intent(in)           :: vGn_n                     ! van Genutchen "n" parameter (-)
 real(dp),intent(in)           :: vGn_m                     ! van Genutchen "m" parameter (-)
 real(dp),intent(in)           :: theta_sat                 ! soil porosity (-)
 real(dp),intent(in)           :: theta_res                 ! soil residual volumetric water content (-)
 real(dp),intent(in)           :: kAnisotropic              ! anisotropy factor for lateral hydraulic conductivity (-)
 real(dp),intent(in)           :: zScale_TOPMODEL           ! scale factor for TOPMODEL-ish baseflow parameterization (m)
 real(dp),intent(in)           :: bpar_VIC                  ! b-parameter in the VIC surface runoff parameterization (-)
 ! output: fluxes at layer interfaces and surface runoff
 real(dp),intent(out)          :: iLayerLiqFlux(0:)         ! liquid fluxes at layer interfaces (m s-1)
 real(dp),intent(out)          :: scalarSurfaceRunoff       ! surface runoff (m s-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                       ! error code
 character(*),intent(out)      :: message                   ! error message
 ! local variables
 real(dp)                      :: satArea                   ! saturated area (-)
 real(dp)                      :: infUnsat                  ! infiltration over unsaturated areas (m s-1)
 real(dp)                      :: psiBot                    ! matric head in the lowest layer (m)
 real(dp)                      :: cflux                     ! capillary flux (m s-1)
 integer(i4b)                  :: iSoil                     ! index of soil layer
 real(dp)                      :: zWater                    ! depth to the water table (m)
 ! initialize error control
 err=0; message="iLayer_liq/"

 ! -------------------------------------------------------------------------------------------------------------------------------
 ! 1) compute fluxes at the upper boundary -- positive downwards
 ! -------------------------------------------------------------------------------------------------------------------------------
 select case(bc_upper)

 ! *****
 ! head condition
 case(prescribedHead)
  scalarSurfaceRunoff = 0._dp
  select case(ixRichards)
   case(moisture); cflux = -iLayerDiffuse(0)*(mLayervolFracLiqTrial(1) - upperBoundTheta) / (mLayerDepth(1)*0.5_dp)
   case(mixdform); cflux = -iLayerHydCond(0)*(mLayerMatricHeadTrial(1) - upperBoundHead) / (mLayerDepth(1)*0.5_dp)
  end select
  iLayerLiqFlux(0) = cflux + iLayerHydCond(0)

 ! *****
 ! flux condition
 case(liquidFlux)
  ! compute the surface runoff (m s-1)
  select case(ixRichards)
   case(moisture); satArea = satArea_liquid(mLayerVolFracLiqTrial(1),mLayerVolFracIceTrial(1),theta_sat,bpar_VIC)
   case(mixdform); satArea = satArea_matric(mLayerMatricHeadTrial(1),mLayerVolFracIceTrial(1),theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m,bpar_VIC)
  endselect
  infUnsat= min(iLayerSatHydCond(0),scalarRainPlusMelt) * iceImpedeFac(1)   ! infiltration over unsaturated areas (m s-1)
  scalarSurfaceRunoff = satArea*scalarRainPlusMelt + (1._dp - satArea)*(scalarRainPlusMelt - infUnsat)
  ! compute the flux at the upper boundary
  iLayerLiqFlux(0) = scalarRainPlusMelt - scalarSurfaceRunoff

 ! ***** error check
 case default
  err=20; message=trim(message)//'unknown upper boundary condition for soil hydrology'; return
 endselect

 ! -------------------------------------------------------------------------------------------------------------------------------
 ! 2) compute fluxes at the lower boundary -- positive downwards
 ! -------------------------------------------------------------------------------------------------------------------------------
 select case(bc_lower)

 ! ***** coupled to groundwater
 case(groundwaterCouple)
  if(ixRichards/=moisture)then; err=20; message=trim(message)//"require moisture-based form of Richards' eqn when using coupled gw"; return; endif
  ! compute the capillary flux
  psiBot = matricHead(mLayervolFracLiqTrial(nLevels),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
  cflux  = -iLayerHydCond(nLevels)*(0._dp - psiBot) / (scalarWaterTableDepth - mLayerHeight(nLevels))
  ! compute the total flux
  iLayerLiqFlux(nLevels) = cflux + iLayerHydCond(nLevels)

 ! ***** head condition
 case(prescribedHead)
  select case(ixRichards)
   case(moisture); cflux = -iLayerDiffuse(nLevels)*(lowerBoundTheta - mLayervolFracLiqTrial(nLevels)) / (mLayerDepth(nLevels)*0.5_dp)
   case(mixdform); cflux = -iLayerHydCond(nLevels)*(lowerBoundHead  - mLayerMatricHeadTrial(nLevels)) / (mLayerDepth(nLevels)*0.5_dp)
  endselect
  iLayerLiqFlux(nLevels) = cflux + iLayerHydCond(nLevels)

 ! *****
 ! function of matric head in the bottom layer
 case(funcBottomHead)
  select case(ixRichards)
   case(moisture); zWater = mLayerHeight(nLevels) - matricHead(mLayerVolFracLiqTrial(nLevels),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   case(mixdform); zWater = mLayerHeight(nLevels) - mLayerMatricHeadTrial(nLevels)
  endselect
  iLayerLiqFlux(nLevels) = kAnisotropic*iLayerSatHydCond(0) * exp(-zWater/zScale_TOPMODEL)

 ! *****
 ! free drainage
 case(freeDrainage); iLayerLiqFlux(nLevels) = iLayerHydCond(nLevels)

 ! ***** error check
 case default
  err=20; message=trim(message)//'unknown lower boundary condition for soil hydrology'; return
 endselect ! (type of boundary condition)

 ! -------------------------------------------------------------------------------------------------------------------------------
 ! 3) compute fluxes within the domain -- positive downwards
 ! -------------------------------------------------------------------------------------------------------------------------------
 do iSoil=1,nLevels-1 ! iSoil=0 is the upper boundary, so flux at iSoil=1 is bottom of the top layer
  ! compute the capillary flux (negative sign means positive downwards)
  select case(ixRichards)
   case(moisture); cflux = -iLayerDiffuse(iSoil)*(mLayervolFracLiqTrial(iSoil+1) - mLayervolFracLiqTrial(iSoil)) / &
                                                 (mLayerHeight(iSoil+1) - mLayerHeight(iSoil))
   case(mixdform); cflux = -iLayerHydCond(iSoil)*(mLayerMatricHeadTrial(iSoil+1) - mLayerMatricHeadTrial(iSoil)) / &
                                                 (mLayerHeight(iSoil+1) - mLayerHeight(iSoil))
  end select
  ! compute the total flux (add gravity flux, positive downwards)
  iLayerLiqFlux(iSoil) = cflux + iLayerHydCond(iSoil)
 end do ! looping through layers within the domain

 !print*, 'iLayerLiqFlux(0:1) = ', iLayerLiqFlux(0:1)

 end subroutine iLayer_liq

 
  
 ! ************************************************************************************************
 ! private subroutine: compute derivative in fluxes at layer interfaces w.r.t.
 !                      volumetric liquid water content in the layer above and the layer below
 ! ************************************************************************************************
 subroutine groundwatr(&
                       dt,                                              & ! input:  time step (seconds) 
                       ! input (state variables)
                       scalarAquiferStorage,                            & ! input: aquifer storage at the start of the time step (m)
                       scalarAquiferStorageTrial,                       & ! input: trial value of aquifer storage (m)
                       ! input (fluxes)
                       scalarAquiferRcharge,                            & ! input: aquifer recharge (m s-1)
                       scalarAquiferTranspire,                          & ! input: aquifer transpiration (m s-1)
                       ! input (available depth)
                       mLayerVolFracLiqIter,                            & ! input: volumetric fraction of liquid water after itertations (-)
                       mLayerVolFracIceIter,                            & ! input: volumetric fraction of ice after itertations (-)
                       iLayerHeight,                                    & ! input: height of each interface (m), NOTE: include height of soil (iLayer=nSnow)
                       mLayerDepth,                                     & ! input: depth of each soil layer (m)
                       ! input (diagnostic variables)
                       k_surf,                                          & ! input: hydraulic conductivity at the surface (m s-1)
                       ! input (parameters)
                       theta_sat,                                       & ! input: soil porosity (-)
                       specificYield,                                   & ! input: specific yield (-)
                       kAnisotropic,                                    & ! input: anisotropy factor for lateral hydraulic conductivity (-)
                       zScale_TOPMODEL,                                 & ! input: scale factor for TOPMODEL-ish baseflow parameterization (m)
                       ! output (states)
                       scalarAquiferStorageNew,                         & ! output: aquifer storage at the end of the time step (m)
                       ! output (diagnostic variables and error control)
                       scalarWaterTableDepth,                           & ! output: water table depth at the start/end of the time step (m)
                       scalarAquiferBaseflow,                           & ! output: baseflow from the aquifer (m s-1)
                       err,message)                                       ! output: error control
 ! compute change in aquifer storage over the time step
 implicit none
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! input (time step)
 real(dp),intent(in)                 :: dt                        ! time step (seconds)
 ! input (state variables)
 real(dp),intent(in)                 :: scalarAquiferStorage      ! aquifer storage at the start of the time step (m)
 real(dp),intent(in)                 :: scalarAquiferStorageTrial ! trial value of aquifer storage (m)
 ! input (fluxes)
 real(dp),intent(in)                 :: scalarAquiferRcharge      ! aquifer recharge averaged over the time step (m s-1)
 real(dp),intent(in)                 :: scalarAquiferTranspire    ! aquifer transpiration averaged over the time step (m s-1)
 ! input (available depth)
 real(dp),intent(in)                 :: mLayerVolFracLiqIter(:)   ! volumetric fraction of liquid water after itertations (-)
 real(dp),intent(in)                 :: mLayerVolFracIceIter(:)   ! volumetric fraction of ice after itertations (-)
 real(dp),intent(in)                 :: iLayerHeight(0:)          ! height at each layer interface (m)
 real(dp),intent(in)                 :: mLayerDepth(:)            ! depth of each soil layer (m)
 ! input (diagnostic variables)
 real(dp),intent(in)                 :: k_surf                    ! hydraulic conductivity at the surface (m s-1) 
 ! input (parameters)
 real(dp),intent(in)                 :: theta_sat                 ! total porosity (-)
 real(dp),intent(in)                 :: specificYield             ! fraction of water volume drained by gravity in an unconfined aquifer (-)
 real(dp),intent(in)                 :: kAnisotropic              ! anisotropy factor for lateral hydraulic conductivity (-)
 real(dp),intent(in)                 :: zScale_TOPMODEL           ! scale factor for TOPMODEL-ish baseflow parameterization (m)
 ! output (states)
 real(dp),intent(out)                :: scalarAquiferStorageNew   ! aquifer storage at the end of the time step (m)
 ! output (diagnostic variables and error control)
 real(dp),intent(out)                :: scalarWaterTableDepth     ! water table depth at the start/end of the time step (m)
 real(dp),intent(out)                :: scalarAquiferBaseflow     ! baseflow from the aquifer (m s-1)
 integer(i4b),intent(out)            :: err                       ! error code
 character(*),intent(out)            :: message                   ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! algorithmic control parameters
 integer(i4b),parameter              :: maxiter=10                ! maximum number of iterations
 real(dp),parameter                  :: incTol= 1.e-6_dp          ! convergence tolerance for the iteration increment (m)
 real(dp),parameter                  :: resTol= 1.e-6_dp          ! convergence tolerance for the residual (m)
 ! define local variables
 character(len=256)                  :: cmessage                  ! error message for downwind routine
 integer(i4b)                        :: iter                      ! iteration index
 integer(i4b)                        :: iSoil                     ! index of each soil layer
 integer(i4b)                        :: nSoil                     ! number of soil layers
 integer(i4b)                        :: nUnsat                    ! number of unsaturated layers
 integer(i4b)                        :: nUnsat_init               ! number of unsaturated layers at the start of the iterations
 real(dp)                            :: baseflowMax               ! maximum baseflow rate from the aquifer (m s-1)
 real(dp),dimension(0:size(mLayerVolFracLiqIter)) :: waterReq2FillPore  ! water required to fill pore space (m)
 real(dp)                            :: drainablePorosity         ! drainable porosity (-)
 real(dp)                            :: aquiferStorageTrial       ! trial value of aquifer storage (m) -- local variable to preserve input
 real(dp)                            :: zWaterTrial               ! depth of the water table (m)
 real(dp)                            :: dBaseflow_dStorage        ! derivative in the baseflow term w.r.t. aquifer storage (s-1)
 real(dp)                            :: netFlux                   ! recharge minus baseflow (m s-1)
 real(dp)                            :: residual                  ! residual in aquifer storage (m)
 real(dp)                            :: dStorage                  ! iteration increment for aquifer storage (m)
 integer(i4b)                        :: idiff                     ! trial integers used to find the bracket
 integer(i4b),parameter              :: maxdiff=100               ! maximum trial values to find the bracket
 real(dp),parameter                  :: xdiff=0.00001_dp          ! increments used to find the bracket
 real(dp)                            :: xbeg,xend                 ! lower and upper points in the bracket
 real(dp)                            :: rbeg,rend                 ! residuals for the lower and upper points in the bracket
 integer(i4b)                        :: itest                     ! trial integers used to test along a line
 integer(i4b),parameter              :: ntest=100000              ! number of values to test along a line
 real(dp)                            :: xinc                      ! increment in storage to test along a line
 real(dp)                            :: xtry,rtry                 ! trial value and residual
 real(dp)                            :: zWaterTrial_init          ! initial value of the water table depth

 ! initialize error control
 err=0; message="groundwatr/"

 ! compute the maximum baseflow rate
 baseflowMax = kAnisotropic*k_surf

 ! compute the number of unsaturated layers, the drainable porosity, and the depth to the water table
 call waterTablePosition(mLayerVolFracLiqIter,       & ! intent(in):  volumetric liquid water content in each soil layer (-)
                         mLayerVolFracIceIter,       & ! intent(in):  volumetric ice content in each soil layer (-)
                         aquiferStorageTrial,        & ! intent(in):  aquifer storage (m)
                         waterReq2FillPore,          & ! intent(in):  water required to fill pore space up to the level of each layer interface (m(
                         iLayerHeight,               & ! intent(in):  height of each interface (m)
                         mLayerDepth,                & ! intent(in):  depth of each soil layer (m)
                         theta_sat,                  & ! intent(in):  soil porosity (-)
                         specificYield,              & ! intent(in):  fraction of water volume drained by gravity in an unconfined aquifer (-)
                         nUnsat_init,                & ! intent(out): number of unsaturated layers at the start of the iterations
                         drainablePorosity,          & ! intent(out): drainable porosity (-)
                         zWaterTrial,                & ! intent(out): water table depth (m)
                         err,cmessage)                 ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! -----
 ! calculate the residual...
 ! -------------------------
 ! calculate the baseflow term (m s-1)
 scalarAquiferBaseflow = baseflowMax*exp(-zWaterTrial/zScale_TOPMODEL)
 ! compute the net flux (m s-1)
 netFlux = (scalarAquiferRcharge + scalarAquiferTranspire) - scalarAquiferBaseflow
 ! compute the residual in aquifer storage (m)
 residual = aquiferStorageTrial - (scalarAquiferStorage + netFlux*dt)




 end subroutine groundwatr


 ! ************************************************************************************************
 ! private subroutine: compute derivative in fluxes at layer interfaces w.r.t.
 !                      volumetric liquid water content in the layer above and the layer below
 ! ************************************************************************************************
 subroutine dq_dVLiquid(&
                        ! input: model control
                        dMethod,               & ! intent(in): method used to compute derivatives (analytical or numerical)
                        bc_upper,bc_lower,     & ! intent(in): index defining the type of boundary conditions
                        ! input: model coordinate variables
                        mLayerDepth,           & ! intent(in): depth of the layer (m)
                        mLayerHeight,          & ! intent(in): height of the layer mid-point (m)
                        ! input: state and diagnostic variables
                        mLayerVolFracLiqTrial, & ! intent(in): volumetric liquid water content (-)
                        mLayerVolFracIceTrial, & ! intent(in): volumetric fraction of ice (-)
                        mLayerdPsi_dTheta,     & ! intent(in): derivative in the soil water characteritic w.r.t. theta (m)
                        ! input: height of the lower boundary and the depth to the water table
                        lowerBoundaryHeight,   & ! intent(in): height of the lower boundary
                        scalarWaterTableDepth, & ! intent(in): water table depth (m)
                        ! input: boundary conditions for volumetric liqquid water content
                        lowerBoundTheta,       & ! intent(in): lower boundary condition (-)
                        upperBoundTheta,       & ! intent(in): upper boundary condition (-)
                        ! input: flux at the upper boundary
                        scalarRainPlusMelt,    & ! intent(in): rain plus melt (m s-1)
                        ! input: transmittance
                        mLayerSatHydCond,      & ! intent(in): saturated hydraulic conductivity in each layer (m s-1)
                        iLayerSatHydCond,      & ! intent(in): saturated hydraulic conductivity at layer interfaces (m s-1)
                        mLayerHydCond,         & ! intent(in): hydraulic conductivity in each layer (m s-1)
                        mLayerDiffuse,         & ! intent(in): hydraulic diffusivity in each layer (m s-1)
                        iLayerHydCond,         & ! intent(in): hydraulic conductivity at layer interfaces (m s-1)
                        iLayerDiffuse,         & ! intent(in): hydraulic diffusivity at layer interfaces (m s-1)
                        iceImpedeFac,          & ! intent(in): ice impedence factor in each layer (-)
                        ! input: soil parameters
                        vGn_alpha,             & ! intent(in): van Genutchen "alpha" parameter (m-1)
                        vGn_n,                 & ! intent(in): van Genutchen "n" parameter (-)
                        VGn_m,                 & ! intent(in): van Genutchen "m" parameter (-)
                        theta_sat,             & ! intent(in): soil porosity (-)
                        theta_res,             & ! intent(in): soil residual volumetric water content (-)
                        kAnisotropic,          & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                        zScale_TOPMODEL,       & ! intent(in): scale factor for TOPMODEL-ish baseflow parameterization (m)
                        bpar_VIC,              & ! intent(in): b-parameter in the VIC surface runoff parameterization (-)
                        f_impede,              & ! intent(in): ice impedence factor (-)
                        ! output
                        dq_dVolLiqAbove,       & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer above (m s-1)
                        dq_dVolLiqBelow,       & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer below (m s-1)
                        err,message)             ! intent(out): error control
 ! compute derivative in fluxes at layer interfaces w.r.t. volumetric liquid water content head in the layer above and the layer below
 USE soil_utils_module,only:darcyFlux_liquid ! compute Darcy's flux
 USE soil_utils_module,only:matricHead       ! compute matric head as a function of volumetric liquid water content (m)
 USE soil_utils_module,only:hydCond_liq      ! compute hydraulic conductivity as a function of volumetric liquid water content (m s-1)
 USE soil_utils_module,only:dHydCond_dLiq    ! compute derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
 USE soil_utils_module,only:satArea_liquid   ! compute saturated area as a function of volumetric liquid water content (-)
 USE soil_utils_module,only:dPsi_dTheta      ! compute derivative in soil water characteristic (m)
 USE soil_utils_module,only:dPsi_dTheta2     ! compute derivative in dPsi_dTheta (m)
 implicit none
 ! input: model control
 integer(i4b),intent(in)       :: dMethod                    ! method used to compute derivatives (analytical or numerical)
 integer(i4b),intent(in)       :: bc_upper,bc_lower          ! index defining the type of boundary conditions (neumann or diriclet)
 ! model coordinate variables
 real(dp),intent(in)           :: mLayerDepth(:)             ! depth of the layer (m)
 real(dp),intent(in),target    :: mLayerHeight(:)            ! height of the layer mid-point (m)
 ! input: state and diagnostic variables
 real(dp),intent(in),target    :: mLayerVolFracLiqTrial(:)   ! volumetric liquid water content (-)
 real(dp),intent(in),target    :: mLayerVolFracIceTrial(:)   ! volumetric fraction of ice in each layer (-)
 real(dp),intent(in)           :: mLayerdPsi_dTheta(:)       ! derivative in the soil water characteritic w.r.t. theta (m)
 ! input: height of the lower boundary and the depth to the water table
 real(dp),intent(in)           :: lowerBoundaryHeight        ! height of the lower boundary (m)
 real(dp),intent(in)           :: scalarWaterTableDepth      ! water table depth (m)
 ! input: diriclet boundary conditions for volumetric liquid water content
 real(dp),intent(in)           :: upperBoundTheta            ! upper boundary condition for volumetric liquid water content (-)
 real(dp),intent(in)           :: lowerBoundTheta            ! lower boundary condition for volumetric liquid water content (-)
 ! input: flux at the upper boundary
 real(dp),intent(in)           :: scalarRainPlusMelt         ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
 ! input: transmittance
 real(dp),intent(in),target    :: mLayerSatHydCond(:)        ! saturated hydraulic conductivity at the mid-point of each layer (m s-1)
 real(dp),intent(in)           :: iLayerSatHydCond(0:)       ! saturated hydraulic conductivity at the interface of each layer (m s-1)
 real(dp),intent(in)           :: mLayerHydCond(:)           ! hydraulic conductivity in each layer (m s-1)
 real(dp),intent(in)           :: mLayerDiffuse(:)           ! hydraulic diffusivity in each layer (m s-1)
 real(dp),intent(in)           :: iLayerHydCond(0:)          ! hydraulic conductivity at layer interfaces (m s-1)
 real(dp),intent(in)           :: iLayerDiffuse(0:)          ! hydraulic diffusivity at layer interfaces (m s-1)
 real(dp),intent(in)           :: iceImpedeFac(:)            ! ice impedence factor in each layer (-)
 ! input: soil parameters
 real(dp),intent(in)           :: vGn_alpha                  ! van Genutchen "alpha" parameter (m-1)
 real(dp),intent(in)           :: vGn_n                      ! van Genutchen "n" parameter (-)
 real(dp),intent(in)           :: vGn_m                      ! van Genutchen "m" parameter (-)
 real(dp),intent(in)           :: theta_sat                  ! soil porosity (-)
 real(dp),intent(in)           :: theta_res                  ! soil residual volumetric water content (-)
 real(dp),intent(in)           :: kAnisotropic               ! anisotropy factor for lateral hydraulic conductivity (-)
 real(dp),intent(in)           :: zScale_TOPMODEL            ! scale factor for TOPMODEL-ish baseflow parameterization (m)
 real(dp),intent(in)           :: bpar_VIC                   ! b-parameter in the VIC surface runoff parameterization (-)
 real(dp),intent(in)           :: f_impede                   ! ice impedence factor (-)
 ! output
 real(dp),intent(out)          :: dq_dVolLiqAbove(0:)        ! derivatives in the flux w.r.t. volumetric liquid water content in the layer above (m s-1)
 real(dp),intent(out)          :: dq_dVolLiqBelow(0:)        ! derivatives in the flux w.r.t. volumetric liquid water content in the layer below (m s-1)
 integer(i4b),intent(out)      :: err                        ! error code
 character(*),intent(out)      :: message                    ! error message
 ! local variables
 real(dp),dimension(size(mLayerVolFracLiqTrial))  :: dHydCond_dVolLiq  ! derivative in hydraulic conductivity at layer mid-points w.r.t. volumetric liquid water content
 real(dp),dimension(size(mLayerVolFracLiqTrial))  :: dDiffuse_dVolLiq  ! derivative in hydraulic diffusivity at layer mid-points w.r.t. volumetric liquid water content
 real(dp)                      :: dHydCondIface_dVolLiqAbove ! derivative in hydraulic conductivity at layer interface w.r.t. volumetric liquid water content in layer above
 real(dp)                      :: dHydCondIface_dVolLiqBelow ! derivative in hydraulic conductivity at layer interface w.r.t. volumetric liquid water content in layer below
 real(dp)                      :: dDiffuseIface_dVolLiqAbove ! derivative in hydraulic diffusivity at layer interface w.r.t. volumetric liquid water content in layer above
 real(dp)                      :: dDiffuseIface_dVolLiqBelow ! derivative in hydraulic diffusivity at layer interface w.r.t. volumetric liquid water content in layer below
 real(dp)                      :: dPsi_dTheta2a              ! derivative in dPsi_dTheta (analytical)
 real(dp)                      :: dPsi_dTheta2n              ! derivative in dPsi_dTheta (numerical)
 real(dp)                      :: fracCap                    ! fraction of pore space filled with liquid water and ice (-)
 real(dp)                      :: fInfArea                   ! area of the landscape where water infiltrates (-)
 real(dp)                      :: dInfArea                   ! derivative in the infiltrating area w.r.t. matric head (m-1)
 real(dp)                      :: fInfRate                   ! infiltration rate at the surface (m s-1)
 real(dp)                      :: dInfRate                   ! derivative in the infiltration rate w.r.t matric head (s-1)
 real(dp)                      :: func0,func1                ! function evaluations used to compute numerical derivatives (general)
 real(dp)                      :: bottomHead                 ! matric head in the lowest soil layer (m)
 real(dp)                      :: dHeight                    ! height difference between water table and bottom soil layer (m) 
 real(dp)                      :: spat0,spat1                ! spatial operators (-)
 real(dp)                      :: dLiq,dz                    ! spatial differenes in volumetric liquid water content and height
 real(dp)                      :: hCnd0,hCnd1                ! hydraulic conductivity in a given layer (m s-1)
 real(dp),pointer              :: iceAbove,iceBelow          ! volumetric ice content in the layer above and the layer below
 real(dp),pointer              :: liqAbove,liqBelow          ! volumetric liquid water content in the layer above and the layer below
 real(dp),pointer              :: sHcAbove,sHcBelow          ! saturated hydraulic conductivity in the layer above and the layer below
 real(dp),pointer              :: zAbove,zBelow              ! height of the layer mid-point in the layer above and the layer below
 real(dp)                      :: flux0,flux1,flux2          ! used to test the numerical derivatives
 integer(i4b)                  :: iLayer                     ! index of layer
 ! initialize error control
 err=0; message="dq_dVLiquid/"

 ! -----------------------------------------------------------------------------------------------------------------------------
 ! ***** compute the derivative in hydraulic conductivity (m s-1) and diffusivity (m2 s-1) at layer mid-points
 !        w.r.t. volumetric liquid water content (s-1)
 ! -----------------------------------------------------------------------------------------------------------------------------
 do iLayer=1,nLevels
  select case(dMethod)
  case(analytical)
   ! compute derivative in hydraulic conductivity (m s-1)
   dHydCond_dVolLiq(iLayer) = dHydCond_dLiq(mLayerVolFracLiqTrial(iLayer),mLayerSatHydCond(iLayer),theta_res,theta_sat,vGn_m,iceImpedeFac(iLayer),.true.)  ! analytical
   ! compute derivative in dPsi_dTheta (m)
   dPsi_dTheta2a = dPsi_dTheta2(mLayerVolFracLiqTrial(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m,.true.)   ! analytical
   !dPsi_dTheta2n = dPsi_dTheta2(mLayerVolFracLiqTrial(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m,.false.)  ! numerical
   !print*, iLayer, dPsi_dTheta2a, dPsi_dTheta2n
   ! compute derivative in hydraulic diffusivity (m2 s-1)
   dDiffuse_dVolLiq(iLayer) = dHydCond_dVolLiq(iLayer)*mLayerdPsi_dTheta(iLayer) + mLayerHydCond(iLayer)*dPsi_dTheta2a
  case(numerical)
   ! compute derivative in hydraulic conductivity (m s-1)
   dHydCond_dVolLiq(iLayer) = dHydCond_dLiq(mLayerVolFracLiqTrial(iLayer),mLayerSatHydCond(iLayer),theta_res,theta_sat,vGn_m,iceImpedeFac(iLayer),.false.) ! numerical
   ! compute derivative in hydraulic diffusivity (m2 s-1)
   func0 = hydCond_liq(mLayerVolFracLiqTrial(iLayer),   mLayerSatHydCond(iLayer),theta_res,theta_sat,vGn_m) * dPsi_dTheta(mLayerVolFracLiqTrial(iLayer),   vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   func1 = hydCond_liq(mLayerVolFracLiqTrial(iLayer)+dx,mLayerSatHydCond(iLayer),theta_res,theta_sat,vGn_m) * dPsi_dTheta(mLayerVolFracLiqTrial(iLayer)+dx,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   dDiffuse_dVolLiq(iLayer) = (func1 - func0)/dx
  case default
   err=10; message=trim(message)//"unknown option to compute derivative in hydraulic conductivity at layer mid-points"; return
  end select
 end do

 ! compute the derivative in flux at layer interfaces w.r.t. volumetric liquid water content
 do iLayer=0,nLevels  ! loop through interfaces

  ! -----------------------------------------------------------------------------------------------------------------------------
  ! ***** the upper boundary
  ! -----------------------------------------------------------------------------------------------------------------------------
  if(iLayer==0)then  ! (upper boundary)

   dq_dVolLiqAbove(iLayer) = valueMissing  ! don't expect this to be used, so deliberately set to a ridiculous value to cause problems

   select case(bc_upper)

   ! *****
   ! head boundary
   case(prescribedHead)      ! head boundary
    ! derivatives in the flux w.r.t. volumetric liquid water content (NOTE: hydraulic diffusivity and conductivity are constant over the step)
    if(dMethod==analytical)then
     dq_dVolLiqBelow(iLayer) = -iLayerDiffuse(0)/(mLayerDepth(1)/2._dp)
    ! compute numerical derivatives
    else
     flux0 = -iLayerDiffuse(0)*( mLayerVolFracLiqTrial(1)     - upperBoundTheta) / (mLayerDepth(1)*0.5_dp) + iLayerHydCond(0)
     flux1 = -iLayerDiffuse(0)*((mLayerVolFracLiqTrial(1)+dx) - upperBoundTheta) / (mLayerDepth(1)*0.5_dp) + iLayerHydCond(0)
     dq_dVolLiqBelow(iLayer) = (flux1 - flux0)/dx
    endif

   ! *****
   ! flux boundary
   case(liquidFlux)          ! flux boundary 
    ! * compute analytical derivatives in the flux w.r.t. volumetric liquid water content
    if(dMethod==analytical)then
     ! compute the function and derivative of the infiltrating area
     fracCap  = min((mLayerVolFracLiqTrial(1) + mLayerVolFracIceTrial(1))/theta_sat, 1._dp - abs(dx))  ! (fraction of capacity filled with liquid water and ice)
     fInfArea = (1._dp - fracCap)**bpar_VIC                                         ! area of the landscape where water infiltrates (-)
     dInfArea = (-1._dp/theta_sat) * bpar_VIC*(1._dp - fracCap)**(bpar_VIC - 1._dp) ! derivative in the infiltrating area w.r.t. volumetric liquid water content (-)
     ! compute the function and derivative of the infiltration rate
     fInfRate = min(iLayerSatHydCond(0),scalarRainPlusMelt)*iceImpedeFac(1)
     dInfRate = 0._dp
     ! ...finally, compute analytical derivatives
     ! NOTE, only first term is used currently, but keep in this form to make life simpler once infiltration depends on volumetric liquid water content
     dq_dVolliqBelow(iLayer) = dInfArea*fInfRate + fInfArea*dInfRate
    ! * compute numerical derivatives
    else
     flux0 = (1._dp - satArea_liquid(mLayerVolFracLiqTrial(1),   mLayerVolFracIceTrial(1),theta_sat,bpar_VIC)) * min(iLayerSatHydCond(0),scalarRainPlusMelt)*iceImpedeFac(1)
     flux1 = (1._dp - satArea_liquid(mLayerVolFracLiqTrial(1)+dx,mLayerVolFracIceTrial(1),theta_sat,bpar_VIC)) * min(iLayerSatHydCond(0),scalarRainPlusMelt)*iceImpedeFac(1)
     dq_dVolliqBelow(iLayer) = (flux1 - flux0)/dx
    endif

   ! *****
   ! error check
   case default
    err=20; message=trim(message)//'unknown upper boundary condition'; return
   end select

  ! -----------------------------------------------------------------------------------------------------------------------------
  ! ***** the lower boundary
  ! -----------------------------------------------------------------------------------------------------------------------------
  elseif(iLayer==nLevels)then  ! (lower boundary)

   dq_dVolLiqBelow(iLayer) = valueMissing  ! don't expect this to be used, so deliberately set to a ridiculous value to cause problems

   select case(bc_lower)

   ! *****
   ! coupled to groundwater
   case(groundwaterCouple)      ! coupled with groundwater
    ! compute intermediate variables
    dHeight = scalarWaterTableDepth - mLayerHeight(nLevels) ! height difference (m)
    if(dHeight < 0._dp)then; err=20; message=trim(message)//'water table is above the mid-point of the lowest soil layer'; return; endif
    ! compute analytical derivatives in the flux w.r.t. volumetric liquid water content
    if(dMethod==analytical)then
     bottomHead  = matricHead(mLayerVolFracLiqTrial(nLevels),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) ! matric head (m)
     dq_dVolLiqAbove(iLayer) = (dHydCond_dVolLiq(iLayer)*bottomHead + mLayerHydCond(iLayer)*mLayerdPsi_dTheta(iLayer))/dHeight + dHydCond_dVolLiq(iLayer)
    ! compute numerical derivatives
    else
     ! (spatial operators)
     spat0 = (0._dp - matricHead(mLayerVolFracLiqTrial(nLevels),   vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m))/dHeight
     spat1 = (0._dp - matricHead(mLayerVolFracLiqTrial(nLevels)+dx,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m))/dHeight
     ! (hydraulic conductivity)
     hCnd0 = hydCond_liq(mLayerVolFracLiqTrial(nLevels),   iLayerSatHydCond(nLevels),theta_res,theta_sat,vGn_m)*iceImpedeFac(nLevels)
     hCnd1 = hydCond_liq(mLayerVolFracLiqTrial(nLevels)+dx,iLayerSatHydCond(nLevels),theta_res,theta_sat,vGn_m)*iceImpedeFac(nLevels)
     ! (fluxes)
     flux0 = -hCnd0*spat0 + hCnd0 
     flux1 = -hCnd1*spat1 + hCnd1
     ! (derivative)
     dq_dVolLiqAbove(iLayer) = (flux1 - flux0)/dx
    endif

   ! *****
   ! head boundary
   case(prescribedHead)      ! head boundary
    ! derivatives in the flux w.r.t. volumetric liquid water content
    if(dMethod==analytical)then
     dq_dVolLiqAbove(iLayer) = iLayerDiffuse(nLevels)/(mLayerDepth(nLevels)/2._dp)
    ! compute numerical derivatives
    else
     flux0 = -iLayerDiffuse(nLevels)*(lowerBoundTheta -  mLayerVolFracLiqTrial(nLevels)    ) / (mLayerDepth(nLevels)*0.5_dp) + iLayerHydCond(nLevels)
     flux1 = -iLayerDiffuse(nLevels)*(lowerBoundTheta - (mLayerVolFracLiqTrial(nLevels)+dx)) / (mLayerDepth(nLevels)*0.5_dp) + iLayerHydCond(nLevels)
     dq_dVolLiqAbove(iLayer) = (flux1 - flux0)/dx
    endif

   ! *****
   ! free draining (analytical vs. numerical switch mage when computing dHydCond_dVolLiq)
   case(freeDrainage); dq_dVolLiqAbove(iLayer) = dHydCond_dVolLiq(iLayer)

   ! *****
   ! function of matric head in the bottom soil layer
   case(funcBottomHead)

    bottomHead = matricHead(mLayerVolFracLiqTrial(nLevels),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
    ! compute analytical derivatives
    if(dMethod==analytical)then
     dq_dVolLiqAbove(iLayer) = kAnisotropic*iLayerSatHydCond(0) * mLayerdPsi_dTheta(nLevels)*exp(-(mLayerHeight(nLevels) - bottomHead)/zScale_TOPMODEL)/zScale_TOPMODEL
    ! compute numerical derivarives
    else
     flux0 = kAnisotropic*iLayerSatHydCond(0) * exp(-(mLayerHeight(nLevels) -  bottomHead    )/zScale_TOPMODEL)
     flux1 = kAnisotropic*iLayerSatHydCond(0) * exp(-(mLayerHeight(nLevels) - (bottomHead+dx))/zScale_TOPMODEL)
     dq_dVolLiqAbove(iLayer) = (flux1 - flux0)/dx
    endif

   ! *****
   ! error check
   case default;
    err=20; message=trim(message)//'unknown lower boundary condition'; return
   endselect ; ! (lower boundary condition)

  ! -----------------------------------------------------------------------------------------------------------------------------
  ! ***** internal layers
  ! -----------------------------------------------------------------------------------------------------------------------------
  else

   if(dMethod==analytical)then
    ! derivatives in hydraulic conductivity at the layer interface (m s-1)
    dHydCondIface_dVolLiqAbove = dHydCond_dVolLiq(iLayer)  *mLayerHydCond(iLayer+1) * 0.5_dp*(mLayerHydCond(iLayer)*mLayerHydCond(iLayer+1))**(-0.5_dp)
    dHydCondIface_dVolLiqBelow = dHydCond_dVolLiq(iLayer+1)*mLayerHydCond(iLayer)   * 0.5_dp*(mLayerHydCond(iLayer)*mLayerHydCond(iLayer+1))**(-0.5_dp)
    ! derivatives in hydraulic diffusivity at the layer interface (m2 s-1)
    dDiffuseIface_dVolLiqAbove = dDiffuse_dVolLiq(iLayer)  *mLayerDiffuse(iLayer+1) * 0.5_dp*(mLayerDiffuse(iLayer)*mLayerDiffuse(iLayer+1))**(-0.5_dp)
    dDiffuseIface_dVolLiqBelow = dDiffuse_dVolLiq(iLayer+1)*mLayerDiffuse(iLayer)   * 0.5_dp*(mLayerDiffuse(iLayer)*mLayerDiffuse(iLayer+1))**(-0.5_dp)
    ! spatial differences in volumetric liquid water content and height
    dLiq  = mLayerVolFracLiqTrial(iLayer+1) - mLayerVolFracLiqTrial(iLayer)
    dz    = mLayerHeight(iLayer+1) - mLayerHeight(iLayer)
    ! derivatives in the flux w.r.t. volumetric liquid water content
    dq_dVolLiqAbove(iLayer) = -dDiffuseIface_dVolLiqAbove*dLiq/dz + iLayerDiffuse(iLayer)/dz + dHydCondIface_dVolLiqAbove
    dq_dVolLiqBelow(iLayer) = -dDiffuseIface_dVolLiqBelow*dLiq/dz - iLayerDiffuse(iLayer)/dz + dHydCondIface_dVolLiqBelow
   else
    ! get some short-cut variables (save typing)
    iceAbove => mLayerVolFracIceTrial(iLayer)
    iceBelow => mLayerVolFracIceTrial(iLayer+1)
    liqAbove => mLayerVolFracLiqTrial(iLayer)
    liqBelow => mLayerVolFracLiqTrial(iLayer+1)
    sHcAbove => mLayerSatHydCond(iLayer)
    sHcBelow => mLayerSatHydCond(iLayer+1)
    zAbove   => mLayerHeight(iLayer)
    zBelow   => mLayerHeight(iLayer+1)
    ! compute numerical derivatives -- note, DarcyFlux computes new hydraulic conductivity and diffusivity
    flux0 = darcyFlux_liquid(iceAbove,iceBelow,liqAbove,   liqBelow,   sHcAbove,sHcBelow,zAbove,zBelow,theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m,f_impede)
    flux1 = darcyFlux_liquid(iceAbove,iceBelow,liqAbove+dx,liqBelow,   sHcAbove,sHcBelow,zAbove,zBelow,theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m,f_impede)
    flux2 = darcyFlux_liquid(iceAbove,iceBelow,liqAbove   ,liqBelow+dx,sHcAbove,sHcBelow,zAbove,zBelow,theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m,f_impede)
    dq_dVolLiqAbove(iLayer) = (flux1 - flux0)/dx
    dq_dVolLiqBelow(iLayer) = (flux2 - flux0)/dx
   endif  ! case for the numerical method

  endif  ! type of layer (upper, internal, or lower)

  ! print results
  !if(iLayer<5) write(*,'(2(i4,1x),2(e20.10,1x))') dMethod, iLayer, dq_dVolLiqAbove(iLayer), dq_dVolLiqBelow(iLayer)

 end do  ! looping through layers
 end subroutine dq_dVLiquid

  
 ! ************************************************************************************************
 ! private subroutine: compute derivative in fluxes at layer interfaces w.r.t.
 !                      matric head in the layer above and the layer below
 ! ************************************************************************************************
  subroutine dq_dMatHead(&
                        ! input: model control
                        dMethod,               & ! intent(in): method used to compute derivatives (analytical or numerical)
                        bc_upper,bc_lower,     & ! intent(in): index defining the type of boundary conditions
                        ! input: model coordinate variables
                        mLayerDepth,           & ! intent(in): depth of the layer (m)
                        mLayerHeight,          & ! intent(in): height of the layer mid-point (m)
                        ! input: state variables
                        mLayerMatricHeadTrial, & ! intent(in): matric head (m)
                        mLayerVolFracIceTrial, & ! intent(in): volumetric fraction of ice (-)
                        mLayerdTheta_dPsi,     & ! intent(in): derivative in the soil water characteritic w.r.t. psi (m-1)
                        ! input: boundary conditions for matric head
                        lowerBoundHead,        & ! intent(in): lower boundary condition (m)
                        upperBoundHead,        & ! intent(in): upper boundary condition (m)
                        ! input: flux at the upper boundary
                        scalarRainPlusMelt,    & ! intent(in): rain plus melt (m s-1)
                        ! input: transmittance
                        mLayerSatHydCond,      & ! intent(in): saturated hydraulic conductivity in each layer (m s-1)
                        iLayerSatHydCond,      & ! intent(in): saturated hydraulic conductivity at layer interfaces (m s-1)
                        mLayerHydCond,         & ! intent(in): hydraulic conductivity in each layer (m s-1)
                        iLayerHydCond,         & ! intent(in): hydraulic conductivity at layer interfaces (m s-1)
                        iceImpedeFac,          & ! intent(in): ice impedence factor in each layer (-)
                        ! input: soil parameters
                        vGn_alpha,             & ! intent(in): van Genutchen "alpha" parameter (m-1)
                        vGn_n,                 & ! intent(in): van Genutchen "n" parameter (-)
                        VGn_m,                 & ! intent(in): van Genutchen "m" parameter (-)
                        theta_sat,             & ! intent(in): soil porosity (-)
                        theta_res,             & ! intent(in): soil residual volumetric water content (-)
                        kAnisotropic,          & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                        zScale_TOPMODEL,       & ! intent(in): scale factor for TOPMODEL-ish baseflow parameterization (m)
                        bpar_VIC,              & ! intent(in): b-parameter in the VIC surface runoff parameterization (-)
                        f_impede,              & ! intent(in): ice impedence factor (-)
                        ! output
                        dq_dMatricAbove,       & ! intent(out): derivatives in the flux w.r.t. matric head in the layer above (s-1)
                        dq_dMatricBelow,       & ! intent(out): derivatives in the flux w.r.t. matric head in the layer below (s-1)
                        err,message)             ! intent(out): error control
 ! compute derivative in fluxes at layer interfaces w.r.t. matric head in the layer above and the layer below
 USE soil_utils_module,only:hydCond_psi      ! compute hydraulic conductivity as a function of matric head
 USE soil_utils_module,only:dHydCond_dPsi    ! compute derivative in hydraulic conductivity w.r.t. matric head
 USE soil_utils_module,only:satArea_matric   ! compute saturated area as a function of matric head (-)
 USE soil_utils_module,only:volFracLiq       ! compute volumetric fraction of liquid water (-)
 USE soil_utils_module,only:darcyFlux_matric ! compute Darcy's flux (m s-1)
 implicit none
 ! input: model control
 integer(i4b),intent(in)       :: dMethod                    ! method used to compute derivatives (analytical or numerical)
 integer(i4b),intent(in)       :: bc_upper,bc_lower          ! index defining the type of boundary conditions (neumann or diriclet)
 ! model coordinate variables
 real(dp),intent(in)           :: mLayerDepth(:)             ! depth of the layer (m)
 real(dp),intent(in),target    :: mLayerHeight(:)            ! height of the layer mid-point (m)
 ! input: state variables
 real(dp),intent(in),target    :: mLayerMatricHeadTrial(:)   ! matric head in each layer (m)
 real(dp),intent(in),target    :: mLayerVolFracIceTrial(:)   ! volumetric fraction of ice in each layer (-)
 real(dp),intent(in)           :: mLayerdTheta_dPsi(:)       ! derivative in the soil water characteritic w.r.t. psi (m-1)
 ! input: diriclet boundary conditions for matric head
 real(dp),intent(in)           :: upperBoundHead             ! upper boundary condition for matric head (m)
 real(dp),intent(in)           :: lowerBoundHead             ! lower boundary condition for matric head (m)
 ! input: flux at the upper boundary
 real(dp),intent(in)           :: scalarRainPlusMelt         ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
 ! input: transmittance
 real(dp),intent(in),target    :: mLayerSatHydCond(:)        ! saturated hydraulic conductivity at the mid-point of each layer (m s-1)
 real(dp),intent(in)           :: iLayerSatHydCond(0:)       ! saturated hydraulic conductivity at the interface of each layer (m s-1)
 real(dp),intent(in)           :: mLayerHydCond(:)           ! hydraulic conductivity in each layer (m s-1)
 real(dp),intent(in)           :: iLayerHydCond(0:)          ! hydraulic conductivity at layer interfaces (m s-1)
 real(dp),intent(in)           :: iceImpedeFac(:)            ! ice impedence factor in each layer (-)
 ! input: soil parameters
 real(dp),intent(in)           :: vGn_alpha                  ! van Genutchen "alpha" parameter (m-1)
 real(dp),intent(in)           :: vGn_n                      ! van Genutchen "n" parameter (-)
 real(dp),intent(in)           :: vGn_m                      ! van Genutchen "m" parameter (-)
 real(dp),intent(in)           :: theta_sat                  ! soil porosity (-)
 real(dp),intent(in)           :: theta_res                  ! soil residual volumetric water content (-)
 real(dp),intent(in)           :: kAnisotropic               ! anisotropy factor for lateral hydraulic conductivity (-)
 real(dp),intent(in)           :: zScale_TOPMODEL            ! scale factor for TOPMODEL-ish baseflow parameterization (m)
 real(dp),intent(in)           :: bpar_VIC                   ! b-parameter in the VIC surface runoff parameterization (-)
 real(dp),intent(in)           :: f_impede                   ! ice impedence factor (-)
 ! output
 real(dp),intent(out)          :: dq_dMatricAbove(0:)        ! derivatives in the flux w.r.t. matric head in the layer above (s-1)
 real(dp),intent(out)          :: dq_dMatricBelow(0:)        ! derivatives in the flux w.r.t. matric head in the layer below (s-1)
 integer(i4b),intent(out)      :: err                        ! error code
 character(*),intent(out)      :: message                    ! error message
 ! local variables
 real(dp),dimension(size(mLayerMatricHeadTrial))  :: dHydCond_dMatric    ! derivative in hydraulic conductivity at layer mid-points w.r.t. matric head
 real(dp)                      :: dHydCondIface_dMatricAbove ! derivative in hydraulic conductivity at layer interface w.r.t. matric head in layer above
 real(dp)                      :: dHydCondIface_dMatricBelow ! derivative in hydraulic conductivity at layer interface w.r.t. matric head in layer below
 real(dp)                      :: vFracLiq                   ! volumetric fraction of liquid water (-)
 real(dp)                      :: fracCap                    ! fraction of pore space filled with liquid water and ice (-)
 real(dp)                      :: fInfArea                   ! area of the landscape where water infiltrates (-)
 real(dp)                      :: dInfArea                   ! derivative in the infiltrating area w.r.t. matric head (m-1)
 real(dp)                      :: fInfRate                   ! infiltration rate at the surface (m s-1)
 real(dp)                      :: dInfRate                   ! derivative in the infiltration rate w.r.t matric head (s-1)
 real(dp)                      :: hCnd0,hCnd1                ! hydraulic conductivity in the surface layer (m s-1)
 real(dp)                      :: dPsi,dz                    ! spatial differenes in matric head and height
 real(dp),pointer              :: iceAbove,iceBelow          ! volumetric ice content in the layer above and the layer below
 real(dp),pointer              :: psiAbove,psiBelow          ! matric head in the layer above and the layer below
 real(dp),pointer              :: sHcAbove,sHcBelow          ! saturated hydraulic conductivity in the layer above and the layer below
 real(dp),pointer              :: zAbove,zBelow              ! height of the layer mid-point in the layer above and the layer below
 real(dp)                      :: flux0,flux1,flux2          ! used to test the numerical derivatives
 real(dp)                      :: infil                      ! infiltration at the surface (m s-1)
 integer(i4b)                  :: iLayer                     ! layer index
 ! initialize error control
 err=0; message="dq_dMatHead/"

 ! -----------------------------------------------------------------------------------------------------------------------------
 ! ***** compute the derivative in hydraulic conductivity at layer mid-points w.r.t. matric head (s-1)
 ! -----------------------------------------------------------------------------------------------------------------------------
 do iLayer=1,nLevels
  select case(dMethod)
  case(analytical)
   dHydCond_dMatric(iLayer) = dHydCond_dPsi(mLayerMatricHeadTrial(iLayer),mLayerSatHydCond(iLayer),vGn_alpha,vGn_n,vGn_m,iceImpedeFac(iLayer),.true.)  ! analytical
  case(numerical)
   dHydCond_dMatric(iLayer) = dHydCond_dPsi(mLayerMatricHeadTrial(iLayer),mLayerSatHydCond(iLayer),vGn_alpha,vGn_n,vGn_m,iceImpedeFac(iLayer),.false.) ! numerical
  case default
   err=10; message=trim(message)//"unknown option to compute derivative in hydraulic conductivity at layer mid-points w.r.t. matric head"; return
  end select
  !if(iLayer < 5) print*, 'dHydCond = ', iLayer, dMethod, dHydCond_dMatric(iLayer), mLayerHydCond(iLayer), iceImpedeFac(iLayer)
 end do

 ! compute the derivative in flux at layer interfaces w.r.t. matric head
 do iLayer=0,nLevels  ! loop through interfaces
  
  ! -----------------------------------------------------------------------------------------------------------------------------
  ! ***** the upper boundary
  ! -----------------------------------------------------------------------------------------------------------------------------
  if(iLayer==0)then  ! (upper boundary)

   dq_dMatricAbove(iLayer) = valueMissing  ! don't expect this to be used, so deliberately set to a ridiculous value to cause problems

   select case(bc_upper)

   ! *****
   ! coupled to groundwater
   case(groundwaterCouple)   ! coupled with groundwater
    err=20; message=trim(message)//"cannot couple mixed form of Richards' eqn with moving groundwater boundary"; return

   ! *****
   ! head boundary
   case(prescribedHead)      ! head boundary
    ! derivatives in the flux w.r.t. matric head
    if(dMethod==analytical)then
     dq_dMatricBelow(iLayer) = -iLayerHydCond(0)/(mLayerDepth(1)/2._dp)
    ! compute numerical derivatives
    else
     flux0 = -iLayerHydCond(0)*( mLayerMatricHeadTrial(1)     - upperBoundHead) / (mLayerDepth(1)*0.5_dp) + iLayerHydCond(0)
     flux1 = -iLayerHydCond(0)*((mLayerMatricHeadTrial(1)+dx) - upperBoundHead) / (mLayerDepth(1)*0.5_dp) + iLayerHydCond(0)
     dq_dMatricBelow(iLayer) = (flux1 - flux0)/dx
    endif

   ! *****
   ! flux boundary
   case(liquidFlux)          ! flux boundary 
    ! * compute derivatives in the flux w.r.t. matric head
    if(dMethod==analytical)then
     ! compute the function and derivative of the infiltrating area
     vFracLiq = volFracLiq(mLayerMatricHeadTrial(1),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)     ! fraction of liquid water 
     fracCap  = min((vFracLiq + mLayerVolFracIceTrial(1))/theta_sat, 1._dp - epsilon(dx))          ! fraction of capacity filled with liquid water and ice
     fInfArea = (1._dp - fracCap)**bpar_VIC                                                        ! area of the landscape where water infiltrates (-)
     dInfArea = (-mLayerdTheta_dPsi(1)/theta_sat) * bpar_VIC*(1._dp - fracCap)**(bpar_VIC - 1._dp) ! derivative in the infiltrating area w.r.t. matric head (m-1)
     ! compute the function and derivative of the infiltration rate
     fInfRate = min(iLayerSatHydCond(0),scalarRainPlusMelt)*iceImpedeFac(1)  ! (m s-1)
     dInfRate = 0._dp  ! (s-1)
     ! ...finally, compute analytical derivatives
     ! NOTE, only first term is used currently, but keep in this form to make life simpler once infiltration depends on matric head
     dq_dMatricBelow(iLayer) = dInfArea*fInfRate + fInfArea*dInfRate
    ! * compute numerical derivatives
    else
     infil = min(iLayerSatHydCond(0),scalarRainPlusMelt)*iceImpedeFac(1)  ! infiltration at the surface
     flux0 = (1._dp - satArea_matric(mLayerMatricHeadTrial(1),   mLayerVolFracIceTrial(1),theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m,bpar_VIC)) * infil 
     flux1 = (1._dp - satArea_matric(mLayerMatricHeadTrial(1)+dx,mLayerVolFracIceTrial(1),theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m,bpar_VIC)) * infil
     dq_dMatricBelow(iLayer) = (flux1 - flux0)/dx
    endif

   ! *****
   ! error check
   case default;
    err=20; message=trim(message)//'unknown upper boundary condition'; return
   endselect ; ! (lower boundary condition)

  ! -----------------------------------------------------------------------------------------------------------------------------
  ! ***** the lower boundary
  ! -----------------------------------------------------------------------------------------------------------------------------
  elseif(iLayer==nLevels)then  ! (lower boundary)

   dq_dMatricBelow(iLayer) = valueMissing ! don't expect this to be used, so deliberately set to a ridiculous value to cause problems

   select case(bc_lower)

   ! *****
   ! head boundary
   case(prescribedHead)      ! head boundary
    ! derivatives in the flux w.r.t. matric head
    if(dMethod==analytical)then
     dq_dMatricAbove(iLayer) = iLayerHydCond(nLevels)/(mLayerDepth(nLevels)/2._dp)
    ! compute numerical derivatives
    else
     flux0 = -iLayerHydCond(nLevels)*(lowerBoundHead -  mLayerMatricHeadTrial(nLevels)    ) / (mLayerDepth(nLevels)*0.5_dp) + iLayerHydCond(nLevels)
     flux1 = -iLayerHydCond(nLevels)*(lowerBoundHead - (mLayerMatricHeadTrial(nLevels)+dx)) / (mLayerDepth(nLevels)*0.5_dp) + iLayerHydCond(nLevels)
     dq_dMatricAbove(iLayer) = (flux1 - flux0)/dx
    endif

   ! *****
   ! free draining (analytical vs. numerical switch mage when computing dHydCond_dVolLiq)
   case(freeDrainage); dq_dMatricAbove(iLayer) = dHydCond_dMatric(iLayer)

   ! *****
   ! function of matric head in the bottom soil layer
   case(funcBottomHead)
    ! compute analytical derivatives
    if(dMethod==analytical)then
     dq_dMatricAbove(iLayer) = kAnisotropic*iLayerSatHydCond(0) * exp(mLayerMatricHeadTrial(nLevels)/zScale_TOPMODEL - mLayerHeight(nLevels)/zScale_TOPMODEL)/zScale_TOPMODEL
    ! compute numerical derivatives
    else
     flux0 = kAnisotropic*iLayerSatHydCond(0) * exp(-(mLayerHeight(nLevels) -  mLayerMatricHeadTrial(nLevels)    )/zScale_TOPMODEL)
     flux1 = kAnisotropic*iLayerSatHydCond(0) * exp(-(mLayerHeight(nLevels) - (mLayerMatricHeadTrial(nLevels)+dx))/zScale_TOPMODEL)
     dq_dMatricAbove(iLayer) = (flux1 - flux0)/dx
    endif

   ! *****
   ! error check
   case default;
    err=20; message=trim(message)//'unknown lower boundary condition'; return
   endselect ; ! (lower boundary condition)

  ! -----------------------------------------------------------------------------------------------------------------------------
  ! ***** internal layers
  ! -----------------------------------------------------------------------------------------------------------------------------
  else

   if(dMethod==analytical)then
    ! derivatives in hydraulic conductivity
    dHydCondIface_dMatricAbove = dHydCond_dMatric(iLayer)  *mLayerHydCond(iLayer+1) * 0.5_dp*(mLayerHydCond(iLayer)*mLayerHydCond(iLayer+1))**(-0.5_dp)
    dHydCondIface_dMatricBelow = dHydCond_dMatric(iLayer+1)*mLayerHydCond(iLayer)   * 0.5_dp*(mLayerHydCond(iLayer)*mLayerHydCond(iLayer+1))**(-0.5_dp)
    ! spatial differences in matric head and height
    dPsi  = mLayerMatricHeadTrial(iLayer+1) - mLayerMatricHeadTrial(iLayer)
    dz    = mLayerHeight(iLayer+1) - mLayerHeight(iLayer)
    ! derivatives in the flux w.r.t. matric head
    dq_dMatricAbove(iLayer) = -dHydCondIface_dMatricAbove*dPsi/dz + iLayerHydCond(iLayer)/dz + dHydCondIface_dMatricAbove
    dq_dMatricBelow(iLayer) = -dHydCondIface_dMatricBelow*dPsi/dz - iLayerHydCond(iLayer)/dz + dHydCondIface_dMatricBelow
   else
    ! get some short-cut variables (save typing)
    iceAbove => mLayerVolFracIceTrial(iLayer)
    iceBelow => mLayerVolFracIceTrial(iLayer+1)
    psiAbove => mLayerMatricHeadTrial(iLayer)
    psiBelow => mLayerMatricHeadTrial(iLayer+1)
    sHcAbove => mLayerSatHydCond(iLayer)
    sHcBelow => mLayerSatHydCond(iLayer+1)
    zAbove   => mLayerHeight(iLayer)
    zBelow   => mLayerHeight(iLayer+1)
    ! compute numerical derivatives
    flux0 = darcyFlux_matric(iceAbove,iceBelow,psiAbove,   psiBelow,   sHcAbove,sHcBelow,zAbove,zBelow,theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m,f_impede)
    flux1 = darcyFlux_matric(iceAbove,iceBelow,psiAbove+dx,psiBelow,   sHcAbove,sHcBelow,zAbove,zBelow,theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m,f_impede)
    flux2 = darcyFlux_matric(iceAbove,iceBelow,psiAbove,   psiBelow+dx,sHcAbove,sHcBelow,zAbove,zBelow,theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m,f_impede)
    dq_dMatricAbove(iLayer) = (flux1 - flux0)/dx
    dq_dMatricBelow(iLayer) = (flux2 - flux0)/dx
   endif  ! case for the numerical method

  endif  ! type of layer (upper, internal, or lower)

  ! print results
  !if(iLayer<5) write(*,'(2(i4,1x),2(e20.10,1x))') dMethod, iLayer, dq_dMatricAbove(iLayer), dq_dMatricBelow(iLayer)

 end do  ! looping through layers
 end subroutine dq_dMatHead



end module soilHydrol_module

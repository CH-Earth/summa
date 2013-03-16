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
 equilWaterTable,            & ! equilibrium water table
 pseudoWaterTable,           & ! pseudo water table
 bigBucket,                  & ! a big bucket (lumped aquifer model)
 noExplicit,                 & ! no explicit groundwater parameterization
 ! look-up values for the choice of boundary conditions for hydrology
 prescribedHead,             & ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
 funcBottomHead,             & ! function of matric head in the lower-most layer
 freeDrainage,               & ! free drainage
 liquidFlux,                 & ! liquid water flux
 zeroFlux                      ! zero flux
! -----------------------------------------------------------------------------------------------------------
implicit none
private
public::soilHydrol
public::eWaterTable
! number of layers
integer(i4b)           :: nSoil                     ! number of soil layers
integer(i4b)           :: nSnow                     ! number of snow layers
integer(i4b)           :: nLevels                   ! total number of soil layers to examine
! constant parameters
real(dp),parameter     :: valueMissing=-9999._dp    ! missing value parameter
real(dp),parameter     :: verySmall=epsilon(1.0_dp) ! a very small number (used to avoid divide by zero)
real(dp),parameter     :: dx=1.e-8_dp               ! finite difference increment
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
                        mpar_data%var(iLookPARAM%aquiferScaleFactor),                 & ! intent(in): scaling factor for aquifer storage in the big bucket (m)
                        mpar_data%var(iLookPARAM%bucketBaseflowExp),                  & ! intent(in): baseflow exponent for the big bucket (-)
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

                        ! factors limiting transpiration (from vegFlux routine)
                        mvar_data%var(iLookMVAR%mLayerRootDensity)%dat,                & ! intent(in): root density in each layer (-)
                        mvar_data%var(iLookMVAR%scalarAquiferRootFrac)%dat(1),         & ! intent(in): fraction of roots below the lowest soil layer (-)
                        mvar_data%var(iLookMVAR%scalarTranspireLim)%dat(1),            & ! intent(in): weighted average of the transpiration limiting factor (-)
                        mvar_data%var(iLookMVAR%mLayerTranspireLim)%dat,               & ! intent(in): transpiration limiting factor in each layer (-)
                        mvar_data%var(iLookMVAR%scalarTranspireLimAqfr)%dat(1),        & ! intent(in): transpiration limiting factor for the aquifer (-)

                        ! evaporation/transpiration fluxes (from vegFlux routine)
                        mvar_data%var(iLookMVAR%scalarCanopyTranspiration)%dat(1),    & ! intent(in): canopy transpiration (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarGroundEvaporation)%dat(1),      & ! intent(in): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarGroundSublimation)%dat(1),      & ! intent(in): ground sublimation/frost -- below canopy or non-vegetated (kg m-2 s-1)

                        ! initial fluxes (intent inout)
                        mvar_data%var(iLookMVAR%iLayerInitLiqFluxSoil)%dat(0:nLevels),& ! intent(inout): liquid flux at layer interfaces at the start of the time step (m s-1)
                        mvar_data%var(iLookMVAR%mLayerInitEjectWater)%dat(1:nLevels), & ! intent(inout): water ejected from each soil layer at the start-of-step (m s-1)
                        mvar_data%var(iLookMVAR%mLayerInitTranspire)%dat(1:nLevels),  & ! intent(inout): transpiration loss from each soil layer at the start of the time step (m s-1)
                        mvar_data%var(iLookMVAR%mLayerInitBaseflow)%dat(1:nLevels),   & ! intent(inout): baseflow from each soil layer at the start-of-step (m s-1)
                        mvar_data%var(iLookMVAR%scalarInitAquiferTranspire)%dat(1),   & ! intent(inout): transpiration loss from the aquifer at the start-of-step (m s-1)
                        mvar_data%var(iLookMVAR%scalarInitAquiferRecharge)%dat(1),    & ! intent(inout): (scalar) recharge to the aquifer at the start-of-step (m s-1)
                        mvar_data%var(iLookMVAR%scalarInitAquiferBaseflow)%dat(1),    & ! intent(inout): (scalar) baseflow from the aquifer at the start-of-step (m s-1)

                        ! diagnostic scalar variables
                        mvar_data%var(iLookMVAR%scalarRainPlusMelt)%dat(1),           & ! intent(out): (scalar) rain plus melt (m s-1)
                        mvar_data%var(iLookMVAR%scalarSurfaceRunoff)%dat(1),          & ! intent(out): (scalar) surface runoff (m s-1)
                        mvar_data%var(iLookMVAR%scalarWaterTableDepth)%dat(1),        & ! intent(out): (scalar) water table depth (m)
                        mvar_data%var(iLookMVAR%scalarAquiferTranspire)%dat(1),       & ! intent(out): transpiration loss from the aquifer (m s-1)
                        mvar_data%var(iLookMVAR%scalarAquiferRecharge)%dat(1),        & ! intent(out): (scalar) recharge to the aquifer at the end-of-step (m s-1)
                        mvar_data%var(iLookMVAR%scalarAquiferBaseflow)%dat(1),        & ! intent(out): (scalar) baseflow from the aquifer at the end-of-step (m s-1)

                        ! model diagnostic variables
                        mvar_data%var(iLookMVAR%mLayerdTheta_dPsi)%dat(1:nLevels),    & ! intent(out): derivative in the soil water characteristic w.r.t. psi (m-1)
                        mvar_data%var(iLookMVAR%mLayerdPsi_dTheta)%dat(1:nLevels),    & ! intent(out): derivative in the soil water characteristic w.r.t. theta (m)
                        mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat(0:nLevels),    & ! intent(out): liquid flux at layer interfaces at the end of the time step (m s-1)
                        mvar_data%var(iLookMVAR%mLayerEjectWater)%dat(1:nLevels),     & ! intent(out): water ejected from each soil layer (m s-1)
                        mvar_data%var(iLookMVAR%mLayerTranspire)%dat(1:nLevels),      & ! intent(out): transpiration loss from each soil layer (m s-1)
                        mvar_data%var(iLookMVAR%mLayerBaseflow)%dat(1:nLevels),       & ! intent(out): baseflow from each soil layer (m s-1)

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
                              dt,                          & ! intent(in): time step (seconds)
                              iter,                        & ! intent(in): current iteration count
                              mLayerMatricHeadIter,        & ! intent(in): matric head in each layer at the current iteration (m)
                              mLayerVolFracLiqIter,        & ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                              mLayerVolFracIceIter,        & ! intent(in): volumetric fraction of ice at the current iteration (-)
                              scalarAquiferStorageIter,    & ! intent(in): aquifer storage (m)

                              ! model decisions
                              ixDerivMethod,               & ! intent(in): choice of method used to compute derivative
                              ixRichards,                  & ! intent(in): choice of the form of Richards' equation
                              hc_Profile,                  & ! intent(in): index defining the option for the hydraulic conductivity profile
                              ixGroundwater,               & ! intent(in): choice of groundwater parameterization
                              ixBcUpperSoilHydrology,      & ! intent(in): choice of upper boundary condition for soil hydrology
                              ixBcLowerSoilHydrology,      & ! intent(in): choice of upper boundary condition for soil hydrology

                              ! model coordinate variables -- NOTE: use of ibeg and iend 
                              mLayerDepth,                 & ! intent(in): depth of the layer (m)
                              mLayerHeight,                & ! intent(in): height of the layer mid-point (m)
                              iLayerHeight,                & ! intent(in): height of the layer interfaces (m)

                              ! boundary conditions for matric head
                              lowerBoundHead,              & ! intent(in): lower boundary condition (m)
                              upperBoundHead,              & ! intent(in): upper boundary condition (m)

                              ! boundary conditions for volumetric liqquid water content
                              lowerBoundTheta,             & ! intent(in): lower boundary condition (-)
                              upperBoundTheta,             & ! intent(in): upper boundary condition (-)

                              ! model forcing
                              scalarRainfall,              & ! intent(in): computed rainfall rate (kg m-2 s-1)
                              scalarLiqFluxSnow,           & ! intent(in): liquid flux from the base of the snowpack (m s-1)

                              ! general model parameters
                              wimplicit,                   & ! intent(in): weight assigned to start-of-step fluxes (-)

                              ! soil parameters
                              vGn_alpha,                   & ! intent(in): van Genutchen "alpha" parameter (m-1)
                              vGn_n,                       & ! intent(in): van Genutchen "n" parameter (-)
                              VGn_m,                       & ! intent(in): van Genutchen "m" parameter (-)
                              theta_sat,                   & ! intent(in): soil porosity (-)
                              theta_res,                   & ! intent(in): soil residual volumetric water content (-)
                              kAnisotropic,                & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                              zScale_TOPMODEL,             & ! intent(in): scale factor for TOPMODEL-ish baseflow parameterization (m)
                              bpar_VIC,                    & ! intent(in): b-parameter in the VIC surface runoff parameterization (-)
                              specificYield,               & ! intent(in): specific yield (-)
                              specificStorage,             & ! intent(in): specific storage coefficient (m-1)
                              aquiferScaleFactor,          & ! intent(in): scaling factor for aquifer storage in the big bucket (m)
                              bucketBaseflowExp,           & ! intent(in): baseflow exponent for the big bucket (-)
                              f_impede,                    & ! intent(in): ice impedence factor (-)

                              ! model state variables -- NOTE: use of ibeg and iend
                              mLayerVolFracIce,            & ! intent(in): volumetric fraction of ice in each layer (-)
                              mLayerVolFracLiq,            & ! intent(in): volumetric fraction of liquid water in each layer (-)
                              mLayerMatricHead,            & ! intent(in): matric head in each layer (m)
                              scalarAquiferStorage,        & ! intent(in): (scalar)    ! aquifer storage at the start of the time step (m)
                              scalarSfcMeltPond,           & ! intent(in): (scalar)    ! ponded water caused by melt of the "snow without a layer" (kg m-2)

                              ! saturated hydraulic conductivity in each layer
                              mLayerSatHydCond,            & ! intent(in): saturated hydraulic conductivity at the mid-point of each layer (m s-1)
                              iLayerSatHydCond,            & ! intent(in): saturated hydraulic conductivity at the interface of each layer (m s-1)

                              ! factors limiting transpiration (from vegFlux routine)
                              mLayerRootDensity,           & ! intent(in): root density in each layer (-)
                              scalarAquiferRootFrac,       & ! intent(in): fraction of roots below the lowest soil layer (-)
                              scalarTranspireLim,          & ! intent(in): weighted average of the transpiration limiting factor (-)
                              mLayerTranspireLim,          & ! intent(in): transpiration limiting factor in each layer (-)
                              scalarTranspireLimAqfr,      & ! intent(in): transpiration limiting factor for the aquifer (-)

                              ! evaporation/transpiration fluxes (from vegFlux routine)
                              scalarCanopyTranspiration,   & ! intent(in): canopy transpiration (W m-2)
                              scalarGroundEvaporation,     & ! intent(in): ground evaporation/condensation -- below canopy or non-vegetated (W m-2)
                              scalarGroundSublimation,     & ! intent(in): ground sublimation/frost -- below canopy or non-vegetated (W m-2)

                              ! initial fluxes (intent inout)
                              iLayerInitLiqFluxSoil,       & ! intent(inout): liquid flux at layer interfaces at the start of the time step (m s-1)
                              mLayerInitEjectWater,        & ! intent(inout): water ejected from each soil layer at the start-of-step (m s-1)
                              mLayerInitTranspire,         & ! intent(inout): transpiration loss from each soil layer at the start of the time step (m s-1)
                              mLayerInitBaseflow,          & ! intent(inout): baseflow from each soil layer at the start-of-step (m s-1)
                              scalarInitAquiferTranspire,  & ! intent(inout): transpiration loss from the aquifer at the start-of-step (m s-1)
                              scalarInitAquiferRecharge,   & ! intent(inout): recharge to the aquifer at the start-of-step (m s-1)
                              scalarInitAquiferBaseflow,   & ! intent(inout): baseflow from the aquifer at the start-of-step (m s-1)

                              ! diagnostic scalar variables
                              scalarRainPlusMelt,          & ! intent(out): rain plus melt (m s-1)
                              scalarSurfaceRunoff,         & ! intent(out): surface runoff (m s-1)
                              scalarWaterTableDepth,       & ! intent(out): water table depth (m)
                              scalarAquiferTranspire,      & ! intent(out): transpiration loss from the aquifer (m s-1)
                              scalarAquiferRecharge,       & ! intent(out): recharge to the aquifer at the end-of-step (m s-1)
                              scalarAquiferBaseflow,       & ! intent(out): baseflow from the aquifer at the end-of-step (m s-1)

                              ! model diagnostic variables
                              mLayerdTheta_dPsi,           & ! intent(out): derivative in the soil water characteristic w.r.t. psi (m-1)
                              mLayerdPsi_dTheta,           & ! intent(out): derivative in the soil water characteristic w.r.t. theta (m)
                              iLayerLiqFluxSoil,           & ! intent(out): liquid flux at layer interfaces at the end of the time step (m s-1)
                              mLayerEjectWater,            & ! intent(out): water ejected from each soil layer (m s-1)
                              mLayerTranspire,             & ! intent(out): transpiration loss from each soil layer (m s-1)
                              mLayerBaseflow,              & ! intent(out): baseflow from each soil layer (m s-1)

                              ! output variables from the soilHydrol routine
                              mLayerMatricHeadNew,         & ! intent(out): matric head in each layer at the next iteration (m)
                              mLayerVolFracLiqNew,         & ! intent(out): volumetric fraction of liquid water at the next iteration (-)
                              scalarAquiferStorageNew,     & ! intent(out): aquifer storage (m)
                              err,message)                   ! intent(out): error control
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
 real(dp),intent(in)              :: dt                           ! time step (seconds)
 integer(i4b),intent(in)          :: iter                         ! iteration index
 real(dp),intent(in)              :: mLayerMatricHeadIter(:)      ! matric head in each layer at the current iteration (m)
 real(dp),intent(in)              :: mLayerVolFracLiqIter(:)      ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),intent(in)              :: mLayerVolFracIceIter(:)      ! volumetric fraction of ice at the current iteration (-)
 real(dp),intent(in)              :: scalarAquiferStorageIter     ! aquifer storage (m)
 ! model decisions
 integer(i4b),intent(in)          :: ixDerivMethod                ! choice of method used to compute derivative
 integer(i4b),intent(in)          :: ixRichards                   ! choice of the form of Richards' equation
 integer(i4b),intent(in)          :: hc_profile                   ! choice of type of hydraulic conductivity profile
 integer(i4b),intent(in)          :: ixGroundwater                ! choice of groundwater parameterization
 integer(i4b),intent(in)          :: ixBcUpperSoilHydrology       ! choice of upper boundary condition for soil hydrology
 integer(i4b),intent(in)          :: ixBcLowerSoilHydrology       ! choice of upper boundary condition for soil hydrology
 ! model coordinate variables
 real(dp),intent(in)              :: mLayerDepth(:)               ! depth of the layer (m)
 real(dp),intent(in)              :: mLayerHeight(:)              ! height of the layer mid-point (m)
 real(dp),intent(in)              :: iLayerHeight(0:)             ! height of the layer interfaces (m)
 ! diriclet boundary conditions for matric head
 real(dp),intent(in)              :: upperBoundHead               ! upper boundary condition for matric head (m)
 real(dp),intent(in)              :: lowerBoundHead               ! lower boundary condition for matric head (m)
 ! diriclet boundary conditions for volumetric liquid water content
 real(dp),intent(in)              :: upperBoundTheta              ! upper boundary condition for volumetric liquid water content (-)
 real(dp),intent(in)              :: lowerBoundTheta              ! lower boundary condition for volumetric liquid water content (-)
 ! model forcing
 real(dp),intent(in)              :: scalarRainfall               ! computed rainfall rate (kg m-2 s-1)
 real(dp),intent(in)              :: scalarLiqFluxSnow            ! liquid flux from the base of the snowpack (m s-1)
 ! general model parameters
 real(dp),intent(in)              :: wimplicit                    ! weight assigned to start-of-step fluxes (-)
 ! soil parameters
 real(dp),intent(in)              :: vGn_alpha                    ! van Genutchen "alpha" parameter (m-1)
 real(dp),intent(in)              :: vGn_n                        ! van Genutchen "n" parameter (-)
 real(dp),intent(in)              :: vGn_m                        ! van Genutchen "m" parameter (-)
 real(dp),intent(in)              :: theta_sat                    ! soil porosity (-)
 real(dp),intent(in)              :: theta_res                    ! soil residual volumetric water content (-)
 real(dp),intent(in)              :: kAnisotropic                 ! anisotropy factor for lateral hydraulic conductivity (-)
 real(dp),intent(in)              :: zScale_TOPMODEL              ! scale factor for TOPMODEL-ish baseflow parameterization (m)
 real(dp),intent(in)              :: bpar_VIC                     ! b-parameter in the VIC surface runoff parameterization (-)
 real(dp),intent(in)              :: specificYield                ! specific yield (-)
 real(dp),intent(in)              :: specificStorage              ! specific storage coefficient (m-1)
 real(dp),intent(in)              :: aquiferScaleFactor           ! scaling factor for aquifer storage in the big bucket (m)
 real(dp),intent(in)              :: bucketBaseflowExp            ! baseflow exponent for the big bucket (-)
 real(dp),intent(in)              :: f_impede                     ! ice impedence factor (-)
 ! state variables
 real(dp),intent(in)              :: mLayerVolFracIce(:)          ! volumetric fraction of ice at the start of the time step (-)
 real(dp),intent(in)              :: mLayerVolFracLiq(:)          ! volumetric fraction of liquid water at the start of the time step (-)
 real(dp),intent(in)              :: mLayerMatricHead(:)          ! matric head in each layer at the start of the time step (m)
 real(dp),intent(in)              :: scalarAquiferStorage         ! aquifer storage at the start of the time step (m)
 real(dp),intent(in)              :: scalarSfcMeltPond            ! ponded water caused by melt of the "snow without a layer" (kg m-2)
 ! saturated hydraulic conductivity
 real(dp),intent(in)              :: mLayerSatHydCond(:)          ! saturated hydraulic conductivity at the mid-point of each layer (m s-1)
 real(dp),intent(in)              :: iLayerSatHydCond(0:)         ! saturated hydraulic conductivity at the interface of each layer (m s-1)
 ! factors limiting transpiration (from vegFlux routine)
 real(dp),intent(in)              :: mLayerRootDensity(:)         ! root density in each layer (-)
 real(dp),intent(in)              :: scalarAquiferRootFrac        ! fraction of roots below the lowest soil layer (-)
 real(dp),intent(in)              :: scalarTranspireLim           ! weighted average of the transpiration limiting factor (-)
 real(dp),intent(in)              :: mLayerTranspireLim(:)        ! transpiration limiting factor in each layer (-)
 real(dp),intent(in)              :: scalarTranspireLimAqfr       ! transpiration limiting factor for the aquifer (-)
 ! evaporation/transpiration fluxes (from vegFlux routine)
 real(dp),intent(in)              :: scalarCanopyTranspiration    ! canopy transpiration (kg m-2 s-1)
 real(dp),intent(in)              :: scalarGroundEvaporation      ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
 real(dp),intent(in)              :: scalarGroundSublimation      ! ground sublimation/frost -- below canopy or non-vegetated (kg m-2 s-1)
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! ***** input/output variables -- start-of-step fluxes
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 real(dp),intent(inout)           :: iLayerInitLiqFluxSoil(0:)    ! liquid flux at soil layer interfaces at the start of the time step (m s-1)
 real(dp),intent(inout)           :: mLayerInitEjectWater(:)      ! water ejected from each soil layer at the start-of-step (m s-1)
 real(dp),intent(inout)           :: mLayerInitTranspire(:)       ! transpiration loss from each soil layer at the start of the time step (m s-1)
 real(dp),intent(inout)           :: mLayerInitBaseflow(:)        ! baseflow from each soil layer at the start-of-step (m s-1)
 real(dp),intent(inout)           :: scalarInitAquiferTranspire   ! transpiration loss from the aquifer at the start-of-step (m s-1)
 real(dp),intent(inout)           :: scalarInitAquiferRecharge    ! recharge to the aquifer at the start-of-step (m s-1)
 real(dp),intent(inout)           :: scalarInitAquiferBaseflow    ! baseflow from the aquifer at the start-of-step (m s-1)
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! ***** output variables
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! diagnostic scalar variables
 real(dp),intent(out)             :: scalarRainPlusMelt           ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
 real(dp),intent(out)             :: scalarSurfaceRunoff          ! surface runoff (m s-1)
 real(dp),intent(out)             :: scalarWaterTableDepth        ! water table depth (m)
 real(dp),intent(out)             :: scalarAquiferTranspire       ! transpiration loss from the aquifer at the start-of-step (m s-1)
 real(dp),intent(out)             :: scalarAquiferRecharge        ! recharge to the aquifer at the end-of-step (m s-1)
 real(dp),intent(out)             :: scalarAquiferBaseflow        ! baseflow from the aquifer at the end-of-step (m s-1)
 ! diagnostic variables for each layer
 real(dp),intent(out)             :: mLayerdTheta_dPsi(:)         ! derivative in the soil water characteristic w.r.t. psi (m-1)
 real(dp),intent(out)             :: mLayerdPsi_dTheta(:)         ! derivative in the soil water characteristic w.r.t. theta (m)
 real(dp),intent(out)             :: iLayerLiqFluxSoil(0:)        ! liquid flux at soil layer interfaces at the end of the time step (m s-1)
 real(dp),intent(out)             :: mLayerEjectWater(:)          ! water ejected from each soil layer (m s-1)
 real(dp),intent(out)             :: mLayerTranspire(:)           ! transpiration loss from each soil layer (m s-1)
 real(dp),intent(out)             :: mLayerBaseflow(:)            ! baseflow from each soil layer (m s-1)
 ! output variables from the soilHydrol routine
 real(dp),intent(out)             :: mLayerMatricHeadNew(:)       ! matric head in each layer at the next iteration (m)
 real(dp),intent(out)             :: mLayerVolFracLiqNew(:)       ! volumetric fraction of liquid water at the next iteration (m)
 real(dp),intent(out)             :: scalarAquiferStorageNew      ! aquifer storage (m)
 integer(i4b),intent(out)         :: err                          ! error code
 character(*),intent(out)         :: message                      ! error message
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! ***** local variables
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables (general)
 character(LEN=256)               :: cmessage                 ! error message of downwind routine
 logical(lgt)                     :: printflag                ! flag to print crap to the screen
 integer(i4b)                     :: iLayer                   ! layer index
 logical(lgt),parameter           :: calcJacobian=.false.     ! flag to compute the Jacobian matrix
 real(dp),parameter               :: maxDepthPond=1._dp       ! max depth of ponding allowed on the soil surface
 real(dp)                         :: availPorosity            ! pore space available to fill
 real(dp)                         :: availLiqWater            ! liquid water available to drain (negative) 
 real(dp)                         :: aquiferStorageTest       ! trial value for aquifer storage (m)
 real(dp),dimension(nLevels)      :: mLayerTranspireFrac      ! fraction of transpiration allocated to each soil layer (-)
 real(dp)                         :: aquiferTranspireFrac     ! fraction of transpiration allocated to the aquifer (-)
 ! check need to compute initial fluxes
 integer(i4b)                     :: nFlux                    ! index for flux calculation
 integer(i4b)                     :: ifluxInit                ! starting index for flux calculations (0 or 1)
 real(dp)                         :: maxdiffMatric(1)         ! used to check if we are starting on the first iteration
 ! flux derivatives
 real(dp),dimension(0:nLevels)    :: dq_dStateAbove           ! change in the flux in layer interfaces w.r.t. state variable in the layer above
 real(dp),dimension(0:nLevels)    :: dq_dStateBelow           ! change in the flux in layer interfaces w.r.t. state variable in the layer below
 real(dp),dimension(nLevels)      :: mLayerEjectWaterDeriv    ! derivative in water ejected from each soil layer (m s-1)
 real(dp)                         :: scalarAquiferRechargeDeriv ! derivative in the recharge flux w.r.t. aquifer storage (s-1)
 real(dp)                         :: scalarAquiferBaseflowDeriv ! derivative in the baseflow flux w.r.t. aquifer storage (s-1)
 ! residuals
 real(dp),dimension(nLevels)      :: mLayerVolFracLiqResidual ! residual in volumetric liquid water content (-)
 real(dp)                         :: scalarAquiferResidual    ! residual in aquifer storage (m) 
 ! tri-diagonal solution
 integer(i4b)                     :: nState                   ! number of state variables
 real(dp)                         :: wtim                     ! weighted time (s)
 real(dp),allocatable             :: d_m1(:)                  ! sub-diagonal elements of the tridiagonal system (m-1)
 real(dp),allocatable             :: diag(:)                  ! diagonal elements of the tridiagonal system (m-1)
 real(dp),allocatable             :: d_p1(:)                  ! super-diagonal elements of the tridiagonal system (m-1)
 real(dp),allocatable             :: rVec(:)                  ! right-hand-side vector (- or m)
 real(dp),allocatable             :: sInc(:)                  ! state increment (- or m)
 real(dp),dimension(nLevels)      :: mLayerMatricHeadDiff     ! iteration increment for matric head (m)
 real(dp),dimension(nLevels)      :: mLayerVolFracLiqDiff     ! iteration increment for volumetric fraction of liquid water (m)
 real(dp)                         :: scalarAquiferStorageDiff ! iteration increment for aquifer storage (m)
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

 ! define upper boundary fluxes (m s-1)
 if(ixBcUpperSoilHydrology==liquidFlux)then
  if(nSnow==0) scalarRainPlusMelt = scalarRainfall/iden_water + (scalarSfcMeltPond/dt)/iden_water  ! rainfall plus melt of the snow without a layer (convert to m s-1)
  if(nSnow>0)  scalarRainPlusMelt = scalarLiqFluxSnow                                              ! liquid water flux from the base of the snowpack (m s-1)
 endif
 if(ixBcUpperSoilHydrology==prescribedHead)then
  scalarRainPlusMelt = 0._dp
 endif

 ! compute the fraction of transpiration loss from each snow-soil layer, and the aquifer
 if(scalarTranspireLim > epsilon(scalarTranspireLim))then ! (transpiration may be non-zero even if the soil moisture limiting factor is zero)
  mLayerTranspireFrac(:) = mLayerRootDensity(:)*mLayerTranspireLim(:)/scalarTranspireLim
  aquiferTranspireFrac   = scalarAquiferRootFrac*scalarTranspireLimAqfr/scalarTranspireLim
 else ! (possible for there to be non-zero conductance and therefore transpiration in this case)
  mLayerTranspireFrac(:) = mLayerRootDensity(:)
  aquiferTranspireFrac   = scalarAquiferRootFrac
 endif
 ! check that the sums are OK
 if(abs(1._dp - sum(mLayerTranspireFrac)+aquiferTranspireFrac) > 1.e-8_dp)then
  message=trim(message)//'problem allocating transpiration flux to soil layers and the aquifer'
  err=20; return
 endif

 ! compute transpiration loss from each soil layer, and the aquifer (kg m-2 s-1 --> m s-1)
 mLayerTranspire        = mLayerTranspireFrac(:)*scalarCanopyTranspiration/iden_water
 scalarAquiferTranspire = aquiferTranspireFrac*scalarCanopyTranspiration/iden_water

 ! include ground evaporation in the top soil layer (if no snow)
 if(nSnow == 0)then
  mLayerTranspire(1)     = mLayerTranspire(1) + (scalarGroundEvaporation + scalarGroundSublimation)/iden_water
 endif

 ! initialize fluxes at the start of the time step
 if(iter == 1)then
  mLayerInitTranspire        = mLayerTranspire
  scalarInitAquiferTranspire = scalarAquiferTranspire
 endif

 ! get the number of state variables, and allocate space for the tri-diagonal matrix
 select case(ixGroundwater)
  case(pseudoWaterTable,bigBucket); nState = nLevels+1
  case(noExplicit,equilWaterTable); nState = nLevels
  case default; err=20; message=trim(message)//'unknown groundwater parameterization'; return
 end select
 allocate(d_m1(nState-1),diag(nState),d_p1(nState-1),rVec(nState),sInc(nState),stat=err)
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
                    .false.,                       & ! intent(in): flag to indicate if derivatives are desired 
                    ! input: trial state variables  (NOTE: use vectors from the start of the step)
                    mLayerMatricHead,              & ! intent(in): matric head (m)
                    mLayerVolFracLiq,              & ! intent(in): volumetric fraction of liquid water (-)
                    mLayerVolFracIce,              & ! intent(in): volumetric fraction of ice (-)
                    scalarAquiferStorage,          & ! intent(in): aquifer storage at the start of the step (m)
                    ! output: derivative in the soil water characteristic
                    mLayerdPsi_dTheta,             & ! intent(out): derivative in the soil water characteristic
                    mLayerdTheta_dPsi,             & ! intent(out): derivative in the soil water characteristic
                    ! output: diagnostic variables
                    scalarSurfaceRunoff,           & ! intent(out): surface runoff (m s-1)
                    scalarWaterTableDepth,         & ! intent(inout): water table depth (m) -- inout because use last value to intialize iterations for the equilbrium water table
                    ! output: fluxes for the extended state vector (soil layers plus aquifer)
                    iLayerInitLiqFluxSoil,         & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                    mLayerInitEjectWater,          & ! intent(out): water ejected because pore volume is close to capacity (m s-1)
                    mLayerInitBaseflow,            & ! intent(out): baseflow from each soil layer (m s-1)
                    ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                    dq_dStateAbove,                & ! intent(out): derivatives in the flux w.r.t. state variable in the layer above (m s-1 or s-1)
                    dq_dStateBelow,                & ! intent(out): derivatives in the flux w.r.t. state variable in the layer below (m s-1 or s-1)
                    mLayerEjectWaterDeriv,         & ! intent(out): derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
                    ! output: fluxes and derivatives for the aquifer
                    scalarInitAquiferRecharge,     & ! intent(out): recharge to the aquifer (m s-1)
                    scalarInitAquiferBaseflow,     & ! intent(out): total baseflow from the aquifer (m s-1)
                    scalarAquiferRechargeDeriv,    & ! intent(out): derivative in recharge to the aquifer w.r.t. aquifer storage (s-1)
                    scalarAquiferBaseflowDeriv,    & ! intent(out): derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
                    ! output: error control
                    err,cmessage)                    ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  
  ! compute fluxes
  ! NOTE: use trial state variables
  else
   call computeFlux(&
                    ! input: model control variables
                    .true.,                        & ! intent(in): flag to indicate if derivatives are desired 
                    ! input: trial state variables  (NOTE: use vectors from the start of the step)
                    mLayerMatricHeadIter,          & ! intent(in): matric head (m)
                    mLayerVolFracLiqIter,          & ! intent(in): volumetric fraction of liquid water (-)
                    mLayerVolFracIceIter,          & ! intent(in): volumetric fraction of ice (-)
                    scalarAquiferStorageIter,      & ! intent(in): aquifer storage at the start of the step (m)
                    ! output: derivative in the soil water characteristic
                    mLayerdPsi_dTheta,             & ! intent(out): derivative in the soil water characteristic
                    mLayerdTheta_dPsi,             & ! intent(out): derivative in the soil water characteristic
                    ! output: diagnostic variables
                    scalarSurfaceRunoff,           & ! intent(out): surface runoff (m s-1)
                    scalarWaterTableDepth,         & ! intent(inout): water table depth (m) -- inout because use last value to intialize iterations for the equilbrium water table
                    ! output: fluxes for the extended state vector (soil layers plus aquifer)
                    iLayerLiqFluxSoil,             & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                    mLayerEjectWater,              & ! intent(out): water ejected because pore volume is close to capacity (m s-1)
                    mLayerBaseflow,                & ! intent(out): baseflow from each soil layer (m s-1)
                    ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric liquid water -- in the layer above and layer below (m s-1 or s-1)
                    dq_dStateAbove,                & ! intent(out): derivatives in the flux w.r.t. state variable in the layer above (m s-1 or s-1)
                    dq_dStateBelow,                & ! intent(out): derivatives in the flux w.r.t. state variable in the layer below (m s-1 or s-1)
                    mLayerEjectWaterDeriv,         & ! intent(out): derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
                    ! output: fluxes and derivatives for the aquifer
                    scalarAquiferRecharge,         & ! intent(out): recharge to the aquifer (m s-1)
                    scalarAquiferBaseflow,         & ! intent(out): total baseflow from the aquifer (m s-1)
                    scalarAquiferRechargeDeriv,    & ! intent(out): derivative in recharge to the aquifer w.r.t. aquifer storage (s-1)
                    scalarAquiferBaseflowDeriv,    & ! intent(out): derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
                    ! output: error control
                    err,cmessage)                    ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! ***** assign initial fluxes, if necessary
   if(iter==1)then
    if(nFlux==ifluxInit)then   ! (handles special case where there is no separate computation for start-of-step fluxes)
     mLayerInitBaseflow        = mLayerBaseflow
     mLayerInitEjectWater      = mLayerEjectWater
     iLayerInitLiqFluxSoil     = iLayerLiqFluxSoil
     scalarInitAquiferBaseflow = scalarAquiferBaseflow
     scalarInitAquiferRecharge = scalarAquiferRecharge
    endif  ! (if computing initial fluxes)
   endif  ! (if the first iteration)

  endif ! (if computing standard fluxes)

 end do  ! looping through initial vectors

 !write(*,'(a,1x,i4,1x,e20.10)') 'ixDerivMethod, scalarAquiferRecharge = ', ixDerivMethod, scalarAquiferRecharge
 !write(*,'(a,1x,i4,1x,e20.10)') 'ixDerivMethod, scalarAquiferBaseflow = ', ixDerivMethod, scalarAquiferBaseflow
 !write(*,'(a,1x,i4,1x,e20.10)') 'ixDerivMethod, scalarAquiferRechargeDeriv = ', ixDerivMethod, scalarAquiferRechargeDeriv
 !write(*,'(a,1x,i4,1x,e20.10)') 'ixDerivMethod, scalarAquiferBaseflowDeriv = ', ixDerivMethod, scalarAquiferBaseflowDeriv
 !print*, 'specific Yield = ', specificYield

 ! *****
 ! compute the residual vector
 call liqResidual(&
                  ! input: control variables
                  dt,                         & ! intent(in): length of the time step (s)
                  wimplicit,                  & ! intent(in): weight assigned to the start-of-step
                  ! input: coordinate variables
                  mLayerDepth,                & ! intent(in): depth of each layer (m)
                  ! input: model parameters for the compressibility term
                  theta_sat,                  & ! intent(in): porosity (-)
                  specificStorage,            & ! intent(in): specific storage (m-1)
                  ! input: initial flux vectors (start of the time step)
                  iLayerInitLiqFluxSoil,      & ! intent(in): initial liquid flux at layer interfaces (m s-1)
                  mLayerInitEjectWater,       & ! intent(in): initial water ejected from each layer because volumetric liquid water is approaching porosity (m s-1)
                  mLayerInitBaseflow,         & ! intent(in): initial baseflow from each layer (m s-1)
                  mLayerInitTranspire,        & ! intent(in): initial transpiration from each layer -- from energy routine (m s-1)
                  ! input: trial flux vectors
                  iLayerLiqFluxSoil,          & ! intent(in): trial liquid flux at layer interfaces (m s-1)
                  mLayerEjectWater,           & ! intent(in): trial water ejected from each layer because volumetric liquid water is approaching porosity (m s-1)
                  mLayerBaseflow,             & ! intent(in): trial baseflow from each layer (m s-1)
                  mLayerTranspire,            & ! intent(in): trial transpiration from each layer -- from energy routine (m s-1)
                  ! input: initial state vectors (start of the time step)
                  mLayerMatricHead,           & ! intent(in): initial matric head (m)
                  mLayerVolFracLiq,           & ! intent(in): initial volumetric liquid water content (-)
                  mLayerVolFracIce,           & ! intent(in): initial volumetric ice content (-)
                  ! input: trial state vectors
                  mLayerMatricHeadIter,       & ! intent(in): trial matric head (m)
                  mLayerVolFracLiqIter,       & ! intent(in): trial volumetric liquid water content (-)
                  mLayerVolFracIceIter,       & ! intent(in): trial volumetric ice content (-)
                  ! output: residual vector (-)
                  mLayerVolFracLiqResidual,   & ! intent(out): residual vector for the volumetric fraction of liquid water (-)
                  ! output: error control
                  err,cmessage)                 ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 !print*, 'mLayerVolFracLiqResidual = ', mLayerVolFracLiqResidual

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
                  scalarInitAquiferRecharge,  & ! intent(in): recharge to the aquifer        (m s-1)
                  scalarInitAquiferBaseflow,  & ! intent(in): baseflow from the aquifer      (m s-1)
                  ! input: end-of-step fluxes
                  scalarAquiferTranspire,     & ! intent(in): transpiration from the aquifer (m s-1)
                  scalarAquiferRecharge,      & ! intent(in): recharge to the aquifer        (m s-1)
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
 d_m1(1:nLevels-1) = (wtim/mLayerDepth(2:nLevels  ))*(-dq_dStateAbove(1:nLevels-1) )
 d_p1(1:nLevels-1) = (wtim/mLayerDepth(1:nLevels-1))*( dq_dStateBelow(1:nLevels-1) )
 ! (get diagonal elements)
 select case(ixRichards)
  case(moisture); diag(1:nLevels) = (wtim/mLayerDepth(1:nLevels))*(-dq_dStateBelow(0:nLevels-1) + dq_dStateAbove(1:nLevels) + mLayerEjectWaterDeriv(1:nLevels) ) + 1._dp
  case(mixdform); diag(1:nLevels) = (wtim/mLayerDepth(1:nLevels))*(-dq_dStateBelow(0:nLevels-1) + dq_dStateAbove(1:nLevels) + mLayerEjectWaterDeriv(1:nLevels) ) + mLayerdTheta_dPsi(1:nLevels)
  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
 endselect
 ! (populate the volumetric liquid water content portion of the residual vector)
 rVec(1:nLevels) = mLayerVolFracLiqResidual

 ! *****
 ! populate the tri-diagnonal matrices for the aquifer
 if(nState == nLevels+1)then
  ! (compute the off-diagonal elements)
  d_p1(nLevels)   = (wtim/mLayerDepth(nLevels))*( dq_dStateBelow(nLevels) )  ! change in drainage flux w.r.t change in the aquifer storage (m-1)
  d_m1(nLevels)   =  wtim                      *(-dq_dStateAbove(nLevels) )  ! change in drainage flux w.r.t. change in soil moisture in the lowest unsaturated node (m)
  ! (compute the diagonal)
  diag(nLevels+1) = wtim * (-scalarAquiferRechargeDeriv + scalarAquiferBaseflowDeriv) + 1._dp  ! flux derivatives (s-1), then diagonal here is dimensionless
  ! (add aquifer storage to the residual vector)
  rVec(nLevels+1) = scalarAquiferResidual
 else
  ! check that there us no aquifer transpiration, recharge, or baseflow
  if(abs(scalarAquiferBaseflow) > tiny(dt) .or. scalarAquiferRecharge > tiny(dt) .or. abs(scalarAquiferTranspire) > tiny(dt))then
   err=20; message=trim(message)//'expect there to be no aquifer baseflow or aquifer transpiration when there is no explicit gw'; return
  endif
 endif  ! (if there is an aquifer) 

 ! *****
 ! test the Jacobian (if desired)
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


 ! *****
 ! solve the tri-diagonal system of equations
 call tridag(d_m1,diag,d_p1,-rVec,sInc,err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 ! *****
 ! extract the iteration increments for the soil layers
 select case(ixRichards)
  case(moisture); mLayerVolFracLiqDiff(1:nLevels) = sInc(1:nLevels)
  case(mixdform); mLayerMatricHeadDiff(1:nLevels) = sInc(1:nLevels)
  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
 endselect

 ! *****
 ! extract iteration increment for the aquifer
 if(nState==nLevels+1) scalarAquiferStorageDiff = sInc(nLevels+1)

 ! *****
 ! check that the soil moisture does not exceed constraints
 if(ixRichards == moisture)then
  do iLayer=1,nLevels
   availPorosity = theta_sat - mLayerVolFracLiqIter(iLayer) - mLayerVolFracIceIter(iLayer)*(iden_ice/iden_water)
   availLiqWater = theta_res - mLayerVolFracLiqIter(iLayer)
   if(mLayerVolFracLiqDiff(iLayer) > availPorosity) mLayerVolFracLiqDiff(iLayer) = availPorosity/2._dp
   if(mLayerVolFracLiqDiff(iLayer) < availLiqWater) mLayerVolFracLiqDiff(iLayer) = availLiqWater/2._dp
  end do
 endif

 ! *****
 ! update the fluxes and state variables for the soil layers
 do iLayer=0,nLevels
  select case(ixRichards)
   ! ***** moisture-based form of Richards' equation
   case(moisture)
    ! (update the fluxes)
    if(iLayer==0)then;           iLayerLiqFluxSoil(iLayer) = iLayerLiqFluxSoil(iLayer) + dq_dStateBelow(iLayer)*mLayerVolFracLiqDiff(iLayer+1)
    elseif(iLayer==nLevels)then; iLayerLiqFluxSoil(iLayer) = iLayerLiqFluxSoil(iLayer) + dq_dStateAbove(iLayer)*mLayerVolFracLiqDiff(iLayer)
    else;                        iLayerLiqFluxSoil(iLayer) = iLayerLiqFluxSoil(iLayer) + dq_dStateAbove(iLayer)*mLayerVolFracLiqDiff(iLayer) &
                                                                                       + dq_dStateBelow(iLayer)*mLayerVolFracLiqDiff(iLayer+1)
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
    if(iLayer==0)then;           iLayerLiqFluxSoil(iLayer) = iLayerLiqFluxSoil(iLayer) + dq_dStateBelow(iLayer)*mLayerMatricHeadDiff(iLayer+1)
    elseif(iLayer==nLevels)then; iLayerLiqFluxSoil(iLayer) = iLayerLiqFluxSoil(iLayer) + dq_dStateAbove(iLayer)*mLayerMatricHeadDiff(iLayer)
    else;                        iLayerLiqFluxSoil(iLayer) = iLayerLiqFluxSoil(iLayer) + dq_dStateAbove(iLayer)*mLayerMatricHeadDiff(iLayer) &
                                                                                       + dq_dStateBelow(iLayer)*mLayerMatricHeadDiff(iLayer+1)
    endif

    if(iLayer > 0) mLayerEjectWater(iLayer) = mLayerEjectWater(iLayer) + mLayerEjectWaterDeriv(iLayer)*mLayerMatricHeadDiff(iLayer)
    ! (update the state variables)
    if(iLayer > 0)then
     mLayerMatricHeadNew(iLayer) = mLayerMatricHeadIter(iLayer) + mLayerMatricHeadDiff(iLayer)
     mLayerVolFracLiqNew(iLayer) = volFracLiq(mLayerMatricHeadNew(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
     ! check for ponding on the soil surface
     if(mLayerMatricHeadNew(1) > maxDepthPond)then
      print*, 'mLayerMatricHeadNew(1:5) = ', mLayerMatricHeadNew(1:5)
      print*, 'iLayerLiqFluxSoil(0:1) = ', iLayerLiqFluxSoil(0:1)
      err=-20; message=trim(message)//'excessive ponding on the soil surface'; return
     endif
     !if(iLayer < 5) write(*,'(2(a,1x,3(e20.10,1x)))') 'VolFracLiq = ', mLayerVolFracLiqNew(iLayer), mLayerVolFracLiqIter(iLayer), mLayerVolFracLiqNew(iLayer) - mLayerVolFracLiqIter(iLayer), &
     !                                                 'matricHead = ', mLayerMatricHeadNew(iLayer), mLayerMatricHeadIter(iLayer), mLayerMatricHeadDiff(iLayer)
    endif
   ! ***** unknown form of Richards' equation
   case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
  endselect
 end do  ! (loop through layers)


 ! *****
 ! update the aquifer storage and fluxes
 if(nState==nLevels+1)then
  ! (update storage)
  scalarAquiferStorageNew = scalarAquiferStorageIter + scalarAquiferStorageDiff
  ! (update fluxes)
  scalarAquiferRecharge = scalarAquiferRecharge + scalarAquiferRechargeDeriv*scalarAquiferStorageDiff  ! recharge to the aquifer   (m s-1)
  scalarAquiferBaseflow = scalarAquiferBaseflow + scalarAquiferBaseflowDeriv*scalarAquiferStorageDiff  ! baseflow from the aquifer (m s-1)
 else
  scalarAquiferStorageNew = 0._dp
  scalarAquiferRecharge   = 0._dp
  scalarAquiferBaseflow   = 0._dp
 endif


 ! *****
 ! deallocate space for the tri-diagonal vectors
 deallocate(d_m1,diag,d_p1,rVec,sInc,stat=err)
 if(err/=0)then; err=20; message=trim(message)//'unable to deallocate space for the tri-diagonal vectors'; return; endif


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
                         deriv_desired,                & ! intent(in): flag indicating if derivatives are desired
                         ! input: trial state variables (NOTE: use "trial" instead of "iter" to avoid using the "iter" vectors defined in main subroutine)
                         mLayerMatricHeadTrial,        & ! intent(in): matric head (m)
                         mLayerVolFracLiqTrial,        & ! intent(in): volumetric fraction of liquid water (-)
                         mLayerVolFracIceTrial,        & ! intent(in): volumetric fraction of ice (-)
                         scalarAquiferStorageTrial,    & ! intent(in): aquifer storage at the start of the step (m)
                         ! output: derivative in the soil water characteristic
                         mLayerdPsi_dTheta,            & ! intent(out): derivative in the soil water characteristic
                         mLayerdTheta_dPsi,            & ! intent(out): derivative in the soil water characteristic
                         ! output: diagnostic variables
                         scalarSurfaceRunoff,          & ! intent(out): surface runoff (m s-1)
                         scalarWaterTableDepth,        & ! intent(inout): water table depth (m) -- inout because use last value to intialize iterations for the equilbrium water table
                         ! output: fluxes at layer interfaces
                         iLayerLiqFluxSoil,            & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                         mLayerEjectWater,             & ! intent(out): water ejected because pore volume is close to capacity (m s-1)
                         mLayerBaseflow,               & ! intent(out): baseflow from each soil layer (m s-1)
                         ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                         dq_dStateAbove,               & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer above (m s-1)
                         dq_dStateBelow,               & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer below (m s-1)
                         mLayerEjectWaterDeriv,        & ! intent(out): derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
                         ! output: fluxes and derivatives for the aquifer storage
                         scalarAquiferRecharge,        & ! intent(out): recharge to the aquifer (m s-1)
                         scalarAquiferBaseflow,        & ! intent(out): total baseflow from the aquifer (m s-1)
                         scalarAquiferRechargeDeriv,   & ! intent(out): derivative in recharge to the aquifer w.r.t. aquifer storage (s-1)
                         scalarAquiferBaseflowDeriv,   & ! intent(out): derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
                         ! output: error control
                         err,message)                    ! intent(out): error control
  USE soil_utils_module,only:hydCond_psi  ! compute hydraulic conductivity as a function of matric head
  USE soil_utils_module,only:hydCond_liq  ! compute hydraulic conductivity as a function of volumetric liquid water content
  implicit none
  ! input: model control
  logical(lgt),intent(in)          :: deriv_desired                 ! flag indicating if derivatives are desired
  ! trial model state variables
  real(dp),intent(in)              :: mLayerMatricHeadTrial(:)      ! matric head in each layer at the current iteration (m)
  real(dp),intent(in)              :: mLayerVolFracLiqTrial(:)      ! volumetric fraction of liquid water at the current iteration (-)
  real(dp),intent(in)              :: mLayerVolFracIceTrial(:)      ! volumetric fraction of ice at the current iteration (-)
  real(dp),intent(in)              :: scalarAquiferStorageTrial     ! aquifer storage at the current iteration (m)
  ! output: derivative in the soil water characteristic
  real(dp),intent(out)             :: mLayerdPsi_dTheta(:)          ! derivative in the soil water characteristic
  real(dp),intent(out)             :: mLayerdTheta_dPsi(:)          ! derivative in the soil water characteristic
  ! output: diagnostic variables
  real(dp),intent(out)             :: scalarSurfaceRunoff           ! surface runoff (m s-1)
  real(dp),intent(inout)           :: scalarWaterTableDepth         ! water table depth (m)  -- inout because use last value to intialize iterations for the equilbrium water table
  ! output: liquid fluxes at layer interfaces
  real(dp),intent(out)             :: iLayerLiqFluxSoil(0:)         ! liquid flux at soil layer interfaces (m s-1)
  real(dp),intent(out)             :: mLayerEjectWater(:)           ! water ejected because pore volume is close to capacity (m s-1)
  real(dp),intent(out)             :: mLayerBaseflow(:)             ! baseflow from each soil layer (m s-1)
  ! output: derivatives in fluxes w.r.t. state variables in the layer above and layer below (m s-1)
  real(dp),intent(out)             :: dq_dStateAbove(0:)            ! change in the flux in layer interfaces w.r.t. state variables in the layer above
  real(dp),intent(out)             :: dq_dStateBelow(0:)            ! change in the flux in layer interfaces w.r.t. state variables in the layer below
  real(dp),intent(out)             :: mLayerEjectWaterDeriv(:)      ! derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
  ! output: fluxes and derivatives for the aquifer storage
  real(dp),intent(out)             :: scalarAquiferRecharge         ! recharge to the aquifer (m s-1)
  real(dp),intent(out)             :: scalarAquiferBaseflow         ! total baseflow from the aquifer (m s-1)
  real(dp),intent(out)             :: scalarAquiferRechargeDeriv    ! derivative in recharge to the aquifer w.r.t. aquifer storage (s-1)
  real(dp),intent(out)             :: scalarAquiferBaseflowDeriv    ! derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
  ! output: error control
  integer(i4b),intent(out)         :: err                           ! error code
  character(*),intent(out)         :: message                       ! error message
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! local variables: general
  character(LEN=256)               :: cmessage                      ! error message of downwind routine
  logical(lgt)                     :: desireAnal                    ! flag to identify if analytical derivatives are desired
  logical(lgt)                     :: computeDrainage               ! flag to identify if there is a need to compute drainage from the bottom of the soil profile
  integer(i4b)                     :: nUnsat                        ! number of unsaturated soil layers
  integer(i4b)                     :: ixWaterTable                  ! index of layer that contains the water table
  integer(i4b)                     :: iSoil                         ! index of soil layer
  integer(i4b)                     :: itry                          ! index of different flux calculations
  integer(i4b)                     :: nFlux                         ! number of flux calculations required (>1 = numerical derivatives with one-sided finite differences)
  integer(i4b),parameter           :: unperturbed=0                 ! named variable to identify the case of unperturbed state variables
  integer(i4b),parameter           :: perturbState=1                ! named variable to identify the case where we perturb the state in the current layer
  integer(i4b),parameter           :: perturbStateAbove=2           ! named variable to identify the case where we perturb the state layer above
  integer(i4b),parameter           :: perturbStateBelow=3           ! named variable to identify the case where we perturb the state layer below
  integer(i4b)                     :: ixPerturb                     ! index of element in 2-element vector to perturb
  integer(i4b)                     :: ixOriginal                    ! index of perturbed element in the original vector
  real(dp)                         :: scalardPsi_dTheta             ! derivative in soil water characteristix, used for perturbations when computing numerical derivatives
  real(dp)                         :: totalWaterDeficit             ! water required to bring entire soil profile to saturation (m)
  real(dp)                         :: sumDepthWeightCond            ! summation of the depth-weighted hydraulic conductivity (m2 s-1)
  real(dp)                         :: totalBaseflow                 ! total baseflow from the soil profile (m s-1)
  real(dp)                         :: dq_dWaterTable                ! derivative in flux w.r.t. water table depth (s-1) 
  ! local variables: scalar trial values
  real(dp)                         :: scalarVolFracLiqTrial         ! trial value of volumetric liquid water content (-)
  real(dp)                         :: scalarMatricHeadTrial         ! trial value of matric head (m)
  real(dp)                         :: scalarWaterTableDepthTrial    ! trial value of depth to the water table (m)
  real(dp)                         :: scalarHydCondTrial            ! trial value of hydraulic conductivity (m s-1)
  ! local variables: vector trial values (2-element vectors)
  real(dp),dimension(2)            :: vectorVolFracLiqTrial         ! trial value of volumetric liquid water content (-)
  real(dp),dimension(2)            :: vectorMatricHeadTrial         ! trial value of matric head (m)
  real(dp),dimension(2)            :: vectorHydCondTrial            ! trial value of hydraulic conductivity (m s-1)
  real(dp),dimension(2)            :: vectorDiffuseTrial            ! trial value of hydraulic diffusivity (m2 s-1)
  ! local variables: temporary variables used to store fluxes (used to compute numerical derivatives)
  real(dp)                         :: scalarFlux                    ! vertical flux (m s-1)
  real(dp)                         :: scalarFlux_dState             ! vertical flux with perturbation to the current state (m s-1)
  real(dp)                         :: scalarFlux_dStateAbove        ! vertical flux with perturbation to the state above (m s-1)
  real(dp)                         :: scalarFlux_dStateBelow        ! vertical flux with perturbation to the state below (m s-1)
  ! local variables: transmittance
  real(dp),dimension(nLevels)      :: iceImpedeFac                  ! ice impedence factor at layer mid-points (-)
  real(dp),dimension(nLevels)      :: mLayerHydCond                 ! hydraulic conductivity at layer mid-point (m s-1)
  real(dp),dimension(nLevels)      :: mLayerDiffuse                 ! diffusivity at layer mid-point (m2 s-1)
  real(dp),dimension(0:nLevels)    :: iLayerHydCond                 ! hydraulic conductivity at layer interface (m s-1)
  real(dp),dimension(0:nLevels)    :: iLayerDiffuse                 ! diffusivity at layer interface (m2 s-1)
  real(dp),dimension(nLevels)      :: dHydCond_dVolLiq              ! derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
  real(dp),dimension(nLevels)      :: dDiffuse_dVolLiq              ! derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
  real(dp),dimension(nLevels)      :: dHydCond_dMatric              ! derivative in hydraulic conductivity w.r.t matric head (m s-1)
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

  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! compute the depth to the water table (m)
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  select case(ixGroundwater)
   ! equilibrium water table
   case(equilWaterTable)
    totalWaterDeficit = sum( (theta_sat - (mLayerVolFracLiqTrial(:) + mLayerVolFracIceTrial(:)) ) * mLayerDepth(:) )  ! water required to bring entire soil profile to saturation (m)
    call eWaterTable(&
                     totalWaterDeficit,          & ! input: total water required to bring entire soil profile to saturation (m)
                     scalarWaterTableDepth,      & ! input/output: trial value for water table depth (m)
                     err,cmessage)                 ! output: error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    !if(scalarWaterTableDepth >= iLayerHeight(nLevels))then
    ! err=20; message=trim(message)//'equilibrium water table cannot be below the depth of the soil profile [add more soil layers]'; return
    !endif
   ! pseudo water table
   case(pseudoWaterTable) 
    scalarWaterTableDepth = iLayerHeight(nLevels) - scalarAquiferStorageTrial/specificYield
    if(scalarAquiferStorageTrial < 0._dp)then
     err=20; message=trim(message)//'pseudo water table cannot be below the depth of the soil profile [add more soil layers]'; return
    endif
   ! no explicit water table
   case(bigBucket,noExplicit)
    scalarWaterTableDepth = valueMissing
   case default; err=20; message=trim(message)//'unknown groundwater option'; return
  end select

  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! identify the index of the layer that contains the water table
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  select case(ixGroundwater)
   ! explicit water table
   case(equilWaterTable,pseudoWaterTable)
    ! (compute the number of unsaturated layers)
    nUnSat = count(iLayerHeight(1:nSoil) <= scalarWaterTableDepth)  ! number of unsaturated layers
    !if(nUnSat == nSoil)then; err=20; message=trim(message)//'water table is not within the soil profile'; return; endif
    ! (identify soil layer that contains the water table)
    ixWaterTable = nUnsat + 1
   ! no explicit water table
   case(bigBucket,noExplicit)
    ixWaterTable = 0
   case default; err=20; message=trim(message)//'unknown groundwater option'; return
  end select


  ! *************************************************************************************************************************************************
  ! *************************************************************************************************************************************************

  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! compute diagnostic variables at the nodes throughout the soil profile
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
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
  !write(*,'(a,i4,1x,e20.10)') 'nUnsat, mLayerHydCond(nUnsat) = ', nUnsat, mLayerHydCond(nUnsat)


  ! *************************************************************************************************************************************************
  ! *************************************************************************************************************************************************

  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! compute infiltraton at the surface and its derivative w.r.t. mass in the upper soil layer
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! either one or multiple flux calls, depending on if using analytical or numerical derivatives
  do itry=nFlux,0,-1  ! (work backwards to ensure all computed fluxes come from the un-perturbed case)

   ! =====
   ! get input state variables...
   ! ============================
   ! identify the type of perturbation
   select case(itry)

    ! skip undesired perturbations
    case(perturbStateAbove); cycle  ! cannot perturb state above (does not exist) -- so keep cycling
    case(perturbState); cycle       ! perturbing the layer below the flux at the top interface

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
   !print*, 'scalarSurfaceRunoff, iLayerLiqFluxSoil(0) = ', scalarSurfaceRunoff, iLayerLiqFluxSoil(0)

   ! get copies of surface flux to compute numerical derivatives
   if(deriv_desired .and. ixDerivMethod==numerical)then
    select case(itry)
     case(unperturbed);       scalarFlux             = iLayerLiqFluxSoil(0)
     case(perturbStateBelow); scalarFlux_dStateBelow = iLayerLiqFluxSoil(0)
     case default; err=10; message=trim(message)//"unknown perturbation"; return
    end select
   endif

  end do  ! (looping through different flux calculations -- one or multiple calls depending if desire for numerical or analytical derivatives)

  ! set derivative w.r.t. state above to missing (does not exist)
  dq_dStateAbove(0) = valueMissing

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
   do itry=nFlux,0,-1  ! (work backwards to ensure all computed fluxes come from the un-perturbed case)


    ! =====
    ! determine layer to perturb
    ! ============================
    select case(itry)
     ! skip undesired perturbations
     case(perturbState); cycle       ! perturbing the layers above and below the flux at the interface
     ! identify the index for the perturbation
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

    ! compute total vertical flux, to compute derivatives
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
   !if(iLayer==1) write(*,'(a)') 'ixDerivMethod, iLayer, dq_dStateBelow(iLayer), dq_dStateAbove(iLayer), mLayerVolFracLiqTrial(iLayer+1), mLayerEjectWater(iLayer+1), mLayerEjectWaterDeriv(iLayer+1) = '
   !if(iLayer< 5) write(*,'(2(i4,1x),10(e15.7,1x))') ixDerivMethod, iLayer, dq_dStateBelow(iLayer), dq_dStateAbove(iLayer), mLayerVolFracLiqTrial(iLayer+1), mLayerEjectWater(iLayer+1), mLayerEjectWaterDeriv(iLayer+1)
   !if(iLayer==1) write(*,'(a)') 'ixDerivMethod, iLayer+1, mLayerVolFracLiqTrial(iLayer+1), mLayerVolFracIceTrial(iLayer+1), mLayerEjectWater(iLayer+1)'
   !if(iLayer< 5) write(*,'(2(i4,1x),10(e15.7,1x))') ixDerivMethod, iLayer+1, mLayerVolFracLiqTrial(iLayer+1), mLayerVolFracIceTrial(iLayer+1), mLayerEjectWater(iLayer+1)


  end do  ! (looping through soil layers)

  ! *************************************************************************************************************************************************
  ! *************************************************************************************************************************************************

  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! * compute water ejected from each layer when the total water content in the layer is approaching porosity
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! loop through all soil layers
  do iLayer=1,nLevels

   ! check if there is a need to eject water
   if(ixRichards==moisture .or. & ! (moisture-based form of RE)
     (ixRichards==mixdform .and. mLayerVolFracIceTrial(iLayer)>0._dp) .or. &  ! (mixed form of RE when ice is present)
      iLayer==1)then ! (top soil layer)

    ! either one or multiple flux calls, depending on if using analytical or numerical derivatives
    do itry=nFlux,0,-1  ! (work backwards to ensure all computed fluxes come from the un-perturbed case)
 
     ! =====
     ! get input state variables...
     ! ============================
     ! identify the type of perturbation
     select case(itry)
 
      ! do not perturb states above/below -- keep cycling
      case(perturbStateAbove); cycle
      case(perturbStateBelow); cycle
 
      ! un-perturbed case
      case(unperturbed)
       scalarVolFracLiqTrial = mLayerVolFracLiqTrial(iLayer)
       scalarMatricHeadTrial = mLayerMatricHeadTrial(iLayer)
 
      ! perturb soil state (one-sided finite differences)
      case(perturbState)
       ! (perturbation depends on the form of Richards' equation)
       select case(ixRichards)
        case(moisture)
         scalarVolFracLiqTrial = mLayerVolFracLiqTrial(iLayer) + dx
         scalarMatricHeadTrial = mLayerMatricHeadTrial(iLayer)
        case(mixdform)
         scalarVolFracLiqTrial = mLayerVolFracLiqTrial(iLayer)
         scalarMatricHeadTrial = mLayerMatricHeadTrial(iLayer) + dx
        case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
      end select ! (form of Richards' equation)
  
      ! check for an unknown perturbation 
      case default; err=10; message=trim(message)//"unknown perturbation"; return
  
     end select ! (type of perturbation)
 
     ! =====
     ! compute water ejected from the current layer...
     ! ===============================================
     call ejectWater(&
                     ! input: model control
                     desireAnal,                       & ! intent(in): flag indicating if derivatives are desired
                     ixRichards,                       & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                     ! input: state variables
                     scalarMatricHeadTrial,            & ! intent(in): matric head in each layer (m)
                     scalarVolFracLiqTrial,            & ! intent(in): volumetric liquid water content (-)
                     mLayerVolFracIceTrial(iLayer),    & ! intent(in): volumetric ice content in each layer (-)
                     ! input: derivative in the soil water characteristic w.r.t. psi
                     mLayerdTheta_dPsi(iLayer),        & ! intent(in): derivative in the soil water characteristic w.r.t. psi (m-1)
                     ! input: soil parameters
                     vGn_alpha,                        & ! intent(in): van Genutchen "alpha" parameter (m-1)
                     vGn_n,                            & ! intent(in): van Genutchen "n" parameter (-)
                     VGn_m,                            & ! intent(in): van Genutchen "m" parameter (-)
                     theta_sat,                        & ! intent(in): soil porosity (-)
                     theta_res,                        & ! intent(in): soil residual volumetric water content (-)
                     ! input: saturated hydraulic conductivity at the surface
                     iLayerSatHydCond(0),              & ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
                     ! output: ejected water flux and derivative
                     mLayerEjectWater(iLayer),         & ! intent(out): water ejected because pore volume is filled (m s-1)
                     mLayerEjectWaterDeriv(iLayer),    & ! intent(out): derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
                     ! output: error control
                     err,cmessage)                       ! intent(out): error control
     if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
     !print*, itry, iLayer, mLayerEjectWater(iLayer)
 
     ! get copies of surface flux to compute numerical derivatives
     if(deriv_desired .and. ixDerivMethod==numerical)then
      select case(itry)
       case(unperturbed);  scalarFlux        = mLayerEjectWater(iLayer)
       case(perturbState); scalarFlux_dState = mLayerEjectWater(iLayer)
       case default; err=10; message=trim(message)//"unknown perturbation"; return
      end select
     endif
 
    end do  ! (looping through different flux calculations -- one or multiple calls depending if desire for numerical or analytical derivatives)

    ! compute numerical derivatives
    if(deriv_desired .and. ixDerivMethod==numerical)then
     mLayerEjectWaterDeriv(iLayer) = (scalarFlux_dState - scalarFlux)/dx ! change in ejected flux w.r.t. change in the volumetric soil moisture or matric head in the current soil layer (m s-1, or s-1)
    endif

   ! not computing ejected water
   else  
    mLayerEjectWater(iLayer) = 0._dp         ! water ejected because pore volume is filled (m s-1)
    mLayerEjectWaterDeriv(iLayer) = 0._dp    ! derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
   endif

  end do  ! looping through soil layers

  ! *************************************************************************************************************************************************
  ! *************************************************************************************************************************************************

  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! * compute drainage flux from the bottom of the soil profile, and its derivative
  ! -------------------------------------------------------------------------------------------------------------------------------------------------

  ! =====
  ! check if there is a need to compute drainage...
  ! ===============================================
  select case(ixGroundWater)
   ! explicit water table
   case(equilWaterTable,pseudoWaterTable)
    if(scalarWaterTableDepth > iLayerHeight(nLevels))then
     computeDrainage = .true.
    else
     computeDrainage = .false.  ! -- no drainage from the bottom of the soil profile when water table is within the soil profile (baseflow instead)
    endif
   ! all other cases
   case default; computeDrainage = .true.
  end select

  ! NOTE: for the case of an explicit water table, only compute drainage flux if the water table is below the soil profile
  if(.not.computeDrainage)then  ! zero-flux lower boundary
   iLayerLiqFluxSoil(nLevels) = 0._dp       ! drainage flux (m s-1)
   dq_dStateAbove(nLevels)    = 0._dp       ! change in drainage flux w.r.t. change in state in lowest unsaturated node (m s-1 or s-1)
   dq_dStateBelow(nLevels)    = 0._dp       ! change in drainage flux w.r.t. change in the aquifer storage (s-1)

  ! ** now, proceed with the calculations
  else

   ! either one or multiple flux calls, depending on if using analytical or numerical derivatives
   do itry=nFlux,0,-1  ! (work backwards to ensure all computed fluxes come from the un-perturbed case)
 
    ! =====
    ! get input state variables...
    ! ============================
    ! identify the type of perturbation
    select case(itry)
 
     ! skip undesired perturbations
     case(perturbStateBelow); cycle   ! only perturb soil state at this time (perhaps perturb aquifer state later)
     case(perturbState); cycle        ! here pertubing the state above the flux at the interface
 
     ! un-perturbed case
     case(unperturbed)
      scalarVolFracLiqTrial      = mLayerVolFracLiqTrial(nLevels)
      scalarMatricHeadTrial      = mLayerMatricHeadTrial(nLevels)
 
     ! perturb soil state (one-sided finite differences)
     case(perturbStateAbove)
      select case(ixRichards)  ! (perturbation depends on the form of Richards' equation)
       case(moisture)
        scalarVolFracLiqTrial = mLayerVolFracLiqTrial(nLevels) + dx
        scalarMatricHeadTrial = mLayerMatricHeadTrial(nLevels)
       case(mixdform)
        scalarVolFracLiqTrial = mLayerVolFracLiqTrial(nLevels)
        scalarMatricHeadTrial = mLayerMatricHeadTrial(nLevels) + dx
       case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
      end select ! (form of Richards' equation)
 
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
    call qDrainFlux(&
                    ! input: model control
                    desireAnal,                      & ! intent(in): flag indicating if derivatives are desired
                    ixRichards,                      & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                    hc_profile,                      & ! intent(in): index defining the decrease of hydraulic conductivity with depth
                    ixBcLowerSoilHydrology,          & ! intent(in): index defining the type of boundary conditions
                    ! input: state variables
                    scalarMatricHeadTrial,           & ! intent(in): matric head in the lowest unsaturated node (m)
                    scalarVolFracLiqTrial,           & ! intent(in): volumetric liquid water content the lowest unsaturated node (-)
                    ! input: model coordinate variables
                    mLayerDepth(nLevels),            & ! intent(in): depth of the lowest unsaturated soil layer (m)
                    mLayerHeight(nLevels),           & ! intent(in): height of the lowest unsaturated soil node (m)
                    ! input: boundary conditions
                    lowerBoundHead,                  & ! intent(in): lower boundary condition (m)
                    lowerBoundTheta,                 & ! intent(in): lower boundary condition (-)
                    ! input: water table depth
                    scalarWaterTableDepth,           & ! intent(in): depth to the water table (m)
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
                    specificYield,                   & ! intent(in): specific yield (-)
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
      case(perturbStateBelow); err=20; message=trim(message)//'lower state should never be perturbed when computing drainage do not expect to get here'; return
      case default; err=10; message=trim(message)//"unknown perturbation"; return
     end select
    endif
 
   end do  ! (looping through different flux calculations -- one or multiple calls depending if desire for numerical or analytical derivatives)
 
   ! compute numerical derivatives
   ! NOTE: drainage derivatives w.r.t. state below are *actually* w.r.t. water table depth, so need to be corrected for aquifer storage
   !       (note also negative sign to account for inverse relationship between water table depth and aquifer storage)
   if(deriv_desired .and. ixDerivMethod==numerical)then
    dq_dStateAbove(nLevels) = (scalarFlux_dStateAbove - scalarFlux)/dx    ! change in drainage flux w.r.t. change in state in lowest unsaturated node (m s-1 or s-1)
    dq_dStateBelow(nLevels) = 0._dp  ! keep this here in case we want to couple some day....
   endif
 
   ! print progress
   !if(deriv_desired)then
   ! print*, 'ixDerivMethod, dq_dStateAbove(nLevels) = ', ixDerivMethod, dq_dStateAbove(nLevels)
   ! print*, 'ixDerivMethod, dq_dStateBelow(nLevels) = ', ixDerivMethod, dq_dStateBelow(nLevels)
   !endif

  endif  ! (if need to compute drainage flux from the bottom of the soil profile)

 
  ! *************************************************************************************************************************************************
  ! *************************************************************************************************************************************************
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! * compute the recharge to the water table and its derivative
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  select case(ixGroundwater)

   ! =====
   ! using the pseudo water table for the groundwater parameterization...
   ! ====================================================================
   case(pseudoWaterTable)
    scalarAquiferRecharge      = iLayerLiqFluxSoil(nUnsat)   ! recharge = drainage flux from the bottom of the soil profile (m s-1)
    scalarAquiferRechargeDeriv = 0._dp                       ! recharge does not depend on aquifer storage

!    ! either one or multiple flux calls, depending on if using analytical or numerical derivatives
!    do itry=nFlux,0,-1  ! (work backwards to ensure all computed fluxes come from the un-perturbed case)
! 
!     ! * perturb state variables
!     select case(itry)  ! (identify the type of perturbation)
!      ! (skip undesired perturbations)
!      case(perturbStateAbove); cycle  ! (only perturb recharge flux w.r.t. water table depth)
!      case(perturbState); cycle       ! (only perturb recharge flux w.r.t. water table depth)
!      ! (perturb away, baby)
!      case(unperturbed);       scalarWaterTableDepthTrial = scalarWaterTableDepth       ! un-perturbed case
!      case(perturbStateBelow); scalarWaterTableDepthTrial = scalarWaterTableDepth + dx  ! perturb aquifer (one-sided finite differences)
!      case default; err=10; message=trim(message)//"unknown perturbation"; return
!     end select ! (type of perturbation)
! 
!     ! * retrieve the matric head of the lowest unsaturated node
!     select case(ixRichards)
!      case(moisture); scalarMatricHeadTrial = matricHead(mLayerVolFracLiqTrial(nUnsat),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
!      case(mixdform); scalarMatricHeadTrial = mLayerMatricHeadTrial(nUnsat) 
!     end select
!     !print*, 'mLayerVolFracLiqTrial(nUnsat) = ', mLayerVolFracLiqTrial(nUnsat)
!     
!     ! * compute the recharge to the aquifer (m s-1) and its derivative w.r.t. water table depth (s-1)
!     call rechargeWT(&
!                     ! input: model control
!                     desireAnal,                & ! intent(in): flag indicating if derivatives are desired
!                     hc_profile,                & ! intent(in): index defining the decrease of hydraulic conductivity with depth
!                     ! input: state variables
!                     scalarMatricHeadTrial,     & ! intent(in): matric head in the lowest unsaturated node (m)
!                     scalarWaterTableDepthTrial,& ! intent(in): depth to the water table (m)
!                     ! input: diagnostic variables and parameters
!                     mLayerHeight(nUnsat),      & ! intent(in): height of the lowest unsaturated soil node (m)
!                     mLayerHydCond(nUnsat),     & ! intent(in): hydraulic conductivity at the node itself (m s-1)
!                     iLayerSatHydCond(0),       & ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
!                     zScale_TOPMODEL,           & ! intent(in): TOPMODEL scaling factor (m)
!                     ! output: recharge flux and its derivative w.r.t. aquifer storage
!                     scalarAquiferRecharge,     & ! intent(out): recharge flux (m s-1)
!                     dq_dWaterTable,            & ! intent(out): change in recharge flux w.r.t. change in the water table depth (s-1)
!                     ! output: error control
!                     err,cmessage)                ! intent(out): error control
!     if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
!     if(deriv_desired)then
!      !write(*,'(a,3(i4,1x),e20.10)') 'itry, unperturbed, perturbStateBelow, scalarAquiferRecharge = ',&
!      !                                itry, unperturbed, perturbStateBelow, scalarAquiferRecharge
!     endif 
!
!     ! * get copies of recharge flux to compute derivatives
!     if(deriv_desired .and. ixDerivMethod==numerical)then
!      select case(itry)
!       case(unperturbed);       scalarFlux             = scalarAquiferRecharge
!       case(perturbStateBelow); scalarFlux_dStateBelow = scalarAquiferRecharge
!       case(perturbStateAbove,perturbState); err=10; message=trim(message)//'only perturb water table depth when computing recharge flux -- should not get here'; return
!       case default; err=10; message=trim(message)//'unknown perturbation'; return
!      end select
!     endif
! 
!    end do  ! (looping through different flux calculations -- one or multiple calls depending if desire for numerical or analytical derivatives)
! 
!    ! compute derivatives
!    ! NOTE: recharge derivatives w.r.t. state below are *actually* w.r.t. water table depth, so need to be corrected for aquifer storage
!    if(deriv_desired)then
!     select case(ixDerivMethod)
!      case(numerical);  scalarAquiferRechargeDeriv = -(scalarFlux_dStateBelow - scalarFlux)/dx/specificYield ! change in drainage flux w.r.t. change in the aquifer storage (s-1)
!      case(analytical); scalarAquiferRechargeDeriv = -dq_dWaterTable/specificYield
!      case default; err=10; message=trim(message)//'unknown numerical method'; return
!     end select
!     !write(*,'(a,e20.10)') 'dq_dWaterTable = ', dq_dWaterTable
!     !write(*,'(a,e20.10)') 'scalarAquiferRechargeDeriv*specificYield = ', scalarAquiferRechargeDeriv*specificYield
!    else
!     scalarAquiferRechargeDeriv = valueMissing
!    endif

   ! =====
   ! using the big bucket...
   ! =======================
   case(bigBucket)
    scalarAquiferRecharge      = iLayerLiqFluxSoil(nLevels)  ! recharge = drainage flux from the bottom of the soil profile (m s-1)
    scalarAquiferRechargeDeriv = 0._dp                       ! recharge does not depend on aquifer storage

   ! =====
   ! no explicit aquifer...
   ! ======================
   case(equilWaterTable,noExplicit)
    scalarAquiferRecharge      = 0._dp 
    scalarAquiferRechargeDeriv = 0._dp

   ! =====
   ! error checking...
   ! =================
   case default; err=20; message=trim(message)//'unknown groundwater parameterization'; return

  end select  ; ! (choice of groundwater parameterization)




  ! *************************************************************************************************************************************************
  ! *************************************************************************************************************************************************
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! * compute the baseflow flux and its derivative
  ! -------------------------------------------------------------------------------------------------------------------------------------------------

  ! select groundwater options
  select case(ixGroundwater)

   ! ------------------------------------------------------------------------------------------------------------------------------------------------
   ! *****
   ! case of no explicit deep groundwater
   case(noExplicit)
    ! (baseflow flux from each soil layer)
    mLayerBaseflow = 0._dp
    ! (baseflow flux from the aquifer)
    scalarAquiferBaseflow      = 0._dp  ! baseflow from the aquifer (m s-1)
    scalarAquiferBaseflowDeriv = 0._dp  ! derivative in baseflow w.r.t. aquifer storage (s-1)

   ! ------------------------------------------------------------------------------------------------------------------------------------------------
   ! *****
   ! pseudo water table and equilibrium water table
   case(pseudoWaterTable,equilWaterTable)

    ! either one or multiple flux calls, depending on if using analytical or numerical derivatives
    do itry=nFlux,0,-1  ! (work backwards to ensure all computed fluxes come from the un-perturbed case)

     ! only perturb states for the pseudoWaterTable
     if(ixGroundwater == equilWaterTable)then
      if(itry /= unperturbed) cycle
     endif

     ! * identify the type of perturbation
     select case(itry)
      ! (skip undersired perturbations)
      case(perturbStateBelow); cycle   ! assume baseflow at "bottom" of aquifer, so perturb w.r.t. state above
      case(perturbState); cycle        ! assume baseflow at "bottom" of aquifer, so perturb w.r.t. state above
      ! (perturb)
      case(unperturbed);       scalarWaterTableDepthTrial = scalarWaterTableDepth
      case(perturbStateAbove); scalarWaterTableDepthTrial = scalarWaterTableDepth + dx
      case default; err=10; message=trim(message)//"unknown perturbation"; return
     end select ! (type of perturbation)

     ! only compute baseflow if water table within soil profile
     if(scalarWaterTableDepthTrial < iLayerHeight(nLevels))then
      ! * compute the total baseflow (m s-1)
      call QBtopmodel(&
                      deriv_desired,               & ! input: flag indicating if derivatives are desired
                      scalarWaterTableDepthTrial,  & ! input: water table depth (m)
                      iLayerSatHydCond(0),         & ! input: saturated hydraulic conductivity at the surface (m s-1)
                      kAnisotropic,                & ! input: anisotropy factor for lateral hydraulic conductivity (-)
                      zScale_TOPMODEL,             & ! input: scale factor for TOPMODEL-ish baseflow parameterization (m)
                      totalBaseflow,               & ! output: total baseflow (m s-1)
                      dq_dWaterTable,              & ! output: derivative in baseflow flux w.r.t. water table depth (m s-1)
                      err,cmessage)                  ! output: error control
      if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
     else
      totalBaseflow  = 0._dp
      dq_dWaterTable = 0._dp
     endif

     ! * get copies of baseflow flux to compute derivatives
     if(deriv_desired .and. ixDerivMethod==numerical)then
      select case(itry)
       case(unperturbed);       scalarFlux             = totalBaseflow
       case(perturbStateAbove); scalarFlux_dStateAbove = totalBaseflow
       case(perturbStateBelow,perturbState); err=10; message=trim(message)//'only perturb water table depth when computing baseflow flux -- should not get here'; return
       case default; err=10; message=trim(message)//'unknown perturbation'; return
      end select
     endif

    end do  ! (looping through different flux calculations -- one or multiple calls depending if desire for numerical or analytical derivatives)

    ! * compute derivatives
    ! NOTE: baseflow derivatives w.r.t. state below are *actually* w.r.t. water table depth, so need to be corrected for aquifer storage
    if(deriv_desired .and. ixGroundwater==pseudoWaterTable)then
     select case(ixDerivMethod)
      case(numerical);  scalarAquiferBaseflowDeriv = -(scalarFlux_dStateAbove - scalarFlux)/dx/specificYield ! change in baseflow flux w.r.t. change in the aquifer storage (s-1)
      case(analytical); scalarAquiferBaseflowDeriv = -dq_dWaterTable/specificYield
      case default; err=10; message=trim(message)//'unknown numerical method'; return
     end select
    else
     scalarAquiferBaseflowDeriv = valueMissing
    endif

    ! * compute the baseflow from each soil layer (m s-1)
    sumDepthWeightCond  = sum(mLayerHydCond(ixWaterTable:nSoil)*mLayerDepth(ixWaterTable:nSoil))
    mLayerBaseflow(1:nUnsat)           = 0._dp
    if(ixWaterTable <=nSoil)then
     mLayerBaseflow(ixWaterTable:nSoil) = totalBaseflow * mLayerHydCond(ixWaterTable:nSoil)*mLayerDepth(ixWaterTable:nSoil)/sumDepthWeightCond
    endif

    ! * save aquifer baseflow
    select case(ixGroundWater)
     case(pseudoWaterTable); scalarAquiferBaseflow = totalBaseflow
     case(equilWaterTable);  scalarAquiferBaseflow = 0._dp    ! do not keep track of the aquifer storage -- water table depth is diagnostic in this case
     case default; err=20; message=trim(message)//'unknown groundwater parameterization (expect only pseudoWaterTable and equilWaterTable at this point)'
    end select

   ! ------------------------------------------------------------------------------------------------------------------------------------------------
   ! ***** the big bucket
   case(bigBucket)

    ! either one or multiple flux calls, depending on if using analytical or numerical derivatives
    do itry=nFlux,0,-1  ! (work backwards to ensure all computed fluxes come from the un-perturbed case)

     ! * identify the type of perturbation
     select case(itry)
      ! (skip undersired perturbations)
      case(perturbStateBelow); cycle   ! assume baseflow at "bottom" of aquifer, so perturb w.r.t. state above
      case(perturbState); cycle        ! assume baseflow at "bottom" of aquifer, so perturb w.r.t. state above
      ! (perturb)
      case(unperturbed);       aquiferStorageTest = scalarAquiferStorageTrial
      case(perturbStateAbove); aquiferStorageTest = scalarAquiferStorageTrial + dx
      case default; err=10; message=trim(message)//"unknown perturbation"; return
     end select ! (type of perturbation)

     ! * compute the total baseflow (m s-1)
     call QBbigbuckt(&
                     deriv_desired,               & ! input: flag indicating if derivatives are desired
                     aquiferStorageTest,          & ! input: trial value of aquifer storage (m)
                     aquiferScaleFactor,          & ! input: scaling factor for aquifer storage in the big bucket (m)
                     iLayerSatHydCond(0),         & ! input: saturated hydraulic conductivity at the surface (m s-1)
                     kAnisotropic,                & ! input: anisotropy factor for lateral hydraulic conductivity (-)
                     bucketBaseflowExp,           & ! input: exponent in bucket baseflow parameterization
                     scalarAquiferBaseflow,       & ! output: total baseflow (m s-1)
                     scalarAquiferBaseflowDeriv,  & ! output: derivative in baseflow flux w.r.t. water table depth (m s-1)
                     err,cmessage)                  ! output: error control
     if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

     ! * get copies of baseflow flux to compute derivatives
     if(deriv_desired .and. ixDerivMethod==numerical)then
      select case(itry)
       case(unperturbed);       scalarFlux             = scalarAquiferBaseflow
       case(perturbStateAbove); scalarFlux_dStateAbove = scalarAquiferBaseflow
       case(perturbStateBelow,perturbState); err=10; message=trim(message)//'only perturb water table depth when computing baseflow flux -- should not get here'; return
       case default; err=10; message=trim(message)//'unknown perturbation'; return
      end select
     endif

    end do  ! (multiple flux calls for computing  numerical derivatives

    ! * compute derivatives
    ! NOTE: baseflow derivatives w.r.t. state below are *actually* w.r.t. water table depth, so need to be corrected for aquifer storage
    if(deriv_desired)then
     if(ixDerivMethod==numerical) scalarAquiferBaseflowDeriv = (scalarFlux_dStateAbove - scalarFlux)/dx
    else
     scalarAquiferBaseflowDeriv = valueMissing
    endif
    !print*, 'scalarAquiferBaseflowDeriv = ', scalarAquiferBaseflowDeriv
    !pause

    ! (baseflow flux from each soil layer)
    mLayerBaseflow = 0._dp

   ! ------------------------------------------------------------------------------------------------------------------------------------------------
   ! ***** error check
   case default
    err=20; message=trim(message)//'cannot identify option for groundwater parameterization'; return

  end select  ! (groundwater options)
  !print*, 'total baseflow = ', sum(mLayerBaseflow)


 
  ! *************************************************************************************************************************************************
  ! *************************************************************************************************************************************************

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
  ! local: fluxes used to calculate residual
  real(dp),dimension(nLevels)      :: mLayerTempEjectWater        ! water ejected because pore volume is close to capacity (m s-1)
  real(dp),dimension(nLevels)      :: mLayerTempBaseflow          ! baseflow flux from each soil layer (m s-1)
  real(dp),dimension(0:nLevels)    :: iLayerTempLiqFluxSoil       ! liquid fluxes at layer interfaces at the start of the time step (m s-1)
  real(dp)                         :: scalarTempAquiferRecharge   ! aquifer recharge flux (m s-1)
  real(dp)                         :: scalarTempAquiferBaseflow   ! aquifer baseflow flux (m s-1)
  real(dp)                         :: scalarTempAquiferRechargeDeriv ! derivative in aquifer recharge flux (m s-1)
  real(dp)                         :: scalarTempAquiferBaseflowDeriv ! derivative in aquifer baseflow flux (m s-1)
  ! Jacobian matrix (used for testing)
  integer(i4b),parameter           :: minLayer=1                  ! minimum layer to test/print
  integer(i4b),parameter           :: maxLayer=100                ! maximum layer to test/print
  real(dp)                         :: eps                         ! finite difference increment
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

  ! initialize the water table depth (m)
  local_scalarWaterTableDepth = scalarWaterTableDepth

  ! loop through desired layers
  do ijac=nState,1,-1

   ! perturb states
   if(iJac > nLevels)then  ! (can ONLY happen if aquifer storage)
    eps = -dx*specificYield
    scalarAquiferStorageTrial = scalarAquiferStorageInput + eps
   else
    eps = dx
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
                    local_scalarWaterTableDepth, & ! intent(inout): water table depth (m) -- inout because used to initialize iterations for equilibrium water table                
                    ! output: fluxes
                    iLayerTempLiqFluxSoil,       & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                    mLayerTempEjectWater,        & ! intent(out): water ejected because pore volume is close to capacity (m s-1)
                    mLayerTempBaseflow,          & ! intent(out): baseflow flux from each soil layer (m s-1)
                    ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                    local_dq_dStateAbove,        & ! intent(out): derivatives in the flux w.r.t. state variable in the layer above (m s-1 or s-1)
                    local_dq_dStateBelow,        & ! intent(out): derivatives in the flux w.r.t. state variable in the layer below (m s-1 or s-1)
                    local_mLayerEjectWaterDeriv, & ! intent(out): derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
                    ! output: fluxes and derivatives for the aquifer storage
                    scalarTempAquiferRecharge,     & ! intent(out): recharge to the aquifer (m s-1)
                    scalarTempAquiferBaseflow,     & ! intent(out): total baseflow from the aquifer (m s-1)
                    scalarTempAquiferRechargeDeriv,& ! intent(out): derivative in recharge to the aquifer w.r.t. aquifer storage (s-1)
                    scalarTempAquiferBaseflowDeriv,& ! intent(out): derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
                    ! output: error control
                    err,cmessage)               ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   !if(iJac == nState)then
   ! write(*,'(a,e20.10)') 'eps = ', eps
   ! write(*,'(a,e20.10)') 'scalarWaterTableDepth = ', scalarWaterTableDepth
   ! write(*,'(a,e20.10)') 'local_scalarWaterTableDepth = ', local_scalarWaterTableDepth
   ! write(*,'(a,e20.10)') 'dWT = ', local_scalarWaterTableDepth - scalarWaterTableDepth
   ! write(*,'(a,e20.10)') 'scalarTempAquiferRecharge = ', scalarTempAquiferRecharge
   ! write(*,'(a,e20.10)') 'scalarTempAquiferBaseflow = ', scalarTempAquiferBaseflow
   ! write(*,'(a,e20.10)') 'scalarAquiferRecharge = ', scalarAquiferRecharge
   ! write(*,'(a,e20.10)') 'scalarAquiferBaseflow = ', scalarAquiferBaseflow
   ! write(*,'(a,e20.10)') 'recharge flux derivative = ', (scalarTempAquiferRecharge - scalarAquiferRecharge)/eps
   ! pause
   !endif

   ! *****
   ! compute the residual vector
   call liqResidual(&
                    ! control variables
                    dt,                        & ! intent(in): length of the time step (s)
                    wimplicit,                 & ! intent(in):weight assigned to the start-of-step
                    ! coordinate variables
                    mLayerDepth,               & ! intent(in): depth of each layer (m)
                    ! input: model parameters for the compressibility term
                    theta_sat,                 & ! intent(in): porosity (-)
                    specificStorage,           & ! intent(in): specific storage (m-1)
                    ! initial flux vectors (start of the time step)
                    iLayerInitLiqFluxSoil,     & ! intent(in): initial liquid flux at layer interfaces (m s-1)
                    mLayerInitEjectWater,      & ! intent(in): initial water ejected from each layer because volumetric liquid water is approaching porosity (m s-1)
                    mLayerInitBaseflow,        & ! intent(in): initial baseflow from each layer (m s-1)
                    mLayerInitTranspire,       & ! intent(in): initial transpiration from each layer -- from energy routine (m s-1)
                    ! trial flux vectors
                    iLayerTempLiqFluxSoil,     & ! intent(in): trial liquid flux at layer interfaces (m s-1)
                    mLayerTempEjectWater,      & ! intent(in): trial water ejected from each layer because volumetric liquid water is approaching porosity (m s-1)
                    mLayerTempBaseflow,        & ! intent(in): trial baseflow from each layer (m s-1)
                    mLayerTranspire,           & ! intent(in): trial transpiration from each layer -- from energy routine (m s-1)
                    ! initial state vectors (start of the time step)
                    mLayerMatricHead,          & ! intent(in): initial matric head (m)
                    mLayerVolFracLiq,          & ! intent(in): initial volumetric liquid water content (-)
                    mLayerVolFracIce,          & ! intent(in): initial volumetric ice content (-)
                    ! trial state vectors
                    mLayerMatricHeadTrial,     & ! intent(in): trial matric head (m)
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
                    scalarInitAquiferRecharge,  & ! intent(in): recharge to the aquifer        (m s-1)
                    scalarInitAquiferBaseflow,  & ! intent(in): baseflow from the aquifer      (m s-1)
                    ! input: end-of-step fluxes
                    scalarAquiferTranspire,     & ! intent(in): transpiration from the aquifer (m s-1)
                    scalarTempAquiferRecharge,  & ! intent(in): recharge to the aquifer        (m s-1)
                    scalarTempAquiferBaseflow,  & ! intent(in): baseflow from the aquifer      (m s-1)
                    ! output: aquifer residual
                    fAquifer,                   & ! intent(out): aquifer residual (m)
                    ! output: error control
                    err,cmessage)                 ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   !if(ijac == nState)then
   ! write(*,'(a,2(e20.10,1x))') 'scalarAquiferResidual, faquifer, dg/dS = ', scalarAquiferResidual, fAquifer, (-(fAquifer - scalarAquiferResidual)/eps - 1._dp)/dt
   ! pause
   !endif

   ! compute Jacobian
   if(ixGroundwater==bigBucket .or. ixGroundwater==pseudoWaterTable)then
    jmat(:,ijac) = ( (/ftest(:),fAquifer/) - (/mLayerVolFracLiqResidual(:),scalarAquiferResidual/) ) / eps
   else
    jmat(:,ijac) = (ftest(:) - mLayerVolFracLiqResidual(:) ) / eps
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
  do iJac=minLayer,min(maxLayer,nState)
   ! print the Jacobian
   if(iJac==1)then;          write(*,'(2(i4,1x),2(a,1x,3(f30.10,1x)))') ixDerivMethod, iJac, 'test hyd Jacobian', (/valueMissing,      jmat(iJac,1:iJac+1)/), '--> tri-diag = ', valueMissing, diag(iJac), d_p1(iJac)
   elseif(iJac==nState)then; write(*,'(2(i4,1x),2(a,1x,3(f30.10,1x)))') ixDerivMethod, iJac, 'test hyd Jacobian', (/jmat(iJac,iJac-1:nState), valueMissing/), '--> tri-diag = ', d_m1(iJac-1), diag(iJac), valueMissing
   else;                     write(*,'(2(i4,1x),2(a,1x,3(f30.10,1x)))') ixDerivMethod, iJac, 'test hyd Jacobian', (/jmat(iJac,iJac-1:iJac+1)              /), '--> tri-diag = ', d_m1(iJac-1), diag(iJac), d_p1(iJac) 
   endif
  end do
  !pause 'testing Jacobian' 

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
 ! public subroutine: compute the equilibrium water table depth
 ! ---------------------------------------------------------------------------------------------------------------------------
 subroutine eWaterTable(&
                        totalWaterDeficit,       & ! intent(in): total water required to bring entire soil profile to saturation (m)
                        zWaterTrial,             & ! intent(in): trial value for water table depth (m)
                        err,message)               ! intent(out): error control
 USE soil_utils_module,only:satDeficit  ! function to integrate
 USE integr8func_module,only:qromb      ! Romberg integration
 implicit none
 ! dummy variables
 real(dp),intent(in)        :: totalWaterDeficit   ! total water required to bring entire soil profile to saturation (m)
 real(dp),intent(inout)     :: zWaterTrial         ! trial value for water table depth (m)
 integer(i4b),intent(out)   :: err                 ! error code
 character(*),intent(out)   :: message             ! error message
 ! local variables
 integer(i4b)               :: iter                ! iteration index
 integer(i4b),parameter     :: maxiter=20          ! maximum number of iterations
 real(dp),parameter         :: finc=1.e-9_dp       ! finite-difference increment (m)
 real(dp),parameter         :: Xtol=1.e-8_dp       ! convergence tolerance (m)
 real(dp)                   :: f0,f1               ! function evaluations in iterative solution (m)
 real(dp)                   :: dfdx                ! derivative in function evaluations (-)
 real(dp)                   :: zInc                ! iteration increment (m)
 ! ---------------------------------------------------------------------------------------------------------------------------
 err=0; message='eWaterTable/' 
 ! iterate until convergence
 do iter=1,maxiter
  ! compute function evaluations (residuals) -- integrate sat deficit from water table depth to the surface
  f0 = (-qromb(satDeficit,0._dp,-zWaterTrial) - totalWaterDeficit)
  f1 = (-qromb(satDeficit,0._dp,-zWaterTrial-finc) - totalWaterDeficit)
  ! compute derivative and iteration increment
  dfdx = (f0 - f1)/finc
  zInc = f0/dfdx
  ! compute new value
  zWaterTrial = zWaterTrial + zInc
  if(zWaterTrial < Xtol) zWaterTrial = Xtol
  !write(*,'(a,10(f16.12,1x))') 'zWaterTrial, f0, f1, zInc = ', zWaterTrial, f0, f1, zInc
  ! check for convergence
  if(abs(f0) < Xtol) return
  if(iter==maxiter)then; err=20; message=trim(message)//'failed to converge'; return; endif
 enddo  ! (iterating)
 end subroutine eWaterTable


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
                       err,message)                 ! intent(out): error control
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables (named variables to provide index of 2-element vectors)
 integer(i4b),parameter        :: ixUpper=1                   ! index of upper node in the 2-element vectors
 integer(i4b),parameter        :: ixLower=2                   ! index of lower node in the 2-element vectors
 ! local variables (Darcy flux)
 real(dp)                      :: dPsi                        ! spatial difference in matric head (m)
 real(dp)                      :: dLiq                        ! spatial difference in volumetric liquid water (-)
 real(dp)                      :: dz                          ! spatial difference in layer mid-points (m)
 real(dp)                      :: cflux                       ! capillary flux (m s-1)
 ! local variables (derivative in Darcy's flux)
 real(dp)                      :: dHydCondIface_dVolLiqAbove  ! derivative in hydraulic conductivity at layer interface w.r.t. volumetric liquid water content in layer above
 real(dp)                      :: dHydCondIface_dVolLiqBelow  ! derivative in hydraulic conductivity at layer interface w.r.t. volumetric liquid water content in layer below
 real(dp)                      :: dDiffuseIface_dVolLiqAbove  ! derivative in hydraulic diffusivity at layer interface w.r.t. volumetric liquid water content in layer above
 real(dp)                      :: dDiffuseIface_dVolLiqBelow  ! derivative in hydraulic diffusivity at layer interface w.r.t. volumetric liquid water content in layer below
 real(dp)                      :: dHydCondIface_dMatricAbove  ! derivative in hydraulic conductivity at layer interface w.r.t. matric head in layer above
 real(dp)                      :: dHydCondIface_dMatricBelow  ! derivative in hydraulic conductivity at layer interface w.r.t. matric head in layer below
 ! ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="iLayerFlux/"

 ! *****
 ! compute the vertical flux of liquid water
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
  case default; err=10; message=trim(message)//"unable to identify option for Richards' equation"; return
 end select
 ! compute the total flux (add gravity flux, positive downwards)
 iLayerLiqFluxSoil = cflux + iLayerHydCond

 ! ** compute the derivatives
 if(deriv_desired)then
  select case(ixRichards)  ! (form of Richards' equation)
   case(moisture)
    ! derivatives in hydraulic conductivity at the layer interface (m s-1)
    dHydCondIface_dVolLiqAbove = dHydCond_dVolLiq(ixUpper)*nodeHydCondTrial(ixLower) * 0.5_dp/max(iLayerHydCond,verySmall)
    dHydCondIface_dVolLiqBelow = dHydCond_dVolLiq(ixLower)*nodeHydCondTrial(ixUpper) * 0.5_dp/max(iLayerHydCond,verySmall)
    ! derivatives in hydraulic diffusivity at the layer interface (m2 s-1)
    dDiffuseIface_dVolLiqAbove = dDiffuse_dVolLiq(ixUpper)*nodeDiffuseTrial(ixLower) * 0.5_dp/max(iLayerDiffuse,verySmall)
    dDiffuseIface_dVolLiqBelow = dDiffuse_dVolLiq(ixLower)*nodeDiffuseTrial(ixUpper) * 0.5_dp/max(iLayerDiffuse,verySmall)
    ! derivatives in the flux w.r.t. volumetric liquid water content
    dq_dStateAbove = -dDiffuseIface_dVolLiqAbove*dLiq/dz + iLayerDiffuse/dz + dHydCondIface_dVolLiqAbove
    dq_dStateBelow = -dDiffuseIface_dVolLiqBelow*dLiq/dz - iLayerDiffuse/dz + dHydCondIface_dVolLiqBelow
   case(mixdform)
    ! derivatives in hydraulic conductivity
    dHydCondIface_dMatricAbove = dHydCond_dMatric(ixUpper)*nodeHydCondTrial(ixLower) * 0.5_dp/max(iLayerHydCond,verySmall)
    dHydCondIface_dMatricBelow = dHydCond_dMatric(ixLower)*nodeHydCondTrial(ixUpper) * 0.5_dp/max(iLayerHydCond,verySmall)
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
 ! private subroutine: compute special case recharge to the water table, and its derivatives
 ! ************************************************************************************************
 subroutine rechargeWT(&
                       ! input: model control
                       deriv_desired,             & ! intent(in): flag indicating if derivatives are desired
                       hc_profile,                & ! intent(in): index defining the decrease of hydraulic conductivity with depth
                       ! input: state variables
                       nodeMatricHead,            & ! intent(in): matric head in the lowest unsaturated node (m)
                       scalarWaterTableDepth,     & ! intent(in): depth to the water table (m)
                       ! input: diagnostic variables and parameters
                       nodeHeight,                & ! intent(in): height of the lowest unsaturated soil node (m)
                       nodeHydCond,               & ! intent(in): hydraulic conductivity at the node itself (m s-1)
                       surfaceSatHydCond,         & ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
                       zScale_TOPMODEL,           & ! intent(in): TOPMODEL scaling factor (m)
                       ! output: recharge flux and its derivative w.r.t. aquifer storage
                       scalarRecharge,            & ! intent(out): recharge flux (m s-1)
                       dq_dWaterTable,            & ! intent(out): change in recharge flux w.r.t. change in the water table depth (s-1)
                       ! output: error control
                       err,message)                 ! intent(out): error control
 ! compute recharge to the water table, and its derivative w.r.t. volumetric liquid water content in the layer above and aquifer storage
 implicit none
 ! -----------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 logical(lgt),intent(in)       :: deriv_desired             ! flag to indicate if derivatives are desired
 integer(i4b),intent(in)       :: hc_profile                ! index defining the decrease of hydraulic conductivity with depth
 ! input: state and diagnostic variables
 real(dp),intent(in)           :: nodeMatricHead            ! matric head in the lowest unsaturated node (m)
 real(dp),intent(in)           :: scalarWaterTableDepth     ! depth to the water table (m)
 ! input: model coordinate variables
 real(dp),intent(in)           :: nodeHeight                ! height of the lowest unsaturated soil node (m)
 ! input: diagnostic variables and parameters
 real(dp),intent(in)           :: nodeHydCond               ! hydraulic conductivity at the node itself (m s-1)
 real(dp),intent(in)           :: surfaceSatHydCond         ! saturated hydraulic conductivity at the surface (m s-1)
 real(dp),intent(in)           :: zScale_TOPMODEL           ! scale factor for TOPMODEL-ish baseflow parameterization (m)
 ! -----------------------------------------------------------------------------------------------------------------------------
 ! output: recharge flux and its derivative w.r.t. aquifer storage
 real(dp),intent(out)          :: scalarRecharge            ! recharge flux (m s-1)
 real(dp),intent(out)          :: dq_dWaterTable            ! change in recharge flux w.r.t. change in the water table depth (s-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                       ! error code
 character(*),intent(out)      :: message                   ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------
 ! local variables: compute flux when coupled to groundwater
 real(dp)                      :: hc_depth                  ! hydraulic conductivity at the depth of the water table (m s-1)
 real(dp)                      :: iHydCond                  ! hydraulic conductivity at the midpoint between the node and the water table (m s-1)
 real(dp)                      :: dHeight                   ! height difference between unsaturated node and water table (m) 
 real(dp)                      :: spGrad                    ! spatial gradient in matric head (-)
 real(dp)                      :: cflux                     ! capillary flux (m s-1)
 ! local variables: compute derivative when coupled to groundwater
 real(dp)                      :: dhc_dzwt                  ! change in hyd cond w.r.t. change in water table depth (s-1)
 real(dp)                      :: dhci_dhcz                 ! change in hydraulic conductivity at the "interface" w.r.t. change in the hydraulic conductivity at the water table depth (-)
 real(dp)                      :: dhci_dzwt                 ! change in hydraulic conductivity at the "interface" w.r.t. change in water table depth (s-1)
 real(dp)                      :: dSpGrad_dzwt              ! change in the spatial derivative w.r.t change in water table depth (m-1)
 real(dp)                      :: dCflux_dzwt               ! change in capillary flux w.r.t. change in the water table depth -- product rule (s-1)
 ! -----------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="rechargeWT/"

 ! error check
 if(scalarWaterTableDepth < nodeHeight)then; err=20; message=trim(message)//'water table is above the mid-point of the soil layer'; return; endif

 ! compute the hydraulic conductivity at the water table depth (m s-1)
 select case(hc_profile)  ! (determine type of hydraulic conductivity profile)
  case(constant);                        hc_depth = surfaceSatHydCond      ! constant, so hyd cond at the depth of the water table = surfaceSatHydCond (m s-1)
  case(exp_profile);                     hc_depth = surfaceSatHydCond * exp(-scalarWaterTableDepth/zScale_TOPMODEL)  ! hyd cond at the depth of the water table (m s-1)
  case(powerLaw_profile,linear_profile); message=trim(message)//'hydraulic conductivity profile not implemented yet for "powerLaw" and "linear" options'; err=10; return
  case default;                          message=trim(message)//'unknown hydraulic conductivity profile'; err=10; return
 end select

 ! get the hydraulic conductivity at the mid-point between (1) the lowest unsaturated node and (2) the water table depth
 iHydCond = (nodeHydCond * hc_depth)**0.5_dp ! NOTE: this is at height = 0.5*(nodeHeight + scalarWaterTableDepth)

 ! compute the capillary flux
 dHeight = scalarWaterTableDepth - nodeHeight
 spGrad  = -(0._dp - nodeMatricHead) / dHeight    ! spatial gradient (-)
 cflux   = iHydCond*spGrad                        ! capillary flux (m s-1)

 ! compute the total flux
 scalarRecharge = cflux + iHydCond
 write(*,'(a,10(e20.10,1x))') 'scalarWaterTableDepth, hc_depth, nodeHydCond, iHydCond, cflux, scalarRecharge = ',&
                               scalarWaterTableDepth, hc_depth, nodeHydCond, iHydCond, cflux, scalarRecharge
 pause

 ! ** compute change in recharge flux w.r.t. change in the aquifer storage
 if(deriv_desired)then

  ! compute change in hyd cond w.r.t. change in water table depth (s-1)
  select case(hc_profile)  ! (determine type of hydraulic conductivity profile)
   case(constant);                        dhc_dzwt = 0._dp
   case(exp_profile);                     dhc_dzwt = -hc_depth/zScale_TOPMODEL
   case(powerLaw_profile,linear_profile); message=trim(message)//'hydraulic conductivity profile not implemented yet for "powerLaw" and "linear" options'; err=10; return
   case default;                          message=trim(message)//'unknown hydraulic conductivity profile'; err=10; return
  end select

  ! compute the change in hydraulic conductivity at the "interface" w.r.t. change in the hydraulic conductivity at the water table depth (-)
  dhci_dhcz    = nodeHydCond/(iHydCond*2._dp)
  ! compute the change in hydraulic conductivity at the "interface" w.r.t. change in water table depth -- chain rule (s-1)
  dhci_dzwt    = dhc_dzwt*dhci_dhcz
  ! compute the change in the spatial derivative w.r.t change in water table depth (m-1)
  dSpGrad_dzwt = -nodeMatricHead/(dHeight**2._dp)
  ! define derivative in capillary flux w.r.t. water table depth -- product rule (s-1)
  dCflux_dzwt  = dhci_dzwt*spGrad + iHydCond*dSpGrad_dzwt

  ! compute final derivative (s-1)
  dq_dWaterTable = dCflux_dzwt + dhci_dzwt

 else     ! (do not desire derivatives)
  dq_dWaterTable = valueMissing
 endif    ! (if desire derivatives)

 end subroutine rechargeWT




 ! ************************************************************************************************
 ! private subroutine: compute the drainage flux from the bottom of the soil profile and its derivative
 ! ************************************************************************************************
 subroutine qDrainFlux(&
                       ! input: model control
                       deriv_desired,             & ! intent(in): flag indicating if derivatives are desired
                       ixRichards,                & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                       hc_profile,                & ! intent(in): index defining the decrease of hydraulic conductivity with depth
                       bc_lower,                  & ! intent(in): index defining the type of boundary conditions
                       ! input: state variables
                       nodeMatricHead,            & ! intent(in): matric head in the lowest unsaturated node (m)
                       nodeVolFracLiq,            & ! intent(in): volumetric liquid water content the lowest unsaturated node (-)
                       ! input: model coordinate variables
                       nodeDepth,                 & ! intent(in): depth of the lowest unsaturated soil layer (m)
                       nodeHeight,                & ! intent(in): height of the lowest unsaturated soil node (m)
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
                       specificYield,             & ! intent(in): drainable porosity (-)
                       ! output: hydraulic conductivity and diffusivity at the surface
                       bottomHydCond,             & ! intent(out): hydraulic conductivity at the bottom of the unsatuarted zone (m s-1)
                       bottomDiffuse,             & ! intent(out): hydraulic diffusivity at the bottom of the unsatuarted zone (m2 s-1)
                       ! output: drainage flux from the bottom of the soil profile
                       scalarDrainage,            & ! intent(out): drainage flux from the bottom of the soil profile (m s-1)
                       ! output: derivatives in drainage flux
                       dq_dStateUnsat,            & ! intent(out): change in drainage flux w.r.t. change in state variable in lowest unsaturated node (m s-1 or s-1)
                       dq_dAquifer,               & ! intent(out): change in drainage flux w.r.t. change in the aquifer storage (s-1)
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
 ! input: model coordinate variables
 real(dp),intent(in)           :: nodeDepth                 ! depth of the lowest unsaturated soil layer (m)
 real(dp),intent(in)           :: nodeHeight                ! height of the lowest unsaturated soil node (m)
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
 real(dp),intent(in)           :: specificYield             ! specific yield (-)
 ! -----------------------------------------------------------------------------------------------------------------------------
 ! output: hydraulic conductivity at the bottom of the unsaturated zone
 real(dp),intent(out)          :: bottomHydCond             ! hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
 real(dp),intent(out)          :: bottomDiffuse             ! hydraulic diffusivity at the bottom of the unsatuarted zone (m2 s-1)
 ! output: drainage flux from the bottom of the soil profile
 real(dp),intent(out)          :: scalarDrainage            ! drainage flux from the bottom of the soil profile (m s-1)
 ! output: derivatives in drainage flux
 real(dp),intent(out)          :: dq_dStateUnsat            ! change in drainage flux w.r.t. change in state variable in lowest unsaturated node (m s-1 or s-1)
 real(dp),intent(out)          :: dq_dAquifer               ! change in drainage flux w.r.t. change in the aquifer storage (s-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                       ! error code
 character(*),intent(out)      :: message                   ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------
 ! local variables
 real(dp)                      :: zWater                    ! effective water table depth (m)
 real(dp)                      :: nodePsi                   ! matric head in the lowest unsaturated node (m)
 real(dp)                      :: cflux                     ! capillary flux (m s-1)
 ! -----------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="qDrainFlux/"

 ! determine lower boundary condition
 select case(bc_lower)

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
   scalarDrainage = cflux + bottomHydCond

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
   scalarDrainage = kAnisotropic*surfaceSatHydCond * exp(-zWater/zScale_TOPMODEL)

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
   scalarDrainage = nodeHydCond*kAnisotropic

   ! compute derivatives
   if(deriv_desired)then
    select case(ixRichards)  ! (form of Richards' equation)
     case(moisture); dq_dStateUnsat = dHydCond_dVolLiq*kAnisotropic
     case(mixdform); dq_dStateUnsat = dHydCond_dMatric*kAnisotropic
     case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
    end select
    ! no dependency on aquifer storage
    dq_dAquifer = 0._dp
   else     ! (do not desire derivatives)
    dq_dStateUnsat = valueMissing
    dq_dAquifer    = valueMissing
   endif


  ! ---------------------------------------------------------------------------------------------
  ! * zero flux
  ! ---------------------------------------------------------------------------------------------
  case(zeroFlux)
   scalarDrainage = 0._dp
   if(deriv_desired)then
    dq_dStateUnsat = 0._dp
    dq_dAquifer    = 0._dp
   else
    dq_dStateUnsat = valueMissing
    dq_dAquifer    = valueMissing
   endif

  ! ---------------------------------------------------------------------------------------------
  ! * error check
  ! ---------------------------------------------------------------------------------------------
  case default; err=20; message=trim(message)//'unknown lower boundary condition for soil hydrology'; return

 endselect ! (type of boundary condition)

 end subroutine qDrainFlux




 ! ************************************************************************************************
 ! private subroutine: compute baseflow flux using a topmodel-type approach
 ! ************************************************************************************************
 subroutine QBtopmodel(&
                       deriv_desired,               & ! input: flag indicating if derivatives are desired
                       scalarWaterTableDepth,       & ! input: water table depth (m)
                       k_surf,                      & ! input: saturated hydraulic conductivity at the surface (m s-1)
                       kAnisotropic,                & ! input: anisotropy factor for lateral hydraulic conductivity (-)
                       zScale_TOPMODEL,             & ! input: scale factor for TOPMODEL-ish baseflow parameterization (m)
                       scalarBaseflow,              & ! output: baseflow from each soil layer (m s-1)
                       scalarBaseflowDeriv,         & ! output: derivative in baseflow flux w.r.t. water table depth (m s-1)
                       err,message)                   ! output: error control
 ! model decision structures
 USE data_struc,only:model_decisions        ! model decision structure
 USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure
 implicit none
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! input
 logical(lgt),intent(in)   :: deriv_desired           ! flag to indicate if derivatives are desired
 real(dp),intent(in)       :: scalarWaterTableDepth   ! trial value of water table depth (m)
 real(dp),intent(in)       :: k_surf                  ! saturated hydraulic conductivity at the surface (m s-1) 
 real(dp),intent(in)       :: kAnisotropic            ! anisotropy factor for lateral hydraulic conductivity (-)
 real(dp),intent(in)       :: zScale_TOPMODEL         ! scale factor for TOPMODEL-ish baseflow parameterization (m)
 ! output
 real(dp),intent(out)      :: scalarBaseflow          ! baseflow (m s-1)
 real(dp),intent(out)      :: scalarBaseflowDeriv     ! derivative in baseflow flux w.r.t. water table depth (m s-1)
 integer(i4b),intent(out)  :: err                     ! error code
 character(*),intent(out)  :: message                 ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='QBtopmodel/'
 ! compute the baseflow (m s-1)
 select case(model_decisions(iLookDECISIONS%hc_profile)%iDecision)
  ! (exponential transmissivity profile)
  case(exp_profile)
   scalarBaseflow = k_surf*kAnisotropic*exp(-scalarWaterTableDepth/zScale_TOPMODEL)
   if(deriv_desired)then
    scalarBaseflowDeriv = -scalarBaseflow/zScale_TOPMODEL
   else
    scalarBaseflowDeriv = valueMissing
   endif
  ! (linear transmissivity profile)
  case(linear_profile)
   scalarBaseflow = k_surf*kAnisotropic*(1._dp - scalarWaterTableDepth/zScale_TOPMODEL)
   if(deriv_desired)then
    scalarBaseflowDeriv = -k_surf*kAnisotropic/zScale_TOPMODEL
   else
    scalarBaseflowDeriv = valueMissing
   endif
  ! (constant transmissivity profile)
  case(constant)
   scalarBaseflow = k_surf*kAnisotropic*exp(-scalarWaterTableDepth/zScale_TOPMODEL)
   if(deriv_desired)then
    scalarBaseflowDeriv = -scalarBaseflow/zScale_TOPMODEL
   else
    scalarBaseflowDeriv = valueMissing
   endif
  ! (power-law transmissivity profile)
  case(powerLaw_profile)
   message=trim(message)//"hydraulic conductivity profile not implemented yet [option="//trim(model_decisions(iLookDECISIONS%hc_profile)%cDecision)//"]"
   err=20; return
  ! (unknown transmissivity profile)
  case default
   message=trim(message)//"unknown hydraulic conductivity profile [option="//trim(model_decisions(iLookDECISIONS%hc_profile)%cDecision)//"]"
   err=20; return
 end select
 end subroutine QBtopmodel




 ! ************************************************************************************************
 ! private subroutine: compute baseflow flux using a conceptual bog bucket
 ! ************************************************************************************************
 subroutine QBbigbuckt(&
                       deriv_desired,               &    ! input: flag indicating if derivatives are desired
                       scalarAquiferStorageTrial,   &    ! input: trial value of aquifer storage (m)
                       aquiferScaleFactor,          &    ! input: scaling factor for aquifer storage in the big bucket (m)
                       k_surf,                      &    ! input: saturated hydraulic conductivity at the surface (m s-1)
                       kAnisotropic,                &    ! input: anisotropy factor for lateral hydraulic conductivity (-)
                       bucketBaseflowExp,           &    ! input: exponent in bucket baseflow parameterization
                       scalarAquiferBaseflow,       &    ! output: total baseflow (m s-1)
                       scalarAquiferBaseflowDeriv,  &    ! output: derivative in baseflow flux w.r.t. aquifer storage (s-1)
                       err,message)                      ! output: error control
 implicit none
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! input
 logical(lgt),intent(in)   :: deriv_desired              ! flag to indicate if derivatives are desired
 real(dp),intent(in)       :: scalarAquiferStorageTrial  ! trial value of water table depth (m)
 real(dp),intent(in)       :: aquiferScaleFactor         ! scaling factor for aquifer storage in the big bucket (m)
 real(dp),intent(in)       :: k_surf                     ! saturated hydraulic conductivity at the surface (m s-1) 
 real(dp),intent(in)       :: kAnisotropic               ! anisotropy factor for lateral hydraulic conductivity (-)
 real(dp),intent(in)       :: bucketBaseflowExp          ! exponent in bucket baseflow parameterization
 ! output
 real(dp),intent(out)      :: scalarAquiferBaseflow      ! baseflow flux (m s-1)
 real(dp),intent(out)      :: scalarAquiferBaseflowDeriv ! derivative in baseflow flux w.r.t. water table depth (s-1)
 integer(i4b),intent(out)  :: err                        ! error code
 character(*),intent(out)  :: message                    ! error message
 ! local
 real(dp)                  :: scaledStorage              ! scaled storage (-)
 real(dp)                  :: maxBaseflowRate            ! maximum baseflow rate (m s-1)
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='QBbigbuckt/'
 ! get temporary variables
 scaledStorage    = scalarAquiferStorageTrial/aquiferScaleFactor
 maxBaseflowRate  = kAnisotropic*k_surf
 ! compute baseflow flux (m s-1)
 scalarAquiferBaseflow = maxBaseflowRate*scaledStorage**bucketBaseflowExp
 ! compute derivative in baseflow flux w.r.t. aquifer storage (s-1)
 if(deriv_desired)then
  scalarAquiferBaseflowDeriv = (maxBaseflowRate/aquiferScaleFactor)*bucketBaseflowExp*scaledStorage**(bucketBaseflowExp - 1._dp)
 else
  scalarAquiferBaseflowDeriv = valueMissing
 endif
 end subroutine QBbigbuckt




 ! ************************************************************************************************
 ! private subroutine: compute flux of water ejected because close to exceeding porosity
 ! ************************************************************************************************
 subroutine ejectWater(&
                       ! input: model control
                       deriv_desired,         & ! intent(in): flag determining if the derivative is desired
                       ixRichards,            & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                       ! input: state variable and diagnostic variables
                       scalarMatricHeadTrial, & ! intent(in): matric head in each layer (m)
                       scalarVolFracLiqTrial, & ! intent(in): volumetric liquid water content (-)
                       scalarVolFracIceTrial, & ! intent(in): volumetric ice content in each layer (-)
                       ! input: derivative in the soil water characteristic w.r.t. psi
                       scalardTheta_dPsi,     & ! intent(in): derivative in the soil water characteristic w.r.t. psi (m-1)
                       ! input: soil parameters
                       vGn_alpha,             & ! intent(in): van Genutchen "alpha" parameter (m-1)
                       vGn_n,                 & ! intent(in): van Genutchen "n" parameter (-)
                       VGn_m,                 & ! intent(in): van Genutchen "m" parameter (-)
                       theta_sat,             & ! intent(in): soil porosity (-)
                       theta_res,             & ! intent(in): soil residual volumetric water content (-)
                       ! input: saturated hydraulic conductivity
                       scalarSatHydCond,      & ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
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
 integer(i4b),intent(in)       :: ixRichards                ! index defining the form of Richards' equation (moisture or mixdform)
 ! input: state variables
 real(dp),intent(in)           :: scalarMatricHeadTrial     ! matric head in each layer (m)
 real(dp),intent(in)           :: scalarVolFracLiqTrial     ! volumetric liquid water content in each layer (-)
 real(dp),intent(in)           :: scalarVolFracIceTrial     ! volumetric ice content in each layer (-)
 ! input: derivative in the soil water characteristic w.r.t. psi
 real(dp),intent(in)           :: scalardTheta_dPsi         ! derivative in the soil water characteristic w.r.t. psi (m-1)
 ! input: soil parameters
 real(dp),intent(in)           :: vGn_alpha                 ! van Genutchen "alpha" parameter (m-1)
 real(dp),intent(in)           :: vGn_n                     ! van Genutchen "n" parameter (-)
 real(dp),intent(in)           :: vGn_m                     ! van Genutchen "m" parameter (-)
 real(dp),intent(in)           :: theta_sat                 ! soil porosity (-)
 real(dp),intent(in)           :: theta_res                 ! soil residual volumetric water content (-)
 ! input: saturated hydraulic conductivity
 real(dp),intent(in)           :: scalarSatHydCond          ! saturated hydraulic conductivity at the layer interface (m s-1)
 ! output: ejected water flux and derivative
 real(dp),intent(out)          :: scalarEjectWater          ! water ejected because pore volume is filled (m s-1)
 real(dp),intent(out)          :: scalarEjectWaterDeriv     ! derivative in ejected water flux (m s-1 [moisture form] s-1 [mixed form])
 ! output: error control
 integer(i4b),intent(out)      :: err                       ! error code
 character(*),intent(out)      :: message                   ! error message
 ! local variables
 real(dp),parameter            :: supersatScale=0.01_dp    ! scaling factor for the water ejection function (-)
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
 !print*, 'supersatThresh = ', supersatThresh
 !pause

 ! define minimum value for calculations
 fracMin = -supersatScale*log(1._dp/fSmall - 1._dp) + supersatThresh

 ! calculate the fraction of pore space filled with liquid water and ice (-)
 select case(ixRichards)
  case(moisture); fracCap = (scalarVolFracLiqTrial + scalarVolFracIceTrial)/theta_sat
  case(mixdform); fracCap = (volFracLiq(scalarMatricHeadTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) + scalarVolFracIceTrial)/theta_sat
  case default; err=10; message=trim(message)//"unable to identify option for Richards' equation"; return
 endselect

 if(fracCap > fracMin)then  ! (check if the fractional capacity is greater than the minimum value)

  ! ***** calculate the upward flux associated with water ejected (m s-1)
  expFunc = exp((supersatThresh - fracCap)/supersatScale)
  scalarEjectWater = scalarSatHydCond/(1._dp + expFunc)
  if(scalarEjectWater < 0._dp)then; message=trim(message)//'ejected water < 0'; err=20; return; endif
  !write(*,'(a,10(e20.10,1x))') 'scalarVolFracLiqTrial, scalarVolFracIceTrial, fracCap, scalarEjectWater = ', &
  !                              scalarVolFracLiqTrial, scalarVolFracIceTrial, fracCap, scalarEjectWater

  ! ***** calculate the derivative in upward ejected flux (m s-1 [moisture form] or s-1 [mixed form])
  if(deriv_desired)then
   select case(ixRichards)
    case(moisture); scalarEjectWaterDeriv = (expFunc/(theta_sat*supersatScale)) * (1._dp + expFunc)**(-2._dp) * scalarSatHydCond
    case(mixdform); scalarEjectWaterDeriv = scalardTheta_dPsi * (expFunc/(theta_sat*supersatScale)) * (1._dp + expFunc)**(-2._dp) * scalarSatHydCond
    case default; err=10; message=trim(message)//"unable to identify option for Richards' equation"; return
   end select
  else  ! (if derivative is not desired)
   scalarEjectWaterDeriv = valueMissing
  endif

 else  ! (if fraction of pore space filled with liquid water and ice is less than the minimum value required for lateral flow)
  scalarEjectWater = 0._dp
  scalarEjectWaterDeriv  = 0._dp
 endif ! (if fraction of pore space filled with liquid water and ice is less than the minimum value required for lateral flow)

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
                        ! input: model parameters for the compressibility term
                        theta_sat,                       & ! intent(in): porosity (-)
                        specificStorage,                 & ! intent(in): specific storage (m-1)
                        ! initial flux vectors (start of the time step)
                        iLayerInitLiqFluxSoil,           & ! intent(in): initial liquid flux at layer interfaces (m s-1)
                        mLayerInitEjectWater,            & ! intent(in): initial water ejected from each layer because volumetric liquid water is approaching porosity (m s-1)
                        mLayerInitBaseflow,              & ! intent(in): initial baseflow from each layer (m s-1)
                        mLayerInitTranspire,             & ! intent(in): initial transpiration from each layer -- from energy routine (m s-1)
                        ! trial flux vectors
                        iLayerTrialLiqFluxSoil,          & ! intent(in): trial liquid flux at layer interfaces (m s-1)
                        mLayerTrialEjectWater,           & ! intent(in): trial water ejected from each layer because volumetric liquid water is approaching porosity (m s-1)
                        mLayerTrialBaseflow,             & ! intent(in): trial baseflow from each layer (m s-1)
                        mLayerTrialTranspire,            & ! intent(in): trial transpiration from each layer -- from energy routine (m s-1)
                        ! initial state vectors (start of the time step)
                        mLayerInitMatricHead,            & ! intent(in): initial matric head (m)
                        mLayerInitVolFracLiq,            & ! intent(in): initial volumetric liquid water content (-)
                        mLayerInitVolFracIce,            & ! intent(in): initial volumetric ice content (-)
                        ! trial state vectors
                        mLayerTrialMatricHead,           & ! intent(in): trial matric head (m)
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
 ! model parameters for the compressibility term
 real(dp),intent(in)          :: theta_sat                 ! intent(in): porosity (-)
 real(dp),intent(in)          :: specificStorage           ! intent(in): specific storage (m-1)
 ! initial flux vectors (start of the time step)
 real(dp),intent(in)          :: iLayerInitLiqFluxSoil(0:) ! intent(in): initial liquid flux at layer interfaces (m s-1)
 real(dp),intent(in)          :: mLayerInitEjectWater(:)   ! intent(in): initial water ejected from each layer because volumetric liquid water is approaching porosity (m s-1)
 real(dp),intent(in)          :: mLayerInitBaseflow(:)     ! intent(in): initial baseflow from each layer (m s-1)
 real(dp),intent(in)          :: mLayerInitTranspire(:)    ! intent(in): initial transpiration from each layer -- from energy routine (m s-1)
 ! trial flux vectors
 real(dp),intent(in)          :: iLayerTrialLiqFluxSoil(0:) ! intent(in): trial liquid flux at layer interfaces (m s-1)
 real(dp),intent(in)          :: mLayerTrialEjectWater(:)  ! intent(in): trial water ejected from each layer because volumetric liquid water is approaching porosity (m s-1)
 real(dp),intent(in)          :: mLayerTrialBaseflow(:)    ! intent(in): trial baseflow from each layer (m s-1)
 real(dp),intent(in)          :: mLayerTrialTranspire(:)   ! intent(in): trial transpiration from each layer -- from energy routine (m s-1)
 ! initial state vectors (start of the time step)
 real(dp),intent(in)          :: mLayerInitMatricHead(:)   ! intent(in): initial matric head (m)
 real(dp),intent(in)          :: mLayerInitVolFracLiq(:)   ! intent(in): initial volumetric liquid water content (-)
 real(dp),intent(in)          :: mLayerInitVolFracIce(:)   ! intent(in): initial volumetric ice content (-)
 ! trial state vectors
 real(dp),intent(in)          :: mLayerTrialMatricHead(:)  ! intent(in): trial matric head (m)
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
 real(dp)                     :: mBase                     ! overall contribution to volumetric liquid water from baseflow (-)
 real(dp)                     :: mPhse                     ! overall contribution to volumetric liquid water from phase change (-) 
 real(dp)                     :: compressibility           ! compressibility term (-)
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
  ! compressibility (-)
  compressibility = 0._dp !(specificStorage*mLayerTrialVolFracLiq(iLayer)/theta_sat) * (mLayerInitMatricHead(iLayer) - mLayerTrialMatricHead(iLayer))
  ! transpiration (-)
  mEvap = wimplicit*mLayerInitTranspire(iLayer)*dt_dz + (1._dp - wimplicit)*mLayerTrialTranspire(iLayer)*dt_dz
  ! ejected water (-)
  mEjct = wimplicit*mLayerInitEjectWater(iLayer)*dt_dz + (1._dp - wimplicit)*mLayerTrialEjectWater(iLayer)*dt_dz
  ! baseflow (-)
  mBase = wimplicit*mLayerInitBaseflow(iLayer)*dt_dz + (1._dp - wimplicit)*mLayerTrialBaseflow(iLayer)*dt_dz
  ! phase change (-)
  mPhse = (iden_ice/iden_water)*(mLayerTrialVolFracIce(iLayer)-mLayerInitVolFracIce(iLayer))
  ! residual (-)
  residualVec(iLayer) = mLayerTrialVolFracLiq(iLayer) - (mLayerInitVolFracLiq(iLayer) + mFlux + mEvap - mEjct - mBase - mPhse - compressibility)

  ! print progress
  !if(iLayer==1)    write(*,'(a)') 'iLayer, residualVec(iLayer), mLayerTrialVolFracLiq(iLayer), mLayerInitVolFracLiq(iLayer), bottom flux, mFlux, mEvap, mEjct, mBase, mPhse'
  !if(iLayer < 5) write(*,'(i4,1x,10(e20.10,1x))') iLayer, residualVec(iLayer), mLayerTrialVolFracLiq(iLayer), mLayerInitVolFracLiq(iLayer), iLayerTrialLiqFluxSoil(iLayer)*dt_dz, mFlux, mEvap, mEjct, mBase, mPhse
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
                        scalarInitAquiferRecharge, & ! intent(in): recharge to the aquifer        (m s-1)
                        scalarInitAquiferBaseflow, & ! intent(in): baseflow from the aquifer      (m s-1)
                        ! input: end-of-step fluxes
                        scalarAquiferTranspire,    & ! intent(in): transpiration from the aquifer (m s-1)
                        scalarAquiferRecharge,     & ! intent(in): recharge to the aquifer        (m s-1)
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
 real(dp),intent(in)          :: scalarInitAquiferRecharge   ! aquifer recharge averaged over the time step (m s-1)
 real(dp),intent(in)          :: scalarInitAquiferBaseflow   ! baseflow from the aquifer (m s-1) 
 ! input: end-of-step fluxes
 real(dp),intent(in)          :: scalarAquiferTranspire      ! aquifer transpiration averaged over the time step (m s-1)
 real(dp),intent(in)          :: scalarAquiferRecharge       ! aquifer recharge averaged over the time step (m s-1)
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
 netFlux0 = (scalarInitAquiferRecharge + scalarInitAquiferTranspire) - scalarInitAquiferBaseflow
 ! compute the flux at the end of the step
 netFlux1 = (scalarAquiferRecharge + scalarAquiferTranspire) - scalarAquiferBaseflow
 ! compute the mean net flux
 netFlux  = wimplicit*netFlux0 + (1._dp - wimplicit)*netFlux1
 ! compute the residual
 scalarAquiferResidual = scalarAquiferStorageTrial - (scalarAquiferStorage + netFlux*dt)
 
 !write(*,'(a)') '***** computing the aquifer residual'
 !write(*,'(a,e20.10)') 'dt = ', dt
 !write(*,'(a,e20.10)') 'scalarAquiferStorage = ', scalarAquiferStorage
 !write(*,'(a,e20.10)') 'scalarAquiferStorageTrial = ', scalarAquiferStorageTrial
 !write(*,'(a,2(e20.10,1x))') 'transpire = ', scalarInitAquiferTranspire, scalarAquiferTranspire
 !write(*,'(a,2(e20.10,1x))') 'recharge = ',  scalarInitAquiferRecharge, scalarAquiferRecharge
 !write(*,'(a,2(e20.10,1x))') 'baseflow = ', scalarInitAquiferBaseflow, scalarAquiferBaseflow
 !write(*,'(a,3(e20.10,1x))') 'netFlux0, netFlux1, netFlux = ', netFlux0, netFlux1, netFlux
 !write(*,'(a,e20.10)') 'trial state = ', scalarAquiferStorage + netFlux*dt
 !write(*,'(a,e20.10)') 'scalarAquiferResidual = ', scalarAquiferResidual
 !write(*,'(a)') '*******************************************************************************************************'
 end subroutine gw_residual




end module soilHydrol_module

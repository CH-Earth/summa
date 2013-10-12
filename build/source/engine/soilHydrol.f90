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
 powerLaw_profile,           & ! power-law profile
 ! look-up values for the choice of groundwater parameterization
 qbaseTopmodel,              & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                  & ! a big bucket (lumped aquifer model)
 noExplicit,                 & ! no explicit groundwater parameterization
 ! look-up values for the choice of boundary conditions for hydrology
 prescribedHead,             & ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
 funcBottomHead,             & ! function of matric head in the lower-most layer
 freeDrainage,               & ! free drainage
 liquidFlux,                 & ! liquid water flux
 zeroFlux,                   & ! zero flux
 ! look-up values for the choice of groundwater representation (local-column, or single-basin)
 localColumn,                & ! separate groundwater representation in each local soil column
 singleBasin                   ! single groundwater store over the entire basin
! -----------------------------------------------------------------------------------------------------------
implicit none
private
public::soilHydrol
! number of layers
integer(i4b)           :: nSoil                     ! number of soil layers
integer(i4b)           :: nSnow                     ! number of snow layers
integer(i4b)           :: nLevels                   ! total number of soil layers to examine
! constant parameters
real(dp),parameter     :: valueMissing=-9999._dp    ! missing value parameter
real(dp),parameter     :: verySmall=1.e-12_dp       ! a very small number (used to avoid divide by zero)
real(dp),parameter     :: dx=1.e-8_dp               ! finite difference increment
contains


 ! ************************************************************************************************
 ! new subroutine: compute change in volumetric liquid water content (or matric head) over the time step
 ! ************************************************************************************************
 subroutine soilHydrol(&
                       ! input
                       dt,                        & ! intent(in): time step (seconds)
                       iter,                      & ! intent(in): iteration index
                       ixBcUpperSoilHydrology,    & ! intent(in): choice of upper boundary condition for soil hydrology
                       upperBoundHead,            & ! intent(in): upper boundary condition for matric head (m)
                       upperBoundTheta,           & ! intent(in): upper boundary condition for volumetric liquid water content (-)
                       scalarCanopyTranspiration, & ! intent(in): canopy transpiration (kg m-2 s-1)
                       scalarGroundEvaporation,   & ! intent(in): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
                       mLayerMatricHeadIter,      & ! intent(in): matric head in each layer at the current iteration (m)
                       mLayerVolFracLiqIter,      & ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                       mLayerVolFracIceIter,      & ! intent(in): volumetric fraction of ice at the current iteration (-)
                       scalarAquiferStorageIter,  & ! intent(in): aquifer storage (m)
                       mLayerMatricIncrOld,       & ! intent(in): iteration increment for matric head of soil layers in the previous iteration (m)
                       mLayerLiquidIncrOld,       & ! intent(in): iteration increment for volumetric liquid water content of soil layers in the previous iteration (-) 
                       ! output
                       mLayerMatricHeadNew,       & ! intent(out): matric head in each layer at the next iteration (m)
                       mLayerVolFracLiqNew,       & ! intent(out): volumetric fraction of liquid water at the next iteration (-)
                       scalarAquiferStorageNew,   & ! intent(out): aquifer storage (m)
                       scalarSurfaceInfiltration, & ! intent(out): surface infiltration rate (m s-1)
                       err,message)
 ! model decision structures
 USE data_struc,only:model_decisions        ! model decision structure
 USE var_lookup,only:iLookDECISIONS         ! named variables for elements of the decision structure
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:mpar_data,forc_data,attr_data,mvar_data,indx_data    ! data structures
 USE var_lookup,only:iLookPARAM,iLookFORCE,iLookATTR,iLookMVAR,iLookINDEX ! named variables for structure elements
 ! compute change in soil moisture or matric head over the time step
 implicit none
 ! input
 real(dp),intent(in)           :: dt                        ! time step (seconds)
 integer(i4b),intent(in)       :: iter                      ! iteration index
 integer(i4b),intent(in)       :: ixBcUpperSoilHydrology    ! choice of upper boundary condition for soil hydrology
 real(dp),intent(in)           :: upperBoundHead            ! upper boundary condition for matric head (m)
 real(dp),intent(in)           :: upperBoundTheta           ! upper boundary condition for volumetric liquid water content (-)
 real(dp),intent(in)           :: scalarCanopyTranspiration ! canopy transpiration (kg m-2 s-1)
 real(dp),intent(in)           :: scalarGroundEvaporation   ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
 real(dp),intent(in)           :: mLayerMatricHeadIter(:)   ! matric head in each layer at the current iteration (m)
 real(dp),intent(in)           :: mLayerVolFracLiqIter(:)   ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),intent(in)           :: mLayerVolFracIceIter(:)   ! volumetric fraction of ice at the current iteration (-)
 real(dp),intent(in)           :: scalarAquiferStorageIter  ! aquifer storage (m)
 real(dp),intent(in)           :: mLayerMatricIncrOld(:)    ! iteration increment for matric head of soil layers in the previous iteration (m)
 real(dp),intent(in)           :: mLayerLiquidIncrOld(:)    ! iteration increment for volumetric liquid water content of soil layers in the previous iteration (-)
 ! output
 real(dp),intent(out)          :: mLayerMatricHeadNew(:)    ! matric head in each layer at the next iteration (m)
 real(dp),intent(out)          :: mLayerVolFracLiqNew(:)    ! volumetric fraction of liquid water at the next iteration (-)
 real(dp),intent(out)          :: scalarAquiferStorageNew   ! aquifer storage (m)
 real(dp),intent(out)          :: scalarSurfaceInfiltration ! surface infiltration rate (m s-1)
 integer(i4b),intent(out)      :: err                       ! error code
 character(*),intent(out)      :: message                   ! error message
 ! internal
 integer(i4b)                  :: ibeg,iend                 ! start and end indices of the soil layers in concatanated snow-soil vector
 integer(i4b)                  :: local_ixGroundwater       ! local representation of groundwater (gw decision applies to regional scale, if regional gw used)
 character(LEN=256)            :: cmessage                  ! error message of downwind routine
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

 ! modify the groundwater representation for this single-column implementation
 select case(model_decisions(iLookDECISIONS%spatial_gw)%iDecision)
  case(singleBasin); local_ixGroundwater = noExplicit ! force no explicit representation of groundwater at the local scale
  case(localColumn); local_ixGroundwater = model_decisions(iLookDECISIONS%groundwatr)%iDecision ! go with the specified decision
  case default; err=20; message=trim(message)//'unable to identify spatial representation of groundwater'; return
 end select ! (modify the groundwater representation for this single-column implementation) 



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
                        ixBcUpperSoilHydrology,                                       & ! intent(in): choice of upper boundary condition for soil hydrology
                        upperBoundHead,                                               & ! intent(in): upper boundary condition for matric head (m)
                        upperBoundTheta,                                              & ! intent(in): upper boundary condition for volumetric liquid water content (-)
                        scalarCanopyTranspiration,                                    & ! intent(in): canopy transpiration (kg m-2 s-1)
                        scalarGroundEvaporation,                                      & ! intent(in): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
                        mLayerMatricHeadIter,                                         & ! intent(in): matric head in each layer at the current iteration (m)
                        mLayerVolFracLiqIter,                                         & ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                        mLayerVolFracIceIter,                                         & ! intent(in): volumetric fraction of ice at the current iteration (-)
                        scalarAquiferStorageIter,                                     & ! intent(in): aquifer storage (m)
                        mLayerMatricIncrOld,                                          & ! intent(in): iteration increment for matric head of soil layers in the previous iteration (m)
                        mLayerLiquidIncrOld,                                          & ! intent(in): iteration increment for volumetric liquid water content of soil layers in the previous iteration (-) 

                        ! named variables for model decisions
                        model_decisions(iLookDECISIONS%fDerivMeth)%iDecision,         & ! intent(in): method used to calculate flux derivatives 
                        model_decisions(iLookDECISIONS%f_Richards)%iDecision,         & ! intent(in): form of Richards' equation
                        model_decisions(iLookDECISIONS%hc_Profile)%iDecision,         & ! intent(in): option for the hydraulic conductivity profile
                        local_ixGroundwater,                                          & ! intent(in): groundwater parameterization
                        model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision,         & ! intent(in): lower boundary conditions for soil hydrology
                        model_decisions(iLookDECISIONS%spatial_gw)%iDecision,         & ! intent(in): spatial representation of groundwater (local-column or single-basin)

                        ! model coordinate variables -- NOTE: use of ibeg and iend 
                        mvar_data%var(iLookMVAR%mLayerDepth)%dat(ibeg:iend),          & ! intent(in): depth of the layer (m)
                        mvar_data%var(iLookMVAR%mLayerHeight)%dat(ibeg:iend),         & ! intent(in): height of the layer mid-point (m)
                        mvar_data%var(iLookMVAR%iLayerHeight)%dat(ibeg-1:iend),       & ! intent(in): height of the layer interfaces (m)

                        ! lower boundary conditions
                        mpar_data%var(iLookPARAM%lowerBoundHead),                     & ! intent(in): lower boundary condition for matric head (m)
                        mpar_data%var(iLookPARAM%lowerBoundTheta),                    & ! intent(in): lower boundary condition for volumetric liquid water content (-)

                        ! model forcing
                        mvar_data%var(iLookMVAR%scalarThroughfallRain)%dat(1),        & ! intent(in): computed throughfall rate (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarCanopyLiqDrainage)%dat(1),      & ! intent(in): computed drainage of liquid water (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%iLayerLiqFluxSnow)%dat(nSnow),        & ! intent(in): liquid flux from the base of the snowpack (m s-1)

                        ! column inflow
                        mvar_data%var(iLookMVAR%mLayerColumnInflow)%dat(1:nLevels),    & ! total inflow to each layer in the soil column (m3 s-1)

                        ! local attributes
                        attr_data%var(iLookATTR%HRUarea),                             & ! intent(in): HRU area (m2)
                        attr_data%var(iLookATTR%tan_slope),                           & ! intent(in): tan water table slope, taken as tan local ground surface slope (-)
                        attr_data%var(iLookATTR%contourLength),                       & ! intent(in): length of contour at downslope edge of HRU (m)
 
                        ! general model parameters
                        mpar_data%var(iLookPARAM%wimplicit),                          & ! intent(in): weight assigned to start-of-step fluxes (-)

                        ! soil parameters
                        mpar_data%var(iLookPARAM%vGn_alpha),                          & ! intent(in): van Genutchen "alpha" parameter (m-1)
                        mpar_data%var(iLookPARAM%vGn_n),                              & ! intent(in): van Genutchen "n" parameter (-)
                        mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1),                  & ! intent(in): van Genutchen "m" parameter (-)
                        mpar_data%var(iLookPARAM%mpExp),                              & ! intent(in): empirical exponent in macropore flow equation (-)
                        mpar_data%var(iLookPARAM%theta_mp),                           & ! intent(in): volumetric liquid water content when macropore flow begins (-)
                        mpar_data%var(iLookPARAM%theta_sat),                          & ! intent(in): soil porosity (-)
                        mpar_data%var(iLookPARAM%theta_res),                          & ! intent(in): soil residual volumetric water content (-)
                        mpar_data%var(iLookPARAM%fieldCapacity),                      & ! intent(in): field capacity (-)
                        mpar_data%var(iLookPARAM%kAnisotropic),                       & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                        mpar_data%var(iLookPARAM%zScale_TOPMODEL),                    & ! intent(in): TOPMODEL scaling factor (m)
                        mpar_data%var(iLookPARAM%qSurfScale),                         & ! intent(in): scaling factor in the surface runoff parameterization (-)
                        mpar_data%var(iLookPARAM%specificYield),                      & ! intent(in): specific yield (-)
                        mpar_data%var(iLookPARAM%specificStorage),                    & ! intent(in): specific storage coefficient (m-1)
                        mpar_data%var(iLookPARAM%aquiferScaleFactor),                 & ! intent(in): scaling factor for aquifer storage in the big bucket (m)
                        mpar_data%var(iLookPARAM%aquiferBaseflowExp),                 & ! intent(in): baseflow exponent for the big bucket (-)
                        mpar_data%var(iLookPARAM%f_impede),                           & ! intent(in): ice impedence factor (-)

                        ! model state variables -- NOTE: use of ibeg and iend
                        mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(ibeg:iend),     & ! intent(in): volumetric fraction of ice in each layer (-)
                        mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(ibeg:iend),     & ! intent(in): volumetric fraction of liquid water in each layer (-)
                        mvar_data%var(iLookMVAR%mLayerMatricHead)%dat(1:nLevels),     & ! intent(in): matric head in each layer (m)
                        mvar_data%var(iLookMVAR%scalarAquiferStorage)%dat(1),         & ! intent(in): (scalar) aquifer storage at the start of the time step (m)
                        mvar_data%var(iLookMVAR%scalarSfcMeltPond)%dat(1),            & ! intent(in): (scalar) ponded water caused by melt of the "snow without a layer" (kg m-2)

                        ! saturated hydraulic conductivity
                        mvar_data%var(iLookMVAR%mLayerSatHydCondMP)%dat,              & ! intent(in): saturated hydraulic conductivity of macropores at the mid-point of each layer (m s-1)
                        mvar_data%var(iLookMVAR%mLayerSatHydCond)%dat,                & ! intent(in): saturated hydraulic conductivity at the mid-point of each layer (m s-1)
                        mvar_data%var(iLookMVAR%iLayerSatHydCond)%dat,                & ! intent(in): saturated hydraulic conductivity at the interface of each layer (m s-1)

                        ! factors limiting transpiration (from vegFlux routine)
                        mvar_data%var(iLookMVAR%mLayerRootDensity)%dat,               & ! intent(in): root density in each layer (-)
                        mvar_data%var(iLookMVAR%scalarAquiferRootFrac)%dat(1),        & ! intent(in): fraction of roots below the lowest soil layer (-)
                        mvar_data%var(iLookMVAR%scalarTranspireLim)%dat(1),           & ! intent(in): weighted average of the transpiration limiting factor (-)
                        mvar_data%var(iLookMVAR%mLayerTranspireLim)%dat,              & ! intent(in): transpiration limiting factor in each layer (-)
                        mvar_data%var(iLookMVAR%scalarTranspireLimAqfr)%dat(1),       & ! intent(in): transpiration limiting factor for the aquifer (-)

                        ! initial fluxes (intent inout)
                        mvar_data%var(iLookMVAR%iLayerInitLiqFluxSoil)%dat(0:nLevels),& ! intent(inout): liquid flux at layer interfaces at the start of the time step (m s-1)
                        mvar_data%var(iLookMVAR%mLayerInitQMacropore)%dat(1:nLevels), & ! intent(inout): liquid water flux to/from the macropore at the start-of-step (m s-1)
                        mvar_data%var(iLookMVAR%mLayerInitTranspire)%dat(1:nLevels),  & ! intent(inout): transpiration loss from each soil layer at the start of the time step (m s-1)
                        mvar_data%var(iLookMVAR%mLayerInitBaseflow)%dat(1:nLevels),   & ! intent(inout): baseflow from each soil layer at the start-of-step (m s-1)
                        mvar_data%var(iLookMVAR%scalarInitAquiferTranspire)%dat(1),   & ! intent(inout): transpiration loss from the aquifer at the start-of-step (m s-1)
                        mvar_data%var(iLookMVAR%scalarInitAquiferRecharge)%dat(1),    & ! intent(inout): (scalar) recharge to the aquifer at the start-of-step (m s-1)
                        mvar_data%var(iLookMVAR%scalarInitAquiferBaseflow)%dat(1),    & ! intent(inout): (scalar) baseflow from the aquifer at the start-of-step (m s-1)

                        ! diagnostic scalar variables
                        mvar_data%var(iLookMVAR%scalarRainPlusMelt)%dat(1),           & ! intent(out): (scalar) rain plus melt (m s-1)
                        mvar_data%var(iLookMVAR%scalarSurfaceRunoff)%dat(1),          & ! intent(out): (scalar) surface runoff (m s-1)
                        mvar_data%var(iLookMVAR%scalarAquiferTranspire)%dat(1),       & ! intent(out): (scalar) transpiration loss from the aquifer (m s-1)
                        mvar_data%var(iLookMVAR%scalarAquiferRecharge)%dat(1),        & ! intent(out): (scalar) recharge to the aquifer at the end-of-step (m s-1)
                        mvar_data%var(iLookMVAR%scalarAquiferBaseflow)%dat(1),        & ! intent(out): (scalar) baseflow from the aquifer at the end-of-step (m s-1)

                        ! model diagnostic variables
                        mvar_data%var(iLookMVAR%mLayerdTheta_dPsi)%dat(1:nLevels),    & ! intent(out): derivative in the soil water characteristic w.r.t. psi (m-1)
                        mvar_data%var(iLookMVAR%mLayerdPsi_dTheta)%dat(1:nLevels),    & ! intent(out): derivative in the soil water characteristic w.r.t. theta (m)
                        mvar_data%var(iLookMVAR%mLayerColumnOutflow)%dat(1:nLevels),  & ! intent(out): total outflow from each layer of the soil column (m3 s-1)
                        mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat(0:nLevels),    & ! intent(out): liquid flux at layer interfaces at the end of the time step (m s-1)
                        mvar_data%var(iLookMVAR%mLayerQMacropore)%dat(1:nLevels),     & ! intent(out): liquid water flux to/from the macropore (m s-1)
                        mvar_data%var(iLookMVAR%mLayerTranspire)%dat(1:nLevels),      & ! intent(out): transpiration loss from each soil layer (m s-1)
                        mvar_data%var(iLookMVAR%mLayerBaseflow)%dat(1:nLevels),       & ! intent(inout): baseflow from each soil layer (m s-1)

                        ! output variables from the soilHydrol routine  -- NOTE: variables are already sized appropriately
                        mLayerMatricHeadNew,                                          & ! intent(out): matric head in each layer at the next iteration (m)
                        mLayerVolFracLiqNew,                                          & ! intent(out): volumetric fraction of liquid water at the next iteration (-)
                        scalarAquiferStorageNew,                                      & ! intent(out): aquifer storage (m)
                        scalarSurfaceInfiltration,                                    & ! intent(out): surface infiltration rate (m s-1)
                        err,cmessage)                                                   ! intent(out): error control
 ! check for errors
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! set un-used portion of diagnostic variables to zero
 if(nLevels < nSoil)then
  mvar_data%var(iLookMVAR%mLayerdTheta_dPsi)%dat(    nLevels+1:nSoil) = 0._dp
  mvar_data%var(iLookMVAR%mLayerdPsi_dTheta)%dat(    nLevels+1:nSoil) = 0._dp
  mvar_data%var(iLookMVAR%iLayerInitLiqFluxSoil)%dat(nLevels+1:nSoil) = 0._dp
  mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat(    nLevels+1:nSoil) = 0._dp
  mvar_data%var(iLookMVAR%mLayerInitQMacropore)%dat( nLevels+1:nSoil) = 0._dp
  mvar_data%var(iLookMVAR%mLayerQMacropore)%dat(     nLevels+1:nSoil) = 0._dp
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
                              ixBcUpperSoilHydrology,      & ! intent(in): choice of upper boundary condition for soil hydrology
                              upperBoundHead,              & ! intent(in): upper boundary condition for matric head (m)
                              upperBoundTheta,             & ! intent(in): upper boundary condition for volumetric liquid water content (-)
                              scalarCanopyTranspiration,   & ! intent(in): canopy transpiration (kg m-2 s-1)
                              scalarGroundEvaporation,     & ! intent(in): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
                              mLayerMatricHeadIter,        & ! intent(in): matric head in each layer at the current iteration (m)
                              mLayerVolFracLiqIter,        & ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                              mLayerVolFracIceIter,        & ! intent(in): volumetric fraction of ice at the current iteration (-)
                              scalarAquiferStorageIter,    & ! intent(in): aquifer storage (m)
                              mLayerMatricIncrOld,         & ! intent(in): iteration increment for matric head of soil layers in the previous iteration (m)
                              mLayerLiquidIncrOld,         & ! intent(in): iteration increment for volumetric liquid water content of soil layers in the previous iteration (-) 

                              ! model decisions
                              ixDerivMethod,               & ! intent(in): choice of method used to compute derivative
                              ixRichards,                  & ! intent(in): choice of the form of Richards' equation
                              hc_Profile,                  & ! intent(in): index defining the option for the hydraulic conductivity profile
                              ixGroundwater,               & ! intent(in): choice of groundwater parameterization
                              ixBcLowerSoilHydrology,      & ! intent(in): choice of upper boundary condition for soil hydrology
                              ixSpatialGroundwater,        & ! intent(in): spatial representation of groundwater (local-column or single-basin)

                              ! model coordinate variables -- NOTE: use of ibeg and iend 
                              mLayerDepth,                 & ! intent(in): depth of the layer (m)
                              mLayerHeight,                & ! intent(in): height of the layer mid-point (m)
                              iLayerHeight,                & ! intent(in): height of the layer interfaces (m)

                              ! lower boundary conditions
                              lowerBoundHead,              & ! intent(in): lower boundary condition for matric head (m)
                              lowerBoundTheta,             & ! intent(in): lower boundary condition for volumetric liquid water content (-)

                              ! model forcing
                              scalarThroughfallRain,       & ! intent(in): computed throughfall rate (kg m-2 s-1)
                              scalarCanopyLiqDrainage,     & ! intent(in): computed drainage of liquid water (kg m-2 s-1)
                              scalarLiqFluxSnow,           & ! intent(in): liquid flux from the base of the snowpack (m s-1)

                              ! column inflow
                              mLayerColumnInflow,          & ! total inflow to the soil column (m3 s-1)

                              ! local attributes
                              HRUarea,                     & ! intent(in): HRU area (m2)
                              tan_slope,                   & ! intent(in): tan water table slope, taken as tan local ground surface slope (-)
                              contourLength,               & ! intent(in): length of contour at downslope edge of HRU (m)

                              ! general model parameters
                              wimplicit,                   & ! intent(in): weight assigned to start-of-step fluxes (-)

                              ! soil parameters
                              vGn_alpha,                   & ! intent(in): van Genutchen "alpha" parameter (m-1)
                              vGn_n,                       & ! intent(in): van Genutchen "n" parameter (-)
                              VGn_m,                       & ! intent(in): van Genutchen "m" parameter (-)
                              mpExp,                       & ! intent(in): empirical exponent in macropore flow equation (-)
                              theta_mp,                    & ! intent(in): volumetric liquid water content when macropore flow begins (-)
                              theta_sat,                   & ! intent(in): soil porosity (-)
                              theta_res,                   & ! intent(in): soil residual volumetric water content (-)
                              fieldCapacity,               & ! intent(in): field capacity (-)
                              kAnisotropic,                & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                              zScale_TOPMODEL,             & ! intent(in): TOPMODEL scaling factor (m)
                              qSurfScale,                  & ! intent(in): scaling factor in the surface runoff parameterization (-)
                              specificYield,               & ! intent(in): specific yield (-)
                              specificStorage,             & ! intent(in): specific storage coefficient (m-1)
                              aquiferScaleFactor,          & ! intent(in): scaling factor for aquifer storage in the big bucket (m)
                              aquiferBaseflowExp,          & ! intent(in): baseflow exponent for the big bucket (-)
                              f_impede,                    & ! intent(in): ice impedence factor (-)

                              ! model state variables -- NOTE: use of ibeg and iend
                              mLayerVolFracIce,            & ! intent(in): volumetric fraction of ice in each layer (-)
                              mLayerVolFracLiq,            & ! intent(in): volumetric fraction of liquid water in each layer (-)
                              mLayerMatricHead,            & ! intent(in): matric head in each layer (m)
                              scalarAquiferStorage,        & ! intent(in): (scalar)    ! aquifer storage at the start of the time step (m)
                              scalarSfcMeltPond,           & ! intent(in): (scalar)    ! ponded water caused by melt of the "snow without a layer" (kg m-2)

                              ! saturated hydraulic conductivity in each layer
                              mLayerSatHydCondMP,          & ! intent(in): saturated hydraulic conductivity of macropores at the mid-point of each layer (m s-1)
                              mLayerSatHydCond,            & ! intent(in): saturated hydraulic conductivity at the mid-point of each layer (m s-1)
                              iLayerSatHydCond,            & ! intent(in): saturated hydraulic conductivity at the interface of each layer (m s-1)

                              ! factors limiting transpiration (from vegFlux routine)
                              mLayerRootDensity,           & ! intent(in): root density in each layer (-)
                              scalarAquiferRootFrac,       & ! intent(in): fraction of roots below the lowest soil layer (-)
                              scalarTranspireLim,          & ! intent(in): weighted average of the transpiration limiting factor (-)
                              mLayerTranspireLim,          & ! intent(in): transpiration limiting factor in each layer (-)
                              scalarTranspireLimAqfr,      & ! intent(in): transpiration limiting factor for the aquifer (-)

                              ! initial fluxes (intent inout)
                              iLayerInitLiqFluxSoil,       & ! intent(inout): liquid flux at layer interfaces at the start of the time step (m s-1)
                              mLayerInitQMacropore,        & ! intent(inout): liquid water flux to/from the macropore at the start-of-step (m s-1)
                              mLayerInitTranspire,         & ! intent(inout): transpiration loss from each soil layer at the start of the time step (m s-1)
                              mLayerInitBaseflow,          & ! intent(inout): baseflow from each soil layer at the start-of-step (m s-1)
                              scalarInitAquiferTranspire,  & ! intent(inout): transpiration loss from the aquifer at the start-of-step (m s-1)
                              scalarInitAquiferRecharge,   & ! intent(inout): recharge to the aquifer at the start-of-step (m s-1)
                              scalarInitAquiferBaseflow,   & ! intent(inout): baseflow from the aquifer at the start-of-step (m s-1)

                              ! diagnostic scalar variables (intent out)
                              scalarRainPlusMelt,          & ! intent(out): rain plus melt (m s-1)
                              scalarSurfaceRunoff,         & ! intent(out): surface runoff (m s-1)
                              scalarAquiferTranspire,      & ! intent(out): transpiration loss from the aquifer (m s-1)
                              scalarAquiferRecharge,       & ! intent(out): recharge to the aquifer at the end-of-step (m s-1)
                              scalarAquiferBaseflow,       & ! intent(out): baseflow from the aquifer at the end-of-step (m s-1)

                              ! model diagnostic variables (intent out)
                              mLayerdTheta_dPsi,           & ! intent(out): derivative in the soil water characteristic w.r.t. psi (m-1)
                              mLayerdPsi_dTheta,           & ! intent(out): derivative in the soil water characteristic w.r.t. theta (m)
                              mLayerColumnOutflow,         & ! intent(out): outflow from each layer of the soil column (m3 s-1)
                              iLayerLiqFluxSoil,           & ! intent(out): liquid flux at layer interfaces at the end of the time step (m s-1)
                              mLayerQMacropore,            & ! intent(out): liquid water flux to/from the macropore (m s-1)
                              mLayerTranspire,             & ! intent(out): transpiration loss from each soil layer (m s-1)
                              mLayerBaseflow,              & ! intent(inout): baseflow from each soil layer -- only compute at the start of the step (m s-1)

                              ! output variables from the soilHydrol routine (intent out)
                              mLayerMatricHeadNew,         & ! intent(out): matric head in each layer at the next iteration (m)
                              mLayerVolFracLiqNew,         & ! intent(out): volumetric fraction of liquid water at the next iteration (-)
                              scalarAquiferStorageNew,     & ! intent(out): aquifer storage (m)
                              scalarSurfaceInfiltration,   & ! intent(out): surface infiltration rate (m s-1)
                              err,message)                   ! intent(out): error control
 ! utility modules
 USE soil_utils_module,only:volFracLiq      ! compute volumetric fraction of liquid water
 USE soil_utils_module,only:matricHead      ! compute matric head (m)
 USE soil_utils_module,only:dTheta_dPsi     ! compute derivative of the soil moisture characteristic w.r.t. psi (m-1)
 USE soil_utils_module,only:dPsi_dTheta     ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
 USE groundwatr_module,only:q_baseflow      ! compute baseflow and its derivative
 USE tridagSolv_module,only:tridag          ! solve tridiagonal system of equations
 implicit none
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! ***** input variables
 ! ------------------------------------------------------------------------------------------------------------------------------------------------- 
 ! input variables from the soilHydrol routine
 real(dp),intent(in)              :: dt                           ! time step (seconds)
 integer(i4b),intent(in)          :: iter                         ! iteration index
 integer(i4b),intent(in)          :: ixBcUpperSoilHydrology       ! choice of upper boundary condition for soil hydrology
 real(dp),intent(in)              :: upperBoundHead               ! upper boundary condition for matric head (m)
 real(dp),intent(in)              :: upperBoundTheta              ! upper boundary condition for volumetric liquid water content (-)
 real(dp),intent(in)              :: scalarCanopyTranspiration    ! canopy transpiration (kg m-2 s-1)
 real(dp),intent(in)              :: scalarGroundEvaporation      ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
 real(dp),intent(in)              :: mLayerMatricHeadIter(:)      ! matric head in each layer at the current iteration (m)
 real(dp),intent(in)              :: mLayerVolFracLiqIter(:)      ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),intent(in)              :: mLayerVolFracIceIter(:)      ! volumetric fraction of ice at the current iteration (-)
 real(dp),intent(in)              :: scalarAquiferStorageIter     ! aquifer storage (m)
 real(dp),intent(in)              :: mLayerMatricIncrOld(:)   ! iteration increment for matric head of soil layers in the previous iteration (m)
 real(dp),intent(in)              :: mLayerLiquidIncrOld(:)   ! iteration increment for volumetric liquid water content of soil layers in the previous iteration (-)
 ! model decisions
 integer(i4b),intent(in)          :: ixDerivMethod                ! choice of method used to compute derivative
 integer(i4b),intent(in)          :: ixRichards                   ! choice of the form of Richards' equation
 integer(i4b),intent(in)          :: hc_profile                   ! choice of type of hydraulic conductivity profile
 integer(i4b),intent(in)          :: ixGroundwater                ! choice of groundwater parameterization
 integer(i4b),intent(in)          :: ixBcLowerSoilHydrology       ! choice of upper boundary condition for soil hydrology
 integer(i4b),intent(in)          :: ixSpatialGroundwater         ! choice of the spatial representation of groundwater (local-column or single-basin)
 ! model coordinate variables
 real(dp),intent(in)              :: mLayerDepth(:)               ! depth of the layer (m)
 real(dp),intent(in)              :: mLayerHeight(:)              ! height of the layer mid-point (m)
 real(dp),intent(in)              :: iLayerHeight(0:)             ! height of the layer interfaces (m)
 ! diriclet lower boundary conditions
 real(dp),intent(in)              :: lowerBoundHead               ! lower boundary condition for matric head (m)
 real(dp),intent(in)              :: lowerBoundTheta              ! lower boundary condition for volumetric liquid water content (-)
 ! model forcing
 real(dp),intent(in)              :: scalarThroughfallRain        ! computed throughfall rate (kg m-2 s-1)
 real(dp),intent(in)              :: scalarCanopyLiqDrainage      ! computed drainage of liquid water (kg m-2 s-1)
 real(dp),intent(in)              :: scalarLiqFluxSnow            ! liquid flux from the base of the snowpack (m s-1)
 ! column inflow
 real(dp),intent(in)              :: mLayerColumnInflow(:)        ! total inflow to each layer in the soil column (m3 s-1)
 ! local attributes
 real(dp),intent(in)              :: HRUarea                      ! HRU area (m2)
 real(dp),intent(in)              :: tan_slope                    ! tan water table slope, taken as tan local ground surface slope (-)
 real(dp),intent(in)              :: contourLength                ! length of contour at downslope edge of HRU (m)
 ! general model parameters
 real(dp),intent(in)              :: wimplicit                    ! weight assigned to start-of-step fluxes (-)
 ! soil parameters
 real(dp),intent(in)              :: vGn_alpha                    ! van Genutchen "alpha" parameter (m-1)
 real(dp),intent(in)              :: vGn_n                        ! van Genutchen "n" parameter (-)
 real(dp),intent(in)              :: vGn_m                        ! van Genutchen "m" parameter (-)
 real(dp),intent(in)              :: mpExp                        ! empirical exponent in macropore flow equation (-)
 real(dp),intent(in)              :: theta_mp                     ! volumetric liquid water content when macropore flow begins (-)
 real(dp),intent(in)              :: theta_sat                    ! soil porosity (-)
 real(dp),intent(in)              :: theta_res                    ! soil residual volumetric water content (-)
 real(dp),intent(in)              :: fieldCapacity                ! field capacity (-)
 real(dp),intent(in)              :: kAnisotropic                 ! anisotropy factor for lateral hydraulic conductivity (-)
 real(dp),intent(in)              :: zScale_TOPMODEL              ! TOPMODEL scaling factor (m)
 real(dp),intent(in)              :: qSurfScale                   ! scaling factor in the surface runoff parameterization (-)
 real(dp),intent(in)              :: specificYield                ! specific yield (-)
 real(dp),intent(in)              :: specificStorage              ! specific storage coefficient (m-1)
 real(dp),intent(in)              :: aquiferScaleFactor           ! scaling factor for aquifer storage in the big bucket (m)
 real(dp),intent(in)              :: aquiferBaseflowExp           ! baseflow exponent for the big bucket (-)
 real(dp),intent(in)              :: f_impede                     ! ice impedence factor (-)
 ! state variables
 real(dp),intent(in)              :: mLayerVolFracIce(:)          ! volumetric fraction of ice at the start of the time step (-)
 real(dp),intent(in)              :: mLayerVolFracLiq(:)          ! volumetric fraction of liquid water at the start of the time step (-)
 real(dp),intent(in)              :: mLayerMatricHead(:)          ! matric head in each layer at the start of the time step (m)
 real(dp),intent(in)              :: scalarAquiferStorage         ! aquifer storage at the start of the time step (m)
 real(dp),intent(in)              :: scalarSfcMeltPond            ! ponded water caused by melt of the "snow without a layer" (kg m-2)
 ! saturated hydraulic conductivity
 real(dp),intent(in)              :: mLayerSatHydCondMP(:)        ! saturated hydraulic conductivity of macropores at the mid-point of each layer (m s-1)
 real(dp),intent(in)              :: mLayerSatHydCond(:)          ! saturated hydraulic conductivity at the mid-point of each layer (m s-1)
 real(dp),intent(in)              :: iLayerSatHydCond(0:)         ! saturated hydraulic conductivity at the interface of each layer (m s-1)
 ! factors limiting transpiration (from vegFlux routine)
 real(dp),intent(in)              :: mLayerRootDensity(:)         ! root density in each layer (-)
 real(dp),intent(in)              :: scalarAquiferRootFrac        ! fraction of roots below the lowest soil layer (-)
 real(dp),intent(in)              :: scalarTranspireLim           ! weighted average of the transpiration limiting factor (-)
 real(dp),intent(in)              :: mLayerTranspireLim(:)        ! transpiration limiting factor in each layer (-)
 real(dp),intent(in)              :: scalarTranspireLimAqfr       ! transpiration limiting factor for the aquifer (-)
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! ***** input/output variables -- start-of-step fluxes
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 real(dp),intent(inout)           :: iLayerInitLiqFluxSoil(0:)    ! liquid flux at soil layer interfaces at the start of the time step (m s-1)
 real(dp),intent(inout)           :: mLayerInitQMacropore(:)      ! liquid water flux to/from the macropore at the start-of-step (m s-1)
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
 real(dp),intent(out)             :: scalarAquiferTranspire       ! transpiration loss from the aquifer at the start-of-step (m s-1)
 real(dp),intent(out)             :: scalarAquiferRecharge        ! recharge to the aquifer at the end-of-step (m s-1)
 real(dp),intent(out)             :: scalarAquiferBaseflow        ! baseflow from the aquifer at the end-of-step (m s-1)
 ! diagnostic variables for each layer
 real(dp),intent(out)             :: mLayerdTheta_dPsi(:)         ! derivative in the soil water characteristic w.r.t. psi (m-1)
 real(dp),intent(out)             :: mLayerdPsi_dTheta(:)         ! derivative in the soil water characteristic w.r.t. theta (m)
 real(dp),intent(out)             :: mLayerColumnOutflow(:)       ! outflow from each layer of the soil column (m3 s-1)
 real(dp),intent(out)             :: iLayerLiqFluxSoil(0:)        ! liquid flux at soil layer interfaces at the end of the time step (m s-1)
 real(dp),intent(out)             :: mLayerQMacropore(:)          ! liquid water flux to/from the macropore (m s-1)
 real(dp),intent(out)             :: mLayerTranspire(:)           ! transpiration loss from each soil layer (m s-1)
 real(dp),intent(inout)           :: mLayerBaseflow(:)            ! baseflow from each soil layer -- only compute at the start of the step (m s-1)
 ! output variables from the soilHydrol routine
 real(dp),intent(out)             :: mLayerMatricHeadNew(:)       ! matric head in each layer at the next iteration (m)
 real(dp),intent(out)             :: mLayerVolFracLiqNew(:)       ! volumetric fraction of liquid water at the next iteration (m)
 real(dp),intent(out)             :: scalarAquiferStorageNew      ! aquifer storage (m)
 real(dp),intent(out)             :: scalarSurfaceInfiltration    ! surface infiltration rate (m s-1)
 integer(i4b),intent(out)         :: err                          ! error code
 character(*),intent(out)         :: message                      ! error message
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! ***** local variables
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables (general)
 character(LEN=256)               :: cmessage                 ! error message of downwind routine
 logical(lgt)                     :: printflag                ! flag to print crap to the screen
 integer(i4b)                     :: iLayer                   ! layer index
 logical(lgt),parameter           :: calcJacobian=.true.     ! flag to compute the Jacobian matrix
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
 real(dp),dimension(nLevels)      :: mLayerQMacroporeDeriv    ! derivative in flow to/from the macropore from each soil layer (m s-1)
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
 real(dp),dimension(nLevels)      :: compressibility          ! compressibility term (m-1)
 real(dp),dimension(nLevels)      :: mLayerMatricHeadDiff     ! iteration increment for matric head (m)
 real(dp),dimension(nLevels)      :: mLayerVolFracLiqDiff     ! iteration increment for volumetric fraction of liquid water (m)
 real(dp)                         :: scalarAquiferStorageDiff ! iteration increment for aquifer storage (m)
 ! impose solution constraints
 integer(i4b),dimension(1)        :: iLayerViolated           ! index of layer that was violated
 real(dp),dimension(1)            :: constraintViolation      ! maximum constraint violation
 real(dp)                         :: scaleIncrement           ! scaling factor for the iteration increment
 real(dp)                         :: notUsed_QMacroporeDeriv  ! derivative in macropore flow
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
  if(nSnow==0) scalarRainPlusMelt = (scalarThroughfallRain + scalarCanopyLiqDrainage)/iden_water + (scalarSfcMeltPond/dt)/iden_water  ! rainfall plus melt of the snow without a layer (convert to m s-1)
  if(nSnow>0)  scalarRainPlusMelt = scalarLiqFluxSnow                                              ! liquid water flux from the base of the snowpack (m s-1)
 endif
 if(ixBcUpperSoilHydrology==prescribedHead)then
  scalarRainPlusMelt = 0._dp
 endif
 !print*, 'scalarRainPlusMelt = ', scalarRainPlusMelt

 ! check
 !if(scalarRainPlusMelt*iden_water > 0.0001_dp)then
 ! write(*,'(a,3(f20.10,1x))') 'scalarThroughfallRain, scalarCanopyLiqDrainage, scalarSfcMeltPond/dt = ', scalarThroughfallRain, scalarCanopyLiqDrainage, scalarSfcMeltPond/dt
 !endif

 ! compute the fraction of transpiration loss from each snow-soil layer, and the aquifer
 if(scalarTranspireLim > epsilon(scalarTranspireLim))then ! (transpiration may be non-zero even if the soil moisture limiting factor is zero)
  mLayerTranspireFrac(:) = mLayerRootDensity(:)*mLayerTranspireLim(:)/scalarTranspireLim
  aquiferTranspireFrac   = scalarAquiferRootFrac*scalarTranspireLimAqfr/scalarTranspireLim
 else ! (possible for there to be non-zero conductance and therefore transpiration in this case)
  mLayerTranspireFrac(:) = mLayerRootDensity(:)
  aquiferTranspireFrac   = scalarAquiferRootFrac
 endif
 ! check that the sums are OK
 if(abs(1._dp - (sum(mLayerTranspireFrac)+aquiferTranspireFrac)) > 1.e-8_dp)then
  print*, 'sum(mLayerTranspireFrac) = ', sum(mLayerTranspireFrac)
  print*, 'aquiferTranspireFrac = ', aquiferTranspireFrac
  message=trim(message)//'problem allocating transpiration flux to soil layers and the aquifer'
  err=20; return
 endif

 ! compute transpiration loss from each soil layer, and the aquifer (kg m-2 s-1 --> m s-1)
 mLayerTranspire        = mLayerTranspireFrac(:)*scalarCanopyTranspiration/iden_water
 scalarAquiferTranspire = aquiferTranspireFrac*scalarCanopyTranspiration/iden_water

 ! initialize fluxes at the start of the time step
 if(iter == 1)then
  mLayerInitTranspire         = mLayerTranspire
  scalarInitAquiferTranspire  = scalarAquiferTranspire
 endif

 ! * get the number of state variables
 select case(ixSpatialGroundwater)
  ! separate groundwater representation in each local soil column
  case(localColumn)
   select case(ixGroundwater)
    case(bigBucket);                nState = nLevels+1
    case(noExplicit,qbaseTopmodel); nState = nLevels
    case default; err=20; message=trim(message)//'unknown groundwater parameterization'; return
   end select
  ! single groundwater store over the entire basin
  case(singleBasin)
   nState = nLevels
   if(ixGroundwater /= noExplicit)then
    message=trim(message)//'expect noExplicit representation of groundwater at the local scale if using a single-basin representation of groundwater'
    err=20; return
   endif
  case default; err=20; message=trim(message)//'unable to identify spatial representation of groundwater'; return
 end select ! (spatial representation of groundwater)

 ! allocate space for the tri-diagonal matrix
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
                    ! output: fluxes for the extended state vector (soil layers plus aquifer)
                    iLayerInitLiqFluxSoil,         & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                    mLayerInitQMacropore,          & ! intent(out): liquid water flux to/from the macropores (m s-1)
                    mLayerInitBaseflow,            & ! intent(out): baseflow from each soil layer (m s-1)
                    ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                    dq_dStateAbove,                & ! intent(out): derivatives in the flux w.r.t. state variable in the layer above (m s-1 or s-1)
                    dq_dStateBelow,                & ! intent(out): derivatives in the flux w.r.t. state variable in the layer below (m s-1 or s-1)
                    mLayerQMacroporeDeriv,         & ! intent(out): derivative in liquid water flux to/from the macropores (m s-1 [moisture form] s-1 [mixed form])
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
                    ! output: fluxes for the extended state vector (soil layers plus aquifer)
                    iLayerLiqFluxSoil,             & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                    mLayerQMacropore,              & ! intent(out): liquid water flux to/from the macropores (m s-1)
                    mLayerBaseflow,                & ! intent(inout): baseflow from each soil layer -- only compute at the start of the step (m s-1)
                    ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric liquid water -- in the layer above and layer below (m s-1 or s-1)
                    dq_dStateAbove,                & ! intent(out): derivatives in the flux w.r.t. state variable in the layer above (m s-1 or s-1)
                    dq_dStateBelow,                & ! intent(out): derivatives in the flux w.r.t. state variable in the layer below (m s-1 or s-1)
                    mLayerQMacroporeDeriv,         & ! intent(out): derivative in liquid water flux to/from the macropores (m s-1 [moisture form] s-1 [mixed form])
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
     mLayerInitQMacropore      = mLayerQMacropore
     iLayerInitLiqFluxSoil     = iLayerLiqFluxSoil
     scalarInitAquiferBaseflow = scalarAquiferBaseflow
     scalarInitAquiferRecharge = scalarAquiferRecharge
    endif  ! (if computing initial fluxes)
   endif  ! (if the first iteration)

  endif ! (if computing standard fluxes)

 end do  ! looping through initial vectors

 ! save the surface infiltration rate (m s-1)
 scalarSurfaceInfiltration = iLayerLiqFluxSoil(0)

 ! compute the compressibility term (m-1) 
 compressibility(1:nLevels) = 0._dp !-specificStorage*(mLayerVolFracLiqIter/theta_sat)   ! compressibility term

 !write(*,'(a,1x,i4,1x,e20.10)') 'ixDerivMethod, scalarAquiferRecharge = ', ixDerivMethod, scalarAquiferRecharge
 !write(*,'(a,1x,i4,1x,e20.10)') 'ixDerivMethod, scalarAquiferBaseflow = ', ixDerivMethod, scalarAquiferBaseflow
 !write(*,'(a,1x,i4,1x,e20.10)') 'ixDerivMethod, scalarAquiferRechargeDeriv = ', ixDerivMethod, scalarAquiferRechargeDeriv
 !write(*,'(a,1x,i4,1x,e20.10)') 'ixDerivMethod, scalarAquiferBaseflowDeriv = ', ixDerivMethod, scalarAquiferBaseflowDeriv
 !print*, 'specific Yield = ', specificYield

 ! *************************************************************************************************************************
 ! *************************************************************************************************************************
 ! *************************************************************************************************************************

 ! *****
 ! compute the residual vector
 call liqResidual(&
                  ! input: control variables
                  dt,                         & ! intent(in): length of the time step (s)
                  wimplicit,                  & ! intent(in): weight assigned to the start-of-step
                  iter,                       & ! intent(in): iteration count
                  ! input: coordinate variables
                  mLayerDepth,                & ! intent(in): depth of each layer (m)
                  ! input: model parameters for the compressibility term
                  theta_sat,                  & ! intent(in): porosity (-)
                  specificStorage,            & ! intent(in): specific storage (m-1)
                  ! input: initial flux vectors (start of the time step)
                  iLayerInitLiqFluxSoil,      & ! intent(in): initial liquid flux at layer interfaces (m s-1)
                  mLayerInitQMacropore,       & ! intent(in): initial liquid water flux to the macropores (m s-1)
                  mLayerInitBaseflow,         & ! intent(in): initial baseflow from each layer (m s-1)
                  mLayerInitTranspire,        & ! intent(in): initial transpiration from each layer -- from energy routine (m s-1)
                  ! input: trial flux vectors
                  iLayerLiqFluxSoil,          & ! intent(in): trial liquid flux at layer interfaces (m s-1)
                  mLayerQMacropore,           & ! intent(in): trial liquid water flux to the macropores (m s-1)
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
 !print*, 'iLayerLiqFluxSoil(0:1) = ', iLayerLiqFluxSoil(0:1) 
 !print*, 'mLayerVolFracLiqResidual = ', mLayerVolFracLiqResidual

 ! *****
 ! compute the residual for the groundwater store
 if(nState == nLevels+1)then
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
 endif

 ! *****
 ! populate the tri-diagnonal matrices for the soil layers
 wtim = (1._dp - wimplicit)*dt  ! weighted time
 ! (get off-diagonal elements)
 d_m1(1:nLevels-1) = (wtim/mLayerDepth(2:nLevels  ))*(-dq_dStateAbove(1:nLevels-1) )
 d_p1(1:nLevels-1) = (wtim/mLayerDepth(1:nLevels-1))*( dq_dStateBelow(1:nLevels-1) )
 ! (get diagonal elements)
 select case(ixRichards)
  case(moisture); diag(1:nLevels) = (wtim/mLayerDepth(1:nLevels))*(-dq_dStateBelow(0:nLevels-1) + dq_dStateAbove(1:nLevels) + mLayerQMacroporeDeriv(1:nLevels) ) + 1._dp
  case(mixdform); diag(1:nLevels) = (wtim/mLayerDepth(1:nLevels))*(-dq_dStateBelow(0:nLevels-1) + dq_dStateAbove(1:nLevels) + mLayerQMacroporeDeriv(1:nLevels) ) + mLayerdTheta_dPsi(1:nLevels) + compressibility(1:nLevels)
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
  ! check that there is no aquifer recharge or baseflow
  ! NOTE: can still have aquifer transpiration, if there is regional groundwater
  if(abs(scalarAquiferBaseflow) > tiny(dt) .or. abs(scalarAquiferRecharge) > tiny(dt))then
   err=20; message=trim(message)//'expect there to be no aquifer recharge or baseflow when there is no aquifer'; return
  endif
 endif  ! (if there is an aquifer) 

 ! *****
 ! test the Jacobian (if desired)
 if(calcJacobian)then
  call cmpJacobian(&
                   ! input: model control
                   ixRichards,                   & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                   ixBcLowerSoilHydrology,       & ! intent(in): index defining the type of boundary conditions
                   ! input: diagnostic variables
                   mLayerBaseflow,               & ! intent(in): baseflow from each soil layer at the start of the step (m s-1)
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
 if(err/=0)then
  print*, 'd_m1 = ', d_m1
  print*, 'diag = ', diag
  print*, 'd_p1 = ', d_p1
  message=trim(message)//trim(cmessage)
  err=-20; return
 endif

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

 ! adjust iteration increments in cases where iterations are oscillating
 !if(iter > 5)then
 ! select case(ixRichards)
 !  case(moisture); if(any(mLayerLiquidIncrOld(1:nLevels)*mLayerVolFracLiqDiff(1:nLevels) < 0._dp))then; mLayerVolFracLiqDiff(1:nLevels) = 0.5_dp*mLayerVolFracLiqDiff(1:nLevels); pause; endif
 !  case(mixdform); if(any(mLayerMatricIncrOld(1:nLevels)*mLayerMatricHeadDiff(1:nLevels) < 0._dp))then; mLayerMatricHeadDiff(1:nLevels) = 0.5_dp*mLayerMatricHeadDiff(1:nLevels); pause; endif
 !  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
 ! endselect
 !endif

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

    ! (update macropore flow)
    ! NOTE: this can cause negative flow when matric head increment is negative - decrease time step and try again
    if(iLayer > 0)then
     if(-mLayerQMacroporeDeriv(iLayer)*mLayerVolFracLiqDiff(iLayer) > mLayerQMacropore(iLayer))then
      message=trim(message)//'negative macropore flow associated with linearization of macropore flow w.r.t. volumetric fraction of liquid water'
      err=-20; return
     endif  ! (if linearization provides negative macropore flow)
     mLayerQMacropore(iLayer) = mLayerQMacropore(iLayer) + mLayerQMacroporeDeriv(iLayer)*mLayerVolFracLiqDiff(iLayer)
    endif

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

    ! (update fluxes)
    if(iLayer==0)then;           iLayerLiqFluxSoil(iLayer) = iLayerLiqFluxSoil(iLayer) + dq_dStateBelow(iLayer)*mLayerMatricHeadDiff(iLayer+1)
    elseif(iLayer==nLevels)then; iLayerLiqFluxSoil(iLayer) = iLayerLiqFluxSoil(iLayer) + dq_dStateAbove(iLayer)*mLayerMatricHeadDiff(iLayer)
    else;                        iLayerLiqFluxSoil(iLayer) = iLayerLiqFluxSoil(iLayer) + dq_dStateAbove(iLayer)*mLayerMatricHeadDiff(iLayer) &
                                                                                       + dq_dStateBelow(iLayer)*mLayerMatricHeadDiff(iLayer+1)
    endif

    ! (update macropore flow)
    ! NOTE: this can cause negative flow when matric head increment is negative - decrease time step and try again
    !if(iLayer > 0)then
     !if(-mLayerQMacroporeDeriv(iLayer)*mLayerMatricHeadDiff(iLayer) > mLayerQMacropore(iLayer))then
     ! write(*,'(a,1x,10(e20.10,1x))') 'mLayerVolFracLiqIter(iLayer), mLayerVolFracIceIter(iLayer), mLayerVolFracLiqIter(iLayer)+mLayerVolFracIceIter(iLayer) = ', &
     !                                  mLayerVolFracLiqIter(iLayer), mLayerVolFracIceIter(iLayer), mLayerVolFracLiqIter(iLayer)+mLayerVolFracIceIter(iLayer)
     ! write(*,'(a,1x,10(e20.10,1x))') 'mLayerMatricHeadDiff(iLayer), mLayerQMacroporeDeriv(iLayer), mLayerQMacropore(iLayer) = ', &
     !                                  mLayerMatricHeadDiff(iLayer), mLayerQMacroporeDeriv(iLayer), mLayerQMacropore(iLayer)
     ! ! call qMacropore again to get estimate of macropore flow
     ! call qMacropore(&
     !                 ! input: model control
     !                 .false.,                       & ! intent(in): flag determining if the derivative is desired
     !                 ixRichards,                    & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
     !                 ! input: soil parameters
     !                 vGn_alpha,                     & ! intent(in): van Genutchen "alpha" parameter (m-1)
     !                 vGn_n,                         & ! intent(in): van Genutchen "n" parameter (-)
     !                 VGn_m,                         & ! intent(in): van Genutchen "m" parameter (-)
     !                 theta_sat,                     & ! intent(in): soil porosity (-)
     !                 theta_res,                     & ! intent(in): soil residual volumetric water content (-)
     !                 ! input: state variable and diagnostic variables
     !                 mLayerMatricHeadNew(iLayer),   & ! intent(in): matric head in each layer (m)
     !                 mLayerVolFracLiqNew(iLayer),   & ! intent(in): volumetric liquid water content (-)
     !                 mLayerVolFracIceIter(iLayer),  & ! intent(in): volumetric ice content in each layer (-)
     !                 ! input: fluxes at the layer interfaces
     !                 iLayerLiqFluxSoil(iLayer-1),   & ! intent(in): liquid flux at the top of the layer (m s-1)
     !                 iLayerLiqFluxSoil(iLayer),     & ! intent(in): liquid flux at the bottom of the layer (m s-1)
     !                 ! input: derivatives in liquid fluxes
     !                 dq_dStateBelow(iLayer-1),      & ! intent(in): derivatives in the flux at top of the layer w.r.t. state variable in the layer below (m s-1 or s-1)
     !                 dq_dStateAbove(iLayer),        & ! intent(in): derivatives in the flux at bottom of the layer w.r.t. state variable in the layer above (m s-1 or s-1)
     !                 ! input: derivative in the soil water characteristic w.r.t. psi
     !                 mLayerdTheta_dPsi(iLayer),     & ! intent(in): derivative in the soil water characteristic w.r.t. psi (m-1)
     !                 ! output: water flux to/from macropores, and derivative
     !                 mLayerQMacropore(iLayer),      & ! intent(out): water flux to/from macropores (m s-1)
     !                 notUsed_QMacroporeDeriv,       & ! intent(out): derivative in macropore flow (m s-1 [moisture form] s-1 [mixed form])
     !                 ! output: error control
     !                 err,cmessage)                   ! intent(out): error control
     ! if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
     ! print*, 'mLayerQMacropore(iLayer) = ', mLayerQMacropore(iLayer)
     ! !pause
     !else  ! (if linearization provides negative macropore flow)
     ! mLayerQMacropore(iLayer) = mLayerQMacropore(iLayer) + mLayerQMacroporeDeriv(iLayer)*mLayerMatricHeadDiff(iLayer)
     !endif
    !endif  ! (if layer is > 0)

    ! (update the state variables)
    if(iLayer > 0)then
     mLayerMatricHeadNew(iLayer) = mLayerMatricHeadIter(iLayer) + mLayerMatricHeadDiff(iLayer)
     mLayerVolFracLiqNew(iLayer) = volFracLiq(mLayerMatricHeadNew(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
     if(mLayerVolFracLiqNew(iLayer) < epsilon(mLayerVolFracLiqNew))then; err=20; message=trim(message)//'very dry soil'; return; endif
    endif

   ! ***** unknown form of Richards' equation
   case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
  endselect
 end do  ! (loop through layers)


 ! *****
 ! impose solution constraints -- simplified bi-section method
 if(ixRichards==mixdform)then
  constraintViolation = maxval(mLayerMatricHeadNew - mLayerHeight)
  if(constraintViolation(1) > 0._dp)then  ! at least one layer where matric head is greater than soil depth
   ! scale the iteration increment to the midpoint between previous iteration and layer depth for the layer where constraints were violated
   iLayerViolated       = maxloc(mLayerMatricHeadNew - mLayerHeight)
   scaleIncrement       = 0.5_dp*(mLayerHeight(iLayerViolated(1)) - mLayerMatricHeadIter(iLayerViolated(1)))/mLayerMatricHeadDiff(iLayerViolated(1))
   mLayerMatricHeadDiff = mLayerMatricHeadDiff*scaleIncrement
   ! update model states again
   do iLayer=1,nLevels
    mLayerMatricHeadNew(iLayer) = mLayerMatricHeadIter(iLayer) + mLayerMatricHeadDiff(iLayer)
    mLayerVolFracLiqNew(iLayer) = volFracLiq(mLayerMatricHeadNew(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   end do
  endif  ! if there was a constraint violation
 endif  ! if using the mixed form of Richards' equation

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
                         ! output: fluxes at layer interfaces
                         iLayerLiqFluxSoil,            & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                         mLayerQMacropore,             & ! intent(out): liquid flux to/from the macropores (m s-1)
                         mLayerBaseflow,               & ! intent(inout): baseflow from each soil layer -- only compute at the start of the step (m s-1)
                         ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                         dq_dStateAbove,               & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer above (m s-1)
                         dq_dStateBelow,               & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer below (m s-1)
                         mLayerQMacroporeDeriv,        & ! intent(out): derivatives in liquid flux to/from the macropores (m s-1 [moisture form] s-1 [mixed form])
                         ! output: fluxes and derivatives for the aquifer storage
                         scalarAquiferRecharge,        & ! intent(out): recharge to the aquifer (m s-1)
                         scalarAquiferBaseflow,        & ! intent(out): total baseflow from the aquifer (m s-1)
                         scalarAquiferRechargeDeriv,   & ! intent(out): derivative in recharge to the aquifer w.r.t. aquifer storage (s-1)
                         scalarAquiferBaseflowDeriv,   & ! intent(out): derivative in baseflow from the aquifer w.r.t. aquifer storage (s-1)
                         ! output: error control
                         err,message)                    ! intent(out): error control
  USE soil_utils_module,only:hydCond_psi    ! compute hydraulic conductivity as a function of matric head
  USE soil_utils_module,only:hydCond_liq    ! compute hydraulic conductivity as a function of volumetric liquid water content
  USE soil_utils_module,only:hydCondMP_liq  ! compute hydraulic conductivity of macropores as a function of volumetric liquid water content
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
  ! output: liquid fluxes at layer interfaces
  real(dp),intent(out)             :: iLayerLiqFluxSoil(0:)         ! liquid flux at soil layer interfaces (m s-1)
  real(dp),intent(out)             :: mLayerQMacropore(:)           ! liquid flux to/from the macropores (m s-1)
  real(dp),intent(inout)           :: mLayerBaseflow(:)             ! baseflow from each soil layer -- only compute at the start of the step (m s-1)
  ! output: derivatives in fluxes w.r.t. state variables in the layer above and layer below (m s-1)
  real(dp),intent(out)             :: dq_dStateAbove(0:)            ! change in the flux in layer interfaces w.r.t. state variables in the layer above
  real(dp),intent(out)             :: dq_dStateBelow(0:)            ! change in the flux in layer interfaces w.r.t. state variables in the layer below
  real(dp),intent(out)             :: mLayerQMacroporeDeriv(:)      ! derivative in liquid flux to/from the macropores (m s-1 [moisture form] s-1 [mixed form])
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
  real(dp)                         :: sumDepthAvgCond               ! sum of depth-weighted hydraulic conductivity (m2 s-1)
  real(dp)                         :: aquiferHydCond                ! effective hydraulic conductivity of the aquifer (m s-1)
  real(dp)                         :: subSurfaceStorage             ! total sub-surface storage (m)
  real(dp)                         :: maximumSoilWater              ! maximum water stored in the soil column (m)
  real(dp)                         :: maximumFlowRate               ! maximum flow rate for baseflow from the soil column (m3 s-1)
  integer(i4b)                     :: ixSaturation                  ! index of highest soil layer that is close to saturation
  real(dp)                         :: totalOutflow                  ! total outflow from the entire soil column (m3 s-1)
  real(dp),dimension(nLevels)      :: fracTotalOutflow              ! fraction of outflow apportioned to each layer (-)
  real(dp),parameter               :: fracMinStorage=0.01_dp        ! minimum storage when baseflow occurs (used to avoid numerical problems associated with very dry soil)
  real(dp)                         :: xMinStorage                   ! minimim storage when baseflow can occur
  ! local variables: scalar trial values
  real(dp)                         :: scalarVolFracLiqTrial         ! trial value of volumetric liquid water content (-)
  real(dp)                         :: scalarMatricHeadTrial         ! trial value of matric head (m)
  real(dp)                         :: scalarHydCondTrial            ! trial value of hydraulic conductivity (m s-1)
  real(dp)                         :: scalarHydCondMicro            ! trial value of hydraulic conductivity of micropores (m s-1)
  real(dp)                         :: scalarHydCondMacro            ! trial value of hydraulic conductivity of macropores (m s-1)
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
  ! special case of exceeding porosity
  integer(i4b)                     :: ixFlag_porosityCheck     ! flag to indicate if modified fluxes to handle solution constraints
  logical(lgt),parameter           :: doPorosityFix=.false.    ! flag to do porosity fix
  real(dp)                         :: xDamp                    ! damping parameter (-)
  real(dp)                         :: aPor                     ! available porosity (-)
  real(dp)                         :: qMax                     ! maximum possible flow into layer (m s-1)
  real(dp)                         :: xFlx                     ! total flux into layer (m s-1)
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
  ! compute diagnostic variables at the nodes throughout the soil profile
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  do iSoil=1,nLevels ! (loop through soil layers)
    call diagv_node(&
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
                    mpExp,                           & ! intent(in): empirical exponent in macropore flow equation (-)
                    theta_sat,                       & ! intent(in): soil porosity (-)
                    theta_res,                       & ! intent(in): soil residual volumetric water content (-)
                    theta_mp,                        & ! intent(in): volumetric liquid water content when macropore flow begins (-)
                    f_impede,                        & ! intent(in): ice impedence factor (-)
                    ! input: saturated hydraulic conductivity
                    mLayerSatHydCond(iSoil),         & ! intent(in): saturated hydraulic conductivity at the mid-point of each layer (m s-1)
                    mLayerSatHydCondMP(iSoil),       & ! intent(in): saturated hydraulic conductivity of macropores at the mid-point of each layer (m s-1)
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

  ! set derivative w.r.t. state above to missing (does not exist)
  dq_dStateAbove(0) = valueMissing

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
                   desireAnal,                         & ! intent(in): flag indicating if derivatives are desired
                   ixRichards,                         & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                   ixBcUpperSoilHydrology,             & ! intent(in): index defining the type of boundary conditions (neumann or diriclet)
                   ! input: state variables
                   scalarMatricHeadTrial,              & ! intent(in): matric head in he upper-most soil layer (m)
                   scalarVolFracLiqTrial,              & ! intent(in): volumetric liquid water content the upper-most soil layer (-)
                   mLayerVolFracIceTrial(1),           & ! intent(in): volumetric ice content in the upper-most soil layer (-)
                   ! input: depth of upper-most soil layer (m)
                   mLayerDepth(1),                     & ! intent(in): depth of upper-most soil layer (m)
                   ! input: boundary conditions
                   upperBoundHead,                     & ! intent(in): upper boundary condition (m)
                   upperBoundTheta,                    & ! intent(in): upper boundary condition (-)
                   ! input: flux at the upper boundary
                   scalarRainPlusMelt,                 & ! intent(in): rain plus melt (m s-1)
                   scalarGroundEvaporation/iden_water, & ! intent(in): ground evaporation (m s-1)
                   ! input: derivative in soil water characteristix
                   mLayerdPsi_dTheta(1),               & ! intent(in): derivative of the soil moisture characteristic w.r.t. theta (m)
                   mLayerdTheta_dPsi(1),               & ! intent(in): derivative of the soil moisture characteristic w.r.t. psi (m-1)
                   ! input: transmittance
                   iLayerSatHydCond(0),                & ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
                   iceImpedeFac(1),                    & ! intent(in): ice impedence factor in the upper-most soil layer (-)
                   ! input: soil parameters
                   vGn_alpha,                          & ! intent(in): van Genutchen "alpha" parameter (m-1)
                   vGn_n,                              & ! intent(in): van Genutchen "n" parameter (-)
                   VGn_m,                              & ! intent(in): van Genutchen "m" parameter (-)
                   theta_sat,                          & ! intent(in): soil porosity (-)
                   theta_res,                          & ! intent(in): soil residual volumetric water content (-)
                   qSurfScale,                         & ! intent(in): scaling factor in the surface runoff parameterization (-)
                   ! output: hydraulic conductivity and diffusivity at the surface
                   iLayerHydCond(0),                   & ! intent(out): hydraulic conductivity at the surface (m s-1)
                   iLayerDiffuse(0),                   & ! intent(out): hydraulic diffusivity at the surface (m2 s-1)
                   ! output: fluxes at layer interfaces and surface runoff
                   scalarSurfaceRunoff,                & ! intent(out): surface runoff (m s-1)
                   iLayerLiqFluxSoil(0),               & ! intent(out): surface infiltration (m s-1)
                   ! output: deriavtives in surface infiltration w.r.t. volumetric liquid water (m s-1) and matric head (s-1) in the upper-most soil layer
                   dq_dStateBelow(0),                  & ! intent(out): derivative in surface infiltration w.r.t. state variable in the upper-most soil layer (m s-1 or s-1)
                   ! output: error control
                   err,cmessage)                         ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   !print*, 'scalarRainPlusMelt, scalarSurfaceRunoff, iLayerLiqFluxSoil(0), iLayerHydCond(0) = ', scalarRainPlusMelt, scalarSurfaceRunoff, iLayerLiqFluxSoil(0), iLayerHydCond(0)

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
  !pause

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
    vectorHydCondTrial(1:2)   = mLayerHydCond(iLayer:iLayer+1)
    vectorDiffuseTrial(1:2) = mLayerDiffuse(iLayer:iLayer+1)
    ! make appropriate perturbations
    if(ixPerturb > 0)then
     select case(ixRichards)
      case(moisture)
       scalardPsi_dTheta             = dPsi_dTheta(vectorVolFracLiqTrial(ixPerturb),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
       vectorHydCondTrial(ixPerturb) = hydCond_liq(vectorVolFracLiqTrial(ixPerturb),mLayerSatHydCond(ixOriginal),theta_res,theta_sat,vGn_m) * iceImpedeFac(ixOriginal)
       vectorDiffuseTrial(ixPerturb) = scalardPsi_dTheta * vectorHydCondTrial(ixPerturb)
      case(mixdform)
       scalarVolFracLiqTrial = volFracLiq(vectorMatricHeadTrial(ixPerturb),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
       scalarHydCondMicro    = hydCond_psi(vectorMatricHeadTrial(ixPerturb),mLayerSatHydCond(ixOriginal),vGn_alpha,vGn_n,vGn_m) * iceImpedeFac(ixOriginal)
       scalarHydCondMacro    = hydCondMP_liq(scalarVolFracLiqTrial,theta_sat,theta_mp,mpExp,mLayerSatHydCondMP(ixOriginal),mLayerSatHydCond(ixOriginal))      
       vectorHydCondTrial(ixPerturb) = scalarHydCondMicro + scalarHydCondMacro
      case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
     end select ! (form of Richards' equation)
    endif

    ! =====
    ! compute vertical flux at layer interface and its derivative w.r.t. the state above and the state below...
    ! =========================================================================================================
    call iLayerFlux(&
                    ! input: model control
                    desireAnal,                         & ! intent(in): flag indicating if derivatives are desired
                    ixRichards,                         & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                    ! input: state variables (adjacent layers)
                    vectorMatricHeadTrial,              & ! intent(in): matric head at the soil nodes (m)
                    vectorVolFracLiqTrial,              & ! intent(in): volumetric liquid water content at the soil nodes (-)
                    ! input: model coordinate variables (adjacent layers)
                    mLayerHeight(iLayer:iLayer+1),      & ! intent(in): height of the soil nodes (m)
                    ! input: transmittance (adjacent layers)
                    vectorHydCondTrial,                 & ! intent(in): hydraulic conductivity at the soil nodes (m s-1)
                    vectorDiffuseTrial,                 & ! intent(in): hydraulic diffusivity at the soil nodes (m2 s-1)
                    ! input: transmittance derivatives (adjacent layers)
                    dHydCond_dVolLiq(iLayer:iLayer+1),  & ! intent(in): change in hydraulic conductivity w.r.t. change in volumetric liquid water content (m s-1)
                    dDiffuse_dVolLiq(iLayer:iLayer+1),  & ! intent(in): change in hydraulic diffusivity w.r.t. change in volumetric liquid water content (m2 s-1)
                    dHydCond_dMatric(iLayer:iLayer+1),  & ! intent(in): change in hydraulic conductivity w.r.t. change in matric head (s-1)
                    ! output: tranmsmittance at the layer interface (scalars)
                    iLayerHydCond(iLayer),              & ! intent(out): hydraulic conductivity at the interface between layers (m s-1)
                    iLayerDiffuse(iLayer),              & ! intent(out): hydraulic diffusivity at the interface between layers (m2 s-1)
                    ! output: vertical flux at the layer interface (scalars)
                    iLayerLiqFluxSoil(iLayer),          & ! intent(out): vertical flux of liquid water at the layer interface (m s-1)
                    ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                    dq_dStateAbove(iLayer),             & ! intent(out): derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer above (m s-1 or s-1)
                    dq_dStateBelow(iLayer),             & ! intent(out): derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer below (m s-1 or s-1) 
                    ! output: error control
                    err,cmessage)                         ! intent(out): error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    !if(iLayer < 5)then
    ! write(*,'(a,i4,1x,10(e20.10,1x))') 'after iLayerFlux: iLayer, vectorMatricHeadTrial(1), iLayerLiqFluxSoil(iLayer), dq_dStateAbove(iLayer)*wtim/mLayerDepth(iLayer+1), dq_dStateBelow(iLayer)*wtim/mLayerDepth(iLayer) = ', &
    !                                                       iLayer, vectorMatricHeadTrial(1), iLayerLiqFluxSoil(iLayer), dq_dStateAbove(iLayer)*wtim/mLayerDepth(iLayer+1), dq_dStateBelow(iLayer)*wtim/mLayerDepth(iLayer)
    !endif

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
   if(iLayer==1) write(*,'(a)') 'ixDerivMethod, iLayer, dq_dStateBelow(iLayer), dq_dStateAbove(iLayer), mLayerVolFracLiqTrial(iLayer+1) = '
   if(iLayer< 9) write(*,'(2(i4,1x),10(e15.7,1x))') ixDerivMethod, iLayer, dq_dStateBelow(iLayer), dq_dStateAbove(iLayer), mLayerVolFracLiqTrial(iLayer+1)

  end do  ! (looping through soil layers)

  ! *************************************************************************************************************************************************
  ! *************************************************************************************************************************************************

  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! * compute drainage flux from the bottom of the soil profile, and its derivative
  ! -------------------------------------------------------------------------------------------------------------------------------------------------

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

 
  ! *************************************************************************************************************************************************
  ! *************************************************************************************************************************************************
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! * compute the recharge to the aquifer and its derivative
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  select case(ixGroundwater)

   ! =====
   ! using the big bucket...
   ! =======================
   case(bigBucket)
    scalarAquiferRecharge      = iLayerLiqFluxSoil(nLevels)  ! recharge = drainage flux from the bottom of the soil profile (m s-1)
    scalarAquiferRechargeDeriv = 0._dp                       ! recharge does not depend on aquifer storage

   ! =====
   ! no explicit aquifer...
   ! ======================
   case(qbaseTopmodel,noExplicit)
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
  select case(ixGroundwater)

   ! ------------------------------------------------------------------------------------------------------------------------------------------------
   ! the big bucket
   case(bigBucket)

    ! compute effective hydraulic conductivity
    aquiferHydCond = kAnisotropic*iLayerSatHydCond(0)

    ! compute baseflow
    call q_baseflow(&
                    ! input: model decisions
                    deriv_desired,        & ! intent(in): flag indicating if derivatives are desired
                    ixDerivMethod,        & ! intent(in): choice of method used to compute derivative
                    ! input: effective parameters
                    aquiferHydCond,       & ! intent(in): effective hydraulic conductivity (m s-1)
                    aquiferScaleFactor,   & ! intent(in): scaling factor for aquifer storage (m)
                    aquiferBaseflowExp,   & ! intent(in): exponent in bucket baseflow parameterization (-)
                    ! input: state variables
                    scalarAquiferStorageTrial,   & ! intent(in): aquifer storage (m)
                    ! output
                    scalarAquiferBaseflow,       & ! intent(out): total baseflow (m s-1)
                    scalarAquiferBaseflowDeriv,  & ! intent(out): derivative in baseflow flux w.r.t. water table depth (m s-1)
                    err,cmessage)                  ! intent(out): error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

    ! baseflow from a separate aquifer, so no baseflow loss in each soil layer
    mLayerBaseflow(:) = 0._dp 

   ! ------------------------------------------------------------------------------------------------------------------------------------------------
   ! baseflow using the TOPMODEL-ish parameterization
   case(qbaseTopmodel)

    ! only compute baseflow at the start of the time step (avoid convergence problems)
    if(iter == 1)then

     ! initialize sub-surface storage (m)
     subSurfaceStorage = 0._dp
     ixSaturation      = nLevels+1

     ! compute the "effective" sub-surface storage
     do iLayer=nLevels,1,-1  ! start at the lowest soil layer and work upwards
      if(mLayerVolFracLiqTrial(iLayer) > fieldCapacity)then
       ixSaturation      = iLayer  ! index of saturated layer -- keeps getting over-written as move upwards
       subSurfaceStorage = subSurfaceStorage + (mLayerVolFracLiqTrial(iLayer) - fieldCapacity)*mLayerDepth(iLayer)
       !write(*,'(a,1x,3(f20.10,1x))') 'mLayerVolFracLiqTrial(iLayer), fieldCapacity, subSurfaceStorage = ', mLayerVolFracLiqTrial(iLayer), fieldCapacity, subSurfaceStorage
      else
       exit                          ! (only consider saturated layer at the bottom of the soil profile) 
      endif
     end do  ! (looping through soil layers)

     ! baseflow occurs
     if(ixSaturation <= nLevels)then

      !print*, 'subSurfaceStorage, ixSaturation = ', subSurfaceStorage, ixSaturation
      
      ! compute effective parameters
      maximumSoilWater = iLayerHeight(nSoil)*(theta_sat - fieldCapacity)       ! maximum aquifer storage (m)
      maximumFlowRate  = tan_slope*iLayerSatHydCond(0)*kAnisotropic*maximumSoilWater*contourLength &
                           / ((theta_sat - fieldCapacity)*aquiferBaseflowExp)  ! effective hydraulic conductivity (m3/s)
      !print*, 'tan_slope, iLayerSatHydCond(0), maximumSoilWater, contourLength, aquiferBaseflowExp = ', &
      !         tan_slope, iLayerSatHydCond(0), maximumSoilWater, contourLength, aquiferBaseflowExp
      !print*, 'maximumFlowRate/HRUarea = ', maximumFlowRate/HRUarea
 
      ! compute total outflow from the entire soil column (m3 s-1)
      xMinStorage    = fracMinStorage*maximumSoilWater
      totalOutflow   = maximumFlowRate*((subSurfaceStorage - xMinStorage)/(maximumSoilWater - xMinStorage))**aquiferBaseflowExp
      !print*, 'subSurfaceStorage, maximumSoilWater, maximumFlowRate, totalOutflow = ', subSurfaceStorage, maximumSoilWater, maximumFlowRate*3600._dp, totalOutflow*3600._dp
 
      ! compute the depth-weighted hydraulic conductivity (m2 s-1)
      sumDepthAvgCond = sum(mLayerHydCond(ixSaturation:nLevels)*mLayerDepth(ixSaturation:nLevels))
      if(sumDepthAvgCond < tiny(theta_sat))then; err=20; message=trim(message)//'zero depth-weighted conductivity when baseflow occurs'; return; endif
      !print*, 'sumDepthAvgCond, tiny(theta_sat) = ', sumDepthAvgCond, tiny(theta_sat)
 
      ! compute the fraction of outflow apportioned to each layer (-)
      fracTotalOutflow(ixSaturation:nLevels) = (mLayerHydCond(ixSaturation:nLevels)*mLayerDepth(ixSaturation:nLevels)) / sumDepthAvgCond
      if(ixSaturation > 1) fracTotalOutflow(1:ixSaturation-1) = 0._dp
      !print*, 'fracTotalOutflow, sum(fracTotalOutflow) = ', fracTotalOutflow, sum(fracTotalOutflow)
      if(abs(1._dp - sum(fracTotalOutflow)) > verySmall)then
       write(*,'(a,1x,10(e20.10,1x))') 'fracTotalOutflow, sum(fracTotalOutflow) = ', fracTotalOutflow, sum(fracTotalOutflow)
       print*, 'abs(1._dp - sum(fracTotalOutflow)) = ', abs(1._dp - sum(fracTotalOutflow))
       message=trim(message)//'fraction of baseflow does not sum to 1'
       err=20; return
      endif
    
      ! compute the outflow from each soil layer (m3 s-1)
      mLayerColumnOutflow(1:nLevels) = fracTotalOutflow(1:nLevels)*totalOutflow
      !print*, 'mLayerColumnOutflow(1:nLevels), fracTotalOutflow(1:nLevels), totalOutflow = ', mLayerColumnOutflow(1:nLevels), fracTotalOutflow(1:nLevels), totalOutflow
      !pause

     ! no layers above field capacity
     else
      mLayerColumnOutflow(1:nLevels) = 0._dp
     endif

     ! compute the net baseflow from each soil layer (m s-1)
     mLayerBaseflow(1:nLevels) = (mLayerColumnOutflow(1:nLevels) - mLayerColumnInflow(1:nLevels))/HRUarea
     !print*, 'mLayerBaseflow = ', mLayerBaseflow
 
    endif  ! (if the first iteration)
 
    ! no explicit aquifer (only use soil layers)
    scalarAquiferBaseflow      = 0._dp
    scalarAquiferBaseflowDeriv = 0._dp

   ! ------------------------------------------------------------------------------------------------------------------------------------------------
   ! no explicit baseflow
   case(noExplicit)

    ! no explicit aquifer (only use soil layers)
    scalarAquiferBaseflow      = 0._dp
    scalarAquiferBaseflowDeriv = 0._dp

    ! no baseflow at all
    mLayerBaseflow(:) = 0._dp

   ! ------------------------------------------------------------------------------------------------------------------------------------------------
   ! error check
   case default
    message=trim(message)//'unable to identify baseflow parameterization'
    err=20; return
 
  end select

  ! *************************************************************************************************************************
  ! *************************************************************************************************************************
  ! *************************************************************************************************************************
  ! *************************************************************************************************************************
  ! ***** HANDLE SPECIAL CASE WHERE LIQUID WATER FILLS PARTIALLY FROZEN LAYER BEYOND AVAILABLE POROSITY *********************
  ! *************************************************************************************************************************
  ! *************************************************************************************************************************
  ! *************************************************************************************************************************
  ! *************************************************************************************************************************

  !print*, 'mLayerSatHydCondMP = ', mLayerSatHydCondMP
  !print*, 'mLayerSatHydCond   = ', mLayerSatHydCond
  !print*, 'mpExp, theta_mp    = ', mpExp, theta_mp
  !pause

  ! loop through soil layers
  do iSoil=1,nLevels
   ! check if there is ice in the layer
   !if(mLayerVolFracIceIter(iSoil) > epsilon(dt) .or. ixRichards==moisture)then
   ! ! compute flow from micropores to macropores
   ! call qMacropore(&
   !                 ! input: model control
   !                 deriv_desired,                & ! intent(in): flag determining if the derivative is desired
   !                 ixRichards,                   & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
   !                 ! input: soil parameters
   !                 vGn_alpha,                    & ! intent(in): van Genutchen "alpha" parameter (m-1)
   !                 vGn_n,                        & ! intent(in): van Genutchen "n" parameter (-)
   !                 VGn_m,                        & ! intent(in): van Genutchen "m" parameter (-)
   !                 theta_sat,                    & ! intent(in): soil porosity (-)
   !                 theta_res,                    & ! intent(in): soil residual volumetric water content (-)
   !                 ! input: state variable and diagnostic variables
   !                 mLayerMatricHeadTrial(iSoil), & ! intent(in): matric head in each layer (m)
   !                 mLayerVolFracLiqTrial(iSoil), & ! intent(in): volumetric liquid water content (-)
   !                 mLayerVolFracIceTrial(iSoil), & ! intent(in): volumetric ice content in each layer (-)
   !                 ! input: liquid fluxes at the layer interfaces
   !                 iLayerLiqFluxSoil(iSoil-1),   & ! intent(in): liquid flux at the top of the layer (m s-1)
   !                 iLayerLiqFluxSoil(iSoil),     & ! intent(in): liquid flux at the bottom of the layer (m s-1)
   !                 ! input: derivatives in liquid fluxes
   !                 dq_dStateBelow(iSoil-1),      & ! intent(in): derivatives in the flux w.r.t. state variable in the layer below (m s-1 or s-1)
   !                 dq_dStateAbove(iSoil),        & ! intent(in): derivatives in the flux w.r.t. state variable in the layer above (m s-1 or s-1)
   !                 ! input: derivative in the soil water characteristic w.r.t. psi
   !                 mLayerdTheta_dPsi(iSoil),     & ! intent(in): derivative in the soil water characteristic w.r.t. psi (m-1)
   !                 ! output: water flux to/from macropores, and derivative
   !                 mLayerQMacropore(iSoil),      & ! intent(out): water flux to/from macropores (m s-1)
   !                 mLayerQMacroporeDeriv(iSoil), & ! intent(out): derivative in macropore flow (m s-1 [moisture form] s-1 [mixed form])
   !                 ! output: error control
   !                 err,cmessage)                   ! intent(out): error control
   ! if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   !! no ice and the mixed form of Richards' equation -- simulate saturated flow
   !else
    mLayerQMacropore(iSoil)      = 0._dp
    mLayerQMacroporeDeriv(iSoil) = 0._dp
    if(mLayerVolFracLiqTrial(iSoil) + mLayerVolFracIceTrial(iSoil) > theta_sat)then
     message=trim(message)//'liquid and ice exceeds porosity'
     err=20; return
    endif
   !endif
  end do  ! (looping through soil layers)
 

 ! ! check if desire the porosity fix
 ! if(doPorosityFix)then
 !
 !  ! *****
 !  ! handle special case where liquid flux exceeds available pore space
 !  ! -- truncate fluxes (based on the SHAW model)
 !  ! -- start at the bottom and work upwards
 !  ixFlag_porosityCheck = 0
 !  do iSoil=nLevels,2,-1
 !   ! * only implement fix if layer has ice
 !   if(mLayerVolFracIceIter(iSoil) > epsilon(dt))then
 !
 !    ! (compute maximum flux that can fill available pore space)
 !    qMax  = (theta_sat - mLayerVolFracLiq(iSoil) - mLayerVolFracIce(iSoil))*mLayerDepth(iSoil)/dt ! m s-1
 !    if(qMax < 0._dp) qMax = 0._dp
 !    ! (compute the total flux)
 !    xFlx = (iLayerLiqFluxSoil(iSoil-1) - iLayerLiqFluxSoil(iSoil)) - mLayerBaseflow(iSoil)     + mLayerTranspire(iSoil)
 !    ! (check if the flux exceeds available pore space)
 !    if(xFlx > qMax)then
 !     ixFlag_porosityCheck = 1

 !     ! part 1: downward flux at top of layer
 !     !  -- try and truncate top flux first
 !     if(iLayerLiqFluxSoil(iSoil-1) > 0._dp)then  ! (downward flux at top of layer)
 !      iLayerLiqFluxSoil(iSoil-1) = iLayerLiqFluxSoil(iSoil) + qMax ! maximum possible flux
 !      dq_dStateAbove(iSoil-1)    = 0._dp ! derivatives in the flux w.r.t. state variable in the layer above (m s-1 or s-1)
 !      dq_dStateBelow(iSoil-1)    = 0._dp ! derivatives in the flux w.r.t. state variable in the layer below (m s-1 or s-1)
 !      !print*, 'part1: iSoil, iLayerLiqFluxSoil(iSoil-1), qMax = ', iSoil, iLayerLiqFluxSoil(iSoil-1), qMax
 !      ! (handle case where upward flux at layer bottom > qmax)
 !      if(iLayerLiqFluxSoil(iSoil-1) < 0._dp)then
 !       iLayerLiqFluxSoil(iSoil-1) = 0._dp
 !       iLayerLiqFluxSoil(iSoil)   = -qMax
 !       if(iSoil == nLevels)then; err=20; message=trim(message)//'negative flow at lower boundary (part1) -- check if this is possible'; return; endif
 !      endif  ! if upward flux at layer bottom > qmax

 !     ! part 2: upward flux at top of layer
 !     !  -- we can only exceed porosity here if upward flux at bottom exceeds upward flux at top
 !     !  -- reduce upward flux at bottom
 !     else  ! (upward flux at top of layer)
 !      iLayerLiqFluxSoil(iSoil) = iLayerLiqFluxSoil(iSoil-1) - qMax
 !      dq_dStateAbove(iSoil)    = 0._dp ! derivatives in the flux w.r.t. state variable in the layer above (m s-1 or s-1)
 !      dq_dStateBelow(iSoil)    = 0._dp ! derivatives in the flux w.r.t. state variable in the layer below (m s-1 or s-1)
 !      !print*, 'part2: iSoil, iLayerLiqFluxSoil(iSoil), qMax = ', iSoil, iLayerLiqFluxSoil(iSoil), qMax
 !      if(iLayerLiqFluxSoil(iSoil) > 0._dp)then; err=20; message=trim(message)//'unexpected positive flux while truncating'; return; endif
 !      if(iSoil == nLevels)then; err=20; message=trim(message)//'negative flow at lower boundary (part 2) -- check if this is possible'; return; endif
 !     endif  ! (switch between upward and downward flux at top of layer)

 !    endif  ! (if exceed porosity)
 !   endif  ! (if layer has ice)
 !  end do ! (looping through soil layers)
 

 !  ! *****
 !  ! handle special case of the top soil layer (don't modify the upper flux)
 !  ! (compute maximum flux that can fill available pore space)
 !  if(mLayerVolFracIceIter(1) > epsilon(dt))then
 !   aPor  = theta_sat - mLayerVolFracLiq(1) - mLayerVolFracIce(1)  ! available porosity (-)
 !   qMax  = aPor*mLayerDepth(1)/dt ! m s-1
 !   if(qMax < 0._dp) qMax = 0._dp
 !   !print*, 'aPor, qMax, iLayerLiqFluxSoil(0), iLayerLiqFluxSoil(1) = ', aPor, qMax, iLayerLiqFluxSoil(0), iLayerLiqFluxSoil(1)
 !   ! (compute the total flux)
 !   xFlx = (iLayerLiqFluxSoil(0) - iLayerLiqFluxSoil(1)) - mLayerBaseflow(1) + mLayerTranspire(1)
 !   ! (check if the flux exceeds available pore space)
 !   if(xFlx > qMax)then
 !    ixFlag_porosityCheck = 1
 !    ! NOTE: always downward flux at top of profile
 !    !print*, 'before update: iLayerLiqFluxSoil(1), qMax = ', iLayerLiqFluxSoil(1), qMax
 !    if(-iLayerLiqFluxSoil(1) > qMax)then
 !     iLayerLiqFluxSoil(1) = -qMax ! maximum possible flux
 !     dq_dStateAbove(1)    = 0._dp ! derivatives in the flux w.r.t. state variable in the layer above (m s-1 or s-1)
 !     dq_dStateBelow(1)    = 0._dp ! derivatives in the flux w.r.t. state variable in the layer below (m s-1 or s-1)
 !    endif   ! (if upward flux at the bottom of the top layer exceeds maximum possible flux)
 !    !print*, 'after update: iLayerLiqFluxSoil(1) = ', iLayerLiqFluxSoil(1)
 !   endif   ! (if the net flux exceeds maximum possible flux)
 !  endif   ! (if layer has ice)



 !  ! now move from the top to bottom to make everything honky-dory
 !  if(ixFlag_porosityCheck > 0)then
 !   do iSoil=2,nLevels
 !    ! * only implement fix if layer has ice
 !    if(mLayerVolFracIceIter(iSoil) > epsilon(dt))then
  
 !     ! (compute maximum flux that can fill available pore space)
 !     qMax  = (theta_sat - mLayerVolFracLiq(iSoil) - mLayerVolFracIce(iSoil))*mLayerDepth(iSoil)/dt ! m s-1
 !     if(qMax < 0._dp) qMax = 0._dp
 !     ! (compute the total flux)
 !     xFlx = (iLayerLiqFluxSoil(iSoil-1) - iLayerLiqFluxSoil(iSoil)) - mLayerBaseflow(iSoil)     + mLayerTranspire(iSoil)
 !     ! (check if the flux exceeds available pore space)
 !     if(xFlx > qMax)then

 !      ! part 1: upward flux at bottom of layer
 !      !  -- try and truncate bottom flux first
 !      if(iLayerLiqFluxSoil(iSoil) < 0._dp)then
 !       iLayerLiqFluxSoil(iSoil) = iLayerLiqFluxSoil(iSoil-1) - qMax ! ! maximum possible flux
 !       dq_dStateAbove(iSoil)    = 0._dp ! derivatives in the flux w.r.t. state variable in the layer above (m s-1 or s-1)
 !       dq_dStateBelow(iSoil)    = 0._dp ! derivatives in the flux w.r.t. state variable in the layer below (m s-1 or s-1)
 !       !print*, 'part1: iSoil, iLayerLiqFluxSoil(iSoil), qMax = ', iSoil, iLayerLiqFluxSoil(iSoil), qMax
 !       ! (handle case where downward flux at layer top > qmax)
 !       if(iLayerLiqFluxSoil(iSoil) > 0._dp)then
 !        iLayerLiqFluxSoil(iSoil)   = 0._dp
 !        iLayerLiqFluxSoil(iSoil-1) = qMax
 !       endif  ! if downward flux at layer top > qmax

 !      ! part2: downward flux at bottom of layer
 !      !  -- we can only exceed porosity here if downward flux at top exceeds downward flux at bottom
 !      !  -- reduce downward flux at top
 !      else
 !       iLayerLiqFluxSoil(iSoil-1) = iLayerLiqFluxSoil(iSoil) + qMax
 !       dq_dStateAbove(iSoil-1)    = 0._dp ! derivatives in the flux w.r.t. state variable in the layer above (m s-1 or s-1)
 !       dq_dStateBelow(iSoil-1)    = 0._dp ! derivatives in the flux w.r.t. state variable in the layer below (m s-1 or s-1)
 !       !print*, 'part2: iSoil, iLayerLiqFluxSoil(iSoil-1), qMax = ', iSoil, iLayerLiqFluxSoil(iSoil-1), qMax
 !      endif  ! (switch between downward and upward flux at top of layer)

 !     endif  ! (if exceed porosity)
 !    endif   ! (if layer has ice)
 !   end do ! (looping through soil layers)
 !   !pause 'check porosity'
 !  endif  ! (if previously modified fluxes to keep soil moisture within porosity)

 ! endif  ! if conducting the porosity fix


  end subroutine computeFlux



  ! ************************************************************************************************
  ! internal subroutine: compute the vector of fluxes at layer interfaces
  ! ************************************************************************************************
  subroutine cmpJacobian(&
                         ! input: model control
                         ixRichards,                   & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                         bc_lower,                     & ! intent(in): index defining the type of boundary conditions
                         ! input: diagnostic variables
                         mLayerBaseflowInput,          & ! intent(in): baseflow from each soil layer at the start of the step (m s-1)
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
  ! input: diagnostic variables
  real(dp),intent(in)              :: mLayerBaseflowInput(:)      ! baseflow from each soil layer at the start of the step (m s-1)
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
  real(dp),dimension(nLevels)      :: local_mLayerQMacroporeDeriv ! derivatives in the liquid flux to/from macropores (m s-1 [moisture form] s-1 [mixed form])
  real(dp)                         :: local_scalarSurfaceRunoff   ! surface runoff (m s-1) 
  ! local: fluxes used to calculate residual
  real(dp),dimension(nLevels)      :: mLayerTempQMacropore        ! liquid flux to/from macropores (m s-1)
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
  character(len=64),parameter      :: jFormat='(2(i4,1x),2(a,1x,3(f20.10,1x)))'  ! format string
  ! -------------------------------------------------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message='cmpJacobian/'

  print*, '**************************************'
  print*, '**************************************'
  print*, '*** start of jacobian calculations ***'

  ! initialize baseflow
  mLayerTempBaseflow = mLayerBaseflowInput

  ! get copies of the state vector to perturb
  mLayerMatricHeadTrial     = mLayerMatricHeadInput
  mLayerVolFracLiqTrial     = mLayerVolFracLiqInput
  scalarAquiferStorageTrial = scalarAquiferStorageInput

  ! loop through desired layers
  do iJac=1,nState

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
                    ! output: fluxes
                    iLayerTempLiqFluxSoil,       & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                    mLayerTempQMacropore,        & ! intent(out): liquid flux to/from macropores (m s-1)
                    mLayerTempBaseflow,          & ! intent(inout): baseflow flux from each soil layer (m s-1)
                    ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                    local_dq_dStateAbove,        & ! intent(out): derivatives in the flux w.r.t. state variable in the layer above (m s-1 or s-1)
                    local_dq_dStateBelow,        & ! intent(out): derivatives in the flux w.r.t. state variable in the layer below (m s-1 or s-1)
                    local_mLayerQMacroporeDeriv, & ! intent(out): derivative in liquid flux to/from macropores (m s-1 [moisture form] s-1 [mixed form])
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
                    iter,                      & ! intent(in): iteration count
                    ! coordinate variables
                    mLayerDepth,               & ! intent(in): depth of each layer (m)
                    ! input: model parameters for the compressibility term
                    theta_sat,                 & ! intent(in): porosity (-)
                    specificStorage,           & ! intent(in): specific storage (m-1)
                    ! initial flux vectors (start of the time step)
                    iLayerInitLiqFluxSoil,     & ! intent(in): initial liquid flux at layer interfaces (m s-1)
                    mLayerInitQMacropore,      & ! intent(in): initial liquid flux to/from macropores (m s-1)
                    mLayerInitBaseflow,        & ! intent(in): initial baseflow from each layer (m s-1)
                    mLayerInitTranspire,       & ! intent(in): initial transpiration from each layer -- from energy routine (m s-1)
                    ! trial flux vectors
                    iLayerTempLiqFluxSoil,     & ! intent(in): trial liquid flux at layer interfaces (m s-1)
                    mLayerTempQMacropore,      & ! intent(in): trial liquid flux to/from macropores (m s-1)
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
   if(iJac >  nLevels)then  ! (can ONLY happen if aquifer storage)
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
   endif
   !if(ijac == nState)then
   ! write(*,'(a,2(e20.10,1x))') 'scalarAquiferResidual, faquifer, dg/dS = ', scalarAquiferResidual, fAquifer, (-(fAquifer - scalarAquiferResidual)/eps - 1._dp)/dt
   ! pause
   !endif

   ! compute Jacobian
   if(ixGroundwater==bigBucket)then
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
   if(iJac==1)then;          write(*,jFormat) ixDerivMethod, iJac, 'test hyd Jacobian', (/valueMissing,      jmat(iJac,1:iJac+1)/), '--> tri-diag = ', valueMissing, diag(iJac)-local_mLayerdTheta_dPsi(iJac), d_p1(iJac)
   elseif(iJac==nState)then; write(*,jFormat) ixDerivMethod, iJac, 'test hyd Jacobian', (/jmat(iJac,iJac-1:nState), valueMissing/), '--> tri-diag = ', d_m1(iJac-1), diag(iJac)-local_mLayerdTheta_dPsi(iJac), valueMissing
   else;                     write(*,jFormat) ixDerivMethod, iJac, 'test hyd Jacobian', (/jmat(iJac,iJac-1:iJac+1)              /), '--> tri-diag = ', d_m1(iJac-1), diag(iJac)-local_mLayerdTheta_dPsi(iJac), d_p1(iJac) 
   endif
  end do
  pause 'testing Jacobian' 

  end subroutine cmpJacobian

 end subroutine soilHydrol_muster



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
                       mpExp,                 & ! intent(in): empirical exponent in macropore flow equation (-)
                       theta_sat,             & ! intent(in): soil porosity (-)
                       theta_res,             & ! intent(in): soil residual volumetric water content (-)
                       theta_mp,              & ! intent(in): volumetric liquid water content when macropore flow begins (-)
                       f_impede,              & ! intent(in): ice impedence factor (-)
                       ! input: saturated hydraulic conductivity
                       scalarSatHydCond,      & ! intent(in): saturated hydraulic conductivity at the mid-point of a given layer (m s-1)
                       scalarSatHydCondMP,    & ! intent(in): saturated hydraulic conductivity of macropores at the mid-point of a given layer (m s-1)
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
 USE soil_utils_module,only:volFracLiq      ! compute volumetric fraction of liquid water as a function of matric head
 USE soil_utils_module,only:hydCond_psi     ! compute hydraulic conductivity as a function of matric head
 USE soil_utils_module,only:hydCond_liq     ! compute hydraulic conductivity as a function of volumetric liquid water content
 USE soil_utils_module,only:hydCondMP_liq   ! compute hydraulic conductivity of macropores as a function of volumetric liquid water content
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
 real(dp),intent(in)           :: mpExp                     ! empirical exponent in macropore flow equation (-)
 real(dp),intent(in)           :: theta_sat                 ! soil porosity (-)
 real(dp),intent(in)           :: theta_res                 ! soil residual volumetric water content (-)
 real(dp),intent(in)           :: theta_mp                  ! volumetric liquid water content when macropore flow begins (-)
 real(dp),intent(in)           :: f_impede                  ! ice impedence factor (-)
 ! input: saturated hydraulic conductivity
 real(dp),intent(in)           :: scalarSatHydCond          ! saturated hydraulic conductivity at the mid-point of a given layer (m s-1)
 real(dp),intent(in)           :: scalarSatHydCondMP        ! saturated hydraulic conductivity of macropores at the mid-point of a given layer (m s-1)
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
 real(dp),intent(out)          :: dHydCond_dMatric          ! derivative in hydraulic conductivity w.r.t matric head (s-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                       ! error code
 character(*),intent(out)      :: message                   ! error message
 ! local variables
 real(dp)                      :: scalarHydCondMP           ! hydraulic conductivity of macropores at layer mid-points (m s-1)
 real(dp)                      :: dHydCondMacro_dVolLiq     ! derivative in hydraulic conductivity of macropores w.r.t volumetric liquid water content (m s-1) 
 real(dp)                      :: dHydCondMacro_dMatric     ! derivative in hydraulic conductivity of macropores w.r.t matric head (s-1)
 real(dp)                      :: dHydCondMicro_dMatric     ! derivative in hydraulic conductivity of micropores w.r.t matric head (s-1)
 real(dp)                      :: dPsi_dTheta2a             ! derivative in dPsi_dTheta (analytical)
 real(dp)                      :: dIceImpede_dLiq           ! derivative in ice impedence factor w.r.t. volumetric liquid water content (-)
 real(dp)                      :: hydCond_noIce             ! hydraulic conductivity in the absence of ice (m s-1)
 real(dp)                      :: dK_dLiq__noIce            ! derivative in hydraulic conductivity w.r.t volumetric liquid water content, in the absence of ice (m s-1)
 real(dp)                      :: dK_dPsi__noIce            ! derivative in hydraulic conductivity w.r.t matric head, in the absence of ice (s-1)
 real(dp)                      :: relSatMP                  ! relative saturation of macropores (-)
 !real(dp)                      :: x1,x2                     ! trial values of theta (-)
 !real(dp),parameter            :: dx = 1.e-8_dp             ! finite difference increment (m)
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
   !x1 = volFracLiq(scalarMatricHeadTrial,   vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   !x2 = volFracLiq(scalarMatricHeadTrial+dx,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   !print*, 'scalardTheta_dPsi = ', scalardTheta_dPsi, (x2 - x1)/dx
   scalardPsi_dTheta = valueMissing  ! (deliberately cause problems if this is ever used)
  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
 end select


 ! *****
 ! compute hydraulic conductivity and its derivative in each soil layer

 ! compute the ice impedence factor and its derivative w.r.t. volumetric liquid water content (-)
 call iceImpede(scalarVolFracIceTrial,scalarVolFracLiqTrial,theta_sat,f_impede,deriv_desired, &  ! input
                iceImpedeFac,dIceImpede_dLiq)                                                    ! output

 select case(ixRichards)
  ! ***** moisture-based form of Richards' equation
  case(moisture)
   ! haven't included macropores yet
   err=20; message=trim(message)//'still need to include macropores for the moisture-based form of Richards eqn'; return
   ! compute the hydraulic conductivity (m s-1) and diffusivity (m2 s-1) for a given layer
   hydCond_noIce = hydCond_liq(scalarVolFracLiqTrial,scalarSatHydCond,theta_res,theta_sat,vGn_m)
   scalarHydCond = hydCond_noIce*iceImpedeFac
   scalarDiffuse = scalardPsi_dTheta * scalarHydCond
   ! compute derivative in hydraulic conductivity (m s-1) and hydraulic diffusivity (m2 s-1)
   if(deriv_desired)then
    if(scalarVolFracIceTrial > epsilon(iceImpedeFac))then
     dK_dLiq__noIce   = dHydCond_dLiq(scalarVolFracLiqTrial,scalarSatHydCond,theta_res,theta_sat,vGn_m,.true.)  ! [.true. = analytical]
     dHydCond_dVolLiq = hydCond_noIce*dIceImpede_dLiq + dK_dLiq__noIce*iceImpedeFac
    else
     dHydCond_dVolLiq = dHydCond_dLiq(scalarVolFracLiqTrial,scalarSatHydCond,theta_res,theta_sat,vGn_m,.true.)
    endif
    dPsi_dTheta2a    = dPsi_dTheta2(scalarVolFracLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m,.true.)   ! [.true. = analytical] compute derivative in dPsi_dTheta (m)
    dDiffuse_dVolLiq = dHydCond_dVolLiq*scalardPsi_dTheta + scalarHydCond*dPsi_dTheta2a
    dHydCond_dMatric = valueMissing ! not used, so cause problems
   endif

  ! ***** mixed form of Richards' equation -- just compute hydraulic condictivity
  case(mixdform)
   ! compute the hydraulic conductivity (m s-1) and diffusivity (m2 s-1) for a given layer
   hydCond_noIce = hydCond_psi(scalarMatricHeadTrial,scalarSatHydCond,vGn_alpha,vGn_n,vGn_m)
   scalarDiffuse = valueMissing ! not used, so cause problems
   ! compute the hydraulic conductivity of macropores (m s-1)
   scalarHydCondMP = hydCondMP_liq(scalarvolFracLiqTrial,theta_sat,theta_mp,mpExp,scalarSatHydCondMP,scalarSatHydCond)
   scalarHydCond   = hydCond_noIce*iceImpedeFac + scalarHydCondMP
   ! compute derivative in hydraulic conductivity (m s-1)
   if(deriv_desired)then
    ! (compute derivative for macropores)
    if(scalarvolFracLiqTrial > theta_mp)then
     relSatMP              = (scalarvolFracLiqTrial - theta_mp)/(theta_sat - theta_mp) 
     dHydCondMacro_dVolLiq = ((scalarSatHydCondMP - scalarSatHydCond)/(theta_sat - theta_mp))*mpExp*(relSatMP**(mpExp - 1._dp))
     dHydCondMacro_dMatric = scalardTheta_dPsi*dHydCondMacro_dVolLiq
    else
     dHydCondMacro_dVolLiq = 0._dp
     dHydCondMacro_dMatric = 0._dp
    endif
    ! (compute derivative for micropores)
    if(scalarVolFracIceTrial > epsilon(iceImpedeFac))then
     dK_dPsi__noIce        = dHydCond_dPsi(scalarMatricHeadTrial,scalarSatHydCond,vGn_alpha,vGn_n,vGn_m,.true.)  ! analytical
     dHydCondMicro_dMatric = hydCond_noIce*dIceImpede_dLiq*scalardTheta_dPsi + dK_dPsi__noIce*iceImpedeFac
    else
     dHydCondMicro_dMatric = dHydCond_dPsi(scalarMatricHeadTrial,scalarSatHydCond,vGn_alpha,vGn_n,vGn_m,.true.)
    endif
    ! (combine derivatives)
    dHydCond_dMatric = dHydCondMicro_dMatric + dHydCondMacro_dMatric
    ! (set values that are not used to missing)
    dHydCond_dVolLiq = valueMissing ! not used, so cause problems
    dDiffuse_dVolLiq = valueMissing ! not used, so cause problems
   endif

  case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return

 endselect

 ! if derivatives are not desired, then set values to missing
 if(.not.deriv_desired)then
  dHydCond_dVolLiq   = valueMissing ! not used, so cause problems
  dDiffuse_dVolLiq   = valueMissing ! not used, so cause problems
  dHydCond_dMatric   = valueMissing ! not used, so cause problems
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
                       scalarMatricHead,          & ! intent(in): matric head in the upper-most soil layer (m)
                       scalarVolFracLiq,          & ! intent(in): volumetric liquid water content the upper-most soil layer (-)
                       scalarVolFracIce,          & ! intent(in): volumetric ice content in the upper-most soil layer (-)
                       ! input: depth of upper-most soil layer (m)
                       upperLayerDepth,           & ! intent(in): depth of upper-most soil layer (m)
                       ! input: boundary conditions
                       upperBoundHead,            & ! intent(in): upper boundary condition (m)
                       upperBoundTheta,           & ! intent(in): upper boundary condition (-)
                       ! input: flux at the upper boundary
                       scalarRainPlusMelt,        & ! intent(in): rain plus melt (m s-1)
                       scalarGroundEvaporation,   & ! intent(in): ground evaporation (m s-1)
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
                       qSurfScale,                & ! intent(in): scaling factor in the surface runoff parameterization (-)
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
 real(dp),intent(in)           :: scalarGroundEvaporation   ! ground evaporation (m s-1)
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
 real(dp),intent(in)           :: qSurfScale                ! scaling factor in the surface runoff parameterization (-)
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
 real(dp)                      :: fInfRaw                   ! infiltrating area before imposing solution constraints (-)
 real(dp)                      :: fInfArea                  ! area of the landscape where water infiltrates (-)
 real(dp)                      :: xMaxInfr                  ! maximum infiltration rate (m s-1)
 real(dp)                      :: uForcing                  ! forcing at the upper boundary, rain + melt + soil evaporation (m s-1)
 real(dp)                      :: fInfRate                  ! infiltration rate at the surface (m s-1)
 real(dp)                      :: dInfRaw                   ! derivative in the infiltrating area before imposing solution constraints (-)
 real(dp)                      :: dInfArea                  ! derivative in the infiltrating area w.r.t. matric head (m-1)
 real(dp)                      :: dInfRate                  ! derivative in the infiltration rate w.r.t matric head (s-1)
 real(dp),parameter            :: maxFracCap=0.999_dp       ! maximum fraction capacity -- used to avoid numerical problems associated with an enormous derivative
 real(dp),parameter            :: scaleFactor=0.00001_dp    ! scale factor for the smoothing function (-)
 logical(lgt)                  :: validCorrection           ! flag to denote if imposing solution constraints was appropriate
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

   ! check scaling factor is reasonable
   if(qSurfScale < 1._dp .or. qSurfScale > 100._dp)then
    message=trim(message)//'surface runoff scaling factor "qSurfScale" must be between 1 and 100'
    err=20; return
   endif

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
   fracCap  = (vFracLiq + scalarVolFracIce)/(theta_sat*maxFracCap)    ! fraction of capacity filled with liquid water and ice
   fInfRaw  = 1._dp - exp(-qSurfScale*(1._dp - fracCap))                ! infiltrating area -- allowed to violate solution constraints
   fInfArea = 0.5_dp*(fInfRaw + sqrt(fInfRaw**2._dp + scaleFactor))   ! infiltrating area -- constrained
   if((vFracLiq + scalarVolFracIce)/theta_sat > 1._dp - epsilon(theta_sat)) fInfArea=0._dp ! ensure no infiltration when storage exceeds capacity
   ! make adjustments for very small storages
   if(fInfArea < 1._dp)then
    validCorrection = .true.
   else
    validCorrection = .false.
    fInfArea = 1._dp
   endif
   ! compute the maximum infiltration rate -- assume head of zero at the surface (m s-1)
   select case(ixRichards)
    case(moisture); err=20; message=trim(message)//"infiltration not yet implemented for the moisture form of Richards' equation"; return
    case(mixdform); xMaxInfr = -surfaceSatHydCond*scalarMatricHead/(upperLayerDepth/2._dp) + surfaceSatHydCond
    case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
   endselect
   ! compute the total forcing at the upper boundary (m s-1)
   uForcing = scalarRainPlusMelt+scalarGroundEvaporation
   ! compute the rate of infiltration over the non-saturated area (m s-1)
   fInfRate = min(uForcing,xMaxInfr)
   ! compute the flux at the upper boundary
   scalarSurfaceInfiltration = fInfArea*fInfRate
   ! compute surface runoff (m s-1)
   scalarSurfaceRunoff = uForcing - scalarSurfaceInfiltration
   ! print progress
   print*, 'scalarRainPlusMelt = ', scalarRainPlusMelt
   print*, 'scalarGroundEvaporation = ', scalarGroundEvaporation
   print*, '(vFracLiq + scalarVolFracIce)/theta_sat = ', (vFracLiq + scalarVolFracIce)/theta_sat
   if(uForcing*iden_ice > 0.0001_dp)then
    write(*,'(a,1x,10(f20.10,1x))') 'qSurfScale, fInfArea, fInfRate, xMaxInfr, fracCap = ', qSurfScale, fInfArea, fInfRate, xMaxInfr, fracCap
    write(*,'(a,1x,5(e20.10,1x))') 'scalarMatricHead, surfaceSatHydCond, scalarRainPlusMelt, scalarGroundEvaporation, scalarSurfaceInfiltration = ', &
                                    scalarMatricHead, surfaceSatHydCond, scalarRainPlusMelt, scalarGroundEvaporation, scalarSurfaceInfiltration
   endif

   ! compute analytical derivatives (product rule)
   if(deriv_desired)then
    ! compute the derivative in the infiltration rate
    select case(ixRichards)
     case(moisture); err=20; message=trim(message)//"infiltration not yet implemented for the moisture form of Richards' equation"; return
     case(mixdform)
      if(uForcing > xMaxInfr)then
       dInfRate = -surfaceSatHydCond/(upperLayerDepth/2._dp)
      else
       dInfRate = 0._dp  ! NOTE: infiltration rate does not depend on state variables in the upper-most layer
      endif
     case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
    endselect
    ! compute the derivative in the infiltrating area
    if(validCorrection)then
     dInfRaw  = (-qSurfScale/(theta_sat*maxFracCap))*exp(-qSurfScale*(1._dp - fracCap))
    else
     dInfRaw  = 0._dp
    endif
    ! modify derivative in infiltrating area to account for solution constraints and form of Richards' equation
    select case(ixRichards)
     case(moisture); dInfArea =                   0.5_dp*(dInfRaw + dInfRaw*fInfRaw / sqrt(fInfRaw**2._dp + scaleFactor))
     case(mixdform); dInfArea = scalardTheta_dPsi*0.5_dp*(dInfRaw + dInfRaw*fInfRaw / sqrt(fInfRaw**2._dp + scaleFactor))
     case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
    end select
    ! put everything together -- product rule
    dq_dState = dInfArea*fInfRate + fInfArea*dInfRate
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
 logical(lgt),parameter        :: useGeometric=.false.        ! switch between the arithmetic and geometric mean
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
 ! compute the hydraulic conductivity at the interface
 if(useGeometric)then
  iLayerHydCond   = (nodeHydCondTrial(ixLower)   * nodeHydCondTrial(ixUpper))**0.5_dp
 else
  iLayerHydCond   = (nodeHydCondTrial(ixLower)   + nodeHydCondTrial(ixUpper))*0.5_dp
 endif
 !write(*,'(a,1x,5(e20.10,1x))') 'in iLayerFlux: iLayerHydCond, iLayerHydCondMP = ', iLayerHydCond, iLayerHydCondMP
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
    ! still need to implement arithmetric mean for the moisture-based form
    if(.not.useGeometric)then
     message=trim(message)//'only currently implemented for geometric mean -- change local flag'
     err=20; return
    endif
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
    if(useGeometric)then
     dHydCondIface_dMatricAbove = dHydCond_dMatric(ixUpper)*nodeHydCondTrial(ixLower) * 0.5_dp/max(iLayerHydCond,verySmall)
     dHydCondIface_dMatricBelow = dHydCond_dMatric(ixLower)*nodeHydCondTrial(ixUpper) * 0.5_dp/max(iLayerHydCond,verySmall)
    else
     dHydCondIface_dMatricAbove = dHydCond_dMatric(ixUpper)/2._dp
     dHydCondIface_dMatricBelow = dHydCond_dMatric(ixLower)/2._dp
    endif
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
 ! private subroutine: compute flux of water between micropores and macropores
 ! ************************************************************************************************
 subroutine qMacropore(&
                       ! input: model control
                       deriv_desired,         & ! intent(in): flag determining if the derivative is desired
                       ixRichards,            & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                       ! input: soil parameters
                       vGn_alpha,             & ! intent(in): van Genutchen "alpha" parameter (m-1)
                       vGn_n,                 & ! intent(in): van Genutchen "n" parameter (-)
                       VGn_m,                 & ! intent(in): van Genutchen "m" parameter (-)
                       theta_sat,             & ! intent(in): soil porosity (-)
                       theta_res,             & ! intent(in): soil residual volumetric water content (-)
                       ! input: state variable and diagnostic variables
                       scalarMatricHeadTrial, & ! intent(in): matric head in each layer (m)
                       scalarVolFracLiqTrial, & ! intent(in): volumetric liquid water content (-)
                       scalarVolFracIceTrial, & ! intent(in): volumetric ice content in each layer (-)
                       ! input: fluxes at the layer interfaces
                       scalarLiqFluxUpper,    & ! intent(in): liquid flux at the top of the layer (m s-1)
                       scalarLiqFluxLower,    & ! intent(in): liquid flux at the bottom of the layer (m s-1)
                       ! input: derivatives in liquid fluxes
                       dqUpper_dStateBelow,   & ! intent(in): derivatives in the flux at top of the layer w.r.t. state variable in the layer below (m s-1 or s-1)
                       dqLower_dStateAbove,   & ! intent(in): derivatives in the flux at bottom of the layer w.r.t. state variable in the layer above (m s-1 or s-1)
                       ! input: derivative in the soil water characteristic w.r.t. psi
                       scalardTheta_dPsi,     & ! intent(in): derivative in the soil water characteristic w.r.t. psi (m-1)
                       ! output: water flux to/from macropores, and derivative
                       scalarQMacropore,      & ! intent(out): water flux to/from macropores (m s-1)
                       scalarQMacroporeDeriv, & ! intent(out): derivative in macropore flow (m s-1 [moisture form] s-1 [mixed form])
                       ! output: error control
                       err,message)             ! intent(out): error control
 USE soil_utils_module,only:volFracLiq          ! compute volumetric fraction of liquid water (-)
 ! avoid super-saturation, by computing flux of water due to displacement (m s-1)
 implicit none
 ! input: model control
 logical(lgt),intent(in)       :: deriv_desired             ! flag determining if the derivative is desired
 integer(i4b),intent(in)       :: ixRichards                ! index defining the form of Richards' equation (moisture or mixdform)
 ! input: soil parameters
 real(dp),intent(in)           :: vGn_alpha                 ! van Genutchen "alpha" parameter (m-1)
 real(dp),intent(in)           :: vGn_n                     ! van Genutchen "n" parameter (-)
 real(dp),intent(in)           :: vGn_m                     ! van Genutchen "m" parameter (-)
 real(dp),intent(in)           :: theta_sat                 ! soil porosity (-)
 real(dp),intent(in)           :: theta_res                 ! soil residual volumetric water content (-)
 ! input: state variables
 real(dp),intent(in)           :: scalarMatricHeadTrial     ! matric head in each layer (m)
 real(dp),intent(in)           :: scalarVolFracLiqTrial     ! volumetric liquid water content in each layer (-)
 real(dp),intent(in)           :: scalarVolFracIceTrial     ! volumetric ice content in each layer (-)
 ! input: fluxes at the layer interfaces
 real(dp),intent(in)           :: scalarLiqFluxUpper        ! liquid flux at the top of the layer (m s-1)
 real(dp),intent(in)           :: scalarLiqFluxLower        ! liquid flux at the bottom of the layer (m s-1)
 ! input: derivatives in liquid fluxes
 real(dp),intent(in)           :: dqUpper_dStateBelow       ! derivatives in the flux at top of the layer w.r.t. state variable in the layer below (m s-1 or s-1)
 real(dp),intent(in)           :: dqLower_dStateAbove       ! derivatives in the flux at bottom of the layer w.r.t. state variable in the layer above (m s-1 or s-1)
 ! input: derivative in the soil water characteristic w.r.t. psi
 real(dp),intent(in)           :: scalardTheta_dPsi         ! derivative in the soil water characteristic w.r.t. psi (m-1)
 ! output: liquid flux to/from macropores and derivative
 real(dp),intent(out)          :: scalarQMacropore          ! liquid flux to/from macropores (m s-1)
 real(dp),intent(out)          :: scalarQMacroporeDeriv     ! derivative in liquid flux to/from macropores (m s-1 [moisture form] s-1 [mixed form])
 ! output: error control
 integer(i4b),intent(out)      :: err                       ! error code
 character(*),intent(out)      :: message                   ! error message
 ! local variables
 real(dp),parameter            :: supersatScale=0.001_dp    ! scaling factor for the logistic function (-)
 real(dp),parameter            :: xMatch = 0.999_dp         ! point where x-value and function value match (-)
 real(dp),parameter            :: fSmall = epsilon(xMatch)  ! smallest possible value to test
 real(dp)                      :: supersatThresh            ! threshold in super-saturation function (-)
 real(dp)                      :: fracMin                   ! minimum fraction of pore space required to be filled in order for lateral flow to occur (-)
 real(dp)                      :: fracCap                   ! fraction of pore space filled with liquid water and ice (-)
 real(dp)                      :: expFunc                   ! exponential function used as part of the flux calculation (-)
 real(dp)                      :: dLogFunc_dLiq             ! derivative in the logistic function w.r.t. volumetric liquid water content (-)
 ! initialize error control
 err=0; message="qMacropore/"

 ! check for an early return -- only interested in cases that may exceed porosity
 if(scalarLiqFluxUpper < scalarLiqFluxLower)then
  scalarQMacropore      = 0._dp
  scalarQMacroporeDeriv = 0._dp
  return
 endif

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

  ! ***** calculate the flux from micropores to macropores (m s-1)
  expFunc = exp((supersatThresh - fracCap)/supersatScale)
  scalarQMacropore = (scalarLiqFluxUpper - scalarLiqFluxLower)/(1._dp + expFunc)
  if(scalarQMacropore< 0._dp)then
   print*, 'scalarQMacropore = ', scalarQMacropore
   print*, 'scalarLiqFluxUpper = ', scalarLiqFluxUpper
   print*, 'scalarLiqFluxLower = ', scalarLiqFluxLower
   print*, 'expFunc = ', expFunc
   message=trim(message)//'liquid flux from micropores to macropores is < 0'
   err=20; return
  endif
  !write(*,'(a,4(e20.10,1x))') 'fracCap, scalarLiqFluxUpper, scalarLiqFluxLower, 1._dp/(1._dp + expFunc) = ', &
  !                             fracCap, scalarLiqFluxUpper, scalarLiqFluxLower, 1._dp/(1._dp + expFunc)

  ! ***** calculate the derivative in the flux from micropores to macropores (m s-1)
  if(deriv_desired)then
   ! calculate the derivative in the logistic function w.r.t. volumetric liquid water content (-)
   dLogFunc_dLiq = (expFunc/(theta_sat*supersatScale)) * (1._dp + expFunc)**(-2._dp)
   ! combine the logistic function with the flux derivatives (m s-1)
   select case(ixRichards)
    case(moisture); scalarQMacroporeDeriv = (dqUpper_dStateBelow/(1._dp + expFunc) + scalarLiqFluxUpper*dLogFunc_dLiq) - &
                                            (dqLower_dStateAbove/(1._dp + expFunc) + scalarLiqFluxLower*dLogFunc_dLiq)
    case(mixdform); scalarQMacroporeDeriv = (dqUpper_dStateBelow/(1._dp + expFunc) + scalarLiqFluxUpper*dLogFunc_dLiq*scalardTheta_dPsi) - &
                                            (dqLower_dStateAbove/(1._dp + expFunc) + scalarLiqFluxLower*dLogFunc_dLiq*scalardTheta_dPsi)
    case default; err=10; message=trim(message)//"unable to identify option for Richards' equation"; return
   end select
  else  ! (if derivative is not desired)
   scalarQMacroporeDeriv = valueMissing
  endif

 else  ! (if fraction of pore space filled with liquid water and ice is less than the minimum value required for lateral flow)
  scalarQMacropore      = 0._dp
  scalarQMacroporeDeriv = 0._dp
 endif ! (if fraction of pore space filled with liquid water and ice is less than the minimum value required for lateral flow)

 end subroutine qMacropore




 ! ************************************************************************************************
 ! private subroutine: compute residual in the water balance
 ! ************************************************************************************************
 subroutine liqResidual(&
                        ! control variables
                        dt,                              & ! intent(in): length of the time step (s)
                        wimplicit,                       & ! intent(in): weight assigned to the start-of-step
                        iter,                            & ! intent(in): iteration count
                        ! coordinate variables
                        mLayerDepth,                     & ! depth of each layer (m)
                        ! input: model parameters for the compressibility term
                        theta_sat,                       & ! intent(in): porosity (-)
                        specificStorage,                 & ! intent(in): specific storage (m-1)
                        ! initial flux vectors (start of the time step)
                        iLayerInitLiqFluxSoil,           & ! intent(in): initial liquid flux at layer interfaces (m s-1)
                        mLayerInitQMacropore,            & ! intent(in): initial liquid flux to/from macropores (m s-1)
                        mLayerInitBaseflow,              & ! intent(in): initial baseflow from each layer (m s-1)
                        mLayerInitTranspire,             & ! intent(in): initial transpiration from each layer -- from energy routine (m s-1)
                        ! trial flux vectors
                        iLayerTrialLiqFluxSoil,          & ! intent(in): trial liquid flux at layer interfaces (m s-1)
                        mLayerTrialQMacropore,           & ! intent(in): trial liquid flux to/from macropores (m s-1)
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
 integer(i4b),intent(in)      :: iter                      ! iteration count
 ! coordinate variables
 real(dp),intent(in)          :: mLayerDepth(:)            ! depth of each layer (m)
 ! model parameters for the compressibility term
 real(dp),intent(in)          :: theta_sat                 ! intent(in): porosity (-)
 real(dp),intent(in)          :: specificStorage           ! intent(in): specific storage (m-1)
 ! initial flux vectors (start of the time step)
 real(dp),intent(in)          :: iLayerInitLiqFluxSoil(0:) ! intent(in): initial liquid flux at layer interfaces (m s-1)
 real(dp),intent(in)          :: mLayerInitQMacropore(:)   ! intent(in): initial liquid flux to/from macropores (m s-1)
 real(dp),intent(in)          :: mLayerInitBaseflow(:)     ! intent(in): initial baseflow from each layer (m s-1)
 real(dp),intent(in)          :: mLayerInitTranspire(:)    ! intent(in): initial transpiration from each layer -- from energy routine (m s-1)
 ! trial flux vectors
 real(dp),intent(in)          :: iLayerTrialLiqFluxSoil(0:) ! intent(in): trial liquid flux at layer interfaces (m s-1)
 real(dp),intent(in)          :: mLayerTrialQMacropore(:)  ! intent(in): trial liquid flux to/from macropores (m s-1)
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
 real(dp)                     :: mPore                     ! overall contribution to volumetric liquid water from flow to/from macropores (-) 
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
  compressibility = 0._dp !(specificStorage*mLayerTrialVolFracLiq(iLayer)/theta_sat) * (mLayerTrialMatricHead(iLayer) - mLayerInitMatricHead(iLayer))
  ! transpiration (-)
  mEvap = wimplicit*mLayerInitTranspire(iLayer)*dt_dz + (1._dp - wimplicit)*mLayerTrialTranspire(iLayer)*dt_dz
  ! flow to/from macropores (-)
  mPore = wimplicit*mLayerInitQMacropore(iLayer)*dt_dz + (1._dp - wimplicit)*mLayerTrialQMacropore(iLayer)*dt_dz
  ! baseflow (-)
  mBase = wimplicit*mLayerInitBaseflow(iLayer)*dt_dz + (1._dp - wimplicit)*mLayerTrialBaseflow(iLayer)*dt_dz
  ! phase change (-)
  mPhse = (iden_ice/iden_water)*(mLayerTrialVolFracIce(iLayer)-mLayerInitVolFracIce(iLayer))
  ! residual (-)
  residualVec(iLayer) = mLayerTrialVolFracLiq(iLayer) - (mLayerInitVolFracLiq(iLayer) + mFlux + mEvap - mPore - mBase - mPhse - compressibility)

  ! print progress
  !if(iter > 30)then
   if(iLayer==1)  write(*,'(a)') 'iter, iLayer, residualVec(iLayer), mLayerTrialVolFracLiq(iLayer), mLayerInitVolFracLiq(iLayer), top flux, mFlux, mEvap, mPore, mBase, mPhse'
   if(iLayer < 10) write(*,'(2(i4,1x),10(e20.10,1x))') iter, iLayer, residualVec(iLayer), mLayerTrialVolFracLiq(iLayer), mLayerInitVolFracLiq(iLayer), iLayerTrialLiqFluxSoil(iLayer-1), mFlux, mEvap, mPore, mBase, mPhse
  !endif

 end do  ! (looping through soil layers)

 end subroutine liqResidual




 ! ************************************************************************************************
 ! private subroutine: compute residual in the groundwater balance
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

 err=20; message=trim(message)//'in gw residual -- should not be here'; return
 
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

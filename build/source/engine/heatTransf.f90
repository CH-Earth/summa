module heatTransf_module
USE nrtype
! physical constants
USE multiconst,only:&
                    sb,          & ! Stefan Boltzman constant      (W m-2 K-4)
                    Em_Sno,      & ! emissivity of snow            (-)
                    Cp_air,      & ! specific heat of air          (J kg-1 K-1)
                    Cp_water,    & ! specifric heat of water       (J kg-1 K-1)
                    LH_fus,      & ! latent heat of fusion         (J kg-1)
                    LH_vap,      & ! latent heat of vaporization   (J kg-1)
                    LH_sub,      & ! latent heat of sublimation    (J kg-1)
                    Tfreeze,     & ! freezing point of pure water  (K)
                    iden_air,    & ! intrinsic density of air      (kg m-3)
                    iden_ice,    & ! intrinsic density of ice      (kg m-3)
                    iden_water     ! intrinsic density of water    (kg m-3)
! named variables for snow and soil
USE data_struc,only:ix_soil,ix_snow                        ! names variables for snow and soil
! provide access to look-up values for model decisions
USE mDecisions_module,only:      &
 ! look-up values for the numerical method
 iterative,                      & ! iterative
 nonIterative,                   & ! non-iterative
 iterSurfEnergyBal,              & ! iterate only on the surface energy balance
 ! look-up values for method used to compute derivative
 numerical,                      & ! numerical solution
 analytical,                     & ! analytical solution
 ! look-up values for choice of boundary conditions for thermodynamics
 prescribedTemp,                 & ! prescribed temperature
 energyFlux,                     & ! energy flux
 zeroFlux,                       & ! zero flux
 ! look-up values for choice of boundary conditions for soil hydrology
 prescribedHead                    ! prescribed head
! -------------------------------------------------------------------------------------------------
implicit none
private
public::heatTransf
! global parameters
real(dp),parameter            :: RHsurf=1._dp             ! relative humidity of the surface (-)
real(dp),parameter            :: dx=1.e-8_dp              ! finite difference increment (K)
real(dp),parameter            :: valueMissing=-9999._dp   ! missing value parameter
real(dp),parameter            :: missingValue_belowAbsoluteZero=-9999._dp ! missing value for canopy temperature (must be < 0) 
! number of soil and snow layers
integer(i4b)                  :: nSoil                    ! number of soil layers
integer(i4b)                  :: nSnow                    ! number of snow layers
integer(i4b)                  :: nLayers                  ! total number of layers
contains


 ! ************************************************************************************************
 ! new subroutine: compute change in temperature over the time step
 ! ************************************************************************************************
 subroutine heatTransf(&

                       ! input
                       dt,                               & ! intent(in): time step (seconds)
                       iter,                             & ! intent(in): current iteration count
                       firstSubstep,                     & ! intent(in): flag to indicate if we are processing the first sub-step
                       computeVegFlux,                   & ! intent(in): flag to indicate if we computing fluxes ovser vegetation (.false. means veg is buried with snow)
                       scalarCanairTempIter,             & ! intent(in): trial temperature of the canopy air space (K)
                       scalarCanopyTempIter,             & ! intent(in): trial temperature of the vegetation canopy at the current iteration (K)
                       scalarCanopyIceIter,              & ! intent(in): trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
                       scalarCanopyLiqIter,              & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                       mLayerTempIter,                   & ! intent(in): trial temperature of each model layer at the current iteration (K)
                       mLayerVolFracIceIter,             & ! intent(in): trial volumetric fraction of ice at the current iteration (-)
                       mLayerVolFracLiqIter,             & ! intent(in): trial volumetric fraction of liquid water at the current iteration (-)
                       mLayerMatricHeadIter,             & ! intent(in): trial matric head at the current iteration (m)
                       canopyTempIncrOld,                & ! intent(in): previous iteration increment in canopy temperature (K)
                       mLayerTempIncrOld,                & ! intent(in): previous iteration increment in temperature of the snow-soil vector (K)

                       ! input/output variables from heatTransf subroutine: canopy air space variables
                       scalarVP_CanopyAir,               & ! intent(inout): trial vapor pressure of the canopy air space (Pa)
                       scalarCanopyStabilityCorrection,  & ! intent(inout): stability correction for the canopy (-)
                       scalarGroundStabilityCorrection,  & ! intent(inout): stability correction for the ground surface (-)

                       ! output
                       scalarCanopyTranspiration,        & ! intent(out): canopy transpiration (kg m-2 s-1)
                       scalarCanopyEvaporation,          & ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
                       scalarGroundEvaporation,          & ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
                       scalarCanairTempDiff,             & ! intent(out): iteration increment for temperature of the canopy air space (K)
                       scalarCanopyTempDiff,             & ! intent(out): iteration increment for temperature of the vegetation canopy (K)
                       mLayerTempDiff,                   & ! intent(out): iteration increment for temperature of the snow-soil system (K)
                       scalarCanairTempNew,              & ! intent(out): new temperature of the canopy air space (K)
                       scalarCanopyTempNew,              & ! intent(out): new temperature of the vegetation canopy (K)
                       scalarCanopyIceNew,               & ! intent(out): mass of ice on the canopy (kg m-2) 
                       scalarCanopyLiqNew,               & ! intent(out): mass of liquid water on the canopy (kg m-2)
                       mLayerTempNew,                    & ! intent(out): new temperature of each model layer (K)
                       mLayerVolFracIceNew,              & ! intent(out): new volumetric fraction of ice (-)
                       mLayerVolFracLiqNew,              & ! intent(out): new volumetric fraction of liquid water (-)
                       mLayerMatricHeadNew,              & ! intent(out): new matric head (m)

                       ! error control
                       err,message)                        ! intent(out): error control
 ! model decisions
 USE data_struc,only:model_decisions                       ! model decision structure
 USE var_lookup,only:iLookDECISIONS                        ! named variables for elements of the decision structure
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:attr_data,type_data,mpar_data,forc_data,mvar_data,indx_data    ! data structures
 USE var_lookup,only:iLookATTR,iLookTYPE,iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX ! named variables for structure elements
 ! compute change in temperature over the time step
 implicit none
 ! input
 real(dp),intent(in)           :: dt                            ! time step (seconds)
 integer(i4b),intent(in)       :: iter                          ! iteration count
 logical(lgt),intent(in)       :: firstSubStep                  ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(in)       :: computeVegFlux                ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 real(dp),intent(in)           :: scalarCanairTempIter          ! trial temperature of the canopy air space (K)
 real(dp),intent(in)           :: scalarCanopyTempIter          ! trial temperature of the vegetation canopy at the current iteration (K)
 real(dp),intent(in)           :: scalarCanopyIceIter           ! trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(in)           :: scalarCanopyLiqIter           ! trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(in)           :: mLayerTempIter(:)             ! trial temperature of each snow/soil layer at the current iteration (K)
 real(dp),intent(in)           :: mLayerVolFracIceIter(:)       ! trial volumetric fraction of ice at the current iteration (-)
 real(dp),intent(in)           :: mLayerVolFracLiqIter(:)       ! trial volumetric fraction of liquid water at the current iteration (-)
 real(dp),intent(in)           :: mLayerMatricHeadIter(:)       ! trial matric head at the current iteration (m)
 real(dp),intent(in)           :: canopyTempIncrOld             ! iteration increment for temperature of the vegetation canopy in the previous iteration (K)
 real(dp),intent(in)           :: mLayerTempIncrOld(:)          ! iteration increment for temperature of the snow-soil vector in the previous iteration (K)
 ! input/output variables from the heatTransf subroutine: canopy air space variables
 real(dp),intent(inout)        :: scalarVP_CanopyAir            ! trial vapor pressure of the canopy air space (Pa)
 real(dp),intent(inout)        :: scalarCanopyStabilityCorrection ! stability correction for the canopy (-)
 real(dp),intent(inout)        :: scalarGroundStabilityCorrection ! stability correction for the ground surface (-)
 ! output
 real(dp),intent(out)          :: scalarCanopyTranspiration     ! canopy transpiration (kg m-2 s-1)
 real(dp),intent(out)          :: scalarCanopyEvaporation       ! canopy evaporation/condensation (kg m-2 s-1)
 real(dp),intent(out)          :: scalarGroundEvaporation       ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
 real(dp),intent(out)          :: scalarCanairTempDiff          ! iteration increment for temperature of the canopy air space (K)
 real(dp),intent(out)          :: scalarCanopyTempDiff          ! iteration increment for temperature of the vegetation canopy (K)
 real(dp),intent(out)          :: mLayerTempDiff(:)             ! iteration increment for temperature of the snow-soil system (K)
 real(dp),intent(out)          :: scalarCanairTempNew           ! new temperature of the canopy air space (K)
 real(dp),intent(out)          :: scalarCanopyTempNew           ! new temperature of the vegetation canopy (K)
 real(dp),intent(out)          :: scalarCanopyIceNew            ! mass of ice on the canopy (kg m-2) 
 real(dp),intent(out)          :: scalarCanopyLiqNew            ! mass of liquid water on the canopy (kg m-2)
 real(dp),intent(out)          :: mLayerTempNew(:)              ! new temperature of each snow/soil layer (K)
 real(dp),intent(out)          :: mLayerVolFracIceNew(:)        ! new volumetric fraction of ice (-)
 real(dp),intent(out)          :: mLayerVolFracLiqNew(:)        ! new volumetric fraction of liquid water (-)
 real(dp),intent(out)          :: mLayerMatricHeadNew(:)        ! new matric head (m)
 integer(i4b),intent(out)      :: err                           ! error code
 character(*),intent(out)      :: message                       ! error message
 ! internal
 character(LEN=256)            :: cmessage                      ! error message of downwind routine
 ! initialize error control
 err=0; message="heatTransf/"

 ! define number of layers
 nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat == ix_soil)  ! number of soil layers
 nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat == ix_snow)  ! number of snow layers
 nLayers = indx_data%var(iLookINDEX%nLayers)%dat(1)                   ! total number of layers

 ! *****
 ! wrapper for the temperature change sub-routine...
 ! *************************************************
 ! used as a trick to both
 !  1) save typing; and
 !  2) "protect" variables through the use of the intent attribute
 call heatTransf_muster(&

                        ! input variables from heatTransf routine
                        dt,                                                            & ! intent(in): time step (seconds)
                        iter,                                                          & ! intent(in): current iteration count
                        firstSubstep,                                                  & ! intent(in): flag to indicate if we are processing the first sub-step
                        computeVegFlux,                                                & ! intent(in): flag to indicate if we computing fluxes ovser vegetation (.false. means veg is buried with snow)
                        scalarCanairTempIter,                                          & ! intent(in): trial temperature of the canopy air space (K)
                        scalarCanopyTempIter,                                          & ! intent(in): trial temperature of the vegetation canopy at the current iteration (K)
                        scalarCanopyIceIter,                                           & ! intent(in): trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
                        scalarCanopyLiqIter,                                           & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                        mLayerTempIter,                                                & ! intent(in): trial temperature of each snow/soil layer at the current iteration (K)
                        mLayerVolFracIceIter,                                          & ! intent(in): trial volumetric fraction of ice at the current iteration (-)
                        mLayerVolFracLiqIter,                                          & ! intent(in): trial volumetric fraction of liquid water at the current iteration (-)
                        mLayerMatricHeadIter,                                          & ! intent(in): trial matric head at the current iteration (m)
                        canopyTempIncrOld,                                             & ! intent(in): previous iteration increment in canopy temperature (K)
                        mLayerTempIncrOld,                                             & ! intent(in): previous iteration increment in the temperature of the snow-soil vector  (K)

                        ! model decisions
                        model_decisions(iLookDECISIONS%num_method)%iDecision,          & ! intent(in): choice of numerical method
                        model_decisions(iLookDECISIONS%fDerivMeth)%iDecision,          & ! intent(in): method used to calculate flux derivatives
                        model_decisions(iLookDECISIONS%bcLowrTdyn)%iDecision,          & ! intent(in): type of lower boundary condition for thermodynamics

                        ! index variables
                        indx_data%var(iLookINDEX%nLayers)%dat(1),                      & ! intent(in): number of layers
                        indx_data%var(iLookINDEX%layerType)%dat,                       & ! intent(in): layer type (ix_soil or ix_snow)

                        ! physical attributes
                        type_data%var(iLookTYPE%vegTypeIndex),                         & ! intent(in): vegetation type index
                        type_data%var(iLookTYPE%soilTypeIndex),                        & ! intent(in): soil type index

                        ! general model parameters
                        mpar_data%var(iLookPARAM%wimplicit),                           & ! intent(in): weight assigned to start-of-step fluxes (-)
                        mpar_data%var(iLookPARAM%snowfrz_scale),                       & ! intent(in): scaling parameter for the snow freezing curve (K-1)
                        mpar_data%var(iLookPARAM%lowerBoundTemp),                      & ! intent(in): temperature of the lower boundary (K)

                        ! vegetation parameters
                        mpar_data%var(iLookPARAM%heightCanopyTop),                     & ! intent(in): height of top of the vegetation canopy above ground surface (m)
                        mpar_data%var(iLookPARAM%heightCanopyBottom),                  & ! intent(in): height of bottom of the vegetation canopy above ground surface (m)

                        ! soil parameters
                        mpar_data%var(iLookPARAM%vGn_alpha),                           & ! intent(in): van Genutchen "alpha" parameter (m-1)
                        mpar_data%var(iLookPARAM%vGn_n),                               & ! intent(in): van Genutchen "n" parameter (-)
                        mpar_data%var(iLookPARAM%theta_sat),                           & ! intent(in): soil porosity (-)
                        mpar_data%var(iLookPARAM%theta_res),                           & ! intent(in): soil residual volumetric water content (-)
                        mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1),                   & ! intent(in): van Genutchen "m" parameter (-)
                        mvar_data%var(iLookMVAR%scalarKappa)%dat(1),                   & ! intent(in): constant used to express matric head as a function of temperature (m K-1)

                        ! model state variables
                        ! NOTE: start-of-sub-step values -- protected by the intent(in) attribute
                        mvar_data%var(iLookMVAR%scalarCanairTemp)%dat(1),              & ! intent(in): temperature of the canopy air space (K)
                        mvar_data%var(iLookMVAR%scalarCanopyTemp)%dat(1),              & ! intent(in): temperature of the vegetation canopy (K)
                        mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1),               & ! intent(in): mass of ice on the vegetation canopy (kg m-2)
                        mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1),               & ! intent(in): mass of liquid water on the vegetation canopy (kg m-2)
                        mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1),               & ! intent(in): snow depth on the ground surface (m)
                        mvar_data%var(iLookMVAR%mLayerTemp)%dat,                       & ! intent(in): temperature of each layer (K)
                        mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat,                 & ! intent(in): volumetric fraction of ice in each layer (-)
                        mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat,                 & ! intent(in): volumetric fraction of liquid water in each layer (-)

                        ! model cooordinate variables (input)
                        mvar_data%var(iLookMVAR%mLayerDepth)%dat,                      & ! intent(in): depth of each layer (m)
                        mvar_data%var(iLookMVAR%mLayerHeight)%dat,                     & ! intent(in): height at the mid-point of each layer (m)

                        ! model diagnostic variables for vegetation (intent in) -- phenology and thermal properties constant over iterations
                        mvar_data%var(iLookMVAR%scalarLAI)%dat(1),                     & ! intent(in): one-sided leaf area index (m2 m-2)
                        mvar_data%var(iLookMVAR%scalarSAI)%dat(1),                     & ! intent(in): one-sided stem area index (m2 m-2)
                        mvar_data%var(iLookMVAR%scalarExposedLAI)%dat(1),              & ! intent(in): exposed leaf area index after burial by snow (m2 m-2)
                        mvar_data%var(iLookMVAR%scalarExposedSAI)%dat(1),              & ! intent(in): exposed stem area index after burial by snow (m2 m-2)
                        mvar_data%var(iLookMVAR%scalarGrowingSeasonIndex)%dat(1),      & ! intent(in): growing season index (0=off, 1=on)
                        mvar_data%var(iLookMVAR%scalarFoliageNitrogenFactor)%dat(1),   & ! intent(in): foliage nitrogen concentration (1.0 = saturated)
                        mvar_data%var(iLookMVAR%scalarBulkVolHeatCapVeg)%dat(1),       & ! intent(in): bulk volumetric heat capacity of vegetation (J m-3 K-1)

                        ! model diagnostic variables from the hydrology routines (intent in) -- liquid water fluxes at layer interfaces
                        mvar_data%var(iLookMVAR%iLayerLiqFluxSnow)%dat,                & ! intent(in): liquid flux at the interface of each snow layer (m s-1)
                        mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat,                & ! intent(in): liquid flux at the interface of each soil layer (m s-1)

                        ! model diagnostic variables for snow-soil vector (intent in) -- thermal properties constant over iterations
                        mvar_data%var(iLookMVAR%mLayerTcrit)%dat,                      & ! intent(in): critical soil temperature above which all water is unfrozen (K)
                        mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat,               & ! intent(in): volumetric heat capacity in each layer (J m-3 K-1) 
                        mvar_data%var(iLookMVAR%iLayerThermalC)%dat,                   & ! intent(in): thermal conductivity at the interface of each layer (W m-1 K-1)

                        ! model diagnostic variables (output)
                        mvar_data%var(iLookMVAR%mLayerdTheta_dTk)%dat,                 & ! intent(out): derivative in the freezing curve (K-1)
                        mvar_data%var(iLookMVAR%iLayerConductiveFlux)%dat,             & ! intent(out): conductive energy flux at layer interfaces at end of time step (W m-2)
                        mvar_data%var(iLookMVAR%iLayerAdvectiveFlux)%dat,              & ! intent(out): advective energy flux at layer interfaces at end of time step (W m-2)
                        mvar_data%var(iLookMVAR%iLayerInitNrgFlux)%dat,                & ! intent(out): energy flux at layer interfaces at the start of the time step (W m-2)
                        mvar_data%var(iLookMVAR%iLayerNrgFlux)%dat,                    & ! intent(out): energy flux at layer interfaces at the end of the time step (W m-2)

                        ! input/output variables from heatTransf subroutine: canopy air space variables
                        scalarVP_CanopyAir,                                            & ! intent(inout): trial vapor pressure of the canopy air space (Pa)
                        scalarCanopyStabilityCorrection,                               & ! intent(inout): stability correction for the canopy (-)
                        scalarGroundStabilityCorrection,                               & ! intent(inout): stability correction for the ground surface (-)

                        ! output variables from heatTransf subroutine
                        scalarCanopyTranspiration,                                     & ! intent(out): canopy transpiration (kg m-2 s-1)
                        scalarCanopyEvaporation,                                       & ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
                        scalarGroundEvaporation,                                       & ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
                        scalarCanairTempDiff,                                          & ! intent(out): iteration increment for temperature of the canopy air space (K)
                        scalarCanopyTempDiff,                                          & ! intent(out): iteration increment for temperature of the vegetation canopy (K)
                        mLayerTempDiff,                                                & ! intent(out): iteration increment for temperature of the snow-soil system (K) 
                        scalarCanairTempNew,                                           & ! intent(out): new temperature of the canopy air space (K)
                        scalarCanopyTempNew,                                           & ! intent(out): new temperature of the vegetation canopy (K)
                        scalarCanopyIceNew,                                            & ! intent(out): mass of ice on the canopy (kg m-2) 
                        scalarCanopyLiqNew,                                            & ! intent(out): mass of liquid water on the canopy (kg m-2)
                        mLayerTempNew,                                                 & ! intent(out): new temperature of each snow/soil layer (K)
                        mLayerVolFracIceNew,                                           & ! intent(out): new volumetric fraction of ice (-)
                        mLayerVolFracLiqNew,                                           & ! intent(out): new volumetric fraction of liquid water (-)
                        mLayerMatricHeadNew,                                           & ! intent(out): new matric head (m)
                        err,cmessage)                                                    ! intent(out): error control

 ! check for errors
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 end subroutine heatTransf



 
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ***** PRIVATE SUBROUTINES **********************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************
 ! ************************************************************************************************

 ! ************************************************************************************************
 ! private subroutine: wrapper for the temperature change subroutine
 ! ************************************************************************************************
 subroutine heatTransf_muster(&

                              ! input variables from heatTransf routine
                              dt,                           & ! intent(in): time step (seconds)
                              iter,                         & ! intent(in): current iteration count
                              firstSubstep,                 & ! intent(in): flag to indicate if we are processing the first sub-step
                              computeVegFlux,               & ! intent(in): flag to indicate if we computing fluxes ovser vegetation (.false. means veg is buried with snow)
                              scalarCanairTempIter,         & ! intent(in): trial temperature of the canopy air space (K)
                              scalarCanopyTempIter,         & ! intent(in): trial temperature of the vegetation canopy (K)
                              scalarCanopyIceIter,          & ! intent(in): trial mass of ice on the vegetation canopy (kg m-2)
                              scalarCanopyLiqIter,          & ! intent(in): trial mass of liquid water on the vegetation canopy (kg m-2)
                              mLayerTempIter,               & ! intent(in): trial temperature of each snow/soil layer at the current iteration (K)
                              mLayerVolFracIceIter,         & ! intent(in): trial volumetric fraction of ice at the current iteration (-)
                              mLayerVolFracLiqIter,         & ! intent(in): trial volumetric fraction of liquid water at the current iteration (-)
                              mLayerMatricHeadIter,         & ! intent(in): trial matric head at the current iteration (m)
                              canopyTempIncrOld,            & ! intent(in): previous iteration increment in canopy temperature (K)
                              mLayerTempIncrOld,            & ! intent(in): previous iteration increment in temperature of the snow-soil vector (K)

                              ! model decisions
                              num_method,                   & ! intent(in): choice of numerical method
                              fDerivMeth,                   & ! intent(in): method used to calculate flux derivatives
                              bcLowrTdyn,                   & ! intent(in): type of lower boundary condition for thermodynamics

                              ! index variables
                              nLayers,                      & ! intent(in): number of layers
                              layerType,                    & ! intent(in): layer type (ix_soil or ix_snow)

                              ! physical attributes
                              vegTypeIndex,                 & ! intent(in): vegetation type index
                              soilTypeIndex,                & ! intent(in): soil type index

                              ! general model parameters
                              wimplicit,                    & ! intent(in): weight assigned to start-of-step fluxes (-)
                              snowfrz_scale,                & ! intent(in): scaling parameter for the snow freezing curve (K-1)
                              lowerBoundTemp,               & ! intent(in): temperature of the lower boundary (K)

                              ! vegetation parameters
                              heightCanopyTop,              & ! intent(in): height of top of the vegetation canopy above ground surface (m)
                              heightCanopyBottom,           & ! intent(in): height of bottom of the vegetation canopy above ground surface (m)

                              ! soil parameters
                              vGn_alpha,                    & ! intent(in): van Genutchen "alpha" parameter (m-1)
                              vGn_n,                        & ! intent(in): van Genutchen "n" parameter (-)
                              theta_sat,                    & ! intent(in): soil porosity (-)
                              theta_res,                    & ! intent(in): soil residual volumetric water content (-)
                              vGn_m,                        & ! intent(in): van Genutchen "m" parameter (-)
                              kappa,                        & ! intent(in): constant used to express matric head as a function of temperature (m K-1)

                              ! model state variables
                              ! NOTE: start-of-sub-step values -- protected by the intent(in) attribute
                              scalarCanairTemp,             & ! intent(in): temperature of the canopy air space (K)
                              scalarCanopyTemp,             & ! intent(in): temperature of the vegetation canopy (K)
                              scalarCanopyIce,              & ! intent(in): mass of ice on the vegetation canopy (kg m-2)
                              scalarCanopyLiq,              & ! intent(in): mass of liquid water on the vegetation canopy (kg m-2)
                              scalarSnowDepth,              & ! intent(in): snow depth on the ground surface (m)
                              mLayerTemp,                   & ! intent(in): temperature of each snow/soil layer (K)
                              mLayerVolFracIce,             & ! intent(in): volumetric fraction of ice in each layer (-)
                              mLayerVolFracLiq,             & ! intent(in): volumetric fraction of liquid water in each layer (-)

                              ! model cooordinate variables (input)
                              mLayerDepth,                  & ! intent(in): depth of each layer (m)
                              mLayerHeight,                 & ! intent(in): height at the mid-point of each layer (m)

                              ! model diagnostic variables for vegetation (intent in) -- phenology and thermal properties constant over iterations
                              scalarLAI,                    & ! intent(in): one-sided leaf area index (m2 m-2)
                              scalarSAI,                    & ! intent(in): one-sided stem area index (m2 m-2)
                              scalarExposedLAI,             & ! intent(in): exposed leaf area index after burial by snow (m2 m-2)
                              scalarExposedSAI,             & ! intent(in): exposed stem area index after burial by snow (m2 m-2)
                              scalarGrowingSeasonIndex,     & ! intent(in): growing season index (0=off, 1=on)
                              scalarFoliageNitrogenFactor,  & ! intent(in): foliage nitrogen concentration (1.0 = saturated)
                              scalarBulkVolHeatCapVeg,      & ! intent(in): bulk volumetric heat capacity of vegetation (J m-3 K-1)

                              ! model diagnostic variables from the hydrology routines (intent in) -- liquid water fluxes at layer interfaces
                              iLayerLiqFluxSnow,            & ! intent(in): liquid flux at the interface of each snow layer (m s-1)
                              iLayerLiqFluxSoil,            & ! intent(in): liquid flux at the interface of each soil layer (m s-1)

                              ! model diagnostic variables (intent in) -- thermal properties constant over iterations
                              mLayerTcrit,                  & ! intent(in): critical soil temperature above which all water is unfrozen (K)
                              mLayerVolHtCapBulk,           & ! intent(in): volumetric heat capacity in each layer (J m-3 K-1) 
                              iLayerThermalC,               & ! intent(in): thermal conductivity at the interface of each layer (W m-1 K-1)

                              ! model diagnostic variables (output)
                              mLayerdTheta_dTk,             & ! intent(out): derivative in the freezing curve (K-1)
                              iLayerConductiveFlux,         & ! intent(out): conductive energy flux at layer interfaces at end of time step (W m-2)
                              iLayerAdvectiveFlux,          & ! intent(out): advective energy flux at layer interfaces at end of time step (W m-2)
                              iLayerInitNrgFlux,            & ! intent(out): energy flux at layer interfaces at the start of the time step (W m-2)
                              iLayerNrgFlux,                & ! intent(out): energy flux at layer interfaces at the end of the time step (W m-2)

                              ! input/output variables from heatTransf subroutine: canopy air space variables
                              scalarVP_CanopyAir,           & ! intent(inout): trial vapor pressure of the canopy air space (Pa)
                              scalarCanopyStabilityCorrection,& ! intent(inout): stability correction for the canopy (-)
                              scalarGroundStabilityCorrection,& ! intent(inout): stability correction for the ground surface (-)

                              ! output variables from heatTransf subroutine
                              scalarCanopyTranspiration,    & ! intent(out): canopy transpiration (kg m-2 s-1)
                              scalarCanopyEvaporation,      & ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
                              scalarGroundEvaporation,      & ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
                              scalarCanairTempDiff,         & ! intent(out): iteration increment for temperature of the canopy air space (K)
                              scalarCanopyTempDiff,         & ! intent(out): iteration increment for temperature of the vegetation canopy (K)
                              mLayerTempDiff,               & ! intent(out): iteration increment for temperature of the snow-soil system (K)
                              scalarCanairTempNew,          & ! intent(out): new temperature of the canopy air space (K)
                              scalarCanopyTempNew,          & ! intent(out): new temperature of the vegetation canopy (K)
                              scalarCanopyIceNew,           & ! intent(out): mass of ice on the canopy (kg m-2) 
                              scalarCanopyLiqNew,           & ! intent(out): mass of liquid water on the canopy (kg m-2)
                              mLayerTempNew,                & ! intent(out): new temperature of the snow-soil system (K)
                              mLayerVolFracIceNew,          & ! intent(out): new volumetric fraction of ice (-)
                              mLayerVolFracLiqNew,          & ! intent(out): new volumetric fraction of liquid water (-)
                              mLayerMatricHeadNew,          & ! intent(out): new matric head (m)
                              err,message)                    ! intent(out): error control
 ! ---------------------------------------------------------------------------------------------------------------------------------------------------------
 ! utility module
 USE vegNrgFlux_module,only:vegNrgFlux                        ! compute energy fluxes for vegetation and ground surface
 USE phseChange_module,only:phseChange                        ! compute change in phase over the time step
 USE snow_utils_module,only:fracliquid                        ! compute the fraction of liquid water at a given temperature (snow)
 USE snow_utils_module,only:templiquid                        ! compute the temperature at a given fraction of liquid water (snow)
 USE snow_utils_module,only:dFracLiq_dTk                      ! differentiate the freezing curve w.r.t. temperature (snow)
 USE soil_utils_module,only:dTheta_dTk                        ! differentiate the freezing curve w.r.t. temperature (soil)
 USE conv_funcs_module,only:relhm2sphm                        ! compute specific humidity 
 USE matrixSolv_module,only:matrixSolv                        ! solve full matrix
 USE tridagSolv_module,only:tridag                            ! solve tridiagonal system of equations
 ! ---------------------------------------------------------------------------------------------------------------------------------------------------------
 implicit none
 ! input variables from the heatTransf subroutine
 real(dp),intent(in)            :: dt                          ! time step (seconds)
 integer(i4b),intent(in)        :: iter                        ! iteration count
 logical(lgt),intent(in)        :: firstSubStep                ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(in)        :: computeVegFlux              ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 real(dp),intent(in)            :: scalarCanairTempIter        ! trial temperature of the canopy air space (K)
 real(dp),intent(in)            :: scalarCanopyTempIter        ! trial vegetation temperature (K)
 real(dp),intent(in)            :: scalarCanopyIceIter         ! trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(in)            :: scalarCanopyLiqIter         ! trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
 real(dp),intent(in)            :: mLayerTempIter(:)           ! trial temperature at the current iteration (K)
 real(dp),intent(in)            :: mLayerVolFracIceIter(:)     ! volumetric fraction of ice at the current iteration (-)
 real(dp),intent(in)            :: mLayerVolFracLiqIter(:)     ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),intent(in)            :: mLayerMatricHeadIter(:)     ! matric head at the current iteration (m)
 real(dp),intent(in)            :: canopyTempIncrOld           ! iteration increment for temperature of the vegetation canopy in the previous iteration (K)
 real(dp),intent(in)            :: mLayerTempIncrOld(:)        ! iteration increment for temperature of the snow-soil vector in the previous iteration (K)
 ! model decisions
 integer(i4b),intent(in)        :: num_method                  ! choice of numerical method
 integer(i4b),intent(in)        :: fDerivMeth                  ! method used to calculate flux derivatives
 integer(i4b),intent(in)        :: bcLowrTdyn                  ! type of upper boundary condition for thermodynamics
 ! model index variables
 integer(i4b),intent(in)        :: nLayers                     ! number of layers
 integer(i4b),intent(in)        :: layerType(:)                ! type of the layer (ix_soil or ix_snow)
 ! physical attributes
 integer(i4b),intent(in)        :: vegTypeIndex                ! vegetation type index
 integer(i4b),intent(in)        :: soilTypeIndex               ! soil type index
 ! general model parameters
 real(dp),intent(in)            :: wimplicit                   ! weight assigned to start-of-step fluxes (-)
 real(dp),intent(in)            :: snowfrz_scale               ! scaling parameter for the snow freezing curve (K-1)
 real(dp),intent(in)            :: lowerBoundTemp              ! temperature of the lower boundary (K)
 ! vegetation parameters
 real(dp),intent(in)            :: heightCanopyTop             ! height of top of the vegetation canopy above ground surface (m)
 real(dp),intent(in)            :: heightCanopyBottom          ! height of bottom of the vegetation canopy above ground surface (m)
 ! soil parameters
 real(dp),intent(in)            :: vGn_alpha                   ! van Genutchen "alpha" parameter
 real(dp),intent(in)            :: vGn_n                       ! van Genutchen "n" parameter
 real(dp),intent(in)            :: theta_sat                   ! soil porosity (-)
 real(dp),intent(in)            :: theta_res                   ! soil residual volumetric water content (-)
 real(dp),intent(in)            :: vGn_m                       ! van Genutchen "m" parameter (-)
 real(dp),intent(in)            :: kappa                       ! constant used to express matric head as a function of temperature (m K-1)
 ! model state variables
 ! NOTE: protected with the intent(in) attribute
 real(dp),intent(in)            :: scalarCanairTemp            ! temperature of the canopy air space (K)
 real(dp),intent(in)            :: scalarCanopyTemp            ! temperature of the vegetation canopy (K)
 real(dp),intent(in)            :: scalarCanopyIce             ! mass of ice on the vegetation canopy (kg m-2)
 real(dp),intent(in)            :: scalarCanopyLiq             ! mass of liquid water on the vegetation canopy (kg m-2)
 real(dp),intent(in)            :: scalarSnowDepth             ! snow depth on the ground surface (m)
 real(dp),intent(in)            :: mLayerTemp(:)               ! temperature of each layer (K)
 real(dp),intent(in)            :: mLayerVolFracIce(:)         ! volumetric fraction of ice in each layer (-)
 real(dp),intent(in)            :: mLayerVolFracLiq(:)         ! volumetric fraction of liquid water in each layer (-)
 ! model coordinate variables (intent in)
 real(dp),intent(in)            :: mLayerDepth(:)              ! depth of each layer (m)
 real(dp),intent(in)            :: mLayerHeight(:)             ! height at the mid-point of each layer (m)
 ! model diagnostic variables for vegetation (intent in) -- phenology and thermal properties constant over iterations
 real(dp),intent(in)            :: scalarLAI                   ! one-sided leaf area index (m2 m-2)
 real(dp),intent(in)            :: scalarSAI                   ! one-sided stem area index (m2 m-2)
 real(dp),intent(in)            :: scalarExposedLAI            ! exposed leaf area index after burial by snow (m2 m-2)
 real(dp),intent(in)            :: scalarExposedSAI            ! exposed stem area index after burial by snow (m2 m-2)
 real(dp),intent(in)            :: scalarGrowingSeasonIndex    ! growing season index (0=off, 1=on)
 real(dp),intent(in)            :: scalarFoliageNitrogenFactor ! foliage nitrogen concentration (1.0 = saturated)
 real(dp),intent(in)            :: scalarBulkVolHeatCapVeg     ! bulk volumetric heat capacity of vegetation (J m-3 K-1)
 ! model diagnostic variables from the hydrology routines (intent in) -- liquid water fluxes at layer interfaces
 real(dp),intent(in)            :: iLayerLiqFluxSnow(0:)       ! liquid flux at the interface of each snow layer (m s-1)
 real(dp),intent(in)            :: iLayerLiqFluxSoil(0:)       ! liquid flux at the interface of each soil layer (m s-1)
 ! model diagnostic variables (intent in) -- thermal properties constant over iterations
 real(dp),intent(in)            :: mLayerTcrit(:)              ! critical soil temperature above which all water is unfrozen (K)
 real(dp),intent(in)            :: mLayerVolHtCapBulk(:)       ! volumetric heat capacity in each layer (J m-3 K-1) 
 real(dp),intent(in)            :: iLayerThermalC(0:)          ! thermal conductivity at the interface of each layer (W m-1 K-1)
 ! model diagnostic variables (intent out)
 real(dp),intent(out)           :: mLayerdTheta_dTk(:)         ! derivative in the freezing curve (K-1)
 real(dp),intent(out)           :: iLayerConductiveFlux(0:)    ! conductive energy flux at layer interfaces at end of time step (W m-2)
 real(dp),intent(out)           :: iLayerAdvectiveFlux(0:)     ! advective energy flux at layer interfaces at end of time step (W m-2)
 real(dp),intent(out)           :: iLayerInitNrgFlux(0:)       ! energy flux at layer interfaces at the start of the time step (W m-2)
 real(dp),intent(out)           :: iLayerNrgFlux(0:)           ! energy flux at layer interfaces at the end of the time step (W m-2)
 ! input/output variables from the heatTransf subroutine: canopy air space variables
 real(dp),intent(inout)         :: scalarVP_CanopyAir          ! trial vapor pressure of the canopy air space (Pa)
 real(dp),intent(inout)         :: scalarCanopyStabilityCorrection ! stability correction for the canopy (-)
 real(dp),intent(inout)         :: scalarGroundStabilityCorrection ! stability correction for the ground surface (-)
 ! output variables from the heatTransf subroutine
 real(dp),intent(out)           :: scalarCanopyTranspiration   ! canopy transpiration (kg m-2 s-1)
 real(dp),intent(out)           :: scalarCanopyEvaporation     ! canopy evaporation/condensation (kg m-2 s-1)
 real(dp),intent(out)           :: scalarGroundEvaporation     ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
 real(dp),intent(out)           :: scalarCanairTempDiff        ! iteration increment for temperature of the canopy air space (K)
 real(dp),intent(out)           :: scalarCanopyTempDiff        ! iteration increment for temperature of the vegetation canopy (K)
 real(dp),intent(out)           :: mLayerTempDiff(:)           ! iteration increment for temperature of the snow-soil system (K) 
 real(dp),intent(out)           :: scalarCanairTempNew         ! new temperature of the canopy air space (K)
 real(dp),intent(out)           :: scalarCanopyTempNew         ! new temperature of the vegetation canopy (K)
 real(dp),intent(out)           :: scalarCanopyIceNew          ! mass of ice on the canopy (kg m-2) 
 real(dp),intent(out)           :: scalarCanopyLiqNew          ! mass of liquid water on the canopy (kg m-2)
 real(dp),intent(out)           :: mLayerTempNew(:)            ! new temperature of each model layer (K)
 real(dp),intent(out)           :: mLayerVolFracIceNew(:)      ! new volumetric fraction of ice (-)
 real(dp),intent(out)           :: mLayerVolFracLiqNew(:)      ! new volumetric fraction of liquid water (-)
 real(dp),intent(out)           :: mLayerMatricHeadNew(:)      ! new matric head (m)
 integer(i4b),intent(out)       :: err                         ! error code
 character(*),intent(out)       :: message                     ! error message
 ! ---------------------------------------------------------------------------------------------------------------------------------
 ! define local variables
 ! ---------------------------------------------------------------------------------------------------------------------------------
 ! define general local variables
 character(LEN=256)             :: cmessage                   ! error message of downwind routine
 integer(i4b)                   :: iLayer,jLayer              ! index of model layers
 logical(lgt)                   :: printflag                  ! .true. if print progress to the screen
 logical(lgt)                   :: fTranspire                 ! .true. if computing transpiration
 logical(lgt),parameter         :: computeJacobian=.false.    ! .true. if desire to compute the Jacobian matrix (just used for testing)
 real(dp)                       :: canopyDepth                ! depth of the vegetation canopy (m)
 real(dp)                       :: theta                      ! total volumetric water content (liquid plus ice)
 real(dp)                       :: critDiff                   ! temperature difference from critical temperature (K)
 real(dp)                       :: maxdiffTemp(1)             ! maximum difference between temperature input and start-of-step temperature (K)
 real(dp),parameter             :: epsT=1.d-10                ! offset from Tcrit when re-setting iterations at the critical temperature (K)
 integer(i4b)                   :: nUnsat                     ! number of unsaturated layers
 real(dp),dimension(1)          :: amaxIncrement              ! maximum iteration increment
 ! define local variables for the fluxes at vegetation and ground surfaces
 real(dp)                       :: saveTemp_CanopyAir         ! trial temperature of the canopy air space (K)
 real(dp)                       :: saveVP_CanopyAir           ! trial vapor pressure of the canopy air space (Pa)
 real(dp)                       :: saveCanopyStabilityCorrection ! stability correction for the canopy (-)
 real(dp)                       :: saveGroundStabilityCorrection ! stability correction for the ground surface (-)
 real(dp)                       :: canairNetFlux              ! net energy flux for the canopy air space (W m-2)
 real(dp)                       :: canopyNetFlux              ! net energy flux for the vegetation canopy (W m-2)
 real(dp)                       :: groundNetFlux              ! net energy flux for the ground surface (W m-2)
 real(dp)                       :: dCanairNetFlux_dCanairTemp ! derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp)                       :: dCanairNetFlux_dCanopyTemp ! derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp)                       :: dCanairNetFlux_dGroundTemp ! derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
 real(dp)                       :: dCanopyNetFlux_dCanairTemp ! derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp)                       :: dCanopyNetFlux_dCanopyTemp ! derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp)                       :: dCanopyNetFlux_dGroundTemp ! derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
 real(dp)                       :: dGroundNetFlux_dCanairTemp ! derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp)                       :: dGroundNetFlux_dCanopyTemp ! derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp)                       :: dGroundNetFlux_dGroundTemp ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
 ! define fluxes at the start of the sub-step
 real(dp),save                  :: canairNetFluxInit          ! net energy flux at the canopy air space at the start of the substep (W m-2)
 real(dp),save                  :: canopyNetFluxInit          ! net energy flux at the canopy at the start of the substep (W m-2)
 real(dp),dimension(0:nLayers)  :: iLayerNrgFluxInit          ! flux at layer interfaces of the snow-soil system at the start of the substep (W m-2)
 ! define the Jacobian matrices
 integer(i4b)                   :: nState                     ! number of state variables
 integer(i4b),parameter         :: ixCas=1                    ! index of layer for the canopy air space
 integer(i4b),parameter         :: ixVeg=2                    ! index of layer for the vegetation canopy
 integer(i4b),parameter         :: ixSfc=3                    ! index of layer for the top snow or soil layer
 real(dp),allocatable           :: aJac(:,:)                  ! Jacobian matrix
 real(dp),allocatable           :: rVec(:)                    ! residual vector
 real(dp),allocatable           :: xInc(:)                    ! iteration increment
 !real(dp),dimension(nLayers+2,nLayers+2) :: aJac              ! Jacobian matrix
 !real(dp),dimension(nLayers+2)           :: rVec              ! residual vector
 !real(dp),dimension(nLayers+2)           :: xInc              ! iteration increment
 ! define the local variables for the solution
 real(dp)                       :: dTheta_dTkCanopy           ! derivative in fraction liquid water w.r.t. canopy temperature (K-1)
 real(dp),dimension(0:nLayers)  :: dFlux_dTempAbove           ! derivative in flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
 real(dp),dimension(0:nLayers)  :: dFlux_dTempBelow           ! derivative in flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
 real(dp)                       :: nrg0,nrg1                  ! energy content at the start of the time step / current iteration (J m-3)
 real(dp)                       :: flx0,flx1                  ! fluxes at the start of the time step / current iteration (J m-3 s-1)
 real(dp)                       :: phse                       ! phase change term (J m-3)
 real(dp)                       :: wtim                       ! weighted time (s-1)
 real(dp)                       :: canairResidual             ! residual for the energy of the canopy air space (J m-3)
 real(dp)                       :: canopyResidual             ! residual for the energy of the vegetation canopy (J m-3)
 real(dp),dimension(nLayers)    :: mLayerResidual             ! residual for the energy of the snow/soil layers (J m-3)
 ! define tri-diagonal matrix elements for the canopy
 real(dp)                       :: diagCanopy                 ! diagonal element for the canopy
 real(dp)                       :: d_m1Canopy                 ! sub-diagonal element for the derivative in net ground flux w.r.t. canopy temperature (J m-3 K-1)
 real(dp)                       :: d_p1Canopy                 ! super-diagonal element for the derivative in net canopy flux w.r.t. ground temperature (J m-3 K-1)
 ! define tri-diagonal matrix elements for the snow-soil system
 real(dp),dimension(nLayers)    :: diag                       ! diagonal (J m-3 K-1)
 real(dp),dimension(nLayers-1)  :: d_m1                       ! sub-diagonal (J m-3 K-1)
 real(dp),dimension(nLayers-1)  :: d_p1                       ! super-diagonal (J m-3 K-1)
 ! define local variables for phase change
 integer(i4b)                   :: iSnow                      ! index of snow layer
 real(dp)                       :: fLiq                       ! fraction of liquid water on the vegetation canopy (-)
 real(dp)                       :: tWat                       ! total water on the vegetation canopy (kg m-2) 
 real(dp)                       :: aflx                       ! average flux from start and end of time step (J m-3 s-1)
 real(dp)                       :: xres                       ! trial residual (J m-3)
 real(dp)                       :: delT                       ! temperature increment (K)
 real(dp)                       :: dTheta_dT                  ! derivative in volume fraction of total water w.r.t. temperature (K-1) 
 real(dp)                       :: tempTrial                  ! trial temperature value (K)
 real(dp)                       :: tempVolFracLiq             ! temporary value of volumetric fraction of liquid water (-)
 real(dp)                       :: tempVolFracIce             ! temporary value of volumetric fraction of ice (-)
 real(dp)                       :: xmin,xmax                  ! bounds for fraction of liquid water (used in bi-section) 
 real(dp),parameter             :: nrgToler=1._dp             ! energy tolerance (J m-3)
 integer(i4b)                   :: jiter                      ! iteration index when attempting to correct phase change
 integer(i4b),parameter         :: maxiter=50                 ! maximum number of iterations used in phase change correction
 ! ---------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="heatTransf_muster/"

 ! initialize print flag
 printflag=.false.

 ! iterate on the surface temperature
 if(num_method==iterSurfEnergyBal)then
  err=20; message=trim(message)//'option "iterSurfEnergyBal" not implemented yet';return
 endif

 ! define the canopy depth (m)
 canopyDepth = heightCanopyTop - heightCanopyBottom
 if(heightCanopyBottom > heightCanopyTop)then
  err=20; message=trim(message)//'height of the bottom of the canopy > top of the canopy'; return
 endif

 ! define the number of state variables
 if(computeVegFlux)then
  nState = nLayers+2 ! +2 accounts for (1) canopy air space; and (2) vegetation canopy
 else
  nState = nLayers
 endif

 ! get an initial value for the canopy air space variables, to be used in the Jacobian calculations
 ! NOTE: canopy air space values are intent(inout), which muck up derivative calculations if used directly
 saveTemp_CanopyAir = scalarCanairTempIter       ! trial temperature of the canopy air space (K)
 saveVP_CanopyAir   = scalarVP_CanopyAir         ! trial vapor pressure of the canopy air space (Pa)
 saveCanopyStabilityCorrection = scalarCanopyStabilityCorrection    ! stability correction for the canopy (-)
 saveGroundStabilityCorrection = scalarGroundStabilityCorrection    ! stability correction for the ground surface (-)

 ! allocate space for the Jacobian matric
 if(computeVegFlux)then
  allocate(aJac(nLayers+2,nLayers+2),rVec(nLayers+2),xInc(nLayers+2),stat=err)
  if(err/=0)then; err=20; message=trim(message)//'problem allocating space for the Jacobian matrix'; return; endif
 endif

 ! --------------------------------------------------------------------------------------------------------------------------------

 ! --------------------------------------------------------------------------------------------------------------------------------
 ! * COMPUTE FLUXES...
 ! --------------------------------------------------------------------------------------------------------------------------------

 ! check temperatures
 if(computeVegFlux)then
  if(scalarCanopyTempIter < 200._dp)then; message=trim(message)//'canopy temperature is very cold'; err=20; return; endif
 endif
 if(mLayerTempIter(1)    < 200._dp)then; message=trim(message)//'ground temperature is very cold'; err=20; return; endif
 !if(mLayerTempIter(1) < 260._dp)then
 ! print*, 'before flux computations: mLayerTempIter(1) = ', mLayerTempIter(1)
 !endif

 ! ***** compute energy fluxes at vegetation and ground surfaces
 call vegNrgFlux(&
                 ! input
                 dt,                              & ! intent(in): time step (seconds)
                 iter,                            & ! intent(in): iteration index
                 firstSubStep,                    & ! intent(in): flag to indicate if we are processing the first sub-step
                 computeVegFlux,                  & ! intent(in): flag to indicate if we need to compute fluxes over vegetation
                 (iter==1),                       & ! intent(in): flag to indicate if we need to compute shortwave radiation
                 scalarCanairTempIter,            & ! intent(in): trial temperature of the canopy air space (K)
                 scalarCanopyTempIter,            & ! intent(in): trial value of canopy temperature (K)
                 mLayerTempIter(1),               & ! intent(in): trial value of ground temperature (K)
                 scalarCanopyIceIter,             & ! intent(in): trial mass of ice on the vegetation canopy (kg m-2)
                 scalarCanopyLiqIter,             & ! intent(in): trial mass of liquid water on the vegetation canopy (kg m-2)
                 vegTypeIndex,                    & ! intent(in): vegetation type index
                 soilTypeIndex,                   & ! intent(in): soil type index
                 scalarLAI,                       & ! intent(in): one-sided leaf area index (m2 m-2)
                 scalarSAI,                       & ! intent(in): one-sided stem area index (m2 m-2)
                 scalarExposedLAI,                & ! intent(in): exposed leaf area index after burial by snow (m2 m-2)
                 scalarExposedSAI,                & ! intent(in): exposed stem area index after burial by snow (m2 m-2)
                 scalarGrowingSeasonIndex,        & ! intent(in): growing season index (0=off, 1=on)
                 scalarFoliageNitrogenFactor,     & ! intent(in): foliage nitrogen concentration (1.0 = saturated)
                 ! input/output: canopy air space variables
                 scalarVP_CanopyAir,              & ! intent(inout): trial vapor pressure of the canopy air space (Pa)
                 scalarCanopyStabilityCorrection, & ! intent(inout): stability correction for the canopy (-)
                 scalarGroundStabilityCorrection, & ! intent(inout): stability correction for the ground surface (-)
                 ! output: liquid water fluxes associated with evaporation/transpiration
                 scalarCanopyTranspiration,       & ! intent(out): canopy transpiration (kg m-2 s-1)
                 scalarCanopyEvaporation,         & ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
                 scalarGroundEvaporation,         & ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
                 ! output: fluxes
                 canairNetFlux,                   & ! intent(out): net energy flux for the canopy air space (W m-2)
                 canopyNetFlux,                   & ! intent(out): net energy flux for the vegetation canopy (W m-2)
                 groundNetFlux,                   & ! intent(out): net energy flux for the ground surface (W m-2)
                 ! output: flux derivatives
                 dCanairNetFlux_dCanairTemp,      & ! intent(out): derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
                 dCanairNetFlux_dCanopyTemp,      & ! intent(out): derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
                 dCanairNetFlux_dGroundTemp,      & ! intent(out): derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
                 dCanopyNetFlux_dCanairTemp,      & ! intent(out): derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
                 dCanopyNetFlux_dCanopyTemp,      & ! intent(out): derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
                 dCanopyNetFlux_dGroundTemp,      & ! intent(out): derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
                 dGroundNetFlux_dCanairTemp,      & ! intent(out): derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
                 dGroundNetFlux_dCanopyTemp,      & ! intent(out): derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
                 dGroundNetFlux_dGroundTemp,      & ! intent(out): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
                 ! output: error control
                 err,cmessage)                      ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 !write(*,'(a,1x,2(f20.8,1x))') 'scalarCanopyStabilityCorrection, scalarGroundStabilityCorrection = ', &
 !                               scalarCanopyStabilityCorrection, scalarGroundStabilityCorrection

 ! ***** compute fluxes at layer interfaces and their derivatives (J m-2 s-1)
 call iLayer_nrg(&
                 ! (model control variables)
                 fDerivMeth,                    & ! intent(in): method used to compute derivatives (numerical or analytical)
                 ! (input)
                 layerType,                     & ! intent(in): type of each layer
                 mLayerDepth,                   & ! intent(in): depth of each layer (m)
                 mLayerHeight,                  & ! intent(in): height of layer mid-points (m)
                 iLayerThermalC,                & ! intent(in): thermal conductivity at layer interfaces (W m-1)
                 iLayerLiqFluxSnow,             & ! intent(in): liquid flux at the interface of each snow layer (m s-1)
                 iLayerLiqFluxSoil,             & ! intent(in): liquid flux at the interface of each soil layer (m s-1)
                 mLayerTempIter,                & ! intent(in): trial temperature at the current iteration (K)
                 groundNetFlux,                 & ! intent(in): total flux at the ground surface (W m-2)
                 dGroundNetFlux_dGroundTemp,    & ! intent(in): derivative in total ground surface flux w.r.t. ground temperature (W m-2 K-1)
                 lowerBoundTemp,                & ! intent(in): temperature of the lower boundary (K)
                 ! (output)
                 iLayerConductiveFlux,          & ! intent(out): conductive energy flux at layer interfaces (W m-2)
                 iLayerAdvectiveFlux,           & ! intent(out): advective energy flux at layer interfaces (W m-2)
                 iLayerNrgFlux,                 & ! intent(out): energy flux at the layer interfaces (W m-2)
                 dFlux_dTempAbove,              & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
                 dFlux_dTempBelow,              & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
                 err,cmessage)                    ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! --------------------------------------------------------------------------------------------------------------------------------


 ! --------------------------------------------------------------------------------------------------------------------------------
 ! * COMPUTE RESIDUAL VECTOR...
 ! --------------------------------------------------------------------------------------------------------------------------------

 ! ***** assign initial fluxes
 if(iter==1)then
  ! check that the temperature matches the temperature at the start of the step
  maxdiffTemp = maxval(abs(mLayerTempIter - mLayerTemp))
  if(maxdiffTemp(1) > 1.d-8)then; err=20; message=trim(message)//'first guess for temperature must match start-of-step value'; return; endif
  ! assign initial fluxes
  canairNetFluxInit = canairNetFlux
  canopyNetFluxInit = canopyNetFlux
  iLayerNrgFluxInit = iLayerNrgFlux
 endif

 ! ***** compute the residual for the vegetation canopy and the canopy air space
 if(computeVegFlux)then

  ! ***** compute the residual for the vegetation canopy
  ! (compute individual terms) 
  nrg0 = scalarBulkVolHeatCapVeg*scalarCanopyTemp                      ! energy content at the start of the time step (J m-3)
  nrg1 = scalarBulkVolHeatCapVeg*scalarCanopyTempIter                  ! energy content at the current iteration (J m-3)
  flx0 = canopyNetFluxInit/canopyDepth                                 ! flux at the start of the time step (J m-3 s-1)
  flx1 = canopyNetFlux/canopyDepth                                     ! flux at the current iteration (J m-3 s-1)
  phse = LH_fus*(scalarCanopyIceIter - scalarCanopyIce)/canopyDepth    ! phase change term (J m-3)
  ! (compute residuals)
  canopyResidual = nrg1 - (nrg0 + (flx0*wimplicit + flx1*(1._dp - wimplicit))*dt + phse)
  !print*, 'nrg1, nrg0, flx1, canopyResidual = ', nrg1, nrg0, flx1, canopyResidual
  !print*, 'scalarCanopyIceIter, scalarCanopyIce = ', scalarCanopyIceIter, scalarCanopyIce
  !print*, 'nrg0, nrg1, nrg1 - nrg0 = ', nrg0, nrg1, nrg1 - nrg0
  !print*, 'wimplicit, scalarBulkVolHeatCapVeg, canopyResidual = ', wimplicit, scalarBulkVolHeatCapVeg, canopyResidual
  !write(*,'(a,1x,10(f20.10,1x))') 'canopy: scalarCanopyTempIter, scalarCanopyTemp, canopyNetFluxInit, canopyNetFlux, canopyDepth, flx1*dt, phse = ', &
  !                                         scalarCanopyTempIter, scalarCanopyTemp, canopyNetFluxInit, canopyNetFlux, canopyDepth, flx1*dt, phse

  ! ***** compute the residual for the canopy air space
  ! (compute individual terms)
  nrg0 = Cp_air*iden_air*scalarCanairTemp                              ! energy content at the start of the time step (J m-3)
  nrg1 = Cp_air*iden_air*scalarCanairTempIter                          ! energy content at the current iteration (J m-3)
  flx0 = canairNetFluxInit/canopyDepth                                 ! flux at the start of the time step (J m-3 s-1)
  flx1 = canairNetFlux/canopyDepth                                     ! flux at the current iteration (J m-3 s-1)
  ! (compute residuals)
  canairResidual = nrg1 - (nrg0 + (flx0*wimplicit + flx1*(1._dp - wimplicit))*dt)
  !print*, 'nrg1, nrg0, flx1, canairResidual = ', nrg1, nrg0, flx1, canairResidual

 endif

 ! ***** compute the residual vector for all snow/soil layers (J m-3)
 do iLayer=1,nLayers
  ! (compute energy storage)
  nrg0 = mLayerVolHtCapBulk(iLayer)*mLayerTemp(iLayer)                                        ! energy content at the start of the time step (J m-3)
  nrg1 = mLayerVolHtCapBulk(iLayer)*mLayerTempIter(iLayer)                                    ! energy content at the current iteration (J m-3)
  ! (compute energy fluxes)
  flx0 = -(iLayerInitNrgFlux(iLayer) - iLayerInitNrgFlux(iLayer-1))/mLayerDepth(iLayer)       ! flux at the start of the time step (J m-3 s-1)
  flx1 = -(iLayerNrgFlux(iLayer) - iLayerNrgFlux(iLayer-1))/mLayerDepth(iLayer)               ! flux at the current iteration (J m-3 s-1)
  ! (compute phase change)
  phse = LH_fus*iden_ice*(mLayerVolFracIceIter(iLayer) - mLayerVolFracIce(iLayer))            ! phase change term (J m-3)
  ! (compute residuals)
  mLayerResidual(iLayer) = nrg1 - (nrg0 + (flx0*wimplicit + flx1*(1._dp - wimplicit))*dt + phse)
  ! (print progress)
  !if(iLayer==1) write(*,'(a)')                                'iLayer, mLayerResidual(iLayer), mLayerTempIter(iLayer), mLayerVolFracLiqIter(iLayer), mLayerVolFracIceIter(iLayer), iLayerNrgFlux(iLayer-1), iLayerNrgFlux(iLayer), flx1*dt, phse = '
  !write(*,'(i4,1x,e20.10,5x,f13.9,1x,2(f9.5,1x),5(e14.4,1x))') iLayer, mLayerResidual(iLayer), mLayerTempIter(iLayer), mLayerVolFracLiqIter(iLayer), mLayerVolFracIceIter(iLayer), iLayerNrgFlux(iLayer-1), iLayerNrgFlux(iLayer), flx1*dt, phse
 end do


 ! --------------------------------------------------------------------------------------------------------------------------------
 ! * COMPUTE TRIDIAGONAL MATRIX...
 ! --------------------------------------------------------------------------------------------------------------------------------

 ! ***** compute the derivative in the freezing curve w.r.t. temperature (K-1)
 ! * vegetation canopy
 if(computeVegFlux)then
  if(scalarCanopyIceIter > 0._dp)then
   theta = (scalarCanopyIceIter + scalarCanopyLiqIter)/(canopyDepth*iden_water)
   dTheta_dTkCanopy = dFracLiq_dTk(scalarCanopyTempIter,snowfrz_scale)*theta
  else
   dTheta_dTkCanopy = 0._dp
  endif
 else
  dTheta_dTkCanopy = valueMissing
 endif
 ! * all snow-soil layers
 do iLayer=1,nLayers
  select case(layerType(iLayer))
   case(ix_snow) ! (snow layers)
    theta = mLayerVolFracIceIter(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiqIter(iLayer)
    mLayerdTheta_dTk(iLayer) = dFracLiq_dTk(mLayerTempIter(iLayer),snowfrz_scale)*theta
   case(ix_soil) ! (soil layers)
    if(mLayerVolFracIceIter(iLayer)>0._dp)then
     mLayerdTheta_dTk(iLayer) = dTheta_dTk(mLayerTempIter(iLayer),theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m)
    else
     mLayerdTheta_dTk(iLayer) = 0._dp
    endif
   case default; err=40; message=trim(message)//"cannot identify the layer as snow or soil"; return
  endselect
  !write(*,'(a,1x,i4,1x,f20.6)') 'mLayerdTheta_dTk(iLayer)*LH_fus*iden_water = ', iLayer, mLayerdTheta_dTk(iLayer)*LH_fus*iden_water
 end do

 ! compute the weighted time for end-of-step values
 wtim = (1._dp - wimplicit)*dt  ! weighted time

 ! ***** assemble the tri-diagonal matrix for the snow-soil layers
 if(.not.computeVegFlux)then
  diag = (wtim/mLayerDepth)*(-dFlux_dTempBelow(0:nLayers-1) + dFlux_dTempAbove(1:nLayers)) + mLayerVolHtCapBulk + mLayerdTheta_dTk*LH_fus*iden_water
  d_m1 = (wtim/mLayerDepth(2:nLayers  ))*(-dFlux_dTempAbove(1:nLayers-1) )
  d_p1 = (wtim/mLayerDepth(1:nLayers-1))*( dFlux_dTempBelow(1:nLayers-1) )
 endif

 ! --------------------------------------------------------------------------------------------------------------------------------
 ! * COMPUTE JACOBIAN MATRIX...
 ! --------------------------------------------------------------------------------------------------------------------------------

 ! case of exposed vegetation
 if(computeVegFlux)then

  ! ***** assemble the Jacobian matrix for the full system
  ! initialize the Jacobian as all zeros
  aJac(1:nState,1:nState) = 0._dp
  ! define Jacobian matrix for the canopy air space (J m-3 K-1)
  aJac(ixCas,ixCas) = (wtim/canopyDepth)*(-dCanairNetFlux_dCanairTemp) + Cp_air*iden_air
  aJac(ixCas,ixVeg) = (wtim/canopyDepth)*(-dCanairNetFlux_dCanopyTemp)
  aJac(ixCas,ixSfc) = (wtim/canopyDepth)*(-dCanairNetFlux_dGroundTemp)
  ! define Jacobian matrix for the vegetation canopy (J m-3 K-1)
  aJac(ixVeg,ixCas) = (wtim/canopyDepth)*(-dCanopyNetFlux_dCanairTemp)
  aJac(ixVeg,ixVeg) = (wtim/canopyDepth)*(-dCanopyNetFlux_dCanopyTemp) + scalarBulkVolHeatCapVeg + dTheta_dTkCanopy*LH_fus*iden_water
  aJac(ixVeg,ixSfc) = (wtim/canopyDepth)*(-dCanopyNetFlux_dGroundTemp)
  ! define Jacobian matric for the surface (J m-3 K-1)
  aJac(ixSfc,ixCas) = (wtim/mLayerDepth(1))*(-dGroundNetFlux_dCanairTemp)
  aJac(ixSfc,ixVeg) = (wtim/mLayerDepth(1))*(-dGroundNetFlux_dCanopyTemp)
  ! define Jacobian matrix for the snow-soil system
  do iLayer=1,nLayers  ! loop through layers in the snow-soil system
   jLayer = iLayer+2
   aJac(jLayer,jLayer)   = (wtim/mLayerDepth(iLayer))*(-dFlux_dTempBelow(iLayer-1) + dFlux_dTempAbove(iLayer)) + mLayerVolHtCapBulk(iLayer) + mLayerdTheta_dTk(iLayer)*LH_fus*iden_water
   if(iLayer > 1)       aJac(jLayer-1,jLayer) = (wtim/mLayerDepth(iLayer-1))*( dFlux_dTempBelow(iLayer-1) )
   if(iLayer < nLayers) aJac(jLayer+1,jLayer) = (wtim/mLayerDepth(iLayer+1))*(-dFlux_dTempAbove(iLayer  ) )
  end do  ! (looping through layers in the snow-soil system)

  ! compute Jacobian matrix (just used for testing)
  if(computeJacobian)then
   call cmpJacobian(&
                    dt,                           & ! intent(in): time step (seconds)
                    scalarCanairTempIter,         & ! intent(in): trial value of canopy air space temperature (K)
                    scalarCanopyTempIter,         & ! intent(in): trial value of canopy temperature (K)
                    mLayerTempIter,               & ! intent(in): trial temperature at the current iteration (K)
                    err,message)                    ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

 endif   ! (if computing veg flux)

 ! --------------------------------------------------------------------------------------------------------------------------------
 ! * COMPUTE ITERATION INCREMENT...
 ! --------------------------------------------------------------------------------------------------------------------------------

 ! ***** solve the tridiagonal system of equations -- returns mLayerTempDiff
 if(.not.computeVegFlux)then
  call tridag(d_m1,                    & ! intent(in): sub-diagonal elements of the tridiagonal system (J m-3 K-1)
              diag,                    & ! intent(in): diagonal elements of the tridiagonal system (J m-3 K-1)
              d_p1,                    & ! intent(in): super-diagonal elements of the tridiagonal system (J m-3 K-1)
              -mLayerResidual,         & ! intent(in): residual vector (J m-3)
              mLayerTempDiff,          & ! intent(out): temperature increment (K)
              err,cmessage)              ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! ensure that the canopy temperature increments are zero
  scalarCanairTempDiff = 0._dp
  scalarCanopyTempDiff = 0._dp
 endif

 ! ***** solve the full matrix
 if(computeVegFlux)then
  rVec = (/canairResidual,canopyResidual,mLayerResidual/)
  call matrixSolv(aJac,-rVec,xInc,err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! save iteration increments
  scalarCanairTempDiff = xInc(ixCas)
  scalarCanopyTempDiff = xInc(ixVeg)
  mLayerTempDiff(1:nLayers) = xInc(1+ixVeg:nLayers+ixVeg)
 endif

 ! adjust iteration increments in cases where iteration increments are too large
 if(computeVegFlux)then
  if(abs(scalarCanairTempDiff) > 1._dp .or. abs(scalarCanopyTempDiff) > 1._dp .or. any(abs(mLayerTempDiff(1:nLayers)) > 1._dp) )then
   amaxIncrement = maxval(abs((/scalarCanairTempDiff,scalarCanopyTempDiff,mLayerTempDiff(1:nLayers)/)))
   !print*, 'scalarCanopyTempDiff = ', scalarCanopyTempDiff
   !print*, 'mLayerTempDiff(1:nLayers) = ', mLayerTempDiff(1:nLayers)
   scalarCanairTempDiff      = scalarCanairTempDiff/amaxIncrement(1)
   scalarCanopyTempDiff      = scalarCanopyTempDiff/amaxIncrement(1)
   mLayerTempDiff(1:nLayers) = mLayerTempDiff(1:nLayers)/amaxIncrement(1)
   !print*, 'scalarCanopyTempDiff = ', scalarCanopyTempDiff
   !print*, 'mLayerTempDiff(1:nLayers) = ', mLayerTempDiff(1:nLayers)
   !pause ' excessive increment'
  endif
 endif

 ! adjust iteration increment near the freezing point
 ! (use simplified bisection to take smaller steps near freezing)
 if(computeVegFlux)then
  if(scalarCanopyIceIter > 0.01_dp .and. scalarCanopyTempIter + scalarCanopyTempDiff > Tfreeze)then
   scalarCanopyTempDiff = (Tfreeze - scalarCanopyTempIter)*0.5_dp  ! go halfway to the freezing point
  endif
  ! (get the difference from freezing point (K)
  critDiff = Tfreeze - scalarCanopyTempIter
  ! (set temperature close to freezing point when it crosses freezing)
  if(critDiff > 0._dp)then; if(scalarCanopyTempDiff > critDiff) scalarCanopyTempDiff = critDiff + 0.0001_dp  ! below freezing crossing zero --> slightly above freezing
                      else; if(scalarCanopyTempDiff < critDiff) scalarCanopyTempDiff = critDiff - 0.0001_dp  ! above freezing crossing zero --> slightly below freezing
  endif
 endif

 ! adjust iteration increments in cases where iterations are oscillating
 if(iter > 5)then
  if(scalarCanopyTempDiff*canopyTempIncrOld < -0.01_dp .or. any(mLayerTempIncrOld(1:nLayers)*mLayerTempDiff(1:nLayers) < -0.01_dp) )then
   !write(*,'(a,1x,20(f20.10,1x))') 'temperature difference (old) = ', scalarCanopyTempDiff, mLayerTempDiff
   !write(*,'(a,1x,20(f20.10,1x))') 'temperature difference (new) = ', canopyTempIncrOld,    mLayerTempIncrOld
   scalarCanopyTempDiff      = 0.5_dp*scalarCanopyTempDiff
   mLayerTempDiff(1:nLayers) = 0.5_dp*mLayerTempDiff(1:nLayers)
   !pause ' iteration either oscillating or iteration increment is too large: cut iteration increment in half' 
  endif
 endif

 ! adjust del temperature for snow
 if(nSnow>0)then
  do iLayer=1,nSnow
   ! adjust del temperature in cases where snow temperature exceeds Tfreeze -- use bi-section
   if(mLayerTempIter(iLayer) + mLayerTempDiff(iLayer) > Tfreeze)then
    mLayerTempDiff(iLayer) = (Tfreeze-mLayerTempIter(iLayer))*0.5_dp
   endif
   ! check that temperature increment is not too large
   if(abs(mLayerTempDiff(iLayer)) > 10._dp)then; err=-20; message=trim(message)//'temperature increment is > 10K'; return; endif
  end do
 endif

 ! adjust del temperature in cases where soil temperature crosses the critical temperature
 do iLayer=nSnow+1,nLayers
  ! get the difference from the critical temperature (K)
  critDiff = mLayerTcrit(iLayer-nSnow) - mLayerTempIter(iLayer)
  ! set temperature to Tcrit in cases where temperatures cross Tcrit
  if(critDiff > 0._dp)then  ! (mLayerTempIter < Tcrit)
   if(mLayerTempDiff(iLayer) > critDiff) mLayerTempDiff(iLayer) = critDiff + epsT
  else                      ! (mLayerTempIter > Tcrit)
   if(mLayerTempDiff(iLayer) < critDiff) mLayerTempDiff(iLayer) = critDiff - epsT
  endif
 end do

 ! update temperatures
 if(computeVegFlux)then
  scalarCanairTempNew = scalarCanairTempIter + scalarCanairTempDiff
  scalarCanopyTempNew = scalarCanopyTempIter + scalarCanopyTempDiff
 else
  scalarCanairTempNew = missingValue_belowAbsoluteZero
  scalarCanopyTempNew = missingValue_belowAbsoluteZero
 endif
 mLayerTempNew        = mLayerTempIter + mLayerTempDiff

 ! --------------------------------------------------------------------------------------------------------------------------------
 ! * COMPUTE PHASE CHANGE...
 ! --------------------------------------------------------------------------------------------------------------------------------

 ! compute phase change of water in the vegetation canopy
 if(computeVegFlux)then
  ! compute the fraction of liquid water
  fLiq = fracliquid(scalarCanopyTempNew,snowfrz_scale)  ! fraction of liquid water (-)
  tWat = scalarCanopyLiqIter + scalarCanopyIceIter      ! total water (kg m-2)
  scalarCanopyLiqNew = fLiq*tWat                        ! mass of liquid water on the canopy (kg m-2)
  scalarCanopyIceNew = (1._dp - fLiq)*tWat              ! mass of ice on the canopy (kg m-2)
 else
  scalarCanopyLiqNew = scalarCanopyLiqIter
  scalarCanopyIceNew = scalarCanopyIceIter
 endif

 ! try and avoid oscillations associated with phase change
 if(nSnow>0)then
  if(iter>3)then
   do iSnow=1,nSnow
    ! get a trial temperature value
    tempTrial = mLayerTempNew(iSnow)
    ! check if there is liquid water
    fLiq = fracliquid(tempTrial,snowfrz_scale)  ! fraction of liquid water (-)
    if(fLiq > 0.01_dp)then
     ! save fluxes 
     nrg0 = mLayerVolHtCapBulk(iSnow)*mLayerTemp(iSnow)                                      ! energy content at the start of the time step (J m-3)
     flx0 = -(iLayerInitNrgFlux(iSnow) - iLayerInitNrgFlux(iSnow-1))/mLayerDepth(iSnow)      ! flux at the start of the time step (J m-3 s-1)
     flx1 = -(iLayerNrgFlux(iSnow) - iLayerNrgFlux(iSnow-1))/mLayerDepth(iSnow)              ! flux at the current iteration (J m-3 s-1)
     aflx = flx0*wimplicit + flx1*(1._dp - wimplicit)                                        ! average flux from start and end of time step (J m-3 s-1)
     ! get volumetric fraction of the liquid equivalent of total water (-)
     theta = mLayerVolFracIceIter(iSnow)*(iden_ice/iden_water) + mLayerVolFracLiqIter(iSnow)
     ! get initial bounds of fraction of liquid water for bi-section
     xmin = 0._dp
     xmax = 1._dp
     ! iterate
     do jiter=1,maxiter
      ! update temperature
      if(jiter > 1) tempTrial = templiquid(fLiq,snowfrz_scale)
      ! update liquid and ice content
      tempVolFracLiq = fLiq*theta
      tempVolFracIce = (theta - tempVolFracLiq)*(iden_water/iden_ice)
      ! re-compute residuals
      nrg1 = mLayerVolHtCapBulk(iSnow)*tempTrial                                              ! energy content at the current iteration (J m-3)
      phse = LH_fus*iden_ice*(tempVolFracIce - mLayerVolFracIce(iSnow))                       ! phase change term (J m-3)
      xres = nrg1 - (nrg0 + aflx*dt + phse)                                                   ! residual (J m-3)
      !write(*,'(a,1x,2(i4,1x),2(f20.10,1x),2(e20.10,1x),4(f20.10,1x))') 'iSnow, jiter, xmin, xmax, xres, phse, fLiq, tempTrial, tempVolFracIce = ', &
      !                                                                   iSnow, jiter, xmin, xmax, xres, phse, fLiq, tempTrial, tempVolFracIce 
      ! check convergence
      if(abs(xres) < nrgToler)exit
      if(jiter==maxiter)then; err=-20; message=trim(message)//'failed to converge [bi-section phase]'; return; endif
      ! update bounds
      if(xres < 0._dp)then
       xmin = fLiq
      else
       xmax = fLiq
      endif
      ! try and refine residuals
      dTheta_dT = dFracLiq_dTk(tempTrial,snowfrz_scale)*theta ! derivative in volume fraction of total water w.r.t. temperature (K-1)     
      delT      = -xres/(mLayerVolHtCapBulk(iSnow) + LH_fus*iden_water*dTheta_dT)
      fLiq      = fracliquid(tempTrial+delT,snowfrz_scale)
      ! use bi-section if outside bounds
      if(fLiq < xmin .or. fLiq > xmax) fLiq = 0.5_dp*(xmin+xmax)
     end do  ! (iterating)
     ! update temperature (NOTE: constrain temperature)
     if(tempTrial > (mLayerTempNew(iSnow)-1._dp)) mLayerTempNew(iSnow) = tempTrial
     !write(*,'(a,1x,3(f20.10,1x))') 'theta, fLiq, fracliquid(tempTrial,snowfrz_scale) = ', theta, fLiq, fracliquid(tempTrial,snowfrz_scale)
     !write(*,'(a,1x,2(f20.10,1x))') 'tempVolFracLiq, tempVolFracIce = ', tempVolFracLiq, tempVolFracIce
    endif   ! if sufficient liquid water to try and correct phase change
   end do  ! looping through snow layers
  endif   ! if sufficient iterations to try and correct phase change
 endif   ! if snow is present

 ! check ice
 !print*, 'before phase change: mLayerTempNew(1:nSnow+2) = ', mLayerTempNew(1:nSnow+2)
 !print*, 'before phase change: mLayerVolFracIceIter(1:nSnow+2) = ', mLayerVolFracIceIter(1:nSnow+2)

 ! compute phase change for the snow-soil vector
 call phsechange(mLayerTempNew,       & ! intent(in): new temperature vector (K)
                 mLayerMatricHeadIter,& ! intent(in): matric head at the current iteration (m)
                 mLayerVolFracLiqIter,& ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                 mLayerVolFracIceIter,& ! intent(in): volumetric fraction of ice at the current iteration (-)
                 mLayerMatricHeadNew, & ! intent(out): new matric head (m)
                 mLayerVolFracLiqNew, & ! intent(out): new volumetric fraction of liquid water (-)
                 mLayerVolFracIceNew, & ! intent(out): new volumetric fraction of ice (-)
                 err,cmessage)          ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 !write(*,'(a,1x,10(f20.10,1x))') 'after phase change: mLayerVolFracLiqNew(nSnow+1) = ', mLayerVolFracLiqNew(nSnow+1)
 !write(*,'(a,1x,10(f20.10,1x))') 'after phase change: mLayerVolFracIceNew(nSnow+1) = ', mLayerVolFracIceNew(nSnow+1)
 !if(iter > 3) pause

 ! ***** update the fluxes at the layer interfaces
 do iLayer=0,nLayers
  if(iLayer==0)then;           iLayerNrgFlux(iLayer) = iLayerNrgFlux(iLayer) + dFlux_dTempBelow(iLayer)*mLayerTempDiff(iLayer+1)
  elseif(iLayer==nLayers)then; iLayerNrgFlux(iLayer) = iLayerNrgFlux(iLayer) + dFlux_dTempAbove(iLayer)*mLayerTempDiff(iLayer)
  else;                        iLayerNrgFlux(iLayer) = iLayerNrgFlux(iLayer) + dFlux_dTempAbove(iLayer)*mLayerTempDiff(iLayer) &
                                                                             + dFlux_dTempBelow(iLayer)*mLayerTempDiff(iLayer+1)
  endif
 end do ! (looping through layers)

 ! deallocate space for the Jacobian matric
 if(computeVegFlux)then
  deallocate(aJac,rVec,xInc,stat=err)
  if(err/=0)then; err=20; message=trim(message)//'problem deallocating space for the Jacobian matrix'; return; endif
 endif

 ! ====================================================================================================================

 contains


  ! ************************************************************************************************
  ! internal subroutine: compute the Jacobian matrix
  ! ************************************************************************************************
  subroutine cmpJacobian(&
                         dt,                            & ! intent(in): time step (seconds)
                         scalarCanairTempInput,         & ! intent(in): trial value of canopy air space temperature (K)
                         scalarCanopyTempInput,         & ! intent(in): trial value of canopy temperature (K)
                         mLayerTempInput,               & ! intent(in): trial temperature at the current iteration (K)
                         err,message)                     ! intent(out): error control
  USE snow_utils_module,only:fracliquid    ! compute volumetric fraction of liquid water
  USE soil_utils_module,only:volFracLiq    ! compute volumetric fraction of liquid water based on matric head
  ! ----------------------------------------------------------------------------------------------------------
  implicit none
  ! dummy variables
  real(dp),intent(in)           :: dt                     ! time step (seconds)
  real(dp),intent(in)           :: scalarCanairTempInput  ! input value of canopy air space temperature (K)
  real(dp),intent(in)           :: scalarCanopyTempInput  ! input value of canopy temperature (K)
  real(dp),intent(in)           :: mLayerTempInput(:)     ! input temperature at the current iteration (K)
  integer(i4b),intent(out)      :: err                    ! error code
  character(*),intent(out)      :: message                ! error message
  ! ----------------------------------------------------------------------------------------------------------
  ! general local variables
  integer(i4b)                  :: iJac                   ! index of state variable
  integer(i4b)                  :: jStart                 ! start index for Jacobian calculations
  real(dp)                      :: scalarCanairTempTrial  ! trial value of canopy air space temperature (K)
  real(dp)                      :: scalarCanopyTempTrial  ! trial value of canopy temperature (K)
  real(dp),dimension(nLayers)   :: mLayerTempTrial        ! trial temperature at the current iteration (K)
  real(dp),allocatable          :: jMat(:,:)              ! jacobian matrix
  real(dp)                      :: fCanair                ! function test (canopy air space)
  real(dp)                      :: fCanopy                ! function test (canopy)
  real(dp),dimension(nLayers)   :: fVector                ! function test (residual snow-soil vector)
  integer(i4b),parameter        :: ixCas=1                ! index for the canopy air space
  integer(i4b),parameter        :: ixVeg=2                ! index for the vegetation canopy
  ! local variables for the canopy air space variables
  real(dp)                      :: local_Temp_CanopyAir             ! trial temperature of the canopy air space (K)
  real(dp)                      :: local_VP_CanopyAir               ! trial vapor pressure of the canopy air space (Pa)
  real(dp)                      :: local_canopyStabilityCorrection  ! stability correction for the canopy (-)
  real(dp)                      :: local_groundStabilityCorrection  ! stability correction for the ground surface (-)
  ! local variables for the mass fluxes associated with evaporation/transpiration
  real(dp)                      :: local_scalarCanopyTranspiration  ! canopy transpiration (kg m-2 s-1)
  real(dp)                      :: local_scalarCanopyEvaporation    ! canopy evaporation/condensation (kg m-2 s-1)
  real(dp)                      :: local_scalarGroundEvaporation    ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
  ! local variables for the fluxes at vegetation and ground surfaces
  real(dp)                      :: local_canairNetFlux              ! net energy flux for the canopy air space (W m-2)
  real(dp)                      :: local_canopyNetFlux              ! net energy flux for the vegetation canopy (W m-2)
  real(dp)                      :: local_groundNetFlux              ! net energy flux for the ground surface (W m-2)
  real(dp)                      :: local_dCanairNetFlux_dCanairTemp ! derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
  real(dp)                      :: local_dCanairNetFlux_dCanopyTemp ! derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
  real(dp)                      :: local_dCanairNetFlux_dGroundTemp ! derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
  real(dp)                      :: local_dCanopyNetFlux_dCanairTemp ! derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
  real(dp)                      :: local_dCanopyNetFlux_dCanopyTemp ! derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
  real(dp)                      :: local_dCanopyNetFlux_dGroundTemp ! derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
  real(dp)                      :: local_dGroundNetFlux_dCanairTemp ! derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
  real(dp)                      :: local_dGroundNetFlux_dCanopyTemp ! derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
  real(dp)                      :: local_dGroundNetFlux_dGroundTemp ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
  ! local variables for fluxes at layer interfaces
  real(dp),dimension(0:nLayers) :: local_iLayerConductiveFlux       ! conductive energy flux at layer interfaces at end of time step (W m-2)
  real(dp),dimension(0:nLayers) :: local_iLayerAdvectiveFlux        ! advective energy flux at layer interfaces at end of time step (W m-2)
  real(dp),dimension(0:nLayers) :: local_iLayerNrgFlux              ! energy flux at layer interfaces at the end of the time step (W m-2)
  real(dp),dimension(0:nLayers) :: local_dFlux_dTempAbove           ! derivative in flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
  real(dp),dimension(0:nLayers) :: local_dFlux_dTempBelow           ! derivative in flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
  ! local variables for phase change
  real(dp)                      :: vfLiq                  ! volumetric fraction of liquid water (-)
  real(dp)                      :: vfIce                  ! volumetric fraction of ice (-)
  ! ----------------------------------------------------------------------------------------------------------
  ! initialize error control
  err=0; message='cmpJacobian/'

  ! allocate space for the Jacobian matrix
  allocate(jMat(nState,nState), stat=err)
  if(err/=0)then; err=20; message=trim(message)//'problem allocating space for the Jacobian matrix'; return; endif

  ! get copies of the state vector to perturb
  if(computeVegFlux)then
   scalarCanairTempTrial = scalarCanairTempInput
   scalarCanopyTempTrial = scalarCanopyTempInput
  endif
  mLayerTempTrial = mLayerTempInput
 
  ! loop through desired layers
  do ijac=1,nState

   ! re-initialize vapor pressure of the canopy air space
   local_VP_CanopyAir = saveVP_CanopyAir
   !print*, 'local_VP_CanopyAir = ', local_VP_CanopyAir

   ! ***** perturb states
   if(computeVegFlux)then
    if(iJac==ixCas) scalarCanairTempTrial = scalarCanairTempInput + dx
    if(iJac==ixVeg) scalarCanopyTempTrial = scalarCanopyTempInput + dx
   endif
   if(iJac > ixVeg) mLayerTempTrial(iJac-ixVeg) = mLayerTempInput(iJac-ixVeg) + dx
   !write(*,'(a,i4,1x,10(f20.8,1x))') 'iJac, scalarCanairTempTrial, scalarCanopyTempTrial, mLayerTempTrial(1) = ', &
   !                                   iJac, scalarCanairTempTrial, scalarCanopyTempTrial, mLayerTempTrial(1)

   ! re-set the stability corrections
   local_CanopyStabilityCorrection = saveCanopyStabilityCorrection
   local_GroundStabilityCorrection = saveGroundStabilityCorrection
 
   ! ***** compute energy fluxes at vegetation and ground surfaces
   call vegNrgFlux(&
                   ! input
                   dt,                                  & ! intent(in): time step (seconds)
                   iter,                                & ! intent(in): iteration index
                   firstSubStep,                        & ! intent(in): flag to indicate if we are processing the first sub-step
                   computeVegFlux,                      & ! intent(in): flag to indicate if we need to compute fluxes over vegetation
                   .false.,                             & ! intent(in): flag to indicate if we need to compute shortwave radiation
                   scalarCanairTempTrial,               & ! intent(in): trial value of canopy air space temperature (K)
                   scalarCanopyTempTrial,               & ! intent(in): trial value of canopy temperature (K)
                   mLayerTempTrial(1),                  & ! intent(in): trial value of ground temperature (K)
                   scalarCanopyIceIter,                 & ! intent(in): trial mass of ice on the vegetation canopy (kg m-2)
                   scalarCanopyLiqIter,                 & ! intent(in): trial mass of liquid water on the vegetation canopy (kg m-2)
                   vegTypeIndex,                        & ! intent(in): vegetation type index
                   soilTypeIndex,                       & ! intent(in): soil type index
                   scalarLAI,                           & ! intent(in): one-sided leaf area index (m2 m-2)
                   scalarSAI,                           & ! intent(in): one-sided stem area index (m2 m-2)
                   scalarExposedLAI,                    & ! intent(in): exposed leaf area index after burial by snow (m2 m-2)
                   scalarExposedSAI,                    & ! intent(in): exposed stem area index after burial by snow (m2 m-2)
                   scalarGrowingSeasonIndex,            & ! intent(in): growing season index (0=off, 1=on)
                   scalarFoliageNitrogenFactor,         & ! intent(in): foliage nitrogen concentration (1.0 = saturated)
                   ! input/output: canopy air space variables
                   local_VP_CanopyAir,                  & ! intent(inout): trial vapor pressure of the canopy air space (Pa)
                   local_CanopyStabilityCorrection,     & ! intent(inout): stability correction for the canopy (-)
                   local_GroundStabilityCorrection,     & ! intent(inout): stability correction for the ground surface (-)
                   ! output: liquid water fluxes associated with evaporation/transpiration
                   local_scalarCanopyTranspiration,     & ! intent(out): canopy transpiration (kg m-2 s-1)
                   local_scalarCanopyEvaporation,       & ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
                   local_scalarGroundEvaporation,       & ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
                   ! output: fluxes
                   local_canairNetFlux,                 & ! intent(out): net energy flux for the canopy air space (W m-2)
                   local_canopyNetFlux,                 & ! intent(out): net energy flux for the vegetation canopy (W m-2)
                   local_groundNetFlux,                 & ! intent(out): net energy flux for the ground surface (W m-2)
                   ! output: flux derivatives
                   local_dCanairNetFlux_dCanairTemp,    & ! intent(out): derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
                   local_dCanairNetFlux_dCanopyTemp,    & ! intent(out): derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
                   local_dCanairNetFlux_dGroundTemp,    & ! intent(out): derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
                   local_dCanopyNetFlux_dCanairTemp,    & ! intent(out): derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
                   local_dCanopyNetFlux_dCanopyTemp,    & ! intent(out): derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
                   local_dCanopyNetFlux_dGroundTemp,    & ! intent(out): derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
                   local_dGroundNetFlux_dCanairTemp,    & ! intent(out): derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
                   local_dGroundNetFlux_dCanopyTemp,    & ! intent(out): derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
                   local_dGroundNetFlux_dGroundTemp,    & ! intent(out): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
                   ! output: error control
                   err,cmessage)                          ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
   !write(*,'(a,i4,10(f30.8,1x))') 'iJac, scalarCanairTempTrial, scalarCanopyTempTrial, local_canairNetFlux, local_canopyNetFlux, local_groundNetFlux = ', &
   !                                iJac, scalarCanairTempTrial, scalarCanopyTempTrial, local_canairNetFlux, local_canopyNetFlux, local_groundNetFlux
 
   ! ***** compute fluxes at layer interfaces and their derivatives (J m-2 s-1)
   call iLayer_nrg(&
                   ! (model control variables)
                   fDerivMeth,                          & ! intent(in): method used to compute derivatives (numerical or analytical)
                   ! (input)
                   layerType,                           & ! intent(in): type of each layer
                   mLayerDepth,                         & ! intent(in): depth of each layer (m)
                   mLayerHeight,                        & ! intent(in): height of layer mid-points (m)
                   iLayerThermalC,                      & ! intent(in): thermal conductivity at layer interfaces (W m-1)
                   iLayerLiqFluxSnow,                   & ! intent(in): liquid flux at the interface of each snow layer (m s-1)
                   iLayerLiqFluxSoil,                   & ! intent(in): liquid flux at the interface of each soil layer (m s-1)
                   mLayerTempTrial,                     & ! intent(in): trial temperature at the current iteration (K)
                   local_groundNetFlux,                 & ! intent(in): total flux at the ground surface (W m-2)
                   local_dGroundNetFlux_dGroundTemp,    & ! intent(in): derivative in total ground surface flux w.r.t. ground temperature (W m-2 K-1)
                   lowerBoundTemp,                      & ! intent(in): temperature of the lower boundary (K)
                   ! (output)
                   local_iLayerConductiveFlux,          & ! intent(out): conductive energy flux at layer interfaces (W m-2)
                   local_iLayerAdvectiveFlux,           & ! intent(out): advective energy flux at layer interfaces (W m-2)
                   local_iLayerNrgFlux,                 & ! intent(out): energy flux at the layer interfaces (W m-2)
                   local_dFlux_dTempAbove,              & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
                   local_dFlux_dTempBelow,              & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
                   err,cmessage)                          ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 
   ! ***** compute the residual for the vegetation canopy and the canopy air space (J m-3)
   if(computeVegFlux)then
    ! compute the residual for the vegetation canopy 
    nrg0 = scalarBulkVolHeatCapVeg*scalarCanopyTemp                      ! energy content at the start of the time step (J m-3)
    nrg1 = scalarBulkVolHeatCapVeg*scalarCanopyTempTrial                 ! energy content at the current iteration (J m-3)
    flx0 = canopyNetFluxInit/canopyDepth                                 ! flux at the start of the time step (J m-3 s-1)
    flx1 = local_canopyNetFlux/canopyDepth                               ! flux at the current iteration (J m-3 s-1)
    phse = LH_fus*(scalarCanopyIceIter - scalarCanopyIce)/canopyDepth    ! phase change term (J m-3)
    fCanopy = nrg1 - (nrg0 + (flx0*wimplicit + flx1*(1._dp - wimplicit))*dt + phse)
    ! compute the residual for the canopy air space
    nrg0 = Cp_air*iden_air*scalarCanairTemp                              ! energy content at the start of the time step (J m-3)
    nrg1 = Cp_air*iden_air*scalarCanairTempTrial                         ! energy content at the current iteration (J m-3)
    flx0 = canairNetFluxInit/canopyDepth                                 ! flux at the start of the time step (J m-3 s-1)
    flx1 = local_canairNetFlux/canopyDepth                               ! flux at the current iteration (J m-3 s-1)
    ! (compute residuals)
    fCanair = nrg1 - (nrg0 + (flx0*wimplicit + flx1*(1._dp - wimplicit))*dt)
   endif 

   ! ***** compute the residual vector for all snow/soil layers (J m-3)
   do iLayer=1,nLayers
    ! (compute energy storage)
    nrg0 = mLayerVolHtCapBulk(iLayer)*mLayerTemp(iLayer)                                        ! energy content at the start of the time step (J m-3)
    nrg1 = mLayerVolHtCapBulk(iLayer)*mLayerTempTrial(iLayer)                                   ! energy content at the current iteration (J m-3)
    ! (compute energy fluxes)
    flx0 = -(iLayerInitNrgFlux(iLayer) - iLayerInitNrgFlux(iLayer-1))/mLayerDepth(iLayer)       ! flux at the start of the time step (J m-3 s-1)
    flx1 = -(local_iLayerNrgFlux(iLayer) - local_iLayerNrgFlux(iLayer-1))/mLayerDepth(iLayer)   ! flux at the current iteration (J m-3 s-1)
    ! (compute phase change)
    theta = mLayerVolFracIceIter(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiqIter(iLayer)
    select case(layerType(iLayer))
     case(ix_snow); vfLiq = fracliquid(mLayerTempTrial(iLayer),snowfrz_scale)*theta
     case(ix_soil)
      if(mLayerTempTrial(iLayer) < mLayerTcrit(iLayer-nSnow))then
       vfLiq = volFracLiq(kappa*(mLayerTempTrial(iLayer) - Tfreeze),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
      else
       vfLiq = theta
      endif
     case default; err=20; message=trim(message)//'cannot identify layer type'; return
    endselect
    vfIce = (theta - vfLiq)*(iden_water/iden_ice)
    phse  = LH_fus*iden_ice*(vfIce - mLayerVolFracIce(iLayer))            ! phase change term (J m-3)
    !write(*,'(a,10(f20.10,1x))') 'in Jacobian, computing residuals: nrg0, nrg1, flx0*dt, flx1*dt = ', nrg0, nrg1, flx0*dt, flx1*dt, mLayerTempTrial(iLayer)
    ! (compute residuals)
    fVector(iLayer) = nrg1 - (nrg0 + (flx0*wimplicit + flx1*(1._dp - wimplicit))*dt + phse)
   end do
   !write(*,'(a,10(e20.10,1x))') 'fCanair,        fCanopy        = ', fCanair,        fCanopy
   !write(*,'(a,10(e20.10,1x))') 'canairResidual, canopyResidual = ', canairResidual, canopyResidual
   !write(*,'(a,10(e20.10,1x))') 'fVector        = ', fVector
   !write(*,'(a,10(e20.10,1x))') 'mLayerResidual = ', mLayerResidual
 
   ! compute Jacobian
   if(computeVegFlux)then
    jmat(:,ijac) = ( (/fCanair,fCanopy,fVector(:)/) - (/canairResidual,canopyResidual,mLayerResidual(:)/) ) / dx
   else
    jmat(:,ijac) = (fVector(:) - mLayerResidual(:) ) / dx
   endif
   !print*, '=========='
 
   ! set the state back to the input value
   if(computeVegFlux)then
    if(iJac==ixCas) scalarCanairTempTrial = scalarCanairTempInput
    if(iJac==ixVeg) scalarCanopyTempTrial = scalarCanopyTempInput
   endif
   if(iJac > ixVeg) mLayerTempTrial(iJac-ixVeg) = mLayerTempInput(iJac-ixVeg)
 
  end do  ! looping through snow-soil layers

  ! print the Jacobian
  ! print Jacobian
  print*, 'analytical jacobian'
  do iLayer=1,nState
   write(*,'(15(f18.5,1x))') aJac(:,iLayer)
  end do
  print*, 'numerical Jacobian'
  do iLayer=1,nState
   write(*,'(15(f18.5,1x))') jMat(:,iLayer)
  end do
  !pause 'testing Jacobian'

  ! deallocate space for the Jacobian matrix
  deallocate(jMat, stat=err) 
  if(err/=0)then; err=20; message=trim(message)//'problem de-allocating space for the Jacobian matrix'; return; endif

  end subroutine cmpJacobian

  ! ===== end internal subroutines

 end subroutine heatTransf_muster


 ! ************************************************************************************************
 ! private subroutine: compute energy fluxes at layer interfaces, and their derivatives
 ! ************************************************************************************************
 subroutine iLayer_nrg(&
                       ! (model control variables)
                       fDerivMeth,                         & ! intent(in): method used to compute derivatives (numerical or analytical)
                       ! (input)
                       layerType,                          & ! intent(in): type of each layer
                       mLayerDepth,                        & ! intent(in): depth of each layer (m)
                       mLayerHeight,                       & ! intent(in): height of layer mid-points (m)
                       iLayerThermalC,                     & ! intent(in): thermal conductivity at layer interfaces (W m-1)
                       iLayerLiqFluxSnow,                  & ! intent(in): liquid flux at the interface of each snow layer (m s-1)
                       iLayerLiqFluxSoil,                  & ! intent(in): liquid flux at the interface of each soil layer (m s-1)
                       mLayerTempTrial,                    & ! intent(in): trial temperature at the current iteration (K)
                       groundNetFlux,                      & ! intent(in): total flux at the ground surface (W m-2)
                       dGroundNetFlux_dGroundTemp,         & ! intent(in): derivative in total ground surface flux w.r.t. ground temperature (W m-2 K-1)
                       lowerBoundTemp,                     & ! intent(in): temperature of the lower boundary (K)
                       ! (output)
                       iLayerConductiveFlux,               & ! intent(out): conductive energy flux at layer interfaces (W m-2)
                       iLayerAdvectiveFlux,                & ! intent(out): advective energy flux at layer interfaces (W m-2)
                       iLayerNrgFlux,                      & ! intent(out): energy flux at the layer interfaces (W m-2)
                       dFlux_dTempAbove,                   & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (W m-2 K-1)
                       dFlux_dTempBelow,                   & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (W m-2 K-1)
                       err,message)                          ! intent(out): error control
 ! compute derivative in fluxes at layer interfaces w.r.t. temperature in the layer above and the layer below
 implicit none
 ! input
 integer(i4b),intent(in)       :: fDerivMeth                 ! intent(in): method used to calculate derivatives
 integer(i4b),intent(in)       :: layerType(:)               ! intent(in): type of the layer (ix_soil or ix_snow)
 real(dp),intent(in)           :: mLayerDepth(:)             ! intent(in): depth of each layer (m)
 real(dp),intent(in)           :: mLayerHeight(:)            ! intent(in): height of layer mid-points (m)
 real(dp),intent(in)           :: iLayerThermalC(0:)         ! intent(in): thermal conductivity at layer interfaces (W m-1)
 real(dp),intent(in)           :: iLayerLiqFluxSnow(0:)      ! intent(in): liquid flux at the interface of each snow layer (m s-1)
 real(dp),intent(in)           :: iLayerLiqFluxSoil(0:)      ! intent(in): liquid flux at the interface of each soil layer (m s-1)
 real(dp),intent(in)           :: mLayerTempTrial(:)         ! intent(in): trial temperature at the current iteration (K)
 real(dp),intent(in)           :: groundNetFlux              ! intent(in): total flux at the ground surface (W m-2)
 real(dp),intent(in)           :: dGroundNetFlux_dGroundTemp ! intent(in): derivative in total ground surface flux w.r.t. ground temperature (W m-2 K-1)
 real(dp),intent(in)           :: lowerBoundTemp             ! intent(in): temperature of the lower boundary (K)
 ! output
 real(dp),intent(out)          :: iLayerConductiveFlux(0:)   ! intent(out): conductive energy flux at layer interfaces at end of time step (W m-2)
 real(dp),intent(out)          :: iLayerAdvectiveFlux(0:)    ! intent(out): advective energy flux at layer interfaces at end of time step (W m-2)
 real(dp),intent(out)          :: iLayerNrgFlux(0:)          ! intent(out): energy flux at the layer interfaces (W m-2)
 real(dp),intent(out)          :: dFlux_dTempAbove(0:)       ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (J m-2 s-1 K-1)
 real(dp),intent(out)          :: dFlux_dTempBelow(0:)       ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (J m-2 s-1 K-1)
 integer(i4b),intent(out)      :: err                        ! intent(out): error code
 character(*),intent(out)      :: message                    ! intent(out): error message
 ! local variables
 integer(i4b)                  :: iLayer                     ! index of model layers
 real(dp)                      :: qFlux                      ! liquid flux at layer interfaces (m s-1)
 real(dp)                      :: dz                         ! height difference (m)
 real(dp)                      :: flux0,flux1,flux2          ! fluxes used to calculate derivatives (W m-2)
 ! initialize error control
 err=0; message='iLayer_nrg/'

 ! set conductive and advective fluxes to missing in the upper boundary
 ! NOTE: advective flux at the upper boundary is included in the ground heat flux
 iLayerConductiveFlux(0) = valueMissing
 iLayerAdvectiveFlux(0)  = valueMissing

 ! -------------------------------------------------------------------------------------------------------------------------
 ! ***** compute the conductive fluxes at layer interfaces *****
 ! -------------------------------------------------------------------------------------------------------------------------
 do iLayer=1,nLayers
  ! compute fluxes at the lower boundary -- positive downwards
  if(iLayer==nLayers)then; iLayerConductiveFlux(nLayers) = -iLayerThermalC(iLayer)*(lowerBoundTemp - mLayerTempTrial(iLayer))/(mLayerDepth(iLayer)*0.5_dp)
  ! compute fluxes within the domain -- positive downwards
  else;                    iLayerConductiveFlux(iLayer)  = -iLayerThermalC(iLayer)*(mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer)) / &
                                                                                   (mLayerHeight(iLayer+1) - mLayerHeight(iLayer))
  endif ! (the type of layer)
 end do

 ! -------------------------------------------------------------------------------------------------------------------------
 ! ***** compute the advective fluxes at layer interfaces *****
 ! -------------------------------------------------------------------------------------------------------------------------
 do iLayer=1,nLayers
  ! get the liquid flux at layer interfaces
  select case(layerType(iLayer))
   case(ix_snow); qFlux = iLayerLiqFluxSnow(iLayer)
   case(ix_soil); qFlux = iLayerLiqFluxSoil(iLayer-nSnow)
   case default; err=20; message=trim(message)//'unable to identify layer type'; return
  end select
  ! compute fluxes at the lower boundary -- positive downwards
  if(iLayer==nLayers)then
   iLayerAdvectiveFlux(iLayer) = -Cp_water*iden_water*qFlux*(lowerBoundTemp - mLayerTempTrial(iLayer))
  ! compute fluxes within the domain -- positive downwards
  else
   iLayerAdvectiveFlux(iLayer) = -Cp_water*iden_water*qFlux*(mLayerTempTrial(iLayer+1) - mLayerTempTrial(iLayer))
  endif
 end do  ! looping through layers

 ! -------------------------------------------------------------------------------------------------------------------------
 ! ***** compute the total fluxes at layer interfaces *****
 ! -------------------------------------------------------------------------------------------------------------------------
 ! NOTE: ignore advective fluxes for now
 iLayerNrgFlux(0)         = groundNetFlux
 iLayerNrgFlux(1:nLayers) = iLayerConductiveFlux(1:nLayers)

 ! -------------------------------------------------------------------------------------------------------------------------
 ! ***** compute the derivative in fluxes at layer interfaces w.r.t temperature in the layer above and the layer below *****
 ! -------------------------------------------------------------------------------------------------------------------------

 ! initialize un-used elements
 dFlux_dTempBelow(nLayers) = -huge(lowerBoundTemp)  ! don't expect this to be used, so deliberately set to a ridiculous value to cause problems

 ! loop through INTERFACES...
 do iLayer=0,nLayers

  ! ***** the upper boundary -- ** NOTE: dTotalSurfaceFlux_dTemp was computed previously using fDerivMeth
  if(iLayer==0)then  ! (upper boundary)
   dFlux_dTempBelow(iLayer) = dGroundNetFlux_dGroundTemp

  ! ***** the lower boundary
  elseif(iLayer==nLayers)then  ! (lower boundary)
   dz = mLayerDepth(iLayer)*0.5_dp
   if(fDerivMeth==analytical)then    ! ** analytical derivatives
    dFlux_dTempAbove(iLayer) = iLayerThermalC(iLayer)/dz
   else                              ! ** numerical derivatives
    flux0 = -iLayerThermalC(iLayer)*(lowerBoundTemp - (mLayerTempTrial(iLayer)   ))/dz
    flux1 = -iLayerThermalC(iLayer)*(lowerBoundTemp - (mLayerTempTrial(iLayer)+dx))/dz
    dFlux_dTempAbove(iLayer) = (flux1 - flux0)/dx
   endif

  ! ***** internal layers
  else
   dz = (mLayerHeight(iLayer+1) - mLayerHeight(iLayer))
   if(fDerivMeth==analytical)then    ! ** analytical derivatives
    dFlux_dTempAbove(iLayer) =  iLayerThermalC(iLayer)/dz
    dFlux_dTempBelow(iLayer) = -iLayerThermalC(iLayer)/dz
   else                              ! ** numerical derivatives
    flux0 = -iLayerThermalC(iLayer)*( mLayerTempTrial(iLayer+1)     -  mLayerTempTrial(iLayer)    ) / dz
    flux1 = -iLayerThermalC(iLayer)*( mLayerTempTrial(iLayer+1)     - (mLayerTempTrial(iLayer)+dx)) / dz
    flux2 = -iLayerThermalC(iLayer)*((mLayerTempTrial(iLayer+1)+dx) -  mLayerTempTrial(iLayer)    ) / dz
    dFlux_dTempAbove(iLayer) = (flux1 - flux0)/dx
    dFlux_dTempBelow(iLayer) = (flux2 - flux0)/dx
   endif

  endif  ! type of layer (upper, internal, or lower)

 end do  ! (looping through layers)
 
 end subroutine iLayer_nrg



end module heatTransf_module


module soilLiqFlx_module
! -----------------------------------------------------------------------------------------------------------
! numerical recipes data types
USE nrtype
! physical constants
USE multiconst,only:&
                    LH_fus,  & ! latent heat of fusion         (J kg-1)
                    LH_vap,  & ! latent heat of vaporization   (J kg-1)
                    LH_sub,  & ! latent heat of sublimation    (J kg-1)
                    gravity, & ! gravitational acceleteration  (m s-2)
                    Tfreeze, & ! freezing point of pure water  (K)
                    iden_air,& ! intrinsic density of air      (kg m-3)
                    iden_ice,& ! intrinsic density of ice      (kg m-3)
                    iden_water ! intrinsic density of water    (kg m-3)
! provide access to the number of snow and soil layers
USE data_struc,only:&
                    nSnow,   & ! number of snow layers  
                    nSoil,   & ! number of soil layers  
                    nLayers    ! total number of layers
! provide access to layer types
USE data_struc,only:ix_soil,ix_snow  ! named variables for snow and soil
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
 zeroFlux                      ! zero flux
! -----------------------------------------------------------------------------------------------------------
implicit none
private
public::soilLiqFlx
! constant parameters
real(dp),parameter     :: valueMissing=-9999._dp    ! missing value parameter
real(dp),parameter     :: verySmall=1.e-12_dp       ! a very small number (used to avoid divide by zero)
real(dp),parameter     :: dx=1.e-8_dp               ! finite difference increment
contains


 ! ************************************************************************************************
 ! new subroutine: compute liquid water fluxes and their derivatives
 ! ************************************************************************************************
 subroutine soilLiqFlx(&
                       ! input: model control
                       doInfiltrate,                 & ! intent(in): flag to compute infiltration
                       deriv_desired,                & ! intent(in): flag indicating if derivatives are desired
                       ! input: trial state variables
                       mLayerTempTrial,              & ! intent(in): temperature (K)
                       mLayerMatricHeadTrial,        & ! intent(in): matric head (m)
                       mLayerVolFracLiqTrial,        & ! intent(in): volumetric fraction of liquid water (-)
                       mLayerVolFracIceTrial,        & ! intent(in): volumetric fraction of ice (-)
                       ! input: pre-computed derivatives
                       mLayerdTheta_dTk,             & ! intent(in): derivative in volumetric liquid water content w.r.t. temperature (K-1)
                       dPsiLiq_dTemp,                & ! intent(in): derivative in liquid water matric potential w.r.t. temperature (m K-1)
                       ! input: fluxes
                       scalarCanopyTranspiration,    & ! intent(in): canopy transpiration (kg m-2 s-1)
                       scalarGroundEvaporation,      & ! intent(in): ground evaporation (kg m-2 s-1)
                       scalarRainPlusMelt,           & ! intent(in): rain plus melt (m s-1)
                       ! output: diagnostic variables for surface runoff
                       xMaxInfilRate,                & ! intent(inout): maximum infiltration rate (m s-1)
                       scalarInfilArea,              & ! intent(inout): fraction of unfrozen area where water can infiltrate (-)
                       scalarFrozenArea,             & ! intent(inout): fraction of area that is considered impermeable due to soil ice (-)
                       scalarSurfaceRunoff,          & ! intent(out): surface runoff (m s-1)
                       ! output: diagnostic variables for model layers
                       mLayerdTheta_dPsi,            & ! intent(out): derivative in the soil water characteristic w.r.t. psi (m-1)
                       mLayerdPsi_dTheta,            & ! intent(out): derivative in the soil water characteristic w.r.t. theta (m)
                       dHydCond_dMatric,             & ! intent(out): derivative in hydraulic conductivity w.r.t matric head (s-1)
                       ! output: fluxes
                       scalarSurfaceInfiltration,    & ! intent(out): surface infiltration rate (m s-1)
                       iLayerLiqFluxSoil,            & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                       mLayerTranspire,              & ! intent(out): transpiration loss from each soil layer (m s-1)
                       mLayerHydCond,                & ! intent(out): hydraulic conductivity in each soil layer (m s-1)
                       ! output: derivatives in fluxes w.r.t. hydrology state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                       dq_dHydStateAbove,            & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer above (m s-1)
                       dq_dHydStateBelow,            & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer below (m s-1)
                       ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
                       dq_dNrgStateAbove,            & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
                       dq_dNrgStateBelow,            & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
                       ! output: error control
                       err,message)                    ! intent(out): error control
 ! model decisions
 USE data_struc,only:model_decisions                         ! model decision structure
 USE var_lookup,only:iLookDECISIONS                          ! named variables for elements of the decision structure
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:attr_data,type_data,mpar_data,forc_data,mvar_data,indx_data    ! data structures
 USE var_lookup,only:iLookATTR,iLookTYPE,iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX ! named variables for structure elements
 implicit none
 ! input: model control
 logical(lgt),intent(in)          :: doInfiltrate                  ! flag to compute infiltration
 logical(lgt),intent(in)          :: deriv_desired                 ! flag indicating if derivatives are desired
 ! input: trial model state variables
 real(dp),intent(in)              :: mLayerTempTrial(:)            ! temperature in each layer at the current iteration (m)
 real(dp),intent(in)              :: mLayerMatricHeadTrial(:)      ! matric head in each layer at the current iteration (m)
 real(dp),intent(in)              :: mLayerVolFracLiqTrial(:)      ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),intent(in)              :: mLayerVolFracIceTrial(:)      ! volumetric fraction of ice at the current iteration (-)
 ! input: pre-computed derivatves
 real(dp),intent(in)              :: mLayerdTheta_dTk(:)           ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
 real(dp),intent(in)              :: dPsiLiq_dTemp(:)              ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
 ! input: model fluxes
 real(dp),intent(in)              :: scalarCanopyTranspiration     ! canopy transpiration (kg m-2 s-1)
 real(dp),intent(in)              :: scalarGroundEvaporation       ! ground evaporation (kg m-2 s-1)
 real(dp),intent(in)              :: scalarRainPlusMelt            ! rain plus melt (m s-1) 
 ! output: diagnostic variables for surface runoff
 real(dp),intent(inout)           :: xMaxInfilRate                 ! maximum infiltration rate (m s-1)
 real(dp),intent(inout)           :: scalarInfilArea               ! fraction of unfrozen area where water can infiltrate (-)
 real(dp),intent(inout)           :: scalarFrozenArea              ! fraction of area that is considered impermeable due to soil ice (-)
 real(dp),intent(out)             :: scalarSurfaceRunoff           ! surface runoff (m s-1)
 ! output: diagnostic variables for each layer
 real(dp),intent(out)             :: mLayerdTheta_dPsi(:)          ! derivative in the soil water characteristic w.r.t. psi (m-1)
 real(dp),intent(out)             :: mLayerdPsi_dTheta(:)          ! derivative in the soil water characteristic w.r.t. theta (m)
 real(dp),intent(out)             :: dHydCond_dMatric(:)           ! derivative in hydraulic conductivity w.r.t matric head (s-1)
 ! output: liquid fluxes
 real(dp),intent(out)             :: scalarSurfaceInfiltration     ! surface infiltration rate (m s-1)
 real(dp),intent(out)             :: iLayerLiqFluxSoil(0:)         ! liquid flux at soil layer interfaces (m s-1)
 real(dp),intent(out)             :: mLayerTranspire(:)            ! transpiration loss from each soil layer (m s-1)
 real(dp),intent(out)             :: mLayerHydCond(:)              ! hydraulic conductivity in each soil layer (m s-1)
 ! output: derivatives in fluxes w.r.t. state variables in the layer above and layer below (m s-1)
 real(dp),intent(out)             :: dq_dHydStateAbove(0:)         ! derivative in the flux in layer interfaces w.r.t. state variables in the layer above
 real(dp),intent(out)             :: dq_dHydStateBelow(0:)         ! derivative in the flux in layer interfaces w.r.t. state variables in the layer below
 ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
 real(dp),intent(out)             :: dq_dNrgStateAbove(0:)         ! derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
 real(dp),intent(out)             :: dq_dNrgStateBelow(0:)         ! derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
 ! output: error control
 integer(i4b),intent(out)         :: err                           ! error code
 character(*),intent(out)         :: message                       ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b)                     :: ibeg,iend                     ! start and end indices of the soil layers in concatanated snow-soil vector
 character(LEN=256)               :: cmessage                      ! error message of downwind routine
 ! initialize error control
 err=0; message='soilLiqFlx/'

 ! get indices for the data structures
 ibeg = nSnow + 1
 iend = nSnow + nSoil

 ! wrapper routine for liquid fluxes
 call soilLiqFlx_muster(&

                        ! input: model control
                        doInfiltrate,                                                  & ! intent(in): flag to compute infiltration
                        deriv_desired,                                                 & ! intent(in): flag indicating if derivatives are desired

                        ! input: model decisions
                        model_decisions(iLookDECISIONS%fDerivMeth)%iDecision,          & ! intent(in): index of the method used to calculate flux derivatives 
                        model_decisions(iLookDECISIONS%f_Richards)%iDecision,          & ! intent(in): index of the form of Richards' equation
                        model_decisions(iLookDECISIONS%hc_Profile)%iDecision,          & ! intent(in): index of the option for the hydraulic conductivity profile
                        model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision,          & ! intent(in): index of the upper boundary conditions for soil hydrology
                        model_decisions(iLookDECISIONS%bcLowrSoiH)%iDecision,          & ! intent(in): index of the lower boundary conditions for soil hydrology

                        ! input: trial state variables
                        mLayerTempTrial,                                               & ! intent(in): temperature (K)
                        mLayerMatricHeadTrial,                                         & ! intent(in): matric head (m)
                        mLayerVolFracLiqTrial,                                         & ! intent(in): volumetric fraction of liquid water (-)
                        mLayerVolFracIceTrial,                                         & ! intent(in): volumetric fraction of ice (-)

                        ! input: pre-computed derivatives
                        mLayerdTheta_dTk,                                              & ! intent(in): derivative in volumetric liquid water content w.r.t. temperature (K-1)
                        dPsiLiq_dTemp,                                                 & ! intent(in): derivative in liquid water matric potential w.r.t. temperature (m K-1)

                        ! input: fluxes
                        scalarCanopyTranspiration,                                     & ! intent(in): canopy transpiration (kg m-2 s-1)
                        scalarGroundEvaporation,                                       & ! intent(in): ground evaporation (kg m-2 s-1)
                        scalarRainPlusMelt,                                            & ! intent(in): rain plus melt (m s-1)

                        ! input: model coordinate variables -- NOTE: use of ibeg and iend 
                        mvar_data%var(iLookMVAR%mLayerDepth)%dat(ibeg:iend),           & ! intent(in): depth of the layer (m)
                        mvar_data%var(iLookMVAR%mLayerHeight)%dat(ibeg:iend),          & ! intent(in): height of the layer mid-point (m)
                        mvar_data%var(iLookMVAR%iLayerHeight)%dat(ibeg-1:iend),        & ! intent(in): height of the layer interfaces (m)

                        ! input: upper boundary conditions
                        mpar_data%var(iLookPARAM%upperBoundHead),                      & ! intent(in): upper boundary condition for matric head (m)
                        mpar_data%var(iLookPARAM%upperBoundTheta),                     & ! intent(in): upper boundary condition for volumetric liquid water content (-)

                        ! input: lower boundary conditions
                        mpar_data%var(iLookPARAM%lowerBoundHead),                      & ! intent(in): lower boundary condition for matric head (m)
                        mpar_data%var(iLookPARAM%lowerBoundTheta),                     & ! intent(in): lower boundary condition for volumetric liquid water content (-)

                        ! input: soil parameters
                        mpar_data%var(iLookPARAM%vGn_alpha),                           & ! intent(in): van Genutchen "alpha" parameter (m-1)
                        mpar_data%var(iLookPARAM%vGn_n),                               & ! intent(in): van Genutchen "n" parameter (-)
                        mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1),                   & ! intent(in): van Genutchen "m" parameter (-)
                        mpar_data%var(iLookPARAM%mpExp),                               & ! intent(in): empirical exponent in macropore flow equation (-)
                        mpar_data%var(iLookPARAM%theta_mp),                            & ! intent(in): volumetric liquid water content when macropore flow begins (-)
                        mpar_data%var(iLookPARAM%theta_sat),                           & ! intent(in): soil porosity (-)
                        mpar_data%var(iLookPARAM%theta_res),                           & ! intent(in): soil residual volumetric water content (-)
                        mpar_data%var(iLookPARAM%wettingFrontSuction),                 & ! intent(in): Green-Ampt wetting front suction (m)
                        mpar_data%var(iLookPARAM%fieldCapacity),                       & ! intent(in): field capacity (-)
                        mpar_data%var(iLookPARAM%rootingDepth),                        & ! intent(in): rooting depth (m)
                        mpar_data%var(iLookPARAM%kAnisotropic),                        & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                        mpar_data%var(iLookPARAM%zScale_TOPMODEL),                     & ! intent(in): TOPMODEL scaling factor (m)
                        mpar_data%var(iLookPARAM%qSurfScale),                          & ! intent(in): scaling factor in the surface runoff parameterization (-)
                        mpar_data%var(iLookPARAM%specificYield),                       & ! intent(in): specific yield (-)
                        mpar_data%var(iLookPARAM%specificStorage),                     & ! intent(in): specific storage coefficient (m-1)
                        mpar_data%var(iLookPARAM%f_impede),                            & ! intent(in): ice impedence factor (-)
                        mpar_data%var(iLookPARAM%soilIceScale),                        & ! intent(in): scaling factor for depth of soil ice, used to get frozen fraction (m)
                        mpar_data%var(iLookPARAM%soilIceCV),                           & ! intent(in): CV of depth of soil ice, used to get frozen fraction (-)

                        ! input: saturated hydraulic conductivity
                        mvar_data%var(iLookMVAR%mLayerSatHydCondMP)%dat,               & ! intent(in): saturated hydraulic conductivity of macropores at the mid-point of each layer (m s-1)
                        mvar_data%var(iLookMVAR%mLayerSatHydCond)%dat,                 & ! intent(in): saturated hydraulic conductivity at the mid-point of each layer (m s-1)
                        mvar_data%var(iLookMVAR%iLayerSatHydCond)%dat,                 & ! intent(in): saturated hydraulic conductivity at the interface of each layer (m s-1)

                        ! input: factors limiting transpiration (from vegFlux routine)
                        mvar_data%var(iLookMVAR%mLayerRootDensity)%dat,                & ! intent(in): root density in each layer (-)
                        mvar_data%var(iLookMVAR%scalarTranspireLim)%dat(1),            & ! intent(in): weighted average of the transpiration limiting factor (-)
                        mvar_data%var(iLookMVAR%mLayerTranspireLim)%dat,               & ! intent(in): transpiration limiting factor in each layer (-)

                        ! output: diagnostic scalar variables
                        xMaxInfilRate,                                                 & ! intent(inout): maximum infiltration rate (m s-1)
                        scalarInfilArea,                                               & ! intent(inout): fraction of unfrozen area where water can infiltrate (-)
                        scalarFrozenArea,                                              & ! intent(inout): fraction of area that is considered impermeable due to soil ice (-)
                        scalarSurfaceRunoff,                                           & ! intent(out):   surface runoff (m s-1)

                        ! output: diagnostic variables for model layers
                        mLayerdTheta_dPsi,                                             & ! intent(out): derivative in the soil water characteristic w.r.t. psi (m-1)
                        mLayerdPsi_dTheta,                                             & ! intent(out): derivative in the soil water characteristic w.r.t. theta (m)
                        dHydCond_dMatric,                                              & ! intent(out): derivative in hydraulic conductivity w.r.t matric head (s-1)

                        ! output: fluxes
                        scalarSurfaceInfiltration,                                     & ! intent(out): surface infiltration rate (m s-1)
                        iLayerLiqFluxSoil,                                             & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                        mLayerTranspire,                                               & ! intent(out): transpiration loss from each soil layer (m s-1)
                        mLayerHydCond,                                                 & ! intent(out): hydraulic conductivity in each soil layer (m s-1)

                        ! output: derivatives in fluxes w.r.t. hydrology state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                        dq_dHydStateAbove,                                             & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer above (m s-1)
                        dq_dHydStateBelow,                                             & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer below (m s-1)

                        ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
                        dq_dNrgStateAbove,                                             & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
                        dq_dNrgStateBelow,                                             & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)

                        ! output: error control
                        err,cmessage)                                                    ! intent(out): error control

 ! check for errors
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! save information in the data structures
 mvar_data%var(iLookMVAR%mLayerdTheta_dPsi)%dat(1:nSoil) = mLayerdTheta_dPsi        ! derivative in the soil water characteristic w.r.t. psi (m-1)
 mvar_data%var(iLookMVAR%mLayerdPsi_dTheta)%dat(1:nSoil) = mLayerdPsi_dTheta        ! derivative in the soil water characteristic w.r.t. theta (m)

 end subroutine soilLiqFlx

 ! ****************************************************************************************************************************************************************************
 ! ****************************************************************************************************************************************************************************
 ! *** PRIVATE SUBROUTINES ****************************************************************************************************************************************************
 ! ****************************************************************************************************************************************************************************
 ! ****************************************************************************************************************************************************************************

 ! ************************************************************************************************
 ! new subroutine: wrapper routine to compute liquid water fluxes and their derivatives
 ! ************************************************************************************************
 subroutine soilLiqFlx_muster(&
 
                              ! input: model control
                              doInfiltrate,                & ! intent(in): flag to compute infiltration
                              deriv_desired,               & ! intent(in): flag indicating if derivatives are desired

                              ! model decisions
                              ixDerivMethod,               & ! intent(in): choice of method used to compute derivative
                              ixRichards,                  & ! intent(in): choice of the form of Richards' equation
                              hc_Profile,                  & ! intent(in): index defining the option for the hydraulic conductivity profile
                              ixBcUpperSoilHydrology,      & ! intent(in): choice of upper boundary condition for soil hydrology
                              ixBcLowerSoilHydrology,      & ! intent(in): choice of lower boundary condition for soil hydrology

                              ! input: trial state variables
                              mLayerTempTrial,             & ! intent(in): temperature (m)
                              mLayerMatricHeadTrial,       & ! intent(in): matric head (m)
                              mLayerVolFracLiqTrial,       & ! intent(in): volumetric fraction of liquid water (-)
                              mLayerVolFracIceTrial,       & ! intent(in): volumetric fraction of ice (-)

                              ! input: pre-computed derivatives
                              mLayerdTheta_dTk,            & ! intent(in): derivative in volumetric liquid water content w.r.t. temperature (K-1)
                              dPsiLiq_dTemp,               & ! intent(in): derivative in liquid water matric potential w.r.t. temperature (m K-1)

                              ! input: fluxes
                              scalarCanopyTranspiration,   & ! intent(in): canopy transpiration (kg m-2 s-1)
                              scalarGroundEvaporation,     & ! intent(in): ground evaporation (kg m-2 s-1)
                              scalarRainPlusMelt,          & ! intent(in): rain plus melt (m s-1)

                              ! input: model coordinate variables -- NOTE: use of ibeg and iend to restrict attention to soil
                              mLayerDepth,                 & ! intent(in): depth of the layer (m)
                              mLayerHeight,                & ! intent(in): height of the layer mid-point (m)
                              iLayerHeight,                & ! intent(in): height of the layer interfaces (m)

                              ! input: upper boundary conditions
                              upperBoundHead,              & ! intent(in): upper boundary condition for matric head (m)
                              upperBoundTheta,             & ! intent(in): upper boundary condition for volumetric liquid water content (-)

                              ! input: lower boundary conditions
                              lowerBoundHead,              & ! intent(in): lower boundary condition for matric head (m)
                              lowerBoundTheta,             & ! intent(in): lower boundary condition for volumetric liquid water content (-)

                              ! input: soil parameters
                              vGn_alpha,                   & ! intent(in): van Genutchen "alpha" parameter (m-1)
                              vGn_n,                       & ! intent(in): van Genutchen "n" parameter (-)
                              VGn_m,                       & ! intent(in): van Genutchen "m" parameter (-)
                              mpExp,                       & ! intent(in): empirical exponent in macropore flow equation (-)
                              theta_mp,                    & ! intent(in): volumetric liquid water content when macropore flow begins (-)
                              theta_sat,                   & ! intent(in): soil porosity (-)
                              theta_res,                   & ! intent(in): soil residual volumetric water content (-)
                              wettingFrontSuction,         & ! intent(in): Green-Ampt wetting front suction (m)
                              fieldCapacity,               & ! intent(in): field capacity (-)
                              rootingDepth,                & ! intent(in): rooting depth (m)
                              kAnisotropic,                & ! intent(in): anisotropy factor for lateral hydraulic conductivity (-)
                              zScale_TOPMODEL,             & ! intent(in): TOPMODEL scaling factor (m)
                              qSurfScale,                  & ! intent(in): scaling factor in the surface runoff parameterization (-)
                              specificYield,               & ! intent(in): specific yield (-)
                              specificStorage,             & ! intent(in): specific storage coefficient (m-1)
                              f_impede,                    & ! intent(in): ice impedence factor (-)
                              soilIceScale,                & ! intent(in): scaling factor for depth of soil ice, used to get frozen fraction (m)
                              soilIceCV,                   & ! intent(in): CV of depth of soil ice, used to get frozen fraction (-)

                              ! input: saturated hydraulic conductivity in each layer
                              mLayerSatHydCondMP,          & ! intent(in): saturated hydraulic conductivity of macropores at the mid-point of each layer (m s-1)
                              mLayerSatHydCond,            & ! intent(in): saturated hydraulic conductivity at the mid-point of each layer (m s-1)
                              iLayerSatHydCond,            & ! intent(in): saturated hydraulic conductivity at the interface of each layer (m s-1)

                              ! input: factors limiting transpiration (from vegFlux routine)
                              mLayerRootDensity,           & ! intent(in): root density in each layer (-)
                              scalarTranspireLim,          & ! intent(in): weighted average of the transpiration limiting factor (-)
                              mLayerTranspireLim,          & ! intent(in): transpiration limiting factor in each layer (-)

                              ! output: diagnostic variables for surface runoff
                              xMaxInfilRate,               & ! intent(inout): maximum infiltration rate (m s-1)
                              scalarInfilArea,             & ! intent(inout): fraction of unfrozen area where water can infiltrate (-)
                              scalarFrozenArea,            & ! intent(inout): fraction of area that is considered impermeable due to soil ice (-)
                              scalarSurfaceRunoff,         & ! intent(out):   surface runoff (m s-1)

                              ! output: diagnostic variables for model layers
                              mLayerdTheta_dPsi,           & ! intent(out): derivative in the soil water characteristic w.r.t. psi (m-1)
                              mLayerdPsi_dTheta,           & ! intent(out): derivative in the soil water characteristic w.r.t. theta (m)
                              dHydCond_dMatric,            & ! intent(out): derivative in hydraulic conductivity w.r.t matric head (s-1)

                              ! output: fluxes
                              scalarSurfaceInfiltration,   & ! intent(out): surface infiltration rate (m s-1)
                              iLayerLiqFluxSoil,           & ! intent(out): liquid fluxes at layer interfaces (m s-1)
                              mLayerTranspire,             & ! intent(out): transpiration loss from each soil layer (m s-1)
                              mLayerHydCond,               & ! intent(out): hydraulic conductivity in each soil layer (m s-1)

                              ! output: derivatives in fluxes w.r.t. hydrology state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                              dq_dHydStateAbove,           & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer above (m s-1)
                              dq_dHydStateBelow,           & ! intent(out): derivatives in the flux w.r.t. volumetric liquid water content in the layer below (m s-1)

                              ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
                              dq_dNrgStateAbove,           & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
                              dq_dNrgStateBelow,           & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)

                              ! output: error control
                              err,message)                   ! intent(out): error control
 ! utility modules
 USE soil_utils_module,only:volFracLiq      ! compute volumetric fraction of liquid water
 USE soil_utils_module,only:matricHead      ! compute matric head (m)
 USE soil_utils_module,only:dTheta_dPsi     ! compute derivative of the soil moisture characteristic w.r.t. psi (m-1)
 USE soil_utils_module,only:dPsi_dTheta     ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
 USE soil_utils_module,only:hydCond_psi     ! compute hydraulic conductivity as a function of matric head
 USE soil_utils_module,only:hydCond_liq     ! compute hydraulic conductivity as a function of volumetric liquid water content
 USE soil_utils_module,only:hydCondMP_liq   ! compute hydraulic conductivity of macropores as a function of volumetric liquid water content
 implicit none
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! ***** input variables
 ! ------------------------------------------------------------------------------------------------------------------------------------------------- 
 ! input: model control
 logical(lgt),intent(in)          :: doInfiltrate                  ! flag to compute infiltration
 logical(lgt),intent(in)          :: deriv_desired                ! flag indicating if derivatives are desired
 ! input: model decisions
 integer(i4b),intent(in)          :: ixDerivMethod                ! choice of method used to compute derivative
 integer(i4b),intent(in)          :: ixRichards                   ! choice of the form of Richards' equation
 integer(i4b),intent(in)          :: hc_profile                   ! choice of type of hydraulic conductivity profile
 integer(i4b),intent(in)          :: ixBcUpperSoilHydrology       ! choice of upper boundary condition for soil hydrology
 integer(i4b),intent(in)          :: ixBcLowerSoilHydrology       ! choice of lower boundary condition for soil hydrology
 ! input: trial model state variables
 real(dp),intent(in)              :: mLayerTempTrial(:)           ! temperature in each layer at the current iteration (K)
 real(dp),intent(in)              :: mLayerMatricHeadTrial(:)     ! matric head in each layer at the current iteration (m)
 real(dp),intent(in)              :: mLayerVolFracLiqTrial(:)     ! volumetric fraction of liquid water at the current iteration (-)
 real(dp),intent(in)              :: mLayerVolFracIceTrial(:)     ! volumetric fraction of ice at the current iteration (-)
 ! input: pre-computed derivatves
 real(dp),intent(in)              :: mLayerdTheta_dTk(:)          ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
 real(dp),intent(in)              :: dPsiLiq_dTemp(:)             ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
 ! input: model fluxes
 real(dp),intent(in)              :: scalarCanopyTranspiration    ! canopy transpiration (kg m-2 s-1)
 real(dp),intent(in)              :: scalarGroundEvaporation      ! ground evaporation (kg m-2 s-1)
 real(dp),intent(in)              :: scalarRainPlusMelt           ! rain plus melt (m s-1) 
 ! input: model coordinate variables
 real(dp),intent(in)              :: mLayerDepth(:)               ! depth of the layer (m)
 real(dp),intent(in)              :: mLayerHeight(:)              ! height of the layer mid-point (m)
 real(dp),intent(in)              :: iLayerHeight(0:)             ! height of the layer interfaces (m)
 ! input: diriclet upper boundary conditions
 real(dp),intent(in)              :: upperBoundHead               ! upper boundary condition for matric head (m)
 real(dp),intent(in)              :: upperBoundTheta              ! upper boundary condition for volumetric liquid water content (-)
 ! input: diriclet lower boundary conditions
 real(dp),intent(in)              :: lowerBoundHead               ! lower boundary condition for matric head (m)
 real(dp),intent(in)              :: lowerBoundTheta              ! lower boundary condition for volumetric liquid water content (-)
 ! input: soil parameters
 real(dp),intent(in)              :: vGn_alpha                    ! van Genutchen "alpha" parameter (m-1)
 real(dp),intent(in)              :: vGn_n                        ! van Genutchen "n" parameter (-)
 real(dp),intent(in)              :: vGn_m                        ! van Genutchen "m" parameter (-)
 real(dp),intent(in)              :: mpExp                        ! empirical exponent in macropore flow equation (-)
 real(dp),intent(in)              :: theta_mp                     ! volumetric liquid water content when macropore flow begins (-)
 real(dp),intent(in)              :: theta_sat                    ! soil porosity (-)
 real(dp),intent(in)              :: theta_res                    ! soil residual volumetric water content (-)
 real(dp),intent(in)              :: wettingFrontSuction          ! Green-Ampt wetting front suction (m)
 real(dp),intent(in)              :: fieldCapacity                ! field capacity (-)
 real(dp),intent(in)              :: rootingDepth                 ! rooting depth (m)
 real(dp),intent(in)              :: kAnisotropic                 ! anisotropy factor for lateral hydraulic conductivity (-)
 real(dp),intent(in)              :: zScale_TOPMODEL              ! TOPMODEL scaling factor (m)
 real(dp),intent(in)              :: qSurfScale                   ! scaling factor in the surface runoff parameterization (-)
 real(dp),intent(in)              :: specificYield                ! specific yield (-)
 real(dp),intent(in)              :: specificStorage              ! specific storage coefficient (m-1)
 real(dp),intent(in)              :: f_impede                     ! ice impedence factor (-)
 real(dp),intent(in)              :: soilIceScale                 ! scaling factor for depth of soil ice, used to get frozen fraction (m)
 real(dp),intent(in)              :: soilIceCV                    ! CV of depth of soil ice, used to get frozen fraction (-)
 ! input: saturated hydraulic conductivity
 real(dp),intent(in)              :: mLayerSatHydCondMP(:)        ! saturated hydraulic conductivity of macropores at the mid-point of each layer (m s-1)
 real(dp),intent(in)              :: mLayerSatHydCond(:)          ! saturated hydraulic conductivity at the mid-point of each layer (m s-1)
 real(dp),intent(in)              :: iLayerSatHydCond(0:)         ! saturated hydraulic conductivity at the interface of each layer (m s-1)
 ! input: factors limiting transpiration (from vegFlux routine)
 real(dp),intent(in)              :: mLayerRootDensity(:)         ! root density in each layer (-)
 real(dp),intent(in)              :: scalarTranspireLim           ! weighted average of the transpiration limiting factor (-)
 real(dp),intent(in)              :: mLayerTranspireLim(:)        ! transpiration limiting factor in each layer (-)
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! ***** output variables
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! output: diagnostic variables for surface runoff
 real(dp),intent(inout)           :: xMaxInfilRate                ! maximum infiltration rate (m s-1)
 real(dp),intent(inout)           :: scalarInfilArea              ! fraction of unfrozen area where water can infiltrate (-)
 real(dp),intent(inout)           :: scalarFrozenArea             ! fraction of area that is considered impermeable due to soil ice (-)
 real(dp),intent(out)             :: scalarSurfaceRunoff          ! surface runoff (m s-1)
 ! output: diagnostic variables for each layer
 real(dp),intent(out)             :: mLayerdTheta_dPsi(:)         ! derivative in the soil water characteristic w.r.t. psi (m-1)
 real(dp),intent(out)             :: mLayerdPsi_dTheta(:)         ! derivative in the soil water characteristic w.r.t. theta (m)
 real(dp),intent(out)             :: dHydCond_dMatric(:)          ! derivative in hydraulic conductivity w.r.t matric head (s-1)
 ! output: liquid fluxes
 real(dp),intent(out)             :: scalarSurfaceInfiltration    ! surface infiltration rate (m s-1)
 real(dp),intent(out)             :: iLayerLiqFluxSoil(0:)        ! liquid flux at soil layer interfaces (m s-1)
 real(dp),intent(out)             :: mLayerTranspire(:)           ! transpiration loss from each soil layer (m s-1)
 real(dp),intent(out)             :: mLayerHydCond(:)             ! hydraulic conductivity in each soil layer (m s-1)
 ! output: derivatives in fluxes w.r.t. state variables in the layer above and layer below (m s-1)
 real(dp),intent(out)             :: dq_dHydStateAbove(0:)        ! derivative in the flux in layer interfaces w.r.t. state variables in the layer above
 real(dp),intent(out)             :: dq_dHydStateBelow(0:)        ! derivative in the flux in layer interfaces w.r.t. state variables in the layer below
 ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
 real(dp),intent(out)             :: dq_dNrgStateAbove(0:)        ! derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
 real(dp),intent(out)             :: dq_dNrgStateBelow(0:)        ! derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
 ! output: error control
 integer(i4b),intent(out)         :: err                          ! error code
 character(*),intent(out)         :: message                      ! error message
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! ***** local variables
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables: general
 character(LEN=256)               :: cmessage                     ! error message of downwind routine
 logical(lgt)                     :: desireAnal                   ! flag to identify if analytical derivatives are desired
 integer(i4b)                     :: iLayer,iSoil,jSoil           ! index of soil layer
 ! additional variables to compute numerical derivatives
 integer(i4b)                     :: nFlux                        ! number of flux calculations required (>1 = numerical derivatives with one-sided finite differences)
 integer(i4b)                     :: itry                         ! index of different flux calculations
 integer(i4b),parameter           :: unperturbed=0                ! named variable to identify the case of unperturbed state variables
 integer(i4b),parameter           :: perturbState=1               ! named variable to identify the case where we perturb the state in the current layer
 integer(i4b),parameter           :: perturbStateAbove=2          ! named variable to identify the case where we perturb the state layer above
 integer(i4b),parameter           :: perturbStateBelow=3          ! named variable to identify the case where we perturb the state layer below
 integer(i4b)                     :: ixPerturb                    ! index of element in 2-element vector to perturb
 integer(i4b)                     :: ixOriginal                   ! index of perturbed element in the original vector
 real(dp)                         :: scalarVolFracLiqTrial        ! trial value of volumetric liquid water content (-)
 real(dp)                         :: scalarMatricHeadTrial        ! trial value of matric head (m)
 real(dp)                         :: scalarHydCondTrial           ! trial value of hydraulic conductivity (m s-1)
 real(dp)                         :: scalarHydCondMicro           ! trial value of hydraulic conductivity of micropores (m s-1)
 real(dp)                         :: scalarHydCondMacro           ! trial value of hydraulic conductivity of macropores (m s-1)
 real(dp)                         :: scalarFlux                   ! vertical flux (m s-1)
 real(dp)                         :: scalarFlux_dState            ! vertical flux with perturbation to the current state (m s-1)
 real(dp)                         :: scalarFlux_dStateAbove       ! vertical flux with perturbation to the state above (m s-1)
 real(dp)                         :: scalarFlux_dStateBelow       ! vertical flux with perturbation to the state below (m s-1)
 real(dp)                         :: scalarFluxUpper              ! flux at the top of the layer
 real(dp)                         :: scalarFluxLower              ! flux at the bottom of the layer
 real(dp)                         :: scalarFluxUpper_dx           ! flux at the top of the layer after perturbation
 real(dp)                         :: scalarFluxLower_dx           ! flux at the bottom of the layer after perturbation
 ! transpiration sink term
 real(dp),dimension(nSoil)        :: mLayerTranspireFrac          ! fraction of transpiration allocated to each soil layer (-)
 ! diagnostic variables 
 real(dp),dimension(nSoil)        :: iceImpedeFac                 ! ice impedence factor at layer mid-points (-)
 real(dp),dimension(nSoil)        :: mLayerDiffuse                ! diffusivity at layer mid-point (m2 s-1)
 real(dp),dimension(nSoil)        :: dHydCond_dVolLiq             ! derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
 real(dp),dimension(nSoil)        :: dDiffuse_dVolLiq             ! derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
 real(dp),dimension(nSoil)        :: dHydCond_dTemp               ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
 real(dp),dimension(0:nSoil)      :: iLayerHydCond                ! hydraulic conductivity at layer interface (m s-1)
 real(dp),dimension(0:nSoil)      :: iLayerDiffuse                ! diffusivity at layer interface (m2 s-1)
 ! compute surface flux
 integer(i4b)                     :: nRoots                       ! number of soil layers with roots
 integer(i4b)                     :: ixIce                        ! index of the lowest soil layer that contains ice
 ! compute fluxes and derivatives at layer interfaces 
 real(dp),dimension(2)            :: vectorVolFracLiqTrial        ! trial value of volumetric liquid water content (-)
 real(dp),dimension(2)            :: vectorMatricHeadTrial        ! trial value of matric head (m)
 real(dp),dimension(2)            :: vectorHydCondTrial           ! trial value of hydraulic conductivity (m s-1)
 real(dp),dimension(2)            :: vectorDiffuseTrial           ! trial value of hydraulic diffusivity (m2 s-1)
 real(dp)                         :: scalardPsi_dTheta             ! derivative in soil water characteristix, used for perturbations when computing numerical derivatives
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='soilLiqFlx_muster/'

 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! preliminaries
 ! -------------------------------------------------------------------------------------------------------------------------------------------------

 ! define the pethod to compute derivatives
 !print*, 'numerical derivatives = ', (ixDerivMethod==numerical)

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
  nFlux=0  ! compute analytical derivatives
 endif

 ! identify the number of layers that contain roots
 nRoots = count(iLayerHeight(0:nSoil-1) < rootingDepth)

 ! identify lowest soil layer with ice
 ! NOTE: cannot use count because there may be an unfrozen wedge
 ixIce = 0  ! initialize the index of the ice layer (0 means no ice in the soil profile)
 do iLayer=1,nSoil ! (loop through soil layers)
  if(mLayerVolFracIceTrial(iLayer) > verySmall) ixIce = iLayer
 end do
 !if(ixIce==nSoil)then; err=20; message=trim(message)//'ice extends to the bottom of the soil profile'; return; endif

 ! *************************************************************************************************************************************************
 ! *************************************************************************************************************************************************

 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! compute the transpiration sink term
 ! -------------------------------------------------------------------------------------------------------------------------------------------------

 ! compute the fraction of transpiration loss from each soil layer
 if(scalarTranspireLim > epsilon(scalarTranspireLim))then ! (transpiration may be non-zero even if the soil moisture limiting factor is zero)
  mLayerTranspireFrac(:) = mLayerRootDensity(:)*mLayerTranspireLim(:)/scalarTranspireLim
 else ! (possible for there to be non-zero conductance and therefore transpiration in this case)
  mLayerTranspireFrac(:) = mLayerRootDensity(:)
 endif

 ! compute transpiration loss from each soil layer (kg m-2 s-1 --> m s-1)
 mLayerTranspire        = mLayerTranspireFrac(:)*scalarCanopyTranspiration/iden_water
 ! (special case of prescribed head -- no transpiration)
 if(ixBcUpperSoilHydrology==prescribedHead) mLayerTranspire(:) = 0._dp
 !print*, trim(message)//'mLayerTranspire = ', mLayerTranspire
 !print*, trim(message)//'scalarCanopyTranspiration = ', scalarCanopyTranspiration


 ! *************************************************************************************************************************************************
 ! *************************************************************************************************************************************************

 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! compute diagnostic variables at the nodes throughout the soil profile
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 do iSoil=1,nSoil ! (loop through soil layers)
  call diagv_node(&
                  ! input: model control
                  desireAnal,                      & ! intent(in): flag indicating if derivatives are desired
                  ixRichards,                      & ! intent(in): index defining the option for Richards' equation (moisture or mixdform)
                  ! input: state variables
                  mLayerTempTrial(iSoil),          & ! intent(in): temperature (K)
                  mLayerMatricHeadTrial(iSoil),    & ! intent(in): matric head in each layer (m)
                  mLayerVolFracLiqTrial(iSoil),    & ! intent(in): volumetric liquid water content in each soil layer (-)
                  mLayerVolFracIceTrial(iSoil),    & ! intent(in): volumetric ice content in each soil layer (-)
                  ! input: pre-computed deriavatives
                  mLayerdTheta_dTk(iSoil),         & ! intent(in): derivative in volumetric liquid water content w.r.t. temperature (K-1)
                  dPsiLiq_dTemp(iSoil),            & ! intent(in): derivative in liquid water matric potential w.r.t. temperature (m K-1)
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
                  mLayerdPsi_dTheta(iSoil),        & ! intent(out): derivative in the soil water characteristic
                  mLayerdTheta_dPsi(iSoil),        & ! intent(out): derivative in the soil water characteristic
                  ! output: transmittance
                  mLayerHydCond(iSoil),            & ! intent(out): hydraulic conductivity at layer mid-points (m s-1)
                  mLayerDiffuse(iSoil),            & ! intent(out): diffusivity at layer mid-points (m2 s-1)
                  iceImpedeFac(iSoil),             & ! intent(out): ice impedence factor in each layer (-)
                  ! output: transmittance derivatives
                  dHydCond_dVolLiq(iSoil),         & ! intent(out): derivative in hydraulic conductivity w.r.t volumetric liquid water content (m s-1)
                  dDiffuse_dVolLiq(iSoil),         & ! intent(out): derivative in hydraulic diffusivity w.r.t volumetric liquid water content (m2 s-1)
                  dHydCond_dMatric(iSoil),         & ! intent(out): derivative in hydraulic conductivity w.r.t matric head (m s-1)
                  dHydCond_dTemp(iSoil),           & ! intent(out): derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
                  ! output: error control
                  err,cmessage)                      ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 end do  ! (looping through soil layers)

 ! *************************************************************************************************************************************************
 ! *************************************************************************************************************************************************

 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! compute infiltraton at the surface and its derivative w.r.t. mass in the upper soil layer
 ! -------------------------------------------------------------------------------------------------------------------------------------------------

 ! set derivative w.r.t. state above to zero (does not exist)
 dq_dHydStateAbove(0) = 0._dp
 dq_dNrgStateAbove(0) = 0._dp

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
    end select ! (form of Richards' equation
   ! check for an unknown perturbation 
   case default; err=10; message=trim(message)//"unknown perturbation"; return

  end select ! (type of perturbation)

  ! =====
  ! compute surface flux and its derivative...
  ! ==========================================
  call surfaceFlx(&
                  ! input: model control
                  doInfiltrate,                       & ! intent(in): flag indicating if desire to compute infiltration
                  desireAnal,                         & ! intent(in): flag indicating if derivatives are desired
                  ixRichards,                         & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                  ixBcUpperSoilHydrology,             & ! intent(in): index defining the type of boundary conditions (neumann or diriclet)
                  nRoots,                             & ! intent(in): number of layers that contain roots
                  ixIce,                              & ! intent(in): index of lowest ice layer
                  ! input: state variables
                  scalarMatricHeadTrial,              & ! intent(in): matric head in the upper-most soil layer (m)
                  scalarVolFracLiqTrial,              & ! intent(in): volumetric liquid water content the upper-most soil layer (-)
                  mLayerVolFracLiqTrial,              & ! intent(in): volumetric liquid water content in each soil layer (-)
                  mLayerVolFracIceTrial,              & ! intent(in): volumetric ice content in each soil layer (-)
                  ! input: depth of upper-most soil layer (m)
                  mLayerDepth,                        & ! intent(in): depth of each soil layer (m)
                  iLayerHeight,                       & ! intent(in): height at the interface of each layer (m)
                  ! input: boundary conditions
                  upperBoundHead,                     & ! intent(in): upper boundary condition (m)
                  upperBoundTheta,                    & ! intent(in): upper boundary condition (-)
                  ! input: flux at the upper boundary
                  scalarRainPlusMelt,                 & ! intent(in): rain plus melt (m s-1)
                  ! input: derivative in soil water characteristix
                  mLayerdPsi_dTheta(1),               & ! intent(in): derivative of the soil moisture characteristic w.r.t. theta (m)
                  mLayerdTheta_dPsi(1),               & ! intent(in): derivative of the soil moisture characteristic w.r.t. psi (m-1)
                  ! input: transmittance
                  iLayerSatHydCond(0),                & ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
                  dHydCond_dTemp(1),                  & ! intent(in): derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)     
                  iceImpedeFac(1),                    & ! intent(in): ice impedence factor in the upper-most soil layer (-)
                  ! input: soil parameters
                  vGn_alpha,                          & ! intent(in): van Genutchen "alpha" parameter (m-1)
                  vGn_n,                              & ! intent(in): van Genutchen "n" parameter (-)
                  VGn_m,                              & ! intent(in): van Genutchen "m" parameter (-)
                  theta_sat,                          & ! intent(in): soil porosity (-)
                  theta_res,                          & ! intent(in): soil residual volumetric water content (-)
                  qSurfScale,                         & ! intent(in): scaling factor in the surface runoff parameterization (-)
                  zScale_TOPMODEL,                    & ! intent(in): scaling factor used to describe decrease in hydraulic conductivity with depth (m)
                  rootingDepth,                       & ! intent(in): rooting depth (m)
                  wettingFrontSuction,                & ! intent(in): Green-Ampt wetting front suction (m)
                  soilIceScale,                       & ! intent(in): soil ice scaling factor in Gamma distribution used to define frozen area (m)
                  soilIceCV,                          & ! intent(in): soil ice CV in Gamma distribution used to define frozen area (-)
                  ! input-output: hydraulic conductivity and diffusivity at the surface
                  iLayerHydCond(0),                   & ! intent(inout): hydraulic conductivity at the surface (m s-1)
                  iLayerDiffuse(0),                   & ! intent(inout): hydraulic diffusivity at the surface (m2 s-1)
                  ! input-output: fluxes at layer interfaces and surface runoff
                  xMaxInfilRate,                      & ! intent(inout): maximum infiltration rate (m s-1)
                  scalarInfilArea,                    & ! intent(inout): fraction of unfrozen area where water can infiltrate (-)
                  scalarFrozenArea,                   & ! intent(inout): fraction of area that is considered impermeable due to soil ice (-)
                  scalarSurfaceRunoff,                & ! intent(out): surface runoff (m s-1)
                  scalarSurfaceInfiltration,          & ! intent(out): surface infiltration (m s-1)
                  ! input-output: deriavtives in surface infiltration w.r.t. volumetric liquid water (m s-1) and matric head (s-1) in the upper-most soil layer
                  dq_dHydStateBelow(0),               & ! intent(inout): derivative in surface infiltration w.r.t. hydrology state variable in the upper-most soil layer (m s-1 or s-1)
                  dq_dNrgStateBelow(0),               & ! intent(out):   derivative in surface infiltration w.r.t. energy state variable in the upper-most soil layer (m s-1 K-1)
                  ! output: error control
                  err,cmessage)                         ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! include base soil evaporation as the upper boundary flux
  iLayerLiqFluxSoil(0) = scalarGroundEvaporation/iden_water + scalarSurfaceInfiltration

  ! get copies of surface flux to compute numerical derivatives
  if(deriv_desired .and. ixDerivMethod==numerical)then
   select case(itry)
    case(unperturbed);       scalarFlux             = iLayerLiqFluxSoil(0)
    case(perturbStateBelow); scalarFlux_dStateBelow = iLayerLiqFluxSoil(0)
    case default; err=10; message=trim(message)//"unknown perturbation"; return
   end select
  endif

  !write(*,'(a,1x,10(f30.15))') 'scalarRainPlusMelt, scalarSurfaceInfiltration = ', scalarRainPlusMelt, scalarSurfaceInfiltration

 end do  ! (looping through different flux calculations -- one or multiple calls depending if desire for numerical or analytical derivatives)

 ! compute numerical derivatives
 if(deriv_desired .and. ixDerivMethod==numerical)then
  dq_dHydStateBelow(0) = (scalarFlux_dStateBelow - scalarFlux)/dx ! change in surface flux w.r.t. change in the soil moisture in the top soil layer (m s-1)
 endif
 !print*, 'scalarSurfaceInfiltration, iLayerLiqFluxSoil(0) = ', scalarSurfaceInfiltration, iLayerLiqFluxSoil(0)
 !print*, '(ixDerivMethod==numerical), dq_dHydStateBelow(0) = ', (ixDerivMethod==numerical), dq_dHydStateBelow(0)
 !pause

 ! *************************************************************************************************************************************************
 ! *************************************************************************************************************************************************

 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! * compute fluxes and derivatives at layer interfaces...
 ! -------------------------------------------------------------------------------------------------------------------------------------------------
 ! loop through soil layers
 do iLayer=1,nSoil-1

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
                   ! input: temperature derivatives
                   dPsiLiq_dTemp(iLayer:iLayer+1),     & ! intent(in): derivative in liquid water matric potential w.r.t. temperature (m K-1)
                   dHydCond_dTemp(iLayer:iLayer+1),    & ! intent(in): derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)     
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
                   dq_dHydStateAbove(iLayer),          & ! intent(out): derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer above (m s-1 or s-1)
                   dq_dHydStateBelow(iLayer),          & ! intent(out): derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer below (m s-1 or s-1) 
                   ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
                   dq_dNrgStateAbove(iLayer),          & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
                   dq_dNrgStateBelow(iLayer),          & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
                   ! output: error control
                   err,cmessage)                         ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

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
   dq_dHydStateAbove(iLayer) = (scalarFlux_dStateAbove - scalarFlux)/dx    ! change in drainage flux w.r.t. change in the state in the layer below (m s-1 or s-1)
   dq_dHydStateBelow(iLayer) = (scalarFlux_dStateBelow - scalarFlux)/dx    ! change in drainage flux w.r.t. change in the state in the layer below (m s-1 or s-1)
  endif
  
  ! check
  !if(iLayer==6) write(*,'(a,i4,1x,10(e25.15,1x))') 'iLayer, vectorMatricHeadTrial, iLayerHydCond(iLayer), iLayerLiqFluxSoil(iLayer) = ',&
  !                                                  iLayer, vectorMatricHeadTrial, iLayerHydCond(iLayer), iLayerLiqFluxSoil(iLayer)
  !if(iLayer==1) write(*,'(a,i4,1x,L1,1x,2(e15.5,1x))') 'iLayer, (ixDerivMethod==numerical), dq_dHydStateBelow(iLayer-1), dq_dHydStateAbove(iLayer) = ', &
  !                                                      iLayer, (ixDerivMethod==numerical), dq_dHydStateBelow(iLayer-1), dq_dHydStateAbove(iLayer)
  !pause

 end do  ! (looping through soil layers)

 ! add infiltration to the upper-most unfrozen layer
 ! NOTE: this is done here rather than in surface runoff
 !iLayerLiqFluxSoil(ixIce) = iLayerLiqFluxSoil(ixIce) + scalarSurfaceInfiltration

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
    scalarVolFracLiqTrial   = mLayerVolFracLiqTrial(nSoil)
    scalarMatricHeadTrial   = mLayerMatricHeadTrial(nSoil)

   ! perturb soil state (one-sided finite differences)
   case(perturbStateAbove)
    select case(ixRichards)  ! (perturbation depends on the form of Richards' equation)
     case(moisture)
      scalarVolFracLiqTrial = mLayerVolFracLiqTrial(nSoil) + dx
      scalarMatricHeadTrial = mLayerMatricHeadTrial(nSoil)
     case(mixdform)
      scalarVolFracLiqTrial = mLayerVolFracLiqTrial(nSoil)
      scalarMatricHeadTrial = mLayerMatricHeadTrial(nSoil) + dx
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
     case(moisture); scalarHydCondTrial = hydCond_liq(scalarVolFracLiqTrial,mLayerSatHydCond(nSoil),theta_res,theta_sat,vGn_m) * iceImpedeFac(nSoil)
     case(mixdform); scalarHydCondTrial = hydCond_psi(scalarMatricHeadTrial,mLayerSatHydCond(nSoil),vGn_alpha,vGn_n,vGn_m) * iceImpedeFac(nSoil)
    end select

   ! (use un-perturbed value)
   case default
    scalarHydCondTrial = mLayerHydCond(nSoil)        ! hydraulic conductivity at the mid-point of the lowest unsaturated soil layer (m s-1)

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
                  mLayerDepth(nSoil),              & ! intent(in): depth of the lowest unsaturated soil layer (m)
                  mLayerHeight(nSoil),             & ! intent(in): height of the lowest unsaturated soil node (m)
                  ! input: boundary conditions
                  lowerBoundHead,                  & ! intent(in): lower boundary condition (m)
                  lowerBoundTheta,                 & ! intent(in): lower boundary condition (-)
                  ! input: derivative in the soil water characteristic
                  mLayerdPsi_dTheta(nSoil),        & ! intent(in): derivative in the soil water characteristic
                  mLayerdTheta_dPsi(nSoil),        & ! intent(in): derivative in the soil water characteristic
                  ! input: transmittance
                  iLayerSatHydCond(0),             & ! intent(in): saturated hydraulic conductivity at the surface (m s-1)
                  iLayerSatHydCond(nSoil),         & ! intent(in): saturated hydraulic conductivity at the bottom of the unsaturated zone (m s-1)
                  scalarHydCondTrial,              & ! intent(in): hydraulic conductivity at the node itself (m s-1)
                  iceImpedeFac(nSoil),             & ! intent(in): ice impedence factor in the lower-most soil layer (-)
                  ! input: transmittance derivatives
                  dHydCond_dVolLiq(nSoil),         & ! intent(in): derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
                  dHydCond_dMatric(nSoil),         & ! intent(in): derivative in hydraulic conductivity w.r.t. matric head (s-1)
                  dHydCond_dTemp(nSoil),           & ! intent(in): derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
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
                  iLayerHydCond(nSoil),            & ! intent(out): hydraulic conductivity at the bottom of the unsatuarted zone (m s-1)
                  iLayerDiffuse(nSoil),            & ! intent(out): hydraulic diffusivity at the bottom of the unsatuarted zone (m2 s-1)
                  ! output: drainage flux
                  iLayerLiqFluxSoil(nSoil),        & ! intent(out): drainage flux (m s-1)
                  ! output: derivatives in drainage flux
                  dq_dHydStateAbove(nSoil),        & ! intent(out): change in drainage flux w.r.t. change in hydrology state in lowest unsaturated node (m s-1 or s-1)
                  dq_dNrgStateAbove(nSoil),        & ! intent(out): change in drainage flux w.r.t. change in energy state in lowest unsaturated node (m s-1 or s-1)
                  ! output: error control
                  err,cmessage)                ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 
  ! get copies of drainage flux to compute derivatives
  if(deriv_desired .and. ixDerivMethod==numerical)then
   select case(itry)
    case(unperturbed);       scalarFlux             = iLayerLiqFluxSoil(nSoil)
    case(perturbStateAbove); scalarFlux_dStateAbove = iLayerLiqFluxSoil(nSoil)
    case(perturbStateBelow); err=20; message=trim(message)//'lower state should never be perturbed when computing drainage do not expect to get here'; return
    case default; err=10; message=trim(message)//"unknown perturbation"; return
   end select
  endif

 end do  ! (looping through different flux calculations -- one or multiple calls depending if desire for numerical or analytical derivatives)
 
 ! compute numerical derivatives
 ! NOTE: drainage derivatives w.r.t. state below are *actually* w.r.t. water table depth, so need to be corrected for aquifer storage
 !       (note also negative sign to account for inverse relationship between water table depth and aquifer storage)
 if(deriv_desired .and. ixDerivMethod==numerical)then
  dq_dHydStateAbove(nSoil) = (scalarFlux_dStateAbove - scalarFlux)/dx    ! change in drainage flux w.r.t. change in state in lowest unsaturated node (m s-1 or s-1)
 endif

 ! no dependence on the aquifer for drainage
 dq_dHydStateBelow(nSoil) = 0._dp  ! keep this here in case we want to couple some day....
 dq_dNrgStateBelow(nSoil) = 0._dp  ! keep this here in case we want to couple some day....

 ! print drainage
 !print*, 'iLayerLiqFluxSoil(nSoil) = ', iLayerLiqFluxSoil(nSoil)
 
 ! end of drainage section


 ! *****************************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************

 
 end subroutine soilLiqFlx_muster

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
                       scalarTempTrial,       & ! intent(in): temperature (K)
                       scalarMatricHeadTrial, & ! intent(in): matric head in a given layer (m)
                       scalarVolFracLiqTrial, & ! intent(in): volumetric liquid water content in a given soil layer (-)
                       scalarVolFracIceTrial, & ! intent(in): volumetric ice content in a given soil layer (-)
                       ! input: pre-computed deriavatives
                       dTheta_dTk,            & ! intent(in): derivative in volumetric liquid water content w.r.t. temperature (K-1)
                       dPsiLiq_dTemp,         & ! intent(in): derivative in liquid water matric potential w.r.t. temperature (m K-1)
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
                       dHydCond_dTemp,        & ! intent(out): derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
                       ! output: error control
                       err,message)             ! intent(out): error control
 USE soil_utils_module,only:iceImpede           ! compute the ice impedence factor
 USE soil_utils_module,only:volFracLiq          ! compute volumetric fraction of liquid water as a function of matric head
 USE soil_utils_module,only:matricHead          ! compute matric head (m)
 USE soil_utils_module,only:hydCond_psi         ! compute hydraulic conductivity as a function of matric head
 USE soil_utils_module,only:hydCond_liq         ! compute hydraulic conductivity as a function of volumetric liquid water content
 USE soil_utils_module,only:hydCondMP_liq       ! compute hydraulic conductivity of macropores as a function of volumetric liquid water content
 USE soil_utils_module,only:dTheta_dPsi         ! compute derivative of the soil moisture characteristic w.r.t. psi (m-1)
 USE soil_utils_module,only:dPsi_dTheta         ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
 USE soil_utils_module,only:dPsi_dTheta2        ! compute derivative in dPsi_dTheta (m)
 USE soil_utils_module,only:dHydCond_dLiq       ! compute derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
 USE soil_utils_module,only:dHydCond_dPsi       ! compute derivative in hydraulic conductivity w.r.t. matric head (s-1)
 USE soil_utils_module,only:dIceImpede_dTemp    ! compute the derivative in the ice impedance factor w.r.t. temperature (K-1)
 ! compute hydraulic transmittance and derivatives for all layers
 implicit none
 ! input: model control
 logical(lgt),intent(in)       :: deriv_desired             ! flag indicating if derivatives are desired
 integer(i4b),intent(in)       :: ixRichards                ! index defining the option for Richards' equation (moisture or mixdform)
 ! input: state and diagnostic variables
 real(dp),intent(in)           :: scalarTempTrial           ! temperature in each layer (K)
 real(dp),intent(in)           :: scalarMatricHeadTrial     ! matric head in each layer (m)
 real(dp),intent(in)           :: scalarVolFracLiqTrial     ! volumetric fraction of liquid water in a given layer (-)
 real(dp),intent(in)           :: scalarVolFracIceTrial     ! volumetric fraction of ice in a given layer (-)
 ! input: pre-computed deriavatives
 real(dp),intent(in)           :: dTheta_dTk                ! derivative in volumetric liquid water content w.r.t. temperature (K-1)
 real(dp),intent(in)           :: dPsiLiq_dTemp             ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
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
 real(dp),intent(out)          :: dHydCond_dTemp            ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                       ! error code
 character(*),intent(out)      :: message                   ! error message
 ! local variables
 real(dp)                      :: localVolFracLiq           ! local volumetric fraction of liquid water
 real(dp)                      :: scalarHydCondMP           ! hydraulic conductivity of macropores at layer mid-points (m s-1)
 real(dp)                      :: dIceImpede_dT             ! derivative in ice impedance factor w.r.t. temperature (K-1)
 real(dp)                      :: dHydCondMacro_dVolLiq     ! derivative in hydraulic conductivity of macropores w.r.t volumetric liquid water content (m s-1) 
 real(dp)                      :: dHydCondMacro_dMatric     ! derivative in hydraulic conductivity of macropores w.r.t matric head (s-1)
 real(dp)                      :: dHydCondMicro_dMatric     ! derivative in hydraulic conductivity of micropores w.r.t matric head (s-1)
 real(dp)                      :: dHydCondMicro_dTemp       ! derivative in hydraulic conductivity of micropores w.r.t temperature (m s-1 K-1)
 real(dp)                      :: dPsi_dTheta2a             ! derivative in dPsi_dTheta (analytical)
 real(dp)                      :: dIceImpede_dLiq           ! derivative in ice impedence factor w.r.t. volumetric liquid water content (-)
 real(dp)                      :: hydCond_noIce             ! hydraulic conductivity in the absence of ice (m s-1)
 real(dp)                      :: dK_dLiq__noIce            ! derivative in hydraulic conductivity w.r.t volumetric liquid water content, in the absence of ice (m s-1)
 real(dp)                      :: dK_dPsi__noIce            ! derivative in hydraulic conductivity w.r.t matric head, in the absence of ice (s-1)
 real(dp)                      :: relSatMP                  ! relative saturation of macropores (-)
 !real(dp)                      :: xConst,vTheta,volLiq,volIce,x1,x2,d1,d2,effSat,psiLiq,dEff,dPsi,hydCon,hydIce   ! test derivative
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
   !scalardPsi_dTheta = valueMissing  ! (deliberately cause problems if this is ever used)
   scalardPsi_dTheta = dPsi_dTheta(scalarvolFracLiqTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
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
   localVolFracLiq = volFracLiq(scalarMatricHeadTrial,vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m) 
   scalarHydCondMP = hydCondMP_liq(localVolFracLiq,theta_sat,theta_mp,mpExp,scalarSatHydCondMP,scalarSatHydCond)
   scalarHydCond   = hydCond_noIce*iceImpedeFac + scalarHydCondMP
   ! compute derivative in hydraulic conductivity (m s-1)
   if(deriv_desired)then
    ! (compute derivative for macropores)
    if(localVolFracLiq > theta_mp)then
     relSatMP              = (localVolFracLiq - theta_mp)/(theta_sat - theta_mp) 
     dHydCondMacro_dVolLiq = ((scalarSatHydCondMP - scalarSatHydCond)/(theta_sat - theta_mp))*mpExp*(relSatMP**(mpExp - 1._dp))
     dHydCondMacro_dMatric = scalardTheta_dPsi*dHydCondMacro_dVolLiq
    else
     dHydCondMacro_dVolLiq = 0._dp
     dHydCondMacro_dMatric = 0._dp
    endif
    ! (compute derivatives for micropores)
    if(scalarVolFracIceTrial > verySmall)then
     dK_dPsi__noIce        = dHydCond_dPsi(scalarMatricHeadTrial,scalarSatHydCond,vGn_alpha,vGn_n,vGn_m,.true.)  ! analytical
     dHydCondMicro_dTemp   = dPsiLiq_dTemp*dK_dPsi__noIce  ! m s-1 K-1
     dHydCondMicro_dMatric = hydCond_noIce*dIceImpede_dLiq*scalardTheta_dPsi + dK_dPsi__noIce*iceImpedeFac
    else
     dHydCondMicro_dTemp   = 0._dp
     dHydCondMicro_dMatric = dHydCond_dPsi(scalarMatricHeadTrial,scalarSatHydCond,vGn_alpha,vGn_n,vGn_m,.true.)
    endif
    ! (combine derivatives)
    dHydCond_dMatric = dHydCondMicro_dMatric + dHydCondMacro_dMatric
    ! (compute analytical derivative for change in ice impedance factor w.r.t. temperature)
    call dIceImpede_dTemp(scalarVolFracIceTrial, & ! intent(in): trial value of volumetric ice content (-)
                          dTheta_dTk,            & ! intent(in): derivative in volumetric liquid water content w.r.t. temperature (K-1)
                          f_impede,              & ! intent(in): ice impedance parameter (-)
                          dIceImpede_dT          ) ! intent(out): derivative in ice impedance factor w.r.t. temperature (K-1)
    ! (compute derivative in hydraulic conductivity w.r.t. temperature)
    dHydCond_dTemp = hydCond_noIce*dIceImpede_dT + dHydCondMicro_dTemp*iceImpedeFac
    ! (test derivative)
    !xConst = LH_fus/(gravity*Tfreeze)                            ! m K-1 (NOTE: J = kg m2 s-2)     
    !vTheta = scalarVolFracIceTrial + scalarVolFracLiqTrial
    !volLiq = volFracLiq(xConst*(scalarTempTrial+dx - Tfreeze),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
    !volIce = vTheta - volLiq
    !effSat = (volLiq - theta_res) / (theta_sat - volIce - theta_res)
    !psiLiq = matricHead(effSat,vGn_alpha,0._dp,1._dp,vGn_n,vGn_m)  ! use effective saturation, so theta_res=0 and theta_sat=1
    !hydCon = hydCond_psi(psiLiq,scalarSatHydCond,vGn_alpha,vGn_n,vGn_m)
    !call iceImpede(volIce,volLiq,theta_sat,f_impede,deriv_desired,iceImpedeFac,dIceImpede_dLiq)
    !hydIce = hydCon*iceImpedeFac
    !print*, 'test derivative: ', (psiLiq - scalarMatricHeadTrial)/dx, dPsiLiq_dTemp 
    !print*, 'test derivative: ', (hydCon - hydCond_noIce)/dx, dHydCondMicro_dTemp
    !print*, 'test derivative: ', (hydIce - scalarHydCond)/dx, dHydCond_dTemp 
    !pause
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
                       doInfiltration,            & ! intent(in): flag indicating if desire to compute infiltration
                       deriv_desired,             & ! intent(in): flag indicating if derivatives are desired
                       ixRichards,                & ! intent(in): index defining the form of Richards' equation (moisture or mixdform)
                       bc_upper,                  & ! intent(in): index defining the type of boundary conditions (neumann or diriclet)
                       nRoots,                    & ! intent(in): number of layers that contain roots
                       ixIce,                     & ! intent(in): index of lowest ice layer
                       ! input: state variables
                       scalarMatricHead,          & ! intent(in): matric head in the upper-most soil layer (m)
                       scalarVolFracLiq,          & ! intent(in): volumetric liquid water content in the upper-most soil layer (-)
                       mLayerVolFracLiq,          & ! intent(in): volumetric liquid water content in each soil layer (-)
                       mLayerVolFracIce,          & ! intent(in): volumetric ice content in each soil layer (-)
                       ! input: depth of upper-most soil layer (m)
                       mLayerDepth,               & ! intent(in): depth of each soil layer (m)
                       iLayerHeight,              & ! intent(in): height at the interface of each layer (m)
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
                       dHydCond_dTemp,            & ! intent(in): derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
                       iceImpedeFac,              & ! intent(in): ice impedence factor in the upper-most soil layer (-)
                       ! input: soil parameters
                       vGn_alpha,                 & ! intent(in): van Genutchen "alpha" parameter (m-1)
                       vGn_n,                     & ! intent(in): van Genutchen "n" parameter (-)
                       VGn_m,                     & ! intent(in): van Genutchen "m" parameter (-)
                       theta_sat,                 & ! intent(in): soil porosity (-)
                       theta_res,                 & ! intent(in): soil residual volumetric water content (-)
                       qSurfScale,                & ! intent(in): scaling factor in the surface runoff parameterization (-)
                       zScale_TOPMODEL,           & ! intent(in): scaling factor used to describe decrease in hydraulic conductivity with depth (m)
                       rootingDepth,              & ! intent(in): rooting depth (m)
                       wettingFrontSuction,       & ! intent(in): Green-Ampt wetting front suction (m)
                       soilIceScale,              & ! intent(in): soil ice scaling factor in Gamma distribution used to define frozen area (m)
                       soilIceCV,                 & ! intent(in): soil ice CV in Gamma distribution used to define frozen area (-)
                       ! input-output: hydraulic conductivity and diffusivity at the surface
                       surfaceHydCond,            & ! intent(inout): hydraulic conductivity at the surface (m s-1)
                       surfaceDiffuse,            & ! intent(inout): hydraulic diffusivity at the surface (m2 s-1)
                       ! input-output: fluxes at layer interfaces and surface runoff
                       xMaxInfilRate,             & ! intent(inout): maximum infiltration rate (m s-1)
                       scalarInfilArea,           & ! intent(inout): fraction of unfrozen area where water can infiltrate (-)
                       scalarFrozenArea,          & ! intent(inout): fraction of area that is considered impermeable due to soil ice (-)
                       scalarSurfaceRunoff,       & ! intent(out): surface runoff (m s-1)
                       scalarSurfaceInfiltration, & ! intent(out): surface infiltration (m s-1)
                       ! input-output: deriavtives in surface infiltration w.r.t. volumetric liquid water (m s-1) and matric head (s-1) in the upper-most soil layer
                       dq_dHydState,              & ! intent(inout): derivative in surface infiltration w.r.t. state variable in the upper-most soil layer (m s-1 or s-1)
                       dq_dNrgState,              & ! intent(out):   derivative in surface infiltration w.r.t. energy state variable in the upper-most soil layer (m s-1 K-1)
                       ! output: error control
                       err,message)                 ! intent(out): error control
 USE soil_utils_module,only:volFracLiq      ! compute volumetric fraction of liquid water as a function of matric head (-)
 USE soil_utils_module,only:hydCond_psi     ! compute hydraulic conductivity as a function of matric head (m s-1)
 USE soil_utils_module,only:hydCond_liq     ! compute hydraulic conductivity as a function of volumetric liquid water content (m s-1)
 USE soil_utils_module,only:dPsi_dTheta     ! compute derivative of the soil moisture characteristic w.r.t. theta (m)
 USE soil_utils_module,only:gammp           ! compute the cumulative probabilty based on the Gamma distribution
 ! compute infiltraton at the surface and its derivative w.r.t. mass in the upper soil layer
 implicit none
 ! -----------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 logical(lgt),intent(in)       :: doInfiltration            ! flag indicating if desire to compute infiltration
 logical(lgt),intent(in)       :: deriv_desired             ! flag to indicate if derivatives are desired
 integer(i4b),intent(in)       :: bc_upper                  ! index defining the type of boundary conditions
 integer(i4b),intent(in)       :: ixRichards                ! index defining the option for Richards' equation (moisture or mixdform)
 integer(i4b),intent(in)       :: nRoots                    ! number of layers that contain roots
 integer(i4b),intent(in)       :: ixIce                     ! index of lowest ice layer
 ! input: state and diagnostic variables
 real(dp),intent(in)           :: scalarMatricHead          ! matric head in the upper-most soil layer (m)
 real(dp),intent(in)           :: scalarVolFracLiq          ! volumetric liquid water content in the upper-most soil layer (-)
 real(dp),intent(in)           :: mLayerVolFracLiq(:)       ! volumetric liquid water content in each soil layer (-)
 real(dp),intent(in)           :: mLayerVolFracIce(:)       ! volumetric ice content in each soil layer (-)
 ! input: depth of upper-most soil layer (m)
 real(dp),intent(in)           :: mLayerDepth(:)            ! depth of upper-most soil layer (m)
 real(dp),intent(in)           :: iLayerHeight(0:)          ! height at the interface of each layer (m)
 ! input: diriclet boundary conditions
 real(dp),intent(in)           :: upperBoundHead            ! upper boundary condition for matric head (m)
 real(dp),intent(in)           :: upperBoundTheta           ! upper boundary condition for volumetric liquid water content (-)
 ! input: flux at the upper boundary
 real(dp),intent(in)           :: scalarRainPlusMelt        ! rain plus melt, used as input to the soil zone before computing surface runoff (m s-1)
 ! input: derivative in soil water characteristix
 real(dp),intent(in)           :: scalardPsi_dTheta         ! derivative of the soil moisture characteristic w.r.t. theta (m)
 real(dp),intent(in)           :: scalardTheta_dPsi         ! derivative of the soil moisture characteristic w.r.t. psi (m-1)
 ! input: transmittance
 real(dp),intent(in)           :: surfaceSatHydCond         ! saturated hydraulic conductivity at the surface (m s-1)
 real(dp),intent(in)           :: dHydCond_dTemp            ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
 real(dp),intent(in)           :: iceImpedeFac              ! ice impedence factor in the upper-most soil layer (-)
 ! input: soil parameters
 real(dp),intent(in)           :: vGn_alpha                 ! van Genutchen "alpha" parameter (m-1)
 real(dp),intent(in)           :: vGn_n                     ! van Genutchen "n" parameter (-)
 real(dp),intent(in)           :: vGn_m                     ! van Genutchen "m" parameter (-)
 real(dp),intent(in)           :: theta_sat                 ! soil porosity (-)
 real(dp),intent(in)           :: theta_res                 ! soil residual volumetric water content (-)
 real(dp),intent(in)           :: qSurfScale                ! scaling factor in the surface runoff parameterization (-)
 real(dp),intent(in)           :: zScale_TOPMODEL           ! scaling factor used to describe decrease in hydraulic conductivity with depth (m)
 real(dp),intent(in)           :: rootingDepth              ! rooting depth (m)
 real(dp),intent(in)           :: wettingFrontSuction       ! Green-Ampt wetting front suction (m)
 real(dp),intent(in)           :: soilIceScale              ! soil ice scaling factor in Gamma distribution used to define frozen area (m)
 real(dp),intent(in)           :: soilIceCV                 ! soil ice CV in Gamma distribution used to define frozen area (-)
 ! -----------------------------------------------------------------------------------------------------------------------------
 ! input-output: hydraulic conductivity and diffusivity at the surface
 ! NOTE: intent(inout) because infiltration may only be computed for the first iteration
 real(dp),intent(inout)        :: surfaceHydCond            ! hydraulic conductivity (m s-1)
 real(dp),intent(inout)        :: surfaceDiffuse            ! hydraulic diffusivity at the surface (m
 ! output: surface runoff and infiltration flux (m s-1)
 real(dp),intent(inout)        :: xMaxInfilRate             ! maximum infiltration rate (m s-1)
 real(dp),intent(inout)        :: scalarInfilArea           ! fraction of unfrozen area where water can infiltrate (-)
 real(dp),intent(inout)        :: scalarFrozenArea          ! fraction of area that is considered impermeable due to soil ice (-)
 real(dp),intent(out)          :: scalarSurfaceRunoff       ! surface runoff (m s-1)
 real(dp),intent(out)          :: scalarSurfaceInfiltration ! surface infiltration (m s-1)
 ! output: deriavtives in surface infiltration w.r.t. volumetric liquid water (m s-1) and matric head (s-1) in the upper-most soil layer
 real(dp),intent(out)          :: dq_dHydState              ! derivative in surface infiltration w.r.t. state variable in the upper-most soil layer (m s-1 or s-1)
 real(dp),intent(out)          :: dq_dNrgState              ! derivative in surface infiltration w.r.t. energy state variable in the upper-most soil layer (m s-1 K-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                       ! error code
 character(*),intent(out)      :: message                   ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------
 ! local variables
 ! (general)
 integer(i4b)                  :: iLayer                    ! index of soil layer
 ! (head boundary condition)
 real(dp)                      :: cFlux                     ! capillary flux (m s-1)
 real(dp)                      :: vFracLiq                  ! volumetric fraction of liquid water
 real(dp)                      :: surfaceInfiltration1      ! perturbed value of surface infiltration (used to compute the numerical derivative) 
 real(dp)                      :: dNum                      ! numerical derivative 
 ! (simplified Green-Ampt infiltration)
 real(dp)                      :: rootZoneLiq               ! depth of liquid water in the root zone (m)
 real(dp)                      :: rootZoneIce               ! depth of ice in the root zone (m)
 real(dp)                      :: availCapacity             ! available storage capacity in the root zone (m)
 real(dp)                      :: depthWettingFront         ! depth to the wetting front (m)
 real(dp)                      :: hydCondWettingFront       ! hydraulic conductivity at the wetting front (m s-1) 
 ! (saturated area associated with variable storage capacity)
 real(dp)                      :: fracCap                   ! fraction of pore space filled with liquid water and ice (-)
 real(dp)                      :: fInfRaw                   ! infiltrating area before imposing solution constraints (-)
 real(dp),parameter            :: maxFracCap=0.995_dp       ! maximum fraction capacity -- used to avoid numerical problems associated with an enormous derivative
 real(dp),parameter            :: scaleFactor=0.000001_dp   ! scale factor for the smoothing function (-)
 ! (fraction of impermeable area associated with frozen ground)
 real(dp)                      :: alpha                     ! shape parameter in the Gamma distribution
 real(dp)                      :: xLimg                     ! upper limit of the integral
 ! initialize error control
 err=0; message="surfaceFlx/"

 ! compute derivative in the energy state
 ! NOTE: revisit the need to do this
 dq_dNrgState = 0._dp

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
     cflux = -surfaceDiffuse*(scalarVolFracLiq - upperBoundTheta) / (mLayerDepth(1)*0.5_dp)
    case(mixdform)
     ! compute the hydraulic conductivity and diffusivity at the boundary
     surfaceHydCond = hydCond_psi(upperBoundHead,surfaceSatHydCond,vGn_alpha,vGn_n,vGn_m) * iceImpedeFac
     surfaceDiffuse = valueMissing
     ! compute the capillary flux
     cflux = -surfaceHydCond*(scalarMatricHead - upperBoundHead) / (mLayerDepth(1)*0.5_dp)
    case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
   end select  ! (form of Richards' eqn)
   ! compute the total flux
   scalarSurfaceInfiltration = cflux + surfaceHydCond
   ! compute the derivative
   if(deriv_desired)then
    ! compute the hydrology derivative
    select case(ixRichards)  ! (form of Richards' equation)
     case(moisture); dq_dHydState = -surfaceDiffuse/(mLayerDepth(1)/2._dp)
     case(mixdform); dq_dHydState = -surfaceHydCond/(mLayerDepth(1)/2._dp)
     case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
    end select
    ! compute the energy derivative
    dq_dNrgState = -(dHydCond_dTemp/2._dp)*(scalarMatricHead - upperBoundHead)/(mLayerDepth(1)*0.5_dp) + dHydCond_dTemp/2._dp
    ! compute the numerical derivative
    !cflux = -surfaceHydCond*((scalarMatricHead+dx) - upperBoundHead) / (mLayerDepth(1)*0.5_dp)
    !surfaceInfiltration1 = cflux + surfaceHydCond
    !dNum  = (surfaceInfiltration1 - scalarSurfaceInfiltration)/dx
   else
    dq_dHydState = 0._dp
    dNum         = 0._dp
   endif
   !write(*,'(a,1x,10(e30.20,1x))') 'scalarMatricHead, scalarSurfaceInfiltration, dq_dHydState, dNum = ', &
   !                                 scalarMatricHead, scalarSurfaceInfiltration, dq_dHydState, dNum

  ! *****
  ! flux condition
  case(liquidFlux)

   ! force infiltration to be constant over the iterations
   if(doInfiltration)then

    ! define the storage in the root zone (m)
    rootZoneLiq = 0._dp
    rootZoneIce = 0._dp
    ! (process layers where the roots extend to the bottom of the layer)
    if(nRoots > 1)then
     do iLayer=1,nRoots-1
      rootZoneLiq = rootZoneLiq + mLayerVolFracLiq(iLayer)*mLayerDepth(iLayer)
      rootZoneIce = rootZoneIce + mLayerVolFracIce(iLayer)*mLayerDepth(iLayer)
     enddo
    endif
    if(rootingDepth < iLayerHeight(nRoots-1))then; err=20; message=trim(message)//'problem with definition of nRoots'; return; endif
    ! (process layers where the roots end in the current layer)
    rootZoneLiq = rootZoneLiq + mLayerVolFracLiq(nRoots)*(rootingDepth - iLayerHeight(nRoots-1))
    rootZoneIce = rootZoneIce + mLayerVolFracIce(nRoots)*(rootingDepth - iLayerHeight(nRoots-1))

    ! define available capacity to hold water (m)
    availCapacity = theta_sat*rootingDepth - rootZoneIce
    if(rootZoneLiq > availCapacity)then; err=20; message=trim(message)//'liquid water in the root zone exceeds capacity'; return; endif

    ! define the depth to the wetting front (m)
    depthWettingFront = (rootZoneLiq/availCapacity)*rootingDepth

    ! define the hydraulic conductivity at depth=depthWettingFront (m s-1)
    hydCondWettingFront =  surfaceSatHydCond * ( (1._dp - depthWettingFront/sum(mLayerDepth))**(zScale_TOPMODEL - 1._dp) )
 
    ! define the maximum infiltration rate (m s-1)
    xMaxInfilRate = hydCondWettingFront*( (wettingFrontSuction + depthWettingFront)/depthWettingFront )  ! maximum infiltration rate (m s-1)
    !write(*,'(a,1x,f9.3,1x,10(e20.10,1x))') 'depthWettingFront, surfaceSatHydCond, hydCondWettingFront, xMaxInfilRate = ', depthWettingFront, surfaceSatHydCond, hydCondWettingFront, xMaxInfilRate

    ! define the infiltrating area for the non-frozen part of the cell/basin
    fracCap         = rootZoneLiq/(maxFracCap*availCapacity)                              ! fraction of available root zone filled with water
    fInfRaw         = 1._dp - exp(-qSurfScale*(1._dp - fracCap))                          ! infiltrating area -- allowed to violate solution constraints
    scalarInfilArea = min(0.5_dp*(fInfRaw + sqrt(fInfRaw**2._dp + scaleFactor)), 1._dp)   ! infiltrating area -- constrained
    !print*, 'scalarInfilArea = ', scalarInfilArea

    ! check to ensure we are not infiltrating into a fully saturated column
    if(sum(mLayerVolFracLiq(ixIce+1:nRoots)*mLayerDepth(ixIce+1:nRoots)) > 0.9999_dp*theta_sat*sum(mLayerDepth(ixIce+1:nRoots))) scalarInfilArea=0._dp
    !print*, 'ixIce, nRoots, scalarInfilArea = ', ixIce, nRoots, scalarInfilArea
    !print*, 'sum(mLayerVolFracLiq(ixIce+1:nRoots)*mLayerDepth(ixIce+1:nRoots)) = ', sum(mLayerVolFracLiq(ixIce+1:nRoots)*mLayerDepth(ixIce+1:nRoots))
    !print*, 'theta_sat*sum(mLayerDepth(ixIce+1:nRoots)) = ', theta_sat*sum(mLayerDepth(ixIce+1:nRoots))

    ! define the impermeable area due to frozen ground
    if(rootZoneIce > tiny(rootZoneIce))then  ! (avoid divide by zero)
     alpha            = 1._dp/(soilIceCV**2._dp)        ! shape parameter in the Gamma distribution
     xLimg            = alpha*soilIceScale/rootZoneIce  ! upper limit of the integral
     scalarFrozenArea = 1._dp - gammp(alpha,xLimg)      ! fraction of frozen area
    else
     scalarFrozenArea = 0._dp
    endif
    !print*, 'scalarFrozenArea, rootZoneIce = ', scalarFrozenArea, rootZoneIce

   endif ! (if desire to compute infiltration)

   ! compute infiltration (m s-1)
   scalarSurfaceInfiltration = (1._dp - scalarFrozenArea)*scalarInfilArea*min(scalarRainPlusMelt,xMaxInfilRate)
   !print*, 'scalarSurfaceInfiltration = ', scalarSurfaceInfiltration
   !print*, '(1._dp - scalarFrozenArea), (1._dp - scalarFrozenArea)*scalarInfilArea = ', (1._dp - scalarFrozenArea), (1._dp - scalarFrozenArea)*scalarInfilArea

   ! compute surface runoff (m s-1)
   scalarSurfaceRunoff = scalarRainPlusMelt - scalarSurfaceInfiltration

   ! set surface hydraulic conductivity and diffusivity to missing (not used for flux condition)
   surfaceHydCond = valueMissing
   surfaceDiffuse = valueMissing

   ! set numerical derivative to zero
   ! NOTE 1: Depends on multiple soil layers and does not jive with the current tridiagonal matrix
   ! NOTE 2: Need to define the derivative at every call, because intent(out)
   dq_dHydState = 0._dp
   dq_dNrgState = 0._dp

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
                       ! input: temperature derivatives
                       dPsiLiq_dTemp,             & ! intent(in): derivative in liquid water matric potential w.r.t. temperature (m K-1)
                       dHydCond_dTemp,            & ! intent(in): derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
                       ! input: transmittance (adjacent layers)
                       nodeHydCondTrial,          & ! intent(in): hydraulic conductivity at the soil nodes (m s-1)
                       nodeDiffuseTrial,          & ! intent(in): hydraulic diffusivity at the soil nodes (m2 s-1)
                       ! input: transmittance derivatives (adjacent layers)
                       dHydCond_dVolLiq,          & ! intent(in): derivative in hydraulic conductivity w.r.t. change in volumetric liquid water content (m s-1)
                       dDiffuse_dVolLiq,          & ! intent(in): derivative in hydraulic diffusivity w.r.t. change in volumetric liquid water content (m2 s-1)
                       dHydCond_dMatric,          & ! intent(in): derivative in hydraulic conductivity w.r.t. change in matric head (s-1)
                       ! output: tranmsmittance at the layer interface (scalars)
                       iLayerHydCond,             & ! intent(out): hydraulic conductivity at the interface between layers (m s-1)
                       iLayerDiffuse,             & ! intent(out): hydraulic diffusivity at the interface between layers (m2 s-1)
                       ! output: vertical flux at the layer interface (scalars)
                       iLayerLiqFluxSoil,         & ! intent(out): vertical flux of liquid water at the layer interface (m s-1)
                       ! output: derivatives in fluxes w.r.t. state variables -- matric head or volumetric lquid water -- in the layer above and layer below (m s-1 or s-1)
                       dq_dHydStateAbove,         & ! intent(out): derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer above (m s-1 or s-1)
                       dq_dHydStateBelow,         & ! intent(out): derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer below (m s-1 or s-1) 
                       ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
                       dq_dNrgStateAbove,         & ! intent(out): derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
                       dq_dNrgStateBelow,         & ! intent(out): derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
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
 ! input: temperature derivatives
 real(dp),intent(in)           :: dPsiLiq_dTemp(:)            ! derivative in liquid water matric potential w.r.t. temperature (m K-1)
 real(dp),intent(in)           :: dHydCond_dTemp(:)           ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
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
 real(dp),intent(out)          :: dq_dHydStateAbove           ! derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer above (m s-1 or s-1)
 real(dp),intent(out)          :: dq_dHydStateBelow           ! derivatives in the flux w.r.t. matric head or volumetric lquid water in the layer below (m s-1 or s-1) 
 ! output: derivatives in fluxes w.r.t. energy state variables -- now just temperature -- in the layer above and layer below (m s-1 K-1)
 real(dp),intent(out)          :: dq_dNrgStateAbove           ! derivatives in the flux w.r.t. temperature in the layer above (m s-1 K-1)
 real(dp),intent(out)          :: dq_dNrgStateBelow           ! derivatives in the flux w.r.t. temperature in the layer below (m s-1 K-1)
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
 !write(*,'(a,1x,10(e20.10,1x))') 'iLayerLiqFluxSoil, dPsi, dz, cflux, iLayerHydCond = ', &
 !                                 iLayerLiqFluxSoil, dPsi, dz, cflux, iLayerHydCond

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
    dq_dHydStateAbove = -dDiffuseIface_dVolLiqAbove*dLiq/dz + iLayerDiffuse/dz + dHydCondIface_dVolLiqAbove
    dq_dHydStateBelow = -dDiffuseIface_dVolLiqBelow*dLiq/dz - iLayerDiffuse/dz + dHydCondIface_dVolLiqBelow
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
    dq_dHydStateAbove = -dHydCondIface_dMatricAbove*dPsi/dz + iLayerHydCond/dz + dHydCondIface_dMatricAbove
    dq_dHydStateBelow = -dHydCondIface_dMatricBelow*dPsi/dz - iLayerHydCond/dz + dHydCondIface_dMatricBelow
    ! derivative in the flux w.r.t. temperature
    dq_dNrgStateAbove = -(dHydCond_dTemp(ixUpper)/2._dp)*dPsi/dz + iLayerHydCond*dPsiLiq_dTemp(ixUpper)/dz + dHydCond_dTemp(ixUpper)/2._dp
    dq_dNrgStateBelow = -(dHydCond_dTemp(ixLower)/2._dp)*dPsi/dz - iLayerHydCond*dPsiLiq_dTemp(ixLower)/dz + dHydCond_dTemp(ixLower)/2._dp
   case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
  end select
 else
  dq_dHydStateAbove = valueMissing
  dq_dHydStateBelow = valueMissing
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
                       iceImpedeFac,              & ! intent(in): ice impedence factor in the lower-most soil layer (-)
                       ! input: transmittance derivatives
                       dHydCond_dVolLiq,          & ! intent(in): derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
                       dHydCond_dMatric,          & ! intent(in): derivative in hydraulic conductivity w.r.t. matric head (s-1)
                       dHydCond_dTemp,            & ! intent(in): derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
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
                       dq_dHydStateUnsat,         & ! intent(out): change in drainage flux w.r.t. change in hydrology state variable in lowest unsaturated node (m s-1 or s-1)
                       dq_dNrgStateUnsat,         & ! intent(out): change in drainage flux w.r.t. change in energy state variable in lowest unsaturated node (m s-1 K-1)
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
 real(dp),intent(in)           :: iceImpedeFac              ! ice impedence factor in the upper-most soil layer (-)
 ! input: transmittance derivatives
 real(dp),intent(in)           :: dHydCond_dVolLiq          ! derivative in hydraulic conductivity w.r.t. volumetric liquid water content (m s-1)
 real(dp),intent(in)           :: dHydCond_dMatric          ! derivative in hydraulic conductivity w.r.t. matric head (s-1)
 real(dp),intent(in)           :: dHydCond_dTemp            ! derivative in hydraulic conductivity w.r.t temperature (m s-1 K-1)
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
 real(dp),intent(out)          :: dq_dHydStateUnsat         ! change in drainage flux w.r.t. change in state variable in lowest unsaturated node (m s-1 or s-1)
 real(dp),intent(out)          :: dq_dNrgStateUnsat         ! change in drainage flux w.r.t. change in energy state variable in lowest unsaturated node (m s-1 K-1)
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
    ! hydrology derivatives
    select case(ixRichards)  ! (form of Richards' equation)
     case(moisture); dq_dHydStateUnsat = bottomDiffuse/(nodeDepth/2._dp)
     case(mixdform); dq_dHydStateUnsat = bottomHydCond/(nodeDepth/2._dp)
     case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
    end select
    ! energy derivatives
    err=20; message=trim(message)//"not yet implemented energy derivatives"; return
   else     ! (do not desire derivatives)
    dq_dHydStateUnsat = valueMissing
    dq_dNrgStateUnsat = valueMissing
   endif

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
    ! hydrology derivatives
    select case(ixRichards)  ! (form of Richards' equation)
     case(moisture); dq_dHydStateUnsat = kAnisotropic*surfaceSatHydCond * node__dPsi_dTheta*exp(-zWater/zScale_TOPMODEL)/zScale_TOPMODEL
     case(mixdform); dq_dHydStateUnsat = kAnisotropic*surfaceSatHydCond * exp(-zWater/zScale_TOPMODEL)/zScale_TOPMODEL
     case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
    end select
    ! energy derivatives
    err=20; message=trim(message)//"not yet implemented energy derivatives"; return
   else     ! (do not desire derivatives)
    dq_dHydStateUnsat = valueMissing
    dq_dNrgStateUnsat = valueMissing
   endif

  ! ---------------------------------------------------------------------------------------------
  ! * free drainage
  ! ---------------------------------------------------------------------------------------------
  case(freeDrainage)

   ! compute flux
   scalarDrainage = nodeHydCond*kAnisotropic

   ! compute derivatives
   if(deriv_desired)then
    ! hydrology derivatives
    select case(ixRichards)  ! (form of Richards' equation)
     case(moisture); dq_dHydStateUnsat = dHydCond_dVolLiq*kAnisotropic
     case(mixdform); dq_dHydStateUnsat = dHydCond_dMatric*kAnisotropic
     case default; err=10; message=trim(message)//"unknown form of Richards' equation"; return
    end select
    ! energy derivatives
    dq_dNrgStateUnsat = dHydCond_dTemp*kAnisotropic
   else     ! (do not desire derivatives)
    dq_dHydStateUnsat = valueMissing
    dq_dNrgStateUnsat = valueMissing
   endif


  ! ---------------------------------------------------------------------------------------------
  ! * zero flux
  ! ---------------------------------------------------------------------------------------------
  case(zeroFlux)
   scalarDrainage = 0._dp
   if(deriv_desired)then
    dq_dHydStateUnsat = 0._dp
    dq_dNrgStateUnsat = 0._dp
   else
    dq_dHydStateUnsat = valueMissing
    dq_dNrgStateUnsat = valueMissing
   endif

  ! ---------------------------------------------------------------------------------------------
  ! * error check
  ! ---------------------------------------------------------------------------------------------
  case default; err=20; message=trim(message)//'unknown lower boundary condition for soil hydrology'; return

 endselect ! (type of boundary condition)

 end subroutine qDrainFlux


 ! *******************************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************************

end module soilLiqFlx_module

module vegNrgFlux_module
USE nrtype
! constants
USE multiconst,only:gravity    ! acceleration of gravity              (m s-2)
USE multiconst,only:vkc        ! von Karman's constant                (-)
USE multiconst,only:w_ratio    ! molecular ratio water to dry air     (-)
USE multiconst,only:R_wv       ! gas constant for water vapor         (Pa K-1 m3 kg-1; J kg-1 K-1)
USE multiconst,only:Cp_air     ! specific heat of air                 (J kg-1 K-1)
USE multiconst,only:Cp_ice     ! specific heat of ice                 (J kg-1 K-1)
USE multiconst,only:Cp_soil    ! specific heat of soil                (J kg-1 K-1)
USE multiconst,only:Cp_water   ! specific heat of liquid water        (J kg-1 K-1)
USE multiconst,only:Tfreeze    ! temperature at freezing              (K)
USE multiconst,only:LH_fus     ! latent heat of fusion                (J kg-1)
USE multiconst,only:LH_vap     ! latent heat of vaporization          (J kg-1)
USE multiconst,only:LH_sub     ! latent heat of sublimation           (J kg-1)
USE multiconst,only:sb         ! Stefan Boltzman constant             (W m-2 K-4)
USE multiconst,only:iden_air   ! intrinsic density of air             (kg m-3)
USE multiconst,only:iden_ice   ! intrinsic density of ice             (kg m-3)
USE multiconst,only:iden_water ! intrinsic density of liquid water    (kg m-3)
! look-up values for method used to compute derivative
USE mDecisions_module,only:  &
 numerical,                  & ! numerical solution
 analytical                    ! analytical solution
! look-up values for choice of stability function
USE mDecisions_module,only:  &
 standard,                   & ! standard MO similarity, a la Anderson (1976) 
 louisInversePower,          & ! Louis (1979) inverse power function
 mahrtExponential              ! Mahrt (1987) exponential
! named variables for snow and soil
USE data_struc,only:ix_soil,ix_snow            ! named variables for snow and soil
! -------------------------------------------------------------------------------------------------
implicit none
private
public::vegNrgFlux
! dimensions
integer(i4b)                  :: nSoil         ! number of soil layers
integer(i4b)                  :: nSnow         ! number of snow layers
integer(i4b)                  :: nLayers       ! total number of layers
integer(i4b),parameter        :: nBands=2      ! number of spectral bands for shortwave radiation
! named variables
integer(i4b),parameter        :: ist     = 1   ! Surface type:  IST=1 => soil;  IST=2 => lake
integer(i4b),parameter        :: isc     = 4   ! Soil color type
integer(i4b),parameter        :: ice     = 0   ! Surface type:  ICE=0 => soil;  ICE=1 => sea-ice
! spatial indices
integer(i4b),parameter        :: iLoc    = 1   ! i-location
integer(i4b),parameter        :: jLoc    = 1   ! j-location
! algorithmic parameters
real(dp),parameter     :: verySmall=epsilon(1._dp) ! a very small number
real(dp),parameter     :: mpe=1.e-6_dp         ! prevents overflow error if division by zero 
real(dp),parameter     :: dx=1.e-6_dp          ! finite difference increment
contains

 ! ************************************************************************************************
 ! new subroutine: muster program to compute energy fluxes at vegetation and ground surfaces
 ! ************************************************************************************************
 subroutine vegNrgFlux(&
                       ! input
                       dt,                                      & ! intent(in): time step (seconds)
                       iter,                                    & ! intent(in): iteration index
                       firstSubStep,                            & ! intent(in): flag to indicate if we are processing the first sub-step
                       canopyTempTrial,                         & ! intent(in): trial value of canopy temperature (K)
                       groundTempTrial,                         & ! intent(in): trial value of ground temperature (K)
                       canopyLiqTrial,                          & ! intent(in): trial value of mass of liquid water on the vegetation canopy (kg m-2)
                       canopyIceTrial,                          & ! intent(in): trial value of mass of ice on the vegetation canopy (kg m-2)
                       ! output
                       canopyNetFlux,                           & ! intent(out): net energy flux for the vegetation canopy (W m-2)
                       groundNetFlux,                           & ! intent(out): net energy flux for the ground surface (W m-2)
                       dCanopyNetFlux_dCanopyTemp,              & ! intent(out): derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
                       dCanopyNetFlux_dGroundTemp,              & ! intent(out): derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
                       dGroundNetFlux_dCanopyTemp,              & ! intent(out): derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
                       dGroundNetFlux_dGroundTemp,              & ! intent(out): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
                       err,message)                               ! intent(out): error control
 ! model decisions
 USE data_struc,only:model_decisions                              ! model decision structure
 USE var_lookup,only:iLookDECISIONS                               ! named variables for elements of the decision structure
 ! common variables
 USE data_struc,only:urbanVegCategory                             ! vegetation category for urban areas
 USE data_struc,only:fracJulday                                   ! fractional julian days since the start of year
 USE data_struc,only:yearLength                                   ! number of days in the current year
 ! model variables, parameters, etc.
 USE data_struc,only:time_data,type_data,attr_data,forc_data,mpar_data,mvar_data,indx_data     ! data structures
 USE var_lookup,only:iLookTIME,iLookTYPE,iLookATTR,iLookFORCE,iLookPARAM,iLookMVAR,iLookINDEX  ! named variables for structure elements
 ! compute energy and mass fluxes for vegetation
 implicit none
 ! input
 real(dp),intent(in)           :: dt                              ! time step (seconds)
 integer(i4b),intent(in)       :: iter                            ! iteration index
 logical(i4b),intent(in)       :: firstSubStep                    ! flag to indicate if we are processing the first sub-step
 real(dp),intent(in)           :: canopyTempTrial                 ! trial value of canopy temperature (K)
 real(dp),intent(in)           :: groundTempTrial                 ! trial value of ground temperature (K)
 real(dp),intent(in)           :: canopyLiqTrial                  ! trial value of mass of liquid water on the vegetation canopy (kg m-2)
 real(dp),intent(in)           :: canopyIceTrial                  ! trial value of mass of ice on the vegetation canopy (kg m-2)
 ! output
 real(dp),intent(out)          :: canopyNetFlux                   ! net energy flux for the vegetation canopy (W m-2)
 real(dp),intent(out)          :: groundNetFlux                   ! net energy flux for the ground surface (W m-2) 
 real(dp),intent(out)          :: dCanopyNetFlux_dCanopyTemp      ! derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp),intent(out)          :: dCanopyNetFlux_dGroundTemp      ! derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
 real(dp),intent(out)          :: dGroundNetFlux_dCanopyTemp      ! derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp),intent(out)          :: dGroundNetFlux_dGroundTemp      ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
 integer(i4b),intent(out)      :: err                             ! error code
 character(*),intent(out)      :: message                         ! error message
 ! local
 character(LEN=256)            :: cmessage                        ! error message of downwind routine
 integer(i4b)                  :: nLayersRoots                    ! number of soil layers that contain roots
 ! initialize error control
 err=0; message="vegNrgFlux_muster/"

 ! define number of layers
 nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat == ix_soil)  ! number of soil layers
 nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat == ix_snow)  ! number of snow layers
 nLayers = indx_data%var(iLookINDEX%nLayers)%dat(1)                   ! total number of layers

 ! preliminaries
 if(iter==1 .and. firstSubStep)then
  ! compute the number of layers with roots
  nLayersRoots = count(mvar_data%var(iLookMVAR%iLayerHeight)%dat(nSnow:nLayers-1) < mpar_data%var(iLookPARAM%rootingDepth)-verySmall)
  if(nLayersRoots == 0)then; err=20; message=trim(message)//'no roots within the soil profile'; return; endif
  ! compute the temperature of the root zone (K)
  mvar_data%var(iLookMVAR%scalarRootZoneTemp)%dat(1) = sum(mvar_data%var(iLookMVAR%mLayerTemp)%dat(nSnow+1:nSnow+nLayersRoots)) / real(nLayersRoots, kind(dp))
 endif

 ! subroutine to compute energy fluxes at vegetation and ground surfaces
 ! NOTE: separate call to use data structure elements, and ensure variables are used as they are intended (input, input/output, output)
 call vegNrgFlux_muster(&

                        ! input
                        dt,                                                            & ! intent(in): time step (seconds)
                        iter,                                                          & ! intent(in): iteration index
                        firstSubStep,                                                  & ! intent(in): flag to indicate if we are processing the first sub-step
                        canopyTempTrial,                                               & ! intent(in): trial value of canopy temperature (K)
                        groundTempTrial,                                               & ! intent(in): trial value of ground temperature (K)
                        canopyLiqTrial,                                                & ! intent(in): trial value of mass of liquid water on the vegetation canopy (kg m-2)
                        canopyIceTrial,                                                & ! intent(in): trial value of mass of ice on the vegetation canopy (kg m-2)

                        ! model control -- intent(in)
                        model_decisions(iLookDECISIONS%fDerivMeth)%iDecision,          & ! intent(in): choice of method to compute derivatives
                        model_decisions(iLookDECISIONS%astability)%iDecision,          & ! intent(in): choice of stability function
                        model_decisions(iLookDECISIONS%soilStress)%iDecision,          & ! intent(in): choice of function for the soil moisture control on stomatal resistance
                        model_decisions(iLookDECISIONS%groundwatr)%iDecision,          & ! intent(in): groundwater parameterization
                        model_decisions(iLookDECISIONS%stomResist)%iDecision,          & ! intent(in): choice of function for stomatal resistance

                        ! time -- intent(in)
                        yearLength,                                                    & ! intent(in): number of days in the current year
                        fracJulday,                                                    & ! intent(in): fractional julian days since the start of year

                        ! physical attributes -- intent(in)
                        attr_data%var(iLookATTR%mHeight),                              & ! intent(in): measurement height (m)
                        attr_data%var(iLookATTR%latitude),                             & ! intent(in): latitude (degrees north)
                        type_data%var(iLookTYPE%vegTypeIndex),                         & ! intent(in): vegetation type index
                        urbanVegCategory,                                              & ! intent(in): vegetation category for urban areas

                        ! model parameters (phenology) -- intent(in)
                        mpar_data%var(iLookPARAM%canopyHeight),                        & ! intent(in): height of the vegetation canopy (m)
                        mpar_data%var(iLookPARAM%maxCanopyLiquid),                     & ! intent(in): maximum storage of liquid water on the vegetation canopy (kg m-2)
                        mpar_data%var(iLookPARAM%maxCanopyIce),                        & ! intent(in): maximum storage of ice on the vegetation canopy (kg m-2)

                        ! model parameters (aerodynamic resistance) -- intent(in)
                        mpar_data%var(iLookPARAM%z0Snow),                              & ! intent(in): roughness length of snow (m)
                        mpar_data%var(iLookPARAM%z0Soil),                              & ! intent(in): roughness length of soil (m)
                        mpar_data%var(iLookPARAM%z0Canopy),                            & ! intent(in): roughness length of the canopy (m)
                        mpar_data%var(iLookPARAM%critRichNumber),                      & ! intent(in): critical value for the bulk Richardson number where turbulence ceases (-)
                        mpar_data%var(iLookPARAM%Louis79_bparam),                      & ! intent(in): parameter in Louis (1979) stability function
                        mpar_data%var(iLookPARAM%Louis79_cStar),                       & ! intent(in): parameter in Louis (1979) stability function
                        mpar_data%var(iLookPARAM%Mahrt87_eScale),                      & ! intent(in): exponential scaling factor in the Mahrt (1987) stability function
                        mpar_data%var(iLookPARAM%windReductionParam),                  & ! intent(in): canopy wind reduction parameter (-)                   
                        mpar_data%var(iLookPARAM%leafExchangeCoeff),                   & ! intent(in): turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
                        mpar_data%var(iLookPARAM%leafDimension),                       & ! intent(in): characteristic leaf dimension (m)

                        ! model parameters (soil stress) -- intent(in)
                        mpar_data%var(iLookPARAM%theta_sat),                           & ! intent(in): soil porosity (-)
                        mpar_data%var(iLookPARAM%plantWiltPsi),                        & ! intent(in): matric head at wilting point (m)
                        mpar_data%var(iLookPARAM%soilStressParam),                     & ! intent(in): parameter in the exponential soil stress function (-)
                        mpar_data%var(iLookPARAM%critSoilWilting),                     & ! intent(in): critical vol. liq. water content when plants are wilting (-)
                        mpar_data%var(iLookPARAM%critSoilTranspire),                   & ! intent(in): critical vol. liq. water content when transpiration is limited (-)
                        mpar_data%var(iLookPARAM%critAquiferTranspire),                & ! intent(in): critical aquifer storage value when transpiration is limited (m)
       
                        ! forcing at the upper boundary -- intent(in)
                        forc_data%var(iLookFORCE%airtemp),                             & ! intent(in): air temperature at some height above the surface (K)
                        forc_data%var(iLookFORCE%windspd),                             & ! intent(in): wind speed at some height above the surface (m s-1)
                        forc_data%var(iLookFORCE%airpres),                             & ! intent(in): air pressure at some height above the surface (Pa)
                        forc_data%var(iLookFORCE%LWRadAtm),                            & ! intent(in): downwelling longwave radiation at the upper boundary (W m-2)
                        mvar_data%var(iLookMVAR%scalarVPair)%dat(1),                   & ! intent(in): vapor pressure at some height above the surface (Pa)
                        mvar_data%var(iLookMVAR%scalarO2air)%dat(1),                   & ! intent(in): atmospheric o2 concentration (Pa)
                        mvar_data%var(iLookMVAR%scalarCO2air)%dat(1),                  & ! intent(in): atmospheric co2 concentration (Pa)
                        mvar_data%var(iLookMVAR%scalarSnowfall)%dat(1),                & ! intent(in): snowfall rate (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarCosZenith)%dat(1),               & ! intent(in): cosine of the solar zenith angle (0-1)
                        mvar_data%var(iLookMVAR%spectralIncomingDirect)%dat(1:nBands), & ! intent(in): incoming direct solar radiation in each wave band (w m-2)
                        mvar_data%var(iLookMVAR%spectralIncomingDiffuse)%dat(1:nBands),& ! intent(in): incoming diffuse solar radiation in each wave band (w m-2)

                        ! water storage -- intent(in)
                        ! NOTE: soil stress only computed at the start of the substep (iter==1)
                        mvar_data%var(iLookMVAR%scalarSWE)%dat(1),                     & ! intent(in): snow water equivalent on the ground (kg m-2)
                        mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1),               & ! intent(in): snow depth on the ground surface (m)
                        mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(nSnow+1:nLayers),& ! intent(in): volumetric fraction of liquid water in each soil layer (-)
                        mvar_data%var(iLookMVAR%mLayerMatricHead)%dat,                 & ! intent(in): matric head in each layer (m)
                        mvar_data%var(iLookMVAR%scalarAquiferStorage)%dat(1),          & ! intent(in): aquifer storage (m)

                        ! vegetation phenology -- intent(inout) because only called at the start of the sub-step
                        mvar_data%var(iLookMVAR%scalarCanopyWetFraction)%dat(1),       & ! intent(inout): fraction of canopy that is wet
                        mvar_data%var(iLookMVAR%scalarFoliageNitrogenFactor)%dat(1),   & ! intent(inout): foliage nitrogen concentration (1.0 = saturated)
                        mvar_data%var(iLookMVAR%scalarRootZoneTemp)%dat(1),            & ! intent(inout): root zone temperature
                        mvar_data%var(iLookMVAR%scalarLAI)%dat(1),                     & ! intent(inout): one-sided leaf area index (m2 m-2)
                        mvar_data%var(iLookMVAR%scalarSAI)%dat(1),                     & ! intent(inout): one-sided stem area index (m2 m-2)
                        mvar_data%var(iLookMVAR%scalarExposedLAI)%dat(1),              & ! intent(inout): exposed leaf area index after burial by snow (m2 m-2)
                        mvar_data%var(iLookMVAR%scalarExposedSAI)%dat(1),              & ! intent(inout): exposed stem area index after burial by snow (m2 m-2)
                        mvar_data%var(iLookMVAR%scalarGrowingSeasonIndex)%dat(1),      & ! intent(inout): growing season index (0=off, 1=on)

                        ! shortwave radiation fluxes -- intent(inout) because only called at the start of the sub-step
                        mvar_data%var(iLookMVAR%scalarAlbedo)%dat(1),                  & ! intent(inout): snow albedo (-)
                        mvar_data%var(iLookMVAR%scalarSnowAge)%dat(1),                 & ! intent(inout): non-dimensional snow age (-)
                        mvar_data%var(iLookMVAR%scalarCanopySunlitFraction)%dat(1),    & ! intent(inout): sunlit fraction of canopy (-)
                        mvar_data%var(iLookMVAR%scalarCanopySunlitLAI)%dat(1),         & ! intent(inout): sunlit leaf area (-)
                        mvar_data%var(iLookMVAR%scalarCanopyShadedLAI)%dat(1),         & ! intent(inout): shaded leaf area (-)
                        mvar_data%var(iLookMVAR%scalarCanopySunlitPAR)%dat(1),         & ! intent(inout): average absorbed par for sunlit leaves (w m-2)
                        mvar_data%var(iLookMVAR%scalarCanopyShadedPAR)%dat(1),         & ! intent(inout): average absorbed par for shaded leaves (w m-2)
                        mvar_data%var(iLookMVAR%scalarCanopyAbsorbedSolar)%dat(1),     & ! intent(inout): solar radiation absorbed by canopy (W m-2)
                        mvar_data%var(iLookMVAR%scalarGroundAbsorbedSolar)%dat(1),     & ! intent(inout): solar radiation absorbed by ground (W m-2)
                        mvar_data%var(iLookMVAR%scalarTotalReflectedSolar)%dat(1),     & ! intent(inout): total reflected solar radiation (W m-2)
                        mvar_data%var(iLookMVAR%scalarTotalAbsorbedSolar)%dat(1),      & ! intent(inout): total absorbed solar radiation (W m-2)
                        mvar_data%var(iLookMVAR%scalarCanopyReflectedSolar)%dat(1),    & ! intent(inout): solar radiation reflected from the canopy (W m-2)
                        mvar_data%var(iLookMVAR%scalarGroundReflectedSolar)%dat(1),    & ! intent(inout): solar radiation reflected from the ground (W m-2) 
                        mvar_data%var(iLookMVAR%scalarBetweenCanopyGapFraction)%dat(1),& ! intent(inout): between canopy gap fraction for beam (-)
                        mvar_data%var(iLookMVAR%scalarWithinCanopyGapFraction)%dat(1), & ! intent(inout): within canopy gap fraction for beam (-)

                        ! longwave radiation fluxes -- intent(out) because called in every step
                        mvar_data%var(iLookMVAR%scalarLWRadCanopy)%dat(1),             & ! intent(out): longwave radiation emitted from the canopy (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWRadGround)%dat(1),             & ! intent(out): longwave radiation emitted at the ground surface (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWRadUbound2Canopy)%dat(1),      & ! intent(out): downward atmospheric longwave radiation absorbed by the canopy (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWRadUbound2Ground)%dat(1),      & ! intent(out): downward atmospheric longwave radiation absorbed by the ground (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWRadUbound2Ubound)%dat(1),      & ! intent(out): atmospheric radiation reflected by the ground and lost thru upper boundary (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWRadCanopy2Ubound)%dat(1),      & ! intent(out): longwave radiation emitted from canopy lost thru upper boundary (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWRadCanopy2Ground)%dat(1),      & ! intent(out): longwave radiation emitted from canopy absorbed by the ground (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWRadCanopy2Canopy)%dat(1),      & ! intent(out): canopy longwave reflected from ground and absorbed by the canopy (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWRadGround2Ubound)%dat(1),      & ! intent(out): longwave radiation emitted from ground lost thru upper boundary (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWRadGround2Canopy)%dat(1),      & ! intent(out): longwave radiation emitted from ground and absorbed by the canopy (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWNetCanopy)%dat(1),             & ! intent(out): net longwave radiation at the canopy (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWNetGround)%dat(1),             & ! intent(out): net longwave radiation at the ground surface (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWNetUbound)%dat(1),             & ! intent(out): net longwave radiation at the upper boundary (W m-2)

                        ! aerodynamic resistance -- intent(out) because called in every step
                        mvar_data%var(iLookMVAR%scalarWindReductionFactor)%dat(1),     & ! intent(out): canopy wind reduction factor (-)
                        mvar_data%var(iLookMVAR%scalarZeroPlaneDisplacement)%dat(1),   & ! intent(out): zero plane displacement (m) 
                        mvar_data%var(iLookMVAR%scalarSfc2AtmExchangeCoeff)%dat(1),    & ! intent(out): surface-atmosphere turbulent exchange coefficient (-)
                        mvar_data%var(iLookMVAR%scalarEddyDiffusCanopyTop)%dat(1),     & ! intent(out): eddy diffusivity for heat at the top of the canopy (m2 s-1)
                        mvar_data%var(iLookMVAR%scalarWindspdCanopyTop)%dat(1),        & ! intent(out): windspeed at the top of the canopy (m s-1)
                        mvar_data%var(iLookMVAR%scalarLeafResistance)%dat(1),          & ! intent(out): mean leaf boundary layer resistance per unit leaf area (s m-1)
                        mvar_data%var(iLookMVAR%scalarGroundResistance)%dat(1),        & ! intent(out): below canopy aerodynamic resistance (s m-1) 
                        mvar_data%var(iLookMVAR%scalarCanopyResistance)%dat(1),        & ! intent(out): above canopy aerodynamic resistance (s m-1)

                        ! soil resistance -- intent(in) and intent(inout) because only called at the first iteration
                        mvar_data%var(iLookMVAR%mLayerRootDensity)%dat,                & ! intent(in): root density in each layer (-)
                        mvar_data%var(iLookMVAR%scalarAquiferRootFrac)%dat(1),         & ! intent(in): fraction of roots below the lowest soil layer (-)
                        mvar_data%var(iLookMVAR%scalarTranspireLim)%dat(1),            & ! intent(inout): weighted average of the transpiration limiting factor (-)
                        mvar_data%var(iLookMVAR%mLayerTranspireLim)%dat,               & ! intent(inout): transpiration limiting factor in each layer (-)
                        mvar_data%var(iLookMVAR%scalarTranspireLimAqfr)%dat(1),        & ! intent(inout): transpiration limiting factor for the aquifer (-)

                        ! stomatal resistance -- intent(inout) because only called at the first iteration
                        mvar_data%var(iLookMVAR%scalarStomResistSunlit)%dat(1),        & ! intent(inout): stomatal resistance for sunlit leaves (s m-1)
                        mvar_data%var(iLookMVAR%scalarStomResistShaded)%dat(1),        & ! intent(inout): stomatal resistance for shaded leaves (s m-1)
                        mvar_data%var(iLookMVAR%scalarPhotosynthesisSunlit)%dat(1),    & ! intent(inout): sunlit photosynthesis (umolco2 m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarPhotosynthesisShaded)%dat(1),    & ! intent(inout): shaded photosynthesis (umolco2 m-2 s-1)

                        ! turbulent heat fluxes
                        mvar_data%var(iLookMVAR%scalarPsycConstCanopy)%dat(1),         & ! intent(inout): psychrometric constant for the vegetation canopy (Pa K-1)
                        mvar_data%var(iLookMVAR%scalarPsycConstGround)%dat(1),         & ! intent(inout): psychrometric constant for the ground surface (Pa K-1)
                        mvar_data%var(iLookMVAR%scalarTemp_CanopyAir)%dat(1),          & ! intent(inout): temperature of the canopy air space (K)
                        mvar_data%var(iLookMVAR%scalarVP_CanopyAir)%dat(1),            & ! intent(inout): vapor pressure of the canopy air space (Pa)
                        mvar_data%var(iLookMVAR%scalarSatVP_CanopyTemp)%dat(1),        & ! intent(out): saturation vapor pressure at the temperature of the vegetation canopy (Pa)
                        mvar_data%var(iLookMVAR%scalarSatVP_GroundTemp)%dat(1),        & ! intent(out): saturation vapor pressure at the temperature of the ground surface (Pa)
                        mvar_data%var(iLookMVAR%scalarSoilRelHumidity)%dat(1),         & ! intent(out): relative humidity in the soil pores [0-1]
                        mvar_data%var(iLookMVAR%scalarSoilResistance)%dat(1),          & ! intent(out): resistance from the soil (s m-1)
                        mvar_data%var(iLookMVAR%scalarSenHeatCanopy)%dat(1),           & ! intent(out): sensible heat flux from the canopy to the canopy air space (W m-2)
                        mvar_data%var(iLookMVAR%scalarLatHeatCanopyEvap)%dat(1),       & ! intent(out): latent heat flux associated with evaporation from the canopy to the canopy air space (W m-2)
                        mvar_data%var(iLookMVAR%scalarLatHeatCanopyTrans)%dat(1),      & ! intent(out): latent heat flux associated with transpiration from the canopy to the canopy air space (W m-2)
                        mvar_data%var(iLookMVAR%scalarSenHeatGround)%dat(1),           & ! intent(out): sensible heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
                        mvar_data%var(iLookMVAR%scalarLatHeatGround)%dat(1),           & ! intent(out): latent heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)

                        ! output
                        canopyNetFlux,                                                 & ! intent(out): net energy flux for the vegetation canopy (W m-2)
                        groundNetFlux,                                                 & ! intent(out): net energy flux for the ground surface (W m-2)
                        dCanopyNetFlux_dCanopyTemp,                                    & ! intent(out): derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
                        dCanopyNetFlux_dGroundTemp,                                    & ! intent(out): derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
                        dGroundNetFlux_dCanopyTemp,                                    & ! intent(out): derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
                        dGroundNetFlux_dGroundTemp,                                    & ! intent(out): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
                        err,cmessage)                                                    ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


 end subroutine vegNrgFlux





 ! **********************************************************************************************************************************************************************************
 ! **********************************************************************************************************************************************************************************
 ! **********************************************************************************************************************************************************************************
 ! ***** LOCAL SUBROUTINE vegNrgFlux: COMPUTE ENERGY FLUXES AT VEGETATION AND GROUND SURFACES
 ! **********************************************************************************************************************************************************************************
 ! **********************************************************************************************************************************************************************************
 ! **********************************************************************************************************************************************************************************

 subroutine vegNrgFlux_muster(&

                              ! input
                              dt,                            & ! intent(in): time step (seconds)
                              iter,                          & ! intent(in): iteration index
                              firstSubStep,                  & ! intent(in): flag to indicate if we are processing the first sub-step
                              canopyTempTrial,               & ! intent(in): trial value of canopy temperature (K)
                              groundTempTrial,               & ! intent(in): trial value of ground temperature (K)
                              canopyLiqTrial,                & ! intent(in): trial value of mass of liquid water on the vegetation canopy (kg m-2)
                              canopyIceTrial,                & ! intent(in): trial value of mass of ice on the vegetation canopy (kg m-2)

                              ! model control -- intent(in)
                              ix_fDerivMeth,                 & ! intent(in): choice of method to compute derivatives
                              ix_astability,                 & ! intent(in): choice of stability function
                              ix_soilStress,                 & ! intent(in): choice of function for the soil moisture control on stomatal resistance
                              ix_groundwatr,                 & ! intent(in): groundwater parameterization
                              ix_stomResist,                 & ! intent(in): choice of function for stomatal resistance

                              ! time -- intent(in)
                              yearLength,                    & ! intent(in): number of days in the current year
                              fracJulday,                    & ! intent(in): fractional julian days since the start of year

                              ! physical attributes -- intent(in)
                              mHeight,                       & ! intent(in): measurement height (m)
                              latitude,                      & ! intent(in): latitude (degrees north)
                              vegTypeIndex,                  & ! intent(in): vegetation type index
                              urbanVegCategory,              & ! intent(in): vegetation category for urban areas

                              ! model parameters (phenology) -- intent(in)
                              canopyHeight,                  & ! intent(in): height of the vegetation canopy (m)
                              maxCanopyLiquid,               & ! intent(in): maximum storage of liquid water on the vegetation canopy (kg m-2)
                              maxCanopyIce,                  & ! intent(in): maximum storage of ice on the vegetation canopy (kg m-2)

                              ! model parameters (aerodynamic resistance) -- intent(in)
                              z0Snow,                        & ! intent(in): roughness length of snow (m)
                              z0Soil,                        & ! intent(in): roughness length of soil (m)
                              z0Canopy,                      & ! intent(in): roughness length of the canopy (m)
                              critRichNumber,                & ! intent(in): critical value for the bulk Richardson number where turbulence ceases (-)
                              Louis79_bparam,                & ! intent(in): parameter in Louis (1979) stability function
                              Louis79_cStar,                 & ! intent(in): parameter in Louis (1979) stability function
                              Mahrt87_eScale,                & ! intent(in): exponential scaling factor in the Mahrt (1987) stability function
                              windReductionParam,            & ! intent(in): canopy wind reduction parameter (-)                   
                              leafExchangeCoeff,             & ! intent(in): turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
                              leafDimension,                 & ! intent(in): characteristic leaf dimension (m)

                              ! model parameters (soil stress) -- intent(in)
                              theta_sat,                     & ! intent(in): soil porosity (-)
                              plantWiltPsi,                  & ! intent(in): matric head at wilting point (m)
                              soilStressParam,               & ! intent(in): parameter in the exponential soil stress function (-)
                              critSoilWilting,               & ! intent(in): critical vol. liq. water content when plants are wilting (-)
                              critSoilTranspire,             & ! intent(in): critical vol. liq. water content when transpiration is limited (-)
                              critAquiferTranspire,          & ! intent(in): critical aquifer storage value when transpiration is limited (m)

                              ! forcing at the upper boundary -- intent(in)
                              airtemp,                       & ! intent(in): air temperature at some height above the surface (K)
                              windspd,                       & ! intent(in): wind speed at some height above the surface (m s-1)
                              airpres,                       & ! intent(in): air pressure at some height above the surface (Pa)
                              LWRadAtm,                      & ! intent(in): downwelling longwave radiation at the upper boundary (W m-2)
                              scalarVPair,                   & ! intent(in): vapor pressure at some height above the surface (Pa)
                              scalarO2air,                   & ! intent(in): atmospheric o2 concentration (Pa)
                              scalarCO2air,                  & ! intent(in): atmospheric co2 concentration (Pa)
                              scalarSnowfall,                & ! intent(in): snowfall rate (kg m-2 s-1)
                              scalarCosZenith,               & ! intent(in): cosine of the solar zenith angle (0-1)
                              spectralIncomingDirect,        & ! intent(in): incoming direct solar radiation in each wave band (w m-2)
                              spectralIncomingDiffuse,       & ! intent(in): incoming diffuse solar radiation in each wave band (w m-2)

                              ! water storage -- intent(in)
                              ! NOTE: soil stress only computed at the start of the substep (iter==1)
                              scalarSWE,                     & ! intent(in): snow water equivalent on the ground (kg m-2)
                              scalarSnowDepth,               & ! intent(in): snow depth on the ground surface (m)
                              mLayerVolFracLiq,              & ! intent(in): volumetric fraction of liquid water in each soil layer (-)
                              mLayerMatricHead,              & ! intent(in): matric head in each layer (m)
                              scalarAquiferStorage,          & ! intent(in): aquifer storage (m)

                              ! vegetation phenology -- intent(inout) because only called at the start of the sub-step
                              scalarCanopyWetFraction,       & ! intent(inout): fraction of canopy that is wet
                              scalarFoliageNitrogenFactor,   & ! intent(inout): foliage nitrogen concentration (1.0 = saturated)
                              scalarRootZoneTemp,            & ! intent(inout): root zone temperature
                              scalarLAI,                     & ! intent(inout): one-sided leaf area index (m2 m-2)
                              scalarSAI,                     & ! intent(inout): one-sided stem area index (m2 m-2)
                              scalarExposedLAI,              & ! intent(inout): exposed leaf area index after burial by snow (m2 m-2)
                              scalarExposedSAI,              & ! intent(inout): exposed stem area index after burial by snow (m2 m-2)
                              scalarGrowingSeasonIndex,      & ! intent(inout): growing season index (0=off, 1=on)

                              ! shortwave radiation fluxes -- intent(inout) because only called at the start of the sub-step
                              scalarAlbedo,                  & ! intent(inout): snow albedo (-)
                              scalarSnowAge,                 & ! intent(inout): non-dimensional snow age (-)
                              scalarCanopySunlitFraction,    & ! intent(inout): sunlit fraction of canopy (-)
                              scalarCanopySunlitLAI,         & ! intent(inout): sunlit leaf area (-)
                              scalarCanopyShadedLAI,         & ! intent(inout): shaded leaf area (-)
                              scalarCanopySunlitPAR,         & ! intent(inout): average absorbed par for sunlit leaves (w m-2)
                              scalarCanopyShadedPAR,         & ! intent(inout): average absorbed par for shaded leaves (w m-2)
                              scalarCanopyAbsorbedSolar,     & ! intent(inout): solar radiation absorbed by canopy (W m-2)
                              scalarGroundAbsorbedSolar,     & ! intent(inout): solar radiation absorbed by ground (W m-2)
                              scalarTotalReflectedSolar,     & ! intent(inout): total reflected solar radiation (W m-2)
                              scalarTotalAbsorbedSolar,      & ! intent(inout): total absorbed solar radiation (W m-2)
                              scalarCanopyReflectedSolar,    & ! intent(inout): solar radiation reflected from the canopy (W m-2)
                              scalarGroundReflectedSolar,    & ! intent(inout): solar radiation reflected from the ground (W m-2) 
                              scalarBetweenCanopyGapFraction,& ! intent(inout): between canopy gap fraction for beam (-)
                              scalarWithinCanopyGapFraction, & ! intent(inout): within canopy gap fraction for beam (-)

                              ! longwave radiation fluxes -- intent(out) because called in every step
                              scalarLWRadCanopy,             & ! intent(out): longwave radiation emitted from the canopy (W m-2)
                              scalarLWRadGround,             & ! intent(out): longwave radiation emitted at the ground surface (W m-2)
                              scalarLWRadUbound2Canopy,      & ! intent(out): downward atmospheric longwave radiation absorbed by the canopy (W m-2)
                              scalarLWRadUbound2Ground,      & ! intent(out): downward atmospheric longwave radiation absorbed by the ground (W m-2)
                              scalarLWRadUbound2Ubound,      & ! intent(out): atmospheric radiation reflected by the ground and lost thru upper boundary (W m-2)
                              scalarLWRadCanopy2Ubound,      & ! intent(out): longwave radiation emitted from canopy lost thru upper boundary (W m-2)
                              scalarLWRadCanopy2Ground,      & ! intent(out): longwave radiation emitted from canopy absorbed by the ground (W m-2)
                              scalarLWRadCanopy2Canopy,      & ! intent(out): canopy longwave reflected from ground and absorbed by the canopy (W m-2)
                              scalarLWRadGround2Ubound,      & ! intent(out): longwave radiation emitted from ground lost thru upper boundary (W m-2)
                              scalarLWRadGround2Canopy,      & ! intent(out): longwave radiation emitted from ground and absorbed by the canopy (W m-2)
                              scalarLWNetCanopy,             & ! intent(out): net longwave radiation at the canopy (W m-2)
                              scalarLWNetGround,             & ! intent(out): net longwave radiation at the ground surface (W m-2)
                              scalarLWNetUbound,             & ! intent(out): net longwave radiation at the upper boundary (W m-2)

                              ! aerodynamic resistance -- intent(out) because called in every step
                              scalarWindReductionFactor,     & ! intent(out): canopy wind reduction factor (-)
                              scalarZeroPlaneDisplacement,   & ! intent(out): zero plane displacement (m) 
                              scalarSfc2AtmExchangeCoeff,    & ! intent(out): surface-atmosphere turbulent exchange coefficient (-)
                              scalarEddyDiffusCanopyTop,     & ! intent(out): eddy diffusivity for heat at the top of the canopy (m2 s-1)
                              scalarWindspdCanopyTop,        & ! intent(out): windspeed at the top of the canopy (m s-1)
                              scalarLeafResistance,          & ! intent(out): mean leaf boundary layer resistance per unit leaf area (s m-1)
                              scalarGroundResistance,        & ! intent(out): below canopy aerodynamic resistance (s m-1) 
                              scalarCanopyResistance,        & ! intent(out): above canopy aerodynamic resistance (s m-1)

                              ! soil resistance -- intent(in) and intent(inout) because only called at the first iteration
                              mLayerRootDensity,             & ! intent(in): root density in each layer (-)
                              scalarAquiferRootFrac,         & ! intent(in): fraction of roots below the lowest soil layer (-)
                              scalarTranspireLim,            & ! intent(inout): weighted average of the transpiration limiting factor (-)
                              mLayerTranspireLim,            & ! intent(inout): transpiration limiting factor in each layer (-)
                              scalarTranspireLimAqfr,        & ! intent(inout): transpiration limiting factor for the aquifer (-)

                              ! stomatal resistance -- intent(inout) because only called at the first iteration
                              scalarStomResistSunlit,        & ! intent(inout): stomatal resistance for sunlit leaves (s m-1)
                              scalarStomResistShaded,        & ! intent(inout): stomatal resistance for shaded leaves (s m-1)
                              scalarPhotosynthesisSunlit,    & ! intent(inout): sunlit photosynthesis (umolco2 m-2 s-1)
                              scalarPhotosynthesisShaded,    & ! intent(inout): shaded photosynthesis (umolco2 m-2 s-1)

                              ! turbulent heat fluxes
                              scalarPsycConstCanopy,         & ! intent(inout): psychrometric constant for the vegetation canopy (Pa K-1)
                              scalarPsycConstGround,         & ! intent(inout): psychrometric constant for the ground surface (Pa K-1)
                              scalarTemp_CanopyAir,          & ! intent(inout): temperature of the canopy air space (K)
                              scalarVP_CanopyAir,            & ! intent(inout): vapor pressure of the canopy air space (Pa)
                              scalarSatVP_canopyTemp,        & ! intent(out): saturation vapor pressure at the temperature of the vegetation canopy (Pa)
                              scalarSatVP_groundTemp,        & ! intent(out): saturation vapor pressure at the temperature of the ground surface (Pa)
                              scalarSoilRelHumidity,         & ! intent(out): relative humidity in the soil pores [0-1]
                              scalarSoilResistance,          & ! intent(out): resistance from the soil (s m-1)
                              scalarSenHeatCanopy,           & ! intent(out): sensible heat flux from the canopy to the canopy air space (W m-2)
                              scalarLatHeatCanopyEvap,       & ! intent(out): latent heat flux associated with evaporation from the canopy to the canopy air space (W m-2)
                              scalarLatHeatCanopyTrans,      & ! intent(out): latent heat flux associated with transpiration from the canopy to the canopy air space (W m-2)
                              scalarSenHeatGround,           & ! intent(out): sensible heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
                              scalarLatHeatGround,           & ! intent(out): latent heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)

                              ! output
                              canopyNetFlux,                 & ! intent(out): net energy flux for the vegetation canopy (W m-2)
                              groundNetFlux,                 & ! intent(out): net energy flux for the ground surface (W m-2)
                              dCanopyNetFlux_dCanopyTemp,    & ! intent(out): derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
                              dCanopyNetFlux_dGroundTemp,    & ! intent(out): derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
                              dGroundNetFlux_dCanopyTemp,    & ! intent(out): derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
                              dGroundNetFlux_dGroundTemp,    & ! intent(out): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
                              err,message                    ) ! intent(out): error control
 ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! conversion functions
 USE conv_funcs_module,only:satVapPress                           ! function to compute the saturated vapor pressure (Pa)
 USE conv_funcs_module,only:psychometric                          ! function to compute the psychrometric constant (Pa K-1)
 ! Noah-MP modules
 USE NOAHMP_ROUTINES,only:phenology                               ! compute vegetation phenology
 USE NOAHMP_ROUTINES,only:radiation                               ! compute radiation fluxes
 ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! input/output
 ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! input
 real(dp),intent(in)            :: dt                             ! time step (seconds)
 integer(i4b),intent(in)        :: iter                           ! iteration index
 logical(lgt),intent(in)        :: firstSubStep                   ! flag to indicate if we are processing the first sub-step
 real(dp),intent(in)            :: canopyTempTrial                ! trial value of canopy temperature (K)
 real(dp),intent(in)            :: groundTempTrial                ! trial value of ground temperature (K)
 real(dp),intent(in)            :: canopyLiqTrial                 ! trial value of mass of liquid water on the vegetation canopy (kg m-2)
 real(dp),intent(in)            :: canopyIceTrial                 ! trial value of mass of ice on the vegetation canopy (kg m-2)

 ! model control -- intent(in)
 integer(i4b),intent(in)        :: ix_fDerivMeth                  ! choice of method to compute derivatives
 integer(i4b),intent(in)        :: ix_astability                  ! choice of stability function
 integer(i4b),intent(in)        :: ix_soilStress                  ! choice of function for the soil moisture control on stomatal resistance
 integer(i4b),intent(in)        :: ix_groundwatr                  ! groundwater parameterization
 integer(i4b),intent(in)        :: ix_stomResist                  ! choice of function for stomatal resistance

 ! time -- intent(in)
 integer(i4b),intent(in)        :: yearLength                     ! number of days in the current year
 real(dp),intent(in)            :: fracJulday                     ! fractional julian days since the start of year

 ! physical attributes -- intent(in)
 real(dp),intent(in)            :: mHeight                        ! measurement height (m)
 real(dp),intent(in)            :: latitude                       ! latitude (degrees north)
 integer(i4b),intent(in)        :: vegTypeIndex                   ! vegetation type index
 integer(i4b),intent(in)        :: urbanVegCategory               ! vegetation category for urban areas

 ! model parameters (phenology) -- intent(in)
 real(dp),intent(in)            :: canopyHeight                   ! height of the vegetation canopy (m)
 real(dp),intent(in)            :: maxCanopyLiquid                ! maximum storage of liquid water on the vegetation canopy (kg m-2)
 real(dp),intent(in)            :: maxCanopyIce                   ! maximum storage of ice on the vegetation canopy (kg m-2)

 ! model parameters (aerodynamic resistance) -- intent(in)
 real(dp),intent(in)            :: z0Snow                         ! roughness length of snow (m)
 real(dp),intent(in)            :: z0Soil                         ! roughness length of soil (m)
 real(dp),intent(in)            :: z0Canopy                       ! roughness length of the canopy (m)
 real(dp),intent(in)            :: critRichNumber                 ! critical value for the bulk Richardson number where turbulence ceases (-)
 real(dp),intent(in)            :: Louis79_bparam                 ! parameter in Louis (1979) stability function
 real(dp),intent(in)            :: Louis79_cStar                  ! parameter in Louis (1979) stability function
 real(dp),intent(in)            :: Mahrt87_eScale                 ! exponential scaling factor in the Mahrt (1987) stability function
 real(dp),intent(in)            :: windReductionParam             ! canopy wind reduction parameter (-)                   
 real(dp),intent(in)            :: leafExchangeCoeff              ! turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
 real(dp),intent(in)            :: leafDimension                  ! characteristic leaf dimension (m)

 ! model parameters (soil stress) -- intent(in)
 real(dp),intent(in)            :: theta_sat                      ! soil porosity (-)
 real(dp),intent(in)            :: plantWiltPsi                   ! matric head at wilting point (m)
 real(dp),intent(in)            :: soilStressParam                ! parameter in the exponential soil stress function (-)
 real(dp),intent(in)            :: critSoilWilting                ! critical vol. liq. water content when plants are wilting (-)
 real(dp),intent(in)            :: critSoilTranspire              ! critical vol. liq. water content when transpiration is limited (-)
 real(dp),intent(in)            :: critAquiferTranspire           ! critical aquifer storage value when transpiration is limited (m)

 ! forcing at the upper boundary -- intent(in)
 real(dp),intent(in)            :: airtemp                        ! air temperature at some height above the surface (K)
 real(dp),intent(in)            :: windspd                        ! wind speed at some height above the surface (m s-1)
 real(dp),intent(in)            :: airpres                        ! air pressure at some height above the surface (Pa)
 real(dp),intent(in)            :: LWRadAtm                       ! downwelling longwave radiation at the upper boundary (W m-2)
 real(dp),intent(in)            :: scalarVPair                    ! vapor pressure at some height above the surface (Pa)
 real(dp),intent(in)            :: scalarO2air                    ! atmospheric o2 concentration (Pa)
 real(dp),intent(in)            :: scalarCO2air                   ! atmospheric co2 concentration (Pa)
 real(dp),intent(in)            :: scalarSnowfall                 ! snowfall rate (kg m-2 s-1)
 real(dp),intent(in)            :: scalarCosZenith                ! cosine of the solar zenith angle (0-1)
 real(dp),intent(in)            :: spectralIncomingDirect(:)      ! incoming direct solar radiation in each wave band (w m-2)
 real(dp),intent(in)            :: spectralIncomingDiffuse(:)     ! incoming diffuse solar radiation in each wave band (w m-2)

 ! water storage -- intent(in)
 ! NOTE: soil stress only computed at the start of the substep (iter==1)
 real(dp),intent(in)            :: scalarSWE                      ! snow water equivalent on the ground (kg m-2)
 real(dp),intent(in)            :: scalarSnowDepth                ! snow depth on the ground surface (m)
 real(dp),intent(in)            :: mLayerVolFracLiq(:)            ! volumetric fraction of liquid water in each soil layer (-)
 real(dp),intent(in)            :: mLayerMatricHead(:)            ! matric head in each layer (m)
 real(dp),intent(in)            :: scalarAquiferStorage           ! aquifer storage (m)

 ! vegetation phenology -- intent(inout) because only called at the start of the sub-step
 real(dp),intent(inout)         :: scalarCanopyWetFraction        ! fraction of canopy that is wet
 real(dp),intent(inout)         :: scalarFoliageNitrogenFactor    ! foliage nitrogen concentration (1.0 = saturated)
 real(dp),intent(inout)         :: scalarRootZoneTemp             ! root zone temperature
 real(dp),intent(inout)         :: scalarLAI                      ! one-sided leaf area index (m2 m-2)
 real(dp),intent(inout)         :: scalarSAI                      ! one-sided stem area index (m2 m-2)
 real(dp),intent(inout)         :: scalarExposedLAI               ! exposed leaf area index after burial by snow (m2 m-2)
 real(dp),intent(inout)         :: scalarExposedSAI               ! exposed stem area index after burial by snow (m2 m-2)
 real(dp),intent(inout)         :: scalarGrowingSeasonIndex       ! growing season index (0=off, 1=on)

 ! shortwave radiation fluxes -- intent(inout) because only called at the start of the sub-step
 real(dp),intent(inout)         :: scalarAlbedo                   ! snow albedo (-)
 real(dp),intent(inout)         :: scalarSnowAge                  ! non-dimensional snow age (-)
 real(dp),intent(inout)         :: scalarCanopySunlitFraction     ! sunlit fraction of canopy (-)
 real(dp),intent(inout)         :: scalarCanopySunlitLAI          ! sunlit leaf area (-)
 real(dp),intent(inout)         :: scalarCanopyShadedLAI          ! shaded leaf area (-)
 real(dp),intent(inout)         :: scalarCanopySunlitPAR          ! average absorbed par for sunlit leaves (w m-2)
 real(dp),intent(inout)         :: scalarCanopyShadedPAR          ! average absorbed par for shaded leaves (w m-2)
 real(dp),intent(inout)         :: scalarCanopyAbsorbedSolar      ! solar radiation absorbed by canopy (W m-2)
 real(dp),intent(inout)         :: scalarGroundAbsorbedSolar      ! solar radiation absorbed by ground (W m-2)
 real(dp),intent(inout)         :: scalarTotalReflectedSolar      ! total reflected solar radiation (W m-2)
 real(dp),intent(inout)         :: scalarTotalAbsorbedSolar       ! total absorbed solar radiation (W m-2)
 real(dp),intent(inout)         :: scalarCanopyReflectedSolar     ! solar radiation reflected from the canopy (W m-2)
 real(dp),intent(inout)         :: scalarGroundReflectedSolar     ! solar radiation reflected from the ground (W m-2) 
 real(dp),intent(inout)         :: scalarBetweenCanopyGapFraction ! between canopy gap fraction for beam (-)
 real(dp),intent(inout)         :: scalarWithinCanopyGapFraction  ! within canopy gap fraction for beam (-)

 ! longwave radiation fluxes -- intent(out) because called in every step
 real(dp),intent(out)           :: scalarLWRadCanopy              ! longwave radiation emitted from the canopy (W m-2)
 real(dp),intent(out)           :: scalarLWRadGround              ! longwave radiation emitted at the ground surface (W m-2)
 real(dp),intent(out)           :: scalarLWRadUbound2Canopy       ! downward atmospheric longwave radiation absorbed by the canopy (W m-2)
 real(dp),intent(out)           :: scalarLWRadUbound2Ground       ! downward atmospheric longwave radiation absorbed by the ground (W m-2)
 real(dp),intent(out)           :: scalarLWRadUbound2Ubound       ! atmospheric radiation reflected by the ground and lost thru upper boundary (W m-2)
 real(dp),intent(out)           :: scalarLWRadCanopy2Ubound       ! longwave radiation emitted from canopy lost thru upper boundary (W m-2)
 real(dp),intent(out)           :: scalarLWRadCanopy2Ground       ! longwave radiation emitted from canopy absorbed by the ground (W m-2)
 real(dp),intent(out)           :: scalarLWRadCanopy2Canopy       ! canopy longwave reflected from ground and absorbed by the canopy (W m-2)
 real(dp),intent(out)           :: scalarLWRadGround2Ubound       ! longwave radiation emitted from ground lost thru upper boundary (W m-2)
 real(dp),intent(out)           :: scalarLWRadGround2Canopy       ! longwave radiation emitted from ground and absorbed by the canopy (W m-2)
 real(dp),intent(out)           :: scalarLWNetCanopy              ! net longwave radiation at the canopy (W m-2)
 real(dp),intent(out)           :: scalarLWNetGround              ! net longwave radiation at the ground surface (W m-2)
 real(dp),intent(out)           :: scalarLWNetUbound              ! net longwave radiation at the upper boundary (W m-2)

 ! aerodynamic resistance -- intent(out) because called in every step
 real(dp),intent(out)           :: scalarWindReductionFactor      ! canopy wind reduction factor (-)
 real(dp),intent(out)           :: scalarZeroPlaneDisplacement    ! zero plane displacement (m) 
 real(dp),intent(out)           :: scalarSfc2AtmExchangeCoeff     ! surface-atmosphere turbulent exchange coefficient (-)
 real(dp),intent(out)           :: scalarEddyDiffusCanopyTop      ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
 real(dp),intent(out)           :: scalarWindspdCanopyTop         ! windspeed at the top of the canopy (m s-1)
 real(dp),intent(out)           :: scalarLeafResistance           ! mean leaf boundary layer resistance per unit leaf area (s m-1)
 real(dp),intent(out)           :: scalarGroundResistance         ! below canopy aerodynamic resistance (s m-1) 
 real(dp),intent(out)           :: scalarCanopyResistance         ! above canopy aerodynamic resistance (s m-1)

 ! soil resistance -- intent(in) and intent(inout) because only called at the first iteration
 real(dp),intent(in)            :: mLayerRootDensity(:)           ! root density in each layer (-)
 real(dp),intent(in)            :: scalarAquiferRootFrac          ! fraction of roots below the lowest soil layer (-)
 real(dp),intent(out)           :: scalarTranspireLim             ! weighted average of the transpiration limiting factor (-)
 real(dp),intent(out)           :: mLayerTranspireLim(:)          ! transpiration limiting factor in each layer (-)
 real(dp),intent(out)           :: scalarTranspireLimAqfr         ! transpiration limiting factor for the aquifer (-)

 ! stomatal resistance -- intent(inout) because only called at the first iteration
 real(dp),intent(inout)         :: scalarStomResistSunlit         ! stomatal resistance for sunlit leaves (s m-1)
 real(dp),intent(inout)         :: scalarStomResistShaded         ! stomatal resistance for shaded leaves (s m-1)
 real(dp),intent(inout)         :: scalarPhotosynthesisSunlit     ! sunlit photosynthesis (umolco2 m-2 s-1)
 real(dp),intent(inout)         :: scalarPhotosynthesisShaded     ! shaded photosynthesis (umolco2 m-2 s-1)

 ! turbulent heat fluxes
 real(dp),intent(inout)         :: scalarPsycConstCanopy          ! psychrometric constant for the vegetation canopy (Pa K-1)
 real(dp),intent(inout)         :: scalarPsycConstGround          ! psychrometric constant for the ground surface (Pa K-1)
 real(dp),intent(inout)         :: scalarTemp_CanopyAir           ! temperature of the canopy air space (K)
 real(dp),intent(inout)         :: scalarVP_CanopyAir             ! vapor pressure of the canopy air space (Pa)
 real(dp),intent(out)           :: scalarSatVP_canopyTemp         ! saturation vapor pressure at the temperature of the vegetation canopy (Pa)
 real(dp),intent(out)           :: scalarSatVP_groundTemp         ! saturation vapor pressure at the temperature of the ground surface (Pa)
 real(dp),intent(out)           :: scalarSoilRelHumidity          ! relative humidity in the soil pores [0-1]
 real(dp),intent(out)           :: scalarSoilResistance           ! resistance from the soil (s m-1)
 real(dp),intent(out)           :: scalarSenHeatCanopy            ! sensible heat flux from the canopy to the canopy air space (W m-2)
 real(dp),intent(out)           :: scalarLatHeatCanopyEvap        ! latent heat flux associated with evaporation from the canopy to the canopy air space (W m-2)
 real(dp),intent(out)           :: scalarLatHeatCanopyTrans       ! latent heat flux associated with transpiration from the canopy to the canopy air space (W m-2)
 real(dp),intent(out)           :: scalarSenHeatGround            ! sensible heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
 real(dp),intent(out)           :: scalarLatHeatGround            ! latent heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)

 ! output
 real(dp),intent(out)           :: canopyNetFlux                  ! net energy flux for the vegetation canopy (W m-2)
 real(dp),intent(out)           :: groundNetFlux                  ! net energy flux for the ground surface (W m-2)
 real(dp),intent(out)           :: dCanopyNetFlux_dCanopyTemp     ! derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp),intent(out)           :: dCanopyNetFlux_dGroundTemp     ! derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
 real(dp),intent(out)           :: dGroundNetFlux_dCanopyTemp     ! derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp),intent(out)           :: dGroundNetFlux_dGroundTemp     ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
 integer(i4b),intent(out)       :: err                            ! error code
 character(*),intent(out)       :: message                        ! error message

 ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! local (general)
 character(LEN=256)            :: cmessage                        ! error message of downwind routine
 real(dp)                      :: snowmassPlusNewsnow             ! sum of snow mass and new snowfall (kg m-2 [mm])
 real(dp)                      :: fracSnow                        ! snow cover fraction (0-1)
 real(dp)                      :: VAI                             ! vegetation area index (m2 m-2)
 real(dp)                      :: exposedVAI                      ! "exposed" vegetation area index (m2 m-2)
 real(dp)                      :: greenVegFraction                ! green vegetation fraction (0-1) 
 real(dp)                      :: relativeCanopyWater             ! water stored on vegetation canopy, expressed as a fraction of maximum storage (-)
 ! local (compute numerical derivatives)
 integer(i4b),parameter        :: unperturbed=1                   ! named variable to identify the case of unperturbed state variables
 integer(i4b),parameter        :: perturbStateCanopy=2            ! named variable to identify the case where we perturb the canopy temperature
 integer(i4b),parameter        :: perturbStateGround=3            ! named variable to identify the case where we perturb the ground temperature
 integer(i4b)                  :: itry                            ! index of flux evaluation
 integer(i4b)                  :: nFlux                           ! number of flux evaluations
 real(dp)                      :: canopyTemp                      ! value of canopy temperature used in flux calculations (may be perturbed)
 real(dp)                      :: groundTemp                      ! value of ground temperature used in flux calculations (may be perturbed)
 real(dp)                      :: try0,try1,try2                  ! trial values to evaluate specific derivatives (testing only)
 ! local (phenology)
 real(dp)                      :: notUsed_canopyHeight            ! for some reason the Noah-MP phenology routines output canopy height
 ! local (saturation vapor pressure of veg)
 real(dp)                      :: TV_celcius                      ! vegetaion temperature (C)
 real(dp)                      :: TG_celcius                      ! ground temperature (C)
 real(dp)                      :: dSVPCanopy_dCanopyTemp          ! derivative in canopy saturated vapor pressure w.r.t. vegetation temperature (Pa/K)
 real(dp)                      :: dSVPGround_dGroundTemp          ! derivative in ground saturated vapor pressure w.r.t. ground temperature (Pa/K)
 ! local (longwave radiation)
 real(dp)                      :: canopyEmissivity                ! effective emissivity of the canopy (-)
 real(dp)                      :: groundEmissivity                ! emissivity of the ground surface (-)
 real(dp),parameter            :: soilEmissivity=0.98_dp          ! emmisivity of the soil (-)
 real(dp),parameter            :: snowEmissivity=0.99_dp          ! emissivity of snow (-)
 real(dp)                      :: dLWNetCanopy_dTCanopy           ! derivative in net canopy radiation w.r.t. canopy temperature (W m-2 K-1)
 real(dp)                      :: dLWNetGround_dTGround           ! derivative in net ground radiation w.r.t. ground temperature (W m-2 K-1)
 real(dp)                      :: dLWNetCanopy_dTGround           ! derivative in net canopy radiation w.r.t. ground temperature (W m-2 K-1)
 real(dp)                      :: dLWNetGround_dTCanopy           ! derivative in net ground radiation w.r.t. canopy temperature (W m-2 K-1)
 ! local (turbulent heat transfer)
 real(dp)                      :: z0Ground                        ! roughness length of the ground (ground below the canopy or non-vegetated surface) (m)
 real(dp)                      :: soilEvapFactor                  ! soil water control on evaporation from non-vegetated surfaces
 real(dp)                      :: soilRelHumidity_noSnow          ! relative humidity in the soil pores [0-1]
 real(dp)                      :: dGroundResistance_dTCanopy      ! derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
 real(dp)                      :: dCanopyResistance_dTCanopy      ! derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
 real(dp)                      :: dLeafResistance_dTCanopy        ! derivative in leaf resistance w.r.t. canopy temperature (s m-1 K-1)
 real(dp)                      :: turbFluxCanopy                  ! total turbulent heat fluxes from the canopy to the canopy air space (W m-2)
 real(dp)                      :: turbFluxGround                  ! total turbulent heat fluxes from the ground to the canopy air space (W m-2)
 ! local (turbulent heat transfer -- compute numerical derivatives)
 ! (temporary scalar resistances when states are perturbed)
 real(dp)                      :: trialLeafResistance             ! mean leaf boundary layer resistance per unit leaf area (s m-1)
 real(dp)                      :: trialGroundResistance           ! below canopy aerodynamic resistance (s m-1) 
 real(dp)                      :: trialCanopyResistance           ! above canopy aerodynamic resistance (s m-1)
 real(dp)                      :: notUsed_WindReductionFactor     ! canopy wind reduction factor (-)
 real(dp)                      :: notUsed_ZeroPlaneDisplacement   ! zero plane displacement (m) 
 real(dp)                      :: notUsed_Sfc2AtmExchangeCoeff    ! surface-atmosphere turbulent exchange coefficient (-)
 real(dp)                      :: notUsed_EddyDiffusCanopyTop     ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
 real(dp)                      :: notUsed_WindspdCanopyTop        ! windspeed at the top of the canopy (m s-1)
 real(dp)                      :: notUsed_dGroundResistance_dTCanopy  ! derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
 real(dp)                      :: notUsed_dCanopyResistance_dTCanopy  ! derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
 real(dp)                      :: notUsed_dLeafResistance_dTCanopy    ! derivative in leaf resistance w.r.t. canopy temperature (s m-1 K-1)
 ! (fluxes after perturbations in model states)
 real(dp)                      :: turbFluxCanopy_dStateCanopy     ! total turbulent heat fluxes from the canopy to the canopy air space, after canopy temperature is perturbed (W m-2)
 real(dp)                      :: turbFluxGround_dStateCanopy     ! total turbulent heat fluxes from the ground to the canopy air space, after canopy temperature is perturbed (W m-2)
 real(dp)                      :: turbFluxCanopy_dStateGround     ! total turbulent heat fluxes from the canopy to the canopy air space, after ground temperature is perturbed (W m-2)
 real(dp)                      :: turbFluxGround_dStateGround     ! total turbulent heat fluxes from the ground to the canopy air space, after ground temperature is perturbed (W m-2)
 ! (flux derivatives)
 real(dp)                      :: dTurbFluxCanopy_dTCanopy        ! derivative in net canopy turbulent fluxes w.r.t. canopy temperature (W m-2 K-1) 
 real(dp)                      :: dTurbFluxGround_dTCanopy        ! derivative in net ground turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
 real(dp)                      :: dTurbFluxCanopy_dTGround        ! derivative in net canopy turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
 real(dp)                      :: dTurbFluxGround_dTGround        ! derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
 ! initialize error control
 err=0; message="vegNrgFlux_muster/"

 ! initialize variables to compute stomatal resistance
 if(iter==1 .and. firstSubStep)then
  ! vapor pressure in the canopy air space initialized as vapor pressure of air above the vegetation canopy
  scalarVP_CanopyAir = scalarVPair
 endif

 ! *******************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************
 ! ***** PHENOLOGY  **************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************

 ! set psychrometric constant for canopy and ground surface (Pa/K)
 ! NOTE: variables are constant over the full time step, to simplify relating energy and mass fluxes
 if(iter==1 .and. firstSubStep)then
  scalarPsycConstCanopy = psychometric(canopyTempTrial,airpres)
  scalarPsycConstGround = psychometric(groundTempTrial,airpres)
 endif

 ! variables are constant over the SUBSTEP
 ! NOTE: recomputing phenology every sub-step accounts for changes in exposed vegetation associated with changes in snow depth
 if(iter==1)then

  ! define the foliage nitrogen factor
  scalarFoliageNitrogenFactor = 1._dp  ! foliage nitrogen concentration (1.0 = saturated)

  ! compute vegetation phenology
  call phenology(&
                 ! input
                 vegTypeIndex,                       & ! intent(in): vegetation type index
                 urbanVegCategory,                   & ! intent(in): vegetation category for urban areas
                 scalarSnowDepth,                    & ! intent(in): snow depth (m)
                 canopyTempTrial,                    & ! intent(in): vegetation temperature (K)
                 latitude,                           & ! intent(in): latitude (degrees north)
                 yearLength,                         & ! intent(in): number of days in the current year
                 fracJulday,                         & ! intent(in): fractional julian days since the start of year
                 scalarLAI,                          & ! intent(inout): one-sided leaf area index (m2 m-2)
                 scalarSAI,                          & ! intent(inout): one-sided stem area index (m2 m-2)
                 scalarRootZoneTemp,                 & ! intent(in): average temperature of the root zone (K)
                 ! output
                 notUsed_canopyHeight,               & ! intent(out): height of the top of the canopy layer (m)
                 scalarExposedLAI,                   & ! intent(out): exposed leaf area index after burial by snow (m2 m-2)
                 scalarExposedSAI,                   & ! intent(out): exposed stem area index after burial by snow (m2 m-2)
                 scalarGrowingSeasonIndex            ) ! intent(out): growing season index (0=off, 1=on)

 endif  ! (variables constant over the SUBSTEP)

 ! compute the snow cover fraction
 if(scalarSWE > 0._dp)then
  fracSnow = 1._dp
 else
  fracSnow = 0._dp
 endif

 ! compute the roughness length of the ground (ground below the canopy or non-vegetated surface)
 z0Ground = z0soil*(1._dp - fracSnow) + z0Snow*fracSnow     ! roughness length (m)

 ! compute the total vegetation area index (leaf plus stem)
 VAI        = scalarLAI + scalarSAI  ! vegetation area index
 exposedVAI = scalarExposedLAI + scalarExposedSAI  !  exposed vegetation area index

 ! compute the green vegetation fraction (used in radiation module -- set to 1 to trun off semi-tile approach)
 greenVegFraction = 1._dp

 ! compute emissivity of the canopy and ground surface (-)
 canopyEmissivity = 1._dp - exp(-exposedVAI)                                     ! effective emissivity of the canopy (-)
 groundEmissivity = fracSnow*snowEmissivity + (1._dp - fracSnow)*soilEmissivity  ! emissivity of the ground surface (-)


 ! *******************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************
 ! ***** SHORTWAVE RADIATION *****************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************

 ! shortwave radiation is constant over the SUBSTEP
 ! NOTE: canopy radiation does depend on canopy temperature (to distinguish between snow covered canopy, for albedo calculations) and this dependency
 !       is not included on the assumption that albedo changes over the sub-step are small compared to other fluxes
 if(iter==1)then

  ! compute the fraction of canopy that is wet
  if(canopyIceTrial > 0._dp)then
   relativeCanopyWater = (canopyIceTrial + canopyLiqTrial) / (maxCanopyIce*exposedVAI)
  else
   relativeCanopyWater = canopyLiqTrial / (maxCanopyLiquid*exposedVAI)
  endif
  scalarCanopyWetFraction = min(relativeCanopyWater, 1._dp)*0.666667_dp

  ! compute the sum of snow mass and new snowfall (kg m-2 [mm])
  snowmassPlusNewsnow = scalarSWE + scalarSnowfall*dt

  ! compute canopy shortwave radiation fluxes
  ! (unchanged Noah-MP routine)
  call radiation(&
                 ! input
                 vegTypeIndex,                       & ! intent(in): vegetation type index
                 ist, isc, ice,                      & ! intent(in): indices to define surface type, soil color, and ice type (constant)
                 nSoil,                              & ! intent(in): number of soil layers               
                 scalarSWE,                          & ! intent(in): snow water equivalent (kg m-2 [mm])
                 snowmassPlusNewsnow,                & ! intent(in): sum of snow mass and new snowfall (kg m-2 [mm])
                 dt,                                 & ! intent(in): time step (s)
                 scalarCosZenith,                    & ! intent(in): cosine of the solar zenith angle (0-1)
                 scalarSnowDepth*1000._dp,           & ! intent(in): snow depth (mm)
                 groundTempTrial,                    & ! intent(in): ground temperature (K)
                 canopyTempTrial,                    & ! intent(in): vegetation temperature (K)
                 fracSnow,                           & ! intent(in): snow cover fraction (0-1)
                 scalarSnowfall,                     & ! intent(in): snowfall (kg m-2 s-1 [mm/s])
                 scalarCanopyWetFraction,            & ! intent(in): fraction of canopy that is wet
                 scalarExposedLAI,                   & ! intent(in): exposed leaf area index after burial by snow (m2 m-2)
                 scalarExposedSAI,                   & ! intent(in): exposed stem area index after burial by snow (m2 m-2)          
                 mLayerVolFracLiq(1:nSoil),          & ! intent(in): volumetric fraction of liquid water in each soil layer (-)
                 spectralIncomingDirect(1:nBands),   & ! intent(in): incoming direct solar radiation in each wave band (w m-2)
                 spectralIncomingDiffuse(1:nBands),  & ! intent(in): incoming diffuse solar radiation in each wave band (w m-2)
                 greenVegFraction,                   & ! intent(in): green vegetation fraction (0-1)
                 iLoc, jLoc,                         & ! intent(in): spatial location indices      
                 ! output
                 scalarAlbedo,                       & ! intent(inout): snow albedo (-)
                 scalarSnowAge,                      & ! intent(inout): non-dimensional snow age (-)
                 scalarCanopySunlitFraction,         & ! intent(out): sunlit fraction of canopy (-)
                 scalarCanopySunlitLAI,              & ! intent(out): sunlit leaf area (-)
                 scalarCanopyShadedLAI,              & ! intent(out): shaded leaf area (-)
                 scalarCanopySunlitPAR,              & ! intent(out): average absorbed par for sunlit leaves (w m-2)
                 scalarCanopyShadedPAR,              & ! intent(out): average absorbed par for shaded leaves (w m-2)
                 scalarCanopyAbsorbedSolar,          & ! intent(out): solar radiation absorbed by canopy (W m-2)
                 scalarGroundAbsorbedSolar,          & ! intent(out): solar radiation absorbed by ground (W m-2)
                 scalarTotalReflectedSolar,          & ! intent(out): total reflected solar radiation (W m-2)
                 scalarTotalAbsorbedSolar,           & ! intent(out): total absorbed solar radiation (W m-2)
                 scalarCanopyReflectedSolar,         & ! intent(out): solar radiation reflected from the canopy (W m-2)
                 scalarGroundReflectedSolar,         & ! intent(out): solar radiation reflected from the ground (W m-2) 
                 scalarBetweenCanopyGapFraction,     & ! intent(out): between canopy gap fraction for beam (-)
                 scalarWithinCanopyGapFraction       ) ! intent(out): within canopy gap fraction for beam (-)
  !print*, 'average absorbed par for sunlit leaves (w m-2) = ', scalarCanopySunlitPAR
  !print*, 'average absorbed par for shaded leaves (w m-2) = ', scalarCanopyShadedPAR
  !print*, 'solar radiation absorbed by canopy (W m-2) = ', scalarCanopyAbsorbedSolar
  !print*, 'solar radiation absorbed by ground (W m-2) = ', scalarGroundAbsorbedSolar

 endif  ! (shortwave radiation is constant over the SUBSTEP)

 ! *******************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************
 ! ***** AERODYNAMIC RESISTANCE *****************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************

 ! NOTE: compute for all iterations

 ! compute aerodynamic resistances
 ! Refs: Choudhury and Monteith (4-layer model for heat budget of homogenous surfaces; QJRMS, 1988)
 !       Niu and Yang (Canopy effects on snow processes; JGR, 2004)
 !       Mahat et al. (Below-canopy turbulence in a snowmelt model, WRR, 2012)
 call aeroResist(&
                 ! input: model control
                 (ix_fDerivMeth ==analytical),       & ! intent(in): logical flag if would like to compute analytical derivaties
                 ix_astability,                      & ! intent(in): choice of stability function
                 ! input: above-canopy forcing data
                 mHeight,                            & ! intent(in): measurement height (m)
                 airtemp,                            & ! intent(in): air temperature at some height above the surface (K)
                 windspd,                            & ! intent(in): wind speed at some height above the surface (m s-1)
                 ! input: canopy and ground temperature
                 canopyTempTrial,                    & ! intent(in): temperature of the vegetation canopy (K)
                 groundTempTrial,                    & ! intent(in): temperature of the ground surface (K)
                 ! input: diagnostic variables
                 exposedVAI,                         & ! intent(in): exposed vegetation area index -- leaf plus stem (m2 m-2)
                 scalarSnowDepth,                    & ! intent(in): snow depth (m)
                 ! input: parameters
                 z0Ground,                           & ! intent(in): roughness length of the ground (below canopy or non-vegetated surface [snow]) (m)
                 z0Canopy,                           & ! intent(in): roughness length of the canopy (m)
                 critRichNumber,                     & ! intent(in): critical value for the bulk Richardson number where turbulence ceases (-)
                 Louis79_bparam,                     & ! intent(in): parameter in Louis (1979) stability function
                 Louis79_cStar,                      & ! intent(in): parameter in Louis (1979) stability function
                 Mahrt87_eScale,                     & ! intent(in): exponential scaling factor in the Mahrt (1987) stability function
                 windReductionParam,                 & ! intent(in): canopy wind reduction parameter (-)                   
                 leafExchangeCoeff,                  & ! intent(in): turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
                 leafDimension,                      & ! intent(in): characteristic leaf dimension (m)
                 canopyHeight,                       & ! intent(in): canopy height (m) 
                 ! output: scalar resistances
                 scalarWindReductionFactor,          & ! intent(out): canopy wind reduction factor (-)
                 scalarZeroPlaneDisplacement,        & ! intent(out): zero plane displacement (m) 
                 scalarSfc2AtmExchangeCoeff,         & ! intent(out): surface-atmosphere turbulent exchange coefficient (-)
                 scalarEddyDiffusCanopyTop,          & ! intent(out): eddy diffusivity for heat at the top of the canopy (m2 s-1)
                 scalarWindspdCanopyTop,             & ! intent(out): windspeed at the top of the canopy (m s-1)
                 scalarLeafResistance,               & ! intent(out): mean leaf boundary layer resistance per unit leaf area (s m-1)
                 scalarGroundResistance,             & ! intent(out): below canopy aerodynamic resistance (s m-1) 
                 scalarCanopyResistance,             & ! intent(out): above canopy aerodynamic resistance (s m-1)
                 ! output: derivatives in scalar resistances
                 dGroundResistance_dTCanopy,         & ! intent(out): derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
                 dCanopyResistance_dTCanopy,         & ! intent(out): derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
                 dLeafResistance_dTCanopy,           & ! intent(out): derivative in leaf resistance w.r.t. canopy temperature (s m-1 K-1)
                 ! output: error control
                 err,cmessage                        ) ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 !print*,         scalarLeafResistance,    & ! mean leaf boundary layer resistance per unit leaf area (s m-1)
 !                scalarGroundResistance,  & ! below canopy aerodynamic resistance (s m-1) 
 !                scalarCanopyResistance,  & ! above canopy aerodynamic resistance (s m-1)
 !                '(leaf, ground, canopy)'

 ! *******************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************
 ! ***** STOMATAL RESISTANCE *****************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************

 ! stomatal resistance is constant over the SUBSTEP
 ! NOTE: This is a simplification, as stomatal resistance does depend on canopy temperature
 !       This "short-cut" made because:
 !         (1) computations are expensive;
 !         (2) derivative calculations are rather complex (iterations within the Ball-Berry routine); and
 !         (3) stomatal resistance does not change rapidly
 if(iter == 1)then

  ! compute the saturation vapor pressure for vegetation temperature
  TV_celcius = canopyTempTrial - Tfreeze
  call satVapPress(TV_celcius, scalarSatVP_CanopyTemp, dSVPCanopy_dCanopyTemp)

  ! compute soil moisture factor controlling stomatal resistance
  call soilResist(&
                  ! input (model decisions)
                  ix_soilStress,                     & ! intent(in): choice of function for the soil moisture control on stomatal resistance
                  ix_groundwatr,                     & ! intent(in): groundwater parameterization
                  ! input (state variables)
                  mLayerMatricHead(1:nSoil),         & ! intent(in): matric head in each layer (m)
                  mLayerVolFracLiq(1:nSoil),         & ! intent(in): volumetric fraction of liquid water in each layer (-)
                  scalarAquiferStorage,              & ! intent(in): aquifer storage (m)
                  ! input (diagnostic variables)
                  mLayerRootDensity(1:nSoil),        & ! intent(in): root density in each layer (-)
                  scalarAquiferRootFrac,             & ! intent(in): fraction of roots below the lowest soil layer (-)
                  ! input (parameters)
                  plantWiltPsi,                      & ! intent(in): matric head at wilting point (m)
                  soilStressParam,                   & ! intent(in): parameter in the exponential soil stress function (-)
                  critSoilWilting,                   & ! intent(in): critical vol. liq. water content when plants are wilting (-)
                  critSoilTranspire,                 & ! intent(in): critical vol. liq. water content when transpiration is limited (-)
                  critAquiferTranspire,              & ! intent(in): critical aquifer storage value when transpiration is limited (m)
                  ! output
                  scalarTranspireLim,                & ! intent(out): weighted average of the transpiration limiting factor (-)
                  mLayerTranspireLim(1:nSoil),       & ! intent(out): transpiration limiting factor in each layer (-)
                  scalarTranspireLimAqfr,            & ! intent(out): transpiration limiting factor for the aquifer (-)
                  err,cmessage                       ) ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  !print*, 'weighted average of the soil moiture factor controlling stomatal resistance (-) = ', scalarTranspireLim

  ! compute stomatal resistance (wrapper sound the Noah-MP routines)
  ! NOTE: canopy temperature and canopy air vapor pressure are from the previous time step
  call stomResist(&
                  ! input (model decisions)
                  ix_stomResist,                     & ! intent(in): choice of function for stomatal resistance
                  ! input (local attributes)
                  vegTypeIndex,                      & ! intent(in): vegetation type index
                  iLoc, jLoc,                        & ! intent(in): spatial location indices      
                  ! input (forcing)
                  airtemp,                           & ! intent(in): air temperature at some height above the surface (K)
                  airpres,                           & ! intent(in): air pressure at some height above the surface (Pa)
                  scalarO2air,                       & ! intent(in): atmospheric o2 concentration (Pa)
                  scalarCO2air,                      & ! intent(in): atmospheric co2 concentration (Pa)
                  scalarCanopySunlitPAR,             & ! intent(in): average absorbed par for sunlit leaves (w m-2)
                  scalarCanopyShadedPAR,             & ! intent(in): average absorbed par for shaded leaves (w m-2)
                  ! input (state and diagnostic variables)
                  scalarGrowingSeasonIndex,          & ! intent(in): growing season index (0=off, 1=on)
                  scalarFoliageNitrogenFactor,       & ! intent(in): foliage nitrogen concentration (1=saturated)
                  scalarTranspireLim,                & ! intent(in): weighted average of the soil moiture factor controlling stomatal resistance (-)
                  scalarLeafResistance,              & ! intent(in): leaf boundary layer resistance (s m-1)
                  canopyTempTrial,                   & ! intent(in): temperature of the vegetation canopy (K)
                  scalarSatVP_CanopyTemp,            & ! intent(in): saturation vapor pressure at the temperature of the veg canopy (Pa)
                  scalarVP_CanopyAir,                & ! intent(in): canopy air vapor pressure (Pa)
                  ! output
                  scalarStomResistSunlit,            & ! intent(out): stomatal resistance for sunlit leaves (s m-1)
                  scalarStomResistShaded,            & ! intent(out): stomatal resistance for shaded leaves (s m-1)
                  scalarPhotosynthesisSunlit,        & ! intent(out): sunlit photosynthesis (umolco2 m-2 s-1)
                  scalarPhotosynthesisShaded,        & ! intent(out): shaded photosynthesis (umolco2 m-2 s-1)
                  err,message                                                    ) ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 endif  ! (if the first iteration in a given sub-step)


 ! *******************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************
 ! ***** LONGWAVE RADIATION  *****************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************

 ! compute canopy longwave radiation balance
 call longwaveBal(&
                  ! input: model control
                  ix_fDerivMeth,                     & ! intent(in): method used to calculate flux derivatives 
                  ! input: canopy and ground temperature
                  canopyTempTrial,                   & ! intent(in): temperature of the vegetation canopy (K)
                  groundTempTrial,                   & ! intent(in): temperature of the ground surface (K)
                  ! input: canopy and ground emissivity
                  canopyEmissivity,                  & ! intent(in): canopy emissivity (-)
                  groundEmissivity,                  & ! intent(in): ground emissivity (-)
                  ! input: forcing
                  LWRadAtm,                          & ! intent(in): downwelling longwave radiation at the upper boundary (W m-2)
                  ! output: emitted radiation from the canopy and ground
                  scalarLWRadCanopy,                 & ! intent(out): longwave radiation emitted from the canopy (W m-2)
                  scalarLWRadGround,                 & ! intent(out): longwave radiation emitted at the ground surface (W m-2)
                  ! output: individual fluxes
                  scalarLWRadUbound2Canopy,          & ! intent(out): downward atmospheric longwave radiation absorbed by the canopy (W m-2)
                  scalarLWRadUbound2Ground,          & ! intent(out): downward atmospheric longwave radiation absorbed by the ground (W m-2)
                  scalarLWRadUbound2Ubound,          & ! intent(out): atmospheric radiation reflected by the ground and lost thru upper boundary (W m-2)
                  scalarLWRadCanopy2Ubound,          & ! intent(out): longwave radiation emitted from canopy lost thru upper boundary (W m-2)
                  scalarLWRadCanopy2Ground,          & ! intent(out): longwave radiation emitted from canopy absorbed by the ground (W m-2)
                  scalarLWRadCanopy2Canopy,          & ! intent(out): canopy longwave reflected from ground and absorbed by the canopy (W m-2)
                  scalarLWRadGround2Ubound,          & ! intent(out): longwave radiation emitted from ground lost thru upper boundary (W m-2)
                  scalarLWRadGround2Canopy,          & ! intent(out): longwave radiation emitted from ground and absorbed by the canopy (W m-2)
                  ! output: net fluxes
                  scalarLWNetCanopy,                 & ! intent(out): net longwave radiation at the canopy (W m-2)
                  scalarLWNetGround,                 & ! intent(out): net longwave radiation at the ground surface (W m-2)
                  scalarLWNetUbound,                 & ! intent(out): net longwave radiation at the upper boundary (W m-2)
                  ! output: flux derivatives
                  dLWNetCanopy_dTCanopy,             & ! intent(out): derivative in net canopy radiation w.r.t. canopy temperature (W m-2 K-1)
                  dLWNetGround_dTGround,             & ! intent(out): derivative in net ground radiation w.r.t. ground temperature (W m-2 K-1)
                  dLWNetCanopy_dTGround,             & ! intent(out): derivative in net canopy radiation w.r.t. ground temperature (W m-2 K-1)
                  dLWNetGround_dTCanopy,             & ! intent(out): derivative in net ground radiation w.r.t. canopy temperature (W m-2 K-1)
                  ! output: error control
                  err,cmessage                       ) ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


 ! *******************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************
 ! ***** TURBULENT HEAT FLUXES  **************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************

 ! check the need to compute numerical derivatives
 if(ix_fDerivMeth == numerical)then
  nFlux=3  ! compute the derivatives using one-sided finite differences
 else
  nFlux=1  ! compute analytical derivatives
 endif

 ! either one or multiple flux calls, depending on if using analytical or numerical derivatives
 do itry=nFlux,1,-1  ! (work backwards to ensure all computed fluxes come from the un-perturbed case)

  ! -------------------------------------------------------------------------------------
  ! state perturbations for numerical deriavtives with one-sided finite differences
  ! note: no perturbations performed using analytical derivatives (nFlux=1)
  ! -------------------------------------------------------------------------------------

  ! identify the type of perturbation
  select case(itry)

   ! un-perturbed case
   case(unperturbed)
    canopyTemp = canopyTempTrial 
    groundTemp = groundTempTrial

   ! perturb canopy temperature
   case(perturbStateCanopy)
    canopyTemp = canopyTempTrial + dx
    groundTemp = groundTempTrial

   ! perturb ground temperature
   case(perturbStateGround)
    canopyTemp = canopyTempTrial
    groundTemp = groundTempTrial + dx

   ! check for an unknown perturbation 
   case default; err=10; message=trim(message)//"unknown perturbation"; return

  end select ! (type of perturbation)

  ! compute the saturation vapor pressure for vegetation temperature
  ! NOTE: saturated vapor pressure derivatives don't seem that accurate....
  TV_celcius = canopyTemp - Tfreeze
  call satVapPress(TV_celcius, scalarSatVP_CanopyTemp, dSVPCanopy_dCanopyTemp)

  ! compute the saturation vapor pressure for ground temperature
  ! NOTE: saturated vapor pressure derivatives don't seem that accurate....
  TG_celcius = groundTemp - Tfreeze
  call satVapPress(TG_celcius, scalarSatVP_GroundTemp, dSVPGround_dGroundTemp)


  ! -------------------------------------------------------------------------------------
  ! calculation block (unperturbed fluxes returned [computed last])
  ! -------------------------------------------------------------------------------------

  ! re-compute aerodynamic resistances for perturbed cases
  ! NOTE: unperturbed fluxes computed earlier, and not over-written
  if(itry /= unperturbed)then
   call aeroResist(&
                   ! input: model control
                   .false.,                             & ! intent(in): logical flag if would like to compute analytical derivaties
                   ix_astability,                       & ! intent(in): choice of stability function
                   ! input: above-canopy forcing data
                   mHeight,                             & ! intent(in): measurement height (m)
                   airtemp,                             & ! intent(in): air temperature at some height above the surface (K)
                   windspd,                             & ! intent(in): wind speed at some height above the surface (m s-1)
                   ! input: canopy and ground temperature
                   canopyTemp,                          & ! intent(in): canopy temperature (K)
                   groundTemp,                          & ! intent(in): ground temperature (K)
                   ! input: diagnostic variables
                   exposedVAI,                          & ! intent(in): exposed vegetation area index -- leaf plus stem (m2 m-2)
                   scalarSnowDepth,                     & ! intent(in): snow depth (m)
                   ! input: parameters
                   z0Ground,                            & ! intent(in): roughness length of the ground (below canopy or non-vegetated surface [snow]) (m)
                   z0Canopy,                            & ! intent(in): roughness length of the canopy (m)
                   critRichNumber,                      & ! intent(in): critical value for the bulk Richardson number where turbulence ceases (-)
                   Louis79_bparam,                      & ! intent(in): parameter in Louis (1979) stability function
                   Louis79_cStar,                       & ! intent(in): parameter in Louis (1979) stability function
                   Mahrt87_eScale,                      & ! intent(in): exponential scaling factor in the Mahrt (1987) stability function
                   windReductionParam,                  & ! intent(in): canopy wind reduction parameter (-)                   
                   leafExchangeCoeff,                   & ! intent(in): turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
                   leafDimension,                       & ! intent(in): characteristic leaf dimension (m)
                   canopyHeight,                        & ! intent(in): canopy height (m) 
                   ! output: scalar resistances
                   notUsed_WindReductionFactor,         & ! intent(out): canopy wind reduction factor (-)
                   notUsed_ZeroPlaneDisplacement,       & ! intent(out): zero plane displacement (m) 
                   notUsed_Sfc2AtmExchangeCoeff,        & ! intent(out): surface-atmosphere turbulent exchange coefficient (-)
                   notUsed_EddyDiffusCanopyTop,         & ! intent(out): eddy diffusivity for heat at the top of the canopy (m2 s-1)
                   notUsed_WindspdCanopyTop,            & ! intent(out): windspeed at the top of the canopy (m s-1)
                   trialLeafResistance,                 & ! intent(out): mean leaf boundary layer resistance per unit leaf area (s m-1)
                   trialGroundResistance,               & ! intent(out): below canopy aerodynamic resistance (s m-1) 
                   trialCanopyResistance,               & ! intent(out): above canopy aerodynamic resistance (s m-1)
                   ! output: derivatives in scalar resistances
                   notUsed_dGroundResistance_dTCanopy,  & ! intent(out): derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
                   notUsed_dCanopyResistance_dTCanopy,  & ! intent(out): derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
                   notUsed_dLeafResistance_dTCanopy,    & ! intent(out): derivative in leaf resistance w.r.t. canopy temperature (s m-1 K-1)
                   ! output: error control
                   err,cmessage                          ) ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! assign scalar resistances for un-perturbed cases
  else
   trialLeafResistance   = scalarLeafResistance
   trialGroundResistance = scalarGroundResistance
   trialCanopyResistance = scalarCanopyResistance

  endif  ! (re-computing resistances for perturbed cases)
  !print*, 'trialLeafResistance = ', trialLeafResistance
  !print*, 'trialGroundResistance = ', trialGroundResistance
  !print*, 'trialCanopyResistance = ', trialCanopyResistance

  ! compute the relative humidity in the top soil layer and the resistance at the ground surface
  ! (soil water evaporation factor [0-1])
  soilEvapFactor = mLayerVolFracLiq(1)/theta_sat
  ! (resistance from the soil [s m-1])
  scalarSoilResistance = fracSnow*1._dp + (1._dp - fracSnow)*exp(8.25_dp - 6.0_dp*soilEvapFactor)
  ! (relative humidity in the soil pores [0-1])
  soilRelHumidity_noSnow = exp( (mLayerMatricHead(1)*gravity) / (groundTemp*R_wv) )
  scalarSoilRelHumidity  = fracSnow*1._dp + (1._dp - fracSnow)*soilRelHumidity_noSnow

  ! compute turbulent heat fluxes
  call turbFluxes(&
                  ! input: model control
                  dt,                                   & ! intent(in): model time step (seconds)
                  ix_fDerivMeth,                        & ! intent(in): method used to calculate flux derivatives
                  ! input: above-canopy forcing data
                  airtemp,                              & ! intent(in): air temperature at some height above the surface (K)
                  scalarVPair,                          & ! intent(in): vapor pressure of the air above the vegetation canopy (Pa)
                  ! input: psychometric constant for canopy and ground
                  scalarPsycConstCanopy,                & ! intent(in): psychrometric constant for the vegetation canopy (Pa/K)
                  scalarPsycConstGround,                & ! intent(in): psychrometric constant for the ground surface (Pa/K)
                  ! input: canopy liquid and canopy ice (used as a solution constraint)
                  canopyLiqTrial,                       & ! intent(in): mass of liquid water on the vegetation canopy (kg m-2)
                  canopyIceTrial,                       & ! intent(in): mass of ice water on the vegetation canopy (kg m-2)
                  ! input: canopy/ground temperature and saturated vapor pressure
                  canopyTemp,                           & ! intent(in): canopy temperature (K)
                  groundTemp,                           & ! intent(in): ground temperature (K)
                  scalarSatVP_CanopyTemp,               & ! intent(in): saturation vapor pressure at the temperature of the veg canopy (Pa)
                  scalarSatVP_GroundTemp,               & ! intent(in): saturation vapor pressure at the temperature of the ground (Pa)
                  dSVPCanopy_dCanopyTemp,               & ! intent(in): derivative in canopy saturation vapor pressure w.r.t. canopy temperature (Pa K-1)
                  dSVPGround_dGroundTemp,               & ! intent(in): derivative in ground saturation vapor pressure w.r.t. ground temperature (Pa K-1)
                  ! input: diagnostic variables
                  exposedVAI,                           & ! intent(in): exposed vegetation area index -- leaf plus stem (m2 m-2)
                  scalarCanopyWetFraction,              & ! intent(in): fraction of canopy that is wet [0-1]
                  scalarCanopySunlitLAI,                & ! intent(in): sunlit leaf area (-)
                  scalarCanopyShadedLAI,                & ! intent(in): shaded leaf area (-)
                  scalarSoilRelHumidity,                & ! intent(in): relative humidity in the soil pores [0-1]
                  scalarSoilResistance,                 & ! intent(in): resistance from the soil (s m-1)
                  trialLeafResistance,                  & ! intent(in): mean leaf boundary layer resistance per unit leaf area (s m-1)
                  trialGroundResistance,                & ! intent(in): below canopy aerodynamic resistance (s m-1) 
                  trialCanopyResistance,                & ! intent(in): above canopy aerodynamic resistance (s m-1)
                  scalarStomResistSunlit,               & ! intent(in): stomatal resistance for sunlit leaves (s m-1)
                  scalarStomResistShaded,               & ! intent(in): stomatal resistance for shaded leaves (s m-1)
                  ! input: derivatives in scalar resistances
                  dGroundResistance_dTCanopy,           & ! intent(in): derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
                  dCanopyResistance_dTCanopy,           & ! intent(in): derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
                  dLeafResistance_dTCanopy,             & ! intent(in): derivative in leaf resistance w.r.t. canopy temperature (s m-1 K-1)
                  ! output: canopy air space variables
                  scalarTemp_CanopyAir,                 & ! intent(out): temperature of the canopy air space (K)
                  scalarVP_CanopyAir,                   & ! intent(out): vapor pressure of the canopy air space (Pa)
                  ! output: fluxes from the vegetation canopy
                  scalarSenHeatCanopy,                  & ! intent(out): sensible heat flux from the canopy to the canopy air space (W m-2)
                  scalarLatHeatCanopyEvap,              & ! intent(out): latent heat flux associated with evaporation from the canopy to the canopy air space (W m-2)
                  scalarLatHeatCanopyTrans,             & ! intent(out): latent heat flux associated with transpiration from the canopy to the canopy air space (W m-2)
                  ! output: fluxes from non-vegetated surfaces (ground surface below vegetation, bare ground, or snow covered vegetation)
                  scalarSenHeatGround,                  & ! intent(out): sensible heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
                  scalarLatHeatGround,                  & ! intent(out): latent heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
                  ! output: net fluxes
                  turbFluxCanopy,                       & ! intent(out): net longwave radiation at the canopy (W m-2)
                  turbFluxGround,                       & ! intent(out): net longwave radiation at the ground surface (W m-2)
                  ! output: flux derivatives
                  dTurbFluxCanopy_dTCanopy,             & ! intent(out): derivative in net canopy turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
                  dTurbFluxGround_dTCanopy,             & ! intent(out): derivative in net ground turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
                  dTurbFluxCanopy_dTGround,             & ! intent(out): derivative in net canopy turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
                  dTurbFluxGround_dTGround,             & ! intent(out): derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
                  ! output: error control
                  err,cmessage                          ) ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! save perturbed fluxes
  if(ix_fDerivMeth == numerical)then
   select case(itry) ! (select type of perturbation)
    case(unperturbed)
     exit
    case(perturbStateCanopy)
     turbFluxCanopy_dStateCanopy = turbFluxCanopy         ! total turbulent heat fluxes from the canopy to the canopy air space (W m-2)
     turbFluxGround_dStateCanopy = turbFluxGround         ! total turbulent heat fluxes from the ground to the canopy air space (W m-2)
    case(perturbStateGround)
     turbFluxCanopy_dStateGround = turbFluxCanopy         ! total turbulent heat fluxes from the canopy to the canopy air space (W m-2)
     turbFluxGround_dStateGround = turbFluxGround         ! total turbulent heat fluxes from the ground to the canopy air space (W m-2)
    case default; err=10; message=trim(message)//"unknown perturbation"; return
   end select ! (type of perturbation)
  endif ! (if numerical)

 end do  ! (looping through different flux perturbations)

 ! compute numerical derivatives
 if(ix_fDerivMeth == numerical)then
  ! derivatives w.r.t. canopy temperature
  dTurbFluxCanopy_dTCanopy = (turbFluxCanopy_dStateCanopy - turbFluxCanopy) / dx  ! derivative in net canopy turbulent fluxes w.r.t. canopy temperature (W m-2 K-1) 
  dTurbFluxGround_dTCanopy = (turbFluxGround_dStateCanopy - turbFluxGround) / dx  ! derivative in net ground turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
  ! derivatives w.r.t. ground temperature
  dTurbFluxCanopy_dTGround = (turbFluxCanopy_dStateGround - turbFluxCanopy) / dx  ! derivative in net canopy turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
  dTurbFluxGround_dTGround = (turbFluxGround_dStateGround - turbFluxGround) / dx  ! derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
 endif

 ! test
 print*, (ix_fDerivMeth == numerical)
 print*, 'dTurbFluxCanopy_dTCanopy = ', dTurbFluxCanopy_dTCanopy
 print*, 'dTurbFluxGround_dTCanopy = ', dTurbFluxGround_dTCanopy
 print*, 'dTurbFluxCanopy_dTGround = ', dTurbFluxCanopy_dTGround
 print*, 'dTurbFluxGround_dTGround = ', dTurbFluxGround_dTGround


 ! *******************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************
 ! ***** AND STITCH EVERYTHING TOGETHER  *****************************************************************************************************************************
 ! *******************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************

 ! compute net fluxes at the canopy and ground surface
 canopyNetFlux = scalarCanopyAbsorbedSolar + scalarLWNetCanopy + turbFluxCanopy
 groundNetFlux = scalarGroundAbsorbedSolar + scalarLWNetGround + turbFluxGround

 ! compute the derivatives
 dCanopyNetFlux_dCanopyTemp = dLWNetCanopy_dTCanopy + dTurbFluxCanopy_dTCanopy
 dGroundNetFlux_dCanopyTemp = dLWNetGround_dTCanopy + dTurbFluxGround_dTCanopy 
 dCanopyNetFlux_dGroundTemp = dLWNetCanopy_dTGround + dTurbFluxCanopy_dTGround
 dGroundNetFlux_dGroundTemp = dLWNetGround_dTGround + dTurbFluxGround_dTGround

 end subroutine vegNrgFlux_muster


 ! *********************************************************************************************************************************************************************
 ! *********************************************************************************************************************************************************************
 ! *********************************************************************************************************************************************************************
 ! ***** PRIVATE SUBROUTINES *******************************************************************************************************************************************
 ! *********************************************************************************************************************************************************************
 ! *********************************************************************************************************************************************************************
 ! *********************************************************************************************************************************************************************


 ! ***************************************************************************************
 ! private subroutine: compute longwave radiation balance at the canopy and ground surface
 ! ***************************************************************************************
 subroutine longwaveBal(&
                        ! input: model control
                        ixDerivMethod,                  & ! intent(in): choice of method used to compute derivative (analytical or numerical)
                        ! input: canopy and ground temperature
                        canopyTemp,                     & ! intent(in): canopy temperature (K)
                        groundTemp,                     & ! intent(in): ground temperature (K)
                        ! input: canopy and ground emissivity
                        emc,                            & ! intent(in): canopy emissivity (-)
                        emg,                            & ! intent(in): ground emissivity (-)
                        ! input: forcing
                        LWRadUbound,                    & ! intent(in): downwelling longwave radiation at the upper boundary (W m-2)
                        ! output: sources
                        LWRadCanopy,                    & ! intent(out): longwave radiation emitted from the canopy (W m-2)
                        LWRadGround,                    & ! intent(out): longwave radiation emitted at the ground surface (W m-2)
                        ! output: individual fluxes
                        LWRadUbound2Canopy,             & ! intent(out): downward atmospheric longwave radiation absorbed by the canopy (W m-2)
                        LWRadUbound2Ground,             & ! intent(out): downward atmospheric longwave radiation absorbed by the ground (W m-2)
                        LWRadUbound2Ubound,             & ! intent(out): atmospheric radiation reflected by the ground and lost thru upper boundary (W m-2)
                        LWRadCanopy2Ubound,             & ! intent(out): longwave radiation emitted from canopy lost thru upper boundary (W m-2)
                        LWRadCanopy2Ground,             & ! intent(out): longwave radiation emitted from canopy absorbed by the ground (W m-2)
                        LWRadCanopy2Canopy,             & ! intent(out): canopy longwave reflected from ground and absorbed by the canopy (W m-2)
                        LWRadGround2Ubound,             & ! intent(out): longwave radiation emitted from ground lost thru upper boundary (W m-2)
                        LWRadGround2Canopy,             & ! intent(out): longwave radiation emitted from ground and absorbed by the canopy (W m-2)
                        ! output: net fluxes
                        LWNetCanopy,                    & ! intent(out): net longwave radiation at the canopy (W m-2)
                        LWNetGround,                    & ! intent(out): net longwave radiation at the ground surface (W m-2)
                        LWNetUbound,                    & ! intent(out): net longwave radiation at the upper boundary (W m-2)
                        ! output: flux derivatives
                        dLWNetCanopy_dTCanopy,          & ! intent(out): derivative in net canopy radiation w.r.t. canopy temperature (W m-2 K-1)
                        dLWNetGround_dTGround,          & ! intent(out): derivative in net ground radiation w.r.t. ground temperature (W m-2 K-1)
                        dLWNetCanopy_dTGround,          & ! intent(out): derivative in net canopy radiation w.r.t. ground temperature (W m-2 K-1)
                        dLWNetGround_dTCanopy,          & ! intent(out): derivative in net ground radiation w.r.t. canopy temperature (W m-2 K-1)
                        ! output: error control
                        err,message                     ) ! intent(out): error control
 ! -----------------------------------------------------------------------------------------------------------------------------------------------
 implicit none
 ! input: model control
 integer(i4b),intent(in)       :: ixDerivMethod            ! choice of method used to compute derivative (analytical or numerical)
 ! input: canopy and ground temperature
 real(dp),intent(in)           :: canopyTemp               ! canopy temperature (K)
 real(dp),intent(in)           :: groundTemp               ! ground temperature (K)
 ! input: canopy and ground emissivity
 real(dp),intent(in)           :: emc                      ! canopy emissivity (-)
 real(dp),intent(in)           :: emg                      ! ground emissivity (-)
 ! input: forcing
 real(dp),intent(in)           :: LWRadUbound              ! downwelling longwave radiation at the upper boundary (W m-2)
 ! output: sources
 real(dp),intent(out)          :: LWRadCanopy              ! longwave radiation emitted from the canopy (W m-2)
 real(dp),intent(out)          :: LWRadGround              ! longwave radiation emitted at the ground surface (W m-2)
 ! output: individual fluxes
 real(dp),intent(out)          :: LWRadUbound2Canopy       ! downward atmospheric longwave radiation absorbed by the canopy (W m-2)
 real(dp),intent(out)          :: LWRadUbound2Ground       ! downward atmospheric longwave radiation absorbed by the ground (W m-2)
 real(dp),intent(out)          :: LWRadUbound2Ubound       ! atmospheric radiation reflected by the ground and lost thru upper boundary (W m-2)
 real(dp),intent(out)          :: LWRadCanopy2Ubound       ! longwave radiation emitted from canopy lost thru upper boundary (W m-2)
 real(dp),intent(out)          :: LWRadCanopy2Ground       ! longwave radiation emitted from canopy absorbed by the ground (W m-2)
 real(dp),intent(out)          :: LWRadCanopy2Canopy       ! canopy longwave reflected from ground and absorbed by the canopy (W m-2)
 real(dp),intent(out)          :: LWRadGround2Ubound       ! longwave radiation emitted from ground lost thru upper boundary (W m-2)
 real(dp),intent(out)          :: LWRadGround2Canopy       ! longwave radiation emitted from ground and absorbed by the canopy (W m-2)
 ! output: net fluxes
 real(dp),intent(out)          :: LWNetCanopy              ! net longwave radiation at the canopy (W m-2)
 real(dp),intent(out)          :: LWNetGround              ! net longwave radiation at the ground surface (W m-2)
 real(dp),intent(out)          :: LWNetUbound              ! net longwave radiation at the upper boundary (W m-2)
 ! output: flux derivatives
 real(dp),intent(out)          :: dLWNetCanopy_dTCanopy    ! derivative in net canopy radiation w.r.t. canopy temperature (W m-2 K-1)
 real(dp),intent(out)          :: dLWNetGround_dTGround    ! derivative in net ground radiation w.r.t. ground temperature (W m-2 K-1)
 real(dp),intent(out)          :: dLWNetCanopy_dTGround    ! derivative in net canopy radiation w.r.t. ground temperature (W m-2 K-1)
 real(dp),intent(out)          :: dLWNetGround_dTCanopy    ! derivative in net ground radiation w.r.t. canopy temperature (W m-2 K-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 integer(i4b),parameter        :: unperturbed=1            ! named variable to identify the case of unperturbed state variables
 integer(i4b),parameter        :: perturbStateCanopy=2     ! named variable to identify the case where we perturb the canopy temperature
 integer(i4b),parameter        :: perturbStateGround=3     ! named variable to identify the case where we perturb the ground temperature
 integer(i4b)                  :: itry                     ! index of flux evaluation
 integer(i4b)                  :: nFlux                    ! number of flux evaluations
 real(dp)                      :: TCan                     ! value of canopy temperature used in flux calculations (may be perturbed)
 real(dp)                      :: TGnd                     ! value of ground temperature used in flux calculations (may be perturbed)
 real(dp)                      :: fluxBalance              ! check energy closure (W m-2)
 real(dp),parameter            :: fluxTolerance=1.e-12_dp  ! tolerance for energy closure (W m-2)
 real(dp)                      :: dLWRadCanopy_dTCanopy    ! derivative in emitted radiation at the canopy w.r.t. canopy temperature
 real(dp)                      :: dLWRadGround_dTGround    ! derivative in emitted radiation at the ground w.r.t. ground temperature
 real(dp)                      :: LWNetCanopy_dStateCanopy ! net lw canopy flux after perturbation in canopy temperature
 real(dp)                      :: LWNetGround_dStateCanopy ! net lw ground flux after perturbation in canopy temperature
 real(dp)                      :: LWNetCanopy_dStateGround ! net lw canopy flux after perturbation in ground temperature
 real(dp)                      :: LWNetGround_dStateGround ! net lw ground flux after perturbation in ground temperature
 ! -----------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='longwaveBal/'

 ! check the need to compute numerical derivatives
 if(ixDerivMethod==numerical)then
  nFlux=3  ! compute the derivatives using one-sided finite differences
 else
  nFlux=1  ! compute analytical derivatives
 endif

 ! either one or multiple flux calls, depending on if using analytical or numerical derivatives
 do itry=nFlux,1,-1  ! (work backwards to ensure all computed fluxes come from the un-perturbed case)

  ! -------------------------------------------------------------------------------------
  ! state perturbations for numerical deriavtives with one-sided finite differences
  ! note: no perturbations performed using analytical derivatives (nFlux=1)
  ! -------------------------------------------------------------------------------------

  ! identify the type of perturbation
  select case(itry)

   ! un-perturbed case
   case(unperturbed)
    TCan = canopyTemp
    TGnd = groundTemp

   ! perturb canopy temperature
   case(perturbStateCanopy)
    TCan = canopyTemp + dx 
    TGnd = groundTemp

   ! perturb ground temperature
   case(perturbStateGround)
    TCan = canopyTemp
    TGnd = groundTemp + dx

   ! check for an unknown perturbation 
   case default; err=10; message=trim(message)//"unknown perturbation"; return

  end select ! (type of perturbation)

  ! -------------------------------------------------------------------------------------
  ! calculation block (unperturbed fluxes returned [computed last])
  ! -------------------------------------------------------------------------------------

  ! compute longwave fluxes from canopy and the ground
  LWRadCanopy = emc*sb*TCan**4._dp                                           ! longwave radiation emitted from the canopy (W m-2)
  LWRadGround = emg*sb*TGnd**4._dp                                           ! longwave radiation emitted at the ground surface (W m-2)
  
  ! compute fluxes originating from the atmosphere
  LWRadUbound2Canopy = (emc + (1._dp - emc)*(1._dp - emg)*emc)*LWRadUbound   ! downward atmospheric longwave radiation absorbed by the canopy (W m-2) 
  LWRadUbound2Ground = (1._dp - emc)*emg*LWRadUbound                         ! downward atmospheric longwave radiation absorbed by the ground (W m-2)
  LWRadUbound2Ubound = (1._dp - emc)*(1._dp - emg)*(1._dp - emc)*LWRadUbound ! atmospheric radiation reflected by the ground and lost thru upper boundary (W m-2)

  ! compute fluxes originating from the canopy
  LWRadCanopy2Ubound = (1._dp + (1._dp - emc)*(1._dp - emg))*LWRadCanopy     ! longwave radiation emitted from canopy lost thru upper boundary (W m-2) 
  LWRadCanopy2Ground = emg*LWRadCanopy                                       ! longwave radiation emitted from canopy absorbed by the ground (W m-2)
  LWRadCanopy2Canopy = emc*(1._dp - emg)*LWRadCanopy                         ! canopy longwave reflected from ground and absorbed by the canopy (W m-2)

  ! compute fluxes originating from the ground surface
  LWRadGround2Ubound = (1._dp - emc)*LWRadGround                             ! longwave radiation emitted from ground lost thru upper boundary (W m-2) 
  LWRadGround2Canopy = emc*LWRadGround                                       ! longwave radiation emitted from ground and absorbed by the canopy (W m-2)

  ! compute net longwave radiation (W m-2)
  LWNetCanopy = LWRadUbound2Canopy + LWRadGround2Canopy + LWRadCanopy2Canopy - 2._dp*LWRadCanopy  ! canopy
  LWNetGround = LWRadUbound2Ground + LWRadCanopy2Ground - LWRadGround                             ! ground surface
  LWNetUbound = LWRadUbound - LWRadUbound2Ubound - LWRadCanopy2Ubound - LWRadGround2Ubound                             ! upper boundary

  ! check the flux balance
  fluxBalance = LWNetUbound - (LWNetCanopy + LWNetGround)
  if(abs(fluxBalance) > fluxTolerance)then
   print*, 'fluxBalance = ', fluxBalance
   print*, 'emg, emc = ', emg, emc
   print*, 'LWRadUbound = ', LWRadUbound
   print*, 'LWRadCanopy = ', LWRadCanopy
   print*, 'LWRadGround = ', LWRadGround
   print*, 'LWRadUbound2Canopy = ', LWRadUbound2Canopy
   print*, 'LWRadUbound2Ground = ', LWRadUbound2Ground
   print*, 'LWRadUbound2Ubound = ', LWRadUbound2Ubound
   print*, 'LWRadCanopy2Ubound = ', LWRadCanopy2Ubound
   print*, 'LWRadCanopy2Ground = ', LWRadCanopy2Ground
   print*, 'LWRadCanopy2Canopy = ', LWRadCanopy2Canopy
   print*, 'LWRadGround2Ubound = ', LWRadGround2Ubound
   print*, 'LWRadGround2Canopy = ', LWRadGround2Canopy
   print*, 'LWNetCanopy = ', LWNetCanopy
   print*, 'LWNetGround = ', LWNetGround
   print*, 'LWNetUbound = ', LWNetUbound
   message=trim(message)//'flux imbalance'
   err=20; return
  endif

  ! --------------------------------------------------------------------------------------
  ! save perturbed fluxes to calculate numerical derivatives (one-sided finite difference)
  ! --------------------------------------------------------------------------------------
  if(ixDerivMethod==numerical)then
   select case(itry) ! (select type of perturbation)
    case(unperturbed); exit
    case(perturbStateCanopy)
     LWNetCanopy_dStateCanopy = LWNetCanopy
     LWNetGround_dStateCanopy = LWNetGround
    case(perturbStateGround)
     LWNetCanopy_dStateGround = LWNetCanopy
     LWNetGround_dStateGround = LWNetGround
    case default; err=10; message=trim(message)//"unknown perturbation"; return
   end select ! (type of perturbation)
  endif ! (if numerical)

 end do  ! looping through different perturbations

 ! -------------------------------------------------------------------------------------
 ! compute derivatives
 ! -------------------------------------------------------------------------------------
 select case(ixDerivMethod)

  ! ***** analytical derivatives
  case(analytical)
   ! compute initial derivatives
   dLWRadCanopy_dTCanopy = 4._dp*emc*sb*TCan**3._dp
   dLWRadGround_dTGround = 4._dp*emg*sb*TGnd**3._dp
   ! compute analytical derivatives
   dLWNetCanopy_dTCanopy = (emc*(1._dp - emg) - 2._dp)*dLWRadCanopy_dTCanopy ! derivative in net canopy radiation w.r.t. canopy temperature (W m-2 K-1)
   dLWNetGround_dTGround = dLWRadGround_dTGround      ! derivative in net ground radiation w.r.t. ground temperature (W m-2 K-1)
   dLWNetCanopy_dTGround = emc*dLWRadGround_dTGround  ! derivative in net canopy radiation w.r.t. ground temperature (W m-2 K-1)
   dLWNetGround_dTCanopy = emg*dLWRadCanopy_dTCanopy  ! derivative in net ground radiation w.r.t. canopy temperature (W m-2 K-1)

  ! ***** numerical derivatives
  case(numerical)
   ! compute numerical derivatives (one-sided finite differences)
   dLWNetCanopy_dTCanopy = (LWNetCanopy_dStateCanopy - LWNetCanopy)/dx  ! derivative in net canopy radiation w.r.t. canopy temperature (W m-2 K-1)
   dLWNetGround_dTGround = (LWNetGround_dStateGround - LWNetGround)/dx  ! derivative in net ground radiation w.r.t. ground temperature (W m-2 K-1)
   dLWNetCanopy_dTGround = (LWNetCanopy_dStateGround - LWNetCanopy)/dx  ! derivative in net canopy radiation w.r.t. ground temperature (W m-2 K-1)
   dLWNetGround_dTCanopy = (LWNetGround_dStateCanopy - LWNetGround)/dx  ! derivative in net ground radiation w.r.t. canopy temperature (W m-2 K-1)

  ! ***** error check
  case default; err=10; message=trim(message)//"unknown method to calculate derivatives"; return

 end select ! (type of method to calculate derivatives)

 end subroutine longwaveBal


 ! ********************************************************************************
 ! private subroutine: compute aerodynamic resistances
 ! ********************************************************************************
 subroutine aeroResist(&
                       ! input: model control
                       derivDesired,                  & ! intent(in): flag to indicate if derivatives are desired
                       ixStability,                   & ! intent(in): choice of stability function
                       ! input: above-canopy forcing data
                       mHeight,                       & ! intent(in): measurement height (m)
                       airtemp,                       & ! intent(in): air temperature at some height above the surface (K)
                       windspd,                       & ! intent(in): wind speed at some height above the surface (m s-1)
                       ! input: canopy and ground temperature
                       canopyTemp,                    & ! intent(in): canopy temperature (K)
                       groundTemp,                    & ! intent(in): ground temperature (K)
                       ! input: diagnostic variables
                       exposedVAI,                    & ! intent(in): exposed vegetation area index -- leaf plus stem (m2 m-2)
                       snowDepth,                     & ! intent(in): snow depth (m)
                       ! input: parameters
                       z0Ground,                      & ! intent(in): roughness length of the ground (below canopy or non-vegetated surface [snow]) (m)
                       z0Canopy,                      & ! intent(in): roughness length of the canopy (m)
                       critRichNumber,                & ! intent(in): critical value for the bulk Richardson number where turbulence ceases (-)
                       Louis79_bparam,                & ! intent(in): parameter in Louis (1979) stability function
                       Louis79_cStar,                 & ! intent(in): parameter in Louis (1979) stability function
                       Mahrt87_eScale,                & ! intent(in): exponential scaling factor in the Mahrt (1987) stability function
                       windReductionParam,            & ! intent(in): canopy wind reduction parameter (-)                   
                       leafExchangeCoeff,             & ! intent(in): turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
                       leafDimension,                 & ! intent(in): characteristic leaf dimension (m)
                       canopyHeight,                  & ! intent(in): canopy height (m) 
                       ! output: scalar resistances
                       windReductionFactor,           & ! intent(out): canopy wind reduction factor (-)
                       zeroPlaneDisplacement,         & ! intent(out): zero plane displacement (m) 
                       sfc2AtmExchangeCoeff,          & ! intent(out): surface-atmosphere turbulent exchange coefficient (-)
                       eddyDiffusCanopyTop,           & ! intent(out): eddy diffusivity for heat at the top of the canopy (m2 s-1)
                       windspdCanopyTop,              & ! intent(out): windspeed at the top of the canopy (m s-1)
                       leafResistance,                & ! intent(out): mean leaf boundary layer resistance per unit leaf area (s m-1)
                       groundResistance,              & ! intent(out): below canopy aerodynamic resistance (s m-1) 
                       canopyResistance,              & ! intent(out): above canopy aerodynamic resistance (s m-1)
                       ! output: derivatives in scalar resistances
                       dGroundResistance_dTCanopy,    & ! intent(out): derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
                       dCanopyResistance_dTCanopy,    & ! intent(out): derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
                       dLeafResistance_dTCanopy,      & ! intent(out): derivative in leaf resistance w.r.t. canopy temperature (s m-1 K-1)
                       ! output: error control
                       err,message                    ) ! intent(out): error control
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! compute aerodynamic resistances
 ! Refs: Choudhury and Monteith (4-layer model for heat budget of homogenous surfaces; QJRMS, 1988)
 !       Niu and Yang (Canopy effects on snow processes; JGR, 2004)
 !       Mahat et al. (Below-canopy turbulence in a snowmelt model, WRR, 2012)
 implicit none
 ! input: model control
 logical(lgt),intent(in)       :: derivDesired             ! flag to indicate if derivatives are desired
 integer(i4b),intent(in)       :: ixStability              ! choice of stability function
 ! input: above-canopy forcing data
 real(dp),intent(in)           :: mHeight                  ! measurement height (m)
 real(dp),intent(in)           :: airtemp                  ! air temperature at some height above the surface (K)
 real(dp),intent(in)           :: windspd                  ! wind speed at some height above the surface (m s-1)
 ! input: canopy and ground temperature
 real(dp),intent(in)           :: canopyTemp               ! canopy temperature (K)
 real(dp),intent(in)           :: groundTemp               ! ground temperature (K)
 ! input: diagnostic variables
 real(dp),intent(in)           :: exposedVAI               ! exposed vegetation area index -- leaf plus stem (m2 m-2)
 real(dp),intent(in)           :: snowDepth                ! snow depth (m)
 ! input: parameters
 real(dp),intent(in)           :: z0Ground                 ! roughness length of the ground (below canopy or non-vegetated surface [snow]) (m)
 real(dp),intent(in)           :: z0Canopy                 ! roughness length of the canopy (m)
 real(dp),intent(in)           :: critRichNumber           ! critical value for the bulk Richardson number where turbulence ceases (-)
 real(dp),intent(in)           :: Louis79_bparam           ! parameter in Louis (1979) stability function
 real(dp),intent(in)           :: Louis79_cStar            ! parameter in Louis (1979) stability function
 real(dp),intent(in)           :: Mahrt87_eScale           ! exponential scaling factor in the Mahrt (1987) stability function
 real(dp),intent(in)           :: windReductionParam       ! canopy wind reduction parameter (-)                   
 real(dp),intent(in)           :: leafExchangeCoeff        ! turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
 real(dp),intent(in)           :: leafDimension            ! characteristic leaf dimension (m)
 real(dp),intent(in)           :: canopyHeight             ! canopy height (m) 
 ! output: scalar resistances
 real(dp),intent(out)          :: windReductionFactor      ! canopy wind reduction factor (-)
 real(dp),intent(out)          :: zeroPlaneDisplacement    ! zero plane displacement (m) 
 real(dp),intent(out)          :: sfc2AtmExchangeCoeff     ! surface-atmosphere turbulent exchange coefficient (-)
 real(dp),intent(out)          :: eddyDiffusCanopyTop      ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
 real(dp),intent(out)          :: windspdCanopyTop         ! windspeed at the top of the canopy (m s-1)
 real(dp),intent(out)          :: leafResistance           ! mean leaf boundary layer resistance per unit leaf area (s m-1)
 real(dp),intent(out)          :: groundResistance         ! below canopy aerodynamic resistance (s m-1) 
 real(dp),intent(out)          :: canopyResistance         ! above canopy aerodynamic resistance (s m-1)
 ! output: derivatives in scalar resistances
 real(dp),intent(out)          :: dGroundResistance_dTCanopy    ! derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
 real(dp),intent(out)          :: dCanopyResistance_dTCanopy    ! derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
 real(dp),intent(out)          :: dLeafResistance_dTCanopy      ! derivative in leaf resistance w.r.t. canopy temperature (s m-1 K-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 character(LEN=256)            :: cmessage                 ! error message of downwind routine
 real(dp),parameter            :: oneThird=1._dp/3._dp     ! 1/3
 real(dp),parameter            :: twoThirds=2._dp/3._dp    ! 2/3
 real(dp)                      :: funcLAI                  ! temporary variable to calculate zero plane displacement for the canopy
 real(dp)                      :: fracCanopyHeight         ! zero plane displacement expressed as a fraction of canopy height
 real(dp)                      :: z0                       ! roughness length (m)
 real(dp)                      :: sfcTemp                  ! surface temperature -- "surface" is either canopy or ground (K)
 real(dp)                      :: windConvFactor           ! factor to convert friction velocity to wind speed at top of canopy (-)
 real(dp)                      :: dFV_dT                   ! derivative in friction velocity w.r.t. canopy temperature
 real(dp)                      :: dUC_dT                   ! derivative in wind speed at canopy top w.r.t. canopy temperature
 real(dp)                      :: dLC_dT                   ! derivative in canopy-average leaf conductance w.r.t. canopy temperature
 real(dp)                      :: dED_dT                   ! derivative in eddy diffusivity at the top of the canopy w.r.t. canopy temperature
 real(dp)                      :: tmp1,tmp2                ! temporary variables used in calculation of ground resistance
 real(dp)                      :: dSfc2AtmExCoef_dTemp     ! derivative in surface-atmosphere exchange coeffcient w.r.t. temperature (K-1)
 real(dp)                      :: frictionVelocity         ! friction velocity (m s-1)
 real(dp)                      :: singleLeafConductance    ! leaf boundary layer conductance (m s-1) 
 real(dp)                      :: canopyLeafConductance    ! leaf boundary layer conductance -- scaled up to the canopy (m s-1)
 real(dp)                      :: leaf2CanopyScaleFactor   ! factor to scale from the leaf to the canopy [m s-(1/2)]
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='aeroResist/'

 ! check that measurement height is above the top of the canopy
 if(mHeight < canopyHeight)then
  err=20; message=trim(message)//'measurement height is below the top of the canopy'; return
 endif

 ! ***** identify zero plane displacement, roughness length, and surface temperature for the canopy (m)
 if(exposedVAI > 0._dp) then ! (if vegetation is exposed)
  ! compute the zero plane displacement for the canopy (m)
  funcLAI          = sqrt(7.5_dp*exposedVAI)
  fracCanopyHeight = -(1._dp - exp(-funcLAI))/funcLAI + 1._dp
  zeroPlaneDisplacement = fracCanopyHeight*canopyHeight
  if(zeroPlaneDisplacement < snowDepth) zeroPlaneDisplacement = snowDepth
  ! assign roughness length to the canopy roughness (m)
  z0 = z0Canopy
  ! assign surface temperature to the canopy temperature (K)
  sfcTemp = canopyTemp

 ! ***** identify zero plane displacement, roughness length, and surface temperature for non-vegetated surfaces
 else
  zeroPlaneDisplacement = snowDepth                    ! zero plane displacement (m)
  z0 = z0Ground                                        ! roughness length (m)
  sfcTemp = groundTemp                                 ! surface temperature (K)
 endif

 ! compute the surface-atmosphere turbulent exchange coefficient (-)
 call turbExCoef(&
                 ! input
                 derivDesired,                   & ! input: logical flag to compute analytical derivatives
                 ixStability,                    & ! input: choice of stability function
                 ! input: forcing data, diagnostic and state variables
                 mHeight,                        & ! input: measurement height (m)
                 zeroPlaneDisplacement,          & ! input: zero plane displacement (m)
                 z0,                             & ! input: roughness length (m)
                 airTemp,                        & ! input: air temperature (K)
                 sfcTemp,                        & ! input: trial value of surface temperature -- "surface" is either canopy or ground (K)
                 windspd,                        & ! input: wind speed (m s-1)
                 ! input: stability parameters
                 critRichNumber,                 & ! input: critical value for the bulk Richardson number where turbulence ceases (-)
                 Louis79_bparam,                 & ! input: parameter in Louis (1979) stability function
                 Louis79_cStar,                  & ! input: parameter in Louis (1979) stability function
                 Mahrt87_eScale,                 & ! input: exponential scaling factor in the Mahrt (1987) stability function
                 ! output
                 sfc2AtmExchangeCoeff,           & ! output: surface-atmosphere turbulent exchange coefficient (-)
                 dSfc2AtmExCoef_dTemp,           & ! output: derivative in surface-atmosphere exchange coeffcient w.r.t. temperature (K-1)
                 err, cmessage                   ) ! output: error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! ***** compute the above-canopy, below-canopy, and leaf resistances (s m-1)
 if(exposedVAI > 0._dp) then ! (if vegetation is exposed)

  ! compute the friction velocity (m s-1)
  frictionVelocity = windspd * sqrt(sfc2AtmExchangeCoeff)

  ! compute the above-canopy resistance (s m-1)
  canopyResistance = 1._dp/(sfc2AtmExchangeCoeff*windspd)

  ! compute windspeed at the top of the canopy (m s-1)
  windConvFactor   = log((canopyHeight - zeroPlaneDisplacement)/z0Canopy)/vkc
  windspdCanopyTop = frictionVelocity*windConvFactor

  ! compute the windspeed reduction
  ! Refs: Norman et al. (Ag. Forest Met., 1995) -- citing Goudriaan (1977 manuscript "crop micrometeorology: a simulation study", Wageningen).
  windReductionFactor = windReductionParam * exposedVAI**twoThirds * canopyHeight**oneThird / leafDimension**oneThird

  ! compute the leaf boundary layer resistance (s m-1)
  singleLeafConductance  = sqrt(windspdCanopyTop/leafDimension)
  leaf2CanopyScaleFactor = (2._dp*leafExchangeCoeff/windReductionFactor) * (1._dp - exp(-windReductionFactor/2._dp)) ! factor to scale from the leaf to the canopy
  canopyLeafConductance  = singleLeafConductance*leaf2CanopyScaleFactor
  leafResistance  = 1._dp/(canopyLeafConductance)

  ! compute eddy diffusivity for heat at the top of the canopy (m2 s-1)
  !   Note: use of friction velocity here includes stability adjustments
  eddyDiffusCanopyTop = max(vkc*FrictionVelocity*(canopyHeight - zeroPlaneDisplacement), mpe)  ! (avoid divide by zero)

  ! compute the resistance between the surface and canopy air
  !  assume exponential profile extends from the surface roughness length to the displacement height plus vegetation roughness
  tmp1 = exp(-windReductionFactor* z0Ground/canopyHeight)
  tmp2 = exp(-windReductionFactor*(z0Canopy+zeroPlaneDisplacement)/canopyHeight)
  groundResistance = ( canopyHeight*exp(windReductionFactor) / (windReductionFactor*eddyDiffusCanopyTop) ) * (tmp1 - tmp2)  ! s m-1

  ! * compute analytical derivatives
  if(derivDesired)then
   ! compute derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
   dCanopyResistance_dTCanopy = -dSfc2AtmExCoef_dTemp/(windspd*sfc2AtmExchangeCoeff**2._dp)
   ! compute derivative in leaf resistance w.r.t. canopy temperature (s m-1 K-1)
   dFV_dT = windspd*dSfc2AtmExCoef_dTemp/(sqrt(sfc2AtmExchangeCoeff)*2._dp)                          ! d(frictionVelocity)/d(canopy temperature)
   dUC_dT = windConvFactor*dFV_dT                                                                    ! d(windspdCanopyTop)/d(canopy temperature)
   dLC_dT = leaf2CanopyScaleFactor*dUC_dT/(leafDimension*sqrt(windspdCanopyTop/leafDimension)*2._dp) ! d(canopyLeafConductance)/d(canopy temperature)
   dLeafResistance_dTCanopy = -dLC_dT/(canopyLeafConductance**2._dp)
   ! compute derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
   dED_dT = dFV_dT*vkc*(canopyHeight - zeroPlaneDisplacement)                                        ! d(eddyDiffusCanopyTop)d(canopy temperature)
   dGroundResistance_dTCanopy = -dED_dT*(tmp1 - tmp2)*canopyHeight*exp(windReductionFactor) / (windReductionFactor*eddyDiffusCanopyTop**2._dp)
  ! * numerical derivatives (computed later)
  else
   dGroundResistance_dTCanopy = 0._dp
   dCanopyResistance_dTCanopy = 0._dp
   dLeafResistance_dTCanopy   = 0._dp
  endif

 ! ***** compute resistances for non-vegetated surfaces (e.g., snow)
 else
  canopyResistance = 0._dp
  leafResistance   = 0._dp
  groundResistance = 1._dp / (sfc2AtmExchangeCoeff*windspd)
  ! set canopy derivatives to zero (non-vegetated, remember)
  dCanopyResistance_dTCanopy = 0._dp
  dLeafResistance_dTCanopy   = 0._dp
  ! compute derivatives for ground resistance
  if(derivDesired)then  ! analytical
   dGroundResistance_dTCanopy = -dSfc2AtmExCoef_dTemp/(windspd*sfc2AtmExchangeCoeff**2._dp)
  else                               ! numerical derivatives (computed later)
   dGroundResistance_dTCanopy = 0._dp
  endif

 endif  ! (switch between vegetated and non-vegetated surfaces)

 end subroutine aeroResist


 ! ********************************************************************************
 ! private subroutine: compute soil moisture factor controlling stomatal resistance
 ! ********************************************************************************
 subroutine soilResist(&
                       ! input (model decisions)
                       ixSoilResist,             & ! intent(in): choice of function for the soil moisture control on stomatal resistance
                       ixGroundwater,            & ! intent(in): choice of groundwater representation
                       ! input (state variables)
                       mLayerMatricHead,         & ! intent(in): matric head in each layer (m)
                       mLayerVolFracLiq,         & ! intent(in): volumetric fraction of liquid water in each layer 
                       scalarAquiferStorage,     & ! intent(in): aquifer storage (m)
                       ! input (diagnostic variables)
                       mLayerRootDensity,        & ! intent(in): root density in each layer (-)
                       scalarAquiferRootFrac,    & ! intent(in): fraction of roots below the lowest unsaturated layer (-)
                       ! input (parameters)
                       plantWiltPsi,             & ! intent(in): matric head at wilting point (m)
                       soilStressParam,          & ! intent(in): parameter in the exponential soil stress function (-)
                       critSoilWilting,          & ! intent(in): critical vol. liq. water content when plants are wilting (-)
                       critSoilTranspire,        & ! intent(in): critical vol. liq. water content when transpiration is limited (-)
                       critAquiferTranspire,     & ! intent(in): critical aquifer storage value when transpiration is limited (m)
                       ! output
                       wAvgTranspireLimitFac,    & ! intent(out): weighted average of the transpiration limiting factor (-)
                       mLayerTranspireLimitFac,  & ! intent(out): transpiration limiting factor in each layer (-)
                       aquiferTranspireLimitFac, & ! intent(out): transpiration limiting factor for the aquifer (-)
                       err,message)                ! intent(out): error control
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 USE mDecisions_module, only: NoahType,CLM_type,SiB_Type  ! options for the choice of function for the soil moisture control on stomatal resistance
 USE mDecisions_module, only: bigBucket                   ! named variable that defines the "bigBucket" groundwater parameterization 
 implicit none
 ! input (model decisions)
 integer(i4b),intent(in)       :: ixSoilResist             ! choice of function for the soil moisture control on stomatal resistance
 integer(i4b),intent(in)       :: ixGroundwater            ! choice of groundwater representation
 ! input (variables)
 real(dp),intent(in)           :: mLayerMatricHead(:)      ! matric head in each layer (m)
 real(dp),intent(in)           :: mLayerVolFracLiq(:)      ! volumetric fraction of liquid water in each layer (-)
 real(dp),intent(in)           :: scalarAquiferStorage     ! aquifer storage (m)
 ! input (diagnostic variables)
 real(dp),intent(in)           :: mLayerRootDensity(:)     ! root density in each layer (-)
 real(dp),intent(in)           :: scalarAquiferRootFrac    ! fraction of roots below the lowest unsaturated layer (-)
 ! input (parameters)
 real(dp),intent(in)           :: plantWiltPsi             ! matric head at wilting point (m)
 real(dp),intent(in)           :: soilStressParam          ! parameter in the exponential soil stress function (-)
 real(dp),intent(in)           :: critSoilWilting          ! critical vol. liq. water content when plants are wilting (-)
 real(dp),intent(in)           :: critSoilTranspire        ! critical vol. liq. water content when transpiration is limited (-)
 real(dp),intent(in)           :: critAquiferTranspire     ! critical aquifer storage value when transpiration is limited (m)
 ! output
 real(dp),intent(out)          :: wAvgTranspireLimitFac    ! intent(out): weighted average of the transpiration limiting factor (-)
 real(dp),intent(out)          :: mLayerTranspireLimitFac(:)  ! intent(out): transpiration limiting factor in each layer (-)
 real(dp),intent(out)          :: aquiferTranspireLimitFac ! intent(out): transpiration limiting factor for the aquifer (-)
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! local variables
 real(dp)                      :: gx                       ! stress function for the soil layers
 integer(i4b)                  :: iLayer                   ! index of soil layer
 ! initialize error control
 err=0; message='soilResist/'

 ! ** compute the factor limiting transpiration for each soil layer (-)
 wAvgTranspireLimitFac = 0._dp  ! (initialize the weighted average)
 do iLayer=1,nSoil
  ! compute the soil stress function
  select case(ixSoilResist)    
   case(NoahType)  ! thresholded linear function of volumetric liquid water content
    gx = (mLayerVolFracLiq(iLayer) - critSoilWilting) / (critSoilTranspire - critSoilWilting)
   case(CLM_type)  ! thresholded linear function of matric head
    gx = 1._dp - mLayerMatricHead(iLayer)/plantWiltPsi
   case(SiB_Type)  ! exponential of the log of matric head
    if(mLayerMatricHead(iLayer) < 0._dp)then  ! (unsaturated)
     gx = 1._dp - exp( -soilStressParam * ( log(plantWiltPsi/mLayerMatricHead(iLayer)) ) )
    else ! (saturated)
     gx = 1._dp
    endif
   case default    ! check identified the option
    err=20; message=trim(message)//'cannot identify option for soil resistance'; return
  end select
  ! save the factor the the given layer (ensure between zero and one)
  mLayerTranspireLimitFac(iLayer) = min(1._dp, max(0._dp,gx) )
  ! compute the weighted average (weighted by root density)
  wAvgTranspireLimitFac = wAvgTranspireLimitFac + mLayerTranspireLimitFac(iLayer)*mLayerRootDensity(iLayer)
 end do ! (looping through soil layers)

 ! ** compute the factor limiting evaporation in the aquifer
 if(scalarAquiferRootFrac > 0._dp)then
  ! check that aquifer root fraction is allowed
  if(ixGroundwater /= bigBucket)then
   message=trim(message)//'aquifer evaporation only allowed for the big groundwater bucket -- increase the soil depth to account for roots'
   err=20; return
  endif
  ! compute the factor limiting evaporation for the aquifer
  aquiferTranspireLimitFac = min(scalarAquiferStorage/critAquiferTranspire, 1._dp)
 else  ! (if there are roots in the aquifer)
  aquiferTranspireLimitFac = 0._dp
 endif
 ! compute the weighted average (weighted by root density)
 wAvgTranspireLimitFac = wAvgTranspireLimitFac + aquiferTranspireLimitFac*scalarAquiferRootFrac

 end subroutine soilResist

 ! ********************************************************************************
 ! private subroutine: compute stomatal resistance
 ! ********************************************************************************
 subroutine stomResist(&
                       ! input (model decisions)
                       ixStomResist,                        & ! intent(in): choice of function for stomatal resistance
                       ! input (local attributes)
                       vegTypeIndex,                        & ! intent(in): vegetation type index
                       iLoc, jLoc,                          & ! intent(in): spatial location indices      
                       ! input (forcing)
                       airtemp,                             & ! intent(in): air temperature at some height above the surface (K)
                       airpres,                             & ! intent(in): air pressure at some height above the surface (Pa)
                       scalarO2air,                         & ! intent(in): atmospheric o2 concentration (Pa)
                       scalarCO2air,                        & ! intent(in): atmospheric co2 concentration (Pa)
                       scalarCanopySunlitPAR,               & ! intent(in): average absorbed par for sunlit leaves (w m-2)
                       scalarCanopyShadedPAR,               & ! intent(in): average absorbed par for shaded leaves (w m-2)
                       ! input (state and diagnostic variables)
                       scalarGrowingSeasonIndex,            & ! intent(in): growing season index (0=off, 1=on)
                       scalarFoliageNitrogenFactor,         & ! intent(in): foliage nitrogen concentration (1=saturated)
                       scalarTranspireLim,                  & ! intent(in): weighted average of the soil moiture factor controlling stomatal resistance (-)
                       scalarLeafResistance,                & ! intent(in): leaf boundary layer resistance (s m-1)
                       scalarVegetationTemp,                & ! intent(in): vegetation temperature (K)
                       scalarSatVP_VegTemp,                 & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                       scalarVP_CanopyAir,                  & ! intent(in): canopy air vapor pressure (Pa)
                       ! output
                       scalarStomResistSunlit,              & ! intent(out): stomatal resistance for sunlit leaves (s m-1)
                       scalarStomResistShaded,              & ! intent(out): stomatal resistance for shaded leaves (s m-1)
                       scalarPhotosynthesisSunlit,          & ! intent(out): sunlit photosynthesis (umolco2 m-2 s-1)
                       scalarPhotosynthesisShaded,          & ! intent(out): shaded photosynthesis (umolco2 m-2 s-1)
                       err,message                          ) ! intent(out): error control
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! Modified from Noah-MP
 ! Compute stomatal resistance and photosynthesis using either
 !  1) Ball-Berry
 !  2) Jarvis
 ! See Niu et al. JGR 2011 for more details
 USE mDecisions_module, only: BallBerry,Jarvis                ! options for the choice of function for stomatal resistance
 USE NOAHMP_ROUTINES,only:stomata                             ! compute canopy resistance based on Ball-Berry
 USE NOAHMP_ROUTINES,only:canres                              ! compute canopy resistance based Jarvis
 implicit none
 ! input (model decisions)
 integer(i4b),intent(in)       :: ixStomResist                ! choice of function for stomatal resistance
 ! input (local attributes)
 integer(i4b),intent(in)       :: vegTypeIndex                ! vegetation type index
 integer(i4b),intent(in)       :: iLoc, jLoc                  ! spatial location indices
 ! input (forcing)
 real(dp),intent(in)           :: airtemp                     ! measured air temperature at some height above the surface (K)
 real(dp),intent(in)           :: airpres                     ! measured air pressure at some height above the surface (Pa)
 real(dp),intent(in)           :: scalarO2air                 ! atmospheric o2 concentration (Pa)
 real(dp),intent(in)           :: scalarCO2air                ! atmospheric co2 concentration (Pa)
 real(dp),intent(in),target    :: scalarCanopySunlitPAR       ! average absorbed par for sunlit leaves (w m-2)
 real(dp),intent(in),target    :: scalarCanopyShadedPAR       ! average absorbed par for shaded leaves (w m-2)
 ! input (state and diagnostic variables)
 real(dp),intent(in)           :: scalarGrowingSeasonIndex    ! growing season index (0=off, 1=on)
 real(dp),intent(in)           :: scalarFoliageNitrogenFactor ! foliage nitrogen concentration (1=saturated)
 real(dp),intent(in)           :: scalarTranspireLim          ! weighted average of the soil moiture factor controlling stomatal resistance (-)
 real(dp),intent(in)           :: scalarLeafResistance        ! leaf boundary layer resistance (s m-1)
 real(dp),intent(in)           :: scalarVegetationTemp        ! vegetation temperature (K)
 real(dp),intent(in)           :: scalarSatVP_VegTemp         ! saturation vapor pressure at vegetation temperature (Pa)
 real(dp),intent(in)           :: scalarVP_CanopyAir          ! canopy air vapor pressure (Pa)
 ! output
 real(dp),intent(out)          :: scalarStomResistSunlit      ! stomatal resistance for sunlit leaves (s m-1)
 real(dp),intent(out)          :: scalarStomResistShaded      ! stomatal resistance for shaded leaves (s m-1)
 real(dp),intent(out)          :: scalarPhotosynthesisSunlit  ! sunlit photosynthesis (umolco2 m-2 s-1)
 real(dp),intent(out)          :: scalarPhotosynthesisShaded  ! sunlit photosynthesis (umolco2 m-2 s-1)
 integer(i4b),intent(out)      :: err                         ! error code
 character(*),intent(out)      :: message                     ! error message
 ! local variables
 integer(i4b),parameter        :: ixSunlit=1                  ! named variable for sunlit leaves
 integer(i4b),parameter        :: ixShaded=2                  ! named variable for shaded leaves
 integer(i4b)                  :: iSunShade                   ! index for sunlit/shaded leaves
 real(dp),pointer              :: PAR                         ! average absorbed PAR for sunlit/shaded leaves (w m-2)
 real(dp)                      :: scalarStomResist            ! stomatal resistance for sunlit/shaded leaves (s m-1)
 real(dp)                      :: scalarPhotosynthesis        ! photosynthesis for sunlit/shaded leaves (umolco2 m-2 s-1)
 ! initialize error control
 err=0; message='stomResist/'

 ! loop through sunlit and shaded leaves
 do iSunShade=1,2

  ! get appropriate value for PAR
  select case(iSunShade)
   case(ixSunlit); PAR => scalarCanopySunlitPAR               ! average absorbed par for sunlit leaves (w m-2)
   case(ixShaded); PAR => scalarCanopyShadedPAR               ! average absorbed par for shaded leaves (w m-2)
   case default; err=20; message=trim(message)//'unable to identify case for sunlit/shaded leaves'; return
  end select

  ! identify option for stomatal resistance
  select case(ixStomResist)

   ! Ball-Berry
   case(BallBerry)
   call stomata(&
                ! input
                vegTypeIndex,                       & ! intent(in): vegetation type index
                mpe,                                & ! intent(in): prevents overflow error if division by zero
                PAR,                                & ! intent(in): average absorbed par (w m-2)
                scalarFoliageNitrogenFactor,        & ! intent(in): foliage nitrogen concentration (1=saturated)
                iLoc, jLoc,                         & ! intent(in): spatial location indices      
                scalarVegetationTemp,               & ! intent(in): vegetation temperature (K)
                scalarSatVP_VegTemp,                & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                scalarVP_CanopyAir,                 & ! intent(in): canopy air vapor pressure (Pa)
                airtemp,                            & ! intent(in): air temperature at some height above the surface (K)
                airpres,                            & ! intent(in): air pressure at some height above the surface (Pa)
                scalarO2air,                        & ! intent(in): atmospheric o2 concentration (Pa)
                scalarCO2air,                       & ! intent(in): atmospheric co2 concentration (Pa)
                scalarGrowingSeasonIndex,           & ! intent(in): growing season index (0=off, 1=on)
                scalarTranspireLim,                 & ! intent(in): weighted average of the soil moiture factor controlling stomatal resistance (-)
                scalarLeafResistance,               & ! intent(in): leaf boundary layer resistance (s m-1)
                ! output
                scalarStomResist,                   & ! intent(out): stomatal resistance (s m-1)
                scalarPhotosynthesis                ) ! intent(out): photosynthesis (umolco2 m-2 s-1)

   ! Jarvis
   case(Jarvis)
   call canres(&
                ! input
                PAR,                                & ! intent(in): average absorbed par (w m-2)
                scalarVegetationTemp,               & ! intent(in): vegetation temperature (K)
                scalarTranspireLim,                 & ! intent(in): weighted average of the soil moiture factor controlling stomatal resistance (-)
                scalarVP_CanopyAir,                 & ! intent(in): canopy air vapor pressure (Pa)
                airpres,                            & ! intent(in): air pressure at some height above the surface (Pa)
                ! output
                scalarStomResist,                   & ! intent(out): stomatal resistance (s m-1)
                scalarPhotosynthesis,               & ! intent(out): photosynthesis (umolco2 m-2 s-1)
                ! location indices (input)
                iLoc, jLoc                          ) ! intent(in): spatial location indices      

   ! check identified an option
   case default; err=20; message=trim(message)//'unable to identify case for stomatal resistance'; return

  end select  ! (selecting option for stomatal resistance)

  ! assign output variables
  select case(iSunShade)
   case(ixSunlit)
    scalarStomResistSunlit     = scalarStomResist
    scalarPhotosynthesisSunlit = scalarPhotosynthesis
   case(ixShaded)
    scalarStomResistShaded     = scalarStomResist
    scalarPhotosynthesisShaded = scalarPhotosynthesis
   case default; err=20; message=trim(message)//'unable to identify case for sunlit/shaded leaves'; return
  end select

 end do  ! (looping through sunlit and shaded leaves)

 end subroutine stomResist



 ! ********************************************************************************
 ! private subroutine: compute turbulent heat fluxes
 ! ********************************************************************************
 subroutine turbFluxes(&
                       ! input: model control
                       dt,                            & ! intent(in): model time step (seconds)
                       ixDerivMethod,                 & ! intent(in): choice of method used to compute derivative (analytical or numerical)
                       ! input: above-canopy forcing data
                       airtemp,                       & ! intent(in): air temperature at some height above the surface (K)
                       VPair,                         & ! intent(in): vapor pressure of the air above the vegetation canopy (Pa)
                       ! input: psychometric constant for canopy and ground
                       psycConstCanopy,               & ! psychrometric constant for the vegetation canopy (Pa/K)
                       psycConstGround,               & ! psychrometric constant for the ground surface (Pa/K)
                       ! input: canoopy liquid and canopy ice (used as a solution constraint)
                       canopyLiquid,                  & ! intent(in): mass of liquid water on the vegetation canopy (kg m-2)
                       canopyIce,                     & ! intent(in): mass of ice water on the vegetation canopy (kg m-2)
                       ! input: canopy and ground temperature
                       canopyTemp,                    & ! intent(in): canopy temperature (K)
                       groundTemp,                    & ! intent(in): ground temperature (K)
                       satVP_CanopyTemp,              & ! intent(in): saturation vapor pressure at the temperature of the veg canopy (Pa)
                       satVP_GroundTemp,              & ! intent(in): saturation vapor pressure at the temperature of the ground (Pa)
                       dSVPCanopy_dCanopyTemp,        & ! intent(in): derivative in canopy saturation vapor pressure w.r.t. canopy temperature (Pa K-1)
                       dSVPGround_dGroundTemp,        & ! intent(in): derivative in ground saturation vapor pressure w.r.t. ground temperature (Pa K-1)
                       ! input: diagnostic variables
                       exposedVAI,                    & ! intent(in): exposed vegetation area index -- leaf plus stem (m2 m-2)
                       canopyWetFraction,             & ! intent(in): fraction of canopy that is wet [0-1]
                       canopySunlitLAI,               & ! intent(in): sunlit leaf area (-)
                       canopyShadedLAI,               & ! intent(in): shaded leaf area (-)
                       soilRelHumidity,               & ! intent(in): relative humidity in the soil pores [0-1]
                       soilResistance,                & ! intent(in): resistance from the soil (s m-1)
                       leafResistance,                & ! intent(in): mean leaf boundary layer resistance per unit leaf area (s m-1)
                       groundResistance,              & ! intent(in): below canopy aerodynamic resistance (s m-1) 
                       canopyResistance,              & ! intent(in): above canopy aerodynamic resistance (s m-1)
                       stomResistSunlit,              & ! intent(in): stomatal resistance for sunlit leaves (s m-1)
                       stomResistShaded,              & ! intent(in): stomatal resistance for shaded leaves (s m-1)
                       ! input: derivatives in scalar resistances
                       dGroundResistance_dTCanopy,    & ! intent(in): derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
                       dCanopyResistance_dTCanopy,    & ! intent(in): derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
                       dLeafResistance_dTCanopy,      & ! intent(in): derivative in leaf resistance w.r.t. canopy temperature (s m-1 K-1)
                       ! output: canopy air space variables
                       temp_CanopyAir,                & ! intent(out): temperature of the canopy air space (K)
                       VP_CanopyAir,                  & ! intent(out): vapor pressure of the canopy air space (Pa)
                       ! output: fluxes from the vegetation canopy
                       senHeatCanopy,                 & ! intent(out): sensible heat flux from the canopy to the canopy air space (W m-2)
                       latHeatCanopyEvap,             & ! intent(out): latent heat flux associated with evaporation from the canopy to the canopy air space (W m-2)
                       latHeatCanopyTrans,            & ! intent(out): latent heat flux associated with transpiration from the canopy to the canopy air space (W m-2)
                       ! output: fluxes from non-vegetated surfaces (ground surface below vegetation, bare ground, or snow covered vegetation)
                       senHeatGround,                 & ! intent(out): sensible heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
                       latHeatGround,                 & ! intent(out): latent heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
                       ! output: net fluxes
                       turbFluxCanopy,                & ! intent(out): net longwave radiation at the canopy (W m-2)
                       turbFluxGround,                & ! intent(out): net longwave radiation at the ground surface (W m-2)
                       ! output: flux derivatives
                       dTurbFluxCanopy_dTCanopy,      & ! intent(out): derivative in net canopy turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
                       dTurbFluxGround_dTCanopy,      & ! intent(out): derivative in net ground turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
                       dTurbFluxCanopy_dTGround,      & ! intent(out): derivative in net canopy turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
                       dTurbFluxGround_dTGround,      & ! intent(out): derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
                       ! output: error control
                       err,message                    ) ! intent(out): error control
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 implicit none
 ! input: model control
 real(dp),intent(in)           :: dt                    ! model time step (seconds)
 integer(i4b),intent(in)       :: ixDerivMethod         ! choice of method used to compute derivative (analytical or numerical)
 ! input: above-canopy forcing data
 real(dp),intent(in)           :: airtemp               ! air temperature at some height above the surface (K)
 real(dp),intent(in)           :: VPair                 ! vapor pressure of the air above the vegetation canopy (Pa)
 ! input: psychometric constant for canopy and ground
 real(dp),intent(in)           :: psycConstCanopy       ! psychrometric constant for the vegetation canopy (Pa/K)
 real(dp),intent(in)           :: psycConstGround       ! psychrometric constant for the ground surface (Pa/K)
 ! input: canopy liquid and canopy ice (used as a solution constraint)
 real(dp),intent(in)           :: canopyLiquid          ! mass of liquid water on the vegetation canopy (kg m-2)
 real(dp),intent(in)           :: canopyIce             ! mass of ice water on the vegetation canopy (kg m-2)
 ! input: canopy and ground temperature
 real(dp),intent(in)           :: canopyTemp            ! canopy temperature (K)
 real(dp),intent(in)           :: groundTemp            ! ground temperature (K)
 real(dp),intent(in)           :: satVP_CanopyTemp      ! saturation vapor pressure at the temperature of the veg canopy (Pa)
 real(dp),intent(in)           :: satVP_GroundTemp      ! saturation vapor pressure at the temperature of the ground (Pa)
 real(dp),intent(in)           :: dSVPCanopy_dCanopyTemp  ! derivative in canopy saturation vapor pressure w.r.t. canopy temperature (Pa K-1)
 real(dp),intent(in)           :: dSVPGround_dGroundTemp  ! derivative in ground saturation vapor pressure w.r.t. ground temperature (Pa K-1)
 ! input: diagnostic variables
 real(dp),intent(in)           :: exposedVAI            ! exposed vegetation area index -- leaf plus stem (m2 m-2)
 real(dp),intent(in)           :: canopyWetFraction     ! fraction of canopy that is wet [0-1]
 real(dp),intent(in)           :: canopySunlitLAI       ! sunlit leaf area (-)
 real(dp),intent(in)           :: canopyShadedLAI       ! shaded leaf area (-)
 real(dp),intent(in)           :: soilRelHumidity       ! relative humidity in the soil pores [0-1]
 real(dp),intent(in)           :: soilResistance        ! resistance from the soil (s m-1)
 real(dp),intent(in)           :: leafResistance        ! mean leaf boundary layer resistance per unit leaf area (s m-1)
 real(dp),intent(in)           :: groundResistance      ! below canopy aerodynamic resistance (s m-1) 
 real(dp),intent(in)           :: canopyResistance      ! above canopy aerodynamic resistance (s m-1)
 real(dp),intent(in)           :: stomResistSunlit      ! stomatal resistance for sunlit leaves (s m-1)
 real(dp),intent(in)           :: stomResistShaded      ! stomatal resistance for shaded leaves (s m-1)
 ! input: derivatives in scalar resistances
 real(dp),intent(in)           :: dGroundResistance_dTCanopy      ! derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
 real(dp),intent(in)           :: dCanopyResistance_dTCanopy      ! derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
 real(dp),intent(in)           :: dLeafResistance_dTCanopy        ! derivative in leaf resistance w.r.t. canopy temperature (s m-1 K-1)
 ! output: canopy air space variables
 real(dp),intent(out)          :: temp_CanopyAir        ! temperature of the canopy air space (K)
 real(dp),intent(out)          :: VP_CanopyAir          ! vapor pressure of the canopy air space (Pa)
 ! output: fluxes from the vegetation canopy
 real(dp),intent(out)          :: senHeatCanopy         ! sensible heat flux from the canopy to the canopy air space (W m-2)
 real(dp),intent(out)          :: latHeatCanopyEvap     ! latent heat flux associated with evaporation from the canopy to the canopy air space (W m-2)
 real(dp),intent(out)          :: latHeatCanopyTrans    ! latent heat flux associated with transpiration from the canopy to the canopy air space (W m-2)
 ! output: fluxes from non-vegetated surfaces (ground surface below vegetation, bare ground, or snow covered vegetation)
 real(dp),intent(out)          :: senHeatGround         ! sensible heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
 real(dp),intent(out)          :: latHeatGround         ! latent heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
 ! output: net fluxes
 real(dp),intent(out)          :: turbFluxCanopy        ! net longwave radiation at the canopy (W m-2)
 real(dp),intent(out)          :: turbFluxGround        ! net longwave radiation at the ground surface (W m-2)
 ! output: flux derivatives
 real(dp),intent(out)          :: dTurbFluxCanopy_dTCanopy       ! derivative in net canopy turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
 real(dp),intent(out)          :: dTurbFluxGround_dTCanopy       ! derivative in net ground turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
 real(dp),intent(out)          :: dTurbFluxCanopy_dTGround       ! derivative in net canopy turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
 real(dp),intent(out)          :: dTurbFluxGround_dTGround       ! derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                   ! error code
 character(*),intent(out)      :: message               ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! local variables -- general
 real(dp)                      :: dpart1,dpart2                ! derivatives for different parts of a function
 ! local variables -- "constants"
 real(dp)                      :: volHeatCapacityAir           ! volumetric heat capacity of air (J m-3)
 real(dp)                      :: latHeatConstCanopy           ! latent heat constant for the canopy (J m-3 K-1)
 real(dp)                      :: latHeatConstGround           ! latent heat constant for the ground (J m-3 K-1)
 ! local variables -- conductance
 real(dp)                      :: leafConductance              ! leaf conductance (m s-1)
 real(dp)                      :: canopyConductance            ! canopy conductance (m s-1)
 real(dp)                      :: groundConductanceSH          ! ground conductance for sensible heat (m s-1)
 real(dp)                      :: groundConductanceLH          ! ground conductance for latent heat -- includes soil resistance (m s-1)
 real(dp)                      :: evapConductance              ! conductance for evaporation (m s-1)
 real(dp)                      :: transConductance             ! conductance for transpiration (m s-1)
 real(dp)                      :: totalConductanceSH           ! total conductance for sensible heat (m s-1)
 real(dp)                      :: totalConductanceLH           ! total conductance for latent heat (m s-1)
 ! local variables -- derivatives for energy conductances
 real(dp)                      :: dLeafCond_dCanopyTemp        ! derivative in leaf conductance w.r.t. canopy temperature
 real(dp)                      :: dCanopyCond_dCanopyTemp      ! derivative in canopy conductance w.r.t. canopy temperature
 real(dp)                      :: dGroundCondSH_dCanopyTemp    ! derivative in ground conductance of sensible heat w.r.t. canopy temperature
 ! local variables -- derivatives for mass conductances 
 real(dp)                      :: dTempSunlit                  ! derivative in sunlit conductance w.r.t. canopy temperature
 real(dp)                      :: dTempShaded                  ! derivative in shaded conductance w.r.t. canopy temperature
 real(dp)                      :: dEvapCond_dCanopyTemp        ! derivative in evaporation conductance w.r.t. canopy temperature
 real(dp)                      :: dTransCond_dCanopyTemp       ! derivative in transpiration conductance w.r.t. canopy temperature
 real(dp)                      :: dGroundCondLH_dCanopyTemp    ! derivative in ground conductance w.r.t. canopy temperature
 ! local variables -- derivatives for the canopy air space variables
 real(dp)                      :: fPart_Temp                   ! part of the function for temperature of the canopy air space
 real(dp)                      :: fPart_VP                     ! part of the function for vapor pressure of the canopy air space
 real(dp)                      :: dTempCanopyAir_dTCanopy      ! derivative in the temperature of the canopy air space w.r.t. temperature of the canopy
 real(dp)                      :: dTempCanopyAir_dTGround      ! derivative in the temperature of the canopy air space w.r.t. temperature of the ground
 real(dp)                      :: dVPCanopyAir_dTCanopy        ! derivative in the vapor pressure of the canopy air space w.r.t. temperature of the canopy 
 real(dp)                      :: dVPCanopyAir_dTGround        ! derivative in the vapor pressure of the canopy air space w.r.t. temperature of the ground
 ! local variables -- sensible heat flux derivatives 
 real(dp)                      :: dSenHeatCanopy_dTCanopy      ! derivative in the canopy sensible heat flux w.r.t. canopy temperature
 real(dp)                      :: dSenHeatGround_dTCanopy      ! derivative in the ground sensible heat flux w.r.t. canopy temperature
 real(dp)                      :: dSenHeatCanopy_dTGround      ! derivative in the canopy sensible heat flux w.r.t. ground temperature
 real(dp)                      :: dSenHeatGround_dTGround      ! derivative in the ground sensible heat flux w.r.t. ground temperature
 ! local variables -- latent heat flux derivatives
 real(dp)                      :: dLatHeatCanopyEvap_dTCanopy  ! derivative in the canopy evaporation flux w.r.t. canopy temperature
 real(dp)                      :: dLatHeatCanopyTrans_dTCanopy ! derivative in the canopy transpiration flux w.r.t. canopy temperature
 real(dp)                      :: dLatHeatGround_dTCanopy      ! derivative in the ground latent heat flux w.r.t. canopy temperature
 real(dp)                      :: dLatHeatCanopyEvap_dTGround  ! derivative in the canopy evaporation flux w.r.t. ground temperature
 real(dp)                      :: dLatHeatCanopyTrans_dTGround ! derivative in the canopy transpiration flux w.r.t. ground temperature
 real(dp)                      :: dLatHeatGround_dTGround      ! derivative in the ground latent heat flux w.r.t. ground temperature
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='turbFluxes/'

 ! compute constants
 volHeatCapacityAir = iden_air*cp_air                          ! volumetric heat capacity of air (J m-3)
 latHeatConstCanopy = volHeatCapacityAir/psycConstCanopy       ! latent heat constant for the canopy (J m-3 K-1)
 latHeatConstGround = volHeatCapacityAir/psycConstGround       ! latent heat constant for the ground (J m-3 K-1)

 ! *****
 ! * compute conductances, and derivatives...
 ! ******************************************

 ! compute conductances for sensible heat (m s-1)
 leafConductance     = exposedVAI/leafResistance
 canopyConductance   = 1._dp/canopyResistance
 groundConductanceSH = 1._dp/groundResistance
 totalConductanceSH  = leafConductance + groundConductanceSH + canopyConductance

 ! compute conductances for latent heat (m s-1)
 evapConductance     = canopyWetFraction*leafConductance
 transConductance    = (1._dp - canopyWetFraction) * ( canopySunlitLAI/(leafResistance+stomResistSunlit) + canopyShadedLAI/(leafResistance+stomResistShaded) )
 groundConductanceLH = 1._dp/(groundResistance + soilResistance)
 totalConductanceLH  = evapConductance + transConductance + groundConductanceLH + canopyConductance

 ! * compute derivatives
 ! NOTE: it may be more efficient to compute these derivatives when computing resistances
 if(ixDerivMethod == analytical)then

  ! compute derivatives in individual conductances for sensible heat w.r.t. canopy temperature (m s-1 K-1)
  dLeafCond_dCanopyTemp     = -exposedVAI*dLeafResistance_dTCanopy/leafResistance**2._dp  ! derivative in leaf conductance w.r.t. canopy temperature
  dCanopyCond_dCanopyTemp   = -dCanopyResistance_dTCanopy/canopyResistance**2._dp         ! derivative in canopy conductance w.r.t. canopy temperature
  dGroundCondSH_dCanopyTemp = -dGroundResistance_dTCanopy/groundResistance**2._dp         ! derivative in ground conductance w.r.t. canopy temperature

  ! compute derivatives in individual conductances for latent heat w.r.t. canopy temperature (m s-1 K-1)
  dTempSunlit = -canopySunlitLAI*dLeafResistance_dTCanopy/(leafResistance+stomResistSunlit)**2._dp ! derivative in sunlit conductance w.r.t. canopy temperature
  dTempShaded = -canopyShadedLAI*dLeafResistance_dTCanopy/(leafResistance+stomResistShaded)**2._dp ! derivative in shaded conductance w.r.t. canopy temperature
  dEvapCond_dCanopyTemp     = canopyWetFraction*dLeafCond_dCanopyTemp                              ! derivative in evaporation conductance w.r.t. canopy temperature
  dTransCond_dCanopyTemp    = (1._dp - canopyWetFraction) * (dTempSunlit + dTempShaded)            ! derivative in transpiration conductance w.r.t. canopy temperature
  dGroundCondLH_dCanopyTemp = -dGroundResistance_dTCanopy/(groundResistance+soilResistance)**2._dp ! derivative in ground conductance w.r.t. canopy temperature

 endif ! (if computing analytical derivatives)


 ! *****
 ! * compute canopy air space variables, and derivatives...
 ! ********************************************************

 ! compute the temperature in the canopy air space (K)
 fPart_Temp     = canopyConductance*airtemp + leafConductance*canopyTemp + groundConductanceSH*groundTemp
 temp_CanopyAir = fPart_Temp/totalConductanceSH

 ! compute the vapor pressure in the canopy air space (Pa)
 fPart_VP     = canopyConductance*VPair + (evapConductance + transConductance)*satVP_CanopyTemp + groundConductanceLH*satVP_GroundTemp*soilRelHumidity
 VP_CanopyAir = fPart_VP/totalConductanceLH

 ! * compute derivatives
 if(ixDerivMethod == analytical)then

  ! compute derivative in temperature of the canopy air space w.r.t. canopy temperature (product rule)
  dPart1 = airtemp*dCanopyCond_dCanopyTemp + (leafConductance + canopyTemp*dLeafCond_dCanopyTemp) + groundTemp*dGroundCondSH_dCanopyTemp
  dPart2 = -(dLeafCond_dCanopyTemp + dGroundCondSH_dCanopyTemp + dCanopyCond_dCanopyTemp)/totalConductanceSH**2._dp
  dTempCanopyAir_dTCanopy = dPart1/totalConductanceSH + fPart_Temp*dPart2

  ! compute derivative in temperature of the canopy air space w.r.t. ground temperature
  dTempCanopyAir_dTGround = groundConductanceSH/totalConductanceSH

  ! compute derivative in vapor pressure of the canopy air space w.r.t. canopy temperature (product rule)
  dPart1 = VPair*dCanopyCond_dCanopyTemp + &                                                                                                   ! derivative in canopy-atmosphere conductance w.r.t. canopy temperature
           (dEvapCond_dCanopyTemp + dTransCond_dCanopyTemp)*satVP_CanopyTemp + (evapConductance + transConductance)*dSVPCanopy_dCanopyTemp + & ! derivative in mass conductance from canopy w.r.t. canopy temperature
           dGroundCondLH_dCanopyTemp*satVP_GroundTemp*soilRelHumidity                                                                          ! derivative in mass conductance from ground w.r.t. canopy temperature
  dPart2 = -(dEvapCond_dCanopyTemp + dTransCond_dCanopyTemp + dGroundCondLH_dCanopyTemp + dCanopyCond_dCanopyTemp)/totalConductanceLH**2._dp
  dVPCanopyAir_dTCanopy = dPart1/totalConductanceLH + fPart_VP*dPart2

  ! compute derivative in vapor pressure of the canopy air space w.r.t. ground temperature 
  dVPCanopyAir_dTGround = dSVPGround_dGroundTemp*groundConductanceLH*soilRelHumidity/totalConductanceLH

 endif ! (if computing analytical derivatives)


 ! *****
 ! * compute sensible and latent heat fluxes, and derivatives...
 ! *************************************************************

 ! compute sensible and latent heat fluxes from the canopy to the canopy air space (W m-2)
 senHeatCanopy      = -volHeatCapacityAir*leafConductance*(canopyTemp - temp_CanopyAir)        ! (positive downwards)
 latHeatCanopyEvap  = -latHeatConstCanopy*evapConductance*(satVP_CanopyTemp - VP_CanopyAir)    ! (positive downwards)
 latHeatCanopyTrans = -latHeatConstCanopy*transConductance*(satVP_CanopyTemp - VP_CanopyAir)   ! (positive downwards)

 ! check that energy for canopy evaporation does not exhaust the available water
 ! NOTE: do this here, rather than enforcing solution constraints, because energy and mass solutions may be uncoupled
 if(canopyTemp > Tfreeze)then
  latHeatCanopyEvap = min(latHeatCanopyEvap, (canopyLiquid + canopyIce)*LH_vap/dt)
 else
  latHeatCanopyEvap = min(latHeatCanopyEvap, (canopyLiquid + canopyIce)*LH_sub/dt)
 endif

 ! compute sensible and latent heat fluxes from the ground to the canopy air space (W m-2)
 senHeatGround      = -volHeatCapacityAir*groundConductanceSH*(groundTemp - temp_CanopyAir)                      ! (positive downwards)
 latHeatGround      = -latHeatConstGround*groundConductanceLH*(satVP_GroundTemp*soilRelHumidity - VP_CanopyAir)  ! (positive downwards)

 ! * compute derivatives
 if(ixDerivMethod == analytical)then

  ! differentiate CANOPY fluxes w.r.t. canopy temperature (product rule)
  dSenHeatCanopy_dTCanopy      = (-volHeatCapacityAir*dLeafCond_dCanopyTemp)*(canopyTemp - temp_CanopyAir) + &            ! d(canopy sensible heat flux)/d(canopy temp)
                                 (-volHeatCapacityAir*leafConductance)*(1._dp - dTempCanopyAir_dTCanopy)
  dLatHeatCanopyEvap_dTCanopy  = (-latHeatConstCanopy*dEvapCond_dCanopyTemp)*(satVP_CanopyTemp - VP_CanopyAir) + &        ! d(canopy evaporation flux)/d(canopy temp)
                                 (-latHeatConstCanopy*evapConductance)*(dSVPCanopy_dCanopyTemp - dVPCanopyAir_dTCanopy)
  dLatHeatCanopyTrans_dTCanopy = (-latHeatConstCanopy*dTransCond_dCanopyTemp)*(satVP_CanopyTemp - VP_CanopyAir) + &       ! d(canopy transpiration flux)/d(canopy temp)
                                 (-latHeatConstCanopy*transConductance)*(dSVPCanopy_dCanopyTemp - dVPCanopyAir_dTCanopy)

  ! differentiate CANOPY fluxes w.r.t. ground temperature
  dSenHeatCanopy_dTGround      = volHeatCapacityAir*leafConductance*dTempCanopyAir_dTGround
  dLatHeatCanopyEvap_dTGround  = latHeatConstCanopy*evapConductance*dVPCanopyAir_dTGround
  dLatHeatCanopyTrans_dTGround = latHeatConstCanopy*transConductance*dVPCanopyAir_dTGround


  ! differentiate GROUND fluxes w.r.t. canopy temperature (product rule)
  dSenHeatGround_dTCanopy = (-volHeatCapacityAir*dGroundCondSH_dCanopyTemp)*(groundTemp - temp_CanopyAir) + &                            ! d(ground sensible heat flux)/d(canopy temp)
                            (-volHeatCapacityAir*groundConductanceSH)*(0._dp - dTempCanopyAir_dTCanopy)
  dLatHeatGround_dTCanopy = (-latHeatConstGround*dGroundCondLH_dCanopyTemp)*(satVP_GroundTemp*soilRelHumidity - VP_CanopyAir) + &        ! d(ground latent heat flux)/d(canopy temp)
                            (-latHeatConstCanopy*groundConductanceLH)*(0._dp - dVPCanopyAir_dTCanopy)

  ! differentiate GROUND fluxes w.r.t. ground temperature
  dSenHeatGround_dTGround = -volHeatCapacityAir*groundConductanceSH*(1._dp - dTempCanopyAir_dTGround)
  dLatHeatGround_dTGround = -latHeatConstGround*groundConductanceLH*(dSVPGround_dGroundTemp*soilRelHumidity - dVPCanopyAir_dTGround)

 endif  ! (if computing analytical derivatives)


 ! *****
 ! * compute net turbulent fluxes, and derivatives...
 ! **************************************************

 ! compute net fluxes
 turbFluxCanopy = senHeatCanopy + latHeatCanopyEvap + latHeatCanopyTrans  ! net turbulent flux at the canopy (W m-2)
 turbFluxGround = senHeatGround + latHeatGround                           ! net turbulent flux at the ground surface (W m-2)

  ! * compute derivatives
 if(ixDerivMethod == analytical)then
  dTurbFluxCanopy_dTCanopy = dSenHeatCanopy_dTCanopy + dLatHeatCanopyEvap_dTCanopy + dLatHeatCanopyTrans_dTCanopy  ! derivative in net canopy turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
  dTurbFluxGround_dTCanopy = dSenHeatGround_dTCanopy + dLatHeatGround_dTCanopy                                     ! derivative in net ground turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
  dTurbFluxCanopy_dTGround = dSenHeatCanopy_dTGround + dLatHeatCanopyEvap_dTGround + dLatHeatCanopyTrans_dTGround  ! derivative in net canopy turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
  dTurbFluxGround_dTGround = dSenHeatGround_dTGround + dLatHeatGround_dTGround
 else ! (just make sure we return something)
  dTurbFluxCanopy_dTCanopy = 0._dp 
  dTurbFluxGround_dTCanopy = 0._dp 
  dTurbFluxCanopy_dTGround = 0._dp
  dTurbFluxGround_dTGround = 0._dp
 endif
 
 end subroutine turbFluxes


 ! ***********************************************************************************************************
 ! private subroutine: compute the surface-atmosphere turbulent exchange coefficient (-)
 ! *********************************************************************************************************** 
 subroutine turbExCoef(&
                       ! input: control
                       computeDerivative,              & ! input: logical flag to compute analytical derivatives
                       ixStability,                    & ! input: choice of stability function
                       ! input: forcing data, diagnostic and state variables
                       mHeight,                        & ! input: measurement height (m)
                       zeroPlaneDisplacement,          & ! input: zero plane displacement (m)
                       z0,                             & ! input: roughness length (m)
                       airTemp,                        & ! input: air temperature (K)
                       sfcTemp,                        & ! input: surface temperature (K)
                       windspd,                        & ! input: wind speed (m s-1)
                       ! input: stability parameters
                       critRichNumber,                 & ! input: critical value for the bulk Richardson number where turbulence ceases (-)
                       Louis79_bparam,                 & ! input: parameter in Louis (1979) stability function
                       Louis79_cStar,                  & ! input: parameter in Louis (1979) stability function
                       Mahrt87_eScale,                 & ! input: exponential scaling factor in the Mahrt (1987) stability function
                       ! output
                       sfc2AtmExchangeCoeff,           & ! output: surface-atmosphere turbulent exchange coefficient (-)
                       dSfc2AtmExCoef_dTemp,           & ! output: derivative in surface-atmosphere exchange coeffcient w.r.t. temperature (K-1)
                       err, message                    ) ! output: error control
 implicit none
 ! input: control
 logical(lgt),intent(in)       :: computeDerivative      ! flag to compute the derivative
 integer(i4b),intent(in)       :: ixStability            ! choice of stability function
 ! input: forcing data, diagnostic and state variables
 real(dp),intent(in)           :: mHeight                ! measurement height (m)
 real(dp),intent(in)           :: zeroPlaneDisplacement  ! zero plane displacement (m)
 real(dp),intent(in)           :: z0                     ! roughness length (m)
 real(dp),intent(in)           :: airtemp                ! air temperature (K)
 real(dp),intent(in)           :: sfcTemp                ! surface temperature (K)
 real(dp),intent(in)           :: windspd                ! wind speed (m s-1)
 ! input: stability parameters
 real(dp),intent(in)           :: critRichNumber         ! critical value for the bulk Richardson number where turbulence ceases (-)
 real(dp),intent(in)           :: Louis79_bparam         ! parameter in Louis (1979) stability function
 real(dp),intent(in)           :: Louis79_cStar          ! parameter in Louis (1979) stability function
 real(dp),intent(in)           :: Mahrt87_eScale         ! exponential scaling factor in the Mahrt (1987) stability function
 ! output
 real(dp),intent(out)          :: sfc2AtmExchangeCoeff   ! surface-atmosphere turbulent exchange coefficient (-)
 real(dp),intent(out)          :: dSfc2AtmExCoef_dTemp   ! derivative in surface-atmosphere exchange coeffcient w.r.t. temperature (K-1)
 integer(i4b),intent(out)      :: err                    ! error code
 character(*),intent(out)      :: message                ! error message
 ! local
 real(dp)                      :: ExNeut                 ! turbulent exchange coefficient under neutral conditions (-)
 real(dp)                      :: RiBulk                 ! bulk Richardson number (-)
 real(dp)                      :: dRiBulk_dTemp          ! derivative in the bulk Richardson number w.r.t. temperature (K-1)
 real(dp)                      :: bPrime                 ! scaled "b" parameter for stability calculations in Louis (1979)
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='turbExCoef/'

 ! check that the measurement height is above the displacement height
 if(mHeight < zeroPlaneDisplacement)then; err=20; message=trim(message)//'measurement height is below the displacement height'; return; endif

 ! check that the measurement height is above the roughness length
 if(mHeight < z0)then; err=20; message=trim(message)//'measurement height is below the roughness length'; return; endif

 ! compute turbulent exchange coefficient under conditions of neutral stability
 ExNeut = (vkc**2._dp) / ( log((mHeight - zeroPlaneDisplacement)/z0))**2._dp

 ! compute the bulk Richardson number (-)
 call bulkRichardson(airTemp,sfcTemp,windspd,mHeight,computeDerivative, & ! (input)
                     RiBulk,dRiBulk_dTemp,err,message)                    ! (output)

 ! set derivative to one if not computing it
 if(.not.computeDerivative) dSfc2AtmExCoef_dTemp = 1._dp

 ! ***** process unstable cases
 if(RiBulk<0._dp)then
  ! compute surface-atmosphere exchange coefficient (-)
  sfc2AtmExchangeCoeff = ExNeut * (1._dp - 16._dp*RiBulk)**0.5_dp
  ! compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
  if(computeDerivative)&
   dSfc2AtmExCoef_dTemp = dRiBulk_dTemp * (-16._dp) * 0.5_dp*(1._dp - 16._dp*RiBulk)**(-0.5_dp) * ExNeut
  return
 endif

 ! ***** process stable cases
 select case(ixStability)

  ! ("standard" stability correction, a la Anderson 1976)
  case(standard)
   ! compute surface-atmosphere exchange coefficient (-)   
   if(RiBulk <  critRichNumber) sfc2AtmExchangeCoeff = ExNeut * (1._dp - 5._dp*RiBulk)**2._dp
   if(RiBulk >= critRichNumber) sfc2AtmExchangeCoeff = epsilon(sfc2AtmExchangeCoeff)
   ! compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
   if(computeDerivative)then
    if(RiBulk <  critRichNumber) dSfc2AtmExCoef_dTemp = dRiBulk_dTemp * (-5._dp) * 2._dp*(1._dp - 5._dp*RiBulk) * ExNeut
    if(RiBulk >= critRichNumber) dSfc2AtmExCoef_dTemp = 0._dp
   endif

  ! (Louis 1979)
  case(louisInversePower)
   ! scale the "b" parameter for stable conditions
   bprime = Louis79_bparam/2._dp
   ! compute surface-atmosphere exchange coefficient (-)
   sfc2AtmExchangeCoeff = ExNeut / ( (1._dp + bprime*RiBulk)**2._dp )
   if(sfc2AtmExchangeCoeff < epsilon(sfc2AtmExchangeCoeff)) sfc2AtmExchangeCoeff = epsilon(sfc2AtmExchangeCoeff)
   ! compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
   if(computeDerivative)&
    dSfc2AtmExCoef_dTemp = dRiBulk_dTemp * bprime * (-2._dp)*(1._dp + bprime*RiBulk)**(-3._dp) * ExNeut

  ! (Mahrt 1987)
  case(mahrtExponential)
   ! compute surface-atmosphere exchange coefficient (-)
   sfc2AtmExchangeCoeff = ExNeut * exp(-Mahrt87_eScale * RiBulk)
   if(sfc2AtmExchangeCoeff < epsilon(sfc2AtmExchangeCoeff)) sfc2AtmExchangeCoeff = epsilon(sfc2AtmExchangeCoeff)
   ! compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
   if(computeDerivative)&
    dSfc2AtmExCoef_dTemp = dRiBulk_dTemp * (-Mahrt87_eScale) * exp(-Mahrt87_eScale * RiBulk) * ExNeut

  ! (return error if the stability correction method is not found)
  case default
   err=10; message=trim(message)//"optionNotFound[stability correction]"; return

 endselect
 
 end subroutine turbExCoef


 ! ***********************************************************************************************************
 ! private subroutine: compute bulk Richardson number
 ! *********************************************************************************************************** 
 subroutine bulkRichardson(airtemp,sfcTemp,windspd,mHeight,computeDerivative, & ! (input)
                           RiBulk,dRiBulk_dTemp,err,message)                    ! (output)
 implicit none
 ! input
 real(dp),intent(in)           :: airtemp       ! air temperature (K)
 real(dp),intent(in)           :: sfcTemp       ! surface temperature (K)
 real(dp),intent(in)           :: windspd       ! wind speed (m s-1)
 real(dp),intent(in)           :: mHeight       ! measurement height (m)
 logical(lgt),intent(in)       :: computeDerivative ! flag to compute the derivative
 ! output
 real(dp),intent(out)          :: RiBulk        ! bulk Richardson number (-)
 real(dp),intent(out)          :: dRiBulk_dTemp ! derivative in the bulk Richardson number w.r.t. temperature (K-1)
 integer(i4b),intent(out)      :: err           ! error code
 character(*),intent(out)      :: message       ! error message
 ! local variables
 real(dp)                      :: T_grad        ! gradient in temperature between the atmosphere and surface (K)
 real(dp)                      :: T_mean        ! mean of the atmosphere and surface temperature (K)
 real(dp)                      :: RiMult        ! dimensionless scaling factor (-)
 ! initialize error control
 err=0; message='bulkRichardson/'
 ! compute local variables
 T_grad = airtemp - sfcTemp
 T_mean = 0.5_dp*(airtemp + sfcTemp)
 RiMult = (gravity*mHeight)/(windspd*windspd)
 ! compute the Richardson number
 RiBulk = (T_grad/T_mean) * RiMult
 ! compute the derivative in the Richardson number
 if(computeDerivative)then
  dRiBulk_dTemp = (-1._dp/T_mean + T_grad*(-0.5_dp)*T_mean**(-2._dp) ) * RiMult
 else
  dRiBulk_dTemp = 1._dp
 endif
 end subroutine bulkRichardson




end module vegNrgFlux_module

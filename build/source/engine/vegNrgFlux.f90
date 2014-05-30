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
! look-up values for the choice of parameterization for vegetation roughness length and displacement height
USE mDecisions_module,only:  &
 Raupach_BLM1994,            & ! Raupach (BLM 1994) "Simplified expressions..."
 CM_QJRMS1998,               & ! Choudhury and Monteith (QJRMS 1998) "A four layer model for the heat budget..."
 vegTypeTable                  ! constant parameters dependent on the vegetation type
! look-up values for the choice of canopy shortwave radiation method
USE mDecisions_module,only:  &
 noah_mp,                    & ! full Noah-MP implementation (including albedo)
 CLM_2stream,                & ! CLM 2-stream model (see CLM documentation)
 UEB_2stream,                & ! UEB 2-stream model (Mahat and Tarboton, WRR 2011)
 NL_scatter,                 & ! Simplified method Nijssen and Lettenmaier (JGR 1999)
 BeersLaw                      ! Beer's Law (as implemented in VIC)
! look-up values for the choice of parameterization for canopy emissivity
USE mDecisions_module,only:  &
 simplExp,                   & ! simple exponential function
 difTrans                      ! parameterized as a function of diffuse transmissivity
! look-up values for the choice of canopy wind profile
USE mDecisions_module,only:  &
 exponential,                & ! exponential wind profile extends to the surface
 logBelowCanopy                ! logarithmic profile below the vegetation canopy
! look-up values for the stomatal resistance formulation
USE mDecisions_module,only:  &
 BallBerry,                  & ! Ball-Berry
 Jarvis,                     & ! Jarvis
 simpleResistance              ! simple resistance formulation
! look-up values for choice of stability function
USE mDecisions_module,only:  &
 standard,                   & ! standard MO similarity, a la Anderson (1976) 
 louisInversePower,          & ! Louis (1979) inverse power function
 mahrtExponential              ! Mahrt (1987) exponential
! look-up values for the choice of groundwater representation (local-column, or single-basin)
USE mDecisions_module,only:  &
 localColumn,                & ! separate groundwater representation in each local soil column
 singleBasin                   ! single groundwater store over the entire basin
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
real(dp),parameter     :: missingValue=-9999._dp  ! missing value, used when diagnostic or state variables are undefined
real(dp),parameter     :: verySmall=1.e-6_dp   ! used as an additive constant to check if substantial difference among real numbers
real(dp),parameter     :: mpe=1.e-6_dp         ! prevents overflow error if division by zero 
real(dp),parameter     :: dx=1.e-6_dp          ! finite difference increment
! control
logical(lgt)           :: printflag            ! flag to turn on printing
contains

 ! ************************************************************************************************
 ! new subroutine: muster program to compute energy fluxes at vegetation and ground surfaces
 ! ************************************************************************************************
 subroutine vegNrgFlux(&
                       ! input
                       dt,                                      & ! intent(in): time step (seconds)
                       iter,                                    & ! intent(in): iteration index
                       firstSubStep,                            & ! intent(in): flag to indicate if we are processing the first sub-step
                       computeVegFlux,                          & ! intent(in): flag to indicate if we need to compute fluxes over vegetation
                       computeShortwave,                        & ! intent(in): flag to indicate if we need to compute shortwave radiation
                       canairTempTrial,                         & ! intent(in): trial value of the canopy air space temperature (K)
                       canopyTempTrial,                         & ! intent(in): trial value of canopy temperature (K)
                       groundTempTrial,                         & ! intent(in): trial value of ground temperature (K)
                       canopyIceTrial,                          & ! intent(in): trial value of mass of ice on the vegetation canopy (kg m-2)
                       canopyLiqTrial,                          & ! intent(in): trial value of mass of liquid water on the vegetation canopy (kg m-2)
                       vegTypeIndex,                            & ! intent(in): vegetation type index
                       soilTypeIndex,                           & ! intent(in): soil type index
                       scalarLAI,                               & ! intent(in): one-sided leaf area index (m2 m-2)
                       scalarSAI,                               & ! intent(in): one-sided stem area index (m2 m-2)
                       scalarExposedLAI,                        & ! intent(in): exposed leaf area index after burial by snow (m2 m-2)
                       scalarExposedSAI,                        & ! intent(in): exposed stem area index after burial by snow (m2 m-2)
                       scalarGrowingSeasonIndex,                & ! intent(in): growing season index (0=off, 1=on)
                       scalarFoliageNitrogenFactor,             & ! intent(in): foliage nitrogen concentration (1.0 = saturated)

                       ! input/output: canopy air space variables
                       scalarVP_CanopyAir,                      & ! intent(inout): trial vapor pressure of the canopy air space (Pa)
                       scalarCanopyStabilityCorrection,         & ! intent(inout): stability correction for the canopy (-)
                       scalarGroundStabilityCorrection,         & ! intent(inout): stability correction for the ground surface (-)

                       ! output: liquid water fluxes associated with evaporation/transpiration
                       scalarCanopyTranspiration,               & ! intent(out): canopy transpiration (kg m-2 s-1)
                       scalarCanopyEvaporation,                 & ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
                       scalarGroundEvaporation,                 & ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)

                       ! output: fluxes
                       canairNetFlux,                           & ! intent(out): net energy flux for the canopy air space (W m-2)
                       canopyNetFlux,                           & ! intent(out): net energy flux for the vegetation canopy (W m-2)
                       groundNetFlux,                           & ! intent(out): net energy flux for the ground surface (W m-2)

                       ! output: flux derivatives
                       dCanairNetFlux_dCanairTemp,              & ! intent(out): derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
                       dCanairNetFlux_dCanopyTemp,              & ! intent(out): derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
                       dCanairNetFlux_dGroundTemp,              & ! intent(out): derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
                       dCanopyNetFlux_dCanairTemp,              & ! intent(out): derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
                       dCanopyNetFlux_dCanopyTemp,              & ! intent(out): derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
                       dCanopyNetFlux_dGroundTemp,              & ! intent(out): derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
                       dGroundNetFlux_dCanairTemp,              & ! intent(out): derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
                       dGroundNetFlux_dCanopyTemp,              & ! intent(out): derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
                       dGroundNetFlux_dGroundTemp,              & ! intent(out): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)

                       ! output: error control
                       err,message)                               ! intent(out): error control
 ! model decisions
 USE data_struc,only:model_decisions                              ! model decision structure
 USE var_lookup,only:iLookDECISIONS                               ! named variables for elements of the decision structure
 ! model variables, parameters, etc.
 USE data_struc,only:time_data,type_data,attr_data,forc_data,mpar_data,mvar_data,bvar_data,indx_data     ! data structures
 USE var_lookup,only:iLookTIME,iLookTYPE,iLookATTR,iLookFORCE,iLookPARAM,iLookMVAR,iLookBVAR,iLookINDEX  ! named variables for structure elements
 ! compute energy and mass fluxes for vegetation
 implicit none
 ! input
 real(dp),intent(in)           :: dt                              ! time step (seconds)
 integer(i4b),intent(in)       :: iter                            ! iteration index
 logical(lgt),intent(in)       :: firstSubStep                    ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(in)       :: computeVegFlux                  ! flag to indicate if computing fluxes over vegetation
 logical(lgt),intent(in)       :: computeShortwave                ! flag to indicate if need to compute shortwave radiation
 real(dp),intent(in)           :: canairTempTrial                 ! trial value of canopy air space temperature (K)
 real(dp),intent(in)           :: canopyTempTrial                 ! trial value of canopy temperature (K)
 real(dp),intent(in)           :: groundTempTrial                 ! trial value of ground temperature (K)
 real(dp),intent(in)           :: canopyIceTrial                  ! trial value of mass of ice on the vegetation canopy (kg m-2)
 real(dp),intent(in)           :: canopyLiqTrial                  ! trial value of mass of liquid water on the vegetation canopy (kg m-2)
 integer(i4b),intent(in)       :: vegTypeIndex                    ! vegetation type index
 integer(i4b),intent(in)       :: soilTypeIndex                   ! soil type index
 real(dp),intent(in)           :: scalarLAI                       ! one-sided leaf area index (m2 m-2)
 real(dp),intent(in)           :: scalarSAI                       ! one-sided stem area index (m2 m-2)
 real(dp),intent(in)           :: scalarExposedLAI                ! exposed leaf area index after burial by snow (m2 m-2)
 real(dp),intent(in)           :: scalarExposedSAI                ! exposed stem area index after burial by snow (m2 m-2)
 real(dp),intent(in)           :: scalarGrowingSeasonIndex        ! growing season index (0=off, 1=on)
 real(dp),intent(in)           :: scalarFoliageNitrogenFactor     ! foliage nitrogen concentration (1.0 = saturated)
 ! input/output: canopy air space variables
 real(dp),intent(inout)        :: scalarVP_CanopyAir              ! trial vapor pressure of the canopy air space (Pa)
 real(dp),intent(inout)        :: scalarCanopyStabilityCorrection ! stability correction for the canopy (-)
 real(dp),intent(inout)        :: scalarGroundStabilityCorrection ! stability correction for the ground surface (-)
 ! output: liquid water fluxes associated with evaporation/transpiration
 real(dp),intent(out)           :: scalarCanopyTranspiration      ! canopy transpiration (kg m-2 s-1)
 real(dp),intent(out)           :: scalarCanopyEvaporation        ! canopy evaporation/condensation (kg m-2 s-1)
 real(dp),intent(out)           :: scalarGroundEvaporation        ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
 ! output: fluxes
 real(dp),intent(out)           :: canairNetFlux                  ! net energy flux for the canopy air space (W m-2)
 real(dp),intent(out)           :: canopyNetFlux                  ! net energy flux for the vegetation canopy (W m-2)
 real(dp),intent(out)           :: groundNetFlux                  ! net energy flux for the ground surface (W m-2)
 ! output: flux derivatives
 real(dp),intent(out)           :: dCanairNetFlux_dCanairTemp     ! derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp),intent(out)           :: dCanairNetFlux_dCanopyTemp     ! derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp),intent(out)           :: dCanairNetFlux_dGroundTemp     ! derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
 real(dp),intent(out)           :: dCanopyNetFlux_dCanairTemp     ! derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp),intent(out)           :: dCanopyNetFlux_dCanopyTemp     ! derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp),intent(out)           :: dCanopyNetFlux_dGroundTemp     ! derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
 real(dp),intent(out)           :: dGroundNetFlux_dCanairTemp     ! derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp),intent(out)           :: dGroundNetFlux_dCanopyTemp     ! derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp),intent(out)           :: dGroundNetFlux_dGroundTemp     ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                             ! error code
 character(*),intent(out)      :: message                         ! error message
 ! local
 character(LEN=256)            :: cmessage                        ! error message of downwind routine
 real(dp)                      :: scalarAquiferStorage            ! aquifer storage (m) -- can be basin-average, or single column
 ! initialize error control
 err=0; message="vegNrgFlux_muster/"

 ! define number of layers
 nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat == ix_soil)  ! number of soil layers
 nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat == ix_snow)  ! number of snow layers
 nLayers = indx_data%var(iLookINDEX%nLayers)%dat(1)                   ! total number of layers

 ! identify the appropriate groundwater variable
 select case(model_decisions(iLookDECISIONS%spatial_gw)%iDecision)
  case(singleBasin); scalarAquiferStorage = bvar_data%var(iLookBVAR%basin__AquiferStorage)%dat(1)
  case(localColumn); scalarAquiferStorage = mvar_data%var(iLookMVAR%scalarAquiferStorage)%dat(1)
  case default; err=20; message=trim(message)//'unable to identify spatial representation of groundwater'; return
 end select ! (modify the groundwater representation for this single-column implementation) 

 ! subroutine to compute energy fluxes at vegetation and ground surfaces
 ! NOTE: separate call to use data structure elements, and ensure variables are used as they are intended (input, input/output, output)
 call vegNrgFlux_muster(&

                        ! input
                        dt,                                                                & ! intent(in): time step (seconds)
                        iter,                                                              & ! intent(in): iteration index
                        firstSubStep,                                                      & ! intent(in): flag to indicate if we are processing the first sub-step
                        computeVegFlux,                                                    & ! intent(in): flag to indicate if we need to compute fluxes over vegetation
                        computeShortwave,                                                  & ! intent(in): flag to indicate if we need to compute shortwave radiation
                        canairTempTrial,                                                   & ! intent(in): trial value of the canopy air space temperature (K)
                        canopyTempTrial,                                                   & ! intent(in): trial value of canopy temperature (K)
                        groundTempTrial,                                                   & ! intent(in): trial value of ground temperature (K)
                        canopyIceTrial,                                                    & ! intent(in): trial value of mass of ice on the vegetation canopy (kg m-2)
                        canopyLiqTrial,                                                    & ! intent(in): trial value of mass of liquid water on the vegetation canopy (kg m-2)
                        vegTypeIndex,                                                      & ! intent(in): vegetation type index
                        soilTypeIndex,                                                     & ! intent(in): soil type index
                        scalarLAI,                                                         & ! intent(in): one-sided leaf area index (m2 m-2)
                        scalarSAI,                                                         & ! intent(in): one-sided stem area index (m2 m-2)
                        scalarExposedLAI,                                                  & ! intent(in): exposed leaf area index after burial by snow (m2 m-2)
                        scalarExposedSAI,                                                  & ! intent(in): exposed stem area index after burial by snow (m2 m-2)
                        scalarGrowingSeasonIndex,                                          & ! intent(in): growing season index (0=off, 1=on)
                        scalarFoliageNitrogenFactor,                                       & ! intent(in): foliage nitrogen concentration (1.0 = saturated)

                        ! model control -- intent(in)
                        model_decisions(iLookDECISIONS%fDerivMeth)%iDecision,              & ! intent(in): choice of method to compute derivatives
                        model_decisions(iLookDECISIONS%veg_traits)%iDecision,              & ! intent(in): choice of parameterization for vegetation roughness length and displacement height
                        model_decisions(iLookDECISIONS%canopySrad)%iDecision,              & ! intent(in): choice of canopy shortwave radiation method
                        model_decisions(iLookDECISIONS%canopyEmis)%iDecision,              & ! intent(in): choice of parameterization for canopy emissivity
                        model_decisions(iLookDECISIONS%windPrfile)%iDecision,              & ! intent(in): choice of canopy wind profile
                        model_decisions(iLookDECISIONS%astability)%iDecision,              & ! intent(in): choice of stability function
                        model_decisions(iLookDECISIONS%soilStress)%iDecision,              & ! intent(in): choice of function for the soil moisture control on stomatal resistance
                        model_decisions(iLookDECISIONS%groundwatr)%iDecision,              & ! intent(in): choice of groundwater parameterization
                        model_decisions(iLookDECISIONS%stomResist)%iDecision,              & ! intent(in): choice of function for stomatal resistance

                        ! model parameters (phenology) -- intent(in)
                        mpar_data%var(iLookPARAM%heightCanopyTop),                         & ! intent(in): height at the top of the vegetation canopy (m)
                        mpar_data%var(iLookPARAM%heightCanopyBottom),                      & ! intent(in): height at the bottom of the vegetation canopy (m)
                        mvar_data%var(iLookMVAR%scalarCanopyIceMax)%dat(1),                & ! intent(in): maximum interception storage capacity for ice (kg m-2)
                        mvar_data%var(iLookMVAR%scalarCanopyLiqMax)%dat(1),                & ! intent(in): maximum interception storage capacity for liquid water (kg m-2)

                        ! model parameters (aerodynamic resistance) -- intent(in)
                        mpar_data%var(iLookPARAM%z0Snow),                                  & ! intent(in): roughness length of snow (m)
                        mpar_data%var(iLookPARAM%z0Soil),                                  & ! intent(in): roughness length of soil (m)
                        mpar_data%var(iLookPARAM%z0Canopy),                                & ! intent(in): roughness length of the canopy (m)
                        mpar_data%var(iLookPARAM%zpdFraction),                             & ! intent(in): zero plane displacement / canopy height (-)
                        mpar_data%var(iLookPARAM%critRichNumber),                          & ! intent(in): critical value for the bulk Richardson number where turbulence ceases (-)
                        mpar_data%var(iLookPARAM%Louis79_bparam),                          & ! intent(in): parameter in Louis (1979) stability function
                        mpar_data%var(iLookPARAM%Louis79_cStar),                           & ! intent(in): parameter in Louis (1979) stability function
                        mpar_data%var(iLookPARAM%Mahrt87_eScale),                          & ! intent(in): exponential scaling factor in the Mahrt (1987) stability function
                        mpar_data%var(iLookPARAM%windReductionParam),                      & ! intent(in): canopy wind reduction parameter (-)                   
                        mpar_data%var(iLookPARAM%leafExchangeCoeff),                       & ! intent(in): turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
                        mpar_data%var(iLookPARAM%leafDimension),                           & ! intent(in): characteristic leaf dimension (m)

                        ! model parameters (soil stress) -- intent(in)
                        mpar_data%var(iLookPARAM%theta_sat),                               & ! intent(in): soil porosity (-)
                        mpar_data%var(iLookPARAM%theta_res),                               & ! intent(in): residual volumetric liquid water content (-)
                        mpar_data%var(iLookPARAM%plantWiltPsi),                            & ! intent(in): matric head at wilting point (m)
                        mpar_data%var(iLookPARAM%soilStressParam),                         & ! intent(in): parameter in the exponential soil stress function (-)
                        mpar_data%var(iLookPARAM%critSoilWilting),                         & ! intent(in): critical vol. liq. water content when plants are wilting (-)
                        mpar_data%var(iLookPARAM%critSoilTranspire),                       & ! intent(in): critical vol. liq. water content when transpiration is limited (-)
                        mpar_data%var(iLookPARAM%critAquiferTranspire),                    & ! intent(in): critical aquifer storage value when transpiration is limited (m)
                        mpar_data%var(iLookPARAM%minStomatalResistance),                   & ! intent(in): mimimum stomatal resistance (s m-1) 
       
                        ! forcing at the upper boundary -- intent(in)
                        attr_data%var(iLookATTR%mHeight),                                  & ! intent(in): measurement height (m)
                        forc_data%var(iLookFORCE%airtemp),                                 & ! intent(in): air temperature at some height above the surface (K)
                        forc_data%var(iLookFORCE%windspd),                                 & ! intent(in): wind speed at some height above the surface (m s-1)
                        forc_data%var(iLookFORCE%airpres),                                 & ! intent(in): air pressure at some height above the surface (Pa)
                        forc_data%var(iLookFORCE%LWRadAtm),                                & ! intent(in): downwelling longwave radiation at the upper boundary (W m-2)
                        mvar_data%var(iLookMVAR%scalarVPair)%dat(1),                       & ! intent(in): vapor pressure at some height above the surface (Pa)
                        mvar_data%var(iLookMVAR%scalarO2air)%dat(1),                       & ! intent(in): atmospheric o2 concentration (Pa)
                        mvar_data%var(iLookMVAR%scalarCO2air)%dat(1),                      & ! intent(in): atmospheric co2 concentration (Pa)
                        mvar_data%var(iLookMVAR%scalarTwetbulb)%dat(1),                    & ! intent(in): wetbulb temperature (K)
                        mvar_data%var(iLookMVAR%scalarRainfall)%dat(1),                    & ! intent(in): computed rainfall rate (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarSnowfall)%dat(1),                    & ! intent(in): computed snowfall rate (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarThroughfallRain)%dat(1),             & ! intent(in): rainfall through the vegetation canopy (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarThroughfallSnow)%dat(1),             & ! intent(in): snowfall through the vegetation canopy (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarCosZenith)%dat(1),                   & ! intent(in): cosine of the solar zenith angle (0-1)
                        mvar_data%var(iLookMVAR%spectralIncomingDirect)%dat(1:nBands),     & ! intent(in): incoming direct solar radiation in each wave band (w m-2)
                        mvar_data%var(iLookMVAR%spectralIncomingDiffuse)%dat(1:nBands),    & ! intent(in): incoming diffuse solar radiation in each wave band (w m-2)

                        ! surface characteristix -- intent(in)
                        mvar_data%var(iLookMVAR%spectralSnowAlbedoDirect)%dat(1:nBands),   & ! intent(in): direct albedo of snow in each spectral band (-)
                        mvar_data%var(iLookMVAR%spectralSnowAlbedoDiffuse)%dat(1:nBands),  & ! intent(in): diffuse albedo of snow in each spectral band (-)
                        mvar_data%var(iLookMVAR%scalarSnowAlbedo)%dat(1),                  & ! intent(inout): snow albedo (-)
                        mvar_data%var(iLookMVAR%scalarSnowAge)%dat(1),                     & ! intent(inout): non-dimensional snow age (-)

                        ! water storage -- intent(in)
                        ! NOTE: soil stress only computed at the start of the substep (iter==1)
                        mvar_data%var(iLookMVAR%scalarSWE)%dat(1),                         & ! intent(in): snow water equivalent on the ground (kg m-2)
                        mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1),                   & ! intent(in): snow depth on the ground surface (m)
                        mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(nSnow+1:nLayers),    & ! intent(in): volumetric fraction of liquid water in each soil layer (-)
                        mvar_data%var(iLookMVAR%mLayerMatricHead)%dat,                     & ! intent(in): matric head in each layer (m)
                        scalarAquiferStorage,                                              & ! intent(in): aquifer storage (m) -- can be local-column or single-basin

                        ! shortwave radiation fluxes -- intent(inout) because only called at the start of the sub-step
                        mvar_data%var(iLookMVAR%scalarCanopyWetFraction)%dat(1),           & ! intent(inout): fraction of canopy that is wet
                        mvar_data%var(iLookMVAR%scalarGroundSnowFraction)%dat(1),          & ! intent(inout): fraction of ground covered with snow (-)
                        mvar_data%var(iLookMVAR%scalarCanopySunlitFraction)%dat(1),        & ! intent(inout): sunlit fraction of canopy (-)
                        mvar_data%var(iLookMVAR%scalarCanopySunlitLAI)%dat(1),             & ! intent(inout): sunlit leaf area (-)
                        mvar_data%var(iLookMVAR%scalarCanopyShadedLAI)%dat(1),             & ! intent(inout): shaded leaf area (-)
                        mvar_data%var(iLookMVAR%scalarCanopySunlitPAR)%dat(1),             & ! intent(inout): average absorbed par for sunlit leaves (w m-2)
                        mvar_data%var(iLookMVAR%scalarCanopyShadedPAR)%dat(1),             & ! intent(inout): average absorbed par for shaded leaves (w m-2)
                        mvar_data%var(iLookMVAR%spectralBelowCanopyDirect)%dat,            & ! intent(inout): downward direct flux below veg layer for each spectral band  W m-2)
                        mvar_data%var(iLookMVAR%spectralBelowCanopyDiffuse)%dat,           & ! intent(inout): downward diffuse flux below veg layer for each spectral band (W m-2)
                        mvar_data%var(iLookMVAR%scalarBelowCanopySolar)%dat(1),            & ! intent(inout): solar radiation transmitted below the canopy (W m-2)
                        mvar_data%var(iLookMVAR%spectralAlbGndDirect)%dat,                 & ! intent(inout): direct  albedo of underlying surface (1:nBands) (-)
                        mvar_data%var(iLookMVAR%spectralAlbGndDiffuse)%dat,                & ! intent(inout): diffuse albedo of underlying surface (1:nBands) (-)
                        mvar_data%var(iLookMVAR%scalarGroundAlbedo)%dat(1),                & ! intent(inout): albedo of the ground surface (-)
                        mvar_data%var(iLookMVAR%scalarCanopyAbsorbedSolar)%dat(1),         & ! intent(inout): solar radiation absorbed by canopy (W m-2)
                        mvar_data%var(iLookMVAR%scalarGroundAbsorbedSolar)%dat(1),         & ! intent(inout): solar radiation absorbed by ground (W m-2)

                        ! longwave radiation fluxes -- intent(out) because called in every step
                        mvar_data%var(iLookMVAR%scalarCanopyEmissivity)%dat(1),            & ! intent(out): effective emissivity of the canopy (-)
                        mvar_data%var(iLookMVAR%scalarLWRadCanopy)%dat(1),                 & ! intent(out): longwave radiation emitted from the canopy (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWRadGround)%dat(1),                 & ! intent(out): longwave radiation emitted at the ground surface (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWRadUbound2Canopy)%dat(1),          & ! intent(out): downward atmospheric longwave radiation absorbed by the canopy (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWRadUbound2Ground)%dat(1),          & ! intent(out): downward atmospheric longwave radiation absorbed by the ground (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWRadUbound2Ubound)%dat(1),          & ! intent(out): atmospheric radiation reflected by the ground and lost thru upper boundary (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWRadCanopy2Ubound)%dat(1),          & ! intent(out): longwave radiation emitted from canopy lost thru upper boundary (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWRadCanopy2Ground)%dat(1),          & ! intent(out): longwave radiation emitted from canopy absorbed by the ground (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWRadCanopy2Canopy)%dat(1),          & ! intent(out): canopy longwave reflected from ground and absorbed by the canopy (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWRadGround2Ubound)%dat(1),          & ! intent(out): longwave radiation emitted from ground lost thru upper boundary (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWRadGround2Canopy)%dat(1),          & ! intent(out): longwave radiation emitted from ground and absorbed by the canopy (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWNetCanopy)%dat(1),                 & ! intent(out): net longwave radiation at the canopy (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWNetGround)%dat(1),                 & ! intent(out): net longwave radiation at the ground surface (W m-2)
                        mvar_data%var(iLookMVAR%scalarLWNetUbound)%dat(1),                 & ! intent(out): net longwave radiation at the upper boundary (W m-2)

                        ! aerodynamic resistance -- intent(out) because called in every step
                        mvar_data%var(iLookMVAR%scalarZ0Canopy)%dat(1),                    & ! intent(out): roughness length of the canopy (m)
                        mvar_data%var(iLookMVAR%scalarWindReductionFactor)%dat(1),         & ! intent(out): canopy wind reduction factor (-)
                        mvar_data%var(iLookMVAR%scalarZeroPlaneDisplacement)%dat(1),       & ! intent(out): zero plane displacement (m) 
                        mvar_data%var(iLookMVAR%scalarRiBulkCanopy)%dat(1),                & ! intent(out): bulk Richardson number for the canopy (-)
                        mvar_data%var(iLookMVAR%scalarRiBulkGround)%dat(1),                & ! intent(out): bulk Richardson number for the ground surface (-)
                        mvar_data%var(iLookMVAR%scalarEddyDiffusCanopyTop)%dat(1),         & ! intent(out): eddy diffusivity for heat at the top of the canopy (m2 s-1)
                        mvar_data%var(iLookMVAR%scalarFrictionVelocity)%dat(1),            & ! intent(out): friction velocity (m s-1)
                        mvar_data%var(iLookMVAR%scalarWindspdCanopyTop)%dat(1),            & ! intent(out): windspeed at the top of the canopy (m s-1)
                        mvar_data%var(iLookMVAR%scalarWindspdCanopyBottom)%dat(1),         & ! intent(out): windspeed at the height of the bottom of the canopy (m s-1)
                        mvar_data%var(iLookMVAR%scalarLeafResistance)%dat(1),              & ! intent(out): mean leaf boundary layer resistance per unit leaf area (s m-1)
                        mvar_data%var(iLookMVAR%scalarGroundResistance)%dat(1),            & ! intent(out): below canopy aerodynamic resistance (s m-1) 
                        mvar_data%var(iLookMVAR%scalarCanopyResistance)%dat(1),            & ! intent(out): above canopy aerodynamic resistance (s m-1)

                        ! soil resistance -- intent(in) and intent(inout) because only called at the first iteration
                        mvar_data%var(iLookMVAR%mLayerRootDensity)%dat,                    & ! intent(in): root density in each layer (-)
                        mvar_data%var(iLookMVAR%scalarAquiferRootFrac)%dat(1),             & ! intent(in): fraction of roots below the lowest soil layer (-)
                        mvar_data%var(iLookMVAR%scalarTranspireLim)%dat(1),                & ! intent(inout): weighted average of the transpiration limiting factor (-)
                        mvar_data%var(iLookMVAR%mLayerTranspireLim)%dat,                   & ! intent(inout): transpiration limiting factor in each layer (-)
                        mvar_data%var(iLookMVAR%scalarTranspireLimAqfr)%dat(1),            & ! intent(inout): transpiration limiting factor for the aquifer (-)
                        mvar_data%var(iLookMVAR%scalarSoilRelHumidity)%dat(1),             & ! intent(inout): relative humidity in the soil pores [0-1]
                        mvar_data%var(iLookMVAR%scalarSoilResistance)%dat(1),              & ! intent(inout): resistance from the soil (s m-1)

                        ! stomatal resistance -- intent(inout) because only called at the first iteration
                        mvar_data%var(iLookMVAR%scalarStomResistSunlit)%dat(1),            & ! intent(inout): stomatal resistance for sunlit leaves (s m-1)
                        mvar_data%var(iLookMVAR%scalarStomResistShaded)%dat(1),            & ! intent(inout): stomatal resistance for shaded leaves (s m-1)
                        mvar_data%var(iLookMVAR%scalarPhotosynthesisSunlit)%dat(1),        & ! intent(inout): sunlit photosynthesis (umolco2 m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarPhotosynthesisShaded)%dat(1),        & ! intent(inout): shaded photosynthesis (umolco2 m-2 s-1)

                        ! turbulent heat fluxes
                        mvar_data%var(iLookMVAR%scalarLatHeatSubVapCanopy)%dat(1),         & ! intent(inout): latent heat of sublimation/vaporization for the vegetation canopy (J kg-1)
                        mvar_data%var(iLookMVAR%scalarLatHeatSubVapGround)%dat(1),         & ! intent(inout): latent heat of sublimation/vaporization for the ground surface (J kg-1)
                        mvar_data%var(iLookMVAR%scalarSatVP_CanopyTemp)%dat(1),            & ! intent(out): saturation vapor pressure at the temperature of the vegetation canopy (Pa)
                        mvar_data%var(iLookMVAR%scalarSatVP_GroundTemp)%dat(1),            & ! intent(out): saturation vapor pressure at the temperature of the ground surface (Pa)
                        mvar_data%var(iLookMVAR%scalarSenHeatTotal)%dat(1),                & ! intent(out): sensible heat from the canopy air space to the atmosphere (W m-2)
                        mvar_data%var(iLookMVAR%scalarSenHeatCanopy)%dat(1),               & ! intent(out): sensible heat flux from the canopy to the canopy air space (W m-2)
                        mvar_data%var(iLookMVAR%scalarSenHeatGround)%dat(1),               & ! intent(out): sensible heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
                        mvar_data%var(iLookMVAR%scalarLatHeatTotal)%dat(1),                & ! intent(out): latent heat from the canopy air space to the atmosphere (W m-2)
                        mvar_data%var(iLookMVAR%scalarLatHeatCanopyEvap)%dat(1),           & ! intent(out): latent heat flux associated with evaporation from the canopy to the canopy air space (W m-2)
                        mvar_data%var(iLookMVAR%scalarLatHeatCanopyTrans)%dat(1),          & ! intent(out): latent heat flux associated with transpiration from the canopy to the canopy air space (W m-2)
                        mvar_data%var(iLookMVAR%scalarLatHeatGround)%dat(1),               & ! intent(out): latent heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)

                        ! advective heat fluxes
                        mvar_data%var(iLookMVAR%scalarCanopyAdvectiveHeatFlux)%dat(1),     & ! intent(out): heat advected to the canopy surface with rain + snow (W m-2)
                        mvar_data%var(iLookMVAR%scalarGroundAdvectiveHeatFlux)%dat(1),     & ! intent(out): heat advected to the ground surface with throughfall (W m-2)

                        ! mass fluxes 
                        mvar_data%var(iLookMVAR%scalarCanopySublimation)%dat(1),           & ! intent(out): canopy sublimation/frost (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarSnowSublimation)%dat(1),             & ! intent(out): snow sublimation/frost -- below canopy or non-vegetated (kg m-2 s-1)

                        ! input/output: canopy air space variables
                        scalarVP_CanopyAir,                                                & ! intent(inout): trial vapor pressure of the canopy air space (Pa)
                        scalarCanopyStabilityCorrection,                                   & ! intent(inout): stability correction for the canopy (-)
                        scalarGroundStabilityCorrection,                                   & ! intent(inout): stability correction for the ground surface (-)

                        ! output: liquid water fluxes
                        scalarCanopyTranspiration,                                         & ! intent(out): canopy transpiration (kg m-2 s-1)
                        scalarCanopyEvaporation,                                           & ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
                        scalarGroundEvaporation,                                           & ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)

                        ! output: fluxes
                        canairNetFlux,                                                     & ! intent(out): net energy flux for the canopy air space (W m-2)
                        canopyNetFlux,                                                     & ! intent(out): net energy flux for the vegetation canopy (W m-2)
                        groundNetFlux,                                                     & ! intent(out): net energy flux for the ground surface (W m-2)

                        ! output: flux derivatives
                        dCanairNetFlux_dCanairTemp,                                        & ! intent(out): derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
                        dCanairNetFlux_dCanopyTemp,                                        & ! intent(out): derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
                        dCanairNetFlux_dGroundTemp,                                        & ! intent(out): derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
                        dCanopyNetFlux_dCanairTemp,                                        & ! intent(out): derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
                        dCanopyNetFlux_dCanopyTemp,                                        & ! intent(out): derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
                        dCanopyNetFlux_dGroundTemp,                                        & ! intent(out): derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
                        dGroundNetFlux_dCanairTemp,                                        & ! intent(out): derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
                        dGroundNetFlux_dCanopyTemp,                                        & ! intent(out): derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
                        dGroundNetFlux_dGroundTemp,                                        & ! intent(out): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)

                        ! output: error control
                        err,cmessage)                                                        ! intent(out): error control
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
                              dt,                                & ! intent(in): time step (seconds)
                              iter,                              & ! intent(in): iteration index
                              firstSubStep,                      & ! intent(in): flag to indicate if we are processing the first sub-step
                              computeVegFlux,                    & ! intent(in): flag to indicate if we need to compute fluxes over vegetation
                              computeShortwave,                  & ! intent(in): flag to indicate if we need to compute shortwave radiation
                              canairTempTrial,                   & ! intent(in): trial value of the canopy air space temperature(K)
                              canopyTempTrial,                   & ! intent(in): trial value of canopy temperature (K)
                              groundTempTrial,                   & ! intent(in): trial value of ground temperature (K)
                              canopyIceTrial,                    & ! intent(in): trial value of mass of ice on the vegetation canopy (kg m-2)
                              canopyLiqTrial,                    & ! intent(in): trial value of mass of liquid water on the vegetation canopy (kg m-2)
                              vegTypeIndex,                      & ! intent(in): vegetation type index
                              soilTypeIndex,                     & ! intent(in): soil type index
                              scalarLAI,                         & ! intent(in): one-sided leaf area index (m2 m-2)
                              scalarSAI,                         & ! intent(in): one-sided stem area index (m2 m-2)
                              scalarExposedLAI,                  & ! intent(in): exposed leaf area index after burial by snow (m2 m-2)
                              scalarExposedSAI,                  & ! intent(in): exposed stem area index after burial by snow (m2 m-2)
                              scalarGrowingSeasonIndex,          & ! intent(in): growing season index (0=off, 1=on)
                              scalarFoliageNitrogenFactor,       & ! intent(in): foliage nitrogen concentration (1.0 = saturated)

                              ! model control -- intent(in)
                              ix_fDerivMeth,                     & ! intent(in): choice of method to compute derivatives
                              ix_veg_traits,                     & ! intent(in): choice of parameterization for vegetation roughness length and displacement height
                              ix_canopySrad,                     & ! intent(in): choice of canopy shortwave radiation method
                              ix_canopyEmis,                     & ! intent(in): choice of parameterization for canopy emissivity
                              ix_windPrfile,                     & ! intent(in): choice of canopy wind profile
                              ix_astability,                     & ! intent(in): choice of stability function
                              ix_soilStress,                     & ! intent(in): choice of function for the soil moisture control on stomatal resistance
                              ix_groundwatr,                     & ! intent(in): groundwater parameterization
                              ix_stomResist,                     & ! intent(in): choice of function for stomatal resistance

                              ! model parameters (phenology) -- intent(in)
                              heightCanopyTop,                   & ! intent(in): height at the top of the vegetation canopy (m)
                              heightCanopyBottom,                & ! intent(in): height at the bottom of the vegetation canopy (m)
                              scalarCanopyIceMax,                & ! intent(in): maximum interception storage capacity for ice (kg m-2)
                              scalarCanopyLiqMax,                & ! intent(in): maximum interception storage capacity for liquid water (kg m-2)

                              ! model parameters (aerodynamic resistance) -- intent(in)
                              z0Snow,                            & ! intent(in): roughness length of snow (m)
                              z0Soil,                            & ! intent(in): roughness length of soil (m)
                              z0CanopyParam,                     & ! intent(in): roughness length of the canopy (m)
                              zpdFraction,                       & ! intent(in): zero plane displacement / canopy height (-)
                              critRichNumber,                    & ! intent(in): critical value for the bulk Richardson number where turbulence ceases (-)
                              Louis79_bparam,                    & ! intent(in): parameter in Louis (1979) stability function
                              Louis79_cStar,                     & ! intent(in): parameter in Louis (1979) stability function
                              Mahrt87_eScale,                    & ! intent(in): exponential scaling factor in the Mahrt (1987) stability function
                              windReductionParam,                & ! intent(in): canopy wind reduction parameter (-)                   
                              leafExchangeCoeff,                 & ! intent(in): turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
                              leafDimension,                     & ! intent(in): characteristic leaf dimension (m)

                              ! model parameters (soil stress) -- intent(in)
                              theta_sat,                         & ! intent(in): soil porosity (-)
                              theta_res,                         & ! intent(in): residual volumetric liquid water content (-)
                              plantWiltPsi,                      & ! intent(in): matric head at wilting point (m)
                              soilStressParam,                   & ! intent(in): parameter in the exponential soil stress function (-)
                              critSoilWilting,                   & ! intent(in): critical vol. liq. water content when plants are wilting (-)
                              critSoilTranspire,                 & ! intent(in): critical vol. liq. water content when transpiration is limited (-)
                              critAquiferTranspire,              & ! intent(in): critical aquifer storage value when transpiration is limited (m)
                              minStomatalResistance,             & ! intent(in): mimimum stomatal resistance (s m-1) 

                              ! forcing at the upper boundary -- intent(in)
                              mHeight,                           & ! intent(in): measurement height (m)
                              airtemp,                           & ! intent(in): air temperature at some height above the surface (K)
                              windspd,                           & ! intent(in): wind speed at some height above the surface (m s-1)
                              airpres,                           & ! intent(in): air pressure at some height above the surface (Pa)
                              LWRadAtm,                          & ! intent(in): downwelling longwave radiation at the upper boundary (W m-2)
                              scalarVPair,                       & ! intent(in): vapor pressure at some height above the surface (Pa)
                              scalarO2air,                       & ! intent(in): atmospheric o2 concentration (Pa)
                              scalarCO2air,                      & ! intent(in): atmospheric co2 concentration (Pa)
                              scalarTwetbulb,                    & ! intent(in): wetbulb temperature (K)
                              scalarRainfall,                    & ! intent(in): computed rainfall rate (kg m-2 s-1)
                              scalarSnowfall,                    & ! intent(in): computed snowfall rate (kg m-2 s-1)
                              scalarThroughfallRain,             & ! intent(in): rainfall through the vegetation canopy (kg m-2 s-1)
                              scalarThroughfallSnow,             & ! intent(in): snowfall through the vegetation canopy (kg m-2 s-1)
                              scalarCosZenith,                   & ! intent(in): cosine of the solar zenith angle (0-1)
                              spectralIncomingDirect,            & ! intent(in): incoming direct solar radiation in each wave band (w m-2)
                              spectralIncomingDiffuse,           & ! intent(in): incoming diffuse solar radiation in each wave band (w m-2)

                              ! surface characteristix -- intent(in)
                              spectralSnowAlbedoDirect,          & ! intent(in): direct albedo of snow in each spectral band (-)
                              spectralSnowAlbedoDiffuse,         & ! intent(in): diffuse albedo of snow in each spectral band (-)
                              scalarSnowAlbedo,                  & ! intent(inout): snow albedo (-)
                              scalarSnowAge,                     & ! intent(inout): non-dimensional snow age (-)

                              ! water storage -- intent(in)
                              ! NOTE: soil stress for transpiration only computed at the start of the substep (iter==1)
                              scalarSWE,                         & ! intent(in): snow water equivalent on the ground (kg m-2)
                              scalarSnowDepth,                   & ! intent(in): snow depth on the ground surface (m)
                              mLayerVolFracLiq,                  & ! intent(in): volumetric fraction of liquid water in each soil layer (-)
                              mLayerMatricHead,                  & ! intent(in): matric head in each layer (m)
                              scalarAquiferStorage,              & ! intent(in): aquifer storage (m)

                              ! shortwave radiation fluxes -- intent(inout) because only called at the start of the sub-step
                              scalarCanopyWetFraction,           & ! intent(inout): fraction of canopy that is wet
                              scalarGroundSnowFraction,          & ! intent(inout): fraction of ground covered with snow (-)
                              scalarCanopySunlitFraction,        & ! intent(inout): sunlit fraction of canopy (-)
                              scalarCanopySunlitLAI,             & ! intent(inout): sunlit leaf area (-)
                              scalarCanopyShadedLAI,             & ! intent(inout): shaded leaf area (-)
                              scalarCanopySunlitPAR,             & ! intent(inout): average absorbed par for sunlit leaves (w m-2)
                              scalarCanopyShadedPAR,             & ! intent(inout): average absorbed par for shaded leaves (w m-2)
                              spectralBelowCanopyDirect,         & ! intent(inout): downward direct flux below veg layer for each spectral band  W m-2)
                              spectralBelowCanopyDiffuse,        & ! intent(inout): downward diffuse flux below veg layer for each spectral band (W m-2)
                              scalarBelowCanopySolar,            & ! intent(inout): radiation transmitted below the canopy (W m-2)
                              spectralAlbGndDirect,              & ! intent(inout): direct  albedo of underlying surface (1:nBands) (-)
                              spectralAlbGndDiffuse,             & ! intent(inout): diffuse albedo of underlying surface (1:nBands) (-)
                              scalarGroundAlbedo,                & ! intent(inout): albedo of the ground surface (-)
                              scalarCanopyAbsorbedSolar,         & ! intent(inout): solar radiation absorbed by canopy (W m-2)
                              scalarGroundAbsorbedSolar,         & ! intent(inout): solar radiation absorbed by ground (W m-2)

                              ! longwave radiation fluxes -- intent(out) because called in every step
                              scalarCanopyEmissivity,            & ! intent(out): effective emissivity of the canopy (-)
                              scalarLWRadCanopy,                 & ! intent(out): longwave radiation emitted from the canopy (W m-2)
                              scalarLWRadGround,                 & ! intent(out): longwave radiation emitted at the ground surface (W m-2)
                              scalarLWRadUbound2Canopy,          & ! intent(out): downward atmospheric longwave radiation absorbed by the canopy (W m-2)
                              scalarLWRadUbound2Ground,          & ! intent(out): downward atmospheric longwave radiation absorbed by the ground (W m-2)
                              scalarLWRadUbound2Ubound,          & ! intent(out): atmospheric radiation reflected by the ground and lost thru upper boundary (W m-2)
                              scalarLWRadCanopy2Ubound,          & ! intent(out): longwave radiation emitted from canopy lost thru upper boundary (W m-2)
                              scalarLWRadCanopy2Ground,          & ! intent(out): longwave radiation emitted from canopy absorbed by the ground (W m-2)
                              scalarLWRadCanopy2Canopy,          & ! intent(out): canopy longwave reflected from ground and absorbed by the canopy (W m-2)
                              scalarLWRadGround2Ubound,          & ! intent(out): longwave radiation emitted from ground lost thru upper boundary (W m-2)
                              scalarLWRadGround2Canopy,          & ! intent(out): longwave radiation emitted from ground and absorbed by the canopy (W m-2)
                              scalarLWNetCanopy,                 & ! intent(out): net longwave radiation at the canopy (W m-2)
                              scalarLWNetGround,                 & ! intent(out): net longwave radiation at the ground surface (W m-2)
                              scalarLWNetUbound,                 & ! intent(out): net longwave radiation at the upper boundary (W m-2)

                              ! aerodynamic resistance -- intent(out) because called in every step
                              scalarZ0Canopy,                    & ! intent(out): roughness length of the canopy (m)
                              scalarWindReductionFactor,         & ! intent(out): canopy wind reduction factor (-)
                              scalarZeroPlaneDisplacement,       & ! intent(out): zero plane displacement (m) 
                              scalarRiBulkCanopy,                & ! intent(out): bulk Richardson number for the canopy (-)
                              scalarRiBulkGround,                & ! intent(out): bulk Richardson number for the ground surface (-)
                              scalarEddyDiffusCanopyTop,         & ! intent(out): eddy diffusivity for heat at the top of the canopy (m2 s-1)
                              scalarFrictionVelocity,            & ! intent(out): friction velocity (m s-1)
                              scalarWindspdCanopyTop,            & ! intent(out): windspeed at the top of the canopy (m s-1)
                              scalarWindspdCanopyBottom,         & ! intent(out): windspeed at the height of the bottom of the canopy (m s-1)
                              scalarLeafResistance,              & ! intent(out): mean leaf boundary layer resistance per unit leaf area (s m-1)
                              scalarGroundResistance,            & ! intent(out): below canopy aerodynamic resistance (s m-1) 
                              scalarCanopyResistance,            & ! intent(out): above canopy aerodynamic resistance (s m-1)

                              ! soil resistance -- intent(in) and intent(inout) because only called at the first iteration
                              mLayerRootDensity,                 & ! intent(in): root density in each layer (-)
                              scalarAquiferRootFrac,             & ! intent(in): fraction of roots below the lowest soil layer (-)
                              scalarTranspireLim,                & ! intent(inout): weighted average of the transpiration limiting factor (-)
                              mLayerTranspireLim,                & ! intent(inout): transpiration limiting factor in each layer (-)
                              scalarTranspireLimAqfr,            & ! intent(inout): transpiration limiting factor for the aquifer (-)
                              scalarSoilRelHumidity,             & ! intent(inout): relative humidity in the soil pores [0-1]
                              scalarSoilResistance,              & ! intent(inout): resistance from the soil (s m-1)

                              ! stomatal resistance -- intent(inout) because only called at the first iteration
                              scalarStomResistSunlit,            & ! intent(inout): stomatal resistance for sunlit leaves (s m-1)
                              scalarStomResistShaded,            & ! intent(inout): stomatal resistance for shaded leaves (s m-1)
                              scalarPhotosynthesisSunlit,        & ! intent(inout): sunlit photosynthesis (umolco2 m-2 s-1)
                              scalarPhotosynthesisShaded,        & ! intent(inout): shaded photosynthesis (umolco2 m-2 s-1)

                              ! turbulent heat fluxes
                              scalarLatHeatSubVapCanopy,         & ! intent(inout): latent heat of sublimation/vaporization for the vegetation canopy (J kg-1)
                              scalarLatHeatSubVapGround,         & ! intent(inout): latent heat of sublimation/vaporization for the ground surface (J kg-1)
                              scalarSatVP_canopyTemp,            & ! intent(out): saturation vapor pressure at the temperature of the vegetation canopy (Pa)
                              scalarSatVP_groundTemp,            & ! intent(out): saturation vapor pressure at the temperature of the ground surface (Pa)
                              scalarSenHeatTotal,                & ! intent(out): sensible heat from the canopy air space to the atmosphere (W m-2)
                              scalarSenHeatCanopy,               & ! intent(out): sensible heat flux from the canopy to the canopy air space (W m-2)
                              scalarSenHeatGround,               & ! intent(out): sensible heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
                              scalarLatHeatTotal,                & ! intent(out): latent heat from the canopy air space to the atmosphere (W m-2)
                              scalarLatHeatCanopyEvap,           & ! intent(out): latent heat flux associated with evaporation from the canopy to the canopy air space (W m-2)
                              scalarLatHeatCanopyTrans,          & ! intent(out): latent heat flux associated with transpiration from the canopy to the canopy air space (W m-2)
                              scalarLatHeatGround,               & ! intent(out): latent heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)

                              ! advective heat flux
                              scalarCanopyAdvectiveHeatFlux,     & ! intent(out): heat advected to the surface with rain + snow (W m-2)
                              scalarGroundAdvectiveHeatFlux,     & ! intent(out): heat advected to the surface with throughfall (W m-2)

                              ! mass fluxes
                              scalarCanopySublimation,           & ! intent(out): canopy sublimation/frost (kg m-2 s-1)
                              scalarSnowSublimation,             & ! intent(out): snow sublimation/frost -- below canopy or non-vegetated (kg m-2 s-1)

                              ! input/output: canopy air space variables
                              scalarVP_CanopyAir,                & ! intent(inout): trial vapor pressure of the canopy air space (Pa)
                              scalarCanopyStabilityCorrection,   & ! intent(inout): stability correction for the canopy (-)
                              scalarGroundStabilityCorrection,   & ! intent(inout): stability correction for the ground surface (-)

                              ! output: liquid water fluxes
                              scalarCanopyTranspiration,         & ! intent(out): canopy transpiration (kg m-2 s-1)
                              scalarCanopyEvaporation,           & ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
                              scalarGroundEvaporation,           & ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)

                              ! output: fluxes
                              canairNetFlux,                     & ! intent(out): net energy flux for the canopy air space (W m-2)
                              canopyNetFlux,                     & ! intent(out): net energy flux for the vegetation canopy (W m-2)
                              groundNetFlux,                     & ! intent(out): net energy flux for the ground surface (W m-2)

                              ! output: flux derivatives
                              dCanairNetFlux_dCanairTemp,        & ! intent(out): derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
                              dCanairNetFlux_dCanopyTemp,        & ! intent(out): derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
                              dCanairNetFlux_dGroundTemp,        & ! intent(out): derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
                              dCanopyNetFlux_dCanairTemp,        & ! intent(out): derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
                              dCanopyNetFlux_dCanopyTemp,        & ! intent(out): derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
                              dCanopyNetFlux_dGroundTemp,        & ! intent(out): derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
                              dGroundNetFlux_dCanairTemp,        & ! intent(out): derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
                              dGroundNetFlux_dCanopyTemp,        & ! intent(out): derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
                              dGroundNetFlux_dGroundTemp,        & ! intent(out): derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)

                              ! output: error control
                              err,message                        ) ! intent(out): error control
 ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! external subroutines
 USE NOAHMP_ROUTINES,only:radiation                                ! subroutine to calculate albedo and shortwave radiaiton in the canopy
 USE vegSWavRad_module,only:vegSWavRad                             ! subroutine to calculate shortwave radiaiton in the canopy
 ! utilities
 USE expIntegral_module,only:expIntegral                           ! subroutine to calculate the exponential integral
 ! conversion functions
 USE conv_funcs_module,only:satVapPress                            ! function to compute the saturated vapor pressure (Pa)
 USE conv_funcs_module,only:getLatentHeatValue                     ! function to identify latent heat of vaporization/sublimation (J kg-1)
 ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! input/output
 ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! input
 real(dp),intent(in)            :: dt                              ! time step (seconds)
 integer(i4b),intent(in)        :: iter                            ! iteration index
 logical(lgt),intent(in)        :: firstSubStep                    ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(in)        :: computeVegFlux                  ! flag to indicate if computing fluxes over vegetation
 logical(lgt),intent(in)        :: computeShortwave                ! flag to indicate if need to compute shortwave radiation
 real(dp),intent(in)            :: canairTempTrial                 ! trial value of canopy air space temperature (K)
 real(dp),intent(in)            :: canopyTempTrial                 ! trial value of canopy temperature (K)
 real(dp),intent(in)            :: groundTempTrial                 ! trial value of ground temperature (K)
 real(dp),intent(in)            :: canopyIceTrial                  ! trial value of mass of ice on the vegetation canopy (kg m-2)
 real(dp),intent(in)            :: canopyLiqTrial                  ! trial value of mass of liquid water on the vegetation canopy (kg m-2)
 integer(i4b),intent(in)        :: vegTypeIndex                    ! vegetation type index
 integer(i4b),intent(in)        :: soilTypeIndex                   ! soil type index
 real(dp),intent(in)            :: scalarLAI                       ! one-sided leaf area index (m2 m-2)
 real(dp),intent(in)            :: scalarSAI                       ! one-sided stem area index (m2 m-2)
 real(dp),intent(in)            :: scalarExposedLAI                ! exposed leaf area index after burial by snow (m2 m-2)
 real(dp),intent(in)            :: scalarExposedSAI                ! exposed stem area index after burial by snow (m2 m-2)
 real(dp),intent(in)            :: scalarGrowingSeasonIndex        ! growing season index (0=off, 1=on)
 real(dp),intent(in)            :: scalarFoliageNitrogenFactor     ! foliage nitrogen concentration (1.0 = saturated)

 ! model control -- intent(in)
 integer(i4b),intent(in)        :: ix_fDerivMeth                   ! choice of method to compute derivatives
 integer(i4b),intent(in)        :: ix_veg_traits                   ! choice of parameterization for vegetation roughness length and displacement height
 integer(i4b),intent(in)        :: ix_canopySrad                   ! choice of canopy shortwave radiation method
 integer(i4b),intent(in)        :: ix_canopyEmis                   ! choice of parameterization for canopy emissivity
 integer(i4b),intent(in)        :: ix_windPrfile                   ! choice of canopy wind profile
 integer(i4b),intent(in)        :: ix_astability                   ! choice of stability function
 integer(i4b),intent(in)        :: ix_soilStress                   ! choice of function for the soil moisture control on stomatal resistance
 integer(i4b),intent(in)        :: ix_groundwatr                   ! groundwater parameterization
 integer(i4b),intent(in)        :: ix_stomResist                   ! choice of function for stomatal resistance

 ! model parameters (phenology) -- intent(in)
 real(dp),intent(in)            :: heightCanopyTop                 ! height at the top of the vegetation canopy (m)
 real(dp),intent(in)            :: heightCanopyBottom              ! height at the bottom of the vegetation canopy (m)
 real(dp),intent(in)            :: scalarCanopyIceMax              ! maximum interception storage capacity for ice (kg m-2)
 real(dp),intent(in)            :: scalarCanopyLiqMax              ! maximum interception storage capacity for liquid water (kg m-2)

 ! model parameters (aerodynamic resistance) -- intent(in)
 real(dp),intent(in)            :: z0Snow                          ! roughness length of snow (m)
 real(dp),intent(in)            :: z0Soil                          ! roughness length of soil (m)
 real(dp),intent(in)            :: z0CanopyParam                   ! roughness length of the canopy (m)
 real(dp),intent(in)            :: zpdFraction                     ! zero plane displacement / canopy height (-)
 real(dp),intent(in)            :: critRichNumber                  ! critical value for the bulk Richardson number where turbulence ceases (-)
 real(dp),intent(in)            :: Louis79_bparam                  ! parameter in Louis (1979) stability function
 real(dp),intent(in)            :: Louis79_cStar                   ! parameter in Louis (1979) stability function
 real(dp),intent(in)            :: Mahrt87_eScale                  ! exponential scaling factor in the Mahrt (1987) stability function
 real(dp),intent(in)            :: windReductionParam              ! canopy wind reduction parameter (-)                   
 real(dp),intent(in)            :: leafExchangeCoeff               ! turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
 real(dp),intent(in)            :: leafDimension                   ! characteristic leaf dimension (m)

 ! model parameters (soil stress) -- intent(in)
 real(dp),intent(in)            :: theta_sat                       ! soil porosity (-)
 real(dp),intent(in)            :: theta_res                       ! residual volumetric liquid water content (-)
 real(dp),intent(in)            :: plantWiltPsi                    ! matric head at wilting point (m)
 real(dp),intent(in)            :: soilStressParam                 ! parameter in the exponential soil stress function (-)
 real(dp),intent(in)            :: critSoilWilting                 ! critical vol. liq. water content when plants are wilting (-)
 real(dp),intent(in)            :: critSoilTranspire               ! critical vol. liq. water content when transpiration is limited (-)
 real(dp),intent(in)            :: critAquiferTranspire            ! critical aquifer storage value when transpiration is limited (m)
 real(dp),intent(in)            :: minStomatalResistance           ! mimimum stomatal resistance (s m-1) 

 ! forcing at the upper boundary -- intent(in)
 real(dp),intent(in)            :: mHeight                         ! measurement height (m)
 real(dp),intent(in)            :: airtemp                         ! air temperature at some height above the surface (K)
 real(dp),intent(in)            :: windspd                         ! wind speed at some height above the surface (m s-1)
 real(dp),intent(in)            :: airpres                         ! air pressure at some height above the surface (Pa)
 real(dp),intent(in)            :: LWRadAtm                        ! downwelling longwave radiation at the upper boundary (W m-2)
 real(dp),intent(in)            :: scalarVPair                     ! vapor pressure at some height above the surface (Pa)
 real(dp),intent(in)            :: scalarO2air                     ! atmospheric o2 concentration (Pa)
 real(dp),intent(in)            :: scalarCO2air                    ! atmospheric co2 concentration (Pa)
 real(dp),intent(in)            :: scalarTwetbulb                  ! wetbulb temperature (K)
 real(dp),intent(in)            :: scalarRainfall                  ! computed rainfall rate (kg m-2 s-1)
 real(dp),intent(in)            :: scalarSnowfall                  ! computed snowfall rate (kg m-2 s-1)
 real(dp),intent(in)            :: scalarThroughfallRain           ! rainfall through the vegetation canopy (kg m-2 s-1)
 real(dp),intent(in)            :: scalarThroughfallSnow           ! snowfall through the vegetation canopy (kg m-2 s-1)
 real(dp),intent(in)            :: scalarCosZenith                 ! cosine of the solar zenith angle (0-1)
 real(dp),intent(in)            :: spectralIncomingDirect(:)       ! incoming direct solar radiation in each wave band (w m-2)
 real(dp),intent(in)            :: spectralIncomingDiffuse(:)      ! incoming diffuse solar radiation in each wave band (w m-2)

 ! surface characteristix -- intent(in)
 real(dp),intent(in)            :: spectralSnowAlbedoDirect(:)     ! direct albedo of snow in each spectral band (-)
 real(dp),intent(in)            :: spectralSnowAlbedoDiffuse(:)    ! diffuse albedo of snow in each spectral band (-)
 real(dp),intent(inout)         :: scalarSnowAlbedo                ! snow albedo (-)
 real(dp),intent(inout)         :: scalarSnowAge                   ! non-dimensional snow age (-)

 ! water storage -- intent(in)
 ! NOTE: soil stress only computed at the start of the substep (iter==1)
 real(dp),intent(in)            :: scalarSWE                       ! snow water equivalent on the ground (kg m-2)
 real(dp),intent(in)            :: scalarSnowDepth                 ! snow depth on the ground surface (m)
 real(dp),intent(in)            :: mLayerVolFracLiq(:)             ! volumetric fraction of liquid water in each soil layer (-)
 real(dp),intent(in)            :: mLayerMatricHead(:)             ! matric head in each layer (m)
 real(dp),intent(in)            :: scalarAquiferStorage            ! aquifer storage (m)

 ! shortwave radiation fluxes -- intent(inout) because only called at the start of the sub-step
 real(dp),intent(inout)         :: scalarCanopyWetFraction         ! fraction of canopy that is wet
 real(dp),intent(inout)         :: scalarGroundSnowFraction        ! fraction of ground covered with snow (-)
 real(dp),intent(inout)         :: scalarCanopySunlitFraction      ! sunlit fraction of canopy (-)
 real(dp),intent(inout)         :: scalarCanopySunlitLAI           ! sunlit leaf area (-)
 real(dp),intent(inout)         :: scalarCanopyShadedLAI           ! shaded leaf area (-)
 real(dp),intent(inout)         :: scalarCanopySunlitPAR           ! average absorbed par for sunlit leaves (w m-2)
 real(dp),intent(inout)         :: scalarCanopyShadedPAR           ! average absorbed par for shaded leaves (w m-2)
 real(dp),intent(inout)         :: spectralBelowCanopyDirect(:)    ! downward direct flux below veg layer for each spectral band  W m-2)
 real(dp),intent(inout)         :: spectralBelowCanopyDiffuse(:)   ! downward diffuse flux below veg layer for each spectral band (W m-2)
 real(dp),intent(inout)         :: scalarBelowCanopySolar          ! solar radiation transmitted below the canopy (W m-2)
 real(dp),intent(inout)         :: spectralAlbGndDirect(:)         ! direct  albedo of underlying surface (1:nBands) (-)
 real(dp),intent(inout)         :: spectralAlbGndDiffuse(:)        ! diffuse albedo of underlying surface (1:nBands) (-)
 real(dp),intent(inout)         :: scalarGroundAlbedo              ! albedo of the ground surface (-)
 real(dp),intent(inout)         :: scalarCanopyAbsorbedSolar       ! solar radiation absorbed by canopy (W m-2)
 real(dp),intent(inout)         :: scalarGroundAbsorbedSolar       ! solar radiation absorbed by ground (W m-2)

 ! longwave radiation fluxes -- intent(out) because called in every step
 real(dp),intent(out)           :: scalarCanopyEmissivity          ! effective emissivity of the canopy (-)
 real(dp),intent(out)           :: scalarLWRadCanopy               ! longwave radiation emitted from the canopy (W m-2)
 real(dp),intent(out)           :: scalarLWRadGround               ! longwave radiation emitted at the ground surface (W m-2)
 real(dp),intent(out)           :: scalarLWRadUbound2Canopy        ! downward atmospheric longwave radiation absorbed by the canopy (W m-2)
 real(dp),intent(out)           :: scalarLWRadUbound2Ground        ! downward atmospheric longwave radiation absorbed by the ground (W m-2)
 real(dp),intent(out)           :: scalarLWRadUbound2Ubound        ! atmospheric radiation reflected by the ground and lost thru upper boundary (W m-2)
 real(dp),intent(out)           :: scalarLWRadCanopy2Ubound        ! longwave radiation emitted from canopy lost thru upper boundary (W m-2)
 real(dp),intent(out)           :: scalarLWRadCanopy2Ground        ! longwave radiation emitted from canopy absorbed by the ground (W m-2)
 real(dp),intent(out)           :: scalarLWRadCanopy2Canopy        ! canopy longwave reflected from ground and absorbed by the canopy (W m-2)
 real(dp),intent(out)           :: scalarLWRadGround2Ubound        ! longwave radiation emitted from ground lost thru upper boundary (W m-2)
 real(dp),intent(out)           :: scalarLWRadGround2Canopy        ! longwave radiation emitted from ground and absorbed by the canopy (W m-2)
 real(dp),intent(out)           :: scalarLWNetCanopy               ! net longwave radiation at the canopy (W m-2)
 real(dp),intent(out)           :: scalarLWNetGround               ! net longwave radiation at the ground surface (W m-2)
 real(dp),intent(out)           :: scalarLWNetUbound               ! net longwave radiation at the upper boundary (W m-2)

 ! aerodynamic resistance -- intent(out) because called in every step
 real(dp),intent(out)           :: scalarZ0Canopy                  ! roughness length of the canopy (m) 
 real(dp),intent(out)           :: scalarWindReductionFactor       ! canopy wind reduction factor (-)
 real(dp),intent(out)           :: scalarZeroPlaneDisplacement     ! zero plane displacement (m) 
 real(dp),intent(out)           :: scalarRiBulkCanopy              ! bulk Richardson number for the canopy (-)
 real(dp),intent(out)           :: scalarRiBulkGround              ! bulk Richardson number for the ground surface (-)
 real(dp),intent(out)           :: scalarEddyDiffusCanopyTop       ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
 real(dp),intent(out)           :: scalarFrictionVelocity          ! friction velocity (m s-1)
 real(dp),intent(out)           :: scalarWindspdCanopyTop          ! windspeed at the top of the canopy (m s-1)
 real(dp),intent(out)           :: scalarWindspdCanopyBottom       ! windspeed at the height of the bottom of the canopy (m s-1)
 real(dp),intent(out)           :: scalarLeafResistance            ! mean leaf boundary layer resistance per unit leaf area (s m-1)
 real(dp),intent(out)           :: scalarGroundResistance          ! below canopy aerodynamic resistance (s m-1) 
 real(dp),intent(out)           :: scalarCanopyResistance          ! above canopy aerodynamic resistance (s m-1)

 ! soil resistance -- intent(in) and intent(inout) because only called at the first iteration
 real(dp),intent(in)            :: mLayerRootDensity(:)            ! root density in each layer (-)
 real(dp),intent(in)            :: scalarAquiferRootFrac           ! fraction of roots below the lowest soil layer (-)
 real(dp),intent(inout)         :: scalarTranspireLim              ! weighted average of the transpiration limiting factor (-)
 real(dp),intent(inout)         :: mLayerTranspireLim(:)           ! transpiration limiting factor in each layer (-)
 real(dp),intent(inout)         :: scalarTranspireLimAqfr          ! transpiration limiting factor for the aquifer (-)
 real(dp),intent(inout)         :: scalarSoilRelHumidity           ! relative humidity in the soil pores [0-1]
 real(dp),intent(inout)         :: scalarSoilResistance            ! resistance from the soil (s m-1)

 ! stomatal resistance -- intent(inout) because only called at the first iteration
 real(dp),intent(inout)         :: scalarStomResistSunlit          ! stomatal resistance for sunlit leaves (s m-1)
 real(dp),intent(inout)         :: scalarStomResistShaded          ! stomatal resistance for shaded leaves (s m-1)
 real(dp),intent(inout)         :: scalarPhotosynthesisSunlit      ! sunlit photosynthesis (umolco2 m-2 s-1)
 real(dp),intent(inout)         :: scalarPhotosynthesisShaded      ! shaded photosynthesis (umolco2 m-2 s-1)

 ! turbulent heat fluxes
 real(dp),intent(inout)         :: scalarLatHeatSubVapCanopy       ! latent heat of sublimation/vaporization for the vegetation canopy (J kg-1)
 real(dp),intent(inout)         :: scalarLatHeatSubVapGround       ! latent heat of sublimation/vaporization for the ground surface (J kg-1)
 real(dp),intent(out)           :: scalarSatVP_canopyTemp          ! saturation vapor pressure at the temperature of the vegetation canopy (Pa)
 real(dp),intent(out)           :: scalarSatVP_groundTemp          ! saturation vapor pressure at the temperature of the ground surface (Pa)
 real(dp),intent(out)           :: scalarSenHeatTotal              ! sensible heat flux from the canopy air space to the atmosphere (W m-2)
 real(dp),intent(out)           :: scalarSenHeatCanopy             ! sensible heat flux from the canopy to the canopy air space (W m-2)
 real(dp),intent(out)           :: scalarSenHeatGround             ! sensible heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
 real(dp),intent(out)           :: scalarLatHeatTotal              ! latent heat flux from the canopy air space to the atmosphere (W m-2)
 real(dp),intent(out)           :: scalarLatHeatCanopyEvap         ! latent heat flux associated with evaporation from the canopy to the canopy air space (W m-2)
 real(dp),intent(out)           :: scalarLatHeatCanopyTrans        ! latent heat flux associated with transpiration from the canopy to the canopy air space (W m-2)
 real(dp),intent(out)           :: scalarLatHeatGround             ! latent heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)

 ! advective heat flux
 real(dp),intent(out)           :: scalarCanopyAdvectiveHeatFlux   ! heat advected to the canopy surface with rain + snow (W m-2)
 real(dp),intent(out)           :: scalarGroundAdvectiveHeatFlux   ! heat advected to the ground surface with throughfall (W m-2)

 ! mass fluxes
 real(dp),intent(out)           :: scalarCanopySublimation         ! canopy sublimation/frost (kg m-2 s-1)
 real(dp),intent(out)           :: scalarSnowSublimation           ! snow sublimation/frost -- below canopy or non-vegetated (kg m-2 s-1)

 ! input/output: canopy air space variables
 real(dp),intent(inout)         :: scalarVP_CanopyAir              ! vapor pressure of the canopy air space (Pa)
 real(dp),intent(inout)         :: scalarCanopyStabilityCorrection ! stability correction for the canopy (-)
 real(dp),intent(inout)         :: scalarGroundStabilityCorrection ! stability correction for the ground surface (-)

 ! output: mass fluxes associated with evaporation/transpiration
 real(dp),intent(out)           :: scalarCanopyTranspiration      ! canopy transpiration (kg m-2 s-1)
 real(dp),intent(out)           :: scalarCanopyEvaporation        ! canopy evaporation/condensation (kg m-2 s-1)
 real(dp),intent(out)           :: scalarGroundEvaporation        ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)

 ! output: fluxes
 real(dp),intent(out)           :: canairNetFlux                   ! net energy flux for the canopy air space (W m-2)
 real(dp),intent(out)           :: canopyNetFlux                   ! net energy flux for the vegetation canopy (W m-2)
 real(dp),intent(out)           :: groundNetFlux                   ! net energy flux for the ground surface (W m-2)

 ! output: flux derivatives
 real(dp),intent(out)           :: dCanairNetFlux_dCanairTemp      ! derivative in net canopy air space flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp),intent(out)           :: dCanairNetFlux_dCanopyTemp      ! derivative in net canopy air space flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp),intent(out)           :: dCanairNetFlux_dGroundTemp      ! derivative in net canopy air space flux w.r.t. ground temperature (W m-2 K-1)
 real(dp),intent(out)           :: dCanopyNetFlux_dCanairTemp      ! derivative in net canopy flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp),intent(out)           :: dCanopyNetFlux_dCanopyTemp      ! derivative in net canopy flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp),intent(out)           :: dCanopyNetFlux_dGroundTemp      ! derivative in net canopy flux w.r.t. ground temperature (W m-2 K-1)
 real(dp),intent(out)           :: dGroundNetFlux_dCanairTemp      ! derivative in net ground flux w.r.t. canopy air temperature (W m-2 K-1)
 real(dp),intent(out)           :: dGroundNetFlux_dCanopyTemp      ! derivative in net ground flux w.r.t. canopy temperature (W m-2 K-1)
 real(dp),intent(out)           :: dGroundNetFlux_dGroundTemp      ! derivative in net ground flux w.r.t. ground temperature (W m-2 K-1)

 ! output: error control
 integer(i4b),intent(out)       :: err                             ! error code
 character(*),intent(out)       :: message                         ! error message

 ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! local (general)
 character(LEN=256)            :: cmessage                         ! error message of downwind routine
 real(dp)                      :: snowmassPlusNewsnow              ! sum of snow mass and new snowfall (kg m-2 [mm])
 real(dp)                      :: VAI                              ! vegetation area index (m2 m-2)
 real(dp)                      :: exposedVAI                       ! exposed vegetation area index (m2 m-2)
 real(dp),parameter            :: scalarVegFraction=1._dp          ! vegetation fraction (=1 forces no canopy gaps and open areas in radiation routine)
 real(dp)                      :: relativeCanopyWater              ! water stored on vegetation canopy, expressed as a fraction of maximum storage (-)
 real(dp)                      :: xArg                             ! argument used in the smoothing function (-)
 real(dp)                      :: smoothFunc                       ! smoothing function used to improve numerical stability at times with limited water storage (-)
 real(dp),parameter            :: smoothThresh=0.01_dp             ! mid-point of the smoothing function (kg m-2 s-1)
 real(dp),parameter            :: smoothScale=0.001_dp             ! scaling factor for the smoothing function (kg m-2 s-1)
 ! local (compute numerical derivatives)
 integer(i4b),parameter        :: unperturbed=1                    ! named variable to identify the case of unperturbed state variables
 integer(i4b),parameter        :: perturbStateGround=2             ! named variable to identify the case where we perturb the ground temperature
 integer(i4b),parameter        :: perturbStateCanopy=3             ! named variable to identify the case where we perturb the canopy temperature
 integer(i4b),parameter        :: perturbStateCanair=4             ! named variable to identify the case where we perturb the canopy air temperature
 integer(i4b)                  :: itry                             ! index of flux evaluation
 integer(i4b)                  :: nFlux                            ! number of flux evaluations
 real(dp)                      :: groundTemp                       ! value of ground temperature used in flux calculations (may be perturbed)
 real(dp)                      :: canopyTemp                       ! value of canopy temperature used in flux calculations (may be perturbed)
 real(dp)                      :: canairTemp                       ! value of canopy air temperature used in flux calculations (may be perturbed)
 real(dp)                      :: try0,try1,try2                   ! trial values to evaluate specific derivatives (testing only)
 ! local (saturation vapor pressure of veg)
 real(dp)                      :: TV_celcius                       ! vegetaion temperature (C)
 real(dp)                      :: TG_celcius                       ! ground temperature (C)
 real(dp)                      :: dSVPCanopy_dCanopyTemp           ! derivative in canopy saturated vapor pressure w.r.t. vegetation temperature (Pa/K)
 real(dp)                      :: dSVPGround_dGroundTemp           ! derivative in ground saturated vapor pressure w.r.t. ground temperature (Pa/K)
 ! local (shortwave radiation)
 real(dp),parameter            :: vegFraction=1._dp                ! vegetation fraction (=1 forces no canopy gaps and open areas in radiation routine)
 real(dp)                      :: scalarTotalReflectedSolar        ! total reflected solar radiation (W m-2)
 real(dp)                      :: scalarTotalAbsorbedSolar         ! total absorbed solar radiation (W m-2)
 real(dp)                      :: scalarCanopyReflectedSolar       ! solar radiation reflected from the canopy (W m-2)
 real(dp)                      :: scalarGroundReflectedSolar       ! solar radiation reflected from the ground (W m-2) 
 real(dp)                      :: scalarBetweenCanopyGapFraction   ! between canopy gap fraction for beam (-)
 real(dp)                      :: scalarWithinCanopyGapFraction    ! within canopy gap fraction for beam (-)
 ! local (longwave radiation)
 real(dp)                      :: expi                             ! exponential integral
 real(dp)                      :: scaleLAI                         ! scaled LAI (computing diffuse transmissivity)
 real(dp)                      :: diffuseTrans                     ! diffuse transmissivity (-)
 real(dp)                      :: groundEmissivity                 ! emissivity of the ground surface (-)
 real(dp),parameter            :: vegEmissivity=0.98_dp            ! emissivity of vegetation (-)
 real(dp),parameter            :: soilEmissivity=0.98_dp           ! emmisivity of the soil (-)
 real(dp),parameter            :: snowEmissivity=0.99_dp           ! emissivity of snow (-)
 real(dp)                      :: dLWNetCanopy_dTCanopy            ! derivative in net canopy radiation w.r.t. canopy temperature (W m-2 K-1)
 real(dp)                      :: dLWNetGround_dTGround            ! derivative in net ground radiation w.r.t. ground temperature (W m-2 K-1)
 real(dp)                      :: dLWNetCanopy_dTGround            ! derivative in net canopy radiation w.r.t. ground temperature (W m-2 K-1)
 real(dp)                      :: dLWNetGround_dTCanopy            ! derivative in net ground radiation w.r.t. canopy temperature (W m-2 K-1)
 ! local (aerodynamic resistance)
 real(dp)                      :: saveTemp_CanopyAir               ! temperature of the canopy air space (K)
 real(dp)                      :: scalarCanopyStabilityCorrection_old    ! stability correction for the canopy (-)
 real(dp)                      :: scalarGroundStabilityCorrection_old    ! stability correction for the ground surface (-)
 ! local (turbulent heat transfer)
 real(dp)                      :: z0Ground                         ! roughness length of the ground (ground below the canopy or non-vegetated surface) (m)
 real(dp)                      :: soilEvapFactor                   ! soil water control on evaporation from non-vegetated surfaces
 real(dp)                      :: soilRelHumidity_noSnow           ! relative humidity in the soil pores [0-1]
 real(dp)                      :: scalarLeafConductance            ! leaf conductance (m s-1)
 real(dp)                      :: scalarCanopyConductance          ! canopy conductance (m s-1)
 real(dp)                      :: scalarGroundConductanceSH        ! ground conductance for sensible heat (m s-1)
 real(dp)                      :: scalarGroundConductanceLH        ! ground conductance for latent heat -- includes soil resistance (m s-1)
 real(dp)                      :: scalarEvapConductance            ! conductance for evaporation (m s-1)
 real(dp)                      :: scalarTransConductance           ! conductance for transpiration (m s-1)
 real(dp)                      :: scalarTotalConductanceSH         ! total conductance for sensible heat (m s-1)
 real(dp)                      :: scalarTotalConductanceLH         ! total conductance for latent heat (m s-1)
 real(dp)                      :: dGroundResistance_dTGround       ! derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
 real(dp)                      :: dGroundResistance_dTCanopy       ! derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
 real(dp)                      :: dGroundResistance_dTCanair       ! derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
 real(dp)                      :: dCanopyResistance_dTCanopy       ! derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
 real(dp)                      :: dCanopyResistance_dTCanair       ! derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)
 real(dp)                      :: turbFluxCanair                   ! total turbulent heat fluxes exchanged at the canopy air space (W m-2)
 real(dp)                      :: turbFluxCanopy                   ! total turbulent heat fluxes from the canopy to the canopy air space (W m-2)
 real(dp)                      :: turbFluxGround                   ! total turbulent heat fluxes from the ground to the canopy air space (W m-2)
 ! local (turbulent heat transfer -- compute numerical derivatives)
 ! (temporary scalar resistances when states are perturbed)
 real(dp)                      :: trialLeafResistance              ! mean leaf boundary layer resistance per unit leaf area (s m-1)
 real(dp)                      :: trialGroundResistance            ! below canopy aerodynamic resistance (s m-1) 
 real(dp)                      :: trialCanopyResistance            ! above canopy aerodynamic resistance (s m-1)
 real(dp)                      :: notUsed_RiBulkCanopy             ! bulk Richardson number for the canopy (-)
 real(dp)                      :: notUsed_RiBulkGround             ! bulk Richardson number for the ground surface (-)
 real(dp)                      :: notUsed_z0Canopy                 ! roughness length of the vegetation canopy (m)
 real(dp)                      :: notUsed_WindReductionFactor      ! canopy wind reduction factor (-)
 real(dp)                      :: notUsed_ZeroPlaneDisplacement    ! zero plane displacement (m) 
 real(dp)                      :: notUsed_scalarCanopyStabilityCorrection  ! stability correction for the canopy (-)
 real(dp)                      :: notUsed_scalarGroundStabilityCorrection  ! stability correction for the ground surface (-)
 real(dp)                      :: notUsed_EddyDiffusCanopyTop      ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
 real(dp)                      :: notUsed_FrictionVelocity         ! friction velocity (m s-1)
 real(dp)                      :: notUsed_WindspdCanopyTop         ! windspeed at the top of the canopy (m s-1)
 real(dp)                      :: notUsed_WindspdCanopyBottom      ! windspeed at the height of the bottom of the canopy (m s-1)
 real(dp)                      :: notUsed_dGroundResistance_dTGround  ! derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
 real(dp)                      :: notUsed_dGroundResistance_dTCanopy  ! derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
 real(dp)                      :: notUsed_dGroundResistance_dTCanair  ! derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
 real(dp)                      :: notUsed_dCanopyResistance_dTCanopy  ! derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
 real(dp)                      :: notUsed_dCanopyResistance_dTCanair  ! derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)
 ! (fluxes after perturbations in model states)
 real(dp)                      :: turbFluxCanair_dStateCanair      ! turbulent exchange from the canopy air space to the atmosphere, after canopy air temperature is perturbed (W m-2)
 real(dp)                      :: turbFluxCanair_dStateCanopy      ! turbulent exchange from the canopy air space to the atmosphere, after canopy temperature is perturbed (W m-2)
 real(dp)                      :: turbFluxCanair_dStateGround      ! turbulent exchange from the canopy air space to the atmosphere, after ground temperature is perturbed (W m-2)
 real(dp)                      :: turbFluxCanopy_dStateCanair      ! total turbulent heat fluxes from the canopy to the canopy air space, after canopy air temperature is perturbed (W m-2)
 real(dp)                      :: turbFluxCanopy_dStateCanopy      ! total turbulent heat fluxes from the canopy to the canopy air space, after canopy temperature is perturbed (W m-2)
 real(dp)                      :: turbFluxCanopy_dStateGround      ! total turbulent heat fluxes from the canopy to the canopy air space, after ground temperature is perturbed (W m-2)
 real(dp)                      :: turbFluxGround_dStateCanair      ! total turbulent heat fluxes from the ground to the canopy air space, after canopy air temperature is perturbed (W m-2)
 real(dp)                      :: turbFluxGround_dStateCanopy      ! total turbulent heat fluxes from the ground to the canopy air space, after canopy temperature is perturbed (W m-2)
 real(dp)                      :: turbFluxGround_dStateGround      ! total turbulent heat fluxes from the ground to the canopy air space, after ground temperature is perturbed (W m-2)
 ! (flux derivatives)
 real(dp)                      :: dTurbFluxCanair_dTCanair         ! derivative in net canopy air space fluxes w.r.t. canopy air temperature (W m-2 K-1)
 real(dp)                      :: dTurbFluxCanair_dTCanopy         ! derivative in net canopy air space fluxes w.r.t. canopy temperature (W m-2 K-1)
 real(dp)                      :: dTurbFluxCanair_dTGround         ! derivative in net canopy air space fluxes w.r.t. ground temperature (W m-2 K-1)
 real(dp)                      :: dTurbFluxCanopy_dTCanair         ! derivative in net canopy turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1) 
 real(dp)                      :: dTurbFluxCanopy_dTCanopy         ! derivative in net canopy turbulent fluxes w.r.t. canopy temperature (W m-2 K-1) 
 real(dp)                      :: dTurbFluxCanopy_dTGround         ! derivative in net canopy turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
 real(dp)                      :: dTurbFluxGround_dTCanair         ! derivative in net ground turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
 real(dp)                      :: dTurbFluxGround_dTCanopy         ! derivative in net ground turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
 real(dp)                      :: dTurbFluxGround_dTGround         ! derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
 ! initialize error control
 err=0; message="vegNrgFlux_muster/"

 ! initialize printflag
 printflag = .false.

 ! set canopy stability corrections to the previous values
 scalarCanopyStabilityCorrection_old = scalarCanopyStabilityCorrection       ! stability correction for the canopy (-)
 scalarGroundStabilityCorrection_old = scalarGroundStabilityCorrection       ! stability correction for the ground surface (-)

 ! initialize variables to compute stomatal resistance
 if(iter==1 .and. firstSubStep)then
  ! vapor pressure in the canopy air space initialized as vapor pressure of air above the vegetation canopy
  ! NOTE: this is needed for the stomatal resistance calculations
  if(scalarVP_CanopyAir < 0._dp)then
   scalarVP_CanopyAir    = scalarVPair - 1._dp    ! "small" offset used to assist in checking initial derivative calculations
  endif
 endif

 ! *******************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************
 ! ***** PRELIMINARIES  **********************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************

 ! set latent heat of sublimation/vaporization for canopy and ground surface (Pa/K)
 ! NOTE: variables are constant over the substep, to simplify relating energy and mass fluxes
 if(iter==1)then
  scalarLatHeatSubVapCanopy = getLatentHeatValue(canopyTempTrial)
  ! case when there is snow on the ground (EXCLUDE "snow without a layer" -- in this case, evaporate from the soil)
  if(nSnow > 0)then
   if(groundTempTrial > Tfreeze)then; err=20; message=trim(message)//'do not expect ground temperature > 0 when snow is on the ground'; return; endif
   scalarLatHeatSubVapGround = LH_sub  ! sublimation from snow
   scalarGroundSnowFraction  = 1._dp
  ! case when the ground is snow-free
  else
   scalarLatHeatSubVapGround = LH_vap  ! evaporation of water in the soil pores: this occurs even if frozen because of super-cooled water
   scalarGroundSnowFraction  = 0._dp
  endif  ! (if there is snow on the ground)
 endif  ! (if the first iteration)
 
 ! compute the roughness length of the ground (ground below the canopy or non-vegetated surface)
 z0Ground = z0soil*(1._dp - scalarGroundSnowFraction) + z0Snow*scalarGroundSnowFraction     ! roughness length (m)

 ! compute the total vegetation area index (leaf plus stem)
 VAI        = scalarLAI + scalarSAI  ! vegetation area index
 exposedVAI = scalarExposedLAI + scalarExposedSAI  !  exposed vegetation area index

 ! compute emissivity of the canopy (-)
 if(computeVegFlux)then
  select case(ix_canopyEmis)
   ! *** simple exponential function
   case(simplExp)
    scalarCanopyEmissivity = 1._dp - exp(-exposedVAI)                                     ! effective emissivity of the canopy (-)
   ! *** canopy emissivity parameterized as a function of diffuse transmissivity
   case(difTrans)
    ! compute the exponential integral
    scaleLAI = 0.5_dp*exposedVAI  
    call expIntegral(1,scaleLAI,expi,err,cmessage)
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    ! compute diffuse transmissivity (-) 
    diffuseTrans = (1._dp - scaleLAI)*exp(-scaleLAI) + (scaleLAI**2._dp)*expi
    ! compute the canopy emissivity
    scalarCanopyEmissivity = (1._dp - diffuseTrans)*vegEmissivity
   ! *** check we found the correct option
   case default
    err=20; message=trim(message)//'unable to identify option for canopy emissivity'; return
  end select
 endif

 ! ensure canopy longwave fluxes are zero when not computing canopy fluxes
 if(.not.computeVegFlux) scalarCanopyEmissivity=0._dp 

 ! compute emissivity of the ground surface (-)
 groundEmissivity = scalarGroundSnowFraction*snowEmissivity + (1._dp - scalarGroundSnowFraction)*soilEmissivity  ! emissivity of the ground surface (-)
 
 ! compute the fraction of canopy that is wet
 ! NOTE: we either sublimate or evaporate over the entire substep
 if(exposedVAI>0._dp .and. computeVegFlux)then
  if(scalarLatHeatSubVapCanopy > LH_vap+verySmall)then ! sublimation
   xarg = (canopyIceTrial - smoothThresh)/smoothScale
   relativeCanopyWater = canopyIceTrial / scalarCanopyIceMax
  else  ! evaporation
   xArg = (canopyLiqTrial - smoothThresh)/smoothScale
   relativeCanopyWater = canopyLiqTrial / max(scalarCanopyLiqMax,10._dp)
  endif
  !print*, 'relativeCanopyWater, canopyLiqTrial, scalarCanopyLiqMax = ', relativeCanopyWater, canopyLiqTrial, scalarCanopyLiqMax
  scalarCanopyWetFraction = min(max(0._dp, relativeCanopyWater), 1._dp)**0.666667_dp
  ! smooth the function near zero, to improve numerical stability
  if(xArg < 50._dp)then
   if(xArg > -50._dp)then
    smoothFunc = 1._dp / (1._dp + exp(-xarg))
   else
    smoothFunc = 0._dp
   endif
   scalarCanopyWetFraction = scalarCanopyWetFraction*smoothFunc
   !if(smoothFunc < 0.7_dp) pause  ' smoothing canopy wet fraction'
  endif
 else
  scalarCanopyWetFraction = 0._dp
 endif


 ! *******************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************
 ! ***** SHORTWAVE RADIATION *****************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************

 ! shortwave radiation is constant over the SUBSTEP
 ! NOTE: canopy radiation does depend on canopy temperature (to distinguish between snow covered canopy, for albedo calculations) and this dependency
 !       is not included on the assumption that albedo changes over the sub-step are small compared to other fluxes
 if(computeShortwave)then

  ! compute the sum of snow mass and new snowfall (kg m-2 [mm])
  snowmassPlusNewsnow = scalarSWE + scalarSnowfall*dt

  select case(ix_canopySrad)

   ! ***** unchanged Noah-MP routine
   case(noah_mp)

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
                   scalarGroundSnowFraction,           & ! intent(in): snow cover fraction (0-1)
                   scalarSnowfall,                     & ! intent(in): snowfall (kg m-2 s-1 [mm/s])
                   scalarCanopyWetFraction,            & ! intent(in): fraction of canopy that is wet
                   scalarExposedLAI,                   & ! intent(in): exposed leaf area index after burial by snow (m2 m-2)
                   scalarExposedSAI,                   & ! intent(in): exposed stem area index after burial by snow (m2 m-2)          
                   mLayerVolFracLiq(1:nSoil),          & ! intent(in): volumetric fraction of liquid water in each soil layer (-)
                   spectralIncomingDirect(1:nBands),   & ! intent(in): incoming direct solar radiation in each wave band (w m-2)
                   spectralIncomingDiffuse(1:nBands),  & ! intent(in): incoming diffuse solar radiation in each wave band (w m-2)
                   vegFraction,                        & ! intent(in): vegetation fraction (=1 forces no canopy gaps and open areas in radiation routine)
                   iLoc, jLoc,                         & ! intent(in): spatial location indices      
                   ! output
                   scalarSnowAlbedo,                   & ! intent(inout): snow albedo (-)
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


   ! **** all other options 
   case(CLM_2stream,UEB_2stream,NL_scatter,BeersLaw)

    call vegSWavRad(&
                    ! input: model control
                    vegTypeIndex,                                       & ! intent(in): index of vegetation type
                    soilTypeIndex,                                      & ! intent(in): index of soil type
                    computeVegFlux,                                     & ! intent(in): logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                    ix_canopySrad,                                      & ! intent(in): index of method used for transmission of shortwave rad through the canopy
                    ! input: model variables
                    scalarCosZenith,                                    & ! intent(in): cosine of direct zenith angle (0-1)
                    spectralIncomingDirect(1:nBands),                   & ! intent(in): incoming direct solar radiation in each wave band (w m-2)
                    spectralIncomingDiffuse(1:nBands),                  & ! intent(in): incoming diffuse solar radiation in each wave band (w m-2)
                    spectralSnowAlbedoDirect(1:nBands),                 & ! intent(in): direct albedo of snow in each spectral band (-)
                    spectralSnowAlbedoDiffuse(1:nBands),                & ! intent(in): diffuse albedo of snow in each spectral band (-)
                    scalarExposedLAI,                                   & ! intent(in): exposed leaf area index after burial by snow (m2 m-2)
                    scalarExposedSAI,                                   & ! intent(in): exposed stem area index after burial by snow (m2 m-2)
                    scalarVegFraction,                                  & ! intent(in): vegetation fraction (=1 forces no canopy gaps and open areas in radiation routine)
                    scalarCanopyWetFraction,                            & ! intent(in): fraction of lai, sai that is wetted (-)
                    scalarGroundSnowFraction,                           & ! intent(in): fraction of ground that is snow covered (-)
                    mLayerVolFracLiq(1),                                & ! intent(in): volumetric liquid water content in the upper-most soil layer (-)
                    canopyTempTrial,                                    & ! intent(in): canopy temperature (k)
                    ! output
                    spectralBelowCanopyDirect,                          & ! intent(out): downward direct flux below veg layer for each spectral band  W m-2)
                    spectralBelowCanopyDiffuse,                         & ! intent(out): downward diffuse flux below veg layer for each spectral band (W m-2)
                    scalarBelowCanopySolar,                             & ! intent(out): solar radiation transmitted below the canopy (W m-2)
                    spectralAlbGndDirect,                               & ! intent(out): direct  albedo of underlying surface (1:nBands) (-)
                    spectralAlbGndDiffuse,                              & ! intent(out): diffuse albedo of underlying surface (1:nBands) (-)
                    scalarGroundAlbedo,                                 & ! intent(out): albedo of the ground surface (-)
                    scalarCanopyAbsorbedSolar,                          & ! intent(out): solar radiation absorbed by the vegetation canopy (W m-2)
                    scalarGroundAbsorbedSolar,                          & ! intent(out): solar radiation absorbed by the ground (W m-2)
                    scalarCanopySunlitFraction,                         & ! intent(out): sunlit fraction of canopy (-)
                    scalarCanopySunlitLAI,                              & ! intent(out): sunlit leaf area (-)
                    scalarCanopyShadedLAI,                              & ! intent(out): shaded leaf area (-)
                    scalarCanopySunlitPAR,                              & ! intent(out): average absorbed par for sunlit leaves (w m-2)
                    scalarCanopyShadedPAR,                              & ! intent(out): average absorbed par for shaded leaves (w m-2)
                    err,cmessage)                                         ! intent(out): error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   case default; err=20; message=trim(message)//'unable to identify option for canopy sw radiation'; return

  end select ! (option for canopy sw radiation)


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
                 iter,                               & ! intent(in): iteration index
                 computeVegFlux,                     & ! intent(in): logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                 (ix_fDerivMeth == analytical),      & ! intent(in): logical flag if would like to compute analytical derivaties
                 ix_veg_traits,                      & ! intent(in): choice of parameterization for vegetation roughness length and displacement height
                 ix_windPrfile,                      & ! intent(in): choice of canopy wind profile
                 ix_astability,                      & ! intent(in): choice of stability function
                 ! input: above-canopy forcing data
                 mHeight,                            & ! intent(in): measurement height (m)
                 airtemp,                            & ! intent(in): air temperature at some height above the surface (K)
                 windspd,                            & ! intent(in): wind speed at some height above the surface (m s-1)
                 ! input: canopy and ground temperature
                 canairTempTrial,                    & ! intent(in): temperature of the canopy air space (K)
                 canopyTempTrial,                    & ! intent(in): temperature of the vegetation canopy (K)
                 groundTempTrial,                    & ! intent(in): temperature of the ground surface (K)
                 ! input: diagnostic variables
                 exposedVAI,                         & ! intent(in): exposed vegetation area index -- leaf plus stem (m2 m-2)
                 scalarSnowDepth,                    & ! intent(in): snow depth (m)
                 ! input: parameters
                 z0Ground,                           & ! intent(in): roughness length of the ground (below canopy or non-vegetated surface [snow]) (m)
                 z0CanopyParam,                      & ! intent(in): roughness length of the canopy (m)
                 zpdFraction,                        & ! intent(in): zero plane displacement / canopy height (-)
                 critRichNumber,                     & ! intent(in): critical value for the bulk Richardson number where turbulence ceases (-)
                 Louis79_bparam,                     & ! intent(in): parameter in Louis (1979) stability function
                 Louis79_cStar,                      & ! intent(in): parameter in Louis (1979) stability function
                 Mahrt87_eScale,                     & ! intent(in): exponential scaling factor in the Mahrt (1987) stability function
                 windReductionParam,                 & ! intent(in): canopy wind reduction parameter (-)                   
                 leafExchangeCoeff,                  & ! intent(in): turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
                 leafDimension,                      & ! intent(in): characteristic leaf dimension (m)
                 heightCanopyTop,                    & ! intent(in): height at the top of the vegetation canopy (m) 
                 heightCanopyBottom,                 & ! intent(in): height at the bottom of the vegetation canopy (m) 
                 ! input: stability correction from the last iteration
                 scalarCanopyStabilityCorrection_old,& ! intent(in): stability correction for the canopy (-)
                 scalarGroundStabilityCorrection_old,& ! intent(in): stability correction for the ground surface (-)
                 ! output: stability corrections
                 scalarRiBulkCanopy,                 & ! intent(out): bulk Richardson number for the canopy (-)
                 scalarRiBulkGround,                 & ! intent(out): bulk Richardson number for the ground surface (-)
                 scalarCanopyStabilityCorrection,    & ! intent(out): stability correction for the canopy (-)
                 scalarGroundStabilityCorrection,    & ! intent(out): stability correction for the ground surface (-)
                 ! output: scalar resistances
                 scalarZ0Canopy,                     & ! intent(out): roughness length of the canopy (m)
                 scalarWindReductionFactor,          & ! intent(out): canopy wind reduction factor (-)
                 scalarZeroPlaneDisplacement,        & ! intent(out): zero plane displacement (m) 
                 scalarEddyDiffusCanopyTop,          & ! intent(out): eddy diffusivity for heat at the top of the canopy (m2 s-1)
                 scalarFrictionVelocity,             & ! intent(out): friction velocity (m s-1)
                 scalarWindspdCanopyTop,             & ! intent(out): windspeed at the top of the canopy (m s-1)
                 scalarWindspdCanopyBottom,          & ! intent(out): windspeed at the height of the bottom of the canopy (m s-1)
                 scalarLeafResistance,               & ! intent(out): mean leaf boundary layer resistance per unit leaf area (s m-1)
                 scalarGroundResistance,             & ! intent(out): below canopy aerodynamic resistance (s m-1) 
                 scalarCanopyResistance,             & ! intent(out): above canopy aerodynamic resistance (s m-1)
                 ! output: derivatives in scalar resistances
                 dGroundResistance_dTGround,         & ! intent(out): derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
                 dGroundResistance_dTCanopy,         & ! intent(out): derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
                 dGroundResistance_dTCanair,         & ! intent(out): derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
                 dCanopyResistance_dTCanopy,         & ! intent(out): derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
                 dCanopyResistance_dTCanair,         & ! intent(out): derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)
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

  !write(*,'(a,1x,10(f20.10,1x))') 'canopyTempTrial, scalarSatVP_CanopyTemp, scalarVP_CanopyAir = ', &
  !                                 canopyTempTrial, scalarSatVP_CanopyTemp, scalarVP_CanopyAir

  ! compute stomatal resistance
  select case(ix_stomResist)

   case(simpleResistance)
    ! check that we don't divide by zero -- should be set to minimum of tiny in runroutine soilResist
    if(scalarTranspireLim < tiny(dt))then; err=20; message=trim(message)//'soil moisture stress factor is < tiny -- this will cause problems'; return; endif
    ! compute stomatal resistance (assume equal for sunlit and shaded leaves)
    scalarStomResistSunlit = minStomatalResistance/scalarTranspireLim
    scalarStomResistShaded = scalarStomResistSunlit
    ! set photosynthesis to missing (not computed)
    scalarPhotosynthesisSunlit = missingValue
    scalarPhotosynthesisShaded = missingValue 

   ! compute stomatal resistance (wrapper around the Noah-MP routines)
   ! NOTE: canopy air vapor pressure is from the previous time step
   case(BallBerry,Jarvis)
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
                    err,cmessage                       ) ! intent(out): error control
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

   ! error check
   case default; err=20; message=trim(message)//'unable to identify option for stomatal resistance'; return

  endselect  ! (identifying option for stomatal resistance)

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
                  scalarCanopyEmissivity,            & ! intent(in): canopy emissivity (-)
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
  nFlux=4  ! compute the derivatives using one-sided finite differences
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
    groundTemp = groundTempTrial
    canopyTemp = canopyTempTrial 
    canairTemp = canairTempTrial

   ! perturb ground temperature
   case(perturbStateGround)
    groundTemp = groundTempTrial + dx
    canopyTemp = canopyTempTrial
    canairTemp = canairTempTrial

   ! perturb canopy temperature
   case(perturbStateCanopy)
    groundTemp = groundTempTrial
    canopyTemp = canopyTempTrial + dx
    canairTemp = canairTempTrial

   ! perturb canopy air temperature
   case(perturbStateCanair)
    groundTemp = groundTempTrial
    canopyTemp = canopyTempTrial
    canairTemp = canairTempTrial + dx

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
                   iter,                                    & ! intent(in): iteration index
                   computeVegFlux,                          & ! intent(in): logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                   .false.,                                 & ! intent(in): logical flag if would like to compute analytical derivaties
                   ix_veg_traits,                           & ! intent(in): choice of parameterization for vegetation roughness length and displacement height
                   ix_windPrfile,                           & ! intent(in): choice of canopy wind profile
                   ix_astability,                           & ! intent(in): choice of stability function
                   ! input: above-canopy forcing data
                   mHeight,                                 & ! intent(in): measurement height (m)
                   airtemp,                                 & ! intent(in): air temperature at some height above the surface (K)
                   windspd,                                 & ! intent(in): wind speed at some height above the surface (m s-1)
                   ! input: temperature (canopy, ground, canopy air space)
                   canairTemp,                              & ! intent(in): temperature of the canopy air space (K)
                   canopyTemp,                              & ! intent(in): canopy temperature (K)
                   groundTemp,                              & ! intent(in): ground temperature (K)
                   ! input: diagnostic variables
                   exposedVAI,                              & ! intent(in): exposed vegetation area index -- leaf plus stem (m2 m-2)
                   scalarSnowDepth,                         & ! intent(in): snow depth (m)
                   ! input: parameters
                   z0Ground,                                & ! intent(in): roughness length of the ground (below canopy or non-vegetated surface [snow]) (m)
                   z0CanopyParam,                           & ! intent(in): roughness length of the canopy (m)
                   zpdFraction,                             & ! intent(in): zero plane displacement / canopy height (-)
                   critRichNumber,                          & ! intent(in): critical value for the bulk Richardson number where turbulence ceases (-)
                   Louis79_bparam,                          & ! intent(in): parameter in Louis (1979) stability function
                   Louis79_cStar,                           & ! intent(in): parameter in Louis (1979) stability function
                   Mahrt87_eScale,                          & ! intent(in): exponential scaling factor in the Mahrt (1987) stability function
                   windReductionParam,                      & ! intent(in): canopy wind reduction parameter (-)                   
                   leafExchangeCoeff,                       & ! intent(in): turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
                   leafDimension,                           & ! intent(in): characteristic leaf dimension (m)
                   heightCanopyTop,                         & ! intent(in): height at the top of the vegetation canopy (m) 
                   heightCanopyBottom,                      & ! intent(in): height at the bottom of the vegetation canopy (m) 
                   ! input: stability correction from the last iteration
                   scalarCanopyStabilityCorrection_old,     & ! intent(in): stability correction for the canopy (-)
                   scalarGroundStabilityCorrection_old,     & ! intent(in): stability correction for the ground surface (-)
                   ! output: stability corrections
                   notUsed_RiBulkCanopy,                    & ! intent(out): bulk Richardson number for the canopy (-)
                   notUsed_RiBulkGround,                    & ! intent(out): bulk Richardson number for the ground surface (-)
                   notUsed_scalarCanopyStabilityCorrection, & ! intent(out): stability correction for the canopy (-)
                   notUsed_scalarGroundStabilityCorrection, & ! intent(out): stability correction for the ground surface (-)
                   ! output: scalar resistances
                   notUsed_z0Canopy,                        & ! intent(out): roughness length of the canopy (m)
                   notUsed_WindReductionFactor,             & ! intent(out): canopy wind reduction factor (-)
                   notUsed_ZeroPlaneDisplacement,           & ! intent(out): zero plane displacement (m) 
                   notUsed_EddyDiffusCanopyTop,             & ! intent(out): eddy diffusivity for heat at the top of the canopy (m2 s-1)
                   notUsed_FrictionVelocity,                & ! intent(out): friction velocity (m s-1)
                   notUsed_WindspdCanopyTop,                & ! intent(out): windspeed at the top of the canopy (m s-1)
                   notUsed_WindspdCanopyBottom,             & ! intent(out): windspeed at the height of the bottom of the canopy (m s-1)
                   trialLeafResistance,                     & ! intent(out): mean leaf boundary layer resistance per unit leaf area (s m-1)
                   trialGroundResistance,                   & ! intent(out): below canopy aerodynamic resistance (s m-1) 
                   trialCanopyResistance,                   & ! intent(out): above canopy aerodynamic resistance (s m-1)
                   ! output: derivatives in scalar resistances
                   notUsed_dGroundResistance_dTGround,      & ! intent(out): derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
                   notUsed_dGroundResistance_dTCanopy,      & ! intent(out): derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
                   notUsed_dGroundResistance_dTCanair,      & ! intent(out): derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
                   notUsed_dCanopyResistance_dTCanopy,      & ! intent(out): derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
                   notUsed_dCanopyResistance_dTCanair,      & ! intent(out): derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)
                   ! output: error control
                   err,cmessage                             ) ! intent(out): error control
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
  ! NOTE: computations are based on start-of-step values, so only compute at the first iteration
  if(iter == 1)then
   ! (soil water evaporation factor [0-1])
   soilEvapFactor = mLayerVolFracLiq(1)/(theta_sat - theta_res)
   ! (resistance from the soil [s m-1])
   !scalarSoilResistance = scalarGroundSnowFraction*1._dp + (1._dp - scalarGroundSnowFraction)*EXP(8.25_dp - 4.225_dp*soilEvapFactor)  ! Sellers (1992)
   scalarSoilResistance = scalarGroundSnowFraction*0._dp + (1._dp - scalarGroundSnowFraction)*exp(8.25_dp - 6.0_dp*soilEvapFactor)    ! Niu adjustment to decrease resitance for wet soil
   ! (relative humidity in the soil pores [0-1])
   if(mLayerMatricHead(1) > -1.e+5_dp)then  ! avoid problems with numerical precision when soil is very dry
    soilRelHumidity_noSnow = exp( (mLayerMatricHead(1)*gravity) / (groundTemp*R_wv) )
   else
    soilRelHumidity_noSnow = 0._dp
   endif ! (if matric head is very low)
   scalarSoilRelHumidity  = scalarGroundSnowFraction*1._dp + (1._dp - scalarGroundSnowFraction)*soilRelHumidity_noSnow
   !print*, 'mLayerMatricHead(1), scalarSoilRelHumidity = ', mLayerMatricHead(1), scalarSoilRelHumidity
  endif  ! (if the first iteration)

  ! compute turbulent heat fluxes
  call turbFluxes(&
                  ! input: model control
                  dt,                                   & ! intent(in): model time step (seconds)
                  computeVegFlux,                       & ! intent(in): logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                  ix_fDerivMeth,                        & ! intent(in): method used to calculate flux derivatives
                  ! input: above-canopy forcing data
                  airtemp,                              & ! intent(in): air temperature at some height above the surface (K)
                  airpres,                              & ! intent(in): air pressure of the air above the vegetation canopy (Pa)
                  scalarVPair,                          & ! intent(in): vapor pressure of the air above the vegetation canopy (Pa)
                  ! input: latent heat of sublimation/vaporization
                  scalarLatHeatSubVapCanopy,            & ! intent(in): latent heat of sublimation/vaporization for the vegetation canopy (J kg-1)
                  scalarLatHeatSubVapGround,            & ! intent(in): latent heat of sublimation/vaporization for the ground surface (J kg-1)
                  ! input: canopy liquid and canopy ice (used as a solution constraint)
                  canopyLiqTrial,                       & ! intent(in): mass of liquid water on the vegetation canopy (kg m-2)
                  canopyIceTrial,                       & ! intent(in): mass of ice water on the vegetation canopy (kg m-2)
                  ! input: canopy/ground temperature and saturated vapor pressure
                  canairTemp,                           & ! intent(in): temperature of the canopy air space (K)
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
                  dGroundResistance_dTGround,           & ! intent(in): derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
                  dGroundResistance_dTCanopy,           & ! intent(in): derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
                  dGroundResistance_dTCanair,           & ! intent(in): derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
                  dCanopyResistance_dTCanopy,           & ! intent(in): derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
                  dCanopyResistance_dTCanair,           & ! intent(in): derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)
                  ! output: conductances (used to check derivative calculations)
                  scalarLeafConductance,                & ! intent(out): leaf conductance (m s-1)
                  scalarCanopyConductance,              & ! intent(out): canopy conductance (m s-1)
                  scalarGroundConductanceSH,            & ! intent(out): ground conductance for sensible heat (m s-1)
                  scalarGroundConductanceLH,            & ! intent(out): ground conductance for latent heat -- includes soil resistance (m s-1)
                  scalarEvapConductance,                & ! intent(out): conductance for evaporation (m s-1)
                  scalarTransConductance,               & ! intent(out): conductance for transpiration (m s-1)
                  scalarTotalConductanceSH,             & ! intent(out): total conductance for sensible heat (m s-1)
                  scalarTotalConductanceLH,             & ! intent(out): total conductance for latent heat (m s-1)
                  ! output: canopy air space variables
                  scalarVP_CanopyAir,                   & ! intent(out): vapor pressure of the canopy air space (Pa)
                  ! output: fluxes from the vegetation canopy
                  scalarSenHeatCanopy,                  & ! intent(out): sensible heat flux from the canopy to the canopy air space (W m-2)
                  scalarLatHeatCanopyEvap,              & ! intent(out): latent heat flux associated with evaporation from the canopy to the canopy air space (W m-2)
                  scalarLatHeatCanopyTrans,             & ! intent(out): latent heat flux associated with transpiration from the canopy to the canopy air space (W m-2)
                  ! output: fluxes from non-vegetated surfaces (ground surface below vegetation, bare ground, or snow covered vegetation)
                  scalarSenHeatGround,                  & ! intent(out): sensible heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
                  scalarLatHeatGround,                  & ! intent(out): latent heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
                  ! output: total heat fluxes to the atmosphere
                  scalarSenHeatTotal,                   & ! intent(out): total sensible heat flux to the atmosphere (W m-2)
                  scalarLatHeatTotal,                   & ! intent(out): total latent heat flux to the atmosphere (W m-2)
                  ! output: net fluxes
                  turbFluxCanair,                       & ! intent(out): net turbulent heat fluxes at the canopy air space (W m-2)
                  turbFluxCanopy,                       & ! intent(out): net turbulent heat fluxes at the canopy (W m-2)
                  turbFluxGround,                       & ! intent(out): net turbulent heat fluxes at the ground surface (W m-2)
                  ! output: flux derivatives
                  dTurbFluxCanair_dTCanair,             & ! intent(out): derivative in net canopy air space fluxes w.r.t. canopy air temperature (W m-2 K-1)
                  dTurbFluxCanair_dTCanopy,             & ! intent(out): derivative in net canopy air space fluxes w.r.t. canopy temperature (W m-2 K-1)
                  dTurbFluxCanair_dTGround,             & ! intent(out): derivative in net canopy air space fluxes w.r.t. ground temperature (W m-2 K-1)
                  dTurbFluxCanopy_dTCanair,             & ! intent(out): derivative in net canopy turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
                  dTurbFluxCanopy_dTCanopy,             & ! intent(out): derivative in net canopy turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
                  dTurbFluxCanopy_dTGround,             & ! intent(out): derivative in net canopy turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
                  dTurbFluxGround_dTCanair,             & ! intent(out): derivative in net ground turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
                  dTurbFluxGround_dTCanopy,             & ! intent(out): derivative in net ground turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
                  dTurbFluxGround_dTGround,             & ! intent(out): derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
                  ! output: error control
                  err,cmessage                          ) ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  !print*, 'scalarSenHeatTotal = ', scalarSenHeatTotal
  !print*, 'scalarSenHeatCanopy = ', scalarSenHeatCanopy
  !print*, 'scalarLatHeatCanopyEvap = ', scalarLatHeatCanopyEvap
  !print*, 'scalarLatHeatCanopyTrans = ', scalarLatHeatCanopyTrans

  !print*, 'scalarSenHeatGround = ', scalarSenHeatGround
  !print*, 'scalarLatHeatGround = ', scalarLatHeatGround

  !notUsed_scalarCanopyStabilityCorrection  ! stability correction for the canopy (-)
  !notUsed_scalarGroundStabilityCorrection  ! stability correction for the ground surface (-)
  !notUsed_EddyDiffusCanopyTop              ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
  !notUsed_FrictionVelocity                 ! friction velocity (m s-1)
  !notUsed_WindspdCanopyTop                 ! windspeed at the top of the canopy (m s-1)
  !notUsed_WindspdCanopyBottom              ! windspeed at the height of the bottom of the canopy (m s-1)
  !trialLeafResistance                      ! mean leaf boundary layer resistance per unit leaf area (s m-1)
  !trialGroundResistance                    ! below canopy aerodynamic resistance (s m-1) 
  !trialCanopyResistance                    ! above canopy aerodynamic resistance (s m-1)

  ! save perturbed fluxes
  if(ix_fDerivMeth == numerical)then
   select case(itry) ! (select type of perturbation)
    case(unperturbed)
     !try0 = scalarSenHeatCanopy
     exit
    case(perturbStateCanair)
     !try1 = scalarSenHeatCanopy
     turbFluxCanair_dStateCanair = turbFluxCanair         ! turbulent exchange from the canopy air space to the atmosphere (W m-2)
     turbFluxCanopy_dStateCanair = turbFluxCanopy         ! total turbulent heat fluxes from the canopy to the canopy air space (W m-2)
     turbFluxGround_dStateCanair = turbFluxGround         ! total turbulent heat fluxes from the ground to the canopy air space (W m-2)
    case(perturbStateCanopy)
     turbFluxCanair_dStateCanopy = turbFluxCanair         ! turbulent exchange from the canopy air space to the atmosphere (W m-2)
     turbFluxCanopy_dStateCanopy = turbFluxCanopy         ! total turbulent heat fluxes from the canopy to the canopy air space (W m-2)
     turbFluxGround_dStateCanopy = turbFluxGround         ! total turbulent heat fluxes from the ground to the canopy air space (W m-2)
    case(perturbStateGround)
     turbFluxCanair_dStateGround = turbFluxCanair         ! turbulent exchange from the canopy air space to the atmosphere (W m-2)
     turbFluxCanopy_dStateGround = turbFluxCanopy         ! total turbulent heat fluxes from the canopy to the canopy air space (W m-2)
     turbFluxGround_dStateGround = turbFluxGround         ! total turbulent heat fluxes from the ground to the canopy air space (W m-2)
    case default; err=10; message=trim(message)//"unknown perturbation"; return
   end select ! (type of perturbation)
  endif ! (if numerical)

 end do  ! (looping through different flux perturbations)

 ! test derivative
 !if(ix_fDerivMeth == numerical)  print*, 'derivative = ', (ix_fDerivMeth == numerical), (try1 - try0)/dx
 !if(ix_fDerivMeth == analytical) print*, 'derivative = ', (ix_fDerivMeth == numerical), dTurbFluxCanair_dTCanopy 
 !pause

 ! compute numerical derivatives
 if(ix_fDerivMeth == numerical)then
  ! derivatives w.r.t. canopy air temperature
  dTurbFluxCanair_dTCanair = (turbFluxCanair_dStateCanair - turbFluxCanair) / dx  ! derivative in net canopy air space fluxes w.r.t. canopy air temperature (W m-2 K-1)
  dTurbFluxCanopy_dTCanair = (turbFluxCanopy_dStateCanair - turbFluxCanopy) / dx  ! derivative in net canopy turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1) 
  dTurbFluxGround_dTCanair = (turbFluxGround_dStateCanair - turbFluxGround) / dx  ! derivative in net ground turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
  ! derivatives w.r.t. canopy temperature
  dTurbFluxCanair_dTCanopy = (turbFluxCanair_dStateCanopy - turbFluxCanair) / dx  ! derivative in net canopy air space fluxes w.r.t. canopy temperature (W m-2 K-1)
  dTurbFluxCanopy_dTCanopy = (turbFluxCanopy_dStateCanopy - turbFluxCanopy) / dx  ! derivative in net canopy turbulent fluxes w.r.t. canopy temperature (W m-2 K-1) 
  dTurbFluxGround_dTCanopy = (turbFluxGround_dStateCanopy - turbFluxGround) / dx  ! derivative in net ground turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
  ! derivatives w.r.t. ground temperature
  dTurbFluxCanair_dTGround = (turbFluxCanair_dStateGround - turbFluxCanair) / dx  ! derivative in net canopy air space fluxes w.r.t. ground temperature (W m-2 K-1)
  dTurbFluxCanopy_dTGround = (turbFluxCanopy_dStateGround - turbFluxCanopy) / dx  ! derivative in net canopy turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
  dTurbFluxGround_dTGround = (turbFluxGround_dStateGround - turbFluxGround) / dx  ! derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
 endif
 !if(heightCanopyBottom < scalarSnowDepth+z0Ground) pause 'bottom of the canopy is covered'

 ! test
 !print*, (ix_fDerivMeth == numerical)
 !print*, 'dTurbFluxCanair_dTCanair = ', dTurbFluxCanair_dTCanair
 !print*, 'dTurbFluxCanair_dTCanopy = ', dTurbFluxCanair_dTCanopy
 !print*, 'dTurbFluxCanair_dTGround = ', dTurbFluxCanair_dTGround
 !print*, 'dTurbFluxCanopy_dTCanair = ', dTurbFluxCanopy_dTCanair
 !print*, 'dTurbFluxCanopy_dTCanopy = ', dTurbFluxCanopy_dTCanopy
 !print*, 'dTurbFluxCanopy_dTGround = ', dTurbFluxCanopy_dTGround
 !print*, 'dTurbFluxGround_dTCanair = ', dTurbFluxGround_dTCanair
 !print*, 'dTurbFluxGround_dTCanopy = ', dTurbFluxGround_dTCanopy
 !print*, 'dTurbFluxGround_dTGround = ', dTurbFluxGround_dTGround
 !print*, '*****'
 !pause

 ! compute the heat advected with precipitation (W m-2)
 ! NOTE: fluxes are in kg m-2 s-1, so no need to use density of water/ice here
 scalarCanopyAdvectiveHeatFlux = -Cp_water*(scalarRainfall - scalarThroughfallRain)*(canopyTempTrial - scalarTwetbulb) + &
                                 (-Cp_ice)*(scalarSnowfall - scalarThroughfallSnow)*(canopyTempTrial - scalarTwetbulb)
 scalarGroundAdvectiveHeatFlux = -Cp_water*scalarThroughfallRain*(groundTempTrial - scalarTwetbulb)         + &
                                 (-Cp_ice)*scalarThroughfallSnow*(groundTempTrial - scalarTwetbulb)         !+ &
 !                                -Cp_water*scalarCanopyLiqDrainage  *(groundTempTrial - canopyTempTrial) + &
 !                                -Cp_ice  *scalarCanopySnowUnloading*(groundTempTrial - canopyTempTrial)
 !print*, 'scalarRainfall, scalarThroughfallRain, scalarSnowfall, scalarThroughfallSnow = ', scalarRainfall, scalarThroughfallRain, scalarSnowfall, scalarThroughfallSnow
 !print*, 'scalarCanopyAdvectiveHeatFlux, scalarGroundAdvectiveHeatFlux = ', scalarCanopyAdvectiveHeatFlux, scalarGroundAdvectiveHeatFlux

 ! compute the mass flux associated with transpiration and evaporation/sublimation (J m-2 s-1 --> kg m-2 s-1)
 ! NOTE: remove water from the snow on the ground in preference to removing water from the water in soil pores
 !print*, 'scalarLatHeatCanopyTrans = ', scalarLatHeatCanopyTrans
 !print*, 'scalarLatHeatGround      = ', scalarLatHeatGround
 ! (canopy transpiration/sublimation)
 if(scalarLatHeatSubVapCanopy > LH_vap+verySmall)then ! sublimation
  scalarCanopyEvaporation = 0._dp
  if(scalarLatHeatCanopyTrans > 0._dp)then ! flux directed towards the veg
   scalarCanopySublimation   = scalarCanopySublimation + scalarLatHeatCanopyTrans/LH_sub ! frost
   scalarCanopyTranspiration = 0._dp
  else
   scalarCanopySublimation   = scalarLatHeatCanopyEvap/LH_sub
   scalarCanopyTranspiration = scalarLatHeatCanopyTrans/LH_vap  ! transpiration is always vapor
  endif
 ! (canopy transpiration/evaporation)
 else                                                 ! evaporation
  scalarCanopyEvaporation = scalarLatHeatCanopyEvap/LH_vap
  scalarCanopySublimation = 0._dp
  if(scalarLatHeatCanopyTrans > 0._dp)then ! flux directed towards the veg
   scalarCanopyEvaporation   = scalarCanopyEvaporation + scalarLatHeatCanopyTrans/LH_vap
   scalarCanopyTranspiration = 0._dp
  else
   scalarCanopyTranspiration = scalarLatHeatCanopyTrans/LH_vap
  endif
 endif
 ! (ground evaporation/sublimation)
 if(scalarLatHeatSubVapGround > LH_vap+verySmall)then ! sublimation
  ! NOTE: this should only occur when we have formed snow layers, so check
  if(nSnow == 0)then; err=20; message=trim(message)//'only expect snow sublimation when we have formed some snow layers'; return; endif
  scalarGroundEvaporation = 0._dp  ! ground evaporation is zero once the snowpack has formed
  scalarSnowSublimation   = scalarLatHeatGround/LH_sub
 else
  ! NOTE: this should only occur when we have no snow layers, so check
  if(nSnow > 0)then; err=20; message=trim(message)//'only expect ground evaporation when there are no snow layers'; return; endif
  scalarGroundEvaporation = scalarLatHeatGround/LH_vap
  scalarSnowSublimation   = 0._dp  ! no sublimation from snow if no snow layers have formed
 endif
 !print*, 'scalarCanopySublimation, scalarLatHeatCanopyEvap = ', scalarCanopySublimation, scalarLatHeatCanopyEvap


 ! *******************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************
 ! ***** AND STITCH EVERYTHING TOGETHER  *****************************************************************************************************************************
 ! *******************************************************************************************************************************************************************
 ! *******************************************************************************************************************************************************************

 ! compute net fluxes at the canopy and ground surface
 canairNetFlux = turbFluxCanair
 canopyNetFlux = scalarCanopyAbsorbedSolar + scalarLWNetCanopy + turbFluxCanopy + scalarCanopyAdvectiveHeatFlux
 groundNetFlux = scalarGroundAbsorbedSolar + scalarLWNetGround + turbFluxGround + scalarGroundAdvectiveHeatFlux
 !print*, 'canairNetFlux, canopyNetFlux, scalarCanopyAbsorbedSolar,  scalarLWNetCanopy, turbFluxCanopy = ', canairNetFlux, canopyNetFlux, scalarCanopyAbsorbedSolar,  scalarLWNetCanopy, turbFluxCanopy
 !print*, 'groundNetFlux, scalarGroundAbsorbedSolar,  scalarLWNetGround, turbFluxGround = ', groundNetFlux, scalarGroundAbsorbedSolar,  scalarLWNetGround, turbFluxGround


 ! compute the derivatives
 dCanairNetFlux_dCanairTemp = dTurbFluxCanair_dTCanair
 dCanairNetFlux_dCanopyTemp = dTurbFluxCanair_dTCanopy
 dCanairNetFlux_dGroundTemp = dTurbFluxCanair_dTGround
 dCanopyNetFlux_dCanairTemp = dTurbFluxCanopy_dTCanair
 dCanopyNetFlux_dCanopyTemp = dLWNetCanopy_dTCanopy + dTurbFluxCanopy_dTCanopy - Cp_water*(scalarRainfall - scalarThroughfallRain) - Cp_ice*(scalarSnowfall - scalarThroughfallSnow)
 dCanopyNetFlux_dGroundTemp = dLWNetCanopy_dTGround + dTurbFluxCanopy_dTGround 
 dGroundNetFlux_dCanairTemp = dTurbFluxGround_dTCanair
 dGroundNetFlux_dCanopyTemp = dLWNetGround_dTCanopy + dTurbFluxGround_dTCanopy 
 dGroundNetFlux_dGroundTemp = dLWNetGround_dTGround + dTurbFluxGround_dTGround - Cp_water*scalarThroughfallRain - Cp_ice*scalarThroughfallSnow

 !print*, (ix_fDerivMeth == numerical)
 !print*, 'dCanopyNetFlux_dCanopyTemp = ', dCanopyNetFlux_dCanopyTemp
 !print*, 'dGroundNetFlux_dCanopyTemp = ', dGroundNetFlux_dCanopyTemp
 !print*, 'dCanopyNetFlux_dGroundTemp = ', dCanopyNetFlux_dGroundTemp
 !print*, 'dGroundNetFlux_dGroundTemp = ', dGroundNetFlux_dGroundTemp

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
 real(dp),parameter            :: fluxTolerance=1.e-10_dp  ! tolerance for energy closure (W m-2)
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
  ! NOTE: emc should be set to zero when not computing canopy fluxes
 
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

  !print*, 'LWRadCanopy = ', LWRadCanopy
  !print*, 'LWRadGround = ', LWRadGround

  !print*, 'LWNetCanopy = ', LWNetCanopy
  !print*, 'LWNetGround = ', LWNetGround
  !print*, 'LWNetUbound = ', LWNetUbound

  ! check the flux balance
  fluxBalance = LWNetUbound - (LWNetCanopy + LWNetGround)
  if(abs(fluxBalance) > fluxTolerance)then
   print*, 'fluxBalance = ', fluxBalance
   print*, 'emg, emc = ', emg, emc
   print*, 'TCan, TGnd = ', TCan, TGnd
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
   dLWNetGround_dTGround = -dLWRadGround_dTGround     ! derivative in net ground radiation w.r.t. ground temperature (W m-2 K-1)
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
                       iter,                          & ! intent(in): iteration index
                       computeVegFlux,                & ! intent(in): logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                       derivDesired,                  & ! intent(in): flag to indicate if analytical derivatives are desired
                       ixVegTraits,                   & ! intent(in): choice of parameterization for vegetation roughness length and displacement height
                       ixWindProfile,                 & ! intent(in): choice of canopy wind profile
                       ixStability,                   & ! intent(in): choice of stability function
                       ! input: above-canopy forcing data
                       mHeight,                       & ! intent(in): measurement height (m)
                       airtemp,                       & ! intent(in): air temperature at some height above the surface (K)
                       windspd,                       & ! intent(in): wind speed at some height above the surface (m s-1)
                       ! input: temperature (canopy, ground, canopy air space)
                       canairTemp,                    & ! intent(in): temperature of the canopy air space (K)
                       canopyTemp,                    & ! intent(in): canopy temperature (K)
                       groundTemp,                    & ! intent(in): ground temperature (K)
                       ! input: diagnostic variables
                       exposedVAI,                    & ! intent(in): exposed vegetation area index -- leaf plus stem (m2 m-2)
                       snowDepth,                     & ! intent(in): snow depth (m)
                       ! input: parameters
                       z0Ground,                      & ! intent(in): roughness length of the ground (below canopy or non-vegetated surface [snow]) (m)
                       z0CanopyParam,                 & ! intent(in): roughness length of the canopy (m)
                       zpdFraction,                   & ! intent(in): zero plane displacement / canopy height (-)
                       critRichNumber,                & ! intent(in): critical value for the bulk Richardson number where turbulence ceases (-)
                       Louis79_bparam,                & ! intent(in): parameter in Louis (1979) stability function
                       Louis79_cStar,                 & ! intent(in): parameter in Louis (1979) stability function
                       Mahrt87_eScale,                & ! intent(in): exponential scaling factor in the Mahrt (1987) stability function
                       windReductionParam,            & ! intent(in): canopy wind reduction parameter (-)                   
                       leafExchangeCoeff,             & ! intent(in): turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
                       leafDimension,                 & ! intent(in): characteristic leaf dimension (m)
                       heightCanopyTop,               & ! intent(in): height at the top of the vegetation canopy (m) 
                       heightCanopyBottom,            & ! intent(in): height at the bottom of the vegetation canopy (m) 
                       ! input: stability correction from the last iteration
                       canopyStabilityCorrection_old, & ! intent(in): stability correction for the canopy (-)
                       groundStabilityCorrection_old, & ! intent(in): stability correction for the ground surface (-)
                       ! output: stability corrections
                       RiBulkCanopy,                  & ! intent(out): bulk Richardson number for the canopy (-)
                       RiBulkGround,                  & ! intent(out): bulk Richardson number for the ground surface (-)
                       canopyStabilityCorrection,     & ! intent(out): stability correction for the canopy (-)
                       groundStabilityCorrection,     & ! intent(out): stability correction for the ground surface (-)
                       ! output: scalar resistances
                       z0Canopy,                      & ! intent(out): roughness length of the canopy (m)
                       windReductionFactor,           & ! intent(out): canopy wind reduction factor (-)
                       zeroPlaneDisplacement,         & ! intent(out): zero plane displacement (m) 
                       eddyDiffusCanopyTop,           & ! intent(out): eddy diffusivity for heat at the top of the canopy (m2 s-1)
                       frictionVelocity,              & ! intent(out): friction velocity (m s-1)
                       windspdCanopyTop,              & ! intent(out): windspeed at the top of the canopy (m s-1)
                       windspdCanopyBottom,           & ! intent(out): windspeed at the height of the bottom of the canopy (m s-1)
                       leafResistance,                & ! intent(out): mean leaf boundary layer resistance per unit leaf area (s m-1)
                       groundResistance,              & ! intent(out): below canopy aerodynamic resistance (s m-1) 
                       canopyResistance,              & ! intent(out): above canopy aerodynamic resistance (s m-1)
                       ! output: derivatives in scalar resistances
                       dGroundResistance_dTGround,    & ! intent(out): derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
                       dGroundResistance_dTCanopy,    & ! intent(out): derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
                       dGroundResistance_dTCanair,    & ! intent(out): derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
                       dCanopyResistance_dTCanopy,    & ! intent(out): derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
                       dCanopyResistance_dTCanair,    & ! intent(out): derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)
                       ! output: error control
                       err,message                    ) ! intent(out): error control
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! compute aerodynamic resistances
 ! Refs: Choudhury and Monteith (4-layer model for heat budget of homogenous surfaces; QJRMS, 1988)
 !       Niu and Yang (Canopy effects on snow processes; JGR, 2004)
 !       Mahat et al. (Below-canopy turbulence in a snowmelt model, WRR, 2012)
 implicit none
 ! input: model control
 integer(i4b),intent(in)       :: iter                          ! iteration index
 logical(lgt),intent(in)       :: computeVegFlux                ! logical flag to compute vegetation fluxes (.false. if veg buried by snow)
 logical(lgt),intent(in)       :: derivDesired                  ! logical flag to indicate if analytical derivatives are desired
 integer(i4b),intent(in)       :: ixVegTraits                   ! choice of parameterization for vegetation roughness length and displacement height
 integer(i4b),intent(in)       :: ixWindProfile                 ! choice of canopy wind profile
 integer(i4b),intent(in)       :: ixStability                   ! choice of stability function
 ! input: above-canopy forcing data
 real(dp),intent(in)           :: mHeight                       ! measurement height (m)
 real(dp),intent(in)           :: airtemp                       ! air temperature at some height above the surface (K)
 real(dp),intent(in)           :: windspd                       ! wind speed at some height above the surface (m s-1)
 ! input: temperature (canopy, ground, canopy air space)
 real(dp),intent(in)           :: canairTemp                    ! temperature of the canopy air space (K)
 real(dp),intent(in)           :: canopyTemp                    ! canopy temperature (K)
 real(dp),intent(in)           :: groundTemp                    ! ground temperature (K)
 ! input: diagnostic variables
 real(dp),intent(in)           :: exposedVAI                    ! exposed vegetation area index -- leaf plus stem (m2 m-2)
 real(dp),intent(in)           :: snowDepth                     ! snow depth (m)
 ! input: parameters
 real(dp),intent(in)           :: z0Ground                      ! roughness length of the ground (below canopy or non-vegetated surface [snow]) (m)
 real(dp),intent(in)           :: z0CanopyParam                 ! roughness length of the canopy (m)
 real(dp),intent(in)           :: zpdFraction                   ! zero plane displacement / canopy height (-)
 real(dp),intent(in)           :: critRichNumber                ! critical value for the bulk Richardson number where turbulence ceases (-)
 real(dp),intent(in)           :: Louis79_bparam                ! parameter in Louis (1979) stability function
 real(dp),intent(in)           :: Louis79_cStar                 ! parameter in Louis (1979) stability function
 real(dp),intent(in)           :: Mahrt87_eScale                ! exponential scaling factor in the Mahrt (1987) stability function
 real(dp),intent(in)           :: windReductionParam            ! canopy wind reduction parameter (-)                   
 real(dp),intent(in)           :: leafExchangeCoeff             ! turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
 real(dp),intent(in)           :: leafDimension                 ! characteristic leaf dimension (m)
 real(dp),intent(in)           :: heightCanopyTop               ! height at the top of the vegetation canopy (m) 
 real(dp),intent(in)           :: heightCanopyBottom            ! height at the bottom of the vegetation canopy (m) 
 ! input: stability correction from the last iteration
 real(dp),intent(in)           :: canopyStabilityCorrection_old ! stability correction for the canopy (-)
 real(dp),intent(in)           :: groundStabilityCorrection_old ! stability correction for the ground surface (-)
 ! output: stability corrections
 real(dp),intent(out)          :: RiBulkCanopy                  ! bulk Richardson number for the canopy (-)
 real(dp),intent(out)          :: RiBulkGround                  ! bulk Richardson number for the ground surface (-)
 real(dp),intent(out)          :: canopyStabilityCorrection     ! stability correction for the canopy (-)
 real(dp),intent(out)          :: groundStabilityCorrection     ! stability correction for the ground surface (-)
 ! output: scalar resistances
 real(dp),intent(out)          :: z0Canopy                      ! roughness length of the vegetation canopy (m)
 real(dp),intent(out)          :: windReductionFactor           ! canopy wind reduction factor (-)
 real(dp),intent(out)          :: zeroPlaneDisplacement         ! zero plane displacement (m) 
 real(dp),intent(out)          :: eddyDiffusCanopyTop           ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
 real(dp),intent(out)          :: frictionVelocity              ! friction velocity (m s-1)
 real(dp),intent(out)          :: windspdCanopyTop              ! windspeed at the top of the canopy (m s-1)
 real(dp),intent(out)          :: windspdCanopyBottom           ! windspeed at the height of the bottom of the canopy (m s-1)
 real(dp),intent(out)          :: leafResistance                ! mean leaf boundary layer resistance per unit leaf area (s m-1)
 real(dp),intent(out)          :: groundResistance              ! below canopy aerodynamic resistance (s m-1) 
 real(dp),intent(out)          :: canopyResistance              ! above canopy aerodynamic resistance (s m-1)
 ! output: derivatives in scalar resistances
 real(dp),intent(out)          :: dGroundResistance_dTGround    ! derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
 real(dp),intent(out)          :: dGroundResistance_dTCanopy    ! derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
 real(dp),intent(out)          :: dGroundResistance_dTCanair    ! derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
 real(dp),intent(out)          :: dCanopyResistance_dTCanopy    ! derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
 real(dp),intent(out)          :: dCanopyResistance_dTCanair    ! derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                           ! error code
 character(*),intent(out)      :: message                       ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! local variables: general
 character(LEN=256)            :: cmessage                      ! error message of downwind routine
 ! local variables: vegetation roughness and dispalcement height
 real(dp),parameter            :: oneThird=1._dp/3._dp          ! 1/3
 real(dp),parameter            :: twoThirds=2._dp/3._dp         ! 2/3
 real(dp),parameter            :: C_r = 0.3                     ! roughness element drag coefficient (-) from Raupach (BLM, 1994) 
 real(dp),parameter            :: C_s = 0.003_dp                ! substrate surface drag coefficient (-) from Raupach (BLM, 1994)
 real(dp),parameter            :: approxDragCoef_max = 0.3_dp   ! maximum value of the approximate drag coefficient (-) from Raupach (BLM, 1994)
 real(dp),parameter            :: psi_h = 0.193_dp              ! roughness sub-layer influence function (-) from Raupach (BLM, 1994)
 real(dp),parameter            :: c_d1 = 7.5_dp                 ! scaling parameter used to define displacement height (-) from Raupach (BLM, 1994)
 real(dp),parameter            :: cd_CM = 0.2_dp                ! mean drag coefficient for individual leaves (-) from Choudhury and Monteith (QJRMS, 1988)
 real(dp)                      :: funcLAI                       ! temporary variable to calculate zero plane displacement for the canopy
 real(dp)                      :: fracCanopyHeight              ! zero plane displacement expressed as a fraction of canopy height
 real(dp)                      :: approxDragCoef                ! approximate drag coefficient used in the computation of canopy roughness length (-)
 ! local variables: resistance
 real(dp)                      :: canopyExNeut                  ! surface-atmosphere exchange coefficient under neutral conditions (-)
 real(dp)                      :: groundExNeut                  ! surface-atmosphere exchange coefficient under neutral conditions (-)
 real(dp)                      :: sfc2AtmExchangeCoeff_canopy   ! surface-atmosphere exchange coefficient after stability corrections (-)
 real(dp)                      :: groundResistanceNeutral       ! ground resistance under neutral conditions (s m-1)
 real(dp)                      :: windConvFactorTop             ! factor to convert friction velocity to wind speed at top of canopy (-)
 real(dp)                      :: windConvFactorBottom          ! factor to convert wind speed at top of canopy to wind speed at bottom of canopy (-)
 real(dp)                      :: referenceHeight               ! reference height used to compute above-ground windspeed (m)
 real(dp)                      :: heightAboveGround             ! height above the snow surface (m)
 ! local variables: derivatives
 real(dp)                      :: dFV_dT                        ! derivative in friction velocity w.r.t. canopy air temperature
 real(dp)                      :: dED_dT                        ! derivative in eddy diffusivity at the top of the canopy w.r.t. canopy air temperature
 real(dp)                      :: dGR_dT                        ! derivative in neutral ground resistance w.r.t. canopy air temperature
 real(dp)                      :: tmp1,tmp2                     ! temporary variables used in calculation of ground resistance
 real(dp)                      :: dLeafResistance_dStateCanopy         ! leaf resistance after perturbing canopy temperature (s m-1)
 real(dp)                      :: dCanopyResistance_dStateCanopy       ! canopy resistance after perturbing canopy temperature (s m-1)   
 real(dp)                      :: dGroundResistance_dStateCanopy       ! ground resistance after perturbing canopy temperature (s m-1)
 real(dp)                      :: dGroundResistance_dStateGround       ! ground resistance after perturbing ground temperature (s m-1)
 real(dp)                      :: dCanopyStabilityCorrection_dRich     ! derivative in stability correction w.r.t. Richardson number for the canopy (-)
 real(dp)                      :: dGroundStabilityCorrection_dRich     ! derivative in stability correction w.r.t. Richardson number for the ground surface (-)
 real(dp)                      :: dCanopyStabilityCorrection_dAirTemp  ! (not used) derivative in stability correction w.r.t. air temperature (K-1) 
 real(dp)                      :: dGroundStabilityCorrection_dAirTemp  ! (not used) derivative in stability correction w.r.t. air temperature (K-1) 
 real(dp)                      :: dCanopyStabilityCorrection_dCasTemp  ! derivative in canopy stability correction w.r.t. canopy air space temperature (K-1)
 real(dp)                      :: dGroundStabilityCorrection_dCasTemp  ! derivative in ground stability correction w.r.t. canopy air space temperature (K-1)
 real(dp)                      :: dGroundStabilityCorrection_dSfcTemp  ! derivative in ground stability correction w.r.t. surface temperature (K-1)
 real(dp)                      :: singleLeafConductance         ! leaf boundary layer conductance (m s-1) 
 real(dp)                      :: canopyLeafConductance         ! leaf boundary layer conductance -- scaled up to the canopy (m s-1)
 real(dp)                      :: leaf2CanopyScaleFactor        ! factor to scale from the leaf to the canopy [m s-(1/2)]
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='aeroResist/'

 ! check that measurement height is above the top of the canopy
 if(mHeight < heightCanopyTop)then
  err=20; message=trim(message)//'measurement height is below the top of the canopy'; return
 endif

 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! * compute vegetation poperties (could be done at the same time as phenology.. does not have to be in the flux routine!)

 if(computeVegFlux) then ! (if vegetation is exposed)

  ! ***** identify zero plane displacement, roughness length, and surface temperature for the canopy (m)
  select case(ixVegTraits)

   ! Raupach (BLM 1994) "Simplified expressions..."
   case(Raupach_BLM1994)
    ! (compute zero-plane displacement)
    funcLAI          = sqrt(c_d1*exposedVAI)
    fracCanopyHeight = -(1._dp - exp(-funcLAI))/funcLAI + 1._dp
    zeroPlaneDisplacement = fracCanopyHeight*heightCanopyTop
    ! (coupute roughness length of the veg canopy)
    approxDragCoef   = min( sqrt(C_s + C_r*exposedVAI/2._dp), approxDragCoef_max)
    z0Canopy         = (1._dp - fracCanopyHeight) * exp(-vkc*approxDragCoef - psi_h) * heightCanopyTop

   ! Choudhury and Monteith (QJRMS 1998) "A four layer model for the heat budget..."
   case(CM_QJRMS1998)
    funcLAI =  cd_CM*exposedVAI
    zeroPlaneDisplacement = 1.1_dp*heightCanopyTop*log(1._dp + funcLAI**0.25_dp)
    if(funcLAI < 0.2_dp)then
     z0Canopy = z0Ground + 0.3_dp*heightCanopyTop*funcLAI**0.5_dp
    else
     z0Canopy = 0.3_dp*heightCanopyTop*(1._dp - zeroPlaneDisplacement/heightCanopyTop)
    endif

   ! constant parameters dependent on the vegetation type
   case(vegTypeTable)
    zeroPlaneDisplacement = zpdFraction*heightCanopyTop  ! zero-plane displacement (m)
    z0Canopy = z0CanopyParam                             ! roughness length of the veg canopy (m)

   ! check
   case default
    err=10; message=trim(message)//"unknown parameterization for vegetation roughness length and displacement height"; return

  end select  ! vegetation traits (z0, zpd)

  ! correct for snow depth
  if(zeroPlaneDisplacement < snowDepth) zeroPlaneDisplacement = snowDepth

  ! check that everything is consistent
  if(mHeight < zeroPlaneDisplacement)then; err=20; message=trim(message)//'measurement height is below the displacement height'; return; endif
  if(mHeight < z0Canopy)then; err=20; message=trim(message)//'measurement height is below the roughness length'; return; endif

 endif  ! if there is a canopy

 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! * compute resistance for the case where the canopy is exposed

 ! check if vegetation is exposed
 if(computeVegFlux) then 

  ! compute the stability correction for resistance from canopy air space to air above the canopy (-)
  call aStability(&
                  ! input
                  derivDesired,                                     & ! input: logical flag to compute analytical derivatives
                  ixStability,                                      & ! input: choice of stability function
                  ! input: forcing data, diagnostic and state variables
                  mHeight,                                          & ! input: measurement height (m)
                  airTemp,                                          & ! input: air temperature above the canopy (K)
                  canairTemp,                                       & ! input: temperature of the canopy air space (K)
                  windspd,                                          & ! input: wind speed above the canopy (m s-1)
                  ! input: stability parameters
                  critRichNumber,                                   & ! input: critical value for the bulk Richardson number where turbulence ceases (-)
                  Louis79_bparam,                                   & ! input: parameter in Louis (1979) stability function
                  Louis79_cStar,                                    & ! input: parameter in Louis (1979) stability function
                  Mahrt87_eScale,                                   & ! input: exponential scaling factor in the Mahrt (1987) stability function
                  ! output
                  RiBulkCanopy,                                     & ! output: bulk Richardson number (-)
                  canopyStabilityCorrection,                        & ! output: stability correction for turbulent heat fluxes (-)
                  dCanopyStabilityCorrection_dRich,                 & ! output: derivative in stability correction w.r.t. Richardson number for the canopy (-)
                  dCanopyStabilityCorrection_dAirTemp,              & ! output: (not used) derivative in stability correction w.r.t. air temperature (K-1)
                  dCanopyStabilityCorrection_dCasTemp,              & ! output: derivative in stability correction w.r.t. canopy air space temperature (K-1)
                  err, cmessage                                     ) ! output: error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! compute turbulent exchange coefficient (-)
  canopyExNeut = (vkc**2._dp) / ( log((mHeight - zeroPlaneDisplacement)/z0Canopy))**2._dp     ! coefficient under conditions of neutral stability
  sfc2AtmExchangeCoeff_canopy = canopyExNeut*canopyStabilityCorrection                        ! after stability corrections

  ! compute the friction velocity (m s-1)
  frictionVelocity = windspd * sqrt(sfc2AtmExchangeCoeff_canopy)

  ! compute the above-canopy resistance (s m-1)
  canopyResistance = 1._dp/(sfc2AtmExchangeCoeff_canopy*windspd)
  if(canopyResistance < 0._dp)then; err=20; message=trim(message)//'canopy resistance < 0'; return; endif
  !write(*,'(a,10(f20.10,1x))') 'in aeroResist: windspd, canairTemp, canopyExNeut, canopyStabilityCorrection, canopyResistance = ', &
  !                                             windspd, canairTemp, canopyExNeut, canopyStabilityCorrection, canopyResistance

  ! compute windspeed at the top of the canopy (m s-1)
  ! NOTE: stability corrections cancel out
  windConvFactorTop = log((heightCanopyTop - zeroPlaneDisplacement)/z0Canopy) / log((mHeight - zeroPlaneDisplacement)/z0Canopy)
  windspdCanopyTop  = windspd*windConvFactorTop
  !windConvFactorTop = log((heightCanopyTop - zeroPlaneDisplacement)/z0Canopy)/(sqrt(canopyStabilityCorrection)*vkc)
  !windspdCanopyTop  = frictionVelocity*windConvFactorTop

  ! compute the windspeed reduction
  ! Refs: Norman et al. (Ag. Forest Met., 1995) -- citing Goudriaan (1977 manuscript "crop micrometeorology: a simulation study", Wageningen).
  windReductionFactor = windReductionParam * exposedVAI**twoThirds * (heightCanopyTop - heightCanopyBottom)**oneThird / leafDimension**oneThird
  !windReductionFactor = 3._dp

  ! compute windspeed at the bottom of the canopy (m s-1)
  !referenceHeight      = max(heightCanopyBottom, min(0.5_dp, heightCanopyTop))
  referenceHeight      = max(heightCanopyBottom, snowDepth+z0Ground)
  windConvFactorBottom = exp(-windReductionFactor*(1._dp - referenceHeight/heightCanopyTop))
  windspdCanopyBottom  = windspdCanopyTop*windConvFactorBottom
  if(referenceHeight > z0Canopy+zeroPlaneDisplacement)then; err=20; message=trim(message)//'reference height > z0Canopy+zeroPlaneDisplacement'; return; endif

  ! compute the leaf boundary layer resistance (s m-1)
  singleLeafConductance  = leafExchangeCoeff*sqrt(windspdCanopyTop/leafDimension)
  leaf2CanopyScaleFactor = (2._dp/windReductionFactor) * (1._dp - exp(-windReductionFactor/2._dp)) ! factor to scale from the leaf to the canopy
  canopyLeafConductance  = singleLeafConductance*leaf2CanopyScaleFactor
  leafResistance  = 1._dp/(canopyLeafConductance)
  if(leafResistance < 0._dp)then; err=20; message=trim(message)//'leaf resistance < 0'; return; endif

  ! compute eddy diffusivity for heat at the top of the canopy (m2 s-1)
  !   Note: use of friction velocity here includes stability adjustments
  eddyDiffusCanopyTop = max(vkc*FrictionVelocity*(heightCanopyTop - zeroPlaneDisplacement), mpe)  ! (avoid divide by zero)

  ! compute the resistance between the surface and canopy air UNDER NEUTRAL CONDITIONS (s m-1)
  select case(ixWindProfile)
   ! case 1: assume exponential profile extends from the surface roughness length to the displacement height plus vegetation roughness
   case(exponential)
    tmp1 = exp(-windReductionFactor* (snowDepth+z0Ground)/heightCanopyTop)
    tmp2 = exp(-windReductionFactor*(z0Canopy+zeroPlaneDisplacement)/heightCanopyTop)
    groundResistanceNeutral = ( heightCanopyTop*exp(windReductionFactor) / (windReductionFactor*eddyDiffusCanopyTop) ) * (tmp1 - tmp2)   ! s m-1
   ! case 2: logarithmic profile below the canopy
   case(logBelowCanopy)
    tmp1 = exp(-windReductionFactor* referenceHeight/heightCanopyTop)
    tmp2 = exp(-windReductionFactor*(z0Canopy+zeroPlaneDisplacement)/heightCanopyTop)
    if(referenceHeight > heightCanopyBottom)then  ! snow is above the bottom of the canopy -- just use the exponential profile
     groundResistanceNeutral = ( heightCanopyTop*exp(windReductionFactor) / (windReductionFactor*eddyDiffusCanopyTop) ) * (tmp1 - tmp2)   ! s m-1
    else  ! snow is below the bottom of the canopy
     groundResistanceNeutral = ( heightCanopyTop*exp(windReductionFactor) / (windReductionFactor*eddyDiffusCanopyTop) ) * (tmp1 - tmp2) & ! s m-1
                                  + (1._dp/(max(0.1_dp,windspdCanopyBottom)*vkc**2._dp))*(log((referenceHeight - snowDepth)/z0Ground))**2._dp
    endif
   ! check that we identified the option   
   case default
    err=20; message=trim(message)//'cannot identify option for canopy wind profile'; return
   end select

  ! compute the stability correction for resistance from the ground to the canopy air space (-)
  call aStability(&
                  ! input
                  derivDesired,                                     & ! input: logical flag to compute analytical derivatives
                  ixStability,                                      & ! input: choice of stability function
                  ! input: forcing data, diagnostic and state variables
                  referenceHeight,                                  & ! input: reference height of wind within the canopy (m)
                  canairTemp,                                       & ! input: temperature of the canopy air space (K)
                  groundTemp,                                       & ! input: temperature of the ground surface (K)
                  max(0.1_dp,windspdCanopyBottom),                  & ! input: wind speed at the reference height (m s-1)
                  ! input: stability parameters
                  critRichNumber,                                   & ! input: critical value for the bulk Richardson number where turbulence ceases (-)
                  Louis79_bparam,                                   & ! input: parameter in Louis (1979) stability function
                  Louis79_cStar,                                    & ! input: parameter in Louis (1979) stability function
                  Mahrt87_eScale,                                   & ! input: exponential scaling factor in the Mahrt (1987) stability function
                  ! output
                  RiBulkGround,                                     & ! output: bulk Richardson number (-)
                  groundStabilityCorrection,                        & ! output: stability correction for turbulent heat fluxes (-)
                  dGroundStabilityCorrection_dRich,                 & ! output: derivative in stability correction w.r.t. Richardson number for the canopy (-)
                  dGroundStabilityCorrection_dCasTemp,              & ! output: derivative in stability correction w.r.t. canopy air space temperature (K-1)
                  dGroundStabilityCorrection_dSfcTemp,              & ! output: derivative in stability correction w.r.t. surface temperature (K-1)
                  err, cmessage                                     ) ! output: error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! compute the ground resistance
  groundResistance = groundResistanceNeutral / groundStabilityCorrection
  if(groundResistance < 0._dp)then; err=20; message=trim(message)//'ground resistance < 0 [vegetation is present]'; return; endif

 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! * compute resistance for the case without a canopy (bare ground, or canopy completely buried with snow)
 else

  ! no canopy, so set huge resistances (not used)
  canopyResistance = 1.e12_dp   ! not used: huge resistance, so conductance is essentially zero
  leafResistance   = 1.e12_dp   ! not used: huge resistance, so conductance is essentially zero
 
  ! check that measurement height above the ground surface is above the roughness length
  if(mHeight < snowDepth+z0Ground)then; err=20; message=trim(message)//'measurement height < snow depth + roughness length'; return; endif
 
  ! compute the resistance between the surface and canopy air UNDER NEUTRAL CONDITIONS (s m-1)
  groundExNeut = (vkc**2._dp) / ( log((mHeight - snowDepth)/z0Ground)**2._dp) ! turbulent transfer coefficient under conditions of neutral stability (-)
  groundResistanceNeutral = 1._dp / (groundExNeut*windspd)
 
  ! define height above the snow surface
  heightAboveGround  = mHeight - snowDepth
 
  ! check that measurement height above the ground surface is above the roughness length
  if(heightAboveGround < z0Ground)then
   print*, 'z0Ground = ', z0Ground
   print*, 'mHeight  = ', mHeight
   print*, 'snowDepth = ', snowDepth
   print*, 'heightAboveGround = ', heightAboveGround
   message=trim(message)//'height above ground < roughness length [likely due to snow accumulation]'
   err=20; return
  endif

  ! compute ground stability correction
  call aStability(&
                   ! input
                  derivDesired,                                     & ! input: logical flag to compute analytical derivatives
                  ixStability,                                      & ! input: choice of stability function
                  ! input: forcing data, diagnostic and state variables
                  heightAboveGround,                                & ! input: measurement height above the ground surface (m)
                  airtemp,                                          & ! input: temperature above the ground surface (K)
                  groundTemp,                                       & ! input: trial value of surface temperature -- "surface" is either canopy or ground (K)
                  windspd,                                          & ! input: wind speed above the ground surface (m s-1)
                  ! input: stability parameters
                  critRichNumber,                                   & ! input: critical value for the bulk Richardson number where turbulence ceases (-)
                  Louis79_bparam,                                   & ! input: parameter in Louis (1979) stability function
                  Louis79_cStar,                                    & ! input: parameter in Louis (1979) stability function
                  Mahrt87_eScale,                                   & ! input: exponential scaling factor in the Mahrt (1987) stability function
                  ! output
                  RiBulkGround,                                     & ! output: bulk Richardson number (-)
                  groundStabilityCorrection,                        & ! output: stability correction for turbulent heat fluxes (-)
                  dGroundStabilityCorrection_dRich,                 & ! output: derivative in stability correction w.r.t. Richardson number for the ground surface (-)
                  dGroundStabilityCorrection_dAirTemp,              & ! output: (not used) derivative in stability correction w.r.t. air temperature (K-1)
                  dGroundStabilityCorrection_dSfcTemp,              & ! output: derivative in stability correction w.r.t. surface temperature (K-1)
                  err, cmessage                                     ) ! output: error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
 
  ! compute the ground resistance (after stability corrections)
  groundResistance = groundResistanceNeutral/groundStabilityCorrection
  if(groundResistance < 0._dp)then; err=20; message=trim(message)//'ground resistance < 0 [no vegetation]'; return; endif

  ! set all canopy variables to missing (no canopy!)
  z0Canopy                   = missingValue   ! roughness length of the vegetation canopy (m)
  RiBulkCanopy               = missingValue   ! bulk Richardson number for the canopy (-)
  windReductionFactor        = missingValue   ! canopy wind reduction factor (-)
  zeroPlaneDisplacement      = missingValue   ! zero plane displacement (m) 
  canopyStabilityCorrection  = missingValue   ! stability correction for the canopy (-)
  eddyDiffusCanopyTop        = missingValue   ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
  frictionVelocity           = missingValue   ! friction velocity (m s-1)
  windspdCanopyTop           = missingValue   ! windspeed at the top of the canopy (m s-1)
  windspdCanopyBottom        = missingValue   ! windspeed at the height of the bottom of the canopy (m s-1)

 endif  ! (if no canopy)

 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! * compute derivatives
 if(derivDesired)then  ! if analytical derivatives are desired

  ! derivatives for the vegetation canopy
  if(computeVegFlux) then ! (if vegetation is exposed) 
  
   ! ***** compute derivatives w.r.t. canopy temperature
   ! NOTE: derivatives are zero because using canopy air space temperature
   dCanopyResistance_dTCanopy = 0._dp ! derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
   dGroundResistance_dTCanopy = 0._dp ! derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)

   ! ***** compute derivatives w.r.t. ground temperature (s m-1 K-1)
   dGroundResistance_dTGround = -(groundResistanceNeutral*dGroundStabilityCorrection_dSfcTemp)/(groundStabilityCorrection**2._dp)

   ! ***** compute derivatives w.r.t. temperature of the canopy air space (s m-1 K-1)
   ! derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)
   dCanopyResistance_dTCanair = -dCanopyStabilityCorrection_dCasTemp/(windspd*canopyExNeut*canopyStabilityCorrection**2._dp)
   ! derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
   ! (compute derivative in NEUTRAL ground resistance w.r.t. canopy air temperature (s m-1 K-1))
   dFV_dT = windspd*canopyExNeut*dCanopyStabilityCorrection_dCasTemp/(sqrt(sfc2AtmExchangeCoeff_canopy)*2._dp)                ! d(frictionVelocity)/d(canopy air temperature)
   dED_dT = dFV_dT*vkc*(heightCanopyTop - zeroPlaneDisplacement)                                                             ! d(eddyDiffusCanopyTop)d(canopy air temperature)
   dGR_dT = -dED_dT*(tmp1 - tmp2)*heightCanopyTop*exp(windReductionFactor) / (windReductionFactor*eddyDiffusCanopyTop**2._dp) ! d(groundResistanceNeutral)/d(canopy air temperature)
   ! (stitch everything together -- product rule)
   dGroundResistance_dTCanair = dGR_dT/groundStabilityCorrection - groundResistanceNeutral*dGroundStabilityCorrection_dCasTemp/(groundStabilityCorrection**2._dp)
 
  ! ***** compute resistances for non-vegetated surfaces (e.g., snow)
  else
 
   ! set canopy derivatives to zero (non-vegetated, remember)
   dCanopyResistance_dTCanopy = 0._dp
   dGroundResistance_dTCanopy = 0._dp
 
   ! compute derivatives for ground resistance
   dGroundResistance_dTGround = -dGroundStabilityCorrection_dSfcTemp/(windspd*groundExNeut*groundStabilityCorrection**2._dp)
 
  endif  ! (switch between vegetated and non-vegetated surfaces)

 ! * analytical derivatives not desired
 else
  dGroundResistance_dTGround = missingValue
  dGroundResistance_dTCanopy = missingValue
  dCanopyResistance_dTCanopy = missingValue
 endif

 ! test
 !print*, 'dGroundResistance_dTGround = ', dGroundResistance_dTGround
 !print*, 'dGroundResistance_dTCanopy = ', dGroundResistance_dTCanopy
 !print*, 'dCanopyResistance_dTCanopy = ', dCanopyResistance_dTCanopy
 !pause 'in aeroResist'

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
 USE mDecisions_module, only: NoahType,CLM_Type,SiB_Type  ! options for the choice of function for the soil moisture control on stomatal resistance
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
   case(CLM_Type)  ! thresholded linear function of matric head
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
  mLayerTranspireLimitFac(iLayer) = min(1._dp, max(tiny(gx),gx) )
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
                       computeVegFlux,                & ! intent(in): logical flag to compute vegetation fluxes (.false. if veg buried by snow)
                       ixDerivMethod,                 & ! intent(in): choice of method used to compute derivative (analytical or numerical)
                       ! input: above-canopy forcing data
                       airtemp,                       & ! intent(in): air temperature at some height above the surface (K)
                       airpres,                       & ! intent(in): air pressure of the air above the vegetation canopy (Pa)
                       VPair,                         & ! intent(in): vapor pressure of the air above the vegetation canopy (Pa)
                       ! input: latent heat of sublimation/vaporization
                       latHeatSubVapCanopy,           & ! intent(in): latent heat of sublimation/vaporization for the vegetation canopy (J kg-1)
                       latHeatSubVapGround,           & ! intent(in): latent heat of sublimation/vaporization for the ground surface (J kg-1)
                       ! input: canoopy liquid and canopy ice (used as a solution constraint)
                       canopyLiquid,                  & ! intent(in): mass of liquid water on the vegetation canopy (kg m-2)
                       canopyIce,                     & ! intent(in): mass of ice water on the vegetation canopy (kg m-2)
                       ! input: canopy and ground temperature
                       canairTemp,                    & ! intent(in): temperature of the canopy air space (K)
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
                       dGroundResistance_dTGround,    & ! intent(in): derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
                       dGroundResistance_dTCanopy,    & ! intent(in): derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
                       dGroundResistance_dTCanair,    & ! intent(in): derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
                       dCanopyResistance_dTCanopy,    & ! intent(in): derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
                       dCanopyResistance_dTCanair,    & ! intent(in): derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)
                       ! output: conductances (used to check derivative calculations)
                       leafConductance,               & ! intent(out): leaf conductance (m s-1)
                       canopyConductance,             & ! intent(out): canopy conductance (m s-1)
                       groundConductanceSH,           & ! intent(out): ground conductance for sensible heat (m s-1)
                       groundConductanceLH,           & ! intent(out): ground conductance for latent heat -- includes soil resistance (m s-1)
                       evapConductance,               & ! intent(out): conductance for evaporation (m s-1)
                       transConductance,              & ! intent(out): conductance for transpiration (m s-1)
                       totalConductanceSH,            & ! intent(out): total conductance for sensible heat (m s-1)
                       totalConductanceLH,            & ! intent(out): total conductance for latent heat (m s-1)
                       ! output: canopy air space variables
                       VP_CanopyAir,                  & ! intent(out): vapor pressure of the canopy air space (Pa)
                       ! output: fluxes from the vegetation canopy
                       senHeatCanopy,                 & ! intent(out): sensible heat flux from the canopy to the canopy air space (W m-2)
                       latHeatCanopyEvap,             & ! intent(out): latent heat flux associated with evaporation from the canopy to the canopy air space (W m-2)
                       latHeatCanopyTrans,            & ! intent(out): latent heat flux associated with transpiration from the canopy to the canopy air space (W m-2)
                       ! output: fluxes from non-vegetated surfaces (ground surface below vegetation, bare ground, or snow covered vegetation)
                       senHeatGround,                 & ! intent(out): sensible heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
                       latHeatGround,                 & ! intent(out): latent heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
                       ! output: total heat fluxes to the atmosphere
                       senHeatTotal,                  & ! intent(out): total sensible heat flux to the atmosphere (W m-2)
                       latHeatTotal,                  & ! intent(out): total latent heat flux to the atmosphere (W m-2)
                       ! output: net fluxes
                       turbFluxCanair,                & ! intent(out): net turbulent heat fluxes at the canopy air space (W m-2)
                       turbFluxCanopy,                & ! intent(out): net turbulent heat fluxes at the canopy (W m-2)
                       turbFluxGround,                & ! intent(out): net turbulent heat fluxes at the ground surface (W m-2)
                       ! output: flux derivatives
                       dTurbFluxCanair_dTCanair,      & ! intent(out): derivative in net canopy air space fluxes w.r.t. canopy air temperature (W m-2 K-1)
                       dTurbFluxCanair_dTCanopy,      & ! intent(out): derivative in net canopy air space fluxes w.r.t. canopy temperature (W m-2 K-1)
                       dTurbFluxCanair_dTGround,      & ! intent(out): derivative in net canopy air space fluxes w.r.t. ground temperature (W m-2 K-1)
                       dTurbFluxCanopy_dTCanair,      & ! intent(out): derivative in net canopy turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
                       dTurbFluxCanopy_dTCanopy,      & ! intent(out): derivative in net canopy turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
                       dTurbFluxCanopy_dTGround,      & ! intent(out): derivative in net canopy turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
                       dTurbFluxGround_dTCanair,      & ! intent(out): derivative in net ground turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
                       dTurbFluxGround_dTCanopy,      & ! intent(out): derivative in net ground turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
                       dTurbFluxGround_dTGround,      & ! intent(out): derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
                       ! output: error control
                       err,message                    ) ! intent(out): error control
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 implicit none
 ! input: model control
 real(dp),intent(in)           :: dt                    ! model time step (seconds)
 logical(lgt),intent(in)       :: computeVegFlux        ! logical flag to compute vegetation fluxes (.false. if veg buried by snow)
 integer(i4b),intent(in)       :: ixDerivMethod         ! choice of method used to compute derivative (analytical or numerical)
 ! input: above-canopy forcing data
 real(dp),intent(in)           :: airtemp               ! air temperature at some height above the surface (K)
 real(dp),intent(in)           :: airpres               ! air pressure of the air above the vegetation canopy (Pa)
 real(dp),intent(in)           :: VPair                 ! vapor pressure of the air above the vegetation canopy (Pa)
 ! input: latent heat of sublimation/vaporization
 real(dp),intent(in)           :: latHeatSubVapCanopy   ! latent heat of sublimation/vaporization for the vegetation canopy (J kg-1)
 real(dp),intent(in)           :: latHeatSubVapGround   ! latent heat of sublimation/vaporization for the ground surface (J kg-1)
 ! input: canopy liquid and canopy ice (used as a solution constraint)
 real(dp),intent(in)           :: canopyLiquid          ! mass of liquid water on the vegetation canopy (kg m-2)
 real(dp),intent(in)           :: canopyIce             ! mass of ice water on the vegetation canopy (kg m-2)
 ! input: canopy and ground temperature
 real(dp),intent(in)           :: canairTemp            ! temperature of the canopy air space (K)
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
 real(dp),intent(in)            :: dGroundResistance_dTGround       ! derivative in ground resistance w.r.t. ground temperature (s m-1 K-1)
 real(dp),intent(in)            :: dGroundResistance_dTCanopy       ! derivative in ground resistance w.r.t. canopy temperature (s m-1 K-1)
 real(dp),intent(in)            :: dGroundResistance_dTCanair       ! derivative in ground resistance w.r.t. canopy air temperature (s m-1 K-1)
 real(dp),intent(in)            :: dCanopyResistance_dTCanopy       ! derivative in canopy resistance w.r.t. canopy temperature (s m-1 K-1)
 real(dp),intent(in)            :: dCanopyResistance_dTCanair       ! derivative in canopy resistance w.r.t. canopy air temperature (s m-1 K-1)
 ! ---------------------------------------------------------------------------------------------------------------------------------------------------------------
 ! output: conductances -- used to test derivatives
 real(dp),intent(out)          :: leafConductance              ! leaf conductance (m s-1)
 real(dp),intent(out)          :: canopyConductance            ! canopy conductance (m s-1)
 real(dp),intent(out)          :: groundConductanceSH          ! ground conductance for sensible heat (m s-1)
 real(dp),intent(out)          :: groundConductanceLH          ! ground conductance for latent heat -- includes soil resistance (m s-1)
 real(dp),intent(out)          :: evapConductance              ! conductance for evaporation (m s-1)
 real(dp),intent(out)          :: transConductance             ! conductance for transpiration (m s-1)
 real(dp),intent(out)          :: totalConductanceSH           ! total conductance for sensible heat (m s-1)
 real(dp),intent(out)          :: totalConductanceLH           ! total conductance for latent heat (m s-1)
 ! output: canopy air space variables
 real(dp),intent(out)          :: VP_CanopyAir                 ! vapor pressure of the canopy air space (Pa)
 ! output: fluxes from the vegetation canopy
 real(dp),intent(out)          :: senHeatCanopy                ! sensible heat flux from the canopy to the canopy air space (W m-2)
 real(dp),intent(out)          :: latHeatCanopyEvap            ! latent heat flux associated with evaporation from the canopy to the canopy air space (W m-2)
 real(dp),intent(out)          :: latHeatCanopyTrans           ! latent heat flux associated with transpiration from the canopy to the canopy air space (W m-2)
 ! output: fluxes from non-vegetated surfaces (ground surface below vegetation, bare ground, or snow covered vegetation)
 real(dp),intent(out)          :: senHeatGround                ! sensible heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
 real(dp),intent(out)          :: latHeatGround                ! latent heat flux from ground surface below vegetation, bare ground, or snow covered vegetation (W m-2)
 ! output: total heat fluxes to the atmosphere
 real(dp),intent(out)          :: senHeatTotal                 ! total sensible heat flux to the atmosphere (W m-2)
 real(dp),intent(out)          :: latHeatTotal                 ! total latent heat flux to the atmosphere (W m-2)
 ! output: net fluxes
 real(dp),intent(out)          :: turbFluxCanair               ! net turbulent heat fluxes at the canopy air space (W m-2)
 real(dp),intent(out)          :: turbFluxCanopy               ! net turbulent heat fluxes at the canopy (W m-2)
 real(dp),intent(out)          :: turbFluxGround               ! net turbulent heat fluxes at the ground surface (W m-2)
 ! output: flux derivatives
 real(dp),intent(out)          :: dTurbFluxCanair_dTCanair     ! derivative in net canopy air space fluxes w.r.t. canopy air temperature (W m-2 K-1)
 real(dp),intent(out)          :: dTurbFluxCanair_dTCanopy     ! derivative in net canopy air space fluxes w.r.t. canopy temperature (W m-2 K-1)
 real(dp),intent(out)          :: dTurbFluxCanair_dTGround     ! derivative in net canopy air space fluxes w.r.t. ground temperature (W m-2 K-1)
 real(dp),intent(out)          :: dTurbFluxCanopy_dTCanair     ! derivative in net canopy turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
 real(dp),intent(out)          :: dTurbFluxCanopy_dTCanopy     ! derivative in net canopy turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
 real(dp),intent(out)          :: dTurbFluxCanopy_dTGround     ! derivative in net canopy turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
 real(dp),intent(out)          :: dTurbFluxGround_dTCanair     ! derivative in net ground turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
 real(dp),intent(out)          :: dTurbFluxGround_dTCanopy     ! derivative in net ground turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
 real(dp),intent(out)          :: dTurbFluxGround_dTGround     ! derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                   ! error code
 character(*),intent(out)      :: message               ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! local variables -- general
 real(dp)                      :: fpart1,fpart2                ! different parts of a function
 real(dp)                      :: dpart1,dpart2                ! derivatives for different parts of a function
 real(dp),parameter            :: evapSmooth=1._dp             ! smoothing parameter for latent heat (W m-2)
 real(dp)                      :: damping                      ! damping factor used to maintain feasible solution (-)
 real(dp)                      :: maxFlux                      ! maximum possible latent heat flux (W m-2) 
 real(dp)                      :: nrgDiff                      ! initial change in energy associated with latent heat (W m-2)
 real(dp)                      :: corDiff                      ! corrected change in energy after imposing solution constraints (W m-2)
 ! local variables -- "constants"
 real(dp)                      :: volHeatCapacityAir           ! volumetric heat capacity of air (J m-3)
 real(dp)                      :: latentHeatConstant           ! latent heat constant (kg m-3 K-1)
 ! local variables -- derivatives for energy conductances
 real(dp)                      :: dCanopyCond_dCanairTemp      ! derivative in canopy conductance w.r.t. canopy air temperature
 real(dp)                      :: dCanopyCond_dCanopyTemp      ! derivative in canopy conductance w.r.t. canopy temperature
 real(dp)                      :: dGroundCondSH_dCanairTemp    ! derivative in ground conductance of sensible heat w.r.t. canopy air temperature
 real(dp)                      :: dGroundCondSH_dCanopyTemp    ! derivative in ground conductance of sensible heat w.r.t. canopy temperature
 real(dp)                      :: dGroundCondSH_dGroundTemp    ! derivative in ground conductance of sensible heat w.r.t. ground temperature
 ! local variables -- derivatives for mass conductances 
 real(dp)                      :: dGroundCondLH_dCanairTemp    ! derivative in ground conductance w.r.t. canopy air temperature
 real(dp)                      :: dGroundCondLH_dCanopyTemp    ! derivative in ground conductance w.r.t. canopy temperature
 real(dp)                      :: dGroundCondLH_dGroundTemp    ! derivative in ground conductance w.r.t. ground temperature
 ! local variables -- derivatives for the canopy air space variables
 real(dp)                      :: fPart_VP                     ! part of the function for vapor pressure of the canopy air space
 real(dp)                      :: dVPCanopyAir_dTCanair        ! derivative in the vapor pressure of the canopy air space w.r.t. temperature of the canopy air space
 real(dp)                      :: dVPCanopyAir_dTCanopy        ! derivative in the vapor pressure of the canopy air space w.r.t. temperature of the canopy 
 real(dp)                      :: dVPCanopyAir_dTGround        ! derivative in the vapor pressure of the canopy air space w.r.t. temperature of the ground
 ! local variables -- sensible heat flux derivatives 
 real(dp)                      :: dSenHeatTotal_dTCanair       ! derivative in the total sensible heat flux w.r.t. canopy air temperature 
 real(dp)                      :: dSenHeatTotal_dTCanopy       ! derivative in the total sensible heat flux w.r.t. canopy air temperature
 real(dp)                      :: dSenHeatTotal_dTGround       ! derivative in the total sensible heat flux w.r.t. ground temperature
 real(dp)                      :: dSenHeatCanopy_dTCanair      ! derivative in the canopy sensible heat flux w.r.t. canopy air temperature
 real(dp)                      :: dSenHeatCanopy_dTCanopy      ! derivative in the canopy sensible heat flux w.r.t. canopy temperature
 real(dp)                      :: dSenHeatCanopy_dTGround      ! derivative in the canopy sensible heat flux w.r.t. ground temperature
 real(dp)                      :: dSenHeatGround_dTCanair      ! derivative in the ground sensible heat flux w.r.t. canopy air temperature
 real(dp)                      :: dSenHeatGround_dTCanopy      ! derivative in the ground sensible heat flux w.r.t. canopy temperature
 real(dp)                      :: dSenHeatGround_dTGround      ! derivative in the ground sensible heat flux w.r.t. ground temperature
 ! local variables -- latent heat flux derivatives
 real(dp)                      :: dLatHeatCanopyEvap_dTCanair  ! derivative in the canopy evaporation flux w.r.t. canopy air temperature
 real(dp)                      :: dLatHeatCanopyEvap_dTCanopy  ! derivative in the canopy evaporation flux w.r.t. canopy temperature
 real(dp)                      :: dLatHeatCanopyEvap_dTGround  ! derivative in the canopy evaporation flux w.r.t. ground temperature
 real(dp)                      :: dLatHeatCanopyTrans_dTCanair ! derivative in the canopy transpiration flux w.r.t. canopy air temperature
 real(dp)                      :: dLatHeatCanopyTrans_dTCanopy ! derivative in the canopy transpiration flux w.r.t. canopy temperature
 real(dp)                      :: dLatHeatCanopyTrans_dTGround ! derivative in the canopy transpiration flux w.r.t. ground temperature
 real(dp)                      :: dLatHeatGround_dTCanair      ! derivative in the ground latent heat flux w.r.t. canopy air temperature
 real(dp)                      :: dLatHeatGround_dTCanopy      ! derivative in the ground latent heat flux w.r.t. canopy temperature
 real(dp)                      :: dLatHeatGround_dTGround      ! derivative in the ground latent heat flux w.r.t. ground temperature
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='turbFluxes/'

 ! compute constants
 volHeatCapacityAir = iden_air*cp_air           ! volumetric heat capacity of air (J m-3)
 latentHeatConstant = iden_air*w_ratio/airpres  ! latent heat constant for (kg m-3 Pa-1)

 ! *****
 ! * compute conductances, and derivatives...
 ! ******************************************

 ! compute conductances for sensible heat (m s-1)
 if(computeVegFlux)then
  leafConductance    = exposedVAI/leafResistance
  canopyConductance  = 1._dp/canopyResistance
 else
  leafConductance    = 0._dp
  canopyConductance  = 0._dp
 endif
 groundConductanceSH = 1._dp/groundResistance

 ! compute total conductance for sensible heat
 totalConductanceSH  = leafConductance + groundConductanceSH + canopyConductance

 ! compute conductances for latent heat (m s-1)
 if(computeVegFlux)then
  evapConductance    = canopyWetFraction*leafConductance
  transConductance   = (1._dp - canopyWetFraction) * ( canopySunlitLAI/(leafResistance+stomResistSunlit) + canopyShadedLAI/(leafResistance+stomResistShaded) )
  !write(*,'(a,10(f20.10,1x))') 'canopySunlitLAI, canopyShadedLAI, stomResistSunlit, stomResistShaded, leafResistance, canopyWetFraction = ', &
  !                              canopySunlitLAI, canopyShadedLAI, stomResistSunlit, stomResistShaded, leafResistance, canopyWetFraction
 else
  evapConductance    = 0._dp
  transConductance   = 0._dp
 endif
 groundConductanceLH = 1._dp/(groundResistance + soilResistance)  ! NOTE: soilResistance accounts for fractional snow, and =0 when snow cover is 100%
 totalConductanceLH  = evapConductance + transConductance + groundConductanceLH + canopyConductance

 ! * compute derivatives
 ! NOTE: it may be more efficient to compute these derivatives when computing resistances
 if(ixDerivMethod == analytical)then

  ! compute derivatives in individual conductances for sensible heat w.r.t. canopy temperature (m s-1 K-1)
  if(computeVegFlux)then
   dCanopyCond_dCanairTemp   = -dCanopyResistance_dTCanair/canopyResistance**2._dp         ! derivative in canopy conductance w.r.t. canopy air emperature
   dCanopyCond_dCanopyTemp   = -dCanopyResistance_dTCanopy/canopyResistance**2._dp         ! derivative in canopy conductance w.r.t. canopy temperature
   dGroundCondSH_dCanairTemp = -dGroundResistance_dTCanair/groundResistance**2._dp         ! derivative in ground conductance w.r.t. canopy air temperature
   dGroundCondSH_dCanopyTemp = -dGroundResistance_dTCanopy/groundResistance**2._dp         ! derivative in ground conductance w.r.t. canopy temperature
   dGroundCondSH_dGroundTemp = -dGroundResistance_dTGround/groundResistance**2._dp         ! derivative in ground conductance w.r.t. ground temperature
  else
   dCanopyCond_dCanairTemp   = 0._dp  ! derivative in canopy conductance w.r.t. canopy air emperature
   dCanopyCond_dCanopyTemp   = 0._dp  ! derivative in canopy conductance w.r.t. canopy temperature
   dGroundCondSH_dCanairTemp = 0._dp  ! derivative in ground conductance w.r.t. canopy air temperature
   dGroundCondSH_dCanopyTemp = 0._dp  ! derivative in ground conductance w.r.t. canopy temperature
   dGroundCondSH_dGroundTemp = -dGroundResistance_dTGround/groundResistance**2._dp         ! derivative in ground conductance w.r.t. ground temperature 
  endif

  ! compute derivatives in individual conductances for latent heat w.r.t. canopy temperature (m s-1 K-1)
  if(computeVegFlux)then
   dGroundCondLH_dCanairTemp = -dGroundResistance_dTCanair/(groundResistance+soilResistance)**2._dp ! derivative in ground conductance w.r.t. canopy air temperature
   dGroundCondLH_dCanopyTemp = -dGroundResistance_dTCanopy/(groundResistance+soilResistance)**2._dp ! derivative in ground conductance w.r.t. canopy temperature
   dGroundCondLH_dGroundTemp = -dGroundResistance_dTGround/(groundResistance+soilResistance)**2._dp ! derivative in ground conductance w.r.t. ground temperature
  else
   dGroundCondLH_dCanairTemp = 0._dp  ! derivative in ground conductance w.r.t. canopy air temperature
   dGroundCondLH_dCanopyTemp = 0._dp  ! derivative in ground conductance w.r.t. canopy temperature
   dGroundCondLH_dGroundTemp = -dGroundResistance_dTGround/(groundResistance+soilResistance)**2._dp ! derivative in ground conductance w.r.t. ground temperature
  endif

 endif ! (if computing analytical derivatives)

 ! *****
 ! * compute sensible and latent heat fluxes, and derivatives...
 ! *************************************************************

 ! * compute sensible and latent heat fluxes from the canopy to the canopy air space (W m-2)
 if(computeVegFlux)then

  ! compute the vapor pressure in the canopy air space (Pa)
  fPart_VP     = canopyConductance*VPair + (evapConductance + transConductance)*satVP_CanopyTemp + groundConductanceLH*satVP_GroundTemp*soilRelHumidity
  VP_CanopyAir = fPart_VP/totalConductanceLH
  !write(*,'(a,10(f20.10,1x))') 'canopyConductance, evapConductance, transConductance, groundConductanceLH, soilRelHumidity = ', &
  !                              canopyConductance, evapConductance, transConductance, groundConductanceLH, soilRelHumidity

  ! compute sensible heat flux from the canopy air space to the atmosphere
  ! NOTE: canairTemp is a state variable
  senHeatTotal = -volHeatCapacityAir*canopyConductance*(canairTemp - airtemp)
  !print*, 'canairTemp, airtemp, senHeatTotal = ', canairTemp, airtemp, senHeatTotal

  ! compute fluxes
  senHeatCanopy      = -volHeatCapacityAir*leafConductance*(canopyTemp - canairTemp)        ! (positive downwards)
  latHeatCanopyEvap  = -latHeatSubVapCanopy*latentHeatConstant*evapConductance*(satVP_CanopyTemp - VP_CanopyAir)    ! (positive downwards)
  latHeatCanopyTrans =              -LH_vap*latentHeatConstant*transConductance*(satVP_CanopyTemp - VP_CanopyAir)   ! (positive downwards)
  !write(*,'(a,10(f20.10,1x))') 'latHeatCanopyTrans, VP_CanopyAir = ', latHeatCanopyTrans, VP_CanopyAir
  
  ! check that energy for canopy evaporation does not exhaust the available water
  ! NOTE: do this here, rather than enforcing solution constraints, because energy and mass solutions may be uncoupled
  if(latHeatSubVapCanopy > LH_vap+verySmall)then ! (sublimation)
   maxFlux = -canopyIce*LH_sub/dt       ! W m-2
  else ! (evaporation)
   maxFlux = -canopyLiquid*LH_vap/dt    ! W m-2
  endif
  ! NOTE: fluxes are positive downwards
  if(latHeatCanopyEvap < maxFlux) latHeatCanopyEvap = maxFlux
  !write(*,'(a,10(f20.10,1x))') 'maxFlux, latHeatCanopyEvap = ', maxFlux, latHeatCanopyEvap

 ! * no vegetation, so fluxes are zero
 else
  senHeatCanopy      = 0._dp
  latHeatCanopyEvap  = 0._dp
  latHeatCanopyTrans = 0._dp
 endif

 ! compute sensible and latent heat fluxes from the ground to the canopy air space (W m-2)
 if(computeVegFlux)then
  senHeatGround      = -volHeatCapacityAir*groundConductanceSH*(groundTemp - canairTemp)                                          ! (positive downwards)
  latHeatGround      = -latHeatSubVapGround*latentHeatConstant*groundConductanceLH*(satVP_GroundTemp*soilRelHumidity - VP_CanopyAir)  ! (positive downwards)
 else
  senHeatGround      = -volHeatCapacityAir*groundConductanceSH*(groundTemp - airtemp)                                                 ! (positive downwards)
  latHeatGround      = -latHeatSubVapGround*latentHeatConstant*groundConductanceLH*(satVP_GroundTemp*soilRelHumidity - VPair)         ! (positive downwards)
  senHeatTotal       = senHeatGround
 endif

 ! compute latent heat flux from the canopy air space to the atmosphere
 ! NOTE: VP_CanopyAir is a diagnostic variable
 latHeatTotal = latHeatCanopyEvap + latHeatCanopyTrans + latHeatGround

 ! * compute derivatives
 if(ixDerivMethod == analytical)then

  ! differentiate CANOPY fluxes
  if(computeVegFlux)then

   ! compute derivatives of vapor pressure in the canopy air space w.r.t. all state variables
   ! (derivative of vapor pressure in the canopy air space w.r.t. temperature of the canopy air space)
   dPart1 = dCanopyCond_dCanairTemp*VPair + dGroundCondLH_dCanairTemp*satVP_GroundTemp*soilRelHumidity
   dPart2 = -(dCanopyCond_dCanairTemp + dGroundCondLH_dCanairTemp)/(totalConductanceLH**2._dp)
   dVPCanopyAir_dTCanair = dPart1/totalConductanceLH + fPart_VP*dPart2
   ! (derivative of vapor pressure in the canopy air space w.r.t. temperature of the canopy)
   dPart1 = dCanopyCond_dCanopyTemp*VPair + (evapConductance + transConductance)*dSVPCanopy_dCanopyTemp + dGroundCondLH_dCanopyTemp*satVP_GroundTemp*soilRelHumidity
   dPart2 = -(dCanopyCond_dCanopyTemp + dGroundCondLH_dCanopyTemp)/(totalConductanceLH**2._dp)
   dVPCanopyAir_dTCanopy = dPart1/totalConductanceLH + fPart_VP*dPart2
   ! (derivative of vapor pressure in the canopy air space w.r.t. temperature of the ground)
   dPart1 = dGroundCondLH_dGroundTemp*satVP_GroundTemp*soilRelHumidity + groundConductanceLH*dSVPGround_dGroundTemp*soilRelHumidity
   dPart2 = -dGroundCondLH_dGroundTemp/(totalConductanceLH**2._dp)
   dVPCanopyAir_dTGround = dPart1/totalConductanceLH + fPart_VP*dPart2
   !write(*,'(a,3(f20.8,1x))') 'dVPCanopyAir_dTCanair, dVPCanopyAir_dTCanopy, dVPCanopyAir_dTGround                      = ', &
   !                            dVPCanopyAir_dTCanair, dVPCanopyAir_dTCanopy, dVPCanopyAir_dTGround

   ! sensible heat from the canopy to the atmosphere
   dSenHeatTotal_dTCanair       = -volHeatCapacityAir*canopyConductance - volHeatCapacityAir*dCanopyCond_dCanairTemp*(canairTemp - airtemp)
   dSenHeatTotal_dTCanopy       = -volHeatCapacityAir*dCanopyCond_dCanopyTemp*(canairTemp - airtemp)
   dSenHeatTotal_dTGround       = 0._dp
   !write(*,'(a,3(f20.8,1x))') 'dSenHeatTotal_dTCanair, dSenHeatTotal_dTCanopy, dSenHeatTotal_dTGround                   = ', &
   !                            dSenHeatTotal_dTCanair, dSenHeatTotal_dTCanopy, dSenHeatTotal_dTGround

   ! sensible heat from the canopy to the canopy air space
   dSenHeatCanopy_dTCanair      =  volHeatCapacityAir*leafConductance
   dSenHeatCanopy_dTCanopy      = -volHeatCapacityAir*leafConductance
   dSenHeatCanopy_dTGround      = 0._dp
   !write(*,'(a,3(f20.8,1x))') 'dSenHeatCanopy_dTCanair, dSenHeatCanopy_dTCanopy, dSenHeatCanopy_dTGround                = ', &
   !                            dSenHeatCanopy_dTCanair, dSenHeatCanopy_dTCanopy, dSenHeatCanopy_dTGround

   ! sensible heat from the ground to the canopy air space
   dSenHeatGround_dTCanair      = -volHeatCapacityAir*dGroundCondSH_dCanairTemp*(groundTemp - canairTemp) + volHeatCapacityAir*groundConductanceSH
   dSenHeatGround_dTCanopy      = -volHeatCapacityAir*dGroundCondSH_dCanopyTemp*(groundTemp - canairTemp)
   dSenHeatGround_dTGround      = -volHeatCapacityAir*dGroundCondSH_dGroundTemp*(groundTemp - canairTemp) - volHeatCapacityAir*groundConductanceSH
   !write(*,'(a,3(f20.8,1x))') 'dSenHeatGround_dTCanair, dSenHeatGround_dTCanopy, dSenHeatGround_dTGround                = ', &
   !                            dSenHeatGround_dTCanair, dSenHeatGround_dTCanopy, dSenHeatGround_dTGround

   ! latent heat associated with canopy evaporation
   dLatHeatCanopyEvap_dTCanair  = -latHeatSubVapCanopy*latentHeatConstant*evapConductance*(-dVPCanopyAir_dTCanair)
   dLatHeatCanopyEvap_dTCanopy  = -latHeatSubVapCanopy*latentHeatConstant*evapConductance*(dSVPCanopy_dCanopyTemp - dVPCanopyAir_dTCanopy)
   dLatHeatCanopyEvap_dTGround  = -latHeatSubVapCanopy*latentHeatConstant*evapConductance*(-dVPCanopyAir_dTGround)
   !write(*,'(a,3(f20.8,1x))') 'dLatHeatCanopyEvap_dTCanair, dLatHeatCanopyEvap_dTCanopy, dLatHeatCanopyEvap_dTGround    = ', &
   !                            dLatHeatCanopyEvap_dTCanair, dLatHeatCanopyEvap_dTCanopy, dLatHeatCanopyEvap_dTGround

   ! latent heat associated with canopy transpiration
   dLatHeatCanopyTrans_dTCanair = -LH_vap*latentHeatConstant*transConductance*(-dVPCanopyAir_dTCanair)
   dLatHeatCanopyTrans_dTCanopy = -LH_vap*latentHeatConstant*transConductance*(dSVPCanopy_dCanopyTemp - dVPCanopyAir_dTCanopy)
   dLatHeatCanopyTrans_dTGround = -LH_vap*latentHeatConstant*transConductance*(-dVPCanopyAir_dTGround)
   !write(*,'(a,3(f20.8,1x))') 'dLatHeatCanopyTrans_dTCanair, dLatHeatCanopyTrans_dTCanopy, dLatHeatCanopyTrans_dTGround = ', &
   !                            dLatHeatCanopyTrans_dTCanair, dLatHeatCanopyTrans_dTCanopy, dLatHeatCanopyTrans_dTGround

   ! latent heat flux from the ground
   fPart1 = -latHeatSubVapGround*latentHeatConstant*groundConductanceLH       ! function of the first part
   fPart2 = (satVP_GroundTemp*soilRelHumidity - VP_CanopyAir)                 ! function of the second part
   dLatHeatGround_dTCanair      = -latHeatSubVapGround*latentHeatConstant*dGroundCondLH_dCanairTemp*fPart2 - dVPCanopyAir_dTCanair*fPart1
   dLatHeatGround_dTCanopy      = -latHeatSubVapGround*latentHeatConstant*dGroundCondLH_dCanopyTemp*fPart2 - dVPCanopyAir_dTCanopy*fPart1
   dLatHeatGround_dTGround      = -latHeatSubVapGround*latentHeatConstant*dGroundCondLH_dGroundTemp*fPart2 + (dSVPGround_dGroundTemp*soilRelHumidity - dVPCanopyAir_dTGround)*fPart1
   !write(*,'(a,3(f20.8,1x))') 'dLatHeatGround_dTCanair, dLatHeatGround_dTCanopy, dLatHeatGround_dTGround                = ', &
   !                            dLatHeatGround_dTCanair, dLatHeatGround_dTCanopy, dLatHeatGround_dTGround

  else  ! canopy is undefined

   ! set derivatives for canopy fluxes to zero (no canopy, so fluxes are undefined)
   dSenHeatTotal_dTCanair       = 0._dp
   dSenHeatTotal_dTCanopy       = 0._dp
   dSenHeatTotal_dTGround       = 0._dp
   dSenHeatCanopy_dTCanair      = 0._dp
   dSenHeatCanopy_dTCanopy      = 0._dp
   dSenHeatCanopy_dTGround      = 0._dp
   dLatHeatCanopyEvap_dTCanair  = 0._dp
   dLatHeatCanopyEvap_dTCanopy  = 0._dp
   dLatHeatCanopyEvap_dTGround  = 0._dp
   dLatHeatCanopyTrans_dTCanair = 0._dp
   dLatHeatCanopyTrans_dTCanopy = 0._dp
   dLatHeatCanopyTrans_dTGround = 0._dp

   ! set derivatives for ground fluxes w.r.t canopy temperature to zero (no canopy, so fluxes are undefined)
   dSenHeatGround_dTCanair = 0._dp
   dSenHeatGround_dTCanopy = 0._dp
   dLatHeatGround_dTCanair = 0._dp
   dLatHeatGround_dTCanopy = 0._dp

   ! compute derivatives for the ground fluxes w.r.t. ground temperature
   dSenHeatGround_dTGround = (-volHeatCapacityAir*dGroundCondSH_dGroundTemp)*(groundTemp - airtemp) + &                               ! d(ground sensible heat flux)/d(ground temp)
                             (-volHeatCapacityAir*groundConductanceSH)
   dLatHeatGround_dTGround = (-latHeatSubVapGround*latentHeatConstant*dGroundCondLH_dGroundTemp)*(satVP_GroundTemp - VPair) + &       ! d(ground latent heat flux)/d(ground temp)
                             (-latHeatSubVapGround*latentHeatConstant*groundConductanceLH)*dSVPGround_dGroundTemp

  endif   ! (if canopy is defined)

 endif  ! (if computing analytical derivatives)


 ! *****
 ! * compute net turbulent fluxes, and derivatives...
 ! **************************************************

 ! compute net fluxes
 turbFluxCanair = senHeatTotal - senHeatCanopy - senHeatGround            ! net turbulent flux at the canopy air space (W m-2)
 turbFluxCanopy = senHeatCanopy + latHeatCanopyEvap + latHeatCanopyTrans  ! net turbulent flux at the canopy (W m-2)
 turbFluxGround = senHeatGround + latHeatGround                           ! net turbulent flux at the ground surface (W m-2)

  ! * compute derivatives
 if(ixDerivMethod == analytical)then
  dTurbFluxCanair_dTCanair = dSenHeatTotal_dTCanair - dSenHeatCanopy_dTCanair - dSenHeatGround_dTCanair            ! derivative in net canopy air space fluxes w.r.t. canopy air temperature (W m-2 K-1)
  dTurbFluxCanair_dTCanopy = dSenHeatTotal_dTCanopy - dSenHeatCanopy_dTCanopy - dSenHeatGround_dTCanopy            ! derivative in net canopy air space fluxes w.r.t. canopy temperature (W m-2 K-1)
  dTurbFluxCanair_dTGround = dSenHeatTotal_dTGround - dSenHeatCanopy_dTGround - dSenHeatGround_dTGround            ! derivative in net canopy air space fluxes w.r.t. ground temperature (W m-2 K-1)
  dTurbFluxCanopy_dTCanair = dSenHeatCanopy_dTCanair + dLatHeatCanopyEvap_dTCanair + dLatHeatCanopyTrans_dTCanair  ! derivative in net canopy turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
  dTurbFluxCanopy_dTCanopy = dSenHeatCanopy_dTCanopy + dLatHeatCanopyEvap_dTCanopy + dLatHeatCanopyTrans_dTCanopy  ! derivative in net canopy turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
  dTurbFluxCanopy_dTGround = dSenHeatCanopy_dTGround + dLatHeatCanopyEvap_dTGround + dLatHeatCanopyTrans_dTGround  ! derivative in net canopy turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
  dTurbFluxGround_dTCanair = dSenHeatGround_dTCanair + dLatHeatGround_dTCanair                                     ! derivative in net ground turbulent fluxes w.r.t. canopy air temperature (W m-2 K-1)
  dTurbFluxGround_dTCanopy = dSenHeatGround_dTCanopy + dLatHeatGround_dTCanopy                                     ! derivative in net ground turbulent fluxes w.r.t. canopy temperature (W m-2 K-1)
  dTurbFluxGround_dTGround = dSenHeatGround_dTGround + dLatHeatGround_dTGround                                     ! derivative in net ground turbulent fluxes w.r.t. ground temperature (W m-2 K-1)
 else ! (just make sure we return something)
  dTurbFluxCanair_dTCanair = 0._dp 
  dTurbFluxCanair_dTCanopy = 0._dp 
  dTurbFluxCanair_dTGround = 0._dp 
  dTurbFluxCanopy_dTCanair = 0._dp 
  dTurbFluxCanopy_dTCanopy = 0._dp 
  dTurbFluxCanopy_dTGround = 0._dp
  dTurbFluxGround_dTCanair = 0._dp 
  dTurbFluxGround_dTCanopy = 0._dp 
  dTurbFluxGround_dTGround = 0._dp
 endif
 
 end subroutine turbFluxes


 ! ***********************************************************************************************************
 ! private subroutine: compute stability corrections for turbulent heat fluxes (-)
 ! *********************************************************************************************************** 
 subroutine aStability(&
                       ! input: control
                       computeDerivative,              & ! input: logical flag to compute analytical derivatives
                       ixStability,                    & ! input: choice of stability function
                       ! input: forcing data, diagnostic and state variables
                       mHeight,                        & ! input: measurement height (m)
                       airTemp,                        & ! input: air temperature (K)
                       sfcTemp,                        & ! input: surface temperature (K)
                       windspd,                        & ! input: wind speed (m s-1)
                       ! input: stability parameters
                       critRichNumber,                 & ! input: critical value for the bulk Richardson number where turbulence ceases (-)
                       Louis79_bparam,                 & ! input: parameter in Louis (1979) stability function
                       Louis79_cStar,                  & ! input: parameter in Louis (1979) stability function
                       Mahrt87_eScale,                 & ! input: exponential scaling factor in the Mahrt (1987) stability function
                       ! output
                       RiBulk,                         & ! output: bulk Richardson number (-)
                       stabilityCorrection,            & ! output: stability correction for turbulent heat fluxes (-)
                       dStabilityCorrection_dRich,     & ! output: derivative in stability correction w.r.t. Richardson number (-)
                       dStabilityCorrection_dAirTemp,  & ! output: derivative in stability correction w.r.t. temperature (K-1)
                       dStabilityCorrection_dSfcTemp,  & ! output: derivative in stability correction w.r.t. temperature (K-1)
                       err, message                    ) ! output: error control
 implicit none
 ! input: control
 logical(lgt),intent(in)       :: computeDerivative      ! flag to compute the derivative
 integer(i4b),intent(in)       :: ixStability            ! choice of stability function
 ! input: forcing data, diagnostic and state variables
 real(dp),intent(in)           :: mHeight                ! measurement height (m)
 real(dp),intent(in)           :: airtemp                ! air temperature (K)
 real(dp),intent(in)           :: sfcTemp                ! surface temperature (K)
 real(dp),intent(in)           :: windspd                ! wind speed (m s-1)
 ! input: stability parameters
 real(dp),intent(in)           :: critRichNumber         ! critical value for the bulk Richardson number where turbulence ceases (-)
 real(dp),intent(in)           :: Louis79_bparam         ! parameter in Louis (1979) stability function
 real(dp),intent(in)           :: Louis79_cStar          ! parameter in Louis (1979) stability function
 real(dp),intent(in)           :: Mahrt87_eScale         ! exponential scaling factor in the Mahrt (1987) stability function
 ! output
 real(dp),intent(out)          :: RiBulk                 ! bulk Richardson number (-)
 real(dp),intent(out)          :: stabilityCorrection    ! stability correction for turbulent heat fluxes (-) 
 real(dp),intent(out)          :: dStabilityCorrection_dRich    ! derivative in stability correction w.r.t. Richardson number (-)
 real(dp),intent(out)          :: dStabilityCorrection_dAirTemp ! derivative in stability correction w.r.t. air temperature (K-1)
 real(dp),intent(out)          :: dStabilityCorrection_dSfcTemp ! derivative in stability correction w.r.t. surface temperature (K-1)
 integer(i4b),intent(out)      :: err                    ! error code
 character(*),intent(out)      :: message                ! error message
 ! local
 real(dp)                      :: dRiBulk_dAirTemp       ! derivative in the bulk Richardson number w.r.t. air temperature (K-1)
 real(dp)                      :: dRiBulk_dSfcTemp       ! derivative in the bulk Richardson number w.r.t. surface temperature (K-1)
 real(dp)                      :: bPrime                 ! scaled "b" parameter for stability calculations in Louis (1979)
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message='aStability/'

 ! compute the bulk Richardson number (-)
 call bulkRichardson(&
                     ! input
                     airTemp,                        & ! input: air temperature (K)
                     sfcTemp,                        & ! input: surface temperature (K)
                     windspd,                        & ! input: wind speed (m s-1)
                     mHeight,                        & ! input: measurement height (m)
                     computeDerivative,              & ! input: flag to compute the derivative
                     ! output
                     RiBulk,                         & ! output: bulk Richardson number (-)
                     dRiBulk_dAirTemp,               & ! output: derivative in the bulk Richardson number w.r.t. air temperature (K-1)
                     dRiBulk_dSfcTemp,               & ! output: derivative in the bulk Richardson number w.r.t. surface temperature (K-1)
                     err,message)                      ! output: error control

 ! set derivative to one if not computing it
 if(.not.computeDerivative)then
  dStabilityCorrection_dRich    = 1._dp
  dStabilityCorrection_dAirTemp = 1._dp
  dStabilityCorrection_dSfcTemp = 1._dp
 endif

 ! ***** process unstable cases
 if(RiBulk<0._dp)then
  ! compute surface-atmosphere exchange coefficient (-)
  stabilityCorrection = (1._dp - 16._dp*RiBulk)**0.5_dp
  ! compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
  if(computeDerivative)then
   dStabilityCorrection_dRich    = (-16._dp) * 0.5_dp*(1._dp - 16._dp*RiBulk)**(-0.5_dp)
   dStabilityCorrection_dAirTemp = dRiBulk_dAirTemp * dStabilityCorrection_dRich
   dStabilityCorrection_dSfcTemp = dRiBulk_dSfcTemp * dStabilityCorrection_dRich
  endif
  return
 endif

 ! ***** process stable cases
 select case(ixStability)

  ! ("standard" stability correction, a la Anderson 1976)
  case(standard)
   ! compute surface-atmosphere exchange coefficient (-)   
   if(RiBulk <  critRichNumber) stabilityCorrection = (1._dp - 5._dp*RiBulk)**2._dp
   if(RiBulk >= critRichNumber) stabilityCorrection = epsilon(stabilityCorrection)
   ! compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
   if(computeDerivative)then
    if(RiBulk <  critRichNumber) dStabilityCorrection_dRich = (-5._dp) * 2._dp*(1._dp - 5._dp*RiBulk)
    if(RiBulk >= critRichNumber) dStabilityCorrection_dRich = 0._dp
   endif

  ! (Louis 1979)
  case(louisInversePower)
   ! scale the "b" parameter for stable conditions
   bprime = Louis79_bparam/2._dp
   ! compute surface-atmosphere exchange coefficient (-)
   stabilityCorrection = 1._dp / ( (1._dp + bprime*RiBulk)**2._dp )
   if(stabilityCorrection < epsilon(stabilityCorrection)) stabilityCorrection = epsilon(stabilityCorrection)
   ! compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
   if(computeDerivative)then
    dStabilityCorrection_dRich = bprime * (-2._dp)*(1._dp + bprime*RiBulk)**(-3._dp)
   endif

  ! (Mahrt 1987)
  case(mahrtExponential)
   ! compute surface-atmosphere exchange coefficient (-)
   stabilityCorrection = exp(-Mahrt87_eScale * RiBulk)
   if(stabilityCorrection < epsilon(stabilityCorrection)) stabilityCorrection = epsilon(stabilityCorrection)
   ! compute derivative in surface-atmosphere exchange coefficient w.r.t. temperature (K-1)
   if(computeDerivative)then
    dStabilityCorrection_dRich = (-Mahrt87_eScale) * exp(-Mahrt87_eScale * RiBulk)
   endif

  ! (return error if the stability correction method is not found)
  case default
   err=10; message=trim(message)//"optionNotFound[stability correction]"; return

 endselect
 
 ! get the stability correction with respect to air temperature and surface temperature
 ! NOTE: air temperature is used for canopy air temperature, which is a model state variable
 if(computeDerivative)then
  dStabilityCorrection_dAirTemp = dRiBulk_dAirTemp * dStabilityCorrection_dRich
  dStabilityCorrection_dSfcTemp = dRiBulk_dSfcTemp * dStabilityCorrection_dRich
 endif

 end subroutine aStability


 ! ***********************************************************************************************************
 ! private subroutine: compute bulk Richardson number
 ! *********************************************************************************************************** 
 subroutine bulkRichardson(&
                           ! input
                           airTemp,                    & ! input: air temperature (K)
                           sfcTemp,                    & ! input: surface temperature (K)
                           windspd,                    & ! input: wind speed (m s-1)
                           mHeight,                    & ! input: measurement height (m)
                           computeDerivative,          & ! input: flag to compute the derivative
                           ! output
                           RiBulk,                     & ! output: bulk Richardson number (-)
                           dRiBulk_dAirTemp,           & ! output: derivative in the bulk Richardson number w.r.t. air temperature (K-1)
                           dRiBulk_dSfcTemp,           & ! output: derivative in the bulk Richardson number w.r.t. surface temperature (K-1)
                           err,message)                  ! output: error control
 implicit none
 ! input
 real(dp),intent(in)           :: airtemp                ! air temperature (K)
 real(dp),intent(in)           :: sfcTemp                ! surface temperature (K)
 real(dp),intent(in)           :: windspd                ! wind speed (m s-1)
 real(dp),intent(in)           :: mHeight                ! measurement height (m)
 logical(lgt),intent(in)       :: computeDerivative      ! flag to compute the derivative
 ! output
 real(dp),intent(inout)        :: RiBulk                 ! bulk Richardson number (-)
 real(dp),intent(out)          :: dRiBulk_dAirTemp       ! derivative in the bulk Richardson number w.r.t. air temperature (K-1)
 real(dp),intent(out)          :: dRiBulk_dSfcTemp       ! derivative in the bulk Richardson number w.r.t. surface temperature (K-1)
 integer(i4b),intent(out)      :: err                    ! error code
 character(*),intent(out)      :: message                ! error message
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
  dRiBulk_dAirTemp =  RiMult/T_mean - RiMult*T_grad/(0.5_dp*((airtemp + sfcTemp)**2._dp))
  dRiBulk_dSfcTemp = -RiMult/T_mean - RiMult*T_grad/(0.5_dp*((airtemp + sfcTemp)**2._dp))
 else
  dRiBulk_dAirTemp = 1._dp
  dRiBulk_dSfcTemp = 1._dp
 endif
 end subroutine bulkRichardson



 ! ****************** EXPONENTIAL INTEGRAL FUNCTION *****************************************
 ! From UEBVeg
 ! Computes the exponential integral function for the given value
 FUNCTION EXPINT (LAI)
 REAL(DP) LAI
 REAL(DP) EXPINT
 REAL(DP) a0,a1,a2,a3,a4,a5,b1,b2,b3,b4
 IF (LAI.EQ.0)THEN
  EXPINT=1._dp

 ELSEIF (LAI.LE.1.0) THEN
  a0=-.57721566_dp
  a1=.99999193_dp
  a2=-.24991055_dp
  a3=.05519968_dp
  a4=-.00976004_dp
  a5=.00107857_dp

  EXPINT = a0+a1*LAI+a2*LAI**2+a3*LAI**3+a4*LAI**4+a5*LAI**5 - log(LAI)

 ELSE
  a1=8.5733287401_dp
  a2=18.0590169730_dp
  a3=8.6347637343_dp
  a4=.2677737343_dp
  b1=9.5733223454_dp
  b2=25.6329561486_dp
  b3=21.0996530827_dp
  b4=3.9584969228_dp

  EXPINT=(LAI**4+a1*LAI**3+a2*LAI**2+a3*LAI+a4)/ &
      ((LAI**4+b1*LAI**3+b2*LAI**2+b3*LAI+b4)*LAI*exp(LAI))

 END IF
 RETURN
 END FUNCTION EXPINT





end module vegNrgFlux_module

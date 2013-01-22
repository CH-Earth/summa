module noahMP_veg_module
USE nrtype
! constants
USE multiconst,only:gravity    ! acceleration of gravity (m s-2)
USE multiconst,only:vkc        ! von Karman's constant (-)
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
USE data_struc,only:ix_soil,ix_snow        ! named variables for snow and soil
! -------------------------------------------------------------------------------------------------
implicit none
private
public::noahMP_veg
! number of soil and snow layers
integer(i4b)                  :: nSoil     ! number of soil layers
integer(i4b)                  :: nSnow     ! number of snow layers
integer(i4b)                  :: nLayers   ! total number of layers
! algorithmic parameters
real(dp),parameter     :: mpe=1.e-6_dp     ! prevents overflow error if division by zero 
real(dp),parameter     :: dx=1.e-8_dp      ! finite difference increment
contains


 ! ************************************************************************************************
 ! new subroutine: compute energy and mass fluxes for vegetation
 ! ************************************************************************************************
 subroutine noahMP_veg(dt,                                  & ! intent(in): time step (seconds)
                       iter,                                & ! intent(in): iteration index
                       err,message)                           ! intent(out): error control
 ! constants
 USE multiconst,only:Tfreeze                                  ! temperature at freezing              (K)
 ! model decisions
 USE data_struc,only:model_decisions                          ! model decision structure
 USE var_lookup,only:iLookDECISIONS                           ! named variables for elements of the decision structure
 ! model variables, parameters, etc.
 USE data_struc,only:time_data,type_data,attr_data,forc_data,mpar_data,mvar_data,indx_data     ! data structures
 USE var_lookup,only:iLookTIME,iLookTYPE,iLookATTR,iLookFORCE,iLookPARAM,iLookMVAR,iLookINDEX  ! named variables for structure elements
 ! conversion functions
 USE conv_funcs_module,only:satVapPress                       ! compute saturated vapor pressure (Pa)
 ! Noah-MP modules
 USE NOAHMP_ROUTINES,only:phenology                           ! compute vegetation phenology
 USE NOAHMP_ROUTINES,only:radiation                           ! compute radiation fluxes
 ! compute energy and mass fluxes for vegetation
 implicit none
 ! input
 real(dp),intent(in)           :: dt                          ! time step (seconds)
 integer(i4b),intent(in)       :: iter                        ! iteration index
 ! output
 integer(i4b),intent(out)      :: err                         ! error code
 character(*),intent(out)      :: message                     ! error message
 ! local (general)
 character(LEN=256)            :: cmessage                    ! error message of downwind routine
 real(dp),parameter            :: verySmall=epsilon(1._dp)    ! a very small number
 integer(i4b),parameter        :: ist     = 1                 ! Surface type:  IST=1 => soil;  IST=2 => lake
 integer(i4b),parameter        :: isc     = 4                 ! Soil color type
 integer(i4b),parameter        :: ice     = 0                 ! Surface type:  ICE=0 => soil;  ICE=1 => sea-ice
 integer(i4b),parameter        :: nBands  = 2                 ! number of spectral bands for shortwave radiation
 integer(i4b),parameter        :: iLoc    = 1                 ! i-location
 integer(i4b),parameter        :: jLoc    = 1                 ! j-location
 real(dp)                      :: snowmassPlusNewsnow         ! sum of snow mass and new snowfall (kg m-2 [mm])
 real(dp)                      :: fracSnow                    ! snow cover fraction (0-1)
 real(dp)                      :: VAI                         ! vegetation area index (m2 m-2)
 real(dp)                      :: exposedVAI                  ! "exposed" vegetation area index (m2 m-2)
 real(dp)                      :: greenVegFraction            ! green vegetation fraction (0-1) 
 ! local (phenology)
 integer(i4b)                  :: nLayersRoots                ! number of soil layers that contain roots
 ! local (saturation vapor pressure of veg)
 real(dp)                      :: TV_celcius                  ! vegetaion temperature (C)
 real(dp)                      :: dSVP_dTV                    ! derivative in saturated vapor pressure w.r.t. vegetation temperature (Pa/K)
 ! local (longwave radiation)
 real(dp)                      :: canopyEmissivity            ! effective emissivity of the canopy (-)
 real(dp)                      :: groundEmissivity            ! emissivity of the ground surface (-)
 real(dp),parameter            :: soilEmissivity=0.98_dp      ! emmisivity of the soil (-)
 real(dp),parameter            :: snowEmissivity=0.99_dp      ! emissivity of snow (-)
 real(dp)                      :: dLWNetCanopy_dTCanopy       ! derivative in net canopy radiation w.r.t. canopy temperature (W m-2 K-1)
 real(dp)                      :: dLWNetGround_dTGround       ! derivative in net ground radiation w.r.t. ground temperature (W m-2 K-1)
 real(dp)                      :: dLWNetCanopy_dTGround       ! derivative in net canopy radiation w.r.t. ground temperature (W m-2 K-1)
 real(dp)                      :: dLWNetGround_dTCanopy       ! derivative in net ground radiation w.r.t. canopy temperature (W m-2 K-1)
 ! local (turbulent heat transfer)
 real(dp)                      :: z0Ground                    ! roughness length of the ground (ground below the canopy or non-vegetated surface) (m)

 ! initialize error control
 err=0; message="noahMP_veg/"

 ! define number of layers
 nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat == ix_soil)  ! number of soil layers
 nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat == ix_snow)  ! number of snow layers
 nLayers = indx_data%var(iLookINDEX%nLayers)%dat(1)                   ! total number of layers

 ! initialize the canopy temperature and vapor pressure
 if(iter==1)then
  mvar_data%var(iLookMVAR%scalarVP_CanopyAir)%dat(1)   = mvar_data%var(iLookMVAR%scalarVPair)%dat(1)
  mvar_data%var(iLookMVAR%scalarTemp_CanopyAir)%dat(1) = forc_data%var(iLookFORCE%airtemp)
 endif

 ! define the foliage nitrogen factor
 mvar_data%var(iLookMVAR%scalarFoliageNitrogenFactor)%dat(1) = 1._dp  ! foliage nitrogen concentration (1.0 = saturated)

 ! compute the root zone temperature
 nLayersRoots = count(mvar_data%var(iLookMVAR%iLayerHeight)%dat(nSnow:nLayers-1) < mpar_data%var(iLookPARAM%rootingDepth)-verySmall)
 if(nLayersRoots == 0)then; err=20; message=trim(message)//'no roots within the soil profile'; return; endif
 mvar_data%var(iLookMVAR%scalarRootZoneTemp)%dat(1) = sum(mvar_data%var(iLookMVAR%mLayerTemp)%dat(nSnow+1:nSnow+nLayersRoots)) / real(nLayersRoots, kind(dp))

 ! compute the sum of snow mass and new snowfall (kg m-2 [mm])
 snowmassPlusNewsnow = mvar_data%var(iLookMVAR%scalarSWE)%dat(1) + mvar_data%var(iLookMVAR%scalarSnowfall)%dat(1)*dt

 ! compute the snow cover fraction
 if(mvar_data%var(iLookMVAR%scalarSWE)%dat(1) > 0._dp)then
  fracSnow = 1._dp
 else
  fracSnow = 0._dp
 endif

 ! compute the roughness length of the ground (ground below the canopy or non-vegetated surface)
 z0Ground = mpar_data%var(iLookPARAM%z0soil)*(1._dp - fracSnow) + mpar_data%var(iLookPARAM%z0Snow)*fracSnow     ! roughness length (m)

 ! compute vegetation phenology
 call phenology(type_data%var(iLookTYPE%vegTypeIndex),                          & ! intent(in): vegetation type index
                time_data%var(iLookTIME%im),                                    & ! intent(in): month
                time_data%var(iLookTIME%id),                                    & ! intent(in): day
                mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1),                & ! intent(in): snow depth (m)
                mvar_data%var(iLookMVAR%scalarVegetationTemp)%dat(1),           & ! intent(in): vegetation temperature (K)
                attr_data%var(iLookATTR%latitude),                              & ! intent(in): latitude (degrees north)
                mvar_data%var(iLookMVAR%scalarLAI)%dat(1),                      & ! intent(inout): one-sided leaf area index (m2 m-2)
                mvar_data%var(iLookMVAR%scalarSAI)%dat(1),                      & ! intent(inout): one-sided stem area index (m2 m-2)
                mvar_data%var(iLookMVAR%scalarRootZoneTemp)%dat(1),             & ! intent(in): average temperature of the root zone (K)
                mvar_data%var(iLookMVAR%scalarCanopyHeight)%dat(1),             & ! intent(out): height of the top of the canopy layer (m)
                mvar_data%var(iLookMVAR%scalarExposedLAI)%dat(1),               & ! intent(out): exposed leaf area index after burial by snow (m2 m-2)
                mvar_data%var(iLookMVAR%scalarExposedSAI)%dat(1),               & ! intent(out): exposed stem area index after burial by snow (m2 m-2)
                mvar_data%var(iLookMVAR%scalarGrowingSeasonIndex)%dat(1))         ! intent(out): growing season index (0=off, 1=on)

 ! compute the total vegetation area index (leaf plus stem)
 VAI  = mvar_data%var(iLookMVAR%scalarLAI)%dat(1) + mvar_data%var(iLookMVAR%scalarSAI)%dat(1) ! vegetation area index
 exposedVAI = mvar_data%var(iLookMVAR%scalarExposedLAI)%dat(1) +  mvar_data%var(iLookMVAR%scalarExposedSAI)%dat(1)  !  exposed vegetation area index

 ! compute the green vegetation fraction (from Noah-MP)
 greenVegFraction = 1._dp - exp(-0.52_dp*VAI)

 ! compute emissivity of the canopy and ground surface (-)
 canopyEmissivity = 1._dp - exp(-exposedVAI)                                     ! effective emissivity of the canopy (-)
 groundEmissivity = fracSnow*snowEmissivity + (1._dp - fracSnow)*soilEmissivity  ! emissivity of the ground surface (-)

 ! compute canopy shortwave radiation fluxes
 ! (unchanged Noah-MP routine)
 call radiation(&
                ! input
                type_data%var(iLookTYPE%vegTypeIndex),                           & ! intent(in): vegetation type index
                ist, isc, ice,                                                   & ! intent(in): indices to define surface type, soil color, and ice type (constant)
                nSoil,                                                           & ! intent(in): number of soil layers               
                mvar_data%var(iLookMVAR%scalarSWE)%dat(1),                       & ! intent(in): snow water equivalent (kg m-2 [mm])
                snowmassPlusNewsnow,                                             & ! intent(in): sum of snow mass and new snowfall (kg m-2 [mm])
                dt,                                                              & ! intent(in): time step (s)
                mvar_data%var(iLookMVAR%scalarCosZenith)%dat(1),                 & ! intent(in): cosine of the solar zenith angle (0-1)
                mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1)*1000._dp,        & ! intent(in): snow depth (mm)
                mvar_data%var(iLookMVAR%mLayerTemp)%dat(1),                      & ! intent(in): ground temperature (K)
                mvar_data%var(iLookMVAR%scalarVegetationTemp)%dat(1),            & ! intent(in): vegetation temperature (K)
                fracSnow,                                                        & ! intent(in): snow cover fraction (0-1)
                mvar_data%var(iLookMVAR%scalarSnowfall)%dat(1),                  & ! intent(in): snowfall (kg m-2 s-1 [mm/s])
                mvar_data%var(iLookMVAR%scalarCanopyWetFraction)%dat(1),         & ! intent(in): fraction of canopy that is wet
                mvar_data%var(iLookMVAR%scalarExposedLAI)%dat(1),                & ! intent(in): exposed leaf area index after burial by snow (m2 m-2)
                mvar_data%var(iLookMVAR%scalarExposedSAI)%dat(1),                & ! intent(in): exposed stem area index after burial by snow (m2 m-2)          
                mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(nSnow+1:nLayers),  & ! intent(in): volumetric fraction of liquid water in each soil layer (-)
                mvar_data%var(iLookMVAR%spectralIncomingDirect)%dat(1:nBands),   & ! intent(in): incoming direct solar radiation in each wave band (w m-2)
                mvar_data%var(iLookMVAR%spectralIncomingDiffuse)%dat(1:nBands),  & ! intent(in): incoming diffuse solar radiation in each wave band (w m-2)
                greenVegFraction,                                                & ! intent(in): green vegetation fraction (0-1)
                iLoc, jLoc,                                                      & ! intent(in): spatial location indices      
                ! output
                mvar_data%var(iLookMVAR%scalarAlbedo)%dat(1),                    & ! intent(inout): snow albedo (-)
                mvar_data%var(iLookMVAR%scalarCanopySunlitFraction)%dat(1),      & ! intent(out): sunlit fraction of canopy (-)
                mvar_data%var(iLookMVAR%scalarCanopySunlitLAI)%dat(1),           & ! intent(out): sunlit leaf area (-)
                mvar_data%var(iLookMVAR%scalarCanopyShadedLAI)%dat(1),           & ! intent(out): shaded leaf area (-)
                mvar_data%var(iLookMVAR%scalarCanopySunlitPAR)%dat(1),           & ! intent(out): average absorbed par for sunlit leaves (w m-2)
                mvar_data%var(iLookMVAR%scalarCanopyShadedPAR)%dat(1),           & ! intent(out): average absorbed par for shaded leaves (w m-2)
                mvar_data%var(iLookMVAR%scalarCanopyAbsorbedSolar)%dat(1),       & ! intent(out): solar radiation absorbed by canopy (W m-2)
                mvar_data%var(iLookMVAR%scalarGroundAbsorbedSolar)%dat(1),       & ! intent(out): solar radiation absorbed by ground (W m-2)
                mvar_data%var(iLookMVAR%scalarTotalReflectedSolar)%dat(1),       & ! intent(out): total reflected solar radiation (W m-2)
                mvar_data%var(iLookMVAR%scalarTotalAbsorbedSolar)%dat(1),        & ! intent(out): total absorbed solar radiation (W m-2)
                mvar_data%var(iLookMVAR%scalarCanopyReflectedSolar)%dat(1),      & ! intent(out): solar radiation reflected from the canopy (W m-2)
                mvar_data%var(iLookMVAR%scalarGroundReflectedSolar)%dat(1),      & ! intent(out): solar radiation reflected from the ground (W m-2) 
                mvar_data%var(iLookMVAR%scalarBetweenCanopyGapFraction)%dat(1),  & ! intent(out): between canopy gap fraction for beam (-)
                mvar_data%var(iLookMVAR%scalarWithinCanopyGapFraction)%dat(1),   & ! intent(out): within canopy gap fraction for beam (-)
                mvar_data%var(iLookMVAR%scalarTotalCanopyGapFraction)%dat(1)     ) ! intent(out): total canopy gap fraction for beam (-)

  ! compute canopy longwave radiation balance
  call longwaveBal(&
                   ! input: model control
                   model_decisions(iLookDECISIONS%fDerivMeth)%iDecision,         & ! intent(in): method used to calculate flux derivatives 
                   ! input: canopy and ground temperature
                   mvar_data%var(iLookMVAR%scalarVegetationTemp)%dat(1),         & ! intent(in): canopy temperature (K)
                   mvar_data%var(iLookMVAR%mLayerTemp)%dat(1),                   & ! intent(in): ground temperature (K)
                   ! input: canopy and ground emissivity
                   canopyEmissivity,                                             & ! intent(in): canopy emissivity (-)
                   groundEmissivity,                                             & ! intent(in): ground emissivity (-)
                   ! input: forcing
                   forc_data%var(iLookFORCE%LWRadAtm),                           & ! intent(in): downwelling longwave radiation at the upper boundary (W m-2)
                   ! output: emitted radiation from the canopy and ground
                   mvar_data%var(iLookMVAR%scalarLWRadCanopy)%dat(1),            & ! intent(out): longwave radiation emitted from the canopy (W m-2)
                   mvar_data%var(iLookMVAR%scalarLWRadGround)%dat(1),            & ! intent(out): longwave radiation emitted at the ground surface (W m-2)
                   ! output: individual fluxes
                   mvar_data%var(iLookMVAR%scalarLWRadUbound2Canopy)%dat(1),     & ! intent(out): downward atmospheric longwave radiation absorbed by the canopy (W m-2)
                   mvar_data%var(iLookMVAR%scalarLWRadUbound2Ground)%dat(1),     & ! intent(out): downward atmospheric longwave radiation absorbed by the ground (W m-2)
                   mvar_data%var(iLookMVAR%scalarLWRadUbound2Ubound)%dat(1),     & ! intent(out): atmospheric radiation reflected by the ground and lost thru upper boundary (W m-2)
                   mvar_data%var(iLookMVAR%scalarLWRadCanopy2Ubound)%dat(1),     & ! intent(out): longwave radiation emitted from canopy lost thru upper boundary (W m-2)
                   mvar_data%var(iLookMVAR%scalarLWRadCanopy2Ground)%dat(1),     & ! intent(out): longwave radiation emitted from canopy absorbed by the ground (W m-2)
                   mvar_data%var(iLookMVAR%scalarLWRadCanopy2Canopy)%dat(1),     & ! intent(out): canopy longwave reflected from ground and absorbed by the canopy (W m-2)
                   mvar_data%var(iLookMVAR%scalarLWRadGround2Ubound)%dat(1),     & ! intent(out): longwave radiation emitted from ground lost thru upper boundary (W m-2)
                   mvar_data%var(iLookMVAR%scalarLWRadGround2Canopy)%dat(1),     & ! intent(out): longwave radiation emitted from ground and absorbed by the canopy (W m-2)
                   ! output: net fluxes
                   mvar_data%var(iLookMVAR%scalarLWNetCanopy)%dat(1),            & ! intent(out): net longwave radiation at the canopy (W m-2)
                   mvar_data%var(iLookMVAR%scalarLWNetGround)%dat(1),            & ! intent(out): net longwave radiation at the ground surface (W m-2)
                   mvar_data%var(iLookMVAR%scalarLWNetUbound)%dat(1),            & ! intent(out): net longwave radiation at the upper boundary (W m-2)
                   ! output: flux derivatives
                   dLWNetCanopy_dTCanopy,                                        & ! intent(out): derivative in net canopy radiation w.r.t. canopy temperature (W m-2 K-1)
                   dLWNetGround_dTGround,                                        & ! intent(out): derivative in net ground radiation w.r.t. ground temperature (W m-2 K-1)
                   dLWNetCanopy_dTGround,                                        & ! intent(out): derivative in net canopy radiation w.r.t. ground temperature (W m-2 K-1)
                   dLWNetGround_dTCanopy,                                        & ! intent(out): derivative in net ground radiation w.r.t. canopy temperature (W m-2 K-1)
                   ! output: error control
                   err,cmessage                                                  ) ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! compute aerodynamic resistances
  ! Refs: Choudhury and Monteith (4-layer model for heat budget of homogenous surfaces; QJRMS, 1988)
  !       Niu and Yang (Canopy effects on snow processes; JGR, 2004)
  !       Mahat et al. (Below-canopy turbulence in a snowmelt model, WRR, 2012)
  call aeroResist(&
                  ! input: model control
                  model_decisions(iLookDECISIONS%fDerivMeth)%iDecision,          & ! intent(in): method used to calculate flux derivatives
                  model_decisions(iLookDECISIONS%astability)%iDecision,          & ! intent(in): choice of stability function
                  ! input: above-canopy forcing data
                  mpar_data%var(iLookPARAM%mheight),                             & ! intent(in): measurement height (m)
                  forc_data%var(iLookFORCE%airtemp),                             & ! intent(in): air temperature at some height above the surface (K)
                  forc_data%var(iLookFORCE%windspd),                             & ! intent(in): wind speed at some height above the surface (m s-1)
                  ! input: canopy and ground temperature
                  mvar_data%var(iLookMVAR%scalarVegetationTemp)%dat(1),          & ! intent(in): canopy temperature (K)
                  mvar_data%var(iLookMVAR%mLayerTemp)%dat(1),                    & ! intent(in): ground temperature (K)
                  ! input: diagnostic variables
                  exposedVAI,                                                    & ! intent(in): exposed vegetation area index -- leaf plus stem (m2 m-2)
                  mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1),               & ! intent(in): snow depth (m)
                  ! input: parameters
                  z0Ground,                                                      & ! intent(in): roughness length of the ground (below canopy or non-vegetated surface [snow]) (m)
                  mpar_data%var(iLookPARAM%z0Canopy),                            & ! intent(in): roughness length of the canopy (m)
                  mpar_data%var(iLookPARAM%critRichNumber),                      & ! intent(in): critical value for the bulk Richardson number where turbulence ceases (-)
                  mpar_data%var(iLookPARAM%Louis79_bparam),                      & ! intent(in): parameter in Louis (1979) stability function
                  mpar_data%var(iLookPARAM%Louis79_cStar),                       & ! intent(in): parameter in Louis (1979) stability function
                  mpar_data%var(iLookPARAM%Mahrt87_eScale),                      & ! intent(in): exponential scaling factor in the Mahrt (1987) stability function
                  mpar_data%var(iLookPARAM%windReductionFactor),                 & ! intent(in): canopy wind reduction factor (-)                   
                  mpar_data%var(iLookPARAM%leafExchangeCoeff),                   & ! intent(in): turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
                  mpar_data%var(iLookPARAM%leafDimension),                       & ! intent(in): characteristic leaf dimension (m)
                  mpar_data%var(iLookPARAM%canopyHeight),                        & ! intent(in): canopy height (m) 
                  ! output: scalar resistances
                  mvar_data%var(iLookMVAR%scalarZeroPlaneDisplacement)%dat(1),   & ! intent(out): zero plane displacement (m) 
                  mvar_data%var(iLookMVAR%scalarSfc2AtmExchangeCoeff)%dat(1),    & ! intent(out): surface-atmosphere turbulent exchange coefficient (-)
                  mvar_data%var(iLookMVAR%scalarEddyDiffusCanopyTop)%dat(1),     & ! intent(out): eddy diffusivity for heat at the top of the canopy (m2 s-1)
                  mvar_data%var(iLookMVAR%scalarWindspdCanopyTop)%dat(1),        & ! intent(out): windspeed at the top of the canopy (m s-1)
                  mvar_data%var(iLookMVAR%scalarLeafResistance)%dat(1),          & ! intent(out): mean leaf boundary layer resistance per unit leaf area (s m-1)
                  mvar_data%var(iLookMVAR%scalarGroundResistance)%dat(1),        & ! intent(out): below canopy aerodynamic resistance (s m-1) 
                  mvar_data%var(iLookMVAR%scalarCanopyResistance)%dat(1),        & ! intent(out): above canopy aerodynamic resistance (s m-1)
                  ! output: error control
                  err,cmessage                                                   ) ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! compute soil moisture factor controlling stomatal resistance
  call soilResist(&
                  ! input (model decisions)
                  model_decisions(iLookDECISIONS%soilStress)%iDecision,          & ! intent(in): choice of function for the soil moisture control on stomatal resistance
                  model_decisions(iLookDECISIONS%groundwatr)%iDecision,          & ! intent(in): groundwater parameterization
                  ! input (state variables)
                  mvar_data%var(iLookMVAR%mLayerMatricHead)%dat,                 & ! intent(in): matric head in each layer (m)
                  mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(nSnow+1:nLayers),& ! intent(in): volumetric fraction of liquid water in each layer (-)
                  mvar_data%var(iLookMVAR%scalarAquiferStorage)%dat(1),          & ! intent(in): aquifer storage (m)
                  ! input (diagnostic variables)
                  mvar_data%var(iLookMVAR%mLayerRootDensity)%dat,                & ! intent(in): root density in each layer (-)
                  mvar_data%var(iLookMVAR%scalarAquiferRootFrac)%dat(1),         & ! intent(in): fraction of roots below the lowest soil layer (-)
                  ! input (parameters)
                  mpar_data%var(iLookPARAM%plantWiltPsi),                        & ! intent(in): matric head at wilting point (m)
                  mpar_data%var(iLookPARAM%soilStressParam),                     & ! intent(in): parameter in the exponential soil stress function (-)
                  mpar_data%var(iLookPARAM%critSoilWilting),                     & ! intent(in): critical vol. liq. water content when plants are wilting (-)
                  mpar_data%var(iLookPARAM%critSoilTranspire),                   & ! intent(in): critical vol. liq. water content when transpiration is limited (-)
                  mpar_data%var(iLookPARAM%critAquiferTranspire),                & ! intent(in): critical aquifer storage value when transpiration is limited (m)
                  ! output
                  mvar_data%var(iLookMVAR%scalarTranspireLim)%dat(1),            & ! intent(out): weighted average of the transpiration limiting factor (-)
                  mvar_data%var(iLookMVAR%mLayerTranspireLim)%dat,               & ! intent(out): transpiration limiting factor in each layer (-)
                  mvar_data%var(iLookMVAR%scalarTranspireLimAqfr)%dat(1),        & ! intent(out): transpiration limiting factor for the aquifer (-)
                  err,cmessage                                                   ) ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! compute the saturation vapor pressure for vegetation temperature
  TV_celcius = mvar_data%var(iLookMVAR%scalarVegetationTemp)%dat(1) - Tfreeze
  call satVapPress(TV_celcius, mvar_data%var(iLookMVAR%scalarSatVP_VegTemp)%dat(1), dSVP_dTV)

  ! compute stomatal resistance
  call stomResist(&
                  ! input (model decisions)
                  model_decisions(iLookDECISIONS%stomResist)%iDecision,          & ! intent(in): choice of function for stomatal resistance
                  ! input (local attributes)
                  type_data%var(iLookTYPE%vegTypeIndex),                         & ! intent(in): vegetation type index
                  iLoc, jLoc,                                                    & ! intent(in): spatial location indices      
                  ! input (forcing)
                  forc_data%var(iLookFORCE%airtemp),                             & ! intent(in): air temperature at some height above the surface (K)
                  forc_data%var(iLookFORCE%airpres),                             & ! intent(in): air pressure at some height above the surface (Pa)
                  mvar_data%var(iLookMVAR%scalarO2air)%dat(1),                   & ! intent(in): atmospheric o2 concentration (Pa)
                  mvar_data%var(iLookMVAR%scalarCO2air)%dat(1),                  & ! intent(in): atmospheric co2 concentration (Pa)
                  mvar_data%var(iLookMVAR%scalarCanopySunlitPAR)%dat(1),         & ! intent(in): average absorbed par for sunlit leaves (w m-2)
                  mvar_data%var(iLookMVAR%scalarCanopyShadedPAR)%dat(1),         & ! intent(in): average absorbed par for shaded leaves (w m-2)
                  ! input (state and diagnostic variables)
                  mvar_data%var(iLookMVAR%scalarGrowingSeasonIndex)%dat(1),      & ! intent(in): growing season index (0=off, 1=on)
                  mvar_data%var(iLookMVAR%scalarFoliageNitrogenFactor)%dat(1),   & ! intent(in): foliage nitrogen concentration (1=saturated)
                  mvar_data%var(iLookMVAR%scalarTranspireLim)%dat(1),            & ! intent(in): weighted average of the soil moiture factor controlling stomatal resistance (-)
                  mvar_data%var(iLookMVAR%scalarLeafResistance)%dat(1),          & ! intent(in): leaf boundary layer resistance (s m-1)
                  mvar_data%var(iLookMVAR%scalarVegetationTemp)%dat(1),          & ! intent(in): vegetation temperature (K)
                  mvar_data%var(iLookMVAR%scalarSatVP_VegTemp)%dat(1),           & ! intent(in): saturation vapor pressure at vegetation temperature (Pa)
                  mvar_data%var(iLookMVAR%scalarVP_CanopyAir)%dat(1),            & ! intent(in): canopy air vapor pressure (Pa)
                  ! output
                  mvar_data%var(iLookMVAR%scalarStomResistSunlit)%dat(1),        & ! intent(out): stomatal resistance for sunlit leaves (s m-1)
                  mvar_data%var(iLookMVAR%scalarStomResistShaded)%dat(1),        & ! intent(out): stomatal resistance for shaded leaves (s m-1)
                  mvar_data%var(iLookMVAR%scalarPhotosynthesisSunlit)%dat(1),    & ! intent(out): sunlit photosynthesis (umolco2 m-2 s-1)
                  mvar_data%var(iLookMVAR%scalarPhotosynthesisShaded)%dat(1),    & ! intent(out): shaded photosynthesis (umolco2 m-2 s-1)
                  err,message                                                    ) ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  pause 'after stomResist'




 end subroutine noahMP_veg


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
 USE multiconst,only:sb ! Stefan Boltzman constant (W m-2 K-4)
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
                       ixDerivMethod,                 & ! intent(in): choice of method used to compute derivative (analytical or numerical)
                       ixStability,                   & ! intent(in): choice of stability function
                       ! input: above-canopy forcing data
                       mheight,                       & ! intent(in): measurement height (m)
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
                       windReductionFactor,           & ! intent(in): canopy wind reduction factor (-)                   
                       leafExchangeCoeff,             & ! intent(in): turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
                       leafDimension,                 & ! intent(in): characteristic leaf dimension (m)
                       canopyHeight,                  & ! intent(in): canopy height (m) 
                       ! output: scalar resistances
                       zeroPlaneDisplacement,         & ! intent(out): zero plane displacement (m) 
                       sfc2AtmExchangeCoeff,          & ! intent(out): surface-atmosphere turbulent exchange coefficient (-)
                       eddyDiffusCanopyTop,           & ! intent(out): eddy diffusivity for heat at the top of the canopy (m2 s-1)
                       windspdCanopyTop,              & ! intent(out): windspeed at the top of the canopy (m s-1)
                       leafResistance,                & ! intent(out): mean leaf boundary layer resistance per unit leaf area (s m-1)
                       groundResistance,              & ! intent(out): below canopy aerodynamic resistance (s m-1) 
                       canopyResistance,              & ! intent(out): above canopy aerodynamic resistance (s m-1)
                       ! output: error control
                       err,message                    ) ! intent(out): error control
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! compute aerodynamic resistances
 ! Refs: Choudhury and Monteith (4-layer model for heat budget of homogenous surfaces; QJRMS, 1988)
 !       Niu and Yang (Canopy effects on snow processes; JGR, 2004)
 !       Mahat et al. (Below-canopy turbulence in a snowmelt model, WRR, 2012)
 implicit none
 ! input: model control
 integer(i4b),intent(in)       :: ixDerivMethod            ! choice of method used to compute derivative (analytical or numerical)
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
 real(dp),intent(in)           :: windReductionFactor      ! canopy wind reduction factor (-)                   
 real(dp),intent(in)           :: leafExchangeCoeff        ! turbulent exchange coeff between canopy surface and canopy air ( m s-(1/2) )
 real(dp),intent(in)           :: leafDimension            ! characteristic leaf dimension (m)
 real(dp),intent(in)           :: canopyHeight             ! canopy height (m) 
 ! output: scalar resistances
 real(dp),intent(out)          :: zeroPlaneDisplacement    ! zero plane displacement (m) 
 real(dp),intent(out)          :: sfc2AtmExchangeCoeff     ! surface-atmosphere turbulent exchange coefficient (-)
 real(dp),intent(out)          :: eddyDiffusCanopyTop      ! eddy diffusivity for heat at the top of the canopy (m2 s-1)
 real(dp),intent(out)          :: windspdCanopyTop         ! windspeed at the top of the canopy (m s-1)
 real(dp),intent(out)          :: leafResistance           ! mean leaf boundary layer resistance per unit leaf area (s m-1)
 real(dp),intent(out)          :: groundResistance         ! below canopy aerodynamic resistance (s m-1) 
 real(dp),intent(out)          :: canopyResistance         ! above canopy aerodynamic resistance (s m-1)
 ! output: error control
 integer(i4b),intent(out)      :: err                      ! error code
 character(*),intent(out)      :: message                  ! error message
 ! -----------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 character(LEN=256)            :: cmessage                 ! error message of downwind routine
 real(dp)                      :: funcLAI                  ! temporary variable to calculate zero plane displacement for the canopy
 real(dp)                      :: z0                       ! roughness length (m)
 real(dp)                      :: sfcTemp                  ! surface temperature (K)
 real(dp)                      :: tmp1,tmp2                ! temporary variables used in calculation of leaf resistance
 real(dp)                      :: dSfc2AtmExCoef_dTemp     ! derivative in surface-atmosphere exchange coeffcient w.r.t. temperature (K-1)
 real(dp)                      :: frictionVelocity         ! friction velocity (m s-1)
 real(dp)                      :: leafConductance          ! leaf boundary layer conductance (m s-1) 
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
  funcLAI   = sqrt(7.5_dp*exposedVAI)
  zeroPlaneDisplacement = -(1._dp - exp(-funcLAI))/funcLAI + 1._dp
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
                 (ixDerivMethod==analytical),    & ! input: logical flag to compute analytical derivatives
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
                 err, cmessage                   ) ! output: error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


 ! ***** compute the above-canopy, below-canopy, and leaf resistances (s m-1)
 if(exposedVAI > 0._dp) then ! (if vegetation is exposed)

  ! compute the friction velocity (m s-1)
  frictionVelocity = windspd * sqrt(sfc2AtmExchangeCoeff)

  ! compute the above-canopy resistance (s m-1)
  canopyResistance = 1._dp/(sfc2AtmExchangeCoeff*windspd)

  ! compute windspeed at the top of the canopy (m s-1) -- don't use stability corrections to simplify calculation of derivatives)
  windspdCanopyTop = (1._dp/vkc)*log(canopyHeight/z0Canopy)/log(mHeight/z0Canopy)

  ! compute the leaf boundary layer resistance (s m-1)
  leafConductance = (2._dp*leafExchangeCoeff/windReductionFactor) * sqrt(windspdCanopyTop/leafDimension) * (1._dp - exp(-windReductionFactor/2._dp))
  leafResistance  = 1._dp/(leafConductance)

  ! compute eddy diffusivity for heat at the top of the canopy (m2 s-1)
  !   Note: use of friction velocity here includes stability adjustments
  eddyDiffusCanopyTop = max(vkc*FrictionVelocity*(canopyHeight - zeroPlaneDisplacement), mpe)  ! (avoid divide by zero)

  ! compute the resistance between the surface and canopy air
  !  assume exponential profile extends from the surface roughness length to the displacement height plus vegetation roughness
  tmp1 = exp(-windReductionFactor* z0Ground/canopyHeight)
  tmp2 = exp(-windReductionFactor*(z0Canopy+zeroPlaneDisplacement)/canopyHeight)
  groundResistance = ( canopyHeight*exp(windReductionFactor) / (windReductionFactor*eddyDiffusCanopyTop) ) * (tmp1 - tmp2)  ! s m-1

 ! ***** compute resistances for non-vegetated surfaces (e.g., snow)
 else
  canopyResistance = 0._dp
  leafResistance   = 0._dp
  groundResistance = 1._dp / (sfc2AtmExchangeCoeff*windspd)
 endif

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

 end do  ! (looping through sunlit and shaded leaves

 end subroutine stomResist

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
 real(dp),intent(in)           :: mheight                ! measurement height (m)
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

 ! compute the bulk Richardson number (-)
 call bulkRichardson(airTemp,sfcTemp,windspd,mheight,computeDerivative, & ! (input)
                     RiBulk,dRiBulk_dTemp,err,message)                    ! (output)

 ! compute turbulent exchange coefficient under conditions of neutral stability
 ExNeut = (vkc**2._dp) / ( log((mHeight - zeroPlaneDisplacement)/z0))**2._dp

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
 subroutine bulkRichardson(airtemp,sfcTemp,windspd,mheight,computeDerivative, & ! (input)
                           RiBulk,dRiBulk_dTemp,err,message)                    ! (output)
 implicit none
 ! input
 real(dp),intent(in)           :: airtemp       ! air temperature (K)
 real(dp),intent(in)           :: sfcTemp       ! surface temperature (K)
 real(dp),intent(in)           :: windspd       ! wind speed (m s-1)
 real(dp),intent(in)           :: mheight       ! measurement height (m)
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
 RiMult = (gravity*mheight)/(windspd*windspd)
 ! compute the Richardson number
 RiBulk = (T_grad/T_mean) * RiMult
 ! compute the derivative in the Richardson number
 if(computeDerivative)then
  dRiBulk_dTemp = (-1._dp/T_mean + T_grad*(-0.5_dp)*T_mean**(-2._dp) ) * RiMult
 else
  dRiBulk_dTemp = 1._dp
 endif
 end subroutine bulkRichardson




end module noahMP_veg_module

module getDiagVar_module
! data types
USE nrtype
implicit none
private
public::getDiagVar
! algorithmic parameters
real(dp),parameter     :: valueMissing=-9999._dp  ! missing value, used when diagnostic or state variables are undefined
real(dp),parameter     :: verySmall=1.e-6_dp   ! used as an additive constant to check if substantial difference among real numbers
contains

 ! ************************************************************************************************
 ! new subroutine: get diagnostic variables treated as constant over a model substep
 ! ************************************************************************************************
 subroutine getDiagVar(&
                       dt,                          & ! intent(in): time step (s)
                       computeVegFlux,              & ! intent(out): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                       err,message)                   ! intent(out): error control
 ! data structures and named variables
 USE data_struc,only:model_decisions                                                ! model decision structure
 USE data_struc,only:type_data,attr_data,forc_data,mpar_data,mvar_data,indx_data    ! data structures
 USE var_lookup,only:iLookDECISIONS                                                 ! named variables for elements of the decision structure
 USE var_lookup,only:iLookTYPE,iLookATTR,iLookFORCE,iLookPARAM,iLookMVAR,iLookINDEX ! named variables for structure elements
 ! common variables
 USE data_struc,only:urbanVegCategory       ! vegetation category for urban areas
 USE data_struc,only:fracJulday             ! fractional julian days since the start of year
 USE data_struc,only:yearLength             ! number of days in the current year
 implicit none
 real(dp),intent(in)                  :: dt                          ! time step (seconds)
 
 logical(lgt),intent(out)             :: computeVegFlux              ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 integer(i4b),intent(out)             :: err                         ! error code
 character(*),intent(out)             :: message                     ! error message
 ! -------------------------------------------------------------------------------------------------
 ! local
 character(LEN=256)                   :: cmessage                    ! error message of downwind routine
 ! initialize error control
 err=0; message="getDiagVar/"

 ! compute variables that are treated as constant over a model substep...
 call getDiagVar_muster(&
                        ! input: model control
                        dt,                                                        & ! intent(in): time step (s)
                        yearLength,                                                & ! intent(in): number of days in the current year
                        fracJulday,                                                & ! intent(in): fractional julian days since the start of year
                        attr_data%var(iLookATTR%latitude),                         & ! intent(in): latitude
                        model_decisions(iLookDECISIONS%canopySrad)%iDecision,      & ! intent(in): index of method used for canopy sw radiation
                        model_decisions(iLookDECISIONS%alb_method)%iDecision,      & ! intent(in): index of method used for snow albedo
                        ! input: state variables
                        mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1),           & ! intent(in): canopy ice content (kg m-2)
                        mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1),           & ! intent(in): canopy liquid water content (kg m-2)
                        mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat,             & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
                        mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat,             & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)
                        ! input: coordinate variables
                        indx_data%var(iLookINDEX%layerType)%dat,                   & ! intent(in): layer type (ix_soil or ix_snow)
                        mvar_data%var(iLookMVAR%mLayerHeight)%dat,                 & ! intent(in): height at the mid-point of each layer (m)
                        mvar_data%var(iLookMVAR%iLayerHeight)%dat,                 & ! intent(in): height at the interface of each layer (m)
                        ! input: vegetation phenology
                        type_data%var(iLookTYPE%vegTypeIndex),                     & ! intent(in): vegetation type index
                        urbanVegCategory,                                          & ! intent(in): vegetation category for urban areas               
                        mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1),           & ! intent(in): snow depth on the ground surface (m)
                        mvar_data%var(iLookMVAR%scalarCanopyTemp)%dat(1),          & ! intent(in): temperature of the vegetation canopy at the start of the sub-step (K)
                        mvar_data%var(iLookMVAR%scalarRootZoneTemp)%dat(1),        & ! intent(in): root zone temperature (K)
                        mpar_data%var(iLookPARAM%heightCanopyTop),                 & ! intent(in): height of the top of the canopy layer (m)
                        mpar_data%var(iLookPARAM%heightCanopyBottom),              & ! intent(in): height of the bottom of the canopy layer (m)
                        mvar_data%var(iLookMVAR%scalarLAI)%dat(1),                 & ! intent(inout): one-sided leaf area index (m2 m-2)
                        mvar_data%var(iLookMVAR%scalarSAI)%dat(1),                 & ! intent(inout): one-sided stem area index (m2 m-2)
                        ! input: heat capacity and thermal conductivity
                        mpar_data%var(iLookPARAM%specificHeatVeg),                 & ! intent(in): specific heat of vegetation (J kg-1 K-1)
                        mpar_data%var(iLookPARAM%maxMassVegetation),               & ! intent(in): maximum mass of vegetation (kg m-2)
                        mpar_data%var(iLookPARAM%soil_dens_intr),                  & ! intent(in): intrinsic density of soil (kg m-3)
                        mpar_data%var(iLookPARAM%thCond_soil),                     & ! intent(in): thermal conductivity of soil (W m-1 K-1)
                        mpar_data%var(iLookPARAM%theta_sat),                       & ! intent(in): soil porosity (-)
                        mpar_data%var(iLookPARAM%frac_sand),                       & ! intent(in): fraction of sand (-)
                        mpar_data%var(iLookPARAM%frac_silt),                       & ! intent(in): fraction of silt (-)
                        mpar_data%var(iLookPARAM%frac_clay),                       & ! intent(in): fraction of clay (-)
                        ! input: snow albedo
                        mpar_data%var(iLookPARAM%Frad_vis),                        & ! intent(in): fraction of radiation in visible part of spectrum (-)
                        mpar_data%var(iLookPARAM%Frad_direct),                     & ! intent(in): fraction direct solar radiation (-)
                        mpar_data%var(iLookPARAM%albedoMax),                       & ! intent(in): maximum snow albedo for a single spectral band (-)
                        mpar_data%var(iLookPARAM%albedoMinWinter),                 & ! intent(in): minimum snow albedo during winter for a single spectral band (-)
                        mpar_data%var(iLookPARAM%albedoMinSpring),                 & ! intent(in): minimum snow albedo during spring for a single spectral band (-)
                        mpar_data%var(iLookPARAM%albedoMaxVisible),                & ! intent(in): maximum snow albedo in the visible part of the spectrum (-)
                        mpar_data%var(iLookPARAM%albedoMinVisible),                & ! intent(in): minimum snow albedo in the visible part of the spectrum (-)
                        mpar_data%var(iLookPARAM%albedoMaxNearIR),                 & ! intent(in): maximum snow albedo in the near infra-red part of the spectrum (-)
                        mpar_data%var(iLookPARAM%albedoMinNearIR),                 & ! intent(in): minimum snow albedo in the near infra-red part of the spectrum (-)
                        mpar_data%var(iLookPARAM%albedoDecayRate),                 & ! intent(in): albedo decay rate (s)
                        mpar_data%var(iLookPARAM%tempScalGrowth),                  & ! intent(in): temperature scaling factor for grain growth (K-1) 
                        mpar_data%var(iLookPARAM%albedoSootLoad),                  & ! intent(in): soot load factor (-)
                        mpar_data%var(iLookPARAM%albedoRefresh),                   & ! intent(in): critical mass necessary for albedo refreshment (kg m-2)
                        mpar_data%var(iLookPARAM%snowfrz_scale),                   & ! intent(in): scaling parameter for the freezing curve for snow (K-1) 
                        mvar_data%var(iLookMVAR%mLayerTemp)%dat(1),                & ! intent(in): surface temperature
                        mvar_data%var(iLookMVAR%scalarSnowfall)%dat(1),            & ! intent(in): snowfall rate (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarCosZenith)%dat(1),           & ! intent(in): cosine of the zenith angle (-)
                        ! ---------------------------------------------------------------------------------------------------------------------------------
                        ! output: vegetation phenology
                        mvar_data%var(iLookMVAR%scalarExposedLAI)%dat(1),          & ! intent(out): exposed leaf area index after burial by snow (m2 m-2)
                        mvar_data%var(iLookMVAR%scalarExposedSAI)%dat(1),          & ! intent(out): exposed stem area index after burial by snow (m2 m-2)
                        mvar_data%var(iLookMVAR%scalarGrowingSeasonIndex)%dat(1),  & ! intent(out): growing season index (0=off, 1=on)
                        computeVegFlux,                                            & ! intent(out): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                        ! output: heat capacity and thermal conductivity
                        mvar_data%var(iLookMVAR%scalarBulkVolHeatCapVeg)%dat(1),   & ! intent(out): volumetric heat capacity of the vegetation (J m-3 K-1)
                        mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat,           & ! intent(out): volumetric heat capacity in each layer (J m-3 K-1) 
                        mvar_data%var(iLookMVAR%mLayerThermalC)%dat,               & ! intent(out): thermal conductivity at the mid-point of each layer (W m-1 K-1)
                        mvar_data%var(iLookMVAR%iLayerThermalC)%dat,               & ! intent(out): thermal conductivity at the interface of each layer (W m-1 K-1)
                        mvar_data%var(iLookMVAR%mLayerVolFracAir)%dat,             & ! intent(out): volumetric fraction of air in each layer (-)
                        ! output: snow albedo
                        mvar_data%var(iLookMVAR%spectralSnowAlbedoDiffuse)%dat,    & ! intent(inout): diffuse snow albedo in each spectral band (-)
                        mvar_data%var(iLookMVAR%spectralSnowAlbedoDirect)%dat,     & ! intent(inout): direct snow albedo in each spectral band (-)
                        mvar_data%var(iLookMVAR%scalarSnowAlbedo)%dat(1),          & ! intent(inout): snow albedo for the entire spectral band (-)
                        ! output: error control
                        err,cmessage)                                                ! intent(out): error control
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

 end subroutine getDiagVar



 ! *****************************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************
 ! ***** PRIVATE SUBROUTINES ***************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************

 ! ************************************************************************************************
 ! private subroutine: get diagnostic variables treated as constant over a model substep
 ! ************************************************************************************************
 subroutine getDiagVar_muster(&
                              ! input: model control
                              dt,                          & ! intent(in): time step (s)
                              yearLength,                  & ! intent(in): number of days in the current year
                              fracJulday,                  & ! intent(in): fractional julian days since the start of year
                              latitude,                    & ! intent(in): latitude
                              ixCanSWMethod,               & ! intent(in): index of the method used for canopy shortwave radiation
                              ixAlbedoMethod,              & ! intent(in): index of the method used for snow albedo
                              ! input: state variables
                              scalarCanopyIce,             & ! intent(in): canopy ice content (kg m-2)
                              scalarCanopyLiq,             & ! intent(in): canopy liquid water content (kg m-2)
                              mLayerVolFracIce,            & ! intent(in): volumetric fraction of ice in each layer (-)
                              mLayerVolFracLiq,            & ! intent(in): volumetric fraction of liquid water in each layer (-)
                              ! input: coordinate variables
                              layerType,                   & ! intent(in): layer type (ix_soil or ix_snow)
                              mLayerHeight,                & ! intent(in): height at the mid-point of each layer (m)
                              iLayerHeight,                & ! intent(in): height at the interface of each layer (m)
                              ! input: vegetation phenology
                              vegTypeIndex,                & ! intent(in): vegetation type index
                              urbanVegCategory,            & ! intent(in): vegetation category for urban areas               
                              scalarSnowDepth,             & ! intent(in): snow depth on the ground surface (m)
                              scalarCanopyTemp,            & ! intent(in): temperature of the vegetation canopy at the start of the sub-step (K)
                              scalarRootZoneTemp,          & ! intent(in): root zone temperature (K)
                              heightCanopyTop,             & ! intent(in): height of the top of the canopy (m)
                              heightCanopyBottom,          & ! intent(in): height of the bottom of the canopy (m)
                              scalarLAI,                   & ! intent(inout): one-sided leaf area index (m2 m-2)
                              scalarSAI,                   & ! intent(inout): one-sided stem area index (m2 m-2)
                              ! input: heat capacity and thermal conductivity
                              specificHeatVeg,             & ! intent(in): specific heat of vegetation (J kg-1 K-1)
                              maxMassVegetation,           & ! intent(in): maximum mass of vegetation (kg m-2)
                              soil_dens_intr,              & ! intent(in): intrinsic density of soil (kg m-3)
                              thCond_soil,                 & ! intent(in): thermal conductivity of soil (W m-1 K-1)
                              theta_sat,                   & ! intent(in): soil porosity (-)
                              frac_sand,                   & ! intent(in): fraction of sand (-)
                              frac_silt,                   & ! intent(in): fraction of silt (-)
                              frac_clay,                   & ! intent(in): fraction of clay (-)
                              ! input: snow albedo
                              Frad_vis,                    & ! intent(in): fraction of radiation in visible part of spectrum (-)
                              Frad_direct,                 & ! intent(in): fraction direct solar radiation (-)
                              albedoMax,                   & ! intent(in): maximum snow albedo for a single spectral band (-)
                              albedoMinWinter,             & ! intent(in): minimum snow albedo during winter for a single spectral band (-)
                              albedoMinSpring,             & ! intent(in): minimum snow albedo during spring for a single spectral band (-)
                              albedoMaxVisible,            & ! intent(in): maximum snow albedo in the visible part of the spectrum (-)
                              albedoMinVisible,            & ! intent(in): minimum snow albedo in the visible part of the spectrum (-)
                              albedoMaxNearIR,             & ! intent(in): maximum snow albedo in the near infra-red part of the spectrum (-)
                              albedoMinNearIR,             & ! intent(in): minimum snow albedo in the near infra-red part of the spectrum (-)
                              albedoDecayRate,             & ! intent(in): albedo decay rate (s)
                              tempScalGrowth,              & ! intent(in): temperature scaling factor for grain growth (K-1) 
                              albedoSootLoad,              & ! intent(in): soot load factor (-)
                              albedoRefresh,               & ! intent(in): critical mass necessary for albedo refreshment (kg m-2)
                              snowfrz_scale,               & ! intent(in): scaling parameter for the freezing curve for snow (K-1) 
                              scalarSurfTemp,              & ! intent(in): surface temperature
                              scalarSnowfall,              & ! intent(in): snowfall rate (kg m-2 s-1)
                              scalarCosZenith,             & ! intent(in): cosine of the zenith angle (-)
                              ! --------------------------------------------------------------------------------
                              ! output: vegetation phenology
                              scalarExposedLAI,            & ! intent(out): exposed leaf area index after burial by snow (m2 m-2)
                              scalarExposedSAI,            & ! intent(out): exposed stem area index after burial by snow (m2 m-2)
                              scalarGrowingSeasonIndex,    & ! intent(out): growing season index (0=off, 1=on)
                              computeVegFlux,              & ! intent(out): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                              ! output: heat capacity and thermal conductivity
                              scalarBulkVolHeatCapVeg,     & ! intent(out): volumetric heat capacity of the vegetation (J m-3 K-1)
                              mLayerVolHtCapBulk,          & ! intent(out): volumetric heat capacity in each layer (J m-3 K-1) 
                              mLayerThermalC,              & ! intent(out): thermal conductivity at the mid-point of each layer (W m-1 K-1)
                              iLayerThermalC,              & ! intent(out): thermal conductivity at the interface of each layer (W m-1 K-1)
                              mLayerVolFracAir,            & ! intent(out): volumetric fraction of air in each layer (-)
                              ! output: snow albedo
                              spectralSnowAlbedoDiffuse,   & ! intent(inout): diffuse snow albedo in each spectral band (-)
                              spectralSnowAlbedoDirect,    & ! intent(inout): direct snow albedo in each spectral band (-)
                              scalarSnowAlbedo,            & ! intent(inout): snow albedo for the entire spectral band (-)
                              ! output: error control
                              err,message)                   ! intent(out): error control
 ! ----------------------------------------------------------------------------------------------------------------------------------
 USE NOAHMP_ROUTINES,only:phenology         ! determine vegetation phenology
 USE diagn_evar_module,only:diagn_evar      ! compute diagnostic energy variables -- thermal conductivity and heat capacity
 USE snowAlbedo_module,only:snowAlbedo      ! compute snow albedo
 ! look-up values for the choice of canopy shortwave radiation method
 USE mDecisions_module,only:      &
  noah_mp,                        &         ! full Noah-MP implementation (including albedo)
  CLM_2stream,                    &         ! CLM 2-stream model (see CLM documentation)
  UEB_2stream,                    &         ! UEB 2-stream model (Mahat and Tarboton, WRR 2011)
  NL_scatter,                     &         ! Simplified method Nijssen and Lettenmaier (JGR 1999)
  BeersLaw                                  ! Beer's Law (as implemented in VIC)
 ! named variables for snow and soil
 USE data_struc,only:ix_soil,ix_snow
 implicit none
 ! input
 real(dp),intent(in)                  :: dt                           ! time step (seconds)
 integer(i4b),intent(in)              :: yearLength                   ! number of days in the current year
 real(dp),intent(in)                  :: fracJulday                   ! fractional julian days since the start of year
 real(dp),intent(in)                  :: latitude                     ! latitude
 integer(i4b),intent(in)              :: ixCanSWMethod                ! index of method used for canopy sw radiation
 integer(i4b),intent(in)              :: ixAlbedoMethod               ! index of method used for snow albedo
 ! input: state variables
 real(dp),intent(in)                  :: scalarCanopyLiq              ! canopy liquid water content (kg m-2)
 real(dp),intent(in)                  :: scalarCanopyIce              ! canopy ice content (kg m-2)
 real(dp),intent(in)                  :: mLayerVolFracLiq(:)          ! volumetric fraction of liquid water in each layer (-)
 real(dp),intent(in)                  :: mLayerVolFracIce(:)          ! volumetric fraction of ice in each layer (-)
 ! input: coordinate variables 
 integer(i4b),intent(in)              :: layerType(:)                 ! type of the layer (ix_soil or ix_snow)
 real(dp),intent(in)                  :: mLayerHeight(:)              ! height at the mid-point of each layer (m)
 real(dp),intent(in)                  :: iLayerHeight(0:)             ! height at the interface of each layer (m)
 ! input: vegetation phenology
 integer(i4b),intent(in)              :: vegTypeIndex                 ! vegetation type index
 integer(i4b),intent(in)              :: urbanVegCategory             ! vegetation category for urban areas               
 real(dp),intent(in)                  :: scalarSnowDepth              ! snow depth on the ground surface (m)
 real(dp),intent(in)                  :: scalarCanopyTemp             ! temperature of the vegetation canopy at the start of the sub-step (K)
 real(dp),intent(in)                  :: scalarRootZoneTemp           ! root zone temperature (K)
 real(dp),intent(in)                  :: heightCanopyTop              ! height of the top of the canopy (m)
 real(dp),intent(in)                  :: heightCanopyBottom           ! height of the bottom of the canopy (m)
 real(dp),intent(inout)               :: scalarLAI                    ! one-sided leaf area index (m2 m-2)
 real(dp),intent(inout)               :: scalarSAI                    ! one-sided stem area index (m2 m-2)
 ! input: heat capacity and thermal conductivity
 real(dp),intent(in)                  :: specificHeatVeg              ! specific heat of vegetation (J kg-1 K-1)
 real(dp),intent(in)                  :: maxMassVegetation            ! maximum mass of vegetation (kg m-2)
 real(dp),intent(in)                  :: soil_dens_intr               ! intrinsic density of soil (kg m-3)
 real(dp),intent(in)                  :: thCond_soil                  ! thermal conductivity of soil (W m-1 K-1)
 real(dp),intent(in)                  :: theta_sat                    ! soil porosity (-)
 real(dp),intent(in)                  :: frac_sand                    ! fraction of sand (-)
 real(dp),intent(in)                  :: frac_silt                    ! fraction of silt (-)
 real(dp),intent(in)                  :: frac_clay                    ! fraction of clay (-)
 ! input: snow albedo     
 real(dp),intent(in)                  :: Frad_vis                     ! fraction of radiation in visible part of spectrum (-)
 real(dp),intent(in)                  :: Frad_direct                  ! fraction direct solar radiation (-)
 real(dp),intent(in)                  :: albedoMax                    ! maximum snow albedo for a single spectral band (-)
 real(dp),intent(in)                  :: albedoMinWinter              ! minimum snow albedo during winter for a single spectral band (-)
 real(dp),intent(in)                  :: albedoMinSpring              ! minimum snow albedo during spring for a single spectral band (-)
 real(dp),intent(in)                  :: albedoMaxVisible             ! maximum snow albedo in the visible part of the spectrum (-)
 real(dp),intent(in)                  :: albedoMinVisible             ! minimum snow albedo in the visible part of the spectrum (-)
 real(dp),intent(in)                  :: albedoMaxNearIR              ! maximum snow albedo in the near infra-red part of the spectrum (-)
 real(dp),intent(in)                  :: albedoMinNearIR              ! minimum snow albedo in the near infra-red part of the spectrum (-)
 real(dp),intent(in)                  :: albedoDecayRate              ! albedo decay rate (s)
 real(dp),intent(in)                  :: tempScalGrowth               ! temperature scaling factor for grain growth (K-1) 
 real(dp),intent(in)                  :: albedoSootLoad               ! soot load factor (-)
 real(dp),intent(in)                  :: albedoRefresh                ! critical mass necessary for albedo refreshment (kg m-2)
 real(dp),intent(in)                  :: snowfrz_scale                ! scaling parameter for the freezing curve for snow (K-1)
 real(dp),intent(in)                  :: scalarSurfTemp               ! surface temperature
 real(dp),intent(in)                  :: scalarSnowfall               ! snowfall rate (kg m-2 s-1)
 real(dp),intent(in)                  :: scalarCosZenith              ! cosine of the zenith angle (-)
 ! output: vegetation phenology
 real(dp),intent(out)                 :: scalarExposedLAI             ! exposed leaf area index after burial by snow (m2 m-2)
 real(dp),intent(out)                 :: scalarExposedSAI             ! exposed stem area index after burial by snow (m2 m-2)
 real(dp),intent(out)                 :: scalarGrowingSeasonIndex     ! growing season index (0=off, 1=on)
 logical(lgt),intent(out)             :: computeVegFlux               ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! output: heat capacity and thermal conductivity
 real(dp),intent(out)                 :: scalarBulkVolHeatCapVeg      ! volumetric heat capacity of the vegetation (J m-3 K-1)
 real(dp),intent(out)                 :: mLayerVolHtCapBulk(:)        ! volumetric heat capacity in each layer (J m-3 K-1) 
 real(dp),intent(out)                 :: mLayerThermalC(:)            ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
 real(dp),intent(out)                 :: iLayerThermalC(0:)           ! thermal conductivity at the interface of each layer (W m-1 K-1)
 real(dp),intent(out)                 :: mLayerVolFracAir(:)          ! volumetric fraction of air in each layer (-)
 ! output: snow albedo
 real(dp),intent(inout)               :: spectralSnowAlbedoDiffuse(:) ! diffuse snow albedo in each spectral band (-)
 real(dp),intent(inout)               :: spectralSnowAlbedoDirect(:)  ! direct snow albedo in each spectral band (-)
 real(dp),intent(inout)               :: scalarSnowAlbedo             ! snow albedo for the entire spectral band (-)
 ! output: error control
 integer(i4b),intent(out)             :: err                          ! error code
 character(*),intent(out)             :: message                      ! error message
 ! -------------------------------------------------------------------------------------------------
 ! local
 character(LEN=256)                   :: cmessage                     ! error message of downwind routine
 real(dp)                             :: notUsed_heightCanopyTop      ! for some reason the Noah-MP phenology routines output canopy height
 real(dp)                             :: exposedVAI                   ! exposed vegetation area index (LAI + SAI)
 real(dp)                             :: canopyDepth                  ! canopy depth (m)
 real(dp)                             :: heightAboveSnow              ! height top of canopy is above the snow surface (m)
 integer(i4b)                         :: nSnow                        ! number of snow layers
 ! initialize error control
 err=0; message="getDiagVar/"

 ! get the number of snow layers
 nSnow = count(layerType==ix_snow)

 ! * vegetation phenology...
 ! -------------------------

 ! determine vegetation phenology
 ! NOTE: recomputing phenology every sub-step accounts for changes in exposed vegetation associated with changes in snow depth
 call phenology(&
                ! input
                vegTypeIndex,                & ! intent(in): vegetation type index
                urbanVegCategory,            & ! intent(in): vegetation category for urban areas               
                scalarSnowDepth,             & ! intent(in): snow depth on the ground surface (m)
                scalarCanopyTemp,            & ! intent(in): temperature of the vegetation canopy at the start of the sub-step (K)
                latitude,                    & ! intent(in): latitude
                yearLength,                  & ! intent(in): number of days in the current year
                fracJulday,                  & ! intent(in): fractional julian days since the start of year
                scalarLAI,                   & ! intent(inout): one-sided leaf area index (m2 m-2)
                scalarSAI,                   & ! intent(inout): one-sided stem area index (m2 m-2)
                scalarRootZoneTemp,          & ! intent(in): root zone temperature (K)
                ! output
                notUsed_heightCanopyTop,     & ! intent(out): height of the top of the canopy layer (m)
                scalarExposedLAI,            & ! intent(out): exposed leaf area index after burial by snow (m2 m-2)
                scalarExposedSAI,            & ! intent(out): exposed stem area index after burial by snow (m2 m-2)
                scalarGrowingSeasonIndex     ) ! intent(out): growing season index (0=off, 1=on)

 ! determine additional phenological variables
 exposedVAI      = scalarExposedLAI + scalarExposedSAI   ! exposed vegetation area index (m2 m-2)
 canopyDepth     = heightCanopyTop - heightCanopyBottom  ! canopy depth (m)
 heightAboveSnow = heightCanopyTop - scalarSnowDepth     ! height top of canopy is above the snow surface (m)

 ! determine if need to include vegetation in the energy flux routines
 computeVegFlux  = (exposedVAI > 0.05_dp .and. heightAboveSnow > 0.05_dp)


 ! * heat capacity and thermal conductivity...
 ! -------------------------------------------

 ! compute diagnostic energy variables (thermal conductivity and volumetric heat capacity)
 call diagn_evar(&
                 ! input: control variables
                 computeVegFlux,            & ! intent(in): flag to denote if computing the vegetation flux
                 ! input: state variables
                 scalarCanopyIce,           & ! intent(in): mass of ice on the vegetation canopy (kg m-2)
                 scalarCanopyLiq,           & ! intent(in): mass of liquid water on the vegetation canopy (kg m-2)
                 mLayerVolFracIce,          & ! intent(in): volumetric fraction of ice in each layer (-)
                 mLayerVolFracLiq,          & ! intent(in): volumetric fraction of liquid water in each layer (-)
                 ! input: coordinate variables
                 layerType,                 & ! intent(in): layer type (ix_soil or ix_snow)
                 mLayerHeight,              & ! intent(in): height at the mid-point of each layer (m)
                 iLayerHeight,              & ! intent(in): height at the interface of each layer (m)
                 ! input: model parameters
                 canopyDepth,               & ! intent(in): depth of the vegetation canopy (m)
                 specificHeatVeg,           & ! intent(in): specific heat of vegetation (J kg-1 K-1)
                 maxMassVegetation,         & ! intent(in): maximum mass of vegetation (full foliage) (kg m-2)
                 soil_dens_intr,            & ! intent(in): intrinsic density of soil (kg m-3)
                 thCond_soil,               & ! intent(in): thermal conductivity of soil (W m-1 K-1)
                 theta_sat,                 & ! intent(in): soil porosity (-)
                 frac_sand,                 & ! intent(in): fraction of sand (-)
                 frac_silt,                 & ! intent(in): fraction of silt (-)
                 frac_clay,                 & ! intent(in): fraction of clay (-)
                 ! output: diagnostic variables
                 scalarBulkVolHeatCapVeg,   & ! intent(out): bulk volumetric heat capacity of vegetation (J m-3 K-1)
                 mLayerVolHtCapBulk,        & ! intent(out): volumetric heat capacity in each layer (J m-3 K-1) 
                 mLayerThermalC,            & ! intent(out): thermal conductivity at the mid-point of each layer (W m-1 K-1)
                 iLayerThermalC,            & ! intent(out): thermal conductivity at the interface of each layer (W m-1 K-1)
                 mLayerVolFracAir,          & ! intent(out): volumetric fraction of air in each layer (-)
                 ! output: error control
                 err,cmessage)                ! intent(out): error control
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! * compute the snow albedo...
 ! ----------------------------

 ! NOTE: if canopy radiation is noah-mp, then albedo is computed within the Noah-MP radiation routine
 if(ixCanSWMethod /= noah_mp)then
  call snowAlbedo(&
                  ! input: model control
                  dt,                         & ! intent(in): model time step
                  (nSnow > 0),                & ! intent(in): logical flag to denote if snow is present
                  ixAlbedoMethod,             & ! intent(in): index of method used for snow albedo
                  ! input: model variables
                  scalarSnowfall,             & ! intent(in): snowfall rate (kg m-2 s-1)
                  scalarSurfTemp,             & ! intent(in): surface temperature
                  scalarCosZenith,            & ! intent(in): cosine of the zenith angle (-)
                  ! input: model parameters
                  Frad_vis,                   & ! intent(in): fraction of radiation in visible part of spectrum (-)
                  Frad_direct,                & ! intent(in): fraction direct solar radiation (-)
                  albedoMax,                  & ! intent(in): maximum snow albedo for a single spectral band (-)
                  albedoMinWinter,            & ! intent(in): minimum snow albedo during winter for a single spectral band (-)
                  albedoMinSpring,            & ! intent(in): minimum snow albedo during spring for a single spectral band (-)
                  albedoMaxVisible,           & ! intent(in): maximum snow albedo in the visible part of the spectrum (-)
                  albedoMinVisible,           & ! intent(in): minimum snow albedo in the visible part of the spectrum (-)
                  albedoMaxNearIR,            & ! intent(in): maximum snow albedo in the near infra-red part of the spectrum (-)
                  albedoMinNearIR,            & ! intent(in): minimum snow albedo in the near infra-red part of the spectrum (-)
                  albedoDecayRate,            & ! intent(in): albedo decay rate (s)
                  tempScalGrowth,             & ! intent(in): temperature scaling factor for grain growth (K-1) 
                  albedoSootLoad,             & ! intent(in): soot load factor (-)
                  albedoRefresh,              & ! intent(in): critical mass necessary for albedo refreshment (kg m-2)
                  snowfrz_scale,              & ! intent(in): scaling parameter for the freezing curve for snow (K-1)
                  ! output
                  spectralSnowAlbedoDiffuse,  & ! intent(inout): direct snow albedo in each spectral band (-)
                  spectralSnowAlbedoDirect,   & ! intent(inout): direct snow albedo in each spectral band (-)
                  scalarSnowAlbedo,           & ! intent(inout): snow albedo for the entire spectral band (-)
                  err,cmessage)                 ! intent(out): error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
 endif  ! (if NOT using the Noah-MP radiation routine)

 end subroutine getDiagVar_muster

end module getDiagVar_module

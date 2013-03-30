module picardSolv_module
! data types
USE nrtype
! layer types
USE data_struc,only:ix_soil,ix_snow ! named variables for snow and soil
! constants
USE multiconst,only:&
                    gravity,      & ! acceleration of gravity              (m s-2)
                    Tfreeze,      & ! temperature at freezing              (K)
                    LH_fus,       & ! latent heat of fusion                (J kg-1)
                    iden_ice,     & ! intrinsic density of ice             (kg m-3)
                    iden_water      ! intrinsic density of liquid water    (kg m-3)
implicit none
private
public::picardSolv
! number of soil and snow layers
integer(i4b)        :: nSoil        ! number of soil layers
integer(i4b)        :: nSnow        ! number of snow layers
integer(i4b)        :: nLayers      ! total number of layers
integer(i4b)        :: nLevels      ! number of soil layers to use in the soil hydrology routine
contains

 ! ************************************************************************************************
 ! new subroutine: run the coupled energy-mass model for one timestep
 ! ************************************************************************************************
 subroutine picardSolv(&
                       ! input
                       dt,            & ! time step (s)
                       maxiter,       & ! maximum number of iterations
                       firstSubstep,  & ! flag to denote first sub-step
                       ! output
                       niter,         & ! number of iterations taken
                       err,message)     ! error control
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:attr_data,type_data,mpar_data,forc_data,mvar_data,indx_data    ! data structures
 USE var_lookup,only:iLookATTR,iLookTYPE,iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX ! named variables for structure elements
 ! model decision structures
 USE data_struc,only:model_decisions    ! model decision structure
 USE var_lookup,only:iLookDECISIONS     ! named variables for elements of the decision structure
 ! external subroutines
 USE var_derive_module,only:calcHeight  ! module to calculate height at layer interfaces and layer mid-point
 USE snwDensify_module,only:snwDensify  ! module to compute densification of snow
 ! common variables
 USE data_struc,only:urbanVegCategory   ! vegetation category for urban areas
 USE data_struc,only:fracJulday         ! fractional julian days since the start of year
 USE data_struc,only:yearLength         ! number of days in the current year
 implicit none
 ! input
 real(dp),intent(in)                  :: dt                       ! time step (seconds)
 integer(i4b),intent(in)              :: maxiter                  ! maximum number of iterations
 logical(lgt),intent(in)              :: firstSubStep             ! flag to indicate if we are processing the first sub-step
 ! output
 integer(i4b),intent(out)             :: niter                    ! number of iterations
 integer(i4b),intent(out)             :: err                      ! error code
 character(*),intent(out)             :: message                  ! error message
 ! local
 character(LEN=256)                   :: cmessage                 ! error message of downwind routine
 logical(lgt)                         :: computeVegFlux           ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 real(dp)                             :: exposedVAI               ! exposed vegetation area index (m2 m-2)
 real(dp)                             :: scalarCanopyTempNew      ! temperature of the vegetation canopy at the end of the sub-step (K)
 real(dp)                             :: scalarCanopyIceNew       ! mass of ice on the vegetation canopy at the end of the sub-step (kg m-2)
 real(dp)                             :: scalarCanopyLiqNew       ! mass of liquid water on the vegetation canopy at the end of the sub-step (kg m-2)
 real(dp),allocatable                 :: mLayerTempNew(:)         ! temperature of each snow/soil layer at the end of the sub-step (K)
 real(dp),allocatable                 :: mLayerVolFracIceNew(:)   ! volumetric fraction of ice at the end of the sub-step (-)
 real(dp),allocatable                 :: mLayerVolFracLiqNew(:)   ! volumetric fraction of liquid water at the end of the sub-step (-)
 real(dp),allocatable                 :: mLayerMatricHeadNew(:)   ! matric head at the end of the sub-step (m)
 real(dp)                             :: scalarAquiferStorageNew  ! aquifer storage at the end of the sub-step (m)

 ! initialize error control
 err=0; message="picardSolv/"

 ! define the total number of snow+soil layers
 nLayers = indx_data%var(iLookINDEX%nLayers)%dat(1) 

 ! identify the number of snow and soil layers
 nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat == ix_snow)
 nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat == ix_soil)

 ! identify the number of soil layers to use in the soil hydrology routine
 nLevels = nSoil  ! NOTE: always pass the full number of soil layers

 ! *****
 ! preliminaries for the picard solver...
 ! **************************************
 call picardSolv_prelim(&

                        ! input: site attributes
                        attr_data%var(iLookATTR%latitude),                           & ! intent(in): latitude
                        type_data%var(iLookTYPE%vegTypeIndex),                       & ! intent(in): vegetation type index
                        urbanVegCategory,                                            & ! intent(in): vegetation category for urban areas               

                        ! input: coordinate variables
                        indx_data%var(iLookINDEX%layerType)%dat,                     & ! intent(in): layer type (ix_soil or ix_snow)
                        mvar_data%var(iLookMVAR%mLayerHeight)%dat,                   & ! intent(in): height at the mid-point of each layer (m)
                        mvar_data%var(iLookMVAR%iLayerHeight)%dat,                   & ! intent(in): height at the interface of each layer (m)

                        ! input: model state variables
                        mvar_data%var(iLookMVAR%scalarCanopyTemp)%dat(1),            & ! intent(in): temperature of the vegetation canopy at the start of the sub-step (K)
                        mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1),             & ! intent(in): mass of ice on the vegetation canopy at the start of the sub-step (kg m-2)
                        mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1),             & ! intent(in): mass of liquid water on the vegetation canopy at the start of the sub-step (kg m-2)
                        mvar_data%var(iLookMVAR%mLayerTemp)%dat,                     & ! intent(in): temperature of each snow/soil layer at the start of the sub-step (K)
                        mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat,               & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
                        mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat,               & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)

                        ! input: diagnostic variables
                        yearLength,                                                  & ! intent(in): number of days in the current year
                        fracJulday,                                                  & ! intent(in): fractional julian days since the start of year
                        mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1),             & ! intent(in): snow depth on the ground surface (m)
                        mvar_data%var(iLookMVAR%scalarLAI)%dat(1),                   & ! intent(inout): one-sided leaf area index (m2 m-2)
                        mvar_data%var(iLookMVAR%scalarSAI)%dat(1),                   & ! intent(inout): one-sided stem area index (m2 m-2)

                        ! input: vegetation parameters
                        mpar_data%var(iLookPARAM%heightCanopyTop),                   & ! intent(in): height of top of the vegetation canopy above ground surface (m)
                        mpar_data%var(iLookPARAM%heightCanopyBottom),                & ! intent(in): height of bottom of the vegetation canopy above ground surface (m)
                        mpar_data%var(iLookPARAM%specificHeatVeg),                   & ! intent(in): specific heat of vegetation (J kg-1 K-1)
                        mpar_data%var(iLookPARAM%maxMassVegetation),                 & ! intent(in): maximum mass of vegetation (full foliage) (kg m-2)
                        mpar_data%var(iLookPARAM%rootingDepth),                      & ! intent(in): rooting depth of the vegetation (m)

                        ! input: soil parameters
                        mpar_data%var(iLookPARAM%soil_dens_intr),                    & ! intent(in): intrinsic density of soil (kg m-3)
                        mpar_data%var(iLookPARAM%theta_sat),                         & ! intent(in): soil porosity (-)
                        mpar_data%var(iLookPARAM%frac_sand),                         & ! intent(in): fraction of sand (-)
                        mpar_data%var(iLookPARAM%frac_silt),                         & ! intent(in): fraction of silt (-)
                        mpar_data%var(iLookPARAM%frac_clay),                         & ! intent(in): fraction of clay (-)

                        ! output: vegetation phenology
                        mvar_data%var(iLookMVAR%scalarFoliageNitrogenFactor)%dat(1), & ! intent(out): foliage nitrogen factor (0-1), 1=saturated
                        mvar_data%var(iLookMVAR%scalarGrowingSeasonIndex)%dat(1),    & ! intent(out): growing season index (0=off, 1=on)
                        mvar_data%var(iLookMVAR%scalarRootZoneTemp)%dat(1),          & ! intent(out): root zone temperature (K)
                        mvar_data%var(iLookMVAR%scalarExposedLAI)%dat(1),            & ! intent(out): exposed leaf area index after burial by snow (m2 m-2)
                        mvar_data%var(iLookMVAR%scalarExposedSAI)%dat(1),            & ! intent(out): exposed stem area index after burial by snow (m2 m-2)

                        ! output: thermal properties
                        mvar_data%var(iLookMVAR%scalarBulkVolHeatCapVeg)%dat(1),     & ! intent(out): bulk volumetric heat capacity of vegetation (J m-3 K-1)
                        mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat,             & ! intent(out): volumetric heat capacity in each layer (J m-3 K-1) 
                        mvar_data%var(iLookMVAR%mLayerThermalC)%dat,                 & ! intent(out): thermal conductivity at the mid-point of each layer (W m-1 K-1)
                        mvar_data%var(iLookMVAR%iLayerThermalC)%dat,                 & ! intent(out): thermal conductivity at the interface of each layer (W m-1 K-1)
                        mvar_data%var(iLookMVAR%mLayerVolFracAir)%dat,               & ! intent(out): volumetric fraction of air in each layer (-)

                        ! output: error control
                        err,cmessage)                                                  ! intent(out): error control

 ! check for errors
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! determine if need to include vegetation in the energy flux routines
 exposedVAI     = mvar_data%var(iLookMVAR%scalarExposedLAI)%dat(1) + mvar_data%var(iLookMVAR%scalarExposedSAI)%dat(1)
 computeVegFlux = (exposedVAI > 0.01_dp)

 ! get an initial canopy temperature if veg just starts protruding through snow on the ground
 if(computeVegFlux)then
  ! (NOTE: if canopy temperature is below absolute zero then canopy was previously buried by snow)
  if(mvar_data%var(iLookMVAR%scalarCanopyTemp)%dat(1) < 0._dp) mvar_data%var(iLookMVAR%scalarCanopyTemp)%dat(1) = forc_data%var(iLookFORCE%airtemp)
 endif

 ! allocate space for snow-soil vectors
 allocate(mLayerTempNew(nLayers),mLayerVolFracIceNew(nLayers),mLayerVolFracLiqNew(nLayers),mLayerMatricHeadNew(nSoil),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for the state variables in the snow-soil vector'; return; endif


 ! *****
 ! wrapper for the picard solver sub-routine...
 ! ********************************************
 ! used as a trick to both
 !  1) save typing; and
 !  2) "protect" variables through the use of the intent attribute
 call picardSolv_muster(&

                        ! input: input variables from picardSolv
                        dt,                                                          & ! intent(in): time step (s)
                        maxiter,                                                     & ! intent(in): maximum number of iterations
                        firstSubstep,                                                & ! intent(in): flag to denote first sub-step
                        computeVegFlux,                                              & ! intent(in): flag to denote if computing energy flux over vegetation

                        ! input: coordinate variables
                        indx_data%var(iLookINDEX%layerType)%dat,                     & ! intent(in): layer type (ix_soil or ix_snow)
                        mvar_data%var(iLookMVAR%mLayerDepth)%dat,                    & ! intent(in): depth of each layer (m)
                        mvar_data%var(iLookMVAR%mLayerHeight)%dat,                   & ! intent(in): height at the mid-point of each layer (m)
                        mvar_data%var(iLookMVAR%iLayerHeight)%dat,                   & ! intent(in): height at the interface of each layer (m)

                        ! input: state variables at the start of the step
                        mvar_data%var(iLookMVAR%scalarCanopyTemp)%dat(1),            & ! intent(in): temperature of the vegetation canopy at the start of the sub-step (K)
                        mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1),             & ! intent(in): mass of ice on the vegetation canopy at the start of the sub-step (kg m-2)
                        mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1),             & ! intent(in): mass of liquid water on the vegetation canopy at the start of the sub-step (kg m-2)
                        mvar_data%var(iLookMVAR%mLayerTemp)%dat,                     & ! intent(in): temperature of each snow/soil layer at the start of the sub-step (K)
                        mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat,               & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
                        mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat,               & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)
                        mvar_data%var(iLookMVAR%mLayerMatricHead)%dat,               & ! intent(in): matric head at the current iteration start of the sub-step (m)

                        ! input: diagnostic variables
                        mvar_data%var(iLookMVAR%scalarBulkVolHeatCapVeg)%dat(1),     & ! intent(in): bulk volumetric heat capacity of vegetation (J m-3 K-1)
                        mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat,             & ! intent(in): volumetric heat capacity in each layer (J m-3 K-1) 

                        ! input: van Genutchen soil parameters
                        mpar_data%var(iLookPARAM%vGn_alpha),                         & ! intent(in): van Genutchen "alpha" parameter (m-1)
                        mpar_data%var(iLookPARAM%vGn_n),                             & ! intent(in): van Genutchen "n" parameter (-)
                        mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1),                 & ! intent(in): van Genutchen "m" parameter (-)
                        mpar_data%var(iLookPARAM%theta_sat),                         & ! intent(in): soil porosity (-)
                        mpar_data%var(iLookPARAM%theta_res),                         & ! intent(in): soil residual volumetric water content (-)

                        ! input: algorithmic control parameters
                        mpar_data%var(iLookPARAM%relConvTol_liquid),                 & ! intent(in): relative convergence tolerance for vol frac liq water (-)
                        mpar_data%var(iLookPARAM%absConvTol_liquid),                 & ! intent(in): absolute convergence tolerance for vol frac liq water (-)
                        mpar_data%var(iLookPARAM%relConvTol_matric),                 & ! intent(in): relative convergence tolerance for matric head (-)
                        mpar_data%var(iLookPARAM%absConvTol_matric),                 & ! intent(in): absolute convergence tolerance for matric head (m)
                        mpar_data%var(iLookPARAM%relConvTol_energy),                 & ! intent(in): relative convergence tolerance for energy (-)
                        mpar_data%var(iLookPARAM%absConvTol_energy),                 & ! intent(in): absolute convergence tolerance for energy (J m-3)
                        mpar_data%var(iLookPARAM%relConvTol_aquifr),                 & ! intent(in): relative convergence tolerance for aquifer storage (-)
                        mpar_data%var(iLookPARAM%absConvTol_aquifr),                 & ! intent(in): absolute convergence tolerance for aquifer storage (m)

                        ! output: diagnostic variables
                        mvar_data%var(iLookMVAR%scalarTemp_CanopyAir)%dat(1),        & ! intent(inout): trial temperature of the canopy air space (K)
                        mvar_data%var(iLookMVAR%scalarVP_CanopyAir)%dat(1),          & ! intent(inout): trial vapor pressure of the canopy air space (Pa)
                        mvar_data%var(iLookMVAR%mLayerTcrit)%dat,                    & ! intent(out): critical soil temperature where liquid water begins to freeze (K)

                        ! output: model state variables at the end of the step
                        ! NOTE: use intent(out) instead of intent(inout) to protect start-of-step variables
                        scalarCanopyTempNew,                                         & ! intent(out): temperature of the vegetation canopy at the end of the sub-step (K)
                        scalarCanopyIceNew,                                          & ! intent(out): mass of ice on the vegetation canopy at the end of the sub-step (kg m-2)
                        scalarCanopyLiqNew,                                          & ! intent(out): mass of liquid water on the vegetation canopy at the end of the sub-step (kg m-2)
                        mLayerTempNew,                                               & ! intent(out): temperature of each snow/soil layer at the end of the sub-step (K)
                        mLayerVolFracIceNew,                                         & ! intent(out): volumetric fraction of ice at the end of the sub-step (-)
                        mLayerVolFracLiqNew,                                         & ! intent(out): volumetric fraction of liquid water at the end of the sub-step (-)
                        mLayerMatricHeadNew,                                         & ! intent(out): matric head at the end of the sub-step (m)
                        scalarAquiferStorageNew,                                     & ! intent(out): aquifer storage at the end of the sub-step (m)

                        ! output: number of iterations
                        niter,                                                       & ! intent(out): number of iterations

                        ! output: error control
                        err,cmessage)                                                  ! intent(out): error control

 ! check for errors
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif


 ! *****
 ! compute sublimation of snow...
 ! ******************************

 if(nSnow > 0._dp)then ! snow layers exist
  ! compute sublimation of snow -- included in densification
  mLayerVolFracIceNew(1) = mLayerVolFracIceNew(1) + &
                            mvar_data%var(iLookMVAR%scalarSnowSublimation)%dat(1)*dt/(mvar_data%var(iLookMVAR%mLayerDepth)%dat(1)*iden_ice)
  ! check that we did not melt/sublimate all of the ice
  if(mLayerVolFracIceNew(1) < 0._dp)then
   message=trim(message)//'sublimated all of the ice'
   err=-20; return ! (negative error code forces a reduction in the length of the sub-step and another trial)
  endif
 endif

 




 ! *****
 ! compute densification...
 ! ************************
 call snwDensify(&

                 ! intent(in): variables
                 dt,                                                     & ! intent(in) time step (s)
                 mLayerTempNew(1:nSnow),                                 & ! intent(in): temperature of each layer (K)
                 mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1:nSnow), & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)

                 ! intent(in): parameters
                 mpar_data%var(iLookPARAM%densScalGrowth),               & ! intent(in): density scaling factor for grain growth (kg-1 m3)
                 mpar_data%var(iLookPARAM%tempScalGrowth),               & ! intent(in): temperature scaling factor for grain growth (K-1)
                 mpar_data%var(iLookPARAM%grainGrowthRate),              & ! intent(in): rate of grain growth (s-1)
                 mpar_data%var(iLookPARAM%densScalOvrbdn),               & ! intent(in): density scaling factor for overburden pressure (kg-1 m3)
                 mpar_data%var(iLookPARAM%tempScalOvrbdn),               & ! intent(in): temperature scaling factor for overburden pressure (K-1)
                 mpar_data%var(iLookPARAM%base_visc),                    & ! intent(in): viscosity coefficient at T=T_frz and snow density=0 (kg m-2 s)

                 ! intent(inout): state variables
                 mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow),      & ! intent(inout): depth of each layer (m)
                 mLayerVolFracLiqNew(1:nSnow),                           & ! intent(inout):  volumetric fraction of liquid water after itertations (-)
                 mLayerVolFracIceNew(1:nSnow),                           & ! intent(inout):  volumetric fraction of ice after itertations (-)

                 ! output: error control
                 err,cmessage)                                             ! intent(out): error control

 ! check for errors
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! update coordinate variables
 call calcHeight(err,cmessage)
 if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif




 ! *****
 ! water balance check...
 ! **********************
 call wBal_check(&
                        
                 ! input: model control
                 dt,                                                                 & ! intent(in): time step (s)
                 mpar_data%var(iLookPARAM%wimplicit),                                & ! intent(in): weight assigned to start-of-step fluxes (-)
                 mvar_data%var(iLookMVAR%mLayerDepth)%dat,                           & ! intent(in): depth of each layer (m)
                 model_decisions(iLookDECISIONS%groundwatr)%iDecision,               & ! intent(in): choice of groundwater parameterization

                 ! input: state variables at the start of the sub-step
                 mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat,                      & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
                 mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat,                      & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)
                 mvar_data%var(iLookMVAR%scalarAquiferStorage)%dat(1),               & ! intent(in): aquifer storage at the end of the sub-step (m)

                 ! input: state variables after iteration                
                 mLayerVolFracIceNew,                                                & ! intent(in): volumetric fraction of ice at the end of the sub-step (-)
                 mLayerVolFracLiqNew,                                                & ! intent(in): volumetric fraction of liquid water at the end of the sub-step (-)
                 scalarAquiferStorageNew,                                            & ! intent(in): aquifer storage at the end of the sub-step (m)

                 ! input: model fluxes at the *START* of the sub-step
                 mvar_data%var(iLookMVAR%iLayerInitLiqFluxSoil)%dat,                 & ! intent(in): liquid water flux at the interface of each layer at the start of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%mLayerInitEjectWater)%dat,                  & ! intent(in): liquid water ejected from each layer at the start of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%mLayerInitBaseflow)%dat,                    & ! intent(in): baseflow from each layer at the start of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%mLayerInitTranspire)%dat,                   & ! intent(in): transpiration from each layer at the start of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%scalarInitAquiferRecharge)%dat(1),          & ! intent(in): recharge to the aquifer at the start of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%scalarInitAquiferBaseflow)%dat(1),          & ! intent(in): baseflow from the aquifer at the start of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%scalarInitAquiferTranspire)%dat(1),         & ! intent(in): transpiration from the aquifer at the sraer of the sub-step (m s-1)

                 ! input: model fluxes at the *END* of the sub-step
                 mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat,                     & ! intent(in): liquid water flux at the interface of each layer at the end of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%mLayerEjectWater)%dat,                      & ! intent(in): liquid water ejected from each layer at the end of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%mLayerBaseflow)%dat,                        & ! intent(in): baseflow from each layer at the end of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%mLayerTranspire)%dat,                       & ! intent(in): transpiration from each layer at the end of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%scalarAquiferRecharge)%dat(1),              & ! intent(in): recharge to the aquifer at the end of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%scalarAquiferBaseflow)%dat(1),              & ! intent(in): baseflow from the aquifer at the end of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%scalarAquiferTranspire)%dat(1),             & ! intent(in): transpiration from the aquifer at the end of the sub-step (m s-1)

                 ! output: water balance of the soil zone
                 mvar_data%var(iLookMVAR%scalarSoilInflux)%dat(1),                   & ! intent(out): influx at the top of the soil zone (m s-1)
                 mvar_data%var(iLookMVAR%scalarSoilBaseflow)%dat(1),                 & ! intent(out): total baseflow from the soil zone (m s-1)
                 mvar_data%var(iLookMVAR%scalarSoilEjection)%dat(1),                 & ! intent(out): total ejected water from all soil layers (m s-1)
                 mvar_data%var(iLookMVAR%scalarSoilDrainage)%dat(1),                 & ! intent(out): drainage from the bottom of the soil profile (m s-1)
                 mvar_data%var(iLookMVAR%scalarSoilTranspiration)%dat(1),            & ! intent(out): total soil transpiration (m s-1)

                 ! output: water balance check
                 mvar_data%var(iLookMVAR%scalarSoilWatBalError)%dat(1),              & ! intent(out): error in the total soil water balance (kg m-2)
                 mvar_data%var(iLookMVAR%scalarAquiferBalError)%dat(1),              & ! intent(out): error in the aquifer water balance (kg m-2)
                 mvar_data%var(iLookMVAR%scalarTotalSoilLiq)%dat(1),                 & ! intent(out): total mass of liquid water in the soil (kg m-2)
                 mvar_data%var(iLookMVAR%scalarTotalSoilIce)%dat(1),                 & ! intent(out): total mass of ice in the soil (kg m-2)

                 ! output: error control
                 err,cmessage)                                                         ! intent(out): error control

 ! check for errors
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif




 ! *****
 ! update states...
 ! ****************

 ! update the states for the canopy
 mvar_data%var(iLookMVAR%scalarCanopyTemp)%dat(1)  = scalarCanopyTempNew      ! temperature of the vegetation canopy (K)
 mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1)   = scalarCanopyIceNew       ! mass of ice on the vegetation canopy at the current iteration (kg m-2)
 mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1)   = scalarCanopyLiqNew       ! mass of liquid water on the vegetation canopy at the current iteration (kg m-2)

 ! update the states for the snow-soil vector
 mvar_data%var(iLookMVAR%mLayerTemp)%dat           = mLayerTempNew            ! New (after heatTransf)
 mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat     = mLayerVolFracIceNew      ! New (after densification)
 mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat     = mLayerVolFracLiqNew      ! New (after densifcaction)
 mvar_data%var(iLookMVAR%mLayerMatricHead)%dat     = mLayerMatricHeadNew      ! New (after densification)
 mvar_data%var(iLookMVAR%scalarAquiferStorage)%dat = scalarAquiferStorageNew  ! New (after groundwater)

 ! compuute total snow depth and SWE
 if(nSnow>0)then
  mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1) = sum(mvar_data%var(iLookMVAR%mLayerDepth)%dat, mask = indx_data%var(iLookINDEX%layerType)%dat==ix_snow)
  mvar_data%var(iLookMVAR%scalarSWE)%dat(1)       = sum(iden_ice*  mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat(1:nSnow)*mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow)) + &  ! total ice
                                                    sum(iden_water*mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat(1:nSnow)*mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow))
 endif

 ! deallocate space for snow-soil vectors
 deallocate(mLayerTempNew,mLayerVolFracIceNew,mLayerVolFracLiqNew,mLayerMatricHeadNew, stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem deallocating space for the state variables in the snow-soil vector'; return; endif

 end subroutine picardSolv



 ! *****************************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************
 ! ***** PRIVATE SUBROUTINES ***************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************************





 ! *********************************************************************************************************
 ! private subroutine: compute preliminary variables needed by the solver (phenology, thermal properties)
 ! *********************************************************************************************************
 subroutine picardSolv_prelim(&

                              ! input: site attributes
                              latitude,                      & ! intent(in): latitude
                              vegTypeIndex,                  & ! intent(in): vegetation type index
                              urbanVegCategory,              & ! intent(in): vegetation category for urban areas               

                              ! input: coordinate variables
                              layerType,                     & ! intent(in): layer type (ix_soil or ix_snow)
                              mLayerHeight,                  & ! intent(in): height at the mid-point of each layer (m)
                              iLayerHeight,                  & ! intent(in): height at the interface of each layer (m)

                              ! input: model state variables
                              scalarCanopyTemp,              & ! intent(in): temperature of the vegetation canopy at the start of the sub-step (K)
                              scalarCanopyIce,               & ! intent(in): mass of ice on the vegetation canopy at the start of the sub-step (kg m-2)
                              scalarCanopyLiq,               & ! intent(in): mass of liquid water on the vegetation canopy at the start of the sub-step (kg m-2)
                              mLayerTemp,                    & ! intent(in): temperature of each snow/soil layer at the start of the sub-step (K)
                              mLayerVolFracIce,              & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
                              mLayerVolFracLiq,              & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)

                              ! input: diagnostic variables
                              yearLength,                    & ! intent(in): number of days in the current year
                              fracJulday,                    & ! intent(in): fractional julian days since the start of year
                              scalarSnowDepth,               & ! intent(in): snow depth on the ground surface (m)
                              scalarLAI,                     & ! intent(inout): one-sided leaf area index (m2 m-2)
                              scalarSAI,                     & ! intent(inout): one-sided stem area index (m2 m-2)

                              ! input: vegetation parameters
                              heightCanopyTop,               & ! intent(in): height at the top of the veg canopy (m)
                              heightCanopyBottom,            & ! intent(in): height at the bottom of the veg canopy (m)
                              specificHeatVeg,               & ! intent(in): specific heat of vegetation (J kg-1 K-1)
                              maxMassVegetation,             & ! intent(in): maximum mass of vegetation (full foliage) (kg m-2)
                              rootingDepth,                  & ! intent(in): rooting depth of the vegetation (m)

                              ! input: soil parameters
                              soil_dens_intr,                & ! intent(in): intrinsic density of soil (kg m-3)
                              theta_sat,                     & ! intent(in): soil porosity (-)
                              frac_sand,                     & ! intent(in): fraction of sand (-)
                              frac_silt,                     & ! intent(in): fraction of silt (-)
                              frac_clay,                     & ! intent(in): fraction of clay (-)

                              ! output: vegetation phenology
                              scalarFoliageNitrogenFactor,   & ! intent(out): foliage nitrogen factor (0-1), 1=saturated
                              scalarGrowingSeasonIndex,      & ! intent(out): growing season index (0=off, 1=on)
                              scalarRootZoneTemp,            & ! intent(out): root zone temperature (K)
                              scalarExposedLAI,              & ! intent(out): exposed leaf area index after burial by snow (m2 m-2)
                              scalarExposedSAI,              & ! intent(out): exposed stem area index after burial by snow (m2 m-2)

                              ! output: thermal properties
                              scalarBulkVolHeatCapVeg,       & ! intent(out): bulk volumetric heat capacity of vegetation (J m-3 K-1)
                              mLayerVolHtCapBulk,            & ! intent(out): volumetric heat capacity in each layer (J m-3 K-1) 
                              mLayerThermalC,                & ! intent(out): thermal conductivity at the mid-point of each layer (W m-1 K-1)
                              iLayerThermalC,                & ! intent(out): thermal conductivity at the interface of each layer (W m-1 K-1)
                              mLayerVolFracAir,              & ! intent(out): volumetric fraction of air in each layer (-)

                              ! output: error control
                              err,message)                     ! intent(out): error control

 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 ! provide access to subroutines
 USE NOAHMP_ROUTINES,only:phenology             ! determine vegetation phenology
 USE diagn_evar_module,only:diagn_evar          ! compute diagnostic energy variables -- thermal conductivity and heat capacity
 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 implicit none
 ! input: site attributes
 real(dp),intent(in)            :: latitude                    ! latitude
 integer(i4b),intent(in)        :: vegTypeIndex                ! vegetation type index
 integer(i4b),intent(in)        :: urbanVegCategory            ! vegetation category for urban areas
 ! input: coordinate variables 
 integer(i4b),intent(in)        :: layerType(:)                ! type of the layer (ix_soil or ix_snow)
 real(dp),intent(in)            :: mLayerHeight(:)             ! height at the mid-point of each layer (m)
 real(dp),intent(in)            :: iLayerHeight(0:)            ! height at the interface of each layer (m)
 ! input: model state variables
 real(dp),intent(in)            :: scalarCanopyTemp            ! temperature of the vegetation canopy (K)
 real(dp),intent(in)            :: scalarCanopyIce             ! mass of ice on the vegetation canopy (kg m-2)
 real(dp),intent(in)            :: scalarCanopyLiq             ! mass of liquid water on the vegetation canopy (kg m-2)
 real(dp),intent(in)            :: mLayerTemp(:)               ! temperature of each layer (K)
 real(dp),intent(in)            :: mLayerVolFracIce(:)         ! volumetric fraction of ice in each layer (-)
 real(dp),intent(in)            :: mLayerVolFracLiq(:)         ! volumetric fraction of liquid water in each layer (-)
 ! input: diagnostic variables
 integer(i4b),intent(in)        :: yearLength                  ! number of days in the current year
 real(dp),intent(in)            :: fracJulday                  ! fractional julian days since the start of year
 real(dp),intent(in)            :: scalarSnowDepth             ! snow depth on the ground surface (m)
 real(dp),intent(inout)         :: scalarLAI                   ! one-sided leaf area index (m2 m-2)
 real(dp),intent(inout)         :: scalarSAI                   ! one-sided stem area index (m2 m-2)
 ! input: vegetation parameters
 real(dp),intent(in)            :: heightCanopyTop             ! height at the top of the veg canopy
 real(dp),intent(in)            :: heightCanopyBottom          ! height at the bottom of the veg canopy
 real(dp),intent(in)            :: specificHeatVeg             ! specific heat of vegetation (J kg-1 K-1)
 real(dp),intent(in)            :: maxMassVegetation           ! maximum mass of vegetation (full foliage) (kg m-2)
 real(dp),intent(in)            :: rootingDepth                ! rooting depth of the vegetation (m)
 ! input: soil parameters
 real(dp),intent(in)            :: soil_dens_intr              ! intrinsic density of soil (kg m-3)
 real(dp),intent(in)            :: theta_sat                   ! soil porosity (-)
 real(dp),intent(in)            :: frac_sand                   ! fraction of sand (-)
 real(dp),intent(in)            :: frac_silt                   ! fraction of silt (-)
 real(dp),intent(in)            :: frac_clay                   ! fraction of clay (-)
 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 ! output: vegetation phenology
 real(dp),intent(out)           :: scalarFoliageNitrogenFactor ! foliage nitrogen concentration (1.0 = saturated)
 real(dp),intent(out)           :: scalarGrowingSeasonIndex    ! growing season index (0=off, 1=on)
 real(dp),intent(out)           :: scalarRootZoneTemp          ! root zone temperature (K)
 real(dp),intent(out)           :: scalarExposedLAI            ! exposed leaf area index after burial by snow (m2 m-2)
 real(dp),intent(out)           :: scalarExposedSAI            ! exposed stem area index after burial by snow (m2 m-2)
 ! output: thermal properties
 real(dp),intent(out)           :: scalarBulkVolHeatCapVeg     ! bulk volumetric heat capacity of vegetation (J m-3 K-1)
 real(dp),intent(out)           :: mLayerVolHtCapBulk(:)       ! volumetric heat capacity in each layer (J m-3 K-1) 
 real(dp),intent(out)           :: mLayerThermalC(:)           ! thermal conductivity at the mid-point of each layer (W m-1 K-1)
 real(dp),intent(out)           :: iLayerThermalC(0:)          ! thermal conductivity at the interface of each layer (W m-1 K-1)
 real(dp),intent(out)           :: mLayerVolFracAir(:)         ! volumetric fraction of air in each layer (-)
 ! output: error control
 integer(i4b),intent(out)       :: err                         ! error code
 character(*),intent(out)       :: message                     ! error message
 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 ! local variables
 character(len=256)             :: cmessage                    ! error message of downwind routine
 integer(i4b)                   :: nLayersRoots                ! number of soil layers that contain roots
 real(dp),parameter             :: verySmall=epsilon(1._dp)    ! a very small number
 real(dp)                       :: notUsed_heightCanopyTop     ! for some reason the Noah-MP phenology routines output canopy height
 real(dp)                       :: exposedVAI                  ! exposed vegetation area index (LAI + SAI)
 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="picardSolv_prelim/"

 ! compute the root zone temperature (used in vegetation phenology)
 ! (compute the number of layers with roots)
 nLayersRoots = count(iLayerHeight(nSnow:nLayers-1) < rootingDepth-verySmall)
 if(nLayersRoots == 0)then; err=20; message=trim(message)//'no roots within the soil profile'; return; endif
 ! (compute the temperature of the root zone)
 scalarRootZoneTemp = sum(mLayerTemp(nSnow+1:nSnow+nLayersRoots)) / real(nLayersRoots, kind(dp))

 ! define the foliage nitrogen factor
 scalarFoliageNitrogenFactor = 1._dp  ! foliage nitrogen concentration (1.0 = saturated)

 ! determine vegetation phenology
 ! NOTE: recomputing phenology every sub-step accounts for changes in exposed vegetation associated with changes in snow depth
 call phenology(&
                ! input
                vegTypeIndex,               & ! intent(in): vegetation type index
                urbanVegCategory,           & ! intent(in): vegetation category for urban areas               
                scalarSnowDepth,            & ! intent(in): snow depth on the ground surface (m)
                scalarCanopyTemp,           & ! intent(in): temperature of the vegetation canopy (K)
                latitude,                   & ! intent(in): latitude (degrees north)
                yearLength,                 & ! intent(in): number of days in the current year
                fracJulday,                 & ! intent(in): fractional julian days since the start of year
                scalarLAI,                  & ! intent(inout): one-sided leaf area index (m2 m-2)
                scalarSAI,                  & ! intent(inout): one-sided stem area index (m2 m-2)
                scalarRootZoneTemp,         & ! intent(in): root zone temperature
                ! output
                notUsed_heightCanopyTop,    & ! intent(out): height of the top of the canopy layer (m)
                scalarExposedLAI,           & ! intent(out): exposed leaf area index after burial by snow (m2 m-2)
                scalarExposedSAI,           & ! intent(out): exposed stem area index after burial by snow (m2 m-2)
                scalarGrowingSeasonIndex    ) ! intent(out): growing season index (0=off, 1=on)

 ! compute diagnostic energy variables (thermal conductivity and volumetric heat capacity)
 call diagn_evar(&
                 ! input: state variables
                 ! NOTE: using start-of-substep variables
                 scalarCanopyIce,           & ! intent(in): mass of ice on the vegetation canopy (kg m-2)
                 scalarCanopyLiq,           & ! intent(in): mass of liquid water on the vegetation canopy (kg m-2)
                 mLayerVolFracIce,          & ! intent(in): volumetric fraction of ice in each layer (-)
                 mLayerVolFracLiq,          & ! intent(in): volumetric fraction of liquid water in each layer (-)
                 ! input: coordinate variables
                 layerType,                 & ! intent(in): layer type (ix_soil or ix_snow)
                 mLayerHeight,              & ! intent(in): height at the mid-point of each layer (m)
                 iLayerHeight,              & ! intent(in): height at the interface of each layer (m)
                 ! input: model parameters
                 heightCanopyTop,           & ! intent(in): height of top of the vegetation canopy above ground surface (m)
                 heightCanopyBottom,        & ! intent(in): height of bottom of the vegetation canopy above ground surface (m)
                 specificHeatVeg,           & ! intent(in): specific heat of vegetation (J kg-1 K-1)
                 maxMassVegetation,         & ! intent(in): maximum mass of vegetation (full foliage) (kg m-2)
                 soil_dens_intr,            & ! intent(in): intrinsic density of soil (kg m-3)
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

 end subroutine picardSolv_prelim


 ! *****************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************





 ! *********************************************************************************************************
 ! private subroutine: wrapper subroutine for the picard solver itself
 ! *********************************************************************************************************
 subroutine picardSolv_muster(&

                              ! input: input variables from picardSolv
                              dt,                            & ! intent(in): time step (s)
                              maxiter,                       & ! intent(in): maximum number of iterations
                              firstSubstep,                  & ! intent(in): flag to denote first sub-step
                              computeVegFlux,                & ! intent(in): flag to denote if computing energy flux over vegetation

                              ! input: coordinate variables
                              layerType,                     & ! intent(in): layer type (ix_soil or ix_snow)
                              mLayerDepth,                   & ! intent(in): depth of each layer (m)
                              mLayerHeight,                  & ! intent(in): height at the mid-point of each layer (m)
                              iLayerHeight,                  & ! intent(in): height at the interface of each layer (m)

                              ! input: model state variables
                              scalarCanopyTemp,              & ! intent(in): temperature of the vegetation canopy at the start of the sub-step (K)
                              scalarCanopyIce,               & ! intent(in): mass of ice on the vegetation canopy at the start of the sub-step (kg m-2)
                              scalarCanopyLiq,               & ! intent(in): mass of liquid water on the vegetation canopy at the start of the sub-step (kg m-2)
                              mLayerTemp,                    & ! intent(in): temperature of each snow/soil layer at the start of the sub-step (K)
                              mLayerVolFracIce,              & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
                              mLayerVolFracLiq,              & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)
                              mLayerMatricHead,              & ! intent(in): matric head at the start of the sub-step (-)

                              ! input: diagnostic variables
                              scalarBulkVolHeatCapVeg,       & ! intent(in): volumetric heat capacity of vegetation (J m-3 K-1)
                              mLayerVolHtCapBulk,            & ! intent(in): volumetric heat capacity of each layer in the snow-soil vector (J m-3 K-1)

                              ! input: van Genutchen soil parameters
                              vGn_alpha,                     & ! intent(in): van Genutchen "alpha" parameter (m-1)
                              vGn_n,                         & ! intent(in): van Genutchen "n" parameter (-)
                              VGn_m,                         & ! intent(in): van Genutchen "m" parameter (-)
                              theta_sat,                     & ! intent(in): soil porosity (-)
                              theta_res,                     & ! intent(in): soil residual volumetric water content (-)

                              ! input: algorithmic control parameters
                              relConvTol_liquid,             & ! intent(in): relative convergence tolerance for vol frac liq water (-)
                              absConvTol_liquid,             & ! intent(in): absolute convergence tolerance for vol frac liq water (-)
                              relConvTol_matric,             & ! intent(in): relative convergence tolerance for matric head (-)
                              absConvTol_matric,             & ! intent(in): absolute convergence tolerance for matric head (m)
                              relConvTol_energy,             & ! intent(in): relative convergence tolerance for energy (-)
                              absConvTol_energy,             & ! intent(in): absolute convergence tolerance for energy (J m-3)
                              relConvTol_aquifr,             & ! intent(in): relative convergence tolerance for aquifer storage (-)
                              absConvTol_aquifr,             & ! intent(in): absolute convergence tolerance for aquifer storage (m)
                               
                              ! output: diagnostic variables
                              scalarTemp_CanopyAir,          & ! intent(inout): trial temperature of the canopy air space (K)
                              scalarVP_CanopyAir,            & ! intent(inout): trial vapor pressure of the canopy air space (Pa)
                              mLayerTcrit,                   & ! intent(out): critical soil temperature where liquid water begins to freeze (K)

                              ! output: model state variables at the end of the step
                              ! NOTE: use intent(out) instead of intent(inout) to protect start-of-step variables
                              scalarCanopyTempNew,           & ! intent(out): temperature of the vegetation canopy at the end of the sub-step (K)
                              scalarCanopyIceNew,            & ! intent(out): mass of ice on the vegetation canopy at the end of the sub-step (kg m-2)
                              scalarCanopyLiqNew,            & ! intent(out): mass of liquid water on the vegetation canopy at the end of the sub-step (kg m-2)
                              mLayerTempNew,                 & ! intent(out): temperature of each snow/soil layer at the end of the sub-step (K)
                              mLayerVolFracIceNew,           & ! intent(out): volumetric fraction of ice at the end of the sub-step (-)
                              mLayerVolFracLiqNew,           & ! intent(out): volumetric fraction of liquid water at the end of the sub-step (-)
                              mLayerMatricHeadNew,           & ! intent(out): matric head at the end of the sub-step (m)
                              scalarAquiferStorageNew,       & ! intent(out): aquifer storage at the end of the sub-step (m)

                              ! output: number of iterations
                              niter,                         & ! intent(out): number of iterations

                              ! output: error control
                              err,message)                     ! intent(out): error control

 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 ! provide access to subroutines
 USE heatTransf_module,only:heatTransf          ! compute change in temperature over the time step
 USE phseChange_module,only:phseChange          ! compute change in phase over the time step
 USE can_Hydrol_module,only:can_Hydrol          ! compute canopy water balance
 USE snowHydrol_module,only:snowHydrol          ! compute liquid water flow through the snowpack
 USE soilHydrol_module,only:soilHydrol          ! compute change in mass over the time step for the soil
 USE snwDensify_module,only:snwDensify          ! compute densification of snow
 USE soil_utils_module,only:crit_soilT          ! compute the critical temperature above which all water is unfrozen
 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 implicit none
 ! input
 real(dp),intent(in)            :: dt                          ! time step (seconds)
 integer(i4b),intent(in)        :: maxiter                     ! maximum number of iterations
 logical(lgt),intent(in)        :: firstSubStep                ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(in)        :: computeVegFlux              ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! input: coordinate variables 
 integer(i4b),intent(in)        :: layerType(:)                ! type of the layer (ix_soil or ix_snow)
 real(dp),intent(in)            :: mLayerDepth(:)              ! depth of each layer (m)
 real(dp),intent(in)            :: mLayerHeight(:)             ! height at the mid-point of each layer (m)
 real(dp),intent(in)            :: iLayerHeight(0:)            ! height at the interface of each layer (m)
 ! input: model state variables
 real(dp),intent(in)            :: scalarCanopyTemp            ! temperature of the vegetation canopy (K)
 real(dp),intent(in)            :: scalarCanopyIce             ! mass of ice on the vegetation canopy (kg m-2)
 real(dp),intent(in)            :: scalarCanopyLiq             ! mass of liquid water on the vegetation canopy (kg m-2)
 real(dp),intent(in)            :: mLayerTemp(:)               ! temperature of each layer (K)
 real(dp),intent(in)            :: mLayerVolFracIce(:)         ! volumetric fraction of ice in each layer (-)
 real(dp),intent(in)            :: mLayerVolFracLiq(:)         ! volumetric fraction of liquid water in each layer (-)
 real(dp),intent(in)            :: mLayerMatricHead(:)         ! matric head in each layer (-)
 ! input: diagnostic variables
 real(dp),intent(in)            :: scalarBulkVolHeatCapVeg     ! bulk volumetric heat capacity of vegetation (J m-3 K-1)
 real(dp),intent(in)            :: mLayerVolHtCapBulk(:)       ! volumetric heat capacity in each layer (J m-3 K-1) 
 ! input: van Genutchen soil parameters
 real(dp),intent(in)            :: vGn_alpha                   ! van Genutchen "alpha" parameter (m-1)
 real(dp),intent(in)            :: vGn_n                       ! van Genutchen "n" parameter (-)
 real(dp),intent(in)            :: vGn_m                       ! van Genutchen "m" parameter (-)
 real(dp),intent(in)            :: theta_sat                   ! soil porosity (-)
 real(dp),intent(in)            :: theta_res                   ! soil residual volumetric water content (-)
 ! input: algorithmic control parameters
 real(dp),intent(in)            :: relConvTol_liquid           ! relative convergence tolerance for vol frac liq water (-)
 real(dp),intent(in)            :: absConvTol_liquid           ! absolute convergence tolerance for vol frac liq water (-)
 real(dp),intent(in)            :: relConvTol_matric           ! relative convergence tolerance for matric head (-)
 real(dp),intent(in)            :: absConvTol_matric           ! absolute convergence tolerance for matric head (m)
 real(dp),intent(in)            :: relConvTol_energy           ! relative convergence tolerance for energy (-)
 real(dp),intent(in)            :: absConvTol_energy           ! absolute convergence tolerance for energy (J m-3)
 real(dp),intent(in)            :: relConvTol_aquifr           ! relative convergence tolerance for aquifer storage (-)
 real(dp),intent(in)            :: absConvTol_aquifr           ! absolute convergence tolerance for aquifer storage (m)
 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 ! output: diagnostic variables
 real(dp),intent(inout)         :: scalarTemp_CanopyAir        ! trial temperature of the canopy air space (K)
 real(dp),intent(inout)         :: scalarVP_CanopyAir          ! trial vapor pressure of the canopy air space (Pa)
 real(dp),intent(out)           :: mLayerTcrit(:)              ! critical soil temperature where liquid water begins to freeze (K)
 ! output: model state variables at the end of the step
 ! NOTE: use intent(out) instead of intent(inout) to protect start-of-step variables
 real(dp),intent(out)           :: scalarCanopyTempNew         ! temperature of the vegetation canopy at the end of the sub-step (K)
 real(dp),intent(out)           :: scalarCanopyIceNew          ! mass of ice on the vegetation canopy at the end of the sub-step (kg m-2)
 real(dp),intent(out)           :: scalarCanopyLiqNew          ! mass of liquid water on the vegetation canopy at the end of the sub-step (kg m-2)
 real(dp),intent(out)           :: mLayerTempNew(:)            ! temperature of each snow/soil layer at the end of the sub-step (K)
 real(dp),intent(out)           :: mLayerVolFracIceNew(:)      ! volumetric fraction of ice at the end of the sub-step (-)
 real(dp),intent(out)           :: mLayerVolFracLiqNew(:)      ! volumetric fraction of liquid water at the end of the sub-step (-)
 real(dp),intent(out)           :: mLayerMatricHeadNew(:)      ! matric head at the current iteration end of the sub-step (m)
 real(dp),intent(out)           :: scalarAquiferStorageNew  ! aquifer storage at the end of the sub-step (m)
 ! output: number of iterations
 integer(i4b),intent(out)       :: niter                       ! number of iterations
 ! output: error control
 integer(i4b),intent(out)       :: err                         ! error code
 character(*),intent(out)       :: message                     ! error message
 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 ! define general local variables
 character(len=256)             :: cmessage                    ! error message of downwind routine
 logical(lgt),parameter         :: printflag=.false.           ! flag to print results to screen (debugging)
 logical(lgt),parameter         :: freeze_infiltrate=.true.    ! flag to freeze infiltrating water
 integer(i4b),parameter         :: minLayer=1                  ! minimum layer to print
 integer(i4b),parameter         :: maxLayer=4                  ! maximum layer to print
 integer(i4b)                   :: iLayer                      ! loop through model layers
 integer(i4b)                   :: iter                        ! iteration index
 real(dp)                       :: theta                       ! liquid water equivalent of the volumetric fraction of total water, liquid plus ice (-)
 real(dp)                       :: volFrac_water               ! total volumetric fraction of water, liquid water plus ice (-) 
 real(dp),parameter             :: eps=1.e-10_dp               ! small increment used to define ice content at the freezing point
 ! define derivatives in canopy air space variables
 real(dp)                       :: dTempCanopyAir_dTCanopy     ! derivative in the temperature of the canopy air space w.r.t. temperature of the canopy
 real(dp)                       :: dTempCanopyAir_dTGround     ! derivative in the temperature of the canopy air space w.r.t. temperature of the ground
 real(dp)                       :: dVPCanopyAir_dTCanopy       ! derivative in the vapor pressure of the canopy air space w.r.t. temperature of the canopy 
 real(dp)                       :: dVPCanopyAir_dTGround       ! derivative in the vapor pressure of the canopy air space w.r.t. temperature of the ground
 ! define state variables for the vegetation canopy
 real(dp)                       :: scalarCanopyTempIter        ! trial value of temperature of the vegetation canopy (K)
 real(dp)                       :: scalarCanopyIceIter         ! trial value of mass of ice on the vegetation canopy (kg m-2)
 real(dp)                       :: scalarCanopyLiqIter         ! trial value of mass of liquid water on the vegetation canopy (kg m-2)
 ! define state variables for the snow-soil system
 real(dp),dimension(nLayers)    :: mLayerTempIter              ! temperature vector for all snow/soil layers at the current iteration (K)
 real(dp),dimension(nLayers)    :: mLayerVolFracIceIter        ! volumetric fraction of ice for all snow/soil layers at the current iteration (-)
 real(dp),dimension(nLayers)    :: mLayerVolFracLiqIter        ! volumetric fraction of liquid water for all snow/soil layers at the current iteration (-)
 real(dp),dimension(nSoil)      :: mLayerMatricHeadIter        ! matric head for all soil layers at the current iteration (-)
 ! define state variables for the aquifer
 real(dp)                       :: scalarAquiferStorageIter    ! trial value of aquifer storage (m)
 ! define model diagnostic variables
 real(dp),dimension(nLayers)    :: mLayerMeltFreeze            ! energy associated with melting/freezing ice (J m-3)
 real(dp),dimension(nLayers)    :: mLayerInfilFreeze           ! energy associated with infiltrating water (J m-3)
 ! define iteration increments
 real(dp)                       :: scalarCanopyTempIncr        ! iteration increment for temperature of the vegetation canopy (K)
 real(dp),dimension(nLayers)    :: mLayerTempIncr              ! iteration increment for temperature of the snow-soil vector (K)
 real(dp),dimension(nLayers)    :: mLayerLiqIncr               ! iteration increment for volumetric liquid water content in the snow-soil vector (-)
 real(dp),dimension(nSoil)      :: mLayerMatIncr               ! iteration increment for matric head in soil layers (m)
 real(dp),dimension(nLayers)    :: mLayerTempIncrOld           ! iteration increment for temperature of the snow-soil vector in the previous iteration (K)
 ! define error monitoring variables
 real(dp)                       :: canopyTempComponent         ! veg canopy: temperature component of the energy increment (J m-3)
 real(dp),dimension(nLayers)    :: mLayerTempComponent         ! snow-soil vector: temperature component of the energy increment (J m-3)
 real(dp),dimension(nLayers)    :: mLayerPhseComponent         ! snow-soil vector: phase component of the energy increment (J m-3)
 real(dp),dimension(nLayers)    :: mLayerInflComponent         ! snow-soil vector: infiltration component of the energy increment (J m-3)
 real(dp)                       :: canopyNrgIncr               ! change in energy in the vegetation canopy from one iteration to the next (J m-3)
 real(dp),dimension(nLayers)    :: mLayerNrgIncr               ! change in energy in the snow-soil vector from one iteration to the next (J m-3)
 real(dp)                       :: scalarAqiIncr               ! change in aquifer storage from one iteration to the next (m)
 ! define variables to assess iteration convergence
 integer(i4b),dimension(1)      :: liquid_pos                  ! position of maximum error
 integer(i4b),dimension(1)      :: matric_pos                  ! position of maximum error
 integer(i4b),dimension(1)      :: energy_pos                  ! position of maximum error
 real(dp),dimension(1)          :: liquid_max                  ! maximum absolute change in volumetric liquid water content for a given iteration (-)
 real(dp),dimension(1)          :: matric_max                  ! maximum absolute change in matric head for a given iteration (m)
 real(dp),dimension(1)          :: energy_max                  ! maximum absolute change in energy for a given iteration (J m-3)
 real(dp)                       :: aquifr_max                  ! absolute change in aquifer storage for a given iteration (m)
 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="picardSolv_muster/"

 ! initialize number of iterations
 niter=0

 ! *****
 ! * initialize...
 ! ***************

 ! initialize canopy temperature
 scalarCanopyTempIter = scalarCanopyTemp

 ! initialize canopy water
 scalarCanopyIceIter = scalarCanopyIce   ! mass of ice on the vegetation canopy (kg m-2)
 scalarCanopyLiqIter = scalarCanopyLiq   ! mass of liquid water on the vegetation canopy (kg m-2)

 ! initialize layer temperatures
 mLayerTempIter = mLayerTemp

 ! initialize volumetric liquid and ice content
 mLayerVolFracIceIter = mLayerVolFracIce   ! volumetric ice content (-)
 mLayerVolFracLiqIter = mLayerVolFracLiq   ! volumetric liquid water content (-)

 ! initialize matric head
 mLayerMatricHeadIter = mLayerMatricHead   ! matric lead (m)

 ! calculate the critical soil temperature above which all water is unfrozen (K)
 do iLayer=nSnow+1,nSoil
  theta = mLayerVolFracIce(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiq(iLayer)
  mLayerTcrit(iLayer-nSnow) = crit_soilT(theta,theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m)
 end do

 ! initialize the melt/freeze vectors
 mLayerMeltFreeze  = 0._dp
 mLayerInfilFreeze = 0._dp

 !if(any(mLayerVolFracIce*iden_ice > 700._dp)) printflag=.true.

 ! **********************************************************************************************************************
 ! **********************************************************************************************************************
 ! **********************************************************************************************************************
 ! ***** iterate ********************************************************************************************************
 ! **********************************************************************************************************************
 ! **********************************************************************************************************************
 ! **********************************************************************************************************************
 do iter=1,maxiter

  ! *****
  ! * thermodynamics...
  ! *******************

  ! increment number of iterations
  niter=niter+1

  ! compute the temperature and ice content at the next iteration
  call heatTransf(&
                  ! input
                  dt,&                        ! intent(in): time step (seconds)
                  iter,&                      ! intent(in): current iteration count
                  firstSubstep,             & ! intent(in): flag to indicate if we are processing the first sub-step
                  computeVegFlux,           & ! intent(in): flag to indicate if we computing fluxes ovser vegetation (.false. means veg is buried with snow)
                  scalarCanopyTempIter,     & ! intent(in): trial temperature of the vegetation canopy at the current iteration (K)
                  scalarCanopyIceIter,      & ! intent(in): trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
                  scalarCanopyLiqIter,      & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                  mLayerTempIter,           & ! intent(in): trial temperature of each snow/soil layer at the current iteration (K)
                  mLayerVolFracIceIter,     & ! intent(in): trial volumetric fraction of ice in each snow/soil layer at the current iteration (-)
                  mLayerVolFracLiqIter,     & ! intent(in): trial volumetric fraction of liquid water in each snow/soil layer at the current iteration (-)
                  mLayerMatricHeadIter,     & ! intent(in): trial matric head of each snow/soil layer at the current iteration (m)

                  ! input/output variables from heatTransf subroutine: canopy air space variables
                  scalarTemp_CanopyAir,     & ! intent(inout): trial temperature of the canopy air space (K)
                  scalarVP_CanopyAir,       & ! intent(inout): trial vapor pressure of the canopy air space (Pa)
                  dTempCanopyAir_dTCanopy,  & ! intent(inout): derivative in the temperature of the canopy air space w.r.t. temperature of the canopy
                  dTempCanopyAir_dTGround,  & ! intent(inout): derivative in the temperature of the canopy air space w.r.t. temperature of the ground
                  dVPCanopyAir_dTCanopy,    & ! intent(inout): derivative in the vapor pressure of the canopy air space w.r.t. temperature of the canopy 
                  dVPCanopyAir_dTGround,    & ! intent(inout): derivative in the vapor pressure of the canopy air space w.r.t. temperature of the ground

                  ! output
                  scalarCanopyTempIncr,     & ! intent(out): iteration increment for temperature of the vegetation canopy (K)
                  mLayerTempIncr,           & ! intent(out): iteration increment for temperature of the snow-soil system (K)
                  scalarCanopyTempNew,      & ! intent(out): new temperature of the vegetation canopy (K)
                  mLayerTempNew,            & ! intent(out): new temperature each snow/soil layer (K)
                  mLayerVolFracIceNew,      & ! intent(out): new volumetric fraction of ice in each snow/soil layer (-)
                  mLayerVolFracLiqNew,      & ! intent(out): new volumetric fraction of liquid water in each snow/soil layer (-)
                  mLayerMatricHeadNew,      & ! intent(out): new matric head in each snow/soil layer (m)
                  
                  err,cmessage)               ! intent(out): error control
  ! check errors
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! test
  write(*,'(a,1x,10(f12.5,1x))') 'scalarTemp_CanopyAir, scalarVP_CanopyAir, scalarCanopyTempNew, mLayerTempNew(1), scalarCanopyTempIncr, mLayerTempIncr(1) = ', &
                                  scalarTemp_CanopyAir, scalarVP_CanopyAir, scalarCanopyTempNew, mLayerTempNew(1), scalarCanopyTempIncr, mLayerTempIncr(1)
  pause

  ! compute melt/freeze in each layer (kg m-3 s-1) -- melt is negative
  mLayerMeltFreeze = mLayerMeltFreeze + iden_ice*(mLayerVolFracIceNew - mLayerVolFracIceIter)/dt

  ! compute the temperature and phase components of the energy increment (J m-3)
  canopyTempComponent = scalarBulkVolHeatCapVeg*scalarCanopyTempIncr
  mLayerTempComponent = mLayerVolHtCapBulk*mLayerTempIncr
  mLayerPhseComponent = iden_ice*LH_fus*(mLayerVolFracIceNew - mLayerVolFracIceIter)

  ! save the temperature increment
  mLayerTempIncrOld   = mLayerTempIncr

  ! =================================================================================================================================================
  ! =================================================================================================================================================


  ! *****
  ! * hydrology...
  ! **************

  ! account for phase change
  mLayerMatricHeadIter = mLayerMatricHeadNew
  mLayerVolFracLiqIter = mLayerVolFracLiqNew

  ! compute the canopy water balance
  call can_Hydrol(&
                  ! input
                  dt,                       & ! intent(in): time step (seconds)
                  iter,                     & ! intent(in): iteration index
                  scalarCanopyIceIter,      & ! intent(in): trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
                  scalarCanopyLiqIter,      & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                  ! output
                  scalarCanopyIceNew,       & ! intent(out): updated mass of ice on the vegetation canopy at the current iteration (kg m-2)
                  scalarCanopyLiqNew,       & ! intent(out): updated mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                  err,cmessage)               ! intent(out): error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! compute the volumetric liquid water content at the next iteration (note: only use snow vectors)
  ! NOTE: ice not modified in the snow hydrology routines, so can stay as "New"
  if(nSnow > 0)then
   call snowHydrol(dt,                               & ! time step (seconds)
                   iter,                             & ! iteration index
                   mLayerVolFracLiqIter(1:nSnow),    & ! volumetric fraction of liquid water at the current iteration (-)
                   mLayerVolFracIceNew(1:nSnow),     & ! volumetric fraction of ice at the current iteration (-)
                   mLayerVolFracLiqNew(1:nSnow),     & ! volumetric fraction of liquid water at the next iteration (-)
                   err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  ! compute the matric head at the next iteration (note liquid water and ice vectors are defined for all layers)
  ! NOTE: ice not modified in the soil hydrology routines, so can stay as "New"
  call soilHydrol(&
                  ! input
                  dt,&                                          ! time step (seconds)
                  iter,&                                        ! current iteration count
                  mLayerMatricHeadIter(1:nLevels),            & ! matric head in each layer at the current iteration (m)
                  mLayerVolFracLiqIter(nSnow+1:nSnow+nLevels),& ! volumetric fraction of liquid water at the current iteration (-)
                  mLayerVolFracIceNew(nSnow+1:nSnow+nLevels), & ! volumetric fraction of ice at the current iteration (-)
                  scalarAquiferStorageIter,                   & ! aquifer storage (m)
                  !output
                  mLayerMatricHeadNew(1:nLevels),             & ! matric head in each layer at the next iteration (m)
                  mLayerVolFracLiqNew(nSnow+1:nSnow+nLevels), & ! volumetric fraction of liquid water at the next iteration (-)
                  scalarAquiferStorageNew,                    & ! aquifer storage (m)
                  err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! copy across "saturated" layers
  if(nLevels < nSoil)then
   mLayerMatricHeadNew(nLevels+1:nSoil)         = mLayerMatricHeadIter(nLevels+1:nSoil)
   mLayerVolFracLiqNew(nSnow+nLevels+1:nLayers) = mLayerVolFracLiqIter(nSnow+nLevels+1:nLayers)
  endif

  ! compute the iteration increment for the matric head and volumetric fraction of liquid water
  mLayerMatIncr = mLayerMatricHeadNew - mLayerMatricHeadIter 
  mLayerLiqIncr = mLayerVolFracLiqNew - mLayerVolFracLiqIter 

  ! compute the iteration increment for aquifer storage
  scalarAqiIncr = scalarAquiferStorageNew - scalarAquiferStorageIter

  ! update the state variables -- used in phase change below
  mLayerMatricHeadIter = mLayerMatricHeadNew
  mLayerVolFracLiqIter = mLayerVolFracLiqNew
  mLayerVolFracIceIter = mLayerVolFracIceNew

  ! =================================================================================================================================================
  ! =================================================================================================================================================


  ! *****
  ! * phase change.....
  ! *******************

  ! calculate the critical soil temperature above which all water is unfrozen (K)
  do iLayer=nSnow+1,nLayers
   theta = mLayerVolFracIceNew(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiqNew(iLayer)
   mLayerTcrit(iLayer-nSnow) = crit_soilT(theta,theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m)
  end do

  ! option to compute phase change associated with infiltrating liquid water
  if(freeze_infiltrate)then
   ! compute phase change associated with infiltrating liquid water
   call phsechange(mLayerTempNew,        & ! intent(in): new temperature vector (K)
                   mLayerMatricHeadIter, & ! intent(in): matric head at the current iteration (m)
                   mLayerVolFracLiqIter, & ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                   mLayerVolFracIceIter, & ! intent(in): volumetric fraction of ice at the current iteration (-)
                   mLayerMatricHeadNew,  & ! intent(out): new matric head (m)
                   mLayerVolFracLiqNew,  & ! intent(out): new volumetric fraction of liquid water (-)
                   mLayerVolFracIceNew,  & ! intent(out): new volumetric fraction of ice (-)
                   err,cmessage)           ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  else 
   mLayerMatricHeadNew = mLayerMatricHeadIter
   mLayerVolFracLiqNew = mLayerVolFracLiqIter
   mLayerVolFracIceNew = mLayerVolFracIceIter
  endif

  ! compute melt/freeze of infiltrating liquid water in each layer (kg m-3 s-1) -- melt is negative
  mLayerInfilFreeze = mLayerInfilFreeze + iden_ice*(mLayerVolFracIceIter - mLayerVolFracIceNew)/dt

  ! compute increment in energy for the vegetation canopy
  canopyNrgIncr       = canopyTempComponent

  ! compute increment in energy for the snow-soil vector -- do here because of phase change
  mLayerInflComponent = iden_ice*LH_fus*(mLayerVolFracIceNew - mLayerVolFracIceIter)
  mLayerNrgIncr       = mLayerTempComponent - mLayerPhseComponent - mLayerInflComponent

  ! =================================================================================================================================================
  ! =================================================================================================================================================


  ! *****
  ! * get ready for the next iteration.....
  ! ***************************************

  ! update scalar variables
  scalarCanopyTempIter = scalarCanopyTempNew
  scalarAquiferStorageIter = scalarAquiferStorageNew

  ! update the state variables -- used in phase change below
  mLayerTempIter       = mLayerTempNew
  mLayerMatricHeadIter = mLayerMatricHeadNew
  mLayerVolFracLiqIter = mLayerVolFracLiqNew
  mLayerVolFracIceIter = mLayerVolFracIceNew


  ! =================================================================================================================================================
  ! =================================================================================================================================================


  ! *****
  ! * check convergence.....
  ! ************************

  ! non-iterative check (do not expect convergence)
  if(maxiter==1) exit  ! NOTE: exit loop here to avoid return statement with error code

  ! compute maximum iteration increment
  liquid_max = maxval(abs(mLayerLiqIncr))
  matric_max = maxval(abs(mLayerMatIncr))
  energy_max = maxval(abs((/canopyNrgIncr,mLayerNrgIncr/)))
  aquifr_max = abs(scalarAqiIncr)

  ! get position of maximum iteration increment
  liquid_pos = maxloc(abs(mLayerLiqIncr))
  matric_pos = maxloc(abs(mLayerMatIncr))
  energy_pos = maxloc(abs((/canopyNrgIncr,mLayerNrgIncr/)))

  ! test
  !write(*,'(a25,1x,i4,1x,10(e20.3,1x))') 'temperature increment = ', energy_pos, scalarCanopyTempIncr, mLayerTempIncr(minLayer:maxLayer)
  !write(*,'(a25,1x,i4,1x,10(f20.7,1x))') 'energy increment = ', energy_pos, canopyNrgIncr, mLayerNrgIncr(minLayer:maxLayer)
  !pause

  ! convergence check: 
  if(liquid_max(1) < absConvTol_liquid .and. &   ! volumetric fraction of liquid water (-)
     matric_max(1) < absConvTol_matric .and. &   ! matric head (m)
     energy_max(1) < absConvTol_energy .and. &   ! energy (J m-3)
     aquifr_max    < absConvTol_aquifr)      &   ! aquifer storage (m)
   exit

  ! check for lack of convergence
  if(niter==maxiter)then; err=-30; message=trim(message)//'failed to converge'; return; endif

 end do  ! (iterating)


 ! *****************************************************************************************************************************************
 ! *****************************************************************************************************************************************
 ! ***** END OF ITERATIONS *****************************************************************************************************************
 ! *****************************************************************************************************************************************
 ! *****************************************************************************************************************************************
 !pause 'after iterations'

 ! =================================================================================================================================================
 ! =================================================================================================================================================


 ! *****
 ! * basic checks.....
 ! *******************

 ! check that we did not melt all of the ice
 do iLayer=1,nSnow
  if(mLayerVolFracIceNew(iLayer) < 0._dp)then
   message=trim(message)//'melted all of the ice'
   err=-20; return ! (negative error code forces a reduction in the length of the sub-step and another trial)
  endif
 end do

 ! check matric head is not ridiculous
 do iLayer=1,nSoil
  if(mLayerMatricHeadNew(iLayer) > 100._dp)then
   write(message,'(a,i0,a,f9.1,a)')trim(message)//"matric head > 100 [iLayer=",iLayer,"; matricHead=",&
         mLayerMatricHeadNew(iLayer),"]"
   err=20; return
  endif
 end do

 ! ** check that total volumetric water (liquid water plus ice) does not exceed porosity
 do iLayer=nSnow+1,nLayers
  ! compute total volumetric fraction filled with water (liquid plus ice)
  volFrac_water = mLayerVolFracIceNew(iLayer) + mLayerVolFracLiqNew(iLayer)
  ! check if the total volumetric water (liquid water plus ice) exceeds porosity
  if(volFrac_water > theta_sat)then
   write(*,'(a,i4,1x,2(f15.10,1x))') 'iLayer, mLayerVolFracIceNew(iLayer), mLayerVolFracLiqNew(iLayer) = ',&
                                      iLayer, mLayerVolFracIceNew(iLayer), mLayerVolFracLiqNew(iLayer)
   write(message,'(a,i0,a,i0,a)')trim(message)//"(liquid + ice) > porosity [iLayer=",iLayer,"; iSoil=",iLayer-nSnow,"]"
   err=20; return
  endif  ! (if the total volumetric water -- liquid water plus ice -- exceeds porosity)
 end do ! (looping through soil layers)

 ! =================================================================================================================================================
 ! =================================================================================================================================================


 end subroutine picardSolv_muster


 ! *****************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************





 ! *********************************************************************************************************
 ! private subroutine: check the water balance
 ! *********************************************************************************************************
 subroutine wBal_check(&

                       ! input: model control
                       dt,                                    & ! intent(in): time step (s)
                       wimplicit,                             & ! intent(in): weight assigned to start-of-step fluxes (-)
                       mLayerDepth,                           & ! intent(in): depth of each layer (m)
                       ix_groundwatr,                         & ! intent(in): choice of groundwater parameterization

                       ! input: state variables at the start of the sub-step
                       mLayerVolFracIce,                      & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
                       mLayerVolFracLiq,                      & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)
                       scalarAquiferStorage,                  & ! intent(in): aquifer storage at the start of the sub-step (m)

                       ! input: state variables after iteration                
                       mLayerVolFracIceNew,                   & ! intent(in): volumetric fraction of ice at the end of the sub-step (-)
                       mLayerVolFracLiqNew,                   & ! intent(in): volumetric fraction of liquid water at the end of the sub-step (-)
                       scalarAquiferStorageNew,               & ! intent(in): aquifer storage at the end of the sub-step (m)

                       ! input: model fluxes at the *START* of the sub-step
                       iLayerInitLiqFluxSoil,                 & ! intent(in): liquid water flux at the interface of each layer at the start of the sub-step (m s-1)
                       mLayerInitEjectWater,                  & ! intent(in): liquid water ejected from each layer at the start of the sub-step (m s-1)
                       mLayerInitBaseflow,                    & ! intent(in): baseflow from each layer at the start of the sub-step (m s-1)
                       mLayerInitTranspire,                   & ! intent(in): transpiration from each layer at the start of the sub-step (m s-1)
                       scalarInitAquiferRecharge,             & ! intent(in): recharge to the aquifer at the start of the sub-step (m s-1)
                       scalarInitAquiferBaseflow,             & ! intent(in): baseflow from the aquifer at the start of the sub-step (m s-1)
                       scalarInitAquiferTranspire,            & ! intent(in): transpiration from the aquifer at the start of the sub-step (m s-1)

                       ! input: model fluxes at the *END* of the sub-step
                       iLayerLiqFluxSoil,                     & ! intent(in): liquid water flux at the interface of each layer at the end of the sub-step (m s-1)
                       mLayerEjectWater,                      & ! intent(in): liquid water ejected from each layer at the end of the sub-step (m s-1)
                       mLayerBaseflow,                        & ! intent(in): baseflow from each layer at the end of the sub-step (m s-1)
                       mLayerTranspire,                       & ! intent(in): transpiration from each layer at the end of the sub-step (m s-1)
                       scalarAquiferRecharge,                 & ! intent(in): recharge to the aquifer at the end of the sub-step (m s-1)
                       scalarAquiferBaseflow,                 & ! intent(in): baseflow from the aquifer at the end of the sub-step (m s-1)
                       scalarAquiferTranspire,                & ! intent(in): transpiration from the aquifer at the end of the sub-step (m s-1)

                       ! output: water balance of the soil zone
                       scalarSoilInflux,                      & ! intent(out): influx at the top of the soil zone (m s-1)
                       scalarSoilBaseflow,                    & ! intent(out): total baseflow from the soil zone (m s-1)
                       scalarSoilEjection,                    & ! intent(out): total ejected water from all soil layers (m s-1)
                       scalarSoilDrainage,                    & ! intent(out): drainage from the bottom of the soil profile (m s-1)
                       scalarSoilTranspiration,               & ! intent(out): total soil transpiration (m s-1)

                       ! output: water balance check
                       scalarSoilWatBalError,                 & ! intent(out): error in the total soil water balance (kg m-2)
                       scalarAquiferBalError,                 & ! intent(out): error in the aquifer water balance (kg m-2)
                       scalarTotalSoilLiq,                    & ! intent(out): total mass of liquid water in the soil (kg m-2)
                       scalarTotalSoilIce,                    & ! intent(out): total mass of ice in the soil (kg m-2)

                       ! output: error control
                       err,message)                             ! intent(out): error control

 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 ! look-up values for the choice of groundwater parameterization
 USE mDecisions_module,only:       &
  equilWaterTable,                 & ! equilibrium water table
  pseudoWaterTable,                & ! pseudo water table
  bigBucket,                       & ! a big bucket (lumped aquifer model)
  noExplicit                         ! no explicit groundwater parameterization
 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 real(dp),intent(in)            :: dt                           ! time step (seconds)
 real(dp),intent(in)            :: wimplicit                    ! weight assigned to start-of-step fluxes (-)
 real(dp),intent(in)            :: mLayerDepth(:)               ! depth of each layer (m)
 integer(i4b),intent(in)        :: ix_groundwatr                ! choice of groundwater parameterization
 ! input: state variables at the start of the sub-step
 real(dp),intent(in)            :: mLayerVolFracIce(:)          ! volumetric fraction of ice in each layer (-)
 real(dp),intent(in)            :: mLayerVolFracLiq(:)          ! volumetric fraction of liquid water in each layer (-)
 real(dp),intent(in)            :: scalarAquiferStorage         ! aquifer storage at the start of the sub-step (m)
 ! input: state variables after iteration                
 real(dp),intent(in)            :: mLayerVolFracIceNew(:)       ! volumetric fraction of ice at the end of the sub-step (-)
 real(dp),intent(in)            :: mLayerVolFracLiqNew(:)       ! volumetric fraction of liquid water at the end of the sub-step (-)
 real(dp),intent(in)            :: scalarAquiferStorageNew      ! aquifer storage at the end of the sub-step (m)
 ! input: model fluxes at the *START* of the sub-step
 real(dp),intent(in)            :: iLayerInitLiqFluxSoil(0:)    ! liquid water flux at the interface of each layer at the start of the sub-step (m s-1)
 real(dp),intent(in)            :: mLayerInitEjectWater(:)      ! liquid water ejected from each layer at the start of the sub-step (m s-1)
 real(dp),intent(in)            :: mLayerInitBaseflow(:)        ! baseflow from each layer at the start of the sub-step (m s-1)
 real(dp),intent(in)            :: mLayerInitTranspire(:)       ! transpiration from each layer at the start of the sub-step (m s-1)
 real(dp),intent(in)            :: scalarInitAquiferRecharge    ! recharge to the aquifer at the end of the sub-step (m s-1)
 real(dp),intent(in)            :: scalarInitAquiferBaseflow    ! baseflow from the aquifer at the end of the sub-step (m s-1)
 real(dp),intent(in)            :: scalarInitAquiferTranspire   ! transpiration from the aquifer at the end of the sub-step (m s-1)
 ! input: model fluxes at the *END* of the sub-step
 real(dp),intent(in)            :: iLayerLiqFluxSoil(0:)        ! liquid water flux at the interface of each layer at the start of the sub-step (m s-1)
 real(dp),intent(in)            :: mLayerEjectWater(:)          ! liquid water ejected from each layer at the start of the sub-step (m s-1)
 real(dp),intent(in)            :: mLayerBaseflow(:)            ! baseflow from each layer at the start of the sub-step (m s-1)
 real(dp),intent(in)            :: mLayerTranspire(:)           ! transpiration from each layer at the start of the sub-step (m s-1)
 real(dp),intent(in)            :: scalarAquiferRecharge        ! recharge to the aquifer at the end of the sub-step (m s-1)
 real(dp),intent(in)            :: scalarAquiferBaseflow        ! baseflow from the aquifer at the end of the sub-step (m s-1)
 real(dp),intent(in)            :: scalarAquiferTranspire       ! transpiration from the aquifer at the end of the sub-step (m s-1)
 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 ! output: water balance of the soil zone
 real(dp),intent(out)           :: scalarSoilInflux             ! intent(out): influx at the top of the soil zone (m s-1)
 real(dp),intent(out)           :: scalarSoilBaseflow           ! intent(out): total baseflow from the soil zone (m s-1)
 real(dp),intent(out)           :: scalarSoilEjection           ! intent(out): total ejected water from all soil layers (m s-1)
 real(dp),intent(out)           :: scalarSoilDrainage           ! intent(out): drainage from the bottom of the soil profile (m s-1)
 real(dp),intent(out)           :: scalarSoilTranspiration      ! intent(out): total soil transpiration (m s-1)
 ! output: water balance check
 real(dp),intent(out)           :: scalarSoilWatBalError        ! error in the total soil water balance (kg m-2)
 real(dp),intent(out)           :: scalarAquiferBalError        ! error in the aquifer water balance (kg m-2)
 real(dp),intent(out)           :: scalarTotalSoilLiq           ! total mass of liquid water in the soil (kg m-2)
 real(dp),intent(out)           :: scalarTotalSoilIce           ! total mass of ice in the soil (kg m-2)
 ! output: error control
 integer(i4b),intent(out)       :: err                          ! error code
 character(*),intent(out)       :: message                      ! error message
 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 ! local: define balance check variables
 integer(i4b)                   :: iLayer                       ! index of model layers
 real(dp)                       :: totalChange                  ! total change in volumetric liquid water content over a layer 
 real(dp)                       :: phaseChange                  ! change in volumetric liquid water content associated with phase change
 real(dp)                       :: flux_Change                  ! change in volumetric liquid water content associated with fluxes
 real(dp)                       :: evap_Change                  ! change in volumetric liquid water content associated with transpiration
 real(dp)                       :: qbaseChange                  ! change in volumetric liquid water content associated with baseflow
 real(dp)                       :: ejectChange                  ! change in volumetric liquid water content associated with ejection of water
 real(dp)                       :: balanceSoilWater0            ! total soil storage at the start of the step (kg m-2)
 real(dp)                       :: balanceSoilWater1            ! total soil storage at the end of the step (kg m-2)
 real(dp)                       :: balanceSoilInflux            ! input to the soil zone
 real(dp)                       :: balanceSoilBaseflow          ! output from the soil zone
 real(dp)                       :: balanceSoilDrainage          ! output from the soil zone
 real(dp)                       :: balanceSoilEjection          ! output from the soil zone
 real(dp)                       :: balanceSoilTranspiration     ! output from the soil zone
 real(dp)                       :: balanceAquifer0              ! total aquifer storage at the start of the step (kg m-2)
 real(dp)                       :: balanceAquifer1              ! total aquifer storage at the end of the step (kg m-2)
 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="wBal_check/"

 ! compute total soil moisture and ice at the *START* of the step (kg m-2)
 scalarTotalSoilLiq = sum(iden_water*mLayerVolFracLiq(nSnow+1:nLayers)*mLayerDepth(nSnow+1:nLayers))
 scalarTotalSoilIce = sum(iden_ice  *mLayerVolFracIce(nSnow+1:nLayers)*mLayerDepth(nSnow+1:nLayers))

 ! get the total water in the soil (liquid plus ice) at the start of the time step (kg m-2)
 balanceSoilWater0 = scalarTotalSoilLiq + scalarTotalSoilIce

 ! get the total aquifer storage at the start of the time step (kg m-2)
 balanceAquifer0 = scalarAquiferStorage*iden_water

 ! compute total soil moisture and ice at the *END* of the step (kg m-2)
 scalarTotalSoilLiq = sum(iden_water*mLayerVolFracLiqNew(nSnow+1:nLayers)*mLayerDepth(nSnow+1:nLayers))
 scalarTotalSoilIce = sum(iden_ice  *mLayerVolFracIceNew(nSnow+1:nLayers)*mLayerDepth(nSnow+1:nLayers))

 ! get the total water in the soil (liquid plus ice) at the end of the time step (kg m-2)
 balanceSoilWater1 = scalarTotalSoilLiq + scalarTotalSoilIce

 ! get the total aquifer storage at the start of the time step (kg m-2)
 balanceAquifer1 = scalarAquiferStorageNew*iden_water

 ! compute the influx, drainage and ejection of water from the soil profile (m s-1)
 scalarSoilInflux        = (wimplicit*iLayerInitLiqFluxSoil(0)       + (1._dp - wimplicit)*iLayerLiqFluxSoil(0)      )
 scalarSoilBaseflow      = (wimplicit*sum(mLayerInitBaseflow)        + (1._dp - wimplicit)*sum(mLayerBaseflow)       )
 scalarSoilDrainage      = (wimplicit*iLayerInitLiqFluxSoil(nLevels) + (1._dp - wimplicit)*iLayerLiqFluxSoil(nLevels))
 scalarSoilEjection      = (wimplicit*sum(mLayerInitEjectWater)      + (1._dp - wimplicit)*sum(mLayerEjectWater)     )
 scalarSoilTranspiration = (wimplicit*sum(mLayerInitTranspire)       + (1._dp - wimplicit)*sum(mLayerTranspire)) - &
                           (wimplicit*scalarInitAquiferTranspire     + (1._dp - wimplicit)*scalarAquiferTranspire)
 ! check ejected water
 if(scalarSoilEjection < 0._dp)then
  print*, 'mLayerInitEjectWater = ', mLayerInitEjectWater
  print*, 'mLayerEjectWater = ', mLayerEjectWater
  message=trim(message)//'ejected water < 0'; err=20; return
 endif

 ! get the input and output to/from the soil zone (kg m-2)
 balanceSoilInflux        = scalarSoilInflux*iden_water*dt
 balanceSoilBaseflow      = scalarSoilBaseflow*iden_water*dt
 balanceSoilDrainage      = scalarSoilDrainage*iden_water*dt
 balanceSoilEjection      = scalarSoilEjection*iden_water*dt
 balanceSoilTranspiration = scalarSoilTranspiration*iden_water*dt

 ! check the soil water balance
 scalarSoilWatBalError  = balanceSoilWater1 - (balanceSoilWater0 + (balanceSoilInflux + balanceSoilTranspiration - balanceSoilBaseflow - balanceSoilDrainage - balanceSoilEjection) )
 if(abs(scalarSoilWatBalError) > 1.d-3)then
  ! check the balance of each layer
  write(*,'(a)') 'water balance of each layer'
  write(*,'(a)') 'Temp0 (K), Temp1(K), Liq0 (-), Liq1 (-), Ice0 (-), Ice1 (-), totalChange (-), phaseChange (-), flux_Change, evap_Change, ejectChange, ',&
                 'phaseChange+flux_Change+evap_Change+ejectChange, totalChange - (phaseChange+flux_Change+evap_Change-ejectChange)'
  do iLayer=1,nSoil
   totalChange = mLayerVolFracLiqNew(iLayer+nSnow) - mLayerVolFracLiq(iLayer+nSnow) ! total change in volumetric liquid water content
   phaseChange = -(iden_ice/iden_water)*(mLayerVolFracIceNew(iLayer+nSnow) - mLayerVolFracIce(iLayer+nSnow))  ! change in liquid water content associated with freezing
   evap_Change = dt*mLayerTranspire(iLayer)/mLayerDepth(iLayer)
   qbaseChange = dt*mLayerBaseflow(iLayer)/mLayerDepth(iLayer)
   ejectChange = dt*mLayerEjectWater(iLayer)/mLayerDepth(iLayer)
   flux_Change = dt*(iLayerLiqFluxSoil(iLayer-1) - iLayerLiqFluxSoil(iLayer))/mLayerDepth(iLayer) ! change in volumetric liquid water content from the interface fluxes
   write(*,'(i4,1x,2(f15.8,1x),20(e15.5,1x))') iLayer, mLayerVolFracLiq(iLayer+nSnow), mLayerVolFracLiqNew(iLayer+nSnow), &
                                                       mLayerVolFracIce(iLayer+nSnow), mLayerVolFracIceNew(iLayer+nSnow), &
                                                       totalChange, phaseChange, flux_Change, evap_Change, qbaseChange, ejectChange, &
                                                       phaseChange+flux_Change+evap_Change-qbaseChange-ejectChange, totalChange - (phaseChange+flux_Change+evap_Change-qbaseChange-ejectChange)
  end do
  ! print the total water balance
  print*, 'dt = ', dt
  print*, 'balanceSoilWater0 (kg m-2) = ',        balanceSoilWater0
  print*, 'balanceSoilWater1 (kg m-2) = ',        balanceSoilWater1
  print*, 'balanceSoilInflux (kg m-2) = ',        balanceSoilInflux
  print*, 'balanceSoilDrainage (kg m-2) = ',      balanceSoilDrainage
  print*, 'balanceSoilEjection (kg m-2) = ',      balanceSoilEjection
  print*, 'balanceSoilTranspiration (kg m-2) = ', balanceSoilTranspiration, sum(mLayerTranspire)*dt*iden_water
  write(message,'(a,e20.10,a)')trim(message)//"abs(scalarSoilWatBalError) > 1.d-3 [error = ",&
                               scalarSoilWatBalError," ]"
  err=10; return
 endif

 ! check the aquifer water balance
 select case(ix_groundwatr)
  ! no explicit aquifer
  case(noExplicit,equilWaterTable)
   scalarAquiferBalError = 0._dp
  ! explicit aquifer
  case(bigBucket,pseudoWaterTable)
   ! check the aquifer water balance
   scalarAquiferBalError = balanceAquifer1 - (balanceAquifer0 + (iden_water*(wimplicit*scalarInitAquiferTranspire + (1._dp - wimplicit)*scalarAquiferTranspire)*dt) + &
                                                                (iden_water*(wimplicit*scalarInitAquiferRecharge  + (1._dp - wimplicit)*scalarAquiferRecharge) *dt) - &
                                                                (iden_water*(wimplicit*scalarInitAquiferBaseflow  + (1._dp - wimplicit)*scalarAquiferBaseflow) *dt) )
   ! print the terms in the aquifer balance if errors are sufficiently large
   if(abs(scalarAquiferBalError) > 1.d-3)then
    write(*,'(a,f20.10)') 'scalarAquiferBalError  = ', scalarAquiferBalError
    write(*,'(a,f20.10)') 'balanceAquifer1        = ', balanceAquifer1
    write(*,'(a,f20.10)') 'balanceAquifer0        = ', balanceAquifer0
    write(*,'(a,f20.10)') 'scalarAquiferTranspire = ', iden_water*(wimplicit*scalarInitAquiferTranspire + (1._dp - wimplicit)*scalarAquiferTranspire)*dt
    write(*,'(a,f20.10)') 'scalarAquiferRecharge  = ', iden_water*(wimplicit*scalarInitAquiferRecharge  + (1._dp - wimplicit)*scalarAquiferRecharge) *dt
    write(*,'(a,f20.10)') 'scalarAquiferBaseflow  = ', iden_water*(wimplicit*scalarInitAquiferBaseflow  + (1._dp - wimplicit)*scalarAquiferBaseflow) *dt
    write(message,'(a,e20.10,a)')trim(message)//"abs(scalarAquiferBalError) > 1.d-3 [error = ",scalarAquiferBalError," ]"
    err=10; return
   endif
  case default; err=20; message=trim(message)//'unknown groundwater parameterization'; return
 end select ! (selecting groundwater parameterization)

 end subroutine wBal_check


end module picardSolv_module

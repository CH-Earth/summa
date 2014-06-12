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
! look-up values for the choice of groundwater parameterization
USE mDecisions_module,only:  &
 qbaseTopmodel,              & ! TOPMODEL-ish baseflow parameterization
 bigBucket,                  & ! a big bucket (lumped aquifer model)
 noExplicit                    ! no explicit groundwater parameterization
! look-up values for the form of Richards' equation
USE mDecisions_module,only:  &
 moisture,                   & ! moisture-based form of Richards' equation
 mixdform                      ! mixed form of Richards' equation
! look-up values for the choice of boundary conditions for hydrology
USE mDecisions_module,only:  &
 prescribedHead,             & ! prescribed head (volumetric liquid water content for mixed form of Richards' eqn)
 funcBottomHead,             & ! function of matric head in the lower-most layer
 freeDrainage,               & ! free drainage
 liquidFlux,                 & ! liquid water flux
 zeroFlux                      ! zero flux
implicit none
private
public::picardSolv
! number of soil and snow layers
integer(i4b)        :: nSoil        ! number of soil layers
integer(i4b)        :: nSnow        ! number of snow layers
integer(i4b)        :: nLayers      ! total number of layers
integer(i4b)        :: nLevels      ! number of soil layers to use in the soil hydrology routine
! control parameters
real(dp),parameter  :: valueMissing=-9999._dp     ! missing value
real(dp),parameter  :: verySmall=1.e-6_dp         ! a very small number
contains

 ! ************************************************************************************************
 ! new subroutine: run the coupled energy-mass model for one timestep
 ! ************************************************************************************************
 subroutine picardSolv(&
                       ! input
                       dt,            & ! time step (s)
                       maxiter,       & ! maximum number of iterations
                       firstSubstep,  & ! flag to denote first sub-step
                       computeVegFlux,& ! flag to denote if computing energy flux over vegetation
                       ! output
                       niter,         & ! number of iterations taken
                       err,message)     ! error control
 ! model variables, parameters, forcing data, etc.
 USE data_struc,only:attr_data,type_data,mpar_data,forc_data,mvar_data,indx_data    ! data structures
 USE var_lookup,only:iLookATTR,iLookTYPE,iLookPARAM,iLookFORCE,iLookMVAR,iLookINDEX ! named variables for structure elements
 ! model decision structures
 USE data_struc,only:model_decisions    ! model decision structure
 USE var_lookup,only:iLookDECISIONS     ! named variables for elements of the decision structure
 ! snow densification
 USE snwDensify_module,only:snwDensify
 implicit none
 ! input
 real(dp),intent(in)                  :: dt                       ! time step (seconds)
 integer(i4b),intent(in)              :: maxiter                  ! maximum number of iterations
 logical(lgt),intent(in)              :: firstSubStep             ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(in)              :: computeVegFlux           ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! output
 integer(i4b),intent(out)             :: niter                    ! number of iterations
 integer(i4b),intent(out)             :: err                      ! error code
 character(*),intent(out)             :: message                  ! error message
 ! local
 character(LEN=256)                   :: cmessage                 ! error message of downwind routine
 real(dp)                             :: scalarCanairTempNew      ! temperature of the canopy air space at the end of the sub-step (K)
 real(dp)                             :: scalarCanopyTempNew      ! temperature of the vegetation canopy at the end of the sub-step (K)
 real(dp)                             :: scalarCanopyIceNew       ! mass of ice on the vegetation canopy at the end of the sub-step (kg m-2)
 real(dp)                             :: scalarCanopyLiqNew       ! mass of liquid water on the vegetation canopy at the end of the sub-step (kg m-2)
 real(dp),allocatable                 :: mLayerTempNew(:)         ! temperature of each snow/soil layer at the end of the sub-step (K)
 real(dp),allocatable                 :: mLayerVolFracIceNew(:)   ! volumetric fraction of ice at the end of the sub-step (-)
 real(dp),allocatable                 :: mLayerVolFracLiqNew(:)   ! volumetric fraction of liquid water at the end of the sub-step (-)
 real(dp),allocatable                 :: mLayerMatricHeadNew(:)   ! matric head at the end of the sub-step (m)
 real(dp)                             :: scalarAquiferStorageNew  ! aquifer storage at the end of the sub-step (m)
 real(dp)                             :: scalarCanopyWater        ! total storage of wate rin the canopy; liquidwater plus ice (kg m-2)
 real(dp)                             :: volSub                   ! volumetric sublimation (kg m-3)
 real(dp)                             :: testSWE                  ! test SWE (kg m-2) 

 ! initialize error control
 err=0; message="picardSolv/"

 ! define the total number of snow+soil layers
 nLayers = indx_data%var(iLookINDEX%nLayers)%dat(1) 

 ! identify the number of snow and soil layers
 nSnow   = count(indx_data%var(iLookINDEX%layerType)%dat == ix_snow)
 nSoil   = count(indx_data%var(iLookINDEX%layerType)%dat == ix_soil)

 ! identify the number of soil layers to use in the soil hydrology routine
 nLevels = nSoil  ! NOTE: always pass the full number of soil layers

 ! initialize the liquid flux variables (used in the energy routines, which is called before the hydrology routines)
 if(nSnow > 0)&
 mvar_data%var(iLookMVAR%iLayerLiqFluxSnow)%dat(0:nSnow) = 0._dp
 mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat(0:nSoil) = 0._dp

 ! get an initial canopy temperature if veg just starts protruding through snow on the ground
 if(computeVegFlux)then
  ! (NOTE: if canopy temperature is below absolute zero then canopy was previously buried by snow)
  if(mvar_data%var(iLookMVAR%scalarCanopyTemp)%dat(1) < 0._dp .or. &
     mvar_data%var(iLookMVAR%scalarCanairTemp)%dat(1) < 0._dp)then
   ! check there is snow (there really has to be)
   if(nSnow == 0)then
    message=trim(message)//'no snow when canopy temperature or canopy air temperature is undefined -- canopy temps can only be undefined when buried with snow'
    err=20; return
   endif
   ! set canopy temperature to the temperature of the top snow layer + small offset to check derivative calculations
   mvar_data%var(iLookMVAR%scalarCanairTemp)%dat(1) = mvar_data%var(iLookMVAR%mLayerTemp)%dat(1) + 0.1_dp
   mvar_data%var(iLookMVAR%scalarCanopyTemp)%dat(1) = mvar_data%var(iLookMVAR%mLayerTemp)%dat(1) + 0.1_dp 
  endif  ! (if canopy temperature undefined -- means canopy previously buried with snow)
 endif  ! (if computing vegetation fluxes -- canopy exposed)

 ! allocate space for snow-soil vectors
 allocate(mLayerTempNew(nLayers),mLayerVolFracIceNew(nLayers),mLayerVolFracLiqNew(nLayers),mLayerMatricHeadNew(nSoil),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for the state variables in the snow-soil vector'; return; endif

 if(mvar_data%var(iLookMVAR%scalarCanopyIceMax)%dat(1) < 0._dp)then
  message=trim(message)//'maximum canopy ice is less than zero'
  err=20; return
 endif

 !print*, 'before picard_muster, ice max = ', mvar_data%var(iLookMVAR%scalarCanopyIceMax)%dat(1)

 ! *****
 ! wrapper for the picard solver sub-routine...
 ! ********************************************
 ! used as a trick to both
 !  1) save typing; and
 !  2) "protect" variables through the use of the intent attribute
 call picardSolv_muster(&

                        ! input: input variables from picardSolv
                        dt,                                                              & ! intent(in): time step (s)
                        maxiter,                                                         & ! intent(in): maximum number of iterations
                        firstSubstep,                                                    & ! intent(in): flag to denote first sub-step
                        computeVegFlux,                                                  & ! intent(in): flag to denote if computing energy flux over vegetation

                        ! input: upper boundary conditions
                        model_decisions(iLookDECISIONS%bcUpprSoiH)%iDecision,            & ! intent(in): choice of upper boundary conditions for soil hydrology
                        mpar_data%var(iLookPARAM%upperBoundHead),                        & ! intent(in): upper boundary condition for matric head (m)
                        mpar_data%var(iLookPARAM%upperBoundTheta),                       & ! intent(in): upper boundary condition for volumetric liquid water content (-)

                        ! input: coordinate variables
                        indx_data%var(iLookINDEX%layerType)%dat,                         & ! intent(in): layer type (ix_soil or ix_snow)
                        mvar_data%var(iLookMVAR%mLayerDepth)%dat,                        & ! intent(in): depth of each layer (m)
                        mvar_data%var(iLookMVAR%mLayerHeight)%dat,                       & ! intent(in): height at the mid-point of each layer (m)
                        mvar_data%var(iLookMVAR%iLayerHeight)%dat,                       & ! intent(in): height at the interface of each layer (m)

                        ! input: state variables at the start of the step
                        mvar_data%var(iLookMVAR%scalarCanairTemp)%dat(1),                & ! intent(in): temperature of the canopy air space at the start of the sub-step (K)
                        mvar_data%var(iLookMVAR%scalarCanopyTemp)%dat(1),                & ! intent(in): temperature of the vegetation canopy at the start of the sub-step (K)
                        mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1),                 & ! intent(in): mass of ice on the vegetation canopy at the start of the sub-step (kg m-2)
                        mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1),                 & ! intent(in): mass of liquid water on the vegetation canopy at the start of the sub-step (kg m-2)
                        mvar_data%var(iLookMVAR%mLayerTemp)%dat,                         & ! intent(in): temperature of each snow/soil layer at the start of the sub-step (K)
                        mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat,                   & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
                        mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat,                   & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)
                        mvar_data%var(iLookMVAR%mLayerMatricHead)%dat,                   & ! intent(in): matric head at the current iteration start of the sub-step (m)
                        mvar_data%var(iLookMVAR%scalarAquiferStorage)%dat(1),            & ! intent(in): aquifer storage at the start of the sub-step (m)

                        ! input: forcing and diagnostic variables
                        forc_data%var(iLookFORCE%airtemp),                               & ! intent(in): air temperature (K)
                        mvar_data%var(iLookMVAR%scalarBulkVolHeatCapVeg)%dat(1),         & ! intent(in): bulk volumetric heat capacity of vegetation (J m-3 K-1)
                        mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat,                 & ! intent(in): volumetric heat capacity in each layer (J m-3 K-1) 

                        ! input: snow parameters
                        mpar_data%var(iLookPARAM%snowfrz_scale),                         & ! intent(in): scaling parameter for the snow freezing curve (K-1)

                        ! input: van Genutchen soil parameters
                        mpar_data%var(iLookPARAM%vGn_alpha),                             & ! intent(in): van Genutchen "alpha" parameter (m-1)
                        mpar_data%var(iLookPARAM%vGn_n),                                 & ! intent(in): van Genutchen "n" parameter (-)
                        mvar_data%var(iLookMVAR%scalarVGn_m)%dat(1),                     & ! intent(in): van Genutchen "m" parameter (-)
                        mpar_data%var(iLookPARAM%theta_sat),                             & ! intent(in): soil porosity (-)
                        mpar_data%var(iLookPARAM%theta_res),                             & ! intent(in): soil residual volumetric water content (-)

                        ! input: algorithmic control parameters
                        mpar_data%var(iLookPARAM%relConvTol_liquid),                     & ! intent(in): relative convergence tolerance for vol frac liq water (-)
                        mpar_data%var(iLookPARAM%absConvTol_liquid),                     & ! intent(in): absolute convergence tolerance for vol frac liq water (-)
                        mpar_data%var(iLookPARAM%relConvTol_matric),                     & ! intent(in): relative convergence tolerance for matric head (-)
                        mpar_data%var(iLookPARAM%absConvTol_matric),                     & ! intent(in): absolute convergence tolerance for matric head (m)
                        mpar_data%var(iLookPARAM%relConvTol_energy),                     & ! intent(in): relative convergence tolerance for energy (-)
                        mpar_data%var(iLookPARAM%absConvTol_energy),                     & ! intent(in): absolute convergence tolerance for energy (J m-3)
                        mpar_data%var(iLookPARAM%relConvTol_aquifr),                     & ! intent(in): relative convergence tolerance for aquifer storage (-)
                        mpar_data%var(iLookPARAM%absConvTol_aquifr),                     & ! intent(in): absolute convergence tolerance for aquifer storage (m)

                        ! output: diagnostic variables
                        mvar_data%var(iLookMVAR%scalarVP_CanopyAir)%dat(1),              & ! intent(out): vapor pressure of the canopy air space (Pa)
                        mvar_data%var(iLookMVAR%scalarCanopyStabilityCorrection)%dat(1), & ! intent(out): stability correction for the canopy (-)
                        mvar_data%var(iLookMVAR%scalarGroundStabilityCorrection)%dat(1), & ! intent(out): stability correction for the ground surface (-)
                        mvar_data%var(iLookMVAR%mLayerTcrit)%dat,                        & ! intent(out): critical soil temperature where liquid water begins to freeze (K)
                        mvar_data%var(iLookMVAR%scalarCanopyMeltFreeze)%dat(1),          & ! intent(out): melt/freeze of water stored in the canopy (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarCanopyTranspiration)%dat(1),       & ! intent(out): canopy transpiration (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarCanopyEvaporation)%dat(1),         & ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
                        mvar_data%var(iLookMVAR%scalarGroundEvaporation)%dat(1),         & ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)

                        ! output: model state variables at the end of the step
                        ! NOTE: use intent(out) instead of intent(inout) to protect start-of-step variables
                        scalarCanairTempNew,                                             & ! intent(out): temperature of the canopy air space at the end of the sub-step (K)
                        scalarCanopyTempNew,                                             & ! intent(out): temperature of the vegetation canopy at the end of the sub-step (K)
                        scalarCanopyIceNew,                                              & ! intent(out): mass of ice on the vegetation canopy at the end of the sub-step (kg m-2)
                        scalarCanopyLiqNew,                                              & ! intent(out): mass of liquid water on the vegetation canopy at the end of the sub-step (kg m-2)
                        mLayerTempNew,                                                   & ! intent(out): temperature of each snow/soil layer at the end of the sub-step (K)
                        mLayerVolFracIceNew,                                             & ! intent(out): volumetric fraction of ice at the end of the sub-step (-)
                        mLayerVolFracLiqNew,                                             & ! intent(out): volumetric fraction of liquid water at the end of the sub-step (-)
                        mLayerMatricHeadNew,                                             & ! intent(out): matric head at the end of the sub-step (m)
                        scalarAquiferStorageNew,                                         & ! intent(out): aquifer storage at the end of the sub-step (m)
  
                        ! output: number of iterations
                        niter,                                                           & ! intent(out): number of iterations

                        ! output: error control
                        err,cmessage)                                                      ! intent(out): error control

 ! check for errors
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 !print*, 'after picard_muster, ice max = ', mvar_data%var(iLookMVAR%scalarCanopyIceMax)%dat(1)
 !testSWE = sum( (mLayerVolFracLiqNew(1:nSnow)*iden_water + mLayerVolFracIceNew(1:nSnow)*iden_ice)*mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow) )
 !print*, 'after picardSolv_muster: SWE = ', testSWE


 ! *****
 ! compute melt for the snow without a layer...
 ! ********************************************
 call implctMelt(&
                 ! input/output: integrated snowpack properties
                 mvar_data%var(iLookMVAR%scalarSWE)%dat(1),          & ! snow water equivalent (kg m-2)
                 mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1),    & ! snow depth (m)
                 mvar_data%var(iLookMVAR%scalarSfcMeltPond)%dat(1),  & ! surface melt pond (kg m-2)
                 ! input/output: properties of the upper-most soil layer
                 mvar_data%var(iLookMVAR%mLayerTemp)%dat(1),         & ! surface layer temperature (K)
                 mvar_data%var(iLookMVAR%mLayerDepth)%dat(1),        & ! surface layer depth (m)
                 mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat(1), & ! surface layer volumetric heat capacity (J m-3 K-1)
                 ! output: error control
                 err,cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 testSWE = sum( (mLayerVolFracLiqNew(1:nSnow)*iden_water + mLayerVolFracIceNew(1:nSnow)*iden_ice)*mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow) )
 !print*, 'after implctMelt: SWE = ', testSWE

 ! *****
 ! compute change in ice content due to canopy sublimation...
 ! **********************************************************

 if(computeVegFlux)then
  ! check that sublimation does not exceed the available water
  ! NOTE: the sublimation flux is positive downward
  scalarCanopyWater = mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1) + mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1)
  if(-mvar_data%var(iLookMVAR%scalarCanopySublimation)%dat(1)*dt > scalarCanopyWater)then
   message=trim(message)//'canopy sublimation exceeds total water storage on the canopy'
   err=-20; return  ! (negative error code forces a reduction in the length of the sub-step and another trial)
  endif
  ! remove mass of ice on the canopy
  mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1) = mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1) + &
                                                    mvar_data%var(iLookMVAR%scalarCanopySublimation)%dat(1)*dt
  ! if removed all ice, take the remaining sublimation from water
  if(mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1) < 0._dp)then
   mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1) = mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1) + mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1)
   mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1) = 0._dp
  endif
 endif  ! (if computing the vegetation flux)


 ! *****
 ! * compute sublimation from the top snow layer...
 ! ************************************************

 !print*, 'dt = ', dt
 !print*, 'mvar_data%var(iLookMVAR%scalarSnowSublimation)%dat(1) = ', mvar_data%var(iLookMVAR%scalarSnowSublimation)%dat(1)
 !print*, 'mvar_data%var(iLookMVAR%mLayerDepth)%dat(1)           = ', mvar_data%var(iLookMVAR%mLayerDepth)%dat(1)

 ! NOTE: this is done BEFORE densification
 if(nSnow > 0._dp)then ! snow layers exist
  ! volumetric sublimation (kg m-3)
  volSub = dt*mvar_data%var(iLookMVAR%scalarSnowSublimation)%dat(1)/mvar_data%var(iLookMVAR%mLayerDepth)%dat(1)
  !print*, 'volSub, iden_ice = ', volSub, iden_ice
  ! update volumetric fraction of ice (-)
  ! NOTE: fluxes are positive downward
  mLayerVolFracIceNew(1) = mLayerVolFracIceNew(1) + volSub/iden_ice
  ! check that we did not melt/sublimate all of the ice
  if(mLayerVolFracIceNew(1) < 0._dp)then
   message=trim(message)//'sublimated all of the ice from the top snow layer'
   err=-20; return ! (negative error code forces a reduction in the length of the sub-step and another trial)
  endif
 endif

 ! compute SWE
 !testSWE = sum( (mLayerVolFracLiqNew(1:nSnow)*iden_water + mLayerVolFracIceNew(1:nSnow)*iden_ice)*mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow) )
 !print*, 'after sublimation: SWE = ', testSWE

 ! *****
 ! * compute densification of the snowpack...
 ! ******************************************

 ! NOTE: this is done BEFORE the state update, to account for the change in ice content over the time step
 call snwDensify(&
                 ! intent(in): variables
                 dt,                                                     & ! intent(in): time step (s)
                 mvar_data%var(iLookMVAR%mLayerTemp)%dat(1:nSnow),       & ! intent(in): temperature of each layer (K)
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
                 mLayerVolFracLiqNew,                                    & ! intent(inout): volumetric fraction of liquid water at the end of the sub-step (-)
                 mLayerVolFracIceNew,                                    & ! intent(inout): volumetric fraction of ice at the end of the sub-step (-)
                 ! output: error control
                 err,cmessage)                                             ! intent(out): error control

 !testSWE = sum( (mLayerVolFracLiqNew(1:nSnow)*iden_water + mLayerVolFracIceNew(1:nSnow)*iden_ice)*mvar_data%var(iLookMVAR%mLayerDepth)%dat(1:nSnow) )
 !print*, 'after densification: SWE = ', testSWE

 ! *****
 ! water balance check...
 ! **********************
 call wBal_check(&
                        
                 ! input: model control
                 dt,                                                                 & ! intent(in): time step (s)
                 mpar_data%var(iLookPARAM%wimplicit),                                & ! intent(in): weight assigned to start-of-step fluxes (-)
                 mvar_data%var(iLookMVAR%mLayerDepth)%dat,                           & ! intent(in): depth of each layer (m)
                 model_decisions(iLookDECISIONS%groundwatr)%iDecision,               & ! intent(in): choice of groundwater parameterization
                 model_decisions(iLookDECISIONS%spatial_gw)%iDecision,               & ! intent(in): choice of method for the spatial representation of groundwater

                 ! input: soil parameters
                 mpar_data%var(iLookPARAM%theta_sat),                                & ! intent(in): soil porosity (-)
                 mpar_data%var(iLookPARAM%specificStorage),                          & ! intent(in): specific storage (m-1)

                 ! input: state variables at the start of the sub-step
                 mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat,                      & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
                 mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat,                      & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)
                 mvar_data%var(iLookMVAR%mLayerMatricHead)%dat,                      & ! intent(in): matric head at the start of the sub-step (m)
                 mvar_data%var(iLookMVAR%scalarAquiferStorage)%dat(1),               & ! intent(in): aquifer storage at the start of the sub-step (m)

                 ! input: state variables after iteration                
                 mLayerVolFracIceNew,                                                & ! intent(in): volumetric fraction of ice at the end of the sub-step (-)
                 mLayerVolFracLiqNew,                                                & ! intent(in): volumetric fraction of liquid water at the end of the sub-step (-)
                 mLayerMatricHeadNew,                                                & ! intent(in): matric head at the end of the sub-step (-)
                 scalarAquiferStorageNew,                                            & ! intent(in): aquifer storage at the end of the sub-step (m)

                 ! input: model fluxes at the *START* of the sub-step
                 mvar_data%var(iLookMVAR%iLayerInitLiqFluxSoil)%dat,                 & ! intent(in): liquid water flux at the interface of each layer at the start of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%iLayerInitFluxReversal)%dat,                & ! intent(in): flow reversal flux at the interface of each layer at the start of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%mLayerInitBaseflow)%dat,                    & ! intent(in): baseflow from each layer at the start of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%mLayerInitTranspire)%dat,                   & ! intent(in): transpiration from each layer at the start of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%scalarInitAquiferRecharge)%dat(1),          & ! intent(in): recharge to the aquifer at the start of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%scalarInitAquiferBaseflow)%dat(1),          & ! intent(in): baseflow from the aquifer at the start of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%scalarInitAquiferTranspire)%dat(1),         & ! intent(in): transpiration from the aquifer at the start of the sub-step (m s-1)

                 ! input: model fluxes at the *END* of the sub-step
                 mvar_data%var(iLookMVAR%iLayerLiqFluxSoil)%dat,                     & ! intent(in): liquid water flux at the interface of each layer at the end of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%iLayerFluxReversal)%dat,                    & ! intent(in): flow reversal flux at the interface of each layer at the end of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%mLayerBaseflow)%dat,                        & ! intent(in): baseflow from each layer at the end of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%mLayerTranspire)%dat,                       & ! intent(in): transpiration from each layer at the end of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%scalarAquiferRecharge)%dat(1),              & ! intent(in): recharge to the aquifer at the end of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%scalarAquiferBaseflow)%dat(1),              & ! intent(in): baseflow from the aquifer at the end of the sub-step (m s-1)
                 mvar_data%var(iLookMVAR%scalarAquiferTranspire)%dat(1),             & ! intent(in): transpiration from the aquifer at the end of the sub-step (m s-1)

                 ! output: water balance of the soil zone
                 mvar_data%var(iLookMVAR%scalarSoilInflux)%dat(1),                   & ! intent(out): influx at the top of the soil zone (m s-1)
                 mvar_data%var(iLookMVAR%scalarSoilBaseflow)%dat(1),                 & ! intent(out): total baseflow from the soil zone (m s-1)
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
 mvar_data%var(iLookMVAR%scalarCanairTemp)%dat(1)  = scalarCanairTempNew      ! temperature of the canopy air space (K)
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
 ! private subroutine: wrapper subroutine for the picard solver itself
 ! *********************************************************************************************************
 subroutine picardSolv_muster(&

                              ! input: input variables from picardSolv
                              dt,                                & ! intent(in): time step (s)
                              maxiter,                           & ! intent(in): maximum number of iterations
                              firstSubstep,                      & ! intent(in): flag to denote first sub-step
                              computeVegFlux,                    & ! intent(in): flag to denote if computing energy flux over vegetation

                              ! input: upper boundary conditions
                              ixBcUpperSoilHydrology,            & ! intent(in): choice of upper boundary condition for soil hydrology
                              upperBoundHead,                    & ! intent(in): upper boundary condition for matric head (m)
                              upperBoundTheta,                   & ! intent(in): upper boundary condition for volumetric liquid water content (-)

                              ! input: coordinate variables
                              layerType,                         & ! intent(in): layer type (ix_soil or ix_snow)
                              mLayerDepth,                       & ! intent(in): depth of each layer (m)
                              mLayerHeight,                      & ! intent(in): height at the mid-point of each layer (m)
                              iLayerHeight,                      & ! intent(in): height at the interface of each layer (m)

                              ! input: model state variables
                              scalarCanairTemp,                  & ! intent(in): temperature of the canopy air space at the start of the sub-step (K)
                              scalarCanopyTemp,                  & ! intent(in): temperature of the vegetation canopy at the start of the sub-step (K)
                              scalarCanopyIce,                   & ! intent(in): mass of ice on the vegetation canopy at the start of the sub-step (kg m-2)
                              scalarCanopyLiq,                   & ! intent(in): mass of liquid water on the vegetation canopy at the start of the sub-step (kg m-2)
                              mLayerTemp,                        & ! intent(in): temperature of each snow/soil layer at the start of the sub-step (K)
                              mLayerVolFracIce,                  & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
                              mLayerVolFracLiq,                  & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)
                              mLayerMatricHead,                  & ! intent(in): matric head at the start of the sub-step (-)
                              scalarAquiferStorage,              & ! intent(in): aquifer storage at the start of the sub-step (m)

                              ! input: forcing and diagnostic variables
                              airtemp,                           & ! intent(in): air temperature (K)
                              scalarBulkVolHeatCapVeg,           & ! intent(in): volumetric heat capacity of vegetation (J m-3 K-1)
                              mLayerVolHtCapBulk,                & ! intent(in): volumetric heat capacity of each layer in the snow-soil vector (J m-3 K-1)

                              ! input: snow parameters
                              snowfrz_scale,                     & ! intent(in): scaling parameter for the snow freezing curve (K-1)

                              ! input: van Genutchen soil parameters
                              vGn_alpha,                         & ! intent(in): van Genutchen "alpha" parameter (m-1)
                              vGn_n,                             & ! intent(in): van Genutchen "n" parameter (-)
                              VGn_m,                             & ! intent(in): van Genutchen "m" parameter (-)
                              theta_sat,                         & ! intent(in): soil porosity (-)
                              theta_res,                         & ! intent(in): soil residual volumetric water content (-)

                              ! input: algorithmic control parameters
                              relConvTol_liquid,                 & ! intent(in): relative convergence tolerance for vol frac liq water (-)
                              absConvTol_liquid,                 & ! intent(in): absolute convergence tolerance for vol frac liq water (-)
                              relConvTol_matric,                 & ! intent(in): relative convergence tolerance for matric head (-)
                              absConvTol_matric,                 & ! intent(in): absolute convergence tolerance for matric head (m)
                              relConvTol_energy,                 & ! intent(in): relative convergence tolerance for energy (-)
                              absConvTol_energy,                 & ! intent(in): absolute convergence tolerance for energy (J m-3)
                              relConvTol_aquifr,                 & ! intent(in): relative convergence tolerance for aquifer storage (-)
                              absConvTol_aquifr,                 & ! intent(in): absolute convergence tolerance for aquifer storage (m)
                               
                              ! output: diagnostic variables
                              scalarVP_CanopyAir,                & ! intent(out): trial vapor pressure of the canopy air space (Pa)
                              scalarCanopyStabilityCorrection,   & ! intent(out): stability correction for the canopy (-)
                              scalarGroundStabilityCorrection,   & ! intent(out): stability correction for the ground surface (-)
                              mLayerTcrit,                       & ! intent(out): critical soil temperature where liquid water begins to freeze (K)
                              scalarCanopyMeltFreeze,            & ! intent(out): melt/freeze of water stored in the canopy (kg m-2 s-1)
                              scalarCanopyTranspiration,         & ! intent(out): canopy transpiration (kg m-2 s-1)
                              scalarCanopyEvaporation,           & ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
                              scalarGroundEvaporation,           & ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)

                              ! output: model state variables at the end of the step
                              ! NOTE: use intent(out) instead of intent(inout) to protect start-of-step variables
                              scalarCanairTempNew,               & ! intent(out): temperature of the canopy air space at the end of the sub-step (K)
                              scalarCanopyTempNew,               & ! intent(out): temperature of the vegetation canopy at the end of the sub-step (K)
                              scalarCanopyIceNew,                & ! intent(out): mass of ice on the vegetation canopy at the end of the sub-step (kg m-2)
                              scalarCanopyLiqNew,                & ! intent(out): mass of liquid water on the vegetation canopy at the end of the sub-step (kg m-2)
                              mLayerTempNew,                     & ! intent(out): temperature of each snow/soil layer at the end of the sub-step (K)
                              mLayerVolFracIceNew,               & ! intent(out): volumetric fraction of ice at the end of the sub-step (-)
                              mLayerVolFracLiqNew,               & ! intent(out): volumetric fraction of liquid water at the end of the sub-step (-)
                              mLayerMatricHeadNew,               & ! intent(out): matric head at the end of the sub-step (m)
                              scalarAquiferStorageNew,           & ! intent(out): aquifer storage at the end of the sub-step (m)

                              ! output: number of iterations
                              niter,                             & ! intent(out): number of iterations

                              ! output: error control
                              err,message)                         ! intent(out): error control

 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 ! provide access to subroutines
 USE heatTransf_module,only:heatTransf          ! compute change in temperature over the time step
 USE phseChange_module,only:phseChange          ! compute change in phase over the time step
 USE can_Hydrol_module,only:can_Hydrol          ! compute canopy water balance
 USE snowHydrol_module,only:snowHydrol          ! compute liquid water flow through the snowpack
 USE soilHydrol_module,only:soilHydrol          ! compute change in mass over the time step for the soil
 ! utility modules
 USE snow_utils_module,only:fracliquid          ! compute the fraction of liquid water at a given temperature (snow)
 USE soil_utils_module,only:crit_soilT          ! compute the critical temperature above which all water is unfrozen
 USE soil_utils_module,only:volFracLiq          ! compute volumetric fraction of liquid water
 USE soil_utils_module,only:matricHead          ! compute matric head (m)
 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 implicit none
 ! input
 real(dp),intent(in)            :: dt                          ! time step (seconds)
 integer(i4b),intent(in)        :: maxiter                     ! maximum number of iterations
 logical(lgt),intent(in)        :: firstSubStep                ! flag to indicate if we are processing the first sub-step
 logical(lgt),intent(in)        :: computeVegFlux              ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 ! input: upper boundary conditions
 integer(i4b),intent(in)        :: ixBcUpperSoilHydrology      ! choice of upper boundary condition for soil hydrology
 real(dp),intent(in)            :: upperBoundHead              ! upper boundary condition for matric head (m)
 real(dp),intent(in)            :: upperBoundTheta             ! upper boundary condition for volumetric liquid water content (-)
 ! input: coordinate variables 
 integer(i4b),intent(in)        :: layerType(:)                ! type of the layer (ix_soil or ix_snow)
 real(dp),intent(in)            :: mLayerDepth(:)              ! depth of each layer (m)
 real(dp),intent(in)            :: mLayerHeight(:)             ! height at the mid-point of each layer (m)
 real(dp),intent(in)            :: iLayerHeight(0:)            ! height at the interface of each layer (m)
 ! input: model state variables
 real(dp),intent(in)            :: scalarCanairTemp            ! temperature of the canopy air space (K)
 real(dp),intent(in)            :: scalarCanopyTemp            ! temperature of the vegetation canopy (K)
 real(dp),intent(in)            :: scalarCanopyIce             ! mass of ice on the vegetation canopy (kg m-2)
 real(dp),intent(in)            :: scalarCanopyLiq             ! mass of liquid water on the vegetation canopy (kg m-2)
 real(dp),intent(in)            :: mLayerTemp(:)               ! temperature of each layer (K)
 real(dp),intent(in)            :: mLayerVolFracIce(:)         ! volumetric fraction of ice in each layer (-)
 real(dp),intent(in)            :: mLayerVolFracLiq(:)         ! volumetric fraction of liquid water in each layer (-)
 real(dp),intent(in)            :: mLayerMatricHead(:)         ! matric head in each layer (-)
 real(dp),intent(in)            :: scalarAquiferStorage        ! aquifer storage at the start of the sub-step (m)
 ! input: forcing and diagnostic variables
 real(dp),intent(in)            :: airtemp                     ! air temperature (K)
 real(dp),intent(in)            :: scalarBulkVolHeatCapVeg     ! bulk volumetric heat capacity of vegetation (J m-3 K-1)
 real(dp),intent(in)            :: mLayerVolHtCapBulk(:)       ! volumetric heat capacity in each layer (J m-3 K-1) 
 ! input: snow parameters
 real(dp),intent(in)            :: snowfrz_scale               ! scaling parameter for the snow freezing curve (K-1)
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
 real(dp),intent(out)           :: scalarVP_CanopyAir          ! trial vapor pressure of the canopy air space (Pa)
 real(dp),intent(out)           :: scalarCanopyStabilityCorrection ! stability correction for the canopy (-)
 real(dp),intent(out)           :: scalarGroundStabilityCorrection ! stability correction for the ground surface (-)
 real(dp),intent(out)           :: mLayerTcrit(:)              ! critical soil temperature where liquid water begins to freeze (K)
 real(dp),intent(out)           :: scalarCanopyMeltFreeze      ! melt/freeze of water stored in the canopy (kg m-2 s-1)
 real(dp),intent(out)           :: scalarCanopyTranspiration   ! canopy transpiration (kg m-2 s-1)
 real(dp),intent(out)           :: scalarCanopyEvaporation     ! canopy evaporation/condensation (kg m-2 s-1)
 real(dp),intent(out)           :: scalarGroundEvaporation     ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
 ! output: model state variables at the end of the step
 ! NOTE: use intent(out) instead of intent(inout) to protect start-of-step variables
 real(dp),intent(out)           :: scalarCanairTempNew         ! temperature of the canopy air space at the end of the sub-step (K)
 real(dp),intent(out)           :: scalarCanopyTempNew         ! temperature of the vegetation canopy at the end of the sub-step (K)
 real(dp),intent(out)           :: scalarCanopyIceNew          ! mass of ice on the vegetation canopy at the end of the sub-step (kg m-2)
 real(dp),intent(out)           :: scalarCanopyLiqNew          ! mass of liquid water on the vegetation canopy at the end of the sub-step (kg m-2)
 real(dp),intent(out)           :: mLayerTempNew(:)            ! temperature of each snow/soil layer at the end of the sub-step (K)
 real(dp),intent(out)           :: mLayerVolFracIceNew(:)      ! volumetric fraction of ice at the end of the sub-step (-)
 real(dp),intent(out)           :: mLayerVolFracLiqNew(:)      ! volumetric fraction of liquid water at the end of the sub-step (-)
 real(dp),intent(out)           :: mLayerMatricHeadNew(:)      ! matric head at the current iteration end of the sub-step (m)
 real(dp),intent(out)           :: scalarAquiferStorageNew     ! aquifer storage at the end of the sub-step (m)
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
 real(dp)                       :: scalarSurfaceInfiltration   ! infiltration rate (m s-1)
 ! define derivatives in canopy air space variables
 real(dp)                       :: dTempCanopyAir_dTCanopy     ! derivative in the temperature of the canopy air space w.r.t. temperature of the canopy
 real(dp)                       :: dTempCanopyAir_dTGround     ! derivative in the temperature of the canopy air space w.r.t. temperature of the ground
 real(dp)                       :: dVPCanopyAir_dTCanopy       ! derivative in the vapor pressure of the canopy air space w.r.t. temperature of the canopy 
 real(dp)                       :: dVPCanopyAir_dTGround       ! derivative in the vapor pressure of the canopy air space w.r.t. temperature of the ground
 ! define state variables for the vegetation canopy
 real(dp)                       :: scalarCanairTempIter        ! trial value of temperature of the canopy air space (K)
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
 real(dp)                       :: scalarCanairTempIncr        ! iteration increment for temperature of the canopy air space (K)
 real(dp)                       :: scalarCanopyTempIncr        ! iteration increment for temperature of the vegetation canopy (K)
 real(dp),dimension(nLayers)    :: mLayerTempIncr              ! iteration increment for temperature of the snow-soil vector (K)
 real(dp),dimension(nLayers)    :: mLayerLiqIncr               ! iteration increment for volumetric liquid water content in the snow-soil vector (-)
 real(dp),dimension(nSoil)      :: mLayerMatIncr               ! iteration increment for matric head in soil layers (m)
 real(dp)                       :: canopyTempIncrOld           ! iteration increment for temperature of the vegetation canopy in the previous iteration (K)
 real(dp),dimension(nLayers)    :: mLayerTempIncrOld           ! iteration increment for temperature of the snow-soil vector in the previous iteration (K)
 real(dp),dimension(nSoil)      :: mLayerMatricIncrOld         ! iteration increment for matric head of soil layers in the previous iteration (m)
 real(dp),dimension(nSoil)      :: mLayerLiquidIncrOld         ! iteration increment for volumetric liquid water content of soil layers in the previous iteration (-)
 real(dp),dimension(nSoil)      :: mLayerResidual              ! residual in the soil hydrology equations (-)
 ! define solution constraint variables
 integer(i4b),dimension(1)      :: iLayerViolated              ! index of layer that was violated
 real(dp),dimension(1)          :: constraintViolation         ! maximum constraint violation
 real(dp),dimension(nLevels)    :: mLayerMatricHeadDiff        ! iteration increment for matric head (m)
 real(dp)                       :: scaleIncrement              ! scaling factor for the iteration increment
 real(dp)                       :: maxMatric                   ! maximum possible value of matric head (m)
 real(dp)                       :: allowablePressureViolation=0.0_dp ! allowable pressure violation before time step reduction (m)
 real(dp),dimension(nLevels)    :: soilVolFracLiqNew           ! volumetric fraction of liquid water in the soil 
 real(dp),dimension(nLevels)    :: soilVolFracIceNew           ! volumetric fraction of ice in the soil 
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
 real(dp),dimension(1)          :: hydrol_max                  ! maximum residual for soil hydrology for a given iteration (-)
 real(dp),dimension(1)          :: liquid_max                  ! maximum absolute change in volumetric liquid water content for a given iteration (-)
 real(dp),dimension(1)          :: matric_max                  ! maximum absolute change in matric head for a given iteration (m)
 real(dp),dimension(1)          :: energy_max                  ! maximum absolute change in energy for a given iteration (J m-3)
 real(dp)                       :: aquifr_max                  ! absolute change in aquifer storage for a given iteration (m)
 ! define local variables for phase change
 real(dp)                       :: fLiq                        ! fraction of liquid water on the vegetation canopy (-)
 real(dp)                       :: tWat                        ! total water on the vegetation canopy (kg m-2) 
 real(dp)                       :: testSWE                     ! SWE (kg m-2)
 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="picardSolv_muster/"

 ! initialize number of iterations
 niter=0

 ! *****
 ! * initialize...
 ! ***************

 ! initialize stability correction (no stability correction)
 scalarCanopyStabilityCorrection = 1._dp          ! stability correction for the canopy (-)
 scalarGroundStabilityCorrection = 1._dp          ! stability correction for the ground surface (-)

 ! initialize temperature of the vegetation canopy and the canopy air space
 scalarCanairTempIter = scalarCanairTemp
 scalarCanopyTempIter = scalarCanopyTemp

 ! initialize canopy water
 scalarCanopyIceIter = scalarCanopyIce            ! mass of ice on the vegetation canopy (kg m-2)
 scalarCanopyLiqIter = scalarCanopyLiq            ! mass of liquid water on the vegetation canopy (kg m-2)

 ! initialize layer temperatures
 mLayerTempIter = mLayerTemp                      ! temperature (K)

 ! check layer temperature
 if(any(mLayerTemp < Tfreeze-80._dp))then
  message=trim(message)//'layer temperature is very cold'
  err=20; return
 endif

 ! initialize volumetric liquid and ice content
 mLayerVolFracIceIter = mLayerVolFracIce          ! volumetric ice content (-)
 mLayerVolFracLiqIter = mLayerVolFracLiq          ! volumetric liquid water content (-)

 ! initialize matric head
 mLayerMatricHeadIter = mLayerMatricHead          ! matric head (m)

 ! initialize the aquifer storage
 scalarAquiferStorageIter = scalarAquiferStorage  ! aquifer storage (m)

 ! calculate the critical soil temperature above which all water is unfrozen (K)
 do iLayer=nSnow+1,nLayers
  theta = mLayerVolFracIce(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiq(iLayer)
  mLayerTcrit(iLayer-nSnow) = crit_soilT(theta,theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m)
 end do

 ! initialize the melt/freeze vectors
 mLayerMeltFreeze  = 0._dp
 mLayerInfilFreeze = 0._dp

 !if(any(mLayerVolFracIce*iden_ice > 700._dp)) printflag=.true.

 ! calculate phase change of water on the vegetation canopy
 ! NOTE: this is necessary to ensure consistency at the start of the iterations (e.g., snowfall on warm canopy)
 if(computeVegFlux)then
  ! compute the fraction of liquid water
  fLiq = fracliquid(scalarCanopyTemp,snowfrz_scale)  ! fraction of liquid water (-)
  tWat = scalarCanopyLiq + scalarCanopyIce           ! total water (kg m-2)
  scalarCanopyLiqIter = fLiq*tWat                    ! mass of liquid water on the canopy (kg m-2)
  scalarCanopyIceIter = (1._dp - fLiq)*tWat          ! mass of ice on the canopy (kg m-2)
 endif

 ! check SWE
 !if(nSnow > 0)then
 ! testSWE = sum( (iden_ice*mLayerVolFracIceIter(1:nSnow) + iden_water*mLayerVolFracLiqIter(1:nSnow))*mLayerDepth(1:nSnow) )
 ! print*, 'start picardSolv_muster: SWE = ', testSWE
 !endif

 !print*, '*** dt, maxiter = ', dt, maxiter
 !print*, 'computeVegFlux = ', computeVegFlux

 ! **********************************************************************************************************************
 ! **********************************************************************************************************************
 ! **********************************************************************************************************************
 ! ***** iterate ********************************************************************************************************
 ! **********************************************************************************************************************
 ! **********************************************************************************************************************
 ! **********************************************************************************************************************
 do iter=1,maxiter

  print*, '***************************************************************'
  print*, '***** iter = ', iter, '*****'
  print*, '***************************************************************'
  

  ! *****
  ! * thermodynamics...
  ! *******************

  ! increment number of iterations
  niter=niter+1

  ! check SWE
  !if(nSnow > 0)then
  ! testSWE = sum( (iden_ice*mLayerVolFracIceIter(1:nSnow) + iden_water*mLayerVolFracLiqIter(1:nSnow))*mLayerDepth(1:nSnow) )
  ! write(*,'(a,1x,i4,1x,f20.10)') 'iterating picardSolv_muster, start of thermo: SWE = ', iter, testSWE
  !endif

  ! compute the temperature and ice content at the next iteration
  call heatTransf(&
                  ! input
                  dt,&                              ! intent(in): time step (seconds)
                  iter,&                            ! intent(in): current iteration count
                  firstSubstep,                   & ! intent(in): flag to indicate if we are processing the first sub-step
                  computeVegFlux,                 & ! intent(in): flag to indicate if we computing fluxes ovser vegetation (.false. means veg is buried with snow)
                  scalarCanairTempIter,           & ! intent(in): trial temperature of the canopy air space at the current iteration (K)
                  scalarCanopyTempIter,           & ! intent(in): trial temperature of the vegetation canopy at the current iteration (K)
                  scalarCanopyIceIter,            & ! intent(in): trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
                  scalarCanopyLiqIter,            & ! intent(in): trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                  mLayerTempIter,                 & ! intent(in): trial temperature of each snow/soil layer at the current iteration (K)
                  mLayerVolFracIceIter,           & ! intent(in): trial volumetric fraction of ice in each snow/soil layer at the current iteration (-)
                  mLayerVolFracLiqIter,           & ! intent(in): trial volumetric fraction of liquid water in each snow/soil layer at the current iteration (-)
                  mLayerMatricHeadIter,           & ! intent(in): trial matric head of each snow/soil layer at the current iteration (m)
                  canopyTempIncrOld,              & ! intent(in): previous iteration increment in canopy temperature (K)
                  mLayerTempIncrOld,              & ! intent(in): previous iteration increment in temperature of the snow-soil vector (K)

                  ! input/output variables from heatTransf subroutine: canopy air space variables
                  scalarVP_CanopyAir,             & ! intent(inout): trial vapor pressure of the canopy air space (Pa)
                  scalarCanopyStabilityCorrection,& ! intent(inout): stability correction for the canopy (-)
                  scalarGroundStabilityCorrection,& ! intent(inout): stability correction for the ground surface (-)

                  ! output
                  scalarCanopyTranspiration,      & ! intent(out): canopy transpiration (kg m-2 s-1)
                  scalarCanopyEvaporation,        & ! intent(out): canopy evaporation/condensation (kg m-2 s-1)
                  scalarGroundEvaporation,        & ! intent(out): ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
                  scalarCanairTempIncr,           & ! intent(out): iteration increment for temperature of the canopy air space (K)
                  scalarCanopyTempIncr,           & ! intent(out): iteration increment for temperature of the vegetation canopy (K)
                  mLayerTempIncr,                 & ! intent(out): iteration increment for temperature of the snow-soil system (K)
                  scalarCanairTempNew,            & ! intent(out): new temperature of the canopy air space (K)
                  scalarCanopyTempNew,            & ! intent(out): new temperature of the vegetation canopy (K)
                  scalarCanopyIceNew,             & ! intent(out): mass of ice on the canopy (kg m-2) 
                  scalarCanopyLiqNew,             & ! intent(out): mass of liquid water on the canopy (kg m-2)
                  mLayerTempNew,                  & ! intent(out): new temperature each snow/soil layer (K)
                  mLayerMatricHeadNew,            & ! intent(out): after phase change: new matric head (m)
                  mLayerVolFracLiqNew,            & ! intent(out): after phase change: new volumetric fraction of liquid water (-)
                  mLayerVolFracIceNew,            & ! intent(out): after phase change: new volumetric fraction of ice (-)
                   
                  err,cmessage)                     ! intent(out): error control
  ! check errors
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

  ! test
  !write(*,'(a,1x,i4,1x,10(f12.8,1x))') 'in picardSolv: iter, mLayerVolFracLiqIter(1:nSnow+10) = ', iter, mLayerVolFracLiqIter(1:nSnow+5)
  !write(*,'(a,1x,i4,1x,10(f12.8,1x))') 'in picardSolv: iter, mLayerVolFracLiqNew(1:nSnow+10)  = ', iter, mLayerVolFracLiqNew(1:nSnow+5)
  write(*,'(a,1x,i4,1x,10(f15.11,1x))') 'in picardSolv: iter, mLayerVolFracIceIter(1:nSnow+5) = ', iter, mLayerVolFracIceIter(1:nSnow+5)
  write(*,'(a,1x,i4,1x,10(f15.11,1x))') 'in picardSolv: iter, mLayerVolFracIceNew(1:nSnow+5)  = ', iter, mLayerVolFracIceNew(1:nSnow+5)
  !write(*,'(a,1x,i4,1x,10(f12.8,1x))') 'in picardSolv: iter, mLayerVolFracIceIncr(1:nSnow) = ', iter, mLayerVolFracIceNew(1:nSnow) - mLayerVolFracIceIter(1:nSnow)
  !write(*,'(a,1x,i4,1x,10(f15.8,1x))') 'in picardSolv: iter, mLayerTempIncr(1:nSnow+5)       = ', iter, mLayerTempIncr(1:nSnow+5)
  write(*,'(a,1x,i4,1x,10(f15.8,1x))') 'in picardSolv: iter, mLayerTempIter(1:nSnow+5)       = ', iter, mLayerTempIter(1:nSnow+5)
  write(*,'(a,1x,i4,1x,10(f15.8,1x))') 'in picardSolv: iter, mLayerTempNew(1:nSnow+5)        = ', iter, mLayerTempNew(1:nSnow+5)

  write(*,'(a)')  'in picardSolv: iter, airtemp, scalarCanairTempIter, scalarCanopyTempIter, mLayerTempIter(1), scalarVP_CanopyAir, scalarCanairTempIncr, scalarCanopyTempIncr, mLayerTempIncr(1) = '
  write(*,'(i4,1x,10(f12.5,1x))') iter, airtemp, scalarCanairTempIter, scalarCanopyTempIter, mLayerTempIter(1), scalarVP_CanopyAir, scalarCanairTempIncr, scalarCanopyTempIncr, mLayerTempIncr(1)

  !write(*,'(a,1x,i4,1x,10(f12.5,1x))') 'in picardSolv: iter, dt, scalarCanopyIceIter, scalarCanopyLiqIter, scalarCanopyIceNew, scalarCanopyLiqNew = ', &
  !                                                     iter, dt, scalarCanopyIceIter, scalarCanopyLiqIter, scalarCanopyIceNew, scalarCanopyLiqNew

  ! compute melt/freeze in each layer (kg m-3 s-1) -- melt is negative
  mLayerMeltFreeze = mLayerMeltFreeze + iden_ice*(mLayerVolFracIceNew - mLayerVolFracIceIter)/dt

  ! compute the temperature and phase components of the energy increment (J m-3)
  canopyTempComponent = scalarBulkVolHeatCapVeg*scalarCanopyTempIncr
  mLayerTempComponent = mLayerVolHtCapBulk*mLayerTempIncr
  mLayerPhseComponent = iden_ice*LH_fus*(mLayerVolFracIceNew - mLayerVolFracIceIter)

  ! save the temperature increment
  canopyTempIncrOld   = scalarCanopyTempIncr
  mLayerTempIncrOld   = mLayerTempIncr
  !print*, 'in PicardSolv: canopyTempIncrOld = ', canopyTempIncrOld


  ! =================================================================================================================================================
  ! =================================================================================================================================================


  ! *****
  ! * hydrology...
  ! **************

  ! check SWE
  !if(nSnow > 0)then
  ! testSWE = sum( (iden_ice*mLayerVolFracIceNew(1:nSnow) + iden_water*mLayerVolFracLiqIter(1:nSnow))*mLayerDepth(1:nSnow) )
  ! write(*,'(a,1x,i4,1x,f20.10)') 'iterating picardSolv_muster, start hydrology: SWE = ', iter, testSWE
  !endif

  ! get initial estimate of liquid water and ice on the vegetation canopy
  ! NOTE: the initial estimate used here is after phase change in heatTrans
  !scalarCanopyIceIter = scalarCanopyIceNew  ! (can pass IceNew to can_Hydrol because it is not changed)
  scalarCanopyLiqIter = scalarCanopyLiqNew

  ! get initial estimate of state variables for the snow-soil vector
  ! NOTE 1: the initial estimate used here is after phase change in heatTrans
  ! NOTE 2: ice content is not changed in soilHydrol routines so can be passed as new
  mLayerMatricHeadIter = mLayerMatricHeadNew
  mLayerVolFracLiqIter = mLayerVolFracLiqNew

  ! compute the canopy water balance
  call can_Hydrol(&
                  ! input
                  dt,                       & ! time step (seconds)
                  iter,                     & ! iteration index
                  computeVegFlux,           & ! flag to denote if computing energy flux over vegetation
                  scalarCanopyIceNew,       & ! *** from heatTrans (not changed here) *** trial mass of ice on the vegetation canopy at the current iteration (kg m-2)
                  scalarCanopyLiqIter,      & ! trial mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                  scalarCanopyEvaporation,  & ! canopy evaporation/condensation (kg m-2 s-1)
                  ! output
                  scalarCanopyLiqNew,       & ! updated mass of liquid water on the vegetation canopy at the current iteration (kg m-2)
                  err,cmessage)               ! error control
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  !print*, 'after can_Hydrol: scalarCanopyLiqIter = ', scalarCanopyLiqIter
  !print*, 'after can_Hydrol: scalarCanopyLiqNew = ', scalarCanopyLiqNew


  ! compute the volumetric liquid water content at the next iteration (note: only use snow vectors)
  ! NOTE: ice not modified in the snow hydrology routines, so can stay as "New"
  if(nSnow > 0)then
   call snowHydrol(&
                   ! input
                   dt,                               & ! time step (seconds)
                   iter,                             & ! iteration index
                   mLayerVolFracIceNew(1:nSnow),     & ! ***** from heatTrans (not changed here) ***** volumetric fraction of ice at the current iteration (-)
                   mLayerVolFracLiqIter(1:nSnow),    & ! volumetric fraction of liquid water at the current iteration (-)
                   ! output
                   mLayerVolFracLiqNew(1:nSnow),     & ! volumetric fraction of liquid water at the next iteration (-)
                   err,cmessage)
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  !print*, 'after snowHydrol: mLayerVolFracLiqIter(1:nSnow) = ', mLayerVolFracLiqIter(1:nSnow)
  !print*, 'after snowHydrol: mLayerVolFracLiqNew(1:nSnow)  = ', mLayerVolFracLiqNew(1:nSnow)

  !write(*,'(a,10(f20.10,1x))')  'before soilHydrol: mLayerVolFracLiqNew(nSnow+1), mLayerVolFracIceNew(nSnow+1), mLayerVolFracLiqNew(nSnow+1) + mLayerVolFracIceNew(nSnow+1) = ', &
  !                                                  mLayerVolFracLiqNew(nSnow+1), mLayerVolFracIceNew(nSnow+1), mLayerVolFracLiqNew(nSnow+1) + mLayerVolFracIceNew(nSnow+1)

  ! compute the matric head at the next iteration (note liquid water and ice vectors are defined for all layers)
  ! NOTE: ice not modified in the soil hydrology routines, so can stay as "New"
  call soilHydrol(&
                  ! input
                  dt,                                         & ! time step (seconds)
                  iter,                                       & ! current iteration count
                  ixBcUpperSoilHydrology,                     & ! choice of upper boundary condition for soil hydrology
                  upperBoundHead,                             & ! upper boundary condition for matric head (m)
                  upperBoundTheta,                            & ! upper boundary condition for volumetric liquid water content (-)
                  scalarCanopyTranspiration,                  & ! canopy transpiration (kg m-2 s-1)
                  scalarGroundEvaporation,                    & ! ground evaporation/condensation -- below canopy or non-vegetated (kg m-2 s-1)
                  mLayerMatricHeadIter(1:nLevels),            & ! matric head in each layer at the current iteration (m)
                  mLayerVolFracLiqIter(nSnow+1:nSnow+nLevels),& ! volumetric fraction of liquid water at the current iteration (-)
                  mLayerVolFracIceNew(nSnow+1:nSnow+nLevels), & ! *** from heatTrans (not changed here) *** volumetric fraction of ice at the current iteration (-)
                  scalarAquiferStorageIter,                   & ! aquifer storage (m)
                  mLayerMatricIncrOld(1:nLevels),             & ! iteration increment for matric head of soil layers in the previous iteration (m)
                  mLayerLiquidIncrOld(1:nLevels),             & ! iteration increment for volumetric liquid water content of soil layers in the previous iteration (-) 
                  ! output
                  mLayerResidual(1:nLevels),                  & ! residual vector (-)
                  mLayerMatricHeadNew(1:nLevels),             & ! matric head in each layer at the next iteration (m)
                  mLayerVolFracLiqNew(nSnow+1:nSnow+nLevels), & ! volumetric fraction of liquid water at the next iteration (-)
                  scalarAquiferStorageNew,                    & ! aquifer storage (m)
                  scalarSurfaceInfiltration,                  & ! surface infiltration rate (m s-1)
                  err,cmessage)
  if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  ! copy across "saturated" layers
  if(nLevels < nSoil)then
   mLayerResidual(nLevels+1:nSoil)              = 0._dp
   mLayerMatricHeadNew(nLevels+1:nSoil)         = mLayerMatricHeadIter(nLevels+1:nSoil)
   mLayerVolFracLiqNew(nSnow+nLevels+1:nLayers) = mLayerVolFracLiqIter(nSnow+nLevels+1:nLayers)
   message=trim(message)//'not sure if nLevels < nSoil still works'
   err=20; return
  endif

  ! check SWE
  if(nSnow > 0)then
   testSWE = sum( (iden_ice*mLayerVolFracIceNew(1:nSnow) + iden_water*mLayerVolFracLiqNew(1:nSnow))*mLayerDepth(1:nSnow) )
   write(*,'(a,1x,i4,1x,f20.10)') 'iterating picardSolv_muster, after hydrology: SWE = ', iter, testSWE
  endif

  ! save iteration increments
  mLayerMatricIncrOld(1:nSoil) = mLayerMatricHeadNew(1:nSoil) - mLayerMatricHeadIter(1:nSoil)
  mLayerLiquidIncrOld(1:nSoil) = mLayerVolFracLiqNew(nSnow+1:nSnow+nSoil) - mLayerVolFracLiqIter(nSnow+1:nSnow+nSoil)

  ! compute the iteration increment for the matric head and volumetric fraction of liquid water
  mLayerMatIncr = mLayerMatricHeadNew - mLayerMatricHeadIter 
  mLayerLiqIncr = mLayerVolFracLiqNew - mLayerVolFracLiqIter 
  write(*,'(a,10(f20.10,1x))') 'in picardSolv: mLayerMatricHeadIter = ', mLayerMatricHeadIter(1:5)
  write(*,'(a,10(f20.10,1x))') 'in picardSolv: mLayerMatricHeadNew  = ', mLayerMatricHeadNew(1:5)
  !pause
  !print*, 'after soilHydrol'
  !print*, 'mLayerVolFracLiq     = ', mLayerVolFracLiq
  !print*, 'mLayerVolFracLiqIter = ', mLayerVolFracLiqIter
  !print*, 'mLayerVolFracLiqNew  = ', mLayerVolFracLiqNew
  write(*,'(a,10(e20.10,1x))') 'in picardSolv: mLayerMatIncr = ', mLayerMatIncr
  !print*, 'mLayerLiqIncr = ', mLayerLiqIncr
  !pause

  ! compute the iteration increment for aquifer storage
  scalarAqiIncr = scalarAquiferStorageNew - scalarAquiferStorageIter

  ! copy across volumetric liquid and ice content in the soil
  ! NOTE: this is done just to simplify the indices
  soilVolFracIceNew = mLayerVolFracIceNew(nSnow+1:nSnow+nLevels)  ! after phase change in heatTrans
  soilVolFracLiqNew = mLayerVolFracLiqNew(nSnow+1:nSnow+nLevels)  ! after liquid flux in soilHydrol

  ! *****
  ! impose solution constraints to deal with excessive storage -- simplified bi-section method
  ! NOTE: doing this after computing the iteration increment to avoid premature convergence
  constraintViolation = maxval(soilVolFracLiqNew + soilVolFracIceNew*(iden_ice/iden_water))
  if(constraintViolation(1) > theta_sat)then  ! at least one layer where volumetric (liquid + ice) content exceeds soil porosity
   ! print initial solution
   write(*,'(a,5(f20.10,1x))') 'mLayerMatricHeadNew(1:5) = ', mLayerMatricHeadNew(1:5)
   write(*,'(a,5(f20.10,1x))') 'mLayerVolFracIceNew(1:5) = ', soilVolFracIceNew(1:5)
   write(*,'(a,5(f20.10,1x))') 'mLayerVolFracLiqNew(1:5) = ', soilVolFracLiqNew(1:5)
   write(*,'(a,5(f20.10,1x))') 'mLayerVolFracWater(1:5)  = ', soilVolFracIceNew(1:5) + soilVolFracLiqNew(1:5)
   write(*,'(a,5(e20.10,1x))') 'mLayerMatIncr(1:5)       = ', mLayerMatIncr(1:5)
   ! identify the layer that was violated
   iLayerViolated       = maxloc(soilVolFracLiqNew + soilVolFracIceNew*(iden_ice/iden_water))
   ! identify the maximum possible matric head for the layer that was violated
   maxMatric            = matricHead(theta_sat - soilVolFracIceNew(iLayerViolated(1))*(iden_ice/iden_water),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   ! check
   print*, 'iLayerViolated = ', iLayerViolated
   print*, 'maxMatric      = ', maxMatric
   if(maxMatric < mLayerMatricHeadIter(iLayerViolated(1)) .or. mLayerMatIncr(iLayerViolated(1)) < -tiny(dt))then
    message=trim(message)//'expect maximum matric head to be greater than initial matric head, and matric head increment to be positive'
    err=-20; return  ! force time step reduction
   endif
   ! scale the iteration increment to the midpoint between previous iteration and maximum value of matric head
   scaleIncrement       = 0.95_dp*(maxMatric - mLayerMatricHeadIter(iLayerViolated(1))) / mLayerMatIncr(iLayerViolated(1))
   mLayerMatricHeadDiff = mLayerMatIncr*scaleIncrement
   ! print the layer that was violated
   print*, 'mLayerVolFracLiqIter(iLayerViolated(1)+nSnow) = ', mLayerVolFracLiqIter(iLayerViolated(1)+nSnow)
   print*, 'mLayerVolFracLiqNew(iLayerViolated(1)+nSnow)  = ', mLayerVolFracLiqNew(iLayerViolated(1)+nSnow)
   print*, 'theta_sat, soilVolFracIceNew(iLayerViolated(1)), soilVolFracLiqNew(iLayerViolated(1)) = ', &
            theta_sat, soilVolFracIceNew(iLayerViolated(1)), soilVolFracLiqNew(iLayerViolated(1))
   print*, 'mLayerMatricHeadIter(iLayerViolated(1)) = ', mLayerMatricHeadIter(iLayerViolated(1))
   print*, 'scaleIncrement = ', scaleIncrement
   ! update model states again
   do iLayer=1,nLevels
    mLayerMatricHeadNew(iLayer)   = mLayerMatricHeadIter(iLayer) + mLayerMatricHeadDiff(iLayer)
    soilVolFracLiqNew(iLayer)     = volFracLiq(mLayerMatricHeadNew(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
    if(iLayer==1) write(*,'(a)')    'iLayer, mLayerMatricHeadDiff(iLayer), mLayerMatricHeadNew(iLayer), mLayerVolFracLiqNew(iLayer), mLayerVolFracIceIter(iLayer), mLayerVolFracWater(iLayer) = '
    write(*,'(i4,1x,10(f20.10,1x))') iLayer, mLayerMatricHeadDiff(iLayer), mLayerMatricHeadNew(iLayer), mLayerVolFracLiqNew(iLayer), mLayerVolFracIceIter(iLayer), mLayerVolFracIceIter(iLayer) + mLayerVolFracLiqNew(iLayer)
   end do
   if(mLayerMatIncr(iLayerViolated(1)) < 0._dp)then; err=20; message=trim(message)//'expect that matric head increment is positive'; return; endif
   !pause 'constraint violation'
  endif  ! if there was a constraint violation

  ! *****
  ! impose solution constraints to deal with excessive pressure -- simplified bi-section method
  ! NOTE: doing this after computing the iteration increment to avoid premature convergence
  constraintViolation = maxval(mLayerMatricHeadNew - mLayerHeight(nSnow+1:nLevels))
  if(constraintViolation(1) > allowablePressureViolation)then  ! at least one layer where matric head is greater than soil depth
   print*, 'excessive pressure:'
   write(*,'(a,10(f20.10,1x))') 'mLayerMatricHeadIter = ', mLayerMatricHeadIter
   write(*,'(a,10(f20.10,1x))') 'mLayerMatricHeadNew  = ', mLayerMatricHeadNew
   write(*,'(a,10(f20.10,1x))') 'mLayerVolFracIceIter = ', soilVolFracIceNew
   write(*,'(a,10(f20.10,1x))') 'mLayerVolFracLiqNew  = ', soilVolFracLiqNew
   write(*,'(a,10(f20.10,1x))') 'mLayerVolFracWater   = ', soilVolFracIceNew + soilVolFracLiqNew
   write(*,'(a,10(f20.10,1x))') 'mLayerMatIncr        = ', mLayerMatIncr
   ! scale the iteration increment to the midpoint between previous iteration and layer depth for the layer where constraints were violated
   iLayerViolated       = maxloc(mLayerMatricHeadNew - mLayerHeight(nSnow+1:nLevels))
   scaleIncrement       = 0.5_dp*(mLayerHeight(nSnow+iLayerViolated(1)) - mLayerMatricHeadIter(iLayerViolated(1)))/mLayerMatIncr(iLayerViolated(1))
   mLayerMatricHeadDiff = mLayerMatIncr*scaleIncrement
   ! update model states again
   do iLayer=1,nLevels
    mLayerMatricHeadNew(iLayer) = mLayerMatricHeadIter(iLayer) + mLayerMatricHeadDiff(iLayer)
    soilVolFracLiqNew(iLayer)   = volFracLiq(mLayerMatricHeadNew(iLayer),vGn_alpha,theta_res,theta_sat,vGn_n,vGn_m)
   end do
  endif  ! if there was a constraint violation

  !write(*,'(a,10(f20.10,1x))')  'before phase change: mLayerVolFracLiqNew(nSnow+1), mLayerVolFracIceNew(nSnow+1), mLayerVolFracLiqNew(nSnow+1) + mLayerVolFracIceNew(nSnow+1) = ', &
  !                                                    mLayerVolFracLiqNew(nSnow+1), mLayerVolFracIceNew(nSnow+1), mLayerVolFracLiqNew(nSnow+1) + mLayerVolFracIceNew(nSnow+1)

  ! =================================================================================================================================================
  ! =================================================================================================================================================

  !print*, 'scalarCanopyLiqIter = ', scalarCanopyLiqIter
  !print*, 'scalarCanopyIceIter = ', scalarCanopyIceIter


  ! *****
  ! * phase change.....
  ! *******************

  ! update the canopy variables -- used in phase change below
  scalarCanopyLiqIter = scalarCanopyLiqNew

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

  !print*, 'computeVegFlux = ', computeVegFlux
  !print*, 'scalarCanopyLiqNew = ', scalarCanopyLiqNew
  !print*, 'scalarCanopyIceNew = ', scalarCanopyIceNew


  ! check
  if(scalarCanopyTempNew > Tfreeze .and. scalarCanopyIceNew > 0._dp)then
   message=trim(message)//'ice content above zero when temperature is above freezing'
   err=20; return
  endif

  ! calculate the critical soil temperature above which all water is unfrozen (K)
  !do iLayer=nSnow+1,nLayers
  ! theta = mLayerVolFracIceNew(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiqNew(iLayer)
  ! mLayerTcrit(iLayer-nSnow) = crit_soilT(theta,theta_res,theta_sat,vGn_alpha,vGn_n,vGn_m)
  !end do

  ! option to compute phase change associated with infiltrating liquid water
  if(freeze_infiltrate)then

   ! update the variables in the snow-soil vector -- used in phase change below
   ! NOTE: efficiency gains possible by using New as input and Iter as output
   mLayerMatricHeadIter = mLayerMatricHeadNew
   mLayerVolFracLiqIter = mLayerVolFracLiqNew
   mLayerVolFracIceIter = mLayerVolFracIceNew

   ! compute change in volumetric ice content associated with infiltrating liquid water
   call phsechange(&
                   ! input
                   mLayerTempNew,        & ! intent(in): new temperature vector (K)
                   mLayerMatricHeadIter, & ! intent(in): matric head at the current iteration (m)
                   mLayerVolFracLiqIter, & ! intent(in): volumetric fraction of liquid water at the current iteration (-)
                   mLayerVolFracIceIter, & ! intent(in): volumetric fraction of ice at the current iteration (-)
                   ! output
                   mLayerMatricHeadNew,  & ! intent(out): new matric head (m)
                   mLayerVolFracLiqNew,  & ! intent(out): new volumetric fraction of liquid water (-)
                   mLayerVolFracIceNew,  & ! intent(out): new volumetric fraction of ice (-)
                   err,cmessage)           ! intent(out): error control
   if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
  endif

  !write(*,'(a,10(f20.10,1x))')  'after phase change: mLayerVolFracLiqNew(nSnow+1), mLayerVolFracIceNew(nSnow+1), mLayerVolFracLiqNew(nSnow+1) + mLayerVolFracIceNew(nSnow+1) = ', &
  !                                                   mLayerVolFracLiqNew(nSnow+1), mLayerVolFracIceNew(nSnow+1), mLayerVolFracLiqNew(nSnow+1) + mLayerVolFracIceNew(nSnow+1)

  !write(*,'(a,1x,100(f20.10,1x))') 'after phseChange: iLayerHeight(0:nSnow+2) = ', iLayerHeight(0:nSnow+2)
  !write(*,'(a,1x,100(f20.10,1x))') 'after phseChange: mLayerVolFracIceIter(1:nSnow+2) = ', mLayerVolFracIceIter(1:nSnow+2)
  !write(*,'(a,1x,100(f20.10,1x))') 'after phseChange: mLayerVolFracIceNew(1:nSnow+2)  = ', mLayerVolFracIceNew(1:nSnow+2)
  !write(*,'(a,1x,100(f20.10,1x))') 'after phseChange: mLayerVolFracLiqIter(1:nSnow+2) = ', mLayerVolFracLiqIter(1:nSnow+2)
  !write(*,'(a,1x,100(f20.10,1x))') 'after phseChange: mLayerVolFracLiqNew(1:nSnow+2)  = ', mLayerVolFracLiqNew(1:nSnow+2)

  print*, 'after phase change = '
  write(*,'(a,10(f20.10,1x))') 'in picardSolv: mLayerTempNew        = ', mLayerTempNew(1:5)
  write(*,'(a,10(f20.10,1x))') 'in picardSolv: mLayerMatricHeadIter = ', mLayerMatricHeadIter(1:5)
  write(*,'(a,10(f20.10,1x))') 'in picardSolv: mLayerMatricHeadNew  = ', mLayerMatricHeadNew(1:5)

  ! compute total melt/freeze of canopy water
  scalarCanopyMeltFreeze = (scalarCanopyIceNew - scalarCanopyIce)/dt     ! melt/freeze of water stored in the canopy (kg m-2 s-1)

  ! compute melt/freeze of infiltrating liquid water in each layer (kg m-3 s-1) -- melt is negative
  mLayerInfilFreeze = mLayerInfilFreeze + iden_ice*(mLayerVolFracIceIter - mLayerVolFracIceNew)/dt
 
  ! compute increment in energy for the vegetation canopy
  canopyNrgIncr       = canopyTempComponent

  ! compute increment in energy for the snow-soil vector -- do here because of phase change
  mLayerInflComponent = iden_ice*LH_fus*(mLayerVolFracIceNew - mLayerVolFracIceIter)
  mLayerNrgIncr       = mLayerTempComponent - mLayerPhseComponent - mLayerInflComponent

  !print*, '*****'
  !write(*,'(a,1x,100(f20.10,1x))') 'after phseChange: mLayerPhseComponent(1:nSnow) = ', mLayerPhseComponent(1:nSnow)
  !write(*,'(a,1x,100(f20.10,1x))') 'after phseChange: mLayerInflComponent(1:nSnow) = ', mLayerInflComponent(1:nSnow)

  ! check for phase change oscillation
  !if(iter>3 .and. nSnow>0)then
  ! do iLayer=1,nSnow
  !  if(mLayerPhseComponent(iLayer)*mLayerInflComponent(iLayer) < -tiny(dt))then
  !   stop ' phase change oscillation'
  !  endif
  ! end do
  !endif


  ! =================================================================================================================================================
  ! =================================================================================================================================================


  ! *****
  ! * get ready for the next iteration.....
  ! ***************************************

  ! update scalar state variables
  scalarCanairTempIter = scalarCanairTempNew
  scalarCanopyTempIter = scalarCanopyTempNew
  scalarCanopyIceIter  = scalarCanopyIceNew
  scalarCanopyLiqIter  = scalarCanopyLiqNew
  scalarAquiferStorageIter = scalarAquiferStorageNew

  ! update the state variables in the snow-soil vector
  mLayerTempIter       = mLayerTempNew
  ! (NOTE: efficiency gains possible by altering calling sequence for phseChange)
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

  ! compute maximum residual
  hydrol_max = maxval(abs(mLayerResidual))

  ! compute maximum iteration increment
  liquid_max = maxval(abs(mLayerLiqIncr))
  matric_max = maxval(abs(mLayerMatIncr) - (absConvTol_matric + abs(relConvTol_matric*mLayerMatricHeadNew) )  )
  energy_max = maxval(abs((/canopyNrgIncr,mLayerNrgIncr/)))
  aquifr_max = abs(scalarAqiIncr)
  !print*, 'mLayerMatIncr(1:5) = ', mLayerMatIncr(1:5)
  print*, 'liquid_max, matric_max, energy_max = ', liquid_max, matric_max, energy_max

  ! get position of maximum iteration increment
  liquid_pos = maxloc(abs(mLayerLiqIncr))
  matric_pos = maxloc(abs(mLayerMatIncr))
  energy_pos = maxloc(abs((/canopyNrgIncr,mLayerNrgIncr/)))
  print*, 'liquid_pos, matric_pos, energy_pos = ', liquid_pos, matric_pos, energy_pos
  !print*, 'canopyNrgIncr = ', canopyNrgIncr
  !if(computeVegFlux) pause 'canopy is exposed'

  ! test
  !write(*,'(a25,1x,2(i4,1x),10(e20.3,1x))') 'temperature increment = ', iter, energy_pos, scalarCanopyTempIncr, mLayerTempIncr(minLayer:maxLayer)
  !write(*,'(a25,1x,i4,1x,10(f20.7,1x))') 'energy increment = ', energy_pos, canopyNrgIncr, mLayerNrgIncr(nSnow+1)
  !pause

  !print*, 'iterating..., scalarCanopyLiqNew = ', scalarCanopyLiqNew

  ! convergence check: 
  if( liquid_max(1) < absConvTol_liquid .or.                       & ! volumetric fraction of liquid water (-)
     (matric_max(1) < 0._dp .or. hydrol_max(1) < 1.e-6_dp)   .and. & ! matric head (m)
      energy_max(1) < absConvTol_energy .and.                      & ! energy (J m-3)
      aquifr_max    < absConvTol_aquifr)                           & ! aquifer storage (m)
   exit

  !print*, 'mLayerMatIncr(matric_pos(1)), absConvTol_matric + abs(relConvTol_matric*mLayerMatricHeadNew(matric_pos(1))) = ', &
  !         mLayerMatIncr(matric_pos(1)), absConvTol_matric + abs(relConvTol_matric*mLayerMatricHeadNew(matric_pos(1)))
  !pause

  ! check for lack of convergence
  if(niter==maxiter)then
   err=-30; message=trim(message)//'failed to converge'
   !if(dt < 10._dp) err=abs(err)
   return
  endif

  !print*, 'iter, dt = ', iter, dt
  !print*, 'press enter key to continue'
  !read(*,*)   ! same as a pause statement

 end do  ! (iterating)
 !print*, 'after iterations: mLayerVolFracIceNew(1) = ', mLayerVolFracIceNew(1)
 !pause 'after iterations'


 ! *****************************************************************************************************************************************
 ! *****************************************************************************************************************************************
 ! ***** END OF ITERATIONS *****************************************************************************************************************
 ! *****************************************************************************************************************************************
 ! *****************************************************************************************************************************************
 !write(*,'(a)') '================================================================================================================================='
 !write(*,'(a)') '================================================================================================================================='
 !pause 'after iterations'

 ! =================================================================================================================================================
 ! =================================================================================================================================================

 ! 
 !write(*,'(a,1x,f20.10,1x,i5))') 'dt, niter = ', dt, niter
 !write(*,'(a,1x,10(f20.10,1x))') 'scalarCanairTempNew      = ', scalarCanairTempNew
 !write(*,'(a,1x,10(f20.10,1x))') 'scalarCanopyTempNew      = ', scalarCanopyTempNew
 !write(*,'(a,1x,10(f20.10,1x))') 'scalarCanopyIceNew       = ', scalarCanopyIceNew
 !write(*,'(a,1x,10(f20.10,1x))') 'scalarCanopyLiqNew       = ', scalarCanopyLiqNew
 !write(*,'(a,1x,10(f16.10,1x))') 'mLayerTempNew(1:10)       = ', mLayerTempNew(1:10)
 !write(*,'(a,1x,10(f20.10,1x))') 'mLayerVolFracIceNew(1:3) = ', mLayerVolFracIceNew(1:3)
 !write(*,'(a,1x,10(f20.10,1x))') 'mLayerVolFracLiqNew(1:3) = ', mLayerVolFracLiqNew(1:3)
 !write(*,'(a,1x,10(f20.10,1x))') 'mLayerMatricHeadNew(1:3) = ', mLayerMatricHeadNew(1:3)

 !write(*,'(a,1x,i4,1x,20(f11.5,1x))') 'niter, dt, mLayerTempNew(1:10)       = ', niter, dt, mLayerTempNew(1:10)

 ! *****
 ! * basic checks.....
 ! *******************

 ! check SWE
 !if(nSnow > 0)then
 ! testSWE = sum( (iden_ice*mLayerVolFracIceNew(1:nSnow) + iden_water*mLayerVolFracLiqNew(1:nSnow))*mLayerDepth(1:nSnow) )
 ! print*, 'after iterations picardSolv_muster: SWE = ', testSWE
 !endif

 ! check volumetric ice content is less than the intrinsic density of ice
 do iLayer=1,nSnow
  if(mLayerVolFracIceNew(iLayer) > iden_ice/iden_water)then
   message=trim(message)//'volumetric ice content > intrinsic density of ice'
   err=20; return ! (negative error code forces a reduction in the length of the sub-step and another trial)
  endif
 end do

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
  volFrac_water = mLayerVolFracIceNew(iLayer)*(iden_ice/iden_water) + mLayerVolFracLiqNew(iLayer)
  ! check if the total volumetric water (liquid water plus ice) exceeds porosity
  if(volFrac_water > theta_sat)then
   write(*,'(a,i4,1x,3(f15.10,1x))') 'iLayer, volFrac_water, mLayerVolFracIceNew(iLayer), mLayerVolFracLiqNew(iLayer) = ',&
                                      iLayer, volFrac_water, mLayerVolFracIceNew(iLayer), mLayerVolFracLiqNew(iLayer)
   write(message,'(a,i0,a,i0,a)')trim(message)//"(liquid + ice) > porosity [iLayer=",iLayer,"; iSoil=",iLayer-nSnow,"]"
   err=-20; return
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
                       ix_spatialGW,                          & ! intent(in): choice of method for the spatial representation of groundwater

                       ! input: soil parameters
                       theta_sat,                             & ! intent(in): soil porosity (-)
                       specificStorage,                       & ! intent(in): specific storage (m-1)
 
                       ! input: state variables at the start of the sub-step
                       mLayerVolFracIce,                      & ! intent(in): volumetric fraction of ice at the start of the sub-step (-)
                       mLayerVolFracLiq,                      & ! intent(in): volumetric fraction of liquid water at the start of the sub-step (-)
                       mLayerMatricHead,                      & ! intent(in): matric head at the start of the sub-step (-)
                       scalarAquiferStorage,                  & ! intent(in): aquifer storage at the start of the sub-step (m)

                       ! input: state variables after iteration                
                       mLayerVolFracIceNew,                   & ! intent(in): volumetric fraction of ice at the end of the sub-step (-)
                       mLayerVolFracLiqNew,                   & ! intent(in): volumetric fraction of liquid water at the end of the sub-step (-)
                       mLayerMatricHeadNew,                   & ! intent(in): matric head at the end of the sub-step (-)
                       scalarAquiferStorageNew,               & ! intent(in): aquifer storage at the end of the sub-step (m)

                       ! input: model fluxes at the *START* of the sub-step
                       iLayerInitLiqFluxSoil,                 & ! intent(in): liquid water flux at the interface of each layer at the start of the sub-step (m s-1)
                       iLayerInitFluxReversal,                & ! intent(in): flow reversal flux at the interface of each layer at the start of the sub-step (m s-1)
                       mLayerInitBaseflow,                    & ! intent(in): baseflow from each layer at the start of the sub-step (m s-1)
                       mLayerInitTranspire,                   & ! intent(in): transpiration from each layer at the start of the sub-step (m s-1)
                       scalarInitAquiferRecharge,             & ! intent(in): recharge to the aquifer at the start of the sub-step (m s-1)
                       scalarInitAquiferBaseflow,             & ! intent(in): baseflow from the aquifer at the start of the sub-step (m s-1)
                       scalarInitAquiferTranspire,            & ! intent(in): transpiration from the aquifer at the start of the sub-step (m s-1)

                       ! input: model fluxes at the *END* of the sub-step
                       iLayerLiqFluxSoil,                     & ! intent(in): liquid water flux at the interface of each layer at the end of the sub-step (m s-1)
                       iLayerFluxReversal,                    & ! intent(in): flow reversal flux at the interface of each layer at the end of the sub-step (m s-1)
                       mLayerBaseflow,                        & ! intent(in): baseflow from each layer at the end of the sub-step (m s-1)
                       mLayerTranspire,                       & ! intent(in): transpiration from each layer at the end of the sub-step (m s-1)
                       scalarAquiferRecharge,                 & ! intent(in): recharge to the aquifer at the end of the sub-step (m s-1)
                       scalarAquiferBaseflow,                 & ! intent(in): baseflow from the aquifer at the end of the sub-step (m s-1)
                       scalarAquiferTranspire,                & ! intent(in): transpiration from the aquifer at the end of the sub-step (m s-1)

                       ! output: water balance of the soil zone
                       scalarSoilInflux,                      & ! intent(out): influx at the top of the soil zone (m s-1)
                       scalarSoilBaseflow,                    & ! intent(out): total baseflow from the soil zone (m s-1)
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
  qbaseTopmodel,                   & ! TOPMODEL-ish baseflow parameterization
  bigBucket,                       & ! a big bucket (lumped aquifer model)
  noExplicit                         ! no explicit groundwater parameterization
 ! look-up values for the choice of groundwater representation (local-column, or single-basin)
 USE mDecisions_module,only:       &
  localColumn,                     & ! separate groundwater representation in each local soil column
  singleBasin                        ! single groundwater store over the entire basin
 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 ! input: model control
 real(dp),intent(in)            :: dt                           ! time step (seconds)
 real(dp),intent(in)            :: wimplicit                    ! weight assigned to start-of-step fluxes (-)
 real(dp),intent(in)            :: mLayerDepth(:)               ! depth of each layer (m)
 integer(i4b),intent(in)        :: ix_groundwatr                ! choice of groundwater parameterization
 integer(i4b),intent(in)        :: ix_spatialGW                 ! choice of method for the spatial representation of groundwater
 ! input: soil parameters
 real(dp),intent(in)            :: theta_sat                    ! soil porosity (-) 
 real(dp),intent(in)            :: specificStorage              ! specific storage (m-1)
 ! input: state variables at the start of the sub-step
 real(dp),intent(in)            :: mLayerVolFracIce(:)          ! volumetric fraction of ice in each layer (-)
 real(dp),intent(in)            :: mLayerVolFracLiq(:)          ! volumetric fraction of liquid water in each layer (-)
 real(dp),intent(in)            :: mLayerMatricHead(:)          ! matric head in each layer (-)
 real(dp),intent(in)            :: scalarAquiferStorage         ! aquifer storage at the start of the sub-step (m)
 ! input: state variables after iteration                
 real(dp),intent(in)            :: mLayerVolFracIceNew(:)       ! volumetric fraction of ice at the end of the sub-step (-)
 real(dp),intent(in)            :: mLayerVolFracLiqNew(:)       ! volumetric fraction of liquid water at the end of the sub-step (-)
 real(dp),intent(in)            :: mLayerMatricHeadNew(:)       ! matric head at the end of the sub-step (-)
 real(dp),intent(in)            :: scalarAquiferStorageNew      ! aquifer storage at the end of the sub-step (m)
 ! input: model fluxes at the *START* of the sub-step
 real(dp),intent(in)            :: iLayerInitLiqFluxSoil(0:)    ! liquid water flux at the interface of each layer at the start of the sub-step (m s-1)
 real(dp),intent(in)            :: iLayerInitFluxReversal(0:)   ! flow reversal flux at the interface of each layer at the start of the sub-step (m s-1)
 real(dp),intent(in)            :: mLayerInitBaseflow(:)        ! baseflow from each layer at the start of the sub-step (m s-1)
 real(dp),intent(in)            :: mLayerInitTranspire(:)       ! transpiration from each layer at the start of the sub-step (m s-1)
 real(dp),intent(in)            :: scalarInitAquiferRecharge    ! recharge to the aquifer at the end of the sub-step (m s-1)
 real(dp),intent(in)            :: scalarInitAquiferBaseflow    ! baseflow from the aquifer at the end of the sub-step (m s-1)
 real(dp),intent(in)            :: scalarInitAquiferTranspire   ! transpiration from the aquifer at the end of the sub-step (m s-1)
 ! input: model fluxes at the *END* of the sub-step
 real(dp),intent(in)            :: iLayerLiqFluxSoil(0:)        ! liquid water flux at the interface of each layer at the end of the sub-step (m s-1)
 real(dp),intent(in)            :: iLayerFluxReversal(0:)       ! flow reversal flux at the interface of each layer at the end of the sub-step (m s-1)
 real(dp),intent(in)            :: mLayerBaseflow(:)            ! baseflow from each layer at the end of the sub-step (m s-1)
 real(dp),intent(in)            :: mLayerTranspire(:)           ! transpiration from each layer at the end of the sub-step (m s-1)
 real(dp),intent(in)            :: scalarAquiferRecharge        ! recharge to the aquifer at the end of the sub-step (m s-1)
 real(dp),intent(in)            :: scalarAquiferBaseflow        ! baseflow from the aquifer at the end of the sub-step (m s-1)
 real(dp),intent(in)            :: scalarAquiferTranspire       ! transpiration from the aquifer at the end of the sub-step (m s-1)
 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 ! output: water balance of the soil zone
 real(dp),intent(out)           :: scalarSoilInflux             ! intent(out): influx at the top of the soil zone (m s-1)
 real(dp),intent(out)           :: scalarSoilBaseflow           ! intent(out): total baseflow from the soil zone (m s-1)
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
 real(dp)                       :: stor_Change                  ! change in volumetric liquid water content from water released due to change in pressure
 real(dp)                       :: scalarSoilRelease            ! total water released due to change in pressure (kg m-2)
 real(dp)                       :: balanceSoilWater0            ! total soil storage at the start of the step (kg m-2)
 real(dp)                       :: balanceSoilWater1            ! total soil storage at the end of the step (kg m-2)
 real(dp)                       :: balanceSoilInflux            ! input to the soil zone
 real(dp)                       :: balanceSoilBaseflow          ! output from the soil zone
 real(dp)                       :: balanceSoilDrainage          ! output from the soil zone
 real(dp)                       :: balanceSoilTranspiration     ! output from the soil zone
 real(dp)                       :: balanceAquifer0              ! total aquifer storage at the start of the step (kg m-2)
 real(dp)                       :: balanceAquifer1              ! total aquifer storage at the end of the step (kg m-2)
 ! ------------------------------------------------------------------------------------------------------------------------------------------------
 ! initialize error control
 err=0; message="wBal_check/"

 ! compute total soil moisture and ice at the *START* of the step (kg m-2)
 scalarTotalSoilLiq = sum(iden_water*mLayerVolFracLiq(nSnow+1:nLayers)*mLayerDepth(nSnow+1:nLayers))
 scalarTotalSoilIce = sum(iden_water*mLayerVolFracIce(nSnow+1:nLayers)*mLayerDepth(nSnow+1:nLayers))
 !scalarTotalSoilIce = sum(iden_ice  *mLayerVolFracIce(nSnow+1:nLayers)*mLayerDepth(nSnow+1:nLayers))

 ! get the total water in the soil (liquid plus ice) at the start of the time step (kg m-2)
 balanceSoilWater0 = scalarTotalSoilLiq + scalarTotalSoilIce

 ! get the total aquifer storage at the start of the time step (kg m-2)
 balanceAquifer0 = scalarAquiferStorage*iden_water

 ! compute total soil moisture and ice at the *END* of the step (kg m-2)
 scalarTotalSoilLiq = sum(iden_water*mLayerVolFracLiqNew(nSnow+1:nLayers)*mLayerDepth(nSnow+1:nLayers))
 scalarTotalSoilIce = sum(iden_water*mLayerVolFracIceNew(nSnow+1:nLayers)*mLayerDepth(nSnow+1:nLayers))
 !scalarTotalSoilIce = sum(iden_ice  *mLayerVolFracIceNew(nSnow+1:nLayers)*mLayerDepth(nSnow+1:nLayers))

 ! get the total water in the soil (liquid plus ice) at the end of the time step (kg m-2)
 balanceSoilWater1 = scalarTotalSoilLiq + scalarTotalSoilIce

 ! get the total aquifer storage at the start of the time step (kg m-2)
 balanceAquifer1 = scalarAquiferStorageNew*iden_water

 ! compute the influx, drainage of water from the soil profile (m s-1)
 scalarSoilInflux        = (wimplicit*iLayerInitLiqFluxSoil(0)       + (1._dp - wimplicit)*iLayerLiqFluxSoil(0)      ) + &
                           (wimplicit*iLayerInitFluxReversal(0)      + (1._dp - wimplicit)*iLayerFluxReversal(0)     )
 scalarSoilBaseflow      = (wimplicit*sum(mLayerInitBaseflow)        + (1._dp - wimplicit)*sum(mLayerBaseflow)       )
 scalarSoilDrainage      = (wimplicit*iLayerInitLiqFluxSoil(nLevels) + (1._dp - wimplicit)*iLayerLiqFluxSoil(nLevels))
 scalarSoilTranspiration = (wimplicit*sum(mLayerInitTranspire)       + (1._dp - wimplicit)*sum(mLayerTranspire)      ) 
 !print*, 'iLayerLiqFluxSoil(0), scalarSoilInflux = ', iLayerLiqFluxSoil(0), scalarSoilInflux

 ! get the total water released due to change in pressure (kg m-2)
 scalarSoilRelease = sum( mLayerDepth(nSnow+1:nLayers)*(mLayerMatricHeadNew(1:nSoil) - mLayerMatricHead(1:nSoil))*mLayerVolFracLiqNew(nSnow+1:nLayers)*specificStorage/theta_sat )*iden_water

 ! get the input and output to/from the soil zone (kg m-2)
 balanceSoilInflux        = scalarSoilInflux*iden_water*dt
 balanceSoilBaseflow      = scalarSoilBaseflow*iden_water*dt
 balanceSoilDrainage      = scalarSoilDrainage*iden_water*dt
 balanceSoilTranspiration = scalarSoilTranspiration*iden_water*dt

 ! check the soil water balance
 scalarSoilWatBalError  = balanceSoilWater1 - (balanceSoilWater0 + (balanceSoilInflux + balanceSoilTranspiration - balanceSoilBaseflow - balanceSoilDrainage - scalarSoilRelease) )
 if(abs(scalarSoilWatBalError) > 1.d-3)then
  ! check the balance of each layer
  write(*,'(a)') 'water balance of each layer'
  write(*,'(a)') 'Liq0 (-), Liq1 (-), Ice0 (-), Ice1 (-), totalChange (-), phaseChange (-), flux_Change, evap_Change, qbaseChange, stor_Change',&
                 'phaseChange+flux_Change+evap_Change-qbaseChange-stor_change, totalChange - (phaseChange+flux_Change+evap_Change-qbaseChange-stor_change)'
  do iLayer=1,min(10,nSoil)
   totalChange = mLayerVolFracLiqNew(iLayer+nSnow) - mLayerVolFracLiq(iLayer+nSnow) ! total change in volumetric liquid water content
   !phaseChange = -(iden_ice/iden_water)*(mLayerVolFracIceNew(iLayer+nSnow) - mLayerVolFracIce(iLayer+nSnow))  ! change in liquid water content associated with freezing
   phaseChange = -(mLayerVolFracIceNew(iLayer+nSnow) - mLayerVolFracIce(iLayer+nSnow))  ! change in liquid water content associated with freezing (NOTE: no volume expansion)
   evap_Change = dt*mLayerTranspire(iLayer)/mLayerDepth(iLayer)
   qbaseChange = dt*mLayerBaseflow(iLayer)/mLayerDepth(iLayer)
   flux_Change = dt*(iLayerLiqFluxSoil(iLayer-1)  - iLayerLiqFluxSoil(iLayer))/mLayerDepth(iLayer) + & ! change in volumetric liquid water content from the interface fluxes
                 dt*(iLayerFluxReversal(iLayer-1) - iLayerFluxReversal(iLayer))/mLayerDepth(iLayer)     ! change in volumetric liquid water content from the interface fluxes
   stor_Change = (mLayerMatricHeadNew(iLayer) - mLayerMatricHead(iLayer))*specificStorage*mLayerVolFracLiqNew(iLayer+nSnow)/theta_sat

   write(*,'(i4,1x,2(f15.8,1x),20(e15.5,1x))') iLayer, mLayerVolFracLiq(iLayer+nSnow), mLayerVolFracLiqNew(iLayer+nSnow), &
                                                       mLayerVolFracIce(iLayer+nSnow), mLayerVolFracIceNew(iLayer+nSnow), &
                                                       totalChange, phaseChange, flux_Change, evap_Change, qbaseChange, stor_change, &
                                                       phaseChange+flux_Change+evap_Change-qbaseChange-stor_change, &
                                                       totalChange - (phaseChange+flux_Change+evap_Change-qbaseChange-stor_change)
  end do
  ! print the total water balance
  print*, 'dt = ', dt
  print*, 'balanceSoilWater0 (kg m-2) = ',        balanceSoilWater0
  print*, 'balanceSoilWater1 (kg m-2) = ',        balanceSoilWater1
  print*, 'balanceSoilInflux (kg m-2) = ',        balanceSoilInflux
  print*, 'balanceSoilDrainage (kg m-2) = ',      balanceSoilDrainage
  print*, 'balanceSoilBaseflow (kg m-2) = ',      balanceSoilBaseflow
  print*, 'balanceSoilTranspiration (kg m-2) = ', balanceSoilTranspiration, sum(mLayerTranspire)*dt*iden_water
  if(abs(scalarSoilWatBalError) > 1.d-3)then
   write(message,'(a,e20.10,a)')trim(message)//"abs(scalarSoilWatBalError) > 1.d-3 [error = ",&
                                scalarSoilWatBalError," ]"
   err=-10; return
  endif
 endif

 ! check the aquifer water balance
 ! NOTE: only check if the aquifer is defined in the local soil column
 if(ix_spatialGW == localColumn)then
  select case(ix_groundwatr)
   ! no explicit aquifer
   case(noExplicit,qbaseTopmodel)
    scalarAquiferBalError = 0._dp
   ! explicit aquifer
   case(bigBucket)
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
 endif  ! (if aquifer is defined in the local column)

 end subroutine wBal_check


 ! *****************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************
 ! *****************************************************************************************************************************************************





 ! *********************************************************************************************************
 ! private subroutine: compute melt of the "snow without a layer"
 ! *********************************************************************************************************
 subroutine implctMelt(&
                       ! input/output: integrated snowpack properties
                       scalarSWE,         & ! snow water equivalent (kg m-2)
                       scalarSnowDepth,   & ! snow depth (m)
                       scalarSfcMeltPond, & ! surface melt pond (kg m-2)
                       ! input/output: properties of the upper-most soil layer
                       soilTemp,          & ! surface layer temperature (K)
                       soilDepth,         & ! surface layer depth (m)
                       soilHeatcap,       & ! surface layer volumetric heat capacity (J m-3 K-1)
                       ! output: error control
                       err,message        )
 implicit none
 ! input/output: integrated snowpack properties
 real(dp),intent(inout)    :: scalarSWE          ! snow water equivalent (kg m-2)
 real(dp),intent(inout)    :: scalarSnowDepth    ! snow depth (m)
 real(dp),intent(inout)    :: scalarSfcMeltPond  ! surface melt pond (kg m-2)
 ! input/output: properties of the upper-most soil layer
 real(dp),intent(inout)    :: soilTemp           ! surface layer temperature (K)
 real(dp),intent(inout)    :: soilDepth          ! surface layer depth (m)
 real(dp),intent(inout)    :: soilHeatcap        ! surface layer volumetric heat capacity (J m-3 K-1)
 ! output: error control
 integer(i4b),intent(out)  :: err                ! error code
 character(*),intent(out)  :: message            ! error message
 ! local variables
 real(dp)                  :: nrgRequired        ! energy required to melt all the snow (J m-2)
 real(dp)                  :: nrgAvailable       ! energy available to melt the snow (J m-2)
 real(dp)                  :: snwDensity         ! snow density (kg m-3)
 ! initialize error control
 err=0; message='implctMelt/'

 ! check for the special case of "snow without a layer"
 if (nSnow==0 .and. scalarSWE > 0._dp)then
  ! only melt if temperature of the top soil layer is greater than Tfreeze
  if(soilTemp > Tfreeze)then
   ! compute the energy required to melt all the snow (J m-2)
   nrgRequired     = scalarSWE*LH_fus
   ! compute the energy available to melt the snow (J m-2)
   nrgAvailable    = soilHeatcap*(soilTemp - Tfreeze)*soilDepth
   ! compute the snow density (not saved)
   snwDensity      = scalarSWE/scalarSnowDepth
   ! compute the amount of melt, and update SWE (kg m-2)
   if(nrgAvailable > nrgRequired)then
    scalarSfcMeltPond  = scalarSWE
    scalarSWE          = 0._dp
   else
    scalarSfcMeltPond  = nrgAvailable/LH_fus
    scalarSWE          = scalarSWE - scalarSfcMeltPond
   endif
   ! update depth
   scalarSnowDepth = scalarSWE/snwDensity
   ! update temperature of the top soil layer (K)
   soilTemp =  soilTemp - (LH_fus*scalarSfcMeltPond/soilDepth)/soilHeatcap
  else  ! melt is zero if the temperature of the top soil layer is less than Tfreeze
   scalarSfcMeltPond = 0._dp  ! kg m-2
  endif ! (if the temperature of the top soil layer is greater than Tfreeze)
 else  ! melt is zero if the "snow without a layer" does not exist
  scalarSfcMeltPond = 0._dp  ! kg m-2
 endif ! (if the "snow without a layer" exists)

 end subroutine implctMelt





end module picardSolv_module

module coupled_em_module
USE nrtype
implicit none
private
public::coupled_em
! algorithmic parameters
real(dp),parameter     :: valueMissing=-9999._dp  ! missing value, used when diagnostic or state variables are undefined
real(dp),parameter     :: verySmall=1.e-6_dp   ! used as an additive constant to check if substantial difference among real numbers
real(dp),parameter     :: mpe=1.e-6_dp         ! prevents overflow error if division by zero 
real(dp),parameter     :: dx=1.e-6_dp          ! finite difference increment
contains

 ! ************************************************************************************************
 ! new subroutine: run the coupled energy-mass model for one timestep
 ! ************************************************************************************************
 subroutine coupled_em(dt_init,err,message)
 ! data structures and named variables
 USE data_struc,only:data_step                                                      ! time step of forcing data (s)
 USE data_struc,only:model_decisions                                                ! model decision structure
 USE data_struc,only:type_data,attr_data,forc_data,mpar_data,mvar_data,indx_data    ! data structures
 USE var_lookup,only:iLookDECISIONS                                                 ! named variables for elements of the decision structure
 USE data_struc,only:ix_soil,ix_snow                                                ! named variables for snow and soil
 USE var_lookup,only:iLookTYPE,iLookATTR,iLookFORCE,iLookPARAM,iLookMVAR,iLookINDEX ! named variables for structure elements
 ! common variables
 USE data_struc,only:urbanVegCategory       ! vegetation category for urban areas
 USE data_struc,only:fracJulday             ! fractional julian days since the start of year
 USE data_struc,only:yearLength             ! number of days in the current year
 ! subroutines
 USE NOAHMP_ROUTINES,only:phenology         ! determine vegetation phenology
 USE canopySnow_module,only:canopySnow      ! compute interception and unloading of snow from the vegetation canopy
 USE snowAlbedo_module,only:snowAlbedo      ! compute snow albedo
 USE newsnwfall_module,only:newsnwfall      ! compute new snowfall
 USE layerMerge_module,only:layerMerge      ! merge snow layers if they are too thin
 USE layerDivide_module,only:layerDivide    ! sub-divide layers if they are too thick
 USE picardSolv_module,only:picardSolv      ! provide access to the Picard solver
 ! look-up values for the numerical method
 USE mDecisions_module,only:      &
  iterative,                      &         ! iterative
  nonIterative,                   &         ! non-iterative
  iterSurfEnergyBal                         ! iterate only on the surface energy balance
 ! look-up values for the choice of canopy shortwave radiation method
 USE mDecisions_module,only:      &
  noah_mp,                        &         ! full Noah-MP implementation (including albedo)
  CLM_2stream,                    &         ! CLM 2-stream model (see CLM documentation)
  UEB_2stream,                    &         ! UEB 2-stream model (Mahat and Tarboton, WRR 2011)
  NL_scatter,                     &         ! Simplified method Nijssen and Lettenmaier (JGR 1999)
  BeersLaw                                  ! Beer's Law (as implemented in VIC)
 USE multiconst,only:&
                     Tfreeze,     &         ! temperature at freezing              (K)
                     LH_fus,      &         ! latent heat of fusion                (J kg-1)
                     Cp_ice,      &         ! specific heat of ice                 (J kg-1 K-1)
                     Cp_water,    &         ! specific heat of liquid water        (J kg-1 K-1)
                     iden_ice,    &         ! intrinsic density of ice             (kg m-3)
                     iden_water             ! intrinsic density of liquid water    (kg m-3)
 implicit none
 ! define output
 real(dp),intent(inout)               :: dt_init                ! used to initialize the size of the sub-step
 integer(i4b),intent(out)             :: err                    ! error code
 character(*),intent(out)             :: message                ! error message
 ! control the length of the sub-step
 real(dp)                             :: dt                     ! length of time step (seconds)
 real(dp)                             :: dt_sub                 ! length of the sub-step (seconds)
 real(dp)                             :: dt_done                ! length of time step completed (seconds)
 integer(i4b)                         :: nsub                   ! number of sub-steps
 integer(i4b)                         :: niter                  ! number of iterations
 integer(i4b),parameter               :: n_inc=5                ! minimum number of iterations to increase time step
 integer(i4b),parameter               :: n_dec=20               ! maximum number of iterations to decrease time step
 real(dp),parameter                   :: F_inc = 1.25_dp        ! factor used to increase time step
 real(dp),parameter                   :: F_dec = 0.90_dp        ! factor used to decrease time step
 integer(i4b)                         :: maxiter                ! maxiumum number of iterations
 integer(i4b)                         :: iSnow                  ! index for snow layers
 ! local pointers to model forcing data
 real(dp),pointer                     :: scalarRainfall         ! rainfall flux (kg m-2 s-1)
 real(dp),pointer                     :: scalarSnowfall         ! snowfall flux (kg m-2 s-1)
 ! local pointers to model index variables
 integer(i4b),pointer                 :: nSoil                  ! number of soil layers
 integer(i4b),pointer                 :: nSnow                  ! number of snow layers
 integer(i4b),pointer                 :: nLayers                ! number of layers
 integer(i4b),pointer                 :: layerType(:)           ! type of the layer (ix_soil or ix_snow)
 ! local pointers to variables defining melt of the "snow without a layer"
 real(dp),pointer                     :: scalarSWE              ! snow water equivalent (kg m-2)
 real(dp),pointer                     :: scalarSnowDepth        ! snow depth (m)
 real(dp),pointer                     :: scalarSfcMeltPond      ! surface melt pond (kg m-2)
 real(dp),pointer                     :: mLayerVolHtCapBulk(:)  ! volumetric heat capacity (J m-3 K-1)
 ! local pointers to model state variables -- all layers
 real(dp),pointer                     :: mLayerTemp(:)          ! temperature of each layer (K)
 real(dp),pointer                     :: mLayerDepth(:)         ! depth of each layer (m)
 real(dp),pointer                     :: mLayerVolFracIce(:)    ! volumetric fraction of ice in each layer (-)
 real(dp),pointer                     :: mLayerVolFracLiq(:)    ! volumetric fraction of liquid water in each layer (-) 
 ! local pointers to flux variables
 real(dp),pointer                     :: scalarThroughfallSnow       ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp),pointer                     :: scalarThroughfallRain       ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp),pointer                     :: scalarCanopySnowUnloading   ! unloading of snow from the vegetion canopy (kg m-2 s-1)
 real(dp),pointer                     :: scalarCanopyLiqDrainage     ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
 real(dp),pointer                     :: scalarCanopyMeltFreeze      ! melt/freeze of water stored in the canopy (kg m-2 s-1)
 real(dp),pointer                     :: scalarSnowSublimation       ! snow sublimation/frost - below canopy or non-vegetated (kg m-2 s-1)
 real(dp),pointer                     :: scalarGroundEvaporation     ! ground evaporation/condensation - below canopy or non-vegetated (kg m-2 s-1)
 real(dp),pointer                     :: scalarRainPlusMelt          ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
 real(dp),pointer                     :: scalarSurfaceRunoff         ! surface runoff (m s-1) 
 real(dp),pointer                     :: scalarSoilInflux            ! influx of water at the top of the soil profile (m s-1)
 real(dp),pointer                     :: scalarSoilBaseflow          ! total baseflow from throughout the soil profile (m s-1)
 real(dp),pointer                     :: scalarSoilDrainage          ! drainage from the bottom of the soil profile (m s-1)
 real(dp),pointer                     :: scalarAquiferRecharge       ! recharge to the aquifer (m s-1)
 real(dp),pointer                     :: scalarAquiferBaseflow       ! baseflow from the aquifer (m s-1)
 real(dp),pointer                     :: scalarAquiferTranspire      ! transpiration from the aquifer (m s-1)
 real(dp),pointer                     :: mLayerColumnOutflow(:)      ! total outflow from each layer in a given soil column (m3 s-1)
 ! local pointers to timestep-average flux variables
 real(dp),pointer                     :: averageThroughfallSnow      ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp),pointer                     :: averageThroughfallRain      ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
 real(dp),pointer                     :: averageCanopySnowUnloading  ! unloading of snow from the vegetion canopy (kg m-2 s-1)
 real(dp),pointer                     :: averageCanopyLiqDrainage    ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
 real(dp),pointer                     :: averageCanopyMeltFreeze     ! melt/freeze of water stored in the canopy (kg m-2 s-1)
 real(dp),pointer                     :: averageSnowSublimation      ! snow sublimation/frost - below canopy or non-vegetated (kg m-2 s-1)
 real(dp),pointer                     :: averageGroundEvaporation    ! ground evaporation/condensation - below canopy or non-vegetated (kg m-2 s-1)
 real(dp),pointer                     :: averageRainPlusMelt         ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
 real(dp),pointer                     :: averageSurfaceRunoff        ! surface runoff (m s-1) 
 real(dp),pointer                     :: averageSoilInflux           ! influx of water at the top of the soil profile (m s-1)
 real(dp),pointer                     :: averageSoilBaseflow         ! total baseflow from throughout the soil profile (m s-1)
 real(dp),pointer                     :: averageSoilDrainage         ! drainage from the bottom of the soil profile (m s-1)
 real(dp),pointer                     :: averageAquiferRecharge      ! recharge to the aquifer (m s-1)
 real(dp),pointer                     :: averageAquiferBaseflow      ! baseflow from the aquifer (m s-1)
 real(dp),pointer                     :: averageAquiferTranspire     ! transpiration from the aquifer (m s-1)
 real(dp),pointer                     :: averageColumnOutflow(:)     ! outflow from each layer in the soil profile (m3 s-1)
 ! local pointers to algorithmic control parameters
 real(dp),pointer                     :: minstep                ! minimum time step length (s)
 real(dp),pointer                     :: maxstep                ! maximum time step length (s)
 ! define local variables
 character(len=256)                   :: cmessage               ! error message
 real(dp)                             :: exposedVAI             ! exposed vegetation area index (LAI + SAI)
 real(dp)                             :: canopyDepth            ! canopy depth (m)
 real(dp)                             :: heightAboveSnow        ! height top of canopy is above the snow surface (m)
 logical(lgt)                         :: computeVegFlux         ! flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
 integer(i4b)                         :: nLayersRoots           ! number of soil layers that contain roots
 real(dp),parameter                   :: verySmall=epsilon(1._dp)  ! a very small number
 real(dp)                             :: notUsed_heightCanopyTop   ! for some reason the Noah-MP phenology routines output canopy height
 real(dp),dimension(:),allocatable    :: arrTemp                ! temporary array, used for testing
 real(dp)                             :: nrgRequired            ! case of "snow without a layer": energy required to melt all the snow (J m-2)
 real(dp)                             :: nrgAvailable           ! case of "snow without a layer": energy available to melt the snow (J m-2)
 real(dp)                             :: snwDensity             ! case of "snow without a layer": snow density (kg m-3)
 real(dp)                             :: dt_wght                ! weight applied to each sub-step, to compute time step average
 integer(i4b)                         :: iLayer                 ! index of model layers
 ! initialize error control
 err=0; message="coupled_em/"

 ! assign local pointers to the model index structures
 nSoil             => indx_data%var(iLookINDEX%nSoil)%dat(1)     ! number of soil layers
 nSnow             => indx_data%var(iLookINDEX%nSnow)%dat(1)     ! number of snow layers
 nLayers           => indx_data%var(iLookINDEX%nLayers)%dat(1)   ! total number of layers
 layerType         => indx_data%var(iLookINDEX%layerType)%dat    ! layer type (ix_soil or ix_snow)

 ! identify the number of snow and soil layers
 nSnow = count(layerType==ix_snow)
 nSoil = count(layerType==ix_soil)

 ! assign pointers to model diagnostic variables -- surface scalars 
 scalarRainfall    => mvar_data%var(iLookMVAR%scalarRainfall)%dat(1)     ! rainfall flux (kg m-2 s-1)
 scalarSnowfall    => mvar_data%var(iLookMVAR%scalarSnowfall)%dat(1)     ! snowfall flux (kg m-2 s-1)

 ! assign pointers to timestep-average model fluxes
 averageThroughfallSnow     => mvar_data%var(iLookMVAR%averageThroughfallSnow)%dat(1)     ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
 averageThroughfallRain     => mvar_data%var(iLookMVAR%averageThroughfallRain)%dat(1)     ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
 averageCanopySnowUnloading => mvar_data%var(iLookMVAR%averageCanopySnowUnloading)%dat(1) ! unloading of snow from the vegetion canopy (kg m-2 s-1)
 averageCanopyLiqDrainage   => mvar_data%var(iLookMVAR%averageCanopyLiqDrainage)%dat(1)   ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
 averageCanopyMeltFreeze    => mvar_data%var(iLookMVAR%averageCanopyMeltFreeze)%dat(1)    ! melt/freeze of water stored in the canopy (kg m-2 s-1)
 averageSnowSublimation     => mvar_data%var(iLookMVAR%averageSnowSublimation)%dat(1)     ! snow sublimation/frost - below canopy or non-vegetated (kg m-2 s-1)
 averageGroundEvaporation   => mvar_data%var(iLookMVAR%averageGroundEvaporation)%dat(1)   ! ground evaporation/condensation - below canopy or non-vegetated (kg m-2 s-1)
 averageRainPlusMelt        => mvar_data%var(iLookMVAR%averageRainPlusMelt)%dat(1)        ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
 averageSurfaceRunoff       => mvar_data%var(iLookMVAR%averageSurfaceRunoff)%dat(1)       ! surface runoff (m s-1)
 averageSoilInflux          => mvar_data%var(iLookMVAR%averageSoilInflux)%dat(1)          ! influx of water at the top of the soil profile (m s-1)
 averageSoilBaseflow        => mvar_data%var(iLookMVAR%averageSoilBaseflow)%dat(1)        ! total baseflow from throughout the soil profile (m s-1)
 averageSoilDrainage        => mvar_data%var(iLookMVAR%averageSoilDrainage)%dat(1)        ! drainage from the bottom of the soil profile (m s-1)
 averageAquiferRecharge     => mvar_data%var(iLookMVAR%averageAquiferRecharge)%dat(1)     ! recharge to the aquifer (m s-1)
 averageAquiferBaseflow     => mvar_data%var(iLookMVAR%averageAquiferBaseflow)%dat(1)     ! baseflow from the aquifer (m s-1)
 averageAquiferTranspire    => mvar_data%var(iLookMVAR%averageAquiferTranspire)%dat(1)    ! transpiration from the aquifer (m s-1)
 averageColumnOutflow       => mvar_data%var(iLookMVAR%averageColumnOutflow)%dat          ! outflow from each layer in the soil profile (m3 s-1)

 ! assign pointers to algorithmic control parameters
 minstep => mpar_data%var(iLookPARAM%minstep)  ! minimum time step (s)
 maxstep => mpar_data%var(iLookPARAM%maxstep)  ! maximum time step (s)
 !print*, 'minstep, maxstep = ', minstep, maxstep

 ! initialize average fluxes
 averageThroughfallSnow     = 0._dp  ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
 averageThroughfallRain     = 0._dp  ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
 averageCanopySnowUnloading = 0._dp  ! unloading of snow from the vegetion canopy (kg m-2 s-1)
 averageCanopyLiqDrainage   = 0._dp  ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
 averageCanopyMeltFreeze    = 0._dp  ! melt/freeze of water stored in the canopy (kg m-2 s-1)
 averageSurfaceRunoff       = 0._dp  ! surface runoff (m s-1)
 averageSnowSublimation     = 0._dp  ! snow sublimation/frost - below canopy or non-vegetated (kg m-2 s-1)
 averageGroundEvaporation   = 0._dp  ! ground evaporation/condensation - below canopy or non-vegetated (kg m-2 s-1)
 averageRainPlusMelt        = 0._dp  ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
 averageSoilInflux          = 0._dp  ! influx of water at the top of the soil profile (m s-1)
 averageSoilBaseflow        = 0._dp  ! total baseflow from throughout the soil profile (m s-1)
 averageSoilDrainage        = 0._dp  ! drainage from the bottom of the soil profile (m s-1)
 averageAquiferRecharge     = 0._dp  ! recharge to the aquifer (m s-1)
 averageAquiferBaseflow     = 0._dp  ! baseflow from the aquifer (m s-1)
 averageAquiferTranspire    = 0._dp  ! transpiration from the aquifer (m s-1)
 averageColumnOutflow       = 0._dp  ! outflow from each layer in the soil profile (m3 s-1)

 ! get the length of the time step (seconds)
 dt = data_step

 ! identify the maximum number of iterations
 select case(model_decisions(iLookDECISIONS%num_method)%iDecision)
  case(iterative);         maxiter=nint(mpar_data%var(iLookPARAM%maxiter))  ! iterative
  case(nonIterative);      maxiter=1              ! non-iterative
  case(iterSurfEnergyBal); maxiter=1              ! iterate only on the surface energy balance
   err=90; message=trim(message)//'numerical method "iterSurfEnergyBal" is not implemented yet'; return
  case default
   err=10; message=trim(message)//'unknown option for the numerical method'; return
 end select

 ! initialize the length of the sub-step
 dt_sub  = min(dt_init,min(dt,maxstep))
 dt_done = 0._dp

 ! initialize the number of sub-steps
 nsub=0

 ! loop through sub-steps
 do  ! continuous do statement with exit clause (alternative to "while")

  ! increment the number of sub-steps
  nsub = nsub+1

  ! compute the root zone temperature (used in vegetation phenology)
  ! (compute the number of layers with roots)
  nLayersRoots = count(mvar_data%var(iLookMVAR%iLayerHeight)%dat(nSnow:nLayers-1) < mpar_data%var(iLookPARAM%rootingDepth)-verySmall)
  if(nLayersRoots == 0)then; err=20; message=trim(message)//'no roots within the soil profile'; return; endif
  ! (compute the temperature of the root zone)
  mvar_data%var(iLookMVAR%scalarRootZoneTemp)%dat(1) = sum(mvar_data%var(iLookMVAR%mLayerTemp)%dat(nSnow+1:nSnow+nLayersRoots)) / real(nLayersRoots, kind(dp))

  ! define the foliage nitrogen factor
  mvar_data%var(iLookMVAR%scalarFoliageNitrogenFactor)%dat(1) = 1._dp  ! foliage nitrogen concentration (1.0 = saturated)

  ! determine vegetation phenology
  ! NOTE: recomputing phenology every sub-step accounts for changes in exposed vegetation associated with changes in snow depth
  call phenology(&
                 ! input
                 type_data%var(iLookTYPE%vegTypeIndex),                       & ! intent(in): vegetation type index
                 urbanVegCategory,                                            & ! intent(in): vegetation category for urban areas               
                 mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1),             & ! intent(in): snow depth on the ground surface (m)
                 mvar_data%var(iLookMVAR%scalarCanopyTemp)%dat(1),            & ! intent(in): temperature of the vegetation canopy at the start of the sub-step (K)
                 attr_data%var(iLookATTR%latitude),                           & ! intent(in): latitude
                 yearLength,                                                  & ! intent(in): number of days in the current year
                 fracJulday,                                                  & ! intent(in): fractional julian days since the start of year
                 mvar_data%var(iLookMVAR%scalarLAI)%dat(1),                   & ! intent(inout): one-sided leaf area index (m2 m-2)
                 mvar_data%var(iLookMVAR%scalarSAI)%dat(1),                   & ! intent(inout): one-sided stem area index (m2 m-2)
                 mvar_data%var(iLookMVAR%scalarRootZoneTemp)%dat(1),          & ! intent(in): root zone temperature (K)
                 ! output
                 notUsed_heightCanopyTop,                                     & ! intent(out): height of the top of the canopy layer (m)
                 mvar_data%var(iLookMVAR%scalarExposedLAI)%dat(1),            & ! intent(out): exposed leaf area index after burial by snow (m2 m-2)
                 mvar_data%var(iLookMVAR%scalarExposedSAI)%dat(1),            & ! intent(out): exposed stem area index after burial by snow (m2 m-2)
                 mvar_data%var(iLookMVAR%scalarGrowingSeasonIndex)%dat(1)     ) ! intent(out): growing season index (0=off, 1=on)

  ! determine if need to include vegetation in the energy flux routines
  exposedVAI      = mvar_data%var(iLookMVAR%scalarExposedLAI)%dat(1) + mvar_data%var(iLookMVAR%scalarExposedSAI)%dat(1)   ! exposed vegetation area index (m2 m-2)
  heightAboveSnow = mpar_data%var(iLookPARAM%heightCanopyTop) - mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1)           ! height top of canopy is above the snow surface (m)
  computeVegFlux  = (exposedVAI > 0.05_dp .and. heightAboveSnow > 0.05_dp)
  ! check that we have not buried trees
  if(.not.computeVegFlux .and. mpar_data%var(iLookPARAM%heightCanopyTop) > 5._dp)then
   message=trim(message)//'unexpected burial of trees'
   err=20; return
  endif

  ! compute the canopy depth (m)
  canopyDepth = mpar_data%var(iLookPARAM%heightCanopyTop) - mpar_data%var(iLookPARAM%heightCanopyBottom)
  if(mpar_data%var(iLookPARAM%heightCanopyBottom) > mpar_data%var(iLookPARAM%heightCanopyTop))then
   err=20; message=trim(message)//'height of the bottom of the canopy > top of the canopy'; return
  endif

  ! compute the bulk volumetric heat capacity of vegetation (J m-3 K-1)
  if(computeVegFlux)then
   mvar_data%var(iLookMVAR%scalarBulkVolHeatCapVeg)%dat(1)  = mpar_data%var(iLookPARAM%specificHeatVeg)*mpar_data%var(iLookPARAM%maxMassVegetation)/canopyDepth + & ! vegetation component
                                                              Cp_water*mvar_data%var(iLookMVAR%scalarCanopyLiq)%dat(1)/canopyDepth                              + & ! liquid water component
                                                              Cp_ice*mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1)/canopyDepth                                    ! ice component
  else
   mvar_data%var(iLookMVAR%scalarBulkVolHeatCapVeg)%dat(1)  = valueMissing
  endif
 
  ! compute the net change in snow in the vegetation canopy
  call canopySnow(&
                  ! input: model control
                  dt_sub,                                                      & ! intent(in): time step (seconds)
                  computeVegFlux,                                              & ! intent(in): flag to indicate if we are computing fluxes over vegetation (.false. means veg is buried with snow)
                  model_decisions(iLookDECISIONS%snowIncept)%iDecision,        & ! intent(in): choice of option to determine maximum snow interception capacity
                  ! input-output: state variables
                  mvar_data%var(iLookMVAR%scalarCanopyIce)%dat(1),             & ! intent(inout): storage of ice on the vegetation canopy (kg m-2)
                  ! input: diagnostic variables
                  exposedVAI,                                                  & ! intent(in): exposed vegetation area index (m2 m-2)
                  forc_data%var(iLookFORCE%airtemp),                           & ! intent(in): air temperature (K)
                  mvar_data%var(iLookMVAR%scalarSnowfall)%dat(1),              & ! intent(in): computed snowfall rate (kg m-2 s-1)
                  mvar_data%var(iLookMVAR%scalarNewSnowDensity)%dat(1),        & ! intent(in): density of new snow (kg m-3)
                  ! input: parameters
                  mpar_data%var(iLookPARAM%refInterceptCapSnow),               & ! intent(in): reference canopy interception capacity for snow per unit leaf area (kg m-2)
                  mpar_data%var(iLookPARAM%throughfallScaleSnow),              & ! intent(in): scaling factor for throughfall (snow) (-)
                  mpar_data%var(iLookPARAM%snowUnloadingCoeff),                & ! intent(in): time constant for unloading of snow from the forest canopy (s-1)
                  ! output: diagnostic variables
                  mvar_data%var(iLookMVAR%scalarCanopyIceMax)%dat(1),          & ! intent(out): maximum interception storage capacity for ice (kg m-2)
                  mvar_data%var(iLookMVAR%scalarThroughfallSnow)%dat(1),       & ! intent(out): snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
                  mvar_data%var(iLookMVAR%scalarCanopySnowUnloading)%dat(1),   & ! intent(out): unloading of snow from the vegetion canopy (kg m-2 s-1)
                  ! output: error control
                  err,cmessage                                                 ) ! intent(out): error control
  if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif

  ! initialize drainage and throughfall
  if(.not.computeVegFlux)then
   mvar_data%var(iLookMVAR%scalarThroughfallRain)%dat(1)   = mvar_data%var(iLookMVAR%scalarRainfall)%dat(1)
   mvar_data%var(iLookMVAR%scalarCanopyLiqDrainage)%dat(1) = 0._dp
  else
   mvar_data%var(iLookMVAR%scalarThroughfallRain)%dat(1)   = 0._dp
   mvar_data%var(iLookMVAR%scalarCanopyLiqDrainage)%dat(1) = 0._dp
  endif

  ! initialize maximum canopy liquid water (kg m-2)
  mvar_data%var(iLookMVAR%scalarCanopyLiqMax)%dat(1) = mpar_data%var(iLookPARAM%refInterceptCapRain)*exposedVAI

  ! compute the snow albedo
  ! NOTE: if canopy radiation is noah-mp, then albedo is computed within the Noah-MP radiation routine
  if(model_decisions(iLookDECISIONS%canopySrad)%iDecision /= noah_mp)then
   call snowAlbedo(&
                   ! input: model control
                   dt_sub,                                                 & ! intent(in): model time step
                   (nSnow > 0),                                            & ! intent(in): logical flag to denote if snow is present
                   model_decisions(iLookDECISIONS%alb_method)%iDecision,   & ! intent(in): index of method used for snow albedo
                   ! input: model variables
                   mvar_data%var(iLookMVAR%scalarSnowfall)%dat(1),         & ! intent(in): snowfall rate (kg m-2 s-1)
                   mvar_data%var(iLookMVAR%mLayerTemp)%dat(1),             & ! intent(in): surface temperature
                   mvar_data%var(iLookMVAR%scalarCosZenith)%dat(1),        & ! intent(in): cosine of the zenith angle (-)
                   ! input: model parameters
                   mpar_data%var(iLookPARAM%Frad_vis),                     & ! intent(in): fraction of radiation in visible part of spectrum (-)
                   mpar_data%var(iLookPARAM%Frad_direct),                  & ! intent(in): fraction direct solar radiation (-)
                   mpar_data%var(iLookPARAM%albedoMax),                    & ! intent(in): maximum snow albedo for a single spectral band (-)
                   mpar_data%var(iLookPARAM%albedoMinWinter),              & ! intent(in): minimum snow albedo during winter for a single spectral band (-)
                   mpar_data%var(iLookPARAM%albedoMinSpring),              & ! intent(in): minimum snow albedo during spring for a single spectral band (-)
                   mpar_data%var(iLookPARAM%albedoMaxVisible),             & ! intent(in): maximum snow albedo in the visible part of the spectrum (-)
                   mpar_data%var(iLookPARAM%albedoMinVisible),             & ! intent(in): minimum snow albedo in the visible part of the spectrum (-)
                   mpar_data%var(iLookPARAM%albedoMaxNearIR),              & ! intent(in): maximum snow albedo in the near infra-red part of the spectrum (-)
                   mpar_data%var(iLookPARAM%albedoMinNearIR),              & ! intent(in): minimum snow albedo in the near infra-red part of the spectrum (-)
                   mpar_data%var(iLookPARAM%albedoDecayRate),              & ! intent(in): albedo decay rate (s)
                   mpar_data%var(iLookPARAM%tempScalGrowth),               & ! intent(in): temperature scaling factor for grain growth (K-1) 
                   mpar_data%var(iLookPARAM%albedoSootLoad),               & ! intent(in): soot load factor (-)
                   mpar_data%var(iLookPARAM%albedoRefresh),                & ! intent(in): critical mass necessary for albedo refreshment (kg m-2)
                   mpar_data%var(iLookPARAM%snowfrz_scale),                & ! intent(in): scaling parameter for the freezing curve for snow (K-1) 
                   ! output
                   mvar_data%var(iLookMVAR%spectralSnowAlbedoDiffuse)%dat, & ! intent(inout): direct snow albedo in each spectral band (-)
                   mvar_data%var(iLookMVAR%spectralSnowAlbedoDirect)%dat,  & ! intent(inout): direct snow albedo in each spectral band (-)
                   mvar_data%var(iLookMVAR%scalarSnowAlbedo)%dat(1),       & ! intent(inout): snow albedo for the entire spectral band (-)
                   err,cmessage)                                             ! intent(out): error control
   if(err/=0)then; err=20; message=trim(message)//trim(cmessage); return; endif
  endif  ! (if NOT using the Noah-MP radiation routine)

  ! **
  ! NOTE: add new snowfall and layer divide/combine here, as vector length changes

  ! add new snowfall to the snow-soil system
  call newsnwfall(dt_sub,            & ! time step (seconds)
                  err,cmessage)        ! error control
  if(err/=0)then; err=30; message=trim(message)//trim(cmessage); return; endif

  ! divide snow layers if too thick
  call layerDivide(err,cmessage)        ! error control
  if(err/=0)then; err=55; message=trim(message)//trim(cmessage); return; endif

  ! merge snow layers if they are too thin
  call layerMerge(err,cmessage)        ! error control
  if(err/=0)then; err=65; message=trim(message)//trim(cmessage); return; endif

  ! **
  ! NOTE: need to re-assign pointers here as data structures have changed

  ! assign local pointers to the model index structures
  nSoil             => indx_data%var(iLookINDEX%nSoil)%dat(1)     ! number of soil layers
  nSnow             => indx_data%var(iLookINDEX%nSnow)%dat(1)     ! number of snow layers
  nLayers           => indx_data%var(iLookINDEX%nLayers)%dat(1)   ! total number of layers
  layerType         => indx_data%var(iLookINDEX%layerType)%dat    ! layer type (ix_soil or ix_snow)

  ! identify the number of snow and soil layers
  nSnow = count(layerType==ix_snow)
  nSoil = count(layerType==ix_soil)

  ! assign pointers to variables defining melt of the "snow without a layer"
  scalarSWE          => mvar_data%var(iLookMVAR%scalarSWE)%dat(1)         ! snow water equivalent (kg m-2)
  scalarSnowDepth    => mvar_data%var(iLookMVAR%scalarSnowDepth)%dat(1)   ! snow depth (m)
  scalarSfcMeltPond  => mvar_data%var(iLookMVAR%scalarSfcMeltPond)%dat(1) ! surface melt pond (kg m-2)
  mLayerVolHtCapBulk => mvar_data%var(iLookMVAR%mLayerVolHtCapBulk)%dat   ! volumetric heat capacity (J m-3 K-1)

  ! assign pointers to model state variables -- all layers
  mLayerTemp        => mvar_data%var(iLookMVAR%mLayerTemp)%dat                   ! temperature of each layer (K)
  mLayerDepth       => mvar_data%var(iLookMVAR%mLayerDepth)%dat                  ! depth of each layer (m)
  mLayerVolFracIce  => mvar_data%var(iLookMVAR%mLayerVolFracIce)%dat             ! volumetric fraction of ice in each layer (-)
  mLayerVolFracLiq  => mvar_data%var(iLookMVAR%mLayerVolFracLiq)%dat             ! volumetric fraction of liquid water in each layer (-)

  ! assign pointers to the model flux variables
  scalarThroughfallSnow     => mvar_data%var(iLookMVAR%scalarThroughfallSnow)%dat(1)     ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
  scalarThroughfallRain     => mvar_data%var(iLookMVAR%scalarThroughfallRain)%dat(1)     ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
  scalarCanopySnowUnloading => mvar_data%var(iLookMVAR%scalarCanopySnowUnloading)%dat(1) ! unloading of snow from the vegetion canopy (kg m-2 s-1)
  scalarCanopyLiqDrainage   => mvar_data%var(iLookMVAR%scalarCanopyLiqDrainage)%dat(1)   ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
  scalarCanopyMeltFreeze    => mvar_data%var(iLookMVAR%scalarCanopyMeltFreeze)%dat(1)    ! melt/freeze of water stored in the canopy (kg m-2 s-1)
  scalarSnowSublimation     => mvar_data%var(iLookMVAR%scalarSnowSublimation)%dat(1)     ! snow sublimation/frost - below canopy or non-vegetated (kg m-2 s-1)
  scalarGroundEvaporation   => mvar_data%var(iLookMVAR%scalarGroundEvaporation)%dat(1)   ! ground evaporation/condensation - below canopy or non-vegetated (kg m-2 s-1)
  scalarRainPlusMelt        => mvar_data%var(iLookMVAR%scalarRainPlusMelt)%dat(1)        ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  scalarSurfaceRunoff       => mvar_data%var(iLookMVAR%scalarSurfaceRunoff)%dat(1)       ! surface runoff (m s-1)
  scalarSoilInflux          => mvar_data%var(iLookMVAR%scalarSoilInflux)%dat(1)          ! influx of water at the top of the soil profile (m s-1)
  scalarSoilBaseflow        => mvar_data%var(iLookMVAR%scalarSoilBaseflow)%dat(1)        ! total baseflow from throughout the soil profile (m s-1)
  scalarSoilDrainage        => mvar_data%var(iLookMVAR%scalarSoilDrainage)%dat(1)        ! drainage from the bottom of the soil profile (m s-1)
  scalarAquiferRecharge     => mvar_data%var(iLookMVAR%scalarAquiferRecharge)%dat(1)     ! recharge to the aquifer (m s-1)
  scalarAquiferBaseflow     => mvar_data%var(iLookMVAR%scalarAquiferBaseflow)%dat(1)     ! baseflow from the aquifer (m s-1)
  scalarAquiferTranspire    => mvar_data%var(iLookMVAR%scalarAquiferTranspire)%dat(1)    ! transpiration from the aquifer (m s-1)
  mLayerColumnOutflow       => mvar_data%var(iLookMVAR%mLayerColumnOutflow)%dat          ! total outflow from each layer in a given soil column (m3 s-1) 

  ! allocate temporary array
  allocate(arrTemp(nLayers),stat=err)
  if(err/=0)then; err=20; message='problem allocating space for temporary array'; return; endif

  ! save the volumetric fraction of ice
  arrTemp = mLayerDepth

  ! use Picard iteration to solve model equations
  do
   ! get the new solution
   call picardSolv(dt_sub,maxiter,(nsub==1),computeVegFlux,&  ! input
                   niter,err,cmessage)                        ! output
   if(err > 0)then; message=trim(message)//trim(cmessage); return; endif
   !if(err<0)then; print*, trim(message)//trim(cmessage); print*, 'dt_sub, minstep = ', dt_sub, minstep; pause; endif 
   ! exit do loop if all is a-ok
   if(err==0) exit
   ! if not ok, reduce time step and try again
   dt_sub = dt_sub*0.1_dp
   print*, dt_sub, minstep, trim(message)//trim(cmessage)
   !if(dt_sub < 10._dp)then
   ! pause ' dt_sub < 10'
   ! err=20; return
   !endif
   ! check that the step size is still appropriate -- if not, use non-iterative solution
   if(dt_sub < minstep)then
    if(err/=0)then; message=trim(message)//'dt_sub is below the minimum time step'; return; endif
    dt_sub  = minstep
    ! just iterate once
    call picardSolv(dt_sub,1,(nsub==1),computeVegFlux,     &  ! input
                    niter,err,cmessage)                       ! output
    if(err/=0)then; message=trim(message)//trim(cmessage); return; endif
    exit ! exit do loop if all is a-ok
   endif
  end do 

  ! compute melt for the case of "snow without a layer"
  if(nSnow==0 .and. scalarSWE > 0._dp)then
   ! only melt if temperature of the top soil layer is greater than Tfreeze
   if(mLayerTemp(1) > Tfreeze)then
    ! compute the energy required to melt all the snow (J m-2)
    nrgRequired     = scalarSWE*LH_fus
    ! compute the energy available to melt the snow (J m-2)
    nrgAvailable    = mLayerVolHtCapBulk(1)*(mLayerTemp(1) - Tfreeze)*mLayerDepth(1)
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
    mLayerTemp(1)= mLayerTemp(1) - (scalarSfcMeltPond/mLayerDepth(1))/mLayerVolHtCapBulk(1)
   else  ! melt is zero if the temperature of the top soil layer is less than Tfreeze
    scalarSfcMeltPond = 0._dp  ! kg m-2
   endif ! (if the temperature of the top soil layer is greater than Tfreeze)
  else  ! melt is zero if the "snow without a layer" does not exist
   scalarSfcMeltPond = 0._dp  ! kg m-2
  endif ! (if the "snow without a layer" exists)

  ! define weight applied to each sub-step
  dt_wght = dt_sub/dt

  ! increment timestep-average fluxes
  averageThroughfallSnow     = averageThroughfallSnow     + scalarThroughfallSnow     *dt_wght ! snow that reaches the ground without ever touching the canopy (kg m-2 s-1)
  averageThroughfallRain     = averageThroughfallRain     + scalarThroughfallRain     *dt_wght ! rain that reaches the ground without ever touching the canopy (kg m-2 s-1)
  averageCanopySnowUnloading = averageCanopySnowUnloading + scalarCanopySnowUnloading *dt_wght ! unloading of snow from the vegetion canopy (kg m-2 s-1)
  averageCanopyLiqDrainage   = averageCanopyLiqDrainage   + scalarCanopyLiqDrainage   *dt_wght ! drainage of liquid water from the vegetation canopy (kg m-2 s-1)
  averageCanopyMeltFreeze    = averageCanopyMeltFreeze    + scalarCanopyMeltFreeze    *dt_wght ! melt/freeze of water stored in the canopy (kg m-2 s-1)
  averageSnowSublimation     = averageSnowSublimation     + scalarSnowSublimation     *dt_wght ! snow sublimation/frost - below canopy or non-vegetated (kg m-2 s-1)
  averageGroundEvaporation   = averageGroundEvaporation   + scalarGroundEvaporation   *dt_wght ! ground evaporation/condensation - below canopy or non-vegetated (kg m-2 s-1)
  averageRainPlusMelt        = averageRainPlusMelt        + scalarRainPlusMelt        *dt_wght ! rain plus melt, as input to soil before calculating surface runoff (m s-1)
  averageSurfaceRunoff       = averageSurfaceRunoff       + scalarSurfaceRunoff       *dt_wght ! surface runoff (m s-1)
  averageSoilInflux          = averageSoilInflux          + scalarSoilInflux          *dt_wght ! influx of water at the top of the soil profile (m s-1)
  averageSoilBaseflow        = averageSoilBaseflow        + scalarSoilBaseflow        *dt_wght ! total baseflow from throughout the soil profile (m s-1)
  averageSoilDrainage        = averageSoilDrainage        + scalarSoilDrainage        *dt_wght ! drainage from the bottom of the soil profile (m s-1)
  averageAquiferRecharge     = averageAquiferRecharge     + scalarAquiferRecharge     *dt_wght ! recharge to the aquifer (m s-1)
  averageAquiferBaseflow     = averageAquiferBaseflow     + scalarAquiferBaseflow     *dt_wght ! baseflow from the aquifer (m s-1)
  averageAquiferTranspire    = averageAquiferTranspire    + scalarAquiferTranspire    *dt_wght ! transpiration from the aquifer (m s-1)
  averageColumnOutflow       = averageColumnOutflow       + mLayerColumnOutflow       *dt_wght ! outflow from each soil layer in a given soil column (m3 s-1)

  !write(*,'(a,1x,f9.1,1x,10(e20.10,1x))') 'dt_sub, averageSoilBaseflow, scalarSoilBaseflow', &
  !                                         dt_sub, averageSoilBaseflow, scalarSoilBaseflow
  !write(*,'(a,10(e20.10,1x))') 'scalarRainPlusMelt, scalarSurfaceRunoff, scalarSoilInflux = ', &
  !                              scalarRainPlusMelt, scalarSurfaceRunoff, scalarSoilInflux

  ! check that snow depth is decreasing (can only increase in the top layer)
  if(nSnow>1)then
   do iSnow=2,nSnow
    if(mLayerDepth(iSnow) > arrTemp(iSnow)+1.e-8_dp)then
     write(*,'(a,1x,100(f20.10))') 'depth1 = ', arrTemp(1:nSnow)
     write(*,'(a,1x,100(f20.10))') 'depth2 = ', mLayerDepth(1:nSnow)
     write(*,'(a,1x,100(f20.10))') 'diff   = ', mLayerDepth(1:nSnow) - arrTemp(1:nSnow)
     stop 'depth is increasing '
    endif
   end do ! looping thru snow layers
  endif

  ! increment the time step increment
  dt_done = dt_done + dt_sub
  !print*, '***** ', dt_done, dt_sub, niter

  ! modify the length of the time step
  if(niter<n_inc) dt_sub = min(dt_sub*F_inc,maxstep)
  if(niter>n_dec) dt_sub =     dt_sub*F_dec
  if(dt_sub < minstep)then; message=trim(message)//'dt_sub is below the minimum time step'; return; endif

  ! save the time step to initialize the subsequent step
  if(dt_done<dt .or. nsub==1) dt_init = dt_sub
  if(dt_init < 0.00001_dp .and. nsub > 10000) then
   write(message,'(a,f13.10,a,f9.2,a,i0,a)')trim(message)//"dt < 0.00001 and nsub > 10000 [dt=",dt_init,"; dt_done=",&
         dt_done,"; nsub=",nsub,"]"
   err=20; return
  endif

  ! exit do-loop if finished
  if(dt_done>=dt)exit

  ! make sure that we don't exceed the step
  dt_sub = min(dt-dt_done, dt_sub)
  
  ! deallocate temporary array
  deallocate(arrTemp,stat=err)
  if(err/=0)then; err=20; message='problem deallocating space for temporary array'; return; endif

 end do  ! (sub-step loop)
 !pause 'completed time step'

 !print*, 'mvar_data%var(iLookMVAR%averageCanopyLiqDrainage)%dat(1) = ', mvar_data%var(iLookMVAR%averageCanopyLiqDrainage)%dat(1)

 ! save the surface temperature (just to make things easier to visualize)
 mvar_data%var(iLookMVAR%scalarSurfaceTemp)%dat(1) = mvar_data%var(iLookMVAR%mLayerTemp)%dat(1)

 iLayer = nSnow+1
 !print*, 'nsub, mLayerTemp(iLayer), mLayerVolFracIce(iLayer) = ', nsub, mLayerTemp(iLayer), mLayerVolFracIce(iLayer)
 print*, 'nsub = ', nsub
 if(nsub>200)then
  message=trim(message)//'number of sub-steps > 200'
  err=20; return
 endif

 !write(*,'(a)') '==========================================================================================================================='
 !write(*,'(a)') '==========================================================================================================================='
 !write(*,'(a)') '==========================================================================================================================='
 !write(*,'(a)') '==========================================================================================================================='
 !write(*,'(a)') '==========================================================================================================================='

 
 !if(mLayerVolFracIce(iLayer) > 0.5_dp) pause 'ice content in top soil layer is huge...'
 
 end subroutine coupled_em

end module coupled_em_module
